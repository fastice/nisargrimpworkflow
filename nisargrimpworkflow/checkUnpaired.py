#!/usr/bin/env python3
"""Check-only report on NISAR products whose mate was never delivered.

Answers two questions without changing anything:

  1. ON DISK   -- has the missing mate turned up in the archive since the
                  product was filed?  (cheap, no network)
  2. ONLINE    -- is the missing mate available at ASF now, even though it is
                  not on disk?  (needs --online; runs searchASF, which dedups
                  against the archive, then matches the returned granules
                  against the open registry entries)

Read-only by default: the registry is not written unless --update is given.
If no registry exists yet the unpaired set is derived from the archive itself,
so this works against archives that predate the registry.

    checkUnpaired --archive /Volumes/insar4/ian/Data/NISAR/Greenland
    checkUnpaired --archive ... --online --region greenland
"""
import argparse
import collections
import datetime
import os
import re
import subprocess
import sys
import tempfile

from nisargrimpworkflow import unpairedProducts as up

PAIR = {'RUNW': 'ROFF', 'ROFF': 'RUNW'}
# reference acquisition date sits in the 9th underscore field of the name
_DATE = re.compile(r'_(\d{8})T\d{6}_')


def unpairedFromArchive(archiveDir):
    """(presentPath, missingType) for every product in the archive with no mate."""
    have = {}
    for prod in PAIR:
        d = os.path.join(archiveDir, prod)
        have[prod] = ({n for n in os.listdir(d) if n.endswith('.h5')}
                      if os.path.isdir(d) else set())
    out = []
    for prod, mate in PAIR.items():
        for n in have[prod]:
            if n.replace(f'_{prod}_', f'_{mate}_') not in have[mate]:
                out.append((os.path.join(archiveDir, prod, n), mate))
    return out, {k: len(v) for k, v in have.items()}


def acquisitionDate(name):
    m = _DATE.search(name)
    return m.group(1) if m else None


def onlineCheck(archiveDir, openEntries, region, bandwidths, verbose=False,
                firstDate=None, lastDate=None):
    """Ask ASF whether the missing mates exist online.

    Runs one searchASF per product type over the date span of the open
    entries.  searchASF dedups against the archive, so anything it returns is
    online-but-not-on-disk; we then keep only the names we are waiting for.
    """
    wanted = collections.defaultdict(set)       # missingType -> {granule names}
    dates = []
    for name, e in openEntries:
        mate = os.path.basename(name).replace(f'_{e["present"]}_',
                                              f'_{e["missing"]}_', 1)
        wanted[e['missing']].add(mate)
        d = acquisitionDate(name)
        if d:
            dates.append(d)
    if not dates:
        return {}, 'no dated entries to search'
    # The window must span the WHOLE archive, not just the dates of the open
    # entries.  The project reprocesses older tracks on request with the latest
    # processor, and those products appear at ASF unpredictably, long after the
    # acquisition -- a window derived from what we already hold would keep
    # missing them.  Default: earliest acquisition anywhere in the archive
    # (less a month of margin) through tomorrow.
    if firstDate and lastDate:
        first, last = firstDate, lastDate
    else:
        onDisk = []
        for prod in PAIR:
            d = os.path.join(archiveDir, prod)
            if os.path.isdir(d):
                onDisk += [acquisitionDate(n) for n in os.listdir(d)
                           if n.endswith('.h5')]
        onDisk = [d for d in onDisk if d] or dates
        firstDt = (datetime.datetime.strptime(min(onDisk), '%Y%m%d')
                   - datetime.timedelta(days=31))
        first = firstDate or firstDt.strftime('%Y-%m-%d')
        last = lastDate or (datetime.date.today()
                            + datetime.timedelta(days=1)).strftime('%Y-%m-%d')
    print(f'  search window: {first} .. {last}')

    found = {}
    for missingType, names in sorted(wanted.items()):
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, f'urls.{missingType}')
            cmd = ['searchASF', '--sensor', 'NISAR',
                   '--archiveDir', os.path.join(archiveDir, '*', '*'),
                   '--products', missingType,
                   '--minVersion', '5012',
                   f'--{region}', first, last, out]
            if bandwidths:
                cmd += ['--bandwidth'] + list(bandwidths)
            if verbose:
                print('  ' + ' '.join(cmd))
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=3600)
            except Exception as exc:                           # noqa: BLE001
                return found, f'searchASF failed: {type(exc).__name__}: {exc}'
            if r.returncode != 0:
                tail = (r.stderr or r.stdout).strip().splitlines()[-3:]
                return found, f'searchASF rc={r.returncode}: {" | ".join(tail)}'
            # keep the URL, not just the basename -- the worklist has to be
            # directly usable by ariaDownload
            online = {}
            for path in (out, out + '.updated'):
                if os.path.exists(path):
                    for u in open(path):
                        u = u.strip()
                        if u:
                            online[os.path.basename(u)] = u
            hits = {n: online[n] for n in names if n in online}
            found[missingType] = dict(searched=len(names), onlineTotal=len(online),
                                      hits=dict(sorted(hits.items())))
    return found, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--archive', required=True, help='archive dir holding RUNW/ ROFF/')
    ap.add_argument('--online', action='store_true',
                    help='also ask ASF whether the missing mates exist there')
    ap.add_argument('--region', default='greenland',
                    choices=['greenland', 'antarctica'])
    ap.add_argument('--bandwidth', nargs='*', default=None,
                    help='searchASF --bandwidth tokens (e.g. 4000 or 7700)')
    ap.add_argument('--update', action='store_true',
                    help='write the registry back (default: read-only)')
    ap.add_argument('--limit', type=int, default=15, help='examples to print')
    ap.add_argument('--outList', default=None,
                    help='write every online-but-not-downloaded mate to this '
                         'file, one granule per line (a download worklist)')
    ap.add_argument('--firstDate', default=None,
                    help='search window start (YYYY-MM-DD); default is the '
                         'earliest acquisition in the archive less one month')
    ap.add_argument('--lastDate', default=None,
                    help='search window end (YYYY-MM-DD); default tomorrow')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    archive = os.path.abspath(args.archive)
    print(f'archive: {archive}')
    registry = up.load(archive)
    seeded = False
    if not registry.get('entries'):
        pairs, counts = unpairedFromArchive(archive)
        print(f'  no registry yet -- deriving from archive '
              f'(RUNW {counts.get("RUNW", 0)}, ROFF {counts.get("ROFF", 0)})')
        for p, m in pairs:
            up.record(registry, p, m)
        seeded = True

    # --- 1. on-disk recheck -------------------------------------------------
    resolved, stillOpen = up.checkResolved(archive, registry,
                                           verbose=args.verbose)
    entries = registry.get('entries', {})
    openEntries = [(n, e) for n, e in entries.items() if not e.get('resolved')]
    print(f'\nON DISK')
    print(f'  registry entries      : {len(entries)}'
          + ('  (seeded from archive this run)' if seeded else ''))
    print(f'  mate arrived since    : {len(resolved)}')
    print(f'  still missing its mate: {stillOpen}')
    byMissing = collections.Counter(e['missing'] for _, e in openEntries)
    for k, v in sorted(byMissing.items()):
        print(f'      awaiting {k}: {v}')
    if openEntries:
        waits = sorted((e.get('waitDays', 0) for _, e in openEntries), reverse=True)
        print(f'  wait (days)           : max {waits[0]:.1f}, '
              f'median {waits[len(waits) // 2]:.1f}')
    for n in resolved[:args.limit]:
        print(f'    RESOLVED: {n}')

    # --- 2. online check ----------------------------------------------------
    if args.online:
        print(f'\nONLINE (ASF)')
        found, err = onlineCheck(archive, openEntries, args.region,
                                 args.bandwidth, args.verbose,
                                 args.firstDate, args.lastDate)
        if err:
            print(f'  could not complete: {err}')
        for missingType, info in sorted(found.items()):
            print(f'  {missingType}: {len(info["hits"])} of {info["searched"]} '
                  f'missing mates are available at ASF '
                  f'({info["onlineTotal"]} granule(s) online but not on disk)')
            for n in list(info['hits'])[:args.limit]:
                print(f'    AVAILABLE: {n}')
            if len(info['hits']) > args.limit:
                print(f'    ... and {len(info["hits"]) - args.limit} more')
        if args.outList:
            # hits maps granule name -> download URL; the worklist is the URLs
            urls = sorted(u for i in found.values() for u in i['hits'].values())
            with open(args.outList, 'w') as fp:
                fp.write('\n'.join(urls) + ('\n' if urls else ''))
            print(f'  download worklist ({len(urls)} URL(s)): {args.outList}')

    # --- 3. write back only if asked ---------------------------------------
    if args.update:
        path = up.save(archive, registry)
        print(f'\nregistry updated: {path}')
    else:
        print(f'\n(check only -- registry not written; pass --update to save)')
    print(up.summarise(registry))


if __name__ == '__main__':
    main()
