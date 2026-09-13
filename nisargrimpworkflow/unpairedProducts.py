#!/usr/bin/env python3
"""Registry of NISAR L2 products delivered without their mate.

A GrIMP frame needs both a RUNW (unwrapped phase) and its ROFF (offsets).  ASF
does not always deliver the pair together -- during an in-progress cycle the
RUNW routinely lands first -- and until the mate arrives the frame cannot be
processed.  Previously that left an empty `orbit_frame/H5/` directory behind
with nothing recording why, so the gap was invisible until someone noticed the
directory.

This keeps a YAML registry beside the archive naming every product whose mate
is absent, so the gap is explicit and can be re-checked cheaply on each run.

The registry is APPEND-ONLY in the sense that matters: entries already present
are never dropped or rewritten from scratch.  A run may add new entries, and
may mark an existing one resolved, but it never discards the history -- a
resolved entry keeps its record so the delivery lag stays auditable.

    <archiveDir>/unpairedProducts.yaml

    generated: ...                 # last time a run touched the file
    entries:                       # a LIST -- granule names exceed YAML's
      - name: NISAR_L1_PR_RUNW_029_001_A_139_...h5    # 128-char simple-key
        present: RUNW              # limit, which would force the ugly
        missing: ROFF              # "? key : value" explicit form
        cycle: '029'
        track: 1
        direction: A
        frame: 139
        productDate: 2026-09-08T22:06:11    # mtime of the present file
        firstSeen: 2026-09-11T14:02:00      # when we first noticed the gap
        lastChecked: 2026-09-11T14:02:00
        checkCount: 1
        resolved: false
        resolvedAt: null
        waitDays: 0.0                       # productDate -> resolvedAt/now
"""
import datetime
import os
import re
import tempfile

import yaml

REGISTRY_NAME = 'unpairedProducts.yaml'
# ROFF is the mate GrIMP actually needs; RIFG is linked when present but its
# absence never blocks a frame, so it is not tracked here.
TRACKED_PAIR = {'RUNW': 'ROFF', 'ROFF': 'RUNW'}
_NAME = re.compile(r'NISAR_L1_PR_(?P<prod>[A-Z]+)_(?P<cycle>\d+)_(?P<track>\d+)'
                   r'_(?P<dir>[AD])_(?P<frame>\d+)_')


def _now():
    return datetime.datetime.now().replace(microsecond=0).isoformat()


def registryPath(archiveDir):
    return os.path.join(archiveDir, REGISTRY_NAME)


def parseGranule(name):
    """cycle/track/direction/frame from a granule filename, or {} if unparsable."""
    m = _NAME.match(os.path.basename(name))
    if not m:
        return {}
    return dict(product=m.group('prod'), cycle=m.group('cycle'),
                track=int(m.group('track')), direction=m.group('dir'),
                frame=int(m.group('frame')))


def load(archiveDir):
    """Read the registry, or an empty one.  Never raises on a bad file."""
    path = registryPath(archiveDir)
    if not os.path.exists(path):
        return {'entries': {}}
    try:
        with open(path) as fp:
            data = yaml.safe_load(fp) or {}
    except Exception:                                          # noqa: BLE001
        # A corrupt registry must not stop a filing run; it is an audit aid,
        # not an input.  Start fresh in memory but keep the file on disk.
        return {'entries': {}, 'unreadable': True}
    if not isinstance(data, dict):
        return {'entries': {}, 'unreadable': True}
    raw = data.get('entries')
    if isinstance(raw, list):            # current on-disk form
        entries = {e['name']: {k: v for k, v in e.items() if k != 'name'}
                   for e in raw if isinstance(e, dict) and 'name' in e}
    elif isinstance(raw, dict):          # tolerate the earlier keyed form
        entries = raw
    else:
        return {'entries': {}, 'unreadable': True}
    data['entries'] = entries
    return data


def save(archiveDir, registry):
    """Write atomically so a crash mid-write cannot truncate the registry."""
    if registry.get('unreadable'):
        # Refuse to overwrite a file we could not parse -- it may hold history
        # worth recovering by hand.
        return None
    path = registryPath(archiveDir)
    entries = registry.get('entries', {})
    out = {'generated': _now(),
           'entries': [dict(name=n, **entries[n]) for n in sorted(entries)]}
    d = os.path.dirname(path) or '.'
    fd, tmp = tempfile.mkstemp(dir=d, prefix='.unpaired', suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as fp:
            yaml.safe_dump(out, fp, default_flow_style=False, sort_keys=True)
        os.replace(tmp, path)
    except Exception:                                          # noqa: BLE001
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return path


def _productDate(path):
    try:
        return datetime.datetime.fromtimestamp(
            os.path.getmtime(path)).replace(microsecond=0).isoformat()
    except OSError:
        return None


def record(registry, presentPath, missingType):
    """Note that `presentPath` has no `missingType` mate.  Idempotent."""
    name = os.path.basename(presentPath)
    entries = registry.setdefault('entries', {})
    e = entries.get(name)
    if e is None:
        info = parseGranule(name)
        e = dict(present=info.get('product', '?'), missing=missingType,
                 cycle=info.get('cycle'), track=info.get('track'),
                 direction=info.get('direction'), frame=info.get('frame'),
                 path=presentPath, productDate=_productDate(presentPath),
                 firstSeen=_now(), checkCount=0,
                 resolved=False, resolvedAt=None, waitDays=0.0)
        entries[name] = e
    elif e.get('resolved'):
        # It was resolved before and the mate has gone away again -- reopen it
        # rather than silently leaving a stale "resolved".
        e['resolved'] = False
        e['resolvedAt'] = None
        e['reopenedAt'] = _now()
    e['lastChecked'] = _now()
    e['checkCount'] = int(e.get('checkCount', 0)) + 1
    e['waitDays'] = _waitDays(e, None)
    return e


def _waitDays(entry, endIso):
    start = entry.get('productDate') or entry.get('firstSeen')
    if not start:
        return 0.0
    try:
        t0 = datetime.datetime.fromisoformat(start)
        t1 = (datetime.datetime.fromisoformat(endIso) if endIso
              else datetime.datetime.now())
        return round((t1 - t0).total_seconds() / 86400.0, 2)
    except ValueError:
        return 0.0


def matePath(archiveDir, name, missingType):
    """Where the mate of `name` would live."""
    info = parseGranule(name)
    if not info:
        return None
    mate = os.path.basename(name).replace(f'_{info["product"]}_',
                                          f'_{missingType}_', 1)
    return os.path.join(archiveDir, missingType, mate)


def checkResolved(archiveDir, registry, verbose=False):
    """Re-check every open entry; mark the ones whose mate has arrived.

    Returns (resolvedNames, stillOpenCount).
    """
    resolved = []
    open_ = 0
    for name, e in registry.get('entries', {}).items():
        if e.get('resolved'):
            continue
        mate = matePath(archiveDir, name, e.get('missing', ''))
        e['lastChecked'] = _now()
        e['checkCount'] = int(e.get('checkCount', 0)) + 1
        if mate and os.path.exists(mate):
            e['resolved'] = True
            e['resolvedAt'] = _now()
            e['mateDate'] = _productDate(mate)
            e['waitDays'] = _waitDays(e, e['resolvedAt'])
            resolved.append(name)
            if verbose:
                print(f'  mate arrived after {e["waitDays"]} d: {name}')
        else:
            open_ += 1
            e['waitDays'] = _waitDays(e, None)
    return resolved, open_


def summarise(registry):
    entries = registry.get('entries', {})
    open_ = [e for e in entries.values() if not e.get('resolved')]
    if not entries:
        return 'unpaired registry: empty'
    oldest = max((e.get('waitDays', 0) for e in open_), default=0)
    byMissing = {}
    for e in open_:
        byMissing[e.get('missing', '?')] = byMissing.get(e.get('missing', '?'), 0) + 1
    bits = ', '.join(f'{v} awaiting {k}' for k, v in sorted(byMissing.items()))
    return (f'unpaired registry: {len(open_)} open'
            + (f' ({bits}; oldest {oldest:.1f} d)' if open_ else '')
            + f', {len(entries) - len(open_)} resolved')
