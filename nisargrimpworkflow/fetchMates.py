#!/usr/bin/env python3
"""Download the missing mates found by checkUnpaired and file them.

Downloads into `<archiveDir>/mateFetch-<YYYY-MM-DD>/`, which is deliberately
NOT a `<MM-DD-YYYY>` staging folder: `autoupdate.stagingFolders()` only matches
that date format, so the nightly's release logic ignores this directory, while
`searchASF`'s dedup glob (`archiveDir/*/*`) still covers it and will not
re-request granules that are mid-flight here.

Only granules that verify are moved into `RUNW/` or `ROFF/`:
an `.aria2` control file means the download was interrupted, and a zero-length
file means nothing arrived.  A partial left in the product directory would be
worse than a missing one -- it would be released as a corrupt product and
would suppress its own re-download via the dedup glob.
"""
import argparse
import datetime
import os
import shutil
import subprocess
import sys

PRODUCTS = ('RUNW', 'ROFF', 'RIFG')


def verify(path):
    """(ok, why) for a downloaded granule."""
    if not os.path.exists(path):
        return False, 'not downloaded'
    if os.path.exists(path + '.aria2'):
        return False, 'incomplete (.aria2 control file present)'
    if os.path.getsize(path) == 0:
        return False, 'zero length'
    return True, ''


def productDir(archiveDir, name):
    for p in PRODUCTS:
        if f'_{p}_' in name:
            return os.path.join(archiveDir, p)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive', required=True)
    ap.add_argument('--urls', required=True, help='URL list from checkUnpaired --outList')
    ap.add_argument('--attempts', type=int, default=3)
    ap.add_argument('--dryRun', action='store_true')
    args = ap.parse_args()

    archive = os.path.abspath(args.archive)
    urls = [u.strip() for u in open(args.urls) if u.strip()]
    names = [os.path.basename(u) for u in urls]
    stamp = datetime.date.today().strftime('%Y-%m-%d')
    fetchDir = os.path.join(archive, f'mateFetch-{stamp}')
    print(f'archive : {archive}')
    print(f'urls    : {len(urls)}  ->  {fetchDir}')
    if args.dryRun:
        for n in names[:10]:
            print('   ', n)
        return
    os.makedirs(fetchDir, exist_ok=True)

    missing = names
    for attempt in range(1, args.attempts + 1):
        if attempt > 1:
            # ariaDownload skips an existing file, so a partial would never be
            # retried; clear it and its control file first.
            for n in missing:
                for stale in (os.path.join(fetchDir, n),
                              os.path.join(fetchDir, n + '.aria2')):
                    if os.path.exists(stale):
                        os.remove(stale)
        print(f'ariaDownload attempt {attempt}/{args.attempts} '
              f'({len(missing)} remaining)', flush=True)
        r = subprocess.run(['ariaDownload', os.path.abspath(args.urls)],
                           cwd=fetchDir)
        if r.returncode != 0:
            print(f'  ariaDownload rc={r.returncode}')
        missing = [n for n in names if not verify(os.path.join(fetchDir, n))[0]]
        if not missing:
            break

    filed = 0
    for n in names:
        src = os.path.join(fetchDir, n)
        ok, why = verify(src)
        if not ok:
            print(f'  NOT FILED ({why}): {n}')
            continue
        dest = productDir(archive, n)
        if dest is None:
            print(f'  NOT FILED (no product subdir): {n}')
            continue
        os.makedirs(dest, exist_ok=True)
        target = os.path.join(dest, n)
        if os.path.exists(target):
            os.remove(src)
        else:
            shutil.move(src, target)
        filed += 1
    print(f'filed {filed}/{len(names)} into RUNW/ROFF')
    leftovers = os.listdir(fetchDir) if os.path.isdir(fetchDir) else []
    if not leftovers:
        os.rmdir(fetchDir)
        print(f'removed empty {fetchDir}')
    else:
        print(f'{len(leftovers)} file(s) left in {fetchDir} (incomplete)')
    return 0 if filed == len(names) else 1


if __name__ == '__main__':
    sys.exit(main())
