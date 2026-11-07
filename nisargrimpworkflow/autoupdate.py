#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated NISAR search/download/update driver, configured by autoupdate.yaml.

Intended to run as a cron job. All paths and options come from autoupdate.yaml,
which normally sits beside the project's project.yaml.

If autoupdate.yaml sets `notifyEmail: you@example.com`, the run mails that
address the per-run summary whenever it hits errors (a failed step, a failed
SetupNISAR orbit, or a crash). A clean run stays silent, and with no
notifyEmail key no mail is ever sent (opt-in, best effort via the local MTA).
"""
import argparse
import datetime
import glob
import json
import os
import re
import shutil
import socket
import subprocess
import tempfile
import threading
import traceback

import yaml
import utilities as u

from nisargrimpworkflow.FileNISARProducts import (parseFileName,
                                                  orbitFromCycleTrack)

# Only two region choices for now; each maps to a searchASF spatial flag.
regionFlags = {
    'antarctica': '--antarctica',
    'greenland': '--greenland',
}

# Dated staging-folder name format (matches the Download.MM-DD-YYYY basename
# already used for the search URL lists).
STAGING_DATE_FMT = '%m-%d-%Y'


def parseArgs():
    '''
    Handle command line args.
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mAutomated NISAR search/download/update driver, '
        'configured by autoupdate.yaml\033[0m\n\n',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('config', type=str, nargs='?', default='autoupdate.yaml',
                        help='Path to autoupdate.yaml (default: ./autoupdate.yaml)')
    parser.add_argument('--noDownload', action='store_true',
                        help='Skip the ASF search + download and start from the '
                        'dated staging folders already present under archiveDir. '
                        'Useful for re-running the processing chain without '
                        're-downloading, and for testing.')
    return parser.parse_args()


def loadConfig(configFile):
    '''
    Read the autoupdate.yaml config into a dict.
    '''
    with open(configFile) as fp:
        config = yaml.safe_load(fp) or {}
    return config


def regionFlag(region):
    '''
    Map a region name (Antarctica/Greenland) to its searchASF flag.
    '''
    key = str(region).strip().lower()
    if key not in regionFlags:
        raise ValueError(
            f"region must be one of {list(regionFlags)}, got '{region}'")
    return regionFlags[key]


def searchGranules(config, today):
    '''
    Search ASF for new NISAR granules over the configured region/time window,
    writing the URL lists and coverage GeoPackage into archiveDir/searchResults.
    `today` is fixed once at program start (main()) so the search window, URL-list
    name, and download-folder name all agree even if the run crosses midnight.
    '''
    archiveDir = config['archiveDir']
    lookBack = int(config['lookBack'])
    products = granuleProducts(config)
    # Time window: today - lookBack days through today.
    firstDate = (today - datetime.timedelta(days=lookBack)).strftime('%Y-%m-%d')
    lastDate = today.strftime('%Y-%m-%d')
    # Output basename Download.MM-DD-YYYY. URL lists go under
    # archiveDir/searchResults; the coverage GeoPackage under searchResults/gpkg.
    d = today.strftime('%m-%d-%Y')
    searchDir = os.path.join(archiveDir, 'searchResults')
    gpkgDir = os.path.join(searchDir, 'gpkg')
    os.makedirs(gpkgDir, exist_ok=True)
    out = os.path.join(searchDir, f'Download.{d}')
    gpkg = os.path.join(gpkgDir, f'Download.{d}.gpkg')
    # Absolute paths throughout so the call is independent of the working
    # directory (robust under cron); dedup glob matches per-product subdirs.
    archiveGlob = os.path.join(archiveDir, '*', '*')
    bandwidths = searchBandwidths(config)
    #
    command = (['searchASF', '--sensor', 'NISAR', '--archiveDir', archiveGlob,
                '--bandwidth'] + bandwidths
               + ['--minVersion', '5012', regionFlag(config['region']),
                  '--gpkg', gpkg, firstDate, lastDate, out, '--products']
               + products)
    print(' '.join(command))
    subprocess.run(command, check=True)
    return out


def granuleProducts(config):
    '''
    Return the product list from config as a list (accepts a space-separated
    string or a YAML list).
    '''
    products = config['products']
    if isinstance(products, str):
        products = products.split()
    return products


def searchBandwidths(config):
    '''
    Return the searchASF --bandwidth tokens from config. Accepts a single
    value, a space-separated string, or a YAML list. The GrIMP convention
    refers to the 77 MHz mode as "80", so 80 is mapped to 77 for the search.
    '''
    bandwidths = config['bandwidths']
    if isinstance(bandwidths, str):
        bandwidths = bandwidths.split()
    elif not isinstance(bandwidths, (list, tuple)):
        bandwidths = [bandwidths]
    return ['77' if str(bw).strip() == '80' else str(bw).strip()
            for bw in bandwidths]


def masterInputPaths(config, projectDir):
    '''
    Return the makeMaster --inputPaths list from config, or None when the key
    is unset (makeMaster then falls back to --projectDir). Accepts a single
    value, a space-separated string, or a YAML list. Relative entries are
    anchored to projectDir so the config reads the same from any cwd.
    '''
    paths = config.get('masterInputPaths')
    if paths is None:
        return None
    if isinstance(paths, str):
        paths = paths.split()
    elif not isinstance(paths, (list, tuple)):
        paths = [paths]
    return [os.path.normpath(path if os.path.isabs(path)
                             else os.path.join(projectDir, path))
            for path in (str(p).strip() for p in paths)]


def masterOutputPath(config, projectDir):
    '''
    Return the absolute directory makeMaster writes inputFile into. Defaults to
    projectDir/Release/masterInput (makeMaster's own default); masterOutputPath
    in autoupdate.yaml overrides it, anchored to projectDir when relative.
    '''
    path = config.get('masterOutputPath') or 'Release/masterInput'
    if not os.path.isabs(path):
        path = os.path.join(projectDir, path)
    return os.path.normpath(path)


def readUrlList(urlListFile):
    '''
    Return the granule filenames (URL basenames) listed in a searchASF URL
    list file.
    '''
    names = []
    with open(urlListFile) as fp:
        for line in fp:
            url = line.strip()
            if url:
                names.append(url.split('/')[-1])
    return names


def incompleteGranules(downloadDir, names):
    '''
    Return the subset of names not fully downloaded in downloadDir. A granule
    is complete once its file exists and aria2c's .aria2 control file is gone.
    '''
    missing = []
    for name in names:
        path = os.path.join(downloadDir, name)
        if not os.path.exists(path) or os.path.exists(f'{path}.aria2'):
            missing.append(name)
    return missing


def productDir(archiveDir, name, products):
    '''
    Map a granule filename to its per-product archive subdir (e.g. RUNW/ROFF)
    by matching a product token in the filename. Returns None if no match.
    '''
    for product in products:
        if f'_{product}_' in name:
            return os.path.join(archiveDir, product)
    return None


# Products fetched slim by downloadNISARoptimized: only the bands the GrIMP
# workflow reads, float32 rounded to 10 mantissa bits (unwrappedPhase 12),
# same internal paths, ~4-5x smaller than the archive product and about half the
# transfer. slimDownload: false in autoupdate.yaml fetches whole products instead.
SLIM_PRODUCTS = ('RUNW', 'ROFF', 'GUNW', 'GOFF')
# Granules downloaded per run (maxDownloads in autoupdate.yaml). Above the busiest normal
# night (638 on 2026-10-19; median ~280), so it only bites when ASF releases a backlog, which
# is then cleared over several nights, oldest acquisitions first.
DEFAULT_MAX_DOWNLOADS = 800
# Slim-download rate cap in MB/s (downloadLimit in autoupdate.yaml; 0 for none). The tool's own
# default is 19; the nightly job runs overnight, when the whole-file aria2c downloads it replaced
# ran uncapped (~49 MB/s), so it allows more (Ian, 2026-11-06).
DEFAULT_DOWNLOAD_LIMIT = 30
SLIM_MANIFEST = 'slimManifest.jsonl'


def collectSlimManifest(archiveDir, downloadDir):
    '''
    Tidy a staging folder after a downloadNISARoptimized pass. Its manifest and
    (empty) .scratch dir are moved/removed because releaseStaging only moves
    *.h5 and removes a folder once it is empty; the manifest is appended to
    archiveDir/slimManifest.jsonl instead (a top-level file, so the searchASF
    dedup glob archiveDir/*/* never sees it). A granule that failed
    verification is deleted, so it counts as missing: retried this run and,
    failing that, re-listed by tomorrow's search.
    '''
    manifest = os.path.join(downloadDir, SLIM_MANIFEST)
    if os.path.exists(manifest):
        with open(manifest) as fp:
            lines = fp.readlines()
        for line in lines:
            entry = json.loads(line)
            if entry['problems']:
                bad = os.path.join(downloadDir, f'{entry["stem"]}.h5')
                print(f'slim granule failed verification, removed: {entry["stem"]} '
                      f'({"; ".join(entry["problems"])})')
                if os.path.exists(bad):
                    os.remove(bad)
        with open(os.path.join(archiveDir, SLIM_MANIFEST), 'a') as fp:
            fp.writelines(lines)
        os.remove(manifest)
    scratch = os.path.join(downloadDir, '.scratch')
    if os.path.isdir(scratch) and not os.listdir(scratch):
        os.rmdir(scratch)


def downloadPass(config, urls, downloadDir):
    '''
    One download pass over urls into downloadDir: RUNW/ROFF slim through
    downloadNISARoptimized (one call per product, since it takes a single
    --product), anything else whole through ariaDownload.
    '''
    slim = config.get('slimDownload', True)
    groups = {}
    for url in urls:
        name = url.split('/')[-1]
        product = next((p for p in SLIM_PRODUCTS if slim and f'_{p}_' in name), None)
        groups.setdefault(product, []).append(url)
    with tempfile.TemporaryDirectory() as tmp:
        for product, group in groups.items():
            listFile = os.path.join(tmp, f'urls.{product or "whole"}')
            with open(listFile, 'w') as fp:
                fp.write('\n'.join(group) + '\n')
            if product is None:
                subprocess.run(['ariaDownload', listFile], cwd=downloadDir,
                               check=True)
            else:
                # No check=True: exit status 69 only means some granule failed
                # verification, which collectSlimManifest deals with.
                subprocess.run(['downloadNISARoptimized', '--product', product,
                                '--urls', listFile, '--outputDir', downloadDir,
                                '--limit', str(config.get('downloadLimit',
                                                          DEFAULT_DOWNLOAD_LIMIT))])
    collectSlimManifest(config['archiveDir'], downloadDir)


def capDownloads(names, maxDownloads):
    '''
    The first maxDownloads of names, oldest acquisition first, never splitting the
    products of one pair (a RUNW and its ROFF share every field but the product).
    Returns (kept, deferred). Deferred granules are simply not downloaded, so the next
    night's search lists them again.
    '''
    if maxDownloads is None or len(names) <= maxDownloads:
        return list(names), []
    groups = {}
    for name in names:
        f = name.split('_')
        key = '_'.join(f[:3] + f[4:]) if len(f) > 4 else name
        groups.setdefault(key, []).append(name)
    # field 11 of a pair-product name is the reference start time (YYYYMMDDTHHMMSS)
    order = sorted(groups, key=lambda k: (k.split('_')[10] if len(k.split('_')) > 10 else k, k))
    kept, deferred = [], []
    for key in order:
        if not deferred and len(kept) + len(groups[key]) <= maxDownloads:
            kept += groups[key]
        else:
            deferred += groups[key]
    return kept, deferred


def downloadGranules(config, urlListFile, today, maxAttempts=3):
    '''
    Download the granules listed in urlListFile into today's dated staging
    folder archiveDir/<MM-DD-YYYY>, retrying incomplete ones (aria2c sometimes
    times out; re-running resumes/skips as needed). RUNW/ROFF are fetched slim
    (see SLIM_PRODUCTS / downloadPass).

    `today` (a datetime.date) is fixed once at program start, so the folder is
    named for the run's start date even if the (multi-hour) download crosses
    midnight -- the reconcile/release step keys off the newest dated folder, not
    the wall clock, so a run spanning midnight stages and releases correctly.

    The granules stay in the dated folder for the one-day completion buffer
    (see consolidateStaging/releaseStaging) instead of being moved straight
    into the per-product RUNW/ROFF subdirs. The dated folder sits one level
    under archiveDir so the searchASF dedup glob (archiveDir/*/*) still covers
    staged-but-not-yet-released granules and they are not re-downloaded.
    '''
    archiveDir = config['archiveDir']
    downloadDir = os.path.join(archiveDir, today.strftime(STAGING_DATE_FMT))
    os.makedirs(downloadDir, exist_ok=True)
    names = readUrlList(urlListFile)
    with open(urlListFile) as fp:
        urlOf = {line.strip().split('/')[-1]: line.strip() for line in fp if line.strip()}
    names, deferred = capDownloads(names, config.get('maxDownloads', DEFAULT_MAX_DOWNLOADS))
    if deferred:
        print(f'maxDownloads: {len(names)} granule(s) this run, {len(deferred)} deferred '
              'to later runs (oldest acquisitions first)')
    print(f'{len(names)} granule(s) to download into {downloadDir}')
    # Both downloaders skip a granule already on disk; loop until everything is
    # complete, re-fetching only what is still missing.
    missing = names
    for attempt in range(1, maxAttempts + 1):
        if attempt > 1:
            # ariaDownload skips any file that already exists, so a partial
            # .h5 left by a timed-out aria2c would never be retried (or
            # resumed). Remove the partial and its .aria2 control file (and a
            # killed slim repack's .partial) so the retry re-downloads it.
            for name in missing:
                for stale in (os.path.join(downloadDir, name),
                              os.path.join(downloadDir, f'{name}.aria2'),
                              os.path.join(downloadDir, f'{name}.partial')):
                    if os.path.exists(stale):
                        os.remove(stale)
        print(f'download attempt {attempt}/{maxAttempts} '
              f'({len(missing)} remaining)')
        downloadPass(config, [urlOf[n] for n in missing], downloadDir)
        missing = incompleteGranules(downloadDir, names)
        if not missing:
            break
        print(f'{len(missing)} granule(s) still incomplete after attempt '
              f'{attempt}')
    if missing:
        # Remove the leftover partials entirely: a truncated .h5 left on disk
        # would (a) be released into RUNW/ROFF as a corrupt product and (b) be
        # indexed by searchASF's dedup glob so the granule would never be
        # re-downloaded. With the partial gone, tomorrow's search re-lists it.
        for name in missing:
            for stale in (os.path.join(downloadDir, name),
                          os.path.join(downloadDir, f'{name}.aria2'),
                          os.path.join(downloadDir, f'{name}.partial')):
                if os.path.exists(stale):
                    os.remove(stale)
        print(f'WARNING: {len(missing)} granule(s) did not download after '
              f'{maxAttempts} attempts; partials removed so tomorrow\'s '
              'search re-downloads them')
    print(f'{len(names) - len(missing)} granule(s) staged in {downloadDir}')
    return missing


# ---------------------------------------------------------------------------
# One-day completion buffer: dated staging folders under archiveDir.
#
# Each daily run downloads into archiveDir/<today>. A track/cycle pass can be
# split across two daily downloads, so before releasing a day's data into the
# general population (RUNW/ROFF) we hold it one day and merge in any matching
# frames that arrive the next day. "Release" then moves the aged (prior-day)
# folders into the per-product subdirs where FileNISARProducts consumes them.
# ---------------------------------------------------------------------------

def stagingFolders(archiveDir):
    '''
    Return [(dateObj, path), ...] for every archiveDir/<MM-DD-YYYY> staging
    folder, sorted ascending by date. Non-dated dirs (RUNW, ROFF,
    searchResults, gpkg, ...) do not parse as a date and are ignored.
    '''
    folders = []
    for path in glob.glob(os.path.join(archiveDir, '*')):
        if not os.path.isdir(path):
            continue
        try:
            d = datetime.datetime.strptime(os.path.basename(path),
                                           STAGING_DATE_FMT).date()
        except ValueError:
            continue
        folders.append((d, path))
    return sorted(folders, key=lambda x: x[0])


def passKeys(folder):
    '''
    Return {(cycle, track): [filenames]} for the granules in a staging folder,
    keyed by acquisition pass. cycle/track come from parseFileName (leading
    zeros stripped), so they match regardless of filename zero-padding.
    '''
    passes = {}
    for path in glob.glob(os.path.join(folder, '*.h5')):
        name = os.path.basename(path)
        if os.path.exists(f'{path}.aria2'):
            u.mywarning(f'passKeys: {name} is an incomplete download '
                        '(.aria2 control file present); ignoring')
            continue
        # Key accesses stay INSIDE the try: parseFileName zips positionally,
        # so a short/odd name returns a dict missing keys rather than raising.
        try:
            pDict, _ = parseFileName(name)
            key = (pDict['cycle'], pDict['track'])
        except (ValueError, KeyError):
            u.mywarning(f'passKeys: cannot parse {name}; skipping')
            continue
        passes.setdefault(key, []).append(name)
    return passes


def consolidateStaging(archiveDir):
    '''
    Merge today's frames of any pass that also appears in the aged (non-today)
    staging folders into the most-recent prior folder, so that pass is released
    complete this run. Pass identity is (cycle, track) -- NOT the reference
    date, which can flip mid-track at a UTC-midnight boundary.

    Today's frames are matched against the UNION of pass keys across ALL prior
    folders, not just the immediately-preceding day. Normally there is only one
    prior folder (yesterday), but if cron runs were skipped, several dated
    folders can accumulate; releaseStaging then sweeps them all, so today's
    frames of a pass that first appeared in any of them belong with that release.
    '''
    folders = stagingFolders(archiveDir)
    if len(folders) < 2:
        return
    todayFolder = folders[-1][1]
    priorFolders = [path for _, path in folders[:-1]]
    target = folders[-2][1]  # most-recent prior; all priors are released together
    priorPassKeys = set()
    for folder in priorFolders:
        priorPassKeys.update(passKeys(folder).keys())
    todayPasses = passKeys(todayFolder)
    for key, names in todayPasses.items():
        if key not in priorPassKeys:
            continue
        for name in names:
            src = os.path.join(todayFolder, name)
            dst = os.path.join(target, name)
            if os.path.exists(dst):          # idempotent re-run
                continue
            shutil.move(src, dst)
        print(f'consolidateStaging: merged cycle{key[0]}/track{key[1]} '
              f'({len(names)} file(s)) from {os.path.basename(todayFolder)} '
              f'into {os.path.basename(target)}')


def releaseStaging(archiveDir, products):
    '''
    Release every aged staging folder (all but today's newest) into the
    per-product subdirs (RUNW/ROFF) via productDir(), then remove the emptied
    folders. Returns the sorted list of released granule basenames -- the work
    list for the processing chain. Idempotent: files already present at the
    destination are skipped, and an already-emptied prior folder releases
    nothing.
    '''
    folders = stagingFolders(archiveDir)
    released = set()
    priorFolders = [path for _, path in folders[:-1]]  # exclude today (newest)
    for folder in priorFolders:
        for path in glob.glob(os.path.join(folder, '*.h5')):
            name = os.path.basename(path)
            if os.path.exists(f'{path}.aria2'):
                u.mywarning(f'releaseStaging: {name} is an incomplete '
                            'download (.aria2 control file present); '
                            'not releasing it')
                continue
            dest = productDir(archiveDir, name, products)
            if dest is None:
                u.mywarning(f'releaseStaging: no product subdir for {name}; '
                            'leaving in place')
                continue
            os.makedirs(dest, exist_ok=True)
            target = os.path.join(dest, name)
            if os.path.exists(target):
                os.remove(path)              # already released; drop the stage copy
            else:
                shutil.move(path, target)
            released.add(name)
        # Remove the folder if now empty (no leftover non-h5 clutter).
        if not os.listdir(folder):
            os.rmdir(folder)
    if released:
        print(f'releaseStaging: released {len(released)} granule(s) from '
              f'{len(priorFolders)} aged folder(s) into RUNW/ROFF')
    return sorted(released)


def orbitsFromReleased(released, projectDir):
    '''
    Map released granule names to the (trackDir, orbit1) pairs that need
    (re)processing, deriving orbit1 = orbitFromCycleTrack(cycle, track).
    Skips pairs whose track-<track>/<orbit1>_* frame dirs do not exist yet
    (FileNISARProducts must have created them first); warns if the formula's
    directory is unexpectedly absent so the mismatch is visible.
    '''
    seen, result = set(), []
    for name in released:
        # Key accesses inside the try: parseFileName zips positionally, so a
        # short/odd name returns a dict missing keys rather than raising.
        try:
            pDict, _ = parseFileName(name)
            track = pDict['track']
            orbit = orbitFromCycleTrack(pDict['cycle'], track)
        except (ValueError, KeyError):
            u.mywarning(f'orbitsFromReleased: cannot parse {name}; '
                        'it will drive no orbit processing')
            continue
        if (track, orbit) in seen:
            continue
        seen.add((track, orbit))
        trackDir = os.path.join(projectDir, f'track-{track}')
        if glob.glob(os.path.join(trackDir, f'{orbit}_*')):
            result.append((trackDir, orbit))
        else:
            u.mywarning(f'orbitsFromReleased: no {trackDir}/{orbit}_* frame '
                        f'dirs for cycle{pDict["cycle"]}/track{track} '
                        '(FileNISARProducts did not create them?); skipping')
    return result


def yearsFromReleased(released):
    '''
    Return the sorted list of calendar years (as strings) of the released
    granules' reference acquisition dates -- passed to setupNISARTracks --year.
    Union across granules, so a pass straddling New Year contributes both years.
    '''
    years = set()
    for name in released:
        try:
            pDict, _ = parseFileName(name)
            years.add(pDict['date1Sstart'].year)
        except (ValueError, KeyError):
            u.mywarning(f'yearsFromReleased: cannot parse {name}; '
                        'it contributes no --year')
    return [str(y) for y in sorted(years)]


def framesByOrbit(released):
    '''Map (track, orbit) -> sorted list of frame numbers released this run,
    for the run summary. orbit derived from (cycle, track).'''
    m = {}
    for name in released:
        try:
            pDict, _ = parseFileName(name)
            orbit = orbitFromCycleTrack(pDict['cycle'], pDict['track'])
            m.setdefault((pDict['track'], orbit),
                         set()).add(int(pDict['frame']))
        except (ValueError, KeyError):
            continue    # already warned in orbitsFromReleased
    return {k: sorted(v) for k, v in m.items()}


def readPendingGranules(pendingFile):
    '''Return the granule names left unprocessed by a previous run (empty list
    if none). The pending file makes the released->processed handoff
    crash-safe: releaseStaging moves files into RUNW/ROFF and removes the
    dated folders, so if a run died before SetupNISAR those granules would
    otherwise never drive any processing (they exist only in that run's
    in-memory released list).'''
    if not os.path.exists(pendingFile):
        return []
    with open(pendingFile) as fp:
        return [line.strip() for line in fp if line.strip()]


def writePendingGranules(pendingFile, names):
    '''Write the to-be-processed granule list (one basename per line).'''
    os.makedirs(os.path.dirname(pendingFile), exist_ok=True)
    with open(pendingFile, 'w') as fp:
        fp.write('\n'.join(names) + '\n')


def _runSetupNISAR(command, cwd, orbit, track, frames, logPath, results):
    '''Thread worker: run one SetupNISAR orbit subprocess with stdout+stderr
    streamed to logPath as it runs (survives a mid-run crash/reboot, allows
    tail -f, and avoids interleaving across parallel orbits), appending a
    per-orbit result dict to the shared results list (list.append is atomic
    under the GIL, so no lock is needed). On failure the reason recorded is
    the tail of the log file.'''
    try:
        with open(logPath, 'w') as fp:
            proc = subprocess.run(command, cwd=cwd, stdout=fp,
                                  stderr=subprocess.STDOUT)
    except OSError as e:
        # Log file unwritable (bad NFS?) -- still run the orbit, unlogged.
        u.mywarning(f'_runSetupNISAR: could not write {logPath} ({e}); '
                    'running without a log')
        proc = subprocess.run(command, cwd=cwd)
    entry = {'orbit': orbit, 'track': track, 'frames': frames,
             'returncode': proc.returncode, 'log': logPath,
             'status': 'SUCCESS' if proc.returncode == 0 else 'FAILED'}
    if proc.returncode != 0:
        try:
            with open(logPath) as fp:
                tail = fp.read().strip().splitlines()[-25:]
        except OSError:
            tail = []
        entry['reason'] = '\n'.join(tail) or '(no output)'
    results.append(entry)


def mailReport(recipient, subject, body):
    '''
    Mail a report via the local MTA (mail, then mailx). Best effort: a machine
    with no working mailer must not fail the run, so any failure is warned and
    swallowed. Mirrors asfSearchAndDownload.autoupdateS1.mailReport.
    '''
    for mailer in (['mail', '-s', subject, recipient],
                   ['mailx', '-s', subject, recipient]):
        try:
            proc = subprocess.run(mailer, input=body, text=True,
                                  capture_output=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            continue
        if proc.returncode == 0:
            print(f'autoupdate: mailed "{subject}" to {recipient}')
            return True
        u.mywarning(f'mailReport: {mailer[0]} exited {proc.returncode}: '
                    f'{proc.stderr.strip()}')
    u.mywarning(f'mailReport: could not mail "{subject}" to {recipient}: '
                'no working mailer')
    return False


# Error signatures worth surfacing out of a track's stderr: a myerror banner
# (utilities.myerror / its local copies), an argparse rejection, the last line
# of a traceback, or the C binaries' own error() format (clib/standard.c prints
# a bare "error:  <msg>" -- no ANSI colour, no program prefix, so it matches
# none of the others).
_MYERROR_RE = re.compile(r'\x1b\[1;31m\s*\*\*\*\s+(.*?)\s+\*\*\*')
_ARGPARSE_RE = re.compile(r'^(\S+): error: (.+)$')
_TRACEBACK_RE = re.compile(r'^(\w*(?:Error|Exception)): (.+)$')
_CERROR_RE = re.compile(r'^error:\s+(.+)$')

# Keep the mailed summary bounded no matter how badly a run goes.
MAX_ERROR_KINDS = 15
MAX_EXAMPLE_TRACKS = 5

# Per-project state remembering which (error, track) pairs are already known,
# so a long-standing failure stops re-reporting in full every night. Lives in
# dailyDownloadLogs beside the summaries; delete it to re-report everything.
KNOWN_ERRORS_FILE = 'knownTrackErrors.json'


def normalizeError(message):
    '''
    Collapse the run-specific parts of an error message (directory prefixes,
    long digit runs like orbit/frame numbers, parenthesised counts) so the same
    failure hitting 150 tracks aggregates to one line instead of 150
    near-identical ones.

    The parenthesised-count rule is what keeps "tiepoints: Insufficient Number
    (0|1|2|3) of Valid tie points" -- one failure mode, four spellings -- from
    occupying four lines and four separate entries in the known-error state.
    The count is the detail, not the identity; the exact value is still in
    track-N/stderr.
    '''
    message = re.sub(r'(/\S+)+/', '', message)
    message = re.sub(r'\d{3,}', 'N', message)
    message = re.sub(r'\((\d+)\)', '(N)', message)
    return ' '.join(message.split())[:120]


def trackStderrLabel(path):
    ''' track-N/stderr -> "track-N". '''
    return os.path.basename(os.path.dirname(path))


def frameStderrLabel(path):
    ''' track-N/<frame>/motion/azest_stderr -> "track-N/<frame> (azest)". '''
    frameDir = os.path.dirname(os.path.dirname(path))
    tool = os.path.basename(path).replace('_stderr', '').replace('stderr', '')
    label = (f'{os.path.basename(os.path.dirname(frameDir))}/'
             f'{os.path.basename(frameDir)}')
    return f'{label} ({tool})' if tool else label


def collectTrackErrors(projectDir, since=None, scanned=None):
    '''
    Aggregate error lines out of the log files the tie stage writes, keyed by
    normalised message -> [where it happened]:

      track-N/stderr                       the per-track tie chain
      track-N/<frame>/motion/*stderr       azparams/rparams per frame

    The second set is where azest/estr failures land ("svInitAzParams :
    insufficient points to fit"), and nothing was reading it -- a frame whose
    azparams dies leaves a zero-byte az.est.yaml that makes mosaic3d fail
    obscurely much later. Both are rewritten on every run, so restrict to
    files touched since `since` (the run start) -- otherwise a run that
    aborted before the tie stage would report yesterday's errors as today's.
    Each distinct message is counted once per file, however often it repeats.

    If a set is passed as `scanned`, every label whose log WAS read this run is
    added to it -- saveKnownErrors needs that to tell "this track was checked
    and is clean now" from "this track was not touched this run".
    '''
    errors = {}
    sources = [(os.path.join(projectDir, 'track-*', 'stderr'),
                trackStderrLabel),
               (os.path.join(projectDir, 'track-*', '*', 'motion', '*stderr'),
                frameStderrLabel)]
    paths = [(p, label) for pattern, label in sources
             for p in sorted(glob.glob(pattern))]
    for path, label in paths:
        try:
            if since is not None and datetime.datetime.fromtimestamp(
                    os.path.getmtime(path)) < since:
                continue
            with open(path, errors='replace') as fp:
                text = fp.read()
        except OSError:
            continue
        track = label(path)
        if scanned is not None:
            scanned.add(track)
        found = set()
        for match in _MYERROR_RE.finditer(text):
            found.add(normalizeError(match.group(1)))
        for line in text.splitlines():
            match = _CERROR_RE.match(line.strip())
            if match:
                found.add(normalizeError(match.group(1)))
                continue
            for regex in (_ARGPARSE_RE, _TRACEBACK_RE):
                match = regex.match(line.strip())
                if match:
                    found.add(normalizeError(
                        f'{match.group(1)}: {match.group(2)}'))
        for message in found:
            errors.setdefault(message, []).append(track)
    return errors


def loadKnownErrors(logDir):
    '''
    Read the known-error state written by saveKnownErrors:
    {normalised message: {track label: first-seen MM-DD-YYYY}}. A missing or
    unreadable file reads as empty, which makes every error look new -- the
    safe direction (over-report rather than silently swallow).
    '''
    path = os.path.join(logDir, KNOWN_ERRORS_FILE)
    try:
        with open(path) as fp:
            data = json.load(fp)
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    return {message: dict(tracks) for message, tracks in data.items()
            if isinstance(tracks, dict)}


def splitTrackErrors(trackErrors, known):
    '''
    Partition this run's errors into (message, track) pairs never seen before
    and pairs that were already failing. Returns (new, recurring, firstSeen):
    new and recurring map message -> sorted tracks, firstSeen maps message ->
    the earliest date any of its recurring tracks started failing. A message
    with both new and old tracks appears in both, listing only its own tracks
    in each -- a known failure spreading to a new track is news.
    '''
    new, recurring, firstSeen = {}, {}, {}
    for message, tracks in trackErrors.items():
        seen = known.get(message, {})
        for track in sorted(set(tracks)):
            if track in seen:
                recurring.setdefault(message, []).append(track)
            else:
                new.setdefault(message, []).append(track)
        if message in recurring:
            firstSeen[message] = min(seen[t] for t in recurring[message])
    return new, recurring, firstSeen


def saveKnownErrors(logDir, trackErrors, known, scanned, today):
    '''
    Rewrite logDir/knownTrackErrors.json so it matches what is currently
    failing: pairs still failing keep their original first-seen date, pairs
    seen for the first time get today's.

    A pair is forgotten only when that track's log was actually re-read this
    run (`scanned`) and no longer carries the message. Pairs from tracks the
    run never touched are carried forward untouched -- collectTrackErrors only
    reads logs rewritten since the run start, so a track with no new granules
    contributes nothing, and pruning on that absence would make it report
    "new" again the next time it does fail.
    '''
    stamp = today.strftime(STAGING_DATE_FMT)
    current = {message: set(tracks) for message, tracks in trackErrors.items()}
    updated = {}
    for message, tracks in known.items():
        for track, first in tracks.items():
            if track in current.get(message, ()) or track not in scanned:
                updated.setdefault(message, {})[track] = first
    for message, tracks in current.items():
        for track in tracks:
            updated.setdefault(message, {}).setdefault(track, stamp)
    with open(os.path.join(logDir, KNOWN_ERRORS_FILE), 'w') as fp:
        json.dump(updated, fp, indent=1, sort_keys=True)
    return updated


def runHadErrors(summary):
    '''
    True if the run recorded any failure worth emailing: a step with a non-zero
    returncode, a FAILED SetupNISAR orbit, or an uncaught exception. A step
    recorded as None (skipped this run, e.g. no tie year) is not an error.
    '''
    if summary.get('exception'):
        return True
    if any(rc not in (None, 0) for rc in summary['steps'].values()):
        return True
    if any(e['status'] == 'FAILED' for e in summary['orbitResults']):
        return True
    # Per-track errors are always written to the summary file; whether they are
    # worth waking someone up is a config choice (notifyOnTrackErrors). Only
    # errors NOT seen in an earlier run count -- a known long-standing failure
    # (e.g. a frame whose phase never ties but whose offsets are fine) would
    # otherwise mail every night forever and train the reader to ignore it.
    if summary.get('notifyOnTrackErrors') and summary.get('newTrackErrors'):
        return True
    return False


def notifyOnErrors(config, summary, projectDir, body):
    '''
    Email the run summary if the run had errors and a notifyEmail recipient is
    configured. Silent on a clean run and when notifyEmail is unset (opt-in):
    adding `notifyEmail: you@example.com` to autoupdate.yaml turns it on.
    '''
    recipient = config.get('notifyEmail')
    if not recipient or not runHadErrors(summary):
        return
    host = socket.gethostname()
    nFail = sum(1 for e in summary['orbitResults']
                if e['status'] == 'FAILED')
    stepFails = [k for k, rc in summary['steps'].items()
                 if rc not in (None, 0)]
    parts = []
    if nFail:
        parts.append(f'{nFail} orbit(s) failed')
    if stepFails:
        parts.append(f'{len(stepFails)} step(s) failed')
    if summary.get('exception'):
        parts.append('crashed')
    newErrors = summary.get('newTrackErrors') or {}
    if newErrors:
        nTracks = len({t for tracks in newErrors.values() for t in tracks})
        parts.append(f'{len(newErrors)} new error kind(s) on {nTracks} track(s)')
    subject = (f'autoupdate {os.path.basename(projectDir)}: '
               f'{", ".join(parts) or "errors"} on {host}')
    mailReport(recipient, subject, body or '(run summary unavailable)')


def processReleased(config, released, projectDir, pendingFile=None):
    '''
    Steps 4-7: file the released granules, (re)build the affected orbits'
    virtual frames, refresh tie points/velocity for the relevant year(s), and
    update the master input file. A per-run summary is written to
    projectDir/dailyDownloadLogs/ (even if a step fails). pendingFile (the
    crash-recovery work list written by main()) is cleared once the per-orbit
    SetupNISAR stage has run; if this run dies before that, the next run
    re-derives the same work from the pending file.
    '''
    start = datetime.datetime.now()
    logDir = os.path.join(projectDir, 'dailyDownloadLogs')
    stamp = start.strftime('%m-%d-%Y_%H%M%S')
    summary = {'start': start, 'stamp': stamp, 'released': released,
               'steps': {}, 'orbitResults': [], 'orbitLogDir': None,
               'years': [], 'tieTracks': [],
               'trackErrors': {}, 'newTrackErrors': {},
               'recurringTrackErrors': {}, 'trackErrorFirstSeen': {}}
    try:
        # Step 4: file released RUNW/ROFF into track-<N>/<orbit>_<frame>/H5/.
        # archiveDir passed explicitly so the run does not depend on
        # projectDir containing an autoupdate.yaml with the same archiveDir
        # (FileNISARProducts' own fallback re-reads that file by name).
        print('processReleased: FileNISARProducts')
        r = subprocess.run(['FileNISARProducts',
                            '--inputPath', config['archiveDir']],
                           cwd=projectDir)
        summary['steps']['FileNISARProducts'] = r.returncode
        if r.returncode != 0:
            u.mywarning('processReleased: FileNISARProducts failed; '
                        'aborting the rest of this run (the pending-granules '
                        'file is kept, so the next run retries)')
            return

        # Step 5: per-orbit SetupNISAR for the affected orbits. No --new/
        # --overWrite so already-converted frames are reused while a virtual
        # frame extended by late-arriving frames is rebuilt in place (same
        # _0000; see plan). Region-dependent processing choices (e.g.
        # sepIceRock) live in project.yaml and are read by SetupNISAR itself,
        # so setupFlags defaults to empty.
        #
        # Orbits run up to nThreads at a time (config 'nThreads', default 20).
        # Each orbit is a separate SetupNISAR SUBPROCESS, so the mosaic-side
        # OpenMP global-state caveats do not apply across orbits (no shared
        # memory). ompThreads (config) sets the per-orbit OpenMP width; balance
        # it against nThreads so nThreads*ompThreads does not badly
        # oversubscribe the cores.
        setupFlags = config.get('setupFlags', '')
        if isinstance(setupFlags, str):
            setupFlags = setupFlags.split()
        ompArgs = (['--ompThreads', str(config['ompThreads'])]
                   if config.get('ompThreads') is not None else [])
        nThreads = int(config.get('nThreads', 20))
        orbits = orbitsFromReleased(released, projectDir)
        framesMap = framesByOrbit(released)
        orbitLogDir = os.path.join(logDir, f'orbitLogs_{stamp}')
        if orbits:
            os.makedirs(orbitLogDir, exist_ok=True)
            summary['orbitLogDir'] = orbitLogDir
        results = []
        orbitThreads = []
        for trackDir, orbit in orbits:
            track = os.path.basename(trackDir).replace('track-', '')
            frames = framesMap.get((track, orbit), [])
            command = ['SetupNISAR', str(orbit)] + ompArgs + setupFlags
            logPath = os.path.join(orbitLogDir,
                                   f'orbit_{orbit}_track-{track}.log')
            print(f'processReleased: {" ".join(command)}  (cwd={trackDir})')
            orbitThreads.append(threading.Thread(
                target=_runSetupNISAR,
                args=(command, trackDir, orbit, track, frames, logPath,
                      results)))
        u.runMyThreads(orbitThreads, nThreads, 'SetupNISAR orbits')
        summary['orbitResults'] = results
        nFail = sum(1 for e in results if e['status'] == 'FAILED')
        if nFail:
            u.mywarning(f'processReleased: {nFail} of {len(orbits)} '
                        'SetupNISAR orbit(s) failed -- continuing with '
                        'tie/master steps for the rest')
        # The per-orbit stage has run (failures are recorded in the summary,
        # not retried daily), so the crash-recovery work list is done.
        if pendingFile is not None and os.path.exists(pendingFile):
            os.remove(pendingFile)

        # Step 6: tie points / velocity for the year(s) of the new data.
        # Prime first with --copyFiles: tie_plan_header/vel_thumb_plan are only
        # instantiated under that flag, so a track appearing for the FIRST time
        # would otherwise have no tie plan and its ties would silently not run.
        # Idempotent -- existing headers are skipped -- and it exits right
        # after setup, so it costs almost nothing on normal days.
        years = yearsFromReleased(released)
        summary['years'] = years
        if years:
            print(f'processReleased: setupNISARTracks --copyFiles  '
                  f'(cwd={projectDir})')
            r = subprocess.run(['setupNISARTracks', '--copyFiles'],
                               cwd=projectDir)
            summary['steps']['setupNISARTracks --copyFiles'] = r.returncode
            # Tie only the tracks that actually got new data. Ties and
            # baselines are per-track quantities, so re-tying a track with no
            # new frames reproduces the same answer -- across ~170 tracks that
            # was 1.5-2 h of redundant work every night. setupNISARTracks adds
            # back any track whose autoclean mask it rebuilds. With no orbits
            # (nothing filed) the flag is omitted, keeping the old all-tracks
            # pass rather than silently skipping the tie stage; set
            # retieAllTracks: true in autoupdate.yaml to force that always.
            tieTracks = sorted({os.path.basename(td) for td, _ in orbits},
                               key=lambda t: int(t.split('-')[1]))
            restrictTies = bool(tieTracks) and not config.get('retieAllTracks',
                                                              False)
            summary['tieTracks'] = tieTracks if restrictTies else []
            command = ['setupNISARTracks', '--year'] + years
            if restrictTies:
                command += ['--tieTracks'] + tieTracks
            print(f'processReleased: {" ".join(command)}  (cwd={projectDir})')
            r = subprocess.run(command, cwd=projectDir)
            summary['steps']['setupNISARTracks'] = r.returncode

        # Step 7: update the master input file. By default makeMaster scans
        # projectDir and writes projectDir/Release/masterInput/inputFile. A
        # project whose mosaics read a master assembled from several project
        # dirs (e.g. Antarctica40 + Antarctica80 -> Antarctica2026/Release/
        # track-all) points this step at that master instead via
        # masterInputPaths/masterOutputPath in autoupdate.yaml -- otherwise the
        # nightly run refreshes a master nothing downstream reads. makeMaster
        # errors if the output directory's parent is missing, so create it
        # first (idempotent).
        outputPath = masterOutputPath(config, projectDir)
        inputPaths = masterInputPaths(config, projectDir)
        os.makedirs(outputPath, exist_ok=True)
        command = ['makeMaster', '--projectDir', projectDir,
                   '--outputPath', outputPath]
        if inputPaths:
            command += ['--inputPaths'] + inputPaths
        print(f'processReleased: {" ".join(command)}')
        r = subprocess.run(command, cwd=projectDir)
        summary['steps']['makeMaster'] = r.returncode
    except Exception:
        # Record a crash so notifyOnErrors mails even when a step raised rather
        # than returning non-zero; re-raise so cron still logs the traceback.
        summary['exception'] = traceback.format_exc()
        raise
    finally:
        summary['end'] = datetime.datetime.now()
        # Gather per-track errors here rather than after the tie step, so they
        # are still collected when a step above raised.
        try:
            summary['notifyOnTrackErrors'] = bool(
                config.get('notifyOnTrackErrors', True))
            scanned = set()
            summary['trackErrors'] = collectTrackErrors(
                projectDir, since=summary['start'], scanned=scanned)
            # Split into new vs already-known so the summary (and the mail
            # trigger) surface what changed rather than the whole standing
            # backlog again.
            known = loadKnownErrors(logDir)
            (summary['newTrackErrors'],
             summary['recurringTrackErrors'],
             summary['trackErrorFirstSeen']) = splitTrackErrors(
                summary['trackErrors'], known)
            saveKnownErrors(logDir, summary['trackErrors'], known, scanned,
                            summary['start'])
        except Exception as e:                 # never let logging break a run
            u.mywarning(f'processReleased: could not collect track errors: {e}')
        body = None
        try:
            body = writeRunSummary(logDir, summary)
        except Exception as e:                 # never let logging break a run
            u.mywarning(f'processReleased: could not write run summary: {e}')
        try:
            notifyOnErrors(config, summary, projectDir, body)
        except Exception as e:                 # never let mailing break a run
            u.mywarning(f'processReleased: could not send notify email: {e}')


def trackErrorLines(errors, indent='  '):
    '''
    Render {message: [tracks]} as the detailed two-line-per-kind block: a
    count-prefixed message and a few example tracks, worst first, capped so a
    bad night cannot produce an unmailable summary.
    '''
    lines = []
    ranked = sorted(errors.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for message, tracks in ranked[:MAX_ERROR_KINDS]:
        shown = sorted(set(tracks))
        examples = ', '.join(shown[:MAX_EXAMPLE_TRACKS])
        if len(shown) > MAX_EXAMPLE_TRACKS:
            examples += f', +{len(shown) - MAX_EXAMPLE_TRACKS} more'
        lines.append(f'{indent}[{len(shown):4d}] {message}')
        lines.append(f'{indent}       {examples}')
    if len(ranked) > MAX_ERROR_KINDS:
        lines.append(f'{indent}... and {len(ranked) - MAX_ERROR_KINDS} more '
                     'distinct message(s)')
    return lines


def writeRunSummary(logDir, summary):
    '''
    Write a human-readable per-run summary to
    logDir/summary_<MM-DD-YYYY>_<HHMMSS>.txt: a quick overview (orbits
    succeeded/failed, frames, granules, step statuses) followed by a per-orbit
    breakdown with frame counts, status, and -- for failures -- the returncode,
    a pointer to the full per-orbit log, and the tail of that log as the reason.

    Returns the summary text so the caller can mail it (notifyOnErrors).
    '''
    os.makedirs(logDir, exist_ok=True)
    path = os.path.join(logDir, f'summary_{summary["stamp"]}.txt')
    released = summary['released']
    nRUNW = sum('_RUNW_' in n for n in released)
    nROFF = sum('_ROFF_' in n for n in released)
    results = summary['orbitResults']
    nOK = sum(1 for e in results if e['status'] == 'SUCCESS')
    nFail = sum(1 for e in results if e['status'] == 'FAILED')
    nFrames = sum(len(e['frames']) for e in results)
    steps = summary['steps']

    def stepStr(name):
        rc = steps.get(name)
        if rc is None:
            return 'not run'
        return 'OK' if rc == 0 else f'FAILED (rc={rc})'

    lines = ['NISAR autoupdate run summary',
             f'  started : {summary["start"]:%Y-%m-%d %H:%M:%S}',
             f'  finished: {summary["end"]:%Y-%m-%d %H:%M:%S}  '
             f'(duration {summary["end"] - summary["start"]})',
             '',
             'OVERVIEW',
             f'  Orbits: {len(results)} processed -- {nOK} succeeded, '
             f'{nFail} failed',
             f'  Frames in processed orbits: {nFrames}',
             f'  Granules released: {len(released)}  '
             f'({nRUNW} RUNW, {nROFF} ROFF)']
    if summary['years']:
        lines.append(f'  Tie/velocity year(s): {" ".join(summary["years"])}')
    # Only set when the tie stage was restricted to the tracks with new data;
    # absent means every track was tiepointed (the pre-restriction behaviour).
    if summary.get('tieTracks'):
        lines.append(f'  Tracks tiepointed: {len(summary["tieTracks"])} '
                     f'({" ".join(summary["tieTracks"])})')
    lines += [f'  FileNISARProducts : {stepStr("FileNISARProducts")}',
              f'  setupNISARTracks --copyFiles : '
              f'{stepStr("setupNISARTracks --copyFiles")}',
              f'  setupNISARTracks  : {stepStr("setupNISARTracks")}',
              f'  makeMaster        : {stepStr("makeMaster")}']

    # Per-track errors, aggregated. These come from track-*/stderr, which is
    # where a failing tie chain reports and where nothing else was looking --
    # a broken tie stage went unnoticed for 11 nights because every step above
    # still said OK.
    trackErrors = summary.get('trackErrors') or {}
    if trackErrors:
        nTracks = len({t for tracks in trackErrors.values() for t in tracks})
        newErrors = summary.get('newTrackErrors') or {}
        recurring = summary.get('recurringTrackErrors') or {}
        firstSeen = summary.get('trackErrorFirstSeen') or {}
        if not newErrors and not recurring:
            # Split never ran (state unreadable, or an older summary dict):
            # report everything in full rather than showing an empty section.
            newErrors = {message: sorted(set(tracks))
                         for message, tracks in trackErrors.items()}
        lines += ['',
                  f'TRACK ERRORS  ({len(newErrors)} new, {len(recurring)} '
                  f'recurring; {nTracks} track(s) affected)']
        if newErrors:
            lines.append('  NEW')
            lines += trackErrorLines(newErrors)
        # Recurring kinds collapse to one line each: the count and how long
        # they have been failing is all that is actionable once the tracks
        # behind them have been looked at (they are all in the state file).
        if recurring:
            lines.append(f'  RECURRING  (already reported; per-track detail in '
                         f'dailyDownloadLogs/{KNOWN_ERRORS_FILE})')
            ranked = sorted(recurring.items(),
                            key=lambda kv: (-len(kv[1]), kv[0]))
            for message, tracks in ranked[:MAX_ERROR_KINDS]:
                since = firstSeen.get(message)
                sinceStr = f'   since {since}' if since else ''
                lines.append(f'  [{len(tracks):4d}] {message}{sinceStr}')
            if len(ranked) > MAX_ERROR_KINDS:
                lines.append(f'  ... and {len(ranked) - MAX_ERROR_KINDS} more '
                             'recurring message(s)')
        lines.append('  (full text: <projectDir>/track-<N>/stderr)')

    lines += ['', 'ORBITS  (failures first)']
    if not results:
        lines.append('  (none)')
    # Failures first, then by orbit number.
    for e in sorted(results, key=lambda x: (x['status'] == 'SUCCESS',
                                            int(x['orbit']))):
        fr = e['frames']
        frStr = f'{len(fr)} frame(s)' + (f' {fr}' if fr else '')
        lines.append(f'  orbit {e["orbit"]} (track-{e["track"]}): '
                     f'{frStr}  {e["status"]}')
        if e['status'] == 'FAILED':
            lines.append(f'      returncode {e["returncode"]}   '
                         f'log: {e.get("log", "")}')
            for rl in e.get('reason', '').splitlines():
                lines.append(f'      | {rl}')
    if summary['orbitLogDir']:
        lines += ['', f'Full per-orbit logs: {summary["orbitLogDir"]}']
    body = '\n'.join(lines) + '\n'
    with open(path, 'w') as fp:
        fp.write(body)
    print(f'processReleased: wrote run summary {path}')
    return body


def main():
    '''
    Run the automated update driven by autoupdate.yaml.
    '''
    args = parseArgs()
    config = loadConfig(args.config)
    # Fail with a clear message on missing keys (a raw KeyError traceback in
    # cron.log is much harder to diagnose). The search keys are only needed
    # when actually searching/downloading.
    required = ['archiveDir', 'products']
    if not args.noDownload:
        required += ['lookBack', 'bandwidths', 'region']
    missingKeys = [k for k in required if k not in config]
    if missingKeys:
        u.myerror(f'autoupdate: required key(s) {missingKeys} missing '
                  f'from {args.config}')
    # Anchor paths: a relative archiveDir would silently break ariaDownload
    # (run with cwd=downloadDir) and the searchASF dedup glob.
    config['archiveDir'] = os.path.abspath(config['archiveDir'])
    archiveDir = config['archiveDir']
    products = granuleProducts(config)
    projectDir = os.path.abspath(config.get('projectDir') or os.path.dirname(
        os.path.abspath(args.config)))

    # Fix "today" once, at program start, so the search window, URL-list name,
    # and download-folder name all use the run's start date even if the run
    # (a multi-hour download) crosses midnight.
    today = datetime.date.today()
    if not args.noDownload:
        urlListFile = searchGranules(config, today)
        downloadGranules(config, urlListFile, today)

    consolidateStaging(archiveDir)                    # step 2
    released = releaseStaging(archiveDir, products)   # step 3
    # Crash recovery: merge in granules a previous run released but never
    # processed (died before/at SetupNISAR). Re-listing an already-processed
    # granule is harmless -- per-frame products self-skip and the virtual
    # frame just rebuilds cheaply.
    pendingFile = os.path.join(projectDir, 'dailyDownloadLogs',
                               'pendingGranules.txt')
    pending = readPendingGranules(pendingFile)
    if pending:
        print(f'autoupdate: {len(pending)} granule(s) pending from a '
              'previous incomplete run')
    workList = sorted(set(released) | set(pending))
    if not workList:
        print('autoupdate: no aged granules to release this run; done.')
        return
    writePendingGranules(pendingFile, workList)
    processReleased(config, workList, projectDir,
                    pendingFile=pendingFile)          # steps 4-7


if __name__ == '__main__':
    main()
