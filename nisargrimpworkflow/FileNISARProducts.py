#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 08:58:04 2026

@author: ian
"""

import argparse
import os
import sys
import glob
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path
from datetime import datetime

import yaml
import nisarhdf


def parseArgs():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mOrganise downloaded NISAR HDF5 products into '
        'the track-based directory tree expected by SetupNISAR.\033[0m\n\n'
        'Expects products pre-sorted by type under inputPath:\n'
        '  inputPath/RUNW/*.h5\n'
        '  inputPath/ROFF/*.h5\n'
        '  inputPath/RIFG/*.h5\n'
        '  inputPath/RSLC/*.h5  (optional)\n\n'
        'For each RUNW the script (1) creates track-{N}/source/ and symlinks '
        'all companion products there, then (2) reads orbit and frame from the '
        'RUNW HDF5 metadata and creates track-{N}/{orbit1}_{frame}/ symlinks '
        'for the L2 products (RUNW, ROFF, RIFG) needed by SetupNISAR. '
        'When multiple products exist for the same track/frame the one with the '
        'shortest temporal baseline is kept; ties are broken by newest '
        'modification time. Extras go into track-{N}/unfiled/ with a log entry. '
        'Mixed-mode frames are skipped. Run once on a fresh download before '
        'calling SetupNISAR.',
        epilog='Example:\n'
               '  FileNISARProducts /data/nisar/downloads '
               '--outputPath /data/nisar/orbits\n\n'
               'Output layout:\n'
               '  outputPath/\n'
               '    track-64/\n'
               '      source/\n'
               '        NISAR_L1_PR_RUNW_....h5  (symlink)\n'
               '        NISAR_L1_PR_ROFF_....h5  (symlink)\n'
               '        NISAR_L1_PR_RIFG_....h5  (symlink)\n'
               '        NISAR_L1_PR_RSLC_....h5  (symlink, if present)\n'
               '      12345_010/\n'
               '        H5/\n'
               '          NISAR_L1_PR_RUNW_....h5  (symlink)\n'
               '          NISAR_L1_PR_ROFF_....h5  (symlink)\n'
               '          NISAR_L1_PR_RIFG_....h5  (symlink)\n'
               '      unfiled/\n'
               '        NISAR_L1_PR_RUNW_....h5  (symlink, duplicate)\n'
               '        log\n'
               '    track-71/...\n'
               '\nPart of the nisargrimpworkflow package.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputPath', type=str, default=None,
                        help='Root directory with RUNW/, ROFF/, RIFG/, RSLC/ '
                        'product subdirectories. If omitted, archiveDir from '
                        'an autoupdate.yaml in --outputPath is used')
    parser.add_argument('--firstOrbit', type=int, default=1,
                        help='Skip orbits numbered below this value [1]')
    parser.add_argument('--lastOrbit', type=int, default=999999,
                        help='Skip orbits numbered above this value [999999]')
    parser.add_argument('--outputPath', type=str, default='.',
                        help='Root directory for the track-{N}/ output tree '
                        '[current directory]')
    parser.add_argument('--products', nargs='+',
                        choices=['RUNW', 'ROFF', 'RIFG', 'RSLC'],
                        default=None, metavar='TYPE',
                        help='Product types to file (RUNW ROFF RIFG RSLC). '
                        'Default: all types. RUNW is always used to determine '
                        'orbit/frame; specifying other types without RUNW '
                        'adds them alongside existing RUNW symlinks.')
    parser.add_argument('--firstDate', type=str, default=None,
                        metavar='YYYYMMDD',
                        help='Skip products whose reference (first) date is '
                        'before this date [no lower limit]')
    parser.add_argument('--lastDate', type=str, default=None,
                        metavar='YYYYMMDD',
                        help='Skip products whose reference (first) date is '
                        'after this date [no upper limit]')
    parser.add_argument('--reFile', action='store_true', default=False,
                        help='Process all products even if track-{N}/source/ '
                        'symlinks already exist. Use when source/ was created '
                        'by a previous run but the orbit_frame subdirectories '
                        'still need to be built')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Print per-file detail: filenames, orbit info, '
                        'already-filed and unfiled messages. Without this flag '
                        'a progress bar is shown and only totals are reported')
    parser.add_argument('--ignoreMissingROFF', action='store_true', default=False,
                        help='Leave RUNW products in their orbit_frame/H5/ '
                        'directory even when no companion ROFF is present. '
                        'Without this flag, lone RUNWs are moved to '
                        'track-{N}/unfiled/')
    parser.add_argument('--nThreads', type=int, default=8,
                        help='Number of worker processes used to read RUNW '
                        'HDF5 orbit/frame metadata in parallel [8]; use 1 '
                        'for serial processing')

    args = parser.parse_args()
    #
    params = {}
    for param in ['firstOrbit', 'lastOrbit', 'outputPath',
                  'reFile', 'verbose', 'ignoreMissingROFF', 'nThreads']:
        params[param] = getattr(args, param)
        if 'Path' in param and params[param] == '.':
            params[param] = os.getcwd()
    # Resolve inputPath: an explicit CLI argument always wins; otherwise fall
    # back to archiveDir from an autoupdate.yaml sitting in the output tree
    # (the autoupdate workflow runs from the project dir where it lives).
    inputPath = args.inputPath
    if inputPath is None:
        autoYaml = os.path.join(params['outputPath'], 'autoupdate.yaml')
        autoCfg = {}
        if os.path.exists(autoYaml):
            with open(autoYaml) as fp:
                autoCfg = yaml.safe_load(fp) or {}
        inputPath = autoCfg.get('archiveDir')
        if inputPath is None:
            detail = ('present but has no archiveDir key'
                      if os.path.exists(autoYaml) else 'no autoupdate.yaml found')
            parser.error(
                'no input path: pass --inputPath, or set archiveDir in an '
                f'autoupdate.yaml in {params["outputPath"]} ({detail})')
        print(f'Using archiveDir from {autoYaml}: {inputPath}')
    if inputPath == '.':
        inputPath = os.getcwd()
    params['inputPath'] = inputPath
    params['products'] = set(args.products) if args.products else None
    params['firstDate'] = (datetime.strptime(args.firstDate, '%Y%m%d')
                           if args.firstDate else None)
    params['lastDate'] = (datetime.strptime(args.lastDate, '%Y%m%d')
                          if args.lastDate else None)
    return params


def parseFileName(product):
    '''
    Parse a NISAR HDF5 filename into a dictionary of metadata fields.
    Supports RSLC (13 fields) and all other product types (15 fields).
    track, frame, and cycle are returned as str(int(...)) — leading zeros
    stripped. Date fields are returned as datetime objects.
    '''
    h5Name = product.split('/')[-1]
    if 'RSLC' in h5Name:
        pDict = dict(zip(['Sensor', 'Level', 'x', 'productType', 'cycle',
                          'track', 'direction', 'frame', 'bw', 'pol', 'mode',
                          'date1Sstart', 'date1End'],
                         h5Name.split('_')))
    else:
        pDict = dict(zip(['Sensor', 'Level', 'x', 'productType', 'cycle',
                          'track', 'direction', 'frame', '004', 'bw', 'pol',
                          'date1Sstart', 'date1End', 'date2Sstart',
                          'date2End'],
                         h5Name.split('_')))
    for key in pDict:
        if 'date' in key:
            pDict[key] = datetime.strptime(pDict[key], "%Y%m%dT%H%M%S")
        if key in ['track', 'frame', 'cycle']:
            pDict[key] = str(int(pDict[key]))
    return pDict, h5Name


def orbitFromCycleTrack(cycle, track):
    '''
    Return the absolute reference orbit number (orbit1) for a NISAR
    (cycle, track) pass.

    NISAR flies 173 orbits per repeat cycle, so the absolute orbit is a linear
    function of cycle and track: orbit1 = 173*cycle + track + 618. The +618
    constant ties the mission-wide cycle count to the absolute orbit epoch.
    Verified against real filed data (Greenland tracks 1/11/14/42/102 over
    cycles 4-15 and the Antarctic cycle21/track83 -> 4334); a single (cycle,
    track) maps to one orbit even when the pass spans both ascending and
    descending frames across the pole. This reproduces the orbit1 that
    FileNISARProducts reads from the HDF5 and uses for the {orbit1}_{frame}
    directory name, letting callers derive that directory straight from a
    granule filename without opening the HDF5.
    '''
    return 173 * int(cycle) + int(track) + 618


def _coverageSeconds(pDict):
    '''
    Total along-track acquisition span (reference + secondary frame durations,
    in seconds) — a proxy for how much data a product covers. Two products of
    the same pair can be framed to different end times; the larger span is the
    more complete one.
    '''
    return ((pDict['date1End'] - pDict['date1Sstart']).total_seconds() +
            (pDict['date2End'] - pDict['date2Sstart']).total_seconds())


def selectBestRUNW(candidates):
    '''
    Select the best RUNW from a list of candidates for the same track/frame.

    Selection criteria applied in order:
      1. Shortest temporal baseline (|date2Start - date1Start| in days)
      2. Newest file modification time (proxy for processing version)
      3. Longest along-track coverage (later end times = a longer frame) —
         breaks ties between the same pair framed to different end times

    Returns (winner, losers, reasons) where reasons is a dict mapping each
    loser path to a human-readable string explaining why it was not selected.
    '''
    def rankKey(RUNW):
        pDict, _ = parseFileName(RUNW)
        baseline = abs((pDict['date2Sstart'] - pDict['date1Sstart']).days)
        mtime = os.path.getmtime(RUNW)
        # shortest baseline first, newest mtime first, longest coverage first
        return (baseline, -mtime, -_coverageSeconds(pDict))

    ranked = sorted(candidates, key=rankKey)
    winner = ranked[0]
    losers = ranked[1:]

    winner_pDict, _ = parseFileName(winner)
    winner_baseline = abs((winner_pDict['date2Sstart'] -
                           winner_pDict['date1Sstart']).days)
    winner_mtime = os.path.getmtime(winner)
    winner_coverage = _coverageSeconds(winner_pDict)

    reasons = {}
    for loser in losers:
        pDict, _ = parseFileName(loser)
        loser_baseline = abs((pDict['date2Sstart'] - pDict['date1Sstart']).days)
        loser_mtime = os.path.getmtime(loser)
        loser_coverage = _coverageSeconds(pDict)
        if loser_baseline > winner_baseline:
            reasons[loser] = (f'longer temporal baseline '
                              f'({loser_baseline}d vs {winner_baseline}d)')
        elif loser_mtime < winner_mtime:
            reasons[loser] = (
                f'older modification time '
                f'({datetime.fromtimestamp(loser_mtime):%Y-%m-%d %H:%M:%S} '
                f'vs winner {datetime.fromtimestamp(winner_mtime):%Y-%m-%d %H:%M:%S}; '
                f'same temporal baseline of {loser_baseline}d)')
        elif loser_coverage < winner_coverage:
            reasons[loser] = (
                f'shorter along-track coverage ({loser_coverage:.0f}s vs '
                f'winner {winner_coverage:.0f}s; same temporal baseline of '
                f'{loser_baseline}d and modification time)')
        else:
            reasons[loser] = (
                f'identical baseline/mtime/coverage — kept the first by '
                f'filename order')

    return winner, losers, reasons


def getCompanion(RUNW, inputPath, productType):
    '''
    Return the path to a companion ROFF or RIFG product for a given RUNW file,
    or None if the file does not exist.
    Companion files live in inputPath/{productType}/ and have the same name
    as the RUNW with the product-type field substituted.
    '''
    basename = os.path.basename(RUNW).replace('_RUNW_', f'_{productType}_')
    companion = os.path.join(inputPath, productType, basename)
    if os.path.exists(companion):
        return companion
    return None


def findRSLC(RUNW, inputPath):
    '''
    Return a list of RSLC files matching the track, direction, and frame of a
    given RUNW. RSLC has a different (single-date) filename format so cannot
    be derived by simple product-type substitution; instead the track,
    direction, and frame fields are extracted from the RUNW filename (keeping
    the original zero-padded forms) and used to glob the RSLC directory.
    '''
    parts = os.path.basename(RUNW).split('_')
    track, direction, frame = parts[5], parts[6], parts[7]
    return glob.glob(
        f'{inputPath}/RSLC/NISAR_*_{track}_{direction}_{frame}_*.h5')


def isMixedMode(myRUNW):
    '''
    Return True if either the reference or secondary SLC granule name
    contains '_M_', indicating a mixed-mode acquisition.
    '''
    inputs = myRUNW.h5['RUNW']['metadata']['processingInformation']['inputs']
    for key in ['l1ReferenceSlcGranules', 'l1SecondarySlcGranules']:
        granule = inputs[key].asstr()[()].item()
        if '_M_' in granule:
            return True
    return False


def readRUNWMeta(RUNW):
    '''
    Open a RUNW HDF5 and return its orbit/frame metadata as plain values,
    suitable for running in a worker process. Returns a dict with keys
    RUNW, orbit1, orbit2, frame, mixed and error (None on success, else a
    message string with the other fields None). Module-level so it is
    picklable for ProcessPoolExecutor.
    '''
    try:
        myRUNW = nisarhdf.nisarRUNWHDF()
        myRUNW.openHDF(RUNW, noLoadData=True)
        return {'RUNW': RUNW, 'orbit1': myRUNW.referenceOrbit,
                'orbit2': myRUNW.secondaryOrbit, 'frame': myRUNW.frame,
                'mixed': isMixedMode(myRUNW), 'error': None}
    except Exception as e:
        return {'RUNW': RUNW, 'orbit1': None, 'orbit2': None, 'frame': None,
                'mixed': None, 'error': str(e)}


def symlink_file(src_path, dst_path, relative=True, overwrite=False):
    '''
    Create a relative (or absolute) symlink at dst_path pointing to src_path.
    Skips silently if the destination already exists unless overwrite=True.
    Returns True if a new symlink was created, False if skipped.
    '''
    src = Path(src_path)
    dst = Path(dst_path)

    src_abs = src.resolve()
    dst_parent = dst.parent.resolve()

    target = os.path.relpath(src_abs, start=dst_parent) if relative \
        else str(src_abs)

    if overwrite and (dst.exists() or dst.is_symlink()):
        dst.unlink()

    if not (dst.exists() or dst.is_symlink()):
        os.symlink(target, str(dst))
        return True
    return False


def progressBar(i, total):
    '''
    Print an in-place progress bar to stdout.
    Call with i = 0-based index of the current item.
    '''
    bar_width = 40
    filled = int(bar_width * (i + 1) / total)
    bar = '#' * filled + '-' * (bar_width - filled)
    sys.stdout.write(f'\r  [{bar}] {i + 1}/{total}')
    sys.stdout.flush()


# Widest status line written so far, so the next one can pad over any leftover
# characters when it is shorter (e.g. counts shrinking in digit width).
_statusWidth = [0]


def statusLine(msg):
    '''
    Rewrite a single status line in place (carriage return, no scrolling).
    Pads with spaces to clear any longer previous message.
    '''
    text = f'  {msg}'
    pad = max(0, _statusWidth[0] - len(text))
    _statusWidth[0] = len(text)
    sys.stdout.write('\r' + text + ' ' * pad)
    sys.stdout.flush()


def readMetadataParallel(runwList, nThreads, verbose,
                         label='reading orbit/frame metadata'):
    '''
    Read orbit/frame/mixed-mode metadata for a list of RUNW files, in parallel
    across nThreads worker processes (nThreads=1 runs serially, no pool). Each
    HDF5 open is ~0.5-1.5s over NFS, so this is where the wall time goes.
    Returns {RUNW: resultDict}. Shows a single in-place progress line unless
    verbose. Separate processes (not threads) side-step h5py thread-safety.
    '''
    total = len(runwList)
    metaMap = {}
    if not total:
        return metaMap
    nThreads = max(1, nThreads)
    if not verbose:
        statusLine(f'{label} 0/{total}')
    if nThreads > 1:
        ctx = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=nThreads, mp_context=ctx) as ex:
            futures = [ex.submit(readRUNWMeta, r) for r in runwList]
            for k, fut in enumerate(as_completed(futures), 1):
                res = fut.result()
                metaMap[res['RUNW']] = res
                if not verbose:
                    statusLine(f'{label} {k}/{total}')
    else:
        for k, RUNW in enumerate(runwList, 1):
            res = readRUNWMeta(RUNW)
            metaMap[res['RUNW']] = res
            if not verbose:
                statusLine(f'{label} {k}/{total}')
    if not verbose:
        sys.stdout.write('\n')
        sys.stdout.flush()
    return metaMap


def printWarning(msg, verbose, i, total):
    '''
    Print a warning message cleanly. In non-verbose mode the progress bar
    occupies the current line, so a newline is written first; the bar is
    then redrawn on a new line so progress is not lost.
    '''
    if not verbose:
        sys.stdout.write('\n')
    print(f'  Warning: {msg}')
    if not verbose:
        progressBar(i, total)


def moveToUnfiled(srcLink, unfiledDir, logFile, reason, verbose):
    '''
    Move a symlink from its current location to unfiledDir.
    Creates a new symlink in unfiledDir pointing to the same real file,
    removes the original link, and appends a log entry.
    Returns True if the destination symlink was newly created.
    '''
    os.makedirs(unfiledDir, exist_ok=True)
    h5Name = os.path.basename(srcLink)
    destLink = os.path.join(unfiledDir, h5Name)
    created = symlink_file(srcLink, destLink, relative=True, overwrite=False)
    if os.path.islink(srcLink):
        os.unlink(srcLink)
    if created:
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(logFile, 'a') as fp:
            print(f'{ts}  {h5Name}', file=fp)
            print(f'  Reason: {reason}', file=fp)
            print(file=fp)
        if verbose:
            print(f'  Moved to unfiled: {h5Name}')
            print(f'    Reason: {reason}')
    return created


def recoverUnfiled(params, verbose):
    '''
    Re-file RUNWs previously moved to track-N/unfiled/ for a missing companion
    ROFF, once that ROFF has arrived in the archive (e.g. the RUNW downloaded
    before its ROFF). Such an orphan is distinguished from a duplicate loser by
    still having a track-N/source/ symlink — losers never get one. For each
    orphan whose RUNW is still present and whose companion ROFF now exists in
    inputPath/ROFF/, the unfiled/ and source/ links are removed so the main pass
    re-files the now-complete pair. Runs before the skip set is built (dropping
    the source/ link is what lets the main pass re-process the RUNW).
    Returns the number recovered.
    '''
    outputPath = params['outputPath']
    inputPath = params['inputPath']
    nRecovered = 0
    for link in glob.glob(f'{outputPath}/track-*/unfiled/*_RUNW_*.h5'):
        base = os.path.basename(link)
        trackDir = os.path.dirname(os.path.dirname(link))
        sourceLink = os.path.join(trackDir, 'source', base)
        # Duplicate losers have no source/ link — only recover true orphans.
        if not os.path.islink(sourceLink):
            continue
        # The RUNW must still be in the archive and its ROFF must now be too.
        roff = os.path.join(inputPath, 'ROFF', base.replace('_RUNW_', '_ROFF_'))
        if not (os.path.exists(link) and os.path.exists(roff)):
            continue
        os.unlink(link)
        os.unlink(sourceLink)
        nRecovered += 1
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        try:
            with open(os.path.join(trackDir, 'unfiled', 'log'), 'a') as fp:
                print(f'{ts}  {base}', file=fp)
                print('  Recovered: companion ROFF arrived; re-filing', file=fp)
                print(file=fp)
        except OSError:
            pass
        if verbose:
            print(f'  Recovered (ROFF arrived): {base}')
    return nRecovered


def main():
    '''
    Organise a directory of NISAR HDF5 products into the track-based directory
    tree expected by SetupNISAR.

    Pre-processing pass: group RUNWs by (track, frame, reference_date); when
    multiple products share the same reference date AND frame (i.e. they would
    both land in the same orbit1_frame/ directory) select the best (shortest
    temporal baseline, then newest modification time) and file the rest into
    track-{N}/unfiled/ with a log entry explaining the choice.  Products from
    different orbit passes (different reference dates) are never treated as
    duplicates even if they share the same track and frame number.

    Main pass (winners only):
      1. Derive track from filename (no HDF open) for the fast-path skip check.
      2. Create track-{N}/ and track-{N}/source/ if absent.
      3. Symlink RUNW, ROFF, RIFG, and any matching RSLC into source/.
      4. Open the RUNW HDF5 to read referenceOrbit and frame.
      5. Create track-{N}/{orbit1}_{frame}/ and symlink RUNW, ROFF, RIFG there.
    '''
    params = parseArgs()
    verbose = params['verbose']
    #
    RUNWs = sorted(glob.glob(f'{params["inputPath"]}/RUNW/*.h5'))
    print(f'Found {len(RUNWs)} RUNW products')
    #
    if not os.path.exists(params['outputPath']):
        os.mkdir(params['outputPath'])
    #
    # --- Pre-processing: group by (track, frame, refDate), select winners ---
    # Two products are duplicates only if they share the same reference date
    # (i.e. the same orbit1) AND the same frame — meaning they would both land
    # in the same orbit1_frame/ directory.  Products from different orbit passes
    # (different reference dates) are independent and must never be grouped.
    groups = defaultdict(list)
    for RUNW in RUNWs:
        pDict, _ = parseFileName(RUNW)
        refDate = pDict['date1Sstart']
        if params['firstDate'] and refDate.date() < params['firstDate'].date():
            continue
        if params['lastDate'] and refDate.date() > params['lastDate'].date():
            continue
        groups[(pDict['track'], pDict['frame'],
                pDict['date1Sstart'])].append(RUNW)
    #
    # Resolve winners and file losers to unfiled/ using filename + mtime only
    # (no HDF opens here). The winner's absolute orbit — needed only for the
    # unfiled log's "Filed winner" line — is read afterwards, and only for the
    # groups that actually produced a NEW unfiled loser (so a re-run, where all
    # losers are already filed, does zero HDF opens in this pass).
    winners = []
    pendingLogs = []          # (logFile, h5Name, track, frame, refDate, reason, winner)
    winnersNeedingOrbit = []  # winners with >=1 new loser to log
    for (track, frame, _refDate), candidates in sorted(groups.items()):
        if len(candidates) == 1:
            winners.append(candidates[0])
            continue
        #
        winner, losers, reasons = selectBestRUNW(candidates)
        winners.append(winner)
        #
        trackDir = f'{params["outputPath"]}/track-{track}'
        unfiledDir = f'{trackDir}/unfiled'
        os.makedirs(unfiledDir, exist_ok=True)
        logFile = f'{unfiledDir}/log'
        #
        for loser in losers:
            h5Name = os.path.basename(loser)
            destLink = f'{unfiledDir}/{h5Name}'
            isNew = symlink_file(loser, destLink, relative=True, overwrite=False)
            if isNew:
                pendingLogs.append((logFile, h5Name, track, frame, _refDate,
                                    reasons[loser], winner))
                winnersNeedingOrbit.append(winner)
    #
    # Read winner orbits (log detail only) in parallel, once, for the winners
    # that produced new losers.
    winnerMeta = {}
    if winnersNeedingOrbit:
        print(f'Resolving {len(pendingLogs)} duplicate(s)...', flush=True)
        winnerMeta = readMetadataParallel(sorted(set(winnersNeedingOrbit)),
                                          params['nThreads'], verbose,
                                          label='reading duplicate-winner orbits')
    #
    nUnfiled = 0
    for logFile, h5Name, track, frame, refDate, reason, winner in pendingLogs:
        nUnfiled += 1
        wres = winnerMeta.get(winner)
        if wres is not None and wres['error'] is None:
            winnerPath = (f"{wres['orbit1']}_{wres['frame']}/"
                          f"{os.path.basename(winner)}")
        else:
            winnerPath = os.path.basename(winner)
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(logFile, 'a') as fp:
            print(f'{ts}  {h5Name}', file=fp)
            print(f'  Track: {track},  Frame: {frame},  '
                  f'RefDate: {refDate:%Y-%m-%d}', file=fp)
            print(f'  Reason: {reason}', file=fp)
            print(f'  Filed winner: {winnerPath}', file=fp)
            print(file=fp)
        if verbose:
            print(f'  Unfiled {h5Name}')
            print(f'    Reason: {reason}')
    #
    # Preserve original sorted order for winners
    winnersSet = set(winners)
    winners = [r for r in RUNWs if r in winnersSet]
    #
    if nUnfiled:
        print(f'{nUnfiled} duplicate product(s) filed to track-N/unfiled/ '
              f'(see track-N/unfiled/log for details)')
    #
    # --- Main pass: process winners ---
    products = params['products']   # set of requested types, or None = all
    if products is not None and 'RUNW' not in products:
        print(f'  (RUNW used as index only — not filed)')
    # When RUNW is not among the requested products (e.g. --products RSLC) we
    # always proceed so orbit/frame can be read from the HDF5.
    wantRunw = products is None or 'RUNW' in products
    #
    # Recover previously-unfiled orphans whose companion ROFF has since arrived
    # (RUNW filed before its ROFF was downloaded). Must run before the skip set
    # below, since it drops the orphan's source/ link so the main pass re-files
    # it. Only relevant when both RUNW and ROFF are being filed.
    nRecovered = 0
    if products is None or {'RUNW', 'ROFF'}.issubset(products):
        print('Checking unfiled/ for arrived companions...', flush=True)
        nRecovered = recoverUnfiled(params, verbose)
        if nRecovered:
            print(f'  {nRecovered} previously-unpaired product(s) restored '
                  f'(companion arrived) — re-filing', flush=True)
    #
    # Smart skip: scan the already-filed RUNW source symlinks in one glob per
    # track and test membership in memory, instead of a per-winner
    # os.path.exists() network stat (thousands of NFS round trips on a re-run).
    filedRUNW = set()
    if wantRunw and not params['reFile']:
        print('Scanning already-filed products...', flush=True)
        for src in glob.glob(
                f'{params["outputPath"]}/track-*/source/*_RUNW_*.h5'):
            filedRUNW.add(os.path.basename(src))
        print(f'  {len(filedRUNW)} product(s) already filed', flush=True)
    #
    skipped = 0
    toProcess = []
    for RUNW in winners:
        if wantRunw and not params['reFile'] \
                and os.path.basename(RUNW) in filedRUNW:
            if verbose:
                print(f'  Already filed: {os.path.basename(RUNW)}')
            skipped += 1
            continue
        toProcess.append(RUNW)
    #
    # Read each new RUNW's orbit/frame/mixed-mode from its HDF5 in parallel —
    # this network-bound open (~0.5-1.5s each) is the dominant cost. The
    # filesystem mutations (mkdir/symlink) stay serial below to avoid races.
    total = len(toProcess)
    metaMap = readMetadataParallel(toProcess, params['nThreads'], verbose)
    #
    filed_rslcs = set()   # track RSLCs filed via RUNW companion search
    nFiled = 0
    nMixed = 0
    nRange = 0
    nError = 0
    for i, RUNW in enumerate(toProcess):
        h5Name = os.path.basename(RUNW)
        pDict, _ = parseFileName(RUNW)
        track = pDict['track']
        trackDir = f'{params["outputPath"]}/track-{track}'
        sourceDir = f'{trackDir}/source'
        #
        if verbose:
            print(RUNW)
        else:
            statusLine(f'{nFiled} filed, {skipped} skipped (already exist), '
                       f'{total - i} to process')
        #
        # Steps 2-3: create track/source dirs and symlink requested source
        # product types. Done before the orbit-range/mixed-mode checks so a
        # frame skipped there still gets its source/ links (original order).
        os.makedirs(sourceDir, exist_ok=True)
        sourceProducts = []
        if products is None or 'RUNW' in products:
            sourceProducts.append(RUNW)
        for pt in ['ROFF', 'RIFG']:
            if products is None or pt in products:
                c = getCompanion(RUNW, params['inputPath'], pt)
                if c:
                    sourceProducts.append(c)
        if products is None or 'RSLC' in products:
            rslcs = findRSLC(RUNW, params['inputPath'])
            filed_rslcs.update(rslcs)
            sourceProducts.extend(rslcs)
        #
        for product in sourceProducts:
            symlink_file(product,
                         f'{sourceDir}/{os.path.basename(product)}',
                         relative=True, overwrite=False)
        #
        # Steps 4-5: use the pre-read orbit/frame, create orbit_frame dir/links.
        # Non-fatal outcomes are tallied (counters below) and reported in the
        # final summary rather than printed inline, so the status line above
        # stays a single non-scrolling line; --verbose still prints each one.
        res = metaMap[RUNW]
        if res['error'] is not None:
            nError += 1
            if verbose:
                printWarning(f'could not read {h5Name}: {res["error"]} — '
                             f'skipping', verbose, i, total)
            continue
        orbit1 = res['orbit1']
        orbit2 = res['orbit2']
        frame = res['frame']
        #
        if orbit1 < params['firstOrbit'] or orbit1 > params['lastOrbit']:
            nRange += 1
            if verbose:
                printWarning(f'orbit {orbit1} outside --firstOrbit/--lastOrbit '
                             f'range — skipping', verbose, i, total)
            continue
        if res['mixed']:
            nMixed += 1
            if verbose:
                printWarning(f'mixed-mode frame {orbit1}_{frame} — skipping',
                             verbose, i, total)
            continue
        #
        if verbose:
            print(f'  orbit1={orbit1}  orbit2={orbit2}  frame={frame}')
        orbitFrameDir = f'{trackDir}/{orbit1}_{frame}'
        h5Dir = f'{orbitFrameDir}/H5'
        os.makedirs(h5Dir, exist_ok=True)
        #
        # Symlink requested L2 products (RUNW, ROFF, RIFG) into H5/.
        # RSLC is source-only (wrapH5WithVRT does not handle it).
        for pt, product in [('RUNW', RUNW)] + [
                (pt, getCompanion(RUNW, params['inputPath'], pt))
                for pt in ['ROFF', 'RIFG']]:
            if product is not None and (products is None or pt in products):
                symlink_file(product,
                             f'{h5Dir}/{os.path.basename(product)}',
                             relative=True, overwrite=False)
        nFiled += 1
    #
    # Clear the status line and print a single final summary.
    if not verbose and total:
        statusLine(f'{nFiled} filed, {skipped} skipped (already exist), '
                   f'0 to process')
        sys.stdout.write('\n')
        sys.stdout.flush()
    summary = f'Filed {nFiled} product(s)'
    extras = []
    if nRecovered:
        extras.append(f'{nRecovered} previously-unpaired restored')
    if skipped:
        extras.append(f'{skipped} skipped (already exist)')
    if nMixed:
        extras.append(f'{nMixed} mixed-mode')
    if nRange:
        extras.append(f'{nRange} outside orbit range')
    if nError:
        extras.append(f'{nError} unreadable')
    if extras:
        summary += '; ' + ', '.join(extras)
    print(summary)
    #
    # --- Standalone RSLC pass: file RSLCs with no companion RUNW ---
    # findRSLC only runs inside the RUNW loop, so RSLCs without a matching RUNW
    # are never seen above.  Scan the RSLC directory directly and file any that
    # were not already picked up.  These go to track-N/source/ only — no
    # orbit_frame/H5/ because there is no L2 product to process.
    if products is None or 'RSLC' in products:
        all_rslcs = sorted(glob.glob(f'{params["inputPath"]}/RSLC/*.h5'))
        orphans = [r for r in all_rslcs if r not in filed_rslcs]
        if orphans:
            print(f'Filing {len(orphans)} RSLC product(s) with no companion RUNW...')
            for j, rslc in enumerate(orphans):
                if not verbose:
                    progressBar(j, len(orphans))
                pDict, _ = parseFileName(rslc)
                refDate = pDict['date1Sstart']
                if params['firstDate'] and refDate.date() < params['firstDate'].date():
                    continue
                if params['lastDate'] and refDate.date() > params['lastDate'].date():
                    continue
                track = pDict['track']
                trackDir = f'{params["outputPath"]}/track-{track}'
                sourceDir = f'{trackDir}/source'
                if not os.path.exists(trackDir):
                    os.mkdir(trackDir)
                if not os.path.exists(sourceDir):
                    os.mkdir(sourceDir)
                if verbose:
                    print(f'  {rslc}')
                symlink_file(rslc, f'{sourceDir}/{os.path.basename(rslc)}',
                             relative=True, overwrite=params['reFile'])
            if not verbose:
                print()
    #
    # --- Companion check: every orbit_frame dir must have both RUNW and ROFF ---
    # Skip when the user explicitly requested only a subset that does not
    # include both RUNW and ROFF — in that case incomplete pairs are expected.
    if products is not None and not {'RUNW', 'ROFF'}.issubset(products):
        return
    # If only one is present the product cannot be processed by SetupNISAR;
    # move the orphan to track-N/unfiled/ and log why.
    print('Checking RUNW/ROFF companions across the tree...', flush=True)
    nMoved = 0
    for trackEntry in sorted(glob.glob(f'{params["outputPath"]}/track-*')):
        if not os.path.isdir(trackEntry):
            continue
        unfiledDir = f'{trackEntry}/unfiled'
        logFile = f'{unfiledDir}/log'
        for entry in sorted(os.listdir(trackEntry)):
            # orbit_frame directories are named {digits}_{digits}
            parts = entry.split('_')
            if not (len(parts) == 2 and all(p.isdigit() for p in parts)):
                continue
            orbitFrameDir = f'{trackEntry}/{entry}'
            if not os.path.isdir(orbitFrameDir):
                continue
            # One directory read (scandir) instead of two globs of H5/ —
            # halves NFS round trips across the whole output tree.
            runwLinks, roffLinks = [], []
            try:
                with os.scandir(f'{orbitFrameDir}/H5') as it:
                    for e in it:
                        if not e.name.endswith('.h5'):
                            continue
                        if '_RUNW_' in e.name:
                            runwLinks.append(e.path)
                        elif '_ROFF_' in e.name:
                            roffLinks.append(e.path)
            except FileNotFoundError:
                pass
            if runwLinks and not roffLinks and not params['ignoreMissingROFF']:
                for link in runwLinks:
                    reason = (f'missing companion ROFF file '
                              f'(moved from {entry}/)')
                    if moveToUnfiled(link, unfiledDir, logFile,
                                     reason, verbose):
                        nMoved += 1
            elif roffLinks and not runwLinks:
                for link in roffLinks:
                    reason = (f'missing companion RUNW file '
                              f'(moved from {entry}/)')
                    if moveToUnfiled(link, unfiledDir, logFile,
                                     reason, verbose):
                        nMoved += 1
            # Remove the orbit_frame dir if it is now empty
            try:
                if not os.listdir(orbitFrameDir):
                    os.rmdir(orbitFrameDir)
            except OSError:
                pass
    if nMoved:
        print(f'{nMoved} product(s) moved to track-N/unfiled/ '
              f'(missing RUNW or ROFF companion; '
              f'see track-N/unfiled/log for details)')


if __name__ == "__main__":
    main()
