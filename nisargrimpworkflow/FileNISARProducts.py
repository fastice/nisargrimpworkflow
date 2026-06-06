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
from collections import defaultdict
from pathlib import Path
from datetime import datetime

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
    parser.add_argument('inputPath', type=str, nargs=1,
                        help='Root directory with RUNW/, ROFF/, RIFG/, RSLC/ '
                        'product subdirectories')
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

    args = parser.parse_args()
    #
    params = {}
    for param in ['firstOrbit', 'lastOrbit', 'outputPath', 'inputPath',
                  'reFile', 'verbose', 'ignoreMissingROFF']:
        params[param] = getattr(args, param)
        if 'Path' in param:
            if params[param] == '.':
                params[param] = os.getcwd()
    params['inputPath'] = params['inputPath'][0]
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


def selectBestRUNW(candidates):
    '''
    Select the best RUNW from a list of candidates for the same track/frame.

    Selection criteria applied in order:
      1. Shortest temporal baseline (|date2Start - date1Start| in days)
      2. Newest file modification time (proxy for processing version)

    Returns (winner, losers, reasons) where reasons is a dict mapping each
    loser path to a human-readable string explaining why it was not selected.
    '''
    def rankKey(RUNW):
        pDict, _ = parseFileName(RUNW)
        baseline = abs((pDict['date2Sstart'] - pDict['date1Sstart']).days)
        mtime = os.path.getmtime(RUNW)
        return (baseline, -mtime)   # shortest baseline first, newest mtime first

    ranked = sorted(candidates, key=rankKey)
    winner = ranked[0]
    losers = ranked[1:]

    winner_pDict, _ = parseFileName(winner)
    winner_baseline = abs((winner_pDict['date2Sstart'] -
                           winner_pDict['date1Sstart']).days)
    winner_mtime = datetime.fromtimestamp(os.path.getmtime(winner))

    reasons = {}
    for loser in losers:
        pDict, _ = parseFileName(loser)
        loser_baseline = abs((pDict['date2Sstart'] - pDict['date1Sstart']).days)
        if loser_baseline > winner_baseline:
            reasons[loser] = (f'longer temporal baseline '
                              f'({loser_baseline}d vs {winner_baseline}d)')
        else:
            loser_mtime = datetime.fromtimestamp(os.path.getmtime(loser))
            reasons[loser] = (
                f'older modification time ({loser_mtime:%Y-%m-%d %H:%M:%S} '
                f'vs winner {winner_mtime:%Y-%m-%d %H:%M:%S}; '
                f'same temporal baseline of {loser_baseline}d)')

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
    winners = []
    nUnfiled = 0
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
        # Open the winner HDF once to get orbit1/frame for the log message.
        # Duplicates are rare so the extra open is acceptable.
        myWinner = nisarhdf.nisarRUNWHDF()
        myWinner.openHDF(winner, noLoadData=True)
        winnerOrbitFrame = f'{myWinner.referenceOrbit}_{myWinner.frame}'
        winnerPath = f'{winnerOrbitFrame}/{os.path.basename(winner)}'
        #
        for loser in losers:
            h5Name = os.path.basename(loser)
            destLink = f'{unfiledDir}/{h5Name}'
            isNew = symlink_file(loser, destLink, relative=True, overwrite=False)
            if isNew:
                nUnfiled += 1
                ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                with open(logFile, 'a') as fp:
                    print(f'{ts}  {h5Name}', file=fp)
                    print(f'  Track: {track},  Frame: {frame},  '
                          f'RefDate: {_refDate:%Y-%m-%d}', file=fp)
                    print(f'  Reason: {reasons[loser]}', file=fp)
                    print(f'  Filed winner: {winnerPath}', file=fp)
                    print(file=fp)
                if verbose:
                    print(f'  Unfiled {h5Name}')
                    print(f'    Reason: {reasons[loser]}')
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
    total = len(winners)
    skipped = 0
    filed_rslcs = set()   # track RSLCs filed via RUNW companion search
    for i, RUNW in enumerate(winners):
        h5Name = os.path.basename(RUNW)
        #
        # Derive track from filename (fast; no HDF open) for the skip check
        pDict, _ = parseFileName(RUNW)
        track = pDict['track']
        trackDir = f'{params["outputPath"]}/track-{track}'
        sourceDir = f'{trackDir}/source'
        #
        # Fast-path: if the RUNW source symlink already exists and RUNW is one
        # of the requested products, skip unless --reFile is set.
        # When RUNW is not in the requested products (e.g. --products RSLC) we
        # always proceed so orbit/frame can be read from the HDF5.
        wantRunw = products is None or 'RUNW' in products
        if wantRunw and not params['reFile']:
            if os.path.exists(f'{sourceDir}/{h5Name}'):
                if verbose:
                    print(f'  Already filed: {sourceDir}/{h5Name}')
                skipped += 1
                if not verbose:
                    progressBar(i, total)
                continue
        #
        if verbose:
            print(RUNW)
        else:
            progressBar(i, total)
        #
        # Step 2: create track dir
        if not os.path.exists(trackDir):
            os.mkdir(trackDir)
        #
        # Step 3: create source dir and symlink requested product types
        if not os.path.exists(sourceDir):
            os.mkdir(sourceDir)
        #
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
        # Step 4-5: read orbit/frame from HDF, create orbit_frame dir and links
        myRUNW = nisarhdf.nisarRUNWHDF()
        myRUNW.openHDF(RUNW, noLoadData=True)
        orbit1 = myRUNW.referenceOrbit
        orbit2 = myRUNW.secondaryOrbit
        frame = myRUNW.frame
        #
        if orbit1 < params['firstOrbit'] or orbit1 > params['lastOrbit']:
            printWarning(f'orbit {orbit1} outside --firstOrbit/--lastOrbit '
                         f'range — skipping', verbose, i, total)
            continue
        if isMixedMode(myRUNW):
            printWarning(f'mixed-mode frame {orbit1}_{frame} — skipping',
                         verbose, i, total)
            continue
        #
        if verbose:
            print(f'  orbit1={orbit1}  orbit2={orbit2}  frame={frame}')
        orbitFrameDir = f'{trackDir}/{orbit1}_{frame}'
        if not os.path.exists(orbitFrameDir):
            os.mkdir(orbitFrameDir)
        h5Dir = f'{orbitFrameDir}/H5'
        if not os.path.exists(h5Dir):
            os.mkdir(h5Dir)
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
    #
    if not verbose:
        print()  # newline after final progress bar
    if skipped:
        print(f'Skipped {skipped} already-filed RUNW product(s) '
              f'(use --reFile to override, --verbose to list them)')
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
            runwLinks = glob.glob(f'{orbitFrameDir}/H5/*_RUNW_*.h5')
            roffLinks = glob.glob(f'{orbitFrameDir}/H5/*_ROFF_*.h5')
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
