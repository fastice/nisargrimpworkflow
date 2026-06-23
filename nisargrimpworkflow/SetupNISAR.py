#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 12:07:17 2026

@author: ian
"""
import utilities as u
import argparse
import nisarhdf
import time
from nisargrimpworkflow.wrapH5WithVRT import wrapH5sInFrameDir
from nisargrimpworkflow.processFrameIonosphere import processFrameIonosphere
import geojson
import copy
import glob
import numpy as np
from subprocess import run, DEVNULL
import os
import shutil
import yaml
from scipy.interpolate import interp1d
from scipy.ndimage import binary_erosion
import sys
import importlib.resources
from osgeo import gdal, osr
import sarfunc
from collections import Counter

def restrictedFrame(x):
    x = int(x)
    if x < 0.0 or x > 999:
        raise argparse.ArgumentTypeError("Value must be between 0 and 999")
    return x


def parseCommandLine():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mConvert NISAR Level-2 HDF5 products (RUNW, '
        'ROFF) into the binary flat-file and VRT formats used by the GrIMP '
        'velocity mosaic pipeline.\033[0m\n\n'
        'Processes all frames of a single orbit pair and consolidates the '
        'per-frame outputs into a single virtual-frame directory using GDAL '
        'VRT mosaics.  Must be run from the track-<N> directory that contains '
        'the <orbit1>_<frame> subdirectories (e.g. 12345_010/, 12345_020/, ...).',
        epilog='Examples:\n'
               '  # Process all frames for orbit 12345\n'
               '  SetupNISAR 12345\n\n'
               '  # Process a subset of frames\n'
               '  SetupNISAR 12345 --firstFrame 10 --lastFrame 30\n\n'
               '  # Reprocess everything from scratch\n'
               '  SetupNISAR 12345 --overWrite\n\n'
               '  # Reprocess phase products only (keep existing ROFF)\n'
               '  SetupNISAR 12345 --overWritePhase\n\n'
               '  # Phase products only, custom virtual frame name\n'
               '  SetupNISAR 12345 --RUNWOnly --virtualFrame 0001\n\n'
               '  # Include mixed-mode frames, print all subprocess output\n'
               '  SetupNISAR 12345 --allowMixedMode --verbose\n'
               '\nPart of the nisargrimpworkflow package.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('orbit1', type=int,
                        help='Reference orbit number')

    parser.add_argument('--virtualFrame', type=str, default=None,
                        help='Frame suffix for the consolidated virtual-frame '
                        'output directory (e.g. 0000 → <orbit1>_0000/). '
                        'When omitted, the number is assigned automatically '
                        'based on frame content and reference orbits.')
    parser.add_argument('--overWrite', action="store_true",
                        help='Re-run RUNW and ROFF conversion even if output '
                        'products already exist')
    parser.add_argument('--overWritePhase', action="store_true",
                        help='Re-run RUNW (phase) conversion only; leave '
                        'existing ROFF products untouched')
    parser.add_argument('--firstFrame', type=restrictedFrame, default=0,
                        help='Skip frames numbered below this value (0–999) [0]')
    parser.add_argument('--lastFrame', type=restrictedFrame, default=999,
                        help='Skip frames numbered above this value (0–999) [999]')
    parser.add_argument('--allowMixedMode', action="store_true",
                        help='Include mixed-mode frames (SLC granule name '
                        'contains _M_); these are skipped by default')
    parser.add_argument('--RUNWOnly', action="store_true",
                        help='Process only RUNW (phase/coherence/ionosphere) '
                        'products; skip ROFF offset conversion and power images')
    parser.add_argument('--noMask', action='store_true', default=False,
                        help='Do not apply the fast-region mask to layer 3 '
                        'during ROFF offset conversion')
    parser.add_argument('--verbose', action='store_true',
                        help='Print all subprocess output to the terminal '
                        '(default: suppressed to keep output readable)')
    parser.add_argument('--ompThreads', type=int, default=6,
                        help='Number of OpenMP threads passed to siminsar '
                        'via ROFFtoGrimp/RUNWtoGrimp [6]')
    parser.add_argument('--phaseDerivedIonosphere', action='store_true',
                        help='Use the RUNW phase-derived ionosphere screen '
                        '(passes --phaseDerivedIonosphere to RUNWtoGrimp). '
                        'When omitted, estimateIonosphere is run after ROFF '
                        'to estimate the ionosphere from range offsets.')
    parser.add_argument('--outputAll', action='store_true',
                        help='Pass --outputAll to estimateIonosphere: write all '
                        'intermediate bands (5-band VRT + stem-named TIFs). '
                        'Default: write only correctedUnwrappedPhase.tif, '
                        'ionosphereCorrection.tif, ionosphereCorrection.offset.tif '
                        'each with a single-band VRT sidecar.')
    parser.add_argument('--phaseThresh', type=float,
                        default=14 * 3.141592653589793, metavar='RAD',
                        help='Pass --phaseThresh to estimateIonosphere: mask '
                        'correctedUnwrappedPhase where '
                        '|correctedPhase - simPhase| >= RAD radians. '
                        'Screens regions with likely incorrect unwrapping. '
                        'Default: 14π rad.')
    parser.add_argument('--noPhaseThreshPass', action='store_true', default=False,
                        help='Pass --noPhaseThreshPass to estimateIonosphere: '
                        'disable the second-pass iono re-estimation that uses '
                        'the phase-residual mask instead of the velocity mask.')
    parser.add_argument('--sepIceRock', action='store_true', default=False,
                        help='Pass --sepIceRock to estimateIonosphere — restrict the '
                        'per-frame (phase+offset_phase)/2 estimate to ice pixels. '
                        'In globalFillIonosphere(), the excluded rock pixels '
                        '(offset_phase only, an absolute "actual vs simulated '
                        'offset" reference) are used once, at the virtual-frame '
                        'level, to anchor the otherwise-floating ice-derived '
                        'ionosphere field via a single additive constant. '
                        'Enables global fill on its own (combine with '
                        '--noGlobalFillIono to disable the rock-anchoring step). '
                        'Mask auto-detected from offsetSims/offsets.geom.mask.vrt.')
    parser.add_argument('--sigmaAz', type=float, default=None, metavar='PX',
                        help='Azimuth Gaussian sigma for iono smoothing, in RUNW pixels '
                        '(passed to estimateIonosphere --sigma-az). Default: 10 px.')
    parser.add_argument('--sigmaRg', type=float, default=None, metavar='PX',
                        help='Range Gaussian sigma for iono smoothing, in RUNW pixels '
                        '(passed to estimateIonosphere --sigma-rg). Default: 30 px.')
    parser.add_argument('--correlationOnly', action='store_true',
                        help='Extract coherence and geodat files only. '
                        'Skips ROFF conversion, ionosphere estimation, and '
                        'virtual-frame assembly.')
    parser.add_argument('--corrOnly', action='store_true',
                        help='Extract .cor files per real frame and assemble '
                        'only the correlation virtual frame. Skips ROFF '
                        'conversion, ionosphere estimation, and power images.')
    parser.add_argument('--geodatsOnly', action='store_true',
                        help='Re-merge virtual-frame geodats from existing '
                        'per-frame geodats and the existing virtual-frame VRT '
                        '(e.g. to pick up a mergedGeodat() fix) without '
                        'rebuilding any VRTs or re-running RUNW/ROFF/ionosphere '
                        'processing. Per-frame products must already exist.')
    parser.add_argument('--noGlobalFillIono', action='store_true',
                        help='Disable full-swath ionosphere gap fill; use '
                        'per-frame fill only (default: global fill is on)')
    parser.add_argument('--retainIntermediateIono', action='store_true',
                        help='Keep per-frame unfilled and per-frame filled '
                        'offset iono files after global fill (useful for '
                        'debugging or comparison)')
    parser.add_argument('--debugIono', action='store_true',
                        help='Rename all ionosphere intermediates with "debug" '
                        'in their filenames and build two extra virtual-frame '
                        'VRTs: (1) assembled unfilled iono across all frames, '
                        '(2) assembled per-frame-filled (pre-global-fill) offset '
                        'correction.  Files named *.debug.* can be found and '
                        'deleted as a group after inspection.')
    parser.add_argument('--clean', action='store_true',
                        help='Remove all computed output files for this orbit -- '
                        'everything --overWrite would replace (workingDir/, offsetSims/, '
                        'simPhase/, range/azimuth offsets, phase/ionosphere products, '
                        'pairinfo, geodats, *.vrt.stats) -- plus mislabeled files (frame '
                        'directories whose name does not match the expected real-frame '
                        'or virtual-frame pattern, e.g. from a past failed run) and '
                        'stale temp files left by an interrupted run (*.tmp*, *.vsmooth). '
                        'debug/ directories are emptied (not removed; see --cleanDebug). '
                        'Builds and prints the full list, then prompts for confirmation '
                        'unless --noPrompt is set. Exits immediately afterward without '
                        'any further processing.')
    parser.add_argument('--cleanDebug', action='store_true',
                        help='Empty the contents of every debug/ directory for this '
                        'orbit (real-frame and virtual-frame), leaving the empty '
                        'directory itself in place. Same print/prompt-unless-noPrompt '
                        'flow as --clean, and also exits immediately afterward. Can be '
                        'combined with --clean (redundant but harmless, since --clean '
                        'already empties debug/ too).')
    parser.add_argument('-noPrompt', '--noPrompt', action='store_true',
                        help='Skip the confirmation prompt for --clean/--cleanDebug')
    args = parser.parse_args()
    #
    params = {}
    for key in ['overWrite', 'overWritePhase', 'firstFrame', 'lastFrame',
                'orbit1', 'allowMixedMode', 'virtualFrame', 'noMask',
                'verbose', 'RUNWOnly', 'ompThreads', 'phaseDerivedIonosphere',
                'outputAll', 'phaseThresh', 'noPhaseThreshPass', 'sepIceRock',
                'sigmaAz', 'sigmaRg',
                'correlationOnly', 'corrOnly', 'geodatsOnly',
                'retainIntermediateIono', 'debugIono',
                'clean', 'cleanDebug', 'noPrompt']:
        params[key] = getattr(args, key)
    params['globalFillIono'] = not args.noGlobalFillIono and (
        bool(args.phaseDerivedIonosphere) or bool(args.debugIono) or bool(args.sepIceRock))
    #
    if args.verbose:
        params['stdout'], params['stderr'] = None, None
    else:
        params['stdout'], params['stderr'] = DEVNULL, DEVNULL
    #
    return params


def splitFrameGroups(frames):
    '''
    Split a sorted list of frame numbers into contiguous groups.
    A group is a run of consecutive integers (no gaps).
    Returns a list of lists, e.g. [61,62,63,71,72,73] -> [[61,62,63],[71,72,73]].
    '''
    if not frames:
        return []
    groups = []
    currentGroup = [frames[0]]
    for f in frames[1:]:
        if f == currentGroup[-1] + 1:
            currentGroup.append(f)
        else:
            groups.append(currentGroup)
            currentGroup = [f]
    groups.append(currentGroup)
    return groups


def splitGroupsBySecondaryEpoch(groups, frameSecondaryInfo):
    '''
    Further split each frame-number-contiguous group wherever a frame's actual
    secondary date (from frameSecondaryInfo, see getSecondaryOrbit()) disagrees
    with the group's majority secondary date.

    A contiguous run of frame numbers can still span two different secondary
    acquisitions: if the intended-epoch secondary product was missing
    upstream for just one frame, NISAR processing substitutes a different
    repeat cycle for that frame alone, while neighboring frames get the
    intended epoch. Merging such a frame into the rest of the virtual frame
    bakes a real, physical baseline discontinuity into the merged geodat/
    state vectors with no warning -- see [[project_track131_rparams_sigma]].
    Splitting it into its own group lets assignVirtualFrameNumbers() give it
    its own virtual frame instead.
    '''
    splitGroups = []
    for group in groups:
        dates = [frameSecondaryInfo[f][1] for f in group if f in frameSecondaryInfo]
        if not dates:
            splitGroups.append(group)
            continue
        majorityDate = Counter(dates).most_common(1)[0][0]
        subGroup = []
        for f in group:
            info = frameSecondaryInfo.get(f)
            if info is not None and info[1] != majorityDate:
                u.mywarning(
                    f'Frame {f}: secondary date {info[1]} does not match the '
                    f'majority secondary date {majorityDate} for its '
                    f'frame-contiguous group -- splitting it into its own '
                    f'virtual frame instead of merging it in (likely a '
                    f'missing-acquisition fallback to a different repeat '
                    f'cycle for this one frame).')
                if subGroup:
                    splitGroups.append(subGroup)
                    subGroup = []
                splitGroups.append([f])
            else:
                subGroup.append(f)
        if subGroup:
            splitGroups.append(subGroup)
    return splitGroups


def writeFramesList(frameDir, frames):
    '''Write the frame numbers that make up this virtual frame to frames.txt
    so future orbits can use it as a reference for splitting.'''
    with open(f'{frameDir}/frames.txt', 'w') as fp:
        print(' '.join(str(f) for f in sorted(frames)), file=fp)


def readReferenceFrameSets(orbit1, cwd='.'):
    '''
    Scan cwd for *_0???/frames.txt files from other orbits.
    Returns {vfStr: set_of_frame_numbers} (union across all reference orbits).
    '''
    refSets = {}
    for framesFile in sorted(glob.glob(os.path.join(cwd, '*_0???', 'frames.txt'))):
        dirName = os.path.basename(os.path.dirname(framesFile))
        parts = dirName.rsplit('_', 1)
        if len(parts) != 2:
            continue
        try:
            refOrbit = int(parts[0])
            vf = parts[1]
        except ValueError:
            continue
        if refOrbit == orbit1:
            continue
        with open(framesFile) as fp:
            frames = {int(x) for x in fp.read().split() if x.strip().isdigit()}
        if frames:
            refSets.setdefault(vf, set()).update(frames)
    return refSets


def assignVirtualFrameNumbers(groups, orbit1, cwd='.'):
    '''
    Assign a 4-digit virtual frame number to each naturally-detected group.

    Numbering scheme:
      N*10        – group N, canonical (new area, or matches/extends a reference group)
      N*10 + K    – group N, Kth fragment (proper subset of reference group N, K≥1)

    Rules applied per group:
      • Overlaps a reference VF and frames ⊇ reference frames → canonical N*10
      • Overlaps a reference VF and frames ⊂ reference frames → fragment N*10+K
      • No reference overlap → new canonical group M*10, M above all reference groups

    Returns list of (groupFrames, virtualFrameStr) pairs.
    '''
    refSets = readReferenceFrameSets(orbit1, cwd)
    # Canonical numbers (multiples of 10) that already appear in references
    refCanonicals = {(int(vf) // 10) for vf in refSets}
    nextNewGroup = max(refCanonicals) + 1 if refCanonicals else 0

    # Read existing virtual frames for THIS orbit so re-runs reuse the same numbers.
    selfSets = {}
    for framesFile in sorted(glob.glob(os.path.join(cwd, f'{orbit1}_0???', 'frames.txt'))):
        dirName = os.path.basename(os.path.dirname(framesFile))
        vf = dirName.rsplit('_', 1)[1]
        with open(framesFile) as fp:
            existingFrames = {int(x) for x in fp.read().split() if x.strip().isdigit()}
        if existingFrames:
            selfSets[vf] = existingFrames

    assignedNums = set()   # tracks numbers chosen during this run

    assignments = []
    for group in groups:
        groupSet = set(group)

        # On re-runs, prefer reusing this orbit's own existing VF numbers over
        # cross-orbit fragment assignment.  Pick the existing VF with the most
        # frame overlap (must not already be taken by an earlier group this run).
        bestSelfVF, bestSelfOverlap = None, 0
        for vf, existingFrames in selfSets.items():
            overlap = len(groupSet & existingFrames)
            if overlap > bestSelfOverlap and int(vf) not in assignedNums:
                bestSelfOverlap, bestSelfVF = overlap, vf
        if bestSelfVF is not None:
            assignedNums.add(int(bestSelfVF))
            assignments.append((group, bestSelfVF))
            continue

        # Find the reference VF with the greatest frame overlap
        bestRefVF, bestOverlap = None, 0
        for vf, refFrames in refSets.items():
            overlap = len(groupSet & refFrames)
            if overlap > bestOverlap:
                bestOverlap, bestRefVF = overlap, vf

        if bestRefVF is None:
            # No reference match → new geographic group
            while nextNewGroup * 10 in assignedNums:
                nextNewGroup += 1
            vfNum = nextNewGroup * 10
            nextNewGroup += 1
        else:
            baseGroup = int(bestRefVF) // 10
            refFrames = refSets[bestRefVF]
            if groupSet >= refFrames:
                # Matches or extends reference → canonical
                vfNum = baseGroup * 10
            else:
                # Proper subset → fragment
                k = 1
                while baseGroup * 10 + k in assignedNums:
                    k += 1
                vfNum = baseGroup * 10 + k

        assignedNums.add(vfNum)
        assignments.append((group, f'{vfNum:04d}'))

    if assignments:
        print('Virtual frame assignments: '
              + ', '.join(f'{vf}={g}' for g, vf in assignments))
    return assignments


def getFrames(myArgs):
    '''
        Get a sort list of the frames to process
    '''
    dirs = glob.glob(f'{myArgs["orbit1"]}_??') + \
             glob.glob(f'{myArgs["orbit1"]}_1??') + \
             glob.glob(f'{myArgs["orbit1"]}_2??')    
    frames = [int(x.split('_')[-1])
              for x in sorted(dirs)]
    return [x for x in frames
            if myArgs['firstFrame'] <= x <= myArgs['lastFrame']]


# Per-frame-directory computed outputs that --clean/--cleanDebug operate on.
_CLEAN_DIR_PATTERNS = ['workingDir', 'offsetSims', 'simPhase']
_CLEAN_FILE_PATTERNS = [
    'range.offsets*', 'azimuth.offsets*', 'offsets.range-azimuth.vrt',
    '*.nisar.cor*', '*.correctedUnwrappedPhase*', '*.ionosphereCorrection*',
    '*.unwrappedPhase.tif', '*.offsetPhase.tif', '*.uw*', '*.ion*',
    '*.pairinfo', 'geodat*.geojson', '*.vrt.stats',
    # Virtual-frame-level merges of per-real-frame simPhase/offsetSims outputs, written
    # directly at the virtual-frame directory root by custom_buildvrtWithOffsets (not
    # inside a subdirectory, unlike the real-frame versions already caught above via
    # the offsetSims/simPhase rmtree).
    'velSim*', 'maskVel*', 'offsets.geom*', 'offsets.velocity*',
    'frames.txt',
]
_STALE_PATTERNS = ['*.tmp*', '*interp.tmp*', '*.vsmooth']
_CLEAN_BUCKET_LABELS = {'normal': 'normally processed', 'mislabeled': 'mislabeled',
                        'stale': 'stale', 'debug': 'debug'}


def cleanFrames(myArgs, debugOnly=False):
    '''
    Build (but do not remove) the list of computed-output files/directories for this
    orbit, categorized into buckets for --clean / --cleanDebug.

    Frame directories are discovered by globbing {orbit1}_* directly rather than via
    getFrames(), so directories with unexpected/bad frame numbers (e.g. left over from a
    failed run) are caught too -- that's the whole point of --clean. "Well-formed" means
    matching one of the three real-frame patterns getFrames() itself uses (_?? / _1?? /
    _2??) or the four-digit leading-zero virtual-frame pattern used elsewhere in this
    file and in readReferenceFrameSets() (_0???); anything else is "mislabeled".

    Parameters
    ----------
    myArgs : dict
        Must contain 'orbit1'.
    debugOnly : bool, optional
        If True (set when --cleanDebug is given without --clean), only the 'debug'
        bucket is populated -- nothing else is touched. Default False.

    Returns
    -------
    dict
        {'normal': [...], 'mislabeled': [...], 'stale': [...], 'debug': [...]}
        Directories appear as single entries (removed wholesale); debug/ contents are
        listed individually (debug/ itself is never included -- it's emptied, not removed).
    '''
    orbit1 = myArgs['orbit1']
    wellFormedDirs = set(
        glob.glob(f'{orbit1}_??') + glob.glob(f'{orbit1}_1??') +
        glob.glob(f'{orbit1}_2??') + glob.glob(f'{orbit1}_0???'))
    allDirs = sorted(glob.glob(f'{orbit1}_*'))

    buckets = {'normal': [], 'mislabeled': [], 'stale': [], 'debug': []}

    for d in allDirs:
        bucket = 'normal' if d in wellFormedDirs else 'mislabeled'

        debugDir = os.path.join(d, 'debug')
        if os.path.isdir(debugDir):
            buckets['debug'] += sorted(glob.glob(os.path.join(debugDir, '*')))

        if debugOnly:
            continue

        staleMatches = set()
        for pat in _STALE_PATTERNS:
            staleMatches.update(glob.glob(os.path.join(d, pat)))
            staleMatches.update(glob.glob(os.path.join(d, '*', pat)))
        buckets['stale'] += sorted(staleMatches)

        for sub in _CLEAN_DIR_PATTERNS:
            p = os.path.join(d, sub)
            if os.path.isdir(p):
                buckets[bucket].append(p)

        for pat in _CLEAN_FILE_PATTERNS:
            for f in glob.glob(os.path.join(d, pat)):
                if f not in staleMatches:
                    buckets[bucket].append(f)

    return buckets


def _removePath(path):
    if os.path.islink(path) or os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        shutil.rmtree(path)


def confirmAndRemove(buckets, noPrompt):
    '''
    Print the counts (by category) built by cleanFrames(), prompt for confirmation
    (unless noPrompt) with [y/n/l] -- 'l' prints the full categorized list and
    re-prompts -- then remove everything on 'y'. Aborts via u.myerror on 'n'.
    '''
    order = ['normal', 'mislabeled', 'stale', 'debug']
    allItems = []
    for key in order:
        allItems += buckets.get(key, [])

    if not allItems:
        print('\nNothing to clean.')
        return

    counts = [(key, len(buckets[key])) for key in order if buckets.get(key)]
    parts = [f'{n} {_CLEAN_BUCKET_LABELS[key]} files' for key, n in counts]
    msg = 'Removing ' + ', '.join(parts[:-1]) + (
        f', and {parts[-1]}' if len(parts) > 1 else parts[0])

    def printFullList():
        for key in order:
            items = buckets.get(key, [])
            if not items:
                continue
            print(f'\n{_CLEAN_BUCKET_LABELS[key].capitalize()} ({len(items)}):')
            for item in items:
                print(f'  {item}')

    if not noPrompt:
        while True:
            ans = input(f'\n\033[1m{msg} [y/n/l]\033[0m\n')
            if ans.lower() == 'y':
                break
            if ans.lower() == 'n':
                u.myerror('User prompted abort')
            if ans.lower() == 'l':
                printFullList()

    for item in allItems:
        _removePath(item)
    print(f'\nRemoved {len(allItems)} items.')


def getSecondaryOrbit(myArgs):
    '''
    Get first and second orbits.

    Also records each frame's own secondary date AND its own bandwidth/looks
    (read directly from that frame's own RUNW/ROFF product) in
    myArgs['frameSecondaryInfo'] as
    {frame: (secondaryOrbit, secondaryDate, NumberRangeLooks, NumberAzimuthLooks, bandwidth)}.
    A frame-number-contiguous run can still span two different secondary
    acquisitions if the intended-epoch secondary product was missing
    upstream for just one frame and NISAR processing substituted a different
    repeat cycle for it (seen on track-131/3517 frame 38, track-126/3685
    frame 51, track-16/3748 frame 35, track-58/3444 frame 44 -- secondaryOrbit
    was unchanged but secondaryDate jumped a full repeat cycle).
    splitGroupsBySecondaryEpoch() uses the secondaryDate part of this map to
    keep such a frame out of the virtual frame it would otherwise be merged
    into. The same kind of run can also mix a frame acquired at a different
    bandwidth (e.g. 40 vs 77 MHz, hence a different NumberRangeLooks) than its
    neighbors -- the per-group refresh in main() uses the looks/bandwidth part
    of this map so a split-out group gets its own correct pairinfo/sensor YAML
    instead of inheriting whichever frame getSecondaryOrbit() found first.
    '''
    myArgs['frameSecondaryInfo'] = {}
    found = False
    for frame in myArgs['frames']:
        for product in ['RUNW', 'ROFF']:
            frameDir = f'{myArgs["orbit1"]}_{frame}'
            files = glob.glob(f'{frameDir}/H5/NISAR*{product}*.h5')
            if len(files) < 1:
                continue
            myProd = getattr(nisarhdf, f'nisar{product}HDF', None)()
            myProd.openHDF(files[0])
            if myProd.referenceOrbit == myArgs['orbit1']:
                myProd.getRangeBandWidth()
                frameBandwidth = myProd.rangeBandwidth/1e6
                myArgs['frameSecondaryInfo'][frame] = \
                    (myProd.secondaryOrbit, myProd.secondaryDate,
                     myProd.NumberRangeLooks, myProd.NumberAzimuthLooks,
                     frameBandwidth)
                if not found:
                    myArgs['orbit2'] = myProd.secondaryOrbit
                    myArgs['bandwidth'] = frameBandwidth
                    myArgs['NumberAzimuthLooks'] = myProd.NumberAzimuthLooks
                    myArgs['NumberRangeLooks'] = myProd.NumberRangeLooks
                    myArgs['myProd'] = myProd
                    myArgs['datetime'] = myProd.Date
                    myArgs['secondaryDateTime'] = myProd.secondaryDate
                    found = True
                break
    return found


def copy_sensor_yaml(myArgs, dest_dir):
    """Copy sensor YAML into `dest_dir` based on `myArgs['bandwidth']` if present.

    Only copies when the destination file does not already exist.
    """
    if 'bandwidth' not in myArgs or myArgs.get('bandwidth') is None:
        return
    
    try:
        bw = int(myArgs.get('bandwidth'))
    except Exception:
        u.mywarning(f'Invalid bandwidth value: {myArgs.get("bandwidth")}')
        return
    sensor_map = {77: '80', 40: '40', 20: '20'}
    if bw not in sensor_map:
        u.mywarning(f'Unsupported bandwidth {bw}; skipping sensor YAML copy')
        return
    sensor = sensor_map[bw]
    sensors_ref = importlib.resources.files('sarfunc').joinpath('sensors')
    yamlfile = f'NISAR{sensor}.yaml'
    src = str(sensors_ref.joinpath(yamlfile))
    dst = os.path.join(dest_dir, f'sensor.{yamlfile}')
    # Exactly one sensor.*.yaml is expected per directory -- tieScript.py's
    # _get_base_nlooks() trusts whichever one it finds. Remove any other-
    # bandwidth sensor YAML left over from a previous run (e.g. this frame's
    # real bandwidth was redetected, or different-bandwidth image data was
    # later incorporated into this same directory) so it can't be picked up
    # by mistake.
    for stale in glob.glob(os.path.join(dest_dir, 'sensor.NISAR*.yaml')):
        if stale != dst:
            os.remove(stale)
            print(f'Removed stale sensor YAML {stale}')
    copied = False
    if os.path.exists(dst):
        print(f'Sensor YAML already exists at {dst}; skipping copy')
    else:
        try:
            shutil.copyfile(src, dst)
            copied = True
            print(f'Copied sensor YAML {src} -> {dst}')
        except FileNotFoundError:
            u.mywarning(f'Sensor YAML not found: {src}')
        except Exception as e:
            u.mywarning(f'Failed copying sensor YAML: {e}')

    # If the destination YAML exists, update intLooksR/intLooksA from myArgs
    if os.path.exists(dst):
        # Only these keys are expected in myArgs for looks
        range_keys = ['NumberRangeLooks']
        az_keys = ['NumberAzimuthLooks']
        range_val = None
        az_val = None
        for k in range_keys:
            if k in myArgs and myArgs.get(k) is not None:
                try:
                    range_val = int(myArgs.get(k))
                except Exception:
                    range_val = None
                break
        for k in az_keys:
            if k in myArgs and myArgs.get(k) is not None:
                try:
                    az_val = int(myArgs.get(k))
                except Exception:
                    az_val = None
                break

        if range_val is None and az_val is None:
            return

        try:
            with open(dst, 'r') as fp:
                y = yaml.safe_load(fp) or {}
        except Exception as e:
            u.mywarning(f'Unable to read YAML {dst}: {e}')
            return

        modified = False
        if range_val is not None:
            if y.get('intLooksR') != range_val:
                y['intLooksR'] = range_val
                modified = True
        if az_val is not None:
            if y.get('intLooksA') != az_val:
                y['intLooksA'] = az_val
                modified = True

        if modified:
            try:
                with open(dst, 'w') as fp:
                    yaml.safe_dump(y, fp, default_flow_style=False)
                print(f'Updated looks in {dst}: intLooksR={range_val}, intLooksA={az_val}')
            except Exception as e:
                u.mywarning(f'Failed writing YAML {dst}: {e}')


def getProductH5(product,  frame, myArgs):
    '''
    Get product for frame
    '''
    inputDir = myArgs['inputDir']
    productDir = f'{inputDir}/{myArgs["orbit1"]}_{frame}/H5/NISAR*{product}*h5'
    myProduct = glob.glob(productDir)
    if len(myProduct) == 0:
        u.mywarning(f'There are no {product} products in {productDir}')
        return None
    if len(myProduct) > 1:
        u.mywarning(
            f'Warning more than one {product} products in {productDir}')
    return myProduct[0]


def variableSmoothingArgs(myArgs):
    '''
    Build the --minTol/--percentSpeed/--maxTol/--maxSmoothRadius/--smoothNIter/
    --noVariableSmoothing CLI args to pass through to ROFFtoGrimp/RUNWtoGrimp, from values
    already resolved (project.yaml or CLI) into myArgs by main(). Returns [] when the
    feature isn't configured (minTol is None) -- same opt-in convention as ROFFtoGrimp/
    RUNWtoGrimp themselves.
    '''
    if myArgs.get('minTol') is None:
        return []
    args = ['--minTol', str(myArgs['minTol']),
           '--percentSpeed', str(myArgs['percentSpeed']),
           '--maxTol', str(myArgs['maxTol']),
           '--maxSmoothRadius', str(myArgs['maxSmoothRadius']),
           '--smoothNIter', str(myArgs['smoothNIter'])]
    if myArgs.get('noVariableSmoothing'):
        args += ['--noVariableSmoothing']
    return args


def processFrameROFF(frame, myArgs):
    '''
    Run programs to unpack RUNW and ROFF files into GrIMP formats.

    Parameters
    ----------
    RUNW : str
        RUNW hdf path.

    Returns
    -------
    None.

    '''
    ROFFFile = getProductH5('ROFF', frame, myArgs)
    if ROFFFile is None:
        return
    # Use the frame dir (not the H5 subdir) so that ROFFtoGrimp's intermediate
    # dirs (workingDir/, offsetSims/), final binaries, and VRTs all land
    # alongside the RUNW products and geodat files.
    orbit1 = myArgs['orbit1']
    outputDir = f'{myArgs["outputDir"]}/{orbit1}_{frame}'
    #
    rangeOffsets = glob.glob(f'{outputDir}/range.offsets')
    #
    if len(rangeOffsets) == 0 or myArgs['overWrite']:
        # Setup command; pass geodats explicitly so ROFFtoGrimp need not
        # search for them in outputDir (where *.nisar.uw files are absent).
        command = ['ROFFtoGrimp',
                   '--outputDir', outputDir,
                   '--geodat1', myArgs['geodat1'][-1],
                   '--geodat2', myArgs['geodat2'][-1],
                   '--ompThreads', str(myArgs['ompThreads']),
                   ROFFFile]
        if myArgs.get('regionFile'):
            command += ['--regionFile', myArgs['regionFile']]
        if myArgs.get('verticalCorrection'):
            command += ['--verticalCorrection', myArgs['verticalCorrection']]
        command += variableSmoothingArgs(myArgs)
        if myArgs['noMask'] is True:
            command += ['--noMask']
        if myArgs.get('debugIono'):
            command += ['--debugIono']
        if myArgs['verbose'] is True:
            command += ['--verbose'] #+['--mergeOnly']
        # Call command
        run(command,  stderr=myArgs['stderr'], stdout=myArgs['stdout'])


def processFrameRUNW(frame, myArgs):
    '''
    Run programs to unpack RUNW and ROFF files into GrIMP formats.

    Parameters
    ----------
    RUNW : str
        RUNW hdf path.

    Returns
    -------
    None.

    '''
    orbit1, orbit2 = myArgs['orbit1'],  myArgs['orbit2']
    RUNWFile = getProductH5('RUNW', frame, myArgs)
    if RUNWFile is None:
        return True
    myRUNW = nisarhdf.nisarRUNWHDF()
    myRUNW.openHDF(RUNWFile)
    nLooksR = myRUNW.NumberRangeLooks
    nLooksA = myRUNW.NumberAzimuthLooks
    # Geodats are written by RUNWtoGrimp to the RUNW output dir (frame dir),
    # not to the H5 subdir where the HDF5 lives.
    ruNWOutputDir = f'{myArgs["outputDir"]}/{orbit1}_{frame}'
    #
    if not myArgs['allowMixedMode']:
        if isMixedMode(myRUNW):
            print('Skipping mixed mode {RUNWFile}')
            return True
    #
    outputDir = f'{myArgs["outputDir"]}/{orbit1}_{frame}'

    geodat = f'{outputDir}/geodat{nLooksR}x{nLooksA}.geojson'
    if myArgs.get('phaseDerivedIonosphere'):
        # Full phase+ionosphere run: expect both VRT outputs and the geodat
        files = glob.glob(
            f'{outputDir}/{orbit1}_{frame}.{orbit2}_{frame}.*x*.vrt')
        missing = not os.path.exists(geodat) or not all(
            any(key in f for f in files) for key in {'ion.filt', 'uw.interp'})
    else:
        # Default (--noPhase --noIon): geodat is always written
        missing = not os.path.exists(geodat)
    if missing or myArgs['overWrite'] or myArgs['overWritePhase']:
        command = ["RUNWtoGrimp",
                   "--frame", str(frame),
                   "--referenceOrbit", str(orbit1),
                   "--secondaryOrbit", str(orbit2),
                   "--ompThreads", str(myArgs['ompThreads']),
                   "--outputDir", f"./{orbit1}_{frame}",
                   RUNWFile]
        if myArgs.get('regionFile'):
            command += ['--regionFile', myArgs['regionFile']]
        if myArgs.get('verticalCorrection'):
            command += ['--verticalCorrection', myArgs['verticalCorrection']]
        command += variableSmoothingArgs(myArgs)
        if myArgs['verbose'] is True:
            command += ['--verbose']
        if myArgs.get('phaseDerivedIonosphere'):
            command += ['--phaseDerivedIonosphere']
            # Phase output is only produced in this branch (the other branch passes
            # --noPhase), so only here is there anything for the variable smoothing-radius
            # map to smooth -- and only here does RUNWtoGrimp need --simPhase to actually
            # run siminsar and produce phaseSim.smr.vrt for it.
            if myArgs.get('minTol') is not None and not myArgs.get('noVariableSmoothing'):
                command += ['--simPhase']
        else:
            command += ['--noPhase', '--noIon']
        # print(' '.join(command))
        run(command, stderr=myArgs['stderr'], stdout=myArgs['stdout'])
    else:
        print(f'skipping {orbit1}_{frame} since products exist')
    # RUNWtoGrimp can exit early (e.g. unresolvable region/epsg) without
    # producing the geodat; only record this frame as usable — and let
    # ROFF/ionosphere run for it — if it actually wrote one.
    if not os.path.exists(geodat):
        print(f'RUNWtoGrimp produced no geodat for {orbit1}_{frame} — skipping frame')
        return True
    myArgs['geodat1'].append(
        f'{ruNWOutputDir}/geodat{nLooksR}x{nLooksA}.geojson')
    myArgs['geodat2'].append(
        f'{ruNWOutputDir}/geodat{nLooksR}x{nLooksA}.secondary.geojson')
    return False


def isMixedMode(myRUNW):
    '''
    Check for mixed mode
    '''
    inputs = myRUNW.h5['RUNW']['metadata']['processingInformation']['inputs']

    for key in ['l1ReferenceSlcGranules', 'l1SecondarySlcGranules']:
        granule = inputs[key].asstr()[()].item()
        if "_M_" in granule:
            return True
    return False


def mergeCorners(geodatFirst, geodatLast):
    '''
    Merge corners first and last frames for a merged product
    '''
    cornersFirst = getPolygon(geodatFirst)
    cornersLast = getPolygon(geodatLast)
    mergedCorners = copy.deepcopy(cornersFirst)
    passType = geodatFirst['properties']['PassType'].lower().strip()
    if passType == 'ascending':
        mergedCorners['ur'] = cornersLast['ur']
        mergedCorners['ul'] = cornersLast['ul']
    elif passType == 'descending':
        mergedCorners['lr'] = cornersLast['lr']
        mergedCorners['ll'] = cornersLast['ll']
    else:
        raise ValueError(f"Unknown PassType: {passType}")
    return [mergedCorners[x] for x in ['ll', 'lr', 'ur', 'ul', 'll']]


def sortAndUniqueSV(t, x):
    '''
    Sort and remove duplicates.

    Raises if two inputs share the same timestamp but disagree in value --
    genuine duplicates (the same orbit slice reused across overlapping
    frames) are bit-for-bit identical, so any real disagreement means the
    inputs actually come from different secondary acquisitions and should
    not have been merged in the first place (same failure mode as the
    StateVectorInterval check below). See [[project_track131_rparams_sigma]].
    '''
    t = np.array(t, dtype=float)
    x = np.array(x, dtype=float)
    idx = np.argsort(t, kind='stable')
    t = t[idx]
    x = x[idx]
    unique_t, unique_idx, counts = np.unique(
        t, return_index=True, return_counts=True)
    for ut, ui, count in zip(unique_t, unique_idx, counts):
        if count > 1 and not np.allclose(x[ui:ui + count], x[ui], atol=1e-3):
            raise ValueError(
                f'mergeStateVectors: state vectors at time {ut} disagree by '
                f'more than 1e-3 across input frames -- inputs likely come '
                f'from different secondary acquisitions and should not be '
                f'merged.')
    t = unique_t
    x = x[unique_idx]
    return t, x


def mergeStateVectors(geos, geoMerged):
    '''
    Merge state vectors
    '''
    pos, vel, tPos, tVel = [], [], [], []
    dTs = []
    for geo in geos:
        # print(geo)
        props = geos[geo]['properties']
        tS = float(props['TimeOfFirstStateVector'])
        dT = float(props['StateVectorInterval'])
        dTs.append(dT)
        for key in props:
            if 'SV_Pos' in key:
                sn = int(key.split('_')[-1]) - 1
                tPos.append(tS + sn * dT)
                pos.append(props[key])
            elif 'SV_Vel' in key:
                sn = int(key.split('_')[-1]) - 1
                tVel.append(tS + sn * dT)
                vel.append(props[key])
    if not np.allclose(dTs, dTs[0]):
        raise ValueError("StateVectorInterval values are inconsistent")
    dT = float(dTs[0])
    # Sort and eliminate duplicates
    tPos, pos = sortAndUniqueSV(tPos, pos)
    tVel, vel = sortAndUniqueSV(tVel, vel)
    tNew = np.arange(np.min(tPos), np.max(tPos) + dT/2, dT)
    # Interpolate state vectors to new positions
    kind = 'cubic' if len(tVel) >= 4 else 'linear'
    newPos = interp1d(tPos, pos, axis=0, kind=kind)(tNew)
    newVel = interp1d(tVel, vel, axis=0, kind=kind)(tNew)
    # Remove all trace of old state vectors
    keys = list(geoMerged['properties'].keys())
    toRemove = ('NumberOfStateVectors',
                'TimeOfFirstStateVector',
                'StateVectorInterval')
    for key in keys:
        if key.startswith('SV_') or key in toRemove:
            geoMerged['properties'].pop(key, None)
    # Now add the new state vector information
    geoMerged['properties']['NumberOfStateVectors'] = len(tNew)
    geoMerged['properties']['TimeOfFirstStateVector'] = float(tNew[0])
    geoMerged['properties']['StateVectorInterval'] = float(dT)
    #
    for i, p, v in zip(range(1, len(tNew)+1), newPos, newVel):
        geoMerged['properties'][f'SV_Pos_{i}'] = list(p)
        geoMerged['properties'][f'SV_Vel_{i}'] = list(v)


def getPolygon(geodat):
    '''
    Get the polygon as dictionary with corners keyed by location
    '''
    c = [x for x in geodat['geometry']['coordinates'][0]]
    return dict(zip(['ll', 'ul', 'ur', 'lr'], c))


def mergedGeodat(geodatFiles, vrtFile, secondary=False):
    # Save first and last frame
    firstFrame = sorted(list(geodatFiles))[0].split('/')[-2].split('_')[-1]
    #
    # Get merged image size
    image = nisarhdf.readVrtAsXarray(vrtFile)
    band_name = list(image.data_vars)[0]
    na, nr = image[band_name].shape
    #
    # fTime, lTime = image.y.min().item(), image.y.max().item()
    # print('Time', fTime, lTime)
    #
    # Read VRT geotransform to derive accurate range fields
    vrtDs = gdal.Open(vrtFile)
    gt = vrtDs.GetGeoTransform()
    vrtDs = None
    x0, pixelWidth = gt[0], gt[1]
    mlNearRange = x0 + 0.5 * pixelWidth
    mlFarRange = x0 + (nr - 0.5) * pixelWidth
    mlCenterRange = x0 + (nr / 2.0) * pixelWidth
    #
    # Load geodats
    geodats = {}
    lastFrame = None
    for geoFile in geodatFiles:
        orbit, frame = geoFile.split('/')[-2].split('_')
        with open(geoFile) as fp:
            geodats[frame] = geojson.load(fp)
        lastFrame = frame
    # Create the merged geodat as a copy of the first
    geodatMerged = copy.deepcopy(geodats[f'{firstFrame}'])
    # Merge the state vectors
    mergeStateVectors(geodats, geodatMerged)
    # Merge the corner coords
    geodatMerged['geometry']['coordinates'][0] = \
        mergeCorners(geodats[f'{firstFrame}'], geodats[f'{lastFrame}'])
    if not secondary:
        # The VRT's geotransform/dimensions reflect the reference/primary
        # grid, so deriving size and near/far/center range from it is
        # correct here -- this also gives one consistent value when source
        # frames' own values differ slightly.
        geodatMerged['properties']['MLAzimuthSize'] = na
        geodatMerged['properties']['MLRangeSize'] = nr
        geodatMerged['properties']['MLNearRange'] = mlNearRange
        geodatMerged['properties']['MLFarRange'] = mlFarRange
        geodatMerged['properties']['MLCenterRange'] = mlCenterRange
    else:
        # MLAzimuthSize (like MLRangeSize) is a co-registration grid
        # dimension, not an independent physical quantity of the secondary's
        # own acquisition -- RUNW/ROFF always resample the secondary onto the
        # primary's output grid, so for any single real frame primary and
        # secondary report identical MLAzimuthSize (verified directly,
        # track-30/2205_34: both 1045). The merged secondary's MLAzimuthSize
        # must therefore be the same merged total (na) the primary branch
        # uses, not just firstFrame's own individual size -- otherwise
        # azparams's svInitAzParams, which projects a synthetic grid from the
        # primary into the secondary's bounds using MLAzimuthSize, ends up
        # using one frame's size instead of the merged total and finds almost
        # nothing inside those bounds regardless of real data coverage.
        geodatMerged['properties']['MLAzimuthSize'] = na
        # MLNearRange/MLFarRange genuinely are independent physical
        # quantities (slant range to first sample depends on orbit/baseline),
        # so they correctly come from the secondary's own geodats, not the
        # (primary-aligned) output grid -- the deep copy above from
        # geodats[firstFrame] already has the right, mutually-consistent
        # values for them and MLRangeSize -- just check the other frames
        # agree on MLNearRange/MLRangeSize.
        firstNearRange = geodats[f'{firstFrame}']['properties']['MLNearRange']
        firstRangeSize = geodats[f'{firstFrame}']['properties']['MLRangeSize']
        for frame, geo in geodats.items():
            frameNearRange = geo['properties']['MLNearRange']
            frameRangeSize = geo['properties']['MLRangeSize']
            if abs(frameNearRange - firstNearRange) > 1.0:
                u.mywarning(
                    f'mergedGeodat: secondary MLNearRange for frame {frame} '
                    f'({frameNearRange}) differs from frame {firstFrame} '
                    f'({firstNearRange}) by '
                    f'{frameNearRange - firstNearRange:.3f} m -- using '
                    f'frame {firstFrame}\'s value for the merged secondary '
                    f'geodat')
            if frameRangeSize != firstRangeSize:
                u.mywarning(
                    f'mergedGeodat: secondary MLRangeSize for frame {frame} '
                    f'({frameRangeSize}) differs from frame {firstFrame} '
                    f'({firstRangeSize}) -- using frame {firstFrame}\'s '
                    f'value for the merged secondary geodat')
    # Return result
    geojsonString = nisarhdf.formatGeojson(geojson.dumps(geodatMerged))
    # Remove existing file to avoid problems with links
    geodatFileName = os.path.join(os.path.dirname(vrtFile),
                                  os.path.basename(geodatFiles[0]))
    # Save the file
    with open(geodatFileName, 'w') as fpGeojson:
        print(geojsonString, file=fpGeojson)


def writePairInfo(myArgs):
    ''' Write a pairinfo file in the virtualFrame dir '''
    virtualFrame = myArgs['virtualFrame']
    frameDir = f'{myArgs["outputDir"]}/{myArgs["orbit1"]}_{virtualFrame}'
    pairFile = f'{frameDir}/{myArgs["orbit1"]}.{myArgs["orbit2"]}.pairinfo'
    # orbit2 can differ between reruns for a given virtualFrame now that a group's
    # epoch is re-derived per group (see splitGroupsBySecondaryEpoch) -- remove any
    # stale pairinfo left from a previous orbit2 so it doesn't linger alongside the
    # current, correct one.
    for stalePairFile in glob.glob(f'{frameDir}/{myArgs["orbit1"]}.*.pairinfo'):
        if stalePairFile != pairFile:
            os.remove(stalePairFile)
    with open(pairFile, 'w') as fp:
        print(f'{myArgs["orbit1"]}  {myArgs["orbit2"]}  {myArgs["datetime"]}  '
              f'{myArgs["secondaryDateTime"]}  '
              f'{myArgs["NumberRangeLooks"]}  {myArgs["NumberAzimuthLooks"]}',
              file=fp)


def createVirtualFrameRUNW(myArgs):
    '''
    Create a virtual frames using vrts to consolodate frame products

    Parameters
    ----------
    myArgs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    orbit = myArgs['orbit1']
    virtualFrame = myArgs["virtualFrame"]
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
    orbit = myArgs['orbit1']
    #
    if myArgs.get('geodatsOnly'):
        # Re-merge geodats only, using the existing virtual-frame VRT --
        # skip rebuilding any VRTs or touching RUNW/ROFF/ionosphere products.
        existingVrts = glob.glob(f'{frameDir}/*.correctedUnwrappedPhase.vrt')
        if not existingVrts:
            u.mywarning(f'createVirtualFrameRUNW: --geodatsOnly requested but no '
                       f'existing correctedUnwrappedPhase VRT found in {frameDir} '
                       f'-- skipping')
            return
        mergedGeodat(myArgs['geodat1'], existingVrts[0])
        mergedGeodat(myArgs['geodat2'], existingVrts[0], secondary=True)
        return
    #
    virtualVRTs = {}

    # Build velSim and maskVel reference mosaics from per-frame simPhase outputs.
    # These anchor the correctedUnwrappedPhase DC biases to the simulated velocity phase.
    velSimFiles = [f for frame in myArgs['frames']
                   for f in glob.glob(f'{orbit}_{frame}/simPhase/velSim.vrt')]
    maskVelFiles = [f for frame in myArgs['frames']
                    for f in glob.glob(f'{orbit}_{frame}/simPhase/maskVel.vrt')]
    velSimVrt  = f'{frameDir}/velSim.vrt'
    maskVelVrt = f'{frameDir}/maskVel.vrt'
    if velSimFiles:
        print(f'Building {velSimVrt}')
        run(['custom_buildvrtWithOffsets', '--overWrite', velSimVrt] + velSimFiles)
    if maskVelFiles:
        print(f'Building {maskVelVrt}')
        run(['custom_buildvrtWithOffsets', '--overWrite', maskVelVrt] + maskVelFiles)

    # --debugIono: mosaic the per-frame unsmoothed-phase comparison copies (written by
    # estimateIonosphere.py's variable-smoothing step) and the velSim.smr radius map
    # itself, so smoothed vs unsmoothed phase can be compared at the virtual-frame level.
    if myArgs.get('debugIono'):
        debugFrameDir = f'{frameDir}/debug'
        os.makedirs(debugFrameDir, exist_ok=True)
        unsmoothedPhaseFiles = [
            f for frame in myArgs['frames']
            for f in glob.glob(f'{orbit}_{frame}/debug/*.correctedUnwrappedPhase.unsmoothed.vrt')]
        if unsmoothedPhaseFiles:
            unsmoothedPhaseVrt = f'{debugFrameDir}/correctedUnwrappedPhase.unsmoothed.vrt'
            # Same --offsets/--referencePhase/--mask as the real correctedUnwrappedPhase
            # merge above -- without them the per-frame DC/ambiguity bias is never
            # reconciled, so adjacent frames show seams that have nothing to do with
            # smoothing and would make the comparison misleading.
            command = ['custom_buildvrtWithOffsets', '--offsets']
            if os.path.exists(velSimVrt):
                command += ['--referencePhase', velSimVrt]
            if os.path.exists(maskVelVrt):
                command += ['--mask', maskVelVrt]
            command += ['--overWrite', unsmoothedPhaseVrt] + unsmoothedPhaseFiles
            print(f'Building {unsmoothedPhaseVrt}')
            run(command)
        radiusFiles = [f for frame in myArgs['frames']
                      for f in glob.glob(f'{orbit}_{frame}/simPhase/velSim.smr.vrt')]
        if radiusFiles:
            radiusVrt = f'{debugFrameDir}/velSim.smr.vrt'
            print(f'Building {radiusVrt}')
            run(['custom_buildvrtWithOffsets', '--overWrite', radiusVrt] + radiusFiles)

    # (glob suffix, use --offsets, label, save as ionosphereRangeOffsetCorrection)
    products = [
        ('*.correctedUnwrappedPhase.vrt',    True,  'correctedUnwrappedPhase',    False),
        ('*.cor.vrt',                         False, 'cor',                         False),
        ('*.ionosphereCorrection.vrt',        True,  'ionosphereCorrection',        False),
        ('*.ionosphereCorrection.offset.vrt', True,  'ionosphereCorrection.offset', True),
    ]
    if myArgs.get('globalFillIono'):
        # Global fill replaces both the per-frame-filled phase and offset VRTs;
        # skip assembling them here and add the sparse unfilled versions so
        # globalFillIonosphere can do a single full-swath fill pass on each.
        products = [(g, o, l, r) for g, o, l, r in products
                    if l not in ('ionosphereCorrection', 'ionosphereCorrection.offset')]
        products.append(
            ('*.ionosphereCorrectionUnfilled.vrt', True,
             'ionosphereCorrectionUnfilled', False))
        products.append(
            ('*.ionosphereCorrectionUnfilled.offset.vrt', True,
             'ionosphereCorrectionUnfilled.offset', False))
    for globPattern, useOffsets, label, isIonoRangeOffset in products:
        # Find files
        myFiles = []
        for frame in myArgs['frames']:
            myFile = glob.glob(f'{orbit}_{frame}/{globPattern}')
            if len(myFile) == 1:
                myFiles += myFile
            else:
                u.mywarning(f'more than one or missing vrt {myFile}'
                            f'\n\tfor {orbit}_{frame}/{globPattern}')
        if not myFiles:
            u.mywarning(f'createVirtualFrameRUNW: no input VRTs found for '
                        f'{label}, skipping')
            continue
        # Create virtual raster
        command = ['custom_buildvrtWithOffsets']
        if useOffsets:
            command.append('--offsets')
        if label == 'correctedUnwrappedPhase':
            if os.path.exists(velSimVrt):
                command += ['--referencePhase', velSimVrt]
            if os.path.exists(maskVelVrt):
                command += ['--mask', maskVelVrt]
        # Virtual product name — derive from the actual source frame in myFiles[0]
        # (frames[0] may have no output if its processing failed, making myFiles[0]
        # belong to a later frame whose number won't match the replace pattern)
        _srcFrame = os.path.dirname(myFiles[0]).split('_')[-1]
        virtualProduct = os.path.basename(
            myFiles[0].replace(f'_{_srcFrame}.', f'_{virtualFrame}.'))
        virtualVRTs[label] = f'{frameDir}/{virtualProduct}'
        #
        # The ionosphere correction on the offset grid goes into range-offset metadata
        # so the geocoder knows to apply the correction.
        if isIonoRangeOffset:
            myArgs['ionosphereRangeOffsetCorrection'] = virtualProduct
        failFile = f'{virtualVRTs[label]}.fail'
        if os.path.exists(failFile):
            os.remove(failFile)
        # Force a new final file
        command.append('--overWrite')
        command.append(f'{frameDir}/{virtualProduct}')
        command += myFiles
        print(f'Building {frameDir}/{virtualProduct}')
        run(command)
        if not os.path.exists(virtualVRTs[label]):
            open(failFile, 'w').close()
    # Create fixed-name symlink so tieScript process_phase_yaml() can locate
    # the unwrapped phase regardless of NISAR product naming conventions.
    if 'correctedUnwrappedPhase' in virtualVRTs:
        phase_link = f'{frameDir}/phase.uw.vrt'
        target = os.path.basename(virtualVRTs['correctedUnwrappedPhase'])
        if os.path.lexists(phase_link):
            os.remove(phase_link)
        os.symlink(target, phase_link)
    # Now create geodats
    mergedGeodat(myArgs['geodat1'], virtualVRTs['correctedUnwrappedPhase'])
    mergedGeodat(myArgs['geodat2'], virtualVRTs['correctedUnwrappedPhase'],
                secondary=True)


def createVirtualFrameROFF(myArgs):
    '''
    Create a virtual frames using vrts to consolodate frame products

    Parameters
    ----------
    myArgs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    orbit = myArgs['orbit1']
    virtualFrame = myArgs["virtualFrame"]
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
    #
    files = ["azimuth.offsets.vrt",
             "offsetSims/offsets.geom.ll.vrt",
             "offsetSims/offsets.geom.mask.vrt",
             "offsetSims/offsets.geom.vrt",
             "offsetSims/offsets.velocity.ll.vrt",
             "offsetSims/offsets.velocity.mask.vrt",
             "offsets.range-azimuth.vrt",
             "offsetSims/offsets.velocity.vrt",
             "range.offsets.vrt"]
    #
    virtualVRTs = {}
    #
    # --debugIono: mosaic the per-frame unsmoothed-offsets comparison copies (written by
    # ROFFtoGrimp.py's variable-smoothing step) and the offsets.velocity.smr radius map
    # itself, so smoothed vs unsmoothed offsets can be compared at the virtual-frame level.
    if myArgs.get('debugIono'):
        debugFrameDir = f'{frameDir}/debug'
        os.makedirs(debugFrameDir, exist_ok=True)
        # Unsmoothed copies are genuine offset values, so use --offsets (cross-frame
        # radiometric continuity correction) the same way the real range/azimuth.offsets
        # do. The radius map is a per-pixel pixel-count, not a continuous physical
        # quantity -- --offsets's continuity correction would corrupt it, so omit it.
        debugFiles = [
            ('range.offsets.unsmoothed.vrt', 'debug/range.offsets.unsmoothed.vrt', True),
            ('azimuth.offsets.unsmoothed.vrt', 'debug/azimuth.offsets.unsmoothed.vrt', True),
            ('offsets.velocity.smr.vrt', 'offsetSims/offsets.velocity.smr.tif', False),
        ]
        for outName, globPattern, useOffsets in debugFiles:
            myDebugFiles = [f for frame in myArgs['frames']
                            for f in glob.glob(f'{orbit}_{frame}/{globPattern}')]
            if not myDebugFiles:
                continue
            debugVrt = f'{debugFrameDir}/{outName}'
            command = ['custom_buildvrtWithOffsets', '--overWrite']
            if useOffsets:
                command.append('--offsets')
            command.append(debugVrt)
            print(f'Building ROFF debug {debugVrt}')
            run(command + myDebugFiles)
    #
    for myFileType in files:
        # Find files
        myFiles = []
        for frame in myArgs['frames']:
            # print(frame, f'{orbit}_{frame}/*.{key}.vrt')
            myFile = glob.glob(f'{orbit}_{frame}/{myFileType}')
            # print(myFile)
            if len(myFile) == 1:
                myFiles += myFile
            else:
                u.mywarning(f'more than one or missing vrt {myFile}'
                          f'\n\tfor {orbit}_{frame}/{myFileType}')
        if not myFiles:
            u.mywarning(f'createVirtualFrameROFF: no input VRTs found for '
                        f'{myFileType}, skipping')
            continue
        # Create virtual raster
        command = ['custom_buildvrtWithOffsets']
        # Virtual product name — derive from the actual source frame in myFiles[0]
        _srcFrame = os.path.dirname(myFiles[0]).split('_')[-1]
        virtualProduct = os.path.basename(
            myFiles[0].replace(f'_{_srcFrame}.', f'_{virtualFrame}.'))
        #
        vrtFile = f'{frameDir}/{virtualProduct}'
        virtualVRTs[myFileType] = vrtFile
        failFile = f'{vrtFile}.fail'
        if os.path.exists(failFile):
            os.remove(failFile)
        # Force a new final file
        command.append('--overWrite')
        command.append('--offsets')
        command.append(vrtFile)
        command += myFiles
        # print(command)
        print(f'Building ROFF {frameDir}/{virtualProduct}')
        #print(' '.join(command))
        #u.myerror('stop')
        run(command)
        if not os.path.exists(virtualVRTs[myFileType]):
            open(failFile, 'w').close()
    #
    ds = gdal.Open(f'{frameDir}/range.offsets.vrt' , gdal.GA_Update)
    if 'ionosphereRangeOffsetCorrection' in myArgs:
        ds.SetMetadataItem("ionosphereRangeOffsetCorrection",
                       myArgs['ionosphereRangeOffsetCorrection'])
        ds = None  # flush and close


def _readIonosphereParams(myArgs):
    """Return (wavelength, slpSpacing, epsg) from the first available RUNW HDF5."""
    orbit = myArgs['orbit1']
    for frame in myArgs['frames']:
        frameDir = f'{orbit}_{frame}'
        h5Files = (glob.glob(f'{frameDir}/H5/NISAR*RUNW*.h5') or
                   glob.glob(f'{frameDir}/NISAR*RUNW*.h5'))
        if not h5Files:
            continue
        try:
            runw = nisarhdf.nisarRUNWHDF()
            runw.openHDF(h5Files[0], frame=frame)
            return float(runw.Wavelength), float(runw.SLCRangePixelSize), int(runw.epsg)
        except Exception as e:
            u.mywarning(f'_readIonosphereParams: could not read {h5Files[0]}: {e}')
    return None, None, None


def globalFillIonosphere(myArgs, sigmaAz=10.0, sigmaRng=30.0):
    """Globally fill both the RUNW-grid (radians) and offset-grid (SLC pixels)
    iono corrections, re-partition back into per-frame tiles, and build
    plain geometry VRTs.  Avoids frame-boundary discontinuities caused by
    independent per-frame fills.

    Also writes a phase-derived offset correction (RUNW-fill scaled by
    λ/(4π·slp)) to per-frame debug/ dirs for quality evaluation.
    """
    from nisargrimpworkflow.estimateIonosphere import (fill_and_smooth_iono,
                                                       local_consistency_mask,
                                                       write_geotiff,
                                                       write_output_vrt)
    orbit = myArgs['orbit1']
    virtualFrame = myArgs['virtualFrame']
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'

    # --- locate assembled sparse offset-grid iono VRT ---
    unfilledVrts = glob.glob(f'{frameDir}/*.ionosphereCorrectionUnfilled.offset.vrt')
    if not unfilledVrts:
        u.mywarning('globalFillIonosphere: no unfilled offset iono VRT found, skipping')
        return
    unfilledVrt = unfilledVrts[0]

    # --- read sparse offset-grid iono (SLC pixels, DC-corrected across frames) ---
    ds = gdal.Open(unfilledVrt)
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    arr = band.ReadAsArray().astype(np.float64)
    fullGt = ds.GetGeoTransform()
    fullProj = ds.GetProjection()
    ds = None

    valid = np.isfinite(arr)
    if nodata is not None:
        valid &= (arr != nodata)
    arr[~valid] = np.nan

    if not valid.any():
        u.mywarning('globalFillIonosphere: no valid pixels, skipping')
        return

    # --- projection fallback from RUNW HDF ---
    wavelength, slpSpacing, epsg = _readIonosphereParams(myArgs)
    if epsg is not None:
        _sr = osr.SpatialReference()
        _sr.ImportFromEPSG(epsg)
        hdfProjWkt = _sr.ExportToWkt()
    else:
        hdfProjWkt = ''
    projWkt = fullProj or hdfProjWkt

    # --- global fill on full-swath RUNW grid (radians) ---
    memRunwDs = None
    runwProjWkt = ''
    perFrameRunwVrts = []
    unfilledRunwVrts = glob.glob(f'{frameDir}/*.ionosphereCorrectionUnfilled.vrt')
    if unfilledRunwVrts:
        dsR = gdal.Open(unfilledRunwVrts[0])
        bandR = dsR.GetRasterBand(1)
        nodataR = bandR.GetNoDataValue()
        arrR = bandR.ReadAsArray().astype(np.float64)
        runwFullGt = dsR.GetGeoTransform()
        runwFullProj = dsR.GetProjection()
        dsR = None
        validR = np.isfinite(arrR)
        if nodataR is not None:
            validR &= (arrR != nodataR)
        arrR[~validR] = np.nan
        runwProjWkt = runwFullProj or hdfProjWkt
        if validR.any():
            print('globalFillIonosphere: running pyramid fill on full-swath RUNW-grid iono...')
            filledRunw = fill_and_smooth_iono(arrR.astype(np.float32), validR,
                                              sigma=(sigmaAz, sigmaRng),
                                              boundary_mode='nearest')
            nrowsR, ncolsR = filledRunw.shape
            memRunwDs = gdal.GetDriverByName('MEM').Create(
                '', ncolsR, nrowsR, 1, gdal.GDT_Float32)
            memRunwDs.SetGeoTransform(runwFullGt)
            if runwProjWkt:
                memRunwDs.SetProjection(runwProjWkt)
            memRunwDs.GetRasterBand(1).WriteArray(filledRunw)
            memRunwDs.GetRasterBand(1).SetNoDataValue(float('nan'))
            for frame in myArgs['frames']:
                unfilledFrameTifs = glob.glob(
                    f'{orbit}_{frame}/*.ionosphereCorrectionUnfilled.tif')
                if not unfilledFrameTifs:
                    u.mywarning(
                        f'globalFillIonosphere: no unfilled RUNW tif for {orbit}_{frame}')
                    continue
                refTif = unfilledFrameTifs[0]
                refDs = gdal.Open(refTif)
                fGt = refDs.GetGeoTransform()
                fProj = refDs.GetProjection() or runwProjWkt
                fNcols = refDs.RasterXSize
                fNrows = refDs.RasterYSize
                refDs = None
                perFrameRunwTif = refTif.replace(
                    '.ionosphereCorrectionUnfilled.tif',
                    '.ionosphereCorrection.globalFill.tif')
                perFrameRunwVrt = perFrameRunwTif.replace('.tif', '.vrt')
                gdal.Warp(perFrameRunwTif, memRunwDs, options=gdal.WarpOptions(
                    format='GTiff',
                    outputBounds=(fGt[0],
                                  fGt[3] + fNrows * fGt[5],
                                  fGt[0] + fNcols * fGt[1],
                                  fGt[3]),
                    width=fNcols, height=fNrows,
                    dstSRS=fProj,
                    resampleAlg='bilinear',
                    srcNodata=float('nan'), dstNodata=-2.0e9,
                    creationOptions=['COMPRESS=LZW']))
                _ds = gdal.Open(perFrameRunwTif, gdal.GA_Update)
                if _ds is not None:
                    _ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
                    _ds.FlushCache()
                    _ds = None
                write_output_vrt(perFrameRunwVrt, [perFrameRunwTif],
                                 ['ionosphereCorrection'], fGt)
                perFrameRunwVrts.append(perFrameRunwVrt)
                print(f'globalFillIonosphere: wrote {perFrameRunwVrt}')
            if perFrameRunwVrts:
                runwVrtBase = os.path.basename(unfilledRunwVrts[0]).replace(
                    '.ionosphereCorrectionUnfilled.vrt', '.ionosphereCorrection.vrt')
                globalFillRunwVrt = f'{frameDir}/{runwVrtBase}'
                run(['custom_buildvrtWithOffsets', '--overWrite',
                     globalFillRunwVrt] + perFrameRunwVrts,
                    stdout=myArgs.get('stdout', DEVNULL),
                    stderr=myArgs.get('stderr', DEVNULL))
                print(f'globalFillIonosphere: RUNW VRT → {globalFillRunwVrt}')
    else:
        u.mywarning(
            'globalFillIonosphere: no unfilled RUNW iono VRT found, skipping phase fill')

    # --- ice-free pre-fill / sepIceRock ice-anchored bias + rock-truth seeding ---
    if myArgs.get('sepIceRock'):
        # Ice pixels (mask==2): compare each pixel's own measured range offset
        # plus its own ice-derived ionosphere estimate against the
        # velocity-inclusive simulation (offsets.velocity, NOT offsets.geom --
        # ice has real motion, so geom alone would force the bias to absorb the
        # real velocity signal). Averaging this residual over every ice pixel
        # in the virtual frame gives a single DC bias C, applied directly to
        # the still-unfilled `arr` -- no extrapolation/smoothing involved, so a
        # localized bad pixel elsewhere can't skew C through a smoothing kernel.
        #
        # Rock pixels (mask==1, v=0 so offsets.geom==offsets.velocity there)
        # give an independent, absolute truth (offsets.geom - range.offsets),
        # seeded directly into arr's remaining gaps -- never compared against
        # the ice-derived estimate, unlike the previous sIce-vs-fv_full design.
        #
        # Both steps run on the unfilled arr/valid; the single shared
        # fill_and_smooth_iono call below (after this whole if/else) produces
        # the final continuous output from the resulting hybrid array.
        _nfr, _nfc = arr.shape
        _fv_sum = np.zeros((_nfr, _nfc), dtype=np.float64)
        _fv_count = np.zeros((_nfr, _nfc), dtype=np.int32)
        _ice_resid_sum = np.zeros((_nfr, _nfc), dtype=np.float64)
        _ice_resid_count = np.zeros((_nfr, _nfc), dtype=np.int32)
        _n_with_mask = 0
        for frame in myArgs['frames']:
            _fdir = f'{orbit}_{frame}'
            _mask_path = os.path.join(_fdir, 'offsetSims', 'offsets.geom.mask.vrt')
            _geom_paths = (glob.glob(f'{_fdir}/offsetSims/offsets.geom.vrt') +
                           glob.glob(f'{_fdir}/H5/offsetSims/offsets.geom.vrt'))
            _vel_paths = (glob.glob(f'{_fdir}/offsetSims/offsets.velocity.vrt') +
                          glob.glob(f'{_fdir}/H5/offsetSims/offsets.velocity.vrt'))
            _roff_paths = (glob.glob(f'{_fdir}/range.offsets.vrt') +
                           glob.glob(f'{_fdir}/H5/range.offsets.vrt'))

            if (not os.path.exists(_mask_path) or not _geom_paths
                    or not _vel_paths or not _roff_paths):
                u.mywarning(f'globalFillIonosphere: no mask/geom/velocity/roff for '
                           f'{_fdir}, skipping')
                continue

            _mds = gdal.Open(_mask_path)
            _mask = _mds.GetRasterBand(1).ReadAsArray()
            _fgt = _mds.GetGeoTransform()
            _fnr, _fnc = _mds.RasterYSize, _mds.RasterXSize
            _mds = None

            _gds = gdal.Open(_geom_paths[0])
            _geom_r = _gds.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _gnd = _gds.GetRasterBand(1).GetNoDataValue()
            _gds = None

            _vds = gdal.Open(_vel_paths[0])
            _vel_r = _vds.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _vnd = _vds.GetRasterBand(1).GetNoDataValue()
            _vds = None

            _rds = gdal.Open(_roff_paths[0])
            _roff_r = _rds.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _rnd = _rds.GetRasterBand(1).GetNoDataValue()
            _rds = None

            _g_ok = (_geom_r != _gnd) if _gnd is not None else np.ones(_mask.shape, bool)
            _v_ok = np.isfinite(_vel_r)
            if _vnd is not None:
                _v_ok &= (_vel_r != _vnd)
            _r_ok = np.isfinite(_roff_r)
            if _rnd is not None:
                _r_ok &= (_roff_r != _rnd)

            # Rock pixels (value=1); erode 5×5 to pull back from ice/water edges
            _rock = (_mask == 1)
            _rock_eroded = binary_erosion(_rock, structure=np.ones((5, 5), bool))
            _fv = np.where(_rock_eroded & _g_ok & _r_ok,
                           (_geom_r - _roff_r).astype(np.float32), np.nan)

            # Ice pixels (value=2): measured - velocitySim; arr's own ice-derived
            # estimate is added in after remapping to the full-swath grid below.
            _ice = (_mask == 2)
            _iceResid = np.where(_ice & _v_ok & _r_ok,
                                 (_roff_r - _vel_r).astype(np.float32), np.nan)

            # Map per-frame pixel coordinates to full-swath array
            _r0 = round((_fgt[3] - fullGt[3]) / fullGt[5])
            _c0 = round((_fgt[0] - fullGt[0]) / fullGt[1])
            _r1, _c1 = _r0 + _fnr, _c0 + _fnc
            _R0 = max(0, _r0); _R1 = min(_nfr, _r1)
            _C0 = max(0, _c0); _C1 = min(_nfc, _c1)
            if _R1 <= _R0 or _C1 <= _C0:
                continue
            _fr0, _fc0 = _R0 - _r0, _C0 - _c0
            _fr1, _fc1 = _fr0 + _R1 - _R0, _fc0 + _C1 - _C0

            _fv_slice = _fv[_fr0:_fr1, _fc0:_fc1]
            _fv_ok = np.isfinite(_fv_slice)
            _fv_sum[_R0:_R1, _C0:_C1][_fv_ok] += _fv_slice[_fv_ok]
            _fv_count[_R0:_R1, _C0:_C1][_fv_ok] += 1

            # Ice residual + arr's own ice-derived estimate, only where arr is
            # itself already valid (unfilled data -- no extrapolation).
            _iceResid_slice = _iceResid[_fr0:_fr1, _fc0:_fc1]
            _arr_slice = arr[_R0:_R1, _C0:_C1]
            _valid_slice = valid[_R0:_R1, _C0:_C1]
            _iceFull_slice = _iceResid_slice + _arr_slice
            _ice_ok = np.isfinite(_iceFull_slice) & _valid_slice
            _ice_resid_sum[_R0:_R1, _C0:_C1][_ice_ok] += _iceFull_slice[_ice_ok]
            _ice_resid_count[_R0:_R1, _C0:_C1][_ice_ok] += 1

            _n_with_mask += 1
            print(f'globalFillIonosphere: sepIceRock found {int(_fv_ok.sum())} rock px, '
                  f'{int(_ice_ok.sum())} ice px in {_fdir}')

        if _n_with_mask == 0:
            print('\033[1;34mWARNING globalFillIonosphere: offsets.geom.mask.vrt not found '
                  'for any frame — sepIceRock correction skipped\033[0m')
        else:
            rock_mask_full = _fv_count > 0
            fv_full = np.where(rock_mask_full, _fv_sum / np.maximum(_fv_count, 1), np.nan)
            ice_resid_mask_full = _ice_resid_count > 0

            # --- Step 1: ice-anchored bias, directly on unfilled arr, vs
            # offsets.velocity (not offsets.geom -- ice has real motion) ---
            if ice_resid_mask_full.any():
                _resid_ice = (_ice_resid_sum / np.maximum(_ice_resid_count, 1))[ice_resid_mask_full]
                _resid_ice = _resid_ice[np.isfinite(_resid_ice)]
            else:
                _resid_ice = np.array([])
            if len(_resid_ice) > 0:
                cIce = float(np.mean(_resid_ice))
                _dc_m = f'{cIce * slpSpacing:+.3f} m' if slpSpacing else ''
                print(f'globalFillIonosphere: sepIceRock ice-anchored bias '
                      f'C = {cIce:+.4f} SLC px {_dc_m} '
                      f'(from {len(_resid_ice):,} ice px, vs offsets.velocity)')
                arr -= cIce
            else:
                u.mywarning('globalFillIonosphere: sepIceRock found no ice residuals '
                            'vs offsets.velocity — bias correction skipped')

            # --- Step 2: seed rock truth directly into arr's gaps -- never
            # compared against the ice-derived estimate ---
            if rock_mask_full.any():
                # Local-agreement check rather than a flat magnitude cutoff:
                # a single bad/misclassified rock pixel disagrees with its
                # neighbors regardless of magnitude, while a large but
                # spatially-coherent regional discrepancy (every nearby
                # rock pixel agrees) should still be trusted as truth.
                _consistent = local_consistency_mask(fv_full, rock_mask_full)
                _apply = _consistent & ~valid
                arr[_apply] = fv_full[_apply]
                valid |= _apply
                print(f'globalFillIonosphere: sepIceRock seeded '
                      f'{int(_apply.sum()):,} rock px into gaps')
            else:
                print('\033[1;34mWARNING globalFillIonosphere: sepIceRock found no rock '
                      'pixels — rock seeding skipped\033[0m')
    else:
        # --- ice-free pre-fill: offsets.geom − range.offsets at non-ice pixels ---
        _n_with_mask = 0
        for frame in myArgs['frames']:
            _fdir = f'{orbit}_{frame}'
            _mask_path = os.path.join(_fdir, 'offsetSims', 'offsets.geom.mask.vrt')
            _geom_paths = (glob.glob(f'{_fdir}/offsetSims/offsets.geom.vrt') +
                           glob.glob(f'{_fdir}/H5/offsetSims/offsets.geom.vrt'))
            _roff_paths = (glob.glob(f'{_fdir}/range.offsets.vrt') +
                           glob.glob(f'{_fdir}/H5/range.offsets.vrt'))

            if not os.path.exists(_mask_path) or not _geom_paths or not _roff_paths:
                u.mywarning(f'globalFillIonosphere: no mask/geom/roff for {_fdir}, skipping')
                continue

            _mds = gdal.Open(_mask_path)
            _mask = _mds.GetRasterBand(1).ReadAsArray()
            _fgt = _mds.GetGeoTransform()
            _fnr, _fnc = _mds.RasterYSize, _mds.RasterXSize
            _mds = None

            _gds = gdal.Open(_geom_paths[0])
            _geom_r = _gds.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _gnd = _gds.GetRasterBand(1).GetNoDataValue()
            _gds = None

            _rds = gdal.Open(_roff_paths[0])
            _roff_r = _rds.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _rnd = _rds.GetRasterBand(1).GetNoDataValue()
            _rds = None

            # Rock pixels (value=1); erode 5×5 to pull back from ice/water edges
            _rock = (_mask == 1)
            _rock_eroded = binary_erosion(_rock, structure=np.ones((5, 5), bool))
            _g_ok = (_geom_r != _gnd) if _gnd is not None else np.ones(_mask.shape, bool)
            _r_ok = np.isfinite(_roff_r)
            if _rnd is not None:
                _r_ok &= (_roff_r != _rnd)
            _fv = np.where(_rock_eroded & _g_ok & _r_ok,
                           (_geom_r - _roff_r).astype(np.float32), np.nan)

            # Map per-frame pixel coordinates to full-swath array
            _r0 = round((_fgt[3] - fullGt[3]) / fullGt[5])
            _c0 = round((_fgt[0] - fullGt[0]) / fullGt[1])
            _r1, _c1 = _r0 + _fnr, _c0 + _fnc
            _nfr, _nfc = arr.shape
            _R0 = max(0, _r0); _R1 = min(_nfr, _r1)
            _C0 = max(0, _c0); _C1 = min(_nfc, _c1)
            if _R1 <= _R0 or _C1 <= _C0:
                continue
            _fr0, _fc0 = _R0 - _r0, _C0 - _c0
            _fr1, _fc1 = _fr0 + _R1 - _R0, _fc0 + _C1 - _C0

            # DC alignment: arr ice values are zero-meaned per frame (iono_mean was
            # subtracted in estimate_and_correct); _fv is the raw absolute iono.
            # Align rock values to the same DC level as the surrounding ice pixels.
            _fv_rock_vals = _fv[np.isfinite(_fv)]
            _arr_ice_vals = arr[_R0:_R1, _C0:_C1][valid[_R0:_R1, _C0:_C1]]
            if len(_fv_rock_vals) > 0 and len(_arr_ice_vals) > 0:
                _dc = float(np.mean(_fv_rock_vals)) - float(np.mean(_arr_ice_vals))
                _fv = np.where(np.isfinite(_fv), _fv - _dc, np.nan)
                _dc_m = f'{_dc * slpSpacing:+.3f} m' if slpSpacing else ''
                print(f'  DC corrected: {_dc:+.4f} SLC px {_dc_m}')
            elif len(_fv_rock_vals) > 0:
                # No ice pixels in this frame footprint — zero-centre as fallback
                _fv = np.where(np.isfinite(_fv),
                               _fv - float(np.mean(_fv_rock_vals)), np.nan)
                print(f'  No ice in frame footprint — rock values zero-centred')

            # Local-agreement check rather than a flat magnitude cutoff -- see
            # the matching --sepIceRock branch above for the rationale.
            _consistent = local_consistency_mask(_fv, np.isfinite(_fv))
            _fv_slice = _fv[_fr0:_fr1, _fc0:_fc1]
            _apply = (_consistent[_fr0:_fr1, _fc0:_fc1]
                      & ~valid[_R0:_R1, _C0:_C1])
            arr[_R0:_R1, _C0:_C1][_apply] = _fv_slice[_apply]
            valid[_R0:_R1, _C0:_C1] |= _apply
            _n_with_mask += 1
            print(f'globalFillIonosphere: ice-free pre-fill {_apply.sum()} px in {_fdir}')

        if _n_with_mask == 0:
            print('\033[1;34mWARNING globalFillIonosphere: offsets.geom.mask.vrt not found '
                  'for any frame — ice-free pre-fill skipped\033[0m')

    # --- global fill on full-swath offset grid (values already in SLC pixels) ---
    print('globalFillIonosphere: running pyramid fill on full-swath offset-grid iono...')
    filled = fill_and_smooth_iono(arr.astype(np.float32), valid,
                                  sigma=(sigmaAz, sigmaRng),
                                  boundary_mode='nearest')

    # Store in a MEM dataset for per-frame cropping
    nrows, ncols = filled.shape
    memDs = gdal.GetDriverByName('MEM').Create('', ncols, nrows, 1, gdal.GDT_Float32)
    memDs.SetGeoTransform(fullGt)
    if projWkt:
        memDs.SetProjection(projWkt)
    memDs.GetRasterBand(1).WriteArray(filled)
    memDs.GetRasterBand(1).SetNoDataValue(float('nan'))

    # --- re-partition: crop filled surface to each frame's offset grid ---
    perFrameOffsetVrts = []
    for frame in myArgs['frames']:
        frameRoffVrts = glob.glob(f'{orbit}_{frame}/range.offsets.vrt')
        if not frameRoffVrts:
            u.mywarning(f'globalFillIonosphere: no range.offsets.vrt for {orbit}_{frame}')
            continue

        roffDs = gdal.Open(frameRoffVrts[0])
        roffGt = roffDs.GetGeoTransform()
        roffProj = roffDs.GetProjection()
        roffNcols = roffDs.RasterXSize
        roffNrows = roffDs.RasterYSize
        roffDs = None
        if not roffProj:
            roffProj = projWkt

        # Derive output name from the per-frame ionosphereCorrection.offset.vrt
        existingOffVrts = glob.glob(f'{orbit}_{frame}/*.ionosphereCorrection.offset.vrt')
        if not existingOffVrts:
            u.mywarning(f'globalFillIonosphere: no offset vrt for {orbit}_{frame}')
            continue
        perFrameOffsetTif = existingOffVrts[0].replace(
            '.ionosphereCorrection.offset.vrt',
            '.ionosphereCorrection.globalFill.offset.tif')
        perFrameOffsetVrt = perFrameOffsetTif.replace('.tif', '.vrt')

        # Crop filled iono (SLC pixels, full-swath offset grid) to frame's offset grid
        warpOpts = gdal.WarpOptions(
            format='GTiff',
            outputBounds=(roffGt[0],
                          roffGt[3] + roffNrows * roffGt[5],
                          roffGt[0] + roffNcols * roffGt[1],
                          roffGt[3]),
            width=roffNcols, height=roffNrows,
            dstSRS=roffProj,
            resampleAlg='bilinear',
            srcNodata=float('nan'), dstNodata=-2.0e9,
            creationOptions=['COMPRESS=LZW'])
        gdal.Warp(perFrameOffsetTif, memDs, options=warpOpts)
        _ds = gdal.Open(perFrameOffsetTif, gdal.GA_Update)
        if _ds is not None:
            _ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
            _ds.FlushCache()
            _ds = None

        write_output_vrt(perFrameOffsetVrt, [perFrameOffsetTif],
                         ['ionosphereCorrection'], roffGt)
        perFrameOffsetVrts.append(perFrameOffsetVrt)
        print(f'globalFillIonosphere: wrote {perFrameOffsetVrt}')

        # --- recompute debug/range.offsets.corrected with the final,
        # globally-filled (and, with --sepIceRock, rock-anchored) iono
        # correction, replacing the per-frame value written by
        # estimateIonosphere.py before global fill ran ---
        if myArgs.get('debugIono'):
            _slpDs = gdal.Open(frameRoffVrts[0])
            _offsetSlp = _slpDs.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _slpNd = _slpDs.GetRasterBand(1).GetNoDataValue()
            _slpDs = None

            _ionoDs = gdal.Open(perFrameOffsetTif)
            _ionoFinal = _ionoDs.GetRasterBand(1).ReadAsArray().astype(np.float64)
            _ionoNd = _ionoDs.GetRasterBand(1).GetNoDataValue()
            _ionoDs = None

            _slpValid = np.isfinite(_offsetSlp)
            if _slpNd is not None:
                _slpValid &= (_offsetSlp != _slpNd)
            _ionoValid = np.isfinite(_ionoFinal)
            if _ionoNd is not None:
                _ionoValid &= (_ionoFinal != _ionoNd)

            _correctedFinal = np.where(_slpValid & _ionoValid,
                                       (_offsetSlp + _ionoFinal).astype(np.float32),
                                       -2.0e9).astype(np.float32)

            _debugDir = f'{orbit}_{frame}/debug'
            os.makedirs(_debugDir, exist_ok=True)
            _corrTif = os.path.join(_debugDir, 'range.offsets.corrected.tif')
            _corrVrt = os.path.join(_debugDir, 'range.offsets.corrected.vrt')
            write_geotiff(_corrTif, _correctedFinal, roffGt, nodata=-2.0e9, proj=roffProj)
            write_output_vrt(_corrVrt, [_corrTif], ['RangeOffsetsCorrected'], roffGt)
            print(f'globalFillIonosphere: updated final-corrected offsets debug → {_corrVrt}')

    if not perFrameOffsetVrts:
        u.mywarning('globalFillIonosphere: no per-frame tiles written, skipping VRT assembly')
        return

    # --- build plain geometry VRT in virtual-frame dir (no --offsets: biases baked in) ---
    _srcFrame = os.path.dirname(perFrameOffsetVrts[0]).split('_')[-1]
    virtualProduct = os.path.basename(perFrameOffsetVrts[0]).replace(
        f'_{_srcFrame}.', f'_{virtualFrame}.')
    globalFillOffsetVrt = f'{frameDir}/{virtualProduct}'
    run(['custom_buildvrtWithOffsets', '--overWrite',
         globalFillOffsetVrt] + perFrameOffsetVrts,
        stdout=myArgs.get('stdout', DEVNULL),
        stderr=myArgs.get('stderr', DEVNULL))

    # --- stamp metadata on virtual-frame range.offsets.vrt ---
    roffVrt = f'{frameDir}/range.offsets.vrt'
    myArgs['ionosphereRangeOffsetCorrection'] = virtualProduct
    ds = gdal.Open(roffVrt, gdal.GA_Update)
    if ds is not None:
        ds.SetMetadataItem('ionosphereRangeOffsetCorrection', virtualProduct)
        ds = None
    print(f'globalFillIonosphere: virtual frame VRT → {globalFillOffsetVrt}')

    # --- phase-derived offset correction from RUNW global fill → always written to debug/ ---
    if memRunwDs is not None and wavelength and slpSpacing:
        scale = wavelength / (4.0 * np.pi * slpSpacing)
        phaseData = memRunwDs.GetRasterBand(1).ReadAsArray().astype(np.float32)
        scaledData = (phaseData * scale).astype(np.float32)
        nrR2, ncR2 = scaledData.shape
        scaledDs = gdal.GetDriverByName('MEM').Create('', ncR2, nrR2, 1, gdal.GDT_Float32)
        scaledDs.SetGeoTransform(memRunwDs.GetGeoTransform())
        if runwProjWkt:
            scaledDs.SetProjection(runwProjWkt)
        scaledDs.GetRasterBand(1).WriteArray(scaledData)
        scaledDs.GetRasterBand(1).SetNoDataValue(float('nan'))
        phaseBasedOffsetVrts = []
        for frame in myArgs['frames']:
            frameRoffVrts = glob.glob(f'{orbit}_{frame}/range.offsets.vrt')
            if not frameRoffVrts:
                continue
            roffDs = gdal.Open(frameRoffVrts[0])
            rGt = roffDs.GetGeoTransform()
            rProj = roffDs.GetProjection() or projWkt
            rNcols = roffDs.RasterXSize
            rNrows = roffDs.RasterYSize
            roffDs = None
            existingGfVrts = glob.glob(
                f'{orbit}_{frame}/*.ionosphereCorrection.globalFill.offset.vrt')
            if not existingGfVrts:
                continue
            pbDebugDir = f'{orbit}_{frame}/debug'
            os.makedirs(pbDebugDir, exist_ok=True)
            pbTif = os.path.join(pbDebugDir, os.path.basename(existingGfVrts[0]).replace(
                '.ionosphereCorrection.globalFill.offset.vrt',
                '.ionosphereCorrection.phaseBasedGlobalFill.offset.tif'))
            pbVrt = pbTif.replace('.tif', '.vrt')
            gdal.Warp(pbTif, scaledDs, options=gdal.WarpOptions(
                format='GTiff',
                outputBounds=(rGt[0], rGt[3] + rNrows * rGt[5],
                              rGt[0] + rNcols * rGt[1], rGt[3]),
                width=rNcols, height=rNrows,
                dstSRS=rProj,
                resampleAlg='bilinear',
                srcNodata=float('nan'), dstNodata=-2.0e9,
                creationOptions=['COMPRESS=LZW']))
            _ds = gdal.Open(pbTif, gdal.GA_Update)
            if _ds is not None:
                _ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
                _ds.FlushCache()
                _ds = None
            write_output_vrt(pbVrt, [pbTif], ['ionosphereCorrection'], rGt)
            phaseBasedOffsetVrts.append(pbVrt)
            print(f'globalFillIonosphere: wrote phase-derived offset → {pbVrt}')
        if phaseBasedOffsetVrts:
            _srcFrame = os.path.dirname(phaseBasedOffsetVrts[0]).split('_')[-1]
            pbVfBase = os.path.basename(phaseBasedOffsetVrts[0]).replace(
                f'_{_srcFrame}.', f'_{virtualFrame}.')
            vfDebugDir = f'{frameDir}/debug'
            os.makedirs(vfDebugDir, exist_ok=True)
            pbVfVrt = os.path.join(vfDebugDir, pbVfBase)
            run(['custom_buildvrtWithOffsets', '--overWrite',
                 pbVfVrt] + phaseBasedOffsetVrts,
                stdout=myArgs.get('stdout', DEVNULL),
                stderr=myArgs.get('stderr', DEVNULL))
            print(f'globalFillIonosphere: phase-derived offset VRT → {pbVfVrt}')

    # --- debug: move intermediates into per-frame debug/ subdirs, build debug VRTs ---
    if myArgs.get('debugIono'):
        debugUnfilledVrts = []
        debugOffsetUnfilledVrts = []
        debugAmbigVrts = []
        debugCorrectedOffsetVrts = []
        # Build a full-swath MEM dataset from arr (post-ice-free-prefill, pre-pyramid-fill)
        # arr is still float64 with NaN for unfilled; fill_and_smooth_iono worked on a copy.
        _preFillArr = arr.astype(np.float32)
        nrPF, ncPF = _preFillArr.shape
        _preFillDs = gdal.GetDriverByName('MEM').Create('', ncPF, nrPF, 1, gdal.GDT_Float32)
        _preFillDs.SetGeoTransform(fullGt)
        if projWkt:
            _preFillDs.SetProjection(projWkt)
        _preFillDs.GetRasterBand(1).WriteArray(_preFillArr)
        _preFillDs.GetRasterBand(1).SetNoDataValue(float('nan'))
        for frame in myArgs['frames']:
            orbitFrame = f'{orbit}_{frame}'
            debugDir = f'{orbitFrame}/debug'
            os.makedirs(debugDir, exist_ok=True)
            # RUNW-grid unfilled iono (radians) → debug/
            for unfTif in glob.glob(f'{orbitFrame}/*.ionosphereCorrectionUnfilled.tif'):
                destTif = os.path.join(debugDir, os.path.basename(unfTif))
                os.rename(unfTif, destTif)
                _ds = gdal.Open(destTif)
                _gt = _ds.GetGeoTransform()
                _ds = None
                destVrt = destTif.replace('.tif', '.vrt')
                write_output_vrt(destVrt, [destTif],
                                 ['ionosphereCorrectionUnfilled'], _gt)
                debugUnfilledVrts.append(destVrt)
            for oldVrt in glob.glob(f'{orbitFrame}/*.ionosphereCorrectionUnfilled.vrt'):
                os.remove(oldVrt)
            # Offset-grid pre-global-fill iono (sparse + ice-free seeds) → debug/
            # Warp from the full-swath pre-fill MEM dataset rather than moving the
            # original unfilled tif, so the ice-free pre-fill seeds are included.
            for unfOffTif in glob.glob(f'{orbitFrame}/*.ionosphereCorrectionUnfilled.offset.tif'):
                destName = os.path.basename(unfOffTif).replace(
                    '.ionosphereCorrectionUnfilled.offset.tif',
                    '.ionosphereCorrection.preGlobalFill.offset.tif')
                destTif = os.path.join(debugDir, destName)
                _pfDs = gdal.Open(unfOffTif)
                _pfGt = _pfDs.GetGeoTransform()
                _pfProj = _pfDs.GetProjection() or projWkt
                _pfNcols = _pfDs.RasterXSize
                _pfNrows = _pfDs.RasterYSize
                _pfDs = None
                os.remove(unfOffTif)
                gdal.Warp(destTif, _preFillDs, options=gdal.WarpOptions(
                    format='GTiff',
                    outputBounds=(_pfGt[0], _pfGt[3] + _pfNrows * _pfGt[5],
                                  _pfGt[0] + _pfNcols * _pfGt[1], _pfGt[3]),
                    width=_pfNcols, height=_pfNrows,
                    dstSRS=_pfProj,
                    resampleAlg='bilinear',
                    srcNodata=float('nan'), dstNodata=-2.0e9,
                    creationOptions=['COMPRESS=LZW']))
                _ds = gdal.Open(destTif, gdal.GA_Update)
                if _ds is not None:
                    _ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
                    _ds.FlushCache()
                    _ds = None
                destVrt = destTif.replace('.tif', '.vrt')
                write_output_vrt(destVrt, [destTif], ['ionosphereCorrection'], _pfGt)
                debugOffsetUnfilledVrts.append(destVrt)
            for oldVrt in glob.glob(f'{orbitFrame}/*.ionosphereCorrectionUnfilled.offset.vrt'):
                os.remove(oldVrt)
            # Per-frame filled iono offset (superseded by globalFill) → debug/
            # Pattern intentionally excludes *.globalFill.offset.tif
            for offTif in glob.glob(f'{orbitFrame}/*.ionosphereCorrection.offset.tif'):
                destTif = os.path.join(debugDir, os.path.basename(offTif))
                os.rename(offTif, destTif)
                oldVrt = offTif.replace('.tif', '.vrt')
                if os.path.exists(oldVrt):
                    os.rename(oldVrt, os.path.join(debugDir, os.path.basename(oldVrt)))
            # Per-frame filled RUNW iono (superseded by globalFill) → debug/
            # Pattern intentionally excludes *.globalFill.tif
            for ionoTif in glob.glob(f'{orbitFrame}/*.ionosphereCorrection.tif'):
                destTif = os.path.join(debugDir, os.path.basename(ionoTif))
                os.rename(ionoTif, destTif)
                oldVrt = ionoTif.replace('.tif', '.vrt')
                if os.path.exists(oldVrt):
                    os.rename(oldVrt, os.path.join(debugDir, os.path.basename(oldVrt)))
            # Ambiguity-corrected (pre-iono) phase VRTs written by estimateIonosphere
            debugAmbigVrts += sorted(
                glob.glob(f'{debugDir}/*.ambiguityCorrectedUnwrappedPhase.vrt'))
            # Ionosphere-corrected range offsets written by estimateIonosphere
            debugCorrectedOffsetVrts += sorted(
                glob.glob(f'{debugDir}/range.offsets.corrected.vrt'))
        # Remove broken assembled virtual-frame VRTs (per-frame TIFs have moved)
        for pattern in ['*.ionosphereCorrectionUnfilled.vrt',
                        '*.ionosphereCorrectionUnfilled.offset.vrt']:
            for oldVrt in glob.glob(f'{frameDir}/{pattern}'):
                os.remove(oldVrt)
        # Virtual-frame debug dir
        vfDebugDir = f'{frameDir}/debug'
        os.makedirs(vfDebugDir, exist_ok=True)
        # Assemble virtual-frame RUNW-grid unfilled iono VRT in debug/
        if debugUnfilledVrts:
            _srcFrame = os.path.dirname(debugUnfilledVrts[0]).split('_')[-1]
            vfProduct = os.path.basename(debugUnfilledVrts[0]).replace(
                f'_{_srcFrame}.', f'_{virtualFrame}.')
            vfVrt = os.path.join(vfDebugDir, vfProduct)
            run(['custom_buildvrtWithOffsets', '--offsets', '--overWrite',
                 vfVrt] + debugUnfilledVrts,
                stdout=myArgs.get('stdout', DEVNULL),
                stderr=myArgs.get('stderr', DEVNULL))
            print(f'globalFillIonosphere: debug unfilled VRT → {vfVrt}')
        # Assemble virtual-frame pre-global-fill offset VRT in debug/
        if debugOffsetUnfilledVrts:
            _srcFrame = os.path.dirname(debugOffsetUnfilledVrts[0]).split('_')[-1]
            vfProduct = os.path.basename(debugOffsetUnfilledVrts[0]).replace(
                f'_{_srcFrame}.', f'_{virtualFrame}.')
            vfVrt = os.path.join(vfDebugDir, vfProduct)
            run(['custom_buildvrtWithOffsets', '--offsets', '--overWrite',
                 vfVrt] + debugOffsetUnfilledVrts,
                stdout=myArgs.get('stdout', DEVNULL),
                stderr=myArgs.get('stderr', DEVNULL))
            print(f'globalFillIonosphere: debug pre-global-fill offset VRT → {vfVrt}')
        # Assemble virtual-frame ambiguity-corrected phase VRT in debug/
        if debugAmbigVrts:
            _srcFrame = os.path.dirname(debugAmbigVrts[0]).split('_')[-1]
            vfProduct = os.path.basename(debugAmbigVrts[0]).replace(
                f'_{_srcFrame}.', f'_{virtualFrame}.')
            vfVrt = os.path.join(vfDebugDir, vfProduct)
            run(['custom_buildvrtWithOffsets', '--offsets', '--overWrite',
                 vfVrt] + debugAmbigVrts,
                stdout=myArgs.get('stdout', DEVNULL),
                stderr=myArgs.get('stderr', DEVNULL))
            print(f'globalFillIonosphere: debug ambiguity-corrected phase VRT → {vfVrt}')
        # Assemble virtual-frame ionosphere-corrected range offsets VRT in debug/
        if debugCorrectedOffsetVrts:
            vfVrt = os.path.join(vfDebugDir, 'range.offsets.corrected.vrt')
            run(['custom_buildvrtWithOffsets', '--offsets', '--overWrite',
                 vfVrt] + debugCorrectedOffsetVrts,
                stdout=myArgs.get('stdout', DEVNULL),
                stderr=myArgs.get('stderr', DEVNULL))
            print(f'globalFillIonosphere: debug corrected offsets VRT → {vfVrt}')
        _preFillDs = None
    # --- clean up superseded and intermediate files (unless --retainIntermediateIono) ---
    elif not myArgs.get('retainIntermediateIono'):
        for frame in myArgs['frames']:
            for pattern in ['*.ionosphereCorrectionUnfilled.tif',
                            '*.ionosphereCorrectionUnfilled.vrt',
                            '*.ionosphereCorrectionUnfilled.offset.tif',
                            '*.ionosphereCorrectionUnfilled.offset.vrt',
                            '*.ionosphereCorrection.tif',
                            '*.ionosphereCorrection.vrt',
                            '*.ionosphereCorrection.offset.tif',
                            '*.ionosphereCorrection.offset.vrt']:
                for f in glob.glob(f'{orbit}_{frame}/{pattern}'):
                    os.remove(f)
        for pattern in ['*.ionosphereCorrectionUnfilled.vrt',
                        '*.ionosphereCorrectionUnfilled.offset.vrt']:
            for f in glob.glob(f'{frameDir}/{pattern}'):
                os.remove(f)


def createVirtualFrameCorr(myArgs):
    '''Assemble per-frame .cor.vrt files into a single virtual-frame correlation mosaic.'''
    orbit = myArgs['orbit1']
    virtualFrame = myArgs['virtualFrame']
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
    if myArgs.get('geodatsOnly'):
        existingVrts = glob.glob(f'{frameDir}/*.cor.vrt')
        if not existingVrts:
            u.mywarning(f'createVirtualFrameCorr: --geodatsOnly requested but no '
                       f'existing .cor.vrt found in {frameDir} -- skipping')
            return
        mergedGeodat(myArgs['geodat1'], existingVrts[0])
        mergedGeodat(myArgs['geodat2'], existingVrts[0], secondary=True)
        return
    myFiles = []
    for frame in myArgs['frames']:
        myFile = glob.glob(f'{orbit}_{frame}/*.cor.vrt')
        if len(myFile) == 1:
            myFiles += myFile
        else:
            u.mywarning(f'more than one or missing .cor.vrt for '
                        f'{orbit}_{frame}: {myFile}')
    if not myFiles:
        u.mywarning('createVirtualFrameCorr: no .cor.vrt files found')
        return
    _srcFrame = os.path.dirname(myFiles[0]).split('_')[-1]
    virtualProduct = os.path.basename(
        myFiles[0].replace(f'_{_srcFrame}.', f'_{virtualFrame}.'))
    vrtFile = f'{frameDir}/{virtualProduct}'
    failFile = f'{vrtFile}.fail'
    if os.path.exists(failFile):
        os.remove(failFile)
    command = ['custom_buildvrtWithOffsets', '--overWrite', vrtFile] + myFiles
    print(f'Building cor {frameDir}/{virtualProduct}')
    run(command)
    if not os.path.exists(vrtFile):
        open(failFile, 'w').close()
        return
    mergedGeodat(myArgs['geodat1'], vrtFile)
    mergedGeodat(myArgs['geodat2'], vrtFile, secondary=True)


def processFramePow(frame, myArgs):
    '''
    check

    Parameters
    ----------
    frame : TYPE
        DESCRIPTION.
    myArgs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    orbit = myArgs['orbit1']
    print(f'{orbit}_{frame}/P{orbit}_{frame}.*x*.pow')
    powFiles = glob.glob(f'{orbit}_{frame}/P{orbit}_{frame}.*x*.pow')
    if len(powFiles) == 1:
        powFile = powFiles[0]
        myArgs['pow'].append(powFile)
        print(myArgs['pow'])
        nr, na = powFile.split('.')[-2].split('x')
        geodatFile = f'{orbit}_{frame}/geodat{nr}x{na}.geojson'
        myArgs['geodatpow'].append(geodatFile)
        if not os.path.exists(geodatFile):
            u.myerror(f'Mising {geodatFile} for powFile')
        return
    elif len(powFiles) == 1:
        u.myerror(f'To many pow files {powFiles}')    
    
def createVirtualFramePower(myArgs):
    if not myArgs['pow']:
        return
    orbit = myArgs['orbit1']
    virtualFrame = myArgs["virtualFrame"]
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
#    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
    command = ['custom_buildvrtWithOffsets.py']
    # Virtual product name
    myFiles = [x + '.vrt' for x in myArgs['pow']]
    frame = os.path.basename(myArgs['pow'][0]).split('_')[1].split('.')[0]
    virtualProduct = os.path.basename(
        myFiles[0].replace(f'_{frame}.', f'_{virtualFrame}.'))
    virtualVRT = f'{frameDir}/{virtualProduct}'
    failFile = f'{virtualVRT}.fail'
    if os.path.exists(failFile):
        os.remove(failFile)
    # Force a new final file
    command.append('--overWrite')
    command.append(f'{frameDir}/{virtualProduct}')
    command += myFiles
    # print(command)
    print(f'Building {frameDir}/{virtualProduct}')
    run(command)
    #
    mergedGeodat(myArgs['geodatpow'],  virtualVRT)


def main():
    '''
        Set up a series of nisar products for
    '''
    myArgs = parseCommandLine()
    if myArgs.get('clean') or myArgs.get('cleanDebug'):
        # Skip the verbose param dump below -- with processTrack.py calling SetupNISAR
        # once per orbit, that dump (mostly identical across orbits) buries the much
        # shorter clean summary between orbits. Just identify which orbit this is.
        print(f'\n=== Cleaning orbit {myArgs["orbit1"]} ===')
        buckets = cleanFrames(myArgs,
                              debugOnly=myArgs.get('cleanDebug')
                              and not myArgs.get('clean'))
        confirmAndRemove(buckets, myArgs.get('noPrompt'))
        return
    # Read optional regionFile from ../project.yaml
    myArgs['regionFile'] = None
    projectYaml = '../project.yaml'
    if os.path.exists(projectYaml):
        with open(projectYaml) as _fp:
            _proj = yaml.safe_load(_fp) or {}
        myArgs['regionFile'] = _proj.get('regionFile', None)
        if myArgs['regionFile']:
            print(f'regionFile from {projectYaml}: {myArgs["regionFile"]}')
    # Read optional verticalCorrection from the region file, passed through
    # to ROFFtoGrimp/RUNWtoGrimp so simoffsets/siminsar apply it.
    myArgs['verticalCorrection'] = None
    if myArgs['regionFile']:
        regionDef = sarfunc.defaultRegionDefs(None, regionFile=myArgs['regionFile'])
        myArgs['verticalCorrection'] = regionDef.verticalCorrection()
        if myArgs['verticalCorrection']:
            print(f'verticalCorrection from {myArgs["regionFile"]}: '
                  f'{myArgs["verticalCorrection"]}')
    # Read optional variable smoothing-radius map params from ../project.yaml, passed
    # through to ROFFtoGrimp/RUNWtoGrimp (see variableSmoothingArgs()). minTol/percentSpeed/
    # maxTol must be given together; missing all three leaves the feature off (default).
    smoothDefaults = {'minTol': None, 'percentSpeed': None, 'maxTol': None,
                      'maxSmoothRadius': 50, 'smoothNIter': 3,
                      'noVariableSmoothing': False}
    myArgs.update(smoothDefaults)
    if os.path.exists(projectYaml):
        with open(projectYaml) as _fp:
            _proj = yaml.safe_load(_fp) or {}
        smoothFlags = [_proj.get('minTol') is not None,
                      _proj.get('percentSpeed') is not None,
                      _proj.get('maxTol') is not None]
        if any(smoothFlags) and not all(smoothFlags):
            u.myerror('SetupNISAR: project.yaml minTol/percentSpeed/maxTol must be '
                      'given together')
        for key, default in smoothDefaults.items():
            myArgs[key] = _proj.get(key, default)
        if myArgs['minTol'] is not None:
            print(f'variable smoothing-radius map from {projectYaml}: minTol='
                  f'{myArgs["minTol"]} percentSpeed={myArgs["percentSpeed"]} maxTol='
                  f'{myArgs["maxTol"]} noVariableSmoothing='
                  f'{myArgs["noVariableSmoothing"]}')
    # Get list of frames
    myArgs['frames'] = getFrames(myArgs)
    print('Frames: ', myArgs['frames'])
    myArgs['outputDir'] = '.'
    myArgs['inputDir'] = '.'
    myArgs['geodat1'], myArgs['geodat2'] = [], []
    myArgs['pow'], myArgs['geodatpow'] = [], []
    # get second orbit
    haveData = getSecondaryOrbit(myArgs)
    if not haveData:
        print('No InSAR products for orbit: myArgs["orbit1"]')
    else:
        print('orbit2:', myArgs['orbit2'])
    # Move bandwidth-based YAML copy to helper (will run after frameDir exists)
    #
    # Process frames, recording per-frame geodat/pow so groups can be assembled
    # independently when there are gaps in the frame sequence.
    #
    perFrameGeodat = {}   # frame -> (geodat1, geodat2)
    perFramePow    = {}   # frame -> (pow, geodatpow)
    t0 = time.time()
    for frame in myArgs['frames']:
        # Wrap all H5 files in the frame directory with VRT
        wrapH5sInFrameDir(myArgs['orbit1'], frame, verbose=myArgs['verbose'])
        if haveData:
            print(f'Processing Frame {frame}...')
            t_step = time.time()
            print('\tRUNW....', end=' ', flush=True)
            n1_before = len(myArgs['geodat1'])
            mixedMode = processFrameRUNW(frame, myArgs)
            if len(myArgs['geodat1']) > n1_before:
                perFrameGeodat[frame] = (myArgs['geodat1'][-1], myArgs['geodat2'][-1])
            dt = time.time() - t_step
            print(f'{dt:.1f}s  (total {time.time()-t0:.1f}s)')
            if not mixedMode and not myArgs['RUNWOnly'] \
                    and not myArgs.get('correlationOnly') \
                    and not myArgs.get('corrOnly') \
                    and not myArgs.get('geodatsOnly'):
                t_step = time.time()
                print('\tROFF....', end=' ', flush=True)
                processFrameROFF(frame, myArgs)
                dt = time.time() - t_step
                print(f'{dt:.1f}s  (total {time.time()-t0:.1f}s)')
                if not myArgs.get('phaseDerivedIonosphere'):
                    t_step = time.time()
                    print('\tIonosphere (estimateIonosphere)....', end=' ', flush=True)
                    processFrameIonosphere(frame, myArgs, simDir='simPhase')
                    dt = time.time() - t_step
                    print(f'{dt:.1f}s  (total {time.time()-t0:.1f}s)')
        # Need to add check for power mixed mode
        if not myArgs.get('corrOnly') and not myArgs.get('geodatsOnly'):
            np_before = len(myArgs['pow'])
            processFramePow(frame, myArgs)
            if len(myArgs['pow']) > np_before:
                perFramePow[frame] = (myArgs['pow'][-1], myArgs['geodatpow'][-1])

    # Split frames into contiguous groups; each gap creates a new virtual frame.
    allFrames = myArgs['frames']
    groups = splitFrameGroups(allFrames)
    # Further split any group containing a frame whose own secondary date
    # disagrees with its neighbors' (real data gap upstream, not safe to merge).
    groups = splitGroupsBySecondaryEpoch(groups, myArgs.get('frameSecondaryInfo', {}))
    if len(groups) > 1:
        print(f'Frame break(s) detected — {len(groups)} virtual frames will be created')
    if myArgs['virtualFrame'] is not None:
        # Explicit override: assign sequentially from the given base number.
        baseVF = int(myArgs['virtualFrame'])
        groupAssignments = [(g, f'{baseVF + i:04d}') for i, g in enumerate(groups)]
    else:
        groupAssignments = assignVirtualFrameNumbers(groups, myArgs['orbit1'])

    for groupFrames, virtualFrame in groupAssignments:
        print(f'Virtual frame {virtualFrame}: frames {groupFrames}')
        myArgs['virtualFrame'] = virtualFrame
        myArgs['frames']      = groupFrames
        # ionosphereRangeOffsetCorrection is set per-group by globalFillIonosphere()
        # (or, without global fill, by createVirtualFrameRUNW()) -- but it's never
        # unset, so if a previous group set it and this group's own attempt fails
        # (e.g. globalFillIonosphere: no valid pixels) before reaching that point,
        # createVirtualFrameROFF would otherwise stamp this group's range.offsets.vrt
        # with the *previous* group's stale correction filename instead of leaving
        # the metadata correctly absent.
        myArgs.pop('ionosphereRangeOffsetCorrection', None)
        myArgs['geodat1']     = [perFrameGeodat[f][0] for f in groupFrames if f in perFrameGeodat]
        myArgs['geodat2']     = [perFrameGeodat[f][1] for f in groupFrames if f in perFrameGeodat]
        myArgs['pow']         = [perFramePow[f][0]    for f in groupFrames if f in perFramePow]
        myArgs['geodatpow']   = [perFramePow[f][1]    for f in groupFrames if f in perFramePow]
        # orbit2/secondaryDateTime/looks/bandwidth were set globally in
        # getSecondaryOrbit() from whichever frame's H5 was found first across the
        # WHOLE orbit -- refresh them to this group's own values so a split-out
        # group (see splitGroupsBySecondaryEpoch) writes its own correct pairinfo
        # and sensor YAML instead of inheriting another group's epoch/bandwidth
        # (a frame can be split out for a bandwidth/looks mismatch too, e.g.
        # track-58 frame 44 at 40 MHz/13 looks vs its 77 MHz/26-look neighbors).
        groupSecInfo = [myArgs['frameSecondaryInfo'][f] for f in groupFrames
                        if f in myArgs.get('frameSecondaryInfo', {})]
        if groupSecInfo:
            (myArgs['orbit2'], myArgs['secondaryDateTime'],
             myArgs['NumberRangeLooks'], myArgs['NumberAzimuthLooks'],
             myArgs['bandwidth']) = Counter(groupSecInfo).most_common(1)[0][0]

        frameDir = f'{myArgs["outputDir"]}/{myArgs["orbit1"]}_{virtualFrame}'
        if not os.path.exists(frameDir):
            os.mkdir(frameDir)
        writeFramesList(frameDir, groupFrames)

        # Create the virtual frame
        if haveData:
            if myArgs.get('corrOnly') or myArgs.get('correlationOnly'):
                copy_sensor_yaml(myArgs, frameDir)
                createVirtualFrameCorr(myArgs)
                writePairInfo(myArgs)
            else:
                copy_sensor_yaml(myArgs, frameDir)
                createVirtualFrameRUNW(myArgs)
                writePairInfo(myArgs)
                if not myArgs['RUNWOnly'] and not myArgs.get('geodatsOnly'):
                    createVirtualFrameROFF(myArgs)
                    if myArgs.get('globalFillIono'):
                        t_step = time.time()
                        print('\tGlobal iono fill....', end=' ', flush=True)
                        globalFillIonosphere(myArgs,
                                             sigmaAz=myArgs.get('sigmaAz') or 10.0,
                                             sigmaRng=myArgs.get('sigmaRg') or 30.0)
                        print(f'{time.time()-t_step:.1f}s  (total {time.time()-t0:.1f}s)')
        # Create a virtual frame if .pow images exist
        if not myArgs['RUNWOnly'] and not myArgs.get('correlationOnly') \
                and not myArgs.get('corrOnly'):
            createVirtualFramePower(myArgs)


if __name__ == "__main__":
    main()
