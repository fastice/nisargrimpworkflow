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
import sys
import importlib.resources
from osgeo import gdal

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
        'Processes all frames of a single orbit pair found in the current '
        'working directory, then consolidates the per-frame outputs into a '
        'single virtual-frame directory using GDAL VRT mosaics.  Must be run '
        'from the directory containing the <orbit1>_<frame> subdirectories '
        '(e.g. 12345_010/, 12345_020/, ...).',
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
    parser.add_argument('--correlationOnly', action='store_true',
                        help='Extract coherence and geodat files only. '
                        'Skips ROFF conversion, ionosphere estimation, and '
                        'virtual-frame assembly.')
    parser.add_argument('--corrOnly', action='store_true',
                        help='Extract .cor files per real frame and assemble '
                        'only the correlation virtual frame. Skips ROFF '
                        'conversion, ionosphere estimation, and power images.')
    parser.add_argument('--noGlobalFillIono', action='store_true',
                        help='Disable full-swath ionosphere gap fill; use '
                        'per-frame fill only (default: global fill is on)')
    parser.add_argument('--retainIntermediateIono', action='store_true',
                        help='Keep per-frame unfilled and per-frame filled '
                        'offset iono files after global fill (useful for '
                        'debugging or comparison)')
    args = parser.parse_args()
    #
    params = {}
    for key in ['overWrite', 'overWritePhase', 'firstFrame', 'lastFrame',
                'orbit1', 'allowMixedMode', 'virtualFrame', 'noMask',
                'verbose', 'RUNWOnly', 'ompThreads', 'phaseDerivedIonosphere',
                'outputAll', 'phaseThresh', 'correlationOnly', 'corrOnly',
                'retainIntermediateIono']:
        params[key] = getattr(args, key)
    params['globalFillIono'] = not args.noGlobalFillIono and bool(args.phaseDerivedIonosphere)
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


def getSecondaryOrbit(myArgs):
    '''
    Get first and second orbits
    '''
    for frame in myArgs['frames']:
        for product in ['RUNW', 'ROFF']:
            frameDir = f'{myArgs["orbit1"]}_{frame}'
            files = glob.glob(f'{frameDir}/H5/NISAR*{product}*.h5')
            if len(files) < 1:
                continue
            myProd = getattr(nisarhdf, f'nisar{product}HDF', None)()
            myProd.openHDF(files[0])
            if myProd.referenceOrbit == myArgs['orbit1']:
                myArgs['orbit2'] = myProd.secondaryOrbit
                myProd.getRangeBandWidth()
                myArgs['bandwidth'] = myProd.rangeBandwidth/1e6
                myArgs['NumberAzimuthLooks'] = myProd.NumberAzimuthLooks
                myArgs['NumberRangeLooks'] = myProd.NumberRangeLooks
                myArgs['myProd'] = myProd
                myArgs['datetime'] = myProd.Date
                myArgs['secondaryDateTime'] = myProd.secondaryDate
                return True
    return False


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
        if myArgs['noMask'] is True:
            command += ['--noMask']
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
    myArgs['geodat1'].append(
        f'{ruNWOutputDir}/geodat{nLooksR}x{nLooksA}.geojson')
    myArgs['geodat2'].append(
        f'{ruNWOutputDir}/geodat{nLooksR}x{nLooksA}.secondary.geojson')
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
        if myArgs['verbose'] is True:
            command += ['--verbose']
        if myArgs.get('phaseDerivedIonosphere'):
            command += ['--phaseDerivedIonosphere']
        else:
            command += ['--noPhase', '--noIon']
        # print(' '.join(command))
        run(command, stderr=myArgs['stderr'], stdout=myArgs['stdout'])
    else:
        print(f'skipping {orbit1}_{frame} since products exist')
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
    Sort and remove duplicates
    '''
    t = np.array(t, dtype=float)
    x = np.array(x, dtype=float)
    idx = np.argsort(t)
    t = t[idx]
    x = x[idx]
    unique_t, unique_idx = np.unique(t, return_index=True)
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


def mergedGeodat(geodatFiles, vrtFile):
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
    # Set size and range fields from actual VRT dimensions
    geodatMerged['properties']['MLAzimuthSize'] = na
    geodatMerged['properties']['MLRangeSize'] = nr
    geodatMerged['properties']['MLNearRange'] = mlNearRange
    geodatMerged['properties']['MLFarRange'] = mlFarRange
    geodatMerged['properties']['MLCenterRange'] = mlCenterRange
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

    # (glob suffix, use --offsets, label, save as ionosphereRangeOffsetCorrection)
    products = [
        ('*.correctedUnwrappedPhase.vrt',    True,  'correctedUnwrappedPhase',    False),
        ('*.cor.vrt',                         False, 'cor',                         False),
        ('*.ionosphereCorrection.vrt',        True,  'ionosphereCorrection',        False),
        ('*.ionosphereCorrection.offset.vrt', True,  'ionosphereCorrection.offset', True),
    ]
    if myArgs.get('globalFillIono'):
        # Global fill replaces the per-frame filled offset; skip building the
        # assembled offset VRT (its per-frame TIFs will be deleted), and add
        # the unfilled iono for the global-fill assembly step instead.
        products = [(g, o, l, r) for g, o, l, r in products
                    if l != 'ionosphereCorrection.offset']
        products.append(
            ('*.ionosphereCorrectionUnfilled.vrt', True,
             'ionosphereCorrectionUnfilled', False))
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
        # Virtual product name — use first frame (myFiles[0] belongs to frames[0])
        virtualProduct = os.path.basename(
            myFiles[0].replace(f'_{myArgs["frames"][0]}.', f'_{virtualFrame}.'))
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
    mergedGeodat(myArgs['geodat2'], virtualVRTs['correctedUnwrappedPhase'])


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
        # Virtual product name — use first frame (myFiles[0] belongs to frames[0])
        virtualProduct = os.path.basename(
            myFiles[0].replace(f'_{myArgs["frames"][0]}.', f'_{virtualFrame}.'))
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
    """Return (wavelength, slpSpacing) from the first available RUNW HDF5."""
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
            return float(runw.Wavelength), float(runw.SLCRangePixelSize)
        except Exception as e:
            u.mywarning(f'_readIonosphereParams: could not read {h5Files[0]}: {e}')
    return None, None


def globalFillIonosphere(myArgs, sigmaAz=100.0, sigmaRng=30.0):
    """Globally fill the assembled virtual-frame unfilled iono, re-partition back
    into per-frame tiles (in per-frame dirs), and build a plain geometry VRT.

    Unless --retainIntermediateIono is set, removes the per-frame unfilled and
    per-frame filled offset iono files, leaving only the globalFill products.
    """
    from nisargrimpworkflow.estimateIonosphere import (fill_and_smooth_iono,
                                                       write_geotiff,
                                                       write_output_vrt)
    orbit = myArgs['orbit1']
    virtualFrame = myArgs['virtualFrame']
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'

    # --- locate assembled unfilled iono VRT ---
    unfilledVrts = glob.glob(f'{frameDir}/*.ionosphereCorrectionUnfilled.vrt')
    if not unfilledVrts:
        u.mywarning('globalFillIonosphere: no unfilled iono VRT found, skipping')
        return
    unfilledVrt = unfilledVrts[0]

    # --- read consistent (bias-corrected) unfilled iono on full-swath RUNW grid ---
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

    # --- global fill on full swath ---
    print('globalFillIonosphere: running pyramid fill on full-swath unfilled iono...')
    filled = fill_and_smooth_iono(arr.astype(np.float32), valid,
                                  sigma=(sigmaAz, sigmaRng))

    # Hold the globally-filled raster in a GDAL MEM dataset — no temp file needed.
    nrows, ncols = filled.shape
    memDs = gdal.GetDriverByName('MEM').Create('', ncols, nrows, 1, gdal.GDT_Float32)
    memDs.SetGeoTransform(fullGt)
    memDs.SetProjection(fullProj)
    memDs.GetRasterBand(1).WriteArray(filled)
    memDs.GetRasterBand(1).SetNoDataValue(float('nan'))

    # --- wavelength and slpSpacing for radians → SLC pixels conversion ---
    wavelength, slpSpacing = _readIonosphereParams(myArgs)
    if wavelength is None or slpSpacing is None:
        u.mywarning('globalFillIonosphere: could not read wavelength/slpSpacing, skipping')
        return
    scale = -wavelength / (4.0 * np.pi * slpSpacing)

    # --- re-partition: warp to each frame's offset grid, write per-frame tiles ---
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

        # Derive output name from the per-frame ionosphereCorrection.offset.vrt
        existingOffVrts = glob.glob(f'{orbit}_{frame}/*.ionosphereCorrection.offset.vrt')
        if not existingOffVrts:
            u.mywarning(f'globalFillIonosphere: no offset vrt for {orbit}_{frame}')
            continue
        perFrameOffsetTif = existingOffVrts[0].replace(
            '.ionosphereCorrection.offset.vrt',
            '.ionosphereCorrection.globalFill.offset.tif')
        perFrameOffsetVrt = perFrameOffsetTif.replace('.tif', '.vrt')

        # Warp filled iono (radians, full-swath RUNW grid) to frame's offset grid
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

        # Apply radians → SLC pixels scale in-place
        ds = gdal.Open(perFrameOffsetTif, gdal.GA_Update)
        if ds is not None:
            arr2 = ds.GetRasterBand(1).ReadAsArray()
            mask2 = arr2 != -2.0e9
            arr2[mask2] *= scale
            ds.GetRasterBand(1).WriteArray(arr2)
            ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
            ds.FlushCache()
            ds = None

        write_output_vrt(perFrameOffsetVrt, [perFrameOffsetTif],
                         ['ionosphereCorrection'], roffGt)
        perFrameOffsetVrts.append(perFrameOffsetVrt)
        print(f'globalFillIonosphere: wrote {perFrameOffsetVrt}')

    if not perFrameOffsetVrts:
        u.mywarning('globalFillIonosphere: no per-frame tiles written, skipping VRT assembly')
        return

    # --- build plain geometry VRT in virtual-frame dir (no --offsets: biases baked in) ---
    virtualProduct = os.path.basename(perFrameOffsetVrts[0]).replace(
        f'_{myArgs["frames"][0]}.', f'_{virtualFrame}.')
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

    # --- clean up superseded and intermediate files (unless --retainIntermediateIono) ---
    if not myArgs.get('retainIntermediateIono'):
        for frame in myArgs['frames']:
            for pattern in ['*.ionosphereCorrectionUnfilled.tif',
                            '*.ionosphereCorrectionUnfilled.vrt',
                            '*.ionosphereCorrection.offset.tif']:
                for f in glob.glob(f'{orbit}_{frame}/{pattern}'):
                    os.remove(f)


def createVirtualFrameCorr(myArgs):
    '''Assemble per-frame .cor.vrt files into a single virtual-frame correlation mosaic.'''
    orbit = myArgs['orbit1']
    virtualFrame = myArgs['virtualFrame']
    frameDir = f'{myArgs["outputDir"]}/{orbit}_{virtualFrame}'
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
    virtualProduct = os.path.basename(
        myFiles[0].replace(f'_{myArgs["frames"][0]}.', f'_{virtualFrame}.'))
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
    mergedGeodat(myArgs['geodat2'], vrtFile)


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
    for key in myArgs:
        print(key, ':', myArgs[key])
    # Read optional regionFile from ../project.yaml
    myArgs['regionFile'] = None
    projectYaml = '../project.yaml'
    if os.path.exists(projectYaml):
        with open(projectYaml) as _fp:
            _proj = yaml.safe_load(_fp) or {}
        myArgs['regionFile'] = _proj.get('regionFile', None)
        if myArgs['regionFile']:
            print(f'regionFile from {projectYaml}: {myArgs["regionFile"]}')
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
                    and not myArgs.get('corrOnly'):
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
        if not myArgs.get('corrOnly'):
            np_before = len(myArgs['pow'])
            processFramePow(frame, myArgs)
            if len(myArgs['pow']) > np_before:
                perFramePow[frame] = (myArgs['pow'][-1], myArgs['geodatpow'][-1])

    # Split frames into contiguous groups; each gap creates a new virtual frame.
    allFrames = myArgs['frames']
    groups = splitFrameGroups(allFrames)
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
        myArgs['geodat1']     = [perFrameGeodat[f][0] for f in groupFrames if f in perFrameGeodat]
        myArgs['geodat2']     = [perFrameGeodat[f][1] for f in groupFrames if f in perFrameGeodat]
        myArgs['pow']         = [perFramePow[f][0]    for f in groupFrames if f in perFramePow]
        myArgs['geodatpow']   = [perFramePow[f][1]    for f in groupFrames if f in perFramePow]

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
                if not myArgs['RUNWOnly']:
                    createVirtualFrameROFF(myArgs)
                    if myArgs.get('globalFillIono'):
                        t_step = time.time()
                        print('\tGlobal iono fill....', end=' ', flush=True)
                        globalFillIonosphere(myArgs)
                        print(f'{time.time()-t_step:.1f}s  (total {time.time()-t0:.1f}s)')
        # Create a virtual frame if .pow images exist
        if not myArgs['RUNWOnly'] and not myArgs.get('correlationOnly') \
                and not myArgs.get('corrOnly'):
            createVirtualFramePower(myArgs)


if __name__ == "__main__":
    main()
