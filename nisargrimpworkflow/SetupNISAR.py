#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 12:07:17 2026

@author: ian
"""
import utilities as u
import argparse
import nisarhdf
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
               '  SetupNISAR 12345 --allowMixedMode --verbose\n',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('orbit1', type=int,
                        help='Reference orbit number')

    parser.add_argument('--virtualFrame', type=str, default='0000',
                        help='Frame suffix for the consolidated virtual-frame '
                        'output directory (e.g. 0000 → <orbit1>_0000/) [0000]')
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
    args = parser.parse_args()
    #
    params = {}
    for key in ['overWrite', 'overWritePhase', 'firstFrame', 'lastFrame',
                'orbit1', 'allowMixedMode', 'virtualFrame', 'noMask',
                'verbose', 'RUNWOnly', 'ompThreads', 'phaseDerivedIonosphere',
                'outputAll', 'phaseThresh', 'correlationOnly']:
        params[key] = getattr(args, key)
    #
    if args.verbose:
        params['stdout'], params['stderr'] = None, None
    else:
        params['stdout'], params['stderr'] = DEVNULL, DEVNULL
    #
    return params


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
    # Load godats
    geodats = {}
    lastFrame = None
    nrs = []
    for geoFile in geodatFiles:
        orbit, frame = geoFile.split('/')[-2].split('_')
        with open(geoFile) as fp:
            geodats[frame] = geojson.load(fp)
        nrs.append(geodats[frame]['properties']['MLRangeSize'])
        lastFrame = frame
    # Create the merged geodat as a copy of the first
    geodatMerged = copy.deepcopy(geodats[f'{firstFrame}'])
    # Merge the state vectors
    mergeStateVectors(geodats, geodatMerged)
    # Merge the corner coords
    geodatMerged['geometry']['coordinates'][0] = \
        mergeCorners(geodats[f'{firstFrame}'], geodats[f'{lastFrame}'])
    #
    geodatMerged['properties']['MLAzimuthSize'] = na
    # Assumes all files to be merged have the same number of range lines
    if len(set(nrs)) != 1:
        print(geodatMerged['properties']['MLRangeSize'], nrs)
        print('Warning RANGE SIZE INCONSISTENCY')
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
    # Process frames
    #
    for frame in myArgs['frames']:
        # Wrap all H5 files in the frame directory with VRT
        wrapH5sInFrameDir(myArgs['orbit1'], frame, verbose=myArgs['verbose'])
        if haveData:
            print(f'Processing Frame {frame}...')
            print('\tRUNW....')
            mixedMode = processFrameRUNW(frame, myArgs)
            if not mixedMode and not myArgs['RUNWOnly'] and not myArgs.get('correlationOnly'):
                print('\tROFF....')
                processFrameROFF(frame, myArgs)
                if not myArgs.get('phaseDerivedIonosphere'):
                    print('\tIonosphere (estimateIonosphere)....')
                    processFrameIonosphere(frame, myArgs, simDir='simPhase')
        # Need to add check for power mixed mode
        processFramePow(frame, myArgs)
        
    #
    virtualFrame = myArgs["virtualFrame"]
    frameDir = f'{myArgs["outputDir"]}/{myArgs["orbit1"]}_{virtualFrame}'
    #print(frameDir)
    #u.myerror('stop')
    if not os.path.exists(frameDir):
        os.mkdir(frameDir)
    
    # Create the virtual frame
    if haveData and not myArgs.get('correlationOnly'):
        # Copy sensor YAML into frameDir if bandwidth is present and file missing
        copy_sensor_yaml(myArgs, frameDir)
        createVirtualFrameRUNW(myArgs)
        writePairInfo(myArgs)
        #u.myerror('stop debug')  # DEBUG: exit after first frame
        if not myArgs['RUNWOnly']:
            createVirtualFrameROFF(myArgs)
    # Create a virtual frame if .pow images exist
    if not myArgs['RUNWOnly'] and not myArgs.get('correlationOnly'):
        createVirtualFramePower(myArgs)


if __name__ == "__main__":
    main()
