#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:34:37 2024

@author: ian
"""
import argparse
import os
import time
import utilities as u
import nisarhdf
import sarfunc
from subprocess import call, DEVNULL
import sys
import numpy as np
import yaml
from nisargrimpworkflow.ROFFtoGrimp import updateSimVrtGeotransforms


def parseArgs():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mConvert RUNW to GrIMP formatted products '
        ' \033[0m\n\n',
        epilog='Part of the nisargrimpworkflow package.')

    parser.add_argument('RUNW', type=str, nargs=1,
                        help='RUNW hdf file to convert')
    parser.add_argument('--outputDir', type=str, default=None,
                        help='outputDir for results, defaults to RUNW path')
    parser.add_argument('--referenceXML', type=str, default=None,
                        help='Reference orbit xml')
    parser.add_argument('--secondaryXML', type=str, default=None,
                        help='Secondary orbit xml')
    parser.add_argument('--referenceOrbit', type=int, default=None,
                        help='Reference orbit: obsolete once embedded in hdf '
                        '[None]')
    parser.add_argument('--secondaryOrbit', type=int, default=None,
                        help='Secondary orbit: obsolete once embedded in hdf '
                        '[None]')
    parser.add_argument('--frame', type=int, default=None,
                        help='Frame: obsolete once embedded in hdf [None]')
    parser.add_argument('--simMask', action='store_true',
                        help='Create and apply mask')
    parser.add_argument('--simPhase', action='store_true',
                        help='Simulate phase')
    choices = [os.path.splitext(f)[0]
               for f in os.listdir(
                   os.path.join(os.path.dirname(sarfunc.__file__), 'regions'))
               if f.endswith('.yaml')]
    parser.add_argument('--region', type=str, choices=choices, default=None,
                        help='region [autodetect greenland or antarctica]')
    parser.add_argument('--regionFile', type=str, default=None,
                        help='YAML file with region-specific paths (velMap, DEM, '
                        'mask etc.). Overrides --region.')
    parser.add_argument('--verticalCorrection', type=str, default=None,
                        help='xyDEM grid (m/yr) of submergence/emergence '
                        'rate; passed to siminsar -verticalCorrection '
                        '[None]')
    parser.add_argument('--verbose', action='store_true',
                        help='Redirect all output to terminal for debugging')
    parser.add_argument('--interpThresh', type=int, default=20,
                        help='Maximum size hole to interpolate [20]')
    parser.add_argument('--islandThresh', type=int, default=20,
                        help='Maximum size isolated area to discard [20]')
    parser.add_argument('--ompThreads', type=int, default=4,
                        help='Number of OpenMP threads for siminsar [4]')
    parser.add_argument('--phaseDerivedIonosphere', action='store_true',
                        help='Write ionosphere layers (.ion, .ion.filt, and '
                        'range-correction VRTs). When omitted all ionosphere '
                        'outputs are skipped.')
    parser.add_argument('--noPhase', action='store_true',
                        help='Suppress unwrapped-phase output (.uw, '
                        '.uw.interp, .uw.interp.vrt)')
    parser.add_argument('--noIon', action='store_true',
                        help='Suppress ionosphere output (.ion, .ion.filt, '
                        'and range-correction VRTs). Ignored when '
                        '--phaseDerivedIonosphere is not set.')
    parser.add_argument('--minTol', type=float, default=None,
                        help='Variable smoothing-radius map (additional pass on top of '
                        'the fixed -thresh/-islandThresh interpolation): m/yr floor for '
                        'the adaptive tolerance clip(percentSpeed/100*speed, minTol, '
                        'maxTol). Required together with --percentSpeed/--maxTol to '
                        'enable the map [None]')
    parser.add_argument('--percentSpeed', type=float, default=None,
                        help='Variable smoothing-radius map: percent of local speed (e.g. '
                        '1 = 1%%) used in the adaptive tolerance. Required together with '
                        '--minTol/--maxTol [None]')
    parser.add_argument('--maxTol', type=float, default=None,
                        help='Variable smoothing-radius map: m/yr ceiling for the '
                        'adaptive tolerance. Required together with --minTol/'
                        '--percentSpeed [None]')
    parser.add_argument('--maxSmoothRadius', type=int, default=50,
                        help='Variable smoothing-radius map: sweep cap in single-look '
                        'pixels, clamped to <= 255 (byte output) [50]')
    parser.add_argument('--smoothNIter', type=int, default=3,
                        help='Variable smoothing-radius map: repeated box-filter passes '
                        'per sweep step (Gaussian-ish) [3]')
    parser.add_argument('--noVariableSmoothing', action='store_true', default=False,
                        help='Disable the variable smoothing-radius map even if '
                        '--minTol/--percentSpeed/--maxTol (or project.yaml) supply values')
    args = parser.parse_args()
    smoothFlags = [args.minTol is not None, args.percentSpeed is not None,
                  args.maxTol is not None]
    if any(smoothFlags) and not all(smoothFlags):
        u.myerror('RUNWtoGrimp: --minTol/--percentSpeed/--maxTol must be given together')
    if args.maxSmoothRadius > 255:
        print(f'WARNING: --maxSmoothRadius {args.maxSmoothRadius} exceeds byte range, '
              'clamping to 255')
        args.maxSmoothRadius = 255
    #
    print('...', args.RUNW[0])
    if not os.path.exists(args.RUNW[0]) and 's3' not in args.RUNW[0]:
        u.myerror(f'RUNW file {args.RUNW[0]} does not exist')
    #
    params = {}
    for param in ['outputDir', 'referenceXML', 'secondaryXML',
                  'referenceOrbit', 'secondaryOrbit', 'frame', 'region',
                  'regionFile', 'simMask', 'simPhase', 'interpThresh',
                  'islandThresh', 'ompThreads', 'phaseDerivedIonosphere',
                  'noPhase', 'noIon', 'verticalCorrection',
                  'minTol', 'percentSpeed', 'maxTol', 'maxSmoothRadius',
                  'smoothNIter', 'noVariableSmoothing']:
        params[param] = getattr(args, param)
    # Set ouput
    if args.verbose:
        params['stdout'], params['stderr'] = None, None
    else:
        params['stdout'], params['stderr'] = DEVNULL, DEVNULL
    return args.RUNW[0], params


def simPhase(geodat, params, dT, outputDir='.', ompThreads=4):
    '''
    Simulate phse

    ----------
    geodat : str
        Geodat file name.
    params : dict
        Input params.
    outputDir : str, optional
        Path for output. The default is '.'.
    ompThreads : int, optional
        Number of OpenMP threads for siminsar. The default is 4.

    Returns
    -------
    None.

    '''
    print(f'Simphase {params["simPhase"]}')
    if not params['simPhase']:
        return
    regionDef = sarfunc.defaultRegionDefs(params.get('region'),
                                          regionFile=params.get('regionFile'))
    #
    output = f'{outputDir}/phaseSim'
    # run command. siminsar requires its last 4 tokens to be demFile/displacementFile/
    # sceneFile/outputFile -- all optional flags must come before them, never after.
    args = f"-ompThreads {ompThreads} -velocity -dT {dT}"
    if params.get('verticalCorrection') is not None:
        args += f" -verticalCorrection {params['verticalCorrection']}"
    if params.get('minTol') is not None and not params.get('noVariableSmoothing'):
        args += f" -minTol {params['minTol']} -percentSpeed {params['percentSpeed']} " \
            f"-maxTol {params['maxTol']} -maxSmoothRadius {params['maxSmoothRadius']} " \
            f"-smoothNIter {params['smoothNIter']}"
    args += f" {regionDef.dem()} {regionDef.velMap()} {outputDir}/{geodat} {output}"
    #
    command = 'siminsar'
    u.callMyProg(command, myArgs=args.split(), screen=True)
    return True


def runInterp(outputDir, inputVRT, outputFile, ratThresh=1,
              thresh=20, islandThresh=20,
              stderr=DEVNULL, stdout=DEVNULL, workingDir='workingDir'):
    '''
    Shell call to run interpolator, reading geometry from a VRT.

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    inputVRT : str
        Basename of the VRT file for the image to be interpolated.  intfloat
        reads dimensions and georeferencing from the VRT; the matching binary
        is referenced inside it.
    outputFile : str
        Basename of the binary output file (written via stdout redirect).
    ratThresh : float, optional
        Allows larger skinny holes, best left at default. The default is 1.
    thresh : int, optional
        Fill only holes with area <=thresh. The default is 20.
    islandThresh : int, optional
        Remove isolated areas <= islandThresh pixels. The default is 20.
    stderr : file pointer, optional
        File for stderr output. Use None for terminal. The default is DEVNULL.
    stdout : file pointer, optional
        File for stdout output. Use None for terminal. The default is DEVNULL.
    workingDir : str, optional
        The location for all of the intermediate outputs. The default is
        'workingDir'.

    Returns
    -------
    None.

    '''
    command = (f'intfloat -wdist -inputVRT {outputDir}/{inputVRT} '
               f'-thresh {thresh} -islandThresh {islandThresh} '
               f'> {outputDir}/{outputFile}')
    print(command)
    # Run command
    call(command, shell=True, stderr=stderr, stdout=stdout)


def writePairInfo(myRUNW, outputDir):
    '''
    Write a pairinfo file

    Parameters
    ----------
    myRUN : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    pairInfoFile = \
        f'{outputDir}/{myRUNW.referenceOrbit}.{myRUNW.secondaryOrbit}.pairinfo'
    date1 = myRUNW.datetime.strftime("%Y-%m-%d")
    date2 = myRUNW.secondary.datetime.strftime("%Y-%m-%d")
    with open(pairInfoFile, 'w') as fp:
        print(f'{myRUNW.referenceOrbit} {myRUNW.secondaryOrbit} '
              f'{date1} {date2} '
              f'{myRUNW.NumberRangeLooks} {myRUNW.NumberAzimuthLooks}',
              file=fp)

# --- Usage ---


def interpPhase(outputDir, myRUNW, phaseDerivedIonosphere=False,
                noPhase=False, noIon=False,
                ratThresh=1, thresh=20,
                islandThresh=20, stderr=DEVNULL, stdout=DEVNULL,
                workingDir='workingDir', minTol=None, percentSpeed=None,
                maxTol=None, maxSmoothRadius=50, smoothNIter=3,
                noVariableSmoothing=False):
    '''
    Call intfloat to do minor hole filling interpolation on phase

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    myRUNW : nisarRUNW
        Unwrapped phase instance.
    phaseDerivedIonosphere : bool, optional
        Write ionosphere layers (.ion, .ion.filt, range-correction VRTs).
        The default is False.
    noPhase : bool, optional
        Suppress unwrapped-phase output (.uw, .uw.interp, .uw.interp.vrt).
        The default is False.
    noIon : bool, optional
        Suppress ionosphere output (.ion, .ion.filt, range-correction VRTs).
        Ignored when phaseDerivedIonosphere is False. The default is False.
    ratThresh : float, optional
        Allows larger skinny holes, best left at default. The default is 1.
    thresh : int, optional
        Fill only holes with area <=thresh. The default is 20.
    islandThresh : int, optional
        Remove isolated areas <= islandThresh pixels. The default is 20.
    stderr : file pointer, optional
        File for stderr output. Use None for terminal. The default is DEVNULL.
    stdout : file pointer, optional
        File for stdout output. Use None for terminal. The default is DEVNULL.
    workingDir : str, optional
        The location for all of the intermediate outputs. The default is
        'workingDir'.
    minTol, percentSpeed, maxTol, maxSmoothRadius, smoothNIter, noVariableSmoothing :
        Variable smoothing-radius map, applied to the interpolated phase as an additional
        pass on top of the fixed thresh/islandThresh interpolation above. See simPhase's
        siminsar -minTol/-percentSpeed/-maxTol for the tolerance formula. Requires
        simPhase() to have been run with the same minTol/percentSpeed/maxTol so that
        phaseSim.smr.vrt exists.

    Returns
    -------
    None.

    '''
    kwargs = {'ratThresh': ratThresh, 'thresh': thresh,
              'islandThresh': islandThresh, 'stderr': stderr, 'stdout': stdout,
              'workingDir': workingDir}
    #
    geodat1 = \
        f'geodat{myRUNW.NumberRangeLooks}x{myRUNW.NumberAzimuthLooks}.geojson'
    geodat2 = geodat1.replace('geojson', 'secondary.geojson')
    phaseFile = \
        f'{outputDir}/{myRUNW.referenceOrbit}_{myRUNW.frame}.' \
        f'{myRUNW.secondaryOrbit}_{myRUNW.frame}.' \
        f'{myRUNW.NumberRangeLooks}x{myRUNW.NumberAzimuthLooks}.nisar.uw'
    ionosphereFile = phaseFile.replace('.uw', '.ion')
    ionosphereCleanedFile = phaseFile.replace('.uw', '.ion.filt')
    corrFile = phaseFile.replace('.uw', '.cor')
    #
    # Assemble metadata and geotransform up-front (needed before runInterp
    # so we can write the VRT that intfloat reads for geometry)
    myRUNW.assembleMeta()
    meta = myRUNW.meta.copy()
    meta['ByteOrder'] = 'MSB'
    myGT = myRUNW.getGeoTransform(grimp=True, tiff=False)
    # Radians to meters for ionosphere offset correction
    radiansToMeters = -0.5 * myRUNW.Wavelength / (2.0 * np.pi)
    radiansToPixels = radiansToMeters / myRUNW.SLCRangePixelSize
    #
    # --- Unwrapped phase ---
    if not noPhase:
        # Save the masked binary
        if hasattr(myRUNW, 'maskedUnwrappedPhase'):
            myRUNW.writeData(phaseFile, ['maskedUnwrappedPhase'], grimp=True,
                             tiff=False, byteOrder='MSB', noSuffix=True)
        else:
            myRUNW.writeData(phaseFile, 'unwrappedPhase',
                             grimp=True, byteOrder='MSB', noSuffix=True)
        # Write VRT wrapper for the raw phase binary so intfloat can read
        # dimensions and georeferencing without -nr/-na flags
        uwVrt = os.path.basename(phaseFile) + '.vrt'
        meta['bands'] = ['unwrappedPhase']
        nisarhdf.writeMultiBandVrt(
            f'{outputDir}/{uwVrt}',
            myRUNW.MLRangeSize, myRUNW.MLAzimuthSize,
            [os.path.basename(phaseFile)], ['unwrappedPhase'],
            geoTransform=myGT, tiff=False, metaData=meta,
            byteOrder=meta['ByteOrder'], scales=[None])
    #
    # --- Ionosphere phase screen ---
    if phaseDerivedIonosphere and not noIon:
        myRUNW.writeData(ionosphereFile, ['ionospherePhaseScreen'], grimp=True,
                         tiff=False, byteOrder='MSB', noSuffix=True)
    #
    # Coherence is always written
    myRUNW.writeData(corrFile, ['coherenceMagnitude'], grimp=True,
                     tiff=False, byteOrder='MSB', noSuffix=True)
    # Save geodats
    myRUNW.writeGeodatGeojson(filename=geodat1, path=outputDir,
                              secondary=False)
    myRUNW.writeGeodatGeojson(filename=geodat2, path=outputDir,
                              secondary=True)
    #
    # --- Interpolation (intfloat reads VRT, writes binary to stdout) ---
    interpBasename = os.path.basename(phaseFile).replace('nisar.uw',
                                                         'nisar.uw.interp')
    if not noPhase:
        runInterp(outputDir, uwVrt, interpBasename, **kwargs)
    # runInterp(outputDir, os.path.basename(ionosphereFile) + '.vrt',
    #           interpBasename.replace('.uw', '.ion'), **kwargs)
    #
    # --- VRT wrappers for interp and ionosphere outputs ---
    interpFile = phaseFile + '.interp'
    rangeCorrectionFile = \
        os.path.basename(ionosphereCleanedFile.replace('.ion.filt',
                                                       '.ion.filt.rangeOffset.vrt'))
    rangeCorrectionFileUnfiltered = \
        os.path.basename(ionosphereFile.replace('.ion',
                                                '.ion.unfilt.rangeOffset.vrt'))
    outputPhase = os.path.basename(f'{interpFile}.vrt')
    #
    # Re-apply the connected-component=0 mask to the interpolated binary.
    # intfloat may have filled some cc=0 holes; those pixels should remain noData.
    if not noPhase and hasattr(myRUNW, 'connectedComponents'):
        arr = np.fromfile(interpFile, dtype='>f4').reshape(
            myRUNW.MLAzimuthSize, myRUNW.MLRangeSize)
        arr[myRUNW.connectedComponents < 1] = -2.0e9
        u.writeImage(interpFile, arr, '>f4')
    #
    # Variable smoothing-radius map (--minTol/--percentSpeed/--maxTol): an additional
    # smoothing pass on top of the fixed thresh/islandThresh interpolation above, applied
    # in place to the final interpolated phase before it's wrapped into outputPhase's VRT.
    if not noPhase and minTol is not None and not noVariableSmoothing:
        radiusMap = f'{outputDir}/phaseSim.smr.vrt'
        if not os.path.exists(radiusMap):
            print(f'WARNING: variable smoothing requested but {radiusMap} not found; '
                  'skipping (was simPhase run with the same --minTol/--percentSpeed/'
                  '--maxTol?)')
        else:
            print('Applying variable smoothing-radius map...')
            pixRatio = myRUNW.SLCAzimuthPixelSize / myRUNW.SLCRangePixelSize
            tmpFile = f'{interpFile}.vsmooth'
            command = f'filterfloat -nr {myRUNW.MLRangeSize} -na {myRUNW.MLAzimuthSize} ' \
                f'-radiusMap {radiusMap} -pixRatio {pixRatio} ' \
                f'-nIterations {smoothNIter} -minValue -2.0e9 ' \
                f'< {interpFile} > {tmpFile}'
            print(command)
            call(command, shell=True, stderr=stderr, stdout=stdout)
            os.replace(tmpFile, interpFile)
    #
    vrtInputs, vrtOutputs, vrtDescriptions, scales, noDataValues = \
        [], [], [], [], []
    if not noPhase:
        vrtInputs.append(interpFile)
        vrtOutputs.append(outputPhase)
        vrtDescriptions.append('Phase')
        scales.append(None)
        noDataValues.append(-2.0e9)
    if phaseDerivedIonosphere and not noIon:
        vrtInputs += [ionosphereCleanedFile, ionosphereFile]
        vrtOutputs += [rangeCorrectionFile, rangeCorrectionFileUnfiltered]
        vrtDescriptions += ['rangeOffsetCorrection', 'RangeOffsetCorrection']
        scales += [radiansToPixels, radiansToPixels]
        noDataValues += [-2.0e9, -2.0e9]
    for inF, outF, desc, scale, noData in \
            zip(vrtInputs, vrtOutputs, vrtDescriptions, scales, noDataValues):
        meta['bands'] = [desc]
        nisarhdf.writeMultiBandVrt(f'{outputDir}/{outF}',
                                   myRUNW.MLRangeSize,
                                   myRUNW.MLAzimuthSize,
                                   [inF], [desc],
                                   geoTransform=myGT,
                                   tiff=False, metaData=meta,
                                   byteOrder=meta['ByteOrder'],
                                   noDataValue=noData,
                                   scales=[scale])
    writePairInfo(myRUNW, outputDir)
    #
    if phaseDerivedIonosphere and not noIon \
            and hasattr(myRUNW, 'ionosphereCleaned'):
        myRUNW.writeData(ionosphereCleanedFile, ['ionosphereCleaned'],
                         grimp=True,
                         tiff=False, byteOrder='MSB', noSuffix=True)


def simIceMask(geodat, params, outputDir='.', ompThreads=4):
    '''
    Sim ice mask that can be used to retain phase only in ice covered areas
    Parameters
    ----------
    geodat : str
        Geodat file name.
    params : dict
        Input params.
    outputDir : str, optional
        Path for output. The default is '.'.
    ompThreads : int, optional
        Number of OpenMP threads for siminsar. The default is 4.

    Returns
    -------
    None.

    '''
    # Get the region definition, which contains mask, dem and other info.
    regionDef = sarfunc.defaultRegionDefs(params.get('region'),
                                          regionFile=params.get('regionFile'))
    if not params['simMask'] or regionDef.icemask() is None:
        return False
    maskFile = f"{outputDir}/icemask"
    # run command
    args = f"-ompThreads {ompThreads} -mask {regionDef.dem()} " \
        f"{regionDef.icemask()} {geodat} {maskFile}"
    command = 'siminsar'
    u.callMyProg(command, myArgs=args.split(), screen=True)
    return True


def _epsgFromProjectYaml(projectYaml='../project.yaml'):
    """Return EPSG from the project.yaml region file, or None if unavailable."""
    if not os.path.exists(projectYaml):
        return None
    try:
        with open(projectYaml) as fp:
            proj = yaml.safe_load(fp) or {}
        regionPath = proj.get('regionFile') or proj.get('region')
        if regionPath and os.path.isfile(str(regionPath)):
            return sarfunc.defaultRegionDefs(None, regionFile=str(regionPath)).epsg()
    except Exception:
        pass
    return None


def resolveRegion(myRUNW, params):
    '''
    If region not defined, then determine from epsg (greenland or
                                                         antarctica)

    Parameters
    ----------
    myRUNW : nisarRUNWHDF
        Current unwrapped instance.
    params : dict
        Parameters including region.

    Returns
    -------
    None.

    '''
    # regionFile provides all paths; no need to resolve a named region
    if params.get('regionFile') is not None:
        return
    if params['region'] is not None:
        return
    # Prefer the EPSG declared in the project.yaml region file over the EPSG
    # embedded in the HDF5, which can be wrong for frames with atypical coverage.
    epsg = _epsgFromProjectYaml() or myRUNW.epsg
    if epsg == 3031:
        params['region'] = 'antarctica'
        return
    elif epsg == 3413:
        params['region'] = 'greenland'
        return
    #
    print('Exited because could not resolve region from epsg')
    sys.exit()


def main():
    '''
    This program extracts the unwrapped phase from an RUNW HDF. It kills off
    all but the largest connected component and masks out any areas specified
    in an icemask file (typically bedrock)

    '''
    # Parse command line args
    workingDir = 'workingDir'
    RUNW, params = parseArgs()
    RUNWPath = os.path.dirname(RUNW)
    if len(RUNWPath) == 0:
        RUNWPath = '.'
    if params['outputDir'] is None:
        params['outputDir'] = RUNWPath
    print(params)
    t0 = time.perf_counter()
    # Instantiate RUNW class
    myRUNW = nisarhdf.nisarRUNWHDF(referenceOrbitXML=params['referenceXML'],
                                   secondaryOrbitXML=params['secondaryXML'])

    # Open and read the hdf
    myRUNW.openHDF(RUNW,
                   referenceOrbit=params['referenceOrbit'],
                   secondaryOrbit=params['secondaryOrbit'],
                   frame=params['frame'])
    t1 = time.perf_counter()
    print(f'[timing] openHDF: {t1 - t0:.1f}s')
    #
    if params['phaseDerivedIonosphere']:
        myRUNW.cleanIonosphere()
        t2 = time.perf_counter()
        print(f'[timing] cleanIonosphere: {t2 - t1:.1f}s')
        t1 = t2
    #
    resolveRegion(myRUNW, params)
    #
    geodat = \
        f'geodat{myRUNW.NumberRangeLooks}x{myRUNW.NumberAzimuthLooks}.geojson'
    # Remove all but the largest connected component.
    myRUNW.maskPhase(largest=True)
    t2 = time.perf_counter()
    print(f'[timing] maskPhase: {t2 - t1:.1f}s')
    t1 = t2

    #
    # Apply ice mask if one exists, which removes coastal rocky areas
    if simIceMask(geodat, params, outputDir=params['outputDir'],
                  ompThreads=params['ompThreads']):
        updateSimVrtGeotransforms(
            f'{params["outputDir"]}/icemask*.vrt', myRUNW)
        myRUNW.applyMask(f"{params['outputDir']}/icemask")
    t2 = time.perf_counter()
    print(f'[timing] simIceMask+applyMask: {t2 - t1:.1f}s')
    t1 = t2
    # Simulate the phase
    # print(params)
    simPhase(geodat, params, myRUNW.dT, outputDir=params['outputDir'],
             ompThreads=params['ompThreads'])
    updateSimVrtGeotransforms(
        f'{params["outputDir"]}/phaseSim*.vrt', myRUNW)
    t2 = time.perf_counter()
    print(f'[timing] simPhase: {t2 - t1:.1f}s')
    t1 = t2

    #
    if not os.path.exists(f'{params["outputDir"]}/{workingDir}'):
        os.mkdir(f'{params["outputDir"]}/{workingDir}')
    #
    # Interpolate and save the result as final version
    interpPhase(params['outputDir'], myRUNW,
                phaseDerivedIonosphere=params['phaseDerivedIonosphere'],
                noPhase=params['noPhase'],
                noIon=params['noIon'],
                ratThresh=1,
                thresh=params['interpThresh'],
                islandThresh=params['islandThresh'], stderr=params['stderr'],
                stdout=params['stdout'], workingDir=workingDir,
                minTol=params['minTol'], percentSpeed=params['percentSpeed'],
                maxTol=params['maxTol'], maxSmoothRadius=params['maxSmoothRadius'],
                smoothNIter=params['smoothNIter'],
                noVariableSmoothing=params['noVariableSmoothing'])
    t2 = time.perf_counter()
    print(f'[timing] interpPhase: {t2 - t1:.1f}s')
    print(f'[timing] total: {t2 - t0:.1f}s')
    #


if __name__ == "__main__":
    main()
