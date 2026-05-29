#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:34:37 2024

@author: ian
"""
import argparse
import os
import utilities as u
import nisarhdf
import sarfunc
from subprocess import call, DEVNULL
import sys
import numpy as np


def parseArgs():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mConvert RUNW to GrIMP formatted products '
        ' \033[0m\n\n',)

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
    choices = [x for x in sarfunc.defaultRegionDefs(None).regionsDef]
    parser.add_argument('--region', type=str, choices=choices, default=None,
                        help='region [autodetect greenland or antarctica]')
    parser.add_argument('--verbose', action='store_true',
                        help='Redirect all output to terminal for debugging')
    parser.add_argument('--interpThresh', type=int, default=20,
                        help='Maximum size hole to interpolate [20]')
    parser.add_argument('--islandThresh', type=int, default=20,
                        help='Maximum size isolated area to discard [20]')
    args = parser.parse_args()
    #
    print('...', args.RUNW[0])
    if not os.path.exists(args.RUNW[0]) and 's3' not in args.RUNW[0]:
        u.myerror(f'RUNW file {args.RUNW[0]} does not exist')
    #
    params = {}
    for param in ['outputDir', 'referenceXML', 'secondaryXML',
                  'referenceOrbit', 'secondaryOrbit', 'frame', 'region',
                  'simMask', 'simPhase', 'interpThresh', 'islandThresh']:
        params[param] = getattr(args, param)
    # Set ouput
    if args.verbose:
        params['stdout'], params['stderr'] = None, None
    else:
        params['stdout'], params['stderr'] = DEVNULL, DEVNULL
    return args.RUNW[0], params


def simPhase(geodat, params, dT, outputDir='.'):
    '''
    Simulate phse

    ----------
    geodat : str
        Geodat file name.
    params : dict
        Input params.
    outputDir : str, optional
        Path for output. The default is '.'.

    Returns
    -------
    None.

    '''
    print(f'Simphase {params["simPhase"]}')
    if not params['simPhase']:
        return
    regionDef = sarfunc.defaultRegionDefs(params['region'])
    #
    output = f'{outputDir}/phaseSim'
    # run command
    args = f"-velocity -dT {dT} {regionDef.dem()} {regionDef.velMap()}" \
        f" {outputDir}/{geodat} {output}"
    #
    command = 'siminsar'
    u.callMyProg(command, myArgs=args.split(), screen=True)
    return True


def runInterp(outputDir, inputFile, outputFile, nr, na, ratThresh=1,
              thresh=20, islandThresh=20,
              stderr=DEVNULL, stdout=DEVNULL, workingDir='workingDir'):
    '''
    Shell call to run interpolator
   ----------
    outputDir : str
        The product directory for final products.
    inputFile : str
        Filename to be interpolated
    outputFile : str
        Filename to be interpolated
    nr : int
        Number of range samples
    na : int
        Number of azimuth samples
    ratThresh : int, optional
        Allows larger skinny holes, best left at default. The default is 1.
    thresh : int, optional
        Fill only holes with area <=thresh. The default is 20.
    islandThresh : int, optional
        Remove isolated areas <= islandThresh pixels. The default is 20.
    stderr : file pointer, optional
        File for stdout output. Use None for stdout. The default is DEVNULL.
    stdout : file pointer, optional
        File for stderr output. Use None for stdout. The default is DEVNULL.
    layers : list, optional
        The layers to include. Only special cases need to use this. The
        default is [1, 2, 3].
    workingDir : str, optional
        The location for all of the intermediate outputs. The default is
        'workingDir'.

    Returns
    -------
    None.

    '''
    command = f'intfloat -wdist -nr {nr} -na {na} -thresh {thresh} ' \
        f'-islandThresh {islandThresh} {outputDir}/{inputFile} ' \
        f' > {outputDir}/{outputFile}'
    print(command)
    #
    # Run command
    # executable='/bin/csh',
    call(command, shell=True,  stderr=stderr, stdout=stdout)


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


def interpPhase(outputDir, myRUNW, ratThresh=1, thresh=20,
                islandThresh=20, stderr=DEVNULL, stdout=DEVNULL,
                workingDir='workingDir'):
    '''
    Call intfloat to do minor hole filling interpolation on phase

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    baseName : str
        Root name for offsets (e.g., offsets for offsets.layer1.cull.interp.da)
    myRUNW : nisarRUNW
        Unwrapped phase instance.
    ratThresh : int, optional
        Allows larger skinny holes, best left at default. The default is 1.
    thresh : int, optional
        Fill only holes with area <=thresh. The default is 20.
    islandThresh : int, optional
        Remove isolated areas <= islandThresh pixels. The default is 20.
    stderr : file pointer, optional
        File for stdout output. Use None for stdout. The default is DEVNULL.
    stdout : file pointer, optional
        File for stderr output. Use None for stdout. The default is DEVNULL.
    workingDir : str, optional
        The location for all of the intermediate outputs. The default is
        'workingDir'.

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
    # Save the masked version
    if hasattr(myRUNW, 'maskedUnwrappedPhase'):
        myRUNW.writeData(phaseFile, ['maskedUnwrappedPhase'], grimp=True,
                         tiff=False, byteOrder='MSB', noSuffix=True)
    else:
        myRUNW.writeData(phaseFile, 'unwrappedPhase',
                         grimp=True, byteOrder='MSB', noSuffix=True)
    #
    myRUNW.writeData(ionosphereFile, ['ionospherePhaseScreen'], grimp=True,
                     tiff=False, byteOrder='MSB', noSuffix=True)
    #
    corrFile = phaseFile.replace('.uw', '.cor')
    myRUNW.writeData(corrFile, ['coherenceMagnitude'], grimp=True,
                     tiff=False, byteOrder='MSB', noSuffix=True)
    # save geodats
    myRUNW.writeGeodatGeojson(filename=geodat1, path=outputDir,
                              secondary=False)
    myRUNW.writeGeodatGeojson(filename=geodat2, path=outputDir,
                              secondary=True)
    # Radians to meters for offset correction
    radiansToMeters = -0.5 * myRUNW.Wavelength/(2.0*np.pi)
    radiansToPixels = radiansToMeters / myRUNW.SLCRangePixelSize
    # Run the interpolation
    outputFile = os.path.basename(phaseFile).replace('nisar.uw',
                                                     'nisar.uw.interp')
    runInterp(outputDir, os.path.basename(phaseFile), outputFile,
              myRUNW.MLRangeSize, myRUNW.MLAzimuthSize, **kwargs)
    # runInterp(outputDir, os.path.basename(ionosphereFile),
    #          outputFile.replace('.uw', '.ion'),
    #          myRUNW.MLRangeSize, myRUNW.MLAzimuthSize, **kwargs)
    #
    # Write the correponding vrt
    myRUNW.assembleMeta()
    meta = myRUNW.meta.copy()
   
    meta['ByteOrder'] = 'MSB'
    interpFile = phaseFile+'.interp'
    rangeCorrectionFile = \
        os.path.basename(ionosphereCleanedFile.replace('.ion.filt',
                                                       '.ion.filt.rangeOffset.vrt'))
    rangeCorrectionFileUnfiltered = \
        os.path.basename(ionosphereFile.replace('.ion',
                                                 '.ion.unfilt.rangeOffset.vrt'))
    outputPhase = os.path.basename(f'{interpFile}.vrt')
    scales = [None, radiansToPixels, radiansToPixels]
    for inputFile, outputFile, description, scale in \
        zip([interpFile, ionosphereCleanedFile, ionosphereFile],
            [outputPhase, rangeCorrectionFile, rangeCorrectionFileUnfiltered],
            ['unwrappedPhase', 'rangeOffsetCorrection', 'RangeOffsetCorrection'],
            scales):
        # print(inputFile, outputFile)s
        meta['bands'] = [description]
        myGT = myRUNW.getGeoTransform(grimp=True, tiff=False)
        nisarhdf.writeMultiBandVrt(f'{outputDir}/{outputFile}',
                                   myRUNW.MLRangeSize,
                                   myRUNW.MLAzimuthSize,
                                   [inputFile],
                                   [description],
                                   geoTransform=myGT,
                                   tiff=False, metaData=meta,
                                   byteOrder=meta['ByteOrder'],
                                   scales=[scale])
    writePairInfo(myRUNW, outputDir)
    #
    if hasattr(myRUNW, 'ionosphereCleaned'):
        myRUNW.writeData(ionosphereCleanedFile, ['ionosphereCleaned'],
                         grimp=True,
                         tiff=False, byteOrder='MSB', noSuffix=True)


def simIceMask(geodat, params, outputDir='.'):
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

    Returns
    -------
    None.

    '''
    # Get the region definition, which contains mask, dem and other info.
    regionDef = sarfunc.defaultRegionDefs(params['region'])
    if not params['simMask'] or regionDef.icemask() is None:
        return False
    maskFile = f"{outputDir}/icemask"
    # run command
    args = f"-mask {regionDef.dem()} {regionDef.icemask()} {geodat} {maskFile}"
    command = 'siminsar'
    u.callMyProg(command, myArgs=args.split(), screen=True)
    return True


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
    if params['region'] is not None:
        return
    if myRUNW.epsg == 3031:
        params['region'] = 'antarctica'
        return
    elif myRUNW.epsg == 3413:
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
    # Instantiate RUNW class
    myRUNW = nisarhdf.nisarRUNWHDF(referenceOrbitXML=params['referenceXML'],
                                   secondaryOrbitXML=params['secondaryXML'])

    # Open and read the hdf
    myRUNW.openHDF(RUNW,
                   referenceOrbit=params['referenceOrbit'],
                   secondaryOrbit=params['secondaryOrbit'],
                   frame=params['frame'])
    #
    myRUNW.cleanIonosphere()
    #
    resolveRegion(myRUNW, params)
    #
    geodat = \
        f'geodat{myRUNW.NumberRangeLooks}x{myRUNW.NumberAzimuthLooks}.geojson'
    # Remove all but the largest connected component.
    myRUNW.maskPhase(largest=True)

    #
    # Apply ice mask if one exists, which removes coastal rocky areas
    if simIceMask(geodat, params, outputDir=params['outputDir']):
        myRUNW.applyMask(f"{params['outputDir']}/icemask")
    # Simulate the phase
    # print(params)
    simPhase(geodat, params, myRUNW.dT, outputDir=params['outputDir'])

    #
    if not os.path.exists(f'{params["outputDir"]}/{workingDir}'):
        os.mkdir(f'{params["outputDir"]}/{workingDir}')
    #
    # Interpolate and save the result as final version
    interpPhase(params['outputDir'], myRUNW, ratThresh=1,
                thresh=params['interpThresh'],
                islandThresh=params['islandThresh'], stderr=params['stderr'],
                stdout=params['stdout'], workingDir=workingDir)
    #


if __name__ == "__main__":
    main()
