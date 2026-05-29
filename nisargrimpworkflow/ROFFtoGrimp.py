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
import glob
from subprocess import call, DEVNULL
import threading
from osgeo import gdal
import numpy as np
import rioxarray
import warnings
import sarfunc
import sys
from pathlib import Path


def parseArgs():
    '''
    Handle command line args

    Returns
    -------
    ROFF : str
        Filename for the ROFF HDF product.
    params : dict
        Dictionary with the parameters passed in from the command line.

    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mConvert ROFF to GrIMP formatted offsets '
        '\033[0m\n\n',)
    # default values
    boxSize = 7
    maxA, maxR = 3., 3.
    nGood = 10
    sa, sr = 3, 3
    parser.add_argument('ROFF', type=str, nargs=1,
                        help='ROFF hdf file to convert')
    parser.add_argument('--geodat1', type=str, default=None,
                        help='geodat file [geodatNLRxNLA.geojson]]')
    parser.add_argument('--geodat2', type=str, default=None,
                        help='geodat file [geodatNLRxNLA.secondary.geojson]]')
    parser.add_argument('--DEM', type=str, default=None,
                        help='Dem [default for region]')
    choices = [x for x in sarfunc.defaultRegionDefs(None).regionsDef]
    parser.add_argument('--region', type=str, choices=choices, default=None,
                        help='region [autodetect greenland or antarctica]')
    parser.add_argument('-regionFile', '--regionFile', type=str,
                        default=None, help='Yaml file with locations of '
                        'velMap, DEM etc for simulating offsests [None]')
    parser.add_argument('--correlationThresholds', type=float, nargs=3,
                        default=[0.07, 0.05, 0.025],
                        help='Correlation thresholds for discarding bad'
                        'matches for layers 1, 2, 3 respectively. '
                        'Three values must be specified.')
    parser.add_argument('--outputDir', type=str, default=None,
                        help='OutputDir for results, defaults to ROFF path')
    # box params
    parser.add_argument('--boxSize', type=int, default=boxSize,
                        help='Size of the box surrounding the point over '
                        f'which to calculate stats. [{boxSize}]')
    parser.add_argument('--nGood', type=float, default=nGood,
                        help='Reject points where there are fewer than nGood '
                        f'points in a box centered on the point. [{nGood}]')
    # Max deviations
    parser.add_argument('--maxR', type=float, default=maxR,
                        help='Max deviation from local median of the box for'
                        f' range offsets. [{maxR}]')
    parser.add_argument('--maxA', type=float, default=maxA,
                        help=f'Max deviation from local median of the box for'
                        f' azimuth offsets. [{maxA}]')
    # Smoothing params
    parser.add_argument('--sr', type=int, default=sr,
                        help='Smoothing length in range. For odd values the '
                        'kernel is uniform (e.g., 111) and shaped for even '
                        f'values (e.g., 0.51110.5) [{sr}]')
    #
    parser.add_argument('--sa', type=int, default=sa,
                        help='Smoothing length in azimuth. For odd values the '
                        'kernel is uniform (e.g., 111) and shaped for even '
                        f'values (e.g., 0.51110.5) [{sa}]')

    parser.add_argument('--interpThresh', type=int, default=20,
                        help='Maximum size hole to interpolate [20]')

    parser.add_argument('--islandThresh', type=int, default=20,
                        help='Maximum size isolated area to discard [20]')

    parser.add_argument('--byteOrder', type=str, default="MSB",
                        help='Byte order for outputs [MSB]')

    parser.add_argument('--verbose', action='store_true',
                        help='Display detail from all program calls')
    parser.add_argument('--noMask', action='store_true', default=False,
                        help='Do not apply mask to layer 3 in fast regions')
    parser.add_argument('--mergeOnly', action='store_true', default=False,
                        help='Do not regenerate data, merge only')

    choices = [x for x in sarfunc.defaultRegionDefs(None).regionsDef]

    args = parser.parse_args()
    if not os.path.exists(args.ROFF[0]) and 's3' not in args.ROFF[0]:
        u.myerror(f'ROFF file {args.ROFF[0]} does not exist')
    ROFFPath = os.path.dirname(args.ROFF[0])
    if len(ROFFPath) == 0:
        ROFFPath = '.'
    #
    params = {}
    if args.outputDir is None:
        params['outputDir'] = ROFFPath
    else:
        params['outputDir'] = args.outputDir
    #
    findGeodat(params, args.geodat1, args.geodat2)
    #
    if args.verbose:
        params['stdout'], params['stderr'] = None, None
    else:
        params['stdout'], params['stderr'] = DEVNULL, DEVNULL
    #
    for param in ['correlationThresholds', 'region', 'interpThresh',
                  'islandThresh', 'noMask', 'DEM', 'byteOrder', 'regionFile',
                  'mergeOnly']:
        params[param] = getattr(args, param)
    # Assemble cull params
    params['cullParams'] = {}
    for param in ['islandThresh', 'boxSize', 'nGood', 'maxR',
                  'maxA', 'sr', 'sa']:
        params['cullParams'][param] = getattr(args, param)
    #
    return args.ROFF[0], ROFFPath, params


def findGeodat(params, geodat1, geodat2):
    '''
    Find geodat names if either pass as None.
    Parameters
    ----------
    geodat1 : str
        Geodat1 name, use None to find name.
    geodat2 : str
        Geodat2 name, use None to find name.
    ROFFPath : str, optional
        Path to ROFF product. The default is '.'.

    Returns
    -------
    dict
        Dictionary with the geodat names.

    '''
    try:
        if geodat1 is None:
            tmp = glob.glob(f'{params["outputDir"]}/*.nisar.uw')
            nlr, nla, = os.path.basename(
                tmp[0]).split('.')[2].split('x')
            geodat1 = glob.glob(
                f'{params["outputDir"]}/geodat{nlr}x{nla}.geojson')[0]
        if geodat2 is None:
            geodat2 = geodat1.replace('.geojson', '.secondary.geojson')
    except Exception:
        print(f'geodat1 {geodat1}')
        print(f'geodat2 {geodat2}')
        u.myerror('Cannot find geodat files')
    params['geo1'] = geodat1
    params['geo2'] = geodat2


def setupGeodats(params):
    ''' Create links for geodat files

    Parameters
    ----------

      params : dict
        Dictionary with params passed in at the command line.
    '''
    for key in ['geo1', 'geo2']:
        geodatFile = os.path.basename(params[key])
        geodatPath = os.path.dirname(params[key])
        if not os.path.exists(f'{geodatPath}/workingDir/{geodatFile}'):
            symlink_file(params[key], f'{geodatPath}/workingDir/{geodatFile}')


def callSim(outputDir, baseName, params,
            includeVel=True, stderr=DEVNULL, stdout=DEVNULL,
            workingDir='workingDir'):
    '''
    Execute a shell command to run the offset simulations.
    includeVel : bool
        InlcudeVel : use velocity for the offset simulations. The default is
        True
    See simulateOffsets for other paremeter definitions.
    '''
    byteOrderFlag = {'MSB': '', 'LSB': '-LSB'}[params['byteOrder']]
    if params['regionFile'] is not None:
        regionArg = f'-regionFile {params["regionFile"]}'
    else:
        regionArg = f'-region={params["region"]}'
    command = f'simoffsets {regionArg} {byteOrderFlag} '
    if params['DEM'] is not None:
        command += f'-dem={params["DEM"] } '
    if includeVel is False:
        command += '-noVel '
    command += f'-offsetsDat={outputDir}/{workingDir}/{baseName}.dat '
    command += f'-azOffsets={outputDir}/{workingDir}/{baseName}.da -syncDat '
    geodat1 = os.path.basename(params["geo1"])
    command += f'-geodatFile={outputDir}/{workingDir}/{geodat1} '
    geodat2 = os.path.basename(params["geo2"])
    command += f'-secondGeodatFile={outputDir}/{workingDir}/{geodat2} '
    print(command)
    # , executable='/bin/csh'
    call(command, shell=True, stderr=stderr, stdout=stdout)



def symlink_file(src_path, dst_path, relative=True, overwrite=False):
    """
    Create a symbolic link to a file.

    This function creates a symbolic link at ``dst_path`` pointing to
    ``src_path``. By default the link target is stored as a relative path
    (relative to the destination directory), which makes directory trees
    more portable if they are moved.

    Parameters
    ----------
    src_path : str or pathlib.Path
        Path to the source file that the symbolic link will reference.
    dst_path : str or pathlib.Path
        Path where the symbolic link will be created.
    relative : bool, optional
        If True (default), create the symlink using a path relative to the
        destination directory. If False, use the absolute path to the source.
    overwrite : bool, optional
        If True, remove any existing file or symlink at ``dst_path`` before
        creating the new link. If False (default), an existing file or link
        will be left unchanged.

    Notes
    -----
    The source path is resolved to an absolute path before computing the
    symlink target. Relative links are generated using ``os.path.relpath``
    with respect to the destination directory.

    Returns
    -------
    None
        The function creates the symlink as a side effect.
    """
    src = Path(src_path)
    dst = Path(dst_path)

    src_abs = src.resolve()
    dst_parent = dst.parent.resolve()

    target = os.path.relpath(src_abs,
                             start=dst_parent) if relative else str(src_abs)

    if overwrite and (dst.exists() or dst.is_symlink()):
        dst.unlink()

    if not (dst.exists() or dst.is_symlink()):
        os.symlink(target, str(dst))


def simulateOffsets(outputDir, baseName, params,
                    stderr=DEVNULL, stdout=DEVNULL, workingDir='.'):
    '''
    Issue multithreaded shell calls to simulate offsets

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    baseName : str
        Root name for offsets (e.g., offsets for offsets.da)
    params : dict
        Dictionary with params passed in at the command line.
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
    print('Simulating offsets ...')
    threads = []
    # Geometry
    threads.append(threading.Thread(target=callSim,
                                    args=[outputDir, f'{baseName}.geom',
                                          params],
                                    kwargs={'includeVel': False,
                                            'stderr': stderr,
                                            'stdout': stdout,
                                            'workingDir': workingDir}))

    # Velocity with mask
    threads.append(threading.Thread(target=callSim,
                                    args=[outputDir, f'{baseName}', params],
                                    kwargs={'includeVel': True,
                                            'stderr': stderr,
                                            'stdout': stdout,
                                            'workingDir': workingDir}))
    quiet = False
    if stdout == DEVNULL:
        quiet = True
    u.runMyThreads(threads, 2, 'simoffsets', quiet=quiet)
    #
    offsetFiles = glob.glob(f'{outputDir}/offsets.*')
    # Add links to working dir
    for offsetFile in offsetFiles:
        symlink_file(offsetFile,
                     f'{outputDir}/workingDir/{os.path.basename(offsetFile)}',
                     relative=True, overwrite=False)


def runCull(outputDir, baseLayerName, boxSize=9, maxA=3, maxR=3, nGood=17,
            sa=3, sr=3, islandThresh=None, stderr=DEVNULL, stdout=DEVNULL,
            workingDir='workingDir'):
    '''
    Execute a shell command to run the culler.
    See cull st for paremeter definitions.
    '''

    command = f'cullst  -boxSize {boxSize}  -maxA {maxA}  -maxR {maxR} '
    command += f' -nGood {nGood} -sa {sa}  -sr {sr} '
    if islandThresh is not None:
        command += f'-islandThresh {islandThresh} '
    #
    command += f'{outputDir}/{workingDir}/{baseLayerName} '
    command += f'{outputDir}/{workingDir}/{baseLayerName}.cull'

    # Run command
    # , executable='/bin/csh'
    call(command, shell=True, stderr=stderr, stdout=stdout)


def cullst(outputDir, baseName, boxSize=7, maxA=3, maxR=3, nGood=10,
           sa=3, sr=3, islandThresh=20, layers=[1, 2, 3],
           stderr=DEVNULL, stdout=DEVNULL, workingDir='workingDir'):
    '''
    Call cullst to cull offsets

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    baseName : str
        Root name for offsets (e.g., offsets for offsets.layer1.cull.da)
    boxSize : int, optional
        Size of the box surrounding the point over which to calculate stats.
        The default is 7.
    maxA : int, optional
        Max deviation from local median. The default is 3.
    maxR : int, optional
        max deviation from local median. The default is 3.
    nGood : str, optional
        Reject points where there are few than nGood points in the box. The \
        default is 10.
    sa : int, optional
        Azimuth smoothing length (for odd kernel is uniform (e.g. 3->111)
        for even kernel is shape (e.g.,4->.5111.5). The default is 3.
    sr : int, optional
        Range smoothing length. The default is 3.
    islandThresh : int, optional
        Remove isolated areas <= islandThresh pixels. The default is 20.
    layers : list, optional
        The layers to include. Only special cases need to use this. The
        default is [1, 2, 3].
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
    print('Culling offsets...')
    kwargs = {'boxSize': boxSize, 'maxA': maxA, 'maxR': maxR, 'nGood': nGood,
              'sa': sa, 'sr': sr, 'islandThresh': islandThresh,
              'stderr': stderr, 'stdout': stdout, 'workingDir': workingDir}
    threads = []
    for layer in layers:
        baseLayerName = f'{baseName}.layer{layer}'
        # print(baseLayerName)
        threads.append(threading.Thread(target=runCull,
                                        args=[outputDir, baseLayerName],
                                        kwargs=kwargs))
    quiet = False
    if stdout == DEVNULL:
        quiet = True
    u.runMyThreads(threads, len(layers), 'culling', quiet=quiet)


def writeInterpVrt(newVRTFile, sourceFiles, descriptions, nr, na,
                   byteOrder=None, eType=gdal.GDT_Float32,
                   geoTransform=[-0.5, 1., 0., -0.5, 0., 1.], metaData=None,
                   noDataValue=-2.0e9):
    '''
    Write a vrt for the file. Note sourcefiles and descriptions have
    to be passed in.

    Parameters
    ----------
    newVRTFile : str
        Name for vrt file.
    sourceFiles : list
        Source files to include in the vrt.
    descriptions : list of str
        Descriptions for each layer to be included in vrt.
    nr : int
        Number of range samples
    na : int
        Number of azimuth samples
    byteOrder : str, optional
        Byte order (MSB or LSB). The default is None, which defaults to MSB
        if not in metaData.
    eType : data type, optional
        The data type. The default is gdal.GDT_Float32.
    geoTransform : list, optional
        The geotransform. The default is [-0.5, 1., 0., -0.5, 0., 1.].
    metaData : dict, optional
        Dict with optional metadata. The default is None.
    noDataValue : same as eType, optional
        The no data value. The default is -2.0e9.

    Returns
    -------
    None.

    '''
    # Make sure source files and descriptions are lists and metaData is dict
    if type(sourceFiles) is not list:
        sourceFiles = [sourceFiles]
    if type(descriptions) is not list:
        descriptions = [descriptions]
    if metaData is None:
        metaData = {}
    #
    # Kill any old file
    if os.path.exists(newVRTFile):
        os.remove(newVRTFile)
    # Create VRT
    bands = len(sourceFiles)
    drv = gdal.GetDriverByName("VRT")
    vrt = drv.Create(newVRTFile, nr, na, bands=0, eType=eType)
    vrt.SetGeoTransform(geoTransform)
    # Set the byte order
    #
    if byteOrder is None:
        if "ByteOrder" in metaData:
            byteOrder = metaData["ByteOrder"]
        else:
            byteOrder = "MSB"
            metaData["ByteOrder"] = byteOrder
    else:
        metaData["ByteOrder"] = byteOrder
    #
    vrt.SetMetadata(metaData)
    # Loop to add bands
    for sourceFile, description, bandNumber in \
            zip(sourceFiles, descriptions, range(1, bands + 1)):
        # Setup options
        options = [f"SourceFilename={sourceFile}", "relativeToVRT=1",
                   "subclass=VRTRawRasterBand", f"BYTEORDER={byteOrder}",
                   bytes(0)]
        # add the new band
        vrt.AddBand(eType, options=options)
        # Set band properties
        band = vrt.GetRasterBand(bandNumber)
        band.SetMetadataItem("Description", description)
        band.SetDescription(description)
        band.SetNoDataValue(noDataValue)
    # Close the vrt
    vrt = None


def runInterp(outputDir, inputFile, outputFile, nr, na, ratThresh=1,
              thresh=20, islandThresh=20, byteOrder='MSB',
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
    byteOrder : str, optional
        Set to MSB or LSB for interp output. The default is MSB.
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
    byteOrderFlag = ''
    if byteOrder == 'LSB':
        byteOrderFlag = '-LSB'
    command = f'intfloat -wdist -nr {nr} -na {na} -thresh {thresh} ' \
        f'{byteOrderFlag} ' \
        f'-islandThresh {islandThresh} {outputDir}/{workingDir}/{inputFile} ' \
        f' > {outputDir}/{workingDir}/{outputFile}'
    #
    # Run command
    # , executable='/bin/csh'
    call(command, shell=True, stderr=stderr, stdout=stdout)


def interpOffsets(outputDir, baseName, ratThresh=1, thresh=20, islandThresh=20,
                  stderr=DEVNULL, stdout=DEVNULL, layers=[1, 2, 3],
                  workingDir='workingDir', byteOrder='MSB'):
    '''
    Call intfloat to interpolate offsets to do minor hole filling.

    Parameters
    ----------
    outputDir : str
        The product directory for final products.
    baseName : str
        Root name for offsets (e.g., offsets for offsets.layer1.cull.interp.da)
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
    print('Interpolating...')
    datPath = f'{outputDir}/{workingDir}/{baseName}.layer?.dat'
    datFiles = glob.glob(datPath)
    if len(datFiles) > 0:
        with open(datFiles[0], 'r') as fp:
            line = fp.readlines()[0].split()
            nr, na = int(line[2]), int(line[3])
    else:
        u.myerror(f'Cannot open {datPath})')
    kwargs = {'ratThresh': ratThresh, 'thresh': thresh,
              'islandThresh': islandThresh, 'stderr': stderr, 'stdout': stdout,
              'workingDir': workingDir, 'byteOrder': byteOrder}
    threads = []
    # for historical resasons sr, sa don't have cull in the name
    inputSuffixes = ['cull.dr', 'cull.da', 'sr', 'sa']
    # Make all output have consistent names
    outputSuffixes = ['cull.interp.dr', 'cull.interp.da', 'cull.interp.sr',
                      'cull.interp.sa']
    #
    for inputSuffix, outputSuffix in zip(inputSuffixes, outputSuffixes):
        for layer in layers:
            inputFile = f'{baseName}.layer{layer}.{inputSuffix}'
            outputFile = f'{baseName}.layer{layer}.{outputSuffix}'
            threads.append(
                threading.Thread(target=runInterp,
                                 args=[outputDir, inputFile, outputFile,
                                       nr, na],
                                 kwargs=kwargs))
    # Run the interpolation
    quiet = False
    if stdout == DEVNULL:
        quiet = True
    u.runMyThreads(threads, len(layers) * 4, 'interpolating', quiet=quiet)
    #
    # Write the correponding vrt
    descriptions = ['RangeOffsets', 'AzimuthOffsets', 'RangeSigma',
                    'AzimuthSigma']
    for layer in layers:
        sourceFiles = [f'{baseName}.layer{layer}.cull.interp{suffix}'
                       for suffix in ['.dr', '.da', '.sr', '.sa']]
        #
        vrtFile = \
            f'{outputDir}/{workingDir}/{baseName}.layer{layer}.cull.interp.vrt'
        writeInterpVrt(vrtFile,
                       sourceFiles, descriptions, nr, na,
                       byteOrder=byteOrder, eType=gdal.GDT_Float32,
                       geoTransform=[-0.5, 1., 0., -0.5, 0., 1.],
                       metaData=None,
                       noDataValue=-2.0e9)


def readVRTAndRenameBands(layerVRT, masked=True, nameKey='long_name'):
    '''
    Read a vrt as rioxarray. Rename bands to match namekey

    Parameters
    ----------
    layerVRT : str
        Vrt with offset layers.
    masked : bool, optional
        Masks the data to replace nodata values with nans. The default is True.
    nameKey : str optional
        Key to use for band name (e.g., long_name or Description). The default
        is 'long_name'.

    Returns
    -------
    xarray
        Data in vrt with the bands names using nameKey.

    '''
    myOffsets = rioxarray.open_rasterio(layerVRT,
                                        band_as_variable=True,
                                        masked=True)
    #
    try:
        # extract band names
        originalBandNames = [a for a in myOffsets.data_vars]
        # Get new band names
        newBandNames = [getattr(a[1], nameKey)
                        for a in myOffsets.data_vars.items()]
        # return renamed xarray with newBandNames instead of original
        return myOffsets.rename(dict(zip(originalBandNames, newBandNames)))
    # Handle errors
    except Exception as errMsg:
        u.myerror(f'{errMsg}\n readAndRenameBands: error renaming band '
                  f'check {nameKey} present in vrt')


def readVRTAndAppend(layerVRT, data):
    '''
    Read a vrt and append to list in a dictionary indexed by bandname.
    This is used to read multiple vrts and stack the results.

    Parameters
    ----------
    layerVRT : str
        VRT filename.
    data : dict
        A dict to stack the data (eg. {'band1': [], 'band2': []...}).

    Returns
    -------
    data : dict
        Dictionary with new data appended to lists in dict.

    '''
    # Read the vrt and rename bands
    myOffsets = readVRTAndRenameBands(layerVRT, masked=True)
    # Append data to list for each band
    for band in myOffsets.data_vars:
        data[band].append(myOffsets[band].data)
    return data


def mergeOffsets(outputDir, baseName='NISARoffsets', layers=[1, 2, 3],
                 noData=-2.e9, byteOrder='MSB', simName='offsets',
                 simDir='.'):
    '''
    Merge offsets by averaging the results from the three layers.

    Parameters
    ----------
    outputDir : str
        Path where the product is located.
    baseName : str, optional
        Default root name for the products. The default is 'NISARoffsets'.
    layers : list optional
        List of layers. The default is [1, 2, 3].
    noData : str, optional
        No data value. The default is -2.e9.

    Returns
    -------
    None.
    '''
    print('Merging offsets...')
    np.seterr(all='ignore')
    # Create dict with list for each band to save layer data
    bands = ['RangeOffsets', 'AzimuthOffsets', 'RangeSigma', 'AzimuthSigma']
    data = dict(zip(bands, [[], [], [], []]))
    #
    # Loop over layers
    for layer in layers:
        # Current layer file
        layerVRT = \
            f'{outputDir}/workingDir/{baseName}.layer{layer}.cull.interp.vrt'
        # read the offsets and save result in a list
        data = readVRTAndAppend(layerVRT, data)
    #
    # Compute means
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        # Keep track of valid data
        validData = np.isfinite(np.stack(data['RangeOffsets']))
        # Compute mean of offset layers
        rgMean = np.nanmean(np.stack(data['RangeOffsets']), axis=0)
        azMean = np.nanmean(np.stack(data['AzimuthOffsets']), axis=0)
        # Read geometry
        geomVRT = f'{outputDir}/{simDir}/{simName}.geom.vrt'
        geom = readVRTAndRenameBands(geomVRT, nameKey='Description')
        print(geomVRT)
        # Add geometry back to valid pixels
        validMean = np.isfinite(azMean)
        print(np.sum(validMean), np.nanmean(rgMean[validMean]))
        #rgMean[validMean] += geom.RangeOffsets.data[validMean]
        rgMean[validMean] = (rgMean[validMean] +
                             geom.RangeOffsets.data[validMean])
        azMean[validMean] = (azMean[validMean] +
                             geom.AzimuthOffsets.data[validMean])
        # multiply sigmas by N/sqrt(N) = sqrt(N), set nodata to 1 to not
        # rescale no data value
        sqrtN = np.sqrt(np.sum(validData, axis=0))
        sqrtN[sqrtN == 0] = 1
        rgSigmaMean = np.nanmean(np.stack(data['RangeSigma']),
                                 axis=0) * sqrtN
        azSigmaMean = np.nanmean(np.stack(data['AzimuthSigma']),
                                 axis=0) * sqrtN
        # Use grimp noData value
        for x in [rgMean, azMean, rgSigmaMean, azSigmaMean]:
            x[~validMean] = noData
    # Save data
    filenames = ['range.offsets', 'azimuth.offsets',
                 'range.offsets.sr', 'azimuth.offsets.sa']
    dataType = {'LSB': 'f4', 'MSB': '>f4'}[byteOrder]
    for filename, var in zip(filenames,
                             [rgMean, azMean, rgSigmaMean, azSigmaMean]):
        u.writeImage(f'{outputDir}/{filename}', var, dataType)


def resolveRegion(myROFF, params):
    '''
    If region not defined, then determine from epsg (greenland or
                                                         antarctica)

    Parameters
    ----------
    myROFF: nisarRUNWHDF
        Current unwrapped instance.
    params : dict
        Parameters including region.

    Returns
    -------
    None.

    '''
    if params['region'] is not None:
        return
    if myROFF.epsg == 3031:
        params['region'] = 'antarctica'
        return
    elif myROFF.epsg == 3413:
        params['region'] = 'greenland'
        return
    #
    print('Exited because could not resolve region from epsg')
    sys.exit()


def writeVRTs(myROFF, ROFFPath, params):
    '''
    Write the VRT's for final product

    Parameters
    ----------
    myROFF : nisarROFFHDF
        Offset.
    ROFFPath : str
        Path to ROFF.
    params : dict
        Dictionary of params.

    Returns
    -------
    None.

    '''
    #
    # Assemble metadata
    print('Writing final vrts...')
    metaData = {}
    for key in params:
        if 'cullParams' in key:
            for key1 in params[key]:
                metaData[key1] = params[key][key1]
        else:
            if key not in ['stderr', 'stdout', 'byteOrder']:
                metaData[key] = params[key]
    for var in ['r0', 'a0', 'deltaR', 'deltaA']:
        metaData[var] = getattr(myROFF, var)
    metaData['sigmaStreaks'] = 0.0
    metaData['sigmaRange'] = 0.0
    metaData['geo1'] = os.path.basename(metaData['geo1'])
    metaData['geo2'] = os.path.basename(metaData['geo2'])
    azFiles = ['azimuth.offsets', 'azimuth.offsets.sa']
    azDescriptions = ['AzimuthOffsets', 'AzimuthSigma']
    rgFiles = ['range.offsets', 'range.offsets.sr']
    rgDescriptions = ['RangeOffsets', 'RangeSigma']
    filenames = ['azimuth.offsets.vrt', 'range.offsets.vrt',
                 'offsets.range-azimuth.vrt']
    print(metaData)
    for filename, files, descriptions in \
        zip(filenames,
            [azFiles, rgFiles, rgFiles + azFiles],
            [azDescriptions, rgDescriptions, rgDescriptions + azDescriptions]):

        # Remove any prior file
        if os.path.exists(f'{ROFFPath}/{filename}'):
            os.remove(f'{ROFFPath}/{filename}')
        # Write new file
        writeInterpVrt(f'{ROFFPath}/{filename}',
                       files,
                       descriptions,
                       myROFF.OffsetRangeSize, myROFF.OffsetAzimuthSize,
                       eType=gdal.GDT_Float32,
                       geoTransform=myROFF.getGeoTransform(grimp=True,
                                                           tiff=False),
                       byteOrder=params['byteOrder'],
                       metaData=metaData)


def main():
    '''
    Command line program to clean and merge the offsets from an ROFF product.

    The program breaks out the individual layers and culsl and interpolates
    them using call to stand alone c programs. Thevintermediate products are
    saved in 'workingDir' just below the directory containing the ROFF product.

    '''
    gdal.UseExceptions()
    # Parse command line args
    ROFF, ROFFPath, params = parseArgs()
    print(params)
    #
    # Setup ROFF object and open HDF file
    myROFF = nisarhdf.nisarROFFHDF()
    myROFF.openHDF(ROFF)
    resolveRegion(myROFF, params)
    #
    # Discard outliers based on correlation peak
    myROFF.removeOutlierOffsets('correlationSurfacePeak',
                                thresholds=params['correlationThresholds'])
    #
    if not os.path.exists(f'{params["outputDir"]}/workingDir'):
        os.mkdir(f'{params["outputDir"]}/workingDir')
    #
    # Write initial dat files for simulations
    myROFF.writeOffsetsDatFile(f'{params["outputDir"]}/offsets.dat',
                               geodat1=params['geo1'])
    myROFF.writeOffsetsDatFile(
        f'{params["outputDir"]}/offsets.geom.dat',
        geodat1=params['geo1'])
    #
    # Setup the the geodats
    setupGeodats(params)
    #
    # Simulate the offsets
    print(params)
    print(params['mergeOnly'])
    if not params['mergeOnly']:
        simulateOffsets(params["outputDir"], 'offsets', params,
                        stdout=params['stdout'], stderr=params['stderr'],
                        workingDir='.')
        #
        # Apply any mask files
        if not params['noMask']:
            myROFF.applyMask(f'{params["outputDir"]}/workingDir/offsets.mask.vrt')
    #
    # Save the offsets to layer files
    
        myROFF.writeData(f'{params["outputDir"]}/workingDir/NISARoffsets',
                         bands=['slantRangeOffset',
                                'alongTrackOffset',
                                'correlationSurfacePeak'],
                         tiff=False,
                         byteOrder=params['byteOrder'],
                         grimp=True,
                         saveMatch=True,
                         scaleToPixels=True)
        #
        # Cull the offset layers
        cullst(params["outputDir"], 'NISARoffsets',
               stdout=params['stdout'], stderr=params['stderr'],
               **params['cullParams'])
        #
        # Interpolate the offset layers
        interpOffsets(params["outputDir"], 'NISARoffsets',
                      ratThresh=1,
                      thresh=params['interpThresh'],
                      islandThresh=params['islandThresh'],
                      byteOrder=params['byteOrder'],
                      stdout=params['stdout'], stderr=params['stderr'])
    #
    mergeOffsets(params["outputDir"],
                 baseName='NISARoffsets',
                 byteOrder=params['byteOrder'],
                 simName='offsets')
    #
    writeVRTs(myROFF, params["outputDir"], params)
    #
    # writeVRTsRD(myROFF, ROFFPath, params)


if __name__ == "__main__":
    main()
