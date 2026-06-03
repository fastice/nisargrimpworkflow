#!/usr/bin/env python3
'''
processFrameIonosphere.py

Run estimateIonosphere for a single NISAR frame when --phaseDerivedIonosphere
is not selected.

Intended to be called from SetupNISAR.main() after processFrameROFF completes:

    from nisargrimpworkflow.processFrameIonosphere import processFrameIonosphere
    ...
    if not mixedMode and not myArgs['RUNWOnly']:
        processFrameROFF(frame, myArgs)
        if not myArgs.get('phaseDerivedIonosphere'):
            processFrameIonosphere(frame, myArgs)
'''
import glob
import os
import sys
from subprocess import run, DEVNULL

# Use the running interpreter explicitly so we don't accidentally pick up a
# different Python from PATH (e.g. a conda 3.10 binary when we're on 3.12).
_PYTHON = sys.executable


def _find_file(base, candidates):
    '''Return the first existing path from candidates (relative to base), or None.'''
    for rel in candidates:
        path = os.path.join(base, rel)
        if os.path.exists(path):
            return path
    return None


def processFrameIonosphere(frame, myArgs, simDir='simPhase'):
    '''
    Run estimateIonosphere for a single frame.

    The output stem follows the RUNWtoGrimp naming convention:
        {orbit1}_{frame}.{orbit2}_{frame}.{nLooksR}x{nLooksA}.nisar
    e.g. 1830_35.2003_35.26x16.nisar

    The subprocess is run with cwd=frameDir so that the geodat GeoJSONs are
    accessible via relative paths.  Simulation outputs (velSim, maskVel) are
    written to simDir inside frameDir.

    Parameters
    ----------
    frame : int
        Frame number.
    myArgs : dict
        SetupNISAR argument dictionary.  Required keys:
            orbit1, orbit2, NumberRangeLooks, NumberAzimuthLooks,
            outputDir, stderr, stdout.
        Optional keys: overWrite, overWritePhase, regionFile, verbose.
    simDir : str, optional
        Directory (relative to frameDir) where siminsar outputs (velSim,
        maskVel) are written.  Created if absent.  Default: 'simPhase'.
    '''
    orbit1 = myArgs['orbit1']
    orbit2 = myArgs['orbit2']
    nLooksR = myArgs['NumberRangeLooks']
    nLooksA = myArgs['NumberAzimuthLooks']
    frameDir = os.path.abspath(f'{myArgs["outputDir"]}/{orbit1}_{frame}')

    # Output stem and VRT path (relative to frameDir; estimateIonosphere runs
    # from there via cwd=frameDir)
    stem = f'{orbit1}_{frame}.{orbit2}_{frame}.{nLooksR}x{nLooksA}.nisar'
    outputVRT = f'{stem}.vrt'

    # Skip if output already exists and neither overwrite flag is set
    if myArgs.get('outputAll'):
        skip_path = os.path.join(frameDir, outputVRT)
    else:
        skip_path = os.path.join(frameDir, stem + '.ionosphereCorrection.vrt')
    if (os.path.exists(skip_path)
            and not myArgs.get('overWrite')
            and not myArgs.get('overWritePhase')):
        print(f'  estimateIonosphere: {os.path.basename(skip_path)} exists — skipping')
        return

    # Locate the RUNW HDF5 (H5/ subdir or directly in frameDir)
    RUNWFiles = (glob.glob(f'{frameDir}/H5/NISAR*RUNW*.h5') or
                 glob.glob(f'{frameDir}/NISAR*RUNW*.h5'))
    if not RUNWFiles:
        print(f'  estimateIonosphere: no RUNW HDF5 found under {frameDir}; '
              f'skipping frame {frame}')
        return
    RUNWFile = os.path.abspath(RUNWFiles[0])

    # Range offsets VRT — ROFFtoGrimp now writes to the frame dir; fall back
    # to the H5 subdir in case of older/different layouts
    offsetVRT = _find_file(frameDir,
                           ['range.offsets.vrt',
                            'H5/range.offsets.vrt'])
    if offsetVRT is None:
        print(f'  estimateIonosphere: range.offsets.vrt not found under '
              f'{frameDir}; skipping frame {frame}')
        return

    # Offset geometry VRT (simulated geometric offsets, band "RangeOffsets")
    # ROFFtoGrimp writes this to offsetSims/ inside the output dir
    geomVRT = _find_file(frameDir,
                         ['offsetSims/offsets.geom.vrt',
                          'H5/offsetSims/offsets.geom.vrt',
                          'offsets.geom.vrt'])

    command = [_PYTHON, '-m', 'nisargrimpworkflow.estimateIonosphere',
               RUNWFile,
               offsetVRT,
               outputVRT,
               '--frame', str(frame),
               '--simDir', simDir]

    if geomVRT is not None:
        command += ['--offset-geometry', geomVRT]
    if myArgs.get('regionFile'):
        command += ['--regionFile', myArgs['regionFile']]
    if myArgs.get('overWrite'):
        command += ['--overWrite']
    if myArgs.get('noInterp'):
        command += ['--noInterp']
    if myArgs.get('interpThresh') is not None:
        command += ['--interpThresh', str(myArgs['interpThresh'])]
    if myArgs.get('islandThresh') is not None:
        command += ['--islandThresh', str(myArgs['islandThresh'])]
    if myArgs.get('phaseThresh') is not None:
        command += ['--phaseThresh', str(myArgs['phaseThresh'])]
    if myArgs.get('outputAll'):
        command += ['--outputAll']
    if myArgs.get('verbose'):
        command += ['--verbose']

    print(f'  Running estimateIonosphere for frame {frame}...')
    run(command, cwd=frameDir,
        stderr=myArgs.get('stderr', DEVNULL),
        stdout=myArgs.get('stdout', DEVNULL))
