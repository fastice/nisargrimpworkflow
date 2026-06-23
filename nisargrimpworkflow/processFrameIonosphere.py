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
import re
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
            orbit1, orbit2, outputDir, stderr, stdout.
        Optional keys: overWrite, overWritePhase, regionFile, verbose,
            NumberRangeLooks, NumberAzimuthLooks (fallback only -- see below).
    simDir : str, optional
        Directory (relative to frameDir) where siminsar outputs (velSim,
        maskVel) are written.  Created if absent.  Default: 'simPhase'.
    '''
    orbit1 = myArgs['orbit1']
    orbit2 = myArgs['orbit2']
    frameDir = os.path.abspath(f'{myArgs["outputDir"]}/{orbit1}_{frame}')

    # nLooksR/nLooksA must come from THIS frame's own geodat, not myArgs --
    # myArgs['NumberRangeLooks']/['NumberAzimuthLooks'] are set once globally in
    # getSecondaryOrbit() from whichever frame's H5 was found first, and a frame
    # can have a genuinely different bandwidth (hence different looks) than its
    # neighbors (e.g. track-58 frame 44 at 40 MHz/13 looks vs 77 MHz/26 looks).
    # Using the wrong global value here bakes a mismatched looks suffix into the
    # correctedUnwrappedPhase/ionosphereCorrection filenames, which then breaks
    # downstream geodat lookups (e.g. insarworkflow.tieScript.process_phase_yaml).
    ownGeodats = [f for f in glob.glob(f'{frameDir}/geodat*x*.geojson')
                  if '.secondary.' not in os.path.basename(f)]
    m = (re.match(r'geodat(\d+)x(\d+)\.geojson$', os.path.basename(ownGeodats[0]))
         if ownGeodats else None)
    if m:
        nLooksR, nLooksA = int(m.group(1)), int(m.group(2))
    else:
        nLooksR = myArgs['NumberRangeLooks']
        nLooksA = myArgs['NumberAzimuthLooks']

    # Output stem and VRT path (relative to frameDir; estimateIonosphere runs
    # from there via cwd=frameDir)
    stem = f'{orbit1}_{frame}.{orbit2}_{frame}.{nLooksR}x{nLooksA}.nisar'
    outputVRT = f'{stem}.vrt'

    # Skip if output already exists and neither overwrite flag is set -- unless variable
    # smoothing is newly requested and wasn't produced by an earlier (pre-smoothing) run
    # (mirrors the same check inside run_vel_sim() for velSim.smr).
    if myArgs.get('outputAll'):
        skip_path = os.path.join(frameDir, outputVRT)
    else:
        skip_path = os.path.join(frameDir, stem + '.ionosphereCorrection.vrt')
    smrMissing = (myArgs.get('minTol') is not None
                 and not myArgs.get('noVariableSmoothing')
                 and not os.path.exists(os.path.join(frameDir, simDir, 'velSim.smr.vrt')))
    if (os.path.exists(skip_path)
            and not myArgs.get('overWrite')
            and not myArgs.get('overWritePhase')
            and not smrMissing):
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

    # Ice/rock/water mask for --sepIceRock (0=water, 1=rock, 2=ice)
    iceRockMaskVRT = _find_file(frameDir,
                                ['offsetSims/offsets.geom.mask.vrt',
                                 'H5/offsetSims/offsets.geom.mask.vrt'])

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
    if myArgs.get('verticalCorrection'):
        command += ['--verticalCorrection', myArgs['verticalCorrection']]
    if myArgs.get('overWrite'):
        command += ['--overWrite']
    if myArgs.get('minTol') is not None:
        command += ['--minTol', str(myArgs['minTol']),
                    '--percentSpeed', str(myArgs['percentSpeed']),
                    '--maxTol', str(myArgs['maxTol']),
                    '--maxSmoothRadius', str(myArgs['maxSmoothRadius']),
                    '--smoothNIter', str(myArgs['smoothNIter'])]
        if myArgs.get('noVariableSmoothing'):
            command += ['--noVariableSmoothing']
    if myArgs.get('noInterp'):
        command += ['--noInterp']
    if myArgs.get('interpThresh') is not None:
        command += ['--interpThresh', str(myArgs['interpThresh'])]
    if myArgs.get('islandThresh') is not None:
        command += ['--islandThresh', str(myArgs['islandThresh'])]
    if myArgs.get('phaseThresh') is not None:
        command += ['--phaseThresh', str(myArgs['phaseThresh'])]
    if myArgs.get('sigmaAz') is not None:
        command += ['--sigma-az', str(myArgs['sigmaAz'])]
    if myArgs.get('sigmaRg') is not None:
        command += ['--sigma-rg', str(myArgs['sigmaRg'])]
    if myArgs.get('noPhaseThreshPass'):
        command += ['--noPhaseThreshPass']
    if myArgs.get('sepIceRock'):
        command += ['--sepIceRock']
        if myArgs.get('iceRockMask'):
            command += ['--iceRockMask', myArgs['iceRockMask']]
        elif iceRockMaskVRT is not None:
            command += ['--iceRockMask', iceRockMaskVRT]
        else:
            print(f'  Warning: --sepIceRock set but no offsets.geom.mask.vrt found '
                  f'under {frameDir}; rock pre-seeding will be skipped.')
    if myArgs.get('outputAll'):
        command += ['--outputAll']
    if myArgs.get('debugIono'):
        command += ['--debugIono']
    if myArgs.get('verbose'):
        command += ['--verbose']

    print(f'  Running estimateIonosphere for frame {frame}...')
    run(command, cwd=frameDir,
        stderr=myArgs.get('stderr', DEVNULL),
        stdout=myArgs.get('stdout', DEVNULL))
