#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ian
"""

import argparse
import glob
import os

requiredFiles = [
    ('offsets', ['range.offsets.vrt', 'azimuth.offsets.vrt']),
    ('phase', ['phase.uw.vrt']),
    ('ioncorrection', ['*.ionosphereCorrection.globalFill.offset.vrt',
                       '*.ionosphereCorrection.offset.vrt']),
]


def checkFrameInputsArgs():
    '''Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mScan track-N/*_00NN virtual frames for '
        'missing offsets, phase, or ionosphere-correction products\033[0m\n\n',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('--projectDir', type=str, default='.',
                        help='Root directory containing track-N '
                        'subdirectories [.] (current directory)')
    args = parser.parse_args()
    return args.projectDir


def findAny(frameDir, patterns):
    '''True if any of patterns (literal names or globs) exists in frameDir.'''
    for pattern in patterns:
        if glob.glob(os.path.join(frameDir, pattern)):
            return True
    return False


def checkFrame(frameDir):
    '''Return a list of (label, patterns) missing from this virtual frame.'''
    missing = []
    for label, patterns in requiredFiles:
        if not findAny(frameDir, patterns):
            missing.append(label)
    return missing


def main():
    projectDir = checkFrameInputsArgs()
    frameDirs = sorted(glob.glob(os.path.join(projectDir, 'track-*', '*_00??')))
    missingCounts = {label: 0 for label, _ in requiredFiles}
    nMissing = 0
    nExcluded = 0
    for frameDir in frameDirs:
        if os.path.exists(os.path.join(frameDir, 'Exclude')):
            nExcluded += 1
            continue
        missing = checkFrame(frameDir)
        if missing:
            nMissing += 1
            for label in missing:
                missingCounts[label] += 1
            print(f'{frameDir}: missing {", ".join(missing)}')
    print(f'\n{nMissing} / {len(frameDirs)} virtual frame(s) missing one or '
          f'more required inputs ({nExcluded} excluded via an Exclude marker, '
          f'not counted above)')
    for label, count in missingCounts.items():
        print(f'  {label}: {count} frame(s) missing')


if __name__ == '__main__':
    main()
