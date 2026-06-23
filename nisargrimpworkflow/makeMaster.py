#!/usr/bin/env python3
"""
Assembles a master inputFile from per-track per-year inputFile-N-all20?? files.

For each track-N, finds track-N/tiepoints/track-N/inputFile-N-all20?? and
extracts the long data lines (skipping ; comment lines and pure-number lines).
The result is written to outputPath/inputFile with a '; track-N' header before
each track's block and a ';' separator between consecutive per-year files.

Before writing, verifies that every file path referenced in those data lines
actually exists -- including any ionosphere-correction file embedded in a
range.offsets VRT's ionosphereRangeOffsetCorrection metadata (see
SetupNISAR.globalFillIonosphere) -- and reports (but does not block on) any
that are missing.
"""

import argparse
import glob
import os
import re
import sys

import utilities as u
from osgeo import gdal


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Assemble master inputFile from per-track '
        'inputFile-N-all20?? files.',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('--projectDir', type=str, default='.',
                        help='Root directory containing track-N '
                        'subdirectories [.] (current directory)')
    parser.add_argument('--outputPath', default='Release/masterInput',
                        metavar='PATH',
                        help='Directory to write inputFile (default: '
                        'Release/masterInput, relative to --projectDir '
                        'unless absolute)')
    return parser.parse_args()


def get_track_numbers(project_dir):
    """Return sorted list of track numbers found in the project dir."""
    dirs = glob.glob(os.path.join(project_dir, 'track-*'))
    numbers = []
    for d in dirs:
        m = re.search(r'track-(\d+)$', d)
        if m:
            numbers.append(int(m.group(1)))
    return sorted(numbers)


def find_input_files(project_dir, track_num):
    """Return sorted list of inputFile-N-all20?? paths for a given track number."""
    pattern = os.path.join(project_dir,
                           f'track-{track_num}', 'tiepoints',
                           f'track-{track_num}',
                           f'inputFile-{track_num}-all20??')
    return sorted(glob.glob(pattern))


def is_data_line(line):
    """Return True for long data lines — skip ; comments and pure-number lines."""
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith(';'):
        return False
    # Pure-number lines: integers or whitespace-separated floats/ints
    try:
        [float(tok) for tok in stripped.split()]
        return False
    except ValueError:
        pass
    return True


def extract_data_lines(filepath):
    """Return the data lines from an inputFile."""
    with open(filepath) as f:
        return [line.rstrip() for line in f if is_data_line(line)]


def is_file_token(tok):
    """True for tokens in a data line that are file paths, not numbers or 'None'."""
    if tok.lower() == 'none':
        return False
    try:
        float(tok)
        return False
    except ValueError:
        return True


def token_exists(tok):
    """True if tok exists as given, or as tok + '.vrt' -- inputFile templates
    write bare names like 'range.offsets'/'azimuth.offsets' unconditionally,
    but projects that have moved to VRT-only storage (no bare flat-binary
    file on disk) rely on rparams/azparams (getROffsets.c) auto-detecting
    and falling back to the .vrt sidecar when given the bare name."""
    return os.path.exists(tok) or os.path.exists(tok + '.vrt')


def find_ion_correction(range_offsets_path):
    """If range_offsets_path's .vrt has an ionosphereRangeOffsetCorrection
    metadata tag (written by SetupNISAR.globalFillIonosphere), return the
    resolved path to that ion-correction file -- relative to
    range_offsets_path's own directory, matching how rparams/
    checkForIonosphereCorrection() resolves it. Returns None if there's no
    .vrt or no such tag (nothing to verify, not an error)."""
    vrtPath = range_offsets_path + '.vrt'
    if not os.path.exists(vrtPath):
        return None
    ds = gdal.Open(vrtPath)
    if ds is None:
        return None
    ionName = ds.GetMetadataItem('ionosphereRangeOffsetCorrection')
    ds = None
    if not ionName:
        return None
    return os.path.join(os.path.dirname(range_offsets_path), ionName)


def verify_referenced_files(track_numbers, project_dir):
    """Check that every file referenced in each track's inputFile-N-all20??
    exists, including any ion-correction file embedded in range.offsets'
    VRT metadata. Prints one warning per missing file. Returns the count of
    missing files found."""
    n_missing = 0
    for track_num in track_numbers:
        for filepath in find_input_files(project_dir, track_num):
            for line in extract_data_lines(filepath):
                for tok in line.split():
                    if not is_file_token(tok):
                        continue
                    if not token_exists(tok):
                        u.mywarning(f'track-{track_num} '
                                   f'{os.path.basename(filepath)}: missing {tok}')
                        n_missing += 1
                        continue
                    if os.path.basename(tok) == 'range.offsets':
                        ionPath = find_ion_correction(tok)
                        if ionPath is not None and not os.path.exists(ionPath):
                            u.mywarning(f'track-{track_num} '
                                       f'{os.path.basename(filepath)}: missing '
                                       f'ion correction {ionPath} (referenced by '
                                       f'{tok})')
                            n_missing += 1
    return n_missing


def assemble_master(track_numbers, project_dir):
    """Build output lines for all tracks."""
    output = []
    for track_num in track_numbers:
        input_files = find_input_files(project_dir, track_num)
        if not input_files:
            continue
        output.append(f'; track-{track_num}')
        for i, filepath in enumerate(input_files):
            if i > 0:
                output.append(';')
            lines = extract_data_lines(filepath)
            output.extend(lines)
    return output


def write_output(lines, output_path):
    """Write assembled lines to outputPath/inputFile.

    Exits with an error if the parent directory (e.g. Release) does not exist.
    Creates the output directory itself if it is missing.
    """
    parent = os.path.dirname(output_path)
    if parent and not os.path.isdir(parent):
        print(f'\033[1;31mError: parent directory does not exist: {parent}\033[0m', file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
        print(f'Created {output_path}')
    # Count only true data lines (not the ; track-N headers or ; separators)
    n_data = sum(1 for l in lines if not l.startswith(';'))
    out_file = os.path.join(output_path, 'inputFile')
    with open(out_file, 'w') as f:
        f.write('0 0 0 0 .2 .2\n')
        f.write(';\n')
        f.write(f'{n_data}\n')
        f.write(';\n')
        for i, line in enumerate(lines):
            f.write(line + '\n')
            if i < len(lines) - 1:
                f.write(';\n')
    print(f'Wrote {n_data} data lines to {out_file}')


def main():
    args = parseArgs()
    project_dir = args.projectDir
    output_path = (args.outputPath if os.path.isabs(args.outputPath)
                   else os.path.join(project_dir, args.outputPath))

    track_numbers = get_track_numbers(project_dir)
    print(f'Found tracks: {track_numbers}')

    n_missing = verify_referenced_files(track_numbers, project_dir)
    if n_missing:
        print(f'\033[1;43m{n_missing} referenced file(s) missing — see warnings '
              f'above\033[0m')
    else:
        print('All referenced files verified present.')

    lines = assemble_master(track_numbers, project_dir)
    write_output(lines, output_path)


if __name__ == '__main__':
    main()
