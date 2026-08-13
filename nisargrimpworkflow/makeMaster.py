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
SetupNISAR.globalFillIonosphere). It also verifies that each range baseline's
recorded offsetCorrectionFile agrees with that same VRT metadata -- the exact
consistency mosaic3d/checkForIonosphereCorrection() enforces at runtime, so a
line that would abort the mosaic (e.g. a baseline still pointing at an old
globalFill name after the ionosphere product was regenerated) is caught here
instead. Entries with any missing file or ionosphere mismatch are skipped from
the output.  A summary of how many were skipped and why is printed at the end;
--verbose adds a sorted list of every missing file path / mismatch per category.
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict
from functools import lru_cache

import utilities as u
from osgeo import gdal


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Assemble master inputFile from per-track '
        'inputFile-N-all20?? files.',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('--projectDir', type=str, default='.',
                        help='Root directory for anchoring --outputPath when '
                        'it is a relative path [.]. Ignored for track '
                        'discovery when --inputPaths is given.')
    parser.add_argument('--inputPaths', nargs='+', metavar='PATH',
                        help='One or more directories to search for track-N '
                        'subdirectories (overrides --projectDir for input). '
                        'Example: . ../anotherDir')
    parser.add_argument('--outputPath', default='Release/masterInput',
                        metavar='PATH',
                        help='Directory to write inputFile (default: '
                        'Release/masterInput, relative to --projectDir '
                        'unless absolute)')
    parser.add_argument('--verbose', action='store_true',
                        help='List every missing file path in the skip summary')
    return parser.parse_args()


def get_track_entries(input_paths):
    """Return sorted list of (track_num, base_dir) pairs found across all input paths."""
    entries = []
    for base_dir in input_paths:
        dirs = glob.glob(os.path.join(base_dir, 'track-*'))
        for d in dirs:
            m = re.search(r'track-(\d+)$', d)
            if m:
                entries.append((int(m.group(1)), base_dir))
    return sorted(entries, key=lambda x: x[0])


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


@lru_cache(maxsize=None)
def find_ion_correction(range_offsets_path):
    """If range_offsets_path's .vrt has an ionosphereRangeOffsetCorrection
    metadata tag (written by SetupNISAR.globalFillIonosphere), return the
    resolved path to that ion-correction file -- relative to
    range_offsets_path's own directory, matching how rparams/
    checkForIonosphereCorrection() resolves it. Returns None if there's no
    .vrt or no such tag (nothing to verify, not an error). Cached: the same
    range.offsets is inspected for both file-existence and baseline
    consistency."""
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


def baseline_offset_correction(baseline_path):
    """Return the basename of the range baseline's offsetCorrectionFile entry,
    or None if it records no correction. 'nil' and a missing/blank entry both
    map to None, matching getRParams()/readOffsets.c which treats 'nil' as
    unset (checkForIonosphereCorrection then returns early, no comparison)."""
    try:
        with open(baseline_path) as f:
            for line in f:
                if line.strip().startswith('offsetCorrectionFile:'):
                    val = line.split(':', 1)[1].strip()
                    if not val or val == 'nil':
                        return None
                    return os.path.basename(val)
    except OSError:
        return None
    return None


def check_line_ion_consistency(line):
    """Return a list of mismatch messages when a data line's range baseline
    records an offsetCorrectionFile whose basename disagrees with the
    range.offsets VRT's ionosphereRangeOffsetCorrection metadata -- the same
    inconsistency mosaic3d/checkForIonosphereCorrection() aborts on. Empty
    when consistent or not applicable (no range.offsets, no range baseline,
    or the baseline records no correction)."""
    range_offsets_tok = None
    baseline_tok = None
    for tok in line.split():
        if not is_file_token(tok):
            continue
        base = os.path.basename(tok)
        if base == 'range.offsets':
            range_offsets_tok = tok
        elif base.startswith('rBaseline') and base.endswith('.yaml'):
            baseline_tok = tok
    if range_offsets_tok is None or baseline_tok is None:
        return []
    blName = baseline_offset_correction(baseline_tok)
    if blName is None:
        return []
    ionPath = find_ion_correction(range_offsets_tok)
    vrtName = os.path.basename(ionPath) if ionPath is not None else None
    if blName != vrtName:
        return [f'{baseline_tok}: offsetCorrectionFile {blName} != VRT '
                f'ionosphereRangeOffsetCorrection {vrtName}']
    return []


def check_line_files(line):
    """Return list of missing file paths for a single data line."""
    missing = []
    for tok in line.split():
        if not is_file_token(tok):
            continue
        if not token_exists(tok):
            missing.append(tok)
        elif os.path.basename(tok) == 'range.offsets':
            ionPath = find_ion_correction(tok)
            if ionPath is not None and not os.path.exists(ionPath):
                missing.append(ionPath)
    return missing


def file_type_label(path):
    """Categorise a missing file path into a short human-readable type label."""
    base = os.path.basename(path)
    if 'rBaseline' in base:
        return 'rBaseline'
    if 'range.offsets' in base:
        return 'range.offsets'
    if 'azimuth.offsets' in base:
        return 'azimuth.offsets'
    if 'ionosphere' in base.lower():
        return 'ionosphere correction'
    if base.endswith('.geojson') or 'geodat' in base:
        return 'geodat'
    return base


def orbit_frame_from_line(line):
    """Extract the orbit_frame token (e.g. '4040_0000') from a data line."""
    for tok in line.split():
        if not is_file_token(tok):
            continue
        # Walk up from the file path looking for a ????_???? component
        parts = tok.replace('\\', '/').split('/')
        for part in parts:
            if re.fullmatch(r'\d+_\d+', part):
                return part
    return None


def assemble_master(track_entries, verbose=False):
    """Build output lines for all tracks, skipping entries with missing files.

    track_entries is a list of (track_num, base_dir) pairs as returned by
    get_track_entries().

    Returns (lines, n_total, n_skipped, skip_counts, skip_files,
             tracks_no_inputfile, frames_in_output, n_mismatch, mismatch_msgs).
      n_total:            total data lines found across all inputFiles
      n_skipped:          lines dropped (missing files and/or ion mismatch)
      skip_counts:        Counter {type_label: lines_missing_that_type}
      skip_files:         {type_label: [path, ...]}  (verbose only)
      tracks_no_inputfile: list of track numbers with no inputFile found
      frames_in_output:   set of orbit_frame strings that made it into output
      n_mismatch:         lines dropped for a baseline/VRT ionosphere mismatch
      mismatch_msgs:      list of mismatch description strings
    """
    output = []
    n_total = 0
    n_skipped = 0
    skip_counts = Counter()
    skip_files = defaultdict(list)
    tracks_no_inputfile = []
    frames_in_output = set()
    n_mismatch = 0
    mismatch_msgs = []

    for track_num, base_dir in track_entries:
        input_files = find_input_files(base_dir, track_num)
        if not input_files:
            tracks_no_inputfile.append(track_num)
            continue
        output.append(f'; track-{track_num}')
        for i, filepath in enumerate(input_files):
            if i > 0:
                output.append(';')
            for line in extract_data_lines(filepath):
                n_total += 1
                missing = check_line_files(line)
                mismatches = check_line_ion_consistency(line)
                if missing or mismatches:
                    n_skipped += 1
                    seen_types = set()
                    for m in missing:
                        label = file_type_label(m)
                        if label not in seen_types:
                            skip_counts[label] += 1
                            seen_types.add(label)
                        if verbose:
                            skip_files[label].append(m)
                    if mismatches:
                        n_mismatch += 1
                        mismatch_msgs.extend(mismatches)
                else:
                    frame = orbit_frame_from_line(line)
                    if frame:
                        frames_in_output.add(frame)
                    output.append(line)
    return (output, n_total, n_skipped, skip_counts, skip_files,
            tracks_no_inputfile, frames_in_output, n_mismatch, mismatch_msgs)


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
    input_paths = args.inputPaths if args.inputPaths else [args.projectDir]
    # Anchor a relative --outputPath to --projectDir (default CWD),
    # independently of --inputPaths.
    output_path = (args.outputPath if os.path.isabs(args.outputPath)
                   else os.path.join(args.projectDir, args.outputPath))

    track_entries = get_track_entries(input_paths)
    track_numbers = sorted(set(n for n, _ in track_entries))
    print(f'Found tracks: {track_numbers}')

    motion_dirs = []
    for p in input_paths:
        motion_dirs.extend(glob.glob(os.path.join(p, 'track-*', '*_0???', 'motion')))
    print(f'Found {len(motion_dirs)} track-*/*_0???/motion directories')

    (lines, n_total, n_skipped, skip_counts, skip_files,
     tracks_no_inputfile, frames_in_output, n_mismatch,
     mismatch_msgs) = assemble_master(track_entries, verbose=args.verbose)

    n_motion = len(motion_dirs)
    n_tracks_no_data = len(tracks_no_inputfile)
    n_frames_out = len(frames_in_output)

    print(f'Tracks with no inputFile (tiepoints not yet run): '
          f'{n_tracks_no_data} of {len(track_entries)}')
    if args.verbose and tracks_no_inputfile:
        print(f'  tracks: {tracks_no_inputfile}')
    print(f'Data lines found in inputFiles: {n_total}  '
          f'kept: {n_total - n_skipped}  '
          f'skipped (missing files): {n_skipped}')
    print(f'Unique frames in output: {n_frames_out} of {n_motion} motion dirs')

    if n_skipped:
        print(f'\033[1;43mSkip breakdown:\033[0m')
        for label, count in sorted(skip_counts.items(), key=lambda x: -x[1]):
            print(f'  {count:4d} missing {label}')
            if args.verbose:
                for f in sorted(set(skip_files[label])):
                    print(f'    {f}')

    if n_mismatch:
        print(f'\033[1;41mBaseline/VRT ionosphere mismatches '
              f'(re-fit baseline): {n_mismatch}\033[0m')
        for msg in sorted(set(mismatch_msgs)):
            print(f'    {msg}')

    write_output(lines, output_path)


if __name__ == '__main__':
    main()
