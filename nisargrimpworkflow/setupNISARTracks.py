#!/usr/bin/env python3
import os
import glob
import fnmatch
import re
import shutil
import subprocess
import argparse
import yaml
from concurrent.futures import ThreadPoolExecutor
import utilities as u

PROJECT_DIR = os.getcwd()


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            'Orchestrate NISAR velocity processing across all track directories\n'
            'in a project. Must be run from the project root directory — the one\n'
            'containing project.yaml and all track-* subdirectories. The project\n'
            'directory is derived from the current working directory.\n\n'
            'By default discovers every track-* directory under the project root,\n'
            'ensures tiepoints/ exists and creates velocityStats/ subdirectories\n'
            'derived from existing virtual-frame directories, then runs\n'
            'refreshties.py across all tracks for the requested\n'
            'years. Use --tracks to restrict processing to specific tracks.\n\n'
            'Optional flags enable first-time setup (--copyFiles), product\n'
            'checking (--check), and thumbnail generation (--runVelThumbs).\n\n'
            '--runVelstatsregions recomputes velocity-stats regions, and since that\n'
            'changes each vel_thumb_header\'s image extent, it also automatically\n'
            'chains a forced tiepointing pass (so makeframetie.py regenerates the\n'
            'per-range vel_thumb_plan/vel_thumbs output under the new extent) and a\n'
            'velocityStats.py rebuild -- all three run together, standalone (see\n'
            '--runVelstatsregions --help for the exact sequence).'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('--tracks', nargs='+', metavar='track-N',
                        help='Restrict processing to these track directories '
                             '(e.g. --tracks track-12 track-64). '
                             'Default: all track-* directories under the project root')
    parser.add_argument('--tiesOnly', action='store_true',
                        help='Run refreshties.py without computing velocity mosaics')
    parser.add_argument('--copyFiles', action='store_true',
                        help='Instantiate tie_plan_header and vel_thumb_plan from '
                             'templates/ into every track tiepoints/ directory '
                             '(substituting <TRACK> and <DEM>); skip if already present')
    parser.add_argument('--runVelThumbs', action='store_true',
                        help='Run "vel_thumbs vel_thumb_plan" in each tiepoints/ dir')
    parser.add_argument('--runVelstatsregions', action='store_true',
                        help='Run makevelstatsregions.py in each track dir, which rewrites each '
                             'vel_thumb_header_<range>\'s resolution/extent line (filling in the '
                             '"0 0 0 0" auto-size sentinel with the common domain covering all '
                             'data). Since the per-range vel_thumb_plan_<year>dash<range> file '
                             'vel_thumbs actually reads is only regenerated from that header by '
                             'makeframetie.py, this flag also automatically chains: (1) '
                             'makevelstatsregions.py, (2) a forced tiepointing pass (ties_only=False, '
                             'so makeframetie.py reruns and regenerates the per-range plan + calls '
                             'vel_thumbs itself under the new extent), (3) a "velocityStats.py '
                             '--doRA" rebuild from the freshly regenerated mosaics. Standalone, '
                             'early-return action -- not combined with --runVelThumbs/the default '
                             'tiepointing flow')
    parser.add_argument('--runVelStats', action='store_true',
                        help='Run "velocityStats.py --doRA" in each track\'s velocityStats/ '
                             'dir, using the project\'s region file and velMap to build the '
                             'reference basemap if missing')
    parser.add_argument('--updateAutoclean', action='store_true',
                        help='Force velStats+autocleanNISAR to (re)run for every track, even '
                             'ones whose frames already have an autoclean mask embedded in '
                             'range.offsets.vrt. Default: only run it for tracks missing masks')
    parser.add_argument('--nThreads', type=int, default=32, metavar='N',
                        help='Parallelism for the heavy steps: passed as '
                             '-nThreads to refreshties.py (tracks in parallel) '
                             'and as -threads to autocleanNISAR.py --tracks '
                             '(frame pool) [32]')
    parser.add_argument('--noUpdateVelStats', action='store_true',
                        help='Skip the "velocityStats.py --doRA" rebuild step in this invocation, '
                             'reusing the existing velStats reference. Use to avoid rebuilding '
                             'velStats twice when running --runVelstatsregions and then '
                             '--updateAutoclean back to back: pass it on whichever of the two '
                             'commands should not redo velStats (the other one leaves a fresh '
                             'reference for autocleanNISAR to consume)')
    parser.add_argument('--overWrite', action='store_true',
                        help='Pass --overWrite to refreshties.py to rerun existing products')
    parser.add_argument('--keepVz', action='store_true',
                        help='Pass --keepVz to refreshties.py to retain .vz and .vz.geodat files')
    parser.add_argument('--check', action='store_true',
                        help='Report which track-*/*_000? dirs have range.offsets.vrt '
                             'but lack velocity/mosaicOffsets.vx; print baseline sigmas')
    parser.add_argument('--checkFrames', action='store_true',
                        help='Check virtual-frame directory health: report orphaned dirs '
                             '(suffix does not match framePattern) and current-format dirs '
                             'missing upstream products (range offsets, phase, geodat, pairinfo)')
    parser.add_argument('--year', type=int, nargs='+', default=[2025, 2026], metavar='YYYY',
                        help='One or more years to pass to refreshties.py (default: 2025 2026)')
    parser.add_argument('--noPhase', action='store_true',
                        help='Disable phase+offsets mode (default: NISAR uses --phaseAndOffsets)')
    parser.add_argument('--quadFit', action='store_true',
                        help='Enable the -deltaBQ quadratic baseline correction estimate '
                             '(default: NISAR uses --noQuadFit)')
    parser.add_argument('--useSquint', dest='useSquint', action='store_true', default=None,
                        help='Apply the squint heading correction (mosaic3d -useSquint and '
                             'tiepoints -motion), overriding project.yaml applySquintCorrection '
                             'for this run')
    parser.add_argument('--noUseSquint', dest='useSquint', action='store_false',
                        help='Disable the squint heading correction for this run, overriding '
                             'project.yaml applySquintCorrection')
    return parser.parse_args()


def get_track_dirs(tracks=None):
    if tracks:
        dirs = [os.path.join(PROJECT_DIR, t) for t in tracks]
        missing = [d for d in dirs if not os.path.isdir(d)]
        if missing:
            u.myerror(f'Track directories not found: {missing}')
        return sorted(dirs, key=lambda p: int(re.search(r'track-(\d+)', p).group(1)))
    return sorted(glob.glob(os.path.join(PROJECT_DIR, 'track-*')),
                  key=lambda p: int(re.search(r'track-(\d+)', p).group(1)))


def _load_project_config():
    """Return (proj dict, dem string, velMap string, template_content dict,
    frame_pattern string) read from project.yaml."""
    proj_path = os.path.join(PROJECT_DIR, 'project.yaml')
    proj = {}
    if os.path.exists(proj_path):
        with open(proj_path) as f:
            proj = yaml.safe_load(f) or {}

    dem = ''
    velMap = ''
    region_path = proj.get('region') or proj.get('regionFile', '')
    for ext in ('', '.yaml'):
        rpath = (region_path + ext) if ext else region_path
        if rpath and os.path.exists(rpath):
            with open(rpath) as f:
                regionYaml = yaml.safe_load(f) or {}
            dem = regionYaml.get('dem', '')
            velMap = regionYaml.get('velMap', '')
            break

    templates_dir = os.path.join(PROJECT_DIR, 'templates')
    template_files = {
        'tie_plan_header': proj.get('tie_plan_header_template',
                                    os.path.join(templates_dir, 'tie_plan_header')),
        'vel_thumb_plan':  proj.get('vel_thumb_plan_template',
                                    os.path.join(templates_dir, 'vel_thumb_plan')),
        'vel_thumb_header': proj.get('vel_thumb_header_template',
                                     os.path.join(templates_dir, 'vel_thumb_header')),
    }
    contents = {}
    for name, path in template_files.items():
        if os.path.exists(path):
            with open(path) as f:
                contents[name] = f.read()
        else:
            contents[name] = None
    frame_pattern = proj.get('framePattern', '00??')
    return proj, dem, velMap, contents, frame_pattern


def _apply_substitutions(content, track_num, dem, tiepoint_file):
    return (content.replace('<TRACK>', track_num)
                    .replace('<DEM>', dem)
                    .replace('<TIEFILE>', tiepoint_file))


def _vel_stats_dirs_for_track(track_dir, frame_pattern):
    """Return sorted list of velocityStats dir names derived from existing virtual-frame dirs.

    Groups virtual-frame directories (matching *_{frame_pattern}) by the first
    variable digit (M), giving one velocityStats dir per M-group spanning
    {prefix}{M}0...0 – {prefix}{M}9...9.
    """
    prefix = frame_pattern.split('?')[0]          # e.g. '00' from '00??'
    variable_len = frame_pattern.count('?')        # e.g. 2
    tail_len = variable_len - 1                    # digits after M

    groups = set()
    for d in glob.glob(os.path.join(track_dir, f'*_{frame_pattern}')):
        vf = os.path.basename(d).rsplit('_', 1)[1]  # e.g. '0035'
        m_digit = vf[len(prefix)]                    # e.g. '3'
        start = prefix + m_digit + '0' * tail_len   # e.g. '0030'
        end   = prefix + m_digit + '9' * tail_len   # e.g. '0039'
        groups.add(f'{start}-{end}')
    return sorted(groups)


def setup_track_dirs(track_dirs, copy_files, overwrite=False):
    proj, dem, _, templates, frame_pattern = _load_project_config()
    tiepoint_file = proj.get('tiepointFile', '')

    n_created = 0
    n_skipped = 0

    for track_dir in track_dirs:
        track = os.path.basename(track_dir)
        track_num = track.split('-')[1]
        tpdir = os.path.join(track_dir, 'tiepoints')

        if not os.path.isdir(tpdir):
            os.makedirs(tpdir)
            print(f'Created {tpdir}')
            n_created += 1

        if copy_files:
            for tmpl_name, dest_name in (('tie_plan_header', 'tie_plan_header'),
                                         ('vel_thumb_plan', 'vel_thumb_plan')):
                dest = os.path.join(tpdir, dest_name)
                exists = os.path.exists(dest)
                if exists and not overwrite:
                    n_skipped += 1
                    continue
                if templates[tmpl_name] is None:
                    print(f'WARNING: template {tmpl_name} not found, skipping {dest}')
                    continue
                content = _apply_substitutions(templates[tmpl_name], track_num, dem, tiepoint_file)
                with open(dest, 'w') as f:
                    f.write(content)
                if exists:
                    print(f'Overwrote {dest}')
                else:
                    print(f'Created {dest}')
                    n_created += 1

        for dir_name in _vel_stats_dirs_for_track(track_dir, frame_pattern):
            vel_dir = os.path.join(track_dir, 'velocityStats', dir_name)
            if not os.path.isdir(vel_dir):
                os.makedirs(vel_dir)
                print(f'Created {vel_dir}')
                n_created += 1
            else:
                n_skipped += 1

        if templates['vel_thumb_header'] is not None:
            for vsd in sorted(glob.glob(os.path.join(track_dir, 'velocityStats', '*-*'))):
                frame_range = os.path.basename(vsd)
                x_dash_y = frame_range.replace('-', 'dash')
                dest = os.path.join(tpdir, f'vel_thumb_header_{x_dash_y}')
                if not os.path.exists(dest):
                    content = _apply_substitutions(templates['vel_thumb_header'], track_num, dem, tiepoint_file)
                    with open(dest, 'w') as f:
                        f.write(content)
                    print(f'Created {dest}')
                    n_created += 1
                else:
                    n_skipped += 1

    if n_skipped:
        print(f'Skipped {n_skipped} files/dirs (already exist), created {n_created} new')


def _extract_sigma(filepath):
    """Return sigma from az.est.const, rBaseline.deltabp (old text), or rBaseline.deltabp.yaml."""
    try:
        with open(filepath) as f:
            for line in f:
                if line.startswith('sigma:'):
                    return float(line.split(':')[1].split('#')[0].strip())
                if 'sigma*sqrt(X2/n)=' in line:
                    return float(line.split('=')[1].strip())
    except (OSError, ValueError):
        pass
    return None


def check_products():
    """Find track-*/*_000? dirs with range.offsets.vrt, report missing velocity/mosaicOffsets.vx,
    and print az/range sigmas sorted by range sigma."""
    pattern = os.path.join(PROJECT_DIR, 'track-*', '*_000?', 'range.offsets.vrt')
    vrt_files = sorted(glob.glob(pattern))

    missing = []
    sigma_rows = []
    for vrt in vrt_files:
        proc_dir = os.path.dirname(vrt)
        rel = os.path.relpath(proc_dir, PROJECT_DIR)
        vel_file = os.path.join(proc_dir, 'velocity', 'mosaicOffsets.vx')
        if not os.path.exists(vel_file):
            missing.append(proc_dir)

        az_sigma     = _extract_sigma(os.path.join(proc_dir, 'motion', 'az.est.const'))
        rng_sigma    = (_extract_sigma(os.path.join(proc_dir, 'motion', 'rBaseline.deltabp')) or
                        _extract_sigma(os.path.join(proc_dir, 'motion', 'rBaseline.deltabp.yaml')))
        rng_sigma_ni = _extract_sigma(os.path.join(proc_dir, 'motion', 'rBaseline.deltabp.noIonosphere'))
        if az_sigma is not None or rng_sigma is not None:
            sigma_rows.append((rel, rng_sigma, az_sigma, rng_sigma_ni))

    ok_count = len(vrt_files) - len(missing)
    print(f'\n{ok_count} directories have both range.offsets.vrt and velocity/mosaicOffsets.vx')
    if missing:
        print(f'\n{len(missing)} directories missing velocity/mosaicOffsets.vx:')
        for d in missing:
            print(f'  {d}')
    else:
        print('All directories are complete.')

    BOLD  = '\033[1m'
    RESET = '\033[0m'

    def fmt_sigma(val, bold):
        s = f'{val:.4f}' if val is not None else '   N/A'
        return f'{BOLD}{s}{RESET}' if bold else s

    if sigma_rows:
        sigma_rows.sort(key=lambda r: r[1] if r[1] is not None else float('-inf'), reverse=True)
        has_noion = any(r[3] is not None for r in sigma_rows)
        if has_noion:
            hdr = f'\n{"Directory":<35}  {"Rng sigma":>10}  {"Rng(noIon)":>10}  {"Az sigma":>10}'
            sep = '-' * 72
        else:
            hdr = f'\n{"Directory":<35}  {"Range sigma":>12}  {"Az sigma":>10}'
            sep = '-' * 62
        print(hdr)
        print(sep)
        n_rng_best = 0
        n_rng_ni_best = 0
        for rel, rng, az, rng_ni in sigma_rows:
            if has_noion:
                if rng is not None and rng_ni is not None:
                    rng_bold    = rng   <= rng_ni
                    rng_ni_bold = rng_ni < rng
                    if rng_bold:
                        n_rng_best += 1
                    else:
                        n_rng_ni_best += 1
                else:
                    rng_bold = rng_ni_bold = False
                rng_str    = fmt_sigma(rng,    rng_bold)
                rng_ni_str = fmt_sigma(rng_ni, rng_ni_bold)
                az_str     = fmt_sigma(az,     False)
                print(f'{rel:<35}  {rng_str:>10}  {rng_ni_str:>10}  {az_str:>10}')
            else:
                rng_str = f'{rng:.4f}' if rng is not None else '     N/A'
                az_str  = f'{az:.4f}'  if az  is not None else '     N/A'
                print(f'{rel:<35}  {rng_str:>12}  {az_str:>10}')

        def rss(vals):
            v = [x for x in vals if x is not None]
            return (sum(x*x for x in v) / len(v)) ** 0.5 if v else float('nan')

        rss_rng = rss([r[1] for r in sigma_rows])
        rss_az  = rss([r[2] for r in sigma_rows])

        if has_noion:
            rss_rng_ni = rss([r[3] for r in sigma_rows])
            best_sigmas = []
            for _, rng, _, rng_ni in sigma_rows:
                if rng is not None and rng_ni is not None:
                    best_sigmas.append(min(rng, rng_ni))
                elif rng is not None:
                    best_sigmas.append(rng)
                elif rng_ni is not None:
                    best_sigmas.append(rng_ni)
            rss_best = rss(best_sigmas)
            print(sep)
            print(f'{"Best count":<35}  {n_rng_best:>10}  {n_rng_ni_best:>10}')
            print(f'{"RSS sigma":<35}  {rss_rng:>10.4f}  {rss_rng_ni:>10.4f}  {rss_az:>10.4f}')
            print(f'{"RSS best":<35}  {rss_best:>10.4f}')
        else:
            print(sep)
            print(f'{"RSS sigma":<35}  {rss_rng:>12.4f}  {rss_az:>10.4f}')


def check_frames():
    """Check physical source frames and virtual frames for missing key files.

    Physical frames (suffix does not match framePattern) are checked for the
    files that virtual-frame VRTs reference.  Virtual frames (suffix matches
    framePattern) are checked for the assembled VRT products.
    """
    _, _, _, _, frame_pattern = _load_project_config()
    track_dirs = get_track_dirs()

    # Files that must exist in each physical source frame so the virtual-frame
    # VRTs can reference them.
    physical_checks = [
        ('range offsets',    ['range.offsets.vrt', 'range.offsets']),
        ('azimuth offsets',  ['azimuth.offsets.vrt', 'azimuth.offsets']),
        ('phase',            ['*.correctedUnwrappedPhase.vrt']),
        ('ionosphere',       ['*.ionosphereCorrection.globalFill.offset.vrt',
                              '*.ionosphereCorrection.offset.vrt']),
        ('geodat',           ['geodat*.geojson']),
        ('secondary geodat', ['geodat*.secondary.geojson']),
    ]

    # Files that must exist in each assembled virtual frame.
    virtual_checks = [
        ('range offsets', ['range.offsets.vrt', 'range.offsets']),
        ('phase',         ['phase.uw.vrt']),
        ('geodat',        ['geodat*.geojson']),
        ('pairinfo',      ['*.pairinfo']),
    ]

    phys_incomplete = []   # (rel_path, [missing_label, ...])
    virt_incomplete = []   # (rel_path, [missing_label, ...])
    empty_h5_stubs = []   # (rel_path, abs_path) — frame dirs with only an empty H5/
    frames_txt_stubs = [] # (rel_path, abs_path) — virtual frames with only frames.txt
    n_phys = 0
    n_virt = 0
    n_excluded = 0

    for track_dir in track_dirs:
        for entry in sorted(os.scandir(track_dir), key=lambda e: e.name):
            if not entry.is_dir():
                continue
            parts = entry.name.rsplit('_', 1)
            if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
                continue
            suffix = parts[1]
            frame_dir = entry.path
            rel = os.path.relpath(frame_dir, PROJECT_DIR)

            if os.path.exists(os.path.join(frame_dir, 'Exclude')) or \
                    os.path.exists(os.path.join(frame_dir, 'Exclude.pending')):
                n_excluded += 1
                continue

            # Detect stub: directory whose only visible content is an empty H5/ subdir
            visible = [f for f in os.listdir(frame_dir) if not f.startswith('.')]
            h5_dir = os.path.join(frame_dir, 'H5')
            if visible == ['H5'] and os.path.isdir(h5_dir) and not os.listdir(h5_dir):
                empty_h5_stubs.append((rel, frame_dir))
                continue

            if not fnmatch.fnmatch(suffix, frame_pattern):
                # Physical source frame
                n_phys += 1
                missing = [label
                           for label, pats in physical_checks
                           if not any(glob.glob(os.path.join(frame_dir, p)) for p in pats)]
                if missing:
                    phys_incomplete.append((rel, missing))
            else:
                # Virtual assembled frame — detect frames.txt-only stub
                if visible == ['frames.txt']:
                    frames_txt_stubs.append((rel, frame_dir))
                    continue
                n_virt += 1
                missing = [label
                           for label, pats in virtual_checks
                           if not any(glob.glob(os.path.join(frame_dir, p)) for p in pats)]
                if missing:
                    virt_incomplete.append((rel, missing))

    if empty_h5_stubs:
        print(f'\n{len(empty_h5_stubs)} frame director{"y" if len(empty_h5_stubs) == 1 else "ies"}'
              f' contain only an empty H5/ subdirectory (download stub, no data):')
        for rel, _ in empty_h5_stubs:
            print(f'  {rel}')
        answer = input(f'Remove all {len(empty_h5_stubs)}? [y/N]: ').strip().lower()
        if answer == 'y':
            for rel, abs_path in empty_h5_stubs:
                shutil.rmtree(abs_path)
                print(f'  Removed {rel}')
        else:
            print('  Skipped.')

    if frames_txt_stubs:
        print(f'\n{len(frames_txt_stubs)} virtual frame director{"y" if len(frames_txt_stubs) == 1 else "ies"}'
              f' contain only frames.txt (no assembled products):')
        for rel, _ in frames_txt_stubs:
            print(f'  {rel}')
        answer = input(f'Remove all {len(frames_txt_stubs)}? [y/N]: ').strip().lower()
        if answer == 'y':
            for rel, abs_path in frames_txt_stubs:
                shutil.rmtree(abs_path)
                print(f'  Removed {rel}')
        else:
            print('  Skipped.')

    excl_str = f' ({n_excluded} excluded)' if n_excluded else ''
    print(f'\nframePattern: {frame_pattern!r}')
    print(f'{n_phys} physical source frame(s), {n_virt} virtual frame(s){excl_str}')

    if phys_incomplete:
        print(f'\n{len(phys_incomplete)} INCOMPLETE physical frame(s)'
              f' (missing files that virtual-frame VRTs reference):')
        for r, miss in phys_incomplete:
            print(f'  {r}: missing {", ".join(miss)}')
    else:
        print(f'\nAll {n_phys} physical frames have required source files.')

    if virt_incomplete:
        print(f'\n{len(virt_incomplete)} INCOMPLETE virtual frame(s)'
              f' (missing assembled VRT products):')
        for r, miss in virt_incomplete:
            print(f'  {r}: missing {", ".join(miss)}')
    else:
        print(f'\nAll {n_virt} virtual frames have required assembled products.')


def run_vel_thumbs(track_dirs, over_write=False):
    overWrite = ' --overWrite' if over_write else ''

    def _run(track_dir):
        tpdir = os.path.join(track_dir, 'tiepoints')
        print(f'Running vel_thumbs in {tpdir}')
        subprocess.run(['csh', '-c', f'vel_thumbs vel_thumb_plan{overWrite}'], cwd=tpdir)

    with ThreadPoolExecutor(max_workers=10) as executor:
        executor.map(_run, track_dirs)


def run_velstats_regions(track_dirs):
    for track_dir in track_dirs:
        print(f'Running makevelstatsregions.py in {track_dir}')
        subprocess.run(['csh', '-c', 'makevelstatsregions.py'], cwd=track_dir)


def run_vel_stats(track_dirs, region_path, vel_map, n_threads=None):
    """Run 'velocityStats.py --doRA' in each track's velocityStats/ dir, passing
    the project's region file and velMap so the RA reference basemap is built
    (if missing) from the project's own velocity map rather than a generic
    region default. Tracks are independent (each velocityStats.py works only
    inside its own velocityStats/ dir), so they run in parallel up to n_threads
    (default: all at once) instead of serially."""
    regionArg = f'-regionFile {region_path} ' if region_path else ''
    velMapArg = f'-velmap {vel_map} ' if vel_map else ''
    cmd = f'velocityStats.py --doRA {regionArg}{velMapArg}'

    jobs = []
    for track_dir in track_dirs:
        vsd = os.path.join(track_dir, 'velocityStats')
        if not os.path.isdir(vsd):
            print(f'Skipping {track_dir}: no velocityStats/ dir')
            continue
        jobs.append(vsd)
    if not jobs:
        return

    def _run(vsd):
        print(f'Running: {cmd}  (in {vsd})')
        subprocess.run(['csh', '-c', cmd], cwd=vsd)

    maxWorkers = n_threads if n_threads else len(jobs)
    with ThreadPoolExecutor(max_workers=maxWorkers) as executor:
        list(executor.map(_run, jobs))


def _track_has_autoclean_mask(track_dir, frame_pattern):
    """Return True if every merged frame dir in track_dir already has an
    autoclean mask embedded in its range.offsets.vrt (<MaskBand> block written
    by autocleanNISAR.py:embedMaskBand()). False if any frame is missing one
    (including tracks with no frame dirs at all, or no range.offsets.vrt
    yet) -- those need (re)running."""
    frame_dirs = glob.glob(os.path.join(track_dir, f'*_{frame_pattern}'))
    if not frame_dirs:
        return False
    for frame_dir in frame_dirs:
        vrt = os.path.join(frame_dir, 'range.offsets.vrt')
        if not os.path.exists(vrt):
            return False
        with open(vrt) as f:
            if '<MaskBand' not in f.read():
                return False
    return True


def run_autoclean(track_dirs, frame_pattern, n_threads=None, new=True):
    """Run 'autocleanNISAR.py --tracks' once from the project root over all
    track dirs -- flags outlier range/azimuth offset pixels against the
    velocityStats reference (just refreshed by run_vel_stats) and embeds the
    resulting masks into each frame's range.offsets.vrt/azimuth.offsets.vrt,
    so the tiepointing step that follows (rparams/azparams, both
    -noMask-aware) fits baselines against cleaned data. --tracks pools every
    frame across all tracks into one multithreaded run (previously this
    looped 'autocleanNISAR.py --allFrames' per track serially, so each
    track's straggler frames idled the pool before the next track started).
    When new is True (the default), '--new' is passed so only frames without
    an embedded mask are (re)cleaned -- already-cleaned frames are skipped;
    --updateAutoclean forces new=False to reclean every frame."""
    if not track_dirs:
        return
    # frame_pattern (e.g. "00??") contains csh glob metacharacters --
    # quote it or csh expands it against cwd's filenames before
    # autocleanNISAR.py ever sees it, aborting with "No match." when
    # nothing happens to match (same reason run_refresh_ties() quotes
    # -toRun="...").
    tracks = ' '.join(os.path.basename(d) for d in track_dirs)
    threadsArg = f' -threads {n_threads}' if n_threads else ''
    newArg = ' --new' if new else ''
    cmd = (f'autocleanNISAR.py --tracks {tracks}'
           f' -framePattern "{frame_pattern}"{threadsArg}{newArg}')
    print(f'Running: {cmd}  (in {PROJECT_DIR})')
    subprocess.run(['csh', '-c', cmd], cwd=PROJECT_DIR)


def run_refresh_ties(track_dirs, ties_only, over_write, keep_vz, year,
                      phase_and_offsets=True, no_quad_fit=True, use_yaml=False,
                      use_squint=False, n_threads=None):
    tracks_str = str([os.path.basename(d) for d in track_dirs]).replace(' ', '')
    tiesOnly = '-tiesOnly ' if ties_only else ''
    overWrite = '--overWrite ' if over_write else ''
    keepVz = '--keepVz ' if keep_vz else ''
    phaseAndOffsets = '--phaseAndOffsets ' if phase_and_offsets else ''
    noQuadFit = '--noQuadFit ' if no_quad_fit else ''
    useYaml = '--yaml ' if use_yaml else ''
    useSquint = '--useSquint ' if use_squint else ''
    nThreadsFlag = f'-nThreads {n_threads} ' if n_threads else ''
    years_str = ' '.join(str(y) for y in year)
    cmd = f'refreshties.py {tiesOnly}{overWrite}{keepVz}{phaseAndOffsets}{noQuadFit}{useYaml}{useSquint}{nThreadsFlag}-toRun="{tracks_str}" {years_str} -noPrompt'
    print(f'Running: {cmd}')
    subprocess.run(['csh', '-c', cmd], cwd=PROJECT_DIR)


def main():
    args = parse_args()

    if args.check:
        check_products()
        return

    if args.checkFrames:
        check_frames()
        return

    track_dirs = get_track_dirs(args.tracks)
    proj, _, vel_map, _, frame_pattern = _load_project_config()
    region_path = proj.get('region') or proj.get('regionFile', '')

    # Drop tracks that have no virtual-frame data yet.
    have_data = [d for d in track_dirs
                 if glob.glob(os.path.join(d, f'*_{frame_pattern}'))]
    empty = sorted(set(os.path.basename(d) for d in track_dirs)
                   - set(os.path.basename(d) for d in have_data),
                   key=lambda t: int(t.split('-')[1]))
    if empty:
        print(f'Skipping {len(empty)} track(s) with no virtual-frame data: {empty}')
    track_dirs = have_data

    setup_track_dirs(track_dirs, args.copyFiles, overwrite=args.overWrite)

    if args.copyFiles:
        return

    projectUseSquint = bool(proj.get('applySquintCorrection', False))
    useSquint = projectUseSquint if args.useSquint is None else args.useSquint

    if args.runVelThumbs:
        run_vel_thumbs(track_dirs, over_write=args.overWrite)

    if args.runVelstatsregions:
        run_velstats_regions(track_dirs)
        # makevelstatsregions.py rewrites each vel_thumb_header_<range>'s
        # resolution/extent line to the common domain covering all data
        # (filling in the "0 0 0 0" auto-size sentinel). The per-range plan
        # file vel_thumbs actually reads (vel_thumb_plan_<year>dash<range>)
        # is only regenerated from that header by makeframetie.py -- which
        # also calls vel_thumbs on it -- so force a full tiepointing pass
        # (ties_only=False, so makeframetie.py actually runs) rather than
        # calling vel_thumbs on the generic vel_thumb_plan directly (a
        # separate, whole-track auto-extent file that never carries
        # per-range resolutions).
        run_refresh_ties(track_dirs, False, True, args.keepVz, args.year,
                         phase_and_offsets=not args.noPhase,
                         no_quad_fit=not args.quadFit,
                         use_yaml=True,
                         use_squint=useSquint,
                         n_threads=args.nThreads)
        # velocityStats aggregates from those same source mosaics, so its
        # own reference (velocity.vr/.va/.er/.ea/.navg) is now stale too --
        # rebuild it from the freshly regenerated mosaics (unless the caller
        # will do it in a following --updateAutoclean run, via --noUpdateVelStats).
        if not args.noUpdateVelStats:
            run_vel_stats(track_dirs, region_path, vel_map,
                          n_threads=args.nThreads)
        return

    if args.runVelStats:
        run_vel_stats(track_dirs, region_path, vel_map,
                      n_threads=args.nThreads)
        return

    # velStats+autoclean only for tracks that don't already have an autoclean
    # mask embedded (or all of them, with --updateAutoclean) -- tiepointing
    # below always runs for every track regardless, since it's what actually
    # benefits from the (possibly just-refreshed) masks.
    needs_autoclean = [d for d in track_dirs
                       if args.updateAutoclean or not _track_has_autoclean_mask(d, frame_pattern)]
    if needs_autoclean:
        # --noUpdateVelStats reuses velStats built by a prior --runVelstatsregions
        # run, so autoclean consumes it directly without a redundant rebuild.
        if not args.noUpdateVelStats:
            run_vel_stats(needs_autoclean, region_path, vel_map,
                          n_threads=args.nThreads)
        run_autoclean(needs_autoclean, frame_pattern, n_threads=args.nThreads,
                      new=not args.updateAutoclean)

    run_refresh_ties(track_dirs, args.tiesOnly, args.overWrite, args.keepVz, args.year,
                     phase_and_offsets=not args.noPhase,
                     no_quad_fit=not args.quadFit,
                     use_yaml=True,
                     use_squint=useSquint,
                     n_threads=args.nThreads)


if __name__ == '__main__':
    main()
