#!/usr/bin/env python3
import os
import glob
import re
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
            'checking (--check), thumbnail generation (--runVelThumbs), and\n'
            'velocity-stats-region recomputation (--runVelstatsregions).'
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
                        help='Run makevelstatsregions.py in each track dir')
    parser.add_argument('--overWrite', action='store_true',
                        help='Pass --overWrite to refreshties.py to rerun existing products')
    parser.add_argument('--keepVz', action='store_true',
                        help='Pass --keepVz to refreshties.py to retain .vz and .vz.geodat files')
    parser.add_argument('--check', action='store_true',
                        help='Report which track-*/*_000? dirs have range.offsets.vrt '
                             'but lack velocity/mosaicOffsets.vx; print baseline sigmas')
    parser.add_argument('--year', type=int, nargs='+', default=[2025, 2026], metavar='YYYY',
                        help='One or more years to pass to refreshties.py (default: 2025 2026)')
    parser.add_argument('--noPhase', action='store_true',
                        help='Disable phase+offsets mode (default: NISAR uses --phaseAndOffsets)')
    parser.add_argument('--quadFit', action='store_true',
                        help='Enable the -deltaBQ quadratic baseline correction estimate '
                             '(default: NISAR uses --noQuadFit)')
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
    """Return (proj dict, dem string, template_content dict) read from project.yaml."""
    proj_path = os.path.join(PROJECT_DIR, 'project.yaml')
    proj = {}
    if os.path.exists(proj_path):
        with open(proj_path) as f:
            proj = yaml.safe_load(f) or {}

    dem = ''
    region_path = proj.get('region') or proj.get('regionFile', '')
    for ext in ('', '.yaml'):
        rpath = (region_path + ext) if ext else region_path
        if rpath and os.path.exists(rpath):
            with open(rpath) as f:
                dem = (yaml.safe_load(f) or {}).get('dem', '')
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
    return proj, dem, contents, frame_pattern


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


def setup_track_dirs(track_dirs, copy_files):
    proj, dem, templates, frame_pattern = _load_project_config()
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
                if not os.path.exists(dest):
                    if templates[tmpl_name] is None:
                        print(f'WARNING: template {tmpl_name} not found, skipping {dest}')
                        continue
                    content = _apply_substitutions(templates[tmpl_name], track_num, dem, tiepoint_file)
                    with open(dest, 'w') as f:
                        f.write(content)
                    print(f'Created {dest}')
                    n_created += 1
                else:
                    n_skipped += 1

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
    """Return sigma from az.est.const, rBaseline.deltabp (old text), or rBaseline.deltab.yaml."""
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
                        _extract_sigma(os.path.join(proc_dir, 'motion', 'rBaseline.deltab.yaml')))
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


def run_vel_thumbs(track_dirs):
    def _run(track_dir):
        tpdir = os.path.join(track_dir, 'tiepoints')
        print(f'Running vel_thumbs in {tpdir}')
        subprocess.run(['csh', '-c', 'vel_thumbs vel_thumb_plan'], cwd=tpdir)

    with ThreadPoolExecutor(max_workers=10) as executor:
        executor.map(_run, track_dirs)


def run_velstats_regions(track_dirs):
    for track_dir in track_dirs:
        print(f'Running makevelstatsregions.py in {track_dir}')
        subprocess.run(['csh', '-c', 'makevelstatsregions.py'], cwd=track_dir)


def run_refresh_ties(track_dirs, ties_only, over_write, keep_vz, year,
                      phase_and_offsets=True, no_quad_fit=True, use_yaml=False):
    tracks_str = str([os.path.basename(d) for d in track_dirs]).replace(' ', '')
    tiesOnly = '-tiesOnly ' if ties_only else ''
    overWrite = '--overWrite ' if over_write else ''
    keepVz = '--keepVz ' if keep_vz else ''
    phaseAndOffsets = '--phaseAndOffsets ' if phase_and_offsets else ''
    noQuadFit = '--noQuadFit ' if no_quad_fit else ''
    useYaml = '--yaml ' if use_yaml else ''
    years_str = ' '.join(str(y) for y in year)
    cmd = f'refreshties.py {tiesOnly}{overWrite}{keepVz}{phaseAndOffsets}{noQuadFit}{useYaml}-toRun="{tracks_str}" {years_str} -noPrompt'
    print(f'Running: {cmd}')
    subprocess.run(['csh', '-c', cmd], cwd=PROJECT_DIR)


def main():
    args = parse_args()

    if args.check:
        check_products()
        return

    track_dirs = get_track_dirs(args.tracks)

    setup_track_dirs(track_dirs, args.copyFiles)

    if args.runVelThumbs:
        run_vel_thumbs(track_dirs)

    if args.runVelstatsregions:
        run_velstats_regions(track_dirs)
        return

    run_refresh_ties(track_dirs, args.tiesOnly, args.overWrite, args.keepVz, args.year,
                     phase_and_offsets=not args.noPhase,
                     no_quad_fit=not args.quadFit,
                     use_yaml=True)


if __name__ == '__main__':
    main()
