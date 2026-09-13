#!/usr/bin/env python3
import os
import sys
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
            'years. Use --tracks to restrict processing to specific tracks, or\n'
            '--tieTracks to restrict only the tiepointing pass (leaving the\n'
            'autoclean-mask scan project-wide).\n\n'
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
    parser.add_argument('--tieTracks', nargs='+', metavar='track-N',
                        help='Restrict only the final tiepointing pass to these '
                             'tracks -- the autoclean-mask scan still covers every '
                             'track in --tracks, and any track whose mask is '
                             '(re)built is added to this list automatically (a '
                             'refreshed mask changes the tie fit). Use from a '
                             'nightly update to re-tie just the tracks that got new '
                             'data, instead of all of them. Default: tiepoint every '
                             'track')
    parser.add_argument('--tiesOnly', action='store_true',
                        help='Run refreshties.py without computing velocity mosaics')
    parser.add_argument('--copyFiles', action='store_true',
                        help='Instantiate tie_plan_header and vel_thumb_plan from '
                             'templates/ into every track tiepoints/ directory '
                             '(substituting <TRACK> and <DEM>); skip if already present')
    parser.add_argument('--migrateVelStatsRegions', action='store_true',
                        help='Switch the selected tracks to the track-wide velocityStats '
                             'range (project.yaml velocityStatsRegions: track): move the '
                             'legacy per-frame-number bins to velocityStats/old/ and their '
                             'tiepoints headers/plans to tiepoints/old/, then create the new '
                             'range. Nothing is deleted. Follow with --runVelstatsregions '
                             '--noUpdateVelStats and --updateAutoclean to regenerate.')
    parser.add_argument('--dryRun', action='store_true',
                        help='With --migrateVelStatsRegions: list what would move')
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
                        help='Parallelism for the per-track steps: passed as '
                             '-nThreads to refreshties.py (tracks in parallel) '
                             'and to velocityStats [32]')
    parser.add_argument('--autocleanThreads', type=int, default=None, metavar='N',
                        help='Parallelism for autocleanNISAR.py --tracks (frame pool), '
                             'separately from --nThreads because the two have completely '
                             'different memory profiles: refreshties parallelises cheap '
                             'per-track work, while one autoclean frame job holds ~25 GB '
                             '(a virtual frame is ~3500x24000, so each of the ~10 rasters '
                             'it reads is ~0.6 GB, before numpy temporaries). At 32 that is '
                             '~810 GB and the run is reliably SIGKILLed part-way, losing '
                             'both the remaining frames and every autocleanSummary.yaml. '
                             'Default: autocleanThreads in project.yaml, else 8 (~200 GB)')
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
    parser.add_argument('--useAzIonosphere', dest='useAzIonosphere', action='store_true', default=None,
                        help='Use the azimuth ionosphere correction where one exists '
                             '(azparams keeps the better sigma, mosaic3d applies its '
                             'choice), overriding project.yaml '
                             'applyAzimuthIonosphereCorrection')
    parser.add_argument('--noUseAzIonosphere', dest='useAzIonosphere', action='store_false',
                        help='Never use the azimuth ionosphere correction, overriding '
                             'project.yaml applyAzimuthIonosphereCorrection')
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
    frame_pattern string) read from project.yaml. autocleanThreads is read from
    proj by main() rather than returned here."""
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


def _vel_stats_dirs_for_track(track_dir, frame_pattern, regions='frameBins'):
    """Return sorted list of velocityStats dir names derived from existing virtual-frame dirs.

    regions='frameBins' (default): groups virtual-frame directories (matching
    *_{frame_pattern}) by the first variable digit (M), giving one velocityStats
    dir per M-group spanning {prefix}{M}0...0 – {prefix}{M}9...9.

    regions='track' (project.yaml velocityStatsRegions: track): ONE range per
    track spanning the whole frame pattern ({prefix}0..0 – {prefix}9..9, e.g.
    0000-0099). velocityStats.py then assembles every product of the track on
    the master grid (see velocityStats.assembleTrackStats); the single range
    keeps every '<F1>-<F2>' / '<F1>dash<F2>' parser downstream unchanged.
    """
    prefix = frame_pattern.split('?')[0]          # e.g. '00' from '00??'
    variable_len = frame_pattern.count('?')        # e.g. 2
    tail_len = variable_len - 1                    # digits after M
    if regions == 'track':
        if not glob.glob(os.path.join(track_dir, f'*_{frame_pattern}')):
            return []
        return [f"{prefix}{'0' * variable_len}-{prefix}{'9' * variable_len}"]

    groups = set()
    for d in glob.glob(os.path.join(track_dir, f'*_{frame_pattern}')):
        vf = os.path.basename(d).rsplit('_', 1)[1]  # e.g. '0035'
        m_digit = vf[len(prefix)]                    # e.g. '3'
        start = prefix + m_digit + '0' * tail_len   # e.g. '0030'
        end   = prefix + m_digit + '9' * tail_len   # e.g. '0039'
        groups.add(f'{start}-{end}')
    return sorted(groups)


def migrate_velstats_regions(track_dirs, dry_run=False):
    """Switch tracks from per-frame-number velocityStats bins to the single
    track-wide range (project.yaml velocityStatsRegions: track). Nothing is
    deleted: velocityStats/<F1>-<F2> bins go to velocityStats/old/, and the
    per-bin tiepoints files (vel_thumb_header_<F1>dash<F2>, tie_plan<year>.<F1>dash<F2>,
    vel_thumb_plan_all<year>dash<F1>dash<F2>) go to tiepoints/old/. The next
    setup_track_dirs() call then creates the track-wide range and its header
    (sentinel resolution; velThumbs computes each product's box itself in
    track mode). Re-run the tie stage afterwards so products are regenerated on
    the master grid and velocityStats/autoclean rebuild from them."""
    proj, _, _, _, frame_pattern = _load_project_config()
    regions = str(proj.get('velocityStatsRegions', 'frameBins')).lower()
    if regions != 'track':
        sys.exit('--migrateVelStatsRegions: set velocityStatsRegions: track in project.yaml first')
    keep = set(_vel_stats_dirs_for_track('.', frame_pattern, 'track')
               or [f"{frame_pattern.split('?')[0]}{'0' * frame_pattern.count('?')}-"
                   f"{frame_pattern.split('?')[0]}{'9' * frame_pattern.count('?')}"])
    for track_dir in track_dirs:
        vsd = os.path.join(track_dir, 'velocityStats')
        tpd = os.path.join(track_dir, 'tiepoints')
        bins = [d for d in sorted(glob.glob(os.path.join(vsd, '*-*')))
                if os.path.basename(d) not in keep]
        files = []
        for pat in ('vel_thumb_header_*dash*', 'tie_plan20??.*dash*', 'vel_thumb_plan_all20??dash*'):
            files += [f for f in sorted(glob.glob(os.path.join(tpd, pat)))
                      if not any(f.endswith(k.replace('-', 'dash')) for k in keep)]
        if not bins and not files:
            print(f'{track_dir}: nothing to migrate')
            continue
        print(f'{track_dir}: {len(bins)} bin dir(s) -> velocityStats/old/, '
              f'{len(files)} tiepoints file(s) -> tiepoints/old/'
              f'{"  (dry run)" if dry_run else ""}')
        if dry_run:
            for x in bins + files:
                print(f'    {x}')
            continue
        for d in bins:
            dest = os.path.join(vsd, 'old')
            os.makedirs(dest, exist_ok=True)
            shutil.move(d, os.path.join(dest, os.path.basename(d)))
        for f in files:
            dest = os.path.join(tpd, 'old')
            os.makedirs(dest, exist_ok=True)
            shutil.move(f, os.path.join(dest, os.path.basename(f)))


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

        regions = str(proj.get('velocityStatsRegions', 'frameBins')).lower()
        regions = 'track' if regions == 'track' else 'frameBins'
        if regions == 'track':
            # Legacy per-digit bins would make the range lookup ambiguous
            # (autocleanNISAR scans velocityStats/*-*); refuse until they are
            # moved aside (setupNISARTracks --migrateVelStatsRegions).
            legacy = [d for d in glob.glob(os.path.join(track_dir, 'velocityStats', '*-*'))
                      if os.path.basename(d) not in _vel_stats_dirs_for_track(track_dir, frame_pattern, 'track')]
            if legacy:
                # Not migrated yet: keep its legacy bins (velocityStats still
                # assembles each bin's members on the master grid and
                # autocleanNISAR scans whatever ranges exist), just do not add
                # the track-wide range beside them.
                print(f'{track_dir}: velocityStatsRegions is "track" but legacy frame-bin '
                      f'dirs remain ({[os.path.basename(d) for d in legacy]}); run '
                      f'setupNISARTracks --migrateVelStatsRegions --tracks '
                      f'{os.path.basename(track_dir)} to switch this track')
                continue
        for dir_name in _vel_stats_dirs_for_track(track_dir, frame_pattern, regions):
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
    # Track-mode composites (velocityStatsRegions: track) hold the full
    # regional velocity map plus one 8-Mpx strip of accumulators per process
    # (~12 GB peak measured on Antarctica40/track-82), so cap the parallelism
    # separately from the tie stage's --nThreads; velStatsThreads in
    # project.yaml overrides the default of 8 (~100 GB).
    proj, _, _, _, _ = _load_project_config()
    if str(proj.get('velocityStatsRegions', 'frameBins')).lower() == 'track':
        maxWorkers = min(maxWorkers, int(proj.get('velStatsThreads', 8)))
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


def autoclean_threads(proj, args):
    """Frame-pool size for autocleanNISAR: --autocleanThreads, else project.yaml's
    autocleanThreads, else 8. Kept apart from --nThreads because one autoclean frame job
    holds ~25 GB of rasters -- at 32 threads that is ~810 GB, more than either machine has,
    and the run is SIGKILLed part-way through."""
    if args.autocleanThreads is not None:
        return args.autocleanThreads
    return int(proj.get('autocleanThreads', 8))


def tie_track_dirs(track_dirs, needs_autoclean, tie_tracks):
    """
    Track dirs for the final tiepointing pass.

    Without --tieTracks this is every track in track_dirs -- the historical
    behaviour, where a nightly run re-tied all ~170 tracks even when only a
    handful had new data. With it, only the named tracks, unioned with the
    tracks whose autoclean mask was just (re)built: a refreshed mask changes
    the tie fit, so those have to be redone whether or not they carry new
    data.

    Named tracks that are not in track_dirs (not in this project, or dropped
    for having no virtual-frame data) are reported and skipped rather than
    failing the run -- unlike --tracks, this list is generated upstream from
    released granules, so a stale name should not abort the nightly update.
    """
    if not tie_tracks:
        return track_dirs
    byName = {os.path.basename(d): d for d in track_dirs}
    wanted = set(tie_tracks)
    missing = sorted(wanted - set(byName))
    if missing:
        print(f'--tieTracks: skipping {len(missing)} track(s) not in the '
              f'project or with no virtual-frame data: {missing}')
    keep = (wanted & set(byName)) | {os.path.basename(d)
                                     for d in needs_autoclean}
    return [d for d in track_dirs if os.path.basename(d) in keep]


def run_refresh_ties(track_dirs, ties_only, over_write, keep_vz, year,
                      phase_and_offsets=True, no_quad_fit=True, use_yaml=False,
                      use_squint=False, use_az_ionosphere=False, n_threads=None):
    tracks_str = str([os.path.basename(d) for d in track_dirs]).replace(' ', '')
    tiesOnly = '-tiesOnly ' if ties_only else ''
    overWrite = '--overWrite ' if over_write else ''
    keepVz = '--keepVz ' if keep_vz else ''
    phaseAndOffsets = '--phaseAndOffsets ' if phase_and_offsets else ''
    noQuadFit = '--noQuadFit ' if no_quad_fit else ''
    useYaml = '--yaml ' if use_yaml else ''
    useSquint = '--useSquint ' if use_squint else ''
    useAzIon = '--useAzIonosphere ' if use_az_ionosphere else ''
    nThreadsFlag = f'-nThreads {n_threads} ' if n_threads else ''
    years_str = ' '.join(str(y) for y in year)
    cmd = f'refreshties.py {tiesOnly}{overWrite}{keepVz}{phaseAndOffsets}{noQuadFit}{useYaml}{useSquint}{useAzIon}{nThreadsFlag}-toRun="{tracks_str}" {years_str} -noPrompt'
    print(f'Running: {cmd}')
    r = subprocess.run(['csh', '-c', cmd], cwd=PROJECT_DIR)
    # Surface a failing tie stage rather than continuing to makeMaster as if it
    # had worked -- refreshties now exits non-zero when a track's tie chain
    # fails, and swallowing it here would keep it out of the nightly summary.
    if r.returncode != 0:
        u.mywarning(f'setupNISARTracks: refreshties.py exited '
                    f'{r.returncode} -- see the per-track stderr files')
    return r.returncode


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

    if args.migrateVelStatsRegions:
        migrate_velstats_regions(track_dirs, dry_run=args.dryRun)
        if args.dryRun:
            return

    setup_track_dirs(track_dirs, args.copyFiles, overwrite=args.overWrite)

    if args.copyFiles or args.migrateVelStatsRegions:
        return

    projectUseSquint = bool(proj.get('applySquintCorrection', False))
    useSquint = projectUseSquint if args.useSquint is None else args.useSquint
    projectUseAzIon = bool(proj.get('applyAzimuthIonosphereCorrection', False))
    useAzIonosphere = projectUseAzIon if args.useAzIonosphere is None else args.useAzIonosphere

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
                         use_az_ionosphere=useAzIonosphere,
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
    # mask embedded (or all of them, with --updateAutoclean). The mask scan
    # stays project-wide even under --tieTracks, since a track can need a
    # first-time mask without having new data this run; the tiepointing below
    # then covers every track that got one, because a fresh mask changes the
    # tie fit.
    needs_autoclean = [d for d in track_dirs
                       if args.updateAutoclean or not _track_has_autoclean_mask(d, frame_pattern)]
    if needs_autoclean:
        # --noUpdateVelStats reuses velStats built by a prior --runVelstatsregions
        # run, so autoclean consumes it directly without a redundant rebuild.
        if not args.noUpdateVelStats:
            run_vel_stats(needs_autoclean, region_path, vel_map,
                          n_threads=args.nThreads)
        # Deliberately NOT args.nThreads -- see --autocleanThreads: a frame job here
        # holds ~25 GB, so the tie-stage thread count over-subscribes memory and the
        # run gets killed.
        run_autoclean(needs_autoclean, frame_pattern,
                      n_threads=autoclean_threads(proj, args),
                      new=not args.updateAutoclean)

    tie_dirs = tie_track_dirs(track_dirs, needs_autoclean, args.tieTracks)
    if not tie_dirs:
        print('Nothing to tiepoint: --tieTracks named no track with data, and '
              'no track needed an autoclean mask')
        return 0
    if len(tie_dirs) != len(track_dirs):
        print(f'Tiepointing {len(tie_dirs)} of {len(track_dirs)} track(s) '
              f'(--tieTracks + {len(needs_autoclean)} needing autoclean)')

    return run_refresh_ties(tie_dirs, args.tiesOnly, args.overWrite,
                            args.keepVz, args.year,
                            phase_and_offsets=not args.noPhase,
                            no_quad_fit=not args.quadFit,
                            use_yaml=True,
                            use_squint=useSquint,
                            use_az_ionosphere=useAzIonosphere,
                            n_threads=args.nThreads)


if __name__ == '__main__':
    sys.exit(main() or 0)
