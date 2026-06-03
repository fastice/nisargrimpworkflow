#!/usr/bin/env python3
import os
import glob
import re
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor
import utilities as u

PROJECT_DIR = '/Volumes/insar1/ian/NISAR/realNISAR/greenlandProject'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tiesOnly', action='store_true', help='Run refreshties.py')
    parser.add_argument('--copyFiles', action='store_true', help='Copy template files from track-88 to all tracks')
    parser.add_argument('--runVelThumbs', action='store_true', help='Run vel_thumbs vel_thumb_plan in each tiepoints dir')
    parser.add_argument('--runVelstatsregions', action='store_true', help='Run makevelstatsregions.py in each track dir')
    parser.add_argument('--overWrite', action='store_true', help='Pass --overWrite to refreshties.py to rerun existing products')
    parser.add_argument('--keepVz', action='store_true', help='Pass --keepVz to refreshties.py to retain .vz and .vz.geodat files')
    parser.add_argument('--check', action='store_true', help='Check which track-*/*_000? dirs have range.offsets.vrt but lack velocity/mosaicOffsets.vx')
    parser.add_argument('--year', type=int, nargs='+', default=[2025, 2026], metavar='YYYY', help='One or more years to pass to refreshties.py (default: 2025 2026)')
    return parser.parse_args()


def get_track_dirs():
    return sorted(glob.glob(os.path.join(PROJECT_DIR, 'track-*')),
                  key=lambda p: int(re.search(r'track-(\d+)', p).group(1)))


def setup_track_dirs(track_dirs, copy_files):
    src_header = os.path.join(PROJECT_DIR, 'track-88', 'tiepoints', 'tie_plan_header')
    src_vel_thumb = os.path.join(PROJECT_DIR, 'track-88', 'tiepoints', 'vel_thumb_plan')

    for track_dir in track_dirs:
        track = os.path.basename(track_dir)
        tpdir = os.path.join(track_dir, 'tiepoints')

        if not os.path.isdir(tpdir):
            os.makedirs(tpdir)
            print(f'Created {tpdir}')

        if copy_files:
            dest = os.path.join(tpdir, 'tie_plan_header')
            if not os.path.exists(dest):
                with open(src_header) as f:
                    content = f.read()
                content = content.replace('track-88', track)
                with open(dest, 'w') as f:
                    f.write(content)
                print(f'Created {dest}')
            else:
                print(f'Skipped {dest} (already exists)')

            dest_vel_thumb = os.path.join(tpdir, 'vel_thumb_plan')
            if not os.path.exists(dest_vel_thumb):
                with open(src_vel_thumb) as f:
                    content = f.read()
                content = content.replace('-88', f'-{track.split("-")[1]}')
                with open(dest_vel_thumb, 'w') as f:
                    f.write(content)
                print(f'Created {dest_vel_thumb}')
            else:
                print(f'Skipped {dest_vel_thumb} (already exists)')

        vel_dir = os.path.join(track_dir, 'velocityStats', '0000-0001')
        if not os.path.isdir(vel_dir):
            os.makedirs(vel_dir)
            print(f'Created {vel_dir}')
        else:
            print(f'Skipped {vel_dir} (already exists)')


def _extract_sigma(filepath):
    """Return the sigma*sqrt(X2/n) value from an az.est.const or rBaseline.deltabp file, or None."""
    try:
        with open(filepath) as f:
            for line in f:
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
        rng_sigma    = _extract_sigma(os.path.join(proc_dir, 'motion', 'rBaseline.deltabp'))
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


def run_refresh_ties(track_dirs, ties_only, over_write, keep_vz, year):
    tracks_str = str([os.path.basename(d) for d in track_dirs]).replace(' ', '')
    tiesOnly = '-tiesOnly ' if ties_only else ''
    overWrite = '--overWrite ' if over_write else ''
    keepVz = '--keepVz ' if keep_vz else ''
    years_str = ' '.join(str(y) for y in year)
    cmd = f'refreshties.py {tiesOnly}{overWrite}{keepVz}-toRun="{tracks_str}" {years_str} -noPrompt'
    print(f'Running: {cmd}')
    subprocess.run(['csh', '-c', cmd], cwd=PROJECT_DIR)


def main():
    args = parse_args()

    if args.check:
        check_products()
        return

    track_dirs = get_track_dirs()

    setup_track_dirs(track_dirs, args.copyFiles)

    if args.runVelThumbs:
        run_vel_thumbs(track_dirs)

    if args.runVelstatsregions:
        run_velstats_regions(track_dirs)
        return

    run_refresh_ties(track_dirs, args.tiesOnly, args.overWrite, args.keepVz, args.year)


if __name__ == "__main__":
    main()
