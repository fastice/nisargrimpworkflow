#!/usr/bin/env python3
"""Redo only the full-swath ionosphere gap fill for existing NISAR virtual frames, then rerun the
range-baseline fit (motion/estR -> rparams), without reprocessing RUNW/ROFF.

For each virtual-frame directory (track-N/<orbit>_<vf> with frames.txt and an assembled
*.ionosphereCorrectionUnfilled.offset.tif):
  1. back up everything the refill overwrites into <dir>/ionoFillBackup/<oldVersion>/ (and the
     same under each sub-frame dir) -- the global-fill tiles and VRTs, motion/rBaseline*.yaml;
  2. run SetupNISAR.globalFillIonosphere with --ionoFill (default laplace): new sub-frame tiles,
     the baked virtual-frame correction, its ionoFillVersion tag and, for laplace, the fill-error
     map <stem>.ionosphereCorrection.sigma.tif;
  3. rerun motion/estR so rparams refits the baseline with/without the new correction.
The azimuth correction is not touched: azIonoCorrection reads the raw unfilled screen.

    refillIonosphere --project /path/newGreenlandProject [--tracks 69 131] [--frames 5531_0001]
                     [--ionoFill laplace] [--noEstR] [--dryRun] [--nWorkers 4]

Restore a directory: copy ionoFillBackup/<version>/* back over the frame (and sub-frame) dirs.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor

import yaml
from osgeo import gdal

gdal.UseExceptions()

TILE_PATTERNS = ['*.ionosphereCorrection.globalFill.offset.tif', '*.ionosphereCorrection.globalFill.offset.vrt',
                 '*.ionosphereCorrection.sigma.tif']
MOTION_PATTERNS = ['motion/rBaseline*.yaml', 'motion/rBaseline.deltabp', 'motion/rBaseline', 'motion/rBaseline.quad']


def parseArgs(argv=None):
    p = argparse.ArgumentParser(description='\033[1mRedo the NISAR ionosphere global fill and the range '
                                'baseline fit for existing virtual frames\033[0m')
    p.add_argument('--project', required=True, help='Project root (holds project.yaml and track-* dirs)')
    p.add_argument('--tracks', nargs='*', type=int, default=None, help='Tracks to process (default all)')
    p.add_argument('--frames', nargs='*', default=None, help='Only these virtual frames, e.g. 5531_0001')
    p.add_argument('--ionoFill', choices=['laplace', 'pyramid'], default='laplace')
    p.add_argument('--noEstR', action='store_true', help='Refill only; do not rerun motion/estR')
    p.add_argument('--force', action='store_true', help='Refill even if already at this fill version')
    p.add_argument('--dryRun', action='store_true', help='List what would be done')
    p.add_argument('--nWorkers', type=int, default=4)
    return p.parse_args(argv)


def pairStem(vfDir):
    '''<orbit1>_<vf>.<orbit2>_<vf> of the pair the frame now holds (from its pairinfo).'''
    orbit, vf = os.path.basename(vfDir).split('_')
    hits = glob.glob(os.path.join(vfDir, f'{orbit}.*.pairinfo'))
    if len(hits) != 1:
        return None
    return f"{orbit}_{vf}.{os.path.basename(hits[0]).split('.')[1]}_{vf}"


def checkFrame(vfDir):
    '''(ok, reason): the frame's own unfilled screen exists and its range correction and
    baseline (if any) point at this pair -- not at files left from an earlier pairing.'''
    stem = pairStem(vfDir)
    if stem is None:
        return False, 'noPairinfo'
    if not glob.glob(os.path.join(vfDir, f'{stem}.ionosphereCorrectionUnfilled.offset.tif')):
        return False, 'noScreenForPair'
    tag = gdal.Open(os.path.join(vfDir, 'range.offsets.vrt')).GetMetadataItem('ionosphereRangeOffsetCorrection')
    if tag and not tag.startswith(stem + '.'):
        return False, f'rangeOffsetsTagOtherPair({tag})'
    y = os.path.join(vfDir, 'motion', 'rBaseline.deltabp.yaml')
    if os.path.exists(y):
        for line in open(y):
            if line.startswith('offsetCorrectionFile:'):
                f = line.split()[1]
                if f not in ('nil', 'null') and not f.startswith(stem + '.'):
                    return False, f'baselineUsesOtherPair({f})'
    return True, ''


def currentVersion(vfDir):
    '''ionoFillVersion of the frame's own correction (the one range.offsets.vrt names);
    pyramid-1 if untagged.'''
    tag = gdal.Open(os.path.join(vfDir, 'range.offsets.vrt')).GetMetadataItem('ionosphereRangeOffsetCorrection')
    path = os.path.join(vfDir, tag) if tag else None
    if path is None or not os.path.exists(path):
        stem = pairStem(vfDir)
        hits = sorted(glob.glob(os.path.join(vfDir, f'{stem}.ionosphereCorrection.globalFill.offset.vrt')))
        if not hits:
            return None
        path = hits[0]
    ds = gdal.Open(path)
    v = ds.GetMetadataItem('ionoFillVersion') if ds is not None else None
    return v or 'pyramid-1'


def virtualFrames(project, tracks=None, frames=None):
    out = []
    for tdir in sorted(glob.glob(os.path.join(project, 'track-*'))):
        try:
            track = int(os.path.basename(tdir).split('-')[1])
        except ValueError:
            continue
        if tracks and track not in tracks:
            continue
        for vf in sorted(glob.glob(os.path.join(tdir, '*_0[0-9][0-9][0-9]'))):
            name = os.path.basename(vf)
            if frames and name not in frames:
                continue
            if (os.path.exists(os.path.join(vf, 'frames.txt'))
                    and glob.glob(os.path.join(vf, '*.ionosphereCorrectionUnfilled.offset.tif'))
                    and glob.glob(os.path.join(vf, '*.ionosphereCorrection.globalFill.offset.vrt'))):
                out.append(vf)
    return out


def backup(vfDir, subDirs, version):
    for d, patterns in [(vfDir, TILE_PATTERNS + MOTION_PATTERNS)] + [(s, TILE_PATTERNS) for s in subDirs]:
        dest = os.path.join(d, 'ionoFillBackup', version)
        for pat in patterns:
            for f in glob.glob(os.path.join(d, pat)):
                rel = os.path.relpath(f, d)
                os.makedirs(os.path.dirname(os.path.join(dest, rel)), exist_ok=True)
                if not os.path.exists(os.path.join(dest, rel)):
                    shutil.copy2(f, os.path.join(dest, rel))


def refillOne(job):
    vfDir, project, method, runEstR, force = job
    trackDir, name = os.path.split(vfDir)
    orbit, vf = name.split('_')
    rec = dict(frame=f'{os.path.basename(trackDir)}/{name}', status='ok')
    try:
        from nisargrimpworkflow import SetupNISAR
        from nisargrimpworkflow.estimateIonosphere import IONO_FILL_VERSION
        ok, why = checkFrame(vfDir)
        if not ok:
            return dict(rec, status='skipped', error=why)
        old = currentVersion(vfDir)
        rec['oldVersion'] = old
        if old == IONO_FILL_VERSION[method] and not force:
            return dict(rec, status='upToDate')
        frames = open(os.path.join(vfDir, 'frames.txt')).read().split()
        subDirs = [os.path.join(trackDir, f'{orbit}_{f}') for f in frames]
        backup(vfDir, subDirs, old)
        with open(os.path.join(project, 'project.yaml')) as fp:
            proj = yaml.safe_load(fp) or {}
        myArgs = dict(orbit1=orbit, virtualFrame=vf, frames=frames, outputDir='.',
                      regionFile=proj.get('regionFile') or proj.get('region'),
                      sepIceRock=bool(proj.get('sepIceRock')), ionoFill=method,
                      debugIono=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cwd = os.getcwd()
        os.chdir(trackDir)
        try:
            with open(os.path.join(vfDir, 'refillIonosphere.log'), 'w') as log:
                sys.stdout = log
                SetupNISAR.globalFillIonosphere(myArgs)
        finally:
            sys.stdout = sys.__stdout__
            os.chdir(cwd)
        rec['newVersion'] = currentVersion(vfDir)
        if runEstR and os.path.exists(os.path.join(vfDir, 'motion', 'estR')):
            with open(os.path.join(vfDir, 'motion', 'estR.refill.log'), 'w') as log:
                # -f: skip ~/.cshrc, which would reset PATH and could pick a different rparams
                rc = subprocess.call(['csh', '-f', 'estR'], cwd=os.path.join(vfDir, 'motion'),
                                     stdout=log, stderr=subprocess.STDOUT,
                                     env=dict(os.environ, OPENBLAS_NUM_THREADS='1'))
            rec['estR'] = rc
    except Exception as exc:                                          # noqa: BLE001
        rec.update(status='failed', error=f'{type(exc).__name__}: {exc}',
                   trace=traceback.format_exc(limit=3))
    return rec


def main(argv=None):
    args = parseArgs(argv)
    project = os.path.abspath(args.project)
    vfs = virtualFrames(project, args.tracks, args.frames)
    print(f'{len(vfs)} virtual frames with an assembled ionosphere screen')
    if args.dryRun:
        for v in vfs:
            ok, why = checkFrame(v)
            print(f'  {os.path.relpath(v, project)}  current fill {currentVersion(v)}'
                  + ('' if ok else f'  SKIP: {why}'))
        return 0
    jobs = [(v, project, args.ionoFill, not args.noEstR, args.force) for v in vfs]
    nFail = 0
    with ProcessPoolExecutor(args.nWorkers) as ex:
        for rec in ex.map(refillOne, jobs):
            print(f"{rec['frame']}: {rec['status']} {rec.get('oldVersion', '')} -> "
                  f"{rec.get('newVersion', '')} estR={rec.get('estR', '-')} {rec.get('error', '')}", flush=True)
            nFail += rec['status'] == 'failed'
    return 1 if nFail else 0


if __name__ == '__main__':
    sys.exit(main())
