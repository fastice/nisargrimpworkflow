#!/usr/bin/env python3
"""Generate the azimuth ionosphere offset correction for NISAR frames.

The ionosphere shifts targets in azimuth in proportion to the ALONG-TRACK
gradient of its range delay: a slow-time-linear phase perturbation is a Doppler
shift, and azimuth compression turns a Doppler shift into a misregistration.

    azCorrection [SLC az pix] = gain * kTheory * d(ionCorrection)/d(row)
    kTheory = R * PRF^2 * rSLPixSize / (vEff^2 * deltaA)

This module only PRODUCES the correction raster and names it in the offsets
VRT, exactly as `estimateIonosphere` does for the range screen.  It does not
decide whether to use it: `azparams` fits with and without and records the
winner in `az.est.yaml`, and `mosaic3d -useAzIonosphere` honours that decision.

The screen is rebuilt from the raw unfilled estimate at `smoothScale` times the
smoothing `estimateIonosphere` ships with, because a gradient consumer needs far
more smoothing than a value consumer does -- differentiation multiplies every
wavenumber by k.

The derivative is taken from MEASURED pixels only (normalised convolution) and
the result filled afterwards.  Filling first and then differentiating fabricates
gradients inside the gaps, and the median Greenland screen is only ~66 % real
data; getting this backwards roughly halved the benefit.

    azIonoCorrection --project <dir>                     # every merged frame
    azIonoCorrection --project <dir> --frames track-1/1234_0000
"""
import argparse
import glob
import json
import os
import re
import shutil
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from osgeo import gdal
from scipy.ndimage import gaussian_filter, gaussian_filter1d

from nisargrimpworkflow.estimateIonosphere import _pyramid_fill

gdal.UseExceptions()

# --- defaults, fixed for Greenland ----------------------------------------
GAIN = 0.4                  # multiple of the theoretical scale (see below)
SMOOTH_SCALE = 8.0          # multiple of the shipped screen smoothing
DERIV_SIGMA_ROWS = 2.0      # short: the screen already carries the smoothing
VRT_KEY = 'ionosphereAzimuthOffsetCorrection'
CORRECTION_SUFFIX = '.ionosphereCorrection.azimuthOffset'

# The smoothing estimateIonosphere ships, in pixels of the 26x16 multi-look
# grid, converted to metres with that grid's own pixel sizes.  Note the grid is
# named range-looks first: "26x16" is 26 RANGE looks by 16 AZIMUTH looks.
DEFAULT_SIGMA_AZ_PX = 10.0
DEFAULT_SIGMA_RG_PX = 30.0
SLC_RANGE_M = 1.5614190520833333
SLC_AZIMUTH_M = 4.473391818775636
ML_RANGE_LOOKS = 26
ML_AZIMUTH_LOOKS = 16


def defaultSigmaMetres():
    """The shipped (sigmaAz, sigmaRg) in metres -- (715.74, 1217.91)."""
    return (DEFAULT_SIGMA_AZ_PX * ML_AZIMUTH_LOOKS * SLC_AZIMUTH_M,
            DEFAULT_SIGMA_RG_PX * ML_RANGE_LOOKS * SLC_RANGE_M)


def geodatParams(geodatFile):
    """Geometry needed for the theoretical scale factor."""
    p = json.load(open(geodatFile))['properties']
    nSV = p['NumberOfStateVectors']
    sv = np.array([p[f'SV_Vel_{i}'] for i in range(1, nSV + 1)])
    vSat = np.linalg.norm(sv, axis=1).mean()
    vGround = p['SLCAzimuthPixelSize'] * p['PRF']
    return dict(prf=p['PRF'],
                rSLPixSize=p['SLCRangePixelSize'],
                aSLPixSize=p['SLCAzimuthPixelSize'],
                centerRange=p['MLCenterRange'],
                vEff=np.sqrt(vSat * vGround))


def theoreticalScale(geo, deltaA):
    """kTheory: SLC azimuth pixels per (range pixel per offset-grid row)."""
    return (geo['centerRange'] * geo['prf'] ** 2 * geo['rSLPixSize']
            / (geo['vEff'] ** 2 * deltaA))


def findScreen(frameDir):
    """Locate the ionosphere screen on the offsets grid.

    Two product vintages exist.  Newer frames carry the raw per-pixel estimate
    (`...ionosphereCorrectionUnfilled.offset.tif`), which is what we want: it
    lets the screen be smoothed from scratch to the length the gradient needs.
    Older frames only carry the delivered screen, already filled and smoothed at
    the shipped resolution; for those the extra smoothing needed to reach the
    same total is applied on top, sigmaExtra = sigmaShipped * sqrt(k^2 - 1),
    exact for Gaussians and near-exact here since those screens are gap-free.

    Returns (path, kind, stem) with kind in {'raw', 'smoothed'}.
    """
    raw = glob.glob(os.path.join(frameDir, '*ionosphereCorrectionUnfilled.offset.tif'))
    if raw:
        return (raw[0], 'raw',
                os.path.basename(raw[0]).replace('.ionosphereCorrectionUnfilled.offset.tif', ''))
    for suffix in ('.ionosphereCorrection.globalFill.offset.tif',
                   '.ionosphereCorrection.offset.tif'):
        hits = sorted(glob.glob(os.path.join(frameDir, '*' + suffix)))
        if hits:
            return hits[0], 'smoothed', os.path.basename(hits[0]).replace(suffix, '')
    return None, None, None


def offsetsMeta(frameDir):
    ds = gdal.Open(os.path.join(frameDir, 'azimuth.offsets.vrt'))
    meta = ds.GetMetadata()
    need = ('r0', 'a0', 'deltaR', 'deltaA')
    missing = [k for k in need if k not in meta]
    if missing:
        raise RuntimeError(f'azimuth.offsets.vrt missing {missing}')
    out = dict(na=ds.RasterYSize, nr=ds.RasterXSize,
               r0=float(meta['r0']), a0=float(meta['a0']),
               deltaA=float(meta['deltaA']), deltaR=float(meta['deltaR']))
    del ds
    return out


def geodatName(frameDir):
    hits = sorted(glob.glob(os.path.join(frameDir, 'geodat*.geojson')))
    hits = [h for h in hits if '.secondary.' not in h]
    return hits[0] if hits else None


def maskedDerivative(raw, sigRow, sigCol, derivSigma, weightMin=0.05):
    """Along-track derivative from VALID data only, then gap-filled.

    Normalised convolution (data and mask through the same kernel, then
    divided) so only measured pixels contribute; differentiate that; fill only
    what is left, where no real data lay within the kernel.
    """
    valid = np.isfinite(raw)
    w = valid.astype(np.float64)
    dv = np.where(valid, raw, 0.0)
    num = gaussian_filter(dv, (sigRow, sigCol), mode='nearest', truncate=3.0)
    den = gaussian_filter(w, (sigRow, sigCol), mode='nearest', truncate=3.0)
    supported = den > weightMin
    smooth = np.where(supported, num / np.maximum(den, 1e-12), np.nan)
    smooth = _pyramid_fill(smooth.astype(np.float32), supported).astype(np.float64)
    deriv = gaussian_filter1d(smooth, derivSigma, axis=0, order=1,
                              mode='nearest', truncate=4.0)
    return _pyramid_fill(np.where(supported, deriv, np.nan).astype(np.float32),
                         supported).astype(np.float64)


def placeOnOffsetsGrid(screenPath, screen, frameDir, na, nr):
    """Put a screen smaller than the offsets grid onto that grid.

    A dual-bandwidth partial sub-frame carries no unwrapped phase at 7700 MHz
    and hence no ionosphere estimate, leaving the merged screen short along
    track.  Align by geotransform and leave the uncovered part NaN -- never
    zero: a zero fill puts a step at the coverage edge and the along-track
    derivative reads that step as an enormous gradient.

    Returns (placed, covered, note), or (None, None, why).
    """
    sd = gdal.Open(screenPath)
    od = gdal.Open(os.path.join(frameDir, 'azimuth.offsets.vrt'))
    sg, og = sd.GetGeoTransform(), od.GetGeoTransform()
    del sd, od
    if sg[2] or sg[4] or og[2] or og[4]:
        return None, None, 'rotated/sheared geotransform'
    if (abs(sg[1] - og[1]) > 1e-6 * abs(og[1])
            or abs(sg[5] - og[5]) > 1e-6 * abs(og[5])):
        return None, None, (f'pixel size differs: screen {sg[1]:g}x{sg[5]:g} '
                            f'offsets {og[1]:g}x{og[5]:g}')
    col = (sg[0] - og[0]) / og[1]
    row = (sg[3] - og[3]) / og[5]
    if abs(col - round(col)) > 0.01 or abs(row - round(row)) > 0.01:
        return None, None, f'screen not on the offsets grid (row {row:.3f}, col {col:.3f})'
    row, col = int(round(row)), int(round(col))
    nAs, nRs = screen.shape
    if row < 0 or col < 0 or row + nAs > na or col + nRs > nr:
        return None, None, (f'screen {screen.shape} at ({row}, {col}) does not fit '
                            f'offsets ({na}, {nr})')
    placed = np.full((na, nr), np.nan)
    placed[row:row + nAs, col:col + nRs] = screen
    covered = np.zeros((na, nr), dtype=bool)
    covered[row:row + nAs, col:col + nRs] = True
    return placed, covered, (f'placed {nAs}x{nRs} at row {row}, col {col} of '
                             f'{na}x{nr} ({100.0 * nAs * nRs / (na * nr):.0f}% covered)')


def writeCorrection(templateTif, data, outTif, outVrt, extraMeta):
    src = gdal.Open(templateTif)
    ds = gdal.GetDriverByName('GTiff').Create(
        outTif, src.RasterXSize, src.RasterYSize, 1, gdal.GDT_Float32,
        options=['COMPRESS=DEFLATE', 'PREDICTOR=3', 'ZLEVEL=6', 'BIGTIFF=IF_SAFER'])
    ds.SetGeoTransform(src.GetGeoTransform())
    ds.SetProjection(src.GetProjection())
    band = ds.GetRasterBand(1)
    band.SetDescription(VRT_KEY)
    band.SetNoDataValue(float('nan'))
    band.WriteArray(data.astype(np.float32))
    ds.SetMetadata({k: str(v) for k, v in extraMeta.items()})
    ds.FlushCache()
    del ds, src
    vrt = gdal.Translate(outVrt, outTif, format='VRT')
    del vrt


def addVrtMetadata(vrtPath, entries, backupSuffix='.preAzIono'):
    """Insert/replace dataset-level MDI keys in a VRT, keeping one backup."""
    txt = open(vrtPath).read()
    backup = vrtPath + backupSuffix
    if not os.path.exists(backup):
        shutil.copy(vrtPath, backup)
    for key, value in entries.items():
        line = f'    <MDI key="{key}">{value}</MDI>'
        pat = re.compile(rf'\s*<MDI key="{re.escape(key)}">[^<]*</MDI>')
        if pat.search(txt):
            txt = pat.sub('\n' + line, txt, count=1)
        else:
            m = re.search(r'(<Metadata>\n)', txt)
            if m is None:
                raise RuntimeError(f'{vrtPath}: no dataset <Metadata> block')
            txt = txt[:m.end()] + line + '\n' + txt[m.end():]
    tmp = vrtPath + '.tmp'
    open(tmp, 'w').write(txt)
    os.replace(tmp, vrtPath)


def isStale(corrTif, screenPath, offsetsTif):
    """True when the correction predates the DATA it was derived from.

    Compare against azimuth.offsets.TIF, never the .vrt: stamping the VRT key is
    the last thing processFrame does, so the .vrt is always newer than the
    correction it names and would report every frame stale forever.  The .vrt
    mtime is not evidence the offsets changed at all -- 21 % of Greenland frames
    have a .vrt newer than their az fit purely because this module stamped them,
    with the .tif untouched since July.
    """
    if not os.path.exists(corrTif):
        return True
    tc = os.path.getmtime(corrTif)
    for src in (screenPath, offsetsTif):
        if os.path.exists(src) and os.path.getmtime(src) > tc + 60:
            return True
    return False


def processFrame(frameDir, gain, smoothScale, overwrite, stamp=True):
    rel = frameDir
    rec = dict(frame=os.path.basename(frameDir.rstrip('/')), status='ok',
               gain=gain, smoothScale=smoothScale)
    try:
        screenPath, screenKind, stem = findScreen(frameDir)
        if screenPath is None:
            return dict(rec, status='noScreen')
        rec['screenKind'] = screenKind
        corrTif = os.path.join(frameDir, stem + CORRECTION_SUFFIX + '.tif')
        corrVrt = corrTif.replace('.tif', '.vrt')
        offsetsVrt = os.path.join(frameDir, 'azimuth.offsets.vrt')
        if not os.path.exists(offsetsVrt):
            return dict(rec, status='noOffsets')
        offsetsTif = os.path.join(frameDir, 'azimuth.offsets.tif')
        if not overwrite and not isStale(corrTif, screenPath, offsetsTif):
            rec['status'] = 'upToDate'
            return rec

        gd = geodatName(frameDir)
        if gd is None:
            return dict(rec, status='noGeodat')
        meta = offsetsMeta(frameDir)
        geo = geodatParams(gd)
        kTheory = theoreticalScale(geo, meta['deltaA'])
        rec['kTheory'] = kTheory

        azShip, rgShip = defaultSigmaMetres()
        if screenKind == 'smoothed':
            # the shipped smoothing is already in the product; add only the
            # remainder so the total still comes to smoothScale x shipped
            extra = np.sqrt(max(smoothScale ** 2 - 1.0, 0.0))
            applyAz, applyRg = extra * azShip, extra * rgShip
        else:
            applyAz, applyRg = smoothScale * azShip, smoothScale * rgShip

        ds = gdal.Open(screenPath)
        raw = ds.ReadAsArray().astype(np.float64)
        del ds
        covered = None
        if raw.shape != (meta['na'], meta['nr']):
            raw, covered, note = placeOnOffsetsGrid(screenPath, raw, frameDir,
                                                    meta['na'], meta['nr'])
            if raw is None:
                return dict(rec, status='gridMismatch', note=note)
            rec['screenPlaced'] = note
        validFrac = float(np.mean(np.isfinite(raw)))
        rec['validFrac'] = validFrac
        if validFrac < 0.01:
            return dict(rec, status='emptyScreen')

        dAz = meta['deltaA'] * geo['aSLPixSize']
        dRg = meta['deltaR'] * geo['rSLPixSize']
        deriv = maskedDerivative(raw, applyAz / dAz, applyRg / dRg, DERIV_SIGMA_ROWS)
        correction = gain * kTheory * deriv
        if covered is not None:
            # nothing measured outside the footprint; leave the correction at
            # zero rather than let the fill extrapolate a derivative there --
            # the convention autocleanNISAR uses for the range screen
            correction = np.where(covered, correction, 0.0)

        template = offsetsVrt if covered is not None else screenPath
        writeCorrection(template, correction, corrTif, corrVrt,
                        dict(ionosphereGain=gain, ionosphereSmoothScale=smoothScale,
                             kTheory=f'{kTheory:.4f}', screenKind=screenKind))
        rec['correctionFile'] = os.path.basename(corrVrt)
        if stamp:
            addVrtMetadata(offsetsVrt, {VRT_KEY: os.path.basename(corrVrt)})
    except Exception as exc:                                    # noqa: BLE001
        rec.update(status='failed', error=f'{type(exc).__name__}: {exc}',
                   trace=traceback.format_exc(limit=3))
    return rec


def findFrames(project, pattern):
    hits = sorted(glob.glob(os.path.join(project, pattern)))
    return [h for h in hits if os.path.isdir(h)
            and os.path.exists(os.path.join(h, 'azimuth.offsets.vrt'))]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--project', default='.', help='project directory [.]')
    ap.add_argument('--pattern', default='track-*/[0-9]*_0[0-9][0-9][0-9]',
                    help='frame glob relative to the project directory')
    ap.add_argument('--frames', nargs='*', default=None,
                    help='explicit frame list (relative or absolute), overrides --pattern')
    ap.add_argument('--gain', type=float, default=GAIN)
    ap.add_argument('--smoothScale', type=float, default=SMOOTH_SCALE)
    ap.add_argument('--overwrite', action='store_true',
                    help='recompute even when the correction is newer than its inputs')
    ap.add_argument('--noStamp', action='store_true',
                    help='write the raster but do not name it in azimuth.offsets.vrt')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default=None, help='write a per-frame CSV summary here')
    args = ap.parse_args()

    project = os.path.abspath(args.project)
    if args.frames:
        frames = [f if os.path.isabs(f) else os.path.join(project, f) for f in args.frames]
    else:
        frames = findFrames(project, args.pattern)
    print(f'project {project}')
    print(f'{len(frames)} frames; gain {args.gain}, smoothing {args.smoothScale:g}x, '
          f'stamp={not args.noStamp}', flush=True)

    results, done = [], 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(processFrame, f, args.gain, args.smoothScale,
                            args.overwrite, not args.noStamp): f for f in frames}
        for fut in as_completed(futs):
            rec = fut.result()
            results.append(rec)
            done += 1
            print(f'[{done}/{len(frames)}] {rec["frame"]:20s} {rec["status"]}'
                  + (f'  {rec.get("error", "")}' if rec['status'] == 'failed' else ''),
                  flush=True)

    counts = {}
    for r in results:
        counts[r['status']] = counts.get(r['status'], 0) + 1
    print(f'\n{counts}')
    if args.out:
        import pandas as pd
        pd.DataFrame(results).to_csv(args.out, index=False)
        print(f'wrote {args.out}')
    return 0 if counts.get('failed', 0) == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
