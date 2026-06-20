#!/usr/bin/env python3
"""
estimateIonosphere.py

Estimate ionospheric correction from NISAR RUNW phase and independent range
offsets (VRT).

Physics:
    For positive differential TEC (ΔTEC = TEC_sec − TEC_ref > 0):
      • range-offset group delay increases apparent range:
            offset_m  = +K·ΔTEC/f²  > 0
            offset_phase = −4π/λ · offset_m  < 0
      • phase-delay shortens the phase path (opposite sign to group delay),
        but the −4π/λ factor in the interferometric phase convention makes
        the contribution the same sign:
            phase contribution  = −4π/λ · (K·ΔTEC/f²)  < 0
      Both terms are negative, so they ADD rather than cancel:
            phase + offset_phase = 2 × iono   (iono < 0 when ΔTEC > 0)
    → iono  = smooth((phase + offset_phase) / 2)

    Corrections (consumers ADD each correction term):
      Δφ_correction   = −iono             →  φ_corrected  = φ  + Δφ_correction
      ΔR_correction   = +λ/(4π) · iono   →  ΔR_corrected = ΔR + ΔR_correction  (metres)
      ΔR_pix_correction = ΔR_correction / slp_spacing = +λ/(4π·slp) · iono
    (for ΔTEC > 0: iono < 0, so ΔR_correction < 0 — subtracts iono-inflated range)

Steps:
  1. Load RUNW (unwrapped phase + connected components).
  2. Load range offsets from VRT; resample to RUNW grid via geotransforms.
  3. Convert offsets to phase:  offset_phase = −4π/λ × offset_m.
  4. Estimate smooth iono from (phase + offset_phase) / 2.
  5. Write two GeoTIFFs + assembled multi-band VRT with named bands:
       Band 1  ionosphereCorrection  (radians)
       Band 2  unwrappedPhase        (radians)
     Also writes *.ionosphereCorrection.offset.tif/.vrt on the native ROFF
     grid in SLC pixels (same units as range.offsets).

VRT metadata linkage:
    SetupNISAR.py stamps the key ionosphereRangeOffsetCorrection = <basename>
    on range.offsets.vrt after this script runs.  rparams and mosaic3d read
    that key when they open range.offsets.vrt and use it to locate and load
    the offset-grid correction file (*.ionosphereCorrection.offset.vrt).

Usage:
    estimateIonosphere.py RUNW.h5 offsets.vrt output.vrt [options]
"""

import argparse
import itertools
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import ndimage
from scipy.ndimage import gaussian_filter

try:
    from sarfunc.defaultRegionDefs import defaultRegionDefs
except ImportError:
    sys.exit("sarfunc package not found — add it to PYTHONPATH")

try:
    import nisarhdf
except ImportError:
    sys.exit("nisarhdf package not found — add it to PYTHONPATH")

try:
    from osgeo import gdal, osr
    gdal.UseExceptions()
except ImportError:
    sys.exit("GDAL (osgeo) not available")

from .ambiguityAnalysis import load_vel_sim, ambiguity_table, apply_ambiguity_correction
from .ROFFtoGrimp import updateSimVrtGeotransforms


# ---------------------------------------------------------------------------
# Progress spinner
# ---------------------------------------------------------------------------

class _Spinner:
    def __init__(self, msg):
        self._msg = msg
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._spin, daemon=True)

    def _spin(self):
        for ch in itertools.cycle(r'|/-\\'):
            if self._stop.is_set():
                break
            print(f'\r  {self._msg} {ch}', end='', flush=True)
            time.sleep(0.1)

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        self._thread.join()
        print(f'\r  {self._msg} done.    ')


# ---------------------------------------------------------------------------
# I/O — RUNW  (borrowed verbatim from fixConnectedComponents.py)
# ---------------------------------------------------------------------------

def load_runw(path, frame):
    print(f"Loading RUNW: {path}")
    runw = nisarhdf.nisarRUNWHDF()
    kwargs = {}
    if frame is not None:
        kwargs['frame'] = frame
    runw.openHDF(path, **kwargs)
    return runw


def apply_runw_mask(phase, runw):
    try:
        mask_ds = runw.h5[runw.product][runw.bands][runw.frequency][runw.productType]['mask']
    except KeyError:
        print("  No interferogram mask found — skipping.")
        return phase
    mask = np.asarray(mask_ds)
    n_bad = int((mask == 0).sum())
    n_total = mask.size
    print(f"  Interferogram mask: {n_bad:,} / {n_total:,} pixels masked "
          f"({100.0 * n_bad / n_total:.2f}%)")
    phase[mask == 0] = np.nan
    return phase


# ---------------------------------------------------------------------------
# I/O — offset VRT
# ---------------------------------------------------------------------------

def load_geom_vrt(path):
    """Read the 'RangeOffsets' band from an offset-geometry VRT (SLP units)."""
    print(f"Loading offset geometry VRT: {path}")
    ds = gdal.Open(path)
    if ds is None:
        sys.exit(f"Cannot open offset geometry VRT: {path}")
    gt = ds.GetGeoTransform()
    geom_band = None
    for i in range(1, ds.RasterCount + 1):
        b = ds.GetRasterBand(i)
        if b.GetMetadataItem('Description') == 'RangeOffsets':
            geom_band = b
            break
    if geom_band is None:
        sys.exit(f"No band named 'RangeOffsets' found in {path}")
    data = geom_band.ReadAsArray().astype(np.float64)
    nodata = geom_band.GetNoDataValue()
    ds = None
    data[data < -1e9] = np.nan
    if nodata is not None:
        data[data == nodata] = np.nan
    print(f"  Geometry grid: {data.shape[0]} az × {data.shape[1]} rg")
    return data, gt


def load_offset_vrt(path):
    print(f"Loading offset VRT: {path}")
    ds = gdal.Open(path)
    if ds is None:
        sys.exit(f"Cannot open offset VRT: {path}")
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray().astype(np.float64)
    nodata = band.GetNoDataValue()
    ds = None
    data[data < -1e9] = np.nan
    if nodata is not None:
        data[data == nodata] = np.nan
    print(f"  Offset grid: {data.shape[0]} az × {data.shape[1]} rg")
    return data, gt


# ---------------------------------------------------------------------------
# I/O — mask file
# ---------------------------------------------------------------------------

def load_mask_file(path, runw):
    """Load mask GeoTIFF/VRT; return boolean array on RUNW grid (True = use pixel).

    mask=1 pixels are included in the iono fit; all others are excluded.
    If the mask dimensions match the RUNW grid it is used directly (pixel-aligned),
    which handles masks whose geotransform is in pixel rather than geographic coords.
    Otherwise geotransform-based resampling is applied.
    """
    print(f"Loading mask file: {path}")
    ds = gdal.Open(path)
    if ds is None:
        sys.exit(f"Cannot open mask file: {path}")
    data = ds.GetRasterBand(1).ReadAsArray()
    file_gt = ds.GetGeoTransform()
    ds = None
    print(f"  Mask grid: {data.shape[0]} az × {data.shape[1]} rg  "
          f"(RUNW: {runw.MLAzimuthSize} az × {runw.MLRangeSize} rg)")
    if data.shape == (runw.MLAzimuthSize, runw.MLRangeSize):
        print(f"  Shapes match — applying pixel-for-pixel.")
        mask = data == 1
    else:
        print(f"  Shapes differ — resampling via geotransforms.")
        resampled = resample_vrt_to_runw(data.astype(np.float64), file_gt, runw)
        mask = resampled > 0.5
    n_use = int(mask.sum())
    print(f"  Mask: {n_use:,} / {mask.size:,} pixels will be used for iono fit")
    return mask


def load_ice_rock_mask(path, runw):
    """Load 3-value mask (0=water, 1=rock, 2=ice); return (ice_mask, rock_mask)
    as boolean arrays on the RUNW ML grid.  Nearest-neighbor resampling is used
    to preserve the categorical values."""
    print(f"Loading ice/rock mask: {path}")
    ds = gdal.Open(path)
    if ds is None:
        sys.exit(f"Cannot open ice/rock mask: {path}")
    data = ds.GetRasterBand(1).ReadAsArray()
    file_gt = ds.GetGeoTransform()
    ds = None
    print(f"  Mask grid: {data.shape[0]} az × {data.shape[1]} rg  "
          f"(RUNW: {runw.MLAzimuthSize} az × {runw.MLRangeSize} rg)")
    if data.shape == (runw.MLAzimuthSize, runw.MLRangeSize):
        values = data.astype(np.int32)
    else:
        values = resample_vrt_to_runw(
            data.astype(np.float64), file_gt, runw, method='nearest'
        ).round().astype(np.int32)
    ice_mask = (values == 2)
    rock_mask = (values == 1)
    print(f"  sepIceRock: {int(ice_mask.sum()):,} ice px, "
          f"{int(rock_mask.sum()):,} rock px on RUNW grid")
    return ice_mask, rock_mask


# ---------------------------------------------------------------------------
# Grid resampling — VRT → RUNW via geotransforms
# ---------------------------------------------------------------------------

def resample_vrt_to_runw(offset_m, vrt_gt, runw, method='linear'):
    """Resample VRT data to the RUNW ML grid.

    Both grids share the same geographic CRS; the geotransforms encode the
    pixel-centre-to-coordinate mapping for each.
    `method` is passed to RegularGridInterpolator ('linear' or 'nearest').
    """
    runw_gt = runw.getGeoTransform(tiff=False)  # SAR time axis (positive dy)
    nr = runw.MLAzimuthSize
    nc = runw.MLRangeSize

    x0, dx = runw_gt[0], runw_gt[1]
    y0, dy = runw_gt[3], runw_gt[5]
    x_centres = x0 + (np.arange(nc) + 0.5) * dx   # shape (nc,)
    y_centres = y0 + (np.arange(nr) + 0.5) * dy   # shape (nr,)

    vx0, vdx = vrt_gt[0], vrt_gt[1]
    vy0, vdy = vrt_gt[3], vrt_gt[5]
    vrt_col = (x_centres - vx0) / vdx - 0.5   # fractional VRT col for each RUNW col
    vrt_row = (y_centres - vy0) / vdy - 0.5   # fractional VRT row for each RUNW row

    vrt_nr, vrt_nc = offset_m.shape
    interp = RegularGridInterpolator(
        (np.arange(vrt_nr), np.arange(vrt_nc)),
        offset_m,
        method=method, bounds_error=False, fill_value=np.nan,
    )

    vrt_row_grid = np.tile(vrt_row[:, np.newaxis], (1, nc))   # (nr, nc)
    vrt_col_grid = np.tile(vrt_col[np.newaxis, :], (nr, 1))   # (nr, nc)
    return interp((vrt_row_grid, vrt_col_grid)).astype(np.float32)


def resample_runw_to_vrt(data, runw, vrt_gt, vrt_shape):
    """Bilinearly resample RUNW-grid data back to the VRT offset grid."""
    runw_gt = runw.getGeoTransform(tiff=False)  # SAR time axis (positive dy)
    vrt_nr, vrt_nc = vrt_shape

    vx0, vdx = vrt_gt[0], vrt_gt[1]
    vy0, vdy = vrt_gt[3], vrt_gt[5]
    vrt_x = vx0 + (np.arange(vrt_nc) + 0.5) * vdx   # (vrt_nc,)
    vrt_y = vy0 + (np.arange(vrt_nr) + 0.5) * vdy   # (vrt_nr,)

    x0, dx = runw_gt[0], runw_gt[1]
    y0, dy = runw_gt[3], runw_gt[5]
    runw_col = (vrt_x - x0) / dx - 0.5   # fractional RUNW col for each VRT col
    runw_row = (vrt_y - y0) / dy - 0.5   # fractional RUNW row for each VRT row

    interp = RegularGridInterpolator(
        (np.arange(runw.MLAzimuthSize), np.arange(runw.MLRangeSize)),
        data.astype(np.float64),
        method='linear', bounds_error=False, fill_value=np.nan,
    )

    row_grid = np.tile(runw_row[:, np.newaxis], (1, vrt_nc))   # (vrt_nr, vrt_nc)
    col_grid = np.tile(runw_col[np.newaxis, :], (vrt_nr, 1))   # (vrt_nr, vrt_nc)
    return interp((row_grid, col_grid)).astype(np.float32)


# ---------------------------------------------------------------------------
# Gap fill + smooth  (pyramid fill borrowed from cleanIonosphereCorr.py)
# ---------------------------------------------------------------------------

_FILL_KERNEL = np.array([[0.25, 0.5, 0.25],
                         [0.5,  0.0, 0.5],
                         [0.25, 0.5, 0.25]], dtype=np.float32)
_FILL_KERNEL /= _FILL_KERNEL.sum()


def _nn_fill(arr, valid):
    _, idx = ndimage.distance_transform_edt(~valid, return_indices=True)
    filled = arr.copy()
    inv = ~valid
    filled[inv] = arr[idx[0][inv], idx[1][inv]]
    return filled


def _conv_fill(arr, valid, iterations, boundary_mode='reflect'):
    filled = arr.copy()
    for _ in range(iterations):
        smoothed = ndimage.convolve(filled, _FILL_KERNEL, mode=boundary_mode)
        filled = np.where(valid, arr, smoothed)
    return filled


def _pyramid_fill(arr, valid, levels=4, iterations_per_level=50, boundary_mode='reflect'):
    arrs, valids = [arr], [valid]
    for _ in range(levels - 1):
        arrs.append(ndimage.zoom(arrs[-1], 0.5, order=1))
        valids.append(ndimage.zoom(valids[-1].astype(np.float32), 0.5, order=0) > 0.5)
    filled = _nn_fill(arrs[-1], valids[-1])
    filled = _conv_fill(filled, valids[-1], iterations_per_level * 2,
                        boundary_mode=boundary_mode)
    for lvl in range(levels - 2, -1, -1):
        zoom_factors = [o / c for o, c in zip(arrs[lvl].shape, filled.shape)]
        upsampled = ndimage.zoom(filled, zoom_factors, order=1)
        upsampled = upsampled[:arrs[lvl].shape[0], :arrs[lvl].shape[1]]
        merged = np.where(valids[lvl], arrs[lvl], upsampled)
        filled = _conv_fill(merged, valids[lvl], iterations_per_level,
                            boundary_mode=boundary_mode)
    return filled


def fill_and_smooth_iono(raw_iono, valid, sigma, boundary_mode='reflect'):
    """Fill cc=0 / NaN gaps with pyramid fill, then Gaussian smooth."""
    f32 = raw_iono.astype(np.float32)
    f32[~valid] = np.nan
    if valid.any():
        f32 = _nn_fill(f32, valid)   # seed invalid pixels so zoom never touches NaN
    filled = _pyramid_fill(f32, valid, boundary_mode=boundary_mode)
    return gaussian_filter(filled.astype(np.float64), sigma).astype(np.float32)


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def estimate_and_correct(phase, cc, offset_phase, sigma, user_mask=None):
    raw_iono = (offset_phase.astype(np.float64) + phase.astype(np.float64)) / 2.0
    valid_iono = np.isfinite(raw_iono) & (cc != 0)
    if user_mask is not None:
        valid_iono &= user_mask
    iono_final = fill_and_smooth_iono(raw_iono, valid_iono, sigma)
    iono_mean = float(np.nanmean(iono_final))
    iono_final -= iono_mean
    print(f"  Iono zero-mean adjustment: {iono_mean:+.4f} rad removed.")
    iono_unfilled = raw_iono.astype(np.float32)
    iono_unfilled[~valid_iono] = np.nan
    iono_unfilled -= iono_mean
    phase_final = (phase.astype(np.float64) - iono_final).astype(np.float32)
    return iono_final.astype(np.float32), phase_final, iono_unfilled


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_geotiff(path, data, gt, nodata=float('nan'), proj=None):
    nrows, ncols = data.shape
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(path, ncols, nrows, 1, gdal.GDT_Float32,
                       ['COMPRESS=LZW', 'BIGTIFF=NO'])
    ds.SetGeoTransform(gt)
    if proj:
        ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    arr = data.astype(np.float32)
    if not np.isnan(nodata):
        arr[~np.isfinite(arr)] = nodata
    band.WriteArray(arr)
    ds.FlushCache()
    ds = None
    print(f"Written: {path}")


def write_output_vrt(vrt_path, band_files, band_names, gt):
    src_ds = gdal.Open(band_files[0])
    xsize = src_ds.RasterXSize
    ysize = src_ds.RasterYSize
    proj = src_ds.GetProjection()
    src_ds = None

    driver = gdal.GetDriverByName('VRT')
    vrt_ds = driver.Create(vrt_path, xsize, ysize, len(band_files), gdal.GDT_Float32)
    vrt_ds.SetGeoTransform(gt)
    if proj:
        vrt_ds.SetProjection(proj)

    vrt_dir = os.path.dirname(os.path.abspath(vrt_path))
    for idx, (src_path, name) in enumerate(zip(band_files, band_names), start=1):
        rel = os.path.relpath(os.path.abspath(src_path), vrt_dir)
        band = vrt_ds.GetRasterBand(idx)
        band.SetDescription(name)
        _src = gdal.Open(src_path)
        src_nd = _src.GetRasterBand(1).GetNoDataValue() if _src else None
        _src = None
        band.SetNoDataValue(src_nd if src_nd is not None else float('nan'))
        source_xml = (
            f'<ComplexSource>\n'
            f'    <SourceFilename relativeToVRT="1">{rel}</SourceFilename>\n'
            f'    <SourceBand>1</SourceBand>\n'
            f'    <SrcRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>\n'
            f'    <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>\n'
            f'</ComplexSource>'
        )
        band.SetMetadataItem('source_0', source_xml, 'new_vrt_sources')

    vrt_ds = None
    print(f"Written: {vrt_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Estimate ionospheric correction from NISAR RUNW phase '
                    'and independent range offsets (VRT).',
        epilog='Part of the nisargrimpworkflow package.')
    p.add_argument('runw', help='NISAR RUNW HDF5 file')
    p.add_argument('offset_vrt',
                   help='Range offset VRT (band 1, metres, geographic coords)')
    p.add_argument('output', help='Output VRT path (2-band named bands)')
    p.add_argument('--frame', type=int, default=None, metavar='N')
    p.add_argument('--regionFile', default=None, metavar='YAML',
                   help='Region YAML file for defaultRegionDefs; '
                        'defaults to built-in greenland definition')
    p.add_argument('--verticalCorrection', default=None, metavar='FILE',
                   help='xyDEM grid (m/yr) of submergence/emergence rate; '
                        'passed to the velSim siminsar call as '
                        '-verticalCorrection so its DC anchor matches '
                        'correctedUnwrappedPhase [None]')
    p.add_argument('--overWrite', action='store_true',
                   help='Force siminsar to run even if velSim/maskVel already exists')
    p.add_argument('--velThresh', type=float, default=100.0, metavar='M/YR',
                   help='Velocity threshold for maskVel siminsar call (default: 100 m/yr)')
    p.add_argument('--sigma-az', type=float, default=10.0, metavar='PX',
                   help='Azimuth Gaussian sigma for iono smoothing (default: 10 px)')
    p.add_argument('--sigma-rg', type=float, default=30.0, metavar='PX',
                   help='Range Gaussian sigma for iono smoothing (default: 30 px)')
    p.add_argument('--offset-geometry', default='offsets.geom.vrt', metavar='VRT',
                   help='Offset geometry VRT with RangeOffsets band to subtract '
                        '(default: offsets.geom.vrt)')
    p.add_argument('--maskFile', default=None, metavar='FILE',
                   help='GeoTIFF/VRT mask; only mask=1 pixels used for iono estimation '
                        '(default: maskVel.vrt, created by siminsar if absent)')
    p.add_argument('--sepIceRock', action='store_true', default=False,
                   help='Restrict the per-frame (phase+offset_phase)/2 ionosphere '
                        'estimate to ice pixels only (rock/water excluded, since '
                        'the estimate assumes nonzero velocity). The excluded rock '
                        'pixels provide an absolute "actual vs simulated offset" '
                        'reference that SetupNISAR.globalFillIonosphere() uses '
                        '(once, at the virtual-frame level) to anchor the '
                        'otherwise-floating ice-derived field. Requires '
                        '--iceRockMask or auto-detected '
                        'offsetSims/offsets.geom.mask.vrt (0=water, 1=rock, 2=ice).')
    p.add_argument('--iceRockMask', default=None, metavar='FILE',
                   help='3-value mask VRT/GeoTIFF for --sepIceRock (0=water, 1=rock, 2=ice). '
                        'Default: auto-detect offsetSims/offsets.geom.mask.vrt in simDir.')
    p.add_argument('--simOffsets', default=None, metavar='VRT',
                   help='VRT of simulated range offsets (same format/units as offset_vrt, '
                        'SLC pixels) for per-CC diagnostic columns in the ambiguity table '
                        '(default: offsets.vrt if present, else skipped)')
    p.add_argument('--simDir', default='.', metavar='DIR',
                   help='Directory (relative to cwd) where siminsar outputs '
                        '(velSim, maskVel) are written.  Created if absent. '
                        '(default: current directory)')
    p.add_argument('--noInterp', action='store_true',
                   help='Skip intfloat hole-filling of correctedUnwrappedPhase')
    p.add_argument('--interpThresh', type=int, default=20, metavar='N',
                   help='intfloat: maximum hole area to fill (pixels) [20]')
    p.add_argument('--islandThresh', type=int, default=20, metavar='N',
                   help='intfloat: maximum isolated-island area to remove '
                        'after filling (pixels) [20]')
    p.add_argument('--referencePhase', default=None, metavar='VRT',
                   help='Simulated reference phase VRT (e.g. from SetupNISAR siminsar output)')
    p.add_argument('--mask', default=None, metavar='FILE',
                   help='Mask file (GeoTIFF/VRT) to apply during processing')
    p.add_argument('--phaseThresh', type=float, default=14.0 * np.pi, metavar='RAD',
                   help='Mask correctedUnwrappedPhase where '
                        '|correctedPhase - simPhase| >= THRESH radians. '
                        'Screens regions with likely incorrect unwrapping. '
                        'Default: 14π rad.')
    p.add_argument('--noPhaseThreshPass', action='store_true', default=False,
                   help='Disable the second-pass iono re-estimation using the '
                        'phase-residual mask (default: second pass is ON when '
                        '--phaseThresh is active).')
    p.add_argument('--outputAll', action='store_true',
                   help='Write all intermediate bands (5-band VRT + stem-named TIFs). '
                        'Default: write only correctedUnwrappedPhase.tif, '
                        'ionosphereCorrection.tif, ionosphereCorrection.offset.tif '
                        'each with a single-band VRT sidecar.')
    p.add_argument('--debugIono', action='store_true',
                   help='Write the ambiguity-corrected (pre-iono) unwrapped phase '
                        'to a debug/ subdirectory alongside the main outputs. '
                        'Intended to be set by SetupNISAR --debugIono.')
    p.add_argument('--verbose', action='store_true',
                   help='Print progress messages')
    p.add_argument('--minTol', type=float, default=None,
                   help='Variable smoothing-radius map (additional pass on top of the '
                   'intfloat hole-filling above), applied to correctedUnwrappedPhase: '
                   'm/yr floor for the adaptive tolerance clip(percentSpeed/100*speed, '
                   'minTol, maxTol). Required together with --percentSpeed/--maxTol to '
                   'enable the map [None]')
    p.add_argument('--percentSpeed', type=float, default=None,
                   help='Variable smoothing-radius map: percent of local speed (e.g. '
                   '1 = 1%%) used in the adaptive tolerance. Required together with '
                   '--minTol/--maxTol [None]')
    p.add_argument('--maxTol', type=float, default=None,
                   help='Variable smoothing-radius map: m/yr ceiling for the adaptive '
                   'tolerance. Required together with --minTol/--percentSpeed [None]')
    p.add_argument('--maxSmoothRadius', type=int, default=50,
                   help='Variable smoothing-radius map: sweep cap in single-look pixels, '
                   'clamped to <= 255 (byte output) [50]')
    p.add_argument('--smoothNIter', type=int, default=3,
                   help='Variable smoothing-radius map: repeated box-filter passes per '
                   'sweep step (Gaussian-ish) [3]')
    p.add_argument('--noVariableSmoothing', action='store_true', default=False,
                   help='Disable the variable smoothing-radius map even if '
                   '--minTol/--percentSpeed/--maxTol (or project.yaml) supply values')
    args = p.parse_args()
    smoothFlags = [args.minTol is not None, args.percentSpeed is not None,
                  args.maxTol is not None]
    if any(smoothFlags) and not all(smoothFlags):
        sys.exit('estimateIonosphere: --minTol/--percentSpeed/--maxTol must be given '
                 'together')
    if args.maxSmoothRadius > 255:
        print(f'WARNING: --maxSmoothRadius {args.maxSmoothRadius} exceeds byte range, '
              'clamping to 255')
        args.maxSmoothRadius = 255
    return args


def run_intfloat(phase_tif, phase_vrt, thresh, island_thresh):
    """Hole-fill correctedUnwrappedPhase in-place via intfloat.

    intfloat reads the entire TIF into memory before writing, so
    -inputVRT and -tiff can safely share the same basename:
      read(phase.tif) → GDALClose → interpolate → write(phase.tif + phase.vrt)
    The single-band phase.vrt sidecar is distinct from the 5-band output VRT.
    """
    command = ['intfloat', '-wdist',
               '-inputVRT', phase_tif,
               '-tiff', phase_vrt,
               '-thresh', str(thresh),
               '-islandThresh', str(island_thresh),
               phase_tif]
    print(f"  {' '.join(command)}")
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


def run_vel_sim(runw, region, simDir='.', overwrite=False, verticalCorrection=None,
                smoothParams=None):
    '''
    smoothParams : dict or None
        If given, {'minTol', 'percentSpeed', 'maxTol', 'maxSmoothRadius', 'smoothNIter'}
        passed through to siminsar to additionally produce velSim.smr(.vrt), the variable
        smoothing-radius map applied to correctedUnwrappedPhase later in main().
    '''
    geodat = f'geodat{runw.NumberRangeLooks}x{runw.NumberAzimuthLooks}.geojson'
    velSimPath = os.path.join(simDir, 'velSim')
    # Re-run if velSim is missing, or if a radius map is now wanted but wasn't produced
    # by an earlier (pre-smoothing) run.
    smrMissing = smoothParams is not None and not os.path.exists(velSimPath + '.smr.vrt')
    if os.path.exists(velSimPath) and not overwrite and not smrMissing:
        print(f"  {velSimPath} already exists — skipping")
        return
    os.makedirs(simDir, exist_ok=True)
    dT = runw.dT
    # siminsar requires its last 4 tokens to be demFile/displacementFile/sceneFile/
    # outputFile -- all optional flags must come before them, never after.
    cmd = ['siminsar', '-dT', str(dT), '-velOnly', '-velocity']
    if verticalCorrection is not None:
        cmd += ['-verticalCorrection', verticalCorrection]
    if smoothParams is not None:
        cmd += ['-minTol', str(smoothParams['minTol']),
               '-percentSpeed', str(smoothParams['percentSpeed']),
               '-maxTol', str(smoothParams['maxTol']),
               '-maxSmoothRadius', str(smoothParams['maxSmoothRadius']),
               '-smoothNIter', str(smoothParams['smoothNIter'])]
    cmd += [region.dem(), region.velMap(), geodat, velSimPath]
    with _Spinner('Creating velSim'):
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    updateSimVrtGeotransforms(os.path.join(simDir, 'velSim*.vrt'), runw)


def run_mask_vel(runw, region, vel_thresh, simDir='.', overwrite=False):
    maskVelVrt = os.path.join(simDir, 'maskVel.vrt')
    if os.path.exists(maskVelVrt) and not overwrite:
        print(f"  {maskVelVrt} already exists — skipping")
        return
    os.makedirs(simDir, exist_ok=True)
    geodat = f'geodat{runw.NumberRangeLooks}x{runw.NumberAzimuthLooks}.geojson'
    cmd = ['siminsar', '-velThresh', str(vel_thresh), '-velocity',
           region.dem(), region.velMap(), geodat, os.path.join(simDir, 'maskVel')]
    with _Spinner('Creating maskVel'):
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    updateSimVrtGeotransforms(os.path.join(simDir, 'maskVel*.vrt'), runw)


def main():
    args = parse_args()

    runw = load_runw(args.runw, args.frame)
    simDir = args.simDir

    if args.regionFile:
        region = defaultRegionDefs(None, regionFile=args.regionFile)
    else:
        region = defaultRegionDefs('greenland')
    smoothParams = None
    if args.minTol is not None and not args.noVariableSmoothing:
        smoothParams = {'minTol': args.minTol, 'percentSpeed': args.percentSpeed,
                        'maxTol': args.maxTol, 'maxSmoothRadius': args.maxSmoothRadius,
                        'smoothNIter': args.smoothNIter}
    run_vel_sim(runw, region, simDir=simDir, overwrite=args.overWrite,
               verticalCorrection=args.verticalCorrection, smoothParams=smoothParams)
    run_mask_vel(runw, region, args.velThresh, simDir=simDir, overwrite=args.overWrite)

    phase = runw.unwrappedPhase.astype(np.float32)
    print("Applying RUNW interferogram mask...")
    apply_runw_mask(phase, runw)
    cc = runw.connectedComponents
    phase[cc < 1] = np.nan
    # NISAR L2 fill pattern: edge rows have phase=0.0, cc=1, but coh=NaN.
    # The cc mask above won't catch them; mask on NaN coherence instead.
    coh_fill_mask = ~np.isfinite(runw.coherenceMagnitude.astype(np.float32))
    n_coh_masked = int(coh_fill_mask.sum())
    if n_coh_masked:
        print(f"  NaN-coherence fill mask: {n_coh_masked:,} additional pixels masked.")
    phase[coh_fill_mask] = np.nan
    wavelength = float(runw.Wavelength)

    print(f"\nWavelength: {wavelength:.6f} m")
    print(f"RUNW grid:  {runw.MLAzimuthSize} az × {runw.MLRangeSize} rg")
    print(f"Connected components: {len(np.unique(cc[cc != 0]))} (excluding cc=0)")

    slp_spacing = float(runw.SLCRangePixelSize)
    print(f"SL range pixel spacing: {slp_spacing:.4f} m")

    offset_slp, vrt_gt = load_offset_vrt(args.offset_vrt)
    geom_slp, geom_gt = load_geom_vrt(args.offset_geometry)
    offset_m = (offset_slp - geom_slp) * slp_spacing

    print("Resampling offsets to RUNW grid...")
    offset_resampled = resample_vrt_to_runw(offset_m, vrt_gt, runw)
    offset_phase = offset_resampled * (-4.0 * np.pi / wavelength)

    mask_path = args.maskFile or os.path.join(simDir, 'maskVel.vrt')
    user_mask = load_mask_file(mask_path, runw)

    sigma = (args.sigma_az, args.sigma_rg)
    vel_sim = load_vel_sim(os.path.join(simDir, 'velSim.vrt'))

    # Optional: simulated offsets for per-CC diagnostic columns in ambiguity table.
    # Loaded through the same pipeline as the actual offsets (SLC pixels → phase).
    sim_offset_phase = None
    sim_offsets_path = args.simOffsets or ('offsets.vrt' if os.path.exists('offsets.vrt') else None)
    if sim_offsets_path:
        print(f"Loading simulated offsets: {sim_offsets_path}")
        sim_slp, _ = load_offset_vrt(sim_offsets_path)   # same grid as offset_vrt; use vrt_gt
        sim_offset_m = (sim_slp - geom_slp) * slp_spacing
        sim_offset_resampled = resample_vrt_to_runw(sim_offset_m, vrt_gt, runw)
        sim_offset_phase = sim_offset_resampled * (-4.0 * np.pi / wavelength)

    # ------------------------------------------------------------------
    # --sepIceRock: load ice/rock mask and pre-seed rock pixels
    # ------------------------------------------------------------------
    ice_mask_sr = None
    if args.sepIceRock:
        _mask_path = args.iceRockMask
        if _mask_path is None:
            _mask_path = os.path.join(simDir, 'offsets.geom.mask.vrt')
            if not os.path.exists(_mask_path):
                _mask_path = os.path.join(os.path.dirname(simDir), 'offsetSims',
                                          'offsets.geom.mask.vrt')
        if not os.path.exists(_mask_path):
            sys.exit(f'--sepIceRock: cannot find ice/rock mask at {_mask_path}. '
                     f'Supply --iceRockMask.')
        ice_mask_sr, _rock_mask_sr = load_ice_rock_mask(_mask_path, runw)
        # Restrict ice estimation to velocity mask AND ice pixels.  The
        # rock-anchored DC correction (using the rock pixels excluded here)
        # is applied once, globally, in SetupNISAR.globalFillIonosphere().
        if user_mask is not None:
            user_mask = user_mask & ice_mask_sr
        else:
            user_mask = ice_mask_sr

    # ------------------------------------------------------------------
    # Iterative ambiguity correction + ionosphere estimation
    #
    # After iono correction the iono absorbs ~half of each unwrapping
    # error, so the apparent n halves each iteration.  The loop converges
    # in ceil(log2(N_max)) passes or fewer.
    # ------------------------------------------------------------------
    MAX_PASSES = 8
    n_accumulated = {}

    for pass_num in range(1, MAX_PASSES + 1):
        print(f"\n--- Ionosphere estimation (pass {pass_num}) ---")
        iono, phase_out, iono_unfilled = estimate_and_correct(
            phase, cc, offset_phase, sigma, user_mask=user_mask,
        )

        # Offset residual: iono-corrected actual offset vs simulated offset
        off_resid = None
        if sim_offset_phase is not None:
            # Scale phase residual (radians) to SLC pixel units
            rad_to_pix = wavelength / (4.0 * np.pi * slp_spacing)
            off_resid = ((offset_phase + iono) - sim_offset_phase) * rad_to_pix

        print(f"\nAmbiguity diagnostic table (pass {pass_num}):")
        n_map = ambiguity_table(phase_out, vel_sim, cc,
                                n_accumulated=n_accumulated,
                                offset_resid=off_resid,
                                offset_units='pix')

        if all(n == 0 for n in n_map.values()):
            print(f"\nConverged after {pass_num} pass(es).")
            break

        # Update running total *after* printing the table for this pass
        for k, n in n_map.items():
            n_accumulated[k] = n_accumulated.get(k, 0) + n

        print(f"\nApplying ambiguity corrections (pass {pass_num})...")
        phase = apply_ambiguity_correction(phase, cc, n_map)
    else:
        print(f"\nWarning: did not fully converge after {MAX_PASSES} passes.")
        for k, n in n_map.items():
            n_accumulated[k] = n_accumulated.get(k, 0) + n

    # ------------------------------------------------------------------
    # Phase-residual threshold mask (optional)
    # ------------------------------------------------------------------
    thresh_mask = None
    if args.phaseThresh is not None:
        # Remove the same DC bias that ambiguity_table removes: mean of
        # (phase_out - vel_sim) over the reference (largest) CC.  Without
        # this correction a global offset between the two fields causes the
        # entire image to exceed even a generous threshold.
        _labels, _counts = np.unique(cc[cc != 0], return_counts=True)
        if len(_labels) > 0:
            _ref_cc = int(_labels[np.argmax(_counts)])
            _ref_mask = (cc == _ref_cc) & np.isfinite(phase_out) & np.isfinite(vel_sim)
            _bias = float(np.mean((phase_out - vel_sim)[_ref_mask])) if _ref_mask.any() else 0.0
        else:
            _bias = 0.0
        thresh_mask = np.abs(phase_out - vel_sim - _bias) >= args.phaseThresh
        n_masked = int(np.sum(thresh_mask & np.isfinite(phase_out)))
        print(f"\nPhase threshold mask (|correctedPhase - simPhase - bias| >= "
              f"{args.phaseThresh:.2f} rad, bias={_bias:+.3f} rad): "
              f"{n_masked:,} pixels masked.")
        phase_out = phase_out.copy()
        phase_out[thresh_mask] = np.nan

    # ------------------------------------------------------------------
    # Second pass: re-estimate iono using phaseThresh mask as gating.
    # Pixels that already agreed with vel_sim in Pass 1 (residual <
    # phaseThresh) form the new mask.  This is more principled than the
    # velocity mask: it includes fast areas whose ambiguity was correctly
    # resolved and excludes only pixels with clear residual anomalies.
    # Uses the ambiguity-corrected `phase` from Pass 1 (not iono-subtracted
    # phase_out) so the iono surface is re-estimated fresh from a better set.
    # ------------------------------------------------------------------
    if (thresh_mask is not None and not args.noPhaseThreshPass
            and (cc != 0).any()):
        pass2_mask = ~thresh_mask & np.isfinite(phase_out) & (cc != 0)
        n_pass2 = int(pass2_mask.sum())
        print(f"\n--- Ionosphere estimation (Pass 2: phaseThresh mask, "
              f"{n_pass2:,} pixels) ---")
        if n_pass2 > 0:
            phase2 = phase.copy()
            n_accumulated2 = {}
            for pass_num2 in range(1, MAX_PASSES + 1):
                iono2, phase_out2, iono_unfilled2 = estimate_and_correct(
                    phase2, cc, offset_phase, sigma, user_mask=pass2_mask)
                n_map2 = ambiguity_table(phase_out2, vel_sim, cc,
                                         n_accumulated=n_accumulated2)
                for k, n in n_map2.items():
                    n_accumulated2[k] = n_accumulated2.get(k, 0) + n
                if all(n == 0 for n in n_map2.values()):
                    print(f"  Pass 2 converged after {pass_num2} sub-pass(es).")
                    break
                phase2 = apply_ambiguity_correction(phase2, cc, n_map2)
            else:
                print(f"  Pass 2: did not fully converge after {MAX_PASSES} passes.")
            thresh_mask2 = np.abs(phase_out2 - vel_sim - _bias) >= args.phaseThresh
            n_masked2 = int(np.sum(thresh_mask2 & np.isfinite(phase_out2)))
            phase_out2 = phase_out2.copy()
            phase_out2[thresh_mask2] = np.nan
            print(f"  Pass 2 phaseThresh masked {n_masked2:,} pixels.")
            iono, phase_out, iono_unfilled, thresh_mask = (
                iono2, phase_out2, iono_unfilled2, thresh_mask2)
        else:
            print("  Pass 2: no valid pixels after phaseThresh mask — keeping Pass 1.")

    # ------------------------------------------------------------------
    # Write final (converged) results as output
    # ------------------------------------------------------------------
    gt = runw.getGeoTransform(tiff=False)   # positive dy (SAR time axis)
    stem = os.path.splitext(os.path.abspath(args.output))[0]
    _sr = osr.SpatialReference()
    _sr.ImportFromEPSG(runw.epsg)
    proj_wkt = _sr.ExportToWkt()

    if args.outputAll:
        iono_tif          = stem + '.ionosphereCorrection.tif'
        phase_tif         = stem + '.correctedUnwrappedPhase.tif'
        phase_vrt         = stem + '.correctedUnwrappedPhase.vrt'
        iono_unfilled_tif = stem + '.ionosphereCorrectionUnfilled.tif'
        offset_phase_tif  = stem + '.offsetPhase.tif'
        unwrapped_tif     = stem + '.unwrappedPhase.tif'
        offset_iono_tif   = stem + '.ionosphereCorrection.offset.tif'
    else:
        # Phase products use the full stem; offset product drops the looks/nisar suffix.
        offset_stem     = re.sub(r'\.\d+x\d+\.nisar$', '', stem)
        iono_tif        = stem + '.ionosphereCorrection.tif'
        phase_tif       = stem + '.correctedUnwrappedPhase.tif'
        phase_vrt       = stem + '.correctedUnwrappedPhase.vrt'
        offset_iono_tif = offset_stem + '.ionosphereCorrection.offset.tif'

    write_geotiff(iono_tif, iono, gt, proj=proj_wkt)
    write_geotiff(phase_tif, phase_out, gt, proj=proj_wkt)
    if not args.noInterp:
        print("\nHole-filling correctedUnwrappedPhase with intfloat ...")
        # Interpolate to a temp file rather than overwriting phase_tif in place, so the
        # variable smoothing-radius map (below) can be applied as a genuinely separate
        # additional pass reading the interpolated result and writing the final phase_tif/
        # phase_vrt, instead of two successive in-place rewrites of the same file.
        interpTif = phase_tif.replace('.tif', '.interp.tmp.tif')
        interpVrt = phase_vrt.replace('.vrt', '.interp.tmp.vrt')
        run_intfloat(phase_tif, interpVrt,
                     thresh=args.interpThresh,
                     island_thresh=args.islandThresh)
        # Re-mask cc=0 and phaseThresh pixels that intfloat may have filled,
        # using -2e9 so the C geocoder sees the expected noData sentinel.
        _ds = gdal.Open(interpTif, gdal.GA_Update)
        if _ds is not None:
            _arr = _ds.GetRasterBand(1).ReadAsArray()
            _arr[cc < 1] = -2.0e9
            if thresh_mask is not None:
                _arr[thresh_mask] = -2.0e9
            _ds.GetRasterBand(1).WriteArray(_arr)
            _ds.GetRasterBand(1).SetNoDataValue(-2.0e9)
            _ds.FlushCache()
            _ds = None
        # Variable smoothing-radius map: an additional pass on top of the fixed
        # interpolation above, using the radius map produced alongside velSim by
        # run_vel_sim() (same siminsar -minTol/-percentSpeed/-maxTol sweep).
        radiusMap = os.path.join(simDir, 'velSim.smr.vrt')
        if smoothParams is not None and os.path.exists(radiusMap):
            print("Applying variable smoothing-radius map...")
            pixRatio = runw.SLCAzimuthPixelSize / runw.SLCRangePixelSize
            command = ['filterfloat', '-inputVRT', interpVrt, '-tiff', phase_vrt,
                      '-radiusMap', radiusMap, '-pixRatio', str(pixRatio),
                      '-nIterations', str(args.smoothNIter), '-minValue', '-2.0e9']
            print(f"  {' '.join(command)}")
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
            if args.debugIono:
                debugDir = os.path.join(os.path.dirname(os.path.abspath(phase_tif)), 'debug')
                os.makedirs(debugDir, exist_ok=True)
                unsmoothedTif = os.path.join(
                    debugDir, os.path.basename(phase_tif).replace('.tif', '.unsmoothed.tif'))
                unsmoothedVrt = unsmoothedTif.replace('.tif', '.vrt')
                shutil.move(interpTif, unsmoothedTif)
                write_output_vrt(unsmoothedVrt, [unsmoothedTif], ['Phase'], gt)
                os.remove(interpVrt)
                print(f"  Saved unsmoothed copy (--debugIono) -> {unsmoothedVrt}")
            else:
                os.remove(interpTif)
                os.remove(interpVrt)
        else:
            os.replace(interpTif, phase_tif)
            # Don't just rename interpVrt -> phase_vrt: its embedded
            # SourceFilename still points at the pre-rename interp.tmp.tif
            # name, which no longer exists, leaving a dangling reference.
            # Write a fresh VRT against the renamed phase_tif instead.
            os.remove(interpVrt)
            write_output_vrt(phase_vrt, [phase_tif], ['Phase'], gt)
        # intfloat/filterfloat write no band description; stamp it now so downstream
        # VRT mosaics (custom_buildvrtWithOffsets) propagate the correct name.
        _ds = gdal.Open(phase_vrt, gdal.GA_Update)
        if _ds is not None:
            _ds.GetRasterBand(1).SetDescription('Phase')
            _ds = None
    elif not args.outputAll:
        write_output_vrt(phase_vrt, [phase_tif], ['Phase'], gt)

    if args.outputAll:
        write_geotiff(iono_unfilled_tif, iono_unfilled, gt, proj=proj_wkt)
        write_geotiff(offset_phase_tif, offset_phase.astype(np.float32), gt, proj=proj_wkt)
        write_geotiff(unwrapped_tif, phase.astype(np.float32), gt, proj=proj_wkt)
        band_names = ['ionosphereCorrection', 'correctedUnwrappedPhase',
                      'ionosphereCorrectionUnfilled', 'offsetPhase', 'unwrappedPhase']
        write_output_vrt(args.output,
                         [iono_tif, phase_tif, iono_unfilled_tif,
                          offset_phase_tif, unwrapped_tif],
                         band_names, gt)
    else:
        write_output_vrt(stem + '.ionosphereCorrection.vrt',
                         [iono_tif], ['ionosphereCorrection'], gt)
        # Unfilled iono on RUNW phase grid (radians, sparse) — debug/comparison only
        iono_unfilled_tif = stem + '.ionosphereCorrectionUnfilled.tif'
        iono_unfilled_vrt = stem + '.ionosphereCorrectionUnfilled.vrt'
        iono_unfilled_out = iono_unfilled.copy()
        iono_unfilled_out[~np.isfinite(iono_unfilled_out)] = -2.0e9
        write_geotiff(iono_unfilled_tif, iono_unfilled_out, gt, nodata=-2.0e9,
                      proj=proj_wkt)
        write_output_vrt(iono_unfilled_vrt, [iono_unfilled_tif],
                         ['ionosphereCorrectionUnfilled'], gt)

    # Ambiguity-corrected-only phase (no iono subtraction) — debug mode only
    if args.debugIono:
        debugDir = os.path.join(os.path.dirname(os.path.abspath(args.output)), 'debug')
        os.makedirs(debugDir, exist_ok=True)
        debugBase = os.path.join(debugDir, os.path.basename(stem))
        ambigTif = debugBase + '.ambiguityCorrectedUnwrappedPhase.tif'
        ambigVrt = debugBase + '.ambiguityCorrectedUnwrappedPhase.vrt'
        _ambig = phase.astype(np.float32)
        _ambig[cc < 1] = -2.0e9
        write_geotiff(ambigTif, _ambig, gt, nodata=-2.0e9, proj=proj_wkt)
        write_output_vrt(ambigVrt, [ambigTif], ['ambiguityCorrectedUnwrappedPhase'], gt)

    # Ionosphere correction regridded to the offset VRT grid, in SLC pixels
    # (range.offsets is in SLC pixels; C geocoder applies this correction directly)
    print("\nRegridding iono correction to offset grid (SLC pixels)...")
    scale_rad_to_pix = wavelength / (4.0 * np.pi * slp_spacing)
    iono_offset_pix = resample_runw_to_vrt(iono, runw, vrt_gt, offset_m.shape)
    iono_offset_pix *= scale_rad_to_pix
    write_geotiff(offset_iono_tif, iono_offset_pix, vrt_gt, nodata=-2.0e9,
                  proj=proj_wkt)

    if args.debugIono:
        _debugDir = os.path.join(os.path.dirname(os.path.abspath(args.output)), 'debug')
        os.makedirs(_debugDir, exist_ok=True)
        _offBase = os.path.splitext(os.path.basename(args.offset_vrt))[0]
        corrTif = os.path.join(_debugDir, _offBase + '.corrected.tif')
        corrVrt = os.path.join(_debugDir, _offBase + '.corrected.vrt')
        write_geotiff(corrTif, (offset_slp + iono_offset_pix).astype(np.float32),
                      vrt_gt, nodata=-2.0e9, proj=proj_wkt)
        write_output_vrt(corrVrt, [corrTif], ['RangeOffsetsCorrected'], vrt_gt)

    if args.outputAll:
        write_output_vrt(stem + '.offset.vrt',
                         [offset_iono_tif], ['ionosphereCorrection'], vrt_gt)
    else:
        write_output_vrt(offset_stem + '.ionosphereCorrection.offset.vrt',
                         [offset_iono_tif], ['ionosphereCorrection'], vrt_gt)
        # Unfilled iono on offset grid (SLC pixels, sparse) — input for global fill
        # Resample sparse RUNW-grid iono to offset grid without spreading NaN:
        # fill source NaN with 0, resample, resample valid mask, re-apply.
        print("Regridding unfilled iono to offset grid (SLC pixels, sparse)...")
        _iu = iono_unfilled.copy()
        _valid_runw = np.isfinite(_iu)
        _iu[~_valid_runw] = 0.0
        _iu_offset = resample_runw_to_vrt(_iu.astype(np.float64), runw, vrt_gt, offset_m.shape)
        _valid_offset = resample_runw_to_vrt(
            _valid_runw.astype(np.float64), runw, vrt_gt, offset_m.shape) > 0.5
        iono_unfilled_offset_pix = (_iu_offset * scale_rad_to_pix).astype(np.float32)
        iono_unfilled_offset_pix[~_valid_offset] = -2.0e9
        iono_unfilled_offset_tif = offset_stem + '.ionosphereCorrectionUnfilled.offset.tif'
        iono_unfilled_offset_vrt = offset_stem + '.ionosphereCorrectionUnfilled.offset.vrt'
        write_geotiff(iono_unfilled_offset_tif, iono_unfilled_offset_pix, vrt_gt,
                      nodata=-2.0e9, proj=proj_wkt)
        write_output_vrt(iono_unfilled_offset_vrt, [iono_unfilled_offset_tif],
                         ['ionosphereCorrectionUnfilled'], vrt_gt)


if __name__ == '__main__':
    main()
