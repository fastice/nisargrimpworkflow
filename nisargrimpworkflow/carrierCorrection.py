"""
Rubbersheet carrier correction for NISAR RUNW phase (opt-in, used by estimateIonosphere --carrierCorrection).

A rubbersheeted RUNW is resampled along-track by the measured shift dA_c (pixelOffsets/alongTrackOffset), which
follows ice motion AND the ionosphere-induced shift. Moving along the squint's Doppler carrier adds -k_dc*dA_c,
which cancels the motion term the squint put into the phase but overshoots by -k_dc*dA_iono. The correction

    C = k_dc * (dA_c - S[dA_model]),    k_dc = 2*pi*f_dc / V_groundtrack

removes that error and keeps the RUNW's zero-Doppler geometry. S is the dense-offset window (box average,
sampled on the offset grid and interpolated like dA_c), so the model removes only the motion the rubbersheet
followed. dA_model comes from the simulated offsets (velocity - geometry, AzimuthOffsets band).
See nisarErrors/Documents/carrierCorrectionFrames.md.
"""
import os
import re

import h5py
import numpy as np
from osgeo import gdal
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

RUNW = 'science/LSAR/RUNW'
# RUNWs from processors before this CRID can carry gross rubbersheet patches (~-70 m, ~-200 m; seen in ISCE3
# 0.25.8 products) that need the frame-wide outlier test as well as the local one.
PATCH_CULL_MAX_CRID = 5032


def granuleCrid(path):
    '''Integer processor CRID from a NISAR granule name (..._P05023_N_F_J_001.h5 -> 5023); None if absent.'''
    for tok in os.path.basename(path).split('_'):
        m = re.fullmatch(r'[A-Za-z](\d{5})', tok)
        if m:
            return int(m.group(1))
    return None


def interpAxes(z, t, r, tq, rq, method='linear'):
    '''Interpolate z(t, r) onto the (tq, rq) grid; NaN outside.'''
    it = RegularGridInterpolator((t, r), z, method=method, bounds_error=False, fill_value=np.nan)
    T, R = np.meshgrid(tq, rq, indexing='ij')
    return it((T, R))


def nearestFill(z):
    '''Fill NaNs from the nearest finite value.'''
    ok = np.isfinite(z)
    if ok.all() or not ok.any():
        return z
    _, idx = ndimage.distance_transform_edt(~ok, return_indices=True)
    return z[idx[0], idx[1]]


def cullOffsets(c, size=9, maxDev=5.0, minGlobal=20.0, nMad=6.0, grow=2, patchCull=True):
    '''Reject gross rubbersheet outliers: values more than maxDev m from the local median and, when patchCull
    (older RUNWs, see PATCH_CULL_MAX_CRID), values more than max(minGlobal m, nMad robust sigma) from the frame
    median -- the ~-70 m / ~-200 m patches are larger than the local-median window. Rejected cells, grown by
    `grow` cells to catch patch edges, become NaN.'''
    ok = np.isfinite(c)
    med = ndimage.median_filter(np.where(ok, c, np.nanmedian(c)), size=size)
    bad = ~ok | (np.abs(c - med) > maxDev)
    keep = c[~bad]
    if patchCull and keep.size:
        gMed = np.median(keep)
        gSig = 1.4826 * np.median(np.abs(keep - gMed))
        bad |= np.abs(c - gMed) > max(minGlobal, nMad * gSig)
    bad = ndimage.binary_dilation(bad, iterations=grow)
    out = c.copy()
    out[bad] = np.nan
    return out


def vrtBand(path, name):
    '''Named band of a radar-geometry VRT and its (zeroDopplerTime, slantRange) pixel-centre axes.'''
    ds = gdal.Open(path)
    if ds is None:
        raise FileNotFoundError(path)
    gt = ds.GetGeoTransform()
    for i in range(1, ds.RasterCount + 1):
        b = ds.GetRasterBand(i)
        if b.GetMetadataItem('Description') == name or b.GetDescription() == name:
            a = b.ReadAsArray().astype(np.float64)
            nd = b.GetNoDataValue()
            a[a < -1e9] = np.nan
            if nd is not None and np.isfinite(nd):
                a[a == nd] = np.nan
            t = gt[3] + (np.arange(ds.RasterYSize) + 0.5) * gt[5]
            r = gt[0] + (np.arange(ds.RasterXSize) + 0.5) * gt[1]
            return a, t, r
    raise ValueError(f'No band {name} in {path}')


def denseOffsetWindow(h):
    '''(azimuth, range) dense-offset correlation window in SLC pixels, from the RUNW run configuration.'''
    rc = h[f'{RUNW}/metadata/processingInformation/parameters/runConfigurationContents'][()]
    rc = rc.decode() if isinstance(rc, bytes) else str(rc)
    block = rc[rc.find("'dense_offsets'"):]
    wa = re.search(r"'window_azimuth':\s*(\d+)", block)
    wr = re.search(r"'window_range':\s*(\d+)", block)
    return (int(wa.group(1)) if wa else 64), (int(wr.group(1)) if wr else 64)


def boxMean(z, sizeAz, sizeRg):
    '''NaN-aware box average.'''
    sizeAz, sizeRg = max(1, int(round(sizeAz))), max(1, int(round(sizeRg)))
    ok = np.isfinite(z)
    num = ndimage.uniform_filter(np.where(ok, z, 0.0), (sizeAz, sizeRg))
    den = ndimage.uniform_filter(ok.astype(float), (sizeAz, sizeRg))
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0.25, num / den, np.nan)


def linearWeights(src, dst):
    '''Indices and weights for linear interpolation from axis src to axis dst, clamped at the ends.'''
    x = np.interp(dst, src, np.arange(len(src)), left=0, right=len(src) - 1)
    i0 = np.clip(np.floor(x).astype(np.int64), 0, len(src) - 2)
    return i0, (x - i0).astype(np.float32)


def resample2d(z, t, r, tq, rq):
    '''Separable bilinear resampling of gap-free z(t, r) onto (tq, rq); values beyond the grid are clamped.'''
    ir, wr = linearWeights(r, rq)
    zr = z[:, ir] * (1 - wr) + z[:, ir + 1] * wr
    it, wt = linearWeights(t, tq)
    return zr[it] * (1 - wt)[:, None] + zr[it + 1] * wt[:, None]


def carrierCorrectionPhase(runwFile, velocityVrt, geometryVrt, verbose=True):
    '''
    Carrier correction (radians) on the RUNW interferogram grid, to be ADDED to the unwrapped phase before the
    ionosphere is estimated.

    Parameters
    ----------
    runwFile : str
        NISAR RUNW HDF5.
    velocityVrt, geometryVrt : str
        Simulated offsets with the velocity model and with zero motion (offsetSims/offsets.velocity.vrt and
        offsets.geom.vrt), SLC pixels, radar coordinates (slant range, zero-Doppler time).

    Returns
    -------
    correction : np.ndarray (float32), info : dict
    '''
    with h5py.File(runwFile, 'r') as h:
        sw = h[f'{RUNW}/swaths/frequencyA']
        tq, rq = sw['interferogram/zeroDopplerTime'][()], sw['interferogram/slantRange'][()]
        po = sw['pixelOffsets']
        c = po['HH/alongTrackOffset'][()].astype(np.float64)
        tPo, rPo = po['zeroDopplerTime'][()], po['slantRange'][()]
        m = h[f'{RUNW}/metadata']
        fdc = float(np.nanmedian(m['processingInformation/parameters/common/frequencyA/dopplerCentroid'][()]))
        vg = float(np.nanmedian(m['geolocationGrid/groundTrackVelocity'][()]))
        ref = m['processingInformation/parameters/reference/frequencyA']
        dtSlc, drSlc = float(ref['zeroDopplerTimeSpacing'][()]), float(ref['slantRangeSpacing'][()])
        winAz, winRg = denseOffsetWindow(h)
    kdc = 2 * np.pi * fdc / vg
    # applied rubbersheet shift (m), culled and gap-filled on its own (coarse) dense-offset grid
    c[c == 0] = np.nan                                       # exact zeros are fill
    crid = granuleCrid(runwFile)
    patchCull = crid is None or crid < PATCH_CULL_MAX_CRID
    c = nearestFill(cullOffsets(c, patchCull=patchCull))
    # model along-track motion (m), averaged over the dense-offset window, sampled at the dense-offset points
    vAz, tS, rS = vrtBand(velocityVrt, 'AzimuthOffsets')
    gAz, _, _ = vrtBand(geometryVrt, 'AzimuthOffsets')
    model = ((vAz - gAz) * dtSlc * vg).astype(np.float32)   # SLC azimuth pixels -> ground metres
    model = nearestFill(boxMean(model, winAz * dtSlc / abs(tS[1] - tS[0]), winRg * drSlc / abs(rS[1] - rS[0])))
    modelPo = resample2d(model, tS, rS, tPo, rPo)
    # ionospheric shift on the coarse grid; only the final correction goes to the interferogram grid
    shift = (c - modelPo).astype(np.float32)
    corr = (kdc * resample2d(shift, tPo, rPo, tq, rq)).astype(np.float32)
    info = dict(kdc=kdc, fdc=fdc, vg=vg, window=(winAz, winRg), crid=crid, patchCull=patchCull,
                sdIonoShift=float(np.std(corr) / kdc), sdCorrection=float(np.std(corr)))
    if verbose:
        print(f'Carrier correction: k_dc {kdc:.4f} rad/m (f_dc {fdc:.1f} Hz, V_gt {vg:.1f} m/s), '
              f'window {winAz}x{winRg} SLC px, CRID {crid} patch cull {patchCull}, sd(dA_c - S[dA_model]) {info["sdIonoShift"]:.2f} m, '
              f'sd(correction) {info["sdCorrection"]:.2f} rad')
    return corr, info
