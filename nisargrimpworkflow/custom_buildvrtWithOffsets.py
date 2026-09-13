import os
import glob
import json
import argparse
import numpy as np
import xml.etree.ElementTree as ET
from scipy.interpolate import RegularGridInterpolator

from osgeo import gdal

import utilities as u
from datetime import datetime

gdal.UseExceptions()

# Threshold (single-look range pixels) for flagging a pair of merged frames
# whose measured overlap difference disagrees with the range-gate-gap
# heuristic's prediction -- i.e. the metadata says "no inter-frame range bias"
# but the data says there is one (or vice versa). The empirical overlap
# correction is what the merge applies either way (see buildVrt's overlap
# equations); this threshold only controls when the disagreement is loud
# enough to record in a rangeDiscontinuityFound file for later QC.
RANGE_DISCONTINUITY_THRESH = 1.0


def loadBandToArray(path, band=1):
    """Open a GDAL file and return (float64 NaN-filled array, geotransform)."""
    ds = gdal.Open(path)
    if ds is None:
        u.myerror(f'Cannot open: {path}')
    gt = ds.GetGeoTransform()
    b = ds.GetRasterBand(band)
    data = b.ReadAsArray().astype(np.float64)
    nodata = b.GetNoDataValue()
    ds = None
    data[data < -1e9] = np.nan
    if nodata is not None:
        data[data == nodata] = np.nan
    return data, gt


def resampleToGrid(srcData, srcGt, dstGt, dstShape, method='linear'):
    """Resample srcData (on srcGt) onto the grid defined by (dstGt, dstShape)."""
    dstNr, dstNc = dstShape
    srcNr, srcNc = srcData.shape

    dstX = dstGt[0] + (np.arange(dstNc) + 0.5) * dstGt[1]
    dstY = dstGt[3] + (np.arange(dstNr) + 0.5) * dstGt[5]

    srcCol = (dstX - srcGt[0]) / srcGt[1] - 0.5
    srcRow = (dstY - srcGt[3]) / srcGt[5] - 0.5

    interp = RegularGridInterpolator(
        (np.arange(srcNr), np.arange(srcNc)),
        srcData.astype(np.float64),
        method=method, bounds_error=False, fill_value=np.nan
    )

    rowGrid = np.tile(srcRow[:, np.newaxis], (1, dstNc))
    colGrid = np.tile(srcCol[np.newaxis, :], (dstNr, 1))
    return interp((rowGrid, colGrid))


def getOverlapStats(file1, file2, yOffset1=0.0, yOffset2=0.0,
                    maskData=None, maskGt=None, verbose=False, band=1):
    """Median difference (f1 - f2) in georeferenced overlap.

    yOffset1/yOffset2 are midnight-wrap corrections (seconds) applied to Y origins.
    band selects which band to compare -- different bands of the same file (e.g.
    RangeOffsets vs AzimuthOffsets in offsets.geom.vrt) are different physical
    quantities and can have entirely independent overlap discrepancies, so the
    caller must compute/apply a bias per band rather than assuming band 1 speaks
    for the whole file.
    Returns dict {'median': float, 'n': int} or None if overlap is insufficient.
    """
    ds1 = gdal.Open(file1)
    ds2 = gdal.Open(file2)

    gt1 = ds1.GetGeoTransform()
    gt2 = ds2.GetGeoTransform()

    nodata1 = ds1.GetRasterBand(band).GetNoDataValue()
    nodata2 = ds2.GetRasterBand(band).GetNoDataValue()

    y0_1 = gt1[3] + yOffset1
    y0_2 = gt2[3] + yOffset2

    b1 = [gt1[0], gt1[0] + ds1.RasterXSize * gt1[1],
          min(y0_1, y0_1 + ds1.RasterYSize * gt1[5]),
          max(y0_1, y0_1 + ds1.RasterYSize * gt1[5])]
    b2 = [gt2[0], gt2[0] + ds2.RasterXSize * gt2[1],
          min(y0_2, y0_2 + ds2.RasterYSize * gt2[5]),
          max(y0_2, y0_2 + ds2.RasterYSize * gt2[5])]

    inter = [max(b1[0], b2[0]), min(b1[1], b2[1]),
             max(b1[2], b2[2]), min(b1[3], b2[3])]

    f1 = os.path.basename(file1)
    f2 = os.path.basename(file2)
    if inter[0] >= inter[1] or inter[2] >= inter[3]:
        print(f"  overlap({f1},{f2}) band {band}: no spatial overlap — b1={b1} b2={b2}")
        return None

    def toPx(gt, yOffset, x, y):
        return (int(round((x - gt[0]) / gt[1])),
                int(round((y - (gt[3] + yOffset)) / gt[5])))

    x1, y1 = toPx(gt1, yOffset1, inter[0], inter[3] if gt1[5] < 0 else inter[2])
    x2, y2 = toPx(gt2, yOffset2, inter[0], inter[3] if gt2[5] < 0 else inter[2])

    pw = int(round((inter[1] - inter[0]) / gt1[1]))
    ph = int(round((inter[3] - inter[2]) / abs(gt1[5])))

    if verbose:
        print(f"  overlap({f1},{f2}) band {band}: inter y=[{inter[2]:.2f},{inter[3]:.2f}] "
              f"px1=({x1},{y1}) px2=({x2},{y2}) pw={pw} ph={ph}")

    try:
        arr1 = ds1.GetRasterBand(band).ReadAsArray(x1, y1, pw, ph).astype(np.float32)
        arr2 = ds2.GetRasterBand(band).ReadAsArray(x2, y2, pw, ph).astype(np.float32)

        valid = ((arr1 != nodata1) & (arr2 != nodata2)
                 & np.isfinite(arr1) & np.isfinite(arr2))

        if maskData is not None and maskGt is not None:
            # Native (non-midnight-corrected) gt of the overlap window in file1 coordinates
            overlapGt = [gt1[0] + x1 * gt1[1], gt1[1], 0,
                         gt1[3] + y1 * gt1[5], 0, gt1[5]]
            maskOnOverlap = resampleToGrid(maskData, maskGt, overlapGt,
                                           (ph, pw), method='nearest') > 0.5
            valid = valid & maskOnOverlap

        nValid = int(np.sum(valid))
        if verbose:
            print(f"  overlap({f1},{f2}) band {band}: {nValid} valid pixels")
        if nValid < 100:
            return None

        diffs = arr1[valid] - arr2[valid]
        return {'median': float(np.median(diffs)), 'n': nValid}

    except Exception as e:
        print(f"  overlap({f1},{f2}) band {band}: exception — {e}")
        return None


# A seam needs at least this many valid overlapping pixels before the merge can
# trust the median difference it ties two frames with. Healthy Antarctic seams
# carry 150,000-260,000; getOverlapStats() already refuses under 100. The margin
# between the two is deliberate -- a sliver of coincidental overlap gives a
# median dominated by whatever few pixels happen to be valid there.
MIN_SEAM_OVERLAP_PIXELS = 500


def frameYOffsets(files):
    """Cumulative midnight-wrap correction per file, matching buildVrt()'s own
    pass so seam tests run in the same Y convention as the merge."""
    MIDNIGHT_WRAP = 86400
    yOffsets = []
    cumulativeYOffset = 0
    prevY0 = None
    for f in files:
        ds = gdal.Open(f)
        if ds is None:
            u.myerror(f'Cannot open: {f}')
        y0 = ds.GetGeoTransform()[3]
        if prevY0 is not None and (y0 + cumulativeYOffset) - prevY0 < -MIDNIGHT_WRAP / 2:
            cumulativeYOffset += MIDNIGHT_WRAP
        y0 += cumulativeYOffset
        prevY0 = y0
        yOffsets.append(cumulativeYOffset)
    return yOffsets


# The simulated-geometry raster whose seam step matches each offset band's.
# NISAR products have the geometry removed and this workflow adds a simulation back
# (optionally on a different DEM); the step at a seam lives ENTIRELY in that
# simulation, not in the measured offsets. Measured over 134 seams in 30 virtual
# frames: median |measured jump| 0.035 px, median |measured - simulated| 0.0012 px
# (max 0.0914). Every simulated jump was an exact integer multiple of the PRF --
# a whole second of azimuth time -- with fractional cases (track-80/5542_0000 seam
# 138|139, sim -466.9512) equally well explained.
#
# Tying on the simulation instead of the measured offsets is better on three counts:
# it is analytic so the median is exact rather than data-dependent; it is defined
# across the whole frame footprint, so seams whose measured overlap is all invalid
# pixels (ice, ocean, decorrelation) still tie -- 15 of 44 otherwise-unbridgeable
# seams in a 40-frame sample, at ~258,000 overlap pixels each; and it removes exactly
# the quantity that was injected. Only a true time gap, where the frames do not
# overlap at all, remains untieable.
#
# Range validated the same way 2026-09-10, over 104 seams in 10 virtual frames:
# of the 15 seams jumping more than 1 px (up to 1648 px), median |measured - simulated|
# was 0.0005 px and the worst 0.0029 px. Range jumps are whole range-gate steps
# (1648, 1248, 1168, 1040, 992, 928, 880, 848, 808, 520, 440, 376, 136, 112, 8),
# the range counterpart of azimuth's whole seconds.
SIM_FOR_BAND = {'AzimuthOffsets': 'offsets.geom.da.tif',
                'RangeOffsets': 'offsets.geom.dr.tif'}


def simTieFiles(files, bandDescription):
    """The per-frame simulated-geometry rasters to tie `files` with, or None.

    Returns None unless every frame has the simulation for this band, so a merge
    either ties entirely on simulation or entirely on measured data -- never a mix
    across frames of the same band, which would compare two different quantities.
    """
    simName = SIM_FOR_BAND.get(bandDescription)
    if simName is None:
        return None
    sims = [os.path.join(os.path.dirname(os.path.abspath(f)), 'offsetSims', simName)
            for f in files]
    return sims if all(os.path.exists(s) for s in sims) else None



# Rows read from each side of a true time gap to carry the simulation across it.
# The simulation varies by under 0.15 px along a whole frame (track-167/4072:
# frame medians 565.198 vs 565.034 px, 150 s apart), so its value at the frame
# edge is its value at the seam to ~0.1 px, and the step the merge must remove
# is simply simA(edge) - simB(edge). This is the same quantity getOverlapStats()
# measures where the footprints overlap, evaluated where they do not.
SIM_GAP_EDGE_ROWS = 50


def simGapTie(sim1, sim2, yOffset1=0.0, yOffset2=0.0, edgeRows=SIM_GAP_EDGE_ROWS):
    """Median (sim1 - sim2) across a seam the two footprints do not span.

    Reads the trailing rows of the earlier frame and the leading rows of the
    later one and differences their medians. Returns {'median', 'n'} or None
    when the frames actually overlap (use getOverlapStats) or either edge has
    no valid pixels.

    Why this is legitimate: the seam step lives entirely in the re-added
    simulated geometry, and that geometry is known on both sides of a gap even
    though the measured offsets are not. Before this, a gap seam was pinned
    from the primary-minus-secondary NominalTime slip, which assumes every
    frame carries its simulation. Where a frame's simulation never ran
    (all-zero offsets.geom.da.tif, 157 frames project-wide) that frame sits in
    the RAW domain, the slip predicts a step that is not there, and the pin
    disagrees with the sim ties by the missing constant -- track-167/4072_0010
    went from 2.6 m to 477 m az.est sigma that way. Tying the gap from the
    simulation itself is consistent with the sim ties by construction, so no
    pin and no split is needed.
    """
    ds1 = gdal.Open(sim1)
    ds2 = gdal.Open(sim2)
    if ds1 is None or ds2 is None:
        return None
    gt1 = ds1.GetGeoTransform()
    gt2 = ds2.GetGeoTransform()
    y0_1 = gt1[3] + yOffset1
    y0_2 = gt2[3] + yOffset2
    end1 = y0_1 + ds1.RasterYSize * gt1[5]
    end2 = y0_2 + ds2.RasterYSize * gt2[5]
    lo1, hi1 = min(y0_1, end1), max(y0_1, end1)
    lo2, hi2 = min(y0_2, end2), max(y0_2, end2)
    if lo2 < hi1 and lo1 < hi2:
        return None  # footprints overlap in time; not a gap
    firstIsEarlier = hi1 <= lo2
    # Row index increases with time when gt[5] > 0.
    ascending1 = gt1[5] > 0
    ascending2 = gt2[5] > 0

    def edge(ds, rows, trailing, ascending):
        n = min(rows, ds.RasterYSize)
        # trailing in TIME is the last rows when rows ascend in time, else the first
        useLast = trailing == ascending
        y = ds.RasterYSize - n if useLast else 0
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray(0, y, ds.RasterXSize, n).astype(np.float32)
        nodata = band.GetNoDataValue()
        valid = np.isfinite(arr)
        if nodata is not None:
            valid &= arr != nodata
        return arr[valid]

    e1 = edge(ds1, edgeRows, trailing=firstIsEarlier, ascending=ascending1)
    e2 = edge(ds2, edgeRows, trailing=not firstIsEarlier, ascending=ascending2)
    if e1.size < 100 or e2.size < 100:
        return None
    return {'median': float(np.median(e1) - np.median(e2)),
            'n': int(min(e1.size, e2.size))}


def unbridgeableSeams(files, band=1, minPixels=MIN_SEAM_OVERLAP_PIXELS,
                      maskData=None, maskGt=None, verbose=False):
    """Indices i where files[i] and files[i+1] cannot be tied by measured overlap.

    The merge ties adjacent frames with the median difference of their values in
    the georeferenced overlap. Where that measurement does not exist -- the
    frames do not overlap in time, or they do but every overlapping pixel is
    invalid -- the least-squares solve has to guess the step from geodat dead
    reckoning, and when the guess is wrong the merged product carries a bias
    discontinuity in the middle of a virtual frame. az.est then fits one
    constant plus one ramp across the staircase and reports a residual that is
    not scatter at all.

    Callers use this to split a group at the seam instead, so each side becomes
    its own internally consistent virtual frame. Metadata alone cannot answer
    this: of 406 Antarctic seams the merge could not bridge, only 223 were
    predictable from frame times, the other 183 overlapped geometrically but had
    no valid pixels to measure.
    """
    files = list(files)
    if len(files) < 2:
        return []
    # Ask the question the merge will actually ask: it ties on the simulated
    # geometry where that exists (SIM_FOR_BAND), so a seam whose measured overlap
    # is all invalid pixels is NOT unbridgeable -- the simulation still spans it.
    # Testing the measured data here would split groups the merge can join.
    simFiles = simTieFiles(files, 'AzimuthOffsets')
    yOffsets = frameYOffsets(files)
    simYOffsets = frameYOffsets(simFiles) if simFiles else None
    seams = []
    for i in range(len(files) - 1):
        result = None
        if simFiles:
            result = getOverlapStats(simFiles[i], simFiles[i + 1],
                                     yOffset1=simYOffsets[i], yOffset2=simYOffsets[i + 1],
                                     verbose=verbose, band=1)
        if result is None:
            result = getOverlapStats(files[i], files[i + 1],
                                     yOffset1=yOffsets[i], yOffset2=yOffsets[i + 1],
                                     maskData=maskData, maskGt=maskGt,
                                     verbose=verbose, band=band)
        if (result is None or result['n'] < minPixels) and simFiles:
            # A true time gap the simulation can still carry -- see simGapTie().
            result = simGapTie(simFiles[i], simFiles[i + 1],
                               yOffset1=simYOffsets[i], yOffset2=simYOffsets[i + 1])
        if result is None or result['n'] < minPixels:
            seams.append(i)
    return seams


def getReferenceStats(filepath, refData, refGt, maskData=None, maskGt=None, band=1):
    """Median/std/n of (file - reference) over valid, unmasked pixels.

    band selects which band to compare -- see getOverlapStats() docstring for why
    this matters for multi-band files where different bands are independent
    physical quantities.

    Returns dict {'median': float, 'std': float, 'n': int} or None.
    """
    ds = gdal.Open(filepath)
    fileGt = ds.GetGeoTransform()
    fileNr, fileNc = ds.RasterYSize, ds.RasterXSize
    fileData = ds.GetRasterBand(band).ReadAsArray().astype(np.float64)
    nodata = ds.GetRasterBand(band).GetNoDataValue()
    ds = None

    fileData[fileData < -1e9] = np.nan
    if nodata is not None:
        fileData[fileData == nodata] = np.nan

    refOnFile = resampleToGrid(refData, refGt, fileGt, (fileNr, fileNc))

    valid = np.isfinite(fileData) & np.isfinite(refOnFile)
    if maskData is not None and maskGt is not None:
        maskOnFile = resampleToGrid(maskData, maskGt, fileGt,
                                    (fileNr, fileNc), method='nearest') > 0.5
        valid = valid & maskOnFile

    nValid = int(valid.sum())
    if nValid < 100:
        return None

    diffs = fileData[valid] - refOnFile[valid]
    return {'median': float(np.median(diffs)),
            'std': float(np.std(diffs)),
            'n': nValid}


def _frameLabel(path):
    """Label a per-frame input by its frame directory (.../4248_47/range.offsets.vrt
    and .../4248_47/offsetSims/offsets.geom.vrt both -> '4248_47').

    Every frame's input has the same basename, so summary rows built from the
    basename alone all read 'range.offsets.vrt -> range.offsets.vrt' and cannot
    be told apart -- useless exactly when you are trying to work out which frame
    carried a bias."""
    parts = [p for p in os.path.dirname(os.path.abspath(path)).split(os.sep) if p]
    for part in reversed(parts):
        if '_' in part:
            return part
    return parts[-1] if parts else path


def writeSummary(outputVrt, files, offsetsDict, overlapResultsByBand, refStatsByBand, bandCount):
    """Write bias-correction and fit statistics to a text file beside the VRT.

    offsetsDict is keyed by (file, band); overlapResultsByBand/refStatsByBand are
    keyed by band. Each band gets its own section -- they're independent physical
    quantities (e.g. RangeOffsets vs AzimuthOffsets) with independent biases, see
    getOverlapStats().
    """
    statsPath = outputVrt + '.stats'
    lines = ['# custom_buildvrtWithOffsets summary',
             f'# output: {outputVrt}', '']

    for band in range(1, bandCount + 1):
        overlapResults = overlapResultsByBand.get(band, {})
        refStats = refStatsByBand.get(band, [None] * len(files))

        if overlapResults:
            lines += [f'## Band {band} overlap residuals (after bias correction)',
                      '#  i  file_i -> file_{i+1}  raw_median  corrected_residual  n_pixels'
                      '  tie_source']
            for i, result in overlapResults.items():
                f1 = _frameLabel(files[i])
                f2 = _frameLabel(files[i + 1])
                if result is None:
                    lines.append(f'  {i}  {f1} -> {f2}  NO_OVERLAP')
                else:
                    rawMed = result['median']
                    o1 = offsetsDict.get((files[i], band), 0.0)
                    o2 = offsetsDict.get((files[i + 1], band), 0.0)
                    residual = rawMed - (o2 - o1)
                    lines.append(f'  {i}  {f1} -> {f2}  '
                                  f'{rawMed:+.4f}  {residual:+.4f}  {result["n"]}'
                                  f'  {result.get("source", "measured")}')
            lines.append('')

            # The residuals above only show *relative* continuity; they read
            # +0.0000 whenever the merge is self-consistent, no matter which
            # frame's convention the whole product was anchored to. Record the
            # absolute bias actually applied per frame as well, so an
            # anchored-on-the-wrong-frame merge is visible here instead of
            # needing the solve to be instrumented to find it.
            lines += [f'## Band {band} applied bias (per frame)',
                      '#  frame  applied_bias']
            for f in files:
                lines.append(f'  {_frameLabel(f)}  '
                             f'{offsetsDict.get((f, band), 0.0):+.4f}')
            lines.append('')

        if any(s is not None for s in refStats):
            lines += [f'## Band {band} reference phase fit (per frame)',
                      '#  file  ref_median  std  n_pixels  applied_offset  corrected_residual']
            for f, s in zip(files, refStats):
                fname = os.path.basename(f)
                applied = offsetsDict.get((f, band), 0.0)
                if s is None:
                    lines.append(f'  {fname}  NO_DATA  applied_offset={applied:+.4f}')
                else:
                    residual = s['median'] + applied
                    lines.append(f'  {fname}  {s["median"]:+.4f}  {s["std"]:.4f}  '
                                  f'{s["n"]}  {applied:+.4f}  {residual:+.4f}')
            lines.append('')

    with open(statsPath, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"Written: {statsPath}")


def _updateRangeDiscontinuityFile(outputVrt, sectionLines):
    """Maintain a rangeDiscontinuityFound diagnostic file beside outputVrt.

    The file holds one section per merged product ('## product: <basename>'
    header lines), since several products (offsets.range-azimuth.vrt,
    range.offsets.vrt, ...) merge into the same directory and each is checked
    independently. This rewrites only the current product's section: existing
    sections from other products are preserved, a stale section for this
    product is dropped, and the file itself is removed when no sections remain
    (i.e. a re-merge that no longer disagrees clears its own old flag)."""
    flagPath = os.path.join(os.path.dirname(os.path.abspath(outputVrt)),
                            'rangeDiscontinuityFound')
    product = os.path.basename(outputVrt)
    kept = []
    if os.path.exists(flagPath):
        keep = True
        with open(flagPath) as fh:
            for line in fh:
                if line.startswith('## product: '):
                    keep = line.strip() != f'## product: {product}'
                if keep:
                    kept.append(line.rstrip('\n'))
    while kept and not kept[-1].strip():
        kept.pop()
    if sectionLines:
        if kept:
            kept.append('')
        kept.append(f'## product: {product}')
        kept.extend(sectionLines)
    if kept:
        with open(flagPath, 'w') as fh:
            fh.write('\n'.join(kept) + '\n')
        if sectionLines:
            u.mywarning(f'Range discontinuity flagged -- see {flagPath}')
    elif os.path.exists(flagPath):
        os.remove(flagPath)


def bakeVrtInPlace(vrtPath, skipIfExists=False):
    """Bake vrtPath (a possibly recursive multi-source VRT) into a flat GeoTIFF,
       then replace vrtPath with a minimal single-source VRT wrapper around that
       GeoTIFF -- same path/filename, same metadata/band descriptions, one direct
       file open for downstream readers instead of a recursive chain.

       If skipIfExists, does nothing when the corresponding .tif already exists
       (i.e. vrtPath was already baked) -- used by --bakeOnly to cheaply retrofit
       virtual frames built before this baking step existed, without redoing the
       (much more expensive) sub-frame mosaicking buildVrt() does."""
    tifPath = os.path.splitext(vrtPath)[0] + '.tif'
    if skipIfExists and os.path.exists(tifPath):
        print(f"Skipping {vrtPath} -- {tifPath} already exists")
        return
    # Guard against self-clobber: baking an ALREADY-baked wrapper (whose only
    # source is tifPath itself) makes gdal.Translate truncate tifPath before
    # reading it back through the wrapper -- destroying the data and writing an
    # all-NaN tif. Detect that case and no-op instead.
    if os.path.exists(vrtPath) and os.path.exists(tifPath):
        try:
            srcList = gdal.Open(vrtPath).GetFileList() or []
        except Exception:
            srcList = []
        if any(os.path.abspath(s) == os.path.abspath(tifPath) for s in srcList):
            print(f"Skipping {vrtPath} -- already baked (wrapper sources {tifPath})")
            return
    gdal.Translate(tifPath, vrtPath, format='GTiff')
    gdal.Translate(vrtPath, tifPath, format='VRT')
    print(f"Baked to {os.path.basename(tifPath)}, rewrapped {os.path.basename(vrtPath)}")


def _findGeodat(directory, secondary=False):
    """Find the primary (or, if secondary=True, secondary) geodat*.geojson in
    directory, or None."""
    if secondary:
        candidates = sorted(glob.glob(os.path.join(directory, 'geodat*.secondary.geojson')))
    else:
        candidates = sorted(f for f in glob.glob(os.path.join(directory, 'geodat*.geojson'))
                            if not f.endswith('.secondary.geojson'))
    return candidates[0] if candidates else None


def _readGeodatProperty(geodatPath, key):
    with open(geodatPath) as fh:
        return json.load(fh).get('properties', {}).get(key)


def _rangeGateGap(directory):
    """primary MLNearRange - secondary MLNearRange for the geodat pair in
    directory, or None if either geodat/property is missing. This -- not
    primary MLNearRange alone -- is what actually determines whether a
    frame's RangeOffsets values need a bias correction relative to its
    neighbors: RangeOffsets measures a *relative* shift between primary and
    secondary, so a frame whose primary and secondary near ranges moved by
    the same amount (gap unchanged) needs no correction at all, even if its
    own absolute near range differs from its neighbors'. Only a frame whose
    gap itself differs carries a real value-domain bias."""
    primaryFile = _findGeodat(directory)
    secondaryFile = _findGeodat(directory, secondary=True)
    if primaryFile is None or secondaryFile is None:
        return None
    primary = _readGeodatProperty(primaryFile, 'MLNearRange')
    secondary = _readGeodatProperty(secondaryFile, 'MLNearRange')
    if primary is None or secondary is None:
        return None
    return primary - secondary


def _nominalSeconds(geodatPath):
    """NominalTime ('H:MM:SS.ssssss') as seconds since midnight, or None."""
    nominal = _readGeodatProperty(geodatPath, 'NominalTime') if geodatPath else None
    if nominal is None:
        return None
    fields = str(nominal).split(':')
    if len(fields) != 3:
        return None
    try:
        return int(fields[0]) * 3600 + int(fields[1]) * 60 + float(fields[2])
    except ValueError:
        return None


def _azimuthTimeGap(directory):
    """primary NominalTime - secondary NominalTime, in seconds, for the geodat
    pair in directory, or None if either geodat/property is missing. Azimuth
    counterpart of _rangeGateGap: a frame whose two acquisitions start at the
    same nominal second needs no correction relative to one that does, while an
    integer-second slip between them is exactly PRF single-look azimuth pixels
    of value-domain bias. Verified to predict the measured azimuth overlap step
    on 2139 of 2139 seams across the Greenland project."""
    primary = _nominalSeconds(_findGeodat(directory))
    secondary = _nominalSeconds(_findGeodat(directory, secondary=True))
    if primary is None or secondary is None:
        return None
    gap = primary - secondary
    # Guard the midnight boundary, where the two times can straddle 00:00:00.
    if gap > 43200:
        gap -= 86400
    elif gap < -43200:
        gap += 86400
    return gap


def buildVrt(outputVrt, inputPattern, overWrite, offsets,
             referenceVrt=None, maskVrt=None, verbose=False):

    files = []
    if isinstance(inputPattern, list):
        for p in inputPattern:
            files.extend(glob.glob(p))
    else:
        files = glob.glob(inputPattern)
    files = sorted(list(set(files)))

    if not files:
        print("Error: No files found.")
        return

    # Load reference and mask once (both may be None)
    refData, refGt = loadBandToArray(referenceVrt) if referenceVrt else (None, None)
    maskData, maskGt = loadBandToArray(maskVrt) if maskVrt else (None, None)

    # 1. COMPUTE MIDNIGHT-CORRECTED Y ORIGINS
    firstDs = gdal.Open(files[0])
    gt = firstDs.GetGeoTransform()
    proj = firstDs.GetProjection()
    bandCount = firstDs.RasterCount
    globalMetadata = firstDs.GetMetadata()
    dtype = firstDs.GetRasterBand(1).DataType
    resX = gt[1]
    resY = gt[5]
    # Fall back to the band's "Description" metadata item when GetDescription()
    # is empty: the simoffsets products (offsetSims/offsets.geom.vrt,
    # offsets.velocity.vrt) carry their band name only in metadata, so keying
    # the RangeOffsets dead reckoning off GetDescription() alone silently
    # skipped them -- leaving the simulated offsets in a different convention
    # than the measured range.offsets they are differenced against. The source
    # loop below already applies this same fallback when writing descriptions.
    bandDescriptions = [firstDs.GetRasterBand(b).GetDescription()
                        or firstDs.GetRasterBand(b).GetMetadataItem('Description') or ''
                        for b in range(1, bandCount + 1)]

    # bandDescriptions comes from files[0] alone but keys which band gets the
    # RangeOffsets/AzimuthOffsets dead reckoning, so a frame whose bands are
    # ordered differently would have a range bias applied to an azimuth band --
    # silently, and catastrophically. Verified consistent across every frame of
    # several segments, so this is a guard rather than a fix for a live bug.
    for f in files[1:]:
        ds = gdal.Open(f)
        other = [ds.GetRasterBand(b).GetDescription()
                 or ds.GetRasterBand(b).GetMetadataItem('Description') or ''
                 for b in range(1, min(ds.RasterCount, bandCount) + 1)]
        if ds.RasterCount != bandCount or other != bandDescriptions:
            u.mywarning(f'{os.path.basename(outputVrt)}: band layout of '
                        f'{_frameLabel(f)} ({ds.RasterCount} bands, {other}) differs from '
                        f'{_frameLabel(files[0])} ({bandCount} bands, {bandDescriptions});'
                        f'\n\tper-band bias solving assumes they match')
        ds = None

    # Every source must sit on the same grid: sources are placed by integer
    # pixel offset ((x0 - minX) / resX), so a frame with a different pixel size
    # is not just misplaced, it redefines the mosaic extent. A frame whose
    # upstream simoffsets run failed is left ungeoreferenced (identity
    # geotransform, pixel size 1.0, origin ~0) and would silently stretch the
    # mosaic from the projection origin out to the real data --
    # track-145/4050_39 did exactly that, turning a 1659x19213 raster into
    # 20773x2240205 and writing 694 GB per band before anyone noticed. Fail
    # loudly instead; the caller already handles a missing output via .fail.
    offGrid = []
    for f in files:
        fGt = gdal.Open(f).GetGeoTransform()
        if (abs(fGt[1] - resX) > abs(resX) * 1e-6
                or abs(fGt[5] - resY) > abs(resY) * 1e-6):
            offGrid.append(f'{_frameLabel(f)} (pixel {fGt[1]:g}x{fGt[5]:g}, '
                           f'origin {fGt[0]:.2f},{fGt[3]:.2f})')
    if offGrid:
        u.myerror(f'{os.path.basename(outputVrt)}: {len(offGrid)} of {len(files)} inputs are '
                  f'not on the {resX:g}x{resY:g} grid of {_frameLabel(files[0])} -- refusing to '
                  f'mosaic.\n\tCheck for a failed upstream run (fail.simoffsets.* marker) in:'
                  f'\n\t  ' + '\n\t  '.join(offGrid))

    datasets = []
    allX = []
    allY = []
    yOffsets = []

    MIDNIGHT_WRAP = 86400
    cumulativeYOffset = 0
    prevY0 = None

    for f in files:
        ds = gdal.Open(f)
        lGt = ds.GetGeoTransform()
        xSize = ds.RasterXSize
        ySize = ds.RasterYSize

        x0 = lGt[0]
        y0 = lGt[3]

        if prevY0 is not None and (y0 + cumulativeYOffset) - prevY0 < -MIDNIGHT_WRAP / 2:
            cumulativeYOffset += MIDNIGHT_WRAP
            print(f"Midnight wrap detected at {f}: adding {cumulativeYOffset} to Y origin")
        y0 += cumulativeYOffset
        prevY0 = y0
        yOffsets.append(cumulativeYOffset)

        x1 = x0 + xSize * resX
        y1 = y0 + ySize * resY

        allX.extend([x0, x1])
        allY.extend([y0, y1])
        datasets.append((f, x0, y0, xSize, ySize))

    # 2. COLLECT REFERENCE STATS AND OVERLAP STATS, INDEPENDENTLY PER BAND.
    # Different bands of the same file (e.g. RangeOffsets vs AzimuthOffsets in
    # offsets.geom.vrt, or vx/vy/vz in offsets.velocity.vrt) are different physical
    # quantities that can have entirely independent overlap discrepancies -- a bias
    # computed from band 1 must not be assumed to apply to band 2, so every band gets
    # its own reference/overlap stats and its own least-squares solve below.
    n = len(files)
    # Anchor the no-reference fallback to whichever file is geometrically leftmost
    # (i.e. the file that sets minX/the merged geotransform's origin below), not
    # file 0 in list order -- otherwise the value-bias solve and the geodat's
    # MLNearRange (derived from minX) can be anchored to two different frames,
    # producing a constant mismatch between the two equal to that frame's range offset.
    anchorIdx = min(range(n), key=lambda i: datasets[i][1])
    # RangeOffsets measures a *relative* shift between primary and secondary,
    # so what actually determines whether a frame needs a bias correction is
    # whether its own (primary MLNearRange - secondary MLNearRange) *gap*
    # differs from its neighbors' -- not whether its primary MLNearRange
    # alone differs. A frame whose primary and secondary near ranges both
    # shifted by the same amount (gap unchanged) needs no correction at all,
    # even though its absolute near range looks anomalous; only a frame whose
    # gap itself differs carries a real value-domain bias. (Confirmed with
    # user 2026-07-13 against two real cases: track-41/3254_133 has a fixed
    # secondary near range across all 7 frames while only its own primary
    # differs, so its gap genuinely differs -- a real ~1774m correction is
    # needed. track-5/3564_133 has primary and secondary shift together
    # (both exactly track each other, frame-by-frame) -- gap is unchanged,
    # so no near-range-driven correction should be applied there at all; an
    # earlier primary-only version of this fix wrongly computed a large
    # correction for that case.) This also fixes frames with no (or
    # insufficient) overlap with their neighbor, where the overlap solve has
    # nothing to anchor them and previously left their bias at whatever the
    # least-squares fallback produced (effectively unconstrained).
    #
    # Reference this to whichever input frame's own per-physical-frame geodat
    # matches the *master* merged geodat's MLNearRange (the value downstream
    # tools like rparams actually read as RNear) -- not just anchorIdx (whichever
    # RangeOffsets input happens to be geometrically leftmost among the files
    # being merged here). In practice these usually name the same frame, but
    # matching against the master geodat directly is what actually makes this
    # bias consistent with RNear, rather than relying on that coincidence.
    #
    # Units: MLNearRange differences are in meters, but the RangeOffsets band
    # itself is stored in *single-look* SLC range pixels -- NOT the multilook
    # grid spacing (resX) used for the merged VRT's geotransform. For this
    # product resX happens to be exactly 15x SLCRangePixelSize (the
    # offset-estimation grid step vs. the SLC pixel the offset is measured
    # in) -- dividing by resX instead of SLCRangePixelSize understated the
    # correction by that same factor.
    masterGeodatPath = _findGeodat(os.path.dirname(os.path.abspath(outputVrt)))
    slcRangePixelSize = (_readGeodatProperty(masterGeodatPath, 'SLCRangePixelSize')
                         if masterGeodatPath else None)
    if slcRangePixelSize is None:
        # SLCRangePixelSize is a fixed sensor property, not something only the
        # master geodat can supply -- fall back to any input frame's own
        # geodat rather than silently using resX (the multilook grid spacing,
        # wrong by ~15x for this product) if the master happens to be missing.
        for f in files:
            g = _findGeodat(os.path.dirname(os.path.abspath(f)))
            slcRangePixelSize = _readGeodatProperty(g, 'SLCRangePixelSize') if g else None
            if slcRangePixelSize is not None:
                break
    masterNearRange = (_readGeodatProperty(masterGeodatPath, 'MLNearRange')
                       if masterGeodatPath else None)
    referenceIdx = anchorIdx
    if masterNearRange is not None:
        ownNearRanges = []
        for f in files:
            g = _findGeodat(os.path.dirname(os.path.abspath(f)))
            ownNearRanges.append(_readGeodatProperty(g, 'MLNearRange') if g else None)
        validOwn = [i for i, v in enumerate(ownNearRanges) if v is not None]
        if validOwn:
            referenceIdx = min(validOwn, key=lambda i: abs(ownNearRanges[i] - masterNearRange))
    pixelSize = slcRangePixelSize if slcRangePixelSize else resX
    gaps = [_rangeGateGap(os.path.dirname(os.path.abspath(f))) for f in files]
    # Reference the dead reckoning to the frame with the SMALLEST |range-gate
    # gap|, not to any particular frame's position in the list. A near-zero gap
    # means the primary and secondary range gates line up, i.e. a normally
    # aligned pair; a gap of hundreds/thousands of metres is a misalignment
    # artifact in that one frame. When a merge mixes the two, the aligned
    # frames define the correct zero and the misaligned ones get corrected onto
    # it -- never the other way round.
    #
    # Matching referenceIdx on primary MLNearRange alone did not achieve this:
    # when every frame shares the same primary near range (the normal case) the
    # match is a tie that lands on file 0, which may be exactly the frame whose
    # *secondary* range gate is anomalous.
    rangeAnchorIdx = referenceIdx
    if all(g is not None for g in gaps):
        # The MERGED GEODAT's own gap is the reference. That is the range
        # geometry simoffsets/rparams simulate against, so anchoring to it is
        # exactly what drives rparams' fitted Bp to ~0 -- confirmed both ways:
        # Greenland segments anchored elsewhere fitted Bp == the gap difference
        # (3575 m, 2724 m, ...) and collapsed to ~0 once re-anchored here, while
        # Antarctica segments already anchored here fit Bp ~ 0 today.
        #
        # Do NOT substitute "the frame with the smallest |gap|". That coincides
        # with the master gap throughout Greenland but is NOT the invariant: in
        # Antarctica a large gap is the norm, and the smallest-|gap| frame is
        # often the odd one out. track-145/4742_0001 has six frames at 1249.1 m
        # and one at 0.0 (master 1249.1); track-10/4780_0000 has four at 2848.0
        # and two at 1399.0 (master 2848.0). Anchoring on the smallest there
        # would drag every good frame onto the anomalous one -- this bug in
        # reverse -- and inject a ~1250-1450 m spurious Bp into products that
        # are correct today.
        masterGap = _rangeGateGap(os.path.dirname(os.path.abspath(outputVrt)))
        # Only distrust the merged geodat if its gap corresponds to no frame at
        # all: its primary and secondary MLNearRange are independent min()s
        # across frames, so in principle they can come from different frames and
        # describe a pairing that never existed. Not observed in any project;
        # this is the guard for it.
        if masterGap is not None and any(abs(masterGap - g) / pixelSize
                                         <= RANGE_DISCONTINUITY_THRESH for g in gaps):
            referenceGap = masterGap
        else:
            referenceGap = min(gaps, key=abs)
            if masterGap is not None:
                u.mywarning(f'{os.path.basename(outputVrt)}: merged geodat range-gate gap '
                            f'{masterGap:.3f} m matches no frame '
                            f'(frames have {sorted({round(g, 1) for g in gaps})});'
                            f'\n\tfalling back to the smallest-|gap| frame '
                            f'({referenceGap:.3f} m) to anchor RangeOffsets')
        deadReckonBias = [(referenceGap - g) / pixelSize for g in gaps]
        rangeAnchorIdx = min(range(n), key=lambda i: abs(deadReckonBias[i]))
        # Keep the pre-existing anchor when it is already in the reference
        # convention, so merges with nothing to correct keep the DC they have.
        # Otherwise a merge whose frames all share one gap (nothing to dead
        # reckon) would still be re-pinned from anchorIdx to file 0 and shift
        # wholesale -- Antarctica80/track-86/4683_0010, where every gap is 0.0
        # and the only seam step is a real co-registration discontinuity the
        # overlap solve already removes, would have moved 38.4 px despite
        # fitting Bp = -0.37 m today.
        if (abs(deadReckonBias[anchorIdx])
                <= abs(deadReckonBias[rangeAnchorIdx]) + RANGE_DISCONTINUITY_THRESH):
            rangeAnchorIdx = anchorIdx
    else:
        # Fall back to the old primary-only comparison if secondary geodats
        # aren't available for some reason -- better than no correction at all.
        # Say so: this silently reverts to the pre-fix behaviour that anchored
        # the merge on whichever frame happened to be leftmost, so a missing
        # geodat is enough to reintroduce the misaligned-anchor bug unnoticed.
        u.mywarning(f'{os.path.basename(outputVrt)}: no range-gate gap for '
                    f'{sum(g is None for g in gaps)} of {n} frames '
                    f'(missing geodat pair?) -- falling back to near-range dead '
                    f'reckoning; RangeOffsets anchoring is NOT gap-verified here')
        deadReckonBias = [(datasets[referenceIdx][1] - datasets[i][1]) / pixelSize for i in range(n)]

    # AzimuthOffsets gets the dead-reckoning fallback but NOT a moved anchor.
    # Frames carry an analogous primary-vs-secondary NominalTime slip (an
    # integer second = PRF single-look pixels, which predicts the measured
    # azimuth overlap step on 2139 of 2139 seams project-wide), so it can pin
    # frames the overlap chain cannot reach. The anchor itself must stay at
    # anchorIdx: mergedGeodat() takes the merged NominalTime from the first
    # frame, so the leftmost frame is already the convention azparams simulates
    # against, and azparams duly fits a sub-metre constant. Re-anchoring azimuth
    # onto the smallest-|slip| frame would instead inject 1520 px (~6800 m) of
    # constant into every merge whose first frame is slipped. The range side
    # differs only because the merged geodat's secondary MLNearRange is a min()
    # across frames rather than the first frame's, so it can land in a different
    # convention than frame 0 -- which is the bug rangeAnchorIdx exists to fix.
    azGaps = [_azimuthTimeGap(os.path.dirname(os.path.abspath(f))) for f in files]
    prf = _readGeodatProperty(masterGeodatPath, 'PRF') if masterGeodatPath else None
    if prf is None:
        for f in files:
            g = _findGeodat(os.path.dirname(os.path.abspath(f)))
            prf = _readGeodatProperty(g, 'PRF') if g else None
            if prf is not None:
                break
    azDeadReckonBias = None
    if prf and all(g is not None for g in azGaps):
        azDeadReckonBias = [(azGaps[anchorIdx] - g) * prf for g in azGaps]

    offsetsDict = {}          # (file, band) -> bias
    overlapResultsByBand = {}  # band -> {i: result}
    refStatsByBand = {}        # band -> [stats per file]
    for band in range(1, bandCount + 1):
        isRangeOffsetBand = bandDescriptions[band - 1] == 'RangeOffsets'
        # Bias used to pin frames the overlap chain cannot reach (see below).
        unreachableBias = deadReckonBias if isRangeOffsetBand else (
            azDeadReckonBias if bandDescriptions[band - 1] == 'AzimuthOffsets' else None)
        refStats = [None] * n
        if refData is not None:
            print(f"Computing reference phase statistics (band {band})...")
            for i, f in enumerate(files):
                refStats[i] = getReferenceStats(f, refData, refGt, maskData, maskGt, band=band)
                if refStats[i]:
                    print(f"  {os.path.basename(f)}: "
                          f"ref median={refStats[i]['median']:+.4f}  "
                          f"std={refStats[i]['std']:.4f}  n={refStats[i]['n']}")
                else:
                    print(f"  {os.path.basename(f)}: no valid reference pixels")
        refStatsByBand[band] = refStats

        # 3. Compute overlap stats first (before deciding the anchor system) --
        # the dead-reckoning fallback below needs to know which frames are
        # already connected to an anchor through real overlapping pixel data,
        # so it only fills in frames with no usable overlap instead of
        # competing with every frame's overlap-measured bias.
        overlapResults = {}
        if offsets:
            # Tie on the simulated geometry when it is available for every frame --
            # see SIM_FOR_BAND. No mask: the simulation is analytic and valid across
            # the whole footprint, and masking it away is exactly what leaves a seam
            # unbridgeable. Falls back per seam to the measured overlap.
            simFiles = simTieFiles(files, bandDescriptions[band - 1])
            simYOffsets = frameYOffsets(simFiles) if simFiles else None
            print(f"Calculating overlap offsets (band {band})"
                  f"{' from simulated geometry' if simFiles else ''}...")
            for i in range(n - 1):
                result = None
                if simFiles:
                    result = getOverlapStats(simFiles[i], simFiles[i + 1],
                                             yOffset1=simYOffsets[i],
                                             yOffset2=simYOffsets[i + 1],
                                             verbose=verbose, band=1)
                    if result is not None:
                        result = dict(result, source='sim')
                if result is None:
                    result = getOverlapStats(files[i], files[i + 1],
                                             yOffset1=yOffsets[i],
                                             yOffset2=yOffsets[i + 1],
                                             maskData=maskData, maskGt=maskGt,
                                             verbose=verbose, band=band)
                    if result is not None:
                        result = dict(result, source='measured')
                if result is None and simFiles:
                    # Neither raster spans the seam: a true time gap. The
                    # simulation is still known on both sides, so tie on it
                    # (simGapTie) instead of leaving the frame for the
                    # NominalTime pin below, which cannot see a frame whose
                    # simulation never ran and mis-steps by the missing constant.
                    result = simGapTie(simFiles[i], simFiles[i + 1],
                                       yOffset1=simYOffsets[i],
                                       yOffset2=simYOffsets[i + 1])
                    if result is not None:
                        result = dict(result, source='simgap')
                overlapResults[i] = result
        overlapResultsByBand[band] = overlapResults

        # 4. BUILD A / bVec for this band.
        # If reference provides valid constraints, those anchor the solution absolutely.
        # Otherwise fall back to pinning the geometrically leftmost file (anchorIdx) to zero.
        validRefIndices = [i for i, s in enumerate(refStats) if s is not None]
        if validRefIndices:
            A = []
            bVec = []
            for i in validRefIndices:
                row = [0] * n
                row[i] = 1
                A.append(row)
                bVec.append(-refStats[i]['median'])
            anchored = set(validRefIndices)
        else:
            # For the RangeOffsets band the absolute level has to be tied to the
            # frames that are actually aligned, so pin the dead-reckoned anchor
            # (rangeAnchorIdx) instead of whichever frame is geometrically
            # leftmost. The overlap equations below only enforce *relative*
            # continuity; whichever frame gets pinned sets the DC for the whole
            # merge. Pinning the geometrically leftmost frame meant that when
            # the first frame carried a range-gate anomaly, all the good frames
            # were dragged onto its DC -- which rparams then fits as a huge
            # spurious Bp (track-170/4248_0001 was -2288 SL px,
            # track-83/4161_0000 -1744).
            pinIdx = rangeAnchorIdx if isRangeOffsetBand else anchorIdx
            A = [[0] * n]
            A[0][pinIdx] = 1
            bVec = [deadReckonBias[pinIdx] if isRangeOffsetBand else 0]
            anchored = {pinIdx}

        if unreachableBias is not None:
            # Frames reachable from an anchored frame via an unbroken chain of
            # valid overlaps need nothing extra -- the overlap equations below
            # already constrain them from real pixel data, same as before this
            # fallback existed. Only frames NOT reachable that way (e.g. a
            # partial frame with no usable overlap with its neighbor) get
            # pinned to their dead-reckoned bias instead of being left
            # unconstrained.
            connected = set(anchored)
            changed = True
            while changed:
                changed = False
                for i in range(n - 1):
                    if overlapResults.get(i) is not None:
                        if i in connected and (i + 1) not in connected:
                            connected.add(i + 1)
                            changed = True
                        elif (i + 1) in connected and i not in connected:
                            connected.add(i)
                            changed = True
            # Pin one equation per unreachable COMPONENT, not one per frame.
            # Unreachable frames that still overlap EACH OTHER already have their
            # relative biases fixed by the overlap equations below; giving every
            # member its own pin makes those pins compete with the overlap chain
            # and least squares splits the difference, satisfying neither. That
            # manufactures a staircase that was never in the data. On
            # Antarctica80/track-23/3236_0000 all six raw seam medians were under
            # 0.03 px, yet per-frame pinning applied 42/84/211/549/1436/3758 px
            # -- x0 times 1, 2, 5, 13, 34, 89, the closed-form least-squares
            # answer for a chain pinned at every node (x[i+1] = 3x[i] - x[i-1]).
            # On Antarctica40/track-80/5542_0000 a two-frame component with a
            # real -466.995 px seam got +155.665/-155.665 instead of the -466.995
            # that would have closed it, turning one honest discontinuity into
            # two spurious ones. Pinning the component MEAN lets the overlap
            # equations set the shape and dead reckoning set only the DC it is
            # entitled to. For a one-frame component this is exactly the old row.
            unreachable = [i for i in range(n) if i not in connected]
            while unreachable:
                component = {unreachable[0]}
                changed = True
                while changed:
                    changed = False
                    for i in range(n - 1):
                        if overlapResults.get(i) is not None:
                            if i in component and (i + 1) not in component:
                                component.add(i + 1)
                                changed = True
                            elif (i + 1) in component and i not in component:
                                component.add(i)
                                changed = True
                members = sorted(component)
                row = [0.0] * n
                for i in members:
                    row[i] = 1.0 / len(members)
                A.append(row)
                bVec.append(sum(unreachableBias[i] for i in members) / len(members))
                unreachable = [i for i in unreachable if i not in component]

        if offsets:
            for i in range(n - 1):
                result = overlapResults.get(i)
                if result is not None:
                    row = [0] * n
                    row[i] = -1
                    row[i + 1] = 1
                    A.append(row)
                    bVec.append(result['median'])

        if offsets or refData is not None:
            bandOffsets = dict(zip(files, np.linalg.lstsq(A, bVec, rcond=None)[0]))
        else:
            bandOffsets = {}
        for f in files:
            offsetsDict[(f, band)] = bandOffsets.get(f, 0.0)

    # 3b. RANGE-DISCONTINUITY CHECK (RangeOffsets band only). The range-gate-gap
    # heuristic (deadReckonBias) predicts the inter-frame range bias from geodat
    # metadata; the overlap median measures it from real pixels. When they
    # disagree by more than RANGE_DISCONTINUITY_THRESH SL pixels -- e.g. a real
    # co-registration DC shift between frames that the (identical) near-range
    # metadata is blind to -- the empirical overlap correction in the solve above
    # is what actually removes the discontinuity; here we just record the
    # disagreement in a rangeDiscontinuityFound file beside the output so QC can
    # see which merges needed a data-driven rescue the metadata missed.
    # _updateRangeDiscontinuityFile also CLEARS this product's old section when a
    # re-merge no longer disagrees, so the call must be unconditional.
    discontinuityLines = []
    for band in range(1, bandCount + 1):
        if bandDescriptions[band - 1] != 'RangeOffsets':
            continue
        for i, result in overlapResultsByBand.get(band, {}).items():
            if result is None:
                continue
            predicted = deadReckonBias[i + 1] - deadReckonBias[i]
            empirical = result['median']
            disagreement = empirical - predicted
            if abs(disagreement) <= RANGE_DISCONTINUITY_THRESH:
                continue
            f1, f2 = files[i], files[i + 1]
            discontinuityLines += [
                f'# {datetime.now().isoformat(timespec="seconds")}  '
                f'heuristic vs empirical disagree by {disagreement:+.4f} SL px '
                f'(threshold {RANGE_DISCONTINUITY_THRESH:g}); empirical overlap '
                f'correction applied in merge',
                f'pair: {f1} -> {f2}',
                f'  rangeGateGap: {gaps[i]} m -> {gaps[i + 1]} m '
                f'(heuristic-predicted step {predicted:+.4f} px)',
                f'  overlap median (f1-f2): {empirical:+.4f} px  '
                f'(n={result["n"]} px)',
                f'  applied bias: f1={offsetsDict[(f1, band)]:+.4f} px, '
                f'f2={offsetsDict[(f2, band)]:+.4f} px',
            ]
    _updateRangeDiscontinuityFile(outputVrt, discontinuityLines)

    # 4. SETUP VRT
    minX = min(allX)
    maxX = max(allX)
    minY = min(allY)
    maxY = max(allY)

    outW = int(round((maxX - minX) / resX))
    outH = int(round((maxY - minY) / resY))

    if os.path.exists(outputVrt) and not overWrite:
        u.myerror("File exists.")

    driver = gdal.GetDriverByName("VRT")
    vrtDs = driver.Create(outputVrt, outW, outH, bandCount, dtype)
    vrtDs.SetProjection(proj)
    vrtDs.SetGeoTransform([minX, resX, 0, minY, 0, resY])
    if globalMetadata:
        vrtDs.SetMetadata(globalMetadata)

    # 5. ADD SOURCES
    scaleDict = {}
    offsetDict = {}
    for f, x0, y0, xSize, ySize in datasets:
        offX = int(round((x0 - minX) / resX))
        offY = int(round((y0 - minY) / resY))
        relPath = os.path.relpath(f, os.path.dirname(os.path.abspath(outputVrt)))

        srcDs = gdal.Open(f)
        for b in range(1, bandCount + 1):
            srcBand = srcDs.GetRasterBand(b)
            scaleDict[(f, b)] = srcBand.GetScale() if srcBand.GetScale() is not None else 1.0
            offsetDict[(f, b)] = srcBand.GetOffset() if srcBand.GetOffset() is not None else 0.0

        # Per band, and an explicit None test: band 1's nodata does not
        # necessarily speak for the rest of the file, and `or` would silently
        # turn a legitimate nodata of 0.0 into the -2e9 sentinel. Neither bites
        # the current products (nodata is uniform per file, and the simoffsets
        # products report None, which still lands on -2e9) -- this keeps that
        # true by construction rather than by luck.
        srcNodata = {}
        for b in range(1, bandCount + 1):
            nd = srcDs.GetRasterBand(b).GetNoDataValue()
            srcNodata[b] = -2000000000 if nd is None else nd

        for b in range(1, bandCount + 1):
            vrtBand = vrtDs.GetRasterBand(b)
            srcBand = srcDs.GetRasterBand(b)

            if bandDescriptions[b - 1]:
                vrtBand.SetDescription(bandDescriptions[b - 1])

            mdiDesc = srcBand.GetMetadataItem("Description")
            if mdiDesc:
                vrtBand.SetMetadataItem("Description", mdiDesc)
                vrtBand.SetDescription(mdiDesc)

            sourceXml = (
                f"<ComplexSource>\n"
                f'    <SourceFilename relativeToVRT="1">{relPath}</SourceFilename>\n'
                f"    <SourceBand>{b}</SourceBand>\n"
                f'    <SrcRect xOff="0" yOff="0" xSize="{xSize}" ySize="{ySize}"/>\n'
                f'    <DstRect xOff="{offX}" yOff="{offY}" xSize="{xSize}" ySize="{ySize}"/>\n'
                f"    <NOData>{srcNodata[b]}</NOData>\n"
                f"</ComplexSource>"
            )
            vrtBand.SetMetadataItem(f"source_{hash(f)}_{b}", sourceXml, "new_vrt_sources")

            if f == files[0]:
                vrtBand.SetNoDataValue(float('nan'))
    vrtDs = None

    # 6. POST-PROCESS: inject ScaleOffset/ScaleRatio into VRT XML
    print("Post-processing VRT to inject calculated offsets...")
    tree = ET.parse(outputVrt)
    root = tree.getroot()
    vrtDir = os.path.dirname(os.path.abspath(outputVrt))

    for band in root.findall("VRTRasterBand"):
        for source in band.findall("ComplexSource"):
            fnameElem = source.find("SourceFilename")
            if fnameElem is None:
                continue
            absFname = os.path.normpath(os.path.join(vrtDir, fnameElem.text))
            match = next((f for f in files
                          if os.path.normpath(os.path.abspath(f)) == absFname), None)
            bandIdx = int(band.get("band", 1))
            srcScale = scaleDict.get((match, bandIdx), 1.0)
            srcOffset = offsetDict.get((match, bandIdx), 0.0)
            bias = offsetsDict.get((match, bandIdx), 0.0) * srcScale

            so = ET.SubElement(source, "ScaleOffset")
            so.text = f"{srcOffset + bias:.12f}"
            sr = ET.SubElement(source, "ScaleRatio")
            sr.text = f"{srcScale:.12f}"

    tree.write(outputVrt, encoding="utf-8", xml_declaration=True)
    print(f"Done. Master VRT created: {outputVrt}")

    # 6b. BAKE: this VRT mosaics N sub-frame sources (each itself potentially a VRT),
    # which is expensive for downstream readers to resolve (mosaic3d's GDAL reads were
    # paying for N nested file opens per band). Flatten into a single GeoTIFF and
    # rewrap as a minimal VRT -- see bakeVrtInPlace().
    bakeVrtInPlace(outputVrt)

    # 7. SUMMARY
    if offsets or refData is not None:
        writeSummary(outputVrt, files, offsetsDict, overlapResultsByBand, refStatsByBand, bandCount)


def main():
    parser = argparse.ArgumentParser(
        description='Build a mosaic VRT from multiple single-frame VRTs, '
                    'optionally solving inter-frame DC offsets via overlap '
                    'statistics and/or a reference phase.',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('output', help='Output VRT path (with --bakeOnly: an '
                        'existing VRT to bake in place)')
    parser.add_argument('inputs', nargs='*',
                        help='Input VRT/TIF files (or glob patterns); not used '
                        'with --bakeOnly')
    parser.add_argument('--overWrite', action='store_true',
                        help='Overwrite output VRT if it already exists')
    parser.add_argument('--offsets', action='store_true',
                        help='Solve and apply inter-frame DC offsets via overlap statistics')
    parser.add_argument('--referencePhase', default=None, metavar='VRT',
                        help='Reference phase VRT for absolute bias anchoring '
                             '(e.g. simulated phase from siminsar)')
    parser.add_argument('--mask', default=None, metavar='FILE',
                        help='Mask file (GeoTIFF/VRT); masked pixels excluded '
                             'from bias estimation')
    parser.add_argument('--verbose', action='store_true',
                        help='Print per-overlap detail lines')
    parser.add_argument('--bakeOnly', action='store_true',
                        help='Skip mosaicking; just bake an existing recursive VRT '
                        '("output") into a flat GeoTIFF and rewrap it as a minimal '
                        'VRT. Much faster than rebuilding from sub-frame inputs -- '
                        'use this to retrofit virtual frames built before baking was '
                        'added to buildVrt(). No-op if the .tif already exists. Not '
                        'needed for newly-built virtual frames, which are baked '
                        'automatically.')
    args = parser.parse_args()
    if args.bakeOnly:
        bakeVrtInPlace(args.output, skipIfExists=True)
        return
    if not args.inputs:
        parser.error('inputs are required unless --bakeOnly is set')
    buildVrt(args.output, args.inputs, args.overWrite, args.offsets,
             referenceVrt=args.referencePhase, maskVrt=args.mask,
             verbose=args.verbose)


if __name__ == '__main__':
    main()
