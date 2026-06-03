import os
import glob
import argparse
import numpy as np
import xml.etree.ElementTree as ET
from scipy.interpolate import RegularGridInterpolator

from osgeo import gdal

import utilities as u

gdal.UseExceptions()


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
                    maskData=None, maskGt=None):
    """Median difference (f1 - f2) in georeferenced overlap.

    yOffset1/yOffset2 are midnight-wrap corrections (seconds) applied to Y origins.
    Returns dict {'median': float, 'n': int} or None if overlap is insufficient.
    """
    ds1 = gdal.Open(file1)
    ds2 = gdal.Open(file2)

    gt1 = ds1.GetGeoTransform()
    gt2 = ds2.GetGeoTransform()

    nodata1 = ds1.GetRasterBand(1).GetNoDataValue()
    nodata2 = ds2.GetRasterBand(1).GetNoDataValue()

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
        print(f"  overlap({f1},{f2}): no spatial overlap — b1={b1} b2={b2}")
        return None

    def toPx(gt, yOffset, x, y):
        return (int(round((x - gt[0]) / gt[1])),
                int(round((y - (gt[3] + yOffset)) / gt[5])))

    x1, y1 = toPx(gt1, yOffset1, inter[0], inter[3] if gt1[5] < 0 else inter[2])
    x2, y2 = toPx(gt2, yOffset2, inter[0], inter[3] if gt2[5] < 0 else inter[2])

    pw = int(round((inter[1] - inter[0]) / gt1[1]))
    ph = int(round((inter[3] - inter[2]) / abs(gt1[5])))

    print(f"  overlap({f1},{f2}): inter y=[{inter[2]:.2f},{inter[3]:.2f}] "
          f"px1=({x1},{y1}) px2=({x2},{y2}) pw={pw} ph={ph}")

    try:
        arr1 = ds1.GetRasterBand(1).ReadAsArray(x1, y1, pw, ph).astype(np.float32)
        arr2 = ds2.GetRasterBand(1).ReadAsArray(x2, y2, pw, ph).astype(np.float32)

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
        print(f"  overlap({f1},{f2}): {nValid} valid pixels")
        if nValid < 100:
            return None

        diffs = arr1[valid] - arr2[valid]
        return {'median': float(np.median(diffs)), 'n': nValid}

    except Exception as e:
        print(f"  overlap({f1},{f2}): exception — {e}")
        return None


def getReferenceStats(filepath, refData, refGt, maskData=None, maskGt=None):
    """Median/std/n of (file - reference) over valid, unmasked pixels.

    Returns dict {'median': float, 'std': float, 'n': int} or None.
    """
    ds = gdal.Open(filepath)
    fileGt = ds.GetGeoTransform()
    fileNr, fileNc = ds.RasterYSize, ds.RasterXSize
    fileData = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nodata = ds.GetRasterBand(1).GetNoDataValue()
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


def writeSummary(outputVrt, files, offsetsDict, overlapResults, refStats):
    """Write bias-correction and fit statistics to a text file beside the VRT."""
    statsPath = outputVrt + '.stats'
    lines = ['# custom_buildvrtWithOffsets summary',
             f'# output: {outputVrt}', '']

    if overlapResults:
        lines += ['## Overlap residuals (after bias correction)',
                  '#  i  file_i -> file_{i+1}  raw_median  corrected_residual  n_pixels']
        for i, result in overlapResults.items():
            f1 = os.path.basename(files[i])
            f2 = os.path.basename(files[i + 1])
            if result is None:
                lines.append(f'  {i}  {f1} -> {f2}  NO_OVERLAP')
            else:
                rawMed = result['median']
                o1 = offsetsDict.get(files[i], 0.0)
                o2 = offsetsDict.get(files[i + 1], 0.0)
                residual = rawMed - (o2 - o1)
                lines.append(f'  {i}  {f1} -> {f2}  '
                              f'{rawMed:+.4f}  {residual:+.4f}  {result["n"]}')
        lines.append('')

    if any(s is not None for s in refStats):
        lines += ['## Reference phase fit (per frame)',
                  '#  file  ref_median  std  n_pixels  applied_offset  corrected_residual']
        for f, s in zip(files, refStats):
            fname = os.path.basename(f)
            applied = offsetsDict.get(f, 0.0)
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


def buildVrt(outputVrt, inputPattern, overWrite, offsets,
             referenceVrt=None, maskVrt=None):

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
    bandDescriptions = [firstDs.GetRasterBand(b).GetDescription()
                        for b in range(1, bandCount + 1)]

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

    # 2. COLLECT REFERENCE STATS (before building A so we can choose anchor strategy)
    n = len(files)
    refStats = [None] * n
    if refData is not None:
        print("Computing reference phase statistics...")
        for i, f in enumerate(files):
            refStats[i] = getReferenceStats(f, refData, refGt, maskData, maskGt)
            if refStats[i]:
                print(f"  {os.path.basename(f)}: "
                      f"ref median={refStats[i]['median']:+.4f}  "
                      f"std={refStats[i]['std']:.4f}  n={refStats[i]['n']}")
            else:
                print(f"  {os.path.basename(f)}: no valid reference pixels")

    # 3. BUILD A / bVec
    # If reference provides valid constraints, those anchor the solution absolutely.
    # Otherwise fall back to pinning file 0 to zero.
    validRefIndices = [i for i, s in enumerate(refStats) if s is not None]
    if validRefIndices:
        A = []
        bVec = []
        for i in validRefIndices:
            row = [0] * n
            row[i] = 1
            A.append(row)
            bVec.append(-refStats[i]['median'])
    else:
        A = [[0] * n]
        A[0][0] = 1
        bVec = [0]

    overlapResults = {}
    if offsets:
        print("Calculating overlap offsets...")
        for i in range(n - 1):
            result = getOverlapStats(files[i], files[i + 1],
                                     yOffset1=yOffsets[i],
                                     yOffset2=yOffsets[i + 1],
                                     maskData=maskData, maskGt=maskGt)
            overlapResults[i] = result
            if result is not None:
                row = [0] * n
                row[i] = -1
                row[i + 1] = 1
                A.append(row)
                bVec.append(result['median'])

    if offsets or refData is not None:
        offsetsDict = dict(zip(files, np.linalg.lstsq(A, bVec, rcond=None)[0]))
    else:
        offsetsDict = {}

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

        srcNodata = srcDs.GetRasterBand(1).GetNoDataValue() or -2000000000

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
                f"    <NOData>{srcNodata}</NOData>\n"
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
            match = next((k for k in offsetsDict
                          if os.path.normpath(os.path.abspath(k)) == absFname), None)
            bandIdx = int(band.get("band", 1))
            srcScale = scaleDict.get((match, bandIdx), 1.0)
            srcOffset = offsetDict.get((match, bandIdx), 0.0)
            bias = offsetsDict.get(match, 0.0) * srcScale

            so = ET.SubElement(source, "ScaleOffset")
            so.text = f"{srcOffset + bias:.12f}"
            sr = ET.SubElement(source, "ScaleRatio")
            sr.text = f"{srcScale:.12f}"

    tree.write(outputVrt, encoding="utf-8", xml_declaration=True)
    print(f"Done. Master VRT created: {outputVrt}")

    # 7. SUMMARY
    if offsets or refData is not None:
        writeSummary(outputVrt, files, offsetsDict, overlapResults, refStats)


def main():
    parser = argparse.ArgumentParser(
        description='Build a mosaic VRT from multiple single-frame VRTs, '
                    'optionally solving inter-frame DC offsets via overlap '
                    'statistics and/or a reference phase.')
    parser.add_argument('output', help='Output VRT path')
    parser.add_argument('inputs', nargs='+',
                        help='Input VRT/TIF files (or glob patterns)')
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
    args = parser.parse_args()
    buildVrt(args.output, args.inputs, args.overWrite, args.offsets,
             referenceVrt=args.referencePhase, maskVrt=args.mask)


if __name__ == '__main__':
    main()
