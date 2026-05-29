# ROFFtoGrimp.py

Converts a NISAR ROFF HDF product (range/azimuth speckle-tracking offsets) into GrIMP-formatted binary flat-file products and VRTs ready for input to `mosaic3d`. Processes each of the three correlation layers independently (cull + interpolate), then merges them into a single weighted-average output.

---

## Usage

```
ROFFtoGrimp [options] ROFF
```

`ROFF` is the path to a NISAR ROFF HDF5 file (local path or S3 URI).

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `ROFF` | — | ROFF HDF5 file to convert. |
| `--geodat1 FILE` | auto | Primary geodat file (`geodatNLRxNLA.geojson`). Auto-detected from `*.nisar.uw` files in the output directory. |
| `--geodat2 FILE` | auto | Secondary geodat file (`geodatNLRxNLA.secondary.geojson`). Derived from `--geodat1` if not given. |
| `--DEM FILE` | region default | DEM file for offset simulation. |
| `--region NAME` | auto | Region name (`greenland` or `antarctica`). Auto-detected from the ROFF EPSG code if not given. |
| `--regionFile FILE` | None | YAML file with region-specific paths (velMap, DEM, mask, etc.). Overrides `--region`. |
| `--correlationThresholds F F F` | `0.07 0.05 0.025` | Correlation peak thresholds for layers 1, 2, and 3. Points below the threshold are discarded before culling. |
| `--outputDir DIR` | ROFF directory | Directory for all output files. |
| `--boxSize N` | 7 | Side length of the local statistics box used by `cullst`. |
| `--nGood N` | 10 | Minimum number of valid neighbours in the box; points with fewer are discarded. |
| `--maxR F` | 3.0 | Maximum allowed deviation from the local median for range offsets (pixels). |
| `--maxA F` | 3.0 | Maximum allowed deviation from the local median for azimuth offsets (pixels). |
| `--sr N` | 3 | Range smoothing kernel half-width for `cullst`. Odd → uniform kernel; even → tapered kernel. |
| `--sa N` | 3 | Azimuth smoothing kernel half-width for `cullst`. |
| `--interpThresh N` | 20 | Maximum hole area (pixels) to fill during interpolation. |
| `--islandThresh N` | 20 | Maximum isolated-region area (pixels) to discard. |
| `--byteOrder STR` | `MSB` | Byte order for binary output files (`MSB` or `LSB`). |
| `--noMask` | False | Skip applying the fast-area mask to layer 3. |
| `--mergeOnly` | False | Skip simulation, culling, and interpolation; re-run merge step only. |
| `--verbose` | False | Print all subprocess output to the terminal. |

---

## Processing flow

```
parseArgs()
  └── findGeodat()       — locate geodat1/geodat2 if not supplied

openHDF(ROFF)
resolveRegion()          — determine region from EPSG if not set

removeOutlierOffsets()   — discard pixels below per-layer correlation thresholds

mkdir workingDir/
writeOffsetsDatFile()    — write offsets.dat and offsets.geom.dat

setupGeodats()           — symlink geodat files into workingDir/

simulateOffsets()        — two parallel simoffsets calls:
  │   offsets.geom.*   — geometry-only (no velocity)
  └── offsets.*        — geometry + velocity

applyMask()              — mask fast-moving areas (unless --noMask)

writeData()              — write NISARoffsets.layer{1,2,3}.{dr,da,sr,sa,cor}

cullst()                 — run cullst on each layer in parallel

interpOffsets()          — run intfloat on each layer×band in parallel
                         — write per-layer *.cull.interp.vrt

mergeOffsets()           — nanmean across layers; add back geometry offsets;
                           scale sigmas by sqrt(N); write final binary files

writeVRTs()              — write azimuth.offsets.vrt, range.offsets.vrt,
                           offsets.range-azimuth.vrt
```

---

## Freestanding programs called

### `simoffsets`

Simulates range and azimuth offset fields from a velocity map and DEM. Called twice in parallel threads:

```
simoffsets -region=REGION [-LSB] [-dem=DEM] [-noVel]
    -offsetsDat=outputDir/offsets[.geom].dat
    -azOffsets=outputDir/offsets[.geom].da -syncDat
    -geodatFile=workingDir/geodat1
    -secondGeodatFile=workingDir/geodat2
```

The geometry-only call (with `-noVel`) produces the static offset baseline (`offsets.geom.*`). The full call (with velocity) produces `offsets.*`, which includes glacier motion.

### `cullst`

Spatial outlier filter for speckle-tracker offsets. Invoked once per layer (layers 1–3) in parallel threads. See `mosaicSource` documentation for parameter details.

### `intfloat`

Hole-filling interpolator. Invoked for each layer × component (`.cull.dr`, `.cull.da`, `.sr`, `.sa`) in parallel threads.

---

## Output files

Final products are written to `outputDir/`:

| File | Description |
|------|-------------|
| `range.offsets` | Merged range offset field (MSB float32) |
| `azimuth.offsets` | Merged azimuth offset field (MSB float32) |
| `range.offsets.sr` | Merged range offset sigma (MSB float32) |
| `azimuth.offsets.sa` | Merged azimuth offset sigma (MSB float32) |
| `range.offsets.vrt` | Two-band VRT: RangeOffsets, RangeSigma |
| `azimuth.offsets.vrt` | Two-band VRT: AzimuthOffsets, AzimuthSigma |
| `offsets.range-azimuth.vrt` | Four-band VRT combining all four fields |

Intermediate files in `outputDir/workingDir/`:

| File pattern | Description |
|--------------|-------------|
| `NISARoffsets.layer{1,2,3}.{dr,da,sr,sa}` | Per-layer offset/sigma binaries extracted from ROFF |
| `NISARoffsets.layer{1,2,3}.cull.{dr,da}` | Culled offset binaries |
| `NISARoffsets.layer{1,2,3}.cull.interp.{dr,da,sr,sa}` | Hole-filled offset/sigma binaries |
| `NISARoffsets.layer{1,2,3}.cull.interp.vrt` | Four-band VRT per layer |
| `offsets.{dr,da,sr,sa,vrt}` | Simulated offsets with velocity |
| `offsets.geom.{dr,da,vrt}` | Geometry-only simulated offsets |

---

## Layer merging

After culling and interpolation, the three layers are merged by:

1. Computing `nanmean` of range and azimuth offsets across layers.
2. Adding back the geometry-only simulated offsets (`offsets.geom.vrt`) at valid pixels.
3. Computing per-pixel sigma as `nanmean(sigma) × sqrt(N)`, where `N` is the number of valid layers at each pixel.

No-data pixels use fill value `−2×10⁹`.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `parseArgs()` | Parse command-line arguments; validate ROFF path; assemble `params` dict. |
| `findGeodat(params, geodat1, geodat2)` | Auto-detect geodat filenames from `*.nisar.uw` files in the output directory. |
| `setupGeodats(params)` | Create symlinks to geodat files inside `workingDir/`. |
| `resolveRegion(myROFF, params)` | Determine region (`greenland`/`antarctica`) from ROFF EPSG code. |
| `simulateOffsets(outputDir, baseName, params, ...)` | Spawn two `simoffsets` threads (geometry-only and full velocity). |
| `callSim(outputDir, baseName, params, ...)` | Build and execute a single `simoffsets` shell command. |
| `cullst(outputDir, baseName, ...)` | Spawn one `cullst` thread per layer. |
| `runCull(outputDir, baseLayerName, ...)` | Build and execute a single `cullst` shell command. |
| `interpOffsets(outputDir, baseName, ...)` | Spawn one `intfloat` thread per layer × component; write per-layer VRTs. |
| `runInterp(outputDir, inputFile, outputFile, nr, na, ...)` | Build and execute a single `intfloat` shell command. |
| `mergeOffsets(outputDir, baseName, ...)` | Read per-layer VRTs; nanmean; add geometry; write final binaries. |
| `writeVRTs(myROFF, ROFFPath, params)` | Write the three final GrIMP VRT products with full metadata. |
| `writeInterpVrt(newVRTFile, sourceFiles, descriptions, nr, na, ...)` | Create a GDAL VRT wrapping one or more binary flat-file bands. |
| `readVRTAndRenameBands(layerVRT, ...)` | Open a VRT with rioxarray and rename bands by metadata key. |
| `readVRTAndAppend(layerVRT, data)` | Read a VRT and append each band's array to a list in a dict. |
| `symlink_file(src_path, dst_path, ...)` | Create a relative or absolute symbolic link, optionally overwriting. |
