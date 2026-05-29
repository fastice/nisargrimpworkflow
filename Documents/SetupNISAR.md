# SetupNISAR.py — NISAR HDF5 to GrIMP Workflow Converter

## Overview

`SetupNISAR.py` converts NISAR Level-2 HDF5 products (RUNW, ROFF) into the
binary flat-file formats used by the GrIMP velocity mosaic pipeline.  It
operates on all frames of a single orbit pair found in the current working
directory, then consolidates the per-frame outputs into a single virtual-frame
directory using GDAL VRT mosaics.

The script is located at:

```
nisargrimpworkflow/nisargrimpworkflow/SetupNISAR.py
```

---

## Directory Layout

### Expected input structure

Run from the directory containing the per-frame subdirectories:

```
<workDir>/
├── <orbit1>_<frame>/          e.g. 12345_010/
│   └── NISAR*RUNW*.h5         unwrapped interferogram HDF5
│   └── NISAR*ROFF*.h5         range/azimuth offset HDF5
├── <orbit1>_<frame>/          e.g. 12345_020/
│   └── ...
└── ...
```

Frame directory names must match the pattern `<orbit1>_NNN` where `NNN` is a
0–999 frame number (zero-padded to 2 or 3 digits).

### Output structure

Per-frame products are written into the same `<orbit1>_<frame>/` directories.
A virtual frame directory consolidates all frames:

```
<workDir>/
├── <orbit1>_<frame>/
│   ├── <orbit1>_<frame>.<orbit2>_<frame>.*x*.vrt   phase/offset products
│   ├── geodat<R>x<A>.geojson                        reference geodat
│   ├── geodat<R>x<A>.secondary.geojson              secondary geodat
│   └── range.offsets, azimuth.offsets, ...          ROFF products
│
└── <orbit1>_<virtualFrame>/                         e.g. 12345_0000/
    ├── <orbit1>_0000.<orbit2>_0000.*.uw.interp.vrt  merged unwrapped phase
    ├── <orbit1>_0000.<orbit2>_0000.*.ion.filt.vrt   merged iono-corrected phase
    ├── <orbit1>_0000.<orbit2>_0000.*.cor.vrt         merged coherence
    ├── <orbit1>_0000.<orbit2>_0000.*.ion.filt.rangeOffset.vrt
    ├── <orbit1>_0000.<orbit2>_0000.*.ion.unfilt.rangeOffset.vrt
    ├── range.offsets.vrt                            merged ROFF range offsets
    ├── azimuth.offsets.vrt                          merged ROFF azimuth offsets
    ├── offsets.vrt, offsets.geom.vrt, ...           additional ROFF products
    ├── geodat<R>x<A>.geojson                        merged reference geodat
    ├── geodat<R>x<A>.secondary.geojson              merged secondary geodat
    ├── sensor.NISAR<BW>.yaml                        sensor parameters
    └── <orbit1>.<orbit2>.pairinfo                   orbit pair metadata
```

---

## Usage

```
SetupNISAR.py <orbit1> [options]
```

Must be run from the directory containing the `<orbit1>_<frame>` subdirectories.

### Positional argument

| Argument | Description |
|----------|-------------|
| `orbit1` | Reference orbit number (integer) |

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--virtualFrame NNNN` | `0000` | Name suffix for the consolidated virtual frame directory |
| `--firstFrame N` | 0 | Skip frames numbered below this value (0–999) |
| `--lastFrame N` | 999 | Skip frames numbered above this value (0–999) |
| `--overWrite` | off | Re-run RUNW and ROFF conversion even if output products already exist |
| `--overWritePhase` | off | Re-run RUNW conversion only if phase products already exist |
| `--allowMixedMode` | off | Include frames whose SLC granule names contain `_M_` (mixed-mode acquisitions are skipped by default) |
| `--RUNWOnly` | off | Process only RUNW (phase) products; skip ROFF and power |
| `--noMask` | off | Do not apply the fast-region mask to layer 3 during ROFF conversion |
| `--verbose` | off | Print all subprocess output to terminal (default: suppressed) |

### Examples

```bash
# Process all frames for orbit 12345, consolidate into 12345_0000/
SetupNISAR.py 12345

# Process a subset of frames only
SetupNISAR.py 12345 --firstFrame 10 --lastFrame 30

# Reprocess from scratch
SetupNISAR.py 12345 --overWrite

# Reprocess phase products only (keep existing ROFF)
SetupNISAR.py 12345 --overWritePhase

# Phase products only, custom virtual frame name
SetupNISAR.py 12345 --RUNWOnly --virtualFrame 0001

# Include mixed-mode frames, print all subprocess output
SetupNISAR.py 12345 --allowMixedMode --verbose
```

---

## Processing Pipeline

### 1. Frame discovery

Scans the current directory for subdirectories matching `<orbit1>_NNN` and
builds a sorted list of frame numbers.  Frames outside `--firstFrame` /
`--lastFrame` are excluded.

### 2. Secondary orbit and metadata extraction

Opens the first available RUNW or ROFF HDF5 found across all frames to
extract:
- Secondary orbit number
- Range bandwidth (MHz) → determines which sensor YAML to use (20/40/80 MHz)
- Number of range and azimuth looks
- Reference and secondary acquisition datetimes

### 3. Per-frame RUNW processing (`RUNWtoGrimp.py`)

For each frame, calls `RUNWtoGrimp.py` to convert the NISAR RUNW HDF5 into
GrIMP format VRTs.  Skips if all expected output products already exist (unless
`--overWrite` or `--overWritePhase`).

Output products per frame (as VRTs):

| Product | Description |
|---------|-------------|
| `*.uw.interp.vrt` | Unwrapped, interpolated phase |
| `*.ion.filt.vrt` | Ionosphere-corrected phase |
| `*.cor.vrt` | Interferometric coherence |
| `*.ion.filt.rangeOffset.vrt` | Ionosphere correction as range offset |
| `*.ion.unfilt.rangeOffset.vrt` | Unfiltered ionosphere range offset |

Mixed-mode frames (SLC granule name contains `_M_`) are silently skipped
unless `--allowMixedMode` is set.

Each frame also produces a reference and secondary geodat GeoJSON
(`geodat<R>x<A>.geojson`) recording image geometry for geocoding.

### 4. Per-frame ROFF processing (`ROFFtoGrimp`)

For each frame, calls `ROFFtoGrimp` to convert the NISAR ROFF HDF5 into
GrIMP offset maps.  Skipped if `--RUNWOnly` is set.

Output products per frame:

| Product | Description |
|---------|-------------|
| `range.offsets` | Range pixel offsets (binary) |
| `azimuth.offsets` | Azimuth pixel offsets (binary) |
| `offsets.vrt` | Combined offset VRT |
| `offsets.geom.vrt` | Geometrically corrected offsets |
| `offsets.ll.vrt` | Offsets in lat/lon |
| `offsets.mask.vrt` | Offset quality mask |
| `offsets.range-azimuth.vrt` | Range-azimuth offset pair |

### 5. Per-frame power image collection

Locates `P<orbit>_<frame>.*x*.pow` power image files and their companion
geodat GeoJSONs for later virtual-frame consolidation.

### 6. Virtual frame directory creation

Creates `<orbit1>_<virtualFrame>/` and consolidates all per-frame VRTs into
single mosaic VRTs using `custom_buildvrtWithOffsets.py`.  The VRT mosaic
preserves pixel offsets between frames so geocoding sees a continuous product.

For RUNW products the `ion.filt.rangeOffset` VRT path is recorded as
`ionosphereRangeOffsetCorrection` metadata on the `range.offsets.vrt`, so the
geocoding pipeline applies the ionosphere correction automatically.

### 7. Geodat merging

Merges per-frame GeoJSON geodat files into a single geodat for the virtual
frame by:
- Taking geometry (corners) from the first and last frames and combining them
  according to pass direction (ascending: take far-end corners from last frame;
  descending: take near-end corners)
- Merging and interpolating orbital state vectors (position + velocity) across
  all frames onto a uniform time grid using cubic spline interpolation
- Updating the azimuth size to match the merged VRT dimensions

### 8. Sensor YAML copy and update

Copies the appropriate sensor parameter YAML (`NISAR20.yaml`, `NISAR40.yaml`,
or `NISAR80.yaml`) from the `sarfunc` package into the virtual frame directory
as `sensor.NISAR<BW>.yaml`.  Updates `intLooksR` and `intLooksA` from the
actual HDF5 look counts.  Skips the copy if the file already exists.

### 9. Pair info file

Writes `<orbit1>.<orbit2>.pairinfo` into the virtual frame directory
containing one line:

```
<orbit1>  <orbit2>  <refDateTime>  <secDateTime>  <nRangeLooks>  <nAzLooks>
```

---

## Skip / Overwrite Logic

| Condition | Behaviour |
|-----------|-----------|
| RUNW output VRTs present, no `--overWrite`/`--overWritePhase` | Skip RUNW frame |
| `--overWrite` | Re-run both RUNW and ROFF for all frames |
| `--overWritePhase` | Re-run RUNW only; leave ROFF untouched |
| ROFF `range.offsets` present, no `--overWrite` | Skip ROFF frame |
| Mixed-mode frame, no `--allowMixedMode` | Skip RUNW silently |
| Virtual frame VRT build fails | Writes `<vrtFile>.fail` sentinel; previous `.fail` removed before each attempt |

---

## Output Files Reference

### Virtual frame directory (`<orbit1>_<virtualFrame>/`)

| File | Written by | Description |
|------|-----------|-------------|
| `*.uw.interp.vrt` | `createVirtualFrameRUNW` | Merged unwrapped phase |
| `*.ion.filt.vrt` | `createVirtualFrameRUNW` | Merged iono-corrected phase |
| `*.cor.vrt` | `createVirtualFrameRUNW` | Merged coherence |
| `*.ion.filt.rangeOffset.vrt` | `createVirtualFrameRUNW` | Iono correction as range offset |
| `*.ion.unfilt.rangeOffset.vrt` | `createVirtualFrameRUNW` | Unfiltered iono range offset |
| `range.offsets.vrt` | `createVirtualFrameROFF` | Merged range offsets (with iono metadata) |
| `azimuth.offsets.vrt` | `createVirtualFrameROFF` | Merged azimuth offsets |
| `offsets.vrt` etc. | `createVirtualFrameROFF` | Additional offset products |
| `geodat<R>x<A>.geojson` | `mergedGeodat` | Merged reference image geometry |
| `geodat<R>x<A>.secondary.geojson` | `mergedGeodat` | Merged secondary image geometry |
| `sensor.NISAR<BW>.yaml` | `copy_sensor_yaml` | Sensor parameters with look counts |
| `<orbit1>.<orbit2>.pairinfo` | `writePairInfo` | Orbit pair metadata |
| `*.vrt.fail` | VRT build loop | Sentinel: written if VRT build fails |

---

## Dependencies

| Tool / Package | Role |
|----------------|------|
| `nisarhdf` | Read NISAR HDF5 products; format GeoJSON |
| `RUNWtoGrimp.py` | Convert RUNW HDF5 to GrIMP phase/coherence VRTs |
| `ROFFtoGrimp` | Convert ROFF HDF5 to GrIMP offset products |
| `custom_buildvrtWithOffsets.py` | Build VRT mosaics with inter-frame pixel offsets |
| `sarfunc` | Provides sensor YAML templates (via `importlib.resources`) |
| `gdal` / `osgeo` | Write VRT metadata (ionosphere correction tag) |
| `geojson` | Read/write geodat GeoJSON files |
| `scipy.interpolate` | Cubic spline interpolation of state vectors |
