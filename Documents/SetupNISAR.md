# SetupNISAR — NISAR HDF5 to GrIMP Workflow Orchestrator

## Overview

`SetupNISAR` is the top-level orchestrator for converting a full set of NISAR
Level-2 HDF5 products for one orbit pass into the binary flat-file and VRT
formats consumed by the GrIMP velocity mosaic pipeline.  It loops over every
`<orbit1>_<frame>` directory found in the current working directory, drives the
per-frame conversion tools ([RUNWtoGrimp](RUNWtoGrimp.md),
[ROFFtoGrimp](ROFFtoGrimp.md), and by default
[estimateIonosphere](estimateIonosphere.md)), then consolidates all per-frame
outputs into a single virtual-frame directory using GDAL VRT mosaics.

---

## Directory Layout

### Required input structure

`SetupNISAR` is run from the **track directory** inside a project tree:

```
projectDir/
├── project.yaml               (optional — see Region configuration)
└── track-<N>/                 ← run SetupNISAR from here
    ├── <orbit1>_<frame>/      e.g. 1830_35/
    │   └── H5/
    │       ├── NISAR*RUNW*.h5
    │       ├── NISAR*ROFF*.h5
    │       └── NISAR*RIFG*.h5 (optional)
    ├── <orbit1>_<frame>/
    │   └── ...
    └── <orbit1>_<virtualFrame>/   created by SetupNISAR
```

Frame directory names must match `<orbit1>_NNN` where `NNN` is a 0–999 frame
number (zero-padded to 2 or 3 digits).

### Output structure per frame

Processing fills each frame directory with several subdirectories and products:

```
<orbit1>_<frame>/
├── H5/
│   ├── NISAR*RUNW*.vrt                  — HDF5 band VRT (wrapH5sInFrameDir)
│   ├── NISAR*ROFF*.layer{1,2,3}.vrt     — per-layer offset VRTs
│   └── NISAR*RIFG*.vrt
├── workingDir/                           — ROFFtoGrimp intermediates
│   ├── NISARoffsets.layer{1,2,3}.{dr,da,sr,sa,cc,mt}  — per-layer extracted binaries
│   ├── NISARoffsets.layer{1,2,3}.dat                   — metadata sidecar
│   ├── NISARoffsets.layer{1,2,3}.cull.{dr,da,cc,mt}   — after cullst
│   ├── NISARoffsets.layer{1,2,3}.cull.interp.{dr,da,sr,sa}  — after intfloat
│   ├── NISARoffsets.layer{1,2,3}.{vrt,mt.vrt,cull.vrt,cull.mt.vrt,cull.interp.vrt}
│   └── offsets.{geom,velocity}.{dr,da,dat,lat,lon,mask,poly,simdat}  — simoffsets binaries
├── offsetSims/                           — simoffsets VRTs and geodats (ROFFtoGrimp)
│   ├── offsets.geom.vrt                  — geometry-only simulated offsets
│   ├── offsets.geom.ll.vrt
│   ├── offsets.geom.mask.vrt
│   ├── offsets.velocity.vrt              — geometry + velocity simulated offsets
│   ├── offsets.velocity.ll.vrt
│   ├── offsets.velocity.mask.vrt
│   ├── geodat{R}x{A}.geojson            — symlinked geodats for simoffsets
│   └── geodat{R}x{A}.secondary.geojson
├── simPhase/                             — siminsar outputs (estimateIonosphere)
│   ├── velSim.vrt                        — velocity-derived interferometric phase (rad), no topography
│   └── maskVel.vrt                       — 0 where speed > velocityThreshold, 1 elsewhere
├── range.offsets                         — merged range offset binary
├── range.offsets.sr                      — merged range sigma binary
├── range.offsets.vrt
├── azimuth.offsets                       — merged azimuth offset binary
├── azimuth.offsets.sa                    — merged azimuth sigma binary
├── azimuth.offsets.vrt
├── offsets.range-azimuth.vrt
├── {pair}.nisar.cor                      — coherence binary (always written)
├── {pair}.nisar.cor.vrt
├── {pair}.correctedUnwrappedPhase.tif    — iono/ambiguity-corrected phase (default path)
├── {pair}.correctedUnwrappedPhase.vrt    ← primary phase input to mosaic
├── {pair}.ionosphereCorrection.tif
├── {pair}.ionosphereCorrection.vrt
├── {stem}.ionosphereCorrection.offset.tif  — iono on offset grid, metres
├── {stem}.ionosphereCorrection.offset.vrt
├── geodat{R}x{A}.geojson
├── geodat{R}x{A}.secondary.geojson
└── <orbit1>.<orbit2>.pairinfo
```

`{pair}` = `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.nisar`.
`{stem}` = `{orbit1}_{frame}.{orbit2}_{frame}` (looks and `.nisar` suffix stripped for offset-grid products).

In the `--phaseDerivedIonosphere` path, `RUNWtoGrimp` additionally writes
`{pair}.nisar.uw`, `.nisar.uw.interp`, `.nisar.ion`, and `.nisar.ion.filt`;
these are not produced in the default path.

### Virtual frame directory

```
<orbit1>_<virtualFrame>/                  e.g. 1830_0000/
├── {pair_vf}.correctedUnwrappedPhase.vrt — merged corrected phase
├── {pair_vf}.cor.vrt                     — merged coherence
├── {pair_vf}.ionosphereCorrection.vrt    — merged iono correction
├── {pair_vf}.ionosphereCorrection.offset.vrt
├── range.offsets.vrt                     — merged ROFF range offsets
│                                           (carries ionosphereRangeOffsetCorrection metadata)
├── azimuth.offsets.vrt
├── offsets.geom.vrt
├── offsets.geom.ll.vrt
├── offsets.geom.mask.vrt
├── offsets.velocity.vrt
├── offsets.velocity.ll.vrt
├── offsets.velocity.mask.vrt
├── offsets.range-azimuth.vrt
├── velSim.vrt                            — merged simulated phase (from simPhase/)
├── maskVel.vrt                           — merged velocity mask
├── geodat{R}x{A}.geojson                 — merged reference geometry
├── geodat{R}x{A}.secondary.geojson
├── sensor.NISAR{BW}.yaml                 — sensor parameters
└── <orbit1>.<orbit2>.pairinfo
```

---

## Usage

```
SetupNISAR <orbit1> [options]
```

Run from the directory containing the `<orbit1>_<frame>` subdirectories.

### Positional argument

| Argument | Description |
|----------|-------------|
| `orbit1` | Reference orbit number (integer) |

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--virtualFrame NNNN` | `0000` | Suffix for the consolidated virtual frame directory |
| `--firstFrame N` | 0 | Skip frames numbered below this value |
| `--lastFrame N` | 999 | Skip frames numbered above this value |
| `--overWrite` | off | Re-run all per-frame conversions even if outputs exist |
| `--overWritePhase` | off | Re-run RUNW and ionosphere steps only; leave ROFF untouched |
| `--allowMixedMode` | off | Include mixed-mode frames (SLC granule name contains `_M_`); skipped by default |
| `--RUNWOnly` | off | Process only RUNW (phase/coherence) — skip ROFF, ionosphere, and power |
| `--noMask` | off | Do not apply the fast-region mask to ROFF layer 3 |
| `--ompThreads N` | 6 | OpenMP thread count passed to `RUNWtoGrimp` |
| `--phaseDerivedIonosphere` | off | Use NISAR-embedded ionosphere screen inside `RUNWtoGrimp` instead of running `estimateIonosphere` (see [Two ionosphere paths](#two-ionosphere-paths)) |
| `--outputAll` | off | Pass `--outputAll` to `estimateIonosphere`: write all 5 intermediate bands instead of the standard 3 outputs |
| `--phaseThresh RAD` | 14π | Pass to `estimateIonosphere`: mask `correctedUnwrappedPhase` where \|correctedPhase − simPhase\| ≥ RAD radians — screens regions of likely incorrect unwrapping |
| `--correlationOnly` | off | Extract coherence and geodat files only — skips ROFF conversion, ionosphere estimation, and virtual-frame assembly |
| `--verbose` | off | Print all subprocess output to terminal (default: suppressed) |

### Examples

```bash
# Standard run — all frames for orbit 1830
SetupNISAR 1830

# Re-process only phase and ionosphere products (ROFF already done)
SetupNISAR 1830 --overWritePhase

# Process a subset of frames, tighter phase residual screen
SetupNISAR 1830 --firstFrame 30 --lastFrame 50 --phaseThresh 7.0

# Phase products only, custom virtual frame name
SetupNISAR 1830 --RUNWOnly --virtualFrame 0001

# Include mixed-mode frames, print all subprocess output
SetupNISAR 1830 --allowMixedMode --verbose

# Extract coherence and geodats only — no ROFF, no ionosphere, no virtual frame
SetupNISAR 1830 --correlationOnly
```

---

## Region configuration

On startup `SetupNISAR` looks for `../project.yaml` (i.e. in `projectDir/`).
If found and it contains a `regionFile` key, that YAML path is forwarded to
[ROFFtoGrimp](ROFFtoGrimp.md) and
[estimateIonosphere](estimateIonosphere.md) so they use the correct
region-specific DEM, velocity map, and ice mask without requiring per-frame
`--regionFile` flags.

---

## Processing pipeline

The steps below occur in the order shown.  Steps inside the per-frame loop
repeat for every frame; the virtual-frame steps run once after all frames are
done.

### Step 1 — Frame discovery

Scans the current directory for subdirectories matching `<orbit1>_NNN` and
builds a sorted list of frame numbers.  Frames outside `--firstFrame` /
`--lastFrame` are excluded.

### Step 2 — Secondary orbit and metadata extraction

Opens the first available RUNW or ROFF HDF5 across all frames to extract:

- Secondary orbit number
- Range bandwidth (MHz) → determines sensor YAML (20 / 40 / 80 MHz)
- Number of range and azimuth looks
- Reference and secondary acquisition datetimes

If no HDF5 can be opened, `haveData = False` and all per-frame and
virtual-frame InSAR steps are skipped (power processing still runs if `.pow`
files are present).

### Step 3 — Per-frame loop

For every frame, the following sub-steps execute in order:

#### 3a. HDF5 wrapping (`wrapH5sInFrameDir`)

Wraps all NISAR HDF5 files in `<orbit1>_<frame>/H5/` with GDAL VRTs that
point directly into the HDF5 bands — no data is extracted.  Recognised product
types and their outputs:

| H5 product | VRT written |
|------------|-------------|
| `RUNW` | `<h5root>.vrt` (5 bands: unwrappedPhase, ionospherePhaseScreen, ionospherePhaseScreenUncertainty, coherenceMagnitude, connectedComponents) |
| `ROFF` | `<h5root>.layer{1,2,3}.vrt` (4 bands per layer: slantRangeOffset, alongTrackOffset, correlationSurfacePeak, snr) |
| `RIFG` | `<h5root>.vrt` (2 bands: wrappedInterferogram, coherenceMagnitude) |
| `RSLC` | **Silently skipped** — not yet supported; no VRT is created |

VRTs already present are not overwritten.

#### 3b. RUNW processing ([RUNWtoGrimp](RUNWtoGrimp.md))

Calls `RUNWtoGrimp` to convert the RUNW HDF5 into GrIMP-format files.  The
command differs depending on the ionosphere path (see
[Two ionosphere paths](#two-ionosphere-paths) below).

Mixed-mode frames (SLC granule contains `_M_`) are silently skipped unless
`--allowMixedMode` is set.

Skip logic: the frame is skipped if the expected output files already exist and
neither `--overWrite` nor `--overWritePhase` is set.

**Default path** (`--phaseDerivedIonosphere` not set):

`RUNWtoGrimp` is called with `--noPhase --noIon`, which writes only the
geodat GeoJSONs and coherence products.  Phase and ionosphere outputs are
deferred to [estimateIonosphere](estimateIonosphere.md) in step 3d.

Outputs written by `RUNWtoGrimp --noPhase --noIon`:

| File | Description |
|------|-------------|
| `geodat{R}x{A}.geojson` | Reference image geometry (GeoJSON) |
| `geodat{R}x{A}.secondary.geojson` | Secondary image geometry (GeoJSON) |
| `{pair}.nisar.cor` | Coherence magnitude (binary MSB float32) |
| `{pair}.nisar.cor.vrt` | Coherence VRT |

**Phase-derived ionosphere path** (`--phaseDerivedIonosphere` set):

`RUNWtoGrimp` is called with `--phaseDerivedIonosphere`, which uses the
NISAR-embedded ionosphere phase screen to correct the unwrapped phase without
independent range offsets.  See [RUNWtoGrimp](RUNWtoGrimp.md) for details.

Outputs include `*.uw.interp.vrt`, `*.ion.filt.rangeOffset.vrt`,
`*.ion.unfilt.rangeOffset.vrt`, coherence, and geodats.

#### 3c. ROFF processing ([ROFFtoGrimp](ROFFtoGrimp.md))

Calls `ROFFtoGrimp` on the ROFF HDF5.  Skipped if `--RUNWOnly`, `--correlationOnly`, or
the frame is mixed-mode.

Skip logic: the frame is skipped if `range.offsets` already exists and
`--overWrite` is not set.

Inside `ROFFtoGrimp` the following happen (see [ROFFtoGrimp](ROFFtoGrimp.md)
for full details):

1. Reads range/azimuth offsets from the ROFF HDF5; discards pixels below
   per-layer correlation peak thresholds.
2. Writes `.dat` metadata files (`offsets.dat`, `offsets.geom.dat`) for
   `simoffsets`.
3. Runs `simoffsets` twice in parallel threads:
   - geometry-only (no velocity) → `offsetSims/offsets.geom.*`
   - geometry + velocity → `offsetSims/offsets.velocity.*`
   Both sets of VRTs are stamped with the NISAR zeroDoppler-time / slant-range
   geotransform so they align correctly with RUNW products.
4. Optionally applies a fast-area mask to layer 3 (unless `--noMask`).
5. Writes per-layer binary flat files to `workingDir/`.
6. Runs `cullst` per layer (threaded spatial outlier filter).
7. Runs `intfloat` per layer × component (threaded hole-filling).
8. Merges three layers by nanmean; adds geometry offsets back; writes final
   binaries `range.offsets`, `azimuth.offsets` and their sigma files.
9. Writes final VRTs: `range.offsets.vrt`, `azimuth.offsets.vrt`,
   `offsets.range-azimuth.vrt`.

#### 3d. Ionosphere estimation ([estimateIonosphere](estimateIonosphere.md))

*Default path only* — skipped when `--phaseDerivedIonosphere`, `--RUNWOnly`,
or `--correlationOnly` is set.

Calls `processFrameIonosphere`, which drives `estimateIonosphere` as a
subprocess from the frame directory (`cwd = <orbit1>_<frame>/`).  Input files
are resolved automatically:

| Input | Resolved from |
|-------|--------------|
| RUNW HDF5 | `H5/NISAR*RUNW*.h5` (or frame root) |
| Range offset VRT | `range.offsets.vrt` |
| Geometry offset VRT | `offsetSims/offsets.geom.vrt` |
| Simulation dir | `simPhase/` (created if absent) |

Inside `estimateIonosphere` the following happen (see
[estimateIonosphere](estimateIonosphere.md) for full theory and details):

1. Runs `siminsar` to produce `simPhase/velSim.vrt` (simulated interferometric
   phase in radians) and `simPhase/maskVel.vrt` (velocity threshold mask).
   Both are stamped with the zeroDoppler-time / slant-range geotransform.
   Skipped if the files already exist and `--overWrite` is not set.
2. Loads unwrapped phase from the RUNW HDF5; masks zero connected-component
   pixels.
3. Computes the range-offset-derived phase: `offset_phase = (meas − geom) × (−4π/λ)`.
4. Iteratively estimates and removes the ionospheric phase screen (up to 8
   passes), using `ambiguityAnalysis` to detect and correct integer-cycle
   unwrapping errors per connected component by comparison with `velSim`.
5. If `--phaseThresh RAD` is set (default 14π), masks any pixel where
   |correctedPhase − simPhase| ≥ RAD, screening regions with likely incorrect
   unwrapping.  Masked pixels are set to −2×10⁹ in the final TIF so the C
   geocoder identifies them as noData.
6. Optionally hole-fills `correctedUnwrappedPhase` via `intfloat`
   (unless `--noInterp`).  cc=0 and phaseThresh-masked pixels are re-applied
   after intfloat so they cannot be filled.
7. Regrids the ionosphere correction to the offset VRT coordinate system and
   writes `*.ionosphereCorrection.offset.tif` with −2×10⁹ noData.

Frame is skipped if `*.ionosphereCorrection.vrt` already exists (or the 5-band
output VRT when `--outputAll`) and neither `--overWrite` nor `--overWritePhase`
is set.

Outputs written per frame:

| File | Description |
|------|-------------|
| `{pair}.correctedUnwrappedPhase.tif` | Corrected phase (float32 GeoTIFF, noData = −2×10⁹) |
| `{pair}.correctedUnwrappedPhase.vrt` | Single-band VRT, band description `Phase` |
| `{pair}.ionosphereCorrection.tif` | Ionosphere estimate (radians) |
| `{pair}.ionosphereCorrection.vrt` | Single-band VRT |
| `{pair}.ionosphereCorrection.offset.tif` | Iono correction on offset grid (metres, noData = −2×10⁹) |
| `{pair}.ionosphereCorrection.offset.vrt` | Single-band offset-grid VRT |
| `simPhase/velSim.vrt` | Velocity-derived interferometric phase (rad), no topography; used for ambiguity correction |
| `simPhase/maskVel.vrt` | 0 where speed > velocityThreshold, 1 elsewhere; controls which pixels enter the iono fit |

When `--outputAll` is set, five bands are written to a single multi-band VRT
instead.

#### 3e. Power image collection (`processFramePow`)

Globs for `P<orbit1>_<frame>.*x*.pow` at the frame root.  If exactly one file
is found, its path is appended to `myArgs['pow']` and its companion
`geodat{R}x{A}.geojson` to `myArgs['geodatpow']` for later virtual-frame
assembly.  Runs regardless of `--RUNWOnly`.

### Step 4 — Virtual frame directory

Skipped entirely when `--correlationOnly` is set.

Creates `<orbit1>_<virtualFrame>/` and writes consolidated products.

#### 4a. Sensor YAML (`copy_sensor_yaml`)

Copies `NISAR{80|40|20}.yaml` from the `sarfunc` package into the virtual
frame directory as `sensor.NISAR{BW}.yaml`, and updates `intLooksR` /
`intLooksA` from the HDF5 look counts.  Skipped if the file already exists.

Bandwidth mapping:

| Bandwidth | YAML copied |
|-----------|-------------|
| ~77 MHz (rounds to 77) | `NISAR80.yaml` |
| ~40 MHz | `NISAR40.yaml` |
| ~20 MHz | `NISAR20.yaml` |

#### 4b. RUNW virtual frame (`createVirtualFrameRUNW`)

Calls `custom_buildvrtWithOffsets` for each product type, assembling per-frame
VRTs into a continuous mosaic that preserves inter-frame pixel offsets.

First, any per-frame `simPhase/velSim.vrt` and `simPhase/maskVel.vrt` files
are assembled into `velSim.vrt` and `maskVel.vrt` in the virtual frame
directory.  These are then passed to `custom_buildvrtWithOffsets` as
`--referencePhase` and `--mask` when building the `correctedUnwrappedPhase`
mosaic, so the DC bias of each frame's phase is anchored to the common velocity
simulation before joining.

Products assembled:

| Glob pattern searched per frame | Virtual product | `--offsets` used |
|---------------------------------|-----------------|-----------------|
| `*.correctedUnwrappedPhase.vrt` | `{pair_vf}.correctedUnwrappedPhase.vrt` | yes |
| `*.cor.vrt` | `{pair_vf}.cor.vrt` | no |
| `*.ionosphereCorrection.vrt` | `{pair_vf}.ionosphereCorrection.vrt` | yes |
| `*.ionosphereCorrection.offset.vrt` | `{pair_vf}.ionosphereCorrection.offset.vrt` | yes |

The path to the virtual `ionosphereCorrection.offset.vrt` is stored in
`myArgs['ionosphereRangeOffsetCorrection']` for use in the next step.

Then geodats are merged (`mergedGeodat`):
- Polygon corners: taken from first frame (near end) and last frame (far end),
  reversed for descending passes.
- Orbital state vectors: sorted, deduplicated, and cubic-spline interpolated
  onto a uniform time grid.
- Azimuth size updated to match the merged VRT dimensions.

#### 4c. Pair info file (`writePairInfo`)

Writes `<orbit1>.<orbit2>.pairinfo` into the virtual frame directory:

```
<orbit1>  <orbit2>  <refDateTime>  <secDateTime>  <nRangeLooks>  <nAzLooks>
```

#### 4d. ROFF virtual frame (`createVirtualFrameROFF`)

Skipped if `--RUNWOnly` is set.

Calls `custom_buildvrtWithOffsets --offsets` for each ROFF product type,
assembling them from the per-frame `offsetSims/` and frame-root locations:

| Source (per frame) | Virtual product |
|--------------------|----------------|
| `azimuth.offsets.vrt` | `azimuth.offsets.vrt` |
| `range.offsets.vrt` | `range.offsets.vrt` |
| `offsets.range-azimuth.vrt` | `offsets.range-azimuth.vrt` |
| `offsetSims/offsets.geom.vrt` | `offsets.geom.vrt` |
| `offsetSims/offsets.geom.ll.vrt` | `offsets.geom.ll.vrt` |
| `offsetSims/offsets.geom.mask.vrt` | `offsets.geom.mask.vrt` |
| `offsetSims/offsets.velocity.vrt` | `offsets.velocity.vrt` |
| `offsetSims/offsets.velocity.ll.vrt` | `offsets.velocity.ll.vrt` |
| `offsetSims/offsets.velocity.mask.vrt` | `offsets.velocity.mask.vrt` |

After assembly, the `ionosphereRangeOffsetCorrection` metadata key is written
onto the virtual `range.offsets.vrt` so the geocoding pipeline (`mosaic3d`)
knows to apply the ionosphere correction automatically.

#### 4e. Power virtual frame (`createVirtualFramePower`)

Skipped if `--RUNWOnly` is set or no `.pow` files were collected in step 3e.

Calls `custom_buildvrtWithOffsets.py` to assemble per-frame power images into
a virtual-frame power mosaic, then calls `mergedGeodat` on the power geodats.

---

## Two ionosphere paths

### Default: `estimateIonosphere` (recommended)

`RUNWtoGrimp` writes only coherence and geodats; the ionosphere and corrected
phase are produced by [estimateIonosphere](estimateIonosphere.md) which
combines the RUNW interferometric phase with the independent ROFF range offsets.
This approach is more accurate because it uses an independent measurement of the
ionospheric group delay rather than the NISAR-internal ionosphere estimate.

Requires: ROFF products to exist for the frame (cannot be used with
`--RUNWOnly`).

### Alternative: `--phaseDerivedIonosphere`

`RUNWtoGrimp` is passed `--phaseDerivedIonosphere` and handles the entire phase
pipeline internally, using the ionosphere phase screen embedded in the RUNW
HDF5.  This path can be used when ROFF products are unavailable.  See
[RUNWtoGrimp](RUNWtoGrimp.md) for details of the internal correction.

Note: the virtual-frame assembly in `createVirtualFrameRUNW` expects
`*.correctedUnwrappedPhase.vrt` from each frame; the phaseDerivedIonosphere
path produces `*.uw.interp.vrt` instead.  The virtual-frame VRT build will
report a warning and write a `.fail` sentinel for the
`correctedUnwrappedPhase` product.

---

## Skip and overwrite logic

| Condition | Behaviour |
|-----------|-----------|
| ROFF `range.offsets` present, no `--overWrite` | Skip ROFF frame |
| Phase geodat present in default path, no `--overWrite`/`--overWritePhase` | Skip RUNW frame |
| `*.correctedUnwrappedPhase.vrt` (or 5-band VRT) present, no `--overWrite`/`--overWritePhase` | Skip ionosphere frame |
| `--overWrite` | Re-run RUNW, ROFF, and ionosphere for all frames |
| `--overWritePhase` | Re-run RUNW and ionosphere only |
| `--correlationOnly` | Skip ROFF, ionosphere, and entire virtual-frame step; only RUNW (coherence/geodats) and power collection run |
| Mixed-mode frame, no `--allowMixedMode` | Skip RUNW silently; ROFF still runs |
| VRT build fails | Writes `<vrtFile>.fail` sentinel; any previous `.fail` is removed before each attempt |
| `simPhase/velSim.vrt` present, no `--overWrite` | Skip siminsar velSim regeneration |
| No `.pow` files collected | `createVirtualFramePower` returns immediately |

---

## noData convention

All final products use **−2×10⁹** as the noData sentinel to match what the C
geocoding programs (`mosaic3d`, `intfloat`, etc.) expect.  This applies to:

- `correctedUnwrappedPhase.tif` — cc=0 and phaseThresh-masked pixels
- `ionosphereCorrection.offset.tif` — pixels outside the offset VRT extent
- ROFF binary flat files — unfilled pixels from `intfloat`

Python-internal intermediate GeoTIFFs (e.g. `ionosphereCorrection.tif`) use
NaN noData and are not read by C programs.

---

## Dependencies

| Tool / Package | Called by | Role |
|----------------|-----------|------|
| [`RUNWtoGrimp`](RUNWtoGrimp.md) | `processFrameRUNW` | RUNW HDF5 → coherence, geodats, optionally phase |
| [`ROFFtoGrimp`](ROFFtoGrimp.md) | `processFrameROFF` | ROFF HDF5 → offset binaries and VRTs |
| [`estimateIonosphere`](estimateIonosphere.md) | `processFrameIonosphere` | Estimate and remove ionospheric phase; correct unwrapping ambiguities |
| `simoffsets` | `ROFFtoGrimp` (via `simulateOffsets`) | Simulate range/azimuth offset fields |
| `siminsar` | `estimateIonosphere` (via `run_vel_sim`, `run_mask_vel`) | Simulate interferometric phase and velocity mask |
| `intfloat` | `ROFFtoGrimp`, `estimateIonosphere` | Hole-filling interpolation |
| `cullst` | `ROFFtoGrimp` | Spatial outlier filtering of offsets |
| `custom_buildvrtWithOffsets` | `createVirtualFrameRUNW`, `createVirtualFrameROFF`, `createVirtualFramePower` | Build inter-frame VRT mosaics with pixel-offset alignment |
| `nisarhdf` | `getSecondaryOrbit`, `processFrameRUNW` | NISAR HDF5 reader |
| `sarfunc` | `copy_sensor_yaml` | Sensor YAML templates |
| `gdal` / `osgeo` | `createVirtualFrameROFF` | Write VRT metadata |
| `geojson`, `scipy` | `mergedGeodat` | Geodat GeoJSON merge and state-vector interpolation |
