# CLAUDE.md — nisargrimpworkflow

Converts NISAR Level-2 HDF5 products into GrIMP binary flat-file format for ingestion by the velocity mosaic pipeline. See the [packages CLAUDE.md](../CLAUDE.md) for the full pipeline context.

## Active project directory

Current NISAR Greenland development lives under:
```
/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/
```
This directory contains `project.yaml`, the track subdirectories, and a `templates/` directory with `tie_plan_header` and `vel_thumb_plan` templates. When searching for `project.yaml` or track data, look here first.

## Programs

Scripts in this package:

| Script | Entry point | Role |
|---|---|---|
| `ROFFtoGrimp` | `ROFFtoGrimp.py:main()` | ROFF HDF5 → GrIMP offsets |
| `RUNWtoGrimp` | `RUNWtoGrimp.py:main()` | RUNW HDF5 → GrIMP phase products |
| `SetupNISAR` | `SetupNISAR.py:main()` | Orchestrate per-frame conversion + virtual-frame assembly |
| `FileNISARProducts` | `FileNISARProducts.py:main()` | Organise raw HDF5 downloads into `track-{N}/source/` tree |
| `processTrack` | `processTrack.py:main()` | Run `SetupNISAR` for every orbit in a single track directory |
| `setupNISARTracks` | `setupNISARTracks.py:main()` | Initialise track dirs and refresh tie points across all tracks |
| `buildFrameGpkg` | `buildFrameGpkg.py:main()` | Per-cycle QC GeoPackages of virtual frames (see [doc](Documents/buildFrameGpkg.md)) |
| `buildFrameLayers` | `buildFrameLayers.py:main()` | QGIS `.qlr` layer tree for `buildFrameGpkg`'s output (see [doc](Documents/buildFrameLayers.md)) |
| `estimateIonosphere` | `estimateIonosphere.py:main()` | Offset-based ionosphere estimation + corrected phase (see [doc](Documents/estimateIonosphere.md)) |
| `custom_buildvrtWithOffsets` | `custom_buildvrtWithOffsets.py:main()` | gdalbuildvrt replacement that solves/applies per-frame DC offsets when merging |
| `checkFrameInputs` | `checkFrameInputs.py:main()` | Scan virtual frames for missing offsets/phase/iono products (see [doc](Documents/checkFrameInputs.md)) |
| `makeMaster` | `makeMaster.py:main()` | Assemble master input file lists (see [doc](Documents/makeMaster.md)) |
| `autoupdateNISAR` | `autoupdate.py:main()` | Cron driver: search/download + one-day completion buffer + full processing chain (see [doc](Documents/autoupdate.md)) |

---

## autoupdateNISAR

Daily cron driver. Searches ASF, downloads into a dated staging folder, holds each
day's data one day so late-arriving frames of the same `(cycle, track)` pass can be
merged in, then releases the completed passes and runs the chain:
`FileNISARProducts` → per-orbit `SetupNISAR` (no `--new`, so virtual frames extended
by late frames rebuild in place) → `setupNISARTracks --year` → `makeMaster`. Orbit is
derived from the filename via `orbitFromCycleTrack(cycle, track) = 173*cycle + track
+ 618` (in `FileNISARProducts.py`), so no HDF5 read is needed to find the
`track-<track>/<orbit1>_*` dirs. `--noDownload` runs the chain on already-staged
folders. Full detail in [autoupdate.md](Documents/autoupdate.md).

---

## ROFFtoGrimp

Converts a single NISAR ROFF (range/azimuth offset) HDF5 product into the GrIMP binary format consumed by `mosaic3d`.

### Processing steps (in order)

1. Opens ROFF HDF5 via `nisarhdf.nisarROFFHDF`
2. Discards offsets below correlation peak thresholds (per layer: default 0.07, 0.05, 0.025)
3. Writes `.dat` metadata files (`offsets.velocity.dat`, `offsets.geom.dat`) into `offsetSims/` for use by `simoffsets`
4. Calls `simoffsets` (GIT64 C binary, two threads); outputs land in `offsetSims/`:
   - geometry-only simulation (no velocity): `offsets.geom.*`
   - full simulation (geometry + velocity): `offsets.velocity.*`
5. Optionally applies a mask (`workingDir/offsets.velocity.mask.vrt`) from the simulation to fast-moving areas (layer 3)
6. Writes per-layer binary flat files to `workingDir/`:
   - `NISARoffsets.layer{N}.dr` — range offsets (big-endian float32)
   - `NISARoffsets.layer{N}.da` — azimuth offsets
   - `NISARoffsets.layer{N}.sr` — range sigma
   - `NISARoffsets.layer{N}.sa` — azimuth sigma
7. Calls `cullst` (GIT64 C binary) per layer (threaded) → `*.layer{N}.cull.{dr,da,sr,sa}`
8. Calls `intfloat` (GIT64 C binary) per layer/component (threaded) → `*.layer{N}.cull.interp.{dr,da,sr,sa}`
9. Merges three layers by nanmean, adds geometry back → `range.offsets`, `azimuth.offsets`, `range.offsets.sr`, `azimuth.offsets.sa` (all big-endian float32)
10. Writes final VRTs: `azimuth.offsets.vrt`, `range.offsets.vrt`, `offsets.range-azimuth.vrt`

### Output VRT metadata

The VRTs carry metadata that `mosaic3d` reads via GDAL:
`ByteOrder`, `geo1`, `geo2`, `r0`, `a0`, `deltaR`, `deltaA`, `sigmaStreaks`, `sigmaRange`, `correlationThresholds`, `region`, and cull parameters.

### Key arguments

```
ROFFtoGrimp [--outputDir DIR] [--noMask] [--verbose] [--mergeOnly]
            [--correlationThresholds T1 T2 T3]
            [--boxSize N] [--nGood N] [--maxR F] [--maxA F] [--sr N] [--sa N]
            [--interpThresh N] [--islandThresh N]
            [--geodat1 F] [--geodat2 F] [--DEM F] [--region R] [--regionFile YAML]
            [--verticalCorrection F] [--ompThreads N] [--byteOrder MSB|LSB]
            [--minTol F --percentSpeed F --maxTol F] [--maxSmoothRadius N]
            [--smoothNIter N] [--noVariableSmoothing] [--debugIono] ROFF_HDF5
```
`--mergeOnly` skips simulation/culling/interpolation and only re-runs the final merge step.
The `--minTol/--percentSpeed/--maxTol` trio (all three required together) enables the variable
smoothing-radius map applied on top of the fixed `--sr/--sa` smoothing.

---

## RUNWtoGrimp

Converts a single NISAR RUNW (unwrapped interferogram) HDF5 product into VRT files for the mosaic pipeline.

### Output files (in `outputDir/orbit1_frame/`)

**Default pipeline path** (`SetupNISAR` invokes it with `--noPhase --noIon`): only

- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.nisar.cor` (+ `.vrt`) — coherence
- `geodat{NLR}x{NLA}.geojson` / `geodat{NLR}x{NLA}.secondary.geojson`

The interpolated/corrected phase and ionosphere products are produced downstream by
`estimateIonosphere` (`*.correctedUnwrappedPhase.vrt`, `*.ionosphereCorrection*.vrt`).

**Legacy `--phaseDerivedIonosphere` path** additionally writes:

- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.nisar.uw.interp.vrt` — unwrapped phase (band description `Phase`)
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.nisar.ion.filt.rangeOffset.vrt` — filtered iono as range offset
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.nisar.ion.unfilt.rangeOffset.vrt` — unfiltered iono

The `radiansToPixels = −λ/(4π·slp)` scale here converts the ionosphere phase screen to a
*correction* (correction = −ionosphere), so like the `estimateIonosphere` products it is
ADDed by consumers — the opposite-looking sign vs `estimateIonosphere`'s `+λ/(4π·slp)` is
not a discrepancy (that scale applies to an already-negative iono estimate). See
`Documents/estimateIonosphere.md` "Background and equations".

---

## SetupNISAR

Orchestrates multi-frame conversion for a full orbit pass and assembles per-frame products into a virtual-frame mosaic.

### Directory structure assumed

```
{orbit1}_{frame}/          ← frame suffix unpadded (e.g. 1830_35)
    H5/
        NISAR*RUNW*.h5
        NISAR*ROFF*.h5
        NISAR*RIFG*.h5
```
Multiple frame directories for the same `orbit1`. (HDF5s at the frame-dir root
are still found as a fallback.)

### Processing flow

1. Discovers frame directories (`{orbit1}_{NN}`) matching the orbit number and splits them into contiguous groups (with secondary-epoch/bandwidth splitting)
2. Determines secondary orbit and bandwidth per group (majority vote across the group's frames)
3. For each frame: calls `processFrameRUNW` → `ROFFtoGrimp` → `processFrameROFF`, then `estimateIonosphere` per frame (unless `--phaseDerivedIonosphere`)
4. Copies `sensor.NISAR{bw}.yaml` into the virtual-frame directory and updates `intLooksR`/`intLooksA`
5. Calls `createVirtualFrameRUNW`:
   - Runs `custom_buildvrtWithOffsets.py` per product type: `correctedUnwrappedPhase`, `cor`, `ionosphereCorrection`, `ionosphereCorrection.offset` (with global fill: the `Unfilled` variants instead, filled later by `globalFillIonosphere()`; with `--phaseDerivedIonosphere`: the legacy `uw.interp` / `ion.*.rangeOffset` products)
   - Writes merged geodat GeoJSONs (merging corners and state vectors across frames)
6. Calls `createVirtualFrameROFF`:
   - Runs `custom_buildvrtWithOffsets.py` for all ROFF VRT types
   - Sets `ionosphereRangeOffsetCorrection` metadata on `range.offsets.vrt` so geocoding applies it
7. Writes a `.pairinfo` file: `orbit1 orbit2 date1 date2 NLR NLA`

### Virtual frame naming

`--virtualFrame` defaults to `None`: virtual-frame numbers are assigned automatically per
contiguous frame group (`assignVirtualFrameNumbers`), with `0000` the canonical full group and
higher suffixes for straggler/fragment groups. Individual frames are `orbit1_NN` (unpadded).
Pass `--virtualFrame VVVV` to force a specific suffix.

### Mixed mode

NISAR products where the SLC granule name contains `_M_` are mixed mode. By default these are skipped unless `--allowMixedMode` is set.

### Key arguments

```
SetupNISAR orbit1 [--virtualFrame VVVV] [--firstFrame N] [--lastFrame N]
           [--firstDate YYYY-MM-DD] [--lastDate YYYY-MM-DD]
           [--overWrite] [--overWritePhase] [--allowMixedMode]
           [--RUNWOnly] [--noMask] [--verbose] [--ompThreads N]
           [--phaseDerivedIonosphere] [--sepIceRock] [--noGlobalFillIono]
           [--retainIntermediateIono] [--debugIono]
           [--phaseThresh RAD] [--noPhaseThreshPass] [--outputAll]
           [--sigmaAz PX] [--sigmaRg PX]
           [--correlationOnly | --corrOnly] [--geodatsOnly] [--bakeOnly]
           [--new] [--clean] [--cleanDebug] [--noPrompt]
```
See `SetupNISAR --help` (or [Documents/SetupNISAR.md](Documents/SetupNISAR.md)) for details;
production runs use `--sepIceRock` (which also enables the global ionosphere fill).

---

## Geodat GeoJSON format

The `.geojson` files carry per-image geometry. Key `properties` fields:

- `Date`, `SecondaryDate` — ISO date strings
- `PassType` — `ascending` or `descending`
- `MLRangeSize`, `MLAzimuthSize` — image dimensions after multi-looking
- `NumberRangeLooks`, `NumberAzimuthLooks`
- `NumberOfStateVectors`, `TimeOfFirstStateVector`, `StateVectorInterval`
- `SV_Pos_N`, `SV_Vel_N` — state vectors (ECEF, metres and m/s)

When merging frames, corners are updated to span first-to-last, and state vectors are merged by sorting, deduplicating, and cubic interpolating onto a uniform time grid.

### Squint (residual Doppler) — extracted, merged, and applied as an opt-in `mosaic3d` correction

3-stage rollout (see `~/progs/GIT64/mosaicSource/CLAUDE.md` "Squint (residual Doppler)
sensitivity" for why the underlying analysis concluded no `mosaic3d` fix is currently needed, and
`~/PycharmProjects/packages/nisarErrors/Documents/plotSquintError.md` for the full derivation):
1. **Done** — `RUNWtoGrimp`/`nisarhdf` measures squint per RUNW sub-frame from the
   `geolocationGrid` cube and writes a 6-parameter polynomial fit into each per-frame geodat's
   `squintAnglePolynomial` key (see `nisarhdf/CLAUDE.md`). `None` for the secondary image.
2. **Done** — `SetupNISAR.mergedGeodat()` combines per-sub-frame fits into one polynomial for
   the virtual frame, via the new `mergeSquintAnglePolynomial(geos, geoMerged)` helper (defined
   just above `mergedGeodat()`, called right after the existing
   `MLNearRange`/`MLFarRange`/`MLCenterRange` recompute block). Each sub-frame's polynomial is
   only valid over its own local azimuth window, so the per-frame coefficients can't be averaged
   directly — instead each sub-frame's fitted surface is resampled onto a grid spanning the
   *merged* range bounds (shared swath geometry across sub-frames of one pass) and *that
   sub-frame's own* reconstructed azimuth bounds, every sub-frame's samples are pooled, and a
   single polynomial is refit over the combined domain (the same "redistribute onto one common
   reference" idea `mergeStateVectors()` already uses). A sub-frame's azimuth span isn't stored
   directly in the geodat schema but is exactly reconstructable from fields that are:
   `azimuthSpan = (MLAzimuthSize - 1) * (NumberAzimuthLooks / PRF)`, centered on
   `squintAnglePolynomial['refAzimuthTime']`. No-op (stays `None`) for the secondary image —
   detected via `geoMerged['properties']['squintAnglePolynomial'] is None`, no need to check the
   `secondary` flag explicitly. Verified on real data (track-25, orbit 3757, frames 51-55): merged
   `c0` = 1.4825° sits sensibly among the per-frame values (1.470°-1.501°, a smooth ~0.03°
   frame-to-frame drift matching the analysis doc's documented noise floor); re-evaluating the
   merged fit against each sub-frame's own fitted surface gives residual std ~0.0002°. Also
   writes flat top-level `squintCoefficients`/`squintRefRange`/`squintRefAzimuthTime` fields
   alongside the nested `squintAnglePolynomial` dict — needed because `mosaic3d`'s C-side reader
   (GDAL's OGR GeoJSON driver) can't read a nested object, only flat scalars/lists.
3. **Done** — `mosaic3d` consumes the merged polynomial, phase path only (`make3DMosaic.c`, not
   `make3DOffsets.c`), behind `-useSquint` (default off). See `mosaicSource/CLAUDE.md`'s squint
   section for the C-side implementation and verification (a real `mosaic3d` run on an actual
   overlapping ascending/descending pair confirmed both the rotation magnitude and sign).

**End-to-end plumbing** (opt-in, off by default): `project.yaml` key `applySquintCorrection`
(default `false`) → `setupNISARTracks --useSquint`/`--noUseSquint` (CLI overrides the project
default) → `refreshties.py --useSquint` → `makeframetie.py --useSquint` → `tie_script
--useSquint`, which appends `-useSquint` to the `mosaic3d` call and to both `tiepoints -motion`
call sites. See `mosaicworkflow/CLAUDE.md` and `insarworkflow/CLAUDE.md` for those packages'
squint entries.

## Bandwidth → sensor YAML mapping

| Bandwidth (MHz) | Sensor YAML | Band |
|---|---|---|
| ~77 (rounds to 77 int) | `NISAR80.yaml` | L-band 80 MHz |
| ~40 | `NISAR40.yaml` | L-band 40 MHz |
| ~20 | `NISAR20.yaml` | L-band 20 MHz |

The YAML sets `intLooksR` and `intLooksA` which C programs use for pixel-spacing calculations.

---

## FileNISARProducts

Organises a download directory of NISAR HDF5 products into the two-level track-based directory tree expected by `SetupNISAR`. Run once on a fresh download before calling `SetupNISAR`.

### Expected input structure

```
inputPath/
  RUNW/*.h5
  ROFF/*.h5
  RIFG/*.h5
  RSLC/*.h5   (optional)
```

### Usage

```
FileNISARProducts inputPath [--firstOrbit N] [--lastOrbit N]
                  [--outputPath DIR] [--reFile] [--verbose]
```

- `inputPath` — root directory with `RUNW/`, `ROFF/`, `RIFG/`, `RSLC/` subdirectories
- `--outputPath` — destination root; default is the current working directory
- `--firstOrbit` / `--lastOrbit` — filter on `referenceOrbit` from HDF5 metadata
- `--reFile` — bypass the already-filed check; use when `source/` was built by a prior run but `orbit_frame/` directories still need to be created
- `--verbose` — print a line per skipped file; otherwise only the total count is shown

### What it does (per RUNW)

1. Derives `track` from the RUNW filename (no HDF open) for the fast-path skip check
2. Creates `outputPath/track-{track}/` if absent
3. Creates `outputPath/track-{track}/source/` if absent; symlinks RUNW, ROFF, RIFG, and any matching RSLC there
4. Opens the RUNW HDF5 (`noLoadData=True`) to read `referenceOrbit`, `secondaryOrbit`, `frame`; applies orbit filter; skips mixed-mode frames
5. Creates `outputPath/track-{track}/{orbit1}_{frame}/`; symlinks RUNW, ROFF, RIFG there

### Output layout

```
outputPath/
  track-64/
    source/
      NISAR_L1_PR_RUNW_....h5   (symlink)
      NISAR_L1_PR_ROFF_....h5   (symlink)
      NISAR_L1_PR_RIFG_....h5   (symlink)
      NISAR_L1_PR_RSLC_....h5   (symlink, if present)
    12345_10/                    ← SetupNISAR reads from here (frame unpadded)
      H5/
        NISAR_L1_PR_RUNW_....h5 (symlink)
        NISAR_L1_PR_ROFF_....h5 (symlink)
        NISAR_L1_PR_RIFG_....h5 (symlink)
    12345_20/
      ...
    unfiled/                     ← duplicate products not selected
      NISAR_L1_PR_RUNW_....h5   (symlink)
      log
```

### Skip logic

The fast-path check tests whether `track-{track}/source/{RUNW_basename}` already exists (using filename parsing — no HDF open). If it does, the file is counted as skipped and the HDF is never opened. Pass `--reFile` to bypass this and always process the HDF (step 4–5 still runs, `symlink_file` silently skips any links that already exist).

### Duplicate product handling

Before the main loop, products are grouped by `(track, frame)`. When more than one RUNW maps to the same slot, `selectBestRUNW` picks the winner:

1. **Shortest temporal baseline** — `|date2Start − date1Start|` in days
2. **Newest modification time** — tiebreaker; proxy for processing version

Losers are symlinked into `track-{N}/unfiled/` and a log entry is appended to `track-{N}/unfiled/log` explaining why each was not selected and which product won. Log entries are only written for newly created symlinks, so re-runs do not produce duplicate log lines.

```
track-64/
  unfiled/
    NISAR_L1_PR_RUNW_....h5   (symlink, 24-day pair)
    log                        ← human-readable reason + winner name
```

### Helper functions

**`getCompanion(RUNW, inputPath, productType)`** — returns the path to the matching ROFF or RIFG file (same filename with product type substituted, in `inputPath/{productType}/`), or `None` if absent.

**`findRSLC(RUNW, inputPath)`** — globs `inputPath/RSLC/` for RSLC files matching the track, direction, and frame extracted from the RUNW basename (using the raw zero-padded fields, not the stripped `parseFileName` values). Returns a list (may be empty).

**`selectBestRUNW(candidates)`** — given a list of RUNW paths for the same track/frame, returns `(winner, losers, reasons)`. `reasons` is a dict mapping each loser path to a human-readable string explaining the selection.

**`parseFileName(product)`** — splits the NISAR filename on `_` and maps fields by position. Two layouts:
- RSLC (13 fields): `NISAR_L1_PR_RSLC_{cycle}_{track}_{direction}_{frame}_{bw}_{pol}_{mode}_{date1Start}_{date1End}`
- All others (15 fields): `NISAR_L1_PR_{type}_{cycle}_{track}_{direction}_{frame}_004_{bw}_{pol}_{date1Start}_{date1End}_{date2Start}_{date2End}`

---

## buildFrameGpkg / buildFrameLayers

A read-only QC/inventory pair, independent of the conversion pipeline above —
they never write into `track-N/`. Full detail in
[buildFrameGpkg.md](Documents/buildFrameGpkg.md) and
[buildFrameLayers.md](Documents/buildFrameLayers.md); summary:

1. **`buildFrameGpkg`** walks `track-N/*_0000/` virtual frames under
   `--projectDir` (default `.`) and writes one GeoPackage per cycle
   (`<outputDir>/cycle<NN>.gpkg`, layers `ascending`/`descending`), pulling
   baseline/rBaseline/azimuth sigma + tiepoint counts, frame list, and
   footprint geometry from the sidecar files `SetupNISAR`/tie-point
   processing already produced. Cycle/direction come from `ImageName` via
   the same `parseFileName` used by `FileNISARProducts` above. Any frame
   missing a required input is skipped with a printed reason, not silently
   patched over.
2. **`buildFrameLayers`** reads that GeoPackage directory and writes a QGIS
   `.qlr` with four groups: `Frames` (sigma field switchable via a QGIS
   project variable) and `rBaseline` (fixed 10-class `sigmaRBaseline`
   coloring), both organized `ascending`/`descending` → `Cycle N`; plus two
   flat groups with just `ascending`/`descending` — no per-cycle split —
   built from an OGR VRT Union Layer across all cycles' GeoPackages plus a
   subset filter: `sigmaRBaseline > <thresh>` (`--offsetsSigmaThresh`,
   default `0.5`) and `Bad` (frames with `sigmaRBaseline = -1 AND
   sigmaRBaselineWithoutIon = -1`, i.e. no baseline solution).

Both are general — they work on any directory laid out as `track-N/*_0000/`,
not just one specific project.

`track`, `frame`, `cycle` returned as `str(int(...))` (leading zeros stripped); date fields as `datetime` objects.
