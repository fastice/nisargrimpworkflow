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

---

## ROFFtoGrimp

Converts a single NISAR ROFF (range/azimuth offset) HDF5 product into the GrIMP binary format consumed by `mosaic3d`.

### Processing steps (in order)

1. Opens ROFF HDF5 via `nisarhdf.nisarROFFHDF`
2. Discards offsets below correlation peak thresholds (per layer: default 0.07, 0.05, 0.025)
3. Writes `.dat` metadata files (`offsets.dat`, `offsets.geom.dat`) for use by `simoffsets`
4. Calls `simoffsets` (GIT64 C binary, two threads):
   - geometry-only simulation (no velocity): `offsets.geom.*`
   - full simulation (geometry + velocity): `offsets.*`
5. Optionally applies a mask (`offsets.mask.vrt`) from the simulation to fast-moving areas (layer 3)
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
            [--interpThresh N] [--islandThresh N] ROFF_HDF5
```
`--mergeOnly` skips simulation/culling/interpolation and only re-runs the final merge step.

---

## RUNWtoGrimp

Converts a single NISAR RUNW (unwrapped interferogram) HDF5 product into VRT files for the mosaic pipeline.

### Output files (in `outputDir/orbit1_frame/`)

- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.uw.interp.vrt` — unwrapped phase
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.cor.vrt` — coherence
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.ion.filt.vrt` — ionosphere correction (filtered)
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.ion.filt.rangeOffset.vrt` — iono as range offset
- `{orbit1}_{frame}.{orbit2}_{frame}.{NLR}x{NLA}.ion.unfilt.rangeOffset.vrt` — unfiltered iono

Also writes `geodat{NLR}x{NLA}.geojson` and `geodat{NLR}x{NLA}.secondary.geojson` for both reference and secondary orbits.

---

## SetupNISAR

Orchestrates multi-frame conversion for a full orbit pass and assembles per-frame products into a virtual-frame mosaic.

### Directory structure assumed

```
{orbit1}_{frame}/
    NISAR*RUNW*.h5
    NISAR*ROFF*.h5
```
Multiple frame directories for the same `orbit1`.

### Processing flow

1. Discovers frame directories (`{orbit1}_{NN}`) matching the orbit number
2. Determines secondary orbit and bandwidth from the first RUNW/ROFF found
3. For each frame: calls `processFrameRUNW` → `ROFFtoGrimp` → `processFrameROFF`
4. Copies `sensor.NISAR{bw}.yaml` into the virtual-frame directory and updates `intLooksR`/`intLooksA`
5. Calls `createVirtualFrameRUNW`:
   - Runs `custom_buildvrtWithOffsets.py` for each product type: `uw.interp`, `cor`, `ion.filt`, `ion.filt.rangeOffset`, `ion.unfilt.rangeOffset`
   - Writes merged geodat GeoJSONs (merging corners and state vectors across frames)
6. Calls `createVirtualFrameROFF`:
   - Runs `custom_buildvrtWithOffsets.py` for all ROFF VRT types
   - Sets `ionosphereRangeOffsetCorrection` metadata on `range.offsets.vrt` so geocoding applies it
7. Writes a `.pairinfo` file: `orbit1 orbit2 date1 date2 NLR NLA`

### Virtual frame naming

Frame `0000` (default `--virtualFrame`) is the virtual merged product. Individual frames are `orbit1_NN`. A virtual-frame directory `orbit1_0000/` holds all merged VRTs.

### Mixed mode

NISAR products where the SLC granule name contains `_M_` are mixed mode. By default these are skipped unless `--allowMixedMode` is set.

### Key arguments

```
SetupNISAR orbit1 [--virtualFrame VVVV] [--firstFrame N] [--lastFrame N]
           [--overWrite] [--overWritePhase] [--allowMixedMode]
           [--RUNWOnly] [--noMask] [--verbose]
```

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
    12345_010/                   ← SetupNISAR reads from here
      NISAR_L1_PR_RUNW_....h5   (symlink)
      NISAR_L1_PR_ROFF_....h5   (symlink)
      NISAR_L1_PR_RIFG_....h5   (symlink)
    12345_020/
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

`track`, `frame`, `cycle` returned as `str(int(...))` (leading zeros stripped); date fields as `datetime` objects.
