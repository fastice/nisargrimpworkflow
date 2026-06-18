# nisargrimpworkflow — Processing Workflow Overview

Converts NISAR Level-2 HDF5 products into GrIMP binary/VRT format and assembles them into virtual-frame mosaics ready for `mosaic3d` velocity estimation.

---

## Directory structure

NLR = NumberRangeLooks, NLA = NumberAzimuthLooks (e.g. 26x16).
orbit1/orbit2 are reference/secondary orbit numbers; NN is the frame number.

```
project/                                   ← PROJECT_DIR (set in setupNISARTracks.py)
  project.yaml                             ← project-level configuration (hand-edited)
  templates/                               ← master templates copied into each track
    tie_plan_header                        ← tiepoint plan header template
    vel_thumb_plan                         ← velocity thumbnail plan template

  track-1/
    source/                                ← all HDF5 symlinks (FileNISARProducts)
      NISAR*RUNW*.h5
      NISAR*ROFF*.h5
      NISAR*RIFG*.h5
      NISAR*RSLC*.h5  (if present)

    orbit1_NN/                             ← one dir per frame (FileNISARProducts)
      H5/                                  ← HDF5 symlinks + per-layer VRT wrappers
        NISAR*RUNW*.h5
        NISAR*RUNW*.vrt
        NISAR*ROFF*.h5
        NISAR*ROFF*.layer{1,2,3}.vrt
        NISAR*RIFG*.h5
        NISAR*RIFG*.vrt
      workingDir/                          ← per-layer offset intermediates (ROFFtoGrimp)
        NISARoffsets.layer{N}.{dr,da,sr,sa}           ← raw binary offsets
        NISARoffsets.layer{N}.{cc,mt,mt.vrt}          ← correlation/mask
        NISARoffsets.layer{N}.cull.{dr,da,sr,sa,vrt}  ← after cullst
        NISARoffsets.layer{N}.cull.interp.{dr,da,sr,sa,vrt}  ← after intfloat
      offsetSims/                          ← simoffsets output
        offsets.geom.{dr,da,lat,lon,mask,poly,simdat,dat}
        offsets.geom.{vrt,ll.vrt,mask.vrt}
        offsets.velocity.{dr,da,lat,lon,mask,dat}
        offsets.velocity.{vrt,ll.vrt,mask.vrt}
        geodat{NLR}x{NLA}.{geojson,secondary.geojson}
      simPhase/                            ← RUNWtoGrimp simulation output
        velSim  velSim.simdat  velSim.vrt
        maskVel  maskVel.simdat  maskVel.vrt
      motion/                              ← baseline/velocity estimation (refreshties)
        estR  estR_stderr
        az.est  azest  az.est.const  azest_stderr  az.est.svlinear
        rBaseline  rBaseline.deltabp  rBaseline.quad
        baseline.{NLR}x{NLA}
      range.offsets  range.offsets.sr      ← merged offset flat files (big-endian float32)
      azimuth.offsets  azimuth.offsets.sa
      range.offsets.vrt
      azimuth.offsets.vrt
      offsets.range-azimuth.vrt
      orbit1_NN.orbit2_NN.NLRxNLA.nisar.cor          ← coherence flat file
      orbit1_NN.orbit2_NN.NLRxNLA.nisar.cor.vrt
      orbit1_NN.orbit2_NN.NLRxNLA.nisar.correctedUnwrappedPhase.{tif,vrt}
      orbit1_NN.orbit2_NN.NLRxNLA.nisar.ionosphereCorrection.{tif,vrt}
      orbit1_NN.orbit2_NN.ionosphereCorrection.offset.{tif,vrt}
      geodat{NLR}x{NLA}.geojson
      geodat{NLR}x{NLA}.secondary.geojson
      orbit1.orbit2.pairinfo

    orbit1_0000/                           ← virtual frame (SetupNISAR, merged across frames)
      frames.txt                           ← list of real frames merged here
      sensor.NISAR{bw}.yaml               ← sensor parameters (bandwidth-derived)
      geodat{NLR}x{NLA}.geojson
      geodat{NLR}x{NLA}.secondary.geojson
      orbit1.orbit2.pairinfo
      orbit1_0000.orbit2_0000.NLRxNLA.nisar.cor.vrt
      orbit1_0000.orbit2_0000.NLRxNLA.nisar.correctedUnwrappedPhase.vrt
      orbit1_0000.orbit2_0000.NLRxNLA.nisar.ionosphereCorrection.vrt
      orbit1_0000.orbit2_0000.ionosphereCorrection.offset.vrt
      range.offsets.vrt
      azimuth.offsets.vrt
      offsets.range-azimuth.vrt
      offsets.geom.{vrt,ll.vrt,mask.vrt}
      offsets.velocity.{vrt,ll.vrt,mask.vrt}
      velSim.vrt
      maskVel.vrt
      phase.uw.vrt
      *.vrt.stats                          ← per-VRT statistics sidecar files
      motion/                              ← baseline/velocity estimation (refreshties)
        (same contents as per-frame motion/)

    tiepoints/                             ← created by setupNISARTracks
      tie_plan_header
      vel_thumb_plan
      tie_plan{YYYY}                       ← per-year tiepoint plans
      tie_plan{YYYY}.{VF1}dash{VF2}
      tie_planAll
      refreshTies/                         ← refreshties.py working dir
      logs/
      old/

    velocityStats/
      {VF1}-{VF2}/                         ← one dir per virtual-frame pair
```

---

## Step 1 — Create project and region configuration files

Before any processing, create `project.yaml`, the region YAML file, and the `templates/` directory by hand in the project root.

### `project.yaml`

The top-level project configuration file. Read by `setupNISARTracks`, `refreshties`, `makevelstatsregions`, and `autoclean` to locate data paths and control processing options. Also read by `SetupNISAR` for `regionFile`/`verticalCorrection` and the variable smoothing-radius map params below, which it passes through to `ROFFtoGrimp`/`RUNWtoGrimp` as explicit CLI flags (see their `--minTol`/`--percentSpeed`/`--maxTol`/`--maxSmoothRadius`/`--smoothNIter`/`--noVariableSmoothing`).

```yaml
sensor: NISAR                    # sensor type; controls which code paths refreshties uses
region: /path/to/project/greenland  # path to the region yaml (without .yaml extension)
maxDays: 37                      # maximum temporal baseline (days) for pair selection
header: tiepoints/tie_plan_header  # relative path to the per-track tiepoint plan header
framePattern: '00??'             # glob matching virtual-frame directory suffixes
tiepointFile: /path/to/SentinelTies.vz  # base tiepoint file (GrIMP format)
tie_plan_header_template: /path/to/project/templates/tie_plan_header  # absolute path; default: PROJECT_DIR/templates/tie_plan_header
vel_thumb_plan_template:  /path/to/project/templates/vel_thumb_plan   # absolute path; default: PROJECT_DIR/templates/vel_thumb_plan
vel_thumb_header_template: /path/to/project/templates/vel_thumb_header # absolute path; default: PROJECT_DIR/templates/vel_thumb_header
velocityStatsMode: 'RA'          # 'RA' or 'XY' — velocity component pair for autoclean
velThumbOutput: tiff             # 'tiff' or 'binary' (default) — output format for vel_thumbs
# Variable smoothing-radius map (SetupNISAR -> ROFFtoGrimp/RUNWtoGrimp), an additional
# smoothing pass on top of the fixed -sr/-sa (offsets) / interpolation (phase) above.
# minTol/percentSpeed/maxTol must be given together to enable it; omit all three to leave
# it off (default). Per-pixel tolerance = clip(percentSpeed/100*speed, minTol, maxTol) m/yr.
minTol: 0.25                     # m/yr floor [no default -- feature off unless all three set]
percentSpeed: 1                  # percent of local speed, e.g. 1 = 1% [no default]
maxTol: 10                       # m/yr ceiling [no default]
maxSmoothRadius: 50              # sweep cap, single-look pixels, clamped to <= 255 [50]
smoothNIter: 3                   # repeated box-filter passes per sweep step (Gaussian-ish) [3]
noVariableSmoothing: false       # true disables the map even if minTol/percentSpeed/maxTol are set [false]
```

### Region YAML file (e.g. `greenland.yaml`)

Defines the geospatial parameters for a region. Referenced by the `region` key in `project.yaml`
(the path is stored without the `.yaml` extension; the `.yaml` is appended when loading).
Read by `defaultRegionDefs` in `sarfunc`, which is used by `makevelstatsregions`,
`velocityStats`, and `autoclean`.

Built-in region definitions for `greenland`, `antarctica`, `amundsen`, `columbia`, and `taku`
are bundled with `sarfunc` and can be used directly by name (e.g. `-region greenland`). For a
NISAR project, create a custom copy in the project root so paths can be tailored:

```yaml
epsg: 3413                        # EPSG code for the polar stereographic projection
wktFile: null                     # path to a WKT projection file, or null to derive from epsg
dem: /path/to/dem/dem.270m        # DEM used by tiepoints and mosaic3d
velMap: /path/to/velocityMap.vrt  # reference velocity map for statistics
fastmask: /path/to/greenlandmask.fast   # fast-moving ice mask
mask: /path/to/greenlandmask            # full ice mask
shelfMask: null                         # ice-shelf mask (null if not applicable)
sigmaShape: /path/to/sigmaScale.shp     # shapefile of per-region sigma scale factors
icemask: /path/to/GimpIceMask_90m       # high-resolution ice mask
```

Key fields:
- `epsg` — projection code; 3413 = North Polar Stereographic, 3031 = South Polar Stereographic
- `dem` — elevation model used for tiepoint coordinate conversion and mosaic geometry
- `velMap` — reference velocity mosaic used to seed statistics and detect outliers
- `fastmask` / `mask` / `icemask` — mask layers controlling which pixels are processed
- `sigmaShape` — shapefile that spatially varies the outlier rejection threshold in autoclean
- `shelfMask` — identifies floating ice shelves for tidal correction (null for Greenland)

### `templates/tie_plan_header`

Copied into each track's `tiepoints/tie_plan_header` by `setupNISARTracks --copyFiles`,
with `<TRACK>` substituted by the actual track number. Controls how `makeframetie` finds
tiepoints. Example:

```
extraties = <TIEFILE>
DEM       = <DEM>
track_root= /path/to/project/
tie_dir   = /path/to/project/track-<TRACK>/tiepoints
default_nDays=12
extra_flags = "-tieThresh 500"
```

- `extraties` — additional tiepoint file (VZ format) merged with the base tiepoints
- `DEM` — DEM path substituted at runtime by refreshties
- `track_root` — absolute path to the project root
- `tie_dir` — absolute path to this track's tiepoints directory
- `default_nDays` — default temporal baseline limit for tiepoint matching
- `extra_flags` — extra flags passed to the tiepoint program

### `templates/vel_thumb_plan`

Copied into each track's `tiepoints/vel_thumb_plan` by `setupNISARTracks --copyFiles`,
with `<TRACK>` substituted. Controls velocity thumbnail generation and mosaic assembly.
Example:

```
DEM       = <DEM>
track_root= /path/to/project/track-<TRACK>
tie_dir   = /path/to/project/track-<TRACK>/tiepoints
vel_dir   = /path/to/project/track-<TRACK>/velocity

extra_flags = "-no3d"
resolution  = "0 0 0 0     0.200000     0.200000"

start
track-<TRACK> all2025
track-<TRACK> all2026
end
```

- `resolution` — output pixel spacing in km (range, azimuth; zeros = native)
- `extra_flags` — extra flags for the velocity mosaic program (`-no3d` = 2-D output only)
- `start`/`end` block — list of track/year combinations to include in the mosaic

---

## Step 2 — Organise downloads: `FileNISARProducts`

Run once after downloading a batch of NISAR HDF5 files.

```
FileNISARProducts <download_dir> [--outputPath <project_dir>]
```

- Reads `<download_dir>/RUNW/`, `ROFF/`, `RIFG/`, `RSLC/`
- Creates `track-{N}/source/` (all products) and `track-{N}/{orbit1}_{frame}/` (selected products)
- Symlinks HDF5 files; does not copy
- Selects the best RUNW when multiple cover the same frame (shortest baseline, then newest)
- Losers go to `track-{N}/unfiled/` with an explanatory log entry
- Run with `--reFile` to redo filing without rebuilding `source/`

---

## Step 3 — Convert HDF5 to GrIMP format

### Option A — Single track: `processTrack` (recommended)

Run from anywhere; pass the track directory as the argument.

```
processTrack <track-dir> [--overWrite] [--overWritePhase] [--RUNWOnly]
```

Discovers all `orbit1` values in `<track-dir>` and calls `SetupNISAR <orbit1>` for each.

### Option B — Single orbit pair: `SetupNISAR`

Run from inside the track directory.

```
cd project/track-12
SetupNISAR 3744 [options]
```

#### What SetupNISAR does (per orbit pair)

For each frame directory (`3744_125`, `3744_126`, …):

1. **`RUNWtoGrimp`** — RUNW HDF5 → VRTs:
   - `*.uw.interp.vrt` — unwrapped phase
   - `*.cor.vrt` — coherence
   - `*.ion.filt.vrt`, `*.ion.filt.rangeOffset.vrt` — ionosphere correction
   - `geodat{NLR}x{NLA}.geojson` — image geometry for reference and secondary

2. **`ROFFtoGrimp`** — ROFF HDF5 → GrIMP offset files:
   - Simulates geometry (`simoffsets`), culls (`cullst`), interpolates (`intfloat`)
   - `range.offsets`, `azimuth.offsets` (big-endian float32)
   - `range.offsets.vrt`, `azimuth.offsets.vrt`

3. **Ionosphere** — estimates ionosphere from range offsets (`estimateIonosphere`) unless `--phaseDerivedIonosphere` was set (in which case RUNW step already wrote the correction)

4. **Virtual frame assembly** — merges all frames into `{orbit1}_0000/`:
   - One merged VRT per product type (uw.interp, cor, ion, offsets, …)
   - Merged geodat GeoJSON (corners and state vectors spanning all frames)
   - `pairinfo` file

#### Common modes

| Goal | Command |
|---|---|
| Full pipeline (default) | `SetupNISAR 3744` |
| Phase + coherence only, no offsets | `SetupNISAR 3744 --RUNWOnly` |
| Correlation/coherence only (QC pass) | `SetupNISAR 3744 --corrOnly` |
| Phase-derived ionosphere (skip range-offset iono) | `SetupNISAR 3744 --phaseDerivedIonosphere` |
| Force full redo | `SetupNISAR 3744 --overWrite` |
| Redo phase products only, keep offsets | `SetupNISAR 3744 --overWritePhase` |
| Limit to a frame range | `SetupNISAR 3744 --firstFrame 125 --lastFrame 130` |

---

## Step 4 — Tiepoints and velocity: `setupNISARTracks`

Run from the project root directory (the one containing `project.yaml` and all `track-*` subdirectories).

```
setupNISARTracks [--tracks track-N ...] [--year YYYY [YYYY ...]] [--tiesOnly]
                 [--overWrite] [--keepVz] [--noPhase] [--quadFit] [--runVelThumbs]
                 [--runVelstatsregions] [--check]
```

- Ensures `tiepoints/` exists and creates `velocityStats/{X}-{Y}/` dirs for each virtual-frame group found (based on `framePattern` from `project.yaml`)
- Creates `tiepoints/vel_thumb_header_XdashY` for each existing `velocityStats/X-Y/` dir (skip if already present; template from `vel_thumb_header_template` or `templates/vel_thumb_header`)
- Calls `refreshties.py --phaseAndOffsets --noQuadFit` across all tracks for the requested years
  (default years: 2025 2026)
- `--tracks track-N [track-N ...]` — restrict to specific track directories
- `--tiesOnly` — run `refreshties` without computing velocity mosaics
- `--overWrite` — rerun even if velocity products already exist
- `--keepVz` — retain `.vz` and `.vz.geodat` files
- `--noPhase` — disable phase+offsets mode (default: `--phaseAndOffsets` is passed)
- `--quadFit` — enable the `-deltaBQ` quadratic baseline correction estimate (default: `--noQuadFit` is passed)
- `--runVelThumbs` — run `vel_thumbs vel_thumb_plan` in each track's `tiepoints/` dir
- `--runVelstatsregions` — run `makevelstatsregions.py` in each track dir, then exit
- `--check` — report which virtual-frame dirs have `range.offsets.vrt` but no `velocity/mosaicOffsets.vx`, and print range/azimuth baseline sigmas sorted by range sigma

### First-time setup

```
setupNISARTracks --copyFiles
```

Copies `tie_plan_header` and `vel_thumb_plan` from the templates directory into every track's `tiepoints/` directory, substituting `<TRACK>` and `<DEM>` placeholders. Template paths default to `PROJECT_DIR/templates/` and can be overridden via `project.yaml`. Run once when initialising a new project.

---

## Typical end-to-end sequence

```bash
# 1. Create project.yaml, greenland.yaml (or other region yaml), and templates/ by hand

# 2. File the downloads
FileNISARProducts /downloads/batch1 --outputPath /project

# 3. Convert all orbit pairs in a track
processTrack /project/track-12

# 4. (Repeat step 3 for other tracks, or parallelise manually)

# 5. First-time only: copy tiepoint templates
setupNISARTracks --copyFiles

# 6. Run tiepoints and velocity estimation
setupNISARTracks --year 2025 2026
```
