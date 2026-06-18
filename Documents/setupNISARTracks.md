# setupNISARTracks — Multi-Track Directory Setup and Tie-Point Refresh

## Overview

`setupNISARTracks` initialises and maintains all `track-*` directories under a
project root.  It creates required subdirectories, optionally distributes
template files into every track, and then drives `refreshties.py` across all
tracks.  Auxiliary modes run velocity thumbnails and velocity-stats regions.

Must be run from the project root directory (the one containing `project.yaml`
and all `track-*` subdirectories).  The project directory is derived from the
current working directory.

---

## Usage

```
setupNISARTracks [options]
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--tracks track-N [track-N ...]` | all `track-*` dirs | Restrict processing to the listed track directories (e.g. `--tracks track-12 track-64`) |
| `--copyFiles` | off | Copy `tie_plan_header` and `vel_thumb_plan` from `track-88/tiepoints/` to every other track, substituting the track number |
| `--tiesOnly` | off | Run `refreshties.py` without computing velocity mosaics |
| `--overWrite` | off | Pass `--overWrite` to `refreshties.py` to rerun existing products |
| `--keepVz` | off | Pass `--keepVz` to `refreshties.py` to retain `.vz` and `.vz.geodat` files |
| `--runVelThumbs` | off | Run `vel_thumbs vel_thumb_plan` in every track's `tiepoints/` directory (10 threads) |
| `--runVelstatsregions` | off | Run `makevelstatsregions.py` in every track directory, then exit |
| `--year YYYY [YYYY ...]` | 2025 2026 | Year(s) forwarded to `refreshties.py` |
| `--noPhase` | off | Disable phase+offsets mode (by default NISAR processing uses `--phaseAndOffsets`) |
| `--quadFit` | off | Enable the `-deltaBQ` quadratic baseline correction estimate (by default NISAR processing uses `--noQuadFit`) |
| `--check` | off | Report product completeness and residual sigmas (see [Check mode](#check-mode)); exits without running any other steps |

### Examples

```bash
# First-time setup: create subdirectories and copy template files from track-88
setupNISARTracks --copyFiles

# Refresh tie points for all tracks, years 2025 and 2026
setupNISARTracks

# Refresh tie-point step only, overwriting existing results
setupNISARTracks --tiesOnly --overWrite

# Regenerate velocity thumbnails in all tiepoints directories
setupNISARTracks --runVelThumbs

# Check which track directories are missing velocity products
setupNISARTracks --check
```

---

## Processing steps (default mode)

### Step 1 — Directory creation

For every `track-*` directory found under `PROJECT_DIR`:

- Creates `<track>/tiepoints/` if it does not exist.
- Creates `velocityStats/{X}-{Y}/` directories for each virtual-frame group found
  (derived from `*_{framePattern}` directories and `framePattern` in `project.yaml`).

### Step 2 — Template file distribution (`--copyFiles`)

Reads `PROJECT_DIR/project.yaml` and instantiates two template files into
every track's `tiepoints/` directory.  Template paths come from `project.yaml`;
if a key is absent the corresponding file under `PROJECT_DIR/templates/` is used
as the default.

| `project.yaml` key | Default path | Destination in each tiepoints/ |
|--------------------|-------------|----------------------------------|
| `tie_plan_header_template` | `templates/tie_plan_header` | `tie_plan_header` |
| `vel_thumb_plan_template` | `templates/vel_thumb_plan` | `vel_thumb_plan` |

The placeholders substituted in each template:

| Placeholder | Value |
|-------------|-------|
| `<TRACK>` | Track number (e.g. `1` from `track-1`) |
| `<DEM>` | `dem` field from the region YAML (`region` / `regionFile` in `project.yaml`) |

Files already present are skipped.

### Step 3 — `vel_thumb_header` distribution (always)

After creating any missing `velocityStats/` subdirectories, `setupNISARTracks`
scans all `velocityStats/*-*` directories in the track and creates a
`vel_thumb_header_XdashY` file in `tiepoints/` for each one that does not
already have one.

The template source is `vel_thumb_header_template` from `project.yaml`, falling
back to `templates/vel_thumb_header`.  The same `<TRACK>` and `<DEM>`
substitutions are applied.  Existing files are **never overwritten** — they may
be manually edited later in the processing workflow.

### Step 4 — Velocity thumbnails (`--runVelThumbs`)

Runs `csh -c 'vel_thumbs vel_thumb_plan'` in each track's `tiepoints/`
directory using a 10-thread pool.

### Step 5 — Velocity stats regions (`--runVelstatsregions`)

Runs `makevelstatsregions.py` sequentially in each track directory, then
exits (skipping the `refreshties.py` call).

### Step 6 — Tie-point refresh

Calls `refreshties.py` once for all tracks:

```
refreshties.py [-tiesOnly] [--overWrite] [--keepVz] -toRun="[track-NN,...]" YYYY [YYYY ...] -noPrompt
```

---

## What each flag sets in motion

`main()` (in `setupNISARTracks.py`) is the actual control flow; the tree below follows
it exactly. Most flags don't change *which* programs run — they just change the arguments on
the same `refreshties.py` call — so only the structurally distinct cases get their own figure.
`--check` and `--runVelstatsregions` are the only two that skip `refreshties.py` entirely.

### No flags (default)

```
setupNISARTracks
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    │   ├── create <track>/tiepoints/                       (if missing)
    │   ├── create <track>/velocityStats/{X}-{Y}/           (if missing)
    │   └── write <track>/tiepoints/vel_thumb_header_XdashY (if missing; never overwritten)
    └── run_refresh_ties()
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--copyFiles`

```
setupNISARTracks --copyFiles
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=True)
    │   ├── create <track>/tiepoints/                       (if missing)
    │   ├── copy tie_plan_header → <track>/tiepoints/       (if missing; <TRACK>/<DEM> substituted)
    │   ├── copy vel_thumb_plan  → <track>/tiepoints/       (if missing; <TRACK>/<DEM> substituted)
    │   ├── create <track>/velocityStats/{X}-{Y}/           (if missing)
    │   └── write <track>/tiepoints/vel_thumb_header_XdashY (if missing)
    └── run_refresh_ties()
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--tiesOnly`

```
setupNISARTracks --tiesOnly
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)        (same as default)
    └── run_refresh_ties()
        └── refreshties.py -tiesOnly --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--overWrite`

```
setupNISARTracks --overWrite
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    └── run_refresh_ties()
        └── refreshties.py --overWrite --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--keepVz`

```
setupNISARTracks --keepVz
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    └── run_refresh_ties()
        └── refreshties.py --keepVz --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--runVelThumbs`

Adds a call but does **not** skip the tie-point refresh — both run.

```
setupNISARTracks --runVelThumbs
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── run_vel_thumbs(track_dirs)                    (10-thread pool, one thread/track)
    │   └── csh -c 'vel_thumbs vel_thumb_plan'         (cwd = <track>/tiepoints/)
    └── run_refresh_ties()                             ← still runs, no early return
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--runVelstatsregions`

The one flag (besides `--check`) that **skips** `refreshties.py`.

```
setupNISARTracks --runVelstatsregions
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    └── run_velstats_regions(track_dirs)
        └── csh -c 'makevelstatsregions.py'            (cwd = <track>/, one at a time)
    ⇥ return                                            (refreshties.py never called)
```

### `--check`

A fully separate branch — exits before any directory setup or track processing.

```
setupNISARTracks --check
└── main()
    └── check_products()
        ├── glob track-*/{framePattern}/range.offsets.vrt
        ├── report dirs missing velocity/mosaicOffsets.vx
        └── print range/azimuth sigma table from
              motion/{rBaseline.deltabp, rBaseline.deltabp.noIonosphere, az.est.const}
    ⇥ return                                            (no setup_track_dirs, no refreshties.py)
```

### Flags that only change the `refreshties.py` arguments

These don't alter the call tree — they apply to whichever figure above matches the rest of
your flags:

| Flag | Effect on the `refreshties.py` command line |
|------|------|
| `--year YYYY [YYYY ...]` | Replaces the trailing `2025 2026` with the year(s) given |
| `--noPhase` | Drops `--phaseAndOffsets` |
| `--quadFit` | Drops `--noQuadFit` (i.e. enables the `-deltaBQ` quadratic baseline estimate) |
| `--tracks track-N [...]` | Narrows `track_dirs` used by *every* step in the tree (default: all `track-*`) |

---

## project.yaml reference

`setupNISARTracks` reads `project.yaml` from the project root directory (wherever
the script is run from).  The same file is read by `makeframetie.py`, which
looks two levels up from the `tiepoints/` directory it runs in — i.e. the
same project root.

| Key | Used by | Purpose |
|-----|---------|---------|
| `tie_plan_header_template` | `setupNISARTracks` | Source template for `tie_plan_header`; default: `templates/tie_plan_header` |
| `vel_thumb_plan_template` | `setupNISARTracks` | Source template for `vel_thumb_plan`; default: `templates/vel_thumb_plan` |
| `vel_thumb_header_template` | `setupNISARTracks` | Source template for `vel_thumb_header_XdashY`; default: `templates/vel_thumb_header` |
| `region` / `regionFile` | `setupNISARTracks` | Path to region YAML; its `dem` field replaces `<DEM>` in all templates |
| `framePattern` | `setupNISARTracks`, `makeframetie.py` | Glob pattern for virtual-frame directories (e.g. `000?`); controls which `velocityStats/` subdirs `makeframetie.py` processes |
| `sensor` | `makeframetie.py` | Sensor type (`NISAR`, `Sentinel1`, `TSX`, etc.) |

`makeframetie.py` will fall back to `sensor.yaml` if `project.yaml` is not
found, but will print a warning — rename `sensor.yaml` to `project.yaml`.

---

## Check mode

`--check` scans `track-*/<framePattern>/range.offsets.vrt` across the project
tree (using `framePattern` from `project.yaml`, defaulting to `000?`) and
reports:

1. **Missing products** — directories that have `range.offsets.vrt` but lack
   `velocity/mosaicOffsets.vx`.
2. **Residual sigma table** — for each processed directory, reads
   `sigma*sqrt(X2/n)` from:
   - `motion/rBaseline.deltabp` (range sigma)
   - `motion/rBaseline.deltabp.noIonosphere` (range sigma without ionosphere correction, if present)
   - `motion/az.est.const` (azimuth sigma)

   Rows are sorted by descending range sigma.  When both ionosphere-corrected and
   uncorrected range sigmas are available the better value is printed in **bold**,
   and per-column RSS sigmas are shown at the bottom.

---

## Project root

`setupNISARTracks` derives the project directory from `os.getcwd()` — run it
from the project root directory (the one containing `project.yaml` and the
`track-*/` directories).  If no `track-*` directories are found it exits with
an error.
