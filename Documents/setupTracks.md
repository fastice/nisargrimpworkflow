# setupTracks — Multi-Track Directory Setup and Tie-Point Refresh

## Overview

`setupTracks` initialises and maintains all `track-*` directories under a
fixed project root (`/Volumes/insar1/ian/NISAR/realNISAR/greenlandProject`).
It creates required subdirectories, optionally distributes template files from
`track-88` to every other track, and then drives `refreshties.py` across all
tracks.  Auxiliary modes run velocity thumbnails and velocity-stats regions.

---

## Usage

```
setupTracks [options]
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--copyFiles` | off | Copy `tie_plan_header` and `vel_thumb_plan` from `track-88/tiepoints/` to every other track, substituting the track number |
| `--tiesOnly` | off | Pass `-tiesOnly` to `refreshties.py` (skip non-tie-point steps) |
| `--overWrite` | off | Pass `--overWrite` to `refreshties.py` to rerun existing products |
| `--keepVz` | off | Pass `--keepVz` to `refreshties.py` to retain `.vz` and `.vz.geodat` files |
| `--runVelThumbs` | off | Run `vel_thumbs vel_thumb_plan` in every track's `tiepoints/` directory (10 threads) |
| `--runVelstatsregions` | off | Run `makevelstatsregions.py` in every track directory, then exit |
| `--year YYYY [YYYY ...]` | 2025 2026 | Year(s) forwarded to `refreshties.py` |
| `--check` | off | Report product completeness and residual sigmas (see [Check mode](#check-mode)); exits without running any other steps |

### Examples

```bash
# First-time setup: create subdirectories and copy template files from track-88
setupTracks --copyFiles

# Refresh tie points for all tracks, years 2025 and 2026
setupTracks

# Refresh tie-point step only, overwriting existing results
setupTracks --tiesOnly --overWrite

# Regenerate velocity thumbnails in all tiepoints directories
setupTracks --runVelThumbs

# Check which track directories are missing velocity products
setupTracks --check
```

---

## Processing steps (default mode)

### Step 1 — Directory creation

For every `track-*` directory found under `PROJECT_DIR`:

- Creates `<track>/tiepoints/` if it does not exist.
- Creates `<track>/velocityStats/0000-0001/` if it does not exist.

### Step 2 — Template file distribution (`--copyFiles`)

Copies two template files from `track-88/tiepoints/` to each other track,
substituting the destination track number:

| Source | Destination | Substitution |
|--------|-------------|--------------|
| `track-88/tiepoints/tie_plan_header` | `<track>/tiepoints/tie_plan_header` | `track-88` → `track-NN` |
| `track-88/tiepoints/vel_thumb_plan` | `<track>/tiepoints/vel_thumb_plan` | `-88` → `-NN` |

Files already present are skipped.

### Step 3 — Velocity thumbnails (`--runVelThumbs`)

Runs `csh -c 'vel_thumbs vel_thumb_plan'` in each track's `tiepoints/`
directory using a 10-thread pool.

### Step 4 — Velocity stats regions (`--runVelstatsregions`)

Runs `makevelstatsregions.py` sequentially in each track directory, then
exits (skipping the `refreshties.py` call).

### Step 5 — Tie-point refresh

Calls `refreshties.py` once for all tracks:

```
refreshties.py [-tiesOnly] [--overWrite] [--keepVz] -toRun="[track-NN,...]" YYYY [YYYY ...] -noPrompt
```

---

## Check mode

`--check` scans `track-*/*_000?/range.offsets.vrt` across the project tree and
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

The project root is hardcoded in the source as:

```
PROJECT_DIR = '/Volumes/insar1/ian/NISAR/realNISAR/greenlandProject'
```

All `track-*` directories are discovered relative to this path.
