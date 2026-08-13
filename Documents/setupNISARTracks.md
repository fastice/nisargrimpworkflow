# setupNISARTracks — Multi-Track Directory Setup and Tie-Point Refresh

## Overview

`setupNISARTracks` initialises and maintains all `track-*` directories under a
project root.  It creates required subdirectories, optionally distributes
template files into every track, keeps each track's velocityStats/autoclean
masks up to date, and then drives `refreshties.py` across all tracks.
Auxiliary modes run velocity thumbnails, velocity-stats regions, and two
product-health reports.

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
| `--copyFiles` | off | Instantiate `tie_plan_header` and `vel_thumb_plan` from `project.yaml`-configured templates (default `templates/`) into every track's `tiepoints/` directory, substituting `<TRACK>`/`<DEM>`; exits after this step (skips `refreshties.py` and all other steps) |
| `--tiesOnly` | off | Run `refreshties.py` without computing velocity mosaics |
| `--overWrite` | off | Pass `--overWrite` to `refreshties.py` to rerun existing products, and to `vel_thumbs` when combined with `--runVelThumbs`; combined with `--copyFiles`, also makes the `tie_plan_header`/`vel_thumb_plan` copy overwrite existing files instead of skipping them |
| `--keepVz` | off | Pass `--keepVz` to `refreshties.py` to retain `.vz` and `.vz.geodat` files |
| `--runVelThumbs` | off | Run `vel_thumbs vel_thumb_plan` (`--overWrite` appended if `--overWrite` given) in every track's `tiepoints/` directory (10 threads); does **not** skip the rest of the default flow |
| `--runVelstatsregions` | off | Standalone, early-return action: run `makevelstatsregions.py`, then force a full tiepointing pass (so `makeframetie.py` regenerates the per-range `vel_thumb_plan`/thumbnails under the new extent), then rebuild the velocityStats reference — see [`--runVelstatsregions`](#--runvelstatsregions) below for why all three are chained |
| `--runVelStats` | off | Standalone, early-return action: run `velocityStats.py --doRA` in each track's `velocityStats/` dir, using the project's region file and velMap |
| `--updateAutoclean` | off | Force velStats+`autocleanNISAR.py` to (re)run for every track **and every frame**, even ones whose `range.offsets.vrt` already has an autoclean mask embedded. Default: only run velStats+autoclean for tracks with at least one unmasked frame, and within those, `autocleanNISAR.py --new` cleans only the not-yet-masked frames (already-cleaned frames are left as-is) |
| `--check` | off | Report product completeness and residual sigmas (see [Check mode](#check-mode)); exits without running any other steps |
| `--checkFrames` | off | Report virtual/physical frame-directory health (see [Check-frames mode](#check-frames-mode)); exits without running any other steps |
| `--year YYYY [YYYY ...]` | 2025 2026 | Year(s) forwarded to `refreshties.py` |
| `--noPhase` | off | Disable phase+offsets mode (by default NISAR processing uses `--phaseAndOffsets`) |
| `--quadFit` | off | Enable the `-deltaBQ` quadratic baseline correction estimate (by default NISAR processing uses `--noQuadFit`) |
| `--useSquint` | project.yaml `applySquintCorrection` | Apply the squint heading correction in `mosaic3d` (`-useSquint`) and `tiepoints -motion` (`addMotionCorrections.c`); overrides `project.yaml`'s `applySquintCorrection` for this run |
| `--noUseSquint` | project.yaml `applySquintCorrection` | Disable the squint correction for this run, overriding `project.yaml` |

### Examples

```bash
# First-time setup: create subdirectories and copy template files
setupNISARTracks --copyFiles

# Refresh tie points for all tracks, years 2025 and 2026 -- also runs
# velStats+autoclean first for any track that doesn't have masks yet
setupNISARTracks

# Refresh tie-point step only, overwriting existing results
setupNISARTracks --tiesOnly --overWrite

# Force velStats+autoclean to rerun even for tracks that already have masks
setupNISARTracks --tracks track-106 --updateAutoclean

# Regenerate velocity thumbnails in all tiepoints directories
setupNISARTracks --runVelThumbs

# Recompute velocity-stats regions and everything downstream of that
# (tiepointing + velStats rebuild)
setupNISARTracks --runVelstatsregions

# Check which track directories are missing velocity products
setupNISARTracks --check

# Check virtual/physical frame directory health
setupNISARTracks --checkFrames
```

---

## Processing steps (default mode)

### Step 0 — Track filtering (every flow except `--check`/`--checkFrames`)

Before anything else, `main()` narrows `track_dirs` (from `--tracks`, or all
`track-*` dirs by default) to only those containing at least one
`*_{framePattern}` virtual-frame directory. Tracks with no virtual-frame data
yet are silently dropped from every step below (with a one-line
`Skipping N track(s) with no virtual-frame data: [...]` notice) — they are
**not** treated as an error. `--check`/`--checkFrames` return before this
filtering runs, so they scan whatever `track-*` dirs exist regardless.

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

Files already present are skipped (unless `--overWrite` is also given).

### Step 3 — `vel_thumb_header` distribution (always)

After creating any missing `velocityStats/` subdirectories, `setupNISARTracks`
scans all `velocityStats/*-*` directories in the track and creates a
`vel_thumb_header_XdashY` file in `tiepoints/` for each one that does not
already have one.

The template source is `vel_thumb_header_template` from `project.yaml`, falling
back to `templates/vel_thumb_header`.  The same `<TRACK>` and `<DEM>`
substitutions are applied.  Existing files are **never overwritten** here
(regardless of `--overWrite`) — `makevelstatsregions.py` is what updates an
existing header's `resolution` line in place (see `--runVelstatsregions`
below); a new header still starts from the generic `0 0 0 0` auto-size
sentinel.

### Step 4 — velStats + autoclean (gated on masks; default flow only)

Skipped entirely by `--copyFiles`, `--runVelstatsregions`, and `--runVelStats`
(each of those is its own standalone/early-return action — see below). In the
plain default flow, for each track directory that does **not** already have
every frame's `range.offsets.vrt` carrying an embedded autoclean `<MaskBand>`
(or for every track, if `--updateAutoclean` is given):

1. `velocityStats.py --doRA` — refresh the region reference (`velocity.vr/.va/.er/.ea/.navg`).
2. `autocleanNISAR.py --tracks ... --new` — flag outlier range/azimuth offset
   pixels against that reference and embed the resulting masks into each
   frame's `range.offsets.vrt`/`azimuth.offsets.vrt`. `--new` restricts this
   to frames that don't already carry an embedded mask, so adding a few pairs
   to an already-processed track only cleans the newly-built frames — the
   existing frames' masks are left untouched. `--updateAutoclean` drops
   `--new`, re-cleaning every frame. (Reference caveat: velStats refreshes the
   region reference over the whole track first, so with `--new` the pre-existing
   frames stay cleaned against their prior reference while the new frames use
   the just-refreshed one; use `--updateAutoclean` if you want every frame
   re-cleaned against a single consistent reference.)

Tracks that already have masks are skipped here (no wasted recompute), but
the tiepointing step below always runs for every track regardless — it's
what actually benefits from the (possibly just-refreshed) masks, since
`rparams`/`azparams` honor an embedded mask by default.

### Step 5 — Tie-point refresh

Calls `refreshties.py` once for all tracks:

```
refreshties.py [-tiesOnly] [--overWrite] [--keepVz] --phaseAndOffsets --noQuadFit --yaml [--useSquint] -toRun="[track-NN,...]" YYYY [YYYY ...] -noPrompt
```

---

## What each flag sets in motion

`main()` (in `setupNISARTracks.py`) is the actual control flow; the tree below follows
it exactly. Most flags don't change *which* programs run — they just change the arguments on
the same `refreshties.py` call — so only the structurally distinct cases get their own figure.
`--check`, `--checkFrames`, `--copyFiles`, `--runVelstatsregions`, and `--runVelStats` are the
flags that skip the default velStats/autoclean/tiepointing sequence entirely.

### No flags (default)

```
setupNISARTracks
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    │   ├── create <track>/tiepoints/                       (if missing)
    │   ├── create <track>/velocityStats/{X}-{Y}/           (if missing)
    │   └── write <track>/tiepoints/vel_thumb_header_XdashY (if missing; never overwritten)
    ├── needs_autoclean = tracks where _track_has_autoclean_mask() is False
    │   ├── run_vel_stats(needs_autoclean, region_path, vel_map)
    │   │   └── csh -c 'velocityStats.py --doRA -regionFile <region> -velmap <velMap>'  (per track, cwd = <track>/velocityStats/)
    │   └── run_autoclean(needs_autoclean, frame_pattern, new=True)
    │       └── csh -c 'autocleanNISAR.py --tracks <t..> -framePattern <pat> --new'  (one pooled run, cwd = PROJECT_DIR)
    │           (--new: clean only frames without an embedded mask)
    └── run_refresh_ties()
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--copyFiles`

`--copyFiles` is a standalone, first-time-setup step: `setupNISARTracks` exits
right after `setup_track_dirs()` and never calls `refreshties.py`, the
velStats/autoclean step, `run_vel_thumbs`, or `run_velstats_regions` — even if
those flags are also passed.

```
setupNISARTracks --copyFiles
└── main()
    └── setup_track_dirs(track_dirs, copyFiles=True, overwrite=False)
        ├── create <track>/tiepoints/                       (if missing)
        ├── copy tie_plan_header → <track>/tiepoints/       (if missing; <TRACK>/<DEM> substituted)
        ├── copy vel_thumb_plan  → <track>/tiepoints/       (if missing; <TRACK>/<DEM> substituted)
        ├── create <track>/velocityStats/{X}-{Y}/           (if missing)
        └── write <track>/tiepoints/vel_thumb_header_XdashY (if missing)
    ⇥ return                                                (nothing else called)
```

### `--copyFiles --overWrite`

Combining the two makes the `tie_plan_header`/`vel_thumb_plan` copy
overwrite existing files instead of skipping them (`vel_thumb_header_XdashY`
is still never overwritten by this step — see Step 3). `--overWrite` is not
forwarded to `refreshties.py` here, since `--copyFiles` returns before that
call would happen.

```
setupNISARTracks --copyFiles --overWrite
└── main()
    └── setup_track_dirs(track_dirs, copyFiles=True, overwrite=True)
        ├── create <track>/tiepoints/                       (if missing)
        ├── copy tie_plan_header → <track>/tiepoints/       (always; <TRACK>/<DEM> substituted)
        ├── copy vel_thumb_plan  → <track>/tiepoints/       (always; <TRACK>/<DEM> substituted)
        ├── create <track>/velocityStats/{X}-{Y}/           (if missing)
        └── write <track>/tiepoints/vel_thumb_header_XdashY (if missing; never overwritten)
    ⇥ return                                                (nothing else called)
```

### `--tiesOnly`

Only changes the `refreshties.py` arguments — the velStats/autoclean gating
step still runs first, same as the default flow.

```
setupNISARTracks --tiesOnly
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)        (same as default)
    ├── needs_autoclean gating                                (same as default)
    └── run_refresh_ties()
        └── refreshties.py -tiesOnly --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--overWrite`

```
setupNISARTracks --overWrite
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── needs_autoclean gating                                (same as default; --overWrite does not
    │                                                           force autoclean -- use --updateAutoclean)
    └── run_refresh_ties()
        └── refreshties.py --overWrite --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--keepVz`

```
setupNISARTracks --keepVz
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── needs_autoclean gating                                (same as default)
    └── run_refresh_ties()
        └── refreshties.py --keepVz --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--updateAutoclean`

Forces the velStats+autoclean step for every track, not just ones missing
masks. Everything else (including tiepointing afterward) proceeds as in the
default flow.

```
setupNISARTracks --updateAutoclean
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── needs_autoclean = ALL track_dirs (forced, ignoring _track_has_autoclean_mask())
    │   ├── run_vel_stats(needs_autoclean, region_path, vel_map)
    │   └── run_autoclean(needs_autoclean, frame_pattern, new=False)
    │       └── csh -c 'autocleanNISAR.py --tracks <t..> -framePattern <pat>'  (no --new: reclean EVERY frame)
    └── run_refresh_ties()
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--runVelThumbs`

Adds a call but does **not** skip the rest of the default flow — velStats/
autoclean gating and the tie-point refresh both still run.

```
setupNISARTracks --runVelThumbs
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── run_vel_thumbs(track_dirs, over_write=args.overWrite)   (10-thread pool, one thread/track)
    │   └── csh -c 'vel_thumbs vel_thumb_plan [--overWrite]'    (cwd = <track>/tiepoints/)
    ├── needs_autoclean gating                                  (same as default)
    └── run_refresh_ties()                                       ← still runs, no early return
        └── refreshties.py --phaseAndOffsets --noQuadFit --yaml \
                -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
```

### `--runVelstatsregions`

A standalone, early-return action — but unlike earlier versions of this
script, it does **not** skip tiepointing; it forces it. `makevelstatsregions.py`
rewrites each `vel_thumb_header_<range>`'s `resolution` line to the common
domain covering all data for that range (filling in the `0 0 0 0` auto-size
sentinel). The actual per-range file `vel_thumbs` reads
(`vel_thumb_plan_<year>dash<range>`) is only regenerated from that header by
`makeframetie.py` — which is itself only invoked as part of a tiepointing
pass, and which also calls `vel_thumbs` on the regenerated file. So this flag
chains all three steps needed to make the new extent actually take effect:

```
setupNISARTracks --runVelstatsregions
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    ├── run_velstats_regions(track_dirs)
    │   └── csh -c 'makevelstatsregions.py'                     (cwd = <track>/, one at a time)
    │       rewrites <track>/tiepoints/vel_thumb_header_XdashY's resolution line
    ├── run_refresh_ties(track_dirs, ties_only=False, over_write=True, ...)  (forced, not from args.tiesOnly/args.overWrite)
    │   └── refreshties.py --overWrite [--keepVz] --phaseAndOffsets --noQuadFit --yaml [--useSquint] \
    │           -toRun="[track-1,track-2,...]" 2025 2026 -noPrompt
    │       (this is what actually invokes makeframetie.py, which regenerates
    │        vel_thumb_plan_<year>dash<range> from the fresh header and calls
    │        vel_thumbs on it)
    └── run_vel_stats(track_dirs, region_path, vel_map)
        └── csh -c 'velocityStats.py --doRA -regionFile <region> -velmap <velMap>'  (cwd = <track>/velocityStats/)
    ⇥ return
```

### `--runVelStats`

Standalone, early-return action — refreshes only the velocityStats reference,
nothing else.

```
setupNISARTracks --runVelStats
└── main()
    ├── setup_track_dirs(track_dirs, copyFiles=False)
    └── run_vel_stats(track_dirs, region_path, vel_map)
        └── csh -c 'velocityStats.py --doRA -regionFile <region> -velmap <velMap>'  (cwd = <track>/velocityStats/)
    ⇥ return                                                (refreshties.py never called)
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

### `--checkFrames`

Also a fully separate branch — exits before any directory setup or track processing.

```
setupNISARTracks --checkFrames
└── main()
    └── check_frames()
        ├── for each physical source frame (suffix doesn't match framePattern):
        │     report missing range/azimuth offsets, phase, ionosphere, geodat files
        ├── for each virtual/merged frame (suffix matches framePattern):
        │     report missing range offsets, phase.uw.vrt, geodat, pairinfo
        ├── detect + optionally remove empty-H5-stub and frames.txt-only-stub directories
        └── print frame counts and incomplete-directory lists
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
| `--useSquint` / `--noUseSquint` | Adds/drops `--useSquint`, overriding `project.yaml`'s `applySquintCorrection` default |
| `--tracks track-N [...]` | Narrows `track_dirs` used by *every* step in the tree (default: all `track-*`) |

Note: `--overWrite` is *not* purely a `refreshties.py`-argument flag — it also
threads into `run_vel_thumbs()` (adds `--overWrite` to the `vel_thumbs` call
under `--runVelThumbs`) and into the `--copyFiles` template-copy step.

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
| `region` / `regionFile` | `setupNISARTracks` | Path to region YAML; its `dem` field replaces `<DEM>` in all templates, and its `velMap` field is passed to `velocityStats.py --doRA -velmap` by `run_vel_stats()` (default velStats/autoclean step, `--runVelStats`, `--runVelstatsregions`) |
| `framePattern` | `setupNISARTracks`, `makeframetie.py` | Glob pattern for virtual-frame directories (e.g. `00??`); controls which `velocityStats/` subdirs are processed and which frame dirs `autocleanNISAR.py --allFrames`/`_track_has_autoclean_mask()` scan |
| `sensor` | `makeframetie.py` | Sensor type (`NISAR`, `Sentinel1`, `TSX`, etc.) |
| `applySquintCorrection` | `setupNISARTracks` | Boolean (default `false`); sets the per-project default for the squint heading correction in `mosaic3d` (`-useSquint`) and `tiepoints -motion` (`addMotionCorrections.c`). Overridable per-run by `--useSquint`/`--noUseSquint` CLI flags. See `mosaicSource/CLAUDE.md` "Squint" for the physics and current NISAR-default rationale. |
| `velThumbOutput` | `vel_thumbs`/`insarworkflow.velThumbs` (not `setupNISARTracks` itself) | `'tiff'` selects GeoTIFF (`mosaicOffsets.vrt`/`.vr.tif`/`.va.tif`) output from `mosaic3d`; anything else (including absent) defaults to legacy binary (`mosaicOffsets.vx`/`.vy` + `.geodat`). `makevelstatsregions.py`'s `_boundsFromVelDir()` supports both formats (tries `.vrt` first, falls back to `.vx.geodat`) for exactly this reason — projects/tracks processed before this key was set (or before it was respected) can have a mix of the two formats on disk |
| `velocityStatsMode` | `velocityStats.py`, `autoclean.py`/`autocleanNISAR.py` (not `setupNISARTracks` itself) | `'RA'` or `'XY'` — which velocity component pair (range/azimuth vs. x/y) autoclean and velocityStats operate on |

`makeframetie.py` will fall back to `sensor.yaml` if `project.yaml` is not
found, but will print a warning — rename `sensor.yaml` to `project.yaml`.

---

## Check mode

`--check` scans `track-*/*_000?/range.offsets.vrt` across the project tree —
note this glob is **hardcoded** to `*_000?` in `check_products()`, unlike
`check_frames()`/the rest of the script, which read `framePattern` from
`project.yaml`; a project using a different `framePattern` convention won't
be scanned correctly by `--check` — and reports:

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

## Check-frames mode

`--checkFrames` walks every entry directly under each `track-*` directory
whose name matches `<digits>_<digits>` (skipping any with an `Exclude` file)
and classifies it by whether its suffix matches `framePattern`:

- **Physical source frames** (suffix does *not* match `framePattern`) are
  checked for the files a virtual-frame VRT would reference: range/azimuth
  offsets, phase, ionosphere correction, geodat, secondary geodat.
- **Virtual/merged frames** (suffix *matches* `framePattern`) are checked for
  the assembled products: range offsets, `phase.uw.vrt`, geodat, `*.pairinfo`.

It also detects two kinds of stub directories left behind by interrupted
downloads/merges — a frame dir containing only an empty `H5/` subdirectory,
or a virtual frame containing only `frames.txt` — and, if any are found,
prompts once (`[y/N]`) to remove all of them.

Output: frame counts (physical/virtual/excluded), a `framePattern` echo, and
a list of any incomplete physical or virtual frames with which files are
missing.

---

## Project root

`setupNISARTracks` derives the project directory from `os.getcwd()` — run it
from the project root directory (the one containing `project.yaml` and the
`track-*/` directories).  If no `track-*` directories are found it exits with
an error.
