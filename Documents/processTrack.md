# processTrack — Batch SetupNISAR Runner for a Single Track

## Overview

`processTrack` runs [SetupNISAR](SetupNISAR.md) for every orbit found in a
single track directory.  It discovers all `<orbit>_<frame>` subdirectories,
extracts the unique orbit numbers, and calls `SetupNISAR <orbit>` for each one
(sequentially by default, or `--threads N` orbits concurrently), passing its
flags through to `SetupNISAR` unchanged.

Use this to (re-)process all pairs in one track with a single command, rather
than calling `SetupNISAR` once per orbit manually.

---

## Usage

```
processTrack <track> [options]
```

### Positional argument

| Argument | Description |
|----------|-------------|
| `track` | Path to the track directory (e.g. `track-88`) |

### Options

| Flag | Description |
|------|-------------|
| `--overWrite all\|ORBITS` | Rebuild existing products. **Requires an explicit scope** — see below |
| `--overWritePhase` | Pass `--overWritePhase` — re-run RUNW and ionosphere steps only |
| `--RUNWOnly` | Pass `--RUNWOnly` — skip ROFF and ionosphere; process phase/coherence only |
| `--correlationOnly` | Pass `--correlationOnly` — process correlation products only (no phase/offsets) |
| `--debugIono` | Pass `--debugIono` — debug-named iono intermediates + extra debug VRTs |
| `--sepIceRock` | Pass `--sepIceRock` — ice-only iono estimation + rock-seeded global fill |
| `--geodatsOnly` | Pass `--geodatsOnly` — re-merge virtual-frame geodats only, no reprocessing |
| `--clean` | Pass `--clean` — remove computed outputs for every orbit, then exit |
| `--cleanDebug` | Pass `--cleanDebug` — empty all `debug/` directories, then exit |
| `--noPrompt` | Pass `--noPrompt` — skip the confirmation prompt for `--clean`/`--cleanDebug` |
| `--bakeOnly` | Pass `--bakeOnly` — bake existing virtual-frame VRTs into flat GeoTIFFs, no reprocessing |
| `--new` | **Deprecated, now the default.** Accepted so old command lines keep working; changes nothing |
| `--firstDate YYYY-MM-DD` | Pass `--firstDate` — only process pairs whose reference date is on/after this date |
| `--lastDate YYYY-MM-DD` | Pass `--lastDate` — only process pairs whose reference date is on/before this date |
| `--threads N` | Number of orbits to process concurrently (default 1; separate from `-ompThreads` inside the C binaries) |

### What runs by default, and how to rebuild

**Default is new-only.** Every orbit in the track is visited, and `--new` is passed to each
`SetupNISAR`, so only virtual frames that do not already have a product VRT are built. Existing
products are left alone. This is safe to re-run.

**Rebuilding requires an explicit scope.** `--overWrite` is not a bare switch; it takes either
the word `all` or a list of orbits:

| Form | Effect |
|------|--------|
| *(omitted)* | every orbit visited, new virtual frames only — existing products untouched |
| `--overWrite all` | every orbit in the track rebuilt |
| `--overWrite 5430` | **only** orbit 5430 is processed, and rebuilt in full |
| `--overWrite 5430,5431` | only those orbits, comma or space separated |

With an orbit list the run is *restricted* to those orbits, so a rebuild cannot reach healthy
orbits by accident. A listed orbit is always done in full: every virtual-frame section of it
(`5430_0000`, `5430_0010`, …) is rebuilt, because `SetupNISAR` works a whole orbit at a time.
Orbits named but not present in the track produce a warning and are ignored.

`--firstDate`/`--lastDate` still apply and compose with the above. They do **not** limit which
orbits `processTrack` dispatches — they are passed to each `SetupNISAR`, which then skips orbits
whose reference date falls outside the window (`Skipping orbit NNNN: reference date ... outside
[...]`). Use them to narrow within a selection; use an orbit list to choose the selection.

`processTrack` prints its selection and mode before starting, e.g.

```
processTrack: track-141  orbits 5430  mode: overWrite orbits 5430
processTrack: date window [2026-08-11, 2026-08-11] -- orbits outside it are skipped by SetupNISAR
```

### Examples

```bash
# Default: build only what is missing, leave existing products alone
processTrack track-88

# Rebuild one bad orbit, in full, without touching the rest of the track
processTrack track-141 --overWrite 5430 -noPrompt

# Rebuild two orbits
processTrack track-23 --overWrite 4793,5139 -noPrompt

# Rebuild the entire track (the old bare --overWrite behaviour)
processTrack track-88 --overWrite all

# Re-run phase and ionosphere only
processTrack track-88 --overWritePhase

# Process only phase products, skip offsets entirely
processTrack track-88 --RUNWOnly
```

---

## Processing steps

1. Globs `<track>/*_*` to find all orbit–frame subdirectories.
2. Extracts unique orbit numbers (the part before the first `_` in each directory name).
3. Skips any directory whose orbit contains `tie` (e.g. `tiepoints`).
4. Calls `SetupNISAR <orbit>` with `cwd=<track>` for each orbit in sorted order.

---

## Dependencies

| Tool | Role |
|------|------|
| [`SetupNISAR`](SetupNISAR.md) | Per-orbit NISAR HDF5 → GrIMP conversion orchestrator |
