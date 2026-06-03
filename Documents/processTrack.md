# processTrack — Batch SetupNISAR Runner for a Single Track

## Overview

`processTrack` runs [SetupNISAR](SetupNISAR.md) for every orbit found in a
single track directory.  It discovers all `<orbit>_<frame>` subdirectories,
extracts the unique orbit numbers, and calls `SetupNISAR <orbit>` for each one
in sequence, passing through the `--overWrite`, `--overWritePhase`, and
`--RUNWOnly` flags unchanged.

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
| `--overWrite` | Pass `--overWrite` to every `SetupNISAR` call — re-run all per-frame conversions |
| `--overWritePhase` | Pass `--overWritePhase` — re-run RUNW and ionosphere steps only |
| `--RUNWOnly` | Pass `--RUNWOnly` — skip ROFF and ionosphere; process phase/coherence only |

### Examples

```bash
# Process all orbits in track-88
processTrack track-88

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
