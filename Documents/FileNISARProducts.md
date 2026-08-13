# FileNISARProducts — NISAR HDF5 Product Filer

## Overview

`FileNISARProducts` organises a flat download directory of NISAR HDF5 products
into the track-based directory tree expected by [SetupNISAR](SetupNISAR.md).
It should be run once on each fresh data download before calling `SetupNISAR`.

RUNW products are the primary index: each RUNW determines the track, orbit, and
frame numbers that define where all companion products (ROFF, RIFG, RSLC) are
filed.  RSLC products with no companion RUNW are discovered in a second pass and
filed to `source/` only.

---

## Input layout

Products must be pre-sorted by type under a common root:

```
inputPath/
├── RUNW/
│   └── NISAR_L2_PR_RUNW_*.h5
├── ROFF/
│   └── NISAR_L2_PR_ROFF_*.h5
├── RIFG/
│   └── NISAR_L2_PR_RIFG_*.h5   (optional)
└── RSLC/
    └── NISAR_L1_PR_RSLC_*.h5   (optional)
```

---

## Output layout

```
outputPath/
└── track-{N}/
    ├── source/                         — all symlinks to original files
    │   ├── NISAR_L2_PR_RUNW_....h5
    │   ├── NISAR_L2_PR_ROFF_....h5
    │   ├── NISAR_L2_PR_RIFG_....h5
    │   └── NISAR_L1_PR_RSLC_....h5    — RSLC here, not in H5/
    ├── {orbit1}_{frame}/               e.g. 12345_10/ (frame unpadded)
    │   └── H5/
    │       ├── NISAR_L2_PR_RUNW_....h5  (symlink)
    │       ├── NISAR_L2_PR_ROFF_....h5  (symlink)
    │       └── NISAR_L2_PR_RIFG_....h5  (symlink, if present)
    └── unfiled/                        — duplicates and orphans
        ├── NISAR_L2_PR_RUNW_....h5     (symlink)
        └── log
```

RSLC files are symlinked into `source/` only — `wrapH5sInFrameDir` in
`SetupNISAR` does not yet support RSLC, so they are not placed in `H5/`.

---

## Usage

```
FileNISARProducts [--inputPath PATH] [options]
```

The input directory comes from **either** `--inputPath` **or** the `archiveDir`
key of an `autoupdate.yaml` sitting in `--outputPath` (the autoupdate workflow
runs from the project dir where that file lives). `--inputPath` always wins if
given; if neither source provides a path the program exits with an error naming
both.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--inputPath PATH` | from `autoupdate.yaml` | Root directory with `RUNW/`, `ROFF/`, `RIFG/`, `RSLC/` product subdirectories. If omitted, `archiveDir` from an `autoupdate.yaml` in `--outputPath` is used |
| `--outputPath PATH` | `.` (cwd) | Root directory for the `track-{N}/` output tree |
| `--products TYPE ...` | all | Product types to file: `RUNW ROFF RIFG RSLC`. RUNW is always used as the index (see [RUNW as driver](#runw-as-driver)) |
| `--firstDate YYYYMMDD` | none | Skip products whose reference date is before this date |
| `--lastDate YYYYMMDD` | none | Skip products whose reference date is after this date |
| `--firstOrbit N` | 1 | Skip orbits numbered below this value |
| `--lastOrbit N` | 999999 | Skip orbits numbered above this value |
| `--reFile` | off | Re-process all products even if `source/` symlinks already exist |
| `--verbose` | off | Print per-file detail; without this flag a progress bar is shown |
| `--nThreads N` | 8 | Worker processes used to read RUNW HDF5 orbit/frame metadata in parallel (the dominant cost). Use `1` for serial processing |

### Examples

```bash
# File everything from a download directory into the current directory
FileNISARProducts --inputPath /data/nisar/downloads

# Autoupdate workflow: run from the project dir; inputPath comes from
# autoupdate.yaml's archiveDir (no --inputPath needed)
FileNISARProducts

# File into a specific output tree, only products from 2025
FileNISARProducts --inputPath /data/nisar/downloads --outputPath /data/nisar/orbits \
    --firstDate 20250101 --lastDate 20251231

# File only RSLC products (including ones with no companion RUNW)
FileNISARProducts --inputPath /data/nisar/downloads --products RSLC

# Re-file after adding new downloads to an existing tree
FileNISARProducts --inputPath /data/nisar/downloads --reFile

# Serial run (no worker pool) — e.g. for debugging or a tiny input
FileNISARProducts --inputPath /data/nisar/downloads --nThreads 1
```

---

## Processing steps

### 1. Duplicate resolution (pre-processing pass)

All RUNW files are grouped by `(track, frame, referenceDate)`.  Two products
are duplicates only if they share the same reference date **and** frame — meaning
they would both land in the same `orbit1_frame/` directory.  Products from
different orbit passes (different reference dates) that happen to share track
and frame are independent and are never treated as duplicates.

When a group has more than one candidate, `selectBestRUNW` picks the winner:

1. **Shortest temporal baseline** — `|date2Start − date1Start|` in days
2. **Newest modification time** — proxy for processing version
3. **Longest along-track coverage** — reference + secondary frame durations
   (later end times = a longer frame). Breaks ties between the *same* pair
   framed to different end times: identical start times and baseline, differing
   only in `date1End`/`date2End`. Without this the winner fell to filename sort
   and could keep the shorter (less-coverage) product.

Losers are symlinked into `track-{N}/unfiled/` and a human-readable log entry
is written to `track-{N}/unfiled/log` explaining why each was not selected.

### 2. Main filing pass (winners only)

The pass runs in three stages: (1) build the already-filed skip set and select
the winners still needing work, (2) read those winners' HDF5 orbit/frame
metadata **in parallel** across `--nThreads` worker processes, then (3) do all
filesystem mutations serially. For each winning RUNW:

1. Derive track from filename (no HDF5 open needed) for the skip check.
2. Apply `--firstDate` / `--lastDate` filter on the reference date from the filename.
3. Create `track-{N}/source/` and symlink all requested product types there:
   - RUNW, ROFF, RIFG: companion lookup by filename substitution
   - RSLC: glob `inputPath/RSLC/` by matching track, direction, and frame fields
4. Read the RUNW HDF5 (`noLoadData=True`) for `referenceOrbit`, `frame`, and
   mixed-mode status — done up front in parallel (see below).
5. Apply `--firstOrbit` / `--lastOrbit` filter; skip mixed-mode frames.
6. Create `track-{N}/{orbit1}_{frame}/H5/` and symlink the L2 products (RUNW,
   ROFF, RIFG) there.  RSLC is not symlinked into `H5/`.

**Smart skip**: the already-filed RUNW source symlinks are collected in a single
glob per track into an in-memory set, and each winner is tested against that set
— avoiding a per-winner `os.path.exists()` network stat (thousands of NFS round
trips on a re-run). If the RUNW is already filed and RUNW is in `--products`,
the entry is skipped without opening the HDF5. Use `--reFile` to override. When
RUNW is not in `--products`, the skip is disabled so orbit/frame can still be
extracted from the HDF5 for other products.

**Parallel metadata read**: opening each new RUNW HDF5 over NFS (~0.5–1.5 s) is
the dominant cost, so the winners' orbit/frame/mixed-mode fields are read up
front in a process pool (`--nThreads`, default 8), and only the fast filesystem
operations (mkdir/symlink) run serially afterwards. Separate worker processes
(not threads) keep each HDF5 open independent and side-step h5py thread-safety
caveats. `--nThreads 1` runs the reads serially with no pool. The output tree is
identical regardless of thread count.

### 3. Standalone RSLC pass

After the main loop, `inputPath/RSLC/` is scanned for any RSLC files that were
not found via a companion RUNW (e.g. acquisitions that were never used as a
reference or secondary in a processed interferogram).  Each unfiled RSLC is:

- Date-filtered by its single acquisition date against `--firstDate`/`--lastDate`
- Symlinked into `track-{N}/source/`

No `orbit_frame/H5/` directory is created because there is no L2 product to
associate with it.

### 4. Companion check

After all filing, every `orbit_frame/` directory is checked for the presence of
both RUNW and ROFF.  A directory with only one of the two cannot be processed by
`SetupNISAR` (the default ionosphere path requires both).  The orphan product is
moved to `track-{N}/unfiled/` with a log entry.

This check is skipped when `--products` does not include both `RUNW` and `ROFF`,
since incomplete pairs are expected in that case.

### 0. Recovery of previously-unfiled orphans (runs first)

Before anything else, `track-{N}/unfiled/` is scanned for RUNWs that were moved
there for a **missing companion ROFF** on an earlier run — the common case being
a RUNW that downloaded before its ROFF.  For each, if the companion ROFF has
since arrived in `inputPath/ROFF/`, the orphan's `unfiled/` and `source/` links
are removed (with a `Recovered: companion ROFF arrived` log entry) so the main
pass re-files the now-complete pair into `orbit_frame/H5/`.

A missing-companion orphan is told apart from a **duplicate loser** (§1) by the
fact that it still has a `track-{N}/source/` symlink — losers never get one — so
duplicate losers are never resurrected, even though their ROFF also exists.  This
runs before the skip set is built (dropping the `source/` link is what lets the
main pass re-process the RUNW) and only when both RUNW and ROFF are being filed.

---

## RUNW as driver

RUNW is always the primary index even when not in `--products`.  This is because:

- ROFF and RIFG companion filenames are derived from the RUNW filename by
  simple product-type substitution.
- RSLC companion search matches the RUNW's track/direction/frame fields.
- Orbit number, frame number, and mixed-mode detection all require opening the
  RUNW HDF5.

When RUNW is not in `--products`, the RUNW files are opened for metadata only
and are not symlinked anywhere.  The progress output notes:

```
Found 1913 RUNW products
  (RUNW used as index only — not filed)
```

---

## Duplicate log format

`track-{N}/unfiled/log` contains one entry per unfiled product:

```
2026-05-15 09:23:11  NISAR_L2_PR_RUNW_....h5
  Track: 64,  Frame: 10,  RefDate: 2025-06-01
  Reason: longer temporal baseline (24d vs 12d)
  Filed winner: 12345_10/NISAR_L2_PR_RUNW_....h5
```

---

## Filename format

`parseFileName` supports both product types:

**RUNW / ROFF / RIFG (15 fields):**
```
NISAR_L2_PR_RUNW_{cycle}_{track}_{direction}_{frame}_004_{bw}_{pol}_{date1start}_{date1end}_{date2start}_{date2end}.h5
```

**RSLC (13 fields):**
```
NISAR_L1_PR_RSLC_{cycle}_{track}_{direction}_{frame}_{bw}_{pol}_{mode}_{date1start}_{date1end}.h5
```

`track`, `frame`, and `cycle` are returned as strings with leading zeros
stripped.  Date fields are returned as `datetime` objects.
