# autoupdate — daily incremental NISAR update driver

`autoupdateNISAR` (`autoupdate.py`) is a cron-oriented driver that searches ASF
for new NISAR granules, downloads them, and runs the full GrIMP processing chain
on the newly-completed data. It is configured entirely by an `autoupdate.yaml`
that normally sits beside the project's `project.yaml`.

## Setup checklist

One-time setup to bring a project onto the daily driver. Detail for each item is
in the sections below; this is the ordered "did I do everything" list.

1. **Conda `base` env has the tools.** `autoupdateNISAR`, `searchASF`,
   `ariaDownload`, `SetupNISAR` (and GDAL) must all be importable/on `PATH` in the
   env the wrapper activates (`base`). Install/update `nisargrimpworkflow` into it.
2. **GrIMP C binaries built and installed** to `~/bin/x86_64` (`simoffsets`,
   `cullst`, `siminsar`, `intfloat`, `mosaic3d`). Override the location with
   `BINDIR` if elsewhere.
3. **Earthdata credentials in place** for `ariaDownload` (standard `~/.netrc`
   machine `urs.earthdata.nasa.gov` login). Without them every download fails.
4. **`project.yaml` present** in the project dir and carrying the region-dependent
   processing keys `SetupNISAR` reads — in particular `sepIceRock` (`True` for
   Greenland, `False` for Antarctica). These do **not** go in `autoupdate.yaml`.
5. **Create `autoupdate.yaml` beside `project.yaml`** with the required keys
   (`archiveDir`, `region`, `lookBack`, `products`, `bandwidths`) plus any optional
   tuning (`nThreads`, `ompThreads`, `setupFlags`, `projectDir`). Minimal example:

   ```yaml
   archiveDir: /Volumes/insar4/ian/Antarctica2026/Antarctica40/archive
   region: Antarctica
   lookBack: 3          # search today-3 .. today
   products: RUNW ROFF
   bandwidths: 40       # GrIMP bandwidth token
   nThreads: 20         # SetupNISAR orbit jobs in parallel
   ompThreads: 6        # OpenMP width per orbit; balance nThreads*ompThreads vs cores
   ```

6. **`archiveDir` exists and is writable** (staging + RUNW/ROFF population live
   under it).
7. **Hand-test before cron:** run
   `scripts/runAutoupdateNISAR.sh <projectDir>` once. This exercises the same PATH /
   mount / credential setup cron will use and catches problems early. Check the
   run's `dailyDownloadLogs/summary_<date>.txt`. Use `--noDownload` to dry-run the
   processing chain on pre-staged folders without hitting ASF.
8. **Add the cron line** (see *Running as a cron job*): absolute path to the
   wrapper, a `flock` lock file **inside the project dir**, and a redirect into
   `<projectDir>/dailyDownloadLogs/cron.log`. One line per project.

## Why the one-day completion buffer

A single daily run can capture a track/cycle acquisition **pass** that is only
partially downloaded — more frames of the same pass may arrive on a later day. If
partial passes flowed straight into processing, virtual frames would be built (and
mosaics concatenated) from incomplete data. So each day's download is held one day
in a dated **staging folder**; the next day any late-arriving frames of the same
`(cycle, track)` are merged in, and only then is the (now-complete) pass released
into the general population and processed.

Pass identity is `(cycle, track)` — **not** the reference date, which can flip
mid-track at a UTC-midnight boundary.

## Pipeline (one run)

```
searchGranules          searchASF over [today-lookBack, today] -> URL list
downloadGranules        ariaDownload into archiveDir/<MM-DD-YYYY>/  (today's staging folder)
consolidateStaging      move today's frames of any pass also present in the aged
                        (non-today) folders into the most-recent prior folder
                        (matched against the UNION of all prior folders, so a
                        skipped-cron backlog reconciles as a group)
releaseStaging          move every aged (non-today) staging folder into
                        archiveDir/{RUNW,ROFF}; returns the released granule list
  (if nothing released -> done)
FileNISARProducts       file released RUNW/ROFF into track-<N>/<orbit>_<frame>/H5/
SetupNISAR <orbit>      per affected orbit (derived from the released granules),
                        cwd=track-<N>, WITHOUT --new  (rebuilds virtual frames in
                        place when late frames extend them; see below).
                        Up to nThreads orbits run in parallel; within an orbit the
                        frames run serially and the C binaries (simoffsets/cullst/
                        siminsar/intfloat) are OpenMP-threaded via ompThreads.
setupNISARTracks --copyFiles   prime new tracks: instantiate tie_plan_header/
                        vel_thumb_plan from the project templates (idempotent;
                        only matters the first time a track appears)
setupNISARTracks --year tie points / velocity for the year(s) of the new data
makeMaster              update Release/masterInput/inputFile
```

`--noDownload` skips `searchGranules`/`downloadGranules` and starts at
`consolidateStaging`, running the chain on the dated staging folders already
present. Useful for reprocessing without re-downloading, and for testing (simulate
"days" by pre-creating dated folders). **Note:** the newest dated folder always
plays the "today" role and is withheld from release — so to force manually-staged
data through, either stage it in a folder dated older than an empty folder dated
today, or just wait for the next real run.

## Running as a cron job

cron starts with an almost-empty environment, so it can't just call
`autoupdateNISAR` — the conda env (autoupdateNISAR/searchASF/ariaDownload/SetupNISAR
+ GDAL libs) and the GrIMP C binaries (`~/bin/x86_64`) aren't on `PATH`. The wrapper
`scripts/runAutoupdateNISAR.sh` (in this repo) sets both up, then runs the update for
one project directory. Point cron at it by absolute path:

```cron
# daily 22:00; flock skips a run if the previous is still going; log the console catch-all
0 22 * * * /bin/flock -n /Volumes/insar4/ian/Antarctica2026/Antarctica40/.autoupdate.lock \
  /home/ian/PycharmProjects/packages/nisargrimpworkflow/scripts/runAutoupdateNISAR.sh \
  /Volumes/insar4/ian/Antarctica2026/Antarctica40 \
  >> /Volumes/insar4/ian/Antarctica2026/Antarctica40/dailyDownloadLogs/cron.log 2>&1
```

Add one line per project (Antarctica40, Antarctica80, …) with its own project dir.
Notes:

- **Lock file in the project dir, not `/tmp`.** The project trees live on shared
  storage, so a lock under the project dir is honored regardless of which machine
  the cron runs on (crontabs are per-machine, so a project migrated between
  machines — or accidentally scheduled on two — can't double-run). Each project's
  lock is naturally distinct since it's inside that project dir. Run on the machine
  where the project disk is *local* (e.g. its NFS server) for best throughput; the
  lock still works because it's a normal file on that shared tree.

- The wrapper takes the project dir as its first argument; extra args pass through
  (e.g. append `--noDownload` to process manually-staged data instead of searching).
- `CONDA_BASE` and `BINDIR` env vars override the machine defaults
  (`/home/ian/miniforge3`, `/home/ian/bin/x86_64`).
- It aborts cleanly if `<projectDir>/autoupdate.yaml` is missing (e.g. an NFS mount
  not yet up), so a half-mounted tree doesn't cause a spurious run.
- Test it by hand first — this catches PATH/mount problems before cron does:
  `scripts/runAutoupdateNISAR.sh /Volumes/.../Antarctica40`.
- The `cron.log` redirect is a catch-all for search/download/staging output and
  tracebacks; the per-run `summary_<date>.txt` files (below) are the primary record.

## Run summary logs

After each run, `processReleased` writes a human-readable summary to
`<projectDir>/dailyDownloadLogs/summary_<MM-DD-YYYY>_<HHMMSS>.txt` (the directory
is created if absent; the summary is written even if a step fails). It opens with
an overview — orbits succeeded/failed, frames processed this run, granules released
(RUNW/ROFF split), and the pass/fail status of `FileNISARProducts`,
`setupNISARTracks`, and `makeMaster` — followed by a per-orbit breakdown (failures
first) giving each orbit's track, the frames it processed this run, and its status.
For a failed orbit it also records the return code and the tail of that orbit's
output as the failure reason. Each orbit's full stdout+stderr is captured to
`dailyDownloadLogs/orbitLogs_<MM-DD-YYYY>_<HHMMSS>/orbit_<orbit>_track-<track>.log`
(also avoids the interleaving you'd get from up to `nThreads` orbits printing to
one console at once).

## Orbit derivation (no HDF5 read)

The absolute reference orbit is a linear function of cycle and track:

```
orbit1 = 173*cycle + track + 618        # 173 orbits per NISAR repeat cycle
```

`orbitFromCycleTrack(cycle, track)` (in `FileNISARProducts.py`) computes this, so
the released granule filenames map directly to the `track-<track>/<orbit1>_<frame>/`
directories that `FileNISARProducts` created — no HDF5 read needed. A single
`(cycle, track)` is one orbit even when an Antarctic pass spans both descending and
ascending frames across the pole. `orbitsFromReleased()` warns and skips any
`(cycle, track)` whose `<orbit1>_*` frame dirs are unexpectedly absent.

## Virtual-frame rebuild when late frames arrive

`SetupNISAR` is run **without `--new`**. When late frames extend an existing
virtual frame, this rebuilds it in place rather than fragmenting or skipping:

- `assignVirtualFrameNumbers()` reads the orbit's own `frames.txt` and reassigns
  the grown contiguous group back to the **same** `_0000`;
- per-frame RUNW/ROFF/ionosphere products self-skip for frames already converted,
  so only the new frames are processed;
- virtual-frame assembly is ungated (except by `--new`), so `_0000` is rebuilt from
  all frames and `frames.txt` rewritten.

Example: `_0000 = [135,136,137,138,139]`; frame 140 arrives → `_0000` is rebuilt as
`[135..140]` (no `_0001`, frames 135–139 not reconverted). A *non-contiguous* late
frame (a gap) correctly becomes its own virtual frame instead.

## autoupdate.yaml keys

| Key | Meaning |
|---|---|
| `archiveDir` | Download/staging + general-population (RUNW/ROFF) root |
| `region` | `Antarctica` / `Greenland` → searchASF spatial flag |
| `lookBack` | Search window length in days (today-lookBack … today) |
| `products` | Product types to search/file (e.g. `RUNW ROFF`) |
| `bandwidths` | searchASF `--bandwidth` tokens (GrIMP `80` → 77 MHz) |
| `projectDir` | (optional) project root for the processing steps; default = the config file's directory |
| `setupFlags` | (optional) extra flags forwarded to `SetupNISAR` per orbit; default empty. Region-dependent processing choices are NOT put here — see `sepIceRock` below |
| `nThreads` | (optional) number of `SetupNISAR` orbit jobs to run in parallel in step 5; default 20 |
| `ompThreads` | (optional) OpenMP width passed to each orbit's `SetupNISAR` (`-ompThreads`); default = SetupNISAR's own (6). Balance `nThreads * ompThreads` against core count |

Region-dependent processing choices live in **`project.yaml`**, read by `SetupNISAR`
itself, not in `setupFlags`. In particular `sepIceRock` (the ice-anchored/rock-seeded
ionosphere path + global fill) is a `project.yaml` key — `True` for Greenland, `False`
for Antarctica for now. `SetupNISAR`, `processTrack` (via SetupNISAR), and a standalone
`estimateIonosphere` all honor it; the `--sepIceRock` CLI flag only forces it on.

The driver passes `--inputPath <archiveDir>` to `FileNISARProducts` explicitly, so
the run does not depend on `projectDir` containing an `autoupdate.yaml` with the
same `archiveDir` (standalone `FileNISARProducts` still falls back to that file
when `--inputPath` is omitted).

## Edge cases / notes

- **Group-merge orphan:** if late frames bridge two previously-separate contiguous
  groups into one, the smaller old `_00K` virtual-frame dir is orphaned (its frames
  now live in `_0000`); cleanup of the orphan is not automated.
- **Idempotent re-runs:** consolidate/release moves are per-file skip-if-exists; a
  re-run whose prior folders are already released releases nothing.
- **Skipped cron runs:** normally only two dated folders exist (yesterday + today)
  and the released one is removed. If runs were skipped, several accumulate; the
  next run sweeps every non-today folder into the archive and reconciles today
  against the union of them all, so nothing is stranded.
- **Partial downloads:** a granule with an aria2c `.aria2` control file next to it
  is a partial. The retry loop deletes partials before re-attempting (ariaDownload
  would otherwise skip the existing file rather than resume), deletes any partials
  still incomplete after the last attempt (so tomorrow's search re-downloads them —
  a leftover partial would otherwise be dedup-indexed and never fetched again), and
  `consolidateStaging`/`releaseStaging` both refuse to count or release any `.h5`
  with a sibling `.aria2`.
- **Crash recovery (pending granules):** the released work list is persisted to
  `dailyDownloadLogs/pendingGranules.txt` before processing and cleared once the
  per-orbit SetupNISAR stage has run. If a run dies between release and
  processing (released granules exist only in RUNW/ROFF at that point), the next
  run merges the pending list back in — re-listing an already-processed granule is
  harmless (per-frame products self-skip; the virtual frame rebuilds cheaply).
- **Per-orbit logs stream to disk** as each orbit runs (not buffered), so a
  mid-run crash keeps everything logged so far and `tail -f` works during a run.
- **Stray files:** unparseable `.h5` names in staging folders are warned about and
  skipped (never crash the run); files with no RUNW/ROFF product token stay in
  their dated folder.
- **searchASF dedup covers staging:** the dated folders sit one level under
  `archiveDir`, so the existing `archiveDir/*/*` dedup glob prevents re-downloading
  a staged-but-not-yet-released pass.
