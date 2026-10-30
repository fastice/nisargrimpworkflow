# autoupdate — daily incremental NISAR update driver

`autoupdateNISAR` (`autoupdate.py`) is a cron-oriented driver that searches ASF
for new NISAR granules, downloads them, and runs the full GrIMP processing chain
on the newly-completed data. It is configured entirely by an `autoupdate.yaml`
that normally sits beside the project's `project.yaml`.

## Setup checklist

One-time setup to bring a project onto the daily driver. Detail for each item is
in the sections below; this is the ordered "did I do everything" list.

1. **Conda `base` env has the tools.** `autoupdateNISAR`, `searchASF`,
   `downloadNISARoptimized`, `ariaDownload`, `SetupNISAR` (and GDAL) must all be importable/on `PATH` in the
   env the wrapper activates (`base`). Install/update `nisargrimpworkflow` into it.
2. **GrIMP C binaries built and installed** to `~/bin/x86_64` (`simoffsets`,
   `cullst`, `siminsar`, `intfloat`, `mosaic3d`). Override the location with
   `BINDIR` if elsewhere.
3. **Earthdata credentials in place** for the downloads (standard `~/.netrc`
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
downloadGranules        RUNW/ROFF/GUNW/GOFF slim via downloadNISARoptimized, anything else via
                        ariaDownload, into archiveDir/<MM-DD-YYYY>/  (today's staging folder)
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
setupNISARTracks --year tie points / velocity for the year(s) of the new data,
                        restricted with --tieTracks to the tracks that got new
                        granules (plus any whose autoclean mask was rebuilt)
makeMaster              update Release/masterInput/inputFile (or the master named
                        by masterInputPaths/masterOutputPath)
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

## Which tracks get tiepointed

Step 6 used to run `setupNISARTracks --year` with no track restriction, and
`setupNISARTracks` tiepoints every track it discovers. So a night that filed 8 new
orbits still re-tied all ~170 tracks: with `-nThreads 32` that is a 1.5–2 hour stage,
and for a track with no new frames it recomputes the identical answer, because ties and
baselines are per-track quantities that depend on nothing outside their own track.

Step 6 now passes `--tieTracks` naming the tracks with new granules — the same
`orbitsFromReleased()` list that drove step 5. Two things are deliberately *not*
restricted:

- The **autoclean-mask scan** still covers every track. A track can need a first-time
  mask without having new data this run, and `setupNISARTracks` adds any track whose
  mask it (re)builds back into the tie list, since a fresh mask changes the tie fit.
- **With no orbits at all** (granules released but no `track-<N>/<orbit>_*` dirs — see
  the `orbitsFromReleased` warning), the flag is omitted and every track is tied. The
  conservative default: an anomalous run should not silently skip the tie stage.

The summary's `Tracks tiepointed:` line lists the restricted set; its absence means all
tracks were tied. Set `retieAllTracks: true` in `autoupdate.yaml` to force the old
behaviour without a code change, or run `setupNISARTracks --year <yyyy>` by hand — with
no `--tieTracks`, that still tiepoints everything.

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

Between the overview and the per-orbit breakdown sits a **`TRACK ERRORS`** section
aggregating what the tie stage wrote to two sets of log files:

| path | holds |
|---|---|
| `track-N/stderr` | the per-track tie chain (`maketies` → `makeframetie.py`) |
| `track-N/<frame>/motion/*stderr` | per-frame `azparams`/`rparams` (`azest_stderr`, `estr_stderr`) |

Four error formats are recognised: `myerror` banners, argparse rejections, traceback
final lines, and the C binaries' own `error()` output — `clib/standard.c` prints a bare
`error:  <msg>` with no colour and no program prefix, so it matches none of the others
and needs its own pattern. Messages are normalised (directory prefixes stripped, long
digit runs collapsed, parenthesised counts collapsed to `(N)`) so one failure hitting
every track becomes a single line with a count and a few example tracks, rather than
hundreds of near-identical ones. The parenthesised-count rule is what merges
`tiepoints: Insufficient Number (0|1|2|3) of Valid tie points` — one failure mode, four
spellings — into a single kind; the exact count is still in `track-N/stderr`. Only files modified since the run started are read, so a run
that aborts before the tie stage cannot report yesterday's errors as today's. It is
capped at 15 distinct messages × 5 example tracks to keep the mail bounded.

### Screening out errors you have already seen

Several standing failures are benign and permanent. The commonest is
`tiepoints: Insufficient Number (0) of Valid tie points`: a frame whose phase never
ties, where `tiepoints` writes the `sigma: -1` no-solution sentinel into
`motion/baseline.NxM.yaml`, `mosaic3d` skips the phase for that frame, and the frame
still contributes through its (perfectly good) range offsets. Nothing is broken and
nothing needs doing — but reported in full every night it buries the errors that do.

So the section splits into **`NEW`** and **`RECURRING`**, keyed on the
(normalised message, track) pair:

```
TRACK ERRORS  (2 new, 3 recurring; 40 track(s) affected)
  NEW
  [   1] readgeodat : could not open file
         track-12
  [   1] tiepoints: Insufficient Number (N) of Valid tie points
         track-107
  RECURRING  (already reported; per-track detail in dailyDownloadLogs/knownTrackErrors.json)
  [  19] tiepoints: Insufficient Number (N) of Valid tie points   since 08-24-2026
  [   8] svInitAzParams : insufficient points to fit              since 08-24-2026
  [   6] No output because no data in range, check dates          since 08-25-2026
```

New pairs get the full message + example-tracks block; recurring ones collapse to one
line carrying the count and the date that message was first seen. A known message
spreading to a **new** track appears in both sections — that track is news even though
the message is not.

The state lives in `dailyDownloadLogs/knownTrackErrors.json`
(`{message: {track: first-seen}}`), rewritten every run:

- a pair still failing keeps its original first-seen date;
- a pair is **forgotten** only when that track's log was actually re-read this run and
  no longer carries the message — so it will be reported as new if it ever comes back;
- a pair from a track the run never touched is **carried forward** untouched.
  `collectTrackErrors` only reads logs rewritten since the run started, so a track with
  no new granules contributes nothing this run; pruning on that absence would make it
  report "new" again the next time it genuinely fails.

Delete the file to re-report everything once; an unreadable or corrupt file reads as
empty, which does the same (over-report rather than silently swallow). The first run
after this feature was added therefore reports the whole standing backlog once, and is
quiet from the next run on.

`notifyOnTrackErrors` now triggers the mail on **new** errors only — a night whose only
errors are recurring sends nothing, while the summary file still records both sections
in full. A failed step, a failed orbit, or a crash still mails regardless.

This section is the reason the whole thing exists: nothing was reading either set of
logs, so a tie stage that failed on every track for eleven consecutive nights (see
*Pointing the nightly master at a multi-project master* for the related master-path bug)
produced summaries where every step still said `OK`. Aggregated, that failure is two
lines. The per-frame logs matter for the same reason: when `azparams` dies it leaves a
zero-byte `az.est.yaml`, which `mosaic3d` reads as all zeros and then fails on obscurely
much later — the useful message (`svInitAzParams : insufficient points to fit`, naming
the frame) only ever existed in `motion/azest_stderr`.

Nothing here changes any failure behaviour — every one of these steps already exited
non-zero (C `error()` exits 1, `myerror` now does too). They were simply never reported.
There is still exactly **one** mail per run: `notifyOnErrors` sends the whole summary, so
the new section rides along in it rather than adding another message.

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
| `slimDownload` | (optional, default `true`) fetch RUNW/ROFF (and GUNW/GOFF, if listed in `products`) slim with `downloadNISARoptimized` -- only the bands the GrIMP workflow reads, float32 rounded to 10 mantissa bits (RUNW `unwrappedPhase` 12), same internal paths, ~4-5x smaller and about half the transfer. `false` fetches whole products with `ariaDownload`, e.g. if ASF range requests misbehave |
| `bandwidths` | searchASF `--bandwidth` tokens (GrIMP `80` → 77 MHz) |
| `projectDir` | (optional) project root for the processing steps; default = the config file's directory |
| `setupFlags` | (optional) extra flags forwarded to `SetupNISAR` per orbit; default empty. Region-dependent processing choices are NOT put here — see `sepIceRock` below |
| `nThreads` | (optional) number of `SetupNISAR` orbit jobs to run in parallel in step 5; default 20 |
| `ompThreads` | (optional) OpenMP width passed to each orbit's `SetupNISAR` (`-ompThreads`); default = SetupNISAR's own (6). Balance `nThreads * ompThreads` against core count |
| `retieAllTracks` | (optional, default `false`) tiepoint **every** track in step 6 instead of only the ones with new data. See *Which tracks get tiepointed* below |
| `masterInputPaths` | (optional) directories `makeMaster` scans for `track-N` subdirs (its `--inputPaths`); a single value, a space-separated string, or a YAML list. Relative entries are anchored to `projectDir`. Default: `projectDir` alone |
| `masterOutputPath` | (optional) directory `makeMaster` writes `inputFile` into (its `--outputPath`), anchored to `projectDir` when relative; default `Release/masterInput` |
| `notifyOnTrackErrors` | (optional, default `true`) let aggregated per-track errors trigger the notification mail on their own — **new** ones only, per *Screening out errors you have already seen*. The `TRACK ERRORS` section is written to the summary file either way; set `false` to keep it out of the mail trigger entirely while you work through a known recurring error |
| `notifyEmail` | (optional) address mailed the run summary **only when the run hits errors** (a failed step, a failed `SetupNISAR` orbit, or a crash). Silent on a clean run; with no key set, no mail is ever sent (opt-in). Best effort via the local MTA (`mail`/`mailx`) — a host with no working mailer just logs a warning and continues |

Region-dependent processing choices live in **`project.yaml`**, read by `SetupNISAR`
itself, not in `setupFlags`. In particular `sepIceRock` (the ice-anchored/rock-seeded
ionosphere path + global fill) is a `project.yaml` key — `True` for Greenland, `False`
for Antarctica for now. `SetupNISAR`, `processTrack` (via SetupNISAR), and a standalone
`estimateIonosphere` all honor it; the `--sepIceRock` CLI flag only forces it on.

### Pointing the nightly master at a multi-project master

By default step 7 rewrites `<projectDir>/Release/masterInput/inputFile`. That is the
file the mosaics read only when the project dir *is* the mosaic root, as for
Greenland. Antarctica splits the season by bandwidth (`Antarctica2026/Antarctica40`,
`Antarctica2026/Antarctica80`), each with its own `autoupdate.yaml`, while the
mosaics read one combined master at `Antarctica2026/Release/track-all/inputFile`.
Without these keys the nightly run happily refreshes `Antarctica40`'s own master and
the combined one silently goes stale. In `Antarctica40/autoupdate.yaml`:

```yaml
masterInputPaths:              # scan both bandwidth trees (relative to projectDir)
  - .
  - ../Antarctica80
masterOutputPath: ../Release/track-all
```

List **every** tree the master should cover, not just the ones being updated. No new
80 MHz data is being collected and `Antarctica80` has no `autoupdate.yaml`, but its
already-processed frames are the majority of the combined master (2471 of 4380 data
lines as of 2026-08-24); leaving it out of `masterInputPaths` would quietly drop
them from the mosaic input. A static tree simply re-emits the same lines each night.

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
- **Slim downloads:** `downloadNISARoptimized` writes each granule atomically
  (`.h5.partial`, renamed when complete), so a present `.h5` is complete; a killed
  repack's `.partial` is removed with the aria partials. After each pass
  `collectSlimManifest` deletes any granule that failed verification (so it is
  retried, then re-listed tomorrow), appends the staging folder's
  `slimManifest.jsonl` to `archiveDir/slimManifest.jsonl` and removes an empty
  `.scratch/` -- `releaseStaging` only moves `*.h5` and removes a folder once it is
  empty, so either would otherwise strand the folder.
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
