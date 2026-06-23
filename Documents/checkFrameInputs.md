# checkFrameInputs — Missing Upstream-Product Scan

## Overview

`checkFrameInputs` walks every virtual frame (`track-N/*_00NN/`) under a
project root and checks for three categories of required upstream product:

| Category | Files checked (any one match satisfies it) |
|---|---|
| `offsets` | `range.offsets.vrt`, `azimuth.offsets.vrt` |
| `phase` | `phase.uw.vrt` |
| `ioncorrection` | `*.ionosphereCorrection.globalFill.offset.vrt`, `*.ionosphereCorrection.offset.vrt` |

It's a read-only reporting tool — it never writes into `track-N/`. Not tied
to any specific project — it operates on whatever `track-N/*_00NN/`
directories it finds under `--projectDir` (default: current directory).

A frame missing a category means [SetupNISAR](SetupNISAR.md) didn't (or
couldn't) produce that product for it — usually because the per-frame
RUNW/ROFF conversion failed for that frame, or (for `ioncorrection`
specifically) the underlying coherence was too low for phase unwrapping to
produce any valid pixels at all, so `globalFillIonosphere()` correctly has
nothing to fill (`globalFillIonosphere: no valid pixels, skipping`) and never
writes the final `...ionosphereCorrection.globalFill.offset.vrt` — see the
track-16/3921_0002 case (frame 36: coherence mean 0.07, unwrapped phase
0/5,432,670 valid pixels).

A frame containing an `Exclude` marker file (any content, including empty —
same convention used by `mosaicworkflow`/`insarScripts` tooling, e.g.
`checkinput.py`'s `addExcludeIfNeeded()`) is skipped entirely and counted
separately, not reported as missing — it's been reviewed and intentionally
marked unusable, not silently broken.

## Usage

```bash
cd /path/to/project   # or pass --projectDir
checkFrameInputs
```

Output: one line per frame with any missing category, followed by a summary:

```
./track-98/1408_0020: missing ioncorrection

9 / 401 virtual frame(s) missing one or more required inputs (1 excluded via
an Exclude marker, not counted above)
  offsets: 0 frame(s) missing
  phase: 0 frame(s) missing
  ioncorrection: 9 frame(s) missing
```

## See also

[buildFrameGpkg](buildFrameGpkg.md) is a complementary, finer-grained check —
it looks at the *downstream* baseline/azimuth-estimation sidecar files
(`motion/baseline*.yaml`, `motion/rBaseline*.yaml`, `motion/az.est.const*`)
and their QC values (sigma, tiepoint counts), not the upstream
offsets/phase/ionosphere products this script checks.
