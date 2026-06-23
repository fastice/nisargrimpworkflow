# makeMaster — Assemble the Master mosaic3d inputFile

## Overview

`makeMaster` walks every `track-N/tiepoints/track-N/inputFile-N-all20??`
file under a project root, extracts their data lines (skipping `;` comments
and pure-number lines), and assembles them into a single master `inputFile`
that `mosaic3d` reads for the velocity mosaic.

Before writing the output, it verifies that every file path referenced in
those data lines actually exists — including, for each `range.offsets`
token, any ionosphere-correction file embedded in that frame's
`range.offsets.vrt` `ionosphereRangeOffsetCorrection` metadata (the same
mechanism `SetupNISAR.globalFillIonosphere()` writes and
`rparams`/`checkForIonosphereCorrection()` reads). Missing files are
reported as warnings, one per missing path, plus a summary count — this
does **not** block writing the output (loud, not silent: a track with
missing inputs still gets assembled, so you can see exactly what's
incomplete rather than have the whole run abort).

Originally a one-off script hardcoded to `newGreenlandProject`
(`/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/makeMaster.py`,
with `main()` called unconditionally at import time — moved into the
package, generalized via `--projectDir`, and given a proper
`if __name__ == '__main__':` guard so importing it no longer has side
effects).

## Usage

```bash
cd /path/to/project   # or pass --projectDir
makeMaster                              # writes <projectDir>/Release/masterInput/inputFile
makeMaster --outputPath /some/other/dir # writes /some/other/dir/inputFile directly
```

| Flag | Default | Description |
|---|---|---|
| `--projectDir` | `.` | Root directory containing `track-N` subdirectories |
| `--outputPath` | `Release/masterInput` | Output directory (relative to `--projectDir` unless absolute) |

Output:

```
Found tracks: [1, 11, 14, ...]
	 ---- track-170 inputFile-170-all2026: missing /path/.../range.offsets ----- 
754 referenced file(s) missing — see warnings above
Wrote 377 data lines to Release/masterInput/inputFile
```

or, when everything resolves:

```
Found tracks: [1, 11, 14, ...]
All referenced files verified present.
Wrote 377 data lines to Release/masterInput/inputFile
```

## Bare-name vs `.vrt` resolution

`inputFile-N-all20??` data lines reference `range.offsets`/`azimuth.offsets`
by their legacy bare (no-extension) name unconditionally — the
tie-script template that generates them doesn't check which form actually
exists on disk. Projects that have moved to VRT-only storage (no flat
binary file written at all, just `range.offsets.vrt`) still work, because
`rparams`/`azparams` (`getROffsets.c`) auto-detect and fall back to the
`.vrt` sidecar when given the bare name. `makeMaster`'s existence check
matches this: a token is considered present if either the bare path or
`<path>.vrt` exists.

## A track silently missing from the master inputFile

If a track has zero matches for `inputFile-N-all20??`, it's skipped with no
warning (same behavior as the original script) — check
`tiepoints/track-N/inputFile-N-all20??` exists and is non-empty for that
track if it's unexpectedly absent from the assembled output. A common
cause: per-year tie-plan segment naming that doesn't match the plain
`all20??` glob (e.g. `all2026dash0000dash0009`-style names from a segmented
tie plan) — not yet handled by this script.

## See also

[checkFrameInputs](checkFrameInputs.md) checks a different, earlier stage
of the pipeline (whether `SetupNISAR` produced the raw upstream
offsets/phase/ionosphere products per virtual frame at all) — `makeMaster`
checks whether the *tie-script-assembled* `inputFile` entries that
reference those products resolve correctly.
