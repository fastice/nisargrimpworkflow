# buildFrameGpkg — Virtual-Frame QC/Inventory GeoPackages

## Overview

`buildFrameGpkg` walks every virtual frame (`track-N/*_0000/`) under a project
root and writes one GeoPackage per NISAR cycle, each containing two layers
(`ascending`, `descending`). Each feature is one virtual frame, with its
footprint geometry plus baseline/range-baseline/azimuth QC values pulled from
the sidecar files the rest of the pipeline ([SetupNISAR](SetupNISAR.md) /
[ROFFtoGrimp](ROFFtoGrimp.md) / [RUNWtoGrimp](RUNWtoGrimp.md) / tie-point
processing) already produces in each frame directory.

It's a read-only reporting tool — it never writes into `track-N/`. Pair it
with [buildFrameLayers](buildFrameLayers.md), which turns its output into a
QGIS Layer Definition for browsing/QC in QGIS.

Not tied to any specific project — it operates on whatever `track-N/*_0000/`
directories it finds under `--projectDir` (default: current directory), so it
works for any project laid out this way.

---

## Required inputs per virtual frame

For each `track-N/OOOO_0000/` directory, all of the following are **required**.
If any is missing, unreadable, or internally inconsistent, that frame is
skipped — loudly, with a one-line reason printed to stdout — rather than
silently filled in with a fallback or null. A skipped frame means that frame
still needs (re)processing upstream, not that this tool failed.

| File | Used for |
|---|---|
| `frames.txt` | Space-separated list of physical sub-frame numbers merged into this virtual frame → `frames` column |
| `geodat*.geojson` (excluding `*.secondary.geojson`) | Footprint geometry (lon/lat) and `properties.ImageName`, which is parsed for `cycle`/`direction` |
| `motion/baseline.<R>x<A>.yaml` | Phase baseline: `sigma`, `nTiepoints`, `nTiepointsGiven` |
| `motion/rBaseline.deltabp.yaml` (also accepts the older `rBaseline.deltab.yaml`) | Range baseline: `sigma`, `nTiepointsUsed`, `nTiepointsGiven`, `sigmaWithoutIonCorrection`, `usingIon` |
| `motion/az.est.const.yaml` (preferred) or `motion/az.est.const` (legacy text) | Azimuth estimate: sigma + tiepoint counts. Tries the yaml form first; falls back to parsing the legacy `;`-comment text format if only that exists |
| `*.pairinfo` (e.g. `1830.2003.pairinfo`) | `orbit1`/`orbit2` — cross-checked against the directory name's orbit; a mismatch is treated as an error (skipped), not silently trusted |

`cycle` and `direction` come from `ImageName` via
[`FileNISARProducts.parseFileName`](FileNISARProducts.md), the same NISAR
filename parser the rest of the package uses — not a separate one-off regex.

---

## Output schema

One GeoPackage per cycle, named `cycle<NN>.gpkg` (2-digit zero-padded), each
with layers `ascending` and `descending`. Geometry is reprojected from the
source geodat's WGS84 lon/lat to **EPSG:3413** (NSIDC North Polar
Stereographic) — appropriate for Arctic/Greenland projects; an Antarctic
project would need EPSG:3031 instead (not currently a flag).

| Field | Type | Source |
|---|---|---|
| `track` | Integer | `track-N` directory name |
| `orbit` | Integer | `OOOO` from the `OOOO_0000` directory name |
| `virtualFrameId` | String | The `OOOO_0000` directory name itself |
| `orbit1`, `orbit2` | Integer | `*.pairinfo` |
| `cycle` | Integer | `ImageName` |
| `direction` | String | `ascending`/`descending`, from `ImageName` |
| `frames` | String | Comma-joined sub-frame numbers from `frames.txt` |
| `sigmaBaseline` | Real | `motion/baseline.*.yaml` → `sigma` (radians) |
| `sigmaRBaseline` | Real | `motion/rBaseline.deltabp.yaml` → `sigma` (meters) |
| `sigmaAz` | Real | `motion/az.est.const(.yaml)` → sigma (meters) |
| `sigmaRBaselineWithoutIon` | Real | `rBaseline.deltabp.yaml` → `sigmaWithoutIonCorrection` |
| `usingIonRBaseline` | Integer (0/1) | `rBaseline.deltabp.yaml` → `usingIon` |
| `nTiepointsBaseline`, `nTiepointsGivenBaseline` | Integer | `baseline.*.yaml` |
| `nTiepointsUsedRBaseline`, `nTiepointsGivenRBaseline` | Integer | `rBaseline.deltabp.yaml` |
| `nTiepointsUsedAz`, `nTiepointsGivenAz` | Integer | `az.est.const(.yaml)` |

---

## Usage

```
buildFrameGpkg [--projectDir DIR] [-o OUTPUTDIR]
```

Run [buildFrameLayers](buildFrameLayers.md) afterward, pointed at this
command's output directory, to turn the result into a QGIS `.qlr` layer tree.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--projectDir DIR` | `.` (cwd) | Root directory containing `track-N/` subdirectories |
| `-o, --outputDir DIR` | `<projectDir>/frameInventory` | Where to write `cycle<NN>.gpkg` files |

### Examples

```bash
# Run from inside the project root
cd /path/to/some/Project
buildFrameGpkg

# Run from elsewhere, pointing at the project explicitly
buildFrameGpkg --projectDir /path/to/some/Project -o /tmp/inventory
```

Output ends with a summary, e.g.:
```
Output dir: ./frameInventory  (EPSG:3413)
  cycle07.gpkg: 7 ascending, 9 descending
  ...
311 frame(s) written, 2 of 313 skipped
```

Re-run any time after reprocessing — it always rewrites `cycle<NN>.gpkg` from
scratch (`os.remove` then recreate), so there's no stale-data risk between runs.

The program itself prints a reminder of the exact `buildFrameLayers
--inventoryDir <outputDir>` command to run next.
