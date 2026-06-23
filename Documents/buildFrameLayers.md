# buildFrameLayers — QGIS Layer Definition for Frame QC

## Overview

`buildFrameLayers` reads the per-cycle GeoPackages written by
[buildFrameGpkg](buildFrameGpkg.md) and writes a single QGIS Layer Definition
file (`.qlr`) that organizes them into a ready-to-load layer tree:

```
Frames
├── ascending
│   ├── Cycle 3
│   ├── Cycle 4
│   └── ... one layer per cycle
└── descending
    └── ... same per-cycle layers

rBaseline
├── ascending
│   └── ... same per-cycle layers
└── descending
    └── ... same per-cycle layers

sigmaRBaseline > 0.5
├── ascending   ← all cycles combined into one layer, no per-cycle split
└── descending  ← all cycles combined into one layer, no per-cycle split
```

`Frames` and `rBaseline` both reference the exact same underlying GeoPackage
layers — they're just two different ways of coloring the same data:

- **`Frames`** — every frame is drawn as an outline only (no fill), colored
  continuously by whichever sigma field you choose, via a QGIS **project**
  variable named `sigma_field` (one of `sigmaBaseline`, `sigmaRBaseline`,
  `sigmaAz`; defaults to `sigmaRBaseline`). Set it once under **Project
  Properties → Variables**, and every layer in this group recolors together —
  no reload needed to switch fields. Low sigma = green, high = red
  (ColorBrewer `RdYlGn`), scaled across that field's actual min/max at
  generation time.
- **`rBaseline`** — fixed, discrete 10-class coloring of `sigmaRBaseline`
  specifically, in 0.1 increments (`< 0.1`, `0.1 - 0.2`, … `0.8 - 0.9`,
  `>= 0.9` — the top class is open-ended so it also catches values above 1).
  Since each class has its own static color (no data-defined expression), the
  QGIS legend swatches are accurate for this group, unlike `Frames`.

The third group is structured differently — **not** organized by cycle at
all:

- **`sigmaRBaseline > <thresh>`** (threshold in the group name itself, set
  via `--offsetsSigmaThresh`, default `0.5`) — one flat layer per direction
  (`ascending`, `descending`) flagging every frame, across *all* cycles, where
  `sigmaRBaseline` exceeds the threshold. Built as an
  [OGR VRT Union Layer](https://gdal.org/drivers/vector/vrt.html) combining
  that direction's layer from every `cycle<NN>.gpkg` into one virtual layer
  (written alongside the GeoPackages as `allAscending.vrt`/
  `allDescending.vrt`), with a QGIS subset filter (`"sigmaRBaseline" >
  <thresh>`) applied so only the flagged frames are visible. Fixed solid red
  outline, no per-feature coloring — membership in this group already tells
  the whole story.

Not tied to any specific project — it just needs a directory of
`cycle<NN>.gpkg` files in the layout `buildFrameGpkg` produces. All paths
embedded in the `.qlr` are resolved to absolute paths at generation time, so
the result is valid no matter what directory QGIS is run from.

---

## Usage

```
buildFrameLayers [--projectDir DIR] [--inventoryDir DIR] [-o OUTPUT]
```

Run [buildFrameGpkg](buildFrameGpkg.md) first — this reads the per-cycle
GeoPackages it writes.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--projectDir DIR` | `.` (cwd) | Used only to derive `--inventoryDir`'s default |
| `--inventoryDir DIR` | `<projectDir>/frameInventory` | Directory containing `cycle*.gpkg` (i.e. `buildFrameGpkg`'s output) |
| `-o, --output PATH` | `<inventoryDir>/frames.qlr` | Output `.qlr` path |
| `--offsetsSigmaThresh VAL` | `0.5` | `sigmaRBaseline` threshold for the flat `sigmaRBaseline > VAL` flag group |

### Examples

```bash
# Typical pairing, run from the project root
cd /path/to/some/Project
buildFrameGpkg
buildFrameLayers

# Pointing at an inventory directory built elsewhere
buildFrameLayers --inventoryDir /tmp/inventory -o /tmp/inventory/frames.qlr

# Flag frames above a stricter threshold
buildFrameLayers --offsetsSigmaThresh 0.3
```

### Loading the result in QGIS

**Layer → Add Layer → Add Layer Definition File...** → pick the `.qlr`. Then,
to use the `Frames` group's switchable coloring: **Project Properties →
Variables → Add Variable**, name `sigma_field`, value `sigmaBaseline` /
`sigmaRBaseline` / `sigmaAz`, Apply. Changing that one variable later
recolors every layer in the `Frames` group at once.

Re-run after every `buildFrameGpkg` run — it recomputes the cycle list and
each field's min/max from whatever's currently in `--inventoryDir`, so the
color scale always reflects the current data.
