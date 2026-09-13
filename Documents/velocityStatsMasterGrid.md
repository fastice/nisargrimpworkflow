# Velocity stats on a master grid (`velocityStatsRegions: track`)

Opt-in replacement for the per-frame-number velocityStats bins
(`velocityStats/0000-0009`, `0010-0019`, ...). Default behaviour is unchanged; the
switch lives in the project's `project.yaml`:

```yaml
velocityStatsRegions: track     # absent or frameBins = original bins
velStatsThreads: 8              # optional: parallel velocityStats.py runs in track mode [8]
```

## What changes in track mode

| step | frameBins (default) | track |
|---|---|---|
| range dirs (`setupNISARTracks`) | one `velocityStats/<M>0-<M>9` per tens digit | one `velocityStats/0000-0099` per track; legacy bins must be migrated first (a warning skips the track otherwise) |
| per-product mosaic box (`vel_thumbs`) | the bin header's `resolution` line, frozen when `makevelstatsregions` last ran | computed per product from its own `offsets.velocity.ll.vrt`: true extent, 3 km pad, origin and far edge snapped to integer km, posting from the header (0.2 km) |
| `makevelstatsregions.py` | rewrites each header's box | no-op: reports each product's origin/posting, never touches the header |
| composite (`velocityStats.py --doRA`) | all members must share one grid (first product's); others skipped | members placed on the union master grid by integer pixel offset; the same two-pass findKeepers/sumStats/normStats run over strips of 8 Mpx; output written tile by tile (DEFLATE, sparse); `products.vrt` indexes the pieces |
| comparison (`autocleanNISAR.py`) | reads the whole bin composite | scans `velocityStats/*-*` for the range holding the frame, reads only the window the frame covers (`geoimage.readDataWindow`), and falls back to the reference-only blend when no composite exists (mode-independent) |
| `checkvelsize.py` | equal grids per range | integer-km origin + common posting per product |

The master grid is defined by the region projection (EPSG:3031 for Antarctica), the
0.2 km posting from the thumb header, and integer-km lower-left pixel-centre origins.
Every existing product already satisfies this (all 4812 checked on 2026-09-11), so no
resampling is ever needed.

## Migration of a project

1. Set `velocityStatsRegions: track` in `project.yaml`.
2. `setupNISARTracks --migrateVelStatsRegions [--tracks track-N ...] [--dryRun]`
   moves `velocityStats/<bins>` to `velocityStats/old/` and the per-bin
   `tiepoints/vel_thumb_header_*dash*`, `tie_plan20??.*dash*`,
   `vel_thumb_plan_all20??dash*` to `tiepoints/old/`, then creates the track-wide
   range and its header (sentinel resolution). Nothing is deleted.
3. Regenerate products and stats exactly as `runUpdate` does today:
   `setupNISARTracks --runVelstatsregions --noUpdateVelStats [--tracks ...]` then
   `setupNISARTracks --updateAutoclean [--tracks ...]`.

Tracks not yet migrated keep working: their bins are assembled on the master grid
per bin and autoclean finds them by scanning.

## Also fixed on 2026-09-11 (both modes)

`vel_thumbs` now passes `-legacyPairRange` to mosaic3d. A thumb is a one-image
mosaic; mosaic3d's joint crossing-orbit solver (default since 2026-08-29) is
ill-conditioned for a single image and its condition gate rejected essentially every
pixel, so every thumb built between 2026-08-31 and 2026-09-11 was empty and the
stats/autoclean chain silently ran on nothing. Thumbs built since must be regenerated
(`--runVelstatsregions --noUpdateVelStats`, then `--updateAutoclean`).

## Validation (2026-09-11, Antarctica40/track-82, scratch copy)

- Composite from the existing products: identical to bin mode pixel for pixel
  (navg, vr, va, er) on the two populated bins, except 5 pixels where a product
  from a neighbouring bin now legitimately contributes.
- `readDataWindow` vs full read + interpolation at 1000 random points: max
  difference 1e-9. Autoclean masks bit-identical to the old code on the same
  composite; 5.5 s and 1.4 GB instead of 55 s and 7.8 GB.
- Regenerated 14 products with per-product boxes: all integer-km origins; 55 M
  valid pixels against 26 M in the clipped/empty originals (5025_0000: 1.2 M to
  10.1 M; the three 5544 products from 0 to 3-11 M each).
- Track composite 24075 x 13900 px in 42 strips: 8.3 min, 12.7 GB peak (the
  regional velocity map alone is 1.2 GB plus interpolators), hence the
  `velStatsThreads` cap.

Backups of every touched file carry the suffix `.preMasterGrid-20260911`
(mosaicworkflow: velocityStats, autocleanNISAR, makevelstatsregions, checkvelsize;
insarworkflow: velThumbs; utilities: geoimage; nisargrimpworkflow: setupNISARTracks).
Plan and evidence: `Release/Velocity/gateAbsScoring/masterGridPlan-2026-09-11.md`.
