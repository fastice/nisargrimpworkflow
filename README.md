# nisargrimpworkflow

Command-line tools for converting NISAR Level-2 HDF5 products into the binary
flat-file formats used by the GrIMP velocity mosaic pipeline.

## Programs

**ROFFtoGrimp** — Converts a NISAR ROFF (range/azimuth offset) HDF5 product
into GrIMP-format offset maps and VRT products.

**RUNWtoGrimp** — Converts a NISAR RUNW (unwrapped interferogram) HDF5 product
into GrIMP-format phase, coherence, and ionosphere-correction VRT products.

**SetupNISAR** — Orchestrates per-frame RUNW/ROFF conversion and consolidates
multi-frame products into a virtual-frame VRT mosaic.

**searchASF** — Search the ASF DAAC for Sentinel-1 IW SLC products within a
date range and spatial area.

## Documentation

- [setupNISARTracks](Documents/setupNISARTracks.md) — orchestrate per-track NISAR processing
- [SetupNISAR](Documents/SetupNISAR.md) — per-frame RUNW/ROFF conversion and VRT mosaicking
- [processTrack](Documents/processTrack.md) — per-track processing workflow
- [ROFFtoGrimp](Documents/ROFFtoGrimp.md) — convert NISAR ROFF HDF5 to GrIMP offset maps
- [RUNWtoGrimp](Documents/RUNWtoGrimp.md) — convert NISAR RUNW HDF5 to GrIMP phase/coherence products
- [estimateIonosphere](Documents/estimateIonosphere.md) — ionospheric range-offset correction
- [FileNISARProducts](Documents/FileNISARProducts.md) — NISAR product file management

## Dependencies

- `nisarhdf`
- `sarfunc`
- `utilities`
- `numpy`, `rioxarray`, `gdal`
