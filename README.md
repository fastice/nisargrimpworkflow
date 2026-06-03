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

See `Documents/ROFFtoGrimp.md` and `Documents/RUNWtoGrimp.md` for detailed
usage and output descriptions.

## Dependencies

- `nisarhdf`
- `sarfunc`
- `utilities`
- `numpy`, `rioxarray`, `gdal`
