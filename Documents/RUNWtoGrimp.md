# RUNWtoGrimp.py

Converts a NISAR RUNW HDF product (unwrapped interferogram) into GrIMP-formatted binary flat files and VRTs ready for input to `mosaic3d`. Retains only the largest connected phase component, optionally applies an ice mask, performs minor hole-filling, and writes the ionosphere phase screen alongside the unwrapped phase.

---

## Usage

```
RUNWtoGrimp [options] RUNW
```

`RUNW` is the path to a NISAR RUNW HDF5 file (local path or S3 URI).

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `RUNW` | — | RUNW HDF5 file to convert. |
| `--outputDir DIR` | RUNW directory | Directory for all output files. |
| `--referenceXML FILE` | None | Reference orbit XML file. Passed to the RUNW reader; used when orbit metadata is not embedded in the HDF. |
| `--secondaryXML FILE` | None | Secondary orbit XML file. |
| `--referenceOrbit N` | None | Reference orbit number override. Obsolete once embedded in HDF. |
| `--secondaryOrbit N` | None | Secondary orbit number override. |
| `--frame N` | None | Frame number override. |
| `--simMask` | False | Simulate and apply an ice mask (removes bedrock / non-ice areas). |
| `--simPhase` | False | Simulate interferometric phase from a DEM and velocity map. |
| `--region NAME` | auto | Region name (`greenland` or `antarctica`). Auto-detected from RUNW EPSG if not given. |
| `--interpThresh N` | 20 | Maximum hole area (pixels) to fill during phase interpolation. |
| `--islandThresh N` | 20 | Maximum isolated-region area (pixels) to discard during interpolation. |
| `--verbose` | False | Print all subprocess output to the terminal. |

---

## Processing flow

```
parseArgs()

openHDF(RUNW)            — load unwrapped phase, coherence, ionosphere screen,
                           and image geometry

cleanIonosphere()        — apply ionosphere correction (produces ionosphereCleaned
                           attribute if successful)

resolveRegion()          — determine region from EPSG if not set

maskPhase(largest=True)  — zero-out all but the largest connected phase component

simIceMask()             — if --simMask: run siminsar to create icemask binary,
                           then applyMask() to zero bedrock/non-ice pixels

simPhase()               — if --simPhase: run siminsar to compute simulated phase

mkdir workingDir/

interpPhase()
  ├── writeData(*.nisar.uw)           — masked unwrapped phase (binary MSB float32)
  ├── writeData(*.nisar.ion)          — ionosphere phase screen
  ├── writeData(*.nisar.cor)          — coherence magnitude
  ├── writeGeodatGeojson(geodat1/2)   — reference and secondary geodat files
  ├── runInterp()                     — intfloat hole-fill on *.nisar.uw → *.nisar.uw.interp
  ├── writeMultiBandVrt(*.uw.interp.vrt)           — VRT for hole-filled phase
  ├── writeMultiBandVrt(*.ion.filt.rangeOffset.vrt) — range correction from cleaned ionosphere
  ├── writeMultiBandVrt(*.ion.unfilt.rangeOffset.vrt) — range correction from raw ionosphere
  └── writePairInfo()                 — *.pairinfo text file
  └── writeData(*.nisar.ion.filt)     — cleaned ionosphere (if cleanIonosphere succeeded)
```

---

## Freestanding programs called

### `siminsar`

SAR product simulator. Used in two optional steps:

**Ice mask simulation** (`simIceMask`):
```
siminsar -mask <dem> <icemask> <geodat> <outputDir>/icemask
```
Produces `icemask`, a binary flat file marking non-ice pixels.

**Phase simulation** (`simPhase`):
```
siminsar -velocity -dT <dT> <dem> <velMap> <outputDir>/<geodat> <outputDir>/phaseSim
```
Produces `phaseSim.*` files containing the simulated interferometric phase.

### `intfloat`

Hole-filling interpolator. Called once on the masked unwrapped phase file:

```
intfloat -wdist -nr <nr> -na <na> -thresh <interpThresh>
    -islandThresh <islandThresh> <outputDir>/<phaseFile>
    > <outputDir>/<phaseFile>.interp
```

---

## Output files

All output is written to `outputDir/`. Filenames use the pattern
`{refOrbit}_{frame}.{secOrbit}_{frame}.{NLR}x{NLA}` for orbit/frame metadata.

| File | Description |
|------|-------------|
| `{pair}.nisar.uw` | Masked unwrapped phase (binary MSB float32, radians) |
| `{pair}.nisar.uw.interp` | Hole-filled unwrapped phase (binary MSB float32) |
| `{pair}.nisar.uw.interp.vrt` | VRT wrapping the hole-filled phase; band description: `unwrappedPhase` |
| `{pair}.nisar.cor` | Coherence magnitude (binary MSB float32) |
| `{pair}.nisar.ion` | Ionosphere phase screen (binary MSB float32, radians) |
| `{pair}.nisar.ion.filt` | Cleaned ionosphere phase screen (written only if `cleanIonosphere` succeeded) |
| `{pair}.ion.filt.rangeOffset.vrt` | VRT for range offset correction derived from cleaned ionosphere (metres → pixels) |
| `{pair}.ion.unfilt.rangeOffset.vrt` | VRT for range offset correction derived from raw ionosphere |
| `{refOrbit}.{secOrbit}.pairinfo` | Text file: `refOrbit secOrbit date1 date2 NLR NLA` |
| `geodat{NLR}x{NLA}.geojson` | Reference image geodat (GeoJSON) |
| `geodat{NLR}x{NLA}.secondary.geojson` | Secondary image geodat (GeoJSON) |
| `icemask` | Ice mask binary (only when `--simMask`) |
| `phaseSim.*` | Simulated phase products (only when `--simPhase`) |

---

## Ionosphere correction

`cleanIonosphere()` is called unconditionally on the RUNW object. If it succeeds, the `ionosphereCleaned` attribute is set and written to `*.nisar.ion.filt`. The ionosphere phase screen is also converted to a range offset correction (radians → pixels via `−λ/4π / SLCRangePixelSize`) and written as a VRT for both the cleaned and uncleaned versions.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `parseArgs()` | Parse command-line arguments; validate RUNW path; assemble `params` dict. |
| `resolveRegion(myRUNW, params)` | Determine region (`greenland`/`antarctica`) from RUNW EPSG code. |
| `simIceMask(geodat, params, outputDir)` | Run `siminsar -mask` to produce an ice mask; return `True` if mask was created. |
| `simPhase(geodat, params, dT, outputDir)` | Run `siminsar -velocity` to produce a simulated phase field. |
| `interpPhase(outputDir, myRUNW, ...)` | Write all output files: phase, coherence, ionosphere, geodats, VRTs, pairinfo; run hole-filling. |
| `runInterp(outputDir, inputFile, outputFile, nr, na, ...)` | Build and execute a single `intfloat` shell command. |
| `writePairInfo(myRUNW, outputDir)` | Write the `{ref}.{sec}.pairinfo` text file with orbit, date, and look metadata. |
