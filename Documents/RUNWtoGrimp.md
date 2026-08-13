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
| `--region NAME` | auto | Region name (`greenland` or `antarctica`). Auto-detected from `project.yaml`/RUNW EPSG if not given. |
| `--regionFile YAML` | None | Region YAML with velMap/DEM/mask paths (overrides `--region`). |
| `--verticalCorrection FILE` | None | xyDEM grid (m/yr) of submergence/emergence rate; passed to `siminsar -verticalCorrection`. |
| `--interpThresh N` | 20 | Maximum hole area (pixels) to fill during phase interpolation. |
| `--islandThresh N` | 20 | Maximum isolated-region area (pixels) to discard during interpolation. |
| `--ompThreads N` | 4 | OpenMP threads passed to `siminsar`. |
| `--phaseDerivedIonosphere` | False | **Legacy path**: run `cleanIonosphere()` and write the phase-derived ionosphere rangeOffset VRTs. Default pipeline leaves iono estimation to `estimateIonosphere`. |
| `--noPhase` | False | Skip writing/interpolating the unwrapped phase (SetupNISAR passes this in the default pipeline). |
| `--noIon` | False | Skip the ionosphere products (SetupNISAR passes this in the default pipeline). |
| `--minTol` / `--percentSpeed` / `--maxTol` | None | Variable smoothing-radius map (all three required together); applied to the interpolated phase. |
| `--maxSmoothRadius N` | 50 | Variable smoothing: sweep cap in single-look pixels (≤255). |
| `--smoothNIter N` | 3 | Variable smoothing: box-filter passes per sweep step. |
| `--noVariableSmoothing` | False | Disable the variable smoothing-radius map even if the trio is supplied. |
| `--verbose` | False | Print all subprocess output to the terminal. |

---

## Processing flow

```
parseArgs()

openHDF(RUNW)            — load unwrapped phase, coherence, ionosphere screen,
                           and image geometry

cleanIonosphere()        — only if --phaseDerivedIonosphere: apply ionosphere
                           correction (produces ionosphereCleaned attribute
                           if successful)

resolveRegion()          — determine region from project.yaml/EPSG if not set

writeGeodatGeojson()     — write the reference geodat early so siminsar
                           (simIceMask/simPhase below) can read it on a
                           fresh frame

maskPhase(largest=True)  — zero-out all but the largest connected phase component

simIceMask()             — if --simMask: run siminsar to create icemask binary,
                           then applyMask() to zero bedrock/non-ice pixels

simPhase()               — if --simPhase: run siminsar to compute simulated phase

mkdir workingDir/

interpPhase()
  ├── writeData(*.nisar.uw)           — masked unwrapped phase (skipped with --noPhase)
  ├── writeData(*.nisar.ion)          — ionosphere phase screen (--phaseDerivedIonosphere only)
  ├── writeData(*.nisar.cor)          — coherence magnitude (always)
  ├── writeGeodatGeojson(geodat1/2)   — reference and secondary geodat files (always)
  ├── runInterp()                     — intfloat hole-fill on *.nisar.uw → *.nisar.uw.interp
  │                                     (skipped with --noPhase)
  ├── writeMultiBandVrt(*.uw.interp.vrt)           — VRT for hole-filled phase (band `Phase`)
  ├── writeMultiBandVrt(*.ion.filt.rangeOffset.vrt) — cleaned iono range correction
  │                                     (--phaseDerivedIonosphere only)
  ├── writeMultiBandVrt(*.ion.unfilt.rangeOffset.vrt) — raw iono range correction
  │                                     (--phaseDerivedIonosphere only)
  └── writePairInfo()                 — *.pairinfo text file
  └── writeData(*.nisar.ion.filt)     — cleaned ionosphere (if cleanIonosphere succeeded)
```

In the default pipeline (`SetupNISAR` passes `--noPhase --noIon`), the outputs reduce to
`*.nisar.cor` (+ `.vrt`), the two geodats, the pairinfo, and any `--simMask`/`--simPhase`
products; phase and ionosphere products come later from `estimateIonosphere`.

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
| `{pair}.nisar.uw` | Masked unwrapped phase (binary MSB float32, radians) — not with `--noPhase` |
| `{pair}.nisar.uw.interp` | Hole-filled unwrapped phase (binary MSB float32) — not with `--noPhase` |
| `{pair}.nisar.uw.interp.vrt` | VRT wrapping the hole-filled phase; band description: `Phase` — not with `--noPhase` |
| `{pair}.nisar.cor` (+ `.vrt`) | Coherence magnitude (binary MSB float32) — always |
| `{pair}.nisar.ion` | Ionosphere phase screen (radians) — `--phaseDerivedIonosphere` only |
| `{pair}.nisar.ion.filt` | Cleaned ionosphere phase screen (written only if `cleanIonosphere` succeeded) |
| `{pair}.nisar.ion.filt.rangeOffset.vrt` | Range offset correction from cleaned ionosphere (radians → pixels) — `--phaseDerivedIonosphere` only |
| `{pair}.nisar.ion.unfilt.rangeOffset.vrt` | Range offset correction from raw ionosphere — `--phaseDerivedIonosphere` only |
| `{refOrbit}.{secOrbit}.pairinfo` | Text file: `refOrbit secOrbit date1 date2 NLR NLA` |
| `geodat{NLR}x{NLA}.geojson` | Reference image geodat (GeoJSON) — always |
| `geodat{NLR}x{NLA}.secondary.geojson` | Secondary image geodat (GeoJSON) — always |
| `icemask` | Ice mask binary (only when `--simMask`) |
| `phaseSim.*` | Simulated phase products (only when `--simPhase`) |

---

## Ionosphere correction

`cleanIonosphere()` is called only under `--phaseDerivedIonosphere` (legacy path; the default
pipeline estimates the ionosphere downstream with `estimateIonosphere`). If it succeeds, the
`ionosphereCleaned` attribute is set and written to `*.nisar.ion.filt`. The ionosphere phase
screen is also converted to a range offset correction (radians → pixels via
`−λ/4π / SLCRangePixelSize`) and written as a VRT for both the cleaned and uncleaned versions.

**Sign convention:** the negative scale here is not a discrepancy with
`estimateIonosphere`'s `+λ/(4π·slp)`. This path converts the ionosphere phase screen itself,
and the applied quantity is the *correction* = −ionosphere; `estimateIonosphere` applies its
positive scale to an iono estimate that is already negative (for positive ΔTEC). Both paths
therefore produce a correction that downstream consumers **ADD** to range — see the
"Background and equations" section of [estimateIonosphere.md](estimateIonosphere.md).

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
