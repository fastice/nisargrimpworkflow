# estimateIonosphere.py

Estimates and removes the dispersive ionospheric phase from a NISAR RUNW interferogram by combining the unwrapped phase with independent range offsets. Because the ionosphere shifts the carrier phase and the group delay (range) in opposite directions and by equal amounts, their combination isolates the ionospheric signal. An iterative ambiguity-correction step (using an independent velocity simulation) resolves any integer-cycle (2π) unwrapping errors that survived the initial phase processing before the final ionosphere estimate is formed.

---

## Background and equations

### The dispersive ionosphere

For a transionospheric radar, the ionosphere is dispersive: it changes the phase velocity and group velocity of the signal in opposite senses. For a round-trip SAR acquisition the key results are:

**Phase contribution** (carrier phase advances through the ionosphere, appearing as a shortened path):

$$\Delta\phi_\mathrm{ion} = -\frac{4\pi}{\lambda} \cdot \frac{K \cdot \Delta\mathrm{TEC}}{f^2}$$

**Range offset contribution** (group delay slows the modulation envelope, making the range appear longer):

$$\Delta R_\mathrm{ion} = +\frac{K \cdot \Delta\mathrm{TEC}}{f^2}$$

where K ≈ 40.28 m³ s⁻² is the ionospheric constant, ΔTEC is the differential total electron content between the two acquisition times (electrons m⁻²), f is the radar frequency (Hz), and λ = c/f is the wavelength. The crucial point is that these two effects are equal in magnitude but opposite in sign:

$$\Delta\phi_\mathrm{ion} = -\frac{4\pi}{\lambda} \cdot \Delta R_\mathrm{ion}$$

### Combining phase and offset

The observed interferometric phase is:

$$\phi_\mathrm{obs} = \phi_\mathrm{signal} + \Delta\phi_\mathrm{ion}$$

where φ_signal is the non-dispersive component (surface deformation, topography, troposphere, orbit). The residual range offset after removing the geometry-only simulation (which accounts for topography and orbit) is approximately the ionospheric group-delay component:

$$\delta R \approx \Delta R_\mathrm{ion}$$

Converting this to an equivalent phase using the dispersive relation:

$$\phi_\mathrm{offset} \equiv -\frac{4\pi}{\lambda} \cdot \delta R = \Delta\phi_\mathrm{ion}$$

Note the sign: because phase and group delay have opposite signs, multiplying the (positive) range offset by −4π/λ gives the (negative, for positive ΔTEC) ionospheric phase contribution directly. Substituting into the expression for the observed phase:

$$\phi_\mathrm{obs} + \phi_\mathrm{offset} = (\phi_\mathrm{signal} + \Delta\phi_\mathrm{ion}) + \Delta\phi_\mathrm{ion}
                                           = \phi_\mathrm{signal} + 2\,\Delta\phi_\mathrm{ion}$$

A heavy spatial Gaussian smooth (σ_az ≈ 100 px, σ_rg ≈ 30 px) eliminates the spatially variable signal component, leaving only the slowly-varying ionosphere:

$$\Delta\phi_\mathrm{ion} \approx \mathrm{smooth}\!\left(\frac{\phi_\mathrm{obs} + \phi_\mathrm{offset}}{2}\right)$$

The corrected phase is then:

$$\phi_\mathrm{corrected} = \phi_\mathrm{obs} - \Delta\phi_\mathrm{ion}$$

### Geometry removal

Before combining with the interferometric phase, the geometric contribution to the range offsets (from topography, orbital geometry, and look angle) is subtracted. `simoffsets` produces two simulated offset fields — one geometry-only (`offsets.geom`) and one including a velocity component (`offsets.velocity`). The residual

$$\delta R = (\delta R_\mathrm{meas} - \delta R_\mathrm{geom}) \times \Delta_\mathrm{slp}$$

(where Δ_slp is the SLC range pixel spacing in metres) isolates the ionospheric and noise components.

### Ambiguity correction

Unwrapping errors in the RUNW phase appear as connected components (CCs) whose mean phase is offset from neighbouring components by integer multiples of 2π. After a first-pass ionosphere estimate, the corrected phase is compared with an independent velocity simulation (velSim), which is produced by `siminsar` from a reference velocity map and DEM and is never involved in the ionosphere estimation. The per-CC mean of

$$\Delta\phi_k = \overline{(\phi_\mathrm{corrected} - \phi_\mathrm{velSim})}_k - c$$

(where c is the bias constant from the largest/reference CC) should be near zero for a correctly unwrapped component and near n × 2π for a component with an n-cycle error. The integer n is estimated as `round(Δφ_k / 2π)` and the correction `n × 2π` is subtracted from the original unwrapped phase before the ionosphere estimation is repeated.

Because the ionosphere estimation absorbs approximately half of each unwrapping error (the smooth iono picks up the low-frequency part of the jump), the apparent residual after correction is roughly n/2 rather than n. The loop therefore converges in at most ⌈log₂(N_max)⌉ passes; in practice one or two passes suffice.

---

## Usage

```
estimateIonosphere.py RUNW offsets.vrt output.vrt [options]
```

| Argument | Description |
|----------|-------------|
| `RUNW` | NISAR RUNW HDF5 file |
| `offsets.vrt` | Range offset VRT (band 1, SLC pixels, geographic SAR coordinates) |
| `output.vrt` | Output VRT path (multi-band named bands) |

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--frame N` | None | Frame number override passed to the RUNW reader. |
| `--regionFile YAML` | greenland | Region YAML for `defaultRegionDefs`; provides DEM and velocity map paths for `siminsar`. |
| `--overWrite` | False | Force `siminsar` to regenerate `velSim` and `maskVel` even if they already exist. |
| `--velThresh M/YR` | 100.0 | Velocity threshold (m/yr) for the `maskVel` `siminsar` call; pixels above this speed are excluded from the ionosphere fit. |
| `--sigma-az PX` | 100.0 | Azimuth Gaussian σ for smoothing the raw ionosphere estimate (pixels). |
| `--sigma-rg PX` | 30.0 | Range Gaussian σ (pixels). |
| `--offset-geometry VRT` | `offsets.geom.vrt` | Offset geometry VRT containing the `RangeOffsets` band to subtract from the measured offsets before ionosphere estimation. |
| `--maskFile FILE` | `maskVel.vrt` | GeoTIFF or VRT mask; only mask=1 pixels contribute to the ionosphere fit. Created by `siminsar` if not supplied. |
| `--simOffsets VRT` | `offsets.vrt` if present | VRT of simulated range offsets (SLC pixels) for per-CC diagnostic columns in the ambiguity table. |
| `--simDir DIR` | `.` (cwd) | Directory where `siminsar` outputs (`velSim.vrt`, `maskVel.vrt`) are written; created if absent. |
| `--noInterp` | off | Skip `intfloat` hole-filling of `correctedUnwrappedPhase`. |
| `--interpThresh N` | 20 | Maximum hole area (pixels) that `intfloat` will fill. |
| `--islandThresh N` | 20 | Maximum isolated-island area (pixels) that `intfloat` removes after filling. |
| `--phaseThresh RAD` | 14π | Mask `correctedUnwrappedPhase` where \|correctedPhase − simPhase\| ≥ RAD radians; screens regions of likely incorrect unwrapping. Masked pixels are written as −2×10⁹. |
| `--outputAll` | off | Write all 5 intermediate bands to a multi-band VRT (`ionosphereCorrection`, `correctedUnwrappedPhase`, `ionosphereCorrectionUnfilled`, `offsetPhase`, `unwrappedPhase`). Default: write only the 3 standard outputs as separate single-band VRTs. |
| `--verbose` | off | Print progress messages to stdout. |

---

## Processing flow

```
load_runw()                — open RUNW HDF; read unwrappedPhase, connectedComponents,
                             wavelength, SLC pixel spacing, geotransform

apply_runw_mask()          — zero pixels excluded by the interferogram mask

run_vel_sim()              — siminsar: generate velSim.vrt (radians, RUNW grid)
                             skipped if velSim already exists and --overWrite not set

run_mask_vel()             — siminsar: generate maskVel.vrt (velocity threshold mask)
                             skipped if maskVel.vrt already exists

load_offset_vrt()          — read range offsets (SLC pixels)
load_geom_vrt()            — read geometry RangeOffsets band (SLC pixels)
offset_m = (meas − geom) × SLC_pixel_size
resample_vrt_to_runw()     — bilinear resample offset_m to RUNW multilooked grid
offset_phase = offset_m × (−4π / λ)

load_mask_file()           — load maskVel.vrt (or --maskFile); resample if needed

for pass in 1..MAX_PASSES:
    estimate_and_correct()
      raw_iono = (offset_phase + phase) / 2
      valid = finite(raw_iono) & (cc ≠ 0) & mask
      fill_and_smooth_iono()
        Stage 1: NN seed     — each gap pixel ← nearest valid neighbour
        Stage 2: pyramid fill — coarse-to-fine diffusion (4 levels × 50 iters;
                                coarsest level gets 100 iters)
        Stage 3: Gaussian smooth — σ = (σ_az, σ_rg), default (100, 30) px
      iono_final -= mean(iono_final)    — zero-mean adjustment
      phase_final = phase − iono_final

    ambiguity_table()      — compare phase_final with velSim per CC;
                             compute bias c from reference (largest) CC;
                             return n_map: {cc → nearest integer n}

    if all n == 0: break   — converged

    apply_ambiguity_correction()
                           — subtract n × 2π from each CC of original phase;
                             update phase for next pass

write_geotiff() × 5        — ionosphereCorrection, correctedUnwrappedPhase,
                             ionosphereCorrectionUnfilled, offsetPhase, unwrappedPhase
write_output_vrt()         — 5-band named VRT on RUNW grid

resample_runw_to_vrt()     — regrid iono correction to offset VRT grid (metres)
write_geotiff() + write_output_vrt()
                           — single-band offset-grid VRT (ionosphereCorrection.offset.vrt)
```

---

## Output files

All outputs are written relative to the stem of the `output` argument (e.g. `output = foo/bar.vrt` → stem = `foo/bar`).

### RUNW-grid outputs (SAR time/range coordinates)

| File | Band | Units | Description |
|------|------|-------|-------------|
| `<stem>.vrt` | 1 `ionosphereCorrection` | rad | Smoothed ionospheric phase estimate Δφ_ion (to be subtracted from observed phase). |
| | 2 `correctedUnwrappedPhase` | rad | Ionosphere- and ambiguity-corrected unwrapped phase. |
| | 3 `ionosphereCorrectionUnfilled` | rad | Raw (un-gap-filled) per-pixel ionosphere estimate; NaN where cc=0 or masked. |
| | 4 `offsetPhase` | rad | Range-offset-derived phase: −(4π/λ) × δR; equal to Δφ_ion before smoothing. |
| | 5 `unwrappedPhase` | rad | Original unwrapped phase after ambiguity corrections (before iono removal). |

### Offset-grid output (geographic coordinates of the offset VRT)

| File | Band | Units | Description |
|------|------|-------|-------------|
| `<stem>.offset.vrt` | 1 `ionosphereCorrection` | m | Ionosphere correction regridded to the offset VRT coordinate system, in metres (for direct comparison with range offsets). |

---

## Gap-filling and smoothing

The gap-fill + smooth pipeline (`fill_and_smooth_iono`) is applied to the raw
per-pixel ionosphere estimate before the final Gaussian smooth. Valid pixels are
those where cc ≠ 0 (connected component is not masked), the raw estimate is
finite, and the velocity mask allows the pixel. The pipeline runs in three stages.

### Stage 1 — Nearest-neighbour seed (`_nn_fill`)

Every invalid (gap) pixel is initialised to the value of its nearest valid
neighbour, found via `scipy.ndimage.distance_transform_edt`. This produces a
Voronoi-tessellated seed image with no NaN values.

The seed step is not the final fill: it is used solely to prevent NaN
propagation into the coarser pyramid levels. `scipy.ndimage.zoom` uses bilinear
interpolation (`order=1`), which would convert a single NaN pixel into a
spreading NaN patch at the coarser level. By seeding all gaps first, the
downsampling step sees only finite values.

### Stage 2 — Coarse-to-fine pyramid fill (`_pyramid_fill`)

The image is processed through a 4-level dyadic pyramid:

```
Level 0 — full resolution          (e.g. 1000 × 2000 px)
Level 1 — ½ resolution             (e.g.  500 × 1000 px)
Level 2 — ¼ resolution             (e.g.  250 ×  500 px)
Level 3 — ⅛ resolution (coarsest)  (e.g.  125 ×  250 px)
```

Each level is built with `ndimage.zoom(..., order=1)` (bilinear) for the data
array and `ndimage.zoom(..., order=0) > 0.5` (nearest-neighbour majority vote)
for the valid-pixel mask. The NN seed from Stage 1 ensures no NaN enters
the zoom.

**Coarsest level (Level 3) — double-weight fill:**
1. Re-apply NN fill on the coarse valid mask to guarantee all coarse pixels
   are initialised.
2. Run 100 convolution fill iterations (`iterations_per_level × 2`).

At ⅛ resolution, 100 iterations of a 3×3 kernel propagate information across
gaps of order 100 × pixel_size_coarse ≈ 100 × 8 = 800 full-resolution pixels,
efficiently seeding even the widest data voids.

**Upsample passes (Levels 3→2, 2→1, 1→0) — 50 iterations each:**

For each level `lvl` from 2 down to 0:
1. **Upsample** the current filled result to the dimensions of `arrs[lvl]`
   using bilinear zoom.
2. **Merge** with the original valid pixels at this level:
   ```
   merged[i,j] = arrs[lvl][i,j]   if validids[lvl][i,j]  (original data wins)
               = upsampled[i,j]    otherwise               (propagated fill)
   ```
3. **Refine** with 50 convolution fill iterations, which smooth the boundary
   between original data and upsampled fill values.

**Convolution kernel (`_FILL_KERNEL`):**

The 3×3 kernel used by `_conv_fill` is:

```
┌──────┬──────┬──────┐
│ 0.25 │ 0.50 │ 0.25 │
├──────┼──────┼──────┤
│ 0.50 │ 0.00 │ 0.50 │
├──────┼──────┼──────┤
│ 0.25 │ 0.50 │ 0.25 │
└──────┴──────┴──────┘    (raw weights, normalised by dividing by sum = 2.5)
```

After normalisation, face-adjacent neighbours have weight **0.20** and
corner-diagonal neighbours have weight **0.10**; the centre pixel has weight
**0**. This is a weighted 8-neighbour diffusion kernel — qualitatively similar
to a single step of explicit heat-equation diffusion on a 2-D grid.

At each iteration, only gap pixels are updated; valid pixels always retain
their original values, so the fill never modifies the measured data.

**Why coarse-to-fine?**

A purely fine-resolution diffusion fill with a 3×3 kernel propagates
information at roughly 1 pixel per iteration. Filling a gap of width *W*
pixels would require O(W) iterations — many thousands for ocean-size gaps.
The pyramid approach instead:

1. Represents the gap at ⅛ resolution, where its width is only *W*/8 pixels.
2. Runs 100 fine (coarse-pixel) iterations there — effectively propagating
   across 800 full-resolution pixels of gap.
3. Upsamples the coarse solution as an initial seed for progressively finer
   levels, each needing only 50 iterations to blend the seam between
   upsampled fill and real data.

Total computation is O(N log N) rather than O(N·W) and fills correctly across
the entire image regardless of gap size.

### Stage 3 — Final Gaussian smooth

`scipy.ndimage.gaussian_filter` with σ = (σ_az, σ_rg), nominally **(100, 30)**
pixels, is applied to the fully-filled image. The heavy azimuthal smoothing
(100 px ≈ tens of km at typical NISAR multilook spacing) removes the spatially
variable surface-displacement signal from the raw ionosphere estimate, leaving
only the slowly-varying ionospheric component. The asymmetric kernel reflects
the fact that ionospheric structures are typically elongated in azimuth
(aligned with the orbital track) and shorter-scale in range.

---

## Ambiguity diagnostic table

`ambiguity_table` prints a table with one row per connected component:

| Column | Description |
|--------|-------------|
| `cc` | Connected-component label |
| `npix` | Total pixels in CC |
| `nfit` | Pixels with finite velSim overlap |
| `mean(rad)` | Mean of (corrected_phase − velSim) relative to reference CC bias |
| `median(rad)` | Median of the same residual |
| `sigma(rad)` | Standard deviation of the residual |
| `stderr(rad)` | Standard error of the mean (σ / √n); flagged `!` if > 0.1 cycles |
| `n_mean` | Nearest integer to mean / 2π — number of 2π cycles to remove |
| `n_median` | Nearest integer to median / 2π |
| `total_n` | Cumulative n across all passes (shown from pass 2 onward) |
| `off_mean`, `off_sig` | Mean and σ of iono-corrected offset residual (if `--simOffsets` given) |
| `note` | `*REF*` for reference CC; warnings for uncertain or inconsistent n |

The reference CC is the component with the most pixels. Its mean residual defines the global bias constant c; all other CCs are compared relative to it.

---

## Pipeline integration (`processFrameIonosphere`)

`estimateIonosphere` is normally invoked by `processFrameIonosphere` rather
than directly. That wrapper is called from `SetupNISAR.main()` immediately
after `processFrameROFF` completes for each frame, when
`--phaseDerivedIonosphere` is **not** selected.

```python
from nisargrimpworkflow.processFrameIonosphere import processFrameIonosphere
processFrameIonosphere(frame, myArgs, simDir='simPhase')
```

It resolves the required paths automatically:

| Input | Where it looks |
|-------|----------------|
| RUNW HDF5 | `{outputDir}/{orbit1}_{frame}/H5/NISAR*RUNW*.h5` (or same dir without `H5/`) |
| Range offset VRT | `{frameDir}/range.offsets.vrt` (or `H5/range.offsets.vrt`) |
| Geometry offset VRT | `{frameDir}/offsetSims/offsets.geom.vrt` (or `H5/...`) |

The output VRT is written to
`{frameDir}/{orbit1}_{frame}.{orbit2}_{frame}.{nLooksR}x{nLooksA}.nisar.vrt`.

Simulation outputs (`velSim`, `maskVel`) are written to `simDir` inside
`frameDir` (default `simPhase/`). The frame is skipped if the output VRT
already exists and neither `--overWrite` nor `--overWritePhase` is set.

---

## Dependencies

- `nisarhdf` — RUNW reader and geotransform utilities
- `sarfunc` — `defaultRegionDefs` for DEM and velocity map paths
- `ambiguityAnalysis` — `load_vel_sim`, `ambiguity_table`, `apply_ambiguity_correction` (co-located module)
- `processFrameIonosphere` — pipeline integration wrapper called by `SetupNISAR`
- `osgeo.gdal` — raster I/O and VRT creation
- `scipy.interpolate.RegularGridInterpolator` — bilinear VRT↔RUNW grid resampling
- `scipy.ndimage` — `distance_transform_edt` (NN seed), `convolve` (diffusion fill), `zoom` (pyramid levels)
- `scipy.ndimage.gaussian_filter` — final ionosphere smoothing
- `siminsar` — C binary on PATH; called for `velSim` and `maskVel` generation
