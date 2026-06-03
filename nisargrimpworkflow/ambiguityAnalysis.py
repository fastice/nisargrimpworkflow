"""
ambiguityAnalysis.py

Diagnose and correct integer-cycle (2π) unwrapping errors between connected
components by comparing the ionosphere-corrected phase with an independent
velocity simulation (velSim).

Key idea
--------
After ionosphere correction the corrected phase and velSim should agree up to
a constant offset (orbital/tropospheric bias) plus integer multiples of 2π.
Because velSim was never used in the ionosphere estimation it is immune to
the circularity that caused the previous offset-based approach to fail.

Workflow (called from estimateIonosphere.py)
--------------------------------------------
1. Load velSim.vrt (radians, same RUNW grid as corrected phase).
2. First-pass: ambiguity_table(phase_out1, vel_sim, cc)
   → prints diagnostic table, returns n_map {cc_label: n_mean}.
3. apply_ambiguity_correction(phase_original, cc, n_map)
   → subtracts n × 2π from each component of the *original* unwrapped phase.
4. Second-pass iono estimation on ambiguity-corrected original phase.
5. Second-pass: ambiguity_table(phase_out2, vel_sim, cc)
   → verification; all n values should be 0.
"""

import math
import os
import sys

import numpy as np

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    sys.exit("ambiguityAnalysis: GDAL (osgeo) not available")


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_vel_sim(path: str = 'velSim.vrt') -> np.ndarray:
    """Load a velSim GeoTIFF or VRT and return a float32 array (radians).

    Parameters
    ----------
    path : str
        Path to velSim.vrt (or velSim.tif / velSim binary that GDAL can open).

    Returns
    -------
    np.ndarray, shape (naz, nrg), dtype float32
        Simulated interferometric phase in radians.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    RuntimeError
        If GDAL cannot open the file or band count is wrong.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"velSim not found at '{path}'. "
            "Run estimateIonosphere with --overWrite to regenerate it, "
            "or check that siminsar completed successfully."
        )
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"GDAL could not open '{path}'")
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        data[data == np.float32(nodata)] = np.nan
    ds = None
    return data


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def ambiguity_table(corrected_phase: np.ndarray,
                    vel_sim: np.ndarray,
                    cc: np.ndarray,
                    n_accumulated: dict = None,
                    offset_resid: np.ndarray = None,
                    offset_units: str = 'rad') -> dict:
    """Identify integer-cycle unwrapping errors per connected component.

    Parameters
    ----------
    corrected_phase : np.ndarray, shape (naz, nrg), float32
        Iono-corrected unwrapped phase (``phase_out`` from
        ``estimate_and_correct``), radians.  After iono correction a component
        with a true N-cycle unwrapping error retains approximately N×π of
        apparent error (the iono estimation absorbs the other half).
        Dividing the per-component mean by 2π gives the apparent n, which
        is ~N/2.  The loop in ``main()`` iterates until all n→0.
    vel_sim : np.ndarray, shape (naz, nrg), float32
        Simulated interferometric phase from velSim.vrt, radians.
    cc : np.ndarray, shape (naz, nrg), integer
        Connected component labels (0 = invalid/masked).
    n_accumulated : dict[int, int], optional
        Running total of n corrections already applied in previous calls
        (cc_label → cumulative n).  When supplied, a ``total_n`` column is
        added showing n_accumulated[k] + n_this_call[k].

    Returns
    -------
    dict[int, int]
        Mapping of CC label → n (integer number of 2π cycles to remove from
        the original unwrapped phase).  The reference CC always has n=0.
        Returns {} on early exit.
    """
    TWO_PI = 2.0 * math.pi
    # Flag when stderr exceeds this fraction of 2π (one step in the n grid).
    # 0.1 × 2π ≈ 0.63 rad — uncertain n estimate.
    STDERR_CYCLE_THRESH = 0.1

    # ------------------------------------------------------------------
    # 1. Shape check
    # ------------------------------------------------------------------
    if corrected_phase.shape != vel_sim.shape:
        raise ValueError(
            f"Shape mismatch: corrected_phase {corrected_phase.shape} "
            f"vs vel_sim {vel_sim.shape}. "
            "velSim must be on the same RUNW multilooked grid."
        )
    if corrected_phase.shape != cc.shape:
        raise ValueError(
            f"Shape mismatch: corrected_phase {corrected_phase.shape} "
            f"vs cc {cc.shape}."
        )

    # ------------------------------------------------------------------
    # 2. Find largest connected component as reference
    # ------------------------------------------------------------------
    labels, counts = np.unique(cc[cc != 0], return_counts=True)
    if len(labels) == 0:
        print("  [ambiguity_table] No valid connected components (all cc=0). Skipping.")
        return {}
    ref_cc = int(labels[np.argmax(counts)])

    # ------------------------------------------------------------------
    # 3. Compute the global bias constant c from the reference CC
    #    c = mean(corrected_phase − vel_sim) over ref CC, finite pixels only
    # ------------------------------------------------------------------
    ref_mask = (
        (cc == ref_cc)
        & np.isfinite(corrected_phase)
        & np.isfinite(vel_sim)
    )
    if not ref_mask.any():
        print(
            f"  [ambiguity_table] Reference CC={ref_cc} has no finite overlap "
            "with velSim. Cannot compute bias. Skipping."
        )
        return {}

    diff = (corrected_phase - vel_sim).astype(np.float64)
    c = float(np.mean(diff[ref_mask]))
    diff_adj = diff - c   # reference CC now has mean ≈ 0; others ≈ N×π

    # ------------------------------------------------------------------
    # 4. Per-component statistics
    # ------------------------------------------------------------------
    results = []
    for k, npix in zip(labels, counts):
        mask_k = (
            (cc == k)
            & np.isfinite(corrected_phase)
            & np.isfinite(vel_sim)
        )
        nfit = int(mask_k.sum())
        if nfit == 0:
            results.append(dict(
                cc=int(k), npix=int(npix), nfit=0,
                mean=np.nan, median=np.nan, sigma=np.nan,
                n_mean=0, n_median=0, is_ref=(int(k) == ref_cc),
            ))
            continue

        resid = diff_adj[mask_k]
        mean_k   = float(np.mean(resid))
        median_k = float(np.median(resid))
        sigma_k  = float(np.std(resid))
        stderr_k = sigma_k / math.sqrt(nfit)   # uncertainty of the mean
        n_mean_k   = int(round(mean_k   / TWO_PI))
        n_median_k = int(round(median_k / TWO_PI))

        # Optional: offset residual stats (iono-corrected offset vs simulated)
        if offset_resid is not None:
            off_mask_k = mask_k & np.isfinite(offset_resid)
            if off_mask_k.any():
                off_mean_k  = float(np.mean(offset_resid[off_mask_k]))
                off_sigma_k = float(np.std(offset_resid[off_mask_k]))
            else:
                off_mean_k = off_sigma_k = np.nan
        else:
            off_mean_k = off_sigma_k = np.nan

        results.append(dict(
            cc=int(k), npix=int(npix), nfit=nfit,
            mean=mean_k, median=median_k, sigma=sigma_k, stderr=stderr_k,
            n_mean=n_mean_k, n_median=n_median_k,
            is_ref=(int(k) == ref_cc),
            off_mean=off_mean_k, off_sigma=off_sigma_k,
        ))

    # Sort by decreasing pixel count (largest = ref first, then others)
    results.sort(key=lambda r: -r['npix'])

    # ------------------------------------------------------------------
    # 5. Print table
    # ------------------------------------------------------------------
    c_m = c  # radians; could convert to cycles: c / TWO_PI
    print(f"\nReference CC: {ref_cc}  ({int(ref_mask.sum()):,} finite pixels in velSim overlap)")
    print(f"Global bias c = {c:+.4f} rad  ({c / TWO_PI:+.4f} cycles)  "
          "[mean(corrected_phase − velSim) over reference CC]\n")

    show_total  = n_accumulated is not None
    show_off    = offset_resid is not None
    total_hdr   = f"  {'total_n':>7}" if show_total else ""
    off_hdr     = (f"  {f'off_mean({offset_units})':>13}  {f'off_sig({offset_units})':>12}"
                   if show_off else "")
    hdr = (f"{'cc':>6}  {'npix':>8}  {'nfit':>8}  "
           f"{'mean(rad)':>10}  {'median(rad)':>11}  "
           f"{'sigma(rad)':>10}  {'stderr(rad)':>11}  "
           f"{'n_mean':>7}  {'n_median':>8}{total_hdr}{off_hdr}  {'note'}")
    print(hdr)
    print('-' * len(hdr))

    for r in results:
        total_n   = (n_accumulated or {}).get(r['cc'], 0) + r.get('n_mean', 0)
        total_col = f"  {total_n:+7d}" if show_total else ""

        if r['nfit'] == 0:
            note    = 'NO VELsim OVERLAP'
            off_col = f"  {'--':>13}  {'--':>12}" if show_off else ""
            print(f"{r['cc']:6d}  {r['npix']:8,}  {r['nfit']:8,}  "
                  f"{'--':>10}  {'--':>11}  {'--':>10}  {'--':>11}  "
                  f"{'--':>7}  {'--':>8}"
                  + (f"  {'--':>7}" if show_total else "")
                  + off_col
                  + f"  {note}")
            continue

        stderr_cycles = r['stderr'] / TWO_PI
        stderr_flag   = '!' if stderr_cycles > STDERR_CYCLE_THRESH else ' '
        note_parts = []
        if r['is_ref']:
            note_parts.append('*REF*')
        if stderr_cycles > STDERR_CYCLE_THRESH:
            note_parts.append(f'stderr={stderr_cycles:.3f} cycles: n uncertain')
        if r['n_mean'] != r['n_median']:
            note_parts.append('n_mean≠n_median')
        note = '  '.join(note_parts)

        if show_off and np.isfinite(r.get('off_mean', np.nan)):
            off_col = f"  {r['off_mean']:+13.4f}  {r['off_sigma']:12.4f}"
        elif show_off:
            off_col = f"  {'--':>13}  {'--':>12}"
        else:
            off_col = ""

        print(f"{r['cc']:6d}  {r['npix']:8,}  {r['nfit']:8,}  "
              f"{r['mean']:+10.4f}  {r['median']:+11.4f}  "
              f"{r['sigma']:10.4f}   {r['stderr']:10.4f}{stderr_flag}  "
              f"{r['n_mean']:+7d}  {r['n_median']:+8d}"
              + total_col
              + off_col
              + f"  {note}")

    print()
    n_flagged = sum(
        1 for r in results
        if r['nfit'] > 0 and (r['stderr'] / TWO_PI) > STDERR_CYCLE_THRESH
    )
    n_nonzero = sum(1 for r in results if r['nfit'] > 0 and not r['is_ref']
                    and (r['n_mean'] != 0 or r['n_median'] != 0))
    print(f"Summary: {len(results)} components  |  "
          f"{n_nonzero} with non-zero n  |  "
          f"{n_flagged} with stderr > {STDERR_CYCLE_THRESH} cycles (flagged '!')")
    print()

    # Build and return the n_map: cc_label → n_mean
    n_map = {r['cc']: r['n_mean'] for r in results if r['nfit'] > 0}
    return n_map


# ---------------------------------------------------------------------------
# Correction
# ---------------------------------------------------------------------------

def apply_ambiguity_correction(phase: np.ndarray,
                                cc: np.ndarray,
                                n_map: dict) -> np.ndarray:
    """Apply integer-cycle corrections to the original unwrapped phase.

    Subtracts ``n × 2π`` from every finite pixel in each connected component,
    where ``n`` comes from the n_map returned by ``ambiguity_table``.

    Parameters
    ----------
    phase : np.ndarray, shape (naz, nrg), float32
        Original (pre-iono-correction) unwrapped phase, radians.
        Not modified in place — a copy is returned.
    cc : np.ndarray, shape (naz, nrg), integer
        Connected component labels (0 = invalid/masked).
    n_map : dict[int, int]
        CC label → integer number of 2π cycles to remove (from ambiguity_table).

    Returns
    -------
    np.ndarray, float32
        Phase with ambiguity corrections applied.
    """
    TWO_PI = 2.0 * math.pi
    corrected = phase.copy()
    n_applied = {k: n for k, n in n_map.items() if n != 0}
    if not n_applied:
        print("  No ambiguity corrections needed (all n = 0).")
        return corrected
    for k, n in sorted(n_applied.items()):
        mask = (cc == k) & np.isfinite(phase)
        corrected[mask] -= np.float32(n * TWO_PI)
        print(f"  CC {k:4d}: n = {n:+d}  →  removed {n * TWO_PI:+.4f} rad "
              f"({int(mask.sum()):,} pixels)")
    return corrected
