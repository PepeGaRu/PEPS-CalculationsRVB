#!/usr/bin/env python3
"""Compute physical quantum correlation lengths on exact spin-1 cylinders."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.sparse.linalg import eigs
from rvb.octahedral_spin1_transfer import Spin1OctahedralCylinder
from rvb.octahedral_spin1_width import Spin1OctahedralCylinder23


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    if isinstance(value, np.ndarray):
        return serial(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.generic):
        return value.item()
    return value


def physical_modes(analysis, threshold=1e-18):
    modes = [m for m in analysis["modes"] if abs(m["weight"]) > threshold and m["xi"] is not None]
    return sorted(modes, key=lambda m: -abs(m["ratio"]))


def fits(analysis, windows):
    r = analysis["distance"]
    c = np.real(analysis["correlation"])
    records = []
    for first, last in windows:
        mask = (r >= first) & (r <= last) & (c != 0)
        slope, intercept = np.polyfit(r[mask], np.log(abs(c[mask])), 1)
        records.append({"window": [first, last], "xi": float(-1/slope),
                        "log_amplitude": float(intercept),
                        "max_log_residual": float(np.max(abs(np.log(abs(c[mask]))-(slope*r[mask]+intercept))))})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results/octahedral-spin1-cylinder")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    tube = Spin1OctahedralCylinder()
    longitudinal = tube.analyze("longitudinal_dimer_projector", distances=range(3, 65))
    transverse = tube.analyze("dimer_projector", distances=range(1, 65))
    spin = tube.analyze("spin_z", distances=range(1, 17))
    leading = physical_modes(longitudinal)[0]
    fits_long = fits(longitudinal, ((6, 16), (12, 28), (24, 48)))
    fits_trans = fits(transverse, ((8, 20), (12, 24), (20, 40)))
    sector_spectra = tube.sector_spectra()
    sector_variation = []
    for sector in range(0, 64, 9):
        analysis = tube.analyze("longitudinal_dimer_projector", distances=range(3, 8), sector=sector)
        sector_variation.append({"sector": sector, "lambda0": analysis["lambda0"],
                                 "mean": analysis["mean"], "xi": physical_modes(analysis)[0]["xi"]})
    other_leading = max(abs(item["eigenvalues"][0]) for item in sector_spectra if item["sector"] != 0)
    sector_selection_scale = float(-1/np.log(other_leading/longitudinal["lambda0"]))
    print(f"2x2: longitudinal dimer xi={leading['xi']:.12f}, mean={longitudinal['mean'].real:.12f}", flush=True)

    wide = Spin1OctahedralCylinder23()
    wide_analysis = wide.analyze(distances=range(1, 49))
    wide_other = wide.analyze(distances=range(1, 13), parity=1)
    wide_leading = physical_modes(wide_analysis)[0]
    wide_sectors = wide.diagonal_spectra()
    cross, _, _ = wide.block(parities=(0, 1))
    cross_radius = float(max(abs(eigs(cross, k=3, which="LM", tol=1e-12,
                                     v0=np.ones(len(cross)), return_eigenvectors=False))))
    torus_wide, torus_sectors = wide.finite_length_two_norm()
    torus_rotated = float((tube.transfer @ tube.transfer).multiply(tube.transfer.T).sum())
    rotation_error = abs(torus_wide-torus_rotated)/torus_rotated
    assert rotation_error < 1e-11
    assert abs(wide_analysis["mean"]-wide_other["mean"]) < 1e-11
    print(f"2x3: transverse dimer xi={wide_leading['xi']:.12f}, axis-rotation norm error={rotation_error:.2e}", flush=True)
    elapsed = perf_counter()-start
    result = {
        "state": "Equal coefficients for normalized spin-1 singlet dimer/trimer coverings; axis-ordered epsilon trimers",
        "geometry": "Infinite x, periodic transverse dimensions 2x2 and 2x3",
        "units": "One unit is the spacing between adjacent octahedron centres along x",
        "contraction": "Exact finite-width transfer, no tensor truncation; floating-point eigenanalysis",
        "dominant_mode_weight_threshold": 1e-18,
        "tube_2x2": {"longitudinal_projector": longitudinal, "transverse_projector": transverse,
                     "spin_z": spin, "longitudinal_window_fits": fits_long,
                     "transverse_window_fits": fits_trans,
                     "sector_spectra": sector_spectra, "diagonal_sector_variation": sector_variation,
                     "norm_transfer_nonzeros": tube.transfer.nnz,
                     "other_sector_selection_scale_not_local_xi": sector_selection_scale},
        "tube_2x3": {"transverse_projector": wide_analysis,
                     "other_leading_sector": wide_other, "diagonal_sector_spectra": wide_sectors,
                     "leading_cross_sector_radius": cross_radius},
        "independent_axis_rotation_check": {"spins": 36, "shape": [3, 2, 2],
            "trace_transfer_2x2_cubed": torus_rotated, "trace_transfer_2x3_squared": torus_wide,
            "relative_error": rotation_error, "wide_sector_representatives": torus_sectors},
        "calculation_seconds": elapsed,
        "timing_scope": "Both widths, sector scans, observables and 36-spin axis-rotation check; excludes plotting and tests",
        "bulk_3d_correlation_length": None,
        "bulk_conclusion": "Finite-width quantum lengths are resolved; width dependence is large, so no 3D extrapolation is established.",
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            "rvb/octahedral_spin1.py", "rvb/octahedral_spin1_transfer.py",
            "rvb/octahedral_spin1_width.py", "scripts/run_octahedral_spin1_cylinder.py")},
    }
    (out/"data.json").write_text(json.dumps(serial(result), indent=2, allow_nan=False)+"\n")

    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    r = longitudinal["distance"]
    c = longitudinal["correlation"].real
    for sign, color, label in ((1, "#276aa4", "Positive C(r)"), (-1, "#c36b27", "Negative C(r)")):
        mask = np.sign(c) == sign
        axes[0].semilogy(r[mask], abs(c[mask]), "o", ms=3.8, color=color, label=label)
    asymptote = abs(leading["weight"])*abs(leading["ratio"])**(r-3)
    axes[0].semilogy(r, asymptote, "k-", lw=1.1, alpha=.7, label=f"Exponential: xi = {leading['xi']:.5f}")
    axes[0].set(title="Longitudinal singlet-projector decay", xlabel="Centre separation r along x",
                ylabel="Absolute connected correlation", xlim=(2, 65))
    axes[0].legend(fontsize=8)
    axes[0].text(.06, .08, "Infinite length; 2 × 2 cross-section\nSigns alternate at successive sites",
                 transform=axes[0].transAxes, fontsize=9)
    effective = -1/np.log(abs(c[1:]/c[:-1]))
    axes[1].plot(r[:-1]+.5, effective, "o", color="#276aa4", ms=3)
    axes[1].axhline(leading["xi"], color="black", lw=1, label="Observable-coupled eigenvalue")
    axes[1].set(title="Local decay length stabilizes", xlabel="Separation r + 1/2",
                ylabel="Effective correlation length", ylim=(2.48, 2.50))
    axes[1].legend(fontsize=8, loc="upper right")
    for analysis, shape, color in ((transverse, "2 × 2", "#276aa4"), (wide_analysis, "2 × 3", "#824c9d")):
        distance = analysis["distance"]
        correlation = np.real(analysis["correlation"])
        display = distance <= 28
        mode = physical_modes(analysis)[0]
        axes[2].semilogy(distance[display], abs(correlation[display]), "o-", ms=3, lw=1,
                         color=color, label=f"{shape}: xi = {mode['xi']:.5f}")
    axes[2].set(title="Width check: transverse projector", xlabel="Centre separation r along x",
                ylabel="Absolute connected correlation")
    axes[2].legend(fontsize=9)
    axes[2].text(.05, .08, "Width dependence remains large\nNo 3D bulk length inferred", transform=axes[2].transAxes, fontsize=9)
    for ax in axes:
        ax.grid(alpha=.17)
    fig.suptitle("Spin-1 singlet coverings: measured quantum correlation lengths on infinite cylinders", fontsize=14)
    fig.savefig(out/"correlation-lengths.png", dpi=170, bbox_inches="tight")
    fig.savefig(out/"correlation-lengths.pdf", bbox_inches="tight")

    page = plt.figure(figsize=(8.27, 11.69))
    page.text(.085, .94, "A quantum correlation length", fontsize=23, weight="bold", color="#276aa4")
    page.text(.085, .89, f"Spin-1 singlet dimers and trimers: xi = {leading['xi']:.7f}", fontsize=17)
    sections = [
        (.83, "What was measured", "The connected correlation of physical singlet projectors on\n"
         "longitudinal diameter pairs. Distances count cubic-centre spacings.\n"
         "The covering sum includes nonorthogonal quantum interference."),
        (.69, "Geometry and exact finite-width method", "The cylinder is infinite along x and periodic with 2 × 2 centres\n"
         "across y,z. Its full norm transfer has 4,096 states, separated into\n"
         "64 invariant sectors. No tensor singular values are discarded."),
        (.55, "Resolved exponential envelope", f"lambda_0 = {longitudinal['lambda0']:.10f}\n"
         f"lambda_dimer = {leading['eigenvalue'].real:.10f}\n"
         f"xi = -1/log|lambda_dimer/lambda_0| = {leading['xi']:.10f}\n"
         "The negative eigenvalue gives alternating signs, not a power law."),
        (.38, "Independent checks", "Physical projector insertions reproduce the exact 24-spin check.\n"
         "Two cylinder contractions of the same 36-spin torus agree.\n"
         "Direct correlation propagation, several fit windows and an\n"
         "80-digit scalar calculation agree on the 2 × 2 decay length."),
        (.21, "What remains unresolved", f"The same transverse projector gives xi = {wide_leading['xi']:.7f} on\n"
         "a 2 × 3 cylinder, versus 2.4901024 on 2 × 2. This large width\n"
         "dependence prevents a three-dimensional bulk extrapolation.\n"
         "These are computed quantum cylinder lengths, not a bulk claim."),
    ]
    for y, title, body in sections:
        page.text(.085, y, title, fontsize=13, weight="bold")
        page.text(.085, y-.026, body, fontsize=10.4, va="top", linespacing=1.55)
    page.text(.085, .05, "Reproduce: python scripts/run_octahedral_spin1_cylinder.py\n"
              "Precision check: python scripts/check_spin1_cylinder_precision.py", fontsize=9, color="#68747f")
    pdf_path = out/"Spin1_Quantum_Correlation_Lengths.pdf"
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(page)
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(page); plt.close(fig)

    fit_rows = "\n".join(f"| {v['window'][0]}–{v['window'][1]} | {v['xi']:.12f} |" for v in fits_long)
    report = f"""# Resolved quantum correlation lengths for spin-1 singlet coverings

**The infinite 2x2 cylinder has a physical singlet-projector correlation length
xi = {leading['xi']:.10f} cubic-centre spacings.** This is a new quantum
calculation including covering overlaps, using the equal-coefficient,
normalized-singlet state defined in [the model notes](../../docs/octahedral-spin1-model.md).
The earlier classical zero-tail result is not used to infer this length.

## Measurement and geometry

Octahedron centres have integer coordinates and the physical spin-1 sites
are cubic-bond midpoints. The x direction is infinite. The transverse y,z
directions are periodic with two centres each.

The principal observable is the physical singlet projector P_x(v) on the
sites v-e_x/2 and v+e_x/2, where `P = ((S_i · S_j)^2 - I)/3`.
We measure `C(r) = <P_x(0) P_x(r e_x)> - <P_x>^2`.
The mean is {longitudinal['mean'].real:.12f}. Each transfer step translates
one cubic-centre spacing. The implementation uses an exact three-slice
insertion and stores non-overlapping insertion blocks at r >= 3; shorter
distances are not assumed zero.

![Quantum correlation lengths and width check](correlation-lengths.png)

## The measured tail

The leading observable-coupled eigenvalues are
`lambda_0 = {longitudinal['lambda0']:.12f}` and
`lambda_dimer = {leading['eigenvalue'].real:.12f}`. Consequently
`q = lambda_dimer/lambda_0 = {leading['ratio'].real:.12f}` and
`xi = -1/log|q| = {leading['xi']:.12f}`.
The negative q gives an alternating tail. With the three-slice convention,
`C(r) ~ {leading['weight'].real:.12f} * ({leading['ratio'].real:.12f})^(r-3)`.
Its magnitude decays exponentially; a negative correlator is not a failure
of the decay fit. The asymptotic term couples strongly to this observable.

Independent fits of the actual propagated envelope give:

| Distance window | Fitted xi |
|---|---:|
{fit_rows}

The centered propagation removes the disconnected contribution algebraically.
It does not obtain a tiny tail by subtracting two order-one numbers.
The spectral and direct correlation curves agree within
{longitudinal['spectral_direct_max_abs_difference']:.3g} in absolute value.
The maximum norm-eigenvector residual in the 164-dimensional selected block
is {longitudinal['max_eigen_relative_residual']:.3g} with the recorded normalization.

The transverse diameter projector P_y has the same asymptotic length on this
cylinder, but its leading weight is only about -9.53e-10. Short-distance fits
can instead see a faster transient with xi about 0.426. An independent
eight-dimensional SU(2)-scalar reduction evaluated at 80 digits confirms
the weak tail and the eigenvalue ratio. That precision audit also compares
unrounded input matrices with rationally reconstructed matrices; the length
changes by about 2e-15. See [high-precision data](high-precision.json).

For comparison, the physical single-spin correlation `<S_z(0)S_z(r)>`
on outgoing y sites has asymptotic xi = {physical_modes(spin)[0]['xi']:.10f}.
Spin and dimer lengths describe different operators and need not agree.

## Wider-cylinder check

For the same transverse P_y observable, the exact 2x3 periodic transverse
section gives **xi = {wide_leading['xi']:.10f}**. Its mean is
{wide_analysis['mean'].real:.12f}. The dominant boundary block has 976 states;
it is retained in full. The two leading parity sectors have equal dominant
eigenvalues and identical means/correlations. Their cross-sector radius is
{cross_radius:.6f}, well below lambda_0 = {wide_analysis['lambda0']:.6f}.
All sixteen diagonal sectors were compared. The mixed-sector spectral
radii are bounded by the geometric means of their diagonal-sector radii;
the cross block between the tied leaders was checked explicitly.

The leading transverse modes on 2x3 are close in decay rate, so a finite-window
fit need not equal the ultimate largest-coupled-eigenvalue length exactly.
All coupled modes and direct data are saved. Odd versus even transverse
periods also change the allowed layer patterns. These two widths are not
a converged scaling sequence. **The bulk three-dimensional correlation
length and bulk exponential-versus-algebraic classification remain open.**

## Boundary sectors and a misleading long scale

The 2x2 norm transfer has 4,096 compressed states and {tube.transfer.nnz:,}
nonzero entries. Its 64 invariant bra/ket ownership blocks have at most
164 states. The largest norm eigenvalue occurs uniquely in sector 0.
This is the infinite-x limit of the finite periodic equal-covering sum.

Other sectors have nearly equal norm growth rates. Using the unrestricted
second eigenvalue would produce a scale of approximately
{sector_selection_scale:.2f} cells. It is a sector-selection scale, not this
local dimer correlation length: the physical insertions preserve the
ownership sectors and do not couple the two leading sectors. It also means
that a short finite torus summed over sectors need not approximate the
infinite-cylinder sector mixture. In separately fixed diagonal sectors,
the longitudinal local lengths range from
{min(v['xi'] for v in sector_variation):.7f} to {max(v['xi'] for v in sector_variation):.7f}.
This small sector variation is not an error bar for increasing the width.

## Validation and reproducibility

The transfer exactly retains the rank-eight support of every unobserved
physical-edge norm metric. Observed edges are inserted before this compression,
including the components that lie outside the norm support. Both the
longitudinal and transverse projector constructions reproduce independently
computed physical quantities on the 24-spin torus.

An additional rotation check contracts the same 36-spin 3x2x2 torus in two
directions: `Tr(T_2x2^3)` and `Tr(T_2x3^2)`. Summing all sectors gives
{torus_rotated:.12f} and {torus_wide:.12f}, relative difference
{rotation_error:.3g}. This checks the wider geometry independently of fitting.
Other tests include local spin symmetry, complete covering sums, identity
insertions, a directly constructed twelve-spin state and the dimer-only limit.

```bash
/workspace/peps-environment/.venv/bin/python scripts/run_octahedral_spin1_cylinder.py
/workspace/peps-environment/.venv/bin/python scripts/check_spin1_cylinder_precision.py
```

The main calculation and sector checks took {elapsed:.2f} seconds here,
excluding plotting and tests. Exact finite width means no tensor truncation;
the standard eigenanalysis uses floating-point arithmetic. Numerical precision
within a fixed tube is distinct from the physical uncertainty of its width.
The optional precision check uses mpmath, already installed in this environment
and available as the project's `precision` extra.

- [PDF with the measured lengths and plots](Spin1_Quantum_Correlation_Lengths.pdf)
- [Raw correlations, spectra, fits and source hashes](data.json)
- [80-digit reduced-matrix check](high-precision.json)
- [2x2 transfer implementation](../../rvb/octahedral_spin1_transfer.py)
- [2x3 transfer implementation](../../rvb/octahedral_spin1_width.py)
"""
    (out/"report.md").write_text(report)
    print(f"Calculations and checks: {elapsed:.2f} seconds", flush=True)
    print(f"Wrote {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
