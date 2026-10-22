#!/usr/bin/env python3
"""Plot this project's square-octagon calculations, with source provenance.

No published correlation curves or assumed decay laws supply the plotted data.
Run the CTM, exact-cylinder, and mapping scripts first to regenerate inputs.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.diagnostics import diagnose_decay
from rvb.square_octagon import SquareOctagonOrthogonalDimers


def main():
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlepad": 12,
                         "legend.frameon": False, "savefig.facecolor": "white"})
    out = ROOT / "results/octagon-square"
    input_paths = []

    def read(relative):
        path = out / relative
        input_paths.append(path)
        return json.loads(path.read_text())

    summary = read("ctmrg/summary.json")
    selected = read("ctmrg/" + summary["selected_run"])
    comparison = read("ctmrg/" + summary["comparison_run"])
    orthogonal = read("orthogonal.json")
    floor = summary["absolute_correlation_floor"]
    windows = [(2, 5), (4, 8), (8, 14), (12, 18), (16, 20)]
    key = "triangle_ring_ring_x"
    diagnostic = summary["selected_observables"][key]["diagnostics"]
    fit_window = next(w for w in diagnostic["windows"]
                      if w["requested_window"] == [12., 18.])
    fit = fit_window["models"]
    blue, orange, green, red = "#245e9d", "#ca7925", "#268368", "#b13c53"

    def values(run, name):
        rows = run["observables"][name]["rows"]
        return (np.array([row["distance"] for row in rows]),
                np.array([row["connected_centered"] for row in rows]))

    def save(fig, name):
        fig.savefig(out / (name + ".png"), dpi=175, bbox_inches="tight")
        fig.savefig(out / (name + ".pdf"), bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.1), layout="constrained")
    ax = axes[0]
    curves = [("triangle_ring_ring_x", "Triangle: two ring spins", blue, "o"),
              ("triangle_ring_dummy_x", "Triangle: ring + dummy spin", orange, "s"),
              ("dummy_bond_x", "Singlet bond along dummy link", green, "^")]
    for name, label, color, marker in curves:
        r, c = values(selected, name)
        keep = (r <= 18) & (c > floor)
        ax.semilogy(r[keep], c[keep], marker=marker, markersize=4,
                    lw=1.4, color=color, label=label)
    rr = np.linspace(12, 18, 80)
    cc = np.exp(fit["exponential"]["log_amplitude"]
                - rr / fit["exponential"]["xi"])
    ax.semilogy(rr, cc, "--", color="#222222", lw=1.6,
                label=r"Exponential fit, $r=12\ldots18$: $\xi\approx0.82$")
    ax.set(title="Singlet RVB: rapidly decaying physical bonds",
           xlabel="Separation along x (simulation cells)",
           ylabel=r"$C_B(r)=\langle B_0B_r\rangle-\langle B_0\rangle\langle B_r\rangle$",
           xlim=(.5, 18.5), ylim=(2e-13, .06))
    ax.set_xticks([1, 4, 8, 12, 16, 18])
    ax.legend(loc="upper right", fontsize=9)
    ax.text(.04, .07, r"$B=\mathbf{S}_i\!\cdot\!\mathbf{S}_j$" + "\n"
            + r"CTMRG, $\chi=64$; characteristic scale $\xi\sim0.8$–$0.9$"
            + "\nFinite-distance estimate; cutoff sensitivity remains",
            transform=ax.transAxes, fontsize=9.5, va="bottom")
    ax.grid(alpha=.2, which="major")

    ax = axes[1]
    distance = np.array(orthogonal["distance"])
    for family, color, marker, label in zip(
            orthogonal["families"], (blue, green), ("o", "s"),
            ("Triangle bond, same orientation", "Dummy link, same orientation")):
        keep = distance <= 6
        ax.plot(distance[keep], np.array(family["connected_correlation"])[keep],
                marker=marker, color=color, lw=1.4, markersize=6,
                markerfacecolor="white", label=label)
    # A distinct orientation demonstrates the allowed nearest-cell covariance.
    model = SquareOctagonOrthogonalDimers(16, 4)
    reference = model.geometry.triangle_edges[0][0][0]
    cross_r = list(range(7))
    cross = [model.connected(reference, model.geometry.triangle_edges[r*4][2][0])
             for r in cross_r]
    ax.plot(cross_r, [float(c) for c in cross], "D:", color=orange, markersize=4,
            label="Triangle bonds, opposite orientations")
    ax.axhline(0, color="#888888", lw=.7, zorder=0)
    ax.set(title="Orthogonal dimers: exactly finite range",
           xlabel="Cell displacement along x (simulation cells)",
           ylabel=r"$\langle n_0n_r\rangle-\langle n_0\rangle\langle n_r\rangle$",
           xlim=(-.2, 6.3), ylim=(-.015, .28), xticks=range(7))
    ax.legend(loc="upper right", fontsize=9)
    ax.text(.48, .37, "Exact rational counting\nAll bond pairs in cells farther\n"
            "than nearest neighbors: zero\n\nNo nonzero decay tail to fit",
            transform=ax.transAxes, fontsize=10, va="center")
    ax.grid(alpha=.2)
    fig.suptitle("Your decorated square–octagon model: new calculations", fontsize=15)
    save(fig, "correlations")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), layout="constrained")
    for run, color, marker, label in (
            (selected, blue, "o", r"Projector cutoff $10^{-14}$"),
            (comparison, orange, "s", r"Projector cutoff $10^{-12}$")):
        r, c = values(run, key)
        diag = diagnose_decay(r, c, windows=windows, floor=floor)
        local = [point for point in diag["local"] if point["r_right"] <= 18]
        position = [(point["r_left"] + point["r_right"]) / 2 for point in local]
        axes[0].plot(position, [p["inverse_xi_eff"] for p in local],
                     marker=marker, markersize=4, color=color, lw=1.5, label=label)
        axes[1].plot(position, [p["p_eff"] for p in local],
                     marker=marker, markersize=4, color=color, lw=1.5, label=label)
    mix = fit["power_exponential"]
    r = np.linspace(8, 17, 120)
    # Exact finite-difference slope for the fitted r^-p exp(-r/xi) curve.
    axes[0].plot(r+.5, mix["inverse_xi"] + mix["p"]*np.log1p(1/r), "--",
                 color=green, label=r"$r^{-p}e^{-r/\xi}$ fit: $p\approx0.89$, $\xi\approx0.86$")
    axes[0].set(title="Exponential diagnostic: inverse local length",
                ylabel=r"$\kappa_{\rm eff}(r)=\log[C(r)/C(r+1)]$",
                ylim=(0, 1.85))
    axes[0].legend(fontsize=9, loc="lower right")
    axes[0].text(.04, .56, "Exponential: approaches a positive constant\n"
                 "Pure power law: tends toward zero",
                 transform=axes[0].transAxes, fontsize=10)
    axes[1].set(title="Power-law diagnostic: local exponent",
                ylabel=r"$p_{\rm eff}(r)=\log[C(r)/C(r+1)]\,/\,\log[(r+1)/r]$",
                ylim=(0, 24))
    axes[1].text(.04, .8, "A pure power law would give a plateau\n"
                 "The measured exponent keeps increasing",
                 transform=axes[1].transAxes, fontsize=10)
    for ax in axes:
        ax.set(xlabel="Separation along x (simulation cells)",
               xlim=(1, 18), xticks=[2, 6, 10, 14, 18])
        ax.grid(alpha=.2)
    fig.suptitle("Ring–ring bond: evidence for short exponential-like decay\n"
                 "Finite CTM data do not prove the infinite-distance limit", fontsize=14)
    save(fig, "decay-diagnostics")

    cylinder_rows = []
    fig, ax = plt.subplots(figsize=(7.6, 4.8), layout="constrained")
    for sector, color, marker in (("even", blue, "o"), ("odd", orange, "s")):
        lengths = []
        for width in (2, 3, 4):
            data = read(f"cylinders/width{width}-{sector}.json")
            xi = data["channels"]["ring_ring"]["correlation_length_spectrum"]["xi"]
            lengths.append(xi)
            cylinder_rows.append({"width": width, "sector": sector, "bond_xi": xi})
        ax.plot([2, 3, 4], lengths, marker=marker, color=color, lw=1.6,
                label=f"Exact cylinder, {sector} parity")
    ax.axhspan(.8, .9, color=green, alpha=.12,
               label="Plane CTMRG: descriptive 0.8–0.9 range")
    ax.set(title="Independent exact-cylinder check of bond decay",
           xlabel="Cylinder circumference (simulation cells)",
           ylabel=r"Axial correlation length $\xi$ (simulation cells)",
           xticks=[2, 3, 4], xlim=(1.85, 4.15), ylim=(.5, 1.06))
    ax.legend(loc="upper right", fontsize=10)
    ax.text(.03, .035, "Finite widths and fixed boundary parity sectors\n"
            "Sector differences remain; the shaded range is not a confidence interval",
            transform=ax.transAxes, fontsize=9.5)
    ax.grid(alpha=.2)
    save(fig, "cylinder-lengths")

    provenance = {
        "source": "New calculations in this repository; no digitized paper curves",
        "state": selected["state_definition"],
        "distance_unit": "Simulation-cell translation: one small square per cell; along a diagonal of the dummy-contracted square graph",
        "scalar_decay_description": "Characteristic scale 0.8–0.9 cells; finite-distance estimate, not a precise bulk limit",
        "selected_ctm": summary["selected_run"],
        "fit_window": fit_window,
        "orthogonal_opposite_triangle_pairs": {
            "reference_triangle": "E, ring-ring bond in cell (0,0)",
            "translated_triangle": "W, ring-ring bond in cell (r,0)",
            "distance": cross_r, "exact_connected_fractions": [str(c) for c in cross]},
        "exact_cylinder_lengths": cylinder_rows,
        "absolute_plot_floor": floor,
        "largest_distance_in_plane_figures": 18,
        "input_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in input_paths},
        "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "limitations": [
            "Environment and projector truncation are approximate even when CTM iterations converge.",
            "Relative projector cutoff changes the late correlations by several percent.",
            "Exact-cylinder results have unresolved circumference and sector dependence.",
            "Local physical bond, orthogonal occupation, and nonlocal loop-connectivity observables differ."]}
    (out / "figure-summary.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    print("Wrote correlations, decay-diagnostics, cylinder-lengths (PNG/PDF) and figure-summary.json")


if __name__ == "__main__":
    main()
