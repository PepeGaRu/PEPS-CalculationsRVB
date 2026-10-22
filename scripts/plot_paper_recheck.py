#!/usr/bin/env python3
"""Plot actual saved calculations and separately identified published points."""
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/paper-recheck"
FLOOR = 1e-13
plt.rcParams.update({"font.size": 11, "axes.titlesize": 13,
                     "axes.spines.top": False, "axes.spines.right": False})


def read(path):
    return json.loads((ROOT / path).read_text())


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def published():
    data = read("results/paper-recheck/published_fig11.json")
    fig, ax = plt.subplots(figsize=(8, 5.4), layout="constrained")
    selected = [0., .1, .2, .3, .4, .5, .6]
    colors = plt.cm.viridis(np.linspace(.05, .9, len(selected)))
    for theta, color in zip(selected, colors):
        for dim in (40, 60):
            curve = next(c for c in data["curves"] if c["theta"] == theta
                         and c["boundary_dimension"] == dim)
            r = [p["r"] for p in curve["points"]]
            c = [p["correlation"] for p in curve["points"]]
            ax.semilogy(r, c, "o-" if dim == 60 else "--", color=color,
                        alpha=1 if dim == 60 else .5,
                        label=fr"$\theta={theta:g}$" if dim == 60 else None)
    ax.set(xlim=(.95, 4.15), ylim=(1e-6, .1), xticks=[1, 2, 3, 4],
           xlabel="Separation r (paper's blocked cells)",
           ylabel="Absolute connected correlation",
           title="Published kagome result: Schuch et al., Fig. 11")
    ax.text(.98, .97, "Points recovered from the PDF drawing\nSolid: boundary D=60; dashed: D=40\nPaper operator: z ⊗ z, z=diag(1,−1,1)",
            transform=ax.transAxes, ha="right", va="top", fontsize=10)
    ax.text(.04, .06, "θ=0: singlet RVB\nθ → 1: orthogonal dimers\nClipped points are omitted, not zero.",
            transform=ax.transAxes, fontsize=10)
    ax.legend(loc="center right", ncol=2, fontsize=10)
    ax.grid(alpha=.2)
    fig.get_layout_engine().set(h_pad=.1)
    save(fig, "published-fig11")


def comparison():
    exact = read("results/reference/data.json")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), layout="constrained")
    for chi in (32, 64, 128):
        d = read(f"results/upstream/square_chi{chi}.json")
        rows = d["correlations"]["DD"]
        r = np.array([p["distance"] for p in rows])
        c = np.abs([p["dd_connected"] for p in rows])
        mask = (r % 2 == 0) & (c > FLOOR)
        axes[0].loglog(r[mask], c[mask], ".-", label=f"Singlets, χ={chi}")
    d = exact["square_orthogonal_plane"]
    r = np.array(d["distance"]); c = np.abs(d["connected_correlation"])
    mask = (r >= 2) & (r % 2 == 0)
    axes[0].loglog(r[mask], c[mask], "k--", label="Orthogonal dimers: exact plane")
    axes[0].text(.04, .08, "Orthogonal: asymptotic exponent 2\nSinglets: finite-χ estimate p ≈ 1.30\n(r=6–20); bulk exponent unresolved",
                 transform=axes[0].transAxes, fontsize=10)
    axes[0].set(title="Square lattice: algebraic dimer decay",
                xlabel="Distance (nearest-neighbor spacings)",
                ylabel="|Connected correlation|", ylim=(3e-6, .08))
    axes[0].set_xticks([2, 4, 10, 20, 40, 100]); axes[0].xaxis.set_major_formatter(ScalarFormatter())
    for chi in (18, 32):
        d = read(f"results/upstream/kagome_chi{chi}.json")
        rows = d["correlations"]["DD_intra[(0, 0),(1, 0)]"]
        r = np.array([p["distance"] for p in rows]); c = np.abs([p["dd_connected"] for p in rows])
        mask = (r <= 14) & (c > FLOOR)
        axes[1].semilogy(r[mask], c[mask], "o-", markersize=4, label=f"Singlets, χ={chi}")
    fit = (r >= 8) & (r <= 14) & (c > FLOOR)
    slope, intercept = np.polyfit(r[fit], np.log(c[fit]), 1)
    axes[1].semilogy(r[fit], np.exp(intercept+slope*r[fit]), "k--", label=f"Baseline fit: ξ={-1/slope:.3f} cells")
    axes[1].text(.04, .42, "Projector cutoff 10⁻⁸\nξ remains cutoff dependent",
                 transform=axes[1].transAxes, fontsize=10)
    axes[1].text(.04, .08, "Orthogonal occupation correlations:\nexactly zero for separated bonds\n(zero cannot appear on a log axis)",
                 transform=axes[1].transAxes, fontsize=10)
    axes[1].set(title="Kagome lattice: short dimer decay",
                xlabel="Distance (unit cells; 1 cell = 2 NN spacings)",
                ylabel="|Connected correlation|", ylim=(1e-12, .1))
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=9, loc="upper right")
    fig.suptitle("Our calculations — finite-χ singlets and exact orthogonal dimers", fontsize=13)
    save(fig, "square-kagome-comparison")


def diagnostics():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    audit = read("results/paper-recheck/tail_audit.json")
    best = max(audit["kagome"], key=lambda d: d["chi"])
    for direction, label in (("(1, 0)", "x translation"), ("(0, 1)", "y translation")):
        vals = best["observables"][f"DD_intra[(0, 0),{direction}]"]["local_effective_decay"]
        vals = [v for v in vals if v["distances"][1] <= 14]
        axes[0].plot([np.mean(v["distances"]) for v in vals], [v["xi"] for v in vals], "o-", label=label)
    axes[0].set(title="A short fit can miss the tail", xlabel="Distance (unit cells)",
                ylabel="Local effective decay length (cells)", ylim=(.3, .7))
    axes[0].text(.03, .07, "x approaches ≈0.638 in the resolved window;\ny has not settled before numerical noise.", transform=axes[0].transAxes, fontsize=10)
    data = read("results/reference/data.json")["singlet_cylinders"]
    for sector in ("even", "odd"):
        vals = [d for d in data if d["lattice"] == "kagome" and d["metadata"]["virtual_vacuum_parity"] == sector]
        axes[1].plot([d["width"] for d in vals], [d["correlation_length"] for d in vals], "o-", label=f"{sector} virtual parity")
    axes[1].set(title="Kagome cylinders: keep the sector fixed", xlabel="Circumference (unit cells)",
                ylabel="Bond-coupled transfer length (axial cells)", xticks=[2, 4], ylim=(.25, .85))
    axes[1].text(.03, .06, "Exact contractions of finite cylinders;\nthese widths do not determine the plane limit.", transform=axes[1].transAxes, fontsize=10)
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=9)
    save(fig, "tail-and-sector-checks")


def interpolation():
    files = list((OUT / "ctmrg").glob("theta*_chi32.json"))
    if not files:
        return
    data = sorted((json.loads(p.read_text()) for p in files), key=lambda d: d["theta"])
    fig, ax = plt.subplots(figsize=(8, 5.3), layout="constrained")
    for d in data:
        if not d["ctm_converged"]:
            continue
        rows = d["observables"]["zz_01_y"]["rows"]
        r = np.array([p["distance"] for p in rows]); c = np.abs([p["connected_centered"] for p in rows])
        mask = c > FLOOR
        if d["theta"] < 1:
            ax.semilogy(r[mask], c[mask], "o-", label=fr"$\theta={d['theta']:g}$")
    ax.set(title="New Eq. (12) calculation: physical spin ZZ bonds", xlabel="Separation along y (kagome unit cells)",
           ylabel="|Connected correlation|", ylim=(FLOOR, .1), xlim=(.7, 12.3))
    ax.text(.98, .97, "χ=32 CTMRG, projector cutoff 10⁻⁸\nBond (0,1) in each up triangle; O = Z₀Z₁\nZ = Iflag ⊗ Pauli σz\nNot identified with the paper's Fig. 11 operator",
            transform=ax.transAxes, ha="right", va="top", fontsize=10)
    endpoint = next((d for d in data if d["theta"] == 1 and d["ctm_converged"]), None)
    if endpoint:
        residual = max(abs(p["connected_centered"]) for p in endpoint["observables"]["zz_01_y"]["rows"])
        ax.text(.04, .06, f"θ=1: all displayed separations below {FLOOR:g}\n(maximum residual {residual:.1e}); omitted from log axis",
                transform=ax.transAxes, fontsize=10,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    ax.legend(loc="center right"); ax.grid(alpha=.2)
    save(fig, "physical-interpolation")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--published-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    published()
    if not args.published_only:
        comparison(); diagnostics(); interpolation()
    print(OUT)
