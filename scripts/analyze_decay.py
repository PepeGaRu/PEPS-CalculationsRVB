#!/usr/bin/env python3
"""Apply descriptive decay diagnostics to saved calculations, not paper curves."""
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
from rvb.orthogonal import square_plane_benchmark


def load(relative):
    path = ROOT / relative
    return json.loads(path.read_text()), {
        "path": relative, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def assemble():
    results = {}
    square_windows = [(4, 12), (6, 20), (10, 30), (20, 40)]
    kagome_windows = [(2, 8), (4, 10), (8, 14)]

    def add(key, r, c, *, windows, source, description, parity=False,
            exact_zero=False, metadata=None):
        r, c = np.asarray(r), np.asarray(c)
        keep = r > 0
        r, c = r[keep], c[keep]
        result = diagnose_decay(r, c, windows=windows, floor=1e-13,
                                branch=(r.astype(int) % 2).tolist() if parity else None,
                                exact_zero=exact_zero)
        results[key] = {"description": description, "source": source,
                        "metadata": metadata or {}, "diagnosis": result}

    for chi in (16, 32, 64, 128):
        d, source = load(f"results/upstream/square_chi{chi}.json")
        for observable, column in (("DD", "dd_connected"), ("SS", "ss")):
            rows = d["correlations"][observable]
            add(f"square_{observable}_chi{chi}", [p["distance"] for p in rows],
                [p[column] for p in rows], windows=square_windows, source=source,
                description="Square collinear scalar spin bonds" if observable == "DD" else "Square spin correlation",
                parity=True, metadata={"chi": chi, "ctm_converged": d["ctm_converged"],
                                       "units": "nearest-neighbor spacing", "bulk_converged": False})
    for chi in (18, 32):
        d, source = load(f"results/upstream/kagome_chi{chi}.json")
        for direction, vector in (("x", "(1, 0)"), ("y", "(0, 1)")):
            rows = d["correlations"][f"DD_intra[(0, 0),{vector}]"]
            add(f"kagome_{direction}_chi{chi}", [p["distance"] for p in rows],
                [p["dd_connected"] for p in rows], windows=kagome_windows, source=source,
                description=f"Kagome scalar spin bond, original upstream {direction} convention",
                metadata={"chi": chi, "ctm_converged": d["ctm_converged"],
                          "projector_cutoff": 1e-8, "units": "Bravais unit cell", "bulk_converged": False})
    d, source = load("results/paper-recheck/ctmrg/strict-cutoff/theta0_chi32.json")
    rows = d["observables"]["ss_01_y"]["rows"]
    add("kagome_x_strict_chi32", [p["distance"] for p in rows],
        [p["connected_centered"] for p in rows], windows=kagome_windows, source=source,
        description="Kagome same scalar bond; new block y maps to original upstream x",
        metadata={"chi": 32, "ctm_converged": d["ctm_converged"],
                  "projector_cutoff": d["projector_svd_relative_cutoff"],
                  "units": "Bravais unit cell", "bulk_converged": False,
                  "purpose": "Unconverged cutoff-sensitivity diagnostic, not an improved estimate"})
    exact, source = load("results/reference/data.json")
    d = exact["square_orthogonal_plane"]
    add("square_orthogonal_collinear", d["distance"], d["connected_correlation"],
        windows=[(16, 32), (32, 64), (64, 128)], source=source,
        description="Exact square-plane occupations: collinear horizontal bonds", parity=True)
    # Leading parity cancellation is an independent exact counterexample to
    # fitting all directions/branches with a single exponent.
    transverse = square_plane_benchmark(128, orientation="vertical")
    add("square_orthogonal_transverse", transverse["distance"], transverse["connected_correlation"],
        windows=[(16, 32), (32, 64), (64, 128)],
        source={"method": "New evaluation of exact Kasteleyn quadrature",
                "code_sha256": hashlib.sha256((ROOT / "rvb/orthogonal.py").read_bytes()).hexdigest()},
        description="Vertical occupation bonds separated horizontally; even branch r^-4, odd r^-2",
        parity=True)
    d = exact["kagome_orthogonal"]
    add("kagome_orthogonal", d["distance"], d["connected_correlation"],
        windows=[(1, 4)], source=source, exact_zero=True,
        description="Exact arrow-constraint occupation zeros for supplied nonzero separations")
    return results, transverse


def make_plot(data, output):
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), layout="constrained")
    for chi in (32, 64, 128):
        rows = [r for r in data[f"square_DD_chi{chi}"]["diagnosis"]["local"]
                if r["branch"] == "0" and r["r_right"] <= 40]
        axes[0].semilogx([r["r_geometric"] for r in rows], [r["p_eff"] for r in rows],
                         "o-", ms=3, label=f"Our singlets, χ={chi}")
    rows = [r for r in data["square_orthogonal_collinear"]["diagnosis"]["local"]
            if r["branch"] == "0" and 8 <= r["r_left"] <= 40]
    axes[0].semilogx([r["r_geometric"] for r in rows], [r["p_eff"] for r in rows],
                     "k--", label="Exact orthogonal dimers")
    axes[0].axhspan(1.12, 1.20, color="gray", alpha=.22)
    axes[0].set(title="Square: look for a stable exponent", xlabel="Separation (nearest-neighbor spacings)",
                ylabel="Effective exponent p", xlim=(3, 40), ylim=(.8, 4.3))
    axes[0].set_xticks([4, 6, 10, 20, 40]); axes[0].xaxis.set_major_formatter(ScalarFormatter())
    axes[0].xaxis.set_minor_formatter(NullFormatter())
    axes[0].text(.04, .56, "Gray band: paper's 1.16 ± 0.04\nDifferent displacement from our curves\nUpward drift: finite-χ cutoff", transform=axes[0].transAxes, fontsize=10)
    for key, label, style in (("kagome_x_chi18", "χ=18, cutoff 10⁻⁸", "o-"),
                              ("kagome_x_chi32", "χ=32, cutoff 10⁻⁸", "s-"),
                              ("kagome_x_strict_chi32", "χ=32, cutoff 10⁻¹² (unconverged)", "x--")):
        rows = [r for r in data[key]["diagnosis"]["local"] if r["r_right"] <= 14]
        axes[1].plot([r["r_geometric"] for r in rows], [r["inverse_xi_eff"] for r in rows],
                     style, ms=4, label=label)
    axes[1].set(title="Kagome: look for a stable inverse length", xlabel="Separation (Bravais unit cells)",
                ylabel="Effective inverse length κ = 1/ξ (cells⁻¹)", xlim=(1, 14), ylim=(1.4, 1.85))
    axes[1].text(.04, .06, "A plateau supports exponential decay.\nCutoff dependence still limits precision.\n1 cell = 2 nearest-neighbor spacings.", transform=axes[1].transAxes, fontsize=10)
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=9, loc="upper left")
    fig.suptitle("New diagnostics of our computed correlations — no extracted paper curves", fontsize=13)
    fig.savefig(output / "effective-decay.png", dpi=160, bbox_inches="tight")
    fig.savefig(output / "effective-decay.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(data, output):
    lines = ["# Decay diagnostics of actual calculations", "",
             "The second figure in the earlier answer was Schuch's published Figure 11 redrawn from PDF coordinates. It was not our reproduction. This report analyzes our saved numerical correlations and exact orthogonal benchmarks; no extracted paper curve is fitted.", "",
             "![Effective exponent and inverse length](effective-decay.png)", "",
             "The gray band is only Albuquerque–Alet's reported exponent 1.16(4), not sampled paper data. Their principal displacement differs from the collinear bonds measured here. The colored lines are diagnostics of our numerical data. The stricter kagome curve did not meet its iteration convergence criterion and diagnoses sensitivity rather than supplying a more accurate result.", "",
             "## Descriptive fit windows", "",
             "Fits below use equal weights in log magnitude. Log RMS residuals are fit-shape diagnostics, not independent statistical errors or proof of a bulk decay law. Branch 0/1 denotes even/odd separations. Compare windows and convergence settings, not just residuals.", "",
             "| Dataset | Branch | Actual range | Points | Power exponent | Power log RMS | Exponential length | Exp. log RMS |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for key in ("square_DD_chi64", "square_DD_chi128", "square_orthogonal_collinear",
                "square_orthogonal_transverse", "kagome_x_chi32", "kagome_y_chi32", "kagome_x_strict_chi32"):
        for window in data[key]["diagnosis"]["windows"]:
            if key.startswith("square_DD") and window["branch"] != "0":
                continue
            power = window["models"].get("power_law", {})
            exponential = window["models"].get("exponential", {})
            if power.get("p") is None or exponential.get("xi") is None:
                continue
            used = window["used_distances"]
            lines.append(f"| {key} | {window['branch']} | {used[0]:g}–{used[-1]:g} | {len(used)} | {power['p']:.5g} | {power['log_rms']:.3g} | {exponential['xi']:.5g} | {exponential['log_rms']:.3g} |")
    lines += ["", "## Interpretation and limits", "",
              "- Square orthogonal collinear bonds approach exponent 2. Transverse bonds give exponent 4 on the even branch and 2 on the odd branch: a real leading-term cancellation, not a failure of the dimer model.",
              "- Square singlet curves still move with chi. Their finite-window exponent near 1.30 is not a reproduction of 1.16(4), and an eventual finite-chi exponential tail is not a measured finite bulk dimer length.",
              "- Kagome's approximately constant inverse length supports rapid exponential decay in the resolved range. Projector-cutoff sensitivity and unconverged stricter runs prevent a precise bulk-length claim.",
              "- The exact kagome occupation zeros are supplied as exact combinatorial results. The generic numerical routine does not infer exact finite range from floating-point zeros.",
              "- A mixed power-times-exponential fit is included in JSON with a scale-independent conditioning diagnostic. It has an extra parameter; a lower residual alone cannot establish its physical meaning.", "",
              "The numerical floor is 1e-13 for this float64 dataset. It is a screening threshold, not a bound on CTMRG truncation error. Optional externally estimated absolute errors are supported by the reusable API. Signs, branch labels, excluded points and segment boundaries remain in the output; local derivatives never bridge unresolved points within a branch.", "",
              "See [the diagnosis protocol](../../docs/decay-diagnostics.md), [the new paper audit](../../docs/paper-1009-audit.md), and [the full machine-readable output](data.json). Source hashes and convergence flags are included. No new infinite-plane CTMRG run was required for this analysis.", "",
              "Repeat with `python scripts/analyze_decay.py`. The exact transverse orthogonal benchmark is reevaluated; all singlet inputs come from saved runs. For the next lattice, specify the state/covering superposition, allowed bonds, operator, direction, units and boundaries before fitting a decay."]
    (output / "report.md").write_text("\n".join(lines) + "\n")


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/decay-diagnostics")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data, transverse = assemble()
    payload = {"scope": "Descriptive finite-data diagnostics; no automatic thermodynamic classification",
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "diagnostics_code_sha256": hashlib.sha256((ROOT / "rvb/diagnostics.py").read_bytes()).hexdigest(),
               "datasets": data, "new_exact_transverse_benchmark": transverse}
    (args.output / "data.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    make_plot(data, args.output)
    write_report(data, args.output)
    print(args.output / "report.md")


if __name__ == "__main__":
    main()
