#!/usr/bin/env python3
"""Compact collaborator figure from saved square–octagon spin-half RVB data.

No published points, fresh contractions, or asymptotic extrapolations are used.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/octagon-square/ctmrg"
OUTPUT = ROOT / "results/collaborator-notes"


def main():
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary_path = INPUT / "summary.json"
    summary = json.loads(summary_path.read_text())
    paths = [INPUT / summary["selected_run"], INPUT / summary["comparison_run"]]
    runs = [json.loads(path.read_text()) for path in paths]
    key = "triangle_ring_ring_x"
    item = summary["selected_observables"][key]
    windows = item["diagnostics"]["windows"]
    fitted = next(window for window in windows if window["requested_window"] == [12., 18.])
    early = next(window for window in windows if window["requested_window"] == [8., 14.])
    floor = summary["absolute_correlation_floor"]
    data = []
    for run in runs:
        rows = run["observables"][key]["rows"]
        r = np.array([row["distance"] for row in rows])
        c = np.array([row["connected_centered"] for row in rows])
        keep = (r >= 2) & (r <= 18) & (c > floor)
        r, c = r[keep], c[keep]
        midpoint = .5 * (r[:-1] + r[1:])
        kappa = np.log(c[:-1] / c[1:]) / np.diff(r)
        data.append(dict(distance=r.tolist(), correlation=c.tolist(),
                         midpoint=midpoint.tolist(), inverse_xi_eff=kappa.tolist(),
                         relative_cutoff=run["projector_svd_relative_cutoff"],
                         effective_ranks=run["projector_effective_ranks_last_sweep"]))

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9,
                         "legend.fontsize": 7.1, "axes.labelsize": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                         "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.45))
    fig.subplots_adjust(left=.087, right=.985, bottom=.22, top=.88, wspace=.32)
    colors, markers = ("#17619a", "#d18331"), ("o", "s")
    for row, color, marker in zip(data, colors, markers):
        cutoff = int(round(np.log10(row["relative_cutoff"])))
        label = rf"cutoff $10^{{{cutoff}}}$"
        axes[0].semilogy(row["distance"], row["correlation"], color=color,
                         marker=marker, ms=2.7, lw=.95, label=label)
        axes[1].plot(row["midpoint"], row["inverse_xi_eff"], color=color,
                     marker=marker, ms=2.7, lw=.95)
    fit = fitted["models"]["exponential"]
    rfit = np.linspace(12, 18, 80)
    axes[0].semilogy(rfit, np.exp(fit["log_amplitude"] - rfit / fit["xi"]),
                     "--", color="black", lw=1., label=r"fit: $\xi=0.82$ cells")
    axes[0].axvspan(12, 18, color=".5", alpha=.08, lw=0)
    axes[0].set(title="Physical singlet-bond correlations", xlim=(1.7, 18.4),
                ylim=(1e-13, 4e-3), ylabel=r"$C_B(r)$", yticks=[1e-3, 1e-6, 1e-9, 1e-12])
    axes[0].legend(loc="upper right", frameon=False, handlelength=1.7,
                   labelspacing=.25, borderpad=.15)
    axes[0].text(.035, .09, r"$B=\mathbf{S}_i\!\cdot\!\mathbf{S}_j$"
                 + "\n" + r"PEPS $D=3$; CTMRG $\chi=64$",
                 transform=axes[0].transAxes, fontsize=7.5)
    axes[1].axhspan(1/.9, 1/.8, color="#4c9976", alpha=.13, lw=0)
    axes[1].set(title="Local decay diagnostic", xlim=(1.7, 18.4), ylim=(0, 1.85),
                ylabel=r"$\kappa_{\rm eff}=\log[C_B(r)/C_B(r+1)]$")
    axes[1].text(.04, .49, r"Resolved scale: $\xi\sim0.8$–$0.9$ cells",
                 transform=axes[1].transAxes, fontsize=7.6)
    axes[1].text(.04, .32, r"Exponential: $\kappa_{\rm eff}\to1/\xi>0$" + "\n"
                 + r"Pure power law: $\kappa_{\rm eff}\to0$",
                 transform=axes[1].transAxes, fontsize=7.3, linespacing=1.5)
    axes[1].text(.04, .10, "Finite-window estimate; bulk tail unconverged.",
                 transform=axes[1].transAxes, fontsize=7.2, color=".30")
    for ax in axes:
        ax.set_xticks([2, 6, 10, 14, 18])
        ax.set_xlabel("Separation r (simulation-cell translations)", labelpad=3)
        ax.grid(alpha=.15)
    for extension in ("pdf", "png"):
        fig.savefig(OUTPUT / f"correlations-2d.{extension}", dpi=240)
    plt.close(fig)

    caption = (
        "We measure the connected physical bond covariance "
        "C_B(r)=<B_0 B_r>-<B_0><B_r>, B=S_i·S_j, using a D=3 PEPS and "
        "infinite-plane CTMRG contraction (χ=64). Boundary-contraction and "
        "sector-resolved transfer diagnostics follow the RVB-PEPS framework "
        "of Schuch et al. [1,2]; these are our model’s data. Over r=12–18, "
        "exponential and power-times-exponential fits give ξ≈0.82 and 0.86 "
        "simulation cells. The fitted pure-power exponent drifts from 13.4 "
        "to 18.1. Cutoff changes alter C_B(18) by 8%, so ξ≈0.8–0.9 is "
        "a finite-window decay scale, not a precise bulk asymptote. "
        "One translation follows a diagonal (1,−1) of the contracted square grid."
    )
    references = [
        {"number": 1, "authors": "N. Schuch, D. Poilblanc, J. I. Cirac, D. Pérez-García",
         "title": "Resonating valence bond states in the PEPS formalism",
         "year": 2012, "arxiv": "1203.4816", "url": "https://arxiv.org/abs/1203.4816",
         "method_note": "Fig. 11 uses a boundary iMPS/iTEBD contraction; Fig. 12 separates leading topological transfer sectors. This figure uses CTMRG instead."},
        {"number": 2, "authors": "D. Poilblanc, N. Schuch, D. Pérez-García, J. I. Cirac",
         "title": "Topological and entanglement properties of resonating valence bond wave functions",
         "year": 2012, "arxiv": "1202.0947", "url": "https://arxiv.org/abs/1202.0947",
         "method_note": "Finite-cylinder sector control supplies an independent diagnostic, not a replacement for bulk convergence."},
    ]
    audit = dict(
        observable="triangle ring–ring physical bond covariance, translated along simulation x",
        state="Conventional equal-coefficient fully covered spin-half singlet RVB; no semionic phase",
        unit="One translation between adjacent small-square centres; contracted grid displacement (1,-1)",
        contraction="D=3 PEPS; exact local physical support rank41; CTMRG chi64",
        data=data, floor=floor, floor_is_systematic_error_bound=False,
        fit_window=[12, 18], fitted_models=fitted["models"],
        earlier_window=[8, 14], earlier_models=early["models"],
        cutoff_comparison=item["cutoff_comparison"],
        caption=caption, references=references,
        limits="Finite-distance and cutoff dependence remain; no precise bulk xi or exclusion of a weaker tail.",
        sources=[dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                 for p in [*paths, summary_path]],
    )
    (OUTPUT / "correlations-2d-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    (OUTPUT / "correlations-2d-caption.md").write_text(caption + "\n\n" + "\n".join(
        f'[{r["number"]}] {r["authors"]}, *{r["title"]}* ({r["year"]}), '
        f'[arXiv:{r["arxiv"]}]({r["url"]}).' for r in references) + "\n")


if __name__ == "__main__":
    main()
