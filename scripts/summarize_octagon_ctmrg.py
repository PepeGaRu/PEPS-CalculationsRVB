#!/usr/bin/env python3
"""Summarize the declared square-octagon RVB contractions without bulk extrapolation."""
from pathlib import Path
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.diagnostics import diagnose_decay


def main():
    directory = ROOT / "results/octagon-square/ctmrg"
    runs = {path.name: json.loads(path.read_text()) for path in sorted(directory.glob("chi*.json"))}
    selected_name = "chi64_cutoff1e-14_mult1e-18.json"
    selected = runs[selected_name]
    comparison_name = "chi64_cutoff1e-12_mult1e-18.json"
    comparison = runs[comparison_name]
    windows = [(2, 5), (4, 8), (8, 14), (12, 18), (16, 20)]
    summary = {"selected_run": selected_name, "comparison_run": comparison_name,
               "recommended_precision": "Bond decay length about 0.8 square unit cells; descriptive finite-distance estimate, not a converged bulk precision claim.",
               "absolute_correlation_floor": 1e-13, "runs": {}, "selected_observables": {}}
    lines = ["# Conventional singlet RVB on the decorated square-octagon lattice", "",
             "This calculation uses equal coefficients of all fully covering singlet configurations, with the specified triangle and dummy-link singlet arrows. It has eight physical spin-1/2 degrees of freedom per square cell. An unused dummy link contains an internal singlet; there are no vacancies. No semionic loop sign is included.", "",
             "**Result:** the measured connected spin-bond correlations decay rapidly. Across the late resolved fitting windows, a useful descriptive scale is **about 0.8 square unit cells**. This is not a precision thermodynamic-limit estimate. The effective length still drifts with distance and truncation; a power prefactor multiplying an exponential is allowed.", "",
             "## Environment and truncation checks", "",
             "The CTMRG convergence flag refers to three consecutive changes in all corner spectra below the specified tolerance. It does not establish convergence in environment dimension or truncation. The implementation has both a relative projector cutoff and a separate absolute cutoff inside its multiplet-preserving SVD. Equal results at larger nominal chi can conceal equal retained ranks.", "",
             "| Run | chi | Relative cutoff | Multiplet absolute cutoff | Final projector ranks | Iterations | Corner residual / tolerance | Converged |",
             "|---|---:|---:|---:|---|---:|---|---|"]
    for filename, run in runs.items():
        entry = {"chi": run["chi"], "projector_relative_cutoff": run["projector_svd_relative_cutoff"],
                 "multiplet_absolute_cutoff": run["ctm_algorithm_parameters"]["projector_multiplet_abstol"],
                 "projector_effective_ranks": run["projector_effective_ranks_last_sweep"],
                 "ctm_converged": run["ctm_converged"], "ctm_iterations": run["ctm_iterations"],
                 "ctm_residual": run["ctm_convergence_residual"], "ctm_tolerance": run["ctm_convergence_tolerance"],
                 "physical_dimension": run["physical_dimension_after_exact_compression"],
                 "bulk_means": run["bulk_means"]}
        summary["runs"][filename] = entry
        lines.append(f"| [{filename}]({filename}) | {entry['chi']} | {entry['projector_relative_cutoff']:.0e} | {entry['multiplet_absolute_cutoff']:.0e} | {entry['projector_effective_ranks']} | {entry['ctm_iterations']} | {entry['ctm_residual']:.2e} / {entry['ctm_tolerance']:.0e} | {entry['ctm_converged']} |")
    lines += ["", "The raw-physical run retains all 256 local physical states and uses the exact integer tensor divided by 16. It is an independent numerical precision check of the exact rank-41 physical compression, not an independent physical model. Its stricter iteration criterion was not reached; this audit did not cure small-singular-value conditioning.", "",
              "## Measured observables and fitting windows", "",
              "Triangle E contains `(r0,r3,dE)` = physical indices `(0,3,4)`. Ring-ring and ring-dummy observables are `S_i dot S_j`. The horizontal dummy bond is `S_dE(cell) dot S_dW(cell+x)`; the vertical analogue uses `dN,dS`. Bond-center distance is measured in square Bravais cells. Two-cell dummy-bond insertions use separation r>=2 to avoid overlap.", "",
              "All connected covariances use two one-point functions from the same finite CTM strip. Centered operator insertions are independently compared with direct subtraction. The table uses the retained data with absolute correlation greater than 1e-13. Fits describe positive bond correlations; oscillating single-spin data use separate even/odd branches. Spin data below the floor must not be interpreted as long-distance spin correlations.", "",
              f"Selected contraction: [{selected_name}]({selected_name}).", "",
              "| Observable | Requested window | Used distances | Exponential xi | Exponential log RMS | Power-law exponent | Power-law log RMS | Power × exponential p, xi |",
              "|---|---|---|---:|---:|---:|---:|---|"]
    for name, observable in selected["observables"].items():
        rows = observable["rows"]
        r = np.array([row["distance"] for row in rows])
        c = np.array([row["connected_centered"] for row in rows])
        diagnostic = diagnose_decay(r, c, windows=windows, floor=1e-13,
                                    branch=r % 2 if "spin_z" in name else None)
        previous_rows = {row["distance"]: row for row in comparison["observables"][name]["rows"]}
        changes = {}
        for distance in (14, 18):
            index = list(r).index(distance)
            old = previous_rows[distance]["connected_centered"]
            changes[str(distance)] = {"selected_value": float(c[index]), "comparison_value": old,
                                       "relative_change": float((c[index] - old) / old) if old else None,
                                       "above_floor": bool(abs(c[index]) > 1e-13)}
        subtraction_max = max(abs(row["subtraction_disagreement"]) for row in rows)
        summary["selected_observables"][name] = {"diagnostics": diagnostic,
                                                  "cutoff_comparison": changes,
                                                  "maximum_centered_vs_subtracted_difference": subtraction_max}
        if "spin_z" in name or name.endswith("_y"):
            continue
        for fit in diagnostic["windows"]:
            models = fit["models"]
            exponential, power, mixed = models["exponential"], models["power_law"], models["power_exponential"]
            if exponential["status"] != "ok":
                continue
            mixed_text = f"{mixed['p']:.3g}, {mixed['xi']:.3g}" if mixed["status"] == "ok" and mixed["xi"] else "not identified"
            lines.append(f"| {name} | {fit['requested_window']} | {fit['used_distances']} | {exponential['xi']:.4f} | {exponential['log_rms']:.3g} | {power['p']:.3g} | {power['log_rms']:.3g} | {mixed_text} |")
    lines += ["", "## Cutoff sensitivity at fixed distance", "",
              f"Comparison: selected relative cutoff 1e-14 versus [{comparison_name}]({comparison_name}) at 1e-12, both chi=64 and absolute multiplet cutoff 1e-18. These differences are diagnostic, not rigorous error bars.", "",
              "| Observable | C(14) | Relative change at 14 | C(18) | Relative change at 18 |",
              "|---|---:|---:|---:|---:|"]
    for name, item in summary["selected_observables"].items():
        if "spin_z" in name or name.endswith("_y"):
            continue
        c14, c18 = item["cutoff_comparison"]["14"], item["cutoff_comparison"]["18"]
        lines.append(f"| {name} | {c14['selected_value']:.6g} | {100*c14['relative_change']:.2f}% | {c18['selected_value']:.6g} | {100*c18['relative_change']:.2f}% |")
    lines += ["", "The stronger cutoff changes late amplitudes enough that quoting many digits in a correlation length would be misleading. A short correlation length is supported on the resolved range; neither a mathematically exact exponential asymptote nor a universal exponent is established by these finite-environment contractions.", "",
              "## Reproduction", "", "```bash",
              "/workspace/peps-environment/.venv/bin/python scripts/run_octagon.py --chi 64 --distance 24 --projector-tol 1e-14 --multiplet-abstol 1e-18 --ctm-tol 1e-8 --max-iter 120",
              "/workspace/peps-environment/.venv/bin/python scripts/summarize_octagon_ctmrg.py", "```", "",
              "The pinned peps-torch checkout and CPU torch environment are required. Factorized two-cell correlations are validated against dense small-physical-dimension operators in all four directions and against a polarized product state. Source hashes, package versions, raw correlations, convergence histories, final corner spectra, and all algorithm settings are in each JSON file."]
    (directory / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    (directory / "report.md").write_text("\n".join(lines) + "\n")
    print(directory / "report.md")


if __name__ == "__main__":
    main()
