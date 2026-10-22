#!/usr/bin/env python3
"""Reproduce fixed nearest-neighbor RVB states using unmodified peps-torch.

The worker adds connected bond correlations through upstream contraction APIs.
It does not optimize tensors or edit the upstream checkout.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import runpy
import subprocess
import sys
import time

UPSTREAM_URL = "https://github.com/jurajHasik/peps-torch.git"
UPSTREAM_COMMIT = "239a2e717a567df0ece80a5bce5443a6f78d048b"
ROOT = Path(__file__).resolve().parents[1]
KAGOME_CORRELATION_FLOOR = 1e-13


def worker(lattice: str, source: Path, argv: list[str]) -> None:
    """Run the real upstream example, extending observables without editing it."""
    example = source / "examples" / ("j1j2" if lattice == "square" else "kagome")
    filename = "ctmrg_j1j2_c4v.py" if lattice == "square" else "ctmrg_spin_half_kagome.py"
    sys.path[:0] = [str(source), str(example)]
    sys.argv = [str(example / filename), *argv]
    namespace = runpy.run_path(str(example / filename), run_name="rvb_reference")
    import torch

    if lattice == "square":
        from models.j1j2 import J1J2_C4V_BIPARTITE
        from ctm.one_site_c4v import corrf_c4v
        from groups import su2

        original = J1J2_C4V_BIPARTITE.eval_corrf_DD_H

        def connected(self, state, env, dist, verbosity=0):
            result = original(self, state, env, dist, verbosity=verbosity)
            rot = su2.get_rot_op(self.phys_dim, dtype=self.dtype, device=self.device)
            bond = torch.einsum("ki,kjcb,ca->ijab", rot, self.SS, rot)
            identity = torch.eye(4, dtype=bond.dtype, device=bond.device).reshape(2, 2, 2, 2)
            bond_other = bond.permute(1, 0, 3, 2).contiguous()
            get_bond = lambda r: bond if r % 2 == 0 else bond_other
            contract = corrf_c4v.corrf_2sOH2sOH_E1
            left = contract(state, env, bond, lambda r: identity, dist)
            right = contract(state, env, identity, get_bond, dist)
            result.update(d_left=left, d_right=right, dd_connected=result["dd"] - left * right)
            return result

        J1J2_C4V_BIPARTITE.eval_corrf_DD_H = connected
    else:
        from models.spin_half_kagome import S_HALF_KAGOME
        from ctm.generic import corrf

        original = S_HALF_KAGOME.eval_corrf_SS

        def connected(self, coord, direction, state, env, dist, site=0):
            result = original(self, coord, direction, state, env, dist, site=site)
            if site == 0:
                # Physical axes are the three spins in the blocked down triangle.
                # SSnnId acts on spins (0,1), with identity on spin 2; j1=1, JD=0.
                bond = self.SSnnId.reshape(8, 8)
                identity = torch.eye(8, dtype=bond.dtype, device=bond.device)
                contract = corrf.corrf_1sO1sO
                dd = contract(coord, direction, state, env, bond, lambda r: bond, dist)
                left = contract(coord, direction, state, env, bond, lambda r: identity, dist)
                right = contract(coord, direction, state, env, identity, lambda r: bond, dist)
                print(f"\nDD_intra[{coord},{direction}] r dd d_left d_right dd_connected")
                for i in range(dist):
                    print(i, dd[i].item(), left[i].item(), right[i].item(), (dd[i] - left[i] * right[i]).item())
            return result

        S_HALF_KAGOME.eval_corrf_SS = connected

    with torch.no_grad():
        namespace["main"]()


def number(value: str) -> float:
    z = complex(value)
    if abs(z.imag) > 1e-8:
        raise ValueError(f"Unexpected non-real observable: {value}")
    return float(z.real)


def parse_log(path: Path, lattice: str, chi: int) -> dict:
    lines = path.read_text().splitlines()
    result = {"lattice": lattice, "chi": chi, "correlations": {}, "transfer_spectra": {}}
    for i, line in enumerate(lines):
        iteration = re.fullmatch(r"(\d+), ([\deE.+-]+)", line)
        if iteration:
            result["ctm_last_step"] = int(iteration[1])
            result["ctm_convergence_residual"] = float(iteration[2])
        if line.startswith("FINAL "):
            result["energy_per_site"] = number(line[6:].split(", ")[0])
        if line.startswith("TIMINGS "):
            result["timings"] = line
        if line.startswith(("SS r ", "SS[", "DD r ", "DD_intra[")):
            name, columns = line.split(" r ", 1)
            columns = columns.split(" canonical")[0].split()
            rows = []
            for following in lines[i + 1:]:
                if not re.match(r"^\d+ ", following):
                    break
                fields = following.split()
                row = {key: number(value) for key, value in zip(columns, fields[1:])}
                # Upstream index 0 means adjacent sites, or two nonoverlapping bonds.
                row["distance"] = int(fields[0]) + (2 if name == "DD" else 1)
                rows.append(row)
            result["correlations"][name] = rows
        if line.startswith("spectrum(T)"):
            spectrum = []
            for following in lines[i + 1:]:
                if not re.match(r"^\d+ ", following):
                    break
                _, real, imag = following.split()
                spectrum.append([float(real), float(imag)])
            result["transfer_spectra"][line] = spectrum
    for key, spectrum in result["transfer_spectra"].items():
        ratio = abs(complex(*spectrum[1])) / abs(complex(*spectrum[0]))
        result.setdefault("leading_transfer_length", {})[key] = -1 / math.log(ratio)
    result["distance_units"] = "nearest-neighbor square spacing" if lattice == "square" else "kagome Bravais unit cells (2 nearest-neighbor spacings)"
    result["connected_definition"] = "<D(0)D(r)> - <D(0)><D(r)>; both onepoints contracted in the same finite CTM strip"
    raw_tolerance = next(float(line.split("=", 1)[1]) for line in lines if line.startswith("ctm_conv_tol="))
    check = "rdm2x1" if lattice == "square" else next(line.split("=", 1)[1].strip() for line in lines if line.startswith("CTM_check="))
    tolerance = raw_tolerance * (2 if check == "Partial_energy" else 1)
    result["ctm_convergence_criterion"] = check
    result["ctm_convergence_tolerance"] = tolerance
    result["ctm_converged"] = result.get("ctm_convergence_residual", float("inf")) < tolerance
    result["transfer_length_caveat"] = "Leading transfer gap is not necessarily the dimer channel. Finite chi imposes a spurious finite length on critical square RVB."
    import numpy as np
    for name, rows in result["correlations"].items():
        if not name.startswith("DD"):
            continue
        lower, upper = (6, 20) if lattice == "square" else (2, 8)
        selected = [row for row in rows if lower <= row["distance"] <= upper
                    and (lattice != "square" or row["distance"] % 2 == 0)]
        r = np.array([row["distance"] for row in selected])
        c = np.abs([row["dd_connected"] for row in selected])
        slope, intercept = np.polyfit(np.log(r) if lattice == "square" else r, np.log(c), 1)
        result.setdefault("finite_window_dimer_fits", {})[name] = {
            "window": [lower, upper], "even_distances_only": lattice == "square",
            "power_law_exponent" if lattice == "square" else "exponential_length": float(-slope if lattice == "square" else -1 / slope),
            "caveat": "Descriptive finite-window fit, not an extrapolated infinite-chi result; compare multiple chi and windows."
        }
    return result


def plot_results(results: list[dict], destination: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(destination / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(destination / ".matplotlib/cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    import numpy as np

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    square_chis = sorted(result["chi"] for result in results if result["lattice"] == "square")
    for result in results:
        chi = result["chi"]
        if result["lattice"] == "square":
            color = colors[square_chis.index(chi) % len(colors)]
            rows = result["correlations"]["DD"]
            r = np.array([row["distance"] for row in rows])
            c = np.abs([row["dd_connected"] for row in rows])
            # Show each parity separately: the dimer channel is staggered.
            for parity in (0, 1):
                mask = r % 2 == parity
                axes[0].loglog(r[mask], c[mask], ".-", color=color, alpha=1 if parity == 0 else .45,
                               label=f"chi={chi}" if parity == 0 else None)
            spin = result["correlations"]["SS"]
            spin_visible = [row for row in spin if abs(row["ss"]) >= 1e-12]
            axes[1].semilogy([row["distance"] for row in spin_visible], [abs(row["ss"]) for row in spin_visible], color=color, label=f"chi={chi}")
        else:
            rows = result["correlations"]["DD_intra[(0, 0),(1, 0)]"]
            r = np.array([row["distance"] for row in rows])
            c = np.abs([row["dd_connected"] for row in rows])
            mask = c >= KAGOME_CORRELATION_FLOOR
            axes[2].semilogy(r[mask], c[mask], ".-", label=f"chi={chi}")
    axes[2].axhline(KAGOME_CORRELATION_FLOOR, color=".5", linestyle=":", label="numerical cutoff: 1e-13")
    titles = ["Square RVB: connected scalar-bond correlations", "Square RVB: spin correlations", "Kagome RVB: scalar-bond correlations, x"]
    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("distance (lattice spacings)" if ax is not axes[2] else "distance (Bravais unit cells)")
        ax.set_ylabel("absolute correlation")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
        ax.grid(True, alpha=.2)
    axes[1].set_ylim(1e-12, 1)
    axes[2].set_ylim(KAGOME_CORRELATION_FLOOR / 2, .1)
    axes[0].set_xticks([2, 4, 10, 20, 40, 80])
    axes[0].xaxis.set_major_formatter(ScalarFormatter())
    axes[0].xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle("Fixed D=3 singlet RVB tensors; approximate infinite-plane CTMRG contraction")
    fig.tight_layout()
    fig.savefig(destination / "correlations.png", dpi=180)
    plt.close(fig)


def exponential_window(rows: list[dict], lower: int, upper: int) -> dict:
    """Describe a resolved window without interpreting it as an asymptotic tail."""
    import numpy as np

    selected = [row for row in rows if lower <= row["distance"] <= upper]
    result = {"window": [lower, upper], "xi": None,
              "absolute_numerical_floor": KAGOME_CORRELATION_FLOOR}
    if {row["distance"] for row in selected} != set(range(lower, upper + 1)):
        return {**result, "reason": "requested window not covered by saved data"}
    r = np.array([row["distance"] for row in selected])
    c = np.abs([row["dd_connected"] for row in selected])
    result["minimum_absolute_correlation"] = float(c.min())
    if np.any(c < KAGOME_CORRELATION_FLOOR):
        return {**result, "reason": "window reaches the conservative numerical floor"}
    slope, intercept = np.polyfit(r, np.log(c), 1)
    result["log_rms_residual"] = float(np.sqrt(np.mean((np.log(c) - intercept - slope * r) ** 2)))
    if slope >= 0:
        return {**result, "reason": "window is not decaying"}
    return {**result, "xi": float(-1 / slope),
            "caveat": "Descriptive finite-window fit at fixed chi; no asymptotic or infinite-chi extrapolation."}


def write_report(results: list[dict], destination: Path) -> None:
    import numpy as np

    metadata_path = destination / "manifest.json"
    versions = json.loads(metadata_path.read_text()).get("versions", {}) if metadata_path.exists() else {}
    version_text = ", ".join(f"{name} {versions.get(name, 'unknown')}" for name in ["torch", "numpy", "scipy"])

    lines = [
        "# Independent peps-torch RVB cross-check",
        "",
        "These are **fixed D=3 nearest-neighbor singlet RVB tensors**, contracted on the infinite plane with finite-environment CTMRG. No energy optimization was performed. Finite chi is a contraction approximation; this is not yet a completed infinite-chi extrapolation or a full reproduction of a paper's exponent.",
        "",
        f"Upstream: [{UPSTREAM_URL}]({UPSTREAM_URL}), commit `{UPSTREAM_COMMIT}`. The source checkout is unmodified. CPU execution; {version_text}; exact commands are in `manifest.json`.",
        "",
        "The adapter in `scripts/run_upstream.py` adds connected dimer correlations via upstream contraction functions. D means the physical bond operator S_i dot S_j. It subtracts the product of the two onepoint functions evaluated in the same finite strip as the four-spin correlator. This avoids a spurious offset from mixing different finite-chi contractions.",
        "",
        "These scalar-bond correlations are not the same observable as the z tensor z component-bond correlations plotted in Fig. 11 of Schuch et al., arXiv:1203.4816v3. That paper also studies an interpolation to orthogonal dimers; this report contains only the fixed singlet-RVB endpoint. Agreement of the qualitative decay is not a numerical reproduction of that figure.",
        "",
        "## Contraction convergence",
        "",
        "Square uses the Frobenius change of the nearest-neighbor reduced density matrix. Kagome uses the criterion recorded in the table: `Partial_energy` monitors one triangle energy, while `SingularValue` monitors all four corner spectra. Passing one of these criteria establishes iteration convergence at that chi, not convergence with chi.",
        "",
        "| Lattice | chi | E/site (J=1) | Last step | Criterion | Residual | Tolerance | Passed |",
        "|---|---:|---:|---:|---|---:|---:|---|",
    ]
    for result in results:
        lines.append(f"| {result['lattice']} | {result['chi']} | {result['energy_per_site']:.12f} | {result.get('ctm_last_step', '?')} | {result.get('ctm_convergence_criterion', '?')} | {result.get('ctm_convergence_residual', float('nan')):.3g} | {result.get('ctm_convergence_tolerance', float('nan')):.3g} | {result.get('ctm_converged', False)} |")
    lines += [
        "", "## Square lattice", "",
        "The dimer function is staggered. The table fits absolute connected correlations on **even separations only**, using an unweighted linear fit of log|C| against log r. No presumed literature exponent is imposed on the figure. Different windows and chi must agree before quoting an asymptotic exponent. The square spin correlation is short ranged; it must not be confused with the critical dimer correlation.",
        "",
        "| chi | p, r=4..12 | p, r=6..20 | p, r=10..30 | p, r=20..40 | Leading transfer length |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    summary = {"square": [], "kagome": []}
    previous = None
    for result in results:
        if result["lattice"] != "square":
            continue
        rows = result["correlations"]["DD"]
        r = np.array([row["distance"] for row in rows])
        c = np.array([row["dd_connected"] for row in rows])
        exponents = []
        for lower, upper in [(4, 12), (6, 20), (10, 30), (20, 40)]:
            mask = (r >= lower) & (r <= upper) & (r % 2 == 0)
            exponents.append(float(-np.polyfit(np.log(r[mask]), np.log(np.abs(c[mask])), 1)[0]))
        xi = next(iter(result["leading_transfer_length"].values()))
        lines.append("| " + " | ".join([str(result["chi"]), *[f"{p:.4f}" for p in exponents], f"{xi:.4f}"]) + " |")
        entry = {"chi": result["chi"], "power_law_windows": [[4, 12], [6, 20], [10, 30], [20, 40]], "exponents": exponents, "leading_transfer_length": xi}
        if previous is not None:
            previous_r = {row["distance"]: row["dd_connected"] for row in previous["correlations"]["DD"]}
            entry["relative_curve_change_from_chi"] = previous["chi"]
            entry["relative_curve_changes"] = {str(distance): float(abs(c[r == distance][0] - previous_r[distance]) / abs(c[r == distance][0])) for distance in [4, 8, 12, 20, 30, 40] if distance in r and distance in previous_r}
        summary["square"].append(entry)
        previous = result
    if len(summary["square"]) > 1:
        lines += ["", "Relative change of the even-distance dimer curve on increasing chi, using the newer curve as denominator:", "",
                  "| chi change | r=4 | r=12 | r=20 | r=40 |", "|---|---:|---:|---:|---:|"]
        for entry in summary["square"][1:]:
            changes = entry["relative_curve_changes"]
            lines.append("| " + " | ".join([f"{entry['relative_curve_change_from_chi']} to {entry['chi']}",
                         *[f"{100 * changes[str(distance)]:.2f}%" if str(distance) in changes else "n/a" for distance in [4, 12, 20, 40]]]) + " |")
    lines += [
        "", "The increasing leading transfer length shows that finite chi cuts off the critical state. It is an effective CTMRG scale, not a converged physical correlation length. Long-distance tails at fixed chi are exponential even for a critical state. The fitted exponents above are descriptive, with no statistical error bars; systematic chi and window effects dominate.",
        "", "## Kagome lattice", "",
        "D is the (0,1) bond within the blocked down triangle. Distances are kagome **Bravais unit cells**, equal to two nearest-neighbor spacings. Directions x and y refer to the two blocked lattice directions for that same oriented bond. The 2..8 entries below are descriptive short-window slopes, not established asymptotic correlation lengths. A change with fit window can reflect several transfer modes with different amplitudes.",
        "",
        "| chi | Short-window xi_x, r=2..8 | Short-window xi_y, r=2..8 | Generic transfer length, x |",
        "|---:|---:|---:|---:|",
    ]
    for result in results:
        if result["lattice"] != "kagome":
            continue
        fits, window_diagnostics = [], {}
        for label, direction in [("x", "(1, 0)"), ("y", "(0, 1)")]:
            rows = result["correlations"][f"DD_intra[(0, 0),{direction}]"]
            diagnostics = [exponential_window(rows, lower, upper)
                           for lower, upper in [(2, 8), (4, 10), (6, 12), (8, 14)]]
            fits.append(diagnostics[0]["xi"])
            window_diagnostics[label] = diagnostics
        xi = result["leading_transfer_length"]["spectrum(T)[(0, 0),(1, 0)]"]
        formatted = [f"{value:.6f}" if value is not None else "unresolved" for value in fits]
        lines.append(f"| {result['chi']} | {formatted[0]} | {formatted[1]} | {xi:.6f} |")
        summary["kagome"].append({"chi": result["chi"], "dimer_xi_x": fits[0], "dimer_xi_y": fits[1],
                                  "dimer_xi_window": [2, 8],
                                  "dimer_xi_interpretation": "Descriptive short-window slopes, not asymptotic correlation lengths.",
                                  "tail_window_diagnostics": window_diagnostics,
                                  "leading_transfer_length": xi,
                                  "transfer_length_caveat": "Generic transfer spectrum without observable weights; not a dimer length."})
    if summary["kagome"]:
        lines += ["", "### Tail-window diagnostics", "",
                  "Each fit uses every distance in the stated window. A window is marked unresolved if it is incomplete or contains an absolute correlation below 1e-13. Log RMS measures fit curvature; it is not a statistical uncertainty or a CTMRG error bound.", "",
                  "| chi | Direction | Window | Fitted xi (cells) | Log RMS residual |", "|---:|---|---|---:|---:|"]
        for entry in summary["kagome"]:
            for direction, diagnostics in entry["tail_window_diagnostics"].items():
                for diagnostic in diagnostics:
                    lower, upper = diagnostic["window"]
                    value = f"{diagnostic['xi']:.6f}" if diagnostic["xi"] is not None else "unresolved"
                    residual = f"{diagnostic['log_rms_residual']:.6g}" if "log_rms_residual" in diagnostic else "n/a"
                    lines.append(f"| {entry['chi']} | {direction} | {lower}..{upper} | {value} | {residual} |")
        largest = max(summary["kagome"], key=lambda entry: entry["chi"])
        for direction, diagnostics in largest["tail_window_diagnostics"].items():
            first, previous_tail, tail = diagnostics[0], diagnostics[-2], diagnostics[-1]
            if any(item["xi"] is None for item in [first, previous_tail, tail]):
                continue
            change = abs(tail["xi"] / previous_tail["xi"] - 1)
            assessment = ("The two later windows agree within 0.5%, giving a resolved finite-chi tail scale."
                          if change < .005 else
                          "The later-window slopes are still drifting; an asymptotic tail is not established.")
            lines += ["", f"At the largest loaded chi={largest['chi']}, direction {direction}: "
                      f"xi changes from {first['xi']:.6f} for r=2..8 to {previous_tail['xi']:.6f} for r=6..12 "
                      f"and {tail['xi']:.6f} cells for r=8..14. {assessment}"]
        lines += ["", "The short-window slopes must not be identified with asymptotic lengths. These data do not establish distinct asymptotic x/y lengths. Larger chi and observable-resolved spectra remain necessary for an extrapolated bulk length."]
    lines += [
        "", "The generic transfer gap can belong to a different channel from the dimer operator and is deliberately reported separately. The ten eigenvalues saved in these runs do not reach the ratio associated with the resolved scalar-bond tail, and their dimer weights have not been evaluated. They also do not implement the topological-sector exclusion and circumference extrapolation in Fig. 12 of arXiv:1203.4816v3.",
        "", "Kagome values below 1e-13 are omitted from the logarithmic figure. This conservative numerical cutoff guards against cancellation of O(0.04) terms in the connected subtraction; it is not a CTMRG systematic-error estimate. Upstream's chi=18 kagome test reference is E/site=-0.3931221584692804; this is an implementation regression target at finite chi, not an exact thermodynamic energy.",
        "", "## Repeat the calculation", "",
        "From the repository, after running `bash scripts/setup-upstream.sh`:",
        "", "```bash", "/workspace/peps-environment/.venv/bin/python scripts/run_upstream.py",
        "# More demanding square calculation, retaining earlier result files:",
        "/workspace/peps-environment/.venv/bin/python scripts/run_upstream.py --square-chi 128 --kagome-chi --distance 96 --max-iter 3000",
        "# Stricter kagome convergence:",
        "/workspace/peps-environment/.venv/bin/python scripts/run_upstream.py --square-chi --kagome-chi 18 32 --ctm-tol 1e-10",
        "```", "",
        "The default source lives outside the project at `/workspace/peps-environment/upstream/peps-torch`. A missing source directory is cloned at the pinned revision; an existing source at another revision or with tracked changes is rejected. `--source` and `--output` can select other locations. Optional missing YASTN/ArrayFire/cuQuantum notices do not affect these dense CPU contractions.",
        "", "![Correlation curves](correlations.png)", "",
    ]
    (destination / "report.md").write_text("\n".join(lines))
    (destination / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/workspace/peps-environment/upstream/peps-torch"))
    parser.add_argument("--output", type=Path, default=ROOT / "results/upstream")
    parser.add_argument("--square-chi", type=int, nargs="*", default=[16, 32, 64])
    parser.add_argument("--kagome-chi", type=int, nargs="*", default=[18, 32])
    parser.add_argument("--distance", type=int, default=60)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-iter", type=int, default=3000)
    parser.add_argument("--ctm-tol", type=float, default=1e-8)
    parser.add_argument("--kagome-ctm-check", choices=["Partial_energy", "Energy", "SingularValue"], default="SingularValue")
    parser.add_argument("--worker", choices=["square", "kagome"], help=argparse.SUPPRESS)
    args, remaining = parser.parse_known_args()
    if args.worker:
        worker(args.worker, args.source.resolve(), remaining)
        return
    if remaining:
        parser.error(f"unknown arguments: {remaining}")
    if args.threads < 1 or args.max_iter < 1 or args.ctm_tol <= 0:
        parser.error("threads, max-iter and ctm-tol must be positive")
    if any(chi < 2 for chi in args.square_chi + args.kagome_chi):
        parser.error("chi must be at least 2")
    if args.square_chi and args.distance < 40:
        parser.error("the square fit report requires --distance >= 40")
    if args.kagome_chi and args.distance < 8:
        parser.error("the kagome fit report requires --distance >= 8")
    source = args.source.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        source.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", "--filter=blob:none", UPSTREAM_URL, str(source)], check=True)
        subprocess.run(["git", "-C", str(source), "checkout", "--detach", UPSTREAM_COMMIT], check=True)
    actual = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if actual != UPSTREAM_COMMIT:
        raise SystemExit(f"Expected upstream {UPSTREAM_COMMIT}, found {actual}. Use a fresh --source directory.")
    changes = subprocess.check_output(["git", "-C", str(source), "diff", "--name-only", "HEAD"], text=True)
    if changes.strip():
        raise SystemExit("Upstream has tracked changes; use a clean --source directory for reproducible results.")
    import numpy
    import scipy
    import torch
    manifest = {"upstream_url": UPSTREAM_URL, "upstream_commit": actual,
                "versions": {"python": sys.version, "numpy": numpy.__version__, "scipy": scipy.__version__, "torch": torch.__version__},
                "runs": []}
    if (output / "manifest.json").exists():
        previous = json.loads((output / "manifest.json").read_text())
        if previous.get("upstream_commit") == actual:
            manifest["runs"] = previous.get("runs", [])
    results = []
    environment = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS=str(args.threads), MKL_NUM_THREADS=str(args.threads))
    for lattice, chis in [("square", args.square_chi), ("kagome", args.kagome_chi)]:
        for chi in chis:
            stem = f"{lattice}_chi{chi}"
            command = [sys.executable, "-u", str(Path(__file__).resolve()), "--worker", lattice,
                       "--source", str(source), "--chi", str(chi), "--corrf_r", str(args.distance if lattice == "square" else min(args.distance, 20)),
                       "--top_n", "10", "--CTMARGS_ctm_max_iter", str(args.max_iter), "--CTMARGS_ctm_conv_tol", str(args.ctm_tol), "--omp_cores", str(args.threads),
                       "--out_prefix", str(output / f"{stem}_internal")]
            command += ["--instate", "test-input/RVB_1x1.in", "--j2", "0.0"] if lattice == "square" else [
                "--ansatz", "IPESS", "--instate", "test-input/IPESS_KAGOME_D3_RVB.in", "--bond_dim", "3", "--GLOBALARGS_dtype", "complex128", "--CTM_check", args.kagome_ctm_check]
            print(f"Running {stem}...", flush=True)
            start = time.monotonic()
            with (output / f"{stem}.log").open("w") as log:
                subprocess.run(command, cwd=source, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            elapsed = time.monotonic() - start
            result = parse_log(output / f"{stem}.log", lattice, chi)
            result["elapsed_seconds"] = elapsed
            (output / f"{stem}.json").write_text(json.dumps(result, indent=2) + "\n")
            results.append(result)
            manifest["runs"] = [run for run in manifest["runs"] if run.get("name") != stem]
            manifest["runs"].append({"name": stem, "command": command, "cwd": str(source), "stdout_log": str(output / f"{stem}.log"), "elapsed_seconds": elapsed})
            (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            print(f"  energy/site={result['energy_per_site']:.12f}, {elapsed:.1f}s", flush=True)
    results = [json.loads(path.read_text()) for path in output.glob("*_chi*.json")]
    results.sort(key=lambda result: (result["lattice"], result["chi"]))
    plot_results(results, output)
    write_report(results, output)
    print(f"Saved logs, JSON data and plot to {output}")


if __name__ == "__main__":
    main()
