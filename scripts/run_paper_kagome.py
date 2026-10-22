#!/usr/bin/env python3
"""Paper-informed kagome RVB/interpolation contractions using pinned peps-torch.

The declared physical ZZ operator is (I_flag x sigma_z) on each spin. This
explicit convention is a paper-motivated comparison, not an asserted exact
identification of the paper's z=diag(1,-1,1) observable.
"""
from __future__ import annotations

import argparse
from functools import reduce
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_COMMIT = "239a2e717a567df0ece80a5bce5443a6f78d048b"


def physical_bond_operator(pair: tuple[int, int], observable: str, flagged: bool) -> np.ndarray:
    """Three-spin operator, local basis [flag,spin] when flagged=True."""
    if pair not in ((0, 1), (0, 2), (1, 2)):
        raise ValueError("pair must be one of (0,1), (0,2), (1,2)")
    z = np.diag([1., -1.])
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    paulis = [z] if observable == "zz" else [x, y, z]
    if observable not in ("zz", "ss"):
        raise ValueError("observable must be zz or ss")
    d = 4 if flagged else 2
    result = np.zeros((d ** 3, d ** 3), dtype=complex)
    for pauli in paulis:
        local = np.kron(np.eye(2), pauli) if flagged else pauli
        factors = [local if site in pair else np.eye(d) for site in range(3)]
        result += reduce(np.kron, factors)
    return np.real_if_close(result if observable == "zz" else result / 4)


def finite_window_fits(rows: list[dict]) -> list[dict]:
    """Descriptive exponential fits; no data below explicit numerical threshold."""
    fits = []
    for lower, upper in [(1, 4), (2, 5), (2, 8), (3, 8)]:
        selected = [row for row in rows if lower <= row["distance"] <= upper
                    and abs(row["connected_centered"]) > 1e-13]
        if len(selected) < 3:
            continue
        r = np.array([row["distance"] for row in selected])
        logc = np.log(np.abs([row["connected_centered"] for row in selected]))
        slope, intercept = np.polyfit(r, logc, 1)
        residual = logc - (slope * r + intercept)
        fits.append({"requested_window": [lower, upper], "used_distances": r.tolist(),
                     "xi_unit_cells": float(-1 / slope) if slope < 0 else None,
                     "log_slope": float(slope), "log_intercept": float(intercept),
                     "log_fit_rms": float(np.sqrt(np.mean(residual ** 2))),
                     "absolute_correlation_floor": 1e-13})
    return fits


def calculate(args) -> dict:
    import torch

    sys.path[:0] = [str(ROOT), str(args.source)]
    import config as cfg
    from ctm.generic import corrf, ctmrg, rdm
    from ctm.generic.env import ENV, init_env, ctmrg_conv_specC
    from ipeps.ipeps import IPEPS
    from rvb.peps import kagome_site_tensor

    torch.set_num_threads(args.threads)
    torch.manual_seed(0)
    cfg.global_args.dtype = "float64"
    cfg.global_args.torch_dtype = torch.float64
    cfg.global_args.device = "cpu"
    cfg.ctm_args.ctm_max_iter = args.max_iter
    cfg.ctm_args.ctm_conv_tol = args.ctm_tol
    cfg.ctm_args.ctm_force_dl = True
    cfg.ctm_args.projector_svd_method = "GESDD"
    cfg.ctm_args.projector_svd_reltol = args.projector_tol
    flagged = args.theta != 0 or args.embedded_qutrit
    if flagged:
        from rvb.interpolation import kagome_interpolating_tensor
        local = kagome_interpolating_tensor(args.theta)
    else:
        local = kagome_site_tensor()
    # Our [p,left,right,up,down] -> upstream [p,up,left,down,right].
    a = torch.from_numpy(local.transpose(0, 3, 1, 4, 2).copy()).to(torch.float64)
    a = a / torch.linalg.norm(a)
    state = IPEPS({(0, 0): a}, vertexToSite=lambda coord: (0, 0), lX=1, lY=1)
    env = ENV(args.chi, state)
    init_env(state, env)

    def convergence(state, env, history, ctm_args):
        stop, history = ctmrg_conv_specC(state, env, history, p="inf", ctm_args=ctm_args)
        n = len(history["conv_crit"])
        print(f"CTMRG iteration={n} corner_spectrum_residual={history['conv_crit'][-1]:.6g}", flush=True)
        return stop, history

    start = time.monotonic()
    env, history, ctm_time, check_time = ctmrg.run(state, env, conv_check=convergence)
    density = rdm.rdm1x1((0, 0), state, env)
    identity = torch.eye(a.shape[0], dtype=a.dtype)
    results = {}
    means = {}
    operator_types = ["zz", "ss"] if args.theta == 0 else ["zz"]
    if args.embedded_qutrit:
        operator_types += ["qutrit_first_first", "qutrit_first_second", "qutrit_second_first", "qutrit_second_second"]
    for observable in operator_types:
        for pair in [(0, 1), (0, 2), (1, 2)]:
            if observable.startswith("qutrit_"):
                from rvb.interpolation import local_qutrit_z
                _, first, second = observable.split("_")
                factors = [np.eye(4) for _ in range(3)]
                factors[pair[0]] = local_qutrit_z(first)
                factors[pair[1]] = local_qutrit_z(second)
                matrix = reduce(np.kron, factors)
            else:
                matrix = physical_bond_operator(pair, observable, flagged)
            op = torch.from_numpy(matrix).to(a.dtype)
            bulk_mean = torch.einsum("ij,ji->", density, op).item()
            name = f"{observable}_{pair[0]}{pair[1]}"
            means[name] = bulk_mean
            right_pairs = [(0, 1), (0, 2), (1, 2)] if args.cross_bonds and observable in ("zz", "ss") else [pair]
            for right_pair in right_pairs:
                op_right = torch.from_numpy(physical_bond_operator(right_pair, observable, flagged)).to(a.dtype) if right_pair != pair else op
                right_suffix = f"_to{right_pair[0]}{right_pair[1]}" if right_pair != pair else ""
                for axis, direction in [("x", (1, 0)), ("y", (0, 1))]:
                    def contract(first, second):
                        return corrf.corrf_1sO1sO((0, 0), direction, state, env,
                                                 first, lambda r: second, args.distance - 1)
                    raw = contract(op, op_right)
                    left = contract(op, identity)
                    right = contract(identity, op_right)
                    # Center each insertion near its long-strip onepoint. The remaining
                    # subtraction is then of tiny quantities, avoiding O(1)-O(1) tails.
                    mu_left = left[-1].item()
                    mu_right = right[-1].item()
                    centered = contract(op - mu_left * identity, op_right - mu_right * identity)
                    covariance = centered - (left - mu_left) * (right - mu_right)
                    rows = [{"distance": r + 1, "raw": raw[r].item(),
                             "left_mean": left[r].item(), "right_mean": right[r].item(),
                             "connected_subtraction": (raw[r] - left[r] * right[r]).item(),
                             "centered_raw": centered[r].item(),
                             "connected_centered": covariance[r].item(),
                             "subtraction_disagreement": (covariance[r] - (raw[r] - left[r] * right[r])).item()}
                            for r in range(args.distance)]
                    key = f"{name}{right_suffix}_{axis}"
                    results[key] = {"observable": observable, "pair": list(pair), "right_pair": list(right_pair), "direction": list(direction),
                                    "bulk_rdm_mean": bulk_mean, "centering_left": mu_left,
                                    "centering_right": mu_right, "rows": rows,
                                    "exponential_fits": finite_window_fits(rows)}
                    print(f"OBS {key} mean={bulk_mean:.12g} connected="
                          + ",".join(f"{row['connected_centered']:.6g}" for row in rows), flush=True)
    final = float(history["conv_crit"][-1])
    result = {"lattice": "kagome", "theta": args.theta, "chi": args.chi,
              "max_distance": args.distance, "physical_dimension": a.shape[0],
              "virtual_bond_dimension": 3, "observables": results, "bulk_means": means,
              "ctm_converged": final < args.ctm_tol, "ctm_convergence_residual": final,
              "ctm_convergence_tolerance": args.ctm_tol,
              "ctm_iterations": len(history["conv_crit"]),
              "ctm_convergence_history": [value if np.isfinite(value) else None for value in history["conv_crit"]],
              "ctm_convergence_definition": "maximum Euclidean change of normalized singular spectra across all four corners",
              "projector_svd_relative_cutoff": args.projector_tol,
              "final_corner_spectra": {str(key): value.tolist() for key, value in env.get_spectra().items()},
              "ctm_algorithm_parameters": vars(cfg.ctm_args),
              "ctm_seconds": ctm_time, "convergence_check_seconds": check_time,
              "total_seconds": time.monotonic() - start,
              "geometry": {"blocked_spins": ["a(x,y)", "b(x,y)", "c(x,y)"],
                           "up_triangle": ["a(x,y)", "b(x,y)", "c(x,y)"],
                           "down_triangle": ["a(x,y)", "b(x-1,y)", "c(x,y-1)"],
                           "distance_unit": "Bravais unit cells = 2 nearest-neighbor spacings",
                           "directions": {"x": [1, 0], "y": [0, 1]},
                           "singlet_orientation": "cyclic vertex order of each triangle"},
              "operator_convention": "ZZ=(I_flag tensor sigma_z)_i (I_flag tensor sigma_z)_j, identity on third spin; SS=S_i dot S_j",
              "paper_comparison_caveat": "Paper-motivated comparison using a declared physical ZZ convention; paper Fig.11 z=diag(1,-1,1) physical embedding is not asserted equivalent.",
              "connected_definition": "<O(0)O(r)> - <O(0)><O(r)> with both means from the same finite CTM strip; evaluated again using centered insertions",
              "theta_zero_compression": not flagged,
              "embedded_qutrit_convention": "When requested, first/second denotes J^dagger [z tensor I / I tensor z] J with z=diag(1,-1,1), and J maps |+>flag|s>spin to |s2>, |->flag|s>spin to |2s>; possible embeddings scanned without asserting paper identification.",
              "upstream_commit": UPSTREAM_COMMIT, "tensor_sha256": hashlib.sha256(local.tobytes()).hexdigest(),
              "versions": {"python": sys.version, "numpy": np.__version__, "torch": torch.__version__},
              "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
              "command": [sys.executable, *sys.argv],
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if "ss_01" in means:
        result["mean_intracell_bond_energy"] = sum(means[f"ss_{i}{j}"] for i, j in [(0, 1), (0, 2), (1, 2)]) / 3
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/workspace/peps-environment/upstream/peps-torch"))
    parser.add_argument("--output", type=Path, default=ROOT / "results/paper-recheck/ctmrg")
    parser.add_argument("--theta", type=float, default=0.)
    parser.add_argument("--chi", type=int, default=16)
    parser.add_argument("--distance", type=int, default=12)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--ctm-tol", type=float, default=1e-10)
    parser.add_argument("--projector-tol", type=float, default=1e-8, help="Relative cutoff on CTM projector singular values; lower values need a fresh convergence check")
    parser.add_argument("--embedded-qutrit", action="store_true", help="Also scan four single-qutrit-leg choices; uses uncompressed p=64 even at theta=0")
    parser.add_argument("--cross-bonds", action="store_true", help="Also correlate different intracell bond orientations")
    args = parser.parse_args()
    if not 0 <= args.theta <= 1 or args.chi < 2 or args.distance < 1 or args.threads < 1:
        parser.error("require 0 <= theta <= 1, chi >= 2, distance >= 1, threads >= 1")
    args.source = args.source.resolve()
    actual = subprocess.check_output(["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True).strip()
    if actual != UPSTREAM_COMMIT:
        raise SystemExit(f"Expected upstream {UPSTREAM_COMMIT}, found {actual}")
    changes = subprocess.check_output(["git", "-C", str(args.source), "diff", "--name-only", "HEAD"], text=True)
    if changes.strip():
        raise SystemExit("Upstream checkout has tracked modifications")
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(args.threads))
    import torch
    with torch.no_grad():
        result = calculate(args)
    suffix = "_qutrit" if args.embedded_qutrit else ""
    suffix += "_cross" if args.cross_bonds else ""
    path = args.output / f"theta{args.theta:g}_chi{args.chi}{suffix}.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Saved {path}; CTMRG converged={result['ctm_converged']}; total_seconds={result['total_seconds']:.3f}")


if __name__ == "__main__":
    main()
