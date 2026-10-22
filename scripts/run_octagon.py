#!/usr/bin/env python3
"""Contract the user's fully covered square-octagon triangle RVB in the plane.

The physical space has eight spin-1/2 degrees of freedom per square cell.
Only the exact support of the local tensor is compressed. No dimer covering
orthogonality is assumed. The two spins on an unused dummy link form a singlet.
"""
from __future__ import annotations

import argparse
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


def factorized_pair_apply(corrf, coord, direction, state, env, edge, terms):
    """Apply sum coefficient * A(cell) B(next cell) without a dense pair tensor."""
    next_coord = tuple(c + d for c, d in zip(coord, direction))
    result = None
    for coefficient, first, second in terms:
        value = corrf.apply_TM_1sO(coord, direction, state, env, edge, op=first)
        value = corrf.apply_TM_1sO(next_coord, direction, state, env, value, op=second)
        result = coefficient * value if result is None else result + coefficient * value
    return result


def factorized_pair_correlations(corrf, coord, direction, state, env,
                                 first_terms, second_terms, max_distance):
    """Two-cell bond correlations at center separation r=2,...,max_distance.

    This duplicates the geometry of upstream corrf_2sOH2sOH_E1, but exploits
    the three-term Heisenberg operator instead of SVD of a large dense matrix.
    Operators cannot overlap; the first pair occupies cells 0 and 1.
    """
    import torch

    reverse = tuple(-d for d in direction)
    edge = corrf.get_edge(coord, reverse, state, env)
    identity_terms = [(1., None, None)]
    e1 = factorized_pair_apply(corrf, coord, direction, state, env, edge, first_terms)
    e0 = factorized_pair_apply(corrf, coord, direction, state, env, edge, identity_terms)
    current = tuple(c + 2 * d for c, d in zip(coord, direction))
    values = []
    for distance in range(2, max_distance + 1):
        e12 = factorized_pair_apply(corrf, current, direction, state, env, e1, second_terms)
        e0 = corrf.apply_TM_1sO(current, direction, state, env, e0)
        e1 = corrf.apply_TM_1sO(current, direction, state, env, e1)
        current = tuple(c + d for c, d in zip(current, direction))
        numerator = corrf.apply_edge(current, direction, state, env, e12)
        denominator = corrf.apply_TM_1sO(current, direction, state, env, e0)
        denominator = corrf.apply_edge(current, direction, state, env, denominator)
        values.append(numerator / denominator)
        scale = e0.abs().max()
        e0, e1 = e0 / scale, e1 / scale
    return torch.stack(values)


def connected_rows(raw, left, right, centered, first_distance):
    mu_left, mu_right = left[-1].item(), right[-1].item()
    covariance = centered - (left - mu_left) * (right - mu_right)
    return [{"distance": index + first_distance, "raw": raw[index].item(),
             "left_mean": left[index].item(), "right_mean": right[index].item(),
             "connected_subtraction": (raw[index] - left[index] * right[index]).item(),
             "centered_raw": centered[index].item(),
             "connected_centered": covariance[index].item(),
             "subtraction_disagreement": (covariance[index] - raw[index] + left[index] * right[index]).item()}
            for index in range(len(raw))]


def calculate(args):
    import torch

    sys.path[:0] = [str(ROOT), str(args.source)]
    import config as cfg
    from ctm.generic import corrf, ctmrg, rdm, ctm_projectors
    from ctm.generic.env import ENV, init_env, ctmrg_conv_specC
    from ipeps.ipeps import IPEPS
    from rvb.octagon_peps import (octagon_compressed_tensor, octagon_local_tensor,
                                  projected_bond_operator, projected_spin_operator,
                                  physical_bond_operator, physical_spin_operator)
    from rvb.diagnostics import diagnose_decay

    torch.set_num_threads(args.threads)
    torch.manual_seed(args.seed)
    cfg.global_args.dtype = "float64"
    cfg.global_args.torch_dtype = torch.float64
    cfg.global_args.device = "cpu"
    cfg.ctm_args.ctm_max_iter = args.max_iter
    cfg.ctm_args.ctm_conv_tol = args.ctm_tol
    cfg.ctm_args.ctm_force_dl = True
    cfg.ctm_args.projector_svd_method = "GESDD"
    cfg.ctm_args.projector_svd_reltol = args.projector_tol
    cfg.ctm_args.projector_multiplet_abstol = args.multiplet_abstol

    compression = octagon_compressed_tensor()
    local = octagon_local_tensor() if args.raw_physical else compression.tensor
    a = torch.from_numpy(local.transpose(0, 3, 1, 4, 2).copy()).to(torch.float64)
    # A power-of-two scale preserves the exact integer tensor's zeros in its
    # double layer during the optional uncompressed precision audit.
    power_two_scale = args.raw_physical or args.power_two_scale
    a = a / 16 if power_two_scale else a / torch.linalg.norm(a)
    state = IPEPS({(0, 0): a}, vertexToSite=lambda coord: (0, 0), lX=1, lY=1)
    env = ENV(args.chi, state)
    init_env(state, env)
    projector_ranks = []
    original_svd = ctm_projectors.truncated_svd_gesdd

    def observed_svd(*positional, **keywords):
        result = original_svd(*positional, **keywords)
        singular = result[1]
        projector_ranks.append(int((singular / singular[0] > args.projector_tol).sum().item()))
        return result

    ctm_projectors.truncated_svd_gesdd = observed_svd

    def convergence(state, env, history, ctm_args):
        _, history = ctmrg_conv_specC(state, env, history, p="inf", ctm_args=ctm_args)
        residuals = history["conv_crit"]
        converged = len(residuals) >= 4 and all(value < args.ctm_tol for value in residuals[-3:])
        print(f"CTMRG iteration={len(residuals)} corner_residual={residuals[-1]:.6g} "
              f"projector_ranks={projector_ranks[-4:]}", flush=True)
        return converged or len(residuals) >= args.max_iter, history

    start = time.monotonic()
    try:
        env, history, ctm_time, check_time = ctmrg.run(state, env, conv_check=convergence)
    finally:
        ctm_projectors.truncated_svd_gesdd = original_svd
    density = rdm.rdm1x1((0, 0), state, env)
    identity = torch.eye(a.shape[0], dtype=a.dtype)
    results, means = {}, {}

    def diagnostic(rows):
        r = [row["distance"] for row in rows]
        c = [row["connected_centered"] for row in rows]
        windows = [(2, 6), (3, 8), (5, 10), (8, 14), (12, args.distance)]
        windows = [window for window in windows if window[0] < window[1] <= args.distance]
        return diagnose_decay(r, c, windows=windows or [(1, args.distance)], floor=args.floor,
                              branch=np.asarray(r) % 2)

    def tensor(matrix):
        return torch.from_numpy(np.real_if_close(matrix).copy()).to(a.dtype)

    def bond_operator(i, j):
        return physical_bond_operator(i, j) if args.raw_physical else projected_bond_operator(i, j, compression=compression)

    def spin_operator(site, component):
        return physical_spin_operator(site, component) if args.raw_physical else projected_spin_operator(site, component, compression=compression)

    operators = {"triangle_ring_ring": bond_operator(0, 3),
                 "triangle_ring_dummy": bond_operator(0, 4),
                 "triangle_other_ring_dummy": bond_operator(3, 4),
                 "ring_spin_z": spin_operator(0, "z"),
                 "dummy_spin_z": spin_operator(4, "z")}
    for name, matrix in operators.items():
        op = tensor(matrix)
        means[name] = torch.einsum("ij,ji->", density, op).item()
        for axis, direction in [("x", (1, 0)), ("y", (0, -1))]:
            def contract(first, second):
                return corrf.corrf_1sO1sO((0, 0), direction, state, env, first,
                                         lambda r: second, args.distance - 1)
            raw, left, right = contract(op, op), contract(op, identity), contract(identity, op)
            centered = contract(op - left[-1] * identity, op - right[-1] * identity)
            rows = connected_rows(raw, left, right, centered, 1)
            results[f"{name}_{axis}"] = {"operator": name, "direction": list(direction),
                                          "bulk_rdm_mean": means[name], "rows": rows,
                                          "diagnostics": diagnostic(rows)}
            print(f"OBS {name}_{axis} mean={means[name]:.12g} C="
                  + ",".join(f"{row['connected_centered']:.6g}" for row in rows), flush=True)

    if not args.skip_dummy_bonds:
        for axis, direction, first_site, second_site in [("x", (1, 0), 4, 6), ("y", (0, -1), 5, 7)]:
            spin = {(site, component): tensor(spin_operator(site, component))
                    for site in (first_site, second_site) for component in ("z", "plus", "minus")}
            terms = [(1., spin[first_site, "z"], spin[second_site, "z"]),
                     (.5, spin[first_site, "plus"], spin[second_site, "minus"]),
                     (.5, spin[first_site, "minus"], spin[second_site, "plus"])]
            id_terms = [(1., None, None)]

            def contract_pairs(first, second):
                return factorized_pair_correlations(corrf, (0, 0), direction, state, env,
                                                    first, second, args.distance)
            raw = contract_pairs(terms, terms)
            left, right = contract_pairs(terms, id_terms), contract_pairs(id_terms, terms)
            first_centered = terms + [(-left[-1].item(), None, None)]
            second_centered = terms + [(-right[-1].item(), None, None)]
            centered = contract_pairs(first_centered, second_centered)
            rows = connected_rows(raw, left, right, centered, 2)
            name = f"dummy_bond_{axis}"
            results[name] = {"operator": "S_dummy_forward(cell) dot S_dummy_backward(next_cell)",
                             "direction": list(direction), "dummy_spin_indices": [first_site, second_site],
                             "rows": rows, "diagnostics": diagnostic(rows),
                             "minimum_distance_reason": "Each bond occupies two adjacent cells; r>=2 avoids overlap."}
            means[name] = .5 * (left[-1] + right[-1]).item()
            print(f"OBS {name} mean={means[name]:.12g} C="
                  + ",".join(f"{row['connected_centered']:.6g}" for row in rows), flush=True)

    final = float(history["conv_crit"][-1])
    return {"lattice": "square_octagon_triangle_decoration", "chi": args.chi,
            "state_definition": "Equal amplitudes of all perfect singlet coverings with cyclic triangle singlets and positive E/N-oriented dummy singlets; both dummy halves locally covered or paired across dummy link.",
            "physical_spins_per_cell": 8, "physical_dimension_before_compression": 256,
            "physical_dimension_after_exact_compression": int(a.shape[0]),
            "physical_compression_used": not args.raw_physical,
            "physical_compression_error": 0. if args.raw_physical else float(compression.compression_error),
            "local_tensor_scaling": "power of two: A/16" if power_two_scale else "unit Frobenius norm",
            "physical_singular_values": compression.singular_values.tolist(),
            "virtual_bond_dimension": 3, "max_distance": args.distance,
            "geometry": {"spin_order": ["r0", "r1", "r2", "r3", "dE", "dN", "dW", "dS"],
                         "triangle_at_E": [0, 3, 4],
                         "distance_unit": "square Bravais unit cells", "direction_convention": "physical x/east uses upstream (1,0); physical y/north uses upstream (0,-1), because upstream positive y points down"},
            "observables": results, "bulk_means": means,
            "connected_definition": "<O0 Or>-<O0><Or>, with means from the same CTM strip; evaluated using centered operator insertions as well as direct subtraction",
            "correlation_floor": args.floor,
            "ctm_converged": len(history["conv_crit"]) >= 4 and all(v < args.ctm_tol for v in history["conv_crit"][-3:]),
            "ctm_convergence_residual": final, "ctm_convergence_tolerance": args.ctm_tol,
            "ctm_iterations": len(history["conv_crit"]),
            "ctm_convergence_history": [v if np.isfinite(v) else None for v in history["conv_crit"]],
            "ctm_convergence_definition": "maximum L2 change of normalized singular spectra over all four corners, below tolerance for three consecutive iterations",
            "projector_svd_relative_cutoff": args.projector_tol,
            "projector_effective_ranks_last_sweep": projector_ranks[-4:],
            "projector_effective_ranks_min_max": [min(projector_ranks), max(projector_ranks)],
            "final_corner_spectra": {str(k): v.tolist() for k, v in env.get_spectra().items()},
            "ctm_algorithm_parameters": vars(cfg.ctm_args),
            "ctm_seconds": ctm_time, "convergence_check_seconds": check_time,
            "total_seconds": time.monotonic() - start,
            "upstream_commit": UPSTREAM_COMMIT,
            "tensor_sha256": hashlib.sha256(local.tobytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "versions": {"python": sys.version, "numpy": np.__version__, "torch": torch.__version__},
            "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "command": [sys.executable, *sys.argv]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/workspace/peps-environment/upstream/peps-torch"))
    parser.add_argument("--output", type=Path, default=ROOT / "results/octagon-square/ctmrg")
    parser.add_argument("--chi", type=int, default=16)
    parser.add_argument("--distance", type=int, default=20)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--ctm-tol", type=float, default=1e-10)
    parser.add_argument("--projector-tol", type=float, default=1e-8)
    parser.add_argument("--multiplet-abstol", type=float, default=1e-14,
                        help="Separate upstream absolute singular-value cutoff used in multiplet preservation")
    parser.add_argument("--floor", type=float, default=1e-13)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--skip-dummy-bonds", action="store_true")
    parser.add_argument("--raw-physical", action="store_true", help="Precision audit using the integer p=256 tensor, scaled by 1/16, without physical SVD compression")
    parser.add_argument("--power-two-scale", action="store_true", help="Use A/16 also for the compressed tensor, to isolate scaling from compression in the precision audit")
    args = parser.parse_args()
    if args.chi < 2 or args.distance < 4 or args.threads < 1 or not 0 < args.projector_tol < 1:
        parser.error("require chi>=2, distance>=4, threads>=1, 0<projector-tol<1")
    actual = subprocess.check_output(["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True).strip()
    changes = subprocess.check_output(["git", "-C", str(args.source), "diff", "--name-only", "HEAD"], text=True)
    if actual != UPSTREAM_COMMIT or changes.strip():
        raise SystemExit("Pinned peps-torch source is not the expected unmodified commit")
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(args.threads))
    import torch
    with torch.no_grad():
        result = calculate(args)
    suffix = "_raw" if args.raw_physical else ""
    suffix += "_scale16" if args.power_two_scale and not args.raw_physical else ""
    suffix += f"_mult{args.multiplet_abstol:g}" if args.multiplet_abstol != 1e-14 else ""
    path = args.output / f"chi{args.chi}_cutoff{args.projector_tol:g}{suffix}.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Saved {path}; converged={result['ctm_converged']}; seconds={result['total_seconds']:.3f}")


if __name__ == "__main__":
    main()
