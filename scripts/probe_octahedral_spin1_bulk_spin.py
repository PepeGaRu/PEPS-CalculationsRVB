#!/usr/bin/env python3
"""Finite cubic 3D spin correlations with shared HOTRG norm projectors.

An exploratory finite-chi calculation, never a bulk correlation-length fit.
Physical operator insertions retain full16 support on their bonds. The first
y coarse-graining puts those observed y bonds inside each paired block.
All subsequent projectors come exclusively from the norm tensor network.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import quimb

from rvb.octahedral_spin1 import spin1_generators
from rvb.octahedral_spin1_bulk import build_network


def shared_coarse_grain(norm, impurities, direction, chi, cutoff=1e-12):
    """Contract one HOTRG step with identical norm-derived projectors.

    Original tensor ids and axis orders survive lazy projector insertion.
    Projectors may reindex the original legs, so replacement uses shape and
    tensor id, retaining the new indices in the projected network.
    """
    old_ids = tuple(norm.tensor_map)
    projected = norm.coarse_grain_hotrg(
        direction, max_bond=chi, cutoff=cutoff, canonize=False,
        lazy=True, equalize_norms=False, optimize="greedy")
    outputs = []
    for impurity in impurities:
        result = projected.copy()
        for tid in old_ids:
            source, target = impurity.tensor_map[tid], result.tensor_map[tid]
            if target.shape != source.shape:
                raise ValueError("Mismatched original tensor shape")
            target.modify(data=source.data)
        outputs.append(result)
    for network in (projected, *outputs):
        for site in network.site_tags:
            network.contract_tags_(site, optimize="greedy")
    for result in outputs:
        if set(result.tensor_map) != set(projected.tensor_map):
            raise ValueError("Tensor ids diverged during common projection")
        for tid, tensor in projected.tensor_map.items():
            if tensor.inds != result.tensor_map[tid].inds:
                raise ValueError("Tensor indices diverged during common projection")
    return projected, outputs


def probe(length, chi, components="xyz"):
    started = time.perf_counter()
    distances = list(range(1, length // 2 + 1))
    edges = [((x, 0, length // 2), 1) for x in range(length // 2 + 1)]
    generators = dict(zip("xyz", spin1_generators()))
    sz = generators["z"]
    norm = build_network((length,) * 3, edges)
    impurities = [
        build_network((length,) * 3, edges,
                      {edges[0]: generators[axis], edges[r]: generators[axis]})
        for r in distances for axis in components
    ]
    impurities.append(build_network((length,) * 3, edges, {edges[0]: sz}))
    steps = []
    # Observed y bonds go from y=0 to1 and are internal to the first pair.
    sequence = ["y", "x", "z"]
    while sequence:
        axis = sequence.pop(0)
        if getattr(norm, "L" + axis) <= 2:
            continue
        sequence.append(axis)
        norm, impurities = shared_coarse_grain(norm, impurities, axis, chi)
        step = dict(direction=axis, seconds=time.perf_counter() - started,
                    shape=[norm.Lx, norm.Ly, norm.Lz])
        steps.append(step)
        print(json.dumps(step), flush=True)
        if sum(getattr(norm, "L" + a) > 2 for a in "xyz") <= 1:
            break
    denominator = complex(norm.contract(optimize="greedy"))
    numerators = [complex(network.contract(optimize="greedy"))
                  for network in impurities]
    correlations = [value / denominator for value in numerators]
    return dict(
        geometry="periodic isotropic cubic 3D", L=length,
        physical_spins=3 * length ** 3, chi=chi,
        method="HOTRG with common norm-derived projectors",
        quimb_version=quimb.__version__, cutoff=1e-12, canonize=False,
        optimize="greedy", sequence=["y", "x", "z"],
        equalize_norms=False, observed_bond_dimension=16,
        norm=[denominator.real, denominator.imag],
        distance=distances,
        components=list(components),
        spin_correlations={
            "C" + axis + axis: [[correlations[i * len(components) + j].real,
                                 correlations[i * len(components) + j].imag]
                                for i in range(len(distances))]
            for j, axis in enumerate(components)},
        one_point_Sz=[correlations[-1].real, correlations[-1].imag],
        steps=steps, seconds=time.perf_counter() - started,
        correlation_length=None,
        warning="Unconverged finite-cube diagnostic; no bulk length inferred.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--chi", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--components", choices=("z", "xyz"), default="xyz")
    args = parser.parse_args()
    if args.L < 3 or args.chi < 1:
        parser.error("L must be >=3 and chi positive")
    result = probe(args.L, args.chi, args.components)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
