#!/usr/bin/env python3
"""Exploratory finite-box spin-one overlap sampling, with exact pair contractions.

This samples |<bra|ket>| times a bounded bias that enhances unequal trimer
location fields. Reweight BOTH observables and normalization by 1/bias.
A local spin estimator is <bra|S_i.S_j|ket>/|<bra|ket>|, including its sign.
These runs do not establish equilibration or a thermodynamic correlation length.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rvb.octahedral import OctahedralCoverings
from rvb.octahedral_spin1_overlap import contract_pair, line_flip_mask


def run(length, steps, seed, initial, bias, output):
    if length < 2 or steps <= 0:
        raise ValueError("L must be at least 2 and steps must be positive")
    if not math.isfinite(bias) or bias <= 0:
        raise ValueError("Bias must be positive and finite")
    rng = random.Random(seed)
    model = OctahedralCoverings(length, length, length)
    ket, bra = model.sample(2, seed=seed)
    if initial == "diagonal":
        bra = ket
    elif initial == "trimer":
        ket = bra = model.all_minus_trimer_reference
    elif initial == "dimer":
        if length % 2:
            raise ValueError("The all-dimer initial state requires even L")
        ket = sum((1 << site) for site, (*xyz, axis) in enumerate(model.site_labels)
                  if xyz[axis] % 2)
        bra = ket
    lines = sorted({line_flip_mask(model, axis, cell)
                    for axis in range(3) for cell in range(model.ncells)})
    cubes = [model.cube_flip_mask(cell) for cell in range(model.ncells)]

    def mismatch(a, b):
        delta = a ^ b
        return sum(((delta >> ports[0]) ^ (delta >> ports[1])) & 1
                   for ports in model.ports)

    current = contract_pair(model, ket, bra)
    if not current.numerator:
        raise ValueError("Initial pair has zero overlap; choose --initial diagonal")
    current_bias = bias if mismatch(ket, bra) else 1.
    moves = defaultdict(lambda: [0, 0])
    records, negative_proposals = [], 0
    started = time.perf_counter()
    for step in range(steps):
        u = rng.random()
        kind = ("single_line" if u < .4 else "joint_line" if u < .6
                else "single_cube" if u < .8 else "joint_cube")
        mask = rng.choice(lines if "line" in kind else cubes)
        proposed_ket, proposed_bra = ket, bra
        if kind.startswith("joint"):
            proposed_ket ^= mask
            proposed_bra ^= mask
        elif rng.getrandbits(1):
            proposed_ket ^= mask
        else:
            proposed_bra ^= mask
        proposed_bias = bias if mismatch(proposed_ket, proposed_bra) else 1.
        proposed = contract_pair(model, proposed_ket, proposed_bra)
        negative_proposals += proposed.sign < 0
        moves[kind][0] += 1
        log_ratio = (proposed.log_abs_value - current.log_abs_value
                     + math.log(proposed_bias/current_bias))
        if proposed.numerator and math.log(rng.random()) < log_ratio:
            ket, bra, current = proposed_ket, proposed_bra, proposed
            current_bias = proposed_bias
            moves[kind][1] += 1
        if (step + 1) % 25 == 0:
            trimers_ket = {o.cell for o in model.covering(ket) if o.kind == "trimer"}
            trimers_bra = {o.cell for o in model.covering(bra) if o.kind == "trimer"}
            entry = dict(step=step+1, weight=1/current_bias,
                         log_weight=-math.log(current_bias), sign=current.sign,
                         log_overlap=current.log_abs_value,
                         arrow_difference=(ket ^ bra).bit_count(),
                         trimer_density=(len(trimers_ket)+len(trimers_bra))/(2*model.ncells),
                         trimer_difference=len(trimers_ket ^ trimers_bra),
                         a=hex(ket), b=hex(bra))
            # Covering densities above are mixing diagnostics, not observables.
            cell, axis = rng.randrange(model.ncells), rng.randrange(3)
            first = 3*cell + axis
            for distance in range(1, length//2+1):
                for name, direction in (("parallel", axis), ("perpendicular", (axis+1) % 3)):
                    other = model.cell_index(model.shift(model.cells[cell], direction, distance))
                    insertion = contract_pair(model, ket, bra, sz_sites=(first, 3*other+axis))
                    estimator = (3*insertion.sign*math.exp(insertion.log_abs_value-current.log_abs_value)
                                 if insertion.numerator else 0.)
                    entry[f"spin_{name}_{distance}"] = estimator
            records.append(entry)
        if (step + 1) % 1000 == 0:
            print(json.dumps(dict(step=step+1, seconds=time.perf_counter()-started,
                                  moves=dict(moves), negative_proposals=negative_proposals)), flush=True)
    result = dict(L=length, steps=steps, seed=seed, initial=initial, bias=bias,
                  seconds=time.perf_counter()-started, moves=dict(moves),
                  negative_proposals=negative_proposals, records=records)
    output.mkdir(parents=True, exist_ok=True)
    path = output/f"mcmc_L{length}_{initial}_{seed}_bias{bias:g}.json"
    path.write_text(json.dumps(result, indent=2))
    print(path, flush=True)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=3)
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=10413)
    parser.add_argument("--initial", choices=["diagonal", "independent", "trimer", "dimer"], default="diagonal")
    parser.add_argument("--bias", type=float, default=1.)
    parser.add_argument("--output", type=Path, default=Path("results/octahedral_bulk_overlap"))
    args = parser.parse_args()
    run(args.L, args.steps, args.seed, args.initial, args.bias, args.output)


if __name__ == "__main__":
    main()
