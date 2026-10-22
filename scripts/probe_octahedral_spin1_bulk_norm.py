#!/usr/bin/env python3
"""Probe actual cubic 3D spin-one norm contractions with HOTRG.

This is a truncation diagnostic, not a bulk correlation-length calculation.
The rank-eight local tensor is exact; HOTRG's finite chi is approximate.
The single physical Cartesian-basis amplitude supplies a rigorous norm lower
bound and can reject catastrophic truncations independently of convergence.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
import math
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import quimb
import quimb.tensor as qtn

from rvb.octahedral import OctahedralCoverings
from rvb.octahedral_spin1_transfer import Spin1OctahedralCylinder


def norm_network(length: int):
    """Exact periodic L cubed network; labels preserve parallel L=2 bonds."""
    labels, amplitudes = Spin1OctahedralCylinder()._local_entries()
    local = np.zeros((8,) * 6)
    np.add.at(local, tuple(labels.T), amplitudes)
    if length >= 3:
        # quimb's fill-function axis order is z-,z+,y-,y+,x+,x-.
        # Keep the original probe's ordering so finite-chi degenerate-space
        # choices follow the same numerical path on a matching BLAS build.
        ordered = local.transpose(5, 4, 3, 2, 0, 1)
        return qtn.TN3D_from_fill_fn(
            lambda shape: ordered.copy(), length, length, length, 8, cyclic=True)
    tn = qtn.TensorNetwork3D.new(
        Lx=length, Ly=length, Lz=length, site_tag_id="I{},{},{}",
        x_tag_id="X{}", y_tag_id="Y{}", z_tag_id="Z{}",
    )
    for v in product(range(length), repeat=3):
        inds = []
        for axis in range(3):
            previous = list(v)
            previous[axis] = (previous[axis] - 1) % length
            for origin in (v, tuple(previous)):
                inds.append("e" + "_".join(map(str, (*origin, axis))))
        tags = [tn.site_tag(*v), tn.x_tag(v[0]), tn.y_tag(v[1]), tn.z_tag(v[2])]
        tn.add_tensor(qtn.Tensor(local.copy(), inds=inds, tags=tags))
    return tn


def norm_lower_bound(length: int):
    """Use the positive amplitude with each edge colored by its own axis.

    Each dimer then contributes 1/sqrt(3), each axis-ordered trimer 1/sqrt(6).
    Enumerating distinct trimer occupancy fields is exponentially cheaper
    than enumerating coverings, but restricted to rank <=20 here.
    """
    model = OctahedralCoverings(length, length, length)
    n = model.ncells
    columns = [
        sum(((model._free_mask(masks[0]) >> j) & 1) << i
            for i, masks in enumerate(model.axis_parity_masks))
        for j in range(model.dimension)
    ]
    basis = {}
    for vector in columns:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = vector
                break
            vector ^= basis[pivot]
    rank = len(basis)
    result = dict(
        covering_dimension=model.dimension, trimer_field_rank=rank,
        simple_log10_norm_lower_bound=(
            2 * model.dimension * math.log10(2) - n * math.log10(6)),
    )
    if rank <= 20:
        vectors = list(basis.values())
        vector, previous = (1 << n) - 1, 0
        histogram = Counter({vector.bit_count(): 1})
        for i in range(1, 1 << rank):
            gray = i ^ (i >> 1)
            vector ^= vectors[(gray ^ previous).bit_length() - 1]
            previous = gray
            histogram[vector.bit_count()] += 1
        ratio = 3 ** .75 / math.sqrt(6)
        amplitude = (2 ** (model.dimension - rank) * 3 ** (-.75 * n)
                     * sum(count * ratio ** t for t, count in histogram.items()))
        result.update(
            trimer_histogram=dict(sorted(histogram.items())),
            axis_color_amplitude=amplitude,
            log10_axis_color_norm_lower_bound=2 * math.log10(amplitude),
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--chi", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--force-first-step", action="store_true",
                        help="Force one x coarse-graining even on L=2 benchmark")
    args = parser.parse_args()
    if args.L < 2 or args.chi < 1:
        parser.error("L must be >=2 and chi positive")
    result = dict(
        method="HOTRG", geometry="periodic isotropic cubic 3D",
        L=args.L, physical_spins=3 * args.L ** 3, chi=args.chi,
        cutoff=1e-12, canonize=False, optimize="greedy",
        quimb_version=quimb.__version__, force_first_step=args.force_first_step,
        correlation_length=None,
        warning="Finite-chi norm diagnostic; no bulk correlation length inferred.",
        positivity_certificate=norm_lower_bound(args.L),
    )
    network = norm_network(args.L)
    started = time.perf_counter()
    if args.force_first_step:
        network = network.coarse_grain_hotrg(
            "x", max_bond=args.chi, cutoff=1e-12, strip_exponent=True,
            canonize=False, optimize="greedy")
    value, exponent = network.contract_hotrg(
        max_bond=args.chi, cutoff=1e-12, strip_exponent=True,
        canonize=False, optimize="greedy", progbar=True)
    result.update(seconds=time.perf_counter() - started,
                  norm_sign=float(np.sign(value)),
                  log10_abs_norm=float(exponent + np.log10(abs(value))))
    bound = result["positivity_certificate"]
    minimum = bound.get("log10_axis_color_norm_lower_bound",
                        bound["simple_log10_norm_lower_bound"])
    result["passes_necessary_positivity_bound"] = bool(
        result["norm_sign"] > 0 and result["log10_abs_norm"] >= minimum - 1e-10)
    if args.L == 2:
        exact = 332821096192 / 19683
        approximate = result["norm_sign"] * 10 ** result["log10_abs_norm"]
        result.update(exact_norm=exact, relative_norm_error=approximate / exact - 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
