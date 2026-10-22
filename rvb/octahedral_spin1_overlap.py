"""Exact pair contractions on fully periodic three-dimensional covering lattices.

These functions contract individual normalized spin-one covering states, not
the sum over every covering.  They therefore supply an overlap Monte Carlo
kernel, but do not on their own produce a quantum expectation or bulk length.
Dimers are eliminated as Kronecker deltas.  Remaining connected epsilon
networks are contracted with arbitrary-precision integer arithmetic.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import exp, log

import numpy as np
import opt_einsum as oe

from .octahedral import OctahedralCoverings
from .octahedral_spin1 import levi_civita


class ContractionTooLarge(RuntimeError):
    """The exact contraction would exceed the requested operation budget."""


@dataclass(frozen=True)
class PairContraction:
    numerator: int
    dimers: int
    trimers: int
    component_sizes: tuple[int, ...]
    estimated_operations: int
    largest_intermediate: int

    @property
    def sign(self) -> int:
        return (self.numerator > 0) - (self.numerator < 0)

    @property
    def log_abs_value(self) -> float:
        if not self.numerator:
            return -float("inf")
        return (log(abs(self.numerator)) - self.dimers*log(3)/2
                - self.trimers*log(6)/2)

    @property
    def value(self) -> float:
        return self.sign*exp(self.log_abs_value)


_EPSILON = levi_civita().astype(object)
# S^z = -i J. Two insertions therefore add an overall minus sign.
_J = _EPSILON[2]


@lru_cache(maxsize=65536)
def _component(key, max_operations):
    args = []
    for kind, indices in key:
        args.extend((_EPSILON if kind == 3 else _J, list(indices)))
    path, info = oe.contract_path(*args, [], optimize="greedy")
    if int(info.opt_cost) > max_operations:
        # A better path can make a substantial difference on bulk graphs.
        path, info = oe.contract_path(*args, [], optimize="auto-hq")
    if int(info.opt_cost) > max_operations:
        raise ContractionTooLarge(
            f"Estimated {int(info.opt_cost)} operations exceeds {max_operations}")
    value = int(oe.contract(*args, [], optimize=path))
    return value, int(info.opt_cost), int(info.largest_intermediate)


def contract_pair(model: OctahedralCoverings, ket: int, bra: int, *,
                  sz_sites: tuple[int, int] | None = None,
                  max_operations: int = 100_000_000) -> PairContraction:
    """Return <bra|ket> or <bra|Sz_i Sz_j|ket> exactly up to normalization.

    Integer numerators are exact. ``value`` restores the normalized dimer
    and trimer square roots in floating point.  ``sz_sites`` must contain
    two different physical site IDs; rotational invariance gives
    <S_i dot S_j> = 3 <Sz_i Sz_j> for these singlet states.
    """
    n = model.nsites
    if sz_sites is not None:
        if (len(sz_sites) != 2 or sz_sites[0] == sz_sites[1]
                or any(not isinstance(s, int) or not 0 <= s < n for s in sz_sites)):
            raise ValueError("Provide two distinct valid physical site IDs")
    observed = set(sz_sites or ())
    covers = model.covering(ket), model.covering(bra)
    parent = list(range(2*n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        parent[find(a)] = find(b)

    for site in range(n):
        if site not in observed:
            union(site, n+site)
    nd = nt = 0
    for layer, objects in enumerate(covers):
        for obj in objects:
            if obj.kind == "dimer":
                nd += 1
                union(layer*n+obj.sites[0], layer*n+obj.sites[1])
            else:
                nt += 1
    factors = []
    for layer, objects in enumerate(covers):
        for obj in objects:
            if obj.kind == "trimer":
                factors.append((3, tuple(find(layer*n+s) for s in obj.sites)))
    for site in sz_sites or ():
        factors.append((2, (find(n+site), find(site))))
    incident = {}
    for i, (_, indices) in enumerate(factors):
        for index in indices:
            incident.setdefault(index, []).append(i)
    closed_loops = len({find(s) for s in range(2*n)} - set(incident))
    numerator = (-1 if sz_sites is not None else 1)*3**closed_loops
    done, sizes = set(), []
    operations = largest = 0
    for start in range(len(factors)):
        if start in done:
            continue
        vertices, stack = [], [start]
        done.add(start)
        while stack:
            v = stack.pop()
            vertices.append(v)
            for index in factors[v][1]:
                for w in incident[index]:
                    if w not in done:
                        done.add(w)
                        stack.append(w)
        # Relabel indices to reuse many identical small components.
        labels, key = {}, []
        for v in sorted(vertices):
            kind, indices = factors[v]
            mapped = tuple(labels.setdefault(i, len(labels)) for i in indices)
            key.append((kind, mapped))
        value, cost, space = _component(tuple(key), max_operations)
        numerator *= value
        sizes.append(len(vertices))
        operations += cost
        largest = max(largest, space)
        if not numerator:
            break
    return PairContraction(numerator, nd, nt, tuple(sizes), operations, largest)


def line_flip_mask(model: OctahedralCoverings, axis: int, cell=0) -> int:
    """Flip every ownership arrow on one periodic coordinate line."""
    model._axis(axis)
    coordinate = list(model.cells[model.cell_index(cell)])
    result = 0
    for value in range(model.shape[axis]):
        coordinate[axis] = value
        result |= 1 << model.edge_index(coordinate, axis)
    return result
