"""Spin-one singlet version of the sixteen-state octahedral covering model.

Cartesian spin basis: dimer delta/sqrt(3), trimer epsilon/sqrt(6).
The three trimer sites are always ordered by x,y,z axis, independent of face
signs. Equal coefficients refer to products of these *normalized* objects.

The local tensors define the state on any of the existing covering lattices.
Exact quantum contractions below deliberately specialize to the 2x2x2
periodic cluster. Its 4096 dimer coverings give the same physical state D;
the other 4096 coverings are G T, with G=product_e(1+Swap_e). Consequently
Psi=4096 D+G T. This period-two degeneracy forbids bulk extrapolation.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import sqrt

import numpy as np
import opt_einsum as oe

from .octahedral import LOCAL_STATES


def levi_civita() -> np.ndarray:
    epsilon = np.zeros((3, 3, 3), dtype=np.int64)
    for p in permutations(range(3)):
        epsilon[p] = (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
    return epsilon


def spin1_dimer() -> np.ndarray:
    return np.eye(3)/sqrt(3)


def spin1_trimer() -> np.ndarray:
    return levi_civita()/sqrt(6)


def spin1_generators() -> np.ndarray:
    """Hermitian Cartesian generators S^a_bc=-i epsilon_abc, with hbar=1."""
    return -1j*levi_civita()


def ownership_projection() -> np.ndarray:
    """P[physical, left, right], vacuum=0 and occupied colors=1,2,3."""
    projection = np.zeros((3, 4, 4))
    for color in range(3):
        projection[color, color+1, 0] = 1
        projection[color, 0, color+1] = 1
    return projection


def octahedron_tensor() -> np.ndarray:
    """Six virtual legs ordered +x,-x,+y,-y,+z,-z, each of dimension four."""
    tensor = np.zeros((4,)*6)
    dimer, trimer = spin1_dimer(), spin1_trimer()
    for state in LOCAL_STATES:
        occupied = [p for p in range(6) if state.incoming[p]]
        for colors in product(range(3), repeat=len(occupied)):
            entry = [0]*6
            for p, color in zip(occupied, colors):
                entry[p] = color+1
            if state.is_trimer:
                # occupied ports are sorted by Cartesian axis, as required.
                amplitude = trimer[tuple(colors)]
            else:
                amplitude = 1.0
                for axis in state.dimer_axes:
                    amplitude *= dimer[entry[2*axis]-1, entry[2*axis+1]-1]
            tensor[tuple(entry)] = amplitude
    return tensor


def cube_geometry():
    """Return (cells, edges); edges are (axis, vertex0, vertex1).

    Cells are lexicographic binary Cartesian coordinates. Edges are grouped
    first by axis; vertex0 has coordinate zero along that axis. Each cube
    edge labels a *pair* of physical spins on the periodic 2x2x2 lattice.
    Physical site 3*v+axis belongs to the plus-face reference trimer at v.
    """
    cells = tuple(product(range(2), repeat=3))
    index = {cell: v for v, cell in enumerate(cells)}
    edges = []
    for axis in range(3):
        for cell in cells:
            if cell[axis] == 0:
                other = list(cell)
                other[axis] = 1
                edges.append((axis, index[cell], index[tuple(other)]))
    return cells, tuple(edges)


def cube_contraction(vertex_tensors, edge_tensors):
    """Contract eight rank-three vertex tensors with twelve edge matrices.

    Vertex legs follow x,y,z; edge matrices join the corresponding legs of
    the ordered endpoints returned by cube_geometry(). Works for integer
    unnormalized tensors or normalized floating-point tensors.
    """
    _, edges = cube_geometry()
    if len(vertex_tensors) != 8 or len(edge_tensors) != 12:
        raise ValueError("The cube needs eight vertex tensors and twelve edge tensors")
    exact = all(np.asarray(tensor).dtype.kind in "iuO" for tensor in (*vertex_tensors, *edge_tensors))
    args = []
    for v, tensor in enumerate(vertex_tensors):
        args.extend((np.asarray(tensor, dtype=object) if exact else tensor, [3*v+a for a in range(3)]))
    for (axis, u, v), tensor in zip(edges, edge_tensors):
        args.extend((np.asarray(tensor, dtype=object) if exact else tensor, [3*u+axis, 3*v+axis]))
    args.append([])
    return np.asarray(oe.contract(*args, optimize="greedy")).item()


class Spin1OctahedralTwoCube:
    """Exact rational norms and projector expectations for all 8192 coverings.

    No stochastic sampling, discarded tensor singular values, or Hilbert-
    space vector of length 3^24 is used. Integer tensor contractions are
    evaluated with arbitrary-precision Python integers, then converted to
    fractions after normalized singlet factors are restored.
    """
    def __init__(self):
        self.cells, self.edges = cube_geometry()
        self.dimer_coverings = self.trimer_coverings = 4096
        self.epsilon = levi_civita()
        self.identity = np.eye(3, dtype=np.int64)
        self.symmetrizer = (
            np.einsum("pa,qb->pqab", self.identity, self.identity)
            + np.einsum("pb,qa->pqab", self.identity, self.identity))
        # Pair bra/ket indices at each vertex leg: (a,a'), (b,b'), (c,c').
        self.double_epsilon = np.einsum("abc,def->adbecf", self.epsilon, self.epsilon).reshape((9,)*3)
        self.overlap_edge = (
            np.einsum("ad,be->adbe", self.identity, self.identity)
            + np.einsum("ae,bd->adbe", self.identity, self.identity)).reshape(9, 9)
        self.projector_edge = np.einsum("ab,de->adbe", self.identity, self.identity).reshape(9, 9)

    def reference_overlap(self) -> Fraction:
        numerator = cube_contraction([self.epsilon]*8, [self.identity]*12)
        return Fraction(int(numerator), 3**6*6**4)

    @lru_cache(maxsize=None)
    def _trimer_sum_dimer_expectation(self, edge_ids: tuple[int, ...]) -> Fraction:
        edges = [self.projector_edge if e in edge_ids else self.overlap_edge for e in range(12)]
        numerator = int(cube_contraction([self.double_epsilon]*8, edges))
        return Fraction(4096*numerator, 6**8)*Fraction(2, 3)**len(edge_ids)

    def norm_components(self):
        amplitude = self.dimer_coverings
        cross = amplitude*self.reference_overlap()  # <D|G T>
        tt = self._trimer_sum_dimer_expectation(())
        dd = Fraction(amplitude**2)
        interference = 2*amplitude*cross
        return {"dimer_part": dd, "trimer_part": tt,
                "interference": interference, "total": dd+tt+interference,
                "reference_overlap": self.reference_overlap(),
                "dimer_trimer_sum_overlap": cross}

    def dimer_projector_expectation(self, edge_ids=()) -> Fraction:
        """Product of physical spin-one pair-singlet projectors on distinct pairs."""
        edge_ids = tuple(sorted(set(edge_ids)))
        if any(not isinstance(e, int) or e < 0 or e >= 12 for e in edge_ids):
            raise ValueError("Pair-edge indices must be integers from 0 to 11")
        parts = self.norm_components()
        return (parts["dimer_part"]+parts["interference"]
                +self._trimer_sum_dimer_expectation(edge_ids))/parts["total"]

    def dimer_projector_covariance(self, first: int, second: int) -> Fraction:
        return (self.dimer_projector_expectation((first, second))
                -self.dimer_projector_expectation((first,))*self.dimer_projector_expectation((second,)))

    @lru_cache(maxsize=None)
    def _trimer_projector_branch(self, branches: str, vertices: tuple[int, ...]) -> Fraction:
        """<branch0| product of plus-face trimer projectors |branch1>.

        D denotes the one normalized product of twelve dimers; T denotes
        the unnormalized sum G T_ref of 4096 normalized trimer products.
        The explicit network uses integer epsilon, delta and I+Swap tensors.
        Arbitrary-precision integer arithmetic avoids overflow assumptions.
        """
        observed_sites = {3*v+a for v in vertices for a in range(3)}
        physical, args, counter = {}, [], 24
        denominator = 1
        for layer, branch in enumerate(branches):
            for site in range(24):
                if site not in observed_sites:
                    physical[layer, site] = site
                else:
                    physical[layer, site] = counter
                    counter += 1
            if branch == "D":
                denominator *= 3**6
                for axis, u, v in self.edges:
                    args.extend((self.identity, [physical[layer, 3*u+axis], physical[layer, 3*v+axis]]))
            elif branch == "T":
                denominator *= 6**4
                legs = {}
                for v in range(8):
                    indices = []
                    for axis in range(3):
                        legs[v, axis] = counter
                        indices.append(counter)
                        counter += 1
                    args.extend((self.epsilon, indices))
                for axis, u, v in self.edges:
                    args.extend((self.symmetrizer, [physical[layer, 3*u+axis], physical[layer, 3*v+axis],
                                                   legs[u, axis], legs[v, axis]]))
            else:
                raise ValueError("Branches must be D or T")
        for v in vertices:
            denominator *= 6
            for layer in range(2):
                args.extend((self.epsilon, [physical[layer, 3*v+a] for a in range(3)]))
        for index in range(0, len(args), 2):
            args[index] = np.asarray(args[index], dtype=object)
        args.append([])
        value = np.asarray(oe.contract(*args, optimize="auto-hq")).item()
        return Fraction(int(value), denominator)

    def trimer_projector_expectation(self, vertices=(0,)) -> Fraction:
        """One or two physical plus-face singlet projectors, at cube vertices.

        These reference faces use disjoint physical sites, even at adjacent
        vertices. Products are consequently commuting observables.
        """
        vertices = tuple(sorted(set(vertices)))
        if len(vertices) > 2 or any(not isinstance(v, int) or v < 0 or v >= 8 for v in vertices):
            raise ValueError("Use at most two distinct vertex indices from 0 to 7")
        amplitude = self.dimer_coverings
        numerator = (amplitude**2*self._trimer_projector_branch("DD", vertices)
                     +2*amplitude*self._trimer_projector_branch("DT", vertices)
                     +self._trimer_projector_branch("TT", vertices))
        return numerator/self.norm_components()["total"]

    def trimer_projector_covariance(self, first: int, second: int) -> Fraction:
        return (self.trimer_projector_expectation((first, second))
                -self.trimer_projector_expectation((first,))*self.trimer_projector_expectation((second,)))
