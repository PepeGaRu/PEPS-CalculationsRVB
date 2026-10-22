"""Fully packed dimers on the user's decorated square--octagon lattice.

Each original square-lattice vertex contains four triangles. Their four
shared ring sites and four outward half-link sites are permanent spin-1/2s.
If neither half-link is paired inside its triangle, the two are paired to
one another. Thus ``empty dummy link`` never means an unpaired spin.

Orthogonal covering statistics below are exact affine GF(2) calculations.
They are not the nonorthogonal spin-singlet RVB norm or spin correlations.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from .lattices import Lattice


DIRECTIONS = ("E", "N", "W", "S")


@dataclass(frozen=True)
class SquareOctagonGeometry:
    """Periodic cell convention shared by exact oracles and PEPS code.

    Sites in cell ``c=x*ny+y`` are ``8*c+[r0,r1,r2,r3,dE,dN,dW,dS]``.
    Corner k=E,N,W,S contains the CCW triangle ``(r_k,r_{k-1},d_k)``.
    A link joins dE(x,y) to dW(x+1,y), or dN(x,y) to dS(x,y+1).
    ``triangle_edges[c][k]`` lists its three cyclic oriented edge indices.
    ``dummy_edges[c]`` lists the east and north oriented link edge indices.
    """

    nx: int
    ny: int
    lattice: Lattice
    triangles: tuple[tuple[tuple[int, int, int], ...], ...]
    triangle_edges: tuple[tuple[tuple[int, int, int], ...], ...]
    dummy_edges: tuple[tuple[int, int], ...]
    edge_kinds: tuple[str, ...]

    def site(self, x: int, y: int, slot: int) -> int:
        if not 0 <= slot < 8:
            raise ValueError("slot must lie between 0 and 7")
        return 8 * ((x % self.nx) * self.ny + y % self.ny) + slot


@dataclass(frozen=True)
class DummyContractedSquareGeometry:
    """Loop graph after contracting dummy links, not physical spin sites.

    Vertex ``2*c+axis`` is the midpoint of cell c's east/north dummy edge.
    Edge ``4*c+k`` is the surviving small-square edge k -> k+1, exactly
    the first 4N edge IDs of ``base_loop_pattern``. Period-one geometries
    can have parallel edges; their IDs remain distinct.

    ``endpoint_ring_edges[v]`` gives the two pairs incident at the two
    ORIGINAL endpoints of dummy link v. These specify the inherited
    routing when all four edges at the contracted vertex are occupied.
    """

    lattice: Lattice
    endpoint_ring_edges: tuple[tuple[tuple[int, int], tuple[int, int]], ...]


def dummy_contracted_square_geometry(nx: int, ny: int) -> DummyContractedSquareGeometry:
    """Contract the 2N dummy links into vertices of a rotated square graph.

    Coordinates use the original cell translations. H(x,y)=(x+.5,y),
    V(x,y)=(x,y+.5); (U,V)=(X+Y-.5,Y-X+.5) maps them to integer square
    coordinates. This has 2N vertices and 4N edges, unlike the N-vertex
    projection obtained by collapsing entire small squares.
    """
    if any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1
           for n in (nx, ny)):
        raise ValueError("nx and ny must be positive integers")

    def cell(x, y):
        return (x % nx) * ny + y % ny

    labels = tuple((x, y, axis) for x in range(nx) for y in range(ny)
                   for axis in range(2))
    positions = tuple((x + (.5 if axis == 0 else 0.),
                       y + (.5 if axis == 1 else 0.)) for x, y, axis in labels)
    edges, endpoint_pairs = [], []
    for x in range(nx):
        for y in range(ny):
            c = cell(x, y)
            vertices = (2*c, 2*c+1, 2*cell(x-1, y), 2*cell(x, y-1)+1)
            edges.extend((vertices[k], vertices[(k+1) % 4]) for k in range(4))
            ce, cn = cell(x+1, y), cell(x, y+1)
            endpoint_pairs.extend((
                ((4*c, 4*c+3), (4*ce+2, 4*ce+1)),
                ((4*c+1, 4*c), (4*cn+3, 4*cn+2))))
    return DummyContractedSquareGeometry(
        Lattice("dummy-contracted square loop graph", positions, tuple(edges), labels),
        tuple(endpoint_pairs))


def square_octagon_geometry(nx: int = 2, ny: int = 2) -> SquareOctagonGeometry:
    """Create the 8-spin, 14-bond cell with periodic boundary conditions.

    Period-one tori are valid finite graphs: east and west half-link sites
    remain distinct, as do north and south. Such tori are validation clusters,
    not estimates of the plane. All triangle singlets point CCW, while link
    singlets point east/north; changing these signs can change the RVB state.
    """
    if any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1
           for n in (nx, ny)):
        raise ValueError("nx and ny must be positive integers")
    offsets = ((.17, .17), (-.17, .17), (-.17, -.17), (.17, -.17),
               (.37, 0.), (0., .37), (-.37, 0.), (0., -.37))
    labels = tuple((x, y, slot) for x in range(nx) for y in range(ny)
                   for slot in range(8))
    positions = tuple((x + offsets[s][0], y + offsets[s][1]) for x, y, s in labels)

    def site(x, y, slot):
        return 8 * ((x % nx) * ny + y % ny) + slot

    edges, triangles, triangle_edges, edge_kinds = [], [], [], []
    for x in range(nx):
        for y in range(ny):
            local_triangles, local_edges = [], []
            for k in range(4):
                triangle = (site(x, y, k), site(x, y, (k - 1) % 4),
                            site(x, y, 4 + k))
                local_triangles.append(triangle)
                ids = tuple(range(len(edges), len(edges) + 3))
                local_edges.append(ids)
                edges.extend((triangle[j], triangle[(j + 1) % 3]) for j in range(3))
                edge_kinds.extend(["triangle"] * 3)
            triangles.append(tuple(local_triangles))
            triangle_edges.append(tuple(local_edges))
    dummy_edges = []
    for x in range(nx):
        for y in range(ny):
            dummy_edges.append((len(edges), len(edges) + 1))
            edges.extend(((site(x, y, 4), site(x + 1, y, 6)),
                          (site(x, y, 5), site(x, y + 1, 7))))
            edge_kinds.extend(["dummy"] * 2)
    if len({frozenset(edge) for edge in edges}) != len(edges):
        raise AssertionError("Repeated physical edge in decorated geometry")
    lattice = Lattice("square-octagon decorated", positions, tuple(edges), labels)
    return SquareOctagonGeometry(nx, ny, lattice, tuple(triangles),
                                 tuple(triangle_edges), tuple(dummy_edges), tuple(edge_kinds))


def square_octagon(nx: int = 2, ny: int = 2) -> Lattice:
    """Return the graph alone, suitable for :mod:`rvb.exact` oracles."""
    return square_octagon_geometry(nx, ny).lattice


class SquareOctagonOrthogonalDimers:
    """Uniform perfect matchings, represented without enumerating coverings.

    Variables z[c,k] say that the ring site r_k belongs to triangle k.
    The other triangle at that site consequently has occupation 1-z[c,k].
    Variable t[e]=1 means BOTH half-link spins belong to their triangles;
    t[e]=0 means their internal link singlet is present. At corner k,

        z[c,k] + z[c,k-1] + t[e(k)] = 1 (mod 2).

    A triangle therefore contains zero or two occupied slots. Each solution
    gives one perfect matching. On a connected N-cell torus there are
    2**(2*N+1) matchings. Their projection onto occupied square-lattice links
    t is 2**N-to-one, with two local resolutions at each original vertex.

    In the plane these uniform even subgraphs can be written as domain walls
    of independent fair Ising variables on the dual square lattice. Distinct
    t edges are therefore pairwise independent. Given t, resolutions at
    different cells are independent, and a triangle-bond mean is either
    t[e]/2 or (1-t[e])/2 for its incident dummy link. Consequently dimer
    occupations in non-neighboring cells are exactly uncorrelated. This does
    not imply the same behavior for singlet RVB or loop-connectivity probes.
    """

    def __init__(self, nx: int = 8, ny: int = 8):
        self.geometry = square_octagon_geometry(nx, ny)
        self.nx, self.ny = nx, ny
        self.ncells = nx * ny
        self.nvariables = 6 * self.ncells
        self.lattice = self.geometry.lattice
        self.edges = self.lattice.edges
        self.edge_constraints: dict[int, tuple[tuple[int, int], ...]] = {}
        self._edge_lookup = {frozenset(e): i for i, e in enumerate(self.edges)}
        self.equations = []
        self._corner_variables = []
        for x in range(nx):
            for y in range(ny):
                c = x * ny + y
                corner_vars = []
                for k in range(4):
                    # Each slot is an affine function variable XOR constant.
                    slots = ((4*c+k, 0), (4*c+(k-1) % 4, 1),
                             (self.dummy_variable(x, y, k), 0))
                    corner_vars.append(slots)
                    mask = 0
                    for variable, _ in slots:
                        mask ^= 1 << variable
                    self.equations.append((mask, 1))
                    for j, edge in enumerate(self.geometry.triangle_edges[c][k]):
                        self.edge_constraints[edge] = tuple(
                            (variable, 1 ^ constant)
                            for variable, constant in (slots[j], slots[(j+1) % 3]))
                self._corner_variables.append(tuple(corner_vars))
                for axis, edge in enumerate(self.geometry.dummy_edges[c]):
                    self.edge_constraints[edge] = ((4*self.ncells+2*c+axis, 0),)
        self._parameterize()

    @property
    def covering_count(self) -> int:
        return 1 << self.free_variables

    def dummy_variable(self, x: int, y: int, corner: int) -> int:
        if corner == 2:
            x -= 1
        if corner == 3:
            y -= 1
        return 4*self.ncells + 2*((x % self.nx)*self.ny + y % self.ny) + corner % 2

    def _parameterize(self):
        basis = {}
        for mask, rhs in self.equations:
            while mask:
                pivot = mask.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = (mask, rhs)
                    break
                other, value = basis[pivot]
                mask, rhs = mask ^ other, rhs ^ value
            if not mask and rhs:
                raise ValueError("The covering constraints are inconsistent")
        self._masks, self._constants, self.free_variables = [], [], 0
        for variable in range(self.nvariables):
            if variable not in basis:
                mask, constant = 1 << self.free_variables, 0
                self.free_variables += 1
            else:
                remaining, constant = basis[variable]
                remaining ^= 1 << variable
                mask = 0
                while remaining:
                    bit = remaining & -remaining
                    other = bit.bit_length() - 1
                    mask ^= self._masks[other]
                    constant ^= self._constants[other]
                    remaining ^= bit
            self._masks.append(mask)
            self._constants.append(constant)
        if self.free_variables != 2*self.ncells + 1:
            raise AssertionError("Unexpected cycle-space dimension")

    def probability(self, constraints) -> Fraction:
        """Exact probability of simultaneous variable=value conditions."""
        basis = {}
        for variable, value in constraints:
            if not 0 <= variable < self.nvariables or value not in (0, 1):
                raise ValueError("Invalid binary-variable constraint")
            mask, rhs = self._masks[variable], int(value) ^ self._constants[variable]
            while mask:
                pivot = mask.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = (mask, rhs)
                    break
                other, other_rhs = basis[pivot]
                mask, rhs = mask ^ other, rhs ^ other_rhs
            if not mask and rhs:
                return Fraction(0)
        return Fraction(1, 1 << len(basis))

    def edge_index(self, edge) -> int:
        if isinstance(edge, (int, np.integer)):
            if 0 <= edge < len(self.edges):
                return int(edge)
            raise ValueError("Edge index out of bounds")
        try:
            return self._edge_lookup[frozenset(edge)]
        except KeyError as exc:
            raise ValueError("Not an edge of this geometry") from exc

    def occupation(self, edge) -> Fraction:
        return self.probability(self.edge_constraints[self.edge_index(edge)])

    def joint_occupation(self, first, second) -> Fraction:
        return self.probability(self.edge_constraints[self.edge_index(first)]
                                + self.edge_constraints[self.edge_index(second)])

    def connected(self, first, second) -> Fraction:
        return self.joint_occupation(first, second) - self.occupation(first)*self.occupation(second)

    def bits_from_parameters(self, parameters: int) -> int:
        if not 0 <= parameters < 1 << self.free_variables:
            raise ValueError("Free-variable assignment is out of bounds")
        return sum((((mask & parameters).bit_count() % 2) ^ constant) << i
                   for i, (mask, constant) in enumerate(zip(self._masks, self._constants)))

    def covering_from_bits(self, bits: int) -> tuple[int, ...]:
        if not 0 <= bits < 1 << self.nvariables:
            raise ValueError("Binary assignment is out of bounds")
        if any((bits & mask).bit_count() % 2 != rhs for mask, rhs in self.equations):
            raise ValueError("Assignment violates a triangle covering constraint")
        return tuple(edge for edge, conditions in self.edge_constraints.items()
                     if all((bits >> variable) & 1 == value for variable, value in conditions))

    def bits_from_covering(self, covering) -> int:
        """Inverse map, rejecting partial matchings and repeated edges."""
        covering = tuple(self.edge_index(edge) for edge in covering)
        selected = set(covering)
        degree = [0] * self.lattice.nsites
        for edge in covering:
            for site in self.edges[edge]:
                degree[site] += 1
        if len(selected) != len(covering) or any(d != 1 for d in degree):
            raise ValueError("Not a perfect matching")
        bits = 0
        for c in range(self.ncells):
            for k in range(4):
                site = 8*c+k
                belongs = any(edge in selected and site in self.edges[edge]
                              for edge in self.geometry.triangle_edges[c][k])
                bits |= int(belongs) << (4*c+k)
            for axis, edge in enumerate(self.geometry.dummy_edges[c]):
                bits |= int(edge not in selected) << (4*self.ncells+2*c+axis)
        if set(self.covering_from_bits(bits)) != selected:
            raise ValueError("Matching fails the dummy-link rule")
        return bits

    def coverings(self, max_coverings: int = 100_000):
        if self.covering_count > max_coverings:
            raise ValueError("Covering enumeration limit exceeded")
        return [self.covering_from_bits(self.bits_from_parameters(parameters))
                for parameters in range(self.covering_count)]

    def square_loop_pattern(self, bits: int) -> tuple[int, ...]:
        """Collapse SMALL SQUARES: occupied links ordered cell east,north.

        This 2^N-to-one projection is distinct from contracting DUMMY links.
        It keeps the dummy-edge bits and loses small-square resolutions.
        Degree-four vertices need their local resolution to specify loops.
        For the user's dummy-erasure map use dummy_contracted_loop_pattern.
        """
        self.covering_from_bits(bits)  # validate the supplied assignment
        return tuple((bits >> (4*self.ncells+i)) & 1 for i in range(2*self.ncells))

    def base_loop_pattern(self, bits: int, reference: int | None = None) -> tuple[int, ...]:
        """Reference-relative even subgraph on the trivalent 4.8.8 graph."""
        self.covering_from_bits(bits)
        reference = self.bits_from_parameters(0) if reference is None else reference
        self.covering_from_bits(reference)
        return tuple(((bits ^ reference) >> i) & 1 for i in range(self.nvariables))

    def dummy_contracted_loop_pattern(self, bits: int,
                                      reference: int | None = None) -> tuple[int, ...]:
        """Bijective reference-relative loop map after contracting dummy links.

        The returned 4N bits are the surviving ring edges, ordered cell then
        E->N, N->W, W->S, S->E. Their graph and inherited local routing are
        given by dummy_contracted_square_geometry and dummy_contracted_routing.
        """
        return self.base_loop_pattern(bits, reference)[:4*self.ncells]

    def lift_dummy_contracted_loop_pattern(self, pattern) -> tuple[int, ...]:
        """Restore the unique erased dummy bits of an even square subgraph."""
        pattern = tuple(pattern)
        if len(pattern) != 4*self.ncells or any(bit not in (0, 1) for bit in pattern):
            raise ValueError("Supply 4N binary ring-edge loop occupations")
        geometry = dummy_contracted_square_geometry(self.nx, self.ny)
        dummy = []
        for left, right in geometry.endpoint_ring_edges:
            left_parity = int(pattern[left[0]]) ^ int(pattern[left[1]])
            right_parity = int(pattern[right[0]]) ^ int(pattern[right[1]])
            if left_parity != right_parity:
                raise ValueError("Contracted loop graph has an odd-degree vertex")
            dummy.append(left_parity)
        return tuple(int(bit) for bit in pattern) + tuple(dummy)

    def bits_from_dummy_contracted_loops(self, pattern,
                                        reference: int | None = None) -> int:
        """Invert dummy contraction to the unique covering for a reference."""
        reference = self.bits_from_parameters(0) if reference is None else reference
        self.covering_from_bits(reference)
        loops = self.lift_dummy_contracted_loop_pattern(pattern)
        bits = reference ^ sum(bit << i for i, bit in enumerate(loops))
        self.covering_from_bits(bits)
        return bits

    def dummy_contracted_routing(self, pattern) -> tuple[tuple[tuple[int, int], ...], ...]:
        """Pair occupied ring edges at each contracted vertex without ambiguity.

        Four occupied edges use the two inherited endpoint pairs. There is
        no additional freely selectable resolution at a degree-four node.
        Drawing those pairs as separate arcs preserves individual loops.
        """
        loops = self.lift_dummy_contracted_loop_pattern(pattern)
        geometry = dummy_contracted_square_geometry(self.nx, self.ny)
        routings = []
        for endpoint_pairs in geometry.endpoint_ring_edges:
            active = [edge for pair in endpoint_pairs for edge in pair if loops[edge]]
            if len(active) == 4:
                pairs = endpoint_pairs
            elif len(active) == 2:
                pairs = (tuple(active),)
            else:
                pairs = ()
            routings.append(tuple(sorted(tuple(sorted(pair)) for pair in pairs)))
        return tuple(routings)


def square_octagon_orthogonal_benchmark(nx: int = 16, ny: int = 4,
                                        max_distance: int | None = None) -> dict:
    """Axial exact dimer correlations for a triangle bond and a dummy bond."""
    model = SquareOctagonOrthogonalDimers(nx, ny)
    if max_distance is None:
        max_distance = nx // 2
    if not 0 <= max_distance <= nx // 2:
        raise ValueError("Use distances no larger than half the periodic length")
    distances = list(range(max_distance+1))
    result = {
        "lattice": "decorated square-octagon", "unit_cells": [nx, ny],
        "geometry": "periodic torus", "sites": model.lattice.nsites,
        "model": "equal-amplitude orthogonal fully packed dimers",
        "empty_dummy_rule": "the two permanent spins form an internal dummy-edge singlet",
        "method": "exact affine GF(2) occupation probabilities, no sampling",
        "covering_count_log2": model.free_variables,
        "square_loop_preimage_count_log2": model.ncells,
        "square_loop_projection_convention": "Legacy small-square collapse; retains dummy-link bits, not the dummy-contraction map",
        "dummy_contraction_preimage_count_log2": 0,
        "dummy_contraction_vertices": 2*model.ncells,
        "dummy_contraction_edges": 4*model.ncells,
        "plane_occupation_correlation_range": "same or nearest-neighbor original-square cells only",
        "exponential_correlation_length": None,
        "correlation_length_status": "strictly finite range; no nonzero exponential tail to fit",
        "distance_unit": "simulation-cell translation along x; one cell contains one small square",
        "distance": distances, "families": [],
        "interpretation": (
            "The reported exact zeros are finite-range occupation correlations, "
            "not an exponential fit. These are not spin-singlet RVB correlations "
            "or correlations of loop connectivity. Uniform square even-subgraphs "
            "are not the ordinary close-packed square dimer model."),
    }
    for kind in ("triangle", "dummy"):
        edge_at = (lambda c: model.geometry.triangle_edges[c][0][0]) if kind == "triangle" else (
            lambda c: model.geometry.dummy_edges[c][0])
        reference = edge_at(0)
        targets = [edge_at(x*ny) for x in distances]
        result["families"].append({
            "bond_kind": kind, "reference_edge": list(model.edges[reference]),
            "occupation": float(model.occupation(reference)),
            "connected_correlation": [float(model.connected(reference, e)) for e in targets],
            "nonzero_correlations_with_reference": [
                {"edge_index": e, "edge": list(edge), "connected_correlation": float(value)}
                for e, edge in enumerate(model.edges)
                if (value := model.connected(reference, e)) != 0],
        })
    return result
