"""Exact uniform classical dimer/trimer coverings of corner-sharing octahedra.

Octahedron centres form a periodic simple-cubic graph. Its 3N bond midpoints
are the physical covering sites. At a centre the allowed objects are any
subset of its three opposite-port dimers OR one of its eight face trimers.
Two complementary trimers are excluded. There are no quantum spin states,
singlet-overlap factors, fugacities or Hamiltonian energies in this module.

Every complete, owner-labelled covering has the same probability. All
allowed sectors are included; sampling does not rely on local cube moves.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from numbers import Integral
import random
from typing import Iterable


@dataclass(frozen=True)
class LocalState:
    """An allowed local support; port order is +x,-x,+y,-y,+z,-z."""

    index: int
    name: str
    support: int
    dimer_axes: tuple[int, ...]
    face_signs: tuple[int, int, int] | None

    @property
    def incoming(self) -> tuple[int, ...]:
        return tuple((self.support >> p) & 1 for p in range(6))

    @property
    def is_trimer(self) -> bool:
        return self.face_signs is not None


def _local_states() -> tuple[LocalState, ...]:
    states = []
    for subset in range(8):
        axes = tuple(a for a in range(3) if (subset >> a) & 1)
        support = sum(3 << (2*a) for a in axes)
        name = "+".join("D" + "xyz"[a] for a in axes) or "empty"
        states.append(LocalState(len(states), name, support, axes, None))
    for signs in product((-1, 1), repeat=3):
        support = sum(1 << (2*a + (s == -1)) for a, s in enumerate(signs))
        name = "T(" + ",".join("+" if s == 1 else "-" for s in signs) + ")"
        states.append(LocalState(len(states), name, support, (), signs))
    return tuple(states)


LOCAL_STATES = _local_states()
STATE_BY_SUPPORT = {state.support: state for state in LOCAL_STATES}


@dataclass(frozen=True)
class CoveringObject:
    """An owner-labelled object, including its two or three physical site IDs.

    Ownership remains part of the label on small periodic clusters. In
    particular, size-two tori can have geometrically different chords with
    the same unordered pair of periodic site IDs.
    """

    kind: str
    cell: int
    sites: tuple[int, ...]
    axis: int | None = None
    face_signs: tuple[int, int, int] | None = None


def binary_rank(rows: Iterable[int]) -> int:
    """Rank of binary row bitmasks; no floating-point arithmetic."""
    basis = {}
    for row in rows:
        row = int(row)
        if row < 0:
            raise ValueError("Bitmask rows must be nonnegative")
        while row:
            pivot = row.bit_length()-1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


class OctahedralCoverings:
    """Uniform owner-labelled coverings, solved by binary linear algebra.

    Cell c=(x*ny+y)*nz+z has outgoing cubic edge/site ID 3*c+axis.
    Ownership bit 1 points from c to c+axis: its midpoint is owned by that
    positive endpoint. Consequently a_plus=1-bit(outgoing) and
    a_minus=bit(incoming). Constraints equate the three opposite-port
    parities. The corresponding equations in edge bits are homogeneous.
    """

    def __init__(self, nx: int, ny: int, nz: int):
        dimensions = (nx, ny, nz)
        if any(isinstance(n, bool) or not isinstance(n, Integral) or n < 2
               for n in dimensions):
            raise ValueError("Periodic dimensions must be integers >= 2")
        self.shape = self.nx, self.ny, self.nz = tuple(map(int, dimensions))
        self.ncells = self.nx*self.ny*self.nz
        self.nsites = self.nvariables = 3*self.ncells
        self.cells = tuple(product(range(self.nx), range(self.ny), range(self.nz)))
        self.site_labels = tuple((*v, a) for v in self.cells for a in range(3))
        self.site_positions = tuple(
            tuple(float(v[i]) + (.5 if i == a else 0.) for i in range(3))
            for v in self.cells for a in range(3))
        self.edge_endpoints = tuple(
            (c, self.cell_index(self.shift(v, a, 1)))
            for c, v in enumerate(self.cells) for a in range(3))
        self.ports = tuple(
            tuple(edge for a in range(3)
                  for edge in (3*c+a, self.edge_index(self.shift(v, a, -1), a)))
            for c, v in enumerate(self.cells))
        self.axis_parity_masks = tuple(
            tuple((1 << ports[2*a]) ^ (1 << ports[2*a+1]) for a in range(3))
            for ports in self.ports)
        self.equations = tuple(
            masks[0] ^ masks[a]
            for masks in self.axis_parity_masks for a in (1, 2))
        self._parameterize()

    def cell_index(self, cell) -> int:
        if isinstance(cell, Integral) and not isinstance(cell, bool):
            if 0 <= cell < self.ncells:
                return int(cell)
            raise ValueError("Cell index out of range")
        try:
            x, y, z = cell
        except (TypeError, ValueError) as exc:
            raise ValueError("Cell must be an index or three integer coordinates") from exc
        if any(isinstance(v, bool) or not isinstance(v, Integral) for v in (x, y, z)):
            raise ValueError("Cell coordinates must be integers")
        return ((int(x) % self.nx)*self.ny + int(y) % self.ny)*self.nz + int(z) % self.nz

    @staticmethod
    def shift(cell, axis: int, step: int) -> tuple[int, int, int]:
        shifted = list(cell)
        shifted[axis] += step
        return tuple(shifted)

    def edge_index(self, cell, axis: int) -> int:
        self._axis(axis)
        return 3*self.cell_index(cell)+int(axis)

    @staticmethod
    def _axis(axis):
        if isinstance(axis, bool) or not isinstance(axis, Integral) or axis not in (0, 1, 2):
            raise ValueError("Axis must be 0, 1 or 2")

    def _parameterize(self):
        basis = {}
        for row in self.equations:
            while row:
                pivot = row.bit_length()-1
                if pivot not in basis:
                    basis[pivot] = row
                    break
                row ^= basis[pivot]
        self.rank = len(basis)
        masks, dimension = [], 0
        for variable in range(self.nvariables):
            if variable not in basis:
                masks.append(1 << dimension)
                dimension += 1
                continue
            remaining = basis[variable] ^ (1 << variable)
            free_mask = 0
            while remaining:
                bit = remaining & -remaining
                free_mask ^= masks[bit.bit_length()-1]
                remaining ^= bit
            masks.append(free_mask)
        self.parameter_masks = tuple(masks)
        self.dimension = self.free_variables = dimension

    @property
    def covering_count(self) -> int:
        return 1 << self.dimension

    @property
    def all_minus_trimer_reference(self) -> int:
        return (1 << self.nvariables)-1

    def valid(self, bits: int) -> bool:
        return (isinstance(bits, Integral) and not isinstance(bits, bool)
                and 0 <= bits < 1 << self.nvariables
                and all((int(bits) & row).bit_count() % 2 == 0 for row in self.equations))

    def validate(self, bits: int):
        if not self.valid(bits):
            raise ValueError("Invalid covering arrow assignment")

    def bits_from_parameters(self, parameters: int) -> int:
        if (not isinstance(parameters, Integral) or isinstance(parameters, bool)
                or not 0 <= parameters < 1 << self.dimension):
            raise ValueError("Free-parameter assignment out of range")
        return sum(((mask & int(parameters)).bit_count() % 2) << e
                   for e, mask in enumerate(self.parameter_masks))

    def sample(self, count: int = 1, seed: int | None = None) -> tuple[int, ...]:
        """Independent exact uniform samples, with no equilibration phase."""
        if not isinstance(count, Integral) or isinstance(count, bool) or count < 0:
            raise ValueError("Sample count must be a nonnegative integer")
        rng = random.Random(seed)
        return tuple(self.bits_from_parameters(rng.getrandbits(self.dimension))
                     for _ in range(count))

    def _free_mask(self, mask: int) -> int:
        if (not isinstance(mask, Integral) or isinstance(mask, bool)
                or not 0 <= mask < 1 << self.nvariables):
            raise ValueError("Parity mask out of range")
        remaining, result = int(mask), 0
        while remaining:
            bit = remaining & -remaining
            result ^= self.parameter_masks[bit.bit_length()-1]
            remaining ^= bit
        return result

    def probability_parities(self, constraints: Iterable[tuple[int, int]]) -> Fraction:
        """P(all dot(edge_mask,bits)=rhs mod 2), in the complete uniform space."""
        basis = {}
        for edge_mask, rhs in constraints:
            if rhs not in (0, 1):
                raise ValueError("Parity right hand side must be binary")
            row, rhs = self._free_mask(edge_mask), int(rhs)
            while row:
                pivot = row.bit_length()-1
                if pivot not in basis:
                    basis[pivot] = row, rhs
                    break
                other, other_rhs = basis[pivot]
                row ^= other
                rhs ^= other_rhs
            if not row and rhs:
                return Fraction(0)
        return Fraction(1, 1 << len(basis))

    def probability(self, constraints: Iterable[tuple[int, int]]) -> Fraction:
        """P(all edge variables equal the specified values)."""
        parity_conditions = []
        for edge, value in constraints:
            if (not isinstance(edge, Integral) or isinstance(edge, bool)
                    or not 0 <= edge < self.nvariables):
                raise ValueError("Edge variable out of range")
            parity_conditions.append((1 << int(edge), value))
        return self.probability_parities(parity_conditions)

    def incoming(self, bits: int, cell) -> tuple[int, ...]:
        self.validate(bits)
        ports = self.ports[self.cell_index(cell)]
        return tuple(((int(bits) >> edge) & 1) ^ (port % 2 == 0)
                     for port, edge in enumerate(ports))

    def local_state(self, bits: int, cell) -> LocalState:
        support = sum(int(a) << p for p, a in enumerate(self.incoming(bits, cell)))
        return STATE_BY_SUPPORT[support]

    @staticmethod
    def _state(state: int | LocalState) -> LocalState:
        if isinstance(state, LocalState):
            if state in LOCAL_STATES:
                return state
            raise ValueError("Unknown local state")
        if isinstance(state, Integral) and not isinstance(state, bool) and 0 <= state < 16:
            return LOCAL_STATES[int(state)]
        raise ValueError("Local state must be an index in 0..15 or a LocalState")

    def state_constraints(self, cell, state: int | LocalState):
        state = self._state(state)
        ports = self.ports[self.cell_index(cell)]
        return tuple((edge, int(a) ^ (port % 2 == 0))
                     for port, (edge, a) in enumerate(zip(ports, state.incoming)))

    def state_probability(self, cell, state: int | LocalState) -> Fraction:
        return self.probability(self.state_constraints(cell, state))

    def joint_state_probability(self, cell1, state1, cell2, state2) -> Fraction:
        return self.probability(self.state_constraints(cell1, state1)
                                + self.state_constraints(cell2, state2))

    def connected_state(self, cell1, state1, cell2, state2) -> Fraction:
        return (self.joint_state_probability(cell1, state1, cell2, state2)
                - self.state_probability(cell1, state1)*self.state_probability(cell2, state2))

    def dimer_constraints(self, cell, axis: int):
        self._axis(axis)
        ports = self.ports[self.cell_index(cell)]
        return ((ports[2*axis], 0), (ports[2*axis+1], 1))

    def dimer_probability(self, cell, axis: int) -> Fraction:
        return self.probability(self.dimer_constraints(cell, axis))

    def connected_dimer(self, cell1, axis1, cell2, axis2) -> Fraction:
        return (self.probability(self.dimer_constraints(cell1, axis1)
                                 + self.dimer_constraints(cell2, axis2))
                - self.dimer_probability(cell1, axis1)*self.dimer_probability(cell2, axis2))

    def trimer_constraints(self, cell):
        # a_plus XOR a_minus = 1 means the two edge-ownership bits agree.
        return ((self.axis_parity_masks[self.cell_index(cell)][0], 0),)

    def trimer_probability(self, cell) -> Fraction:
        return self.probability_parities(self.trimer_constraints(cell))

    def connected_trimer(self, cell1, cell2) -> Fraction:
        return (self.probability_parities(self.trimer_constraints(cell1)
                                          + self.trimer_constraints(cell2))
                - self.trimer_probability(cell1)*self.trimer_probability(cell2))

    @staticmethod
    def face_state(signs) -> LocalState:
        signs = tuple(signs)
        if len(signs) != 3 or any(s not in (-1, 1) for s in signs):
            raise ValueError("Face signs must be three values in {-1,+1}")
        return next(state for state in LOCAL_STATES if state.face_signs == signs)

    def face_probability(self, cell, signs) -> Fraction:
        return self.state_probability(cell, self.face_state(signs))

    def connected_face(self, cell1, signs1, cell2, signs2) -> Fraction:
        return self.connected_state(cell1, self.face_state(signs1),
                                    cell2, self.face_state(signs2))

    def covering(self, bits: int) -> tuple[CoveringObject, ...]:
        """Return actual owner-labelled objects and check complete coverage."""
        self.validate(bits)
        bits = int(bits)
        objects, degree = [], [0]*self.nsites
        for cell, ports in enumerate(self.ports):
            # Global validity was checked once above; avoid repeating that
            # O(N) check for every local object in a large sampled covering.
            support = sum((((bits >> edge) & 1) ^ (port % 2 == 0)) << port
                          for port, edge in enumerate(ports))
            state = STATE_BY_SUPPORT[support]
            for axis in state.dimer_axes:
                objects.append(CoveringObject("dimer", cell, ports[2*axis:2*axis+2],
                                              axis=axis))
            if state.is_trimer:
                sites = tuple(ports[p] for p, a in enumerate(state.incoming) if a)
                objects.append(CoveringObject("trimer", cell, sites,
                                              face_signs=state.face_signs))
        for obj in objects:
            for site in obj.sites:
                degree[site] += 1
        if any(d != 1 for d in degree):
            raise AssertionError("Arrow covering did not cover every physical site exactly once")
        return tuple(objects)

    def cage_net(self, bits: int, reference: int | None = None) -> int:
        """Reference-relative cubic-edge flip pattern; same homogeneous code."""
        reference = self.all_minus_trimer_reference if reference is None else reference
        self.validate(bits)
        self.validate(reference)
        return int(bits) ^ int(reference)

    def cube_flip_mask(self, cell) -> int:
        """The twelve edges of a unit cube; a valid local cage-net move."""
        origin = self.cells[self.cell_index(cell)]
        result = 0
        for axis in range(3):
            perpendicular = [a for a in range(3) if a != axis]
            for offsets in product((0, 1), repeat=2):
                corner = list(origin)
                for a, offset in zip(perpendicular, offsets):
                    corner[a] += offset
                result ^= 1 << self.edge_index(corner, axis)
        return result

    @property
    def cube_move_rank(self) -> int:
        return binary_rank(self.cube_flip_mask(c) for c in range(self.ncells))


# Both names describe the same corner-sharing octahedral geometry.
OctochloreCoverings = OctahedralCoverings
