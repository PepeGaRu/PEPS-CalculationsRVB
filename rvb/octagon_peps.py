"""Equal-amplitude singlet PEPS for the triangle-decorated square-octagon graph.

Each square cell contains eight physical spins, ordered
``(r0, r1, r2, r3, dE, dN, dW, dS)``.  Its four triangles have cyclically
oriented edges ``(r_k, r_(k-1), d_k)``.  A dummy edge carries two physical
spins: either they form its intercell singlet or both participate in triangle
singlets.  No spin or vacancy is removed from the Hilbert space.

The four virtual legs, ordered left/right/up/down, have dimension three.
States 0 and 1 carry an intercell singlet; state 2 means that the corresponding
dummy-endpoint spin is paired inside the cell.  All singlets are unnormalized,
which is a common global factor because every covering contains equally many.
"""

from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import numpy as np


EPSILON = np.array([[0., 1.], [-1., 0.]])
TRIANGLES = tuple((k, (k - 1) % 4, 4 + k) for k in range(4))
INTRACELL_EDGES = tuple((t[k], t[(k + 1) % 3])
                        for t in TRIANGLES for k in range(3))
SPIN_LABELS = ("r0", "r1", "r2", "r3", "dE", "dN", "dW", "dS")


def _internal_matchings(remaining):
    """All matchings of the selected spins using triangle edges only."""
    if not remaining:
        yield ()
        return
    first = min(remaining)
    for a, b in INTRACELL_EDGES:
        if first not in (a, b) or a not in remaining or b not in remaining:
            continue
        for rest in _internal_matchings(remaining - {a, b}):
            yield ((a, b),) + rest


@lru_cache(maxsize=1)
def octagon_local_tensor():
    """Return the exact integer tensor with axes ``(spin, left,right,up,down)``.

    The singlet arrows on dummy edges point east and north.  Reversing selected
    arrows generally changes the quantum state; these conventions therefore
    form part of the model definition, not merely the plotting convention.
    """
    tensor = np.zeros((256, 3, 3, 3, 3))
    for virtual in product(range(3), repeat=4):
        left, right, up, down = virtual
        corner_virtual = (right, up, left, down)
        external = tuple(k for k, q in enumerate(corner_virtual) if q != 2)
        if len(external) % 2:
            continue
        remaining = set(range(8)) - {4 + k for k in external}
        for matching in _internal_matchings(remaining):
            for bits in product((0, 1), repeat=len(matching)):
                physical = [0] * 8
                sign = 1
                for (a, b), bit in zip(matching, bits):
                    physical[a], physical[b] = bit, 1 - bit
                    sign *= 1 if bit == 0 else -1
                for k in external:
                    q = corner_virtual[k]
                    if k < 2:  # Outgoing east/north: epsilon[physical, virtual].
                        physical[4 + k] = 1 - q
                        sign *= 1 if q == 1 else -1
                    else:  # Incoming west/south: identity[physical, virtual].
                        physical[4 + k] = q
                index = sum(bit << (7 - i) for i, bit in enumerate(physical))
                tensor[(index, *virtual)] += sign
    tensor.setflags(write=False)
    return tensor


@dataclass(frozen=True)
class OctagonCompression:
    """An exact physical basis reduction; virtual dimensions are unchanged."""

    tensor: np.ndarray
    isometry: np.ndarray
    singular_values: np.ndarray
    rank: int
    compression_error: float

    def project_operator(self, operator):
        """Return U† O U, retaining physical observables in the reduced basis."""
        operator = np.asarray(operator)
        if operator.shape != (256, 256):
            raise ValueError("Physical operator must have shape (256, 256)")
        return self.isometry.conj().T @ operator @ self.isometry


@lru_cache(maxsize=1)
def octagon_compressed_tensor():
    """Reduce only the physical support by an SVD with relative cutoff 1e-12.

    The discarded singular values vanish to floating-point precision.  This
    is not a PEPS/environment truncation and does not approximate the state.
    ``isometry @ tensor.reshape(rank,81)`` reconstructs the original tensor.
    """
    original = octagon_local_tensor().reshape(256, 81)
    u, s, vh = np.linalg.svd(original, full_matrices=False)
    rank = int(np.count_nonzero(s > 1e-12 * s[0]))
    isometry = u[:, :rank]
    reduced = (s[:rank, None] * vh[:rank]).reshape(rank, 3, 3, 3, 3)
    error = float(np.linalg.norm(original - isometry @ reduced.reshape(rank, 81))
                  / np.linalg.norm(original))
    for array in (isometry, reduced, s):
        array.setflags(write=False)
    return OctagonCompression(reduced, isometry, s, rank, error)


def _site_index(site):
    if not isinstance(site, (int, np.integer)) or not 0 <= site < 8:
        raise ValueError("Physical spin index must be an integer in 0..7")
    return int(site)


def physical_spin_operator(site, component="z"):
    """Spin-1/2 operator on one of the eight spins, in the full physical basis."""
    site = _site_index(site)
    components = {
        "z": np.diag([.5, -.5]),
        "plus": np.array([[0., 1.], [0., 0.]]),
        "minus": np.array([[0., 0.], [1., 0.]]),
        "x": np.array([[0., .5], [.5, 0.]]),
        "y": np.array([[0., -.5j], [.5j, 0.]]),
    }
    if component not in components:
        raise ValueError("component must be z, plus, minus, x, or y")
    return np.kron(np.kron(np.eye(2**site), components[component]),
                   np.eye(2**(7 - site)))


def physical_bond_operator(i, j, kind="spin_dot"):
    """Intracell S_i·S_j, Pauli-Z_i Pauli-Z_j, or singlet projector.

    A singlet projector is a quantum spin observable; it is not the classical
    occupation of an edge in a nonorthogonal covering expansion.
    """
    i, j = _site_index(i), _site_index(j)
    if i == j:
        raise ValueError("Choose two distinct physical spins")
    z = physical_spin_operator(i, "z") @ physical_spin_operator(j, "z")
    if kind == "zz":
        return 4 * z
    if kind not in ("spin_dot", "singlet_projector"):
        raise ValueError("kind must be spin_dot, zz, or singlet_projector")
    dot = z + .5 * (physical_spin_operator(i, "plus") @ physical_spin_operator(j, "minus")
                    + physical_spin_operator(i, "minus") @ physical_spin_operator(j, "plus"))
    return .25 * np.eye(256) - dot if kind == "singlet_projector" else dot


def projected_spin_operator(site, component="z", compression=None):
    """Project a single-spin observable into the compressed physical basis."""
    compression = octagon_compressed_tensor() if compression is None else compression
    return compression.project_operator(physical_spin_operator(site, component))


def projected_bond_operator(i, j, kind="spin_dot", compression=None):
    """Project an intracell bond observable into the compressed physical basis."""
    compression = octagon_compressed_tensor() if compression is None else compression
    return compression.project_operator(physical_bond_operator(i, j, kind))
