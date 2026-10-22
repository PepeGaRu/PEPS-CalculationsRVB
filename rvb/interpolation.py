"""Schuch et al., arXiv:1203.4816v3, Eq. (12), on the kagome lattice.

Each physical site has basis ``(+ up, + down, - up, - down)``.  Here +/-
are orthonormal flag labels, not spin labels.  At theta=0 the flag is fixed
and the spin state is the usual RVB state.  At theta=1 the flag records which
of the two adjacent triangles contains the dimer, making coverings orthogonal.
The tensors retain the paper's unnormalized convention.

The paper's Fig. 11 specifies qutrit z=diag(1,-1,1) without unambiguously
specifying its embedding in the four-dimensional site space.  Our default
bond observable is explicitly the physical spin Pauli-Z product; it should
not be silently identified with that figure's observable.  Qutrit diagnostic
operators below make alternative embeddings explicit.
"""
from __future__ import annotations

import numpy as np

from .peps import IDENTITY, SM, SP, SZ, triangle_tensor


def _theta(theta: float) -> float:
    value = float(theta)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("theta must be finite and between zero (RVB) and one (orthogonal)")
    return value


def interpolation_projector(theta: float) -> np.ndarray:
    """Return P_theta[flag*2+spin, first_qutrit, second_qutrit], Eq. (12).

    A spin on the first input gets flag |+>+theta|->; a spin on the
    second gets |+>-theta|->.  Virtual qutrit 2 means no dimer endpoint.
    """
    theta = _theta(theta)
    p = np.zeros((4, 3, 3))
    for spin in range(2):
        p[spin, spin, 2] = p[spin, 2, spin] = 1.
        p[2 + spin, spin, 2] = theta
        p[2 + spin, 2, spin] = -theta
    return p


def kagome_interpolating_tensor(theta: float) -> np.ndarray:
    """Return K[(a,b,c),left,right,up,down] with physical dimension 4**3.

    Blocking and singlet arrows exactly match :func:`rvb.peps.kagome_site_tensor`.
    In this blocking, the up triangle is the second input of site a and the
    first input of sites b,c.  Each flattened physical site is flag*2+spin.
    At theta=0, selecting flag=+ at all three sites recovers the spin-only K.
    """
    t, p = triangle_tensor(), interpolation_projector(theta)
    result = np.einsum("ald,bce,pab,qcr,teu->pqtlrud", t, t, p, p, p, optimize=True)
    return result.reshape(64, 3, 3, 3, 3)


def normalized_location_overlap(theta: float) -> float:
    """Overlap of the two normalized local dimer-location flag states."""
    theta = _theta(theta)
    return (1. - theta**2) / (1. + theta**2)


def paper_loop_parameter(theta: float) -> float:
    """Return x(theta) stated below Eq. (12), q(theta)/sqrt(2).

    This parameter alone does not specify an overlap for arbitrary covering
    pairs: direct Eq. (12) contraction gives spin-covering overlap times
    q(theta)**m, where m counts sites whose occupied triangle changes.
    """
    return normalized_location_overlap(theta) / np.sqrt(2.)


def local_spin_z() -> np.ndarray:
    """I_flag tensor Pauli-Z; eigenvalues +/-1, not spin S_z=+/-1/2."""
    return np.kron(IDENTITY, 2 * SZ)


def qutrit_embedding() -> np.ndarray:
    """Isometry J[3*i+j,flag*2+spin] into the two-qutrit physical space.

    J|+,s>=(|s,2>+|2,s>)/sqrt(2), and
    J|-,s>=(|s,2>-|2,s>)/sqrt(2).  Thus J P_1=sqrt(2) P_perp,
    with P_perp as in the paper's Eq. (3).  This is a fixed, explicit
    convention for all theta, including the RVB endpoint.
    """
    j = np.zeros((9, 4))
    for spin in range(2):
        j[3 * spin + 2, spin] = j[6 + spin, spin] = 1 / np.sqrt(2.)
        j[3 * spin + 2, 2 + spin] = 1 / np.sqrt(2.)
        j[6 + spin, 2 + spin] = -1 / np.sqrt(2.)
    return j


def local_qutrit_z(leg: str = "both") -> np.ndarray:
    """J^dagger z_leg J, for qutrit z=diag(1,-1,1).

    ``first`` and ``second`` act on one qutrit; ``both`` acts with z tensor z
    and equals :func:`local_spin_z`.  The one-qutrit choices are distinct
    observables.  None is asserted to reconstruct Fig. 11's unspecified legs.
    """
    z, identity = np.diag([1., -1., 1.]), np.eye(3)
    if leg == "first":
        operator = np.kron(z, identity)
    elif leg == "second":
        operator = np.kron(identity, z)
    elif leg == "both":
        operator = np.kron(z, z)
    else:
        raise ValueError("leg must be 'first', 'second', or 'both'")
    j = qutrit_embedding()
    return j.T @ operator @ j


def location_projector(leg: str) -> np.ndarray:
    """Project onto first/second occupied-qutrit flag at theta=1.

    These fixed projectors remain valid physical observables at other theta,
    but they no longer exactly read a covering's dimer location there.
    """
    if leg not in {"first", "second"}:
        raise ValueError("leg must be 'first' or 'second'")
    vector = np.array([1., 1. if leg == "first" else -1.]) / np.sqrt(2.)
    return np.kron(np.outer(vector, vector), IDENTITY)


def _embed_pair(first: np.ndarray, second: np.ndarray, endpoints: tuple[int, int]) -> np.ndarray:
    local = [np.eye(4) for _ in range(3)]
    local[endpoints[0]], local[endpoints[1]] = first, second
    return np.kron(np.kron(local[0], local[1]), local[2])


def bond_operator(kind: str = "zz", endpoints: tuple[int, int] = (0, 1)) -> np.ndarray:
    """Physical 64x64 operator for a bond within the blocked up triangle.

    ``zz`` means Pauli-Z_i Pauli-Z_j, acting trivially on the flags.
    ``spin_dot`` means S_i dot S_j with spin one-half normalization.
    ``occupation`` means both endpoints are assigned to this up triangle;
    at theta=1 the dimer constraint makes this exactly its bond occupation.
    At theta<1 this last observable is a flag-projector product, not a
    classical dimer occupation.  Indices 0,1,2 are the three sites a,b,c.
    """
    endpoints = tuple(endpoints)
    if len(endpoints) != 2 or len(set(endpoints)) != 2 or any(i not in (0, 1, 2) for i in endpoints):
        raise ValueError("endpoints must name two distinct sites among 0, 1, 2")
    if kind == "zz":
        return _embed_pair(local_spin_z(), local_spin_z(), endpoints)
    if kind == "spin_dot":
        z, plus, minus = (np.kron(IDENTITY, op) for op in (SZ, SP, SM))
        return (_embed_pair(z, z, endpoints)
                + .5 * _embed_pair(plus, minus, endpoints)
                + .5 * _embed_pair(minus, plus, endpoints))
    if kind == "occupation":
        up_legs = ("second", "first", "first")
        first, second = (location_projector(up_legs[i]) for i in endpoints)
        return _embed_pair(first, second, endpoints)
    raise ValueError("kind must be 'zz', 'spin_dot', or 'occupation'")
