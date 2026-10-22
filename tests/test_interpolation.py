"""Independent finite-torus, local-isometry and singlet-overlap checks."""
import numpy as np
import pytest

from rvb.exact import occupancy_statistics, perfect_matchings, rvb_state
from rvb.interpolation import (
    bond_operator, interpolation_projector, kagome_interpolating_tensor,
    local_qutrit_z, local_spin_z, location_projector,
    normalized_location_overlap, paper_loop_parameter, qutrit_embedding,
)
from rvb.lattices import Lattice, kagome
from rvb.peps import IDENTITY, SM, SP, SZ, column_transfer, double_tensor, kagome_site_tensor


def test_rvb_endpoint_reduces_to_original_spin_tensor_and_observables():
    enlarged, original = kagome_interpolating_tensor(0), kagome_site_tensor()
    embedded = enlarged.reshape(4, 4, 4, 3, 3, 3, 3)
    np.testing.assert_array_equal(embedded[:2, :2, :2].reshape(original.shape), original)
    np.testing.assert_allclose(double_tensor(enlarged), double_tensor(original), atol=1e-13)
    dot = np.kron(np.kron(SZ, SZ) + .5 * (np.kron(SP, SM) + np.kron(SM, SP)), IDENTITY)
    zz = np.kron(np.kron(2 * SZ, 2 * SZ), IDENTITY)
    for kind, operator in (("spin_dot", dot), ("zz", zz)):
        np.testing.assert_allclose(double_tensor(enlarged, bond_operator(kind)),
                                   double_tensor(original, operator), atol=1e-13)


def test_qutrit_embedding_and_orthogonal_endpoint():
    j = qutrit_embedding()
    np.testing.assert_allclose(j.T @ j, np.eye(4), atol=1e-15)
    mapped = j @ interpolation_projector(1).reshape(4, 9)
    expected = np.zeros((9, 9))
    expected[[2, 5, 6, 7], [2, 5, 6, 7]] = np.sqrt(2.)
    np.testing.assert_allclose(mapped, expected, atol=1e-15)
    np.testing.assert_allclose(local_qutrit_z("both"), local_spin_z(), atol=1e-15)
    np.testing.assert_allclose(local_qutrit_z("first") @ local_qutrit_z("second"),
                               local_spin_z(), atol=1e-15)
    first, second = location_projector("first"), location_projector("second")
    np.testing.assert_allclose(first + second, np.eye(4), atol=1e-15)
    np.testing.assert_allclose(first @ second, 0, atol=1e-15)


def _covering_gram_oracle(theta):
    """32 states of 12 spins and a 32x32 flag Gram; never allocate 4**12."""
    lattice = kagome(2, 2, periodic_x=True, periodic_y=True)
    coverings = perfect_matchings(lattice)
    states, flags = [], []
    for covering in coverings:
        states.append(rvb_state(lattice, [covering])[0])
        # Up triangles have all three vertices within one unit cell.
        local_flags = np.empty(lattice.nsites, dtype=int)
        for edge in covering:
            a, b = lattice.edges[edge]
            up_triangle = lattice.labels[a][:2] == lattice.labels[b][:2]
            for site in (a, b):
                # First input is down for a and up for b,c.
                first = up_triangle != (lattice.labels[site][2] == 0)
                local_flags[site] = 1 if first else -1
        flags.append(local_flags)
    states, flags = np.asarray(states), np.asarray(flags)
    changed = np.sum(flags[:, None, :] != flags[None, :, :], axis=-1)
    gram = states @ states.T * normalized_location_overlap(theta) ** changed
    # Eq12 has unnormalized singlets and flags.
    norm = gram.sum() * 2 ** (lattice.nsites // 2) * (1 + theta**2) ** lattice.nsites
    return lattice, coverings, float(norm), gram


@pytest.mark.parametrize("theta", [0., .3, .7, 1.])
def test_small_torus_norm_matches_independent_covering_gram(theta):
    lattice, coverings, expected, gram = _covering_gram_oracle(theta)
    local_norm = double_tensor(kagome_interpolating_tensor(theta))
    transfer = column_transfer([local_norm, local_norm])
    norm = (transfer @ transfer).diagonal().sum()
    assert norm == pytest.approx(expected, rel=3e-13)
    assert len(coverings) == 32
    if theta == 1:
        np.testing.assert_allclose(gram, np.eye(32), atol=1e-15)
        assert norm == pytest.approx(32 * 2 ** (3 * lattice.nsites // 2))


def test_orthogonal_endpoint_occupation_matches_classical_coverings():
    lattice = kagome(2, 2, periodic_x=True, periodic_y=True)
    means, connected = occupancy_statistics(lattice)
    reference = next(i for i, e in enumerate(lattice.edges) if set(e) == {0, 1})
    comparison = next(i for i, e in enumerate(lattice.edges) if set(e) == {6, 7})
    tensor = kagome_interpolating_tensor(1)
    norm = double_tensor(tensor)
    transfer = column_transfer([norm, norm])
    inserted = column_transfer([double_tensor(tensor, bond_operator("occupation")), norm])
    partition = (transfer @ transfer).diagonal().sum()
    mean = (inserted @ transfer).diagonal().sum() / partition
    covariance = (inserted @ inserted).diagonal().sum() / partition - mean**2
    assert mean == pytest.approx(means[reference], abs=1e-13)
    assert covariance == pytest.approx(connected[reference, comparison], abs=1e-13)


@pytest.mark.parametrize("theta", [0., .2, .5, 1.])
def test_six_spin_loop_overlap_directly_follows_eq12(theta):
    # An alternating hexagon: all six sites switch their occupied triangle.
    edges = tuple((i, (i + 1) % 6) for i in range(6))
    lattice = Lattice("loop", tuple((float(i), 0.) for i in range(6)), edges,
                      tuple((i,) for i in range(6)))
    a, _ = rvb_state(lattice, [(0, 2, 4)])
    b, _ = rvb_state(lattice, [(1, 3, 5)])
    p = interpolation_projector(theta)
    flag_a = p[:, 0, 2].reshape(2, 2)[:, 0]
    flag_b = p[:, 2, 0].reshape(2, 2)[:, 0]
    flag_a /= np.linalg.norm(flag_a)
    flag_b /= np.linalg.norm(flag_b)
    overlap = abs(a @ b) * (flag_a @ flag_b)**6
    q = normalized_location_overlap(theta)
    assert overlap == pytest.approx(q**6 / 4, abs=1e-14)
    assert paper_loop_parameter(theta) == pytest.approx(q / np.sqrt(2.))


@pytest.mark.parametrize("theta", [-.1, 1.1, np.nan, np.inf])
def test_interpolation_rejects_invalid_parameters(theta):
    with pytest.raises(ValueError):
        kagome_interpolating_tensor(theta)


def test_operator_arguments_are_unambiguous():
    for kind in ("zz", "spin_dot", "occupation"):
        np.testing.assert_allclose(bond_operator(kind), bond_operator(kind).T, atol=1e-15)
    with pytest.raises(ValueError):
        bond_operator("unknown")
    with pytest.raises(ValueError):
        bond_operator(endpoints=(0, 0))
    with pytest.raises(ValueError):
        local_qutrit_z("spin")
