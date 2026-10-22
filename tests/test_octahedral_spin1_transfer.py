"""Independent finite-size and observable checks of the spin-one cylinder."""
from fractions import Fraction

import numpy as np
import pytest
from scipy import sparse

from rvb.octahedral_spin1 import spin1_generators, spin1_trimer
from rvb.octahedral_spin1_transfer import (
    Spin1OctahedralCylinder,
    full_edge_metric,
    norm_edge_basis,
)


@pytest.fixture(scope="module")
def cylinder():
    return Spin1OctahedralCylinder()


def test_norm_edge_compression_is_exact_but_spin_insertions_need_more_support():
    basis = norm_edge_basis()
    norm = full_edge_metric()
    np.testing.assert_allclose(basis.T @ basis, np.eye(8), atol=3e-16)
    assert np.linalg.matrix_rank(norm, tol=1e-12) == 8
    np.testing.assert_allclose(basis @ (basis.T @ norm @ basis) @ basis.T, norm, atol=1e-15)
    # Compressing an observable's own edge in the norm basis would silently
    # remove physical spin terms. The implementation keeps its internal
    # observable edges uncompressed instead.
    spin = full_edge_metric(spin1_generators()[2])
    assert np.linalg.norm(spin-basis @ (basis.T @ spin @ basis) @ basis.T) > 1


def test_identity_insertion_and_exact_sector_blocks(cylinder):
    norm = cylinder.transfer
    identity = cylinder.observable_transfer("identity")
    assert sparse.linalg.norm(norm-identity) < 2e-11
    sizes = [len(cylinder.sector_indices(s)) for s in range(64)]
    assert sorted(set(sizes)) == [36, 60, 164]
    assert sum(sizes) == 4096
    for observable in (norm, cylinder.observable_transfer("dimer_projector"),
                       cylinder.observable_transfer("spin_z")):
        entries = observable.tocoo()
        assert np.array_equal(cylinder.sector_labels[entries.row], cylinder.sector_labels[entries.col])


def test_two_slice_norm_matches_independent_complete_quantum_covering_sum(cylinder):
    # The previous calculation sums all 8192 normalized spin coverings via
    # independent cube swap networks, without any cylinder edge compression.
    observed = cylinder.transfer.multiply(cylinder.transfer.T).sum()
    expected = Fraction(332821096192, 19683)
    assert observed == pytest.approx(float(expected), rel=3e-14)


def test_two_slice_physical_projector_mean_and_joint_match_exact_cube(cylinder):
    transfer = cylinder.transfer
    projector = cylinder.observable_transfer("dimer_projector")
    norm = transfer.multiply(transfer.T).sum()
    mean = projector.multiply(transfer.T).sum()/norm
    joint = projector.multiply(projector.T).sum()/norm
    assert mean == pytest.approx(float(Fraction(2583378509, 2600164814)), rel=3e-14)
    assert joint == pytest.approx(float(Fraction(99253814, 100006339)), rel=3e-14)


def test_two_slice_spin_correlation_matches_independent_24_spin_network(cylinder):
    # Direct DD/DT/TT contraction of the old 24-spin state, with two physical
    # Cartesian Sz generators at outgoing-y sites (000,y) and (100,y), gives
    # DD=DT=0 and TT=3299920/729 before division by the full norm.
    spin = cylinder.observable_transfer("spin_z")
    norm = cylinder.transfer.multiply(cylinder.transfer.T).sum()
    actual = spin.multiply(spin.T).sum()/norm
    assert actual == pytest.approx(float(Fraction(428355, 1600101424)), rel=3e-13)


def test_formal_one_slice_trace_against_dense_twelve_spin_superposition(cylinder):
    # A one-slice closure identifies +x and -x on each octahedron. Only a
    # trimer is possible, with two equivalent x-port choices per cell.
    # This is a contraction test, not the geometry used for any fitted length.
    state = spin1_trimer()
    for _ in range(3):
        state = np.multiply.outer(state, spin1_trimer())
    cells = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for axis in range(2):
        for vertex, coord in enumerate(cells):
            if coord[axis] == 0:
                other = list(coord)
                other[axis] = 1
                target = cells.index(tuple(other))
                first, second = 3*vertex+axis+1, 3*target+axis+1
                state = state+state.swapaxes(first, second)
    state *= 16
    exact_dense_norm = np.vdot(state, state)
    assert cylinder.transfer.diagonal().sum() == pytest.approx(exact_dense_norm, rel=3e-14)


def test_dimer_only_limit_has_known_norm_and_leading_transfer_scale():
    cylinder = Spin1OctahedralCylinder(trimer_weight=0)
    transfer = cylinder.transfer
    # Four independent longitudinal dimer chains, and four transverse
    # singlet pairs with two identical owner choices each per slice.
    # The latter give amplitude16 and norm factor256 per slice.
    assert transfer.multiply(transfer.T).sum() == pytest.approx(4096**2, rel=3e-14)
    spectra = cylinder.sector_spectra()
    leading = max(abs(record["eigenvalues"][0]) for record in spectra)
    assert leading == pytest.approx(256, rel=3e-14)


def test_observable_tail_agrees_with_centered_propagation_and_window_fit(cylinder):
    result = cylinder.analyze("dimer_projector", range(1, 29))
    assert result["sector"] == 0
    assert result["max_eigen_relative_residual"] < 1e-12
    assert result["spectral_direct_max_abs_difference"] < 1e-16
    values = np.real_if_close(result["correlation"])
    assert np.max(abs(np.imag(values))) < 1e-15
    values = np.real(values)
    relevant = [m for m in result["modes"] if abs(m["weight"]) > 1e-13]
    slow = max(relevant, key=lambda m: m["xi"])
    assert slow["xi"] == pytest.approx(2.490102401687835, rel=2e-12)
    assert slow["weight"].real < 0
    assert slow["ratio"].real < 0
    # Late-window fit uses the directly propagated connected correlation,
    # not the transfer eigenvalues from which the spectral length is defined.
    tail = result["distance"] >= 14
    slope, _ = np.polyfit(result["distance"][tail], np.log(abs(values[tail])), 1)
    assert -1/slope == pytest.approx(slow["xi"], rel=3e-6)
    assert np.all(values[13:-1]*values[14:] < 0)


def test_expanded_longitudinal_edges_and_period_two_projector_match_exact_cube(cylinder):
    raw, incoming, outgoing = cylinder._expanded_longitudinal_slice()
    assert sparse.linalg.norm(outgoing @ raw @ incoming-cylinder.transfer) < 3e-11
    numerator = 0
    for a in range(3):
        for b in range(3):
            elementary = np.zeros((3, 3))
            elementary[a, b] = 1
            connector = sparse.kron(sparse.csr_matrix(full_edge_metric(elementary).T),
                                    sparse.eye(512), format="csr")
            observed = raw @ connector
            # The two physical x bonds of a period-two line are exactly
            # the two spins of its diameter. P0=sum_ab E_ab tensor E_ab/3.
            numerator += observed.multiply(observed.T).sum()/3
    norm = cylinder.transfer.multiply(cylinder.transfer.T).sum()
    assert numerator/norm == pytest.approx(float(Fraction(2583378509, 2600164814)), rel=3e-14)
    norm_connector = sparse.kron(sparse.csr_matrix(full_edge_metric().T), sparse.eye(512), format="csr")
    identity_block = outgoing @ raw @ norm_connector @ raw @ norm_connector @ raw @ incoming
    transfer_cube = cylinder.transfer @ cylinder.transfer @ cylinder.transfer
    assert sparse.linalg.norm(identity_block-transfer_cube)/sparse.linalg.norm(transfer_cube) < 3e-14


def test_longitudinal_projector_has_strong_alternating_exponential_tail(cylinder):
    result = cylinder.analyze("longitudinal_dimer_projector", range(3, 35))
    assert result["insertion_slices"] == 3
    assert result["spectral_direct_max_abs_difference"] < 1e-14
    assert result["mean"] == pytest.approx(.482091798498199, abs=2e-13)
    values = np.real_if_close(result["correlation"])
    assert np.all(values[:-1]*values[1:] < 0)
    assert abs(values[0]) > .06
    relevant = [mode for mode in result["modes"] if abs(mode["weight"]) > 1e-10]
    slow = max(relevant, key=lambda mode: mode["xi"])
    assert slow["xi"] == pytest.approx(2.490102401687835, rel=2e-12)
    assert slow["weight"] == pytest.approx(-.0616426603780744, abs=2e-13)
    for start, stop in ((4, 12), (12, 22), (22, 34)):
        chosen = (result["distance"] >= start) & (result["distance"] <= stop)
        slope, _ = np.polyfit(result["distance"][chosen], np.log(abs(values[chosen])), 1)
        assert -1/slope == pytest.approx(slow["xi"], rel=1e-7)
