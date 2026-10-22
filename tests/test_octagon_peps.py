"""Independent covering-sum checks for the fixed-spin square-octagon PEPS."""

import numpy as np
import pytest

from rvb.exact import bond_expectation, perfect_matchings, rvb_state
from rvb.octagon_peps import (
    octagon_compressed_tensor, octagon_local_tensor, physical_bond_operator,
    physical_spin_operator, projected_bond_operator, projected_spin_operator,
)
from rvb.peps import column_transfer, double_tensor
from rvb.square_octagon import square_octagon


def _contract_small_torus(nx, ny, tensor):
    # Explicit index contractions are independent of the covering enumerator.
    if (nx, ny) == (1, 1):
        return np.einsum("aiijj->a", tensor)
    if (nx, ny) == (1, 2):
        return np.einsum("aiijk,bllkj->ab", tensor, tensor).reshape(-1)
    if (nx, ny) == (2, 1):
        return np.einsum("aijkk,bjill->ab", tensor, tensor).reshape(-1)
    raise ValueError("Only validation tori of one or two cells are supported")


@pytest.mark.parametrize("nx,ny,ncoverings", [(1, 1, 8), (1, 2, 32), (2, 1, 32)])
def test_raw_peps_matches_independent_oriented_covering_sum(nx, ny, ncoverings):
    lattice = square_octagon(nx, ny)
    coverings = perfect_matchings(lattice)
    assert len(coverings) == ncoverings
    state, normalized_singlet_norm = rvb_state(lattice, coverings)
    tensor_state = _contract_small_torus(nx, ny, octagon_local_tensor())
    # The PEPS uses unnormalized singlets. Test amplitudes, not just fidelity:
    # this catches configuration-dependent signs and unintended multiplicities.
    expected = state * np.sqrt(normalized_singlet_norm * 2**(lattice.nsites // 2))
    np.testing.assert_allclose(tensor_state, expected, atol=3e-14)


def test_physical_compression_is_exact_isometry_and_preserves_double_tensor():
    reduced = octagon_compressed_tensor()
    original = octagon_local_tensor()
    assert reduced.rank == 41
    assert reduced.compression_error < 1e-13
    np.testing.assert_allclose(reduced.isometry.T @ reduced.isometry,
                               np.eye(reduced.rank), atol=3e-14)
    np.testing.assert_allclose(reduced.isometry @ reduced.tensor.reshape(41, 81),
                               original.reshape(256, 81), atol=5e-14)
    np.testing.assert_allclose(double_tensor(reduced.tensor), double_tensor(original),
                               atol=2e-12)


@pytest.mark.parametrize("kind", ["spin_dot", "zz", "singlet_projector"])
def test_compression_preserves_actual_bond_insertion(kind):
    reduced = octagon_compressed_tensor()
    full_operator = physical_bond_operator(0, 4, kind=kind)
    operator = projected_bond_operator(0, 4, kind=kind)
    np.testing.assert_allclose(operator, operator.conj().T, atol=2e-14)
    np.testing.assert_allclose(double_tensor(reduced.tensor, operator),
                               double_tensor(octagon_local_tensor(), full_operator),
                               atol=1e-12)


@pytest.mark.parametrize("component", ["z", "plus", "minus", "x", "y"])
def test_compression_preserves_single_spin_insertion(component):
    reduced = octagon_compressed_tensor()
    np.testing.assert_allclose(
        double_tensor(reduced.tensor, projected_spin_operator(4, component)),
        double_tensor(octagon_local_tensor(), physical_spin_operator(4, component)),
        atol=1e-12,
    )


def test_torus_bond_observables_and_singlet_symmetry_match_oracle():
    # The one-cell torus includes physical dummy bonds 4--6 and 5--7.
    oracle, _ = rvb_state(square_octagon(1, 1))
    reduced = octagon_compressed_tensor()
    reduced_state = _contract_small_torus(1, 1, reduced.tensor)
    reduced_state /= np.linalg.norm(reduced_state)
    for edge in ((0, 3), (0, 4), (4, 6), (5, 7)):
        mean = reduced_state @ projected_bond_operator(*edge) @ reduced_state
        assert mean == pytest.approx(bond_expectation(oracle, edge), abs=2e-14)
    total_sz = sum(projected_spin_operator(i, "z") for i in range(8))
    total_sp = sum(projected_spin_operator(i, "plus") for i in range(8))
    np.testing.assert_allclose(total_sz @ reduced_state, 0, atol=2e-14)
    np.testing.assert_allclose(total_sp @ reduced_state, 0, atol=2e-14)


def test_periodic_column_norm_matches_two_cell_covering_sum():
    # Independent double-layer contraction checks the inherited column-axis
    # convention used by the cylinder transfer calculations.
    lattice = square_octagon(2, 1)
    _, norm = rvb_state(lattice)
    transfer = column_transfer([double_tensor(octagon_local_tensor())])
    assert (transfer @ transfer).diagonal().sum() == pytest.approx(
        norm * 2**(lattice.nsites // 2), abs=1e-10)


def test_operator_validation_and_singlet_projector_identity():
    with pytest.raises(ValueError):
        projected_spin_operator(8)
    with pytest.raises(ValueError):
        projected_spin_operator(0, "unknown")
    with pytest.raises(ValueError):
        projected_bond_operator(0, 0)
    with pytest.raises(ValueError):
        projected_bond_operator(0, 4, "occupation")
    p = physical_bond_operator(0, 4, "singlet_projector")
    np.testing.assert_allclose(p @ p, p, atol=1e-15)
