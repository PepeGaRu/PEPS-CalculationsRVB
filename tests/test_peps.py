"""Tests of contraction and observable-channel transfer lengths."""
import numpy as np
import pytest

from rvb.peps import RVBTransfer, column_transfer, rvb_cylinder_benchmark


def test_sparse_column_matches_independent_dense_contraction():
    rng = np.random.default_rng(21)
    tensors = []
    for _ in range(2):
        tensor = np.zeros((9, 9, 9, 9))
        indices = rng.integers(0, 9, size=(150, 4))
        tensor[tuple(indices.T)] = rng.normal(size=150)
        tensors.append(tensor)
    dense = np.einsum("abuv,cdvu->acbd", *tensors).reshape(81, 81)
    np.testing.assert_allclose(column_transfer(tensors).toarray(), dense, atol=1e-14)


def test_kagome_width_two_observable_mode_matches_real_space_decay():
    model = RVBTransfer("kagome", 2)
    result = model.correlations(8)
    spectrum = model.dimer_correlation_length()
    # After the first column the exact thin-cylinder bond correlator has one
    # nonzero exponential mode. This validates spectral weights, normalization,
    # connected subtraction, and distance conventions against direct powers.
    ratios = result.connected[2:] / result.connected[1:-1]
    np.testing.assert_allclose(ratios, spectrum["eigenvalue_ratio_real"], rtol=1e-8)
    assert spectrum["eigenvalue_ratio_imag"] == 0
    assert model.eigen_residual < 1e-10
    assert 0.3 < spectrum["xi"] < 0.4


def test_square_width_scaling_is_not_reported_as_a_bulk_power_law():
    narrow = rvb_cylinder_benchmark("square", 4, 8)
    wider = rvb_cylinder_benchmark("square", 6, 8)
    assert wider["correlation_length"] > narrow["correlation_length"]
    assert wider["correlation_length_spectrum"]["larger_eigenspaces_matched"]
    assert "do not establish" in wider["metadata"]["limitation"]
    assert narrow["eigen_residual"] < 1e-10
    assert wider["eigen_residual"] < 1e-10


@pytest.mark.parametrize("lattice,width", [("square", 3), ("square", 8), ("kagome", 3), ("kagome", 6)])
def test_unsupported_circumferences_fail_before_allocating(lattice, width):
    with pytest.raises(ValueError):
        RVBTransfer(lattice, width)
