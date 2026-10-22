"""Independent topology and sector checks for the wider spin-one tube."""
import numpy as np
import pytest

from rvb.octahedral_spin1_transfer import Spin1OctahedralCylinder
from rvb.octahedral_spin1_width import Spin1OctahedralCylinder23


@pytest.fixture(scope="module")
def wider():
    return Spin1OctahedralCylinder23()


def test_wider_identity_insertion_preserves_the_exact_transfer(wider):
    matrix, indices, pairs = wider.block()
    identity, other_indices, other_pairs = wider.block(observable="identity")
    assert matrix.shape == (976, 976)
    np.testing.assert_array_equal(indices, other_indices)
    for first, second in zip(pairs, other_pairs):
        np.testing.assert_array_equal(first, second)
    np.testing.assert_allclose(identity, matrix, rtol=2e-14, atol=3e-13)


def test_wider_all_sector_norm_agrees_with_axis_rotated_torus(wider):
    # A length-two 2x3 tube and a length-three 2x2 tube describe the same
    # 36-spin torus after permuting Cartesian axes. This tests all sectors,
    # not just the one used to extract an infinite-length correlation.
    norm, representatives = wider.finite_length_two_norm()
    narrow = Spin1OctahedralCylinder().transfer
    rotated = (narrow @ narrow).multiply(narrow.T).sum()
    assert len(representatives) == 80
    assert sum(record["multiplicity"] for record in representatives) == 256
    assert norm == pytest.approx(rotated, rel=3e-14)


def test_tied_dominant_parities_have_the_same_local_correlations(wider):
    first = wider.analyze(distances=range(1, 25), parity=0)
    second = wider.analyze(distances=range(1, 25), parity=1)
    assert first["lambda0"] == pytest.approx(second["lambda0"], rel=3e-13)
    assert first["mean"] == pytest.approx(second["mean"], abs=3e-13)
    np.testing.assert_allclose(first["correlation"], second["correlation"], rtol=3e-9, atol=3e-18)
    for result in (first, second):
        assert result["dominant_right_residual"] < 1e-12
        assert result["dominant_left_residual"] < 1e-12
        assert result["spectral_direct_max_abs_difference"] < 1e-14
        physical = [mode for mode in result["modes"] if abs(mode["weight"]) > 1e-12]
        slowest = sorted(physical, key=lambda mode: mode["xi"], reverse=True)[:2]
        values = np.real_if_close(result["correlation"])
        chosen = result["distance"] >= 12
        slope, _ = np.polyfit(result["distance"][chosen], np.log(abs(values[chosen])), 1)
        # Two negative eigenvalues have nearby decay lengths. At these
        # distances both still contribute, so a single-exponential fit
        # lies between their scales rather than isolating the slowest one.
        assert slowest[1]["xi"] < -1/slope < slowest[0]["xi"]
        np.testing.assert_allclose(result["correlation"], result["spectral_correlation"],
                                   rtol=1e-8, atol=0)
