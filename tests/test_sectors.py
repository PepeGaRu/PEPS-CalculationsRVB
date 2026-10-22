"""Paper-motivated virtual-sector checks for finite kagome cylinders."""
import copy

import numpy as np
import pytest

from rvb.peps import RVBTransfer, rvb_cylinder_benchmark


def test_default_kagome_sector_changes_with_width_and_is_reported():
    # This is why the original unrestricted two-width curve was not a
    # controlled fixed-sector circumference study.
    for width, expected in ((2, "odd"), (4, "even")):
        model = RVBTransfer("kagome", width)
        assert model.selected_sector == expected
        assert model.sector_requested == "dominant"
        assert model.correlations(2).as_dict()["virtual_vacuum_parity"] == expected


@pytest.mark.parametrize("sector,mean_bond", [
    ("even", -0.17783121635129667),
    ("odd", -0.2319146909766659),
])
def test_explicit_sector_has_definite_ket_and_bra_parity(sector, mean_bond):
    model = RVBTransfer("kagome", 2, sector=sector)
    expected = int(sector == "odd")
    # Decode the virtual two-layer basis independently in base three.
    digits = np.array([np.unravel_index(int(index), (3,) * 4)
                       for index in model.active])
    for layer in (0, 1):
        np.testing.assert_array_equal(np.sum(digits[:, layer::2] == 2, axis=1) % 2,
                                      expected)
    assert model.selected_sector == sector
    assert model.mean_bond == pytest.approx(mean_bond, abs=1e-12)
    assert model.eigen_residual < 1e-11


def test_even_sector_alternating_correlation_matches_exact_thin_cylinder_mode():
    model = RVBTransfer("kagome", 2, sector="even")
    values = model.correlations(9).connected
    expected_ratio = -(2 - np.sqrt(3))
    np.testing.assert_allclose(values[2:] / values[1:-1], expected_ratio,
                               rtol=1e-10, atol=1e-12)
    spectrum = model.dimer_correlation_length()
    assert spectrum["xi"] == pytest.approx(1 / np.log(2 + np.sqrt(3)), abs=1e-12)
    assert spectrum["eigenvalue_ratio_real"] == pytest.approx(expected_ratio, abs=1e-12)


def test_explicit_odd_matches_previous_dominant_correlations_and_records_sector():
    original = RVBTransfer("kagome", 2)
    explicit = RVBTransfer("kagome", 2, sector="odd")
    np.testing.assert_allclose(explicit.correlations(8).connected,
                               original.correlations(8).connected, atol=1e-13)
    result = rvb_cylinder_benchmark("kagome", 2, 8, sector="odd")
    assert result["metadata"]["virtual_vacuum_parity"] == "odd"
    assert result["metadata"]["sector_requested"] == "odd"
    assert result["metadata"]["vison_inserted"] is False


def test_mixed_fixed_point_is_rejected_instead_of_assigning_a_sector():
    model = copy.copy(RVBTransfer("kagome", 2))
    ket, bra = model._vacuum_parities(model.active)
    opposite = np.flatnonzero((ket == 0) & (bra == 0))[0]
    model.right = model.right.copy()
    model.right[opposite] += 1
    with pytest.raises(RuntimeError, match="mixes virtual parity"):
        model._identify_kagome_sector()


@pytest.mark.parametrize("lattice,sector", [
    ("kagome", "mixed"), ("square", "even"), ("square", "odd"),
])
def test_unsupported_sector_selection_fails(lattice, sector):
    with pytest.raises(ValueError):
        RVBTransfer(lattice, 2, sector=sector)
