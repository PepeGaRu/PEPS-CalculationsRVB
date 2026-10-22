"""Finite covering-oracle checks of the sector-resolved cylinder machinery."""

import numpy as np
import pytest

from rvb.exact import bond_expectation, perfect_matchings, rvb_state
from rvb.octagon_transfer import (OctagonCylinder, boundary_vacuum_parities,
                                  _full_norm_transfer)
from rvb.square_octagon import square_octagon_geometry


@pytest.mark.parametrize("width", [1, 2])
@pytest.mark.parametrize("sector", ["even", "odd"])
def test_periodic_sector_norm_and_bond_insertion_against_covering_oracle(width, sector):
    # One x cell is closed into a torus. The physical covering itself tells us
    # whether each dummy edge crossing its x cut is absent (=virtual vacuum).
    geometry = square_octagon_geometry(1, width)
    coverings = perfect_matchings(geometry.lattice)
    east_edges = [pair[0] for pair in geometry.dummy_edges]
    parity = int(sector == "odd")
    selected = [c for c in coverings if sum(edge not in c for edge in east_edges) % 2 == parity]
    oracle, norm = rvb_state(geometry.lattice, selected)
    model = OctagonCylinder(width, sector)
    partition = model.transfer.diagonal().sum()
    assert partition == pytest.approx(norm * 2**(geometry.lattice.nsites // 2), abs=2e-11)
    for name, edge in (("ring_ring", (0, 3)), ("ring_dummy_0", (0, 4))):
        mean = model.observable_transfer(name).diagonal().sum() / partition
        assert mean == pytest.approx(bond_expectation(oracle, edge), abs=3e-14)


@pytest.mark.parametrize("width", [2, 3])
def test_norm_and_physical_insertions_preserve_each_virtual_parity(width):
    full = _full_norm_transfer(width)
    ket, bra = boundary_vacuum_parities(np.arange(full.shape[0]), width)
    rows, cols = full.nonzero()
    assert np.all(ket[rows] == ket[cols])
    assert np.all(bra[rows] == bra[cols])
    # Independent reference parity: decode the base-3 ket/bra boundary states.
    for index in range(9**width):
        digits = [(index // 9**j) % 9 for j in range(width)]
        assert ket[index] == sum(d // 3 == 2 for d in digits) % 2
        assert bra[index] == sum(d % 3 == 2 for d in digits) % 2


@pytest.mark.parametrize("sector", ["even", "odd"])
def test_spectral_channel_reconstructs_separated_physical_correlations(sector):
    model = OctagonCylinder(2, sector)
    assert model.eigen_residual < 1e-12
    for name in ("ring_ring", "ring_dummy_0", "spin_ring_z"):
        _, connected = model.correlations(name, max_distance=8)
        spectrum = model.channel_spectrum(name)
        reconstructed = np.zeros(8, dtype=complex)
        for mode in spectrum["coupled_modes"]:
            ratio = complex(mode["eigenvalue_ratio_real"], mode["eigenvalue_ratio_imag"])
            weight = complex(mode["spectral_weight_real"], mode["spectral_weight_imag"])
            reconstructed += weight * ratio**np.arange(8)
        # Exact zero modes may contribute to contact distance1, hence compare2+.
        np.testing.assert_allclose(reconstructed[1:], connected[1:], atol=5e-12)
        leading = spectrum["leading_coupled_mode"]
        assert leading["right_mode_residual"] < 1e-12
        assert leading["left_mode_residual"] < 1e-12
    assert abs(model.correlations("spin_ring_z")[0]) < 1e-13


def test_rejects_uncontrolled_sector_and_excessive_width():
    for width, sector in ((5, "even"), (0, "even"), (2, "dominant")):
        with pytest.raises(ValueError):
            OctagonCylinder(width, sector)
