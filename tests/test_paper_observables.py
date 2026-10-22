"""Independent physical checks of the operators used in paper rechecks."""
import numpy as np
import pytest

from scripts.run_paper_kagome import physical_bond_operator


@pytest.mark.parametrize("pair", [(0, 1), (0, 2), (1, 2)])
def test_singlet_has_correct_bond_observables(pair):
    # Build |singlet_ij> |up_k> directly in the computational basis.
    psi = np.zeros((2, 2, 2))
    up_down = [0, 0, 0]
    up_down[pair[1]] = 1
    down_up = [0, 0, 0]
    down_up[pair[0]] = 1
    psi[tuple(up_down)] = 1 / np.sqrt(2)
    psi[tuple(down_up)] = -1 / np.sqrt(2)
    psi = psi.ravel()
    zz = physical_bond_operator(pair, "zz", flagged=False)
    ss = physical_bond_operator(pair, "ss", flagged=False)
    assert np.allclose(zz @ psi, -psi)
    assert np.allclose(ss @ psi, -.75 * psi)
    # A symmetric pair has scalar eigenvalue +1/4.
    triplet = np.zeros(8)
    triplet[0] = 1
    assert np.allclose(ss @ triplet, .25 * triplet)


@pytest.mark.parametrize("kind", ["zz", "ss"])
def test_flagged_operator_is_identity_on_every_flag(kind):
    unflagged = physical_bond_operator((0, 2), kind, flagged=False)
    flagged = physical_bond_operator((0, 2), kind, flagged=True)
    # All eight flag assignments independently give the identical spin operator;
    # no operator matrix element changes flags.
    matrix = flagged.reshape((2, 2) * 6)
    for flags in np.ndindex(2, 2, 2):
        indices = (flags[0], slice(None), flags[1], slice(None), flags[2], slice(None))
        block = matrix[indices + indices].reshape(8, 8)
        assert np.allclose(block, unflagged)
    assert np.count_nonzero(flagged) == 8 * np.count_nonzero(unflagged)


def test_operator_rejects_ambiguous_name_or_pair():
    with pytest.raises(ValueError):
        physical_bond_operator((0, 1), "dimer", False)
    with pytest.raises(ValueError):
        physical_bond_operator((0, 0), "zz", False)
