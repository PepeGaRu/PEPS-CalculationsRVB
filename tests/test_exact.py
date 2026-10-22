import numpy as np
import pytest
from rvb.lattices import square, kagome
from rvb.exact import perfect_matchings, rvb_state, peps_state, bond_expectation, connected_bond_correlation, occupancy_statistics


def test_square_known_singlet_and_orthogonal_results():
    lattice = square()
    matchings = perfect_matchings(lattice)
    assert len(matchings) == 2
    state, norm = rvb_state(lattice, matchings)
    assert norm == pytest.approx(3.)
    assert bond_expectation(state, (0, 1)) == pytest.approx(-.5)
    assert bond_expectation(state, (0, 3)) == pytest.approx(.25)
    assert connected_bond_correlation(state, (0, 1), (2, 3)) == pytest.approx(3/16)
    mean, connected = occupancy_statistics(lattice, matchings)
    assert np.allclose(mean, .5)
    e0 = next(i for i,e in enumerate(lattice.edges) if set(e)=={0,1})
    e1 = next(i for i,e in enumerate(lattice.edges) if set(e)=={2,3})
    assert connected[e0,e1] == pytest.approx(.25)


@pytest.mark.parametrize("lattice", [square(), kagome(2,2), kagome(2,2,True,True)])
def test_peps_matches_independent_singlet_coverings(lattice):
    oracle, norm = rvb_state(lattice)
    state, peps_norm = peps_state(lattice)
    assert peps_norm == pytest.approx(norm, abs=1e-11)
    assert abs(oracle @ state) == pytest.approx(1., abs=1e-11)
    total_s2 = .75*lattice.nsites + 2*sum(bond_expectation(state,(a,b)) for a in range(lattice.nsites) for b in range(a+1,lattice.nsites))
    assert abs(total_s2) < 1e-10


def test_kagome_geometry_and_exact_covering_count():
    lattice = kagome(2,2,True,True)
    assert lattice.nsites == 12
    assert len(lattice.edges) == 24
    assert all(sum(i in e for e in lattice.edges)==4 for i in range(12))
    assert len(perfect_matchings(lattice)) == 32


def test_oracles_refuse_unphysical_or_excessive_requests():
    with pytest.raises(ValueError): square(3,4,True,False)
    with pytest.raises(ValueError): perfect_matchings(square(), max_coverings=1)
    with pytest.raises(ValueError): rvb_state(square(6,4))
    state, _ = rvb_state(square())
    with pytest.raises(ValueError): connected_bond_correlation(state,(0,1),(1,3))
