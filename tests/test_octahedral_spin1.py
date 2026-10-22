"""Independent spin algebra and dense checks of the spin-one covering state."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product

import numpy as np
import pytest

from rvb.octahedral import OctahedralCoverings
from rvb.octahedral_spin1 import (
    cube_contraction,
    cube_geometry,
    Spin1OctahedralTwoCube,
    spin1_dimer,
    spin1_generators,
    spin1_trimer,
    octahedron_tensor,
    ownership_projection,
)


def independent_generators():
    """Vector-representation spin matrices, written without the model tensors."""
    return np.array([
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
        [[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]],
        [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]],
    ], dtype=complex)


def embed_one(operator, site, nsites):
    result = np.array([[1.]])
    for i in range(nsites):
        result = np.kron(result, operator if i == site else np.eye(3))
    return result


def test_spin_one_singlets_are_normalized_total_spin_zero():
    generators = independent_generators()
    np.testing.assert_allclose(spin1_generators(), generators, atol=1e-14)
    np.testing.assert_allclose(sum(s @ s for s in generators), 2*np.eye(3), atol=1e-14)
    for nsites, tensor in ((2, spin1_dimer()), (3, spin1_trimer())):
        state = tensor.ravel()
        assert np.vdot(state, state) == pytest.approx(1, abs=1e-14)
        for generator in generators:
            total = sum(embed_one(generator, site, nsites) for site in range(nsites))
            np.testing.assert_allclose(total @ state, 0, atol=1e-14)
    np.testing.assert_allclose(spin1_dimer().T, spin1_dimer(), atol=1e-14)
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        np.testing.assert_allclose(spin1_trimer().transpose(perm), sign*spin1_trimer(), atol=1e-14)


def test_singlet_projectors_agree_with_physical_spin_polynomials():
    generators = independent_generators()
    bond = sum(np.kron(s, s) for s in generators)
    dimer_projector = np.outer(spin1_dimer().ravel(), spin1_dimer().ravel())
    np.testing.assert_allclose(dimer_projector, (bond @ bond-np.eye(9))/3, atol=1e-14)
    total_spin_squared = sum(
        np.linalg.matrix_power(sum(embed_one(s, site, 3) for site in range(3)), 2)
        for s in generators)
    identity = np.eye(27)
    trimer_projector = np.outer(spin1_trimer().ravel(), spin1_trimer().ravel())
    polynomial = (identity-total_spin_squared/2) @ (identity-total_spin_squared/6) @ (identity-total_spin_squared/12)
    np.testing.assert_allclose(trimer_projector, polynomial, atol=1e-13)
    assert np.linalg.matrix_rank(trimer_projector, tol=1e-12) == 1
    # A pair within the three-spin singlet is spin one, never a pair singlet.
    trimer = spin1_trimer().ravel()
    assert np.vdot(trimer, np.kron(dimer_projector, np.eye(3)) @ trimer) == pytest.approx(0, abs=1e-14)
    assert np.vdot(trimer, np.kron(bond, np.eye(3)) @ trimer) == pytest.approx(-1, abs=1e-14)


def test_distinct_coverings_have_signed_nonzero_overlaps():
    dimer = spin1_dimer()
    state = np.einsum('ab,cd->abcd', dimer, dimer)
    alternative = np.einsum('ac,bd->abcd', dimer, dimer)
    assert np.vdot(state, alternative) == pytest.approx(1/3, abs=1e-14)
    trimers = np.einsum('abc,def->abcdef', spin1_trimer(), spin1_trimer())
    assert np.vdot(trimers, trimers.transpose(0, 1, 3, 2, 4, 5)) == pytest.approx(1/3, abs=1e-14)
    assert np.vdot(trimers, trimers.transpose(0, 3, 5, 1, 2, 4)) == pytest.approx(-1/3, abs=1e-14)


def test_period_two_coverings_split_into_two_special_families():
    model = OctahedralCoverings(2, 2, 2)
    counts = Counter()
    dimer_pairings = set()
    trimer_assignments = set()
    for parameter in range(model.covering_count):
        bits = model.bits_from_parameters(parameter)
        covering = model.covering(bits)
        kinds = {obj.kind for obj in covering}
        assert len(kinds) == 1
        kind = kinds.pop()
        counts[kind] += 1
        if kind == 'dimer':
            dimer_pairings.add(tuple(sorted(tuple(sorted(obj.sites)) for obj in covering)))
        else:
            trimer_assignments.add(tuple(sorted((obj.cell, obj.sites) for obj in covering)))
    assert counts == {'dimer': 4096, 'trimer': 4096}
    assert len(dimer_pairings) == 1
    assert len(trimer_assignments) == 4096


def independent_cube():
    cells = tuple(product((0, 1), repeat=3))
    ports = np.zeros((8, 3), dtype=int)
    edges = []
    for vertex, coord in enumerate(cells):
        for axis in range(3):
            if coord[axis] == 0:
                other = list(coord)
                other[axis] = 1
                target = cells.index(tuple(other))
                ports[vertex, axis] = ports[target, axis] = len(edges)
                edges.append((vertex, target, axis))
    return cells, tuple(edges), ports


def test_cube_cross_overlap_by_independent_edge_color_enumeration():
    _, _, ports = independent_cube()
    # Enumerate the 3^12 edge colors directly; no tensor contraction is used.
    colors = np.array(np.unravel_index(np.arange(3**12), (3,)*12)).T
    signs = np.ones(len(colors), dtype=int)
    for local in ports:
        a, b, c = colors[:, local].T
        signs *= (a-b)*(b-c)*(c-a)//2
    assert np.count_nonzero(signs) == int(signs.sum()) == 24
    direct = signs.sum()/(3**6*6**4)
    contracted = cube_contraction([spin1_trimer()]*8, [spin1_dimer()]*12)
    assert contracted == pytest.approx(direct, abs=1e-16)
    assert direct == pytest.approx(1/39366, abs=1e-16)


def double_vertex():
    trimer = spin1_trimer()
    return np.einsum('abc,def->adbecf', trimer, trimer).reshape(9, 9, 9)


def reshuffled_pair_operator(operator):
    # Ordinary two-site operator indices are (a,b; a',b'); the vertex legs
    # group a with a' and b with b'.
    return operator.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)


@pytest.mark.parametrize('include_identity', [False, True])
def test_cube_subgraph_against_dense_twelve_spin_swaps(include_identity):
    cells, edges, _ = independent_cube()
    selected_vertices = [i for i, c in enumerate(cells) if c[2] == 0]
    selected_edges = [e for e, (v, w, _) in enumerate(edges)
                      if v in selected_vertices and w in selected_vertices]
    assert len(selected_vertices) == len(selected_edges) == 4
    swap = np.eye(9).reshape(3, 3, 3, 3).transpose(1, 0, 2, 3).reshape(9, 9)
    active_operator = swap + np.eye(9) if include_identity else swap
    own_tensors = [reshuffled_pair_operator(active_operator if e in selected_edges else np.eye(9))
                   for e in range(12)]
    _, implementation_edges = cube_geometry()
    edge_tensors = [own_tensors[edges.index((u, v, a))] for a, u, v in implementation_edges]
    contracted = cube_contraction([double_vertex()]*8, edge_tensors)
    # A 4-vertex cycle involves just four product trimers. Form its entire
    # 3^12-dimensional state and apply swaps as explicit physical-site
    # permutations, independently of the cube contraction and reshuffling.
    state = spin1_trimer()
    for _ in range(3):
        state = np.multiply.outer(state, spin1_trimer())
    transformed = state.copy()
    for e in selected_edges:
        v, w, axis = edges[e]
        first = 3*selected_vertices.index(v)+axis
        second = 3*selected_vertices.index(w)+axis
        swapped = transformed.swapaxes(first, second)
        transformed = transformed+swapped if include_identity else swapped
    expected = np.vdot(state, transformed)
    assert contracted == pytest.approx(expected, abs=2e-13)


def test_all_pair_swaps_factor_into_two_epsilon_cubes():
    swap = np.eye(9).reshape(3, 3, 3, 3).transpose(1, 0, 2, 3).reshape(9, 9)
    contracted = cube_contraction([double_vertex()]*8, [reshuffled_pair_operator(swap)]*12)
    # Swapping every physical pair connects bra-even to ket-odd and vice
    # versa. Each of the two components has 24 positive edge colorings.
    assert contracted == pytest.approx(24**2/6**8, abs=1e-15)


def test_dense_physical_pair_projectors_on_four_trimers():
    cells, edges, _ = independent_cube()
    selected_vertices = [i for i, c in enumerate(cells) if c[2] == 0]
    selected_edges = [e for e, (v, w, _) in enumerate(edges)
                      if v in selected_vertices and w in selected_vertices]
    singlet = spin1_dimer().ravel()
    projector = np.outer(singlet, singlet)
    own_tensors = [reshuffled_pair_operator(projector if e in selected_edges else np.eye(9))
                   for e in range(12)]
    _, implementation_edges = cube_geometry()
    edge_tensors = [own_tensors[edges.index((u, v, a))] for a, u, v in implementation_edges]
    contracted = cube_contraction([double_vertex()]*8, edge_tensors)
    state = spin1_trimer()
    for _ in range(3):
        state = np.multiply.outer(state, spin1_trimer())
    transformed = state.copy()
    for e in selected_edges:
        v, w, axis = edges[e]
        first = 3*selected_vertices.index(v)+axis
        second = 3*selected_vertices.index(w)+axis
        moved = np.moveaxis(transformed, (first, second), (-2, -1))
        shape = moved.shape
        acted = (moved.reshape(-1, 9) @ projector.T).reshape(shape)
        transformed = np.moveaxis(acted, (-2, -1), (first, second))
    assert contracted == pytest.approx(np.vdot(state, transformed), abs=1e-14)


def test_projecting_all_diameter_pairs_gives_dimer_state_fidelity():
    model = Spin1OctahedralTwoCube()
    parts = model.norm_components()
    assert parts['reference_overlap'] == Fraction(1, 39366)
    assert parts['total'] > parts['dimer_part'] > 0
    assert model.dimer_projector_expectation(()) == 1
    # The twelve disjoint rank-one projectors have joint range exactly |D>.
    fidelity = (4096+parts['dimer_trimer_sum_overlap'])**2/parts['total']
    assert model.dimer_projector_expectation(tuple(range(12))) == fidelity
    singles = [model.dimer_projector_expectation((edge,)) for edge in range(12)]
    assert len(set(singles)) == 1
    assert 0 < fidelity <= singles[0] < 1
    assert model.dimer_projector_covariance(0, 0) == singles[0]*(1-singles[0])


def test_trimer_network_empty_insertions_reproduce_independent_norm_network():
    model = Spin1OctahedralTwoCube()
    parts = model.norm_components()
    assert model._trimer_projector_branch('DD', ()) == 1
    assert model._trimer_projector_branch('DT', ()) == parts['dimer_trimer_sum_overlap']
    assert model._trimer_projector_branch('TT', ()) == parts['trimer_part']
    assert model.trimer_projector_expectation(()) == 1


def test_trimer_projectors_on_the_dimer_product_have_maximally_mixed_marginals():
    model = Spin1OctahedralTwoCube()
    assert model._trimer_projector_branch('DD', (0,)) == Fraction(1, 27)
    # Distinct reference faces involve disjoint spins. Even when their
    # centres share one diameter, tracing the other legs gives I/3 on each
    # end of that diameter and the connected expectation is zero.
    for target in (1, 3, 7):
        assert model._trimer_projector_branch('DD', (0, target)) == Fraction(1, 729)
    mean = model.trimer_projector_expectation((0,))
    assert 0 < mean < 1
    assert model.trimer_projector_covariance(0, 0) == mean*(1-mean)


def test_general_octahedron_tensor_has_equal_weight_singlet_supports():
    tensor = octahedron_tensor()
    assert tensor.shape == (4,)*6
    assert np.vdot(tensor, tensor) == pytest.approx(16, abs=1e-13)
    # Check all occupation supports directly against the three pair parities.
    # Each allowed support must carry one normalized spin state.
    for mask in range(64):
        indices = tuple(slice(1, 4) if (mask >> p) & 1 else 0 for p in range(6))
        local = tensor[indices]
        parities = {((mask >> (2*a)) ^ (mask >> (2*a+1))) & 1 for a in range(3)}
        expected_norm = int(len(parities) == 1)
        assert np.vdot(local, local) == pytest.approx(expected_norm, abs=1e-14)
    for physical_generator in independent_generators():
        virtual_generator = np.zeros((4, 4), dtype=complex)
        virtual_generator[1:, 1:] = physical_generator
        total_action = np.zeros_like(tensor, dtype=complex)
        for leg in range(6):
            acted = np.tensordot(virtual_generator, tensor, axes=(1, leg))
            total_action += np.moveaxis(acted, 0, leg)
        np.testing.assert_allclose(total_action, 0, atol=1e-14)


def test_general_ownership_map_intertwines_spin_rotations():
    projection = ownership_projection().reshape(3, 16)
    # Both choices of owner map with unit amplitude to the same physical spin.
    np.testing.assert_allclose(projection @ projection.T, 2*np.eye(3), atol=1e-14)
    for physical_generator in independent_generators():
        virtual_generator = np.zeros((4, 4), dtype=complex)
        virtual_generator[1:, 1:] = physical_generator
        incoming_generator = (np.kron(virtual_generator, np.eye(4))
                              + np.kron(np.eye(4), virtual_generator))
        np.testing.assert_allclose(physical_generator @ projection,
                                   projection @ incoming_generator, atol=1e-14)
