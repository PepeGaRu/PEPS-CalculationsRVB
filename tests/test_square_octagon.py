"""Independent matching enumeration checks the arrow construction."""

from collections import Counter
from fractions import Fraction

import pytest

from rvb.exact import occupancy_statistics, perfect_matchings
from rvb.square_octagon import (SquareOctagonOrthogonalDimers, square_octagon,
                               square_octagon_geometry, square_octagon_orthogonal_benchmark,
                               dummy_contracted_square_geometry)


@pytest.mark.parametrize("nx,ny", [(1, 1), (1, 2), (2, 1), (2, 2)])
def test_geometry_and_independent_covering_bijection(nx, ny):
    model = SquareOctagonOrthogonalDimers(nx, ny)
    lattice = model.lattice
    assert lattice.nsites == 8*nx*ny
    assert len(lattice.edges) == 14*nx*ny
    assert len(set(map(frozenset, lattice.edges))) == len(lattice.edges)
    independent = perfect_matchings(lattice)
    actual = model.coverings()
    assert len(independent) == 2**(2*nx*ny+1)
    assert {frozenset(c) for c in actual} == {frozenset(c) for c in independent}
    for covering in independent:
        bits = model.bits_from_covering(covering)
        assert set(model.covering_from_bits(bits)) == set(covering)
        pattern = model.square_loop_pattern(bits)
        for x in range(nx):
            for y in range(ny):
                occupied = [((bits >> model.dummy_variable(x, y, k)) & 1) for k in range(4)]
                assert sum(occupied) % 2 == 0
        base_loops = model.base_loop_pattern(bits)
        assert all(sum(base_loops[i] for i in range(model.nvariables) if (mask >> i) & 1) % 2 == 0
                   for mask, _ in model.equations)


@pytest.mark.parametrize("nx,ny", [(1, 1), (1, 2), (2, 2)])
def test_exact_marginals_and_joint_probabilities_match_independent_enumeration(nx, ny):
    model = SquareOctagonOrthogonalDimers(nx, ny)
    mean, covariance = occupancy_statistics(model.lattice)
    for a in range(len(model.edges)):
        assert float(model.occupation(a)) == pytest.approx(mean[a], abs=1e-14)
        for b in range(len(model.edges)):
            assert float(model.connected(a, b)) == pytest.approx(covariance[a, b], abs=1e-14)


def test_projection_has_two_resolutions_per_original_square_vertex():
    model = SquareOctagonOrthogonalDimers(2, 2)
    counts = Counter(model.square_loop_pattern(model.bits_from_parameters(p))
                     for p in range(model.covering_count))
    assert len(counts) == 2**(model.ncells+1)
    assert set(counts.values()) == {2**model.ncells}


@pytest.mark.parametrize("nx,ny", [(1, 1), (1, 2), (2, 2)])
def test_dummy_contraction_bijects_coverings_with_all_even_subgraphs(nx, ny):
    model = SquareOctagonOrthogonalDimers(nx, ny)
    graph = dummy_contracted_square_geometry(nx, ny).lattice
    assert graph.nsites == 2*nx*ny
    assert len(graph.edges) == 4*nx*ny
    incidence = [[e for e, endpoints in enumerate(graph.edges) if vertex in endpoints]
                 for vertex in range(graph.nsites)]
    assert all(len(edges) == 4 for edges in incidence)
    # Independent edge-subset enumeration includes every even graph, not just
    # projections generated from the matching-constraint solver.
    even_subgraphs = set()
    for integer in range(1 << len(graph.edges)):
        pattern = tuple((integer >> e) & 1 for e in range(len(graph.edges)))
        if all(sum(pattern[e] for e in edges) % 2 == 0 for edges in incidence):
            even_subgraphs.add(pattern)
    reference = model.bits_from_parameters(model.covering_count-1)
    counts = Counter()
    for p in range(model.covering_count):
        bits = model.bits_from_parameters(p)
        pattern = model.dummy_contracted_loop_pattern(bits, reference)
        counts[pattern] += 1
        assert model.lift_dummy_contracted_loop_pattern(pattern) == model.base_loop_pattern(bits, reference)
        assert model.bits_from_dummy_contracted_loops(pattern, reference) == bits
    assert set(counts) == even_subgraphs
    assert set(counts.values()) == {1}
    assert len(counts) == 2**(2*nx*ny+1)


def test_dummy_contraction_retains_the_fixed_four_leg_routing_and_small_loops():
    model = SquareOctagonOrthogonalDimers(3, 3)
    geometry = dummy_contracted_square_geometry(3, 3)
    pattern = (1,) * (4*model.ncells)
    lifted = model.lift_dummy_contracted_loop_pattern(pattern)
    assert lifted[4*model.ncells:] == (0,) * (2*model.ncells)
    routing = model.dummy_contracted_routing(pattern)
    # All small-square perimeter loops survive. A crossing or different
    # four-leg resolution would join these loops and fail this check.
    adjacency = {edge: set() for edge in range(len(pattern))}
    for v, pairs in enumerate(routing):
        assert set(pairs) == {tuple(sorted(pair)) for pair in geometry.endpoint_ring_edges[v]}
        for a, b in pairs:
            adjacency[a].add(b)
            adjacency[b].add(a)
    unseen = set(adjacency)
    components = []
    while unseen:
        stack = [unseen.pop()]
        component = set(stack)
        while stack:
            for neighbor in adjacency[stack.pop()]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    assert {frozenset(c) for c in components} == {
        frozenset(range(4*c, 4*c+4)) for c in range(model.ncells)}


def test_dummy_contracted_positions_form_the_rotated_square_and_reject_open_ends():
    model = SquareOctagonOrthogonalDimers(3, 3)
    geometry = dummy_contracted_square_geometry(3, 3)
    graph = geometry.lattice
    # Interior small square centered at (1,1); after rotation its perimeter
    # is a unit square and its four nodes are exact dummy-link midpoints.
    c = 1*3+1
    positions = [graph.positions[a] for a, _ in graph.edges[4*c:4*c+4]]
    assert positions == [(1.5, 1.), (1., 1.5), (.5, 1.), (1., .5)]
    rotated = [(x+y-.5, y-x+.5) for x, y in positions]
    assert rotated == [(2., 0.), (2., 1.), (1., 1.), (1., 0.)]
    for invalid in ((), (0,)*35, (2,)*36, (1,)+(0,)*35):
        with pytest.raises(ValueError):
            model.lift_dummy_contracted_loop_pattern(invalid)
    for dimensions in ((0, 2), (2, -1), (True, 2)):
        with pytest.raises(ValueError):
            dummy_contracted_square_geometry(*dimensions)


def test_dummy_empty_means_a_singlet_and_never_an_unpaired_spin():
    model = SquareOctagonOrthogonalDimers(1, 2)
    for covering in model.coverings():
        bits = model.bits_from_covering(covering)
        selected = set(covering)
        for c, links in enumerate(model.geometry.dummy_edges):
            for axis, edge in enumerate(links):
                t = (bits >> (4*model.ncells+2*c+axis)) & 1
                assert (edge in selected) == (t == 0)
        assert len(covering)*2 == model.lattice.nsites


def test_plane_like_torus_has_strictly_local_occupation_correlations():
    model = SquareOctagonOrthogonalDimers(8, 6)
    for e, kind in enumerate(model.geometry.edge_kinds):
        assert model.occupation(e) == (Fraction(1, 4) if kind == "triangle" else Fraction(1, 2))
    # All fourteen cell-zero bond types against every bond belonging to a
    # translated cell at least two cells away. Fractions test exact zeros.
    reference = [e for triangle in model.geometry.triangle_edges[0] for e in triangle]
    reference += list(model.geometry.dummy_edges[0])
    for x in range(model.nx):
        for y in range(model.ny):
            dx, dy = min(x, model.nx-x), min(y, model.ny-y)
            if dx + dy <= 1:
                continue
            c = x*model.ny+y
            targets = [e for triangle in model.geometry.triangle_edges[c] for e in triangle]
            targets += list(model.geometry.dummy_edges[c])
            assert all(model.connected(a, b) == 0 for a in reference for b in targets)


def test_uniform_loop_soup_does_not_get_a_spurious_power_law_fit():
    report = square_octagon_orthogonal_benchmark(8, 4)
    assert report["covering_count_log2"] == 65
    assert report["square_loop_preimage_count_log2"] == 32
    for family in report["families"]:
        assert family["connected_correlation"][0] > 0
        assert family["connected_correlation"][2:] == [0., 0., 0.]
    assert "correlation_length" not in report


def test_invalid_geometries_and_assignments_are_rejected():
    for dimensions in ((0, 2), (-1, 2), (1.5, 2), (True, 2)):
        with pytest.raises(ValueError):
            square_octagon(*dimensions)
    model = SquareOctagonOrthogonalDimers(1, 1)
    with pytest.raises(ValueError):
        model.covering_from_bits(0)
    with pytest.raises(ValueError):
        model.bits_from_covering(())
    with pytest.raises(ValueError):
        model.coverings(max_coverings=1)
    with pytest.raises(ValueError):
        square_octagon_orthogonal_benchmark(4, 4, max_distance=3)
