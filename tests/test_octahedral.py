"""Independent covering enumeration and exact checks for the 16-state model."""
from collections import Counter
from fractions import Fraction
from itertools import product

import pytest

from rvb.octahedral import LOCAL_STATES, OctahedralCoverings


def independent_geometry(shape):
    """Construct physical site sets directly, without using the parity code."""
    nx, ny, nz = shape
    cells = tuple(product(range(nx), range(ny), range(nz)))

    def index(x, y, z):
        return ((x % nx)*ny + y % ny)*nz + z % nz

    ports = []
    objects = []
    for owner, coord in enumerate(cells):
        local = []
        for axis in range(3):
            previous = list(coord)
            previous[axis] -= 1
            local.extend((3*owner+axis, 3*index(*previous)+axis))
        ports.append(tuple(local))
        for axis in range(3):
            objects.append((owner, "dimer", tuple(local[2*axis:2*axis+2])))
        for picks in product((0, 1), repeat=3):
            objects.append((owner, "trimer", tuple(local[2*a+picks[a]] for a in range(3))))
    return cells, tuple(ports), tuple(objects)


def support_from_ownership(bits, ports):
    return sum((((bits >> site) & 1) ^ (p % 2 == 0)) << p
               for p, site in enumerate(ports))


@pytest.fixture(scope="module")
def enumerated_two_cube():
    """Enumerate hypergraph exact covers, not solutions of arrow equations.

    Size-two chords with identical endpoint sets remain distinct objects when
    their owning octahedra differ. Choosing the first uncovered site avoids
    permutation overcounting. The trimer-owner mask excludes the forbidden
    complementary-trimer realization of a fully occupied octahedron.
    """
    cells, ports, objects = independent_geometry((2, 2, 2))
    nsites = 3*len(cells)
    candidates = [[] for _ in range(nsites)]
    for owner, kind, sites in objects:
        occupied = sum(1 << site for site in sites)
        # Each site is a directed cubic bond; its tail has bit 0 ownership.
        ownership = sum(1 << site for site in sites if site // 3 != owner)
        trimer_owner = (1 << owner) if kind == "trimer" else 0
        entry = occupied, ownership, trimer_owner
        for site in sites:
            candidates[site].append(entry)
    counts = Counter()
    full = (1 << nsites)-1

    def visit(occupied, ownership, trimer_owners):
        remaining = full ^ occupied
        if not remaining:
            counts[ownership] += 1
            return
        first = (remaining & -remaining).bit_length()-1
        for sites, arrows, trimer_owner in candidates[first]:
            if occupied & sites or trimer_owners & trimer_owner:
                continue
            visit(occupied | sites, ownership | arrows, trimer_owners | trimer_owner)

    visit(0, 0, 0)
    return counts, ports, objects


def test_catalog_is_exactly_the_sixteen_allowed_supports():
    supports = {state.support for state in LOCAL_STATES}
    parity_allowed = {
        mask for mask in range(64)
        if len({((mask >> (2*a)) & 1) ^ ((mask >> (2*a+1)) & 1)
                for a in range(3)}) == 1}
    assert supports == parity_allowed
    assert len(LOCAL_STATES) == len(supports) == 16
    assert sum(state.is_trimer for state in LOCAL_STATES) == 8
    assert Counter(len(s.dimer_axes) for s in LOCAL_STATES if not s.is_trimer) == {
        0: 1, 1: 3, 2: 3, 3: 1}
    full = next(state for state in LOCAL_STATES if state.support == 63)
    assert full.dimer_axes == (0, 1, 2)
    assert not full.is_trimer
    for state in LOCAL_STATES:
        if state.is_trimer:
            assert all(sum(state.incoming[2*a:2*a+2]) == 1 for a in range(3))


def test_arrow_code_equals_independent_hypergraph_coverings(enumerated_two_cube):
    counts, _, _ = enumerated_two_cube
    model = OctahedralCoverings(2, 2, 2)
    assert sum(counts.values()) == len(counts) == 8192
    assert set(counts.values()) == {1}
    actual = {model.bits_from_parameters(i) for i in range(model.covering_count)}
    assert actual == set(counts)


def test_size_two_keeps_owner_labels_for_coincident_site_sets(enumerated_two_cube):
    _, _, objects = enumerated_two_cube
    multiplicity = Counter((kind, tuple(sorted(sites))) for _, kind, sites in objects)
    assert any(number == 2 for (kind, _), number in multiplicity.items() if kind == "dimer")
    # Losing those owner labels would merge physically different periodic chords.
    assert len(objects) == 8*(3+8)


@pytest.mark.parametrize("target", [1, 4, 7])
def test_joint_probabilities_against_full_cover_enumeration(enumerated_two_cube, target):
    counts, ports, _ = enumerated_two_cube
    model = OctahedralCoverings(2, 2, 2)
    joint = Counter((support_from_ownership(bits, ports[0]),
                     support_from_ownership(bits, ports[target])) for bits in counts)
    for first, second in product(LOCAL_STATES, repeat=2):
        expected = Fraction(joint[first.support, second.support], len(counts))
        assert model.joint_state_probability(0, first, target, second) == expected


@pytest.mark.parametrize("shape", [(2, 2, 2), (3, 3, 3), (4, 4, 4), (3, 4, 5)])
def test_constraint_and_cube_move_ranks(shape):
    model = OctahedralCoverings(*shape)
    n, lengths = model.ncells, sum(shape)
    assert model.nsites == 3*n
    assert model.rank == 2*n-lengths+1
    assert model.dimension == n+lengths-1
    assert model.covering_count == 2**(n+lengths-1)
    assert model.cube_move_rank == n-lengths+2
    assert model.dimension-model.cube_move_rank == 2*lengths-3


def test_sampled_coverings_are_actual_disjoint_local_objects():
    model = OctahedralCoverings(3, 4, 5)
    _, ports, objects = independent_geometry(model.shape)
    allowed = {(owner, kind, frozenset(sites)) for owner, kind, sites in objects}
    for bits in (model.all_minus_trimer_reference, 0, *model.sample(20, seed=739)):
        covering = model.covering(bits)
        site_counts = Counter(site for obj in covering for site in obj.sites)
        assert site_counts == Counter({site: 1 for site in range(model.nsites)})
        assert all((obj.cell, obj.kind, frozenset(obj.sites)) in allowed for obj in covering)
        trimer_counts = Counter(obj.cell for obj in covering if obj.kind == "trimer")
        assert all(number == 1 for number in trimer_counts.values())
        for cell, local_ports in enumerate(ports):
            owned_sites = {site for obj in covering if obj.cell == cell for site in obj.sites}
            observed = sum((site in owned_sites) << p for p, site in enumerate(local_ports))
            assert observed == support_from_ownership(bits, local_ports)
    reference = model.covering(model.all_minus_trimer_reference)
    assert len(reference) == model.ncells
    assert all(obj.kind == "trimer" and obj.face_signs == (-1, -1, -1) for obj in reference)


def test_cube_moves_flip_twelve_edges_and_preserve_coverings():
    model = OctahedralCoverings(3, 4, 5)
    sample = model.sample(1, seed=147)[0]
    for origin in model.cells:
        move = model.cube_flip_mask(origin)
        assert move.bit_count() == 12
        degree = Counter(endpoint for edge, ends in enumerate(model.edge_endpoints)
                         if (move >> edge) & 1 for endpoint in ends)
        assert len(degree) == 8
        assert set(degree.values()) == {3}
        assert model.valid(move)
        assert model.valid(sample ^ move)
        assert model.cage_net(sample ^ move, sample) == move


def test_local_means_and_number_of_covered_sites():
    model = OctahedralCoverings(4, 4, 4)
    for state in LOCAL_STATES:
        assert model.state_probability((0, 0, 0), state) == Fraction(1, 16)
    assert [model.dimer_probability(0, a) for a in range(3)] == [Fraction(1, 4)]*3
    assert model.trimer_probability(0) == Fraction(1, 2)
    for signs in product((-1, 1), repeat=3):
        assert model.face_probability(0, signs) == Fraction(1, 16)
    average_dimers = sum(len(s.dimer_axes)*model.state_probability(0, s) for s in LOCAL_STATES)
    assert average_dimers == Fraction(3, 4)
    assert 2*average_dimers+3*model.trimer_probability(0) == 3


def test_non_neighbor_local_catalogs_are_independent():
    model = OctahedralCoverings(3, 3, 3)
    for cell in model.cells:
        distance = sum(min(q, n-q) for q, n in zip(cell, model.shape))
        if distance <= 1:
            continue
        for first, second in product(LOCAL_STATES, repeat=2):
            assert model.joint_state_probability((0, 0, 0), first, cell, second) == Fraction(1, 256)


def test_nearest_catalogs_only_share_the_one_ownership_bit():
    model = OctahedralCoverings(4, 4, 4)
    for first, second in product(LOCAL_STATES, repeat=2):
        compatible = first.incoming[0] != second.incoming[1]
        expected = Fraction(1, 128) if compatible else Fraction(0)
        assert model.joint_state_probability((0, 0, 0), first, (1, 0, 0), second) == expected


def test_nontrivial_object_covariances_and_finite_support():
    model = OctahedralCoverings(4, 4, 4)
    origin, face = (0, 0, 0), (1, 1, 1)
    assert model.connected_dimer(origin, 0, origin, 0) == Fraction(3, 16)
    assert model.connected_trimer(origin, origin) == Fraction(1, 4)
    assert model.connected_face(origin, face, origin, face) == Fraction(15, 256)
    assert model.connected_dimer(origin, 0, (1, 0, 0), 0) == Fraction(-1, 16)
    assert model.connected_dimer(origin, 0, (0, 1, 0), 0) == 0
    for target in model.cells[1:]:
        distance = sum(min(q, n-q) for q, n in zip(target, model.shape))
        assert model.connected_trimer(origin, target) == 0
        assert model.connected_face(origin, face, target, face) == (Fraction(1, 256) if distance == 1 else 0)
        if distance > 1:
            for first_axis, second_axis in product(range(3), repeat=2):
                assert model.connected_dimer(origin, first_axis, target, second_axis) == 0


def test_short_periodic_cluster_has_extra_correlations():
    model = OctahedralCoverings(2, 2, 2)
    # Every length-two line forces equal trimer bits, hence all eight bits agree.
    assert model.connected_trimer((0, 0, 0), (1, 1, 1)) == Fraction(1, 4)
    assert OctahedralCoverings(3, 3, 3).connected_trimer((0, 0, 0), (1, 1, 1)) == 0


def test_nonlocal_line_parity_survives_zero_two_point_tail():
    model = OctahedralCoverings(3, 4, 5)
    for bits in model.sample(12, seed=251):
        for axis in range(3):
            others = [a for a in range(3) if a != axis]
            for offsets in product(*(range(model.shape[a]) for a in others)):
                parity = 0
                for t in range(model.shape[axis]):
                    cell = [0, 0, 0]
                    cell[axis] = t
                    for a, value in zip(others, offsets):
                        cell[a] = value
                    parity ^= model.local_state(bits, cell).is_trimer
                assert parity == model.shape[axis] % 2


def test_inconsistent_probability_events_and_invalid_covering():
    model = OctahedralCoverings(3, 3, 3)
    assert model.probability(()) == 1
    assert model.probability(((0, 0), (0, 1))) == 0
    assert model.probability_parities(((model.equations[0], 1),)) == 0
    assert not model.valid(1)
    with pytest.raises(ValueError):
        model.covering(1)
    with pytest.raises(ValueError):
        OctahedralCoverings(1, 3, 3)
