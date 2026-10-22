"""Independent finite-graph and analytic checks of orthogonal benchmarks."""

import numpy as np
import pytest

from rvb.orthogonal import (
    KagomeOrthogonalDimers,
    SquareDimerCylinder,
    kagome_orthogonal_benchmark,
    square_plane_benchmark,
)


def perfect_matchings(nsites, edges):
    """Tiny independent exponential oracle, used only on 12/16-site graphs."""
    neighbours = [set() for _ in range(nsites)]
    for i, j in edges:
        neighbours[i].add(j)
        neighbours[j].add(i)

    def visit(remaining, matching):
        if not remaining:
            yield frozenset(matching)
            return
        i = min(remaining)
        for j in sorted(neighbours[i] & remaining):
            yield from visit(remaining - {i, j}, matching + [tuple(sorted((i, j)))])

    return list(visit(set(range(nsites)), []))


def test_square_transfer_against_enumerated_four_by_four_torus():
    width = length = 4
    site = lambda x, y: (x % length) * width + y % width
    edge = lambda a, b: tuple(sorted((a, b)))
    edges = set()
    for x in range(length):
        for y in range(width):
            edges.add(edge(site(x, y), site(x + 1, y)))
            edges.add(edge(site(x, y), site(x, y + 1)))
    all_coverings = perfect_matchings(width * length, edges)
    assert len(all_coverings) == 272
    coverings = [covering for covering in all_coverings
                 if sum((-1) ** y * (edge(site(0, y), site(length - 1, y)) in covering)
                        for y in range(width)) == 0]
    model = SquareDimerCylinder(width)
    transfer = model.matrix.toarray()
    partition = np.trace(np.linalg.matrix_power(transfer, length))
    assert partition == len(coverings)
    horizontal = np.diag(model.states & 1)
    vertical = model.vertical_insertion.toarray()

    for orientation, operator in (("horizontal", horizontal), ("vertical", vertical)):
        def observed_bond(x):
            return (edge(site(x, 0), site(x + 1, 0)) if orientation == "horizontal"
                    else edge(site(x, 0), site(x, 1)))

        for distance in range(1, length):
            exact_joint = np.mean([observed_bond(0) in c and observed_bond(distance) in c
                                   for c in coverings])
            if orientation == "horizontal":
                numerator = np.trace(operator @ np.linalg.matrix_power(transfer, distance)
                                     @ operator @ np.linalg.matrix_power(transfer, length - distance))
            else:
                numerator = np.trace(operator @ np.linalg.matrix_power(transfer, distance - 1)
                                     @ operator @ np.linalg.matrix_power(transfer, length - distance - 1))
            assert numerator / partition == pytest.approx(exact_joint, abs=1e-14)


def test_square_fixed_point_and_exclusion():
    model = SquareDimerCylinder(6)
    matrix = model.matrix
    assert (matrix - matrix.T).nnz == 0
    lam, fixed = model.fixed_point
    assert np.linalg.norm(matrix @ fixed - lam * fixed) < 1e-11
    horizontal, _ = model._observable("horizontal")
    vertical, _ = model._observable("vertical")
    assert 2 * horizontal + 2 * vertical == pytest.approx(1)
    assert model.correlations(1)[1] == pytest.approx(-horizontal**2)
    assert model.correlations(0)[0] == pytest.approx(horizontal * (1 - horizontal))


def test_width_four_known_spectrum():
    model = SquareDimerCylinder(4)
    assert model.fixed_point[0] == pytest.approx(2 + np.sqrt(3))
    assert model.transfer_length()["correlation_length"] == pytest.approx(1 / np.log(2 + np.sqrt(3)))


def test_square_plane_integral_and_power_law():
    result = square_plane_benchmark(128)
    correlation = result["connected_correlation"]
    assert correlation[0] == pytest.approx(3 / 16)
    assert correlation[1] == pytest.approx(-1 / 16)
    # g(3)=2/pi-3/4 follows by substituting t=sin(k) in the integral.
    assert correlation[2] == pytest.approx(-(2 / np.pi - 3 / 4) / 4, abs=1e-13)
    assert correlation[3] == pytest.approx(-(2 / np.pi - 3 / 4) ** 2, abs=1e-13)
    assert result["fitted_even_distance_exponent"] == pytest.approx(2, abs=0.01)
    assert correlation[128] * np.pi**2 * 128**2 == pytest.approx(1, abs=0.001)


def test_square_widening_approaches_plane_and_gap_closes():
    plane = np.array(square_plane_benchmark(3)["connected_correlation"])
    narrow, wide = SquareDimerCylinder(4), SquareDimerCylinder(8)
    assert np.linalg.norm(wide.correlations(3) - plane) < np.linalg.norm(narrow.correlations(3) - plane)
    assert wide.transfer_length()["correlation_length"] > narrow.transfer_length()["correlation_length"]


def test_square_plane_transverse_bonds_and_parity_cancellation():
    result = square_plane_benchmark(128, orientation="vertical")
    correlations = np.array(result["connected_correlation"])
    assert correlations[1] == pytest.approx(1 / 16)
    assert correlations[2] == pytest.approx(-(2 / np.pi - 1 / 2) ** 2 / 4)
    assert result["fitted_odd_distance_exponent"] == pytest.approx(2, abs=0.01)
    assert result["fitted_even_distance_exponent"] == pytest.approx(4, abs=0.01)
    assert -correlations[128] * np.pi**2 * 128**4 == pytest.approx(1, abs=0.001)
    narrow, wide = SquareDimerCylinder(4), SquareDimerCylinder(8)
    assert (np.linalg.norm(wide.correlations(3, "vertical") - correlations[:4])
            < np.linalg.norm(narrow.correlations(3, "vertical") - correlations[:4]))
    assert wide.transfer_length("vertical")["correlation_length"] == pytest.approx(
        wide.transfer_length("horizontal")["correlation_length"])


def test_kagome_arrow_bijection_and_all_pair_probabilities():
    model = KagomeOrthogonalDimers(2, 2)
    coverings = perfect_matchings(model.nsites, model.edges)
    assert len(coverings) == 2**model.free_arrows == 32
    for edge in model.edges:
        assert model.occupation(edge) == np.mean([edge in c for c in coverings]) == 0.25
        for other in model.edges:
            exact = np.mean([edge in c and other in c for c in coverings])
            assert model.joint_occupation(edge, other) == exact


def test_kagome_large_torus_exact_locality():
    model = KagomeOrthogonalDimers(8, 8)
    assert model.free_arrows == 8 * 8 + 1
    reference = (model.site(0, 0, 0), model.site(0, 0, 1))
    assert model.connected(reference, reference) == 3 / 16
    # Site-sharing bonds cannot both be occupied.
    shared = (model.site(0, 0, 0), model.site(0, 0, 2))
    assert model.connected(reference, shared) == -1 / 16
    # Two bonds have a connected correlation only when their host triangles
    # coincide or are adjacent in the underlying honeycomb graph.
    host_triangles = [t for t in model.triangles if set(reference).issubset(t)]
    local_sites = set(host_triangles[0])
    adjacent = [t for t in model.triangles if local_sites.intersection(t)]
    local_edges = {tuple(sorted((t[i], t[j]))) for t in adjacent
                   for i, j in ((0, 1), (1, 2), (2, 0))}
    for other in model.edges:
        if other not in local_edges:
            assert model.connected(reference, other) == 0
    result = kagome_orthogonal_benchmark(8, 8)
    assert result["connected_correlation"][1:] == [0, 0, 0, 0]


@pytest.mark.parametrize("width", [1, 2, 3, 5, 18])
def test_invalid_square_widths(width):
    with pytest.raises(ValueError):
        SquareDimerCylinder(width)


def test_invalid_kagome_parity():
    with pytest.raises(ValueError):
        KagomeOrthogonalDimers(3, 3)
