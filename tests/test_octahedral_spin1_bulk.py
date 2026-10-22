"""Checks for three-dimensional geometry and physical operator support."""
from fractions import Fraction

import numpy as np
import pytest

from rvb.octahedral_spin1 import (
    cube_contraction,
    cube_geometry,
    levi_civita,
    spin1_generators,
)
from rvb.octahedral_spin1_bulk import (
    build_network,
    edge_label,
    local_base,
    replace_observables,
    site_tag,
)


def test_periodic_cube_has_distinct_parallel_bonds_and_no_missing_neighbors():
    network = build_network((2, 2, 2))
    assert network.num_tensors == 8
    assert len(network.ind_map) == 24
    assert not network.outer_inds()
    for v in np.ndindex(2, 2, 2):
        assert network[site_tag(v)].ndim == 6
        for axis in range(3):
            edge = edge_label(v, axis)
            assert len(network.ind_map[edge]) == 2
            assert network.ind_size(edge) == 8


def test_open_cube_norm_matches_its_two_complete_quantum_coverings():
    # Each corner has three exposed spins, which must all be covered locally.
    # The local state is either its exposed-face trimer or all three dimers.
    # Shared internal sites then force the two checkerboard assignments.
    # Each covering is normalized. Their overlap is the epsilon-cube value
    # 24 divided by the 24 dimer and 8 trimer normalization factors.
    actual = build_network((2, 2, 2), cyclic=False).contract(all, optimize="greedy")
    expected = 2 + Fraction(48, 3**12 * 6**4)
    assert actual == pytest.approx(float(expected), rel=1e-13)


def test_open_cube_spin_insertions_match_independent_covering_overlap():
    edges = (((0, 0, 0), 1), ((1, 0, 0), 1))
    sz = spin1_generators()[2]
    observed = build_network((2, 2, 2), edges, {e: sz for e in edges}, cyclic=False)
    actual = observed.contract(all, optimize="greedy")

    # The two internal spins belong to different dimers in each covering,
    # so only the two cross overlaps contribute. At an even cube vertex
    # the ket owns the internal spin, hence the edge matrix is transposed.
    cells, geometry = cube_geometry()
    matrices = []
    for axis, u, _ in geometry:
        matrix = np.eye(3, dtype=complex)
        if (cells[u], axis) in edges:
            matrix = sz.T if sum(cells[u]) % 2 == 0 else sz
        matrices.append(matrix)
    expected = 2 * cube_contraction([levi_civita()] * 8, matrices) / (3**12 * 6**4)
    assert abs(expected) > 1e-9
    assert actual == pytest.approx(expected, rel=1e-12, abs=1e-20)


def test_full_physical_support_preserves_a_dimer_spin_correlation():
    # Isolate one x-axis dimer within the local octahedron. Its two-spin
    # correlation is -2/3. Compressing its observed legs to the norm support
    # would remove the traceless occupied/occupied Sz components entirely.
    tensor = local_base((0, 1))
    for leg in reversed(range(2, 6)):
        vacuum = np.zeros(8)
        vacuum[0 if leg % 2 else 1] = 1 if leg % 2 else 1 / np.sqrt(3)
        tensor = np.tensordot(tensor, vacuum, axes=(leg, 0))
    sz = spin1_generators()[2]
    trace = np.zeros((4, 4))
    trace[1:, 1:] = np.eye(3)
    insertion = np.zeros((4, 4), dtype=complex)
    insertion[1:, 1:] = sz.T
    norm = np.einsum("ij,i,j", tensor, trace.ravel(), trace.ravel())
    numerator = np.einsum("ij,i,j", tensor, insertion.ravel(), insertion.ravel())
    assert norm == pytest.approx(1)
    assert numerator == pytest.approx(-2 / 3)


def test_identity_full_edges_and_replacement_keep_the_same_norm():
    edges = (((0, 0, 0), 1), ((1, 0, 0), 1))
    shape = (2, 2, 2)
    plain = build_network(shape, cyclic=False).contract(all, optimize="greedy")
    full = build_network(shape, edges, cyclic=False)
    assert full.ind_size(edge_label(*edges[0])) == 16
    assert full.contract(all, optimize="greedy") == pytest.approx(plain, rel=1e-13)
    sz = spin1_generators()[2]
    operators = {e: sz for e in edges}
    replaced = replace_observables(full, shape, edges, operators, cyclic=False)
    direct = build_network(shape, edges, operators, cyclic=False)
    assert replaced.contract(all, optimize="greedy") == pytest.approx(
        direct.contract(all, optimize="greedy"), rel=1e-13, abs=1e-20
    )
    assert full.contract(all, optimize="greedy") == pytest.approx(plain, rel=1e-13)


@pytest.mark.parametrize("shape", [(2, 2), (1, 2, 2), (2.5, 2, 2)])
def test_invalid_geometry_is_rejected(shape):
    with pytest.raises(ValueError):
        build_network(shape)


def test_an_operator_cannot_silently_use_the_compressed_norm_support():
    with pytest.raises(ValueError, match="full support"):
        build_network((2, 2, 2), operators={((0, 0, 0), 1): spin1_generators()[2]})
