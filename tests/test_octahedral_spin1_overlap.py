import math

import opt_einsum as oe
import pytest

from rvb.octahedral import OctahedralCoverings, binary_rank
from rvb.octahedral_spin1 import spin1_dimer, spin1_generators, spin1_trimer
from rvb.octahedral_spin1_overlap import contract_pair, line_flip_mask


def test_diagonal_normalization_and_object_spin_correlations():
    model = OctahedralCoverings(3, 3, 3)
    for covering in model.sample(5, seed=812):
        assert contract_pair(model, covering, covering).value == pytest.approx(1)
        objects = model.covering(covering)
        for kind, expected in (("dimer", -2/3), ("trimer", -1/3)):
            obj = next(o for o in objects if o.kind == kind)
            result = contract_pair(model, covering, covering, sz_sites=obj.sites[:2])
            assert result.value == pytest.approx(expected)
        sites = objects[0].sites[0], objects[1].sites[0]
        assert contract_pair(model, covering, covering, sz_sites=sites).numerator == 0


def test_local_cube_reference_overlap_and_hermiticity():
    model = OctahedralCoverings(3, 3, 3)
    a = model.all_minus_trimer_reference
    b = a ^ model.cube_flip_mask(0)
    assert contract_pair(model, a, b).value == pytest.approx(1/39366)
    for sites in (None, (0, 3), (0, 1)):
        ab = contract_pair(model, a, b, sz_sites=sites)
        ba = contract_pair(model, b, a, sz_sites=sites)
        assert ab.value == pytest.approx(ba.value)
        assert math.isfinite(ab.value)


def test_line_and_cube_moves_span_all_coverings():
    for length in (3, 4):
        model = OctahedralCoverings(length, length, length)
        moves = [line_flip_mask(model, a, c)
                 for a in range(3) for c in range(model.ncells)]
        moves += [model.cube_flip_mask(c) for c in range(model.ncells)]
        assert all(model.valid(move) for move in moves)
        assert binary_rank(moves) == model.dimension


def test_winding_line_overlap_formula():
    for length in (3, 4):
        model = OctahedralCoverings(length, length, length)
        for covering in model.sample(4, seed=913):
            for axis in range(3):
                line_cells = {model.cell_index(tuple(i if a == axis else 0 for a in range(3)))
                              for i in range(length)}
                trimers = sum(obj.kind == "trimer" and obj.cell in line_cells
                              for obj in model.covering(covering))
                flipped = covering ^ line_flip_mask(model, axis)
                exact = contract_pair(model, covering, flipped).value
                assert exact == pytest.approx(3**(1-(length+trimers)/2))


def test_invalid_spin_sites_are_rejected():
    model = OctahedralCoverings(3, 3, 3)
    with pytest.raises(ValueError):
        contract_pair(model, 0, 0, sz_sites=(0, 0))


def test_off_diagonal_spin_insertion_against_unreduced_network():
    """Independent raw tensors retain every dimer and physical spin index."""
    model = OctahedralCoverings(3, 3, 3)
    ket = model.all_minus_trimer_reference
    for bra in (ket ^ model.cube_flip_mask(0), ket ^ line_flip_mask(model, 0)):
        for sites in ((0, 1), (0, 3), (0, model.edge_index((1, 0, 0), 0))):
            args = []
            def index(layer, site):
                return site + (model.nsites if layer == 1 and site in sites else 0)
            for layer, covering in enumerate((ket, bra)):
                for obj in model.covering(covering):
                    tensor = spin1_dimer() if obj.kind == "dimer" else spin1_trimer()
                    args.extend((tensor, [index(layer, site) for site in obj.sites]))
            for site in sites:
                args.extend((spin1_generators()[2], [index(1, site), index(0, site)]))
            direct = oe.contract(*args, [], optimize="greedy")
            reduced = contract_pair(model, ket, bra, sz_sites=sites).value
            assert complex(direct).imag == pytest.approx(0)
            assert float(direct.real) == pytest.approx(reduced, abs=1e-13)
