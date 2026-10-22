"""Check observable transport through shared 3D HOTRG projectors."""
import numpy as np
import quimb.tensor as qtn

from scripts.probe_octahedral_spin1_bulk_spin import shared_coarse_grain


def test_untruncated_shared_projectors_match_independent_exact_contraction():
    # Open 3x2x2 rank-two network: chi=4 keeps every paired external bond.
    norm = qtn.TN3D_rand(3, 2, 2, 2, seed=23)
    impurity = norm.copy()
    tid = next(iter(impurity.tensor_map))
    tensor = impurity.tensor_map[tid]
    tensor.modify(data=tensor.data + np.ones(tensor.shape))
    exact_norm = norm.contract(optimize="greedy")
    exact_impurity = impurity.contract(optimize="greedy")
    coarse_norm, coarse_impurities = shared_coarse_grain(
        norm, [impurity], "x", 4, cutoff=0.)
    assert np.isclose(coarse_norm.contract(optimize="greedy"),
                      exact_norm, rtol=1e-12)
    assert np.isclose(coarse_impurities[0].contract(optimize="greedy"),
                      exact_impurity, rtol=1e-12)


def test_truncated_shared_projectors_preserve_linear_observable_scaling():
    norm = qtn.TN3D_rand(3, 2, 2, 2, seed=27)
    impurity = norm.copy()
    tid = next(iter(impurity.tensor_map))
    tensor = impurity.tensor_map[tid]
    tensor.modify(data=3.25 * tensor.data)
    coarse_norm, coarse_impurities = shared_coarse_grain(
        norm, [impurity], "x", 2)
    denominator = coarse_norm.contract(optimize="greedy")
    numerator = coarse_impurities[0].contract(optimize="greedy")
    assert np.isclose(numerator / denominator, 3.25, rtol=1e-12)
