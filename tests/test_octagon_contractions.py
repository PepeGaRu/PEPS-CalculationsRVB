"""Check the factorized two-cell measurements against small dense operators.

These integration tests require the optional pinned peps-torch environment.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

from scripts.run_octagon import factorized_pair_correlations


@pytest.fixture
def upstream():
    torch = pytest.importorskip("torch")
    source = Path("/workspace/peps-environment/upstream/peps-torch")
    if not source.exists():
        pytest.skip("optional pinned peps-torch checkout is not installed")
    sys.path.insert(0, str(source))
    try:
        import config as cfg
        from ctm.generic import corrf
        from ctm.generic.env import ENV, init_env
        from ipeps.ipeps import IPEPS
        cfg.global_args.dtype = "float64"
        cfg.global_args.torch_dtype = torch.float64
        cfg.global_args.device = "cpu"
        torch.set_num_threads(2)
        yield torch, corrf, ENV, init_env, IPEPS
    finally:
        sys.path.remove(str(source))


@pytest.mark.parametrize("direction", [(1, 0), (-1, 0), (0, 1), (0, -1)])
def test_factorized_pair_matches_dense_small_physical_operator(upstream, direction):
    torch, corrf, ENV, init_env, IPEPS = upstream
    torch.manual_seed(173)
    a = torch.randn((2, 2, 2, 2, 2), dtype=torch.float64)
    a /= torch.linalg.norm(a)
    state = IPEPS({(0, 0): a}, vertexToSite=lambda coord: (0, 0), lX=1, lY=1)
    env = ENV(3, state)
    init_env(state, env)
    z = torch.tensor([[.5, 0.], [0., -.5]], dtype=a.dtype)
    plus = torch.tensor([[0., 1.], [0., 0.]], dtype=a.dtype)
    minus = plus.T
    terms = [(1., z, z), (.5, plus, minus), (.5, minus, plus)]
    dense = sum(coef * torch.einsum("ij,kl->ikjl", first, second)
                for coef, first, second in terms)
    actual = factorized_pair_correlations(corrf, (0, 0), direction, state, env,
                                         terms, terms, max_distance=5)
    expected = corrf.corrf_2sOH2sOH_E1((0, 0), direction, state, env,
                                      dense, lambda r: dense, dist=3)
    assert np.allclose(actual.numpy(), expected.numpy(), atol=2e-13, rtol=2e-12)


def test_factorized_pair_in_polarized_product_state(upstream):
    torch, corrf, ENV, init_env, IPEPS = upstream
    a = torch.tensor([1., 0.], dtype=torch.float64).reshape(2, 1, 1, 1, 1)
    state = IPEPS({(0, 0): a}, vertexToSite=lambda coord: (0, 0), lX=1, lY=1)
    env = ENV(1, state)
    init_env(state, env)
    z = torch.tensor([[.5, 0.], [0., -.5]], dtype=a.dtype)
    plus = torch.tensor([[0., 1.], [0., 0.]], dtype=a.dtype)
    terms = [(1., z, z), (.5, plus, plus.T), (.5, plus.T, plus)]
    raw = factorized_pair_correlations(corrf, (0, 0), (1, 0), state, env,
                                      terms, terms, 5)
    centered_terms = terms + [(-.25, None, None)]
    centered = factorized_pair_correlations(corrf, (0, 0), (1, 0), state, env,
                                           centered_terms, centered_terms, 5)
    assert np.allclose(raw.numpy(), 1 / 16, atol=1e-15)
    assert np.allclose(centered.numpy(), 0, atol=1e-15)
