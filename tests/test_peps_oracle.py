"""Independent covering-sum checks of blocked tensors and physical insertions."""

import numpy as np
import pytest
import quimb.tensor as qtn

from rvb.exact import (
    bond_expectation,
    connected_bond_correlation,
    rvb_state,
)
from rvb.lattices import kagome, square
from rvb.peps import (
    IDENTITY,
    SM,
    SP,
    SZ,
    RVBTransfer,
    column_transfer,
    double_tensor,
    kagome_site_tensor,
    square_site_tensor,
)


def test_blocked_kagome_torus_matches_explicit_singlet_coverings():
    lattice = kagome(2, 2, periodic_x=True, periodic_y=True)
    oracle, covering_norm = rvb_state(lattice)
    local = kagome_site_tensor()
    tensors = []
    for x in range(2):
        for y in range(2):
            indices = (
                f"p{x}_{y}", f"h{(x - 1) % 2}_{y}", f"h{x}_{y}",
                f"v{x}_{y}", f"v{x}_{(y - 1) % 2}",
            )
            tensors.append(qtn.Tensor(local, inds=indices))
    order = tuple(f"p{x}_{y}" for x in range(2) for y in range(2))
    state = qtn.TensorNetwork(tensors).contract(
        all, output_inds=order, optimize="greedy"
    ).data.reshape(-1)
    norm = float(state @ state)
    assert norm == pytest.approx(covering_norm * 2 ** (lattice.nsites // 2))
    assert abs(state @ oracle) ** 2 / norm == pytest.approx(1., abs=1e-12)


@pytest.mark.parametrize("name,width", [("square", 4), ("kagome", 2)])
def test_infinite_correlations_include_norm_inactive_boundary_states(name, width):
    model = RVBTransfer(name, width)
    norm_tensor = double_tensor(model.site_tensor)
    transfer = column_transfer([norm_tensor] * width)
    if name == "kagome":
        operator = np.kron(
            np.kron(SZ, SZ) + .5 * (np.kron(SP, SM) + np.kron(SM, SP)),
            IDENTITY,
        )
        inserted = column_transfer(
            [double_tensor(model.site_tensor, operator)] + [norm_tensor] * (width - 1)
        )
    else:
        inserted = None
        for operator, coefficient in ((SZ, -1.), (SP, -.5), (SM, -.5)):
            local = double_tensor(model.site_tensor, operator)
            term = coefficient * column_transfer([local, local] + [norm_tensor] * (width - 2))
            inserted = term if inserted is None else inserted + term
    right, left = np.zeros(model.full_dimension), np.zeros(model.full_dimension)
    right[model.active], left[model.active] = model.right, model.left
    vector = inserted @ right / model.eigenvalue
    full_correlations = []
    for _ in range(2):
        joint = left @ (inserted @ vector) / model.eigenvalue
        full_correlations.append(joint - model.mean_bond ** 2)
        vector = transfer @ vector / model.eigenvalue
    # An operator may create norm-inactive boundary states, which contribute
    # at r=1 even though one norm-transfer step annihilates them for r>=2.
    assert model.correlations(2).connected == pytest.approx(full_correlations, abs=1e-12)


def test_kagome_transfer_insertions_match_exact_torus_observables():
    lattice = kagome(2, 2, periodic_x=True, periodic_y=True)
    oracle, covering_norm = rvb_state(lattice)
    local = kagome_site_tensor()
    norm_tensor = double_tensor(local)
    operator = np.kron(
        np.kron(SZ, SZ) + .5 * (np.kron(SP, SM) + np.kron(SM, SP)),
        IDENTITY,
    )
    transfer = column_transfer([norm_tensor, norm_tensor])
    inserted = column_transfer([double_tensor(local, operator), norm_tensor])
    partition = float((transfer @ transfer).diagonal().sum())
    mean = float((inserted @ transfer).diagonal().sum() / partition)
    joint = float((inserted @ inserted).diagonal().sum() / partition)
    assert partition == pytest.approx(covering_norm * 2 ** (lattice.nsites // 2))
    assert mean == pytest.approx(bond_expectation(oracle, (0, 1)), abs=1e-12)
    # Consecutive axial columns require zero norm-transfer columns between
    # insertions, so this also checks the separation-one convention.
    assert joint - mean ** 2 == pytest.approx(
        connected_bond_correlation(oracle, (0, 1), (6, 7)), abs=1e-12
    )


def test_square_marshall_rotation_matches_explicit_singlet_coverings():
    lattice = square(2, 4, periodic_y=True)
    oracle, covering_norm = rvb_state(lattice)
    local = square_site_tensor()
    # U[s,v] = epsilon[v,s] maps a target-end delta bond to an A->B singlet.
    rotation = np.array([[0., -1.], [1., 0.]])
    tensors = []
    for x in range(2):
        for y in range(4):
            data = np.einsum("sp,plrud->slrud", rotation, local) if (x + y) % 2 else local
            # The two open axial boundaries have their virtual legs in vacuum.
            data = data[:, 2, :, :, :] if x == 0 else data[:, :, 2, :, :]
            indices = (f"p{x}_{y}", f"h{y}", f"v{x}_{y}", f"v{x}_{(y - 1) % 4}")
            tensors.append(qtn.Tensor(data, inds=indices))
    order = tuple(f"p{x}_{y}" for x in range(2) for y in range(4))
    state = qtn.TensorNetwork(tensors).contract(
        all, output_inds=order, optimize="greedy"
    ).data.reshape(-1)
    norm = float(state @ state)
    assert norm == pytest.approx(covering_norm * 2 ** (lattice.nsites // 2))
    assert abs(state @ oracle) ** 2 / norm == pytest.approx(1., abs=1e-12)
