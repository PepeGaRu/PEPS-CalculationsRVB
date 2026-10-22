"""Independent finite-cluster oracles; exponential cost, intended for validation."""

from itertools import product
import numpy as np


def perfect_matchings(lattice, max_coverings=100_000):
    """Enumerate coverings as edge-index tuples, failing instead of truncating."""
    if lattice.nsites % 2:
        return []
    incident = [[] for _ in range(lattice.nsites)]
    for e, (a, b) in enumerate(lattice.edges):
        incident[a].append((b, e))
        incident[b].append((a, e))
    result = []

    def visit(remaining, chosen):
        if not remaining:
            if len(result) >= max_coverings:
                raise ValueError("Covering limit exceeded; use cylinder transfer calculations")
            result.append(tuple(chosen))
            return
        a = min(remaining, key=lambda v: sum(b in remaining for b, _ in incident[v]))
        for b, e in incident[a]:
            if b in remaining:
                visit(remaining - {a, b}, chosen + [e])

    visit(set(range(lattice.nsites)), [])
    return result


def rvb_state(lattice, coverings=None, max_sites=20):
    """Sum normalized oriented singlet products, then normalize the total state.

    The returned second value is the norm squared before final normalization.
    No covering-overlap approximation is made.
    """
    n = lattice.nsites
    if n > max_sites:
        raise ValueError(f"Dense oracle limited to {max_sites} spins; got {n}")
    coverings = perfect_matchings(lattice) if coverings is None else coverings
    state = np.zeros(2**n)
    factor = 2.0 ** (-n / 4)
    for covering in coverings:
        for bits in product((0, 1), repeat=len(covering)):
            index, sign = 0, 1
            for e, bit in zip(covering, bits):
                a, b = lattice.edges[e]
                index |= bit << (n-1-a)
                index |= (1-bit) << (n-1-b)
                sign *= 1 if bit == 0 else -1
            state[index] += factor * sign
    norm = float(state @ state)
    if norm <= 0:
        raise ValueError("No nonzero singlet RVB state on this graph")
    return state / np.sqrt(norm), norm


def spin_bond_action(state, edge):
    """Apply S_i·S_j = (1/2) SWAP_ij - (1/4) I for spin 1/2."""
    n = state.size.bit_length() - 1
    a, b = edge
    tensor = state.reshape((2,) * n)
    return .5 * tensor.swapaxes(a, b).reshape(-1) - .25 * state


def bond_expectation(state, edge):
    return float(np.vdot(state, spin_bond_action(state, edge)).real)


def connected_bond_correlation(state, edge_a, edge_b):
    if set(edge_a) & set(edge_b):
        raise ValueError("Choose disjoint bonds for commuting dimer observables")
    mean_a, mean_b = bond_expectation(state, edge_a), bond_expectation(state, edge_b)
    joint = np.vdot(spin_bond_action(state, edge_a), spin_bond_action(state, edge_b)).real
    return float(joint - mean_a * mean_b)


def occupancy_statistics(lattice, coverings=None):
    coverings = perfect_matchings(lattice) if coverings is None else coverings
    if not coverings:
        raise ValueError("No perfect matchings")
    occupancy = np.zeros((len(coverings), len(lattice.edges)))
    for k, covering in enumerate(coverings):
        occupancy[k, list(covering)] = 1
    mean = occupancy.mean(axis=0)
    connected = occupancy.T @ occupancy / len(coverings) - np.outer(mean, mean)
    return mean, connected


def peps_state(lattice, max_sites=16):
    """Contract a graph D=3 RVB PEPS with quimb; independent of covering oracle."""
    import quimb.tensor as qtn
    if lattice.nsites > max_sites:
        raise ValueError(f"Dense PEPS check limited to {max_sites} spins")
    epsilon = np.array([[0., 1.], [-1., 0.]])
    tensors = []
    for site in range(lattice.nsites):
        incident = [e for e, edge in enumerate(lattice.edges) if site in edge]
        data = np.zeros((2,) + (3,) * len(incident))
        for leg, e in enumerate(incident):
            for virtual in (0, 1):
                indices = [2] * len(incident)
                indices[leg] = virtual
                for spin in (0, 1):
                    value = float(spin == virtual) if site == lattice.edges[e][0] else epsilon[virtual, spin]
                    data[(spin, *indices)] = value * 2**(-.25)
        tensors.append(qtn.Tensor(data, inds=(f"s{site}", *(f"e{e}" for e in incident))))
    result = qtn.TensorNetwork(tensors).contract(all, output_inds=tuple(f"s{i}" for i in range(lattice.nsites)), optimize="greedy")
    state = np.asarray(result.data).reshape(-1)
    norm = float(state @ state)
    if norm <= 0:
        raise ValueError("PEPS contraction produced zero state")
    return state / np.sqrt(norm), norm


def cluster_benchmark(lattice):
    coverings = perfect_matchings(lattice)
    state, norm = rvb_state(lattice, coverings)
    peps, peps_norm = peps_state(lattice)
    fidelity = float(abs(np.vdot(peps, state))**2)
    if not np.isclose(fidelity, 1, atol=1e-11) or not np.isclose(norm, peps_norm, atol=1e-11):
        raise AssertionError("PEPS disagrees with independent covering sum")
    means, connected = occupancy_statistics(lattice, coverings)
    reference = 0
    comparisons = []
    for e, edge in enumerate(lattice.edges):
        if not set(edge) & set(lattice.edges[reference]):
            comparisons.append({"edge_index": e, "edge": list(edge),
                                "orthogonal_connected_occupation": float(connected[reference, e]),
                                "rvb_connected_bond_energy": connected_bond_correlation(state, lattice.edges[reference], edge)})
    total_s2 = .75 * lattice.nsites + 2 * sum(bond_expectation(state, (a,b)) for a in range(lattice.nsites) for b in range(a+1, lattice.nsites))
    return {"lattice": lattice.name, "sites": lattice.nsites, "edges": [list(e) for e in lattice.edges],
            "coverings": len(coverings), "singlet_covering_sum_norm_squared": norm,
            "peps_fidelity_to_covering_sum": fidelity, "total_spin_squared": total_s2,
            "reference_edge": list(lattice.edges[reference]), "orthogonal_mean_occupation": float(means[reference]),
            "rvb_mean_bond_energy": bond_expectation(state, lattice.edges[reference]),
            "correlations": comparisons,
            "scope": "Exact finite cluster; does not establish thermodynamic correlation length or exponent."}
