"""Finite three-dimensional octahedral networks and exact physical impurities.

Every linear dimension is finite here. A bulk correlation length requires
convergence as all dimensions and contraction accuracy increase. Periodic
odd linear sizes force odd trimer parity on their straight lines; even sizes
permit dimer-only coverings. Prefer an even-size sequence for bulk convergence.

No observable is projected onto the rank-eight norm support on its own edge.
A full sixteen-dimensional observed bond joins one raw endpoint to an endpoint
with the physical operator metric absorbed. All other bonds use the complete
rank-eight norm support, which is an exact reduction.
"""
from functools import lru_cache
from itertools import product

import numpy as np
import quimb.tensor as qtn

from .octahedral_spin1 import octahedron_tensor, spin1_generators
from .octahedral_spin1_transfer import (
    EDGE_PERMUTATION,
    EDGE_WEIGHTS,
    full_edge_metric,
)


@lru_cache(None)
def local_base(full_legs=()):
    """Double tensor, with raw sixteen-dimensional indices on ``full_legs``.

    Leg order is +x,-x,+y,-y,+z,-z. Other indices use the exact rank-eight
    norm basis, with the norm metric absorbed on positive-directed legs.
    The cached result should be treated as read-only.
    """
    full_legs = tuple(full_legs)
    single = octahedron_tensor()
    index = np.array(single.nonzero()).T
    value = single[tuple(index.T)]
    ket, bra = index[:, None, :], index[None, :, :]
    compressed = tuple(i for i in range(6) if i not in full_legs)
    valid = np.all(
        ((ket == 0) | (bra == 0) | (ket == bra))[:, :, compressed], axis=-1
    )
    first, second = np.where(valid)
    ket, bra = index[first], index[second]
    labels = np.zeros(ket.shape, dtype=int)
    amplitudes = value[first] * value[second]
    for leg in range(6):
        if leg in full_legs:
            labels[:, leg] = 4 * ket[:, leg] + bra[:, leg]
            continue
        both = (ket[:, leg] > 0) & (bra[:, leg] > 0)
        ket_only = (ket[:, leg] > 0) & (bra[:, leg] == 0)
        bra_only = (ket[:, leg] == 0) & (bra[:, leg] > 0)
        labels[both, leg] = 1
        labels[ket_only, leg] = ket[ket_only, leg] + 1
        labels[bra_only, leg] = bra[bra_only, leg] + 4
        amplitudes[both] /= np.sqrt(3)
        if leg % 2 == 0:
            amplitudes *= EDGE_WEIGHTS[labels[:, leg]]
            labels[:, leg] = EDGE_PERMUTATION[labels[:, leg]]
    result = np.zeros(tuple(16 if i in full_legs else 8 for i in range(6)))
    np.add.at(result, tuple(labels.T), amplitudes)
    return result


def apply_factor(data, leg, factor):
    """Replace raw index a by j using factor[a,j], preserving axis order."""
    result = np.tensordot(data, factor, axes=(leg, 0))
    return np.moveaxis(result, -1, leg)


def shifted(vertex, axis, amount, shape, cyclic):
    result = list(vertex)
    result[axis] += amount
    if cyclic:
        result[axis] %= shape[axis]
    return tuple(result)


def edge_label(vertex, axis):
    """Label an oriented edge, preserving both parallel bonds at period two."""
    return "e" + "_".join(map(str, vertex)) + "_" + str(axis)


def site_tag(vertex):
    return "ORIG_" + "_".join(map(str, vertex))


def build_network(shape, observed_edges=(), operators=None, cyclic=True):
    """Create a closed tensor network for a finite three-dimensional box.

    ``observed_edges`` contains ``(origin_xyz, axis)`` positive-directed
    physical bonds. ``operators`` maps such edges to 3x3 physical operators;
    an omitted operator is the identity. These bonds retain dimension sixteen.

    Periodic geometry has 3*Lx*Ly*Lz physical spins. Open geometry keeps all six
    physical ports at every boundary octahedron, forcing an exposed spin to
    belong to the only existing octahedron. This is a particular boundary
    state, not the same finite model as a periodic box. An open cubic box of
    side L has 3*L**2*(L+1) physical spins.
    """
    shape = tuple(shape)
    if len(shape) != 3 or any(
        not isinstance(s, (int, np.integer)) or s < 2 for s in shape
    ):
        raise ValueError("Three integer dimensions >=2 are required")
    observed_edges = tuple((tuple(v), int(a)) for v, a in observed_edges)
    operators = {} if operators is None else operators
    if len(set(observed_edges)) != len(observed_edges):
        raise ValueError("Observed edges must be distinct")
    if set(operators) - set(observed_edges):
        raise ValueError("Every operator edge must retain full support")
    for vertex, axis in observed_edges:
        if (
            len(vertex) != 3
            or axis not in range(3)
            or any(not 0 <= vertex[i] < shape[i] for i in range(3))
        ):
            raise ValueError("Observed edge outside the lattice")
    for operator in operators.values():
        if np.shape(operator) != (3, 3) or not np.all(np.isfinite(operator)):
            raise ValueError("A physical spin-one operator must be a finite 3 by 3 array")

    full_factors = {v: {} for v in product(*(range(s) for s in shape))}
    for vertex, axis in observed_edges:
        neighbor = shifted(vertex, axis, 1, shape, cyclic)
        if neighbor not in full_factors:
            raise ValueError("Observed edge must be internal")
        full_factors[vertex][2 * axis] = np.eye(16)
        full_factors[neighbor][2 * axis + 1] = full_edge_metric(
            operators.get((vertex, axis))
        ).T

    network = qtn.TensorNetwork3D.new(
        Lx=shape[0], Ly=shape[1], Lz=shape[2],
        site_tag_id="I{},{},{}", x_tag_id="X{}", y_tag_id="Y{}", z_tag_id="Z{}",
    )
    for vertex, factors in full_factors.items():
        data = local_base(tuple(sorted(factors))).copy()
        for leg, factor in factors.items():
            data = apply_factor(data, leg, factor)
        indices = []
        boundaries = []
        for axis in range(3):
            plus = shifted(vertex, axis, 1, shape, False)
            minus = shifted(vertex, axis, -1, shape, False)
            for positive, neighbor in ((True, plus), (False, minus)):
                leg = 2 * axis + int(not positive)
                boundary = not cyclic and not 0 <= neighbor[axis] < shape[axis]
                if boundary:
                    # Positive-leg norm metric is already absorbed; these two
                    # vectors both force the local endpoint to own the spin.
                    trace = np.zeros(8)
                    trace[0 if positive else 1] = 1 if positive else np.sqrt(3)
                    boundaries.append((leg, trace))
                    indices.append(None)
                else:
                    origin = vertex if positive else shifted(vertex, axis, -1, shape, cyclic)
                    indices.append(edge_label(origin, axis))
        # Remove boundary dimensions in descending order to preserve axes.
        for leg, trace in reversed(boundaries):
            data = np.tensordot(data, trace, axes=(leg, 0))
        indices = [i for i in indices if i is not None]
        tags = [
            network.site_tag(*vertex), network.x_tag(vertex[0]),
            network.y_tag(vertex[1]), network.z_tag(vertex[2]), site_tag(vertex),
        ]
        network.add_tensor(qtn.Tensor(data, inds=indices, tags=tags))
    return network


def replace_observables(norm_network, shape, observed_edges, operators, cyclic=True):
    """Replace original tensors while keeping a shared contraction environment.

    Use only when original observable-endpoint tensors survive untouched,
    including their bond indices. This helper checks that condition. It does
    not ensure an approximate environment has converged or retained every
    relevant spin channel.
    """
    impurity = build_network(shape, observed_edges, operators, cyclic=cyclic)
    result = norm_network.copy()
    changed = {shifted(v, axis, 1, shape, cyclic) for v, axis in operators}
    for vertex in changed:
        destination = result[site_tag(vertex)]
        source = impurity[site_tag(vertex)]
        if not isinstance(destination, qtn.Tensor) or destination.inds != source.inds:
            raise ValueError("Observable endpoint was compressed; do not insert blindly")
        destination.modify(data=source.data)
    return result
