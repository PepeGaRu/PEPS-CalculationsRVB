"""Exact 2 by 3 transverse spin-one singlet-covering cylinder.

No tensor truncation is used. The z-ring has length three and the y-ring
length two. Each slice contains six octahedra and eighteen physical spins.
Ownership constraints split the rank-eight double-layer cut into 256
blocks. The diagonal bra/ket blocks have dimension 976. This module keeps
all colors within the requested sector; the scalar physical observable is
inserted before compression on its own y-directed physical sites.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product, combinations_with_replacement
from math import factorial
from time import perf_counter

import numpy as np
from numba import njit
from scipy import linalg, sparse
from scipy.sparse.linalg import eigs

from .octahedral_spin1 import spin1_dimer, spin1_generators
from .octahedral_spin1_transfer import Spin1OctahedralCylinder

# Each compressed edge component has a definite ket/bra ownership flag.
FLAGS = np.array([(0, 0), (1, 1)] + [(1, 0)]*3 + [(0, 1)]*3)


@njit(cache=True)
def _ring_three(indices, dimensions, last, first_ptr, rows, columns, values, second_ptr):
    n = len(indices)
    result = np.zeros((n, n), dtype=last.dtype)
    for outgoing in range(n):
        for incoming in range(n):
            q0 = indices[outgoing, 0]*dimensions[0]+indices[incoming, 0]
            q1 = indices[outgoing, 1]*dimensions[1]+indices[incoming, 1]
            q2 = indices[outgoing, 2]*dimensions[2]+indices[incoming, 2]
            value = 0.0*last[0, 0, 0]
            for e in range(first_ptr[q0], first_ptr[q0+1]):
                i, j = rows[e], columns[e]
                for f in range(second_ptr[q1, j], second_ptr[q1, j+1]):
                    value += values[e]*values[f]*last[q2, columns[f], i]
            result[outgoing, incoming] = value
    return result


def _contract_ring(tensors, indices):
    """Tr(Q0 Q1 Q2) with sparse first and second factors, for each cut pair."""
    first_ptr, rows, columns, values = [0], [], [], []
    for q in tensors[0]:
        r, c = np.nonzero(q)
        rows.extend(r); columns.extend(c); values.extend(q[r, c])
        first_ptr.append(len(rows))
    second_ptr = np.zeros((len(tensors[1]), 65), dtype=np.int64)
    for i, q in enumerate(tensors[1]):
        r, c = np.nonzero(q)
        offset = len(rows)
        rows.extend(r); columns.extend(c); values.extend(q[r, c])
        second_ptr[i] = np.r_[0, np.cumsum(np.bincount(r, minlength=64))]+offset
    dtype = np.result_type(*(q.dtype for q in tensors))
    return _ring_three(np.asarray(indices),
                       np.array([round(len(q)**.5) for q in tensors]),
                       np.asarray(tensors[2], dtype=dtype), np.array(first_ptr),
                       np.array(rows), np.array(columns), np.array(values, dtype=dtype), second_ptr)


class Spin1OctahedralCylinder23:
    """Infinite x and periodic (Ly,Lz)=(2,3), with exact finite-width transfer."""
    transverse_shape = (2, 3)
    sites_per_slice = 18
    boundary_dimension = 8**6

    def __init__(self):
        # Reuse the independently tested local tensors and physical connectors;
        # the cross-section contraction below is a different ring topology.
        self.base = Spin1OctahedralCylinder()

    @lru_cache(maxsize=None)
    def _pair_data(self, ket_difference, bra_difference, observable=None):
        if observable is None:
            operator = None
        elif observable == "identity":
            operator = np.eye(9)
        elif observable == "dimer_projector":
            d = spin1_dimer().ravel()
            operator = np.outer(d, d)
        elif observable == "spin_z":
            operator = np.kron(spin1_generators()[2], np.eye(3))
        else:
            raise ValueError("Unknown physical observable")
        matrix = self.base._pair_matrix(0, operator).tocoo()
        p0, p1 = divmod(matrix.row, 64)
        z0, z1 = divmod(matrix.col, 64)
        r0, l0 = divmod(p0, 8); r1, l1 = divmod(p1, 8)
        zr0, zl0 = divmod(z0, 8); zr1, zl1 = divmod(z1, 8)
        outgoing, incoming = r0*8+r1, l0*8+l1
        z_out, z_in = zr0*8+zr1, zl0*8+zl1
        pairs = np.array([(a, b) for a in range(8) for b in range(8)
                          if np.array_equal(FLAGS[a] ^ FLAGS[b], (ket_difference, bra_difference))])
        lookup = np.full(64, -1, dtype=int)
        lookup[pairs[:, 0]*8+pairs[:, 1]] = np.arange(len(pairs))
        valid = (lookup[outgoing] >= 0) & (lookup[incoming] >= 0)
        row = lookup[outgoing[valid]]*len(pairs)+lookup[incoming[valid]]
        q = sparse.coo_matrix((matrix.data[valid], (row, z_out[valid]*64+z_in[valid])),
                              shape=(len(pairs)**2, 4096)).toarray().reshape(len(pairs)**2, 64, 64)
        return pairs, q

    def block(self, ket_difference=(0, 0, 0), bra_difference=None, parities=(0, 0), observable=None):
        """Return a dense sector transfer and its compact local boundary basis.

        Differences specify ket/bra ownership XOR between the two y sites,
        separately at each z. Parities specify XOR of the three y=0 flags.
        These 3+1 invariants give sixteen sectors in each individual layer.
        """
        if bra_difference is None:
            bra_difference = ket_difference
        ket_difference, bra_difference = tuple(ket_difference), tuple(bra_difference)
        if (len(ket_difference) != 3 or len(bra_difference) != 3 or len(parities) != 2
                or any(v not in (0, 1) for v in (*ket_difference, *bra_difference, *parities))):
            raise ValueError("Differences require three binary values; parities require two")
        data = [self._pair_data(k, b, observable if z == 0 else None)
                for z, (k, b) in enumerate(zip(ket_difference, bra_difference))]
        pairs = [item[0] for item in data]
        indices = np.array(list(product(*(range(len(p)) for p in pairs))))
        bits = sum(FLAGS[p[indices[:, z], 0]] for z, p in enumerate(pairs))
        indices = indices[np.all(bits % 2 == parities, axis=1)]
        matrix = _contract_ring([item[1] for item in data], indices)
        return matrix, indices, pairs

    def diagonal_spectra(self):
        records = []
        for differences in product((0, 1), repeat=3):
            for parity in (0, 1):
                block, _, _ = self.block(differences, parities=(parity, parity))
                values = eigs(sparse.csr_matrix(block), k=6, which="LM", tol=1e-12,
                              v0=np.ones(len(block)), return_eigenvectors=False)
                values = values[np.argsort(-abs(values))]
                records.append({"difference": differences, "parity": parity,
                                "dimension": len(block), "eigenvalues": values})
        return records

    def finite_length_two_norm(self):
        """Trace of the complete transfer squared, summing all256 sectors.

        The z-ring of length3 has the full permutation symmetry of a triangle.
        We contract80 representatives of its256 sector labels and multiply
        by exact permutation multiplicities. This gives an independent
        axis-rotation check against Tr(T_2x2**3), the same36-spin torus.
        """
        total, records = 0.0, []
        types = ((0, 0), (0, 1), (1, 0), (1, 1))
        for columns in combinations_with_replacement(types, 3):
            multiplicity = factorial(3)
            for kind in types:
                multiplicity //= factorial(columns.count(kind))
            for parities in product((0, 1), repeat=2):
                matrix, _, _ = self.block(tuple(c[0] for c in columns),
                                          tuple(c[1] for c in columns), parities)
                value = float(np.einsum("ij,ji->", matrix, matrix))
                total += multiplicity*value
                records.append({"column_types": columns, "parities": parities,
                                "multiplicity": multiplicity, "trace_squared": value})
        return total, records

    def analyze(self, observable="dimer_projector", distances=range(1, 65), parity=0):
        """Full976-state spectrum and centered physical correlations in d=000."""
        distances = tuple(int(r) for r in distances)
        if not distances or min(distances) < 1:
            raise ValueError("Distances must be positive integers")
        matrix, _, _ = self.block(parities=(parity, parity))
        inserted, _, _ = self.block(parities=(parity, parity), observable=observable)
        values, left, right = linalg.eig(matrix, left=True, right=True)
        order = np.argsort(-abs(values)); values = values[order]
        left, right = left[:, order].conj(), right[:, order]
        lambda0 = float(values[0].real)
        r0, l0 = right[:, 0].real, left[:, 0].real
        l0 = l0/(l0 @ r0)
        operator = inserted/lambda0
        mean = l0 @ operator @ r0
        source = operator @ r0-mean*r0
        sink = l0 @ operator-mean*l0
        vector, connected = source.copy(), {}
        normalized = matrix/lambda0
        for distance in range(1, max(distances)+1):
            if distance in distances:
                connected[distance] = sink @ vector
            vector = normalized @ vector
            vector -= r0*(l0 @ vector)
        modes, used = [], set()
        for i, value in enumerate(values):
            if i in used or abs(value) < 1e-10*lambda0:
                continue
            group = np.flatnonzero(abs(values-value) < 1e-10*lambda0)
            used.update(map(int, group))
            if 0 in group:
                continue
            rr, ll = right[:, group], left[:, group]
            overlap = ll.T @ rr
            weight = (sink @ rr) @ np.linalg.pinv(overlap, rcond=1e-11) @ (ll.T @ source)
            ratio = value/lambda0
            residual = max(np.linalg.norm(matrix @ rr-rr*values[group])/np.linalg.norm(rr),
                           np.linalg.norm(matrix.T @ ll-ll*values[group])/np.linalg.norm(ll))/lambda0
            modes.append({"eigenvalue": value, "ratio": ratio, "weight": weight,
                          "multiplicity": len(group), "xi": float(-1/np.log(abs(ratio))),
                          "max_relative_residual": float(residual)})
        spectral = np.array([sum(m["weight"]*m["ratio"]**(r-1) for m in modes) for r in distances])
        direct = np.array([connected[r] for r in distances])
        return {"observable": observable, "parity": parity, "sector_dimension": len(matrix),
                "lambda0": lambda0, "mean": mean, "distance": np.array(distances),
                "correlation": direct, "spectral_correlation": spectral, "modes": modes,
                "dominant_right_residual": float(np.linalg.norm(normalized @ r0-r0)/np.linalg.norm(r0)),
                "dominant_left_residual": float(np.linalg.norm(normalized.T @ l0-l0)/np.linalg.norm(l0)),
                "spectral_direct_max_abs_difference": float(np.max(abs(spectral-direct)))}
