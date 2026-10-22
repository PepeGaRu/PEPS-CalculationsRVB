"""Exact longitudinal transfer for a 2 by 2 periodic spin-1 octahedral tube.

The physical state and normalized singlet phases are those of
:mod:`rvb.octahedral_spin1`. No covering overlaps, tensor singular values,
or virtual sectors are discarded. A shared physical edge has a rank-eight
norm metric on its sixteen bra/ket virtual pairs. Projecting onto its full
support reduces the slice boundary to 8**4=4096 states. Sparse contraction
then gives only about thirty thousand nonzero transfer entries.

The transverse period-two identifications remain; this is an infinite-length
finite-width tube, not a three-dimensional thermodynamic-limit calculation.
Physical observables are inserted on the two y-directed sites of the
periodic line at z=0, before compressing their internal virtual pairs. The
unobserved x-cut edges are compressed exactly. This distinction is essential:
projecting an observable onto the norm support on its own edge loses terms.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import numpy as np
import opt_einsum as oe
from scipy import linalg, sparse

from .octahedral_spin1 import octahedron_tensor, ownership_projection, spin1_dimer, spin1_generators


EDGE_PERMUTATION = np.array([1, 0, 5, 6, 7, 2, 3, 4])
EDGE_WEIGHTS = np.array([np.sqrt(3), np.sqrt(3), 1., 1., 1., 1., 1., 1.])


def norm_edge_basis():
    """Orthonormal columns: 00, trace(occupied)/sqrt3, s0, and 0s."""
    basis = np.zeros((16, 8))
    basis[0, 0] = 1
    for s in range(1, 4):
        basis[4*s+s, 1] = 1/np.sqrt(3)
        basis[4*s, s+1] = 1
        basis[s, s+4] = 1
    return basis


def full_edge_metric(operator=None):
    """Physical one-site operator contracted into the two ownership tensors.

    Rows and columns group ket/bra virtual indices at each endpoint.
    """
    if operator is None:
        operator = np.eye(3)
    projection = ownership_projection()
    return np.einsum("pac,qbd,qp->abcd", projection, projection, operator).reshape(16, 16)


def boundary_sectors():
    """Label bra and ket ownership patterns modulo simultaneous complement.

    The four sites are in lexicographic (y,z) order. Bit significance follows
    their array order, so the first site is the least significant bit here.
    Each individual single-layer orbit is min(pattern,15-pattern) in 0..7.
    """
    indices = np.array(np.unravel_index(np.arange(8**4), (8,)*4)).T
    ket = np.sum(((indices == 1) | ((indices >= 2) & (indices <= 4))) * 2**np.arange(4), axis=1)
    bra = np.sum(((indices == 1) | (indices >= 5)) * 2**np.arange(4), axis=1)
    return np.minimum(bra, 15-bra)*8 + np.minimum(ket, 15-ket)


class Spin1OctahedralCylinder:
    """Exact transfer on an infinite x tube with a periodic 2x2 y,z section.

    ``trimer_weight`` multiplies each normalized trimer amplitude. The user
    state has its default value one; zero supplies an analytic test limit.
    """
    transverse_shape = (2, 2)
    sites_per_slice = 12
    boundary_dimension = 4096

    def __init__(self, trimer_weight=1.0):
        if not np.isfinite(trimer_weight) or trimer_weight < 0:
            raise ValueError("Trimer amplitude must be finite and nonnegative")
        self.trimer_weight = float(trimer_weight)
        tensor = octahedron_tensor()
        occupied = np.count_nonzero(np.indices(tensor.shape), axis=0)
        tensor[occupied == 3] *= self.trimer_weight
        self._single_indices = np.array(tensor.nonzero()).T
        self._single_values = tensor[tuple(self._single_indices.T)]
        self._local_cache = {}
        self._pair_cache = {}
        self.sector_labels = boundary_sectors()
        self.transfer = self._assemble_transfer()

    def _local_entries(self, full_y=False):
        if full_y in self._local_cache:
            return self._local_cache[full_y]
        indices, values = self._single_indices, self._single_values
        ket, bra = indices[:, None, :], indices[None, :, :]
        projected_legs = (0, 1, 4, 5) if full_y else tuple(range(6))
        valid = np.all(((ket == 0) | (bra == 0) | (ket == bra))[:, :, projected_legs], axis=-1)
        first, second = np.where(valid)
        ket, bra = indices[first], indices[second]
        result = np.zeros(ket.shape, dtype=int)
        amplitudes = values[first]*values[second]
        for leg in range(6):
            if full_y and leg in (2, 3):
                result[:, leg] = 4*ket[:, leg]+bra[:, leg]
                continue
            both = (ket[:, leg] > 0) & (bra[:, leg] > 0)
            ket_only = (ket[:, leg] > 0) & (bra[:, leg] == 0)
            bra_only = (ket[:, leg] == 0) & (bra[:, leg] > 0)
            result[both, leg] = 1
            result[ket_only, leg] = ket[ket_only, leg]+1
            result[bra_only, leg] = bra[bra_only, leg]+4
            amplitudes[both] /= np.sqrt(3)
        # Absorb the sparse norm metric on positive-directed legs. On full-y
        # legs the physical two-site connector below supplies both metrics.
        for leg in ((0, 4) if full_y else (0, 2, 4)):
            amplitudes *= EDGE_WEIGHTS[result[:, leg]]
            result[:, leg] = EDGE_PERMUTATION[result[:, leg]]
        self._local_cache[full_y] = result, amplitudes
        return result, amplitudes

    def _cell_matrix(self, y, z, full_y=False):
        indices, amplitudes = self._local_entries(full_y)
        indices = indices.copy()
        if y:
            indices[:, [2, 3]] = indices[:, [3, 2]]
        if z:
            indices[:, [4, 5]] = indices[:, [5, 4]]
        physical = indices[:, 0]*8+indices[:, 1]
        z_pair = indices[:, 4]*8+indices[:, 5]
        y_dimension = 16 if full_y else 8
        y_pair = indices[:, 2]*y_dimension+indices[:, 3]
        result = sparse.coo_matrix((amplitudes, (physical*64+z_pair, y_pair)),
                                   shape=(4096, y_dimension**2)).tocsr()
        result.eliminate_zeros()
        return result

    @staticmethod
    def _two_site_connector(operator):
        """Exact physical operator on the two sites of one y-periodic line."""
        operator = np.asarray(operator)
        if operator.shape != (9, 9):
            raise ValueError("The physical two-spin operator must be 9 by 9")
        p = ownership_projection()
        # Ket local ports u,v and remote t,w; bra ports u',v',t',w'.
        value = oe.contract(p, [0, 4, 8], p, [1, 6, 10], p, [2, 5, 9], p, [3, 7, 11],
                            operator.reshape(3, 3, 3, 3), [2, 3, 0, 1],
                            [4, 5, 6, 7, 8, 9, 10, 11], optimize="optimal").reshape(256, 256)
        return sparse.csr_matrix(value)

    def _pair_matrix(self, z, operator=None):
        if operator is None and z in self._pair_cache:
            return self._pair_cache[z]
        if operator is None:
            matrix = (self._cell_matrix(0, z) @ self._cell_matrix(1, z).T).tocoo()
        else:
            matrix = (self._cell_matrix(0, z, True) @ self._two_site_connector(operator)
                      @ self._cell_matrix(1, z, True).T).tocoo()
        p0, z0 = divmod(matrix.row, 64)
        p1, z1 = divmod(matrix.col, 64)
        result = sparse.coo_matrix((matrix.data, (p0*64+p1, z0*64+z1)), shape=(4096, 4096)).tocsr()
        result.eliminate_zeros()
        if operator is None:
            self._pair_cache[z] = result
        return result

    def _assemble_transfer(self, operator=None):
        matrix = (self._pair_matrix(0, operator) @ self._pair_matrix(1).T).tocoo()
        matrix.eliminate_zeros()
        p00, p10 = divmod(matrix.row, 64)
        p01, p11 = divmod(matrix.col, 64)
        outgoing = np.zeros(matrix.nnz, dtype=int)
        incoming = np.zeros(matrix.nnz, dtype=int)
        for pair in (p00, p01, p10, p11):
            outgoing = outgoing*8+pair//8
            incoming = incoming*8+pair % 8
        result = sparse.coo_matrix((matrix.data, (outgoing, incoming)), shape=(4096, 4096)).tocsr()
        result.eliminate_zeros()
        return result

    @lru_cache(maxsize=None)
    def _expanded_longitudinal_slice(self):
        """One x line uses its full 16-dimensional bra/ket pairs internally.

        Return raw U, incoming basis C_in and outgoing metric C_out so that
        the ordinary one-slice transfer is C_out @ U @ C_in.
        """
        indices, values = self._single_indices, self._single_values
        ket, bra = indices[:, None, :], indices[None, :, :]
        valid = np.all(((ket == 0) | (bra == 0) | (ket == bra))[:, :, (2, 3, 4, 5)], axis=-1)
        first, second = np.where(valid)
        ket, bra = indices[first], indices[second]
        labels = np.zeros(ket.shape, dtype=int)
        amplitudes = values[first]*values[second]
        for leg in range(6):
            if leg in (0, 1):
                labels[:, leg] = 4*ket[:, leg]+bra[:, leg]
                continue
            both = (ket[:, leg] > 0) & (bra[:, leg] > 0)
            ket_only = (ket[:, leg] > 0) & (bra[:, leg] == 0)
            bra_only = (ket[:, leg] == 0) & (bra[:, leg] > 0)
            labels[both, leg] = 1
            labels[ket_only, leg] = ket[ket_only, leg]+1
            labels[bra_only, leg] = bra[bra_only, leg]+4
            amplitudes[both] /= np.sqrt(3)
        for leg in (2, 4):
            amplitudes *= EDGE_WEIGHTS[labels[:, leg]]
            labels[:, leg] = EDGE_PERMUTATION[labels[:, leg]]
        physical = labels[:, 0]*16+labels[:, 1]
        z_pair = labels[:, 4]*8+labels[:, 5]
        y_pair = labels[:, 2]*8+labels[:, 3]
        cell = sparse.coo_matrix((amplitudes, (physical*64+z_pair, y_pair)),
                                 shape=(16384, 64)).tocsr()
        pair = (cell @ self._cell_matrix(1, 0).T).tocoo()
        p0, z0 = divmod(pair.row, 64)
        p1, z1 = divmod(pair.col, 64)
        pair = sparse.coo_matrix((pair.data, (p0*64+p1, z0*64+z1)), shape=(16384, 4096)).tocsr()
        plane = (pair @ self._pair_matrix(1).T).tocoo()
        p00, p10 = divmod(plane.row, 64)
        p01, p11 = divmod(plane.col, 64)
        outgoing, incoming = p00//16, p00 % 16
        for p in (p01, p10, p11):
            outgoing = outgoing*8+p//8
            incoming = incoming*8+p % 8
        raw = sparse.coo_matrix((plane.data, (outgoing, incoming)), shape=(8192, 8192)).tocsr()
        raw.eliminate_zeros()
        basis = norm_edge_basis()
        incoming = sparse.kron(sparse.csr_matrix(basis), sparse.eye(512), format="csr")
        outgoing = sparse.kron(sparse.csr_matrix(basis.T @ full_edge_metric().T),
                               sparse.eye(512), format="csr")
        return raw, incoming, outgoing

    @lru_cache(maxsize=None)
    def longitudinal_projector_transfer(self):
        """Three-slice insertion for the x-diameter singlet projector.

        The two physical x sites lie on either side of the middle slice's
        octahedron at y=z=0. P0=(1/3) sum_ab |a><b| tensor |a><b|. Keep
        both observed bonds' full virtual supports, then compress only the
        two outer unobserved bonds. Products of these three-slice blocks
        represent correlations at centre separation r>=3.
        """
        raw, incoming, outgoing = self._expanded_longitudinal_slice()
        left, right = outgoing @ raw, raw @ incoming
        result = sparse.csr_matrix((4096, 4096))
        for a, b in product(range(3), repeat=2):
            elementary = np.zeros((3, 3))
            elementary[a, b] = 1
            connector = sparse.kron(sparse.csr_matrix(full_edge_metric(elementary).T),
                                    sparse.eye(512), format="csr")
            result += left @ connector @ raw @ connector @ right/3
        result.eliminate_zeros()
        return result

    def physical_transfer(self, operator):
        """Transfer with a chosen two-spin physical operator on y sites at z=0."""
        return self._assemble_transfer(operator)

    @lru_cache(maxsize=None)
    def observable_transfer(self, observable="dimer_projector"):
        if observable == "longitudinal_dimer_projector":
            return self.longitudinal_projector_transfer()
        if observable == "dimer_projector":
            singlet = spin1_dimer().ravel()
            operator = np.outer(singlet, singlet)
        elif observable == "spin_z":
            operator = np.kron(spin1_generators()[2], np.eye(3))
        elif observable == "identity":
            operator = np.eye(9)
        else:
            raise ValueError("Observable must be dimer_projector, longitudinal_dimer_projector, spin_z, or identity")
        return self.physical_transfer(operator)

    def sector_indices(self, sector):
        if not isinstance(sector, (int, np.integer)) or not 0 <= sector < 64:
            raise ValueError("Sector label must be an integer from 0 through 63")
        return np.flatnonzero(self.sector_labels == sector)

    @lru_cache(maxsize=None)
    def sector_spectra(self):
        result = []
        for sector in range(64):
            indices = self.sector_indices(sector)
            block = self.transfer[indices, :][:, indices].toarray()
            values = linalg.eigvals(block)
            values = values[np.argsort(-abs(values))]
            result.append({"sector": sector, "bra_orbit": sector//8, "ket_orbit": sector % 8,
                           "dimension": len(indices), "eigenvalues": values})
        return result

    def dominant_sector(self):
        return max(self.sector_spectra(), key=lambda item: abs(item["eigenvalues"][0]))["sector"]

    def analyze(self, observable="dimer_projector", distances=range(1, 65), sector=None):
        """Observable-resolved eigenmodes and centered transfer propagation.

        C(r)=sum_j w_j (lambda_j/lambda_0)**(r-L), r>=L, with L=1
        for transverse insertions and L=3 for the longitudinal projector.
        The mean term is
        removed algebraically, so the small tail is not obtained by
        subtracting two nearly equal order-one expectation values.
        """
        distances = tuple(int(r) for r in distances)
        insertion_slices = 3 if observable == "longitudinal_dimer_projector" else 1
        if not distances or min(distances) < insertion_slices:
            raise ValueError(f"Use centre separations at least {insertion_slices} for this insertion")
        if sector is None:
            sector = self.dominant_sector()
        indices = self.sector_indices(sector)
        matrix = self.transfer[indices, :][:, indices].toarray()
        values, right = linalg.eig(matrix)
        order = np.argsort(-abs(values))
        values, right = values[order], right[:, order]
        left = linalg.inv(right)
        lambda0 = values[0]
        if abs(lambda0.imag) > 1e-10 or lambda0.real <= 0:
            raise ValueError("Selected sector does not have the required positive leading eigenvalue")
        lambda0 = lambda0.real
        operator = self.observable_transfer(observable)[indices, :][:, indices].toarray()/lambda0**insertion_slices
        transformed = left @ operator @ right
        mean = transformed[0, 0]
        weights = transformed[0, :]*transformed[:, 0]
        weights[0] = 0
        normalized = matrix/lambda0
        r0, l0 = right[:, 0], left[0, :]
        vector = operator @ r0-mean*r0
        row = l0 @ operator-mean*l0
        propagated = {}
        for distance in range(insertion_slices, max(distances)+1):
            if distance in distances:
                propagated[distance] = row @ vector
            vector = normalized @ vector
            vector -= r0*(l0 @ vector)
        residual = np.linalg.norm(matrix @ right-right*values[None, :], axis=0)
        residual /= np.linalg.norm(matrix)*np.linalg.norm(right, axis=0)
        groups = []
        used = set()
        for j, value in enumerate(values):
            if j in used or j == 0:
                continue
            members = np.flatnonzero(abs(values-value) < 1e-10*lambda0)
            members = members[members != 0]
            used.update(map(int, members))
            weight = np.sum(weights[members])
            ratio = value/lambda0
            xi = float(-1/np.log(abs(ratio))) if 0 < abs(ratio) < 1 else None
            groups.append({"eigenvalue": value, "ratio": ratio, "weight": weight,
                           "multiplicity": len(members), "xi": xi,
                           "max_relative_residual": float(np.max(residual[members]))})
        spectral = np.array([np.sum(weights*(values/lambda0)**(r-insertion_slices)) for r in distances])
        direct = np.array([propagated[r] for r in distances])
        return {"observable": observable, "insertion_slices": insertion_slices,
                "sector": sector, "sector_dimension": len(indices),
                "lambda0": lambda0, "mean": mean, "distance": np.array(distances),
                "correlation": direct, "spectral_correlation": spectral,
                "modes": groups, "eigenvector_condition": float(np.linalg.cond(right)),
                "max_eigen_relative_residual": float(np.max(residual)),
                "spectral_direct_max_abs_difference": float(np.max(abs(spectral-direct)))}
