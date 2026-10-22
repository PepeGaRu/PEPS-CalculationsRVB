"""Exact sparse transfer matrices of nearest-neighbour singlet RVB PEPS.

No boundary-MPS truncation is used. The exponentially growing circumference
limits this reference implementation to square widths 4/6 and kagome 2/4.
Distances and correlation lengths are measured in axial unit cells. For the
square lattice the tensor is in the Marshall-rotated physical basis; for kagome
one cell contains three physical spins and two triangle simplex tensors.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh


SZ = np.diag([0.5, -0.5])
SP = np.array([[0.0, 1.0], [0.0, 0.0]])
SM = SP.T
IDENTITY = np.eye(2)


def square_site_tensor() -> np.ndarray:
    """Return A[p,left,right,up,down], virtual 2 meaning no dimer.

    Virtual bonds are delta tensors. The physical B-sublattice rotation i*sigma_y
    turns the represented |00>+|11> bonds into oriented A-to-B singlets (up to an
    overall phase). Observables must be transformed by this same rotation.
    """
    a = np.zeros((2, 3, 3, 3, 3))
    for physical in range(2):
        for leg in range(4):
            virtual = [2] * 4
            virtual[leg] = physical
            a[(physical, *virtual)] = 1.0
    return a


def triangle_tensor() -> np.ndarray:
    """T = epsilon_ijk + |222>, with clockwise/counterclockwise fixed by labels."""
    t = np.zeros((3, 3, 3))
    for perm in permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        t[perm] = (-1.0) ** inversions
    t[2, 2, 2] = 1.0
    return t


def kagome_site_tensor() -> np.ndarray:
    """Return K[(a,b,c),left,right,up,down], physical dimension eight.

    Up triangles are (a_xy,b_xy,c_xy), down triangles
    (a_xy,b_(x-1)y,c_x(y-1)), with singlet arrows following each listed cyclic
    order. P^s_ij = delta_i,s delta_j,2 + delta_i,2 delta_j,s enforces one
    incident dimer per physical site. The index 2 is vacuum, not spin S_z=0.
    """
    t = triangle_tensor()
    p = np.zeros((2, 3, 3))
    for spin in range(2):
        p[spin, spin, 2] = p[spin, 2, spin] = 1.0
    k = np.einsum("ald,bce,pab,qcr,teu->pqtlrud", t, t, p, p, p, optimize=True)
    return k.reshape(8, 3, 3, 3, 3)


def double_tensor(a: np.ndarray, operator: np.ndarray | None = None) -> np.ndarray:
    """Local norm/operator double tensor with axes left,right,up,down (size 9)."""
    if operator is None:
        operator = np.eye(a.shape[0])
    result = np.einsum("plrud,pq,qLRUD->lLrRuUdD", a.conj(), operator, a, optimize=True)
    return np.real_if_close(result.reshape((9,) * 4))


def column_transfer(local_tensors: Sequence[np.ndarray]) -> sparse.csr_matrix:
    """Contract a periodic column exactly, enumerating only nonzero entries.

    Rows/columns encode incoming/outgoing double virtual legs, most significant
    digit at y=0. The circumference has len(local_tensors) cells. Repeated
    configurations are summed by scipy, including exact sign cancellations.
    """
    width = len(local_tensors)
    if width < 1:
        raise ValueError("A column must contain at least one cell")
    transitions = []
    for tensor in local_tensors:
        if tensor.shape != (9, 9, 9, 9):
            raise ValueError("Local double tensors must have shape (9,9,9,9)")
        grouped = [[] for _ in range(9)]
        for left, right, up, down in np.argwhere(tensor != 0):
            grouped[int(down)].append((int(left), int(right), int(up), tensor[left, right, up, down]))
        transitions.append(grouped)
    rows, cols, values = [], [], []

    def visit(y, start, vertical, left, right, value):
        for l, r, u, weight in transitions[y][vertical]:
            if y == width - 1:
                if u == start:
                    rows.append(left * 9 + l)
                    cols.append(right * 9 + r)
                    values.append(value * weight)
            else:
                visit(y + 1, start, u, left * 9 + l, right * 9 + r, value * weight)

    for start in range(9):
        visit(0, start, start, 0, 0, 1.0)
    result = sparse.coo_matrix((values, (rows, cols)), shape=(9**width, 9**width)).tocsr()
    result.eliminate_zeros()
    return result


def _largest_eigenpair(matrix: sparse.csr_matrix, tolerance: float):
    # Fixed v0 makes runs reproducible; the largest real eigenvalue is the
    # physical norm eigenvalue. LM would also select negative/complex sectors.
    if matrix.shape[0] < 5:
        vals, vectors = np.linalg.eig(matrix.toarray())
        index = np.argmax(vals.real)
        value, vector = vals[index], vectors[:, index]
    else:
        vals, vectors = eigs(matrix, k=1, which="LR", tol=tolerance,
                             v0=np.ones(matrix.shape[0]), maxiter=10000)
        value, vector = vals[0], vectors[:, 0]
    if abs(value.imag) > 1e-8 * max(1.0, abs(value.real)):
        raise RuntimeError("The physical leading eigenvalue should be real")
    vector = np.real_if_close(vector, tol=1000).real
    vector /= np.linalg.norm(vector)
    return float(value.real), vector


@dataclass
class CorrelationResult:
    lattice: str
    width: int
    distances: np.ndarray
    connected: np.ndarray
    mean_bond: float
    leading_eigenvalue: float
    eigen_residual: float
    transfer_dimension: int
    transfer_nonzeros: int
    selected_sector: str = "dominant"

    def as_dict(self):
        return {
            "lattice": self.lattice,
            "width_unit_cells": self.width,
            "sites_per_axial_slice": self.width * (3 if self.lattice == "kagome" else 1),
            "observable": "connected correlation of physical nearest-neighbour S_i dot S_j",
            "distance_unit": "axial unit cells",
            "distances": self.distances.tolist(),
            "connected": self.connected.tolist(),
            "mean_bond": self.mean_bond,
            "leading_eigenvalue": self.leading_eigenvalue,
            "eigen_residual": self.eigen_residual,
            "transfer_dimension": self.transfer_dimension,
            "transfer_nonzeros": self.transfer_nonzeros,
            "boundary": "infinite axial length; periodic circumference",
            "virtual_vacuum_parity": self.selected_sector if self.lattice == "kagome" else None,
            "vison_inserted": False,
        }


class RVBTransfer:
    """Infinite-cylinder singlet RVB PEPS and its physical bond correlations.

    Width counts square sites or three-spin kagome unit cells around the
    circumference. Even widths avoid odd-spin unit columns and preserve square
    bipartiteness. Square width two has parallel periodic graph edges; use four
    or six for substantive comparisons. A finite cylinder always has an eventual
    exponential envelope; establishing a bulk power law requires width scaling.

    Kagome ``sector='even'`` or ``'odd'`` fixes ket and bra virtual-vacuum
    parity. The default ``'dominant'`` identifies whichever sector has the
    largest norm eigenvalue; that sector can change with circumference. No
    vison line is inserted in any of these calculations.
    """

    def __init__(self, lattice: str, width: int, tolerance: float = 1e-11,
                 sector: str = "dominant"):
        if lattice not in {"square", "kagome"}:
            raise ValueError("lattice must be square or kagome")
        if width < 2 or width % 2:
            raise ValueError("Use a positive even circumference of at least two")
        if (lattice == "square" and width > 6) or (lattice == "kagome" and width > 4):
            raise ValueError("Exact reference solver supports square widths <=6 and kagome <=4")
        if sector not in {"dominant", "even", "odd"}:
            raise ValueError("sector must be dominant, even, or odd")
        if lattice != "kagome" and sector != "dominant":
            raise ValueError("Explicit even/odd virtual-vacuum parity is supported only for kagome")
        self.lattice, self.width = lattice, width
        self.sector_requested = sector
        self.site_tensor = square_site_tensor() if lattice == "square" else kagome_site_tensor()
        self.local_norm = double_tensor(self.site_tensor)
        full_transfer = column_transfer([self.local_norm] * width)
        # Removing completely empty boundary states is exact and saves memory.
        rows, cols = full_transfer.nonzero()
        self.active = np.union1d(rows, cols)
        if sector != "dominant":
            ket_parity, bra_parity = self._vacuum_parities(self.active)
            parity = int(sector == "odd")
            self.active = self.active[(ket_parity == parity) & (bra_parity == parity)]
        self.full_dimension = full_transfer.shape[0]
        self.transfer = full_transfer[self.active][:, self.active].tocsr()
        self.eigenvalue, self.right = _largest_eigenpair(self.transfer, tolerance)
        left_eigenvalue, self.left = _largest_eigenpair(self.transfer.T.tocsr(), tolerance)
        overlap = self.left @ self.right
        if abs(overlap) < 1e-8:
            raise RuntimeError("Left/right fixed points selected incompatible degenerate sectors")
        self.left /= overlap
        self.selected_sector = self._identify_kagome_sector() if lattice == "kagome" else "dominant"
        self.eigen_residual = max(
            np.linalg.norm(self.transfer @ self.right - self.eigenvalue * self.right)
            / (self.eigenvalue * np.linalg.norm(self.right)),
            np.linalg.norm(self.transfer.T @ self.left - self.eigenvalue * self.left)
            / (self.eigenvalue * np.linalg.norm(self.left)),
            abs(left_eigenvalue - self.eigenvalue) / self.eigenvalue,
        )
        self.bond_transfer = self._bond_transfer()
        self.mean_bond = float(self.left @ (self.bond_transfer @ self.right) / self.eigenvalue)

    def _vacuum_parities(self, indices):
        """Ket/bra parity of virtual vacuum 2 in each flattened double index.

        Each base-nine digit is 3 * ket + bra, as fixed by ``double_tensor``.
        This labels the implemented cylindrical tensor geometry. No vison is
        inserted; a match to a paper's physical boundary convention is separate.
        """
        indices = np.asarray(indices).copy()
        ket = np.zeros(len(indices), dtype=np.int8)
        bra = np.zeros(len(indices), dtype=np.int8)
        for _ in range(self.width):
            digit = indices % 9
            indices //= 9
            ket ^= (digit // 3 == 2)
            bra ^= (digit % 3 == 2)
        return ket, bra

    def _identify_kagome_sector(self):
        """Require a common, unmixed ket/bra parity in both fixed points."""
        ket, bra = self._vacuum_parities(self.active)
        selected = []
        for vector in (self.left, self.right):
            weights = np.array([
                np.linalg.norm(vector[(ket == parity) & (bra == parity)]) ** 2
                / np.linalg.norm(vector) ** 2
                for parity in (0, 1)
            ])
            parity = int(np.argmax(weights))
            if weights[parity] < 1 - 1e-8:
                raise RuntimeError("Kagome fixed point mixes virtual parity sectors or has "
                                   "unequal ket/bra parity; request sector='even' or 'odd'")
            selected.append(parity)
        if selected[0] != selected[1]:
            raise RuntimeError("Left/right fixed points have different kagome virtual parities")
        return "odd" if selected[0] else "even"

    def _bond_transfer(self):
        if self.lattice == "kagome":
            # Physical up-triangle bond (a,b) in the y=0 unit cell.
            operator = np.kron(np.kron(SZ, SZ) + 0.5 * (np.kron(SP, SM) + np.kron(SM, SP)), IDENTITY)
            inserted = double_tensor(self.site_tensor, operator)
            full = column_transfer([inserted] + [self.local_norm] * (self.width - 1))
        else:
            # Bond from y=0 to y=1 in the Marshall basis. Rotating exactly
            # one endpoint maps S.S to -SzSz - (S+S+ + S-S-)/2.
            full = None
            for op, coefficient in ((SZ, -1.0), (SP, -0.5), (SM, -0.5)):
                inserted = double_tensor(self.site_tensor, op)
                term = column_transfer([inserted, inserted] + [self.local_norm] * (self.width - 2))
                full = coefficient * term if full is None else full + coefficient * term
        # A bond insertion can enter states on which the norm transfer is
        # identically zero. Adjacent insertions can return from those states;
        # cropping B itself would lose that contact contribution at distance 1.
        # One intervening norm column annihilates them, so longer distances and
        # nonzero transfer eigenmodes can safely use the reduced active space.
        padded_right = np.zeros(self.full_dimension)
        padded_left = np.zeros(self.full_dimension)
        padded_right[self.active], padded_left[self.active] = self.right, self.left
        inserted_right = full @ padded_right / self.eigenvalue
        inserted_left = full.T @ padded_left / self.eigenvalue
        inactive = np.ones(self.full_dimension, dtype=bool)
        inactive[self.active] = False
        self.inactive_contact_correction = float(inserted_left[inactive] @ inserted_right[inactive])
        return full[self.active][:, self.active].tocsr()

    def correlations(self, max_distance: int = 20) -> CorrelationResult:
        if max_distance < 1:
            raise ValueError("max_distance must be at least one")
        normalized_transfer = self.transfer / self.eigenvalue
        normalized_bond = self.bond_transfer / self.eigenvalue
        # Center at both endpoints to avoid subtracting nearly equal mean^2
        # terms when the connected correlation becomes small.
        vector = normalized_bond @ self.right - self.mean_bond * self.right
        left_insert = normalized_bond.T @ self.left - self.mean_bond * self.left
        connected = []
        for distance in range(1, max_distance + 1):
            value = float(left_insert @ vector)
            if distance == 1:
                value += self.inactive_contact_correction
            connected.append(value)
            vector = normalized_transfer @ vector
            vector -= self.right * (self.left @ vector)
        return CorrelationResult(
            self.lattice, self.width, np.arange(1, max_distance + 1),
            np.asarray(connected), self.mean_bond, self.eigenvalue,
            float(self.eigen_residual), self.full_dimension, self.transfer.nnz,
            self.selected_sector,
        )

    def dimer_correlation_length(self, amplitude_threshold: float = 1e-10):
        """Select the leading eigenmode with nonzero physical dimer weight.

        Topological sectors and spin-carrying eigenmodes can have a longer
        generic transfer length but zero overlap with a local scalar bond.
        Left/right degenerate eigenspaces are paired before evaluating weights.
        A complex pair gives an oscillating exponential envelope. The result is
        a finite-cylinder length, not an infinite-plane extrapolation.
        """
        transfer = self.transfer / self.eigenvalue
        bond = self.bond_transfer / self.eigenvalue
        right_source = bond @ self.right - self.mean_bond * self.right
        left_source = bond.T @ self.left - self.mean_bond * self.left
        dimension = transfer.shape[0]
        symmetric = self.lattice == "square"
        if dimension <= 128:
            if symmetric:
                eigenvalues, right_modes = np.linalg.eigh(transfer.toarray())
                left_values, left_modes = eigenvalues, right_modes
            else:
                eigenvalues, right_modes = np.linalg.eig(transfer.toarray())
                left_values, left_modes = np.linalg.eig(transfer.toarray().T)
        elif symmetric:
            eigenvalues, right_modes = eigsh(transfer, k=24, which="LM", tol=1e-11,
                                             v0=right_source)
            left_values, left_modes = eigenvalues, right_modes
        else:
            eigenvalues, right_modes = eigs(transfer, k=100, which="LM", tol=1e-11,
                                            v0=np.ones(dimension), maxiter=10000)
            left_values, left_modes = eigs(transfer.T, k=100, which="LM", tol=1e-11,
                                           v0=np.ones(dimension), maxiter=10000)
        visited = set()
        spectral_floor = min(abs(eigenvalues))
        unmatched_larger_groups = 0
        for index in np.argsort(-np.abs(eigenvalues)):
            if int(index) in visited:
                continue
            value = eigenvalues[index]
            right_indices = np.flatnonzero(abs(eigenvalues - value) < 1e-7)
            left_indices = np.flatnonzero(abs(left_values - value) < 1e-7)
            visited.update(map(int, right_indices))
            if len(left_indices) != len(right_indices):
                unmatched_larger_groups += 1
                continue
            r_modes, l_modes = right_modes[:, right_indices], left_modes[:, left_indices]
            overlap = l_modes.T @ r_modes
            weight = ((left_source @ r_modes) @ np.linalg.pinv(overlap, rcond=1e-10)
                      @ (l_modes.T @ right_source))
            if abs(weight) <= amplitude_threshold:
                continue
            if abs(value) >= 1 - 1e-8:
                return {"xi": None, "reason": "a nondecaying transfer mode has nonzero connected bond weight",
                        "eigenvalue_ratio_real": float(value.real),
                        "eigenvalue_ratio_imag": float(np.imag(value)),
                        "spectral_weight_magnitude": float(abs(weight))}
            if unmatched_larger_groups:
                return {"xi": None, "reason": "unmatched larger left/right eigenspaces prevent a completeness claim"}
            if dimension > 128 and abs(value) <= spectral_floor * (1 + 1e-4):
                return {"xi": None, "reason": "observable mode lies at the computed spectral cutoff"}
            residual = np.linalg.norm(transfer @ r_modes - r_modes * eigenvalues[right_indices])
            residual /= max(np.linalg.norm(r_modes), 1e-30)
            return {
                "xi": float(-1.0 / np.log(abs(value))),
                "method": "leading transfer eigenmode with nonzero connected physical-bond spectral weight",
                "eigenvalue_ratio_real": float(value.real),
                "eigenvalue_ratio_imag": float(np.imag(value)),
                "spectral_weight_real": float(weight.real),
                "spectral_weight_imag": float(np.imag(weight)),
                "mode_eigen_residual": float(residual),
                "spectral_weight_threshold": amplitude_threshold,
                "examined_modes": len(eigenvalues),
                "smallest_computed_eigenvalue_magnitude": float(spectral_floor),
                "larger_eigenspaces_matched": True,
                "units": "axial unit cells",
            }
        return {"xi": None, "reason": "no coupled decaying mode resolved in computed spectrum"}


def exponential_tail_length(distances, values, relative_floor=1e-9):
    """Fit a finite-cylinder exponential tail, never a bulk critical exponent.

    Uses the last six points above a conservative absolute/relative roundoff
    floor. Oscillations, competing modes and short cylinders can invalidate a
    single exponential; the returned residual exposes a poor fit.
    """
    r, amplitude = np.asarray(distances), np.abs(np.asarray(values))
    threshold = max(1e-13, float(amplitude.max(initial=0)) * relative_floor)
    selected = np.flatnonzero((amplitude > threshold) & (r >= 2))[-6:]
    if len(selected) < 3:
        return {"xi": None, "reason": "too few correlations above numerical floor"}
    x, y = r[selected], np.log(amplitude[selected])
    slope, intercept = np.polyfit(x, y, 1)
    error = np.sqrt(np.mean((y - (intercept + slope * x)) ** 2))
    if slope >= 0:
        return {"xi": None, "reason": "tail is not decaying", "log_rms_residual": float(error)}
    return {"xi": float(-1.0 / slope), "fit_distances": x.tolist(),
            "log_rms_residual": float(error), "method": "finite-cylinder exponential tail fit"}


def rvb_cylinder_benchmark(lattice: str, width: int, max_distance: int = 32,
                           sector: str = "dominant"):
    """JSON-compatible exact singlet-RVB benchmark for the command-line driver."""
    model = RVBTransfer(lattice, width, sector=sector)
    result = model.correlations(max_distance)
    fit = exponential_tail_length(result.distances, result.connected)
    spectral_length = model.dimer_correlation_length()
    return {
        "lattice": lattice,
        "width": width,
        "circumference": width,
        "distance": result.distances.tolist(),
        "connected_correlation": result.connected.tolist(),
        "mean_bond_energy": result.mean_bond,
        "eigen_residual": result.eigen_residual,
        "correlation_length": spectral_length["xi"],
        "correlation_length_spectrum": spectral_length,
        "correlation_length_fit": fit,
        "metadata": {
            "method": "exact sparse D=3 singlet-RVB PEPS transfer matrix",
            "operator": "physical S_i dot S_j; connected correlation between translated disjoint bonds",
            "bond": "circumferential y=0,y=1" if lattice == "square" else "up-triangle a,b at y=0",
            "units": "axial unit cells",
            "circumference_units": "sites" if lattice == "square" else "three-spin kagome unit cells",
            "sector": (f"{model.selected_sector} ket/bra virtual-vacuum parity" if lattice == "kagome"
                       else "dominant norm sector selected by positive largest-real eigenvalue"),
            "sector_requested": sector,
            "virtual_vacuum_parity": model.selected_sector if lattice == "kagome" else None,
            "vison_inserted": False,
            "sector_convention": "Implemented virtual-boundary charge; matching a paper's physical boundary geometry requires a separate mapping.",
            "convergence_residual": result.eigen_residual,
            "leading_eigenvalue": result.leading_eigenvalue,
            "transfer_dimension": result.transfer_dimension,
            "transfer_nonzeros": result.transfer_nonzeros,
            "boundary": "infinite axial direction; periodic circumference",
            "limitation": "Finite widths do not establish a two-dimensional exponent or bulk correlation length.",
        },
    }
