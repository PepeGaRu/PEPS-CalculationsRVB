"""Exact finite-cylinder physical channels of the square-octagon singlet RVB.

There is no CTMRG or boundary-bond truncation.  The circumference remains finite,
so even a critical two-dimensional state can have exponentially decaying
cylinder correlations.  Ket and bra virtual-vacuum parity are fixed separately
to the same requested value; no vison line is inserted.
"""

from functools import lru_cache
from time import perf_counter

import numpy as np
from scipy.linalg import eig
from scipy.sparse.linalg import eigs

from .octagon_peps import (octagon_local_tensor, physical_bond_operator,
                           physical_spin_operator)
from .peps import column_transfer, double_tensor, exponential_tail_length


OBSERVABLES = {
    "ring_ring": {"sites": (0, 3), "kind": "spin_dot"},
    "ring_dummy_0": {"sites": (0, 4), "kind": "spin_dot"},
    "ring_dummy_3": {"sites": (3, 4), "kind": "spin_dot"},
    "spin_ring_z": {"sites": (0,), "kind": "spin_z"},
    "spin_dummy_z": {"sites": (4,), "kind": "spin_z"},
}


def boundary_vacuum_parities(indices, width):
    """Parity of virtual value 2, separately for the two copies of the PEPS."""
    remaining = np.asarray(indices, dtype=np.int64).copy()
    ket, bra = np.zeros(remaining.shape, dtype=np.int8), np.zeros(remaining.shape, dtype=np.int8)
    for _ in range(width):
        digit = remaining % 9
        remaining //= 9
        ket ^= (digit // 3 == 2)
        bra ^= (digit % 3 == 2)
    return ket, bra


@lru_cache(maxsize=4)
def _full_norm_transfer(width):
    return column_transfer([double_tensor(octagon_local_tensor())] * width)


@lru_cache(maxsize=20)
def _full_observable_transfer(width, name):
    specification = OBSERVABLES[name]
    operator = (physical_spin_operator(specification["sites"][0], "z")
                if specification["kind"] == "spin_z" else
                physical_bond_operator(*specification["sites"], kind="spin_dot"))
    tensor = octagon_local_tensor()
    return column_transfer([double_tensor(tensor, operator)]
                           + [double_tensor(tensor)] * (width - 1))


class OctagonCylinder:
    """Infinite x direction and periodic y direction, at a specified parity."""

    def __init__(self, width, sector="even", tolerance=1e-12):
        if isinstance(width, bool) or not isinstance(width, (int, np.integer)) or not 1 <= width <= 4:
            raise ValueError("Exact cylinder width must be an integer from 1 through 4")
        if sector not in ("even", "odd"):
            raise ValueError("Specify the even or odd ket/bra vacuum-parity sector")
        self.width, self.sector = int(width), sector
        full = _full_norm_transfer(self.width)
        ket, bra = boundary_vacuum_parities(np.arange(full.shape[0]), self.width)
        parity = int(sector == "odd")
        self.indices = np.flatnonzero((ket == parity) & (bra == parity))
        self.transfer = full[self.indices][:, self.indices].tocsr()
        # Retain zero rows and columns within this sector. Observable insertions
        # can visit them at separation one, although the norm transfer kills them.
        if self.transfer.shape[0] <= 4:
            vals, vectors = np.linalg.eig(self.transfer.toarray())
            index = np.argmax(vals.real)
            value, right = vals[index], vectors[:, index]
            vals_l, vectors_l = np.linalg.eig(self.transfer.toarray().T)
            index_l = np.argmin(abs(vals_l - value))
            left_value, left = vals_l[index_l], vectors_l[:, index_l]
        else:
            options = dict(k=1, which="LR", tol=tolerance,
                           v0=np.ones(len(self.indices)), maxiter=20000)
            vals, vectors = eigs(self.transfer, **options)
            vals_l, vectors_l = eigs(self.transfer.T, **options)
            value, right = vals[0], vectors[:, 0]
            left_value, left = vals_l[0], vectors_l[:, 0]
        if abs(np.imag(value)) > tolerance * abs(value):
            raise RuntimeError("The dominant norm eigenvalue is not real")
        self.eigenvalue = float(np.real(value))
        self.right, self.left = np.real_if_close(right).real, np.real_if_close(left).real
        overlap = self.left @ self.right
        if abs(overlap) < 1e-10:
            raise RuntimeError("Incompatible left/right dominant eigenspaces")
        self.left /= overlap
        self.normalized = self.transfer / self.eigenvalue
        self.eigen_residual = float(max(
            np.linalg.norm(self.normalized @ self.right - self.right) / np.linalg.norm(self.right),
            np.linalg.norm(self.normalized.T @ self.left - self.left) / np.linalg.norm(self.left),
            abs(left_value - value) / abs(value),
        ))
        self._spectrum = None

    def observable_transfer(self, name):
        if name not in OBSERVABLES:
            raise ValueError(f"Unknown observable {name!r}; choose {tuple(OBSERVABLES)}")
        full = _full_observable_transfer(self.width, name)
        return full[self.indices][:, self.indices].tocsr()

    def sources(self, name):
        inserted = self.observable_transfer(name) / self.eigenvalue
        mean = float(self.left @ (inserted @ self.right))
        return (mean, inserted @ self.right - mean * self.right,
                inserted.T @ self.left - mean * self.left)

    def correlations(self, name, max_distance=32):
        """Compute centered insertions directly, without subtracting two means²."""
        if max_distance < 1:
            raise ValueError("max_distance must be positive")
        mean, vector, left_source = self.sources(name)
        connected = []
        for _ in range(max_distance):
            connected.append(float(left_source @ vector))
            vector = self.normalized @ vector
            vector -= self.right * (self.left @ vector)
        return mean, np.asarray(connected)

    def spectrum(self):
        """All eigenmodes, including left/right vectors; no spectral cutoff."""
        if self._spectrum is None:
            start = perf_counter()
            values, left, right = eig(self.normalized.toarray(), left=True, right=True)
            # scipy gives adjoint-left vectors; below we consistently use transpose.
            self._spectrum = values, left.conj(), right, perf_counter() - start
        return self._spectrum

    def channel_spectrum(self, name, weight_threshold=1e-10):
        """Select complete degenerate eigenspaces coupled to this observable."""
        _, right_source, left_source = self.sources(name)
        values, left, right, seconds = self.spectrum()
        visited, modes = set(), []
        for index in np.argsort(-abs(values)):
            if int(index) in visited:
                continue
            value = values[index]
            # Exact zero eigenmodes can give only distance-one contacts.
            if abs(value) < 1e-9:
                break
            group = np.flatnonzero(abs(values - value) < 1e-8)
            visited.update(map(int, group))
            rr, ll = right[:, group], left[:, group]
            overlap = ll.T @ rr
            weight = ((left_source @ rr) @ np.linalg.pinv(overlap, rcond=1e-11)
                      @ (ll.T @ right_source))
            if abs(weight) <= weight_threshold:
                continue
            right_residual = (np.linalg.norm(self.normalized @ rr - rr * values[group])
                              / np.linalg.norm(rr))
            left_residual = (np.linalg.norm(self.normalized.T @ ll - ll * values[group])
                             / np.linalg.norm(ll))
            magnitude = float(abs(value))
            modes.append({
                "eigenvalue_ratio_real": float(value.real),
                "eigenvalue_ratio_imag": float(value.imag),
                "eigenvalue_ratio_magnitude": magnitude,
                "xi": float(-1 / np.log(magnitude)) if magnitude < 1 - 1e-8 else None,
                "spectral_weight_real": float(weight.real),
                "spectral_weight_imag": float(weight.imag),
                "spectral_weight_magnitude": float(abs(weight)),
                "eigenspace_dimension": len(group),
                "right_mode_residual": float(right_residual),
                "left_mode_residual": float(left_residual),
            })
        return {
            "xi": modes[0]["xi"] if modes else None,
            "leading_coupled_mode": modes[0] if modes else None,
            "coupled_modes": modes,
            "weight_threshold": weight_threshold,
            "examined_modes": len(values),
            "full_sector_spectrum": True,
            "diagonalization_seconds": seconds,
            "method": "largest-magnitude nonzero transfer eigenspace with connected physical-observable weight",
            "limitation": "Exact finite-cylinder channel; not a two-dimensional correlation length",
        }

    def benchmark(self, max_distance=32, names=None):
        names = tuple(OBSERVABLES) if names is None else tuple(names)
        channels = {}
        distances = np.arange(1, max_distance + 1)
        for name in names:
            mean, connected = self.correlations(name, max_distance)
            spectrum = self.channel_spectrum(name)
            reconstructed = np.zeros(max_distance, dtype=complex)
            for mode in spectrum["coupled_modes"]:
                ratio = complex(mode["eigenvalue_ratio_real"], mode["eigenvalue_ratio_imag"])
                weight = complex(mode["spectral_weight_real"], mode["spectral_weight_imag"])
                reconstructed += weight * ratio**np.arange(max_distance)
            spectrum["correlation_reconstruction_max_error_r_ge_2"] = float(
                np.max(abs(reconstructed[1:] - connected[1:]), initial=0))
            channels[name] = {
                "observable": OBSERVABLES[name],
                "mean": mean,
                "distance": distances.tolist(),
                "connected_correlation": connected.tolist(),
                "correlation_length_spectrum": spectrum,
                "correlation_length_fit": exponential_tail_length(distances, connected),
                "correlation_absolute_plot_floor": 1e-13,
                "resolved_distances_above_plot_floor": distances[abs(connected) > 1e-13].tolist(),
            }
        return {
            "lattice": "triangle-decorated square-octagon",
            "state": "equal-coefficient oriented-singlet RVB; no semionic loop sign",
            "width": self.width,
            "sector": self.sector,
            "channels": channels,
            "metadata": {
                "method": "exact D=3 cylinder transfer; complete eigenspectrum within fixed sector",
                "sector": "same virtual-vacuum parity in ket and bra",
                "vison_inserted": False,
                "axial_direction": "x",
                "circumference_direction": "y",
                "distance_units": "original square-lattice unit cells",
                "sites_per_circumference": 8 * self.width,
                "leading_eigenvalue": self.eigenvalue,
                "dominant_left_right_residual": self.eigen_residual,
                "sector_dimension": len(self.indices),
                "transfer_nonzeros": self.transfer.nnz,
                "physical_sites": "r0,r1,r2,r3,dE,dN,dW,dS",
                "limitation": "Finite circumference; no conclusion about bulk criticality from a single width",
            },
        }


def main():
    """Regenerate the sector-resolved data without editing a Python script."""
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", nargs="+", type=int, default=[2, 3, 4])
    parser.add_argument("--sectors", nargs="+", choices=("even", "odd"), default=["even", "odd"])
    parser.add_argument("--max-distance", type=int, default=32)
    parser.add_argument("--output", type=Path, default=Path("results/octagon-square/cylinders"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for width in args.widths:
        for sector in args.sectors:
            result = OctagonCylinder(width, sector).benchmark(args.max_distance)
            destination = args.output / f"width{width}-{sector}.json"
            destination.write_text(json.dumps(result, indent=2) + "\n")
            print(destination, flush=True)


if __name__ == "__main__":
    main()
