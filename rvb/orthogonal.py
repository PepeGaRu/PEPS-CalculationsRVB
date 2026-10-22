"""Exact benchmarks for orthogonal nearest-neighbour dimer states.

Orthogonal coverings have a diagonal norm, so their equal-amplitude state has
the same dimer-occupation statistics as the uniform classical dimer model.
These calculations do not replace the nonorthogonal spin-singlet RVB norm.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Iterable

import numpy as np
from scipy.integrate import quad
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class SquareDimerCylinder:
    """Exact column transfer matrix on an infinite square-lattice cylinder.

    ``width`` is the periodic circumference, in lattice spacings. The matrix
    acts on horizontal dimers crossing a column boundary. We use the zero
    staggered-flux block, which contains the dominant eigenvector for even
    circumference. Each unoccupied site is paired vertically inside a column
    or horizontally to the next column. No length or sampling approximation
    enters the transfer construction.
    """

    width: int

    def __post_init__(self) -> None:
        if isinstance(self.width, bool) or self.width < 4 or self.width % 2:
            raise ValueError("width must be an even integer >= 4")
        if self.width > 16:
            raise ValueError("width > 16 is outside this exact benchmark's resource limit")

    @cached_property
    def states(self) -> np.ndarray:
        return np.array([
            mask for mask in range(1 << self.width)
            if sum(((mask >> y) & 1) * (-1 if y % 2 else 1)
                   for y in range(self.width)) == 0
        ], dtype=np.int64)

    @cached_property
    def _operators(self) -> tuple[csr_matrix, csr_matrix]:
        indices = {int(mask): i for i, mask in enumerate(self.states)}
        rows, columns, values, vertical_values = [], [], [], []
        full = (1 << self.width) - 1

        for incoming, row in indices.items():
            def visit(occupied: int, outgoing: int, vertical: int) -> None:
                if occupied == full:
                    rows.append(row)
                    columns.append(indices[outgoing])
                    values.append(1.0)
                    vertical_values.append(float(vertical & 1))
                    return
                unoccupied = full ^ occupied
                bit = unoccupied & -unoccupied
                y = bit.bit_length() - 1
                visit(occupied | bit, outgoing | bit, vertical)
                if y + 1 < self.width and not occupied & (bit << 1):
                    visit(occupied | bit | (bit << 1), outgoing,
                          vertical | bit)

            # Without the periodic edge (width-1, 0), then with that edge.
            visit(incoming, 0, 0)
            wrap_sites = 1 | (1 << (self.width - 1))
            if not incoming & wrap_sites:
                visit(incoming | wrap_sites, 0, 1 << (self.width - 1))

        shape = (len(self.states), len(self.states))
        matrix = coo_matrix((values, (rows, columns)), shape=shape).tocsr()
        vertical = coo_matrix((vertical_values, (rows, columns)),
                              shape=shape).tocsr()
        vertical.eliminate_zeros()
        return matrix, vertical

    @property
    def matrix(self) -> csr_matrix:
        return self._operators[0]

    @property
    def vertical_insertion(self) -> csr_matrix:
        """Insert occupation of vertical edge (y=0, y=1) in a column."""
        return self._operators[1]

    @cached_property
    def fixed_point(self) -> tuple[float, np.ndarray]:
        # A fixed starting vector makes numerical benchmark runs reproducible.
        values, vectors = eigsh(self.matrix, k=1, which="LA", tol=1e-13,
                                v0=np.ones(len(self.states)))
        vector = vectors[:, 0]
        if vector.sum() < 0:
            vector = -vector
        return float(values[0]), vector / np.linalg.norm(vector)

    def _observable(self, orientation: str) -> tuple[float, np.ndarray]:
        lam, vector = self.fixed_point
        if orientation == "horizontal":
            applied = (self.states & 1) * vector
        elif orientation == "vertical":
            applied = self.vertical_insertion @ vector / lam
        else:
            raise ValueError("orientation must be 'horizontal' or 'vertical'")
        mean = float(vector @ applied)
        return mean, applied - mean * vector

    def correlations(self, max_distance: int,
                     orientation: str = "horizontal") -> np.ndarray:
        """Connected <n(0)n(r)>-<n>² along the infinite cylinder axis.

        Horizontal means a bond parallel to the cylinder axis. Vertical means
        a bond around its circumference. Both bonds have transverse position
        y=0 and the same orientation; distances include r=0.
        """
        if max_distance < 0:
            raise ValueError("max_distance must be nonnegative")
        mean, centered = self._observable(orientation)
        lam, fixed = self.fixed_point
        result = np.empty(max_distance + 1)
        result[0] = mean * (1 - mean)
        propagated = centered.copy()
        first_power = 1 if orientation == "horizontal" else 0
        for distance in range(1, max_distance + 1):
            if distance > 1 or first_power:
                propagated = self.matrix @ propagated / lam
                # Remove the fixed point to avoid a roundoff-induced plateau.
                propagated -= fixed * (fixed @ propagated)
            result[distance] = centered @ propagated
        return result

    def transfer_length(self, orientation: str = "horizontal") -> dict:
        """Longest transfer length with nonzero overlap with the observable.

        Eigenvalues in different flux blocks cannot contribute to these dimer
        correlations. A gap between arbitrary blocks would be a different
        quantity. The returned length is finite on a finite circumference.
        """
        lam, _ = self.fixed_point
        _, centered = self._observable(orientation)
        size = self.matrix.shape[0]
        if size <= 64:
            values, vectors = np.linalg.eigh(self.matrix.toarray())
        else:
            # Degenerate levels may have a vector with zero overlap; retain a
            # generous group of largest-magnitude eigenvectors on both signs.
            values, vectors = eigsh(self.matrix, k=min(32, size - 2),
                                    which="LM", tol=1e-12,
                                    v0=np.arange(1, size + 1, dtype=float))
        weights = np.square(vectors.T @ centered)
        candidates = [i for i in range(len(values))
                      if weights[i] > 1e-12 * max(float(centered @ centered), 1e-30)
                      and abs(values[i]) < lam * (1 - 1e-10)]
        if not candidates:
            raise RuntimeError("No coupled subleading transfer eigenvalue found")
        index = max(candidates, key=lambda i: abs(values[i]))
        ratio = abs(float(values[index])) / lam
        return {"correlation_length": float(-1 / np.log(ratio)),
                "dominant_eigenvalue": lam,
                "coupled_eigenvalue": float(values[index]),
                "eigenvalue_ratio_abs": ratio}


def square_cylinder_benchmark(width: int = 8, max_distance: int = 64,
                              orientation: str = "horizontal") -> dict:
    """A JSON-serializable, exactly contracted orthogonal-dimer benchmark."""
    cylinder = SquareDimerCylinder(width)
    mean, _ = cylinder._observable(orientation)
    return {
        "model": "orthogonal nearest-neighbour dimers",
        "lattice": "square", "geometry": "infinite cylinder",
        "width": width, "orientation": orientation,
        "observable": "connected dimer occupation",
        "method": "exact sparse column transfer matrix, zero staggered flux",
        "transfer_dimension": len(cylinder.states),
        "occupation": mean,
        "distance": list(range(max_distance + 1)),
        "connected_correlation": cylinder.correlations(max_distance, orientation).tolist(),
        **cylinder.transfer_length(orientation),
        "interpretation": (
            "Finite circumference gives exponential long-distance decay. "
            "The infinite square lattice has algebraic dimer correlations; "
            "take a circumference scaling limit before claiming a plane exponent."
        ),
    }


def square_plane_benchmark(max_distance: int = 128,
                           orientation: str = "horizontal") -> dict:
    """Infinite-plane axial dimer correlations by Kasteleyn quadrature.

    Choose the bipartite Kasteleyn matrix with horizontal entries 1 and
    vertical entries i. Its Fourier symbol is 2*cos(kx)+2j*cos(ky). Integrating
    ky gives, for odd n,

      g(n,0) = integral_0^(pi/2) cos(n*k)/sqrt(1+cos(k)**2) dk / pi.

    For horizontal (collinear) bonds, the two-edge determinant gives
    -g(r-1)*g(r+1) at even r and -g(r)**2 at odd r. For vertical bonds
    separated along the horizontal axis, it gives -(g(r-1)+g(r+1))**2/4
    at even r and +g(r)**2 at odd r. The latter follows also from the
    inverse-matrix relation g(r-1)+g(r+1)+2j*g(r,1)=0.

    Boundary integration yields g(n) ~ sin(n*pi/2)/(pi*n), hence horizontal
    C(r) ~ (-1)**r/(pi*r)**2. Vertical C(r) has a positive r**-2 odd branch
    and a negative r**-4 even branch due to a cancellation of leading terms.
    These are numerical evaluations of an exact classical model formula,
    not a PEPS simulation or synthetic points from an assumed power law.
    """
    if max_distance < 1:
        raise ValueError("max_distance must be positive")
    if orientation not in ("horizontal", "vertical"):
        raise ValueError("orientation must be 'horizontal' or 'vertical'")
    green, errors = {}, []
    for n in range(1, max_distance + 2, 2):
        value, error = quad(lambda k: 1 / np.sqrt(1 + np.cos(k) ** 2),
                            0, np.pi / 2, weight="cos", wvar=n,
                            epsabs=1e-13, epsrel=1e-13)
        green[n] = value / np.pi
        errors.append(error / np.pi)
    correlations = [3 / 16]
    for distance in range(1, max_distance + 1):
        if orientation == "vertical" and distance % 2:
            correlations.append(green[distance] ** 2)
        elif orientation == "vertical":
            correlations.append(-(green[distance - 1] + green[distance + 1]) ** 2 / 4)
        elif distance % 2:
            correlations.append(-green[distance] ** 2)
        else:
            correlations.append(-green[distance - 1] * green[distance + 1])
    result = {
        "model": "orthogonal nearest-neighbour dimers",
        "lattice": "square", "geometry": "infinite plane",
        "orientation": orientation, "displacement_direction": "horizontal",
        "occupation": 0.25,
        "observable": "connected dimer occupation",
        "method": "exact Kasteleyn inverse evaluated by one-dimensional quadrature",
        "distance": list(range(max_distance + 1)),
        "connected_correlation": correlations,
        "largest_green_function_quadrature_error": max(errors),
        "asymptotic_exponent": 2.0,
        "asymptotic_exponent_even": 2.0 if orientation == "horizontal" else 4.0,
        "asymptotic_exponent_odd": 2.0,
        "asymptotic_formula": (
            "C(r) ~ (-1)^r / (pi^2 r^2)" if orientation == "horizontal" else
            "C(r odd) ~ +1/(pi^2 r^2); C(r even) ~ -1/(pi^2 r^4)"
        ),
        "interpretation": (
            "The infinite-plane uniform orthogonal dimer model has algebraic "
            "dimer correlations and no finite exponential correlation length. "
            "This exponent is not the exponent of nonorthogonal singlet RVB dimers."
        ),
    }
    if max_distance >= 16:
        for parity, label in ((0, "even"), (1, "odd")):
            fit_distances = np.arange(max(8, max_distance // 4), max_distance + 1)
            fit_distances = fit_distances[fit_distances % 2 == parity]
            exponent = -np.polyfit(np.log(fit_distances),
                                   np.log(np.abs(np.array(correlations)[fit_distances])), 1)[0]
            result[f"fitted_{label}_distance_exponent"] = float(exponent)
            result[f"fit_{label}_distance_range"] = [int(fit_distances[0]), int(fit_distances[-1])]
        # Keep the original key for callers of the horizontal benchmark.
        result["fit_distance_range"] = result["fit_even_distance_range"]
    return result


class KagomeOrthogonalDimers:
    """Exact uniform kagome dimers from the honeycomb arrow representation.

    Kagome sites are edges of a honeycomb graph and triangles are its vertices.
    An arrow points to the triangle containing that site's dimer. A triangle
    has zero or two incoming arrows. This is an affine system over GF(2); each
    solution corresponds to exactly one perfect matching. No sampling is used.

    Site indexing is 3*(x*ny+y)+s with s=(a,b,c). Up triangles are (a[x,y],
    b[x,y],c[x,y]); down triangles are (a[x,y],b[x-1,y],c[x,y-1]).
    The torus must have an even number of unit cells. Each period is >=2 to
    avoid self-identification of edges in this simple graph convention.
    """

    def __init__(self, nx: int = 8, ny: int = 8):
        if nx < 2 or ny < 2 or nx * ny % 2:
            raise ValueError("nx, ny must be >=2 and nx*ny must be even")
        self.nx, self.ny = nx, ny
        self.nsites = 3 * nx * ny
        self.triangles: list[tuple[int, int, int]] = []
        self.edge_constraints: dict[tuple[int, int], tuple[tuple[int, int], ...]] = {}
        equations = []
        for x in range(nx):
            for y in range(ny):
                up = (self.site(x, y, 0), self.site(x, y, 1), self.site(x, y, 2))
                down = (self.site(x, y, 0), self.site(x - 1, y, 1),
                        self.site(x, y - 1, 2))
                for triangle, arrow_value, parity in ((up, 0, 1), (down, 1, 0)):
                    self.triangles.append(triangle)
                    lhs = sum(1 << site for site in triangle)
                    equations.append((lhs, parity))
                    for i, j in ((0, 1), (1, 2), (2, 0)):
                        edge = tuple(sorted((triangle[i], triangle[j])))
                        self.edge_constraints[edge] = (
                            (edge[0], arrow_value), (edge[1], arrow_value))
        self._solve_parities(equations)

    def site(self, x: int, y: int, sublattice: int) -> int:
        return 3 * ((x % self.nx) * self.ny + y % self.ny) + sublattice

    @property
    def edges(self) -> tuple[tuple[int, int], ...]:
        return tuple(sorted(self.edge_constraints))

    def _solve_parities(self, equations: Iterable[tuple[int, int]]) -> None:
        basis: dict[int, tuple[int, int]] = {}
        for lhs, rhs in equations:
            while lhs:
                pivot = lhs.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = (lhs, rhs)
                    break
                other, value = basis[pivot]
                lhs, rhs = lhs ^ other, rhs ^ value
            if not lhs and rhs:
                raise ValueError("No perfect matching satisfies the triangle parities")
        self._masks: list[int] = []
        self._constants: list[int] = []
        self.free_arrows = 0
        for site in range(self.nsites):
            if site not in basis:
                mask, constant = 1 << self.free_arrows, 0
                self.free_arrows += 1
            else:
                lhs, constant = basis[site]
                lhs ^= 1 << site
                mask = 0
                while lhs:
                    bit = lhs & -lhs
                    other = bit.bit_length() - 1
                    mask ^= self._masks[other]
                    constant ^= self._constants[other]
                    lhs ^= bit
            self._masks.append(mask)
            self._constants.append(constant)

    def _probability(self, constraints: Iterable[tuple[int, int]]) -> float:
        basis: dict[int, tuple[int, int]] = {}
        for site, value in constraints:
            lhs, rhs = self._masks[site], value ^ self._constants[site]
            while lhs:
                pivot = lhs.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = (lhs, rhs)
                    break
                other, other_rhs = basis[pivot]
                lhs, rhs = lhs ^ other, rhs ^ other_rhs
            if not lhs and rhs:
                return 0.0
        return 2.0 ** -len(basis)

    def occupation(self, edge: tuple[int, int]) -> float:
        return self._probability(self.edge_constraints[tuple(sorted(edge))])

    def joint_occupation(self, first: tuple[int, int], second: tuple[int, int]) -> float:
        constraints = (self.edge_constraints[tuple(sorted(first))]
                       + self.edge_constraints[tuple(sorted(second))])
        return self._probability(constraints)

    def connected(self, first: tuple[int, int], second: tuple[int, int]) -> float:
        return self.joint_occupation(first, second) - self.occupation(first) * self.occupation(second)


def kagome_orthogonal_benchmark(nx: int = 8, ny: int = 8) -> dict:
    """Exact occupation probabilities and an axial correlation family."""
    model = KagomeOrthogonalDimers(nx, ny)
    reference = (model.site(0, 0, 0), model.site(0, 0, 1))
    distances = list(range(nx // 2 + 1))
    targets = [(model.site(x, 0, 0), model.site(x, 0, 1)) for x in distances]
    nonzero = [{"edge": list(edge), "connected_correlation": value}
               for edge in model.edges
               if (value := model.connected(reference, edge)) != 0]
    return {
        "model": "orthogonal nearest-neighbour dimers",
        "lattice": "kagome", "geometry": "periodic torus",
        "unit_cells": [nx, ny], "sites": model.nsites,
        "method": "exact honeycomb-arrow parity constraints over GF(2)",
        "covering_count_log2": model.free_arrows,
        "occupation": model.occupation(reference),
        "observable": "connected dimer occupation",
        "reference_edge": list(reference),
        "target_edges": [list(edge) for edge in targets],
        "distance_unit": "unit-cell translations along a1",
        "distance": distances,
        "connected_correlation": [model.connected(reference, edge) for edge in targets],
        "all_nonzero_correlations_with_reference": nonzero,
        "interpretation": (
            "Uniform orthogonal kagome dimers have strictly finite-range dimer "
            "correlations. Exact zeros are not an exponential fit and do not "
            "imply the same correlations as the nonorthogonal singlet RVB state."
        ),
    }
