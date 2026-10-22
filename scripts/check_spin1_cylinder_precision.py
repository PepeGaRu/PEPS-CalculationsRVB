#!/usr/bin/env python3
"""Audit the weak dimer-projector cylinder tail in an eight-dimensional SU(2) block.

The full transfer is constructed by the production code. Its spin-zero block
has a basis of two vacuum/trace states and three delta pairings in each of
two four-vector color spaces. Rational reconstruction is checked against
the floating-point construction; mpmath then evaluates the small spectral
weights without subtracting order-one expectation values.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.octahedral_spin1_transfer import Spin1OctahedralCylinder
from rvb.octahedral_spin1 import levi_civita


def scalar_basis(indices):
    """Unnormalized invariant tensors in the dominant ownership sector."""
    colors = np.array(np.unravel_index(indices, (8,)*4)).T
    basis = np.zeros((len(indices), 8))
    basis[np.all(colors == 0, axis=1), 0] = 1
    basis[np.all(colors == 1, axis=1), 1] = 1
    pairings = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
    for offset, first_column in ((2, 2), (5, 5)):
        for k, pairs in enumerate(pairings):
            selected = np.all((colors >= offset) & (colors < offset+3), axis=1)
            for a, b in pairs:
                selected &= colors[:, a] == colors[:, b]
            basis[selected, first_column+k] = 1
    assert np.linalg.matrix_rank(basis) == 8
    return basis


def evaluate(matrix, operator):
    values, right = mp.eig(matrix, left=False, right=True)
    order = sorted(range(len(values)), key=lambda k: abs(values[k]), reverse=True)
    values = [values[k] for k in order]
    right = mp.matrix([[right[i, k] for k in order] for i in range(len(values))])
    left = right**-1
    transformed = left*(operator/values[0])*right
    weights = [transformed[0, k]*transformed[k, 0] for k in range(1, len(values))]
    return values, transformed[0, 0], weights


def vector_basis(indices):
    """J=1, magnetic-component-zero tensors epsilon_abz delta_cd."""
    colors = np.array(np.unravel_index(indices, (8,)*4)).T
    columns = []
    epsilon = levi_civita()
    for offset in (2, 5):
        selected = np.all((colors >= offset) & (colors < offset+3), axis=1)
        local = colors[selected]-offset
        for a, b in combinations(range(4), 2):
            c, d = [k for k in range(4) if k not in (a, b)]
            column = np.zeros(len(indices))
            column[selected] = epsilon[local[:, a], local[:, b], 2]*(local[:, c] == local[:, d])
            columns.append(column)
    basis = np.array(columns).T
    assert np.linalg.matrix_rank(basis) == 12
    return basis


def spin_audit(cylinder, indices, scalar, scalar_left, scalar_matrix, serial):
    """Use polynomial spectral projectors to avoid degenerate eigenvectors."""
    vector = vector_basis(indices)
    vector_left = np.linalg.solve(vector.T@vector, vector.T)
    transfer = cylinder.transfer[indices, :][:, indices].toarray()
    insertion = cylinder.observable_transfer("spin_z")[indices, :][:, indices].toarray()
    raw_transfer = vector_left@transfer@vector
    raw_to_vector = np.real_if_close(-1j*(vector_left@insertion@scalar)).real
    raw_to_scalar = np.real_if_close(-1j*(scalar_left@insertion@vector)).real
    def rational_matrix(raw):
        fractions = [[Fraction(float(value)).limit_denominator(3**12) for value in row] for row in raw]
        numeric = np.array([[float(value) for value in row] for row in fractions])
        assert np.max(abs(raw-numeric)) < 1e-10
        return mp.matrix([[mp.mpf(value.numerator)/value.denominator for value in row] for row in fractions])
    matrix = rational_matrix(raw_transfer)
    to_vector, to_scalar = rational_matrix(raw_to_vector), rational_matrix(raw_to_scalar)
    invariant_residual = float(np.linalg.norm(transfer@vector-vector@raw_transfer))
    assert invariant_residual < 1e-10
    spectrum = np.linalg.eigvals(raw_transfer)
    # These spin-sector eigenvalues are rational. Check the polynomial
    # eigenspace projectors rather than trusting one vector in a degeneracy.
    assert max(abs(spectrum.imag)) < 1e-10
    distinct = sorted({Fraction(float(value.real)).limit_denominator(3**12) for value in spectrum},
                      key=lambda value: -abs(value))
    eigenvalues = [mp.mpf(value.numerator)/value.denominator for value in distinct]
    scalar_values, scalar_right = mp.eig(scalar_matrix, left=False, right=True)
    k0 = max(range(len(scalar_values)), key=lambda k: abs(scalar_values[k]))
    lambda0 = scalar_values[k0]
    scalar_left_mp = scalar_right**-1
    r0, l0 = scalar_right[:, k0], scalar_left_mp[k0, :]
    identity = mp.eye(12)
    projector_sum = mp.zeros(12)
    modes = []
    for j, value in enumerate(eigenvalues):
        projector = mp.eye(12)
        for k, other in enumerate(eigenvalues):
            if j != k:
                projector = projector*(matrix-other*identity)/(value-other)
        projector_sum += projector
        # Each Sz insertion is i times its stored real reduced matrix.
        weight = -(l0*to_scalar*projector*to_vector*r0)[0]/lambda0**2
        modes.append({"eigenvalue": serial(value), "weight": serial(weight),
                      "xi": serial(-1/mp.log(abs(value/lambda0)))})
    assert mp.norm(projector_sum-identity) < mp.mpf("1e-50")
    coupled = [mode for mode in modes if abs(mp.mpf(mode["weight"])) > mp.mpf("1e-30")]
    assert len(coupled) == 2
    return {"vector_subspace_dimension": 12,
            "invariant_subspace_absolute_residual": invariant_residual,
            "method": "Epsilon-delta J=1 basis; rational matrices; polynomial projectors onto degenerate transfer eigenspaces",
            "distinct_modes": modes, "coupled_modes": coupled,
            "asymptotic_xi": coupled[0]["xi"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--output", type=Path, default=ROOT/"results/octahedral-spin1-cylinder")
    args = parser.parse_args()
    if args.digits < 40:
        parser.error("Use at least 40 digits for this precision audit")
    mp.mp.dps = args.digits
    args.output.mkdir(parents=True, exist_ok=True)
    cylinder = Spin1OctahedralCylinder()
    assert cylinder.dominant_sector() == 0
    indices = cylinder.sector_indices(0)
    basis = scalar_basis(indices)
    left_inverse = np.linalg.solve(basis.T@basis, basis.T)
    reduced, raw, records = {}, {}, {}
    for name, full in (("transfer", cylinder.transfer),
                       ("projector", cylinder.observable_transfer("dimer_projector"))):
        block = full[indices, :][:, indices].toarray()
        small = left_inverse@block@basis
        invariant_residual = float(np.linalg.norm(block@basis-basis@small))
        fractions = [[Fraction(float(value)).limit_denominator(3**12) for value in row]
                     for row in small]
        rational_values = np.array([[float(value) for value in row] for row in fractions])
        reconstruction_error = float(np.max(abs(small-rational_values)))
        assert invariant_residual < 1e-10
        assert reconstruction_error < 1e-10
        reduced[name] = mp.matrix([[mp.mpf(value.numerator)/value.denominator for value in row]
                                   for row in fractions])
        raw[name] = mp.matrix([[mp.mpf(float(value)) for value in row] for row in small])
        records[name] = {"matrix_exact_reconstruction": [[str(value) for value in row] for row in fractions],
                         "invariant_subspace_absolute_residual": invariant_residual,
                         "rational_reconstruction_max_abs_error": reconstruction_error}
        (args.output/(name+"-scalar-matrix.txt")).write_text(
            "\n".join(" ".join(str(value) for value in row) for row in fractions)+"\n")
    values, mean, weights = evaluate(reduced["transfer"], reduced["projector"])
    raw_values, raw_mean, raw_weights = evaluate(raw["transfer"], raw["projector"])
    xi = -1/mp.log(abs(values[1]/values[0]))
    raw_xi = -1/mp.log(abs(raw_values[1]/raw_values[0]))
    # The first subleading scalar mode really couples: this is far above
    # both arithmetic residuals and the physical symmetry-zero weights.
    assert abs(weights[0]) > mp.mpf("1e-10")
    assert abs(raw_weights[0]-weights[0]) < mp.mpf("1e-15")
    assert abs(raw_xi-xi) < mp.mpf("1e-10")
    distances = list(range(1, 65))
    correlations = [sum(weights[k-1]*(values[k]/values[0])**(r-1)
                        for k in range(1, len(values))) for r in distances]
    def serial(value):
        assert abs(mp.im(value)) < mp.mpf("1e-30")
        return mp.nstr(mp.re(value), args.digits-10)
    fits = []
    for lo, hi in ((8, 20), (12, 24), (20, 40)):
        xs = list(range(lo, hi+1))
        ys = [mp.log(abs(correlations[r-1])) for r in xs]
        xbar, ybar = mp.mpf(sum(xs))/len(xs), sum(ys)/len(xs)
        slope = sum((x-xbar)*(y-ybar) for x, y in zip(xs, ys))/sum((x-xbar)**2 for x in xs)
        fits.append({"window": [lo, hi], "fitted_envelope_xi": serial(-1/slope)})
    result = {
        "geometry": "Infinite x; periodic 2x2 transverse octahedron centres",
        "observable": "Physical spin-1 singlet projector on the y-diameter pair at z=0",
        "ownership_sector": 0, "sector_dimension": len(indices), "scalar_dimension": 8,
        "arithmetic_digits": args.digits,
        "precision_scope": "80-digit-style evaluation of rationally reconstructed reduced matrices; input construction is also checked without reconstruction. This is not a statement of bulk physical accuracy.",
        "reduction": records,
        "lambda0": serial(values[0]), "mean": serial(mean),
        "eigenvalues": [serial(value) for value in values],
        "connected_mode_weights": [serial(value) for value in weights],
        "asymptotic_xi": serial(xi),
        "raw_float_input_asymptotic_xi": serial(raw_xi),
        "raw_vs_reconstructed_xi_difference": serial(abs(raw_xi-xi)),
        "raw_vs_reconstructed_tail_weight_difference": serial(abs(raw_weights[0]-weights[0])),
        "distance": distances,
        "connected_correlations": [serial(value) for value in correlations],
        "finite_window_envelope_fits": fits,
        "spin_z_independent_audit": spin_audit(cylinder, indices, basis, left_inverse,
                                                reduced["transfer"], serial),
        "source_sha256": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in ("rvb/octahedral_spin1_transfer.py", "scripts/check_spin1_cylinder_precision.py")},
    }
    path = args.output/"high-precision.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print("Asymptotic dimer-projector xi:", serial(xi))
    print("Weak tail weight:", serial(weights[0]))
    print("Raw-input/reconstructed weight difference:", serial(abs(raw_weights[0]-weights[0])))
    print("Wrote", path)


if __name__ == "__main__":
    main()
