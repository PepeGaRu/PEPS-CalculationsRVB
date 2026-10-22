# Conventional singlet RVB on the decorated square-octagon lattice

This calculation uses equal coefficients of all fully covering singlet configurations, with the specified triangle and dummy-link singlet arrows. It has eight physical spin-1/2 degrees of freedom per square cell. An unused dummy link contains an internal singlet; there are no vacancies. No semionic loop sign is included.

**Result:** the measured connected spin-bond correlations decay rapidly. Across the late resolved fitting windows, a useful descriptive scale is **about 0.8 square unit cells**. This is not a precision thermodynamic-limit estimate. The effective length still drifts with distance and truncation; a power prefactor multiplying an exponential is allowed.

## Environment and truncation checks

The CTMRG convergence flag refers to three consecutive changes in all corner spectra below the specified tolerance. It does not establish convergence in environment dimension or truncation. The implementation has both a relative projector cutoff and a separate absolute cutoff inside its multiplet-preserving SVD. Equal results at larger nominal chi can conceal equal retained ranks.

| Run | chi | Relative cutoff | Multiplet absolute cutoff | Final projector ranks | Iterations | Corner residual / tolerance | Converged |
|---|---:|---:|---:|---|---:|---|---|
| [chi16_cutoff1e-08.json](chi16_cutoff1e-08.json) | 16 | 1e-08 | 1e-14 | [7, 7, 7, 7] | 20 | 8.29e-12 / 1e-10 | True |
| [chi16_cutoff1e-12.json](chi16_cutoff1e-12.json) | 16 | 1e-12 | 1e-14 | [14, 14, 14, 14] | 21 | 2.52e-11 / 1e-10 | True |
| [chi32_cutoff1e-08.json](chi32_cutoff1e-08.json) | 32 | 1e-08 | 1e-14 | [7, 7, 7, 7] | 20 | 8.29e-12 / 1e-10 | True |
| [chi32_cutoff1e-12.json](chi32_cutoff1e-12.json) | 32 | 1e-12 | 1e-14 | [30, 30, 30, 30] | 1000 | 5.87e-10 / 1e-10 | False |
| [chi32_cutoff1e-12_raw.json](chi32_cutoff1e-12_raw.json) | 32 | 1e-12 | 1e-14 | [32, 32, 32, 32] | 100 | 1.60e-07 / 1e-10 | False |
| [chi32_cutoff1e-12_scale16.json](chi32_cutoff1e-12_scale16.json) | 32 | 1e-12 | 1e-14 | [30, 30, 30, 30] | 18 | 8.05e-10 / 1e-08 | True |
| [chi64_cutoff1e-12.json](chi64_cutoff1e-12.json) | 64 | 1e-12 | 1e-14 | [30, 30, 30, 30] | 18 | 5.37e-10 / 1e-08 | True |
| [chi64_cutoff1e-12_mult1e-18.json](chi64_cutoff1e-12_mult1e-18.json) | 64 | 1e-12 | 1e-18 | [34, 34, 34, 34] | 18 | 8.75e-10 / 1e-08 | True |
| [chi64_cutoff1e-14.json](chi64_cutoff1e-14.json) | 64 | 1e-14 | 1e-14 | [30, 30, 30, 30] | 18 | 4.27e-09 / 1e-08 | True |
| [chi64_cutoff1e-14_mult1e-18.json](chi64_cutoff1e-14_mult1e-18.json) | 64 | 1e-14 | 1e-18 | [35, 35, 35, 35] | 22 | 5.03e-09 / 1e-08 | True |

The raw-physical run retains all 256 local physical states and uses the exact integer tensor divided by 16. It is an independent numerical precision check of the exact rank-41 physical compression, not an independent physical model. Its stricter iteration criterion was not reached; this audit did not cure small-singular-value conditioning.

## Measured observables and fitting windows

Triangle E contains `(r0,r3,dE)` = physical indices `(0,3,4)`. Ring-ring and ring-dummy observables are `S_i dot S_j`. The horizontal dummy bond is `S_dE(cell) dot S_dW(cell+x)`; the vertical analogue uses `dN,dS`. Bond-center distance is measured in square Bravais cells. Two-cell dummy-bond insertions use separation r>=2 to avoid overlap.

All connected covariances use two one-point functions from the same finite CTM strip. Centered operator insertions are independently compared with direct subtraction. The table uses the retained data with absolute correlation greater than 1e-13. Fits describe positive bond correlations; oscillating single-spin data use separate even/odd branches. Spin data below the floor must not be interpreted as long-distance spin correlations.

Selected contraction: [chi64_cutoff1e-14_mult1e-18.json](chi64_cutoff1e-14_mult1e-18.json).

| Observable | Requested window | Used distances | Exponential xi | Exponential log RMS | Power-law exponent | Power-law log RMS | Power × exponential p, xi |
|---|---|---|---:|---:|---:|---:|---|
| triangle_ring_ring_x | [2.0, 5.0] | [2.0, 3.0, 4.0, 5.0] | 0.6901 | 0.0336 | 4.7 | 0.184 | not identified |
| triangle_ring_ring_x | [4.0, 8.0] | [4.0, 5.0, 6.0, 7.0, 8.0] | 0.7496 | 0.026 | 7.67 | 0.166 | not identified |
| triangle_ring_ring_x | [8.0, 14.0] | [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0] | 0.8001 | 0.0155 | 13.4 | 0.186 | 1.04, 0.867 |
| triangle_ring_ring_x | [12.0, 18.0] | [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0] | 0.8163 | 0.00703 | 18.1 | 0.136 | 0.891, 0.858 |
| triangle_ring_ring_x | [16.0, 20.0] | [16.0, 17.0, 18.0, 19.0, 20.0] | 0.8228 | 0.00187 | 21.8 | 0.0548 | not identified |
| triangle_ring_dummy_x | [2.0, 5.0] | [2.0, 3.0, 4.0, 5.0] | 0.6737 | 0.0386 | 4.82 | 0.184 | not identified |
| triangle_ring_dummy_x | [4.0, 8.0] | [4.0, 5.0, 6.0, 7.0, 8.0] | 0.7405 | 0.0305 | 7.77 | 0.164 | not identified |
| triangle_ring_dummy_x | [8.0, 14.0] | [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0] | 0.7976 | 0.0172 | 13.4 | 0.185 | 1.15, 0.871 |
| triangle_ring_dummy_x | [12.0, 18.0] | [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0] | 0.8153 | 0.00762 | 18.1 | 0.136 | 0.966, 0.861 |
| triangle_ring_dummy_x | [16.0, 20.0] | [16.0, 17.0, 18.0, 19.0] | 0.8214 | 0.00132 | 21.2 | 0.0335 | not identified |
| triangle_other_ring_dummy_x | [2.0, 5.0] | [2.0, 3.0, 4.0, 5.0] | 0.6737 | 0.0386 | 4.82 | 0.184 | not identified |
| triangle_other_ring_dummy_x | [4.0, 8.0] | [4.0, 5.0, 6.0, 7.0, 8.0] | 0.7405 | 0.0305 | 7.77 | 0.164 | not identified |
| triangle_other_ring_dummy_x | [8.0, 14.0] | [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0] | 0.7976 | 0.0172 | 13.4 | 0.185 | 1.15, 0.871 |
| triangle_other_ring_dummy_x | [12.0, 18.0] | [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0] | 0.8153 | 0.00762 | 18.1 | 0.136 | 0.966, 0.861 |
| triangle_other_ring_dummy_x | [16.0, 20.0] | [16.0, 17.0, 18.0, 19.0] | 0.8214 | 0.00132 | 21.2 | 0.0335 | not identified |
| dummy_bond_x | [2.0, 5.0] | [2.0, 3.0, 4.0, 5.0] | 0.6856 | 0.0343 | 4.73 | 0.185 | not identified |
| dummy_bond_x | [4.0, 8.0] | [4.0, 5.0, 6.0, 7.0, 8.0] | 0.7471 | 0.0275 | 7.7 | 0.165 | not identified |
| dummy_bond_x | [8.0, 14.0] | [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0] | 0.7998 | 0.016 | 13.4 | 0.185 | 1.07, 0.868 |
| dummy_bond_x | [12.0, 18.0] | [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0] | 0.8164 | 0.00711 | 18.1 | 0.136 | 0.901, 0.859 |
| dummy_bond_x | [16.0, 20.0] | [16.0, 17.0, 18.0, 19.0, 20.0] | 0.8229 | 0.00187 | 21.8 | 0.0548 | not identified |

## Cutoff sensitivity at fixed distance

Comparison: selected relative cutoff 1e-14 versus [chi64_cutoff1e-12_mult1e-18.json](chi64_cutoff1e-12_mult1e-18.json) at 1e-12, both chi=64 and absolute multiplet cutoff 1e-18. These differences are diagnostic, not rigorous error bars.

| Observable | C(14) | Relative change at 14 | C(18) | Relative change at 18 |
|---|---:|---:|---:|---:|
| triangle_ring_ring_x | 3.05461e-10 | 3.46% | 2.31175e-12 | 8.09% |
| triangle_ring_dummy_x | 1.25908e-10 | 2.92% | 9.48077e-13 | 7.28% |
| triangle_other_ring_dummy_x | 1.25908e-10 | 2.92% | 9.48077e-13 | 7.28% |
| dummy_bond_x | 5.97317e-10 | 2.97% | 4.52289e-12 | 7.36% |

The stronger cutoff changes late amplitudes enough that quoting many digits in a correlation length would be misleading. A short correlation length is supported on the resolved range; neither a mathematically exact exponential asymptote nor a universal exponent is established by these finite-environment contractions.

## Reproduction

```bash
/workspace/peps-environment/.venv/bin/python scripts/run_octagon.py --chi 64 --distance 24 --projector-tol 1e-14 --multiplet-abstol 1e-18 --ctm-tol 1e-8 --max-iter 120
/workspace/peps-environment/.venv/bin/python scripts/summarize_octagon_ctmrg.py
```

The pinned peps-torch checkout and CPU torch environment are required. Factorized two-cell correlations are validated against dense small-physical-dimension operators in all four directions and against a polarized product state. Source hashes, package versions, raw correlations, convergence histories, final corner spectra, and all algorithm settings are in each JSON file.
