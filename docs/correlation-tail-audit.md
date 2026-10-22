# Independent audit of the saved correlation tails

This audit reanalyses the saved `results/upstream/*_chi*.json` outputs. It does not rerun CTMRG and does not turn finite-chi fits into extrapolated bulk results. Complete fit windows, residuals, local slopes, source hashes, and covariance checks are recorded in [`tail_audit.json`](../results/paper-recheck/tail_audit.json).

**Subsequent cutoff check:** new Eq. (12) contractions show that the default
projector cutoff can limit the effective environment rank. Stricter-cutoff
diagnostics change tiny tails and do not all converge. Thus the resolved
baseline x fit near 0.638 cells below remains cutoff dependent, even when
the fit window and nominal chi look stable. See the
[complete paper recheck](../results/paper-recheck/report.md).

## What is actually measured

The adapter measures the connected physical scalar-bond correlation

\[
C_D(r)=\langle (\mathbf S_i\cdot\mathbf S_j)(\mathbf S_k\cdot\mathbf S_l)\rangle-
\langle\mathbf S_i\cdot\mathbf S_j\rangle\langle\mathbf S_k\cdot\mathbf S_l\rangle.
\]

This is not a bond-occupation correlation in a classical/orthogonal dimer ensemble. It is also not the operator used in Fig. 11 of Schuch et al., arXiv:1203.4816v3. Section IV.B.1 on PDF page 7 specifies connected correlations of two `z⊗z` operators, with `z=diag(1,-1,1)`, along the dimer–RVB interpolation. Direct reproduction of that figure requires implementing its operator, its selected bond and direction, and its interpolation. Matching the fixed D=3 RVB tensor and the qualitative exponential decay alone is insufficient.

The scalar and individual spin-component bond correlators need not have identical tails. In an SU(2)-invariant state, the single-bond expectation of `SzSz` equals one third of `S·S`, but this does not imply that a four-spin `SzSz` correlation equals one ninth of the scalar-bond correlation; the former also contains a rank-2 channel.

The square distance convention is correct: the upstream index zero represents two nonoverlapping collinear bonds with centers separated by two lattice spacings. The kagome index zero means neighboring blocked cells, hence distance one. Kagome distances are Bravais-cell translations of length two nearest-neighbor spacings. Thus a fitted length of 0.638 cells corresponds to about 1.28 nearest-neighbor spacings in this direction.

## Kagome: the earlier 0.624 and 0.525 are short-window slopes

The saved chi=32 physical-bond results give:

| Fit window (cells) | xi, x direction | log-fit RMS, x | xi, y direction | log-fit RMS, y |
|---|---:|---:|---:|---:|
| 1–4 | 0.604465 | 0.01882 | 0.452144 | 0.16115 |
| 2–8 | 0.624214 | 0.04451 | 0.525421 | 0.09400 |
| 4–10 | 0.635880 | 0.01061 | 0.547779 | 0.03332 |
| 6–12 | 0.638267 | 0.000962 | 0.559541 | 0.03358 |
| 8–14 | 0.638305 | 0.000477 | 0.573290 | 0.04100 |

Here the reported RMS is the residual in `ln|C|`; it measures lack of a single exponential within a window, not physical uncertainty. The x tail has a convincing finite-chi plateau near 0.638 cells over the reliably resolved distances. Chi=18 gives 0.637302 in the 8–14 window, a 0.16% change. This is better evidence for the resolved x-tail scale than the old 2–8 fit, although larger chi and observable-resolved spectra are still needed for a bulk asymptotic claim.

The y tail has no resolved plateau. The local effective length

\[
\xi_{\mathrm{eff}}(r)=-1/\log |C(r+1)/C(r)|
\]

rises from 0.4916 at 2→3 to 0.5507 at 7→8 and 0.5938 at 13→14. Therefore **0.525 must not be called the asymptotic y correlation length**. Multiple transfer modes or a small amplitude of a slower mode can explain the drift; the available curve does not distinguish them conclusively. There is no justified conclusion of distinct asymptotic x/y lengths.

The generic CTMRG transfer length is 1.02218 cells at chi=32. It is not the dimer correlation length: coupling to the scalar bond was not evaluated. The saved ten eigenvalues stop at magnitude 0.21781, whereas the resolved x tail corresponds to a ratio near `exp(-1/0.6383)=0.20875`. Thus the existing spectrum does not even reach this scalar-tail scale. Multiplicities alone are insufficient to assign physical symmetry sectors.

The transfer ratio in Fig. 12 of arXiv:1203.4816v3 additionally excludes two topological near-degenerate leading eigenvalues and is computed on cylinders, then extrapolated in circumference. It is not directly comparable to the second eigenvalue of the present approximate infinite-plane CTMRG strip.

## Square: finite-chi behavior still limits the exponent

For even separations the saved physical-bond fits give:

| chi | p, 4–12 | p, 6–20 | p, 10–30 | p, 20–40 |
|---:|---:|---:|---:|---:|
| 16 | 1.9262 | 2.6771 | 3.9816 | 5.9583 |
| 32 | 1.4083 | 1.4268 | 1.7184 | 2.5054 |
| 64 | 1.3897 | 1.3751 | 1.4320 | 1.6983 |
| 128 | 1.3539 | 1.3051 | 1.2976 | 1.3102 |

At chi=128, odd separations give 1.2965 for 6–20 and 1.2973 for 10–30, consistent with the staggered even-sector envelope. The apparent power-law interval is real in these data. The 40–80 fit rises to 1.5953, showing the finite-chi cutoff. The curve changes by 8.56% at r=20 and 30.32% at r=40 when chi increases from 64 to 128. Consequently **p≈1.30 remains a descriptive finite-window estimate**, not a reproduced bulk exponent or a contradiction of another numerical estimate. No statistical error bar can replace the remaining systematic chi extrapolation.

The generic leading transfer length increases from 6.21 to 14.00, 23.56, and 49.92 spacings as chi increases. Quoting 49.92 as a physical finite square-RVB correlation length would be misleading. The square spin channel is separately short ranged; its behavior is consistent with the papers and must not be confused with the dimer channel.

## Connected subtraction and numerical floor

Reconstruction of every stored value from `dd - d_left*d_right` agrees exactly at the stored precision. The onepoint factors are contracted with the same strip and normalization as the two-point function, which is the appropriate finite-chi covariance. Replacing them with the separately measured energy-derived mean would introduce artificial constants of approximately `-6.10e-6` at square chi=128 and `-2.24e-9` at kagome chi=32. Those constants would overwhelm genuine long-distance decay. The present subtraction avoids that inconsistency.

The kagome subtraction nevertheless cancels quantities of order 0.04. Float64 cancellation and contraction roundoff eventually dominate values around 1e-16–1e-17. This audit uses `1e-13` as a conservative plotting/fitting floor; it is a numerical guard, not an estimated CTMRG systematic error. In the x 8–14 fit the smallest value is `7.48e-12`; in the y 8–14 fit it is `2.25e-13`. Fits beyond this range should not be used to resolve the drifting y tail from the current subtraction-based data. Centered operator insertions or transfer fixed-point methods would improve tail resolution.

## Required corrections to the presentation

1. Label the 0.624 and 0.525 values as short-window slopes; replace any claim of a definitive kagome correlation length by a qualified statement about the resolved x tail near 0.638 cells and the unresolved y tail.
2. Retain the explicit distinction between physical scalar-bond, component-bond, spin, classical occupation, and generic transfer correlations.
3. Show only above-floor points in logarithmic plots, with the floor and finite-chi limitation stated.
4. To reproduce paper figures numerically, implement the specified observables and sector/geometry definitions before comparing numbers.
5. Keep the square critical conclusion qualitative until a systematic contraction-dimension study establishes a stable exponent.
