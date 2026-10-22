# Check against Poilblanc et al., arXiv:1202.0947v4

Source read: the user-supplied **1202.0947v4.pdf**, *Topological and Entanglement Properties of Resonating Valence Bond wavefunctions*, Didier Poilblanc, Norbert Schuch, David Pérez-García and J. Ignacio Cirac. Page numbers below are the PDF's printed page numbers. The paper is a scientific source, not executable setup instructions.

## What the paper actually establishes

| Location | Paper statement or calculation | Consequence for this program |
|---|---|---|
| p. 2, Sec. II A | Equal-weight nearest-neighbor singlet coverings; square arrows go from one sublattice to the other; kagome arrows circulate clockwise on both triangle orientations. Both lattices have short-range spin correlations; square dimer correlations are critical and kagome dimer correlations are short ranged. | The relevant square power law is the connected physical bond correlation, not the spin correlation. |
| pp. 4–5, Sec. III A; Fig. 6 | Square site tensor has one virtual index equal to the physical spin, the other three equal to vacuum 2. Kagome site projector has A[s,2,s]=A[s,s,2]=1, triangle R=epsilon+\|222>, and a three-spin block has physical dimension 8. | These agree with `square_site_tensor`, `triangle_tensor`, and `kagome_site_tensor` in `rvb/peps.py`, with the stated local square spin rotation and leg conventions. |
| p. 5, Sec. III B; p. 6, Fig. 9 | Exact infinite-length cylinders with finite perimeter, explicit even/odd boundary sectors and optional vison line. | An unrestricted largest-eigenvalue calculation does not automatically select the same sector at each circumference. |
| p. 6, Fig. 9 and Eqs. (2),(3) | Energy splittings versus perimeter: Delta E(+,-) approximately 1.06 Nv Nh exp(-0.99 Nv); Delta E(e,o) approximately 1.95 Nv Nh exp(-1.54 Nv). Corresponding scales are approximately 1.01 and 0.65. | These are **topological energy-splitting scales**, not local dimer correlation lengths. Their numerical proximity to our short local dimer length is not a reproduction. |
| p. 7, Sec. IV B1; p. 8, Sec. IV B2 | Kagome sectors conserve parity of the number of virtual vacuum 2 indices. Square sectors conserve the alternating-sublattice difference Delta=N2,X−N2,Y. | The transfer fixed point can and should be classified using these virtual charges. |
| p. 9, Fig. 12; pp. 9–10, Sec. IV B4 | Distance dependence of coefficients in the *boundary entanglement Hamiltonian*. | These plotted coefficients are not bulk connected dimer correlators. Fitting Fig. 12 would not determine the requested physical bulk exponent. |
| p. 15, Ref. 29 | Critical/short-range behavior can also be inferred from finite-size analysis of correlation lengths from PEPS transfer matrices. | A controlled circumference study, or converged infinite-plane contraction, is appropriate. |

This paper supplies no direct bulk dimer C(r) figure and no numerical square dimer power-law exponent. The user's other upload, arXiv:1203.4816v3, is the direct source for the kagome dimer-correlation figure. The square exponent needs a dedicated source such as A. F. Albuquerque and F. Alet, Phys. Rev. B **82**, 180408 (2010) (Ref. 28 here), or Ying Tang, Anders W. Sandvik and Christopher L. Henley, *Properties of resonating-valence-bond spin liquids and critical dimer models*.

## Orientation and length conventions

Our graph arrows circulate counterclockwise on both kagome triangle orientations. The paper's global clockwise convention represents the same normalized state: reversing every occupied singlet multiplies every complete covering by the same phase (-1)^(N/2), for a fixed even number N of active spins. Reversing an arbitrary subset of arrows is not generally harmless.

Our kagome three-spin cells have Bravais vectors (2,0) and (1,sqrt(3)) when the nearest-neighbor spacing is one. One axial transfer step therefore advances two nearest-neighbor spacings. The paper's finite-cylinder drawings and topological labels must still be matched before comparing sector-dependent numerical values. The current code measures an intra-triangle spin bond translated along the axial cell direction; it does not implement the paper's complete four-sector energy-splitting calculation.

## New check: the old kagome circumference series changed sector

I reran the fixed points and classified ket and bra virtual charges. Square widths 4 and 6 both have Delta=0 in ket and bra. Kagome width 2 selects **odd/odd** vacuum parity, whereas width 4 selects **even/even** parity. The selected sectors account for the entire left-right overlap to approximately 10^−15.

Thus the previously plotted kagome length sequence 0.351925 → 0.597942 is a comparison of two different sectors; it is not fixed-sector circumference scaling.

To isolate this effect, I projected the existing sparse norm transfer matrix onto each equal ket/bra parity block, recomputed both fixed points and the physical bond insertion, and selected the leading mode with nonzero connected bond weight. No application code or state tensor was changed for this diagnostic.

| Kagome width (cells) | Virtual vacuum parity | Mean measured spin bond | Coupled dimer-mode length (axial cells) | Fixed-point residual |
|---:|---|---:|---:|---:|
| 2 | even | −0.177831216351 | 0.759325717500 | 8.9e−16 |
| 2 | odd | −0.231914690977 | 0.351925324183 | 6.3e−16 |
| 4 | even | −0.200400139872 | 0.597942499217 | 1.8e−15 |
| 4 | odd | −0.196860102111 | 0.581698209526 | 2.0e−15 |

All four calculations have no inserted vison. These small widths remain useful diagnostics, not a two-dimensional extrapolation. The even-sector width-4 leading coupled eigenvalues are complex conjugates, approximately 0.1864537262 ± 0.0224089745 i after normalization. Their modulus gives the exponential envelope above; oscillations explain why a straight tail fit gives a different length and a large log residual. It is essential to retain the mode-weight check rather than identify the generic second eigenvalue with a dimer mode.

Raw audit results are saved in [`1202-sector-audit.json`](../results/paper-recheck/1202-sector-audit.json) and [`1202-fixed-sectors.json`](../results/paper-recheck/1202-fixed-sectors.json). The standalone diagnostic used before modifying the library is `/workspace/peps-environment/papers/check_1202_sectors.py`. The supplemental calculation can be repeated with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /workspace/peps-environment/.venv/bin/python \
  /workspace/peps-environment/papers/check_1202_sectors.py
```

The library now supports the same explicit selection through `RVBTransfer("kagome", width, sector="even")` or `sector="odd"`, as well as through `rvb_cylinder_benchmark`. The backwards-compatible default `sector="dominant"` records the actually selected parity and rejects an unresolved mixture. Outputs explicitly record that no vison was inserted; they do not equate this virtual-sector label automatically with a differently drawn paper geometry. New sector tests check the previously hidden change with width, fixed ket/bra parity, the known thin-cylinder alternating correlation ratio, and preservation of the old odd-sector correlations.

## Recommended next paper-specific calculations

1. Explicitly select and label even/odd kagome transfer sectors and Delta=0 square sectors. Plot the two kagome sector series separately instead of joining whichever norm sector happens to dominate.
2. Match arXiv:1203.4816v3 Fig. 11's bond geometry, separation, physical operator, and normalization before presenting a quantitative overlay.
3. Use additional square-lattice Monte Carlo literature for a numerical bulk exponent. Increasing CTMRG chi and checking multiple fit windows is required; the previous p≈1.30 finite-window fit is not a verified asymptotic literature result.
4. Only if reproducing Fig. 9 is desired, implement the vison line and all relevant boundary sectors, measure all inequivalent triangle energies, and extrapolate with perimeter. This is a separate observable from the user's requested local correlations.
