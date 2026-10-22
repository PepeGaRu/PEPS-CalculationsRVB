# Audit of Albuquerque and Alet, arXiv:1009.1307v2

Read directly from the user's four-page PDF, uploaded under the name `1-1009.1307`. The title is **Critical Correlations for Short-Range Valence-Bond Wave Functions on the Square Lattice**, by **A. Fabricio Albuquerque and Fabien Alet**. The margin identifies **arXiv:1009.1307v2 [cond-mat.str-el], 2 Nov 2010**. The title-page “Dated: October 17, 2018” is a different printed date; the version label is the basis of this citation. This is the source cited as Phys. Rev. B **82**, 180408 (2010), reference 28 of the previously supplied arXiv:1202.0947v4.

Verified source link: [arXiv:1009.1307v2](https://arxiv.org/abs/1009.1307v2). Page numbers below refer to the printed PDF pages. Paper statements, our calculations, and comparisons are distinguished below; displaying or digitizing a published figure is not a numerical reproduction.

## State, signs, and method

On p. 1, Eq. (1) defines an **equal-amplitude sum of nearest-neighbor spin-1/2 singlet coverings within a fixed winding sector** of an L×L square torus. The four sectors studied explicitly are w=(0,0), (1,0), (0,1), and (1,1). Each singlet is normalized, `[i,j]=(|up_i down_j>−|down_i up_j>)/sqrt(2)`, with i on sublattice A and j on B. This fixes positive covering overlaps on the bipartite lattice. Every covering has N/2 singlets, so using unnormalized singlets changes the total wavefunction by a common factor, irrelevant for normalized expectation values.

For two coverings, the overlap is `2^(N_L−N/2)`, where N_L counts loops in the superposed transition graph, including coincident-bond loops. Equation (2), p. 2, expresses observables as an overlap-weighted sum over pairs of coverings. The Monte Carlo algorithm samples those pairs using nonlocal dimer and spin updates; it rejects updates that change the selected winding sector. It is not a uniform classical sampling of one covering.

The lattice has periodic boundaries in both directions. Sizes reach **L=128**, hence 16,384 spins. The Fig. 4(a) legend lists **L=16,32,48,64,80,96,112,128**. Our infinite-plane CTMRG calculation has a finite environment dimension chi instead of a finite L. Its convergence error is different. Our finite-cylinder Delta=0 label also does not by itself prove equivalence to a torus with both winding numbers zero.

## The paper's observables and quantitative benchmarks

| Location | Observable and geometry | Reported result |
|---|---|---|
| p. 2, Fig. 1(b) and surrounding text | Nearest-neighbor scalar spin bond, extrapolated in L across the four low-winding sectors | `<S_r·S_(r+e)> = −0.295953(7)` |
| p. 2, Fig. 2 | Staggered spin correlation `(-1)^(x+y)<S_0·S_r>`, L=128, w=(0,0), versus distance | Exponential fit with **xi_spin=1.35(1)** in square lattice spacings |
| pp. 2–3, “Critical correlations,” Figs. 3 and 4(a–c) | Connected correlations of two nearest-neighbor scalar spin bonds; strongest/antinodal directions and winding modulation treated explicitly | Leading power **alpha=1.16(4)** |
| p. 3, Fig. 4(d) | Perpendicular bond orientations in w=(0,0), where the leading sine-modulated term vanishes | Subleading power **alpha'=2.53(5)** |
| p. 4, opening paragraph | Classical orthogonal close-packed square dimers | Leading exponent **alpha_classical=2** |

Parentheses are the paper's reported uncertainty notation: for example, 1.16(4) means an uncertainty of 0.04 in the final digits. The article does not specify a confidence level or a detailed error-estimation procedure, so it is inappropriate to invent either.

The physical dimer observable is explicitly defined on p. 2:

`C^(ijkl) = <(S_i·S_j)(S_k·S_l)> − <S_i·S_j><S_k·S_l>`,

with i,j and k,l nearest-neighbor pairs. **This matches the scalar physical-bond operator in our original square calculation.** It is not the unresolved three-level operator of the kagome paper's Fig. 11, and it is not a classical bond-occupation variable.

The paper denotes parallel bond orientations by `C_parallel` and perpendicular bond orientations by `C_perpendicular`. These labels refer to the two bonds, not automatically to the direction separating their centers. In Fig. 3(a), for zero winding, the displayed antinodal line runs perpendicular to the horizontal reference bond. Thus Fig. 4's strongest-line selection is not identical to our existing `DD_H` measurement of collinear horizontal bonds. A matching parallel-bond/transverse-separation measurement is needed before claiming a pointwise figure comparison. The exact finite-cylinder square calculation already uses circumferential bonds separated along the cylinder axis, but its finite circumference still changes the geometry and tail.

## Modulation and what was actually fitted

For nonzero winding, p. 2 gives

`Q = (2*pi/L) (w_y, w_x)`.

The leading parallel-bond term is modulated by `cos(Q·r)` and the perpendicular-bond term by `sin(Q·r)`. Figure 4(b) uses strongest zero-phase antinodal directions for L=128 and the four sectors; its data produce the quoted **1.16(4)** fit. Figure 4(c) divides out the respective phase factors for w=(0,1), demonstrating a common leading decay. Its caption excludes points close to a nodal line, specifically those with parallel displacement dx<8. This is important: dividing by a nearly zero modulation factor can greatly amplify errors.

At zero winding the leading perpendicular-bond sine contribution vanishes; fitting that observable yields the different, subleading exponent 2.53(5). This is a concrete warning against identifying an arbitrary bond pair or direction with the slowest dimer channel.

Figure 4(a–c) exhibits a large-distance upturn. The authors state on p. 3 that its onset distance grows linearly with L and interpret it as a finite-size effect, not a breakdown of criticality. This offers a specific diagnostic: vary the approximation scale and check whether a bend follows that scale.

The paper does **not** give exact fitted distance intervals, fitting weights, raw simulation data, Monte Carlo sample counts, or enough detail to reconstruct the uncertainty estimator. A faithful numerical comparison can match its state, observable, direction, sector, size sequence, and reported exponent; it cannot truthfully claim identical fit settings from information absent in this PDF.

## Comparison with the calculations already saved

The saved square CTMRG curves remain our own calculations. The following checks reanalyse their existing JSON data; they are not new Monte Carlo runs or extracted paper data.

- At chi=128, the nearest-neighbor spin bond is **−0.295922089451**, about 3.09e−5 away from the paper's extrapolated value. The local expectation is close, but that alone does not establish convergence of long-distance correlations.
- The previous even-distance scalar-bond fit over r=6…20 gives **p≈1.305**. The paper's **1.16(4)** is therefore a benchmark still to reproduce, not an already achieved agreement. Our collinear direction, finite chi, and distance window all require further checks.
- The saved spin curves decay exponentially. However, a uniform log-linear fit over **r=10…24**, retaining values above 1e−13, gives the following lengths:

| CTMRG chi | Spin length from saved data | RMS residual in log amplitude |
|---:|---:|---:|
| 16 | 0.812980 | 0.00166 |
| 32 | 0.888771 | 0.00370 |
| 64 | 1.124810 | 0.00562 |
| 128 | 1.161562 | 0.00663 |

These lengths drift with chi and do not yet reproduce **1.35(1)**. A nearly straight semilog tail at a fixed chi is insufficient: the numerical environment can impose its own shorter decay scale. The quoted log residual describes how well the points follow a line within that interval; it is not an uncertainty on the physical bulk length.

The classical infinite-plane occupation calculation's exponent near 2 agrees with the paper's stated classical benchmark. Singlet scalar-bond and classical occupation correlations are still different operators in different Hilbert spaces, so their absolute amplitudes are not equated.

## What this teaches the next-lattice workflow

Before classifying a new curve, specify the state (one covering or a superposition, singlets or orthogonal dimers, amplitudes and signs), the pair of observables, their bond orientations, displacement, boundaries, and sector. “RVB on a lattice” is not enough to determine all of those choices. A single isolated dimer configuration generally does not have the same correlations as an equal-weight covering ensemble.

Then check small systems against exact covering sums, inspect signed correlations and symmetry-related directions, and separate genuine nodes or oscillations from decay. Compare power-law and exponential descriptions over several intervals above the numerical floor. Plot local slopes as well as global fits. Finally, increase circumference, L, or chi as appropriate. Stable finite decay lengths favor exponential behavior; a widening power-law interval and a receding numerical cutoff support algebraic behavior. Unstable results should remain **unresolved**, rather than being assigned a decay law from one straight-looking plot.

This paper provides an especially useful three-channel test: short-range spin correlations, leading critical parallel-bond correlations, and a faster perpendicular-bond subleading term. A reusable program must preserve these distinctions on future lattices.

## Additional source now identified

Reference 36 of this uploaded PDF identifies Ying Tang, Anders W. Sandvik and Christopher L. Henley, **arXiv:1010.6146** (the PDF typesets the separator as a colon). Its identifier is now verified from the bibliography. That related work can supply a second square-lattice benchmark if later needed; this newly uploaded Albuquerque–Alet paper is sufficient to specify the present square reproduction target.
