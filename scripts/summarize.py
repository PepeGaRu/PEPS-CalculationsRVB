"""Build the overview from saved calculations and the primary-paper audit."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    exact = json.loads((ROOT / "results/reference/data.json").read_text())
    square = json.loads((ROOT / "results/upstream/square_chi128.json").read_text())
    exponent = exact["square_orthogonal_plane"]["fitted_even_distance_exponent"]
    p = square["finite_window_dimer_fits"]["DD"]["power_law_exponent"]
    content = f"""# RVB and dimer/trimer covering calculations

The current target is the **three-dimensional bulk spin correlation length**
of the octahedral spin-1 singlet state. The
[3D calculation status](results/octahedral-spin1-bulk/report.md) records
isotropic-box contraction and exact-overlap tests. **No converged bulk
quantum correlation length is established.** The
[geometry diagram](results/octahedral-spin1-bulk/geometry.pdf) explains why
the earlier narrow cylinders do not supply this result.

The [new spin-1 cylinder calculation](results/octahedral-spin1-cylinder/report.md)
and [PDF with decay plots](results/octahedral-spin1-cylinder/Spin1_Quantum_Correlation_Lengths.pdf)
resolve **xi = 2.4901024** cubic-centre spacings for physical singlet-projector
correlations on an infinite cylinder with 2x2 transverse periodic dimensions.
The tail alternates and its envelope decays exponentially. The same transverse
projector on 2x3 has **xi = 0.6381923**. Both transfers retain all states within
their boundary sectors; their large width dependence prevents a bulk 3D
extrapolation. Independent finite-volume contractions and an 80-digit reduced
matrix calculation validate the result.

The [spin-1 quantum report](results/octahedral-spin1/report.md) and
[PDF with quantum checks and plots](results/octahedral-spin1/Octahedral_Spin1_Quantum_Checks.pdf)
extend the octahedral construction to singlet dimers and singlet trimers.
The [state definition](docs/octahedral-spin1-model.md) fixes the phases and
gives an exact virtual-dimension-four tensor network. Quantum contractions
of all 8,192 coverings are validated on a 24-spin periodic cluster.
Its period-two identification of dimer states is a severe finite-size
effect: the bulk quantum correlation length remains undetermined.

The [octahedral dimer/trimer report](results/octahedral-classical/report.md)
and [PDF with correlation plots](results/octahedral-classical/Octahedral_Classical_First_Results.pdf)
study the newly supplied three-dimensional covering notes. With equal weight
for each complete classical covering, the sixteen local states form a binary
parity model. Exact counting and local occupation probabilities are practical
here: the benchmark reaches 1,728 octahedra and 5,184 covering sites.
Connected correlations of single-octahedron observables vanish beyond
neighboring centres, while nonlocal line constraints remain. Quantum dimer
and trimer states are not assumed. Run
`python scripts/run_octahedral_classical.py` to reproduce this calculation.

The [square–octagon report and plots](results/octagon-square/report.md)
study the specified fully covered spin lattice: both unused dummy endpoints
form an internal singlet. See the [model and loop-map derivation](docs/square-octagon-model.md)
for the precise state. Orthogonal-dimer occupation correlations are
**exactly finite range**, vanishing beyond neighboring simulation cells.
Singlet CTMRG calculations support short-range bond decay with a characteristic
scale of about **0.8–0.9 simulation cells (one small square per cell)**; the finite environment and
contraction checks do not yet determine a precise bulk correlation length.
The [new correlation plots](results/octagon-square/correlations.png) and
[decay diagnostics](results/octagon-square/decay-diagnostics.png) display
our calculations. Run `python scripts/summarize_octagon.py` to regenerate
those figures and their source-data summary from saved data.

The [spin-position diagram](results/octagon-square/spin-lattice.png) and
[dummy-link contraction](results/octagon-square/dummy-contraction.png) show the
intended geometry. Contracting dummy links gives a rotated square graph and
a bijection of even-edge patterns. Our earlier 2^N-loss warning applied to
the different operation of collapsing the small squares. The simulation-cell
axes are diagonal relative to the contracted graph; the numeric correlation
values are unchanged. The earlier benchmarks remain below.

The Python program constructs fixed nearest-neighbor spin-1/2 RVB states,
checks them against exact covering sums, and measures connected bond
correlations. The uploaded **1203.4816v3**, **1202.0947v4**, and
**1009.1307v2** have now been audited. Read the [paper recheck](results/paper-recheck/report.md) for the
corrections and the distinction between our calculations and published data.
The supplied **1407.7773v2** is audited separately for the new lattice's
[reference-arrow loop construction](docs/paper-1407-audit.md).

| State | Calculated behavior | Accuracy and interpretation |
|---|---|---|
| Square orthogonal dimers | Power law, fitted exponent {exponent:.5f}; exact asymptotic exponent 2 | Exact plane quadrature; horizontal bonds separated along their direction |
| Square singlet RVB | Intermediate power law, p≈{p:.3f} at χ=128, even r=6–20 | Finite χ and finite fit window; no converged bulk exponent yet |
| Kagome orthogonal dimers | Distant connected occupation correlations exactly zero | Exact arrow constraints; no positive exponential length to fit |
| Kagome singlet RVB | Rapid exponential bond decay; baseline x-tail fit ≈0.638 cells, r=8–14 | Finite χ and projector-cutoff dependent; precise bulk length unresolved |

The square power law concerns **dimer correlations**, while spin correlations
can be short ranged. A power-law correlator has no finite asymptotic
exponential correlation length. One kagome Bravais cell equals two
nearest-neighbor spacings. The labels x/y refer to a specified oriented bond,
not a claim of lattice anisotropy.

![Our square and kagome calculations](results/paper-recheck/square-kagome-comparison.png)

Singlet curves measure the physical bond operator **Sᵢ·Sⱼ**. Orthogonal curves
measure **occupation nᵢⱼ**. Their amplitudes need not agree: declaring
coverings orthogonal changes the state and removes interference between
different coverings. For occupation measurements it reduces the norm problem
to a uniform classical dimer model.

## Corrections after reading the papers

- Kagome cylinder runs keep **even and odd virtual parity sectors separate**.
  Earlier unrestricted runs changed sector between widths 2 and 4.
- Eq. (12) of 1203.4816v3 is implemented with the physical location flag,
  interpolating from singlets at θ=0 to orthogonal coverings at θ=1. Exact
  covering sums validate the tensor and endpoint orthogonality.
- Figure 11 measures a stated **qutrit z⊗z** operator. Its embedding and bond
  convention remain unmatched to our spin-bond or physical ZZ operators.
  Its published points are redrawn from PDF coordinates and labeled separately.
- The ≈0.65 and ≈1.01 scales in Figure 9 of 1202.0947v4 concern **topological
  energy splittings versus circumference**, not the local dimer length.
- Short r=2–8 fits are no longer presented as asymptotic lengths. We hide
  subtraction-level noise and compare multiple windows. Stricter CTMRG
  projector cutoffs expose additional sensitivity; nominal χ agreement alone
  is not sufficient convergence evidence.

![Published Figure 11, extracted plotted points](results/paper-recheck/published-fig11.png)

These are the authors' **published plotted coordinates**, not their raw
simulation data and not a new reproduction. Points below the plotted range
are omitted, not zero. The D=40/D=60 difference at θ=0, r=4 is about 23%.

## Checks and reproducibility

The validation suite checks explicit singlet sums, graph and
blocked PEPS, spin observables, fixed parity sectors, Eq. (12) interpolation,
orthogonal covering Gram matrices, Kasteleyn formulas, and kagome arrow
probabilities, and the reusable decay diagnostics. Numerical tests validate these constructions, not an
infinite-χ extrapolation. The earlier 32-test core also passed in a fresh
environment installed from hashed dependencies; the expanded suite passed
in the current environment.

- [Paper recheck, new calculations, and limitations](results/paper-recheck/report.md)
- [Octahedral classical dimer/trimer counts, correlations and limitations](results/octahedral-classical/report.md)
- [Octahedral spin-1 quantum construction and finite-cluster checks](results/octahedral-spin1/report.md)
- [Spin-1 quantum correlation lengths on infinite cylinders](results/octahedral-spin1-cylinder/report.md)
- [Decorated square–octagon model, correlations and limitations](results/octagon-square/report.md)
- [Square–octagon covering and loop-map derivation](docs/square-octagon-model.md)
- [Exact cylinder and classical results](results/reference/report.md)
- [Independent peps-torch CTMRG results](results/upstream/report.md)
- [Library and literature review](docs/libraries-and-references.md)
- [Plain-language guide](docs/reading-the-results.md)
- [Exponential versus power-law diagnosis](docs/decay-diagnostics.md)
- [New square-paper audit](docs/paper-1009-audit.md)
- [Run commands and setup](README.md)

Run `bash run.sh --profile reference --output results/reference` for exact
benchmarks, then `python scripts/plot_paper_recheck.py` to redraw the saved
comparison figures. `python scripts/summarize.py` refreshes this overview.

The supplied Albuquerque–Alet paper **1009.1307v2**, *Critical Correlations
for Short-Range Valence-Bond Wave Functions on the Square Lattice*, reports
a leading dimer exponent **1.16(4)** and separate spin length **1.35(1)**.
Our finite-χ estimate near 1.30 does not quantitatively reproduce that
exponent. Although the scalar-bond operator agrees, their main fit also uses
a different displacement direction and fixed-winding finite tori. The
[new diagnostics](results/decay-diagnostics/report.md) track local slopes,
window changes and numerical cutoffs without turning a single fit into a
bulk classification.
"""
    (ROOT / "RESULTS.md").write_text(content)
    print(ROOT / "RESULTS.md")


if __name__ == "__main__":
    main()
