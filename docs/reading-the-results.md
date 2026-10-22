# Reading the calculations without programming

The starting state is a sum over ways to pair neighboring sites. There are
two different versions of this idea:

- In the **spin RVB state**, each pair is a spin singlet. Different pairing
  patterns can describe overlapping quantum states. Their interference must
  be included in the norm and every expectation value.
- In the **orthogonal dimer state**, different pairing patterns are declared
  orthogonal. For an equal-weight state and dimer-occupation measurements,
  this is the same calculation as a uniform classical dimer model.

The PEPS representation encodes the spin RVB sum locally. Its virtual bond
dimension is three: two spin labels and an empty label. A local projector
requires exactly one occupied virtual leg per physical spin. Kagome can be
organized using triangle tensors; three physical spins are grouped into a
single tensor for the transfer calculation.

The new calculation also follows the continuous interpolation in Eq. (12)
of your uploaded Schuch paper, `1203.4816v3.pdf`. It attaches an extra label
to each spin. At parameter theta=0 that label carries no information and
the state is the usual RVB. At theta=1 it records where the dimer lies,
making different coverings orthogonal. Values in between are calculated
from the actual tensor equation; they are not an average of the two
endpoint plots.

## What is a connected correlation?

Measure a bond at one location and another far away. Subtract the product
of their separate means:

`C(r) = <B(0) B(r)> - <B(0)> <B(r)>`.

This subtraction removes the part that would be present even if the two
measurements were independent. For spin RVB, `B = S_i · S_j`. For orthogonal
dimers, replace `B` by the occupation `n`, which is one when a bond is
occupied and zero otherwise.

The occupation variable is not the same physical operator as a singlet
projector. That distinction is why the report compares decay behavior
without identifying the two models' absolute amplitudes.

There is another distinction when reading the supplied paper: its Fig. 11
uses an operator written `z ⊗ z` with three-level `z = diag(1,-1,1)`.
Our original graphs measure `S_i · S_j`. We have added explicitly defined
spin-component and three-level-operator checks, but the exact assignment
of the paper's operator to the selected sites and internal labels is still
unresolved. The original curves should therefore be read as independent
calculations of short-range behavior, not a point-by-point reproduction
of Fig. 11.

## Exponential and power-law decay

An exponential envelope has the form `|C(r)| ≈ A exp(-r/xi)`. On a plot with
distance on a linear axis and correlation on a logarithmic axis, its
envelope is a straight line. `xi` is a decay length, in the distance units
used by the calculation. Oscillations can produce dips beneath the
envelope.

A power law has the form `|C(r)| ≈ A/r^p`. It is straight on a plot with
both axes logarithmic. Such a correlator has no finite asymptotic
exponential decay length. Its exponent `p` and an exponential decay length
are different quantities.

Square RVB's critical behavior is in **dimer correlations**. The separate
spin-spin correlation can still be short ranged. The independent-library
figure shows both, so they are not confused.

## Why the program uses two calculations

The custom code wraps the lattice onto a cylinder. Its length is infinite,
but its circumference is finite. Within that geometry the PEPS transfer
is contracted exactly, with numerical eigensolver error reported. This
gives strong, inspectable tests and avoids choosing an approximate
environment dimension. The circumference still limits what can be inferred
about the two-dimensional plane.

The `peps-torch` calculation instead approximates the infinite plane using
CTMRG. Its parameter `chi` controls the size of the surrounding tensor
environment. Increasing `chi` improves that approximation but takes more
time. The program checks both iteration convergence **at a given chi** and
changes **between chi values**. These are separate checks.

For the critical square state, either finite circumference or finite chi
introduces an eventual exponential cutoff. Seeing that cutoff in one run
does not disprove a power law in the exact two-dimensional state. A useful
check is whether the cutoff recedes as the calculation is improved.

For kagome RVB, the physical dimer decay stays short. For orthogonal kagome
dimers, sufficiently separated occupation correlations are exactly zero.
Fitting an exponential to those zeros would be meaningless.

## What changed after reading your papers

The papers confirm the state tensors and the distinction between square
and kagome behavior. They also exposed an important boundary detail.
A kagome cylinder can occupy different global dimer sectors, identified
in our tensors by even or odd virtual parity. The earlier width-2 result
selected the odd sector, while width 4 selected the even sector. Joining
those points did not give a comparison at fixed sector. The corrected
program selects and labels each sector explicitly; the two series should
be read separately. None of these runs inserts a vison, the extra
topological twist studied in the papers.

The revised analysis also checks several fitting intervals. The earlier
kagome lengths 0.624 and 0.525 cells came from a short interval and are
not established asymptotic lengths. The resolved x-direction scalar-bond
tail gives approximately **0.638 cells**, or **1.28 nearest-neighbor
spacings**, in the tested calculations. The y-direction estimate keeps
drifting with the interval and remains unresolved. The square estimate
near **1.30** describes an apparent power law over a limited interval;
larger surrounding tensor environments are still needed to establish the
bulk exponent.

A further check changed how small CTMRG modes are discarded. It changed
the tiny kagome correlations and some stricter runs did not converge.
The approximate 0.638-cell fit therefore remains dependent on contraction
settings; agreement at two nominal environment sizes is not enough to
assign it a precise physical value. The [paper recheck](../results/paper-recheck/report.md)
records the settings and failed diagnostics.

Do not compare the local decay length directly with the numbers around
0.65 and 1.01 in Fig. 9 of `1202.0947v4.pdf`. Those describe how differences
in **energy between topological sectors** shrink as the cylinder gets
wider. They measure a different quantity. That paper's Fig. 12 concerns
an effective boundary Hamiltonian, also a different observable.

## How to read the files

1. Start with [RESULTS.md](../RESULTS.md) for the current plots and summary.
   The earlier exact benchmarks are in
   [results/reference/report.md](../results/reference/report.md).
2. Read [the tail recheck](correlation-tail-audit.md) alongside
   `results/upstream/report.md` for the independent CTMRG comparison
   and its convergence tables. `results/upstream/kagome_strict/report.md`
   records the tighter kagome convergence check.
3. Consult the JSON files only if you want the numerical values, signs,
   units, eigenvalue diagnostics, fit windows, or exact software versions.
4. Use `README.md` to repeat a calculation. You can also ask the assistant
   to change a parameter or explain a particular curve.

For the paper-specific details, see the
[kagome interpolation audit](paper-1203-audit.md) and
[cylinder-sector audit](paper-1202-audit.md). The subsequently uploaded
**1009.1307v2** has also been read; see the
[square-paper audit](paper-1009-audit.md). It reports dimer exponent
**1.16(4)** and spin correlation length **1.35(1)**. Our existing square
data have not reproduced those numbers quantitatively, and the main dimer
measurement also needs its separation direction matched.

The [decay-diagnosis guide](decay-diagnostics.md) explains the local-slope,
branch and convergence checks now implemented for future lattice studies.
The earlier second figure was a replot of Schuch's published coordinates,
not an independent reproduction. New [diagnostic figures](../results/decay-diagnostics/report.md)
use our computed correlations and exact benchmarks.

The implemented states are fixed nearest-neighbor RVB states, not optimized
ground states of a Hamiltonian. A different lattice, longer-range dimers,
boundary sector, or choice of state amplitudes is a new physical model and
must be checked explicitly.
