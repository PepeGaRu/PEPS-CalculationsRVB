# Deciding whether a correlation is exponential or algebraic

We measure a specified connected correlation,

\[
C(r)=\langle O(0)O(r)\rangle-\langle O(0)\rangle\langle O(r)\rangle.
\]

The operator, bond orientations, displacement direction, distance units,
state amplitudes and boundary sector are part of the question. A fixed
dimer covering, a coherent RVB superposition, and an orthogonal-dimer
ensemble are different states. Dimer and spin correlations in the same
state can have different decay laws.

## Two local checks

For resolved neighboring points on the same oscillation/parity branch,

\[
p_{\rm eff}=-\frac{\Delta\log|C|}{\Delta\log r},\qquad
\kappa_{\rm eff}=-\frac{\Delta\log|C|}{\Delta r}.
\]

| Decay envelope | Effective exponent | Inverse decay length |
|---|---|---|
| \(A r^{-p}\) | Approaches the constant \(p\) | Approaches zero as \(p/r\) |
| \(A e^{-r/\xi}\) | Grows approximately as \(r/\xi\) | Approaches the positive constant \(1/\xi\) |
| \(A r^{-p}e^{-r/\xi}\) | Approximately \(p+r/\xi\) | Approximately \(1/\xi+p/r\) |

For finite point spacing these derivatives are finite differences; the exact
identity is \(p_{\rm eff}=\kappa_{\rm eff}\Delta r/\Delta\log r\).
A stable positive inverse length can resolve an exponential with an
algebraic prefactor. A pure power law has no finite asymptotic exponential
length. Negative local slopes, oscillatory zeros, or a plateau require
inspection; they are not silently converted into positive lengths.

These checks complement semilog and log-log plots. A short curve can look
nearly straight on both axes, and a small fit residual alone cannot decide
the large-distance behavior.

## Numerical protocol

1. Validate the lattice, tensor, singlet signs and measurement on small
   systems. Record which physical distance one transfer step represents.
2. Subtract one-point means consistently with the two-point contraction.
   Where possible, compare direct subtraction and centered insertions.
3. Keep oscillating branches separate. Track signs before taking absolute
   values. Exclude nodes and unresolved values; never connect a slope across
   a point removed for numerical noise. A leading term can vanish by geometry
   or symmetry, exposing a different subleading exponent.
4. Examine both local checks and fit several distance windows to power,
   exponential, and power times exponential envelopes. Fits are descriptive;
   a mixed model has an extra parameter and often cannot separate its power
   and exponential terms over a short window. Report its conditioning.
5. Repeat as cylinder width, system size, environment dimension and
   truncation thresholds improve. Require the correlation curve and fitted
   parameters to stabilize in a growing resolved distance interval. Passing
   an iteration convergence test at one setting is insufficient.
6. Check observable-coupled transfer eigenmodes within a fixed sector.
   \(\xi=-a/\log|\lambda/\lambda_0|\) applies to a mode that actually
   couples to the chosen observable. Leading topological degeneracies or
   other spin channels must not be substituted for the dimer channel.

Finite cylinders and finite CTMRG environments can turn an algebraic plane
correlation into an eventual exponential. A decay length growing with these
cutoffs is evidence of a limitation of the approximation, not by itself a
finite physical length. On periodic finite samples, wrap-around produces
another distortion. One should use a controlled finite-size scaling or a
range well below the system size; the conformal-cylinder chord formula is
not a universal replacement for distance on every two-dimensional torus.

A numerical zero establishes only that a correlation is unresolved.
An exact finite-range result requires independent analytical or exact
combinatorial support, as in the orthogonal kagome occupation benchmark.
A nonzero plateau may indicate order, mixed sectors, incomplete subtraction,
or an unconverged contraction. Do not remove it by an arbitrary fitted offset
and then declare a decay law.

## What the new paper teaches

Albuquerque and Alet, [1009.1307v2](https://arxiv.org/abs/1009.1307v2), use
the same connected scalar spin-bond operator as our original calculations.
They report a leading dimer exponent **1.16(4)**, a subleading exponent
**2.53(5)** for a channel with a vanishing leading term, and a separate
spin length **1.35(1)**. Their finite tori extend to **128×128** with fixed
winding sectors. The winding numbers produce cosine/sine modulations;
dividing by these factors near nodes would amplify noise, and the authors
exclude such points. The large-distance upturn moves with system size.

Their main zero-winding fit uses parallel bonds separated perpendicular to
the bond direction. Our current square CTMRG dimer measurement uses
collinear bonds. The state, operator, geometry and convergence must all be
matched before claiming a direct quantitative reproduction. See the
[paper audit](paper-1009-audit.md).

## Reusable program

`rvb.diagnostics.diagnose_decay` analyzes real, connected correlations with
explicit distance windows, branch labels, numerical floors and optional
absolute error estimates. It returns point masks, constant-sign segments,
local slopes and descriptive fits in JSON-compatible form. It does not
automatically certify a bulk exponent, finite length, or statistically
preferred model. A supplied error scale is used to screen unresolved data,
not to manufacture independent statistical error bars for correlated CTMRG
outputs. Monte Carlo data with uncertainties require a separate analysis
accounting for autocorrelation and covariance.

Run `python scripts/analyze_decay.py` to diagnose the saved calculations and
the exact orthogonal benchmarks. The script writes the source hashes,
parameters, fit tables and actual-data plots to `results/decay-diagnostics/`.
It does not fit extracted Schuch figure points as though they were our data.

For a new lattice we will first specify the allowed bonds and state weights,
singlet versus orthogonal coverings, the measured bonds and directions,
and boundary conditions. The outcome may be exponential, algebraic, ordered,
exactly finite range, or unresolved within available precision. A reliable
program must allow that last outcome.
