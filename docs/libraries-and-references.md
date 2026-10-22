# Libraries, reference states, and a reproduction protocol

This project studies the **equal-amplitude nearest-neighbor spin-1/2 RVB
wavefunction**, and the equal-weight state of **orthogonal dimer coverings**.
Neither state is automatically the ground state of the nearest-neighbor
Heisenberg Hamiltonian. Ground-state optimization is a different calculation.

The two uploaded Schuch/Poilblanc papers and the Albuquerque–Alet square-lattice
paper have now been read and checked against the implementation. Their tensor
formulas support the implemented states, but figure comparisons require
specific observables and boundary conventions. See the
[square-lattice 1009.1307v2 audit](paper-1009-audit.md),
[1203.4816v3 audit](paper-1203-audit.md),
[1202.0947v4 audit](paper-1202-audit.md), and
[correlation-tail audit](correlation-tail-audit.md).

## Which correlation should be measured?

The square-lattice algebraic decay concerns **connected dimer/bond correlations**.
It does not mean that spin-spin correlations are algebraic. The accessible
[European Tensor Network School 2023 tutorial](https://github.com/jurajHasik/ETN-school-2023_peps-tutorial/blob/e0ad72f5c93caafa2c5fb892de42bb8b41ba98af/peps-solutions.ipynb)
explicitly distinguishes rapidly decaying spin correlations from increasingly
slow dimer correlations as the CTMRG environment dimension increases.

For the singlet RVB state use, for example,

\[
B_{ij}=\mathbf S_i\cdot\mathbf S_j,\qquad
C_B(r)=\langle B_0B_r\rangle-\langle B_0\rangle\langle B_r\rangle.
\]

Equivalently one can use the singlet projector
\(P_{ij}^{s}=1/4-\mathbf S_i\cdot\mathbf S_j\), which has the same connected
two-bond correlation. For orthogonal dimers the appropriate observable is the
occupation number \(n_{ij}\in\{0,1\}\), with
\(C_n(r)=\langle n_0n_r\rangle-\langle n_0\rangle\langle n_r\rangle\).
These are analogous observables in different Hilbert spaces, not numerically
identical operators. A physical singlet projector is not a label indicating
which bond occurs in one chosen valence-bond covering.

Figure 11 of Schuch et al., arXiv:1203.4816v3, instead specifies
`z ⊗ z`, with three-level `z = diag(1,-1,1)`. Its precise embedding and
selected legs must be matched before comparing amplitudes or fitted lengths.
The original `S_i · S_j` results are not a direct numerical reproduction of
that figure. Likewise, a correlation between two physical `S_z S_z` bonds is
not generally one ninth of the scalar-bond correlation: the four-spin
operators contain different symmetry channels.

Orthogonal coverings remove the overlaps between different coverings from the
norm. For diagonal occupations their equal-weight quantum state therefore
reduces to the uniform classical close-packed dimer ensemble. Physical singlet
coverings are nonorthogonal; the double-layer network retains their overlaps
and orientation-dependent signs. Replacing them by classical dimers changes
the state, even when it makes the contraction much simpler.

## Library survey

The links below point to project sources inspected during setup. Features are
the maintainers' documented capabilities, not a claim that every package has
been installed or benchmarked here.

| Library | Language | Fit for this project | Practical qualification |
|---|---|---|---|
| [quimb](https://github.com/jcmgray/quimb) | Python | Flexible tensors, arbitrary lattice graphs, PEPS and reference exact calculations; useful for small-system checks and new lattice geometry. | A contraction library does not supply the correct RVB tensor, observable, sector, or convergence study automatically. |
| [peps-torch](https://github.com/jurajHasik/peps-torch/tree/239a2e717a567df0ece80a5bce5443a6f78d048b) | Python/PyTorch | Particularly useful **existing reference implementation**: supplied square and kagome RVB tensors, CTMRG, transfer spectra, square dimer correlations, and regression tests. | Kagome's example prints spin correlations, not bond-bond correlations; a bond insertion is still needed. No energy optimization is required for a fixed RVB state. |
| [PEPSKit.jl](https://github.com/QuantumKitHub/PEPSKit.jl/tree/5f9370e5d753035b7181acd59bb103d0f06261a4) with [TensorKit.jl](https://github.com/QuantumKitHub/TensorKit.jl) | Julia | Strong longer-term choice for infinite PEPS: CTMRG and boundary MPS, explicit correlation-length/transfer-spectrum routines, generic unit cells, and symmetry-aware tensors. | A new language/toolchain and a custom RVB construction are required; this is an alternative backend, not a prerequisite for using the Python program. |
| [variPEPS](https://github.com/variPEPS/variPEPS_Python) | Python | Infinite-PEPS contraction and automatic-differentiation variational calculations; useful if the future task optimizes states. | The project's main workflow is variational optimization. Fixed-state RVB reproduction still needs a defined tensor and measurements. |
| [TeNeS](https://github.com/issp-center-dev/TeNeS) | C++ with Python tools | PEPS/CTM solver with optional MPI/OpenMP for larger calculations. | Requires a C++/CMake/BLAS toolchain; arbitrary fixed RVB inputs need additional integration. |
| [TeNPy](https://github.com/tenpy/tenpy) | Python | Established MPS/DMRG/TEBD workflows and useful independent finite-cylinder comparisons. | Its standard MPS workflows are not a drop-in 2D PEPS/CTMRG implementation. |
| [TensorTrack](https://github.com/quantumghent/TensorTrack) | MATLAB | Symmetry-aware tensor algebra and network contraction. | The README requires MATLAB R2020b+, Parallel Computing Toolbox, Deep Learning Toolbox, and a MEX-compatible C++ compiler. It is a toolbox, not a ready RVB reproduction script. |

For a first usable program, Python avoids a MATLAB license and makes exact
NumPy/SciPy checks straightforward. The included small-cylinder calculations
provide inspectable baselines; `peps-torch` supplies a useful independent
infinite-lattice comparison. PEPSKit/TensorKit is a credible later option if
symmetry handling and larger infinite systems become the priority.

## Existing independent RVB calculations

The inspected `peps-torch` commit is
`239a2e717a567df0ece80a5bce5443a6f78d048b`.
The inspected ETN tutorial commit is
`e0ad72f5c93caafa2c5fb892de42bb8b41ba98af`.
Commands below are reference commands to run **inside that upstream checkout**;
they do not install its dependencies and are not claims of successful runs in
this project. Consult the project's generated results for runs actually made.

### Square singlet RVB

The tutorial uses the supplied fixed state, without optimizing it:

```bash
python examples/j1j2/ctmrg_j1j2_c4v.py \
  --instate test-input/RVB_1x1.in \
  --j2 0.0 --chi 18 --corrf_r 20 --top_n 10 \
  --CTMARGS_ctm_max_iter 1000
```

Source: [square CTMRG example and regression test](https://github.com/jurajHasik/peps-torch/blob/239a2e717a567df0ece80a5bce5443a6f78d048b/examples/j1j2/ctmrg_j1j2_c4v.py).
The `TestRVB` regression test uses `chi=16`, `j2=0.5`, at most 200
iterations and expects energy per site `-0.47684229` with tolerance `1e-8`,
and vanishing magnetization. This is a **J1–J2 energy expectation in the fixed
RVB state**. It must not be compared directly with the `j2=0` command above.

The tutorial warns that higher `chi` can require **more than 1000 CTMRG
iterations**. A small energy change alone is not a guarantee that long-distance
correlations or the transfer gap have converged.

### Kagome singlet RVB

The upstream regression test supplies a D=3 simplex state:

```bash
python examples/kagome/ctmrg_spin_half_kagome.py \
  --ansatz IPESS --instate test-input/IPESS_KAGOME_D3_RVB.in \
  --bond_dim 3 --chi 18 --GLOBALARGS_dtype complex128 \
  --corrf_r 20 --top_n 8
```

Source: [kagome CTMRG example and `TestCtmrg_IPESS_D3_RVB`](https://github.com/jurajHasik/peps-torch/blob/239a2e717a567df0ece80a5bce5443a6f78d048b/examples/kagome/ctmrg_spin_half_kagome.py).
At `J1=1` and `chi=18`, the stored test reference is energy per physical spin
`-0.3931221584692804`, with tolerance `1e-6`. The six nearest-neighbor bond
expectations are approximately `-0.196561`; up/down triangle sums are
approximately `-0.589683`. These are upstream numerical regression values,
not exact analytical constants.

The [point-group-symmetric input](https://github.com/jurajHasik/peps-torch/blob/239a2e717a567df0ece80a5bce5443a6f78d048b/test-input/IPESS_PG_KAGOME_D3_RVB.in)
explicitly stores triangle tensor
\(T_{abc}=\epsilon_{abc}+\delta_{a2}\delta_{b2}\delta_{c2}\), with
\(\epsilon_{012}=+1\). Its site projector has nonzero entries
\(B_{s,s,2}=B_{s,2,s}=1\) for \(s=0,1\).
The ordering of triangle legs fixes the singlet orientations. This ordering
must be kept when moving between a graph description and a blocked tensor.
On the square lattice a consistent bipartite A-to-B singlet convention is a
useful choice; an arbitrary relabeling of singlet arrows can change relative
covering amplitudes.

The uploaded papers confirm these tensor definitions. Our counterclockwise
arrows reverse every singlet relative to the papers' clockwise convention,
giving the same overall phase to every complete covering. This global
reversal leaves normalized observables unchanged.

The kagome example prints spin correlations for all three sublattices and
horizontal/vertical transfer spectra. To measure bond correlations, insert
\(\mathbf S_0\cdot\mathbf S_1\otimes I_2\) into the three-spin physical
block, use the generic one-block correlation routine, and subtract the
product of one-point means. The upstream
[kagome model](https://github.com/jurajHasik/peps-torch/blob/239a2e717a567df0ece80a5bce5443a6f78d048b/models/spin_half_kagome.py)
already constructs this operator as `SSnnId`. Index ordering and the physical
distance per blocked step must be checked against the chosen lattice drawing.

### The paper's RVB-to-orthogonal interpolation

`rvb/interpolation.py` implements the explicit projector in Eq. (12) of
arXiv:1203.4816v3. Each site retains a spin and an additional two-state label
recording the adjacent triangle. At theta=0 the label factors out and the
usual spin RVB state is recovered; at theta=1 different dimer coverings
become orthogonal in the enlarged physical space. Intermediate theta is a
new tensor-network calculation, not a weighted average of the endpoint
correlations.

The implementation follows Eq. (12) directly and is checked against explicit
spin-and-label states. Its observables separately identify physical spin
bonds, spin-component bonds, dimer-location occupation, and specified
three-level embeddings. Matching one of these choices uniquely to Fig. 11
remains unresolved. An apparent limitation of the paper's scalar overlap
shorthand is documented with a six-site-loop check in the
[paper audit](paper-1203-audit.md); the tensor equation is the reference used
by the program.

## How to make the correlation claim convincing

1. **Fix the state and observable.** Record the lattice, singlet arrows,
   boundary conditions, sector, bond orientation, and physical distance units.
   Use connected bond correlations for the square-versus-kagome comparison.
2. **Check small exact cases.** Compare network contractions with explicitly
   enumerated coverings or spin amplitudes. Test normalization and one-point
   observables before attempting a tail fit.
3. **Distinguish the two limits.** An exactly contracted cylinder has a finite
   circumference error; a 2D CTMRG calculation has a finite environment
   dimension `chi` error. Increasing one does not remove the other.
   For kagome, also fix the virtual parity sector: `RVBTransfer` and
   `rvb_cylinder_benchmark` accept `sector="even"` or `sector="odd"`.
   The old unrestricted run selected odd parity at width 2 but even parity
   at width 4, so it was not a controlled fixed-sector circumference series.
   Outputs record the selected parity and absence of a vison insertion;
   these labels do not by themselves establish the paper's physical geometry.
4. **Converge the environment.** For a CTMRG calculation, start with modest
   `chi` and then increase it, for example `18, 27, 36, 54`, while monitoring
   the correlation curve, transfer spectrum, and fixed-point residual. These
   are suggested settings, not guaranteed sufficient values. Keep the RVB
   tensor fixed rather than optimizing its energy.
5. **Identify the observable's transfer channel.** For an eigenvalue with
   nonzero overlap with the observable, a mode has
   \(\xi=-a/\log|\lambda/\lambda_0|\), where `a` is the distance advanced
   by one transfer step. Degenerate leading eigenvalues can encode sectors;
   they do not by themselves prove nondecaying local correlations. A blindly
   selected second eigenvalue need not give the dimer correlation length.
6. **Test both tail descriptions.** Examine semilog and log-log plots, fit
   several distance windows, and separate oscillating/parity components.
   Do not fit values at numerical zero. Check that the inferred length or
   exponent is stable as the circumference or `chi` grows.
7. **Interpret finite-size tails carefully.** Even a critical 2D state can
   show an exponential tail on a finite cylinder or at finite `chi`. Square
   criticality requires a growing correlation scale and stable intermediate
   algebraic behavior as the approximation is improved. Kagome's expected
   short-range behavior requires a stable finite scale. The orthogonal
   kagome dimer model can have strictly vanishing separated occupation
   correlations, so fitting a positive exponential length to zeros is
   inappropriate.

A small-cylinder benchmark or a single low-`chi` curve is useful evidence
that the implementation works; it is not a complete thermodynamic-limit
reproduction or a precise literature exponent/correlation length.

## Primary literature read and further references

The following versions were supplied as PDFs by the user and read directly.
The earlier network-access limitation no longer prevents checking these three
sources.

- Norbert Schuch, Didier Poilblanc, J. Ignacio Cirac, and David Pérez-García,
  **“Resonating valence bond states in the PEPS formalism.”**
  [arXiv:1203.4816v3](https://arxiv.org/abs/1203.4816v3).
  Equations (1)–(3) define kagome tensors and orthogonal dimers; Eq. (12)
  defines their interpolation. Figure 11 plots a specified three-level
  observable, and Fig. 12 treats a transfer-eigenvalue ratio after accounting
  for two leading topological eigenvalues.
  [Detailed source/implementation check](paper-1203-audit.md).
- Didier Poilblanc, Norbert Schuch, David Pérez-García, and J. Ignacio Cirac,
  **“Topological and entanglement properties of resonating valence bond wave
  functions.”**
  [arXiv:1202.0947v4](https://arxiv.org/abs/1202.0947v4).
  Section II A distinguishes short-range spin correlations from critical
  square dimer correlations. Sections III–IV specify tensor and boundary
  sectors. **Fig. 9's lengths near 0.65 and 1.01 describe topological energy
  splittings versus circumference, not local dimer decay.** Figure 12 plots
  boundary-Hamiltonian coefficients rather than bulk dimer correlations.
  [Detailed audit and corrected fixed-sector results](paper-1202-audit.md).
- A. Fabricio Albuquerque and Fabien Alet,
  **“Critical Correlations for Short-Range Valence-Bond Wave Functions on the
  Square Lattice.”**
  [arXiv:1009.1307v2](https://arxiv.org/abs/1009.1307v2), identified in the
  uploaded PDF's margin as 2 November 2010. This is Phys. Rev. B **82**,
  180408 (2010), reference 28 in 1202.0947v4. It directly defines the same
  connected scalar-bond observable used by our square code, and reports a
  leading dimer exponent **1.16(4)**, spin length **1.35(1)**, and classical
  orthogonal-dimer exponent **2**. Its simulations use fixed winding sectors
  on square tori up to L=128; the selected bond-separation direction must
  also be matched. [Detailed source and numerical audit](paper-1009-audit.md).

The Albuquerque–Alet upload supplies the square benchmark previously missing.
Our finite-chi dimer estimate near **1.30** remains unextrapolated and does not
yet reproduce **1.16(4)**. The saved square spin length also drifts with chi;
for example, the r=10…24 fit reaches about **1.162** at chi=128, compared
with the paper's **1.35(1)**. Matching qualitative decay and local energy
does not establish quantitative convergence of either tail. No further
paper is required to define this reproduction target.

Additional useful sources, not yet read in this project, are:

- Ying Tang, Anders W. Sandvik, and Christopher L. Henley,
  **“Properties of resonating-valence-bond spin liquids and critical dimer
  models.”**
  [arXiv:1010.6146](https://arxiv.org/abs/1010.6146), identified in reference
  36 of the uploaded Albuquerque–Alet paper.
  Relevant for distinguishing square singlet-RVB and classical-dimer
  correlation exponents. Do not assume their exponents are equal.
- Grégoire Misguich, Denis Serban, and Vincent Pasquier,
  **“Quantum dimer model on the kagome lattice: Solvable dimer-liquid and
  Ising gauge theory.”**
  [arXiv:cond-mat/0204428](https://arxiv.org/abs/cond-mat/0204428),
  identified by reference 11 of the uploaded 1203.4816v3 paper.
  Relevant for the particularly short-ranged orthogonal kagome dimer case.

For every numerical comparison, record the paper version, figure,
observable, sector, geometry, fit window, and distance convention. Agreement
on short-range versus critical behavior alone does not establish a
point-by-point reproduction.
