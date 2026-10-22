# Spin-1 singlet coverings on the octahedral lattice

This defines a quantum extension of the supplied sixteen-state covering
model. Every physical site at a cubic-bond midpoint carries spin 1. Each
allowed diameter dimer and triangular-face trimer is placed in its
normalized total-spin-zero state. We take equal coefficients for all
complete, owner-labelled coverings, with the phase convention below.
No Hamiltonian or ground-state claim is assumed.

The subsequent [infinite-cylinder calculation](../results/octahedral-spin1-cylinder/report.md)
now resolves physical singlet-projector decay: xi=2.4901024 on a periodic
2x2 transverse section and xi=0.6381923 for the same transverse observable
on 2x3. These are quantum finite-width results. Their width dependence
prevents a three-dimensional bulk extrapolation.

## Local singlets and their phases

Use a Cartesian basis \(\{|x\rangle,|y\rangle,|z\rangle\}\) in which a
spin-1 transforms as a three-dimensional vector. The normalized tensors
are

\[
 d_{ab}=\frac{\delta_{ab}}{\sqrt3},\qquad
 t_{abc}=\frac{\epsilon_{abc}}{\sqrt6},\qquad \epsilon_{xyz}=+1 .
\]

Both are invariant under simultaneous spin rotations, and both have unit
norm. Two spin-1 particles have one singlet; three spin-1 particles also
have exactly one singlet. In the usual magnetic basis the dimer can
equivalently be written
\((|1,-1\rangle-|0,0\rangle+|-1,1\rangle)/\sqrt3\), up to an overall phase
associated with the Cartesian-basis convention.
When converting the full covering sum to a magnetic basis, both local
tensors must be transformed consistently. Independently rephasing an
object type can alter relative coefficients when coverings contain
different numbers of dimers and trimers.

The dimer is symmetric under exchanging its sites. The trimer is
antisymmetric: exchanging two site labels reverses its sign. For every
face with vertices
\[
 v+s_x e_x/2,\quad v+s_y e_y/2,\quad v+s_z e_z/2 ,
 \qquad s_\mu=\pm1 ,
\]
the first, second and third trimer indices refer to the **x, y and z
vertices in that order**, irrespective of the face signs.

This is a specified phase convention for the superposition. Ordering the
vertices by an outward face orientation instead multiplies this trimer
by \(s_xs_ys_z\). With coefficients held equal, such face-dependent sign
changes can change the quantum state. They should not be silently
discarded as an irrelevant overall phase.

For a covering \(C\), let \(|C\rangle\) be the tensor product of its
normalized local singlets. The proposed state is

\[
 |\Psi\rangle=Z^{-1/2}\sum_C |C\rangle,\qquad
 Z=\sum_{C,C'}\langle C'|C\rangle .
\]

The allowed coverings are exactly those of the classical code: any
subset of the three diameter dimers on an octahedron, or one face
trimer; two complementary trimers on one octahedron remain forbidden.
Coverings with different object labels can represent nonorthogonal, or
on exceptionally small periodic clusters even identical, spin states.
Their amplitudes still add in the specified sum.

## An exact tensor-network representation

There is a direct decorated-cubic tensor network with virtual dimension
four. Each of an octahedron's six legs has a vacuum state \(0\) and
three Cartesian spin states \(x,y,z\).

At an octahedron, a rank-six tensor sums the sixteen allowed local
supports. A dimer inserts \(\delta/\sqrt3\) on its opposite legs; a
trimer inserts the axis-ordered \(\epsilon/\sqrt6\) on its face legs.
All unused legs are vacuum. Multiple allowed diameter dimers insert
the product of their normalized tensors.

At every shared physical site, connect its two octahedral legs using

\[
 P^\sigma_{ab}
 =\delta_{a0}\delta_{b\sigma}+\delta_{b0}\delta_{a\sigma},
 \qquad \sigma=x,y,z .
\]

This enforces exactly one owner of each physical spin. Contracting
these tensors gives precisely the unnormalized covering sum above,
including interference between different coverings. Equivalently,
blocking the three outgoing physical sites with each cubic centre
gives a three-dimensional PEPS with three spin-1 sites per block and
virtual dimension four.

The representation is exact; contracting its norm on a large
three-dimensional lattice is a different computational problem from
solving the classical binary constraints. The classical solver does
not include quantum overlaps, and its finite-range occupation
correlations cannot be carried over to this state.

## Physical observables

A classical label saying that a covering contains a particular dimer
is not the same as a physical singlet-projector measurement when
covering states overlap. We should measure specified spin operators
and their connected correlations.

For two spin-1 sites, write \(X=\mathbf S_i\cdot\mathbf S_j\). Their
singlet projector is

\[
 P^{(0)}_{ij}=\frac{X^2-1}{3}.
\]

The eigenvalues of \(X\) in total-spin sectors \(0,1,2\) are
\(-2,-1,+1\), respectively. For three spin-1 sites, set
\(Q=(\mathbf S_i+\mathbf S_j+\mathbf S_k)^2\). The unique trimer
singlet projector is

\[
 P^{(0)}_{ijk}
 =\left(1-\frac Q2\right)
  \left(1-\frac Q6\right)
  \left(1-\frac Q{12}\right)
 =|t\rangle\langle t|.
\]

An isolated singlet dimer has \(\langle\mathbf S_i\cdot\mathbf
S_j\rangle=-2\). Inside an isolated singlet trimer each pair has
\(\langle\mathbf S_i\cdot\mathbf S_j\rangle=-1\), and that pair is in
total spin 1, not a two-site singlet.

Suitable quantities for the full state include spin correlations,
dimer-projector correlations and face-trimer-projector correlations.
A decay law requires enough spatial separations and systematic size
or contraction-convergence checks; a local overlap calculation cannot
supply a bulk correlation length.

## Exact overlap and small-cluster checks

Start from the covering with a negative-face trimer on every
octahedron, and flip the twelve arrows of one elementary cube. On
the affected 24 spins this replaces eight trimers by twelve dimers.
With the axis-order convention,

\[
 \langle C_{\rm dimers}|C_{\rm trimers}\rangle
 =\frac{24}{3^6\,6^4}=\frac1{39366}.
\]

The contraction reduces to the cube graph with an epsilon tensor at
each of its eight vertices. There are 24 proper three-edge-colourings;
all contribute positively with this axis order. The denominator
comes from the twelve normalized dimers and eight normalized trimers.
This provides a direct nonorthogonal-overlap benchmark without
constructing a \(3^{24}\)-component state vector.

A periodic \(2\times2\times2\) cluster is especially exceptional.
Its classical trimer indicators must all agree. Of its 8,192
owner-labelled coverings, 4,096 contain only dimers and 4,096 only
trimers. All 4,096 dimer coverings represent the same spin product
state: opposite geometric chords can have identical physical
endpoints on a period-two line.

Let this normalized twelve-pair product be \(|D\rangle\), and let
\(|T\rangle\) be the eight-trimer negative-face reference. Then the
full unnormalized covering sum on this torus is

\[
 |\Psi_{\rm un}\rangle=4096|D\rangle+G|T\rangle,\qquad
 G=\prod_{e=1}^{12}(I+\mathrm{Swap}_e),
\]

where each swap exchanges the two physical sites in one of the twelve
period-two coordinate lines. This special reduction makes exact
quantum checks possible but also enhances the dimer contribution
through a geometric identification absent on larger lattices.
Results on this cluster must therefore be presented as validation
data, not evidence for a three-dimensional bulk decay law.

## What changes for spin one-half

For three spin-one-half particles,

\[
 \tfrac12\otimes\tfrac12\otimes\tfrac12
 =\tfrac12\oplus\tfrac12\oplus\tfrac32 .
\]

There is no local singlet. The total-spin-one-half sector has
**two distinct copies** of a doublet, so selecting total spin one-half
does not even fix the trimer's internal state. One must choose its
multiplicity degree of freedom, for example through an intermediate
pair coupling or a chirality convention, as well as specify how the
remaining effective spins on different trimers are combined.

Total spin three-halves has one quartet, but it too leaves an
effective spin requiring an inter-trimer state. Simply choosing one
magnetic component produces a different construction and generally
breaks spin-rotation invariance.

Spin 1 is therefore the direct choice for a covering soup whose
individual dimers and trimers are all singlets. Spin one-half is a
reasonable alternative when the intended physics includes residual
effective spins carried by trimers, but it needs extra model data
and is not automatically simpler to calculate.
