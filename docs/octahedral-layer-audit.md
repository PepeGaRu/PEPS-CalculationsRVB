# Octochlore layers and the square–octagon construction

The supplied notes define octahedra at cubic centres `v in Z^3`, with
physical sites `v +/- e_mu/2`. This audit separates a precise geometric
relation from a prospective quantum foliation. It uses that definition,
not the different FCC arrangement of original octahedra in a
tetrahedral–octahedral tiling.

## Geometric statement that can be drawn exactly

Take the octahedra with centres `(m,n,z0)` and look along z. Their equatorial
vertices are `(m +/- 1/2,n)` and `(m,n +/- 1/2)`. Keep only equatorial
perimeter edges, and omit the upper and lower apex sites and their edges.
The resulting graph is a regular square lattice, with a checkerboard
choice of faces corresponding to the octahedra. Adjacent marked squares
share a corner. It is the medial graph of the horizontal cubic-centre
lattice; its vertices are the horizontal physical spin sites.

For an upright drawing use

`(u,v) = ((x+y)/sqrt(2), (y-x)/sqrt(2))`.

The marked equatorial squares are then axis-aligned. Their sides are
actual equatorial octahedron edges, rather than a different square made
from a decoration's triangle bonds.

Resolve each shared corner into two vertices, one belonging to each
adjacent marked square, and join the pair by a new dummy edge. For example,
replace the contact at `(m+1/2,n)` by

`(m+rho,n) --- (m+1-rho,n)`, with `0<rho<1/2`.

Similarly split vertical contacts. The marked faces remain squares and
the complementary faces acquire eight sides: this is exactly the
trivalent square–octagon (4.8.8) graph. The regular-metric choice is
`rho = 1/(2+sqrt(2))`; other rho have the same graph. Contracting every
dummy edge is the inverse operation and recovers the equatorial square
graph. This is the same graph operation as the user's 2D dummy-link map.

Thus a useful accurate sentence is: **A resolved equatorial projection of
one octochlore layer is the square–octagon graph.** A literal top view of
the full layer includes projected apex sites and spokes and is not just
4.8.8; the omission and contact resolution must be indicated.

The figure [layer-relation.pdf](../results/collaborator-notes/layer-relation.pdf)
shows the construction. Its split black vertices are base-graph vertices,
not the physical spin slots of the subsequently triangle-decorated 2D
model. The coordinate operation alone does not split a spin-1 Hilbert space.

## A separate exact planar parity relation

Choose a full reference covering and let n be the difference of ownership
arrows modulo two. These reference-difference cage bits agree at the two
ends of each shared edge. Write

`p_mu(v) = n_(v,+mu) XOR n_(v,-mu)`.

The 16-state ownership rule implies the same homogeneous cage constraint,
`p_x = p_y = p_z`. It implies

`n_(+x) XOR n_(-x) XOR n_(+y) XOR n_(-y) = 0`,

and similarly in the yz and zx planes. Hence the horizontal edges in each
coordinate plane form an even subgraph of that plane's **cubic-centre
square lattice**. This centre graph and the equatorial medial graph above
are both square graphs but have different vertices and different site
densities; a drawing should not silently identify their physical sites.

The layer patterns are not independent: the vertical condition enforces
`n_(+z) XOR n_(-z) = p_x = p_y`. A turn in an xy loop has `p_x=p_y=1`
and therefore needs odd vertical incidence. An isolated xy plaquette loop
with no vertical edges is not an allowed 3D cage net. Two identical
plaquette loops on adjacent z layers, connected by their four vertical
corner edges, form the allowed elementary cube cage. This supplies a
concrete layer coupling and a natural starting point for a foliation study.
Direct enumeration of the sixteen allowed local supports verifies that
their horizontal projections are precisely the eight even four-port
supports, each with two possible vertical completions.

Any finite-plane even-edge pattern can be extended to a 3D allowed pattern
by copying it onto a second plane and connecting matching odd-pair-parity
vertices vertically. On periodic vertical lines the total turn parity must
be even. This extension is a kinematic statement; it is not factorization
into independently chosen 2D layers.

## What remains to establish a quantum relation

The 2D state uses spin-1/2 singlet dimers and doubled dummy endpoints; the
3D state uses physical spin 1 and both two-spin and three-spin singlets.
Their local Hilbert spaces, overlaps, and relative singlet phases differ.
Even the natural representation of a spin 1 as the symmetric subspace of
two spin-1/2 particles does not directly identify the models: the 2D
internal dummy singlet lies in the antisymmetric spin-zero subspace.

It is therefore justified to propose the resolved geometry and planar
parity equations as a way to organize a layer construction. A claim of
quantum foliation additionally requires an explicit map of states or
Hamiltonians, including the interlayer couplings and treatment of the
spin representations. No such equivalence or phase classification follows
from the current projection or from the measured local correlations.

## Suggested brief collaborator wording

“Viewed along a cubic axis, the equatorial skeleton of one octochlore layer
is a square lattice. Splitting each shared equatorial corner into two
vertices linked by a dummy bond produces the square–octagon graph; dummy
contraction reverses this construction. The 3D parity rule also gives even
edge patterns in every coordinate plane, coupled by the vertical parity.
These geometric and kinematic relations motivate a layer construction,
but do not yet identify the spin-1 3D singlet state with stacked spin-1/2
2D RVB states.”
# Historical audit: superseded layer interpretation

The collaborator-note revision uses the intended projection of two distinct
vertical apex sites onto each cubic centre. See
[the corrected layer audit](octahedral-layer-audit-v2.md). The equatorial-perimeter
construction below is a different graph interpretation retained only as a record.
