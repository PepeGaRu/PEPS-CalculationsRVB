# Octochlore covering model: replacement source for collaborator notes

This section incorporates the scientific content of José Garre Rubio,
*Put the objects on octahedra / Octahedral coverings and X-cube cage nets*,
19 October 2026, and separately specifies the subsequent spin-1 extension.

## Geometry and tetrahedral origin

The octochlore lattice consists of corner-sharing octahedra with simple-cubic
centres. For each \(v\in\mathbb Z^3\), let

\[
O_v=\operatorname{conv}\{v\pm\hat x/2,\ v\pm\hat y/2,\ v\pm\hat z/2\}.
\]

The shared corner \(v+\hat\mu/2\) belongs to \(O_v\) and
\(O_{v+\hat\mu}\): it is one physical site, one cubic bond, and has two
possible owners. Centres are graph vertices, not additional physical sites.
A cube of centres has eight centres, twelve bonds and twelve shared sites;
each octahedron has six neighbouring octahedra, one per cubic direction.

The six edge midpoints of a regular tetrahedron form a regular octahedron.
For example the tetrahedron vertices \((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)\)
have edge midpoints \(\pm\hat x,\pm\hat y,\pm\hat z\); a common rescaling
gives the coordinates above. Its eight faces are four old-tetrahedron-face
triangles and four old-vertex-star triangles. Replacing each tetrahedron by
this **new inscribed octahedron** retains the six covering sites and the
specified cubic two-owner incidence construction: only the supporting
polyhedron changes. This is not an exchange of the original cells of the
tetrahedral–octahedral tiling: centres of all its original octahedra form an
FCC lattice, not the required simple-cubic graph.

Figure: `results/collaborator-notes/geometry-origin.pdf` (blue dots are sites,
orange dots are centres; tetrahedron face/vertex-star families have 4 members
each, with one member of each shaded).

## Sixteen allowed local configurations

At an octahedron allow **any subset of its three diameter dimers**, or
**exactly one of its eight triangular-face trimers**. Every physical site is
covered exactly once globally. The cubic-bond arrow points toward the
octahedron whose object covers that site. Thus a diameter consumes an opposite
pair of sites, while a face trimer consumes one site from each axis.

| Objects at one octahedron | Count | Incoming sites |
| --- | ---: | --- |
| Zero, one, two, or three diameter dimers | \(1+3+3+1=8\) | Any union of opposite pairs |
| One face trimer | \(2^3=8\) | One of \(\pm x\), one of \(\pm y\), one of \(\pm z\) |

The full catalogue is
\(\varnothing,D_x,D_y,D_z,D_x+D_y,D_x+D_z,D_y+D_z,D_x+D_y+D_z\)
and \(T_{+++},T_{++-},T_{+-+},T_{+--},T_{-++},T_{-+-},T_{--+},T_{---}\).
Figure: `results/collaborator-notes/octahedral-catalogue.pdf`.
Filled sites mean incoming ownership, gray sites outgoing ownership;
face signs are ordered \((x,y,z)\). These are physical support diagrams,
not reference-dependent cage bits. Diameter chords may cross at a centre,
which is not a covering site.

**Two complementary face trimers remain forbidden.** If included, the
all-six-incoming support would have five realizations: the three diameters,
or one of four opposite-face pairs. Excluding those pairs makes each
allowed support determine its objects uniquely.

## Reference-XOR correspondence with X-cube cage nets

For incoming ownership bits \(a^\pm_{v,\mu}\), exactly the above sixteen
supports satisfy

\[
a^+_{v,x}\oplus a^-_{v,x}
=a^+_{v,y}\oplus a^-_{v,y}
=a^+_{v,z}\oplus a^-_{v,z}.
\]

Fix a complete reference covering \(R\). For another full covering \(C\),
set \(n_e=a_{v,e}(C)\oplus a_{v,e}(R)\). This is independent of which end of
the cubic bond is used, because both ownership bits complement at the other
end. The difference bits obey the same equality of opposite-axis parities;
equivalently there is even parity in each coordinate plane through every
vertex. They therefore form an X-cube cage net. Conversely, reversing the
reference arrows on any such net preserves one owner per site and produces
an allowed support at every octahedron. The unique local-object assignment
then gives a full covering. This proves a bijection, for fixed reference
and compatible boundary conditions. The combinatorial correspondence is
the same as in the tetrahedral representation.

## Specified spin-1 singlet extension

The original geometric support model alone does not define a quantum spin
state. In the extension studied here each site carries spin 1 and every
object is a normalized singlet. In the Cartesian basis,

\[
d_{ab}=\delta_{ab}/\sqrt3,\qquad
t_{abc}=\epsilon_{abc}/\sqrt6,\qquad\epsilon_{xyz}=+1,
\qquad
|\Psi\rangle=Z^{-1/2}\sum_C|C\rangle.
\]

For a face \((v+s_x\hat x/2,v+s_y\hat y/2,v+s_z\hat z/2)\), trimer indices
are ordered \(x,y,z\), independently of the signs. The dimer is symmetric
and the trimer antisymmetric under exchange; changing face phases can
change the superposition. \(|C\rangle\) is the product of normalized
object singlets and \(Z=\sum_{C,C'}\langle C'|C\rangle\), including their
nonorthogonal overlaps. This is an equal-coefficient covering sum, with no
Hamiltonian or ground-state claim. An exact virtual-dimension-four tensor
network represents it. Uniform classical covering correlations do not
determine its quantum correlation length.

## Literature from the original note

1. O. Benton and R. Moessner, “Topological Route to New and Unusual Coulomb
   Spin Liquids,” *Phys. Rev. Lett.* **127**, 107202 (2021),
   [arXiv:2103.10817](https://arxiv.org/abs/2103.10817).
2. M. Stern, M. D. Burke, M. J. P. Gingras, J. Romhányi and K. T. K. Chung,
   “Fracton Spin Liquid and Exotic Frustrated Phases in Ising-like Octochlore
   Magnets” (2026), [arXiv:2603.12313](https://arxiv.org/abs/2603.12313),
   Sections II and V. Their Ising cage-net liquid has **ten** allowed local
   ground-state arrow patterns; it does not use this **sixteen**-state
   covering catalogue. Geometry can be shared without assuming identical
   quantum phases or correlations.

The two geometry/support figures are generated by
`scripts/plot_collaborator_octahedral_catalogue.py`, using the actual sixteen
supports from `rvb/octahedral.py` and verifying unique supports and parity.
