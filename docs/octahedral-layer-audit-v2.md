# Octochlore layer: coincident apex projection and square–octagon opening

This corrects the interpretation in `octahedral-layer-audit.md`. The intended
projection uses the **cubic-centre grid with its edge-midpoint spins**, not the
network of equatorial octahedron perimeters. The two dummy slots arise from two
different vertical physical sites, not from splitting a single physical spin.

## Exact top view

For centres \(v=(m,n,0)\), the six sites of one octahedron are

\[
 v\pm e_x/2,\quad v\pm e_y/2,\quad v\pm e_z/2.
\]

Projection \(P(x,y,z)=(x,y)\) leaves the four equatorial sites at the midpoints
of the horizontal edges of the cubic-centre grid. The two distinct apex sites
both project to \(P(v\pm e_z/2)=(m,n)\). Thus the proposed dummy pair really is
**two physical sites directly above one another**, sharing one location in the
top view. The spin-1 extension changes the on-site Hilbert space, not this
geometric fact.

Choose the checkerboard family of centre-grid plaquettes whose lower-left
corner \((i,j)\) has odd \(i+j\). At each centre, two selected plaquettes meet in
opposite quadrants. The incident horizontal ports divide as follows:

| Centre parity | Upper face frame | Lower face frame |
|---|---|---|
| \(m+n\) even | \(\{y_+,x_-,z_+\}\) (NW) | \(\{y_-,x_+,z_-\}\) (SE) |
| \(m+n\) odd | \(\{y_+,x_+,z_+\}\) (NE) | \(\{y_-,x_-,z_-\}\) (SW) |

Each frame consists of two adjacent equatorial sites and one apex. The frames
are a complementary pair of octahedron faces and use all six local sites once.
In projection their apex sites coincide; resolving the coincidence separates
the two triangular frames joined by the dummy link. Checkerboard labels can be
interchanged without changing the construction.

The grey square in this construction is a plaquette of the **centre grid**.
Its four physical equatorial spins lie at its edge midpoints. It is not the
diamond obtained by joining those four physical spins to one another. With the
cubic x and y axes horizontal and vertical, the relevant grey squares are
upright.

## Continuous drawing deformation

For a selected plaquette with centre \(w=(i+1/2,j+1/2)\) and a corner \(v\), let
\(s(v,w)=\operatorname{sign}(w-v)\), with both components equal to \(\pm1\).
For \(0\leq\delta<1/2\), set

\[
 q(v,w)=v+\delta s(v,w),\qquad
 d(v,w)=v+\frac{\delta}{2}s(v,w).
\]

Here \(q\) is the grey base-graph vertex; \(d\) is the orange physical apex-slot
position in the drawing. Join the four \(q\)'s around each selected plaquette
to make its grey square; put a blue physical site at each grey side midpoint.
Join the two grey vertices belonging to opposite selected plaquettes around
one original centre \(v\) by a dummy edge. The two orange sites lie along this
dummy. At each grey vertex join its orange site to its two adjacent blue sites
with a green dashed triangular frame.

For \(\delta=0\), blue sites are exactly the original horizontal spin positions
and both orange sites at each centre are exactly their common projected
position. Thus the zero-opening limit reproduces the physical top view. For
\(\delta>0\), every selected square stays upright and shrinks. The complementary
checkerboard squares gain a short diagonal at each of their four corners and
become octagons. The resulting trivalent base graph is the square–octagon
(4.8.8) graph. The ordinary grey side has length \(1-2\delta\) and the dummy has
length \(2\sqrt2\delta\), so

\[
 \delta=\frac{1}{2(1+\sqrt2)}
\]

produces equal base-edge lengths. This is a **display deformation**: the actual
physical top projection has coincident apices, not separated in-plane apices.
The blue positions also deform continuously as the grey squares shrink.

For a periodic centre layer with \(N\) centres, there are \(2N\) horizontal
physical sites and \(2N\) apex sites. The opened construction has \(2N\) ordinary
edge spins, \(N\) dummy edges with two spins each, and \(2N\) triangle frames:
the site and frame counts agree.

## What this establishes, and what remains

The construction establishes the user's **geometric site-incidence relation**:
a centre layer viewed along z has two superposed dummy sites, and opening those
pairs gives the decorated square–octagon geometry. No physical equatorial spin
is duplicated. The selected green frames are geometry, not simultaneously
occupied trimers.

It does not by itself identify the 3D quantum covering state with the 2D RVB
state:

- There are eight possible octahedron face trimers in the full 3D model. The
  figure displays a complementary pair of face frames; the other six faces
  remain allowed. In particular, choosing frames for the picture does not
  delete or restrict those local terms.
- The original 3D rule permits one face trimer at a centre, but forbids the
  simultaneous occupation of both complementary faces by trimers.
- The 3D diameter dimers \(D_x,D_y,D_z\) pair opposite octahedron sites. The
  vertical diameter \(D_z\) connects the two apex sites and therefore has the
  same incidence as the 2D internal dummy singlet. The horizontal diameters
  connect opposite ports, whereas 2D triangle singlets pair adjacent ports or
  an adjacent port with its dummy slot. A face trimer occupies three sites.
- The current 3D extension uses spin-1 dimer and trimer singlets; the earlier
  2D RVB construction used spin-1/2 dimers. A reduction must state which local
  Hilbert-space map and which covering amplitudes it uses.
- Consecutive centre layers share physical vertical sites: the upper apex of
  \((m,n,k)\) is the lower apex of \((m,n,k+1)\). A foliation must therefore
  account for these shared degrees of freedom and their covering constraints.
  This picture alone does not furnish independent 2D tensor factors.

A useful next calculation is to write the 3D local support tensor in the
selected adjacent-port basis and examine how the shared vertical leg couples
the planar layers. This can test whether a constrained sector or a local
projection realizes the known 2D tensor. It is not justified to transfer a
2D correlation length to the full 3D state from the geometric picture alone.

## Reproduction

`python scripts/plot_collaborator_layer.py` writes
`results/collaborator-notes/layer-relation.pdf`, its PNG preview, and
`layer-projection-coordinates.json`. The two panels show a local octahedron and
the exact layer top view, connected by a projection arrow. The square–octagon
opening is discussed in the text above rather than drawn as a third panel. Blue denotes
equatorial physical spins, orange the two apex sites, and green dashed lines
the selected geometric frames. The top-view panel uses concentric orange marks
specifically to display the two distinct sites at identical x,y.
