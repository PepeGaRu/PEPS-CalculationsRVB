# Contracting dummy edges: the user's square-loop map

The user's instruction to erase the dummy edges naturally means **contracting
each dummy edge to a vertex** after constructing reference-relative loops on
the square–octagon graph. This differs from the earlier figure, which collapsed
each entire small square and retained the dummy links as coarse edges.

The physical spin lattice and its computed correlations are unchanged by this
clarification. This is a map of loop drawings, not a contraction of physical
spin degrees of freedom or a change to the singlet wavefunction.

## The two square graphs

For N simulation cells, each containing one small square:

| Construction | Square-graph vertices | Square-graph edges | Coverings per reference-relative edge pattern |
|---|---:|---:|---:|
| Contract the dummy edges (user's map) | 2N | 4N surviving small-square edges | 1 |
| Collapse whole small squares (previous auxiliary projection) | N | 2N dummy links | 2^N |

Thus the earlier warning about losing 2^N local choices applies to the second
construction. It does **not** apply to contracting dummy links. The two
square graphs have different vertices and lattice vectors.

## Geometry and a direct inverse

Label a small-square center by c=(x,y), and its base corners E,N,W,S. Contract
E(c)—W(c+ex) to H(c)=(x+1/2,y), and N(c)—S(c+ey) to V(c)=(x,y+1/2).
The surviving small-square perimeter becomes

`H(c) — V(c) — H(c−ex) — V(c−ey) — H(c)`.

These edges form a square lattice rotated by 45 degrees. If (X,Y) are the
display coordinates, `(U,V)=(X+Y−1/2, Y−X+1/2)` gives integer coordinates
on that square lattice. The original simulation translation by ex maps to
the square-lattice displacement (1,−1); ey maps to (1,1). Reported numerical
lengths retain their simulation-cell units and are not silently reassigned
to steps along the contracted square's axes.

Choose a valid reference covering. Its binary arrow differences from a
covering give a degree-zero-or-two subgraph on the trivalent square–octagon
graph. At either endpoint of a dummy edge, let a,b be the two surviving
ring-edge loop bits and d the dummy loop bit. Even parity requires

`d = a XOR b`.

The four surviving edges at the contracted vertex have even total parity,
so the other original endpoint gives the same d. Every even subgraph on
the contracted square graph therefore has one and only one lift back to
the trivalent graph. Relative to the fixed reference, that lift determines
one covering. The cycle-space dimensions agree: both graphs have 2N+1
independent loop bits on a connected torus, including all winding sectors.

## What to draw at a four-leg vertex

If four surviving edges are occupied, d=0: each original endpoint joins its
own pair of ring edges. Draw those as two separate arcs at the contracted
vertex. Their pairing is fixed by which original endpoint they belonged to;
it is not an additional free choice and need not be a crossing. This fixed
routing preserves the individual loops, including loops surrounding a small
square. If only two edges are occupied, their connection is unique, whether
or not the lifted path traverses the dummy edge.

A binary edge picture without visible routing can appear to have touching
loops. The contracted graph plus its fixed endpoint pairing supplies the
information needed to read it correctly.

## Implementation and checks

In `rvb/square_octagon.py`, `dummy_contracted_square_geometry(nx,ny)` returns
the loop graph and the inherited endpoint pairs. Vertex `2*c+axis` is the
east/north dummy midpoint (axis 0/1). Surviving edge `4*c+k` follows the small
square from corner k to k+1, using k=E,N,W,S. These IDs match the first 4N
entries returned by `base_loop_pattern`.

`dummy_contracted_loop_pattern` implements contraction;
`lift_dummy_contracted_loop_pattern` restores erased loop bits;
`bits_from_dummy_contracted_loops` recovers the covering; and
`dummy_contracted_routing` supplies the fixed edge pairings for a drawing.
The older `square_loop_pattern` remains available and explicitly documents
that it implements the different small-square-collapse projection.

The 17 geometry tests pass. New tests independently enumerate **all** even
edge subsets through 2×2 cells, verify the one-to-one covering map and its
inverse for a nontrivial reference, and check that occupied small-square
loops survive four-leg contraction without being spuriously joined.
