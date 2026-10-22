# Diagram and local-map audit for the collaborator note

The source report has nine pages. Its page 6 is headed **“Corrected mapping:
contract the dummy links”** and uses `dummy-contraction.pdf`, produced by
`scripts/plot_octagon_geometry.py`. The older `mapping.pdf` instead collapses
whole small squares and must not be used as the intended map.

## Which square must be upright

The relevant square is the **grey base-graph square** whose four sides carry
the four blue physical spins. Its vertices sit inside the four decorating
triangles and are not physical spins. The smaller polygon joining the blue spins consists of
triangle bases and should not be labelled “square” in the new figure.

Apply the same clockwise 45-degree rotation to every coordinate in the
first two mapping panels:

\[
R(x,y)=\frac1{\sqrt2}(x+y,\;y-x).
\]

In the existing script the grey corners are `0.28*[(1,0),(0,1),(-1,0),(0,-1)]`.
They become the corners of an upright square after R. The blue positions
`(+.14,+.14),(-.14,+.14),(-.14,-.14),(+.14,-.14)` become the midpoints of
its right, top, left, and bottom sides. Rotate singlets, grey edges, orange
slots and triangle outlines together. Do not rotate only the legend or
grey overlay. The final contracted-square panel already uses this rotation,
with an additional uniform scale and translation, so the three panels can
share one consistent orientation.

## The dummy-link table has eight variants

At the two endpoints of one dummy link, label the ordinary slot occupations
`(a,b)` and `(c,d)`; a value 1 means that spin is paired inside the displayed
triangle. Let `t=1` mean both orange endpoint spins are used by their
respective triangles. At each endpoint,

\[
a\oplus b=t=c\oplus d.
\]

There are **four variants in each of the two dummy classes**, not just one
representative per class and not four external variants plus one internal
variant:

| Dummy class | Ordinary slot patterns | Local singlets |
|---|---|---|
| `t=0`: the orange spins form their internal singlet | `0000`, `1100`, `0011`, `1111` | Each triangle independently is empty or contains its ordinary–ordinary singlet. |
| `t=1`: both orange spins pair into their triangles | `1010`, `1001`, `0110`, `0101` | Each orange endpoint independently pairs with one of its two ordinary spins. |

An ordinary spin with local occupation 0 is **not free**: it pairs in its
other adjacent triangle. Show that by a short outgoing singlet stub, or
explain it explicitly in the caption. The existing two-row dummy figure
shows one representative from each class; its internal representative has
both ordinary–ordinary bonds and therefore omits the three other internal
possibilities.

Exact affine counting on a 3x3 torus confirms every one of these eight local
patterns has a global completion and probability 1/8 in the **uniform
orthogonal-covering ensemble**. Do not assign those probabilities to
physical nonorthogonal singlets.

## Reference and square-loop icons

Loops are occupation changes relative to one fixed complete reference
covering, not the singlet bonds themselves. The simple local reference
`t0=0, (a0,b0,c0,d0)=0000` can be extended globally: put an internal singlet
on every dummy link, ordinary–ordinary singlets in the E/W triangles around
each base square, and leave N/S triangles empty. Each vertical dummy link
then has the stated local reference. A horizontal-looking local table is
obtained by rotating that representative.

For this local reference the four surviving loop bits equal `a,b,c,d`.
Place the square-vertex ports as

```
        a=N
         |
 b=W — vertex — c=E
         |
        d=S
```

The inherited endpoint groups are `(N,W)` and `(E,S)`; the **cyclic** order
is `a,c,d,b`, not `a,b,c,d`. The eight loop icons are consequently:

| Bits `abcd` | Contracted square-loop vertex |
|---|---|
| 0000 | Empty |
| 1100 | NW turn |
| 0011 | ES turn |
| 1111 | Two separate arcs, NW and ES |
| 1010 | NE turn |
| 1001 | NS straight |
| 0110 | WE straight |
| 0101 | WS turn |

For a general reference the displayed loop bits must instead be XORed with
its ordinary occupations, and the erased dummy loop bit is `t XOR t0`.
At either original endpoint it is uniquely recovered as the XOR of the
two surviving loop bits. Thus shrinking every dummy edge to a vertex is a
bijection of even-edge configurations. A four-leg vertex retains the fixed
two-arc routing above; there is no additional free crossing choice.

“Cutting/erasing the dummy edge” here means **contracting its endpoints in
the loop graph**, not deleting a physical singlet from the quantum state.
The map preserves the covering labels but does not itself preserve or
compute quantum singlet overlaps.

For N base squares on a connected periodic torus, the contracted graph has
2N vertices, 4N edges, and 2N+1 independent cycle bits. It is distinct from
the older N-vertex graph produced by collapsing each whole base square.
