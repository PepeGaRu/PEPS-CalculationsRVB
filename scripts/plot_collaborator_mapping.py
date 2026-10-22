"""Compact, upright scientific diagrams for the collaborator note.

All positions are constructed from the model. The eight local rows exhaust
the allowed support patterns of two triangles joined by one dummy link.
The global mapping example is the checked covering in the previous note.
"""
from collections import defaultdict
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, PathPatch
from matplotlib.path import Path as MplPath
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.square_octagon import SquareOctagonOrthogonalDimers

OUT = ROOT / "results/collaborator-notes"
BLUE, ORANGE, GREEN, RED = "#246bb1", "#df811c", "#25835d", "#c6374b"
GRAY, PALE = "#86949e", "#dce2e6"
DIR = np.array(((1, 0), (0, 1), (-1, 0), (0, -1)))
OFF = np.array(((.14, .14), (-.14, .14), (-.14, -.14), (.14, -.14),
                (.36, 0), (0, .36), (-.36, 0), (0, -.36)))
ROT = np.array(((1, 1), (-1, 1))) / np.sqrt(2)


def upright(points):
    return np.asarray(points) @ ROT.T


def line(ax, p, q, **kwargs):
    ax.plot(*np.array((p, q)).T, **kwargs)


def clean(ax):
    ax.set_aspect("equal")
    ax.set_axis_off()


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=.03)
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight", pad_inches=.03)
    plt.close(fig)


def geometry():
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.12),
                             gridspec_kw={"width_ratios": (1.3, 1)}, layout="constrained")
    for ax, cells in zip(axes, ([(0, 0), (1, 0), (0, 1), (1, 1)], [(0, 0)])):
        for center in map(np.array, cells):
            corners = upright(center + .28*DIR)
            ax.add_patch(Polygon(corners, facecolor="#f2f4f5", edgecolor=GRAY,
                                 linewidth=1.3, zorder=1))
            for k, direction in enumerate(DIR):
                p, q = upright((center+.28*direction, center+.56*direction))
                line(ax, p, q, color=GRAY, lw=1.05, ls=(0, (3, 2)), zorder=2)
                tri = upright(center + OFF[[k, (k-1) % 4, k+4]])
                ax.add_patch(Polygon(tri, fc="#e9f4ed", ec=GREEN, lw=.65,
                                     linestyle="--", zorder=2))
            # The intended square lies inside the triangle decorations;
            # draw its outline above their fill, so it remains unambiguous.
            ax.add_patch(Polygon(corners, fc="none", ec=GRAY, lw=1.4, zorder=3))
            points = upright(center+OFF)
            for selected, color in ((slice(0, 4), BLUE), (slice(4, 8), ORANGE)):
                ax.scatter(*points[selected].T, c=color, s=14, edgecolor="white",
                           linewidth=.25, zorder=5)
        clean(ax)
    axes[0].set_title("Relevant gray squares are upright", fontsize=9, pad=3)
    axes[1].set_title("One square: edge spins and corner triangles", fontsize=9, pad=3)
    axes[1].text(0, 0, r"$P$", color=GRAY, ha="center", va="center", fontsize=11)
    axes[1].annotate("spin ½ on an edge", xy=upright(OFF[0]), xytext=(.37, .20),
                     fontsize=7.2, ha="center", color=BLUE,
                     arrowprops={"arrowstyle": "-", "lw": .6, "color": BLUE})
    axes[1].annotate("one dummy endpoint", xy=upright(OFF[7]), xytext=(-.36, -.43),
                     fontsize=7.2, ha="center", color=ORANGE,
                     arrowprops={"arrowstyle": "-", "lw": .6, "color": ORANGE})
    axes[1].set_xlim(-.58, .58)
    axes[1].set_ylim(-.5, .5)
    save(fig, "upright-geometry")


def local_cover(ax, bits, t, center=(0, .32), scale=.96):
    # This is the actual right corner of cell (0,0) and left corner of
    # cell (1,0), translated to their midpoint and rotated exactly as in
    # geometry().  The gray square sides therefore remain horizontal /
    # vertical, and their four arms contract to a=N,b=W,c=E,d=S.
    ordinary = upright(np.array((OFF[0], OFF[3], OFF[1]+(1, 0),
                                 OFF[2]+(1, 0))) - (.5, 0))
    dummy = upright(np.array((OFF[4], OFF[6]+(1, 0))) - (.5, 0))
    corners = upright(np.array(((.28, 0), (.72, 0))) - (.5, 0))
    port_directions = np.array(((0, 1), (-1, 0), (1, 0), (0, -1)))
    for i, direction in enumerate(port_directions):
        arm = ordinary[i]-corners[i//2]
        assert np.allclose(arm, np.linalg.norm(arm)*direction)
    cv = lambda p: np.asarray(p)*scale + np.array(center)
    for end, inds in enumerate(((0, 1), (2, 3))):
        tri = cv(np.vstack((ordinary[list(inds)], dummy[end])))
        ax.add_patch(Polygon(tri, fc="#e9f4ed", ec=GREEN, lw=.65,
                             linestyle="--", zorder=1))
        for i in inds:
            # The visible L is part of the relevant upright gray square,
            # not an edge of a decorating triangle.
            tip = ordinary[i]+.075*port_directions[i]
            line(ax, cv(corners[end]), cv(tip), color=GRAY, lw=1.2, zorder=2)
        if t:
            chosen = next(i for i in inds if bits[i])
            line(ax, cv(dummy[end]), cv(ordinary[chosen]), color=GREEN, lw=2, zorder=3)
        elif bits[inds[0]]:
            line(ax, cv(ordinary[inds[0]]), cv(ordinary[inds[1]]), color=GREEN, lw=2, zorder=3)
    line(ax, *cv(corners), color=GRAY, lw=.75, ls=(0, (3, 2)), zorder=2)
    line(ax, *cv(dummy), color=GREEN if not t else GRAY,
         lw=2 if not t else .7, ls="-" if not t else ":", zorder=3)
    for i, p in enumerate(ordinary):
        if not bits[i]:
            delta = .12*port_directions[i]
            line(ax, cv(p), cv(p+delta), color=GREEN, lw=1.3, ls=":", zorder=3)
        label_point = cv(p+.10*port_directions[i]
                         + .075*np.array((-port_directions[i, 1], port_directions[i, 0])))
        ax.text(*label_point, "abcd"[i], fontsize=6.8, ha="center", va="center", color="#4c5964")
    ax.scatter(*cv(ordinary).T, color=BLUE, s=13, edgecolor="white", linewidth=.3, zorder=5)
    ax.scatter(*cv(dummy).T, color=ORANGE, s=13, edgecolor="white", linewidth=.3, zorder=5)


def loop_icon(ax, bits, center=(.02, -.43), scale=.22):
    # Ports are a=N,b=W,c=E,d=S; this is not their cyclic order.
    origin = np.array(center)
    ends = origin + scale*np.array(((0, 1), (-1, 0), (1, 0), (0, -1)))
    for p in ends:
        line(ax, origin, p, color=PALE, lw=1.2, zorder=1)
    active = np.flatnonzero(bits)
    if len(active) == 2:
        p, q = ends[active]
        line(ax, p, origin, color=RED, lw=1.7)
        line(ax, origin, q, color=RED, lw=1.7)
    elif len(active) == 4:
        for pair in ((0, 1), (2, 3)):
            p, q = ends[list(pair)]
            u, v = origin+.25*(p-origin), origin+.25*(q-origin)
            line(ax, p, u, color=RED, lw=1.7)
            line(ax, q, v, color=RED, lw=1.7)
            ax.add_patch(PathPatch(MplPath([u, origin, v],
                              [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]),
                              fc="none", ec=RED, lw=1.7))
    else:
        ax.scatter(*origin, s=8, c=GRAY)
    ax.text(-.43, -.43, "→", fontsize=12, ha="center", va="center", color=GRAY)
    ax.text(.43, -.46, "".join(map(str, bits)), fontsize=7.5, ha="center", va="center")


def covering_table():
    patterns = [((0, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 1, 1)),
                ((1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 0), (0, 1, 0, 1))]
    assert len({bits for row in patterns for bits in row}) == 8
    for t, row in enumerate(patterns):
        for a, b, c, d in row:
            assert (a+b+t) % 2 == 0 and (c+d+t) % 2 == 0
    fig = plt.figure(figsize=(6.7, 2.5), layout="constrained")
    grid = fig.add_gridspec(2, 5, width_ratios=(.95, 1, 1, 1, 1), hspace=.06, wspace=.02)
    for t, row in enumerate(patterns):
        label = fig.add_subplot(grid[t, 0])
        label.text(.02, .62, "INTERNAL" if t == 0 else "EXTERNAL", fontsize=8, weight="bold")
        label.text(.02, .34, "dummy singlet\n$ t=0 $" if t == 0 else "both endpoints used\n$ t=1 $",
                   fontsize=7.5, linespacing=1.6)
        label.set_axis_off()
        for col, bits in enumerate(row, 1):
            ax = fig.add_subplot(grid[t, col])
            local_cover(ax, bits, t)
            loop_icon(ax, bits)
            ax.set(xlim=(-.69, .69), ylim=(-.76, .75))
            clean(ax)
    save(fig, "dummy-covering-table")
    (OUT / "dummy-covering-table.json").write_text(json.dumps({
        "ordinary_port_order": "a=N,b=W,c=E,d=S after contraction",
        "local_geometry": "Two actual adjacent corner triangles, rotated as in the lattice; gray L-shaped sides are upright square edges",
        "reference": "0000 on a representative link; dummy endpoints internally paired",
        "patterns_by_dummy_occupation": patterns,
        "locally_unused_ordinary_spins": "Covered in their other neighboring triangles",
        "erased_dummy_loop_bit": "a XOR b = c XOR d",
        "all_eight_even_vertices_present": True,
    }, indent=2)+"\n")


def mapping():
    model = SquareOctagonOrthogonalDimers(3, 3)
    data = json.loads((ROOT / "results/octagon-square/mapping-example.json").read_text())
    bits = model.bits_from_covering(data["covering_edge_indices"])
    ref = model.bits_from_covering(data["reference_covering_edge_indices"])
    covering = set(model.covering_from_bits(bits))
    loops = model.base_loop_pattern(bits, ref)
    contracted = model.dummy_contracted_loop_pattern(bits, ref)
    assert model.bits_from_dummy_contracted_loops(contracted, ref) == bits
    coords = np.array([(x, y) for x in range(3) for y in range(3)])
    original_positions = np.array([np.array((x, y))+OFF[k] for x, y, k in model.lattice.labels])
    fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.35), layout="constrained")
    for c, center in enumerate(coords):
        corners = upright(center+.28*DIR)
        for ax in axes[:2]:
            ax.add_patch(Polygon(corners, fc="#f3f5f6", ec=GRAY, lw=.75, zorder=1))
        for k in range(4):
            tri = upright(center + OFF[[k, (k-1) % 4, k+4]])
            axes[0].add_patch(Polygon(tri, fc="none", ec=GREEN, lw=.65,
                                      linestyle="--", zorder=2))
            if loops[4*c+k]:
                line(axes[1], corners[k], corners[(k+1) % 4], color=RED, lw=1.8, zorder=4)
        for axis, direction in enumerate(DIR[:2]):
            if max(center+direction) > 2:
                continue
            p, q = upright((center+.28*direction, center+.72*direction))
            active = loops[4*model.ncells+2*c+axis]
            line(axes[0], p, q, color=GRAY, lw=.6, ls=":", zorder=1)
            line(axes[1], p, q, color=RED if active else GRAY,
                 lw=1.8 if active else .65, ls="-" if active else ":", zorder=3)
    for edge, (i, j) in enumerate(model.edges):
        if edge not in covering:
            continue
        p, q = original_positions[[i, j]]
        if np.max(np.abs(p-q)) > 1.5:
            direction = int(np.argmax(abs(q-p)))
            for point, sign in ((p, 1), (q, -1)):
                end = point.copy()
                end[direction] += .14*sign
                line(axes[0], *upright((point, end)), color=GREEN, lw=1.45, zorder=4)
        else:
            line(axes[0], *upright((p, q)), color=GREEN, lw=1.45, zorder=4)
    positions = upright(original_positions)
    for condition, color in ((True, BLUE), (False, ORANGE)):
        selected = np.array([i % 8 < 4 for i in range(len(positions))]) == condition
        axes[0].scatter(*positions[selected].T, s=5.3, color=color, zorder=5)
    incidence, segments = defaultdict(list), {}
    # sqrt(2)*rotation maps the contracted square to a unit grid.
    square = lambda p: np.asarray(p) @ (ROT*np.sqrt(2)).T
    for c, center in enumerate(coords):
        corners = square(center+.5*DIR)
        for k in range(4):
            edge = 4*c+k
            p, q = corners[k], corners[(k+1) % 4]
            segments[edge] = (p, q)
            line(axes[2], p, q, color=PALE, lw=.65)
            if contracted[edge]:
                incidence[tuple(p)].append((q, c, k))
                incidence[tuple(q)].append((p, c, (k+1) % 4))
    for edge, (p, q) in segments.items():
        if not contracted[edge]:
            continue
        u = (q-p)/np.linalg.norm(q-p)
        line(axes[2], p+.13*u if len(incidence[tuple(p)]) == 4 else p,
             q-.13*u if len(incidence[tuple(q)]) == 4 else q, color=RED, lw=1.8)
    for point, entries in incidence.items():
        assert len(entries) in (2, 4)
        if len(entries) != 4:
            continue
        groups = defaultdict(list)
        for other, c, corner in entries:
            groups[c, corner].append(other)
        p = np.array(point)
        for q, r in groups.values():
            a = p+.13*(q-p)/np.linalg.norm(q-p)
            b = p+.13*(r-p)/np.linalg.norm(r-p)
            axes[2].add_patch(PathPatch(MplPath([a, p, b],
                             [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]),
                             fc="none", ec=RED, lw=1.8))
    for ax, title in zip(axes, ("1  Singlet covering", "2  Reference-relative loops",
                                "3  Contract dummy edges")):
        clean(ax)
        ax.set_title(title, fontsize=8.5, pad=3)
        ax.margins(.08)
    save(fig, "corrected-mapping-upright")


def main():
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                         "savefig.facecolor": "white"})
    geometry()
    covering_table()
    mapping()
    print(OUT)


if __name__ == "__main__":
    main()
