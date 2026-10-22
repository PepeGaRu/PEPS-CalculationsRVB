#!/usr/bin/env python3
"""Draw the physical spin graph and the user's dummy-link contraction.

Schematic distances are adjustable; spin IDs, bonds, coverings and loops are
read from the model. These are scientific diagrams, not generated artwork.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from collections import defaultdict

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.square_octagon import SquareOctagonOrthogonalDimers

BLUE = "#246bb1"
ORANGE = "#df811c"
GREEN = "#25835d"
RED = "#c6374b"
GRAY = "#9ca6af"
DIR = np.array(((1, 0), (0, 1), (-1, 0), (0, -1)))
OFF = np.array(((.14, .14), (-.14, .14), (-.14, -.14), (.14, -.14),
                (.36, 0), (0, .36), (-.36, 0), (0, -.36)))


def main():
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon, PathPatch
    from matplotlib.path import Path as MplPath

    plt.rcParams.update({"font.size": 12, "axes.titlepad": 14,
                         "savefig.facecolor": "white", "legend.frameon": False})
    out = ROOT / "results/octagon-square"
    out.mkdir(parents=True, exist_ok=True)

    def line(ax, p, q, **kwargs):
        ax.plot(*np.array((p, q)).T, **kwargs)

    def save(fig, name):
        fig.savefig(out / (name + ".png"), dpi=155, bbox_inches="tight")
        fig.savefig(out / (name + ".svg"), bbox_inches="tight")
        fig.savefig(out / (name + ".pdf"), bbox_inches="tight")
        plt.close(fig)

    def cell(ax, center, triangles=True, labels=False):
        center = np.array(center)
        corners = center + .28*DIR
        for k in range(4):
            line(ax, corners[k], corners[(k+1) % 4], color=GRAY, lw=1.7, zorder=3)
            if triangles:
                tri = center + OFF[[k, (k-1) % 4, 4+k]]
                ax.add_patch(Polygon(tri, facecolor="#e8f3eb", edgecolor=GREEN,
                                     lw=1.1, linestyle="--", zorder=2))
        ax.scatter(*corners.T, s=20, facecolors="white", edgecolors=GRAY,
                   linewidth=1, zorder=4)
        points = center + OFF
        ax.scatter(*points[:4].T, color=BLUE, s=83, zorder=5,
                   edgecolor="white", linewidth=.8)
        ax.scatter(*points[4:].T, color=ORANGE, s=83, zorder=5,
                   edgecolor="white", linewidth=.8)
        if labels:
            names = ["0", "1", "2", "3", "$d_E$", "$d_N$", "$d_W$", "$d_S$"]
            for i, (p, name) in enumerate(zip(points, names)):
                delta = OFF[i] / np.linalg.norm(OFF[i]) * .065
                ax.text(*(p+delta), name, color=BLUE if i < 4 else ORANGE,
                        ha="center", va="center", fontsize=14, zorder=6)

    def decorate(ax):
        ax.set_aspect("equal")
        ax.set_axis_off()

    fig, ax = plt.subplots(figsize=(8.4, 8.0), layout="constrained")
    for x in range(2):
        for y in range(2):
            c = np.array((x, y))
            for k, direction in enumerate(DIR):
                # Outgoing links continue to neighboring cells outside this patch.
                endpoint = c + (.72 if 0 <= x+direction[0] < 2 and
                                0 <= y+direction[1] < 2 else .48)*direction
                line(ax, c+.28*direction, endpoint, color=GRAY, lw=2,
                     linestyle=(0, (4, 3)), zorder=3)
            cell(ax, c)
            ax.text(x, y, "square", color="#68747f", ha="center", va="center", fontsize=10)
    ax.text(.5, .5, "OCTAGON", color="#68747f", ha="center", va="center", fontsize=16)
    ax.annotate("Two spins on this dummy link", xy=(.5, 0), xytext=(.5, -.32),
                ha="center", fontsize=11, arrowprops={"arrowstyle": "->", "color": ORANGE})
    ax.set(xlim=(-.52, 1.52), ylim=(-.55, 1.5),
           title="The lattice used in the calculation\nEvery filled dot is one spin ½")
    decorate(ax)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=BLUE, lw=0, markersize=9, label="One spin on each small-square edge"),
        Line2D([], [], marker="o", color=ORANGE, lw=0, markersize=9, label="Two endpoint spins on each dummy link"),
        Line2D([], [], color=GREEN, linestyle="--", label="Allowed bonds within each triangle"),
        Line2D([], [], color=GRAY, lw=2, label="Underlying square–octagon graph"),
        Line2D([], [], marker="o", color=GRAY, markerfacecolor="white", lw=0,
               markersize=5, label="Hollow nodes: graph vertices, not extra spins")],
        loc="upper center", bbox_to_anchor=(.5, -.015), fontsize=11)
    save(fig, "spin-lattice")

    fig, ax = plt.subplots(figsize=(7.6, 7.4), layout="constrained")
    for k, direction in enumerate(DIR):
        line(ax, .28*direction, .76*direction, color=GRAY, lw=2, linestyle="--", zorder=3)
        partner = .64*direction
        ax.scatter(*partner, color=ORANGE, alpha=.45, s=83, zorder=4)
    cell(ax, (0, 0), labels=True)
    ax.annotate("One triangle:\nspins 0, 3, $d_E$", xy=(.22, 0), xytext=(.46, -.36),
                fontsize=12, ha="center", arrowprops={"arrowstyle": "->", "color": GREEN})
    ax.annotate("The second spin\non the same dummy link", xy=(-.64, 0), xytext=(-.51, -.38),
                fontsize=11, ha="center", arrowprops={"arrowstyle": "->", "color": ORANGE})
    ax.text(0, .64, "", ha="center")
    ax.set(xlim=(-.82, .82), ylim=(-.72, .8),
           title="Zoom: one simulation cell\n4 blue spins + 4 orange endpoint spins = 8 spins")
    ax.text(.5, -.015, "Pale orange dots belong to neighboring cells.\n"
            "Green dashed lines show allowed bonds, not a chosen covering.",
            transform=ax.transAxes, ha="center", fontsize=11, va="top")
    decorate(ax)
    save(fig, "spin-unit-cell")

    # Show both allowed states of one physical dummy link.
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.1), layout="constrained")
    for index, ax in enumerate(axes):
        left = np.array(((-.45, .4), (-.45, -.4), (0, 0)))
        right = np.array(((1.45, .4), (1.45, -.4), (1, 0)))
        for tri in (left, right):
            ax.add_patch(Polygon(tri, facecolor="#e8f3eb", edgecolor=GRAY, lw=1))
            ax.scatter(*tri[:2].T, color=BLUE, s=90, zorder=5)
            ax.scatter(*tri[2], color=ORANGE, s=90, zorder=5)
        line(ax, (0, 0), (1, 0), color=GRAY, lw=1, linestyle="--")
        if index == 0:
            for p, q in (((0, 0), (1, 0)), (left[0], left[1]), (right[0], right[1])):
                line(ax, p, q, color=GREEN, lw=4)
        else:
            for tri, sign in ((left, -1), (right, 1)):
                line(ax, tri[0], tri[2], color=GREEN, lw=4)
                line(ax, tri[1], tri[1]+np.array((.2*sign, -.18)), color=GREEN, lw=3,
                     linestyle=":")
        ax.set_title("Unused by triangles: the orange spins pair with each other" if index == 0 else
                     "Both used by triangles: each orange spin pairs with a blue spin", fontsize=12)
        ax.set(xlim=(-.85, 1.85), ylim=(-.65, .55))
        decorate(ax)
    axes[1].text(.5, -.07, "Dotted green continuations: the remaining blue spins pair in adjacent triangles.",
                 transform=axes[1].transAxes, ha="center", fontsize=10)
    save(fig, "dummy-link-states")

    model = SquareOctagonOrthogonalDimers(3, 3)
    example = json.loads((out / "mapping-example.json").read_text())
    bits = model.bits_from_covering(example["covering_edge_indices"])
    reference = model.bits_from_covering(example["reference_covering_edge_indices"])
    loops = model.base_loop_pattern(bits, reference)
    contracted = model.dummy_contracted_loop_pattern(bits, reference)
    assert model.bits_from_dummy_contracted_loops(contracted, reference) == bits
    covering = set(model.covering_from_bits(bits))
    positions = np.array([np.array((x, y)) + OFF[k]
                          for x, y, k in model.lattice.labels])
    coords = np.array([(x, y) for x in range(3) for y in range(3)])

    def draw_cover(ax):
        for index, (i, j) in enumerate(model.edges):
            p, q = positions[[i, j]]
            kwargs = {"color": GREEN if index in covering else "#d2d8dc",
                      "lw": 3 if index in covering else .7}
            if np.max(np.abs(p-q)) > 1.5:
                delta = q-p
                direction = int(np.argmax(abs(delta)))
                for point, sign in ((p, 1), (q, -1)):
                    end = point.copy(); end[direction] += .18*sign
                    line(ax, point, end, **kwargs)
            else:
                line(ax, p, q, **kwargs)
        for ordinary, color in ((True, BLUE), (False, ORANGE)):
            keep = np.array([i % 8 < 4 for i in range(len(positions))]) == ordinary
            ax.scatter(*positions[keep].T, s=23, color=color, zorder=4)
        ax.set(title="1. One allowed singlet covering", xlim=(-.6, 2.6), ylim=(-.6, 2.6))
        decorate(ax)

    def draw_base(ax):
        for c, center in enumerate(coords):
            corners = center+.28*DIR
            for k in range(4):
                active = loops[4*c+k]
                line(ax, corners[k], corners[(k+1) % 4],
                     color=RED if active else "#d2d8dc", lw=3 if active else 1)
            for axis, direction in enumerate(DIR[:2]):
                if max(center+direction) > 2:
                    continue
                active = loops[4*model.ncells+2*c+axis]
                line(ax, center+.28*direction, center+.72*direction,
                     color=RED if active else GRAY, lw=3 if active else 1,
                     linestyle="-" if active else "--")
        ax.set(title="2. Loops on the square–octagon graph", xlim=(-.6, 2.6), ylim=(-.6, 2.6))
        decorate(ax)

    def square_coords(p):
        return np.array((p[0]+p[1]-.5, p[1]-p[0]+.5))

    def draw_contracted(ax):
        # Draw a periodic patch unwrapped at its boundary. No active dummy
        # link crosses this example's boundary; local closed loops stay intact.
        segments, incidence = {}, defaultdict(list)
        for c, center in enumerate(coords):
            corners = [square_coords(center+.5*d) for d in DIR]
            for k in range(4):
                edge = 4*c+k
                p, q = corners[k], corners[(k+1) % 4]
                segments[edge] = (p, q)
                line(ax, p, q, color="#d2d8dc", lw=1)
                if contracted[edge]:
                    incidence[tuple(p)].append((edge, q, k, c))
                    incidence[tuple(q)].append((edge, p, (k+1) % 4, c))
        for edge, (p, q) in segments.items():
            if not contracted[edge]:
                continue
            u = (q-p)/np.linalg.norm(q-p)
            p0 = p+.14*u if len(incidence[tuple(p)]) == 4 else p
            q0 = q-.14*u if len(incidence[tuple(q)]) == 4 else q
            line(ax, p0, q0, color=RED, lw=3)
        for point, entries in incidence.items():
            assert len(entries) in (2, 4)
            if len(entries) != 4:
                continue
            # At degree four, each original endpoint pairs its own two ring
            # edges. These arcs keep that inherited routing visibly separate.
            groups = defaultdict(list)
            for edge, other, corner, c in entries:
                groups[(c, corner)].append(other)
            p = np.array(point)
            assert sorted(map(len, groups.values())) == [2, 2]
            for endpoints in groups.values():
                q, r = endpoints
                a = p+.14*(q-p)/np.linalg.norm(q-p)
                b = p+.14*(r-p)/np.linalg.norm(r-p)
                ax.add_patch(PathPatch(MplPath([a, p, b], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]),
                                       facecolor="none", edgecolor=RED, lw=3))
        ax.set(title="3. Shrink dummy links → square-lattice loops")
        decorate(ax)

    drawers = [draw_cover, draw_base, draw_contracted]
    names = ["mapping-singlets", "mapping-octagons", "mapping-square"]
    for draw, name in zip(drawers, names):
        fig, ax = plt.subplots(figsize=(7.4, 6.8), layout="constrained")
        draw(ax)
        caption = ("Green segments are singlets; every colored dot is one spin ½."
                   if name == "mapping-singlets" else
                   "Red edges show changes relative to the stated reference covering."
                   if name == "mapping-octagons" else
                   "Dummy links become vertices; small-square edges survive.\n"
                   "Turned 45° for display. Small arcs retain the inherited loop connections.")
        ax.text(.5, -.015, caption, transform=ax.transAxes, ha="center", va="top", fontsize=10)
        save(fig, name)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout="constrained")
    for draw, ax in zip(drawers, axes):
        draw(ax)
        ax.set_title(ax.get_title(), fontsize=12)
    fig.suptitle("Your mapping: covering → square–octagon loops → contract dummy links", fontsize=15)
    fig.supxlabel("The last step is one-to-one, with the inherited loop routing. It does not collapse the small squares.", fontsize=11)
    save(fig, "dummy-contraction")
    (out / "dummy-contraction-example.json").write_text(json.dumps({
        "source_covering": "mapping-example.json",
        "map": "Contract dummy links; retain ring edges, with inherited routing",
        "reference_description": example["reference"],
        "base_loop_bits": loops, "contracted_ring_loop_bits": contracted,
        "reconstructed_bits_equal_original": True,
        "routing_pairs_by_contracted_vertex": model.dummy_contracted_routing(contracted),
        "physical_spin_graph_unchanged": True,
        "coordinates_are_schematic": True}, indent=2) + "\n")
    print("Wrote spin lattice, cell zoom, dummy states, and corrected contraction figures.")


if __name__ == "__main__":
    main()
