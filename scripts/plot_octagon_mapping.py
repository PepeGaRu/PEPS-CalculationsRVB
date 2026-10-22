#!/usr/bin/env python3
"""Draw an actual valid covering and its two loop representations."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.square_octagon import SquareOctagonOrthogonalDimers, square_octagon_orthogonal_benchmark


def configuration(model, active_links):
    """Choose one explicit local resolution for a specified even square graph."""
    bits = 0
    for c, axis in active_links:
        bits |= 1 << (4 * model.ncells + 2 * c + axis)
    for x in range(model.nx):
        for y in range(model.ny):
            c = x * model.ny + y
            z = [1]
            for k in range(1, 4):
                t = (bits >> model.dummy_variable(x, y, k)) & 1
                z.append(1 ^ t ^ z[-1])
            bits |= sum(value << (4 * c + k) for k, value in enumerate(z))
    model.covering_from_bits(bits)
    return bits


def main():
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    out = ROOT / "results/octagon-square"
    out.mkdir(parents=True, exist_ok=True)
    model = SquareOctagonOrthogonalDimers(3, 3)
    links = set()
    # Two plaquettes meeting at one vertex expose the coarse degree-four case.
    for x, y in ((0, 0), (1, 1)):
        links.update(((x * 3 + y, 0), ((x + 1) * 3 + y, 1),
                      (x * 3 + y + 1, 0), (x * 3 + y, 1)))
    bits = configuration(model, links)
    reference = configuration(model, set())
    covering = set(model.covering_from_bits(bits))
    base_loops = model.base_loop_pattern(bits, reference)
    coarse = model.square_loop_pattern(bits)
    positions = np.array(model.lattice.positions)
    directions = np.array(((1, 0), (0, 1), (-1, 0), (0, -1)))
    base_positions = np.array([(x + .34 * dx, y + .34 * dy)
                              for x in range(3) for y in range(3)
                              for dx, dy in directions])
    base_edges = []
    for c in range(9):
        base_edges.extend((4*c+k, 4*c+(k+1) % 4) for k in range(4))
    for x in range(3):
        for y in range(3):
            c = x*3+y
            base_edges.extend(((4*c, 4*(((x+1) % 3)*3+y)+2),
                               (4*c+1, 4*(x*3+(y+1) % 3)+3)))
    degree = [0] * 36
    for active, (a, b) in zip(base_loops, base_edges):
        if active:
            degree[a] += 1; degree[b] += 1
    assert all(d in (0, 2) for d in degree)
    coarse_degree = [0] * 9
    for c in range(9):
        x, y = divmod(c, 3)
        for axis, target in ((0, ((x+1) % 3)*3+y), (1, x*3+(y+1) % 3)):
            if coarse[2*c+axis]:
                coarse_degree[c] += 1; coarse_degree[target] += 1
    assert all(d % 2 == 0 for d in coarse_degree)

    fig, axes = plt.subplots(1, 3, figsize=(13, 5), layout="constrained")
    for index, (a, b) in enumerate(model.edges):
        p, q = positions[[a, b]]
        color = "#18865b" if index in covering else "#c9cdd1"
        width = 2.5 if index in covering else .65
        if np.max(np.abs(q-p)) > 1.5:
            # Periodic boundary singlets are drawn as short outgoing stubs.
            delta = q-p
            axis = int(np.argmax(np.abs(delta)))
            for point, sign in ((p, 1), (q, -1)):
                end = point.copy(); end[axis] += sign * .25
                axes[0].plot(*np.array([point, end]).T, color=color, lw=width)
        else:
            axes[0].plot(*np.array([p, q]).T, color=color, lw=width)
    ordinary = np.array([i % 8 < 4 for i in range(72)])
    axes[0].scatter(*positions[ordinary].T, s=18, color="#3169b2", zorder=3)
    axes[0].scatter(*positions[~ordinary].T, s=18, color="#dc8324", zorder=3)
    axes[0].legend(handles=[Line2D([], [], marker="o", lw=0, color="#3169b2", label="One spin on each square edge"),
                            Line2D([], [], marker="o", lw=0, color="#dc8324", label="Two spins on each dummy link"),
                            Line2D([], [], color="#18865b", lw=2.5, label="One allowed singlet covering")],
                   loc="upper center", bbox_to_anchor=(.5, -.01), fontsize=9)
    for active, (a, b) in zip(base_loops, base_edges):
        p, q = base_positions[[a, b]]
        if np.max(np.abs(q-p)) > 1.5:
            continue
        axes[1].plot(*np.array([p, q]).T, color="#c73c48" if active else "#c9cdd1",
                     lw=3 if active else 1)
    axes[1].scatter(*base_positions.T, s=10, color="#45515c", zorder=3)
    for x in range(3):
        for y in range(3):
            for axis, delta in ((0, (1, 0)), (1, (0, 1))):
                xx, yy = x + delta[0], y + delta[1]
                if xx >= 3 or yy >= 3:
                    continue
                active = coarse[2*(x*3+y)+axis]
                axes[2].plot([x, xx], [y, yy], color="#c73c48" if active else "#c9cdd1",
                             lw=3 if active else 1)
    axes[2].scatter(*np.array([(x, y) for x in range(3) for y in range(3)]).T,
                    s=22, color="#45515c", zorder=3)
    axes[2].scatter([1], [1], s=180, facecolor="white", edgecolor="#c73c48", lw=2, zorder=4)
    axes[2].annotate("Four incident loop edges:\nretain the local resolution", xy=(1, 1),
                     xytext=(.1, -.35), fontsize=10,
                     arrowprops={"arrowstyle": "->", "color": "#45515c"})
    titles = ("1. Your decorated spin lattice", "2. Closed loops on square–octagon", "3. Projection onto the square lattice")
    for ax, title in zip(axes, titles):
        ax.set(title=title, xlim=(-.55, 2.55), ylim=(-.55, 2.55), aspect="equal")
        ax.set_axis_off()
    axes[1].text(.5, -.04, "Flips relative to one stated reference covering\nEvery vertex has loop degree 0 or 2", transform=axes[1].transAxes,
                 ha="center", va="top", fontsize=10)
    fig.suptitle("A checked covering-to-loop map; all physical spins remain paired", fontsize=14)
    fig.savefig(out / "mapping.png", dpi=160, bbox_inches="tight")
    fig.savefig(out / "mapping.pdf", bbox_inches="tight")
    plt.close(fig)
    (out / "mapping-example.json").write_text(json.dumps({
        "geometry": "3x3 periodic cells; boundary singlets drawn as stubs",
        "covering_edge_indices": sorted(covering), "covering_edges": [model.edges[i] for i in sorted(covering)],
        "reference_covering_edge_indices": model.covering_from_bits(reference),
        "reference": "All dummy links internally paired; E/W ring-ring singlets in every cell",
        "base_loop_edges": [list(edge) for active, edge in zip(base_loops, base_edges) if active],
        "base_loop_degrees": degree, "coarse_square_edge_bits": coarse,
        "coarse_square_degrees": coarse_degree,
        "coarse_projection_multiplicity": "2^N local resolutions for each square even-edge pattern",
        "coordinates_are_schematic": True}, indent=2) + "\n")
    exact = square_octagon_orthogonal_benchmark(16, 4)
    (out / "orthogonal.json").write_text(json.dumps(exact, indent=2) + "\n")
    fig, ax = plt.subplots(figsize=(7, 4.3), layout="constrained")
    for family in exact["families"]:
        ax.plot(exact["distance"], family["connected_correlation"], "o-", label=f"{family['bond_kind']} occupation")
    ax.axhline(0, color="gray", lw=.7)
    ax.text(.98, .8, "Exact rational counting\nNo sampling or exponential fit\nDistant occupation correlations are zero",
            ha="right", transform=ax.transAxes, fontsize=10)
    ax.set(title="Orthogonal dimers: exactly finite-range occupation correlations",
           xlabel="Separation along x (original square cells)", ylabel="Connected occupation correlation")
    ax.legend(); ax.grid(alpha=.2)
    fig.savefig(out / "orthogonal-correlations.png", dpi=160, bbox_inches="tight")
    fig.savefig(out / "orthogonal-correlations.pdf", bbox_inches="tight")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
