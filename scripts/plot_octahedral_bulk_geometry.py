#!/usr/bin/env python3
"""Explain finite transverse tubes versus the three-dimensional bulk target.

This figure is a geometry schematic, not a correlation measurement. Periodic
wrapping bonds are omitted for clarity. Physical spins lie on cubic-bond
midpoints; the cubic vertices are the centres of corner-sharing octahedra.
"""
from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#1769a6"
RED = "#b52b33"
GRAY = "#89939e"


def draw_window(ax, shape, box_aspect, observed):
    centres = np.array(list(product(*(range(n) for n in shape))), float)
    spins = []
    for centre in centres:
        for axis in range(3):
            end = centre.copy()
            end[axis] += 1
            if end[axis] >= shape[axis]:
                continue
            ax.plot(*np.stack([centre, end]).T, color=GRAY, lw=.7, alpha=.4)
            spins.append((centre+end)/2)
    spins = np.array(spins)
    ax.scatter(*centres.T, s=9, color="#6a747e", depthshade=False)
    ax.scatter(*spins.T, s=15, color=BLUE, depthshade=False)
    pair = np.array(observed, float)
    ax.plot(*pair.T, color=RED, lw=1.9, linestyle="--", zorder=10)
    ax.scatter(*pair.T, color=RED, edgecolors="white", s=68,
               linewidths=.9, depthshade=False, zorder=11)
    for point, label in zip(pair, ("$i$", "$j$")):
        ax.text(*(point+np.array([0, -.1, .15])), label, color=RED,
                fontsize=12, weight="bold")
    for axis, n in enumerate(shape):
        getattr(ax, f"set_{'xyz'[axis]}lim")(-.3, n-.6)
        getattr(ax, f"set_{'xyz'[axis]}ticks")([])
        getattr(ax, f"set_{'xyz'[axis]}label")("xyz"[axis], labelpad=-10)
        getattr(ax, f"{'xyz'[axis]}axis").pane.fill = False
    ax.grid(False)
    ax.set_box_aspect(box_aspect)
    ax.view_init(elev=22, azim=-57)


def draw_octahedron(ax):
    points = np.vstack([np.eye(3)*.5, -np.eye(3)*.5])
    for i in range(6):
        for j in range(i):
            if not np.isclose(np.dot(points[i], points[j]), 0):
                continue
            ax.plot(*points[[i, j]].T, color=GRAY, alpha=.6, lw=1)
    for a, colour in enumerate((RED, "#4a8b49", "#8f62ae")):
        ax.plot(*points[[a, a+3]].T, color=colour, lw=2.3, alpha=.85)
    ax.scatter(*points.T, color=BLUE, s=75, depthshade=False,
               edgecolors="white", linewidths=.7)
    ax.scatter([0], [0], [0], color="#666666", marker="+", s=70)
    ax.text(.04, .04, .08, "centre\n(no spin)", fontsize=9, color="#4b535b")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=21, azim=-56)
    ax.set_axis_off()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT/"results/octahedral-spin1-bulk")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "pdf.fonttype": 42})
    fig = plt.figure(figsize=(12.4, 8.4), facecolor="white")
    fig.suptitle("Your lattice is three-dimensional. The previous calculation used a narrow tube.",
                 x=.5, y=.969, fontsize=15, weight="bold")
    left = fig.add_axes([.035, .405, .45, .475], projection="3d")
    right = fig.add_axes([.53, .405, .435, .475], projection="3d")
    draw_window(left, (7, 2, 2), (6, 1.5, 1.5),
                [(.5, 0, 0), (3.5, 0, 0)])
    draw_window(right, (4, 4, 4), (1, 1, 1),
                [(.5, 0, 0), (2.5, 0, 0)])
    left.set_title("Previous approximation: finite transverse size", fontsize=12, pad=8)
    right.set_title("Requested target: the bulk 3D lattice", fontsize=12, pad=8)
    left.text2D(.5, -.005,
                "$L_x\\to\\infty$, but $L_y=2$, $L_z=2$\n"
                "Periodic boundaries in $y,z$; a finite window is shown.",
                transform=left.transAxes, ha="center", va="top", fontsize=10)
    right.text2D(.5, -.005,
                 "$L_x,L_y,L_z\\to\\infty$ together\n"
                 "Increase boxes in all three directions; test convergence.",
                 transform=right.transAxes, ha="center", va="top", fontsize=10)
    # Arrows drawn in figure coordinates express continuation of the tube,
    # rather than falsely adding physical bonds at the cut window boundaries.
    left.text2D(.035, .52, "$\\cdots$", transform=left.transAxes,
                fontsize=23, color="#606971")
    left.text2D(.93, .46, "$\\cdots$", transform=left.transAxes,
                fontsize=23, color="#606971")
    legend_handles = [
        Line2D([], [], marker="o", linestyle="", color=BLUE, markersize=6,
               label="Physical spin 1 at a cubic-bond midpoint"),
        Line2D([], [], marker="o", linestyle="", color="#6a747e", markersize=3,
               label="Octahedron centre: no physical spin"),
        Line2D([], [], color=GRAY, lw=1,
               label="Cubic grid: geometric guide, not a singlet bond"),
    ]
    fig.legend(handles=legend_handles, loc="center", bbox_to_anchor=(.5, .325),
               ncol=1, frameon=False, fontsize=10)
    local = fig.add_axes([.045, .04, .23, .24], projection="3d")
    draw_octahedron(local)
    local.set_title("One octahedron: six physical spins", fontsize=10, pad=-2)
    fig.text(.29, .235,
             "Local geometry\n"
             "Coloured diameters show one allowed three-dimer support.\n"
             "A trimer instead occupies one triangular face. Adjacent\n"
             "octahedra share a midpoint spin, which is covered once.",
             fontsize=10, linespacing=1.5, va="top")
    fig.text(.29, .126,
             "The red pair measures a physical spin correlation along $x$:\n"
             "$C(r)=\\langle\\mathbf{S}_i\\!\\cdot\\!\\mathbf{S}_{i+r\\hat{x}}\\rangle$. "
             "The surrounding $y,z$ directions still affect it.",
             fontsize=11, linespacing=1.65, va="top")
    fig.text(.5, .025,
             "Geometry schematic only. Periodic wrapping bonds are omitted. "
             "No converged bulk correlation length is claimed by this figure.",
             fontsize=9, color="#565e66", ha="center")
    for extension in ("png", "pdf"):
        path = args.output/f"geometry.{extension}"
        fig.savefig(path, dpi=190, bbox_inches="tight")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
