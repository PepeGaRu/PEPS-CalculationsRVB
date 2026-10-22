"""Exact top projection and display opening of an octochlore centre layer.

Two distinct vertical apex sites project onto each cubic centre. Selecting a
checkerboard family of centre-lattice plaquettes groups the horizontal ports
into two complementary face frames. Opening their common projected apex
produces the triangle-decorated square--octagon incidence graph. Green frames
are geometry: they do NOT denote simultaneous occupied trimers.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collaborator-notes"
BLUE, ORANGE = "#256da4", "#d78324"
GREY, GREEN = "#9ca3aa", "#458b68"
# Four checkerboard plaquettes surrounding one unselected centre-grid square.
PLAQUETTES = np.array(((-.5, .5), (.5, -.5), (1.5, .5), (.5, 1.5)))
CORNERS = .5 * np.array(((-1, -1), (1, -1), (1, 1), (-1, 1)))
DELTA = 1 / (2 * (1 + np.sqrt(2)))


def polyline(ax, points, **kwargs):
    points = np.asarray(points)
    ax.plot(*points.T, **kwargs)


def closed(points):
    return np.concatenate((points, points[:1]))


def green_frame(ax, points, zorder=3):
    polyline(ax, closed(np.asarray(points)), color=GREEN, lw=1.05,
             ls=(0, (3, 2)), zorder=zorder)


def projected_layer(ax, delta):
    all_v = set()
    blue_points = []
    for w in PLAQUETTES:
        vertices = w + CORNERS
        q = w + (1 - 2 * delta) * CORNERS
        blue = (q + np.roll(q, -1, axis=0)) / 2
        ax.add_patch(Polygon(q, facecolor="#f1f2f3", edgecolor=GREY,
                             lw=1.45, zorder=1))
        blue_points.extend(blue)
        for i, v in enumerate(vertices):
            all_v.add(tuple(v))
            s = np.sign(w - v)
            apex = v + .5 * delta * s
            green_frame(ax, (blue[(i - 1) % 4], blue[i], apex))

    for vt in sorted(all_v):
        v = np.asarray(vt)
        # A selected plaquette above v is NW on even sites, NE on odd sites.
        sx = -1 if round(v.sum()) % 2 == 0 else 1
        s = np.array((sx, 1))
        if delta:
            q_plus, q_minus = v + delta * s, v - delta * s
            polyline(ax, (q_plus, q_minus), color=GREY, lw=1.45, zorder=1)
            # Both slots are shown even where the partner face is outside patch.
            for sign in (-1, 1):
                p = v + sign * .5 * delta * s
                ax.scatter(*p, s=12, c=ORANGE, zorder=6, edgecolors="white", linewidths=.3)
        else:
            # Concentric marks are two physical sites at exactly the same xy.
            ax.scatter(*v, s=26, facecolors="white", edgecolors=ORANGE,
                       linewidths=1.1, zorder=6)
            ax.scatter(*v, s=5.5, c=ORANGE, zorder=7)
    blue_points = np.asarray(blue_points)
    ax.scatter(*blue_points.T, s=12, c=BLUE, zorder=5, edgecolors="white", linewidths=.3)
    ax.set(xlim=(-1.38, 2.38), ylim=(-1.38, 2.38), aspect="equal")
    ax.axis("off")


def octahedron(ax):
    # A transparent affine projection makes both vertical apices unambiguous.
    xyz = np.array(((.5, 0, 0), (0, .5, 0), (-.5, 0, 0),
                    (0, -.5, 0), (0, 0, .5), (0, 0, -.5)))
    projection = np.array(((1., .16), (.48, .52), (0., 1.)))
    p = xyz @ projection
    faces = ((1, 2, 4), (0, 3, 5))
    for f in faces:
        ax.add_patch(Polygon(p[list(f)], facecolor=GREEN, alpha=.055, edgecolor="none"))
    for i in range(4):
        for j in ((i + 1) % 4, 4, 5):
            polyline(ax, (p[i], p[j]), color=GREY, lw=.75, alpha=.85, zorder=1)
    polyline(ax, (p[4], p[5]), color=ORANGE, lw=1.0, ls=":", alpha=.8)
    for f in faces:
        green_frame(ax, p[list(f)])
    ax.scatter(*p[:4].T, s=25, c=BLUE, edgecolors="white", linewidths=.5, zorder=5)
    ax.scatter(*p[4:].T, s=30, c=ORANGE, edgecolors="white", linewidths=.5, zorder=5)
    for i, label, offset in ((4, r"$z_+$", (.06, .025)),
                              (5, r"$z_-$", (.06, -.035)),
                              (0, r"$x_+$", (.035, -.01)),
                              (2, r"$x_-$", (-.15, -.045)),
                              (1, r"$y_+$", (.025, .035)),
                              (3, r"$y_-$", (-.14, -.055))):
        ax.text(*(p[i] + offset), label, fontsize=8.3)
    ax.set(xlim=(-.71, .71), ylim=(-.66, .70), aspect="equal")
    ax.axis("off")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(7.15, 2.6))
    ax0 = fig.add_axes((.10, .18, .28, .67))
    ax1 = fig.add_axes((.58, .16, .32, .72))
    octahedron(ax0)
    projected_layer(ax1, 0.)
    for x, title in ((.24, "(a) One octahedron"),
                     (.74, "(b) Exact layer top view")):
        fig.text(x, .97, title, ha="center", va="top", fontsize=8.8, weight="bold")
    fig.text(.24, .87, "Two complementary face frames", ha="center", fontsize=7.4)
    fig.text(.74, .87, r"Project along $z$: $z_+$ and $z_-$ coincide", ha="center", fontsize=7.5)
    fig.text(.49, .50, r"$\longrightarrow$", ha="center", fontsize=22, color="#555c65")
    fig.text(.49, .66, r"View along $z$", ha="center", fontsize=8, color="#555c65")
    fig.text(.24, .135, "Frames do not mean occupied trimers", ha="center", fontsize=7.4)
    fig.text(.74, .135, "Concentric orange marks = two sites", ha="center", fontsize=7.5)
    fig.text(.5, .052,
             "Blue: equatorial spins     Orange: upper/lower apex spins     Green dashed: selected triangle frames",
             ha="center", fontsize=7.4)
    fig.savefig(OUT / "layer-relation.pdf", bbox_inches="tight", pad_inches=.03)
    fig.savefig(OUT / "layer-relation.png", dpi=240, bbox_inches="tight", pad_inches=.03)
    plt.close(fig)
    metadata = {
        "figure_panels": ["(a) One octahedron", "(b) Exact layer top view"],
        "opening_display": "Discussed in the audit text; not drawn as a figure panel.",
        "centre": "v=(m,n,0)",
        "horizontal_sites": "v +/- e_x/2, v +/- e_y/2",
        "apex_sites": "v +/- e_z/2; both project to xy=v",
        "selected_plaquettes": "lower-left corner parity m+n odd",
        "port_pairs_even_vertex": ["NW,z+", "SE,z-"],
        "port_pairs_odd_vertex": ["NE,z+", "SW,z-"],
        "opened_grey_corner": "q(v,w)=v+delta*sign(w-v)",
        "opened_orange_site": "d(v,w)=v+delta/2*sign(w-v)",
        "opened_blue_site": "midpoint of consecutive q corners of selected square",
        "delta_regular": float(DELTA),
        "zero_delta_limit": "exact projected coordinates for every physical site",
        "scope": "Geometric incidence only; face frames are not simultaneous trimers."
    }
    (OUT / "layer-projection-coordinates.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
