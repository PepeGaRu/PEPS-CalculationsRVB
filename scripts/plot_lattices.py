"""Draw the finite lattice examples and the singlet orientation convention."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rvb.lattices import square, kagome


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, graph in zip(axes, (square(4, 4), kagome(3, 3))):
        for a, b in graph.edges:
            xa, ya = graph.positions[a]
            xb, yb = graph.positions[b]
            ax.plot([xa, xb], [ya, yb], color="0.65", lw=1)
            ax.annotate("", xy=(xa*.35+xb*.65, ya*.35+yb*.65),
                        xytext=(xa*.65+xb*.35, ya*.65+yb*.35),
                        arrowprops={"arrowstyle": "->", "color": "#3264a8", "lw": 1.1})
        x, y = zip(*graph.positions)
        ax.scatter(x, y, s=24, c="#183c64", zorder=3)
        ax.set_aspect("equal")
        ax.set_axis_off()
        convention = "A → B singlet arrows" if graph.name == "square" else "Cyclic singlet arrows on both triangle types"
        ax.set_title(graph.name.capitalize()+" lattice\n"+convention)
    fig.suptitle("Lattice geometry and fixed singlet sign conventions\nEach covering pairs every site once; the RVB state sums over coverings", fontsize=13)
    destination = Path(__file__).resolve().parents[1] / "results/lattices.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=170)
    print(destination)


if __name__ == "__main__":
    main()
