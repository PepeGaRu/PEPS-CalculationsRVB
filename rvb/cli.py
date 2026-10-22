"""One-command, reproducible numerical benchmarks and readable figures."""

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import time

import numpy as np

from .exact import cluster_benchmark
from .lattices import square, kagome
from .orthogonal import square_cylinder_benchmark, square_plane_benchmark, kagome_orthogonal_benchmark


def _run(label, function):
    print(label, flush=True)
    start = time.perf_counter()
    result = function()
    result["wall_seconds"] = time.perf_counter() - start
    print(f"  completed in {result['wall_seconds']:.2f} s", flush=True)
    return result


def make_figures(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    plane = data["square_orthogonal_plane"]
    r = np.asarray(plane["distance"])[1:]
    c = np.abs(plane["connected_correlation"])[1:]
    axes[0, 0].loglog(r, c, label="infinite plane: Kasteleyn quadrature", lw=2)
    axes[0, 0].loglog(r, 1/(np.pi*r)**2, "--", label="asymptotic 1/(π²r²)", color="black", alpha=.65)
    axes[0, 0].set(title="Square · orthogonal dimers (collinear bonds)", xlabel="separation r (lattice spacings)", ylabel="|connected occupation correlation|")
    axes[0, 0].legend(fontsize=8)
    for result in data["singlet_cylinders"]:
        ax = axes[0, 1] if result["lattice"] == "square" else axes[1, 0]
        r = np.asarray(result["distance"])
        c = np.abs(result["connected_correlation"])
        keep = c > max(1e-13, c.max(initial=0)*1e-10)
        plot = ax.loglog if result["lattice"] == "square" else ax.semilogy
        sector = result.get("metadata", {}).get("virtual_vacuum_parity")
        label = f"width {result['width']}" + (f", {sector} sector" if sector else "")
        plot(r[keep], c[keep], "o-", ms=3, label=label)
    axes[0, 1].set(title="Square · singlet RVB (finite cylinders)", xlabel="axial separation r (lattice spacings)", ylabel="|connected spin-bond correlation|")
    axes[1, 0].set(title="Kagome · singlet RVB (finite cylinders)", xlabel="axial separation r (unit cells)", ylabel="|connected spin-bond correlation|")
    for ax in (axes[0, 1], axes[1, 0]):
        ax.legend(fontsize=8)
    orth = data["kagome_orthogonal"]
    axes[1, 1].plot(orth["distance"], orth["connected_correlation"], "o-", color="tab:green")
    axes[1, 1].axhline(0, color="gray", lw=.7)
    axes[1, 1].set(title="Kagome · orthogonal dimers (exact parity solution)", xlabel="axial separation r (unit cells)", ylabel="connected occupation correlation")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Nearest-neighbor RVB and orthogonal dimers\nFinite-cylinder singlet results are not thermodynamic exponent estimates", fontsize=14)
    fig.savefig(output / "correlations.png", dpi=170)
    fig.savefig(output / "correlations.svg")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, lattice in zip(axes, ("square", "kagome")):
        results = [r for r in data["singlet_cylinders"] if r["lattice"] == lattice and r.get("correlation_length") is not None]
        if lattice == "kagome":
            for sector in ("even", "odd"):
                selected = [r for r in results if r["metadata"].get("virtual_vacuum_parity") == sector]
                ax.plot([r["width"] for r in selected], [r["correlation_length"] for r in selected], "o-", label=f"singlet RVB, {sector} sector")
        else:
            ax.plot([r["width"] for r in results], [r["correlation_length"] for r in results], "o-", label="singlet RVB: coupled transfer mode")
        if lattice == "square":
            results = data["square_orthogonal_cylinders"]
            ax.plot([r["width"] for r in results], [r["correlation_length"] for r in results], "s--", label="orthogonal: transfer eigenvalue")
        ax.set(title=lattice.capitalize(), xlabel="circumference (square sites / kagome unit cells)", ylabel="axial decay length (unit cells)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Finite-cylinder lengths: units and fit diagnostics are in data.json")
    fig.savefig(output / "correlation-lengths.png", dpi=170)
    plt.close(fig)


def write_report(data, output):
    lines = ["# RVB correlation benchmark results", "", f"Generated: {data['generated_utc']}", "",
             "These results concern equal-amplitude nearest-neighbor spin-1/2 RVB states and uniform orthogonal dimers, not optimized Heisenberg ground states.", "",
             "## What was calculated", "",
             "- Orthogonal square dimers: numerical quadrature of the exact infinite-plane Kasteleyn solution, plus exact finite-cylinder transfer matrices.",
             "- Orthogonal kagome dimers: exact arrow/parity probabilities on a periodic torus; no sampling noise.",
             "- Singlet RVB: exact sparse PEPS transfer matrices on infinite cylinders of specified finite circumference. No boundary bond-dimension truncation is used.",
             "- Small clusters: independent sums over singlet coverings are compared with quimb PEPS contractions.", "",
             "## Numerical results", "",
             f"The square orthogonal-dimer even-distance fit gives exponent **{data['square_orthogonal_plane']['fitted_even_distance_exponent']:.6f}**, approaching the exact exponent **2** (fit range {data['square_orthogonal_plane']['fit_distance_range']}). This is not the singlet-RVB exponent.", "",
             "| Singlet lattice | Width | Virtual parity | Mean spin bond | Dimer axial length | Eigenvector residual |", "|---|---:|---|---:|---:|---:|"]
    for r in data["singlet_cylinders"]:
        xi = r.get("correlation_length")
        xi_text = f"{xi:.6g}" if xi is not None else "not resolved"
        parity = r.get("metadata", {}).get("virtual_vacuum_parity") or "zero staggered flux"
        lines.append(f"| {r['lattice']} | {r['width']} | {parity} | {r['mean_bond_energy']:.9f} | {xi_text} | {r['eigen_residual']:.2g} |")
    lines += ["", "Lengths count axial unit-cell translations. A kagome unit-cell translation is twice the nearest-neighbor spacing. Width counts square sites or three-site kagome cells. These units must be matched before comparing to a paper.", "",
              "The lengths come from transfer eigenmodes with nonzero coupling to the physical bond observable; spectral weights, solver residuals and separate tail-fit diagnostics are in `data.json`. Circumference scaling is required for bulk statements. In particular, a finite square-cylinder length does not mean the square plane is gapped.", "",
              "The paper audit found that unrestricted kagome fixed points choose odd virtual-vacuum parity at width 2 and even parity at width 4. The revised benchmarks select and compare both parities separately. No vison line is inserted. These labels specify the implemented boundary basis; they do not assert that every paper cylinder boundary convention is identical.", "",
              "Both square cylinder models measure circumferential (vertical) bonds separated along the axis. The infinite-plane orthogonal reference instead uses collinear horizontal bonds; it is explicitly labeled separately. Bond orientation can change amplitudes, oscillations, and cancellations of leading terms.", "",
              "The kagome orthogonal-dimer calculation gives exact zero connected correlation beyond its local support. Exact zeros are not measurements of a small positive exponential decay length.", "",
              "## Validation", "", "| Cluster | Spins | Coverings | PEPS / covering fidelity | Total spin squared |", "|---|---:|---:|---:|---:|"]
    for r in data["clusters"]:
        lines.append(f"| {r['lattice']} | {r['sites']} | {r['coverings']} | {r['peps_fidelity_to_covering_sum']:.12f} | {r['total_spin_squared']:.2g} |")
    lines += ["", "## Interpretation and remaining limits", "",
              "RVB uses physical spin-bond observables B=Si·Sj (equivalently connected singlet-projector correlations, since P=1/4−B). Orthogonal dimers use occupation n=0 or 1. These are different observables in different Hilbert spaces; their absolute amplitudes should not be equated.", "",
              "An algebraic dimer correlator has no finite asymptotic exponential correlation length. Spin correlations and dimer correlations must not be conflated. On square RVB the critical behavior is in the dimer correlations.", "",
              "A full reproduction of a paper's singlet-RVB bulk exponent and quoted lengths additionally needs its exact geometry, sector, observable, units and systematic circumference or CTMRG bond-dimension extrapolation. The finite-cylinder data here are a validated starting calculation, not that extrapolation.", "",
              "See `../upstream` (when generated) for independent peps-torch CTMRG results and `../../docs/libraries-and-references.md` for source and library review.", "",
              "![Correlation curves](correlations.png)", "", "![Cylinder correlation lengths](correlation-lengths.png)", ""]
    (output / "report.md").write_text("\n".join(lines))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Calculate and compare square/kagome RVB and orthogonal dimer correlations.")
    parser.add_argument("--profile", choices=("quick", "reference"), default="quick", help="quick: one singlet width per lattice; reference: two widths")
    parser.add_argument("--output", type=Path, default=Path("results/quick"))
    parser.add_argument("--distance", type=int, default=32)
    args = parser.parse_args(argv)
    if args.distance < 8:
        parser.error("--distance must be at least 8 for useful decay checks")
    from .peps import rvb_cylinder_benchmark
    args.output.mkdir(parents=True, exist_ok=True)
    data = {"generated_utc": datetime.now(timezone.utc).isoformat(), "profile": args.profile,
            "versions": {name: version(name) for name in ("numpy", "scipy", "quimb", "matplotlib")},
            "clusters": [], "singlet_cylinders": [], "square_orthogonal_cylinders": []}
    def checkpoint():
        (args.output / "data.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    for lattice in (square(), kagome(2,2,True,True)):
        data["clusters"].append(_run(f"Exact {lattice.name} cluster check", lambda g=lattice: cluster_benchmark(g)))
        checkpoint()
    data["square_orthogonal_plane"] = _run("Square orthogonal dimers: infinite plane", lambda: square_plane_benchmark(128))
    for width in (4, 6, 8, 10, 12):
        data["square_orthogonal_cylinders"].append(_run(f"Square orthogonal cylinder W={width}", lambda w=width: square_cylinder_benchmark(w, args.distance, orientation="vertical")))
    data["kagome_orthogonal"] = _run("Kagome orthogonal dimers: exact 8×8 torus", lambda: kagome_orthogonal_benchmark(8, 8))
    checkpoint()
    for lattice, widths in (("square", (4,6) if args.profile == "reference" else (4,)), ("kagome", (2,4) if args.profile == "reference" else (2,))):
        for width in widths:
            for sector in (("even", "odd") if lattice == "kagome" else ("dominant",)):
                data["singlet_cylinders"].append(_run(f"{lattice.capitalize()} singlet RVB cylinder W={width}, sector={sector}", lambda l=lattice,w=width,s=sector: rvb_cylinder_benchmark(l, w, args.distance, sector=s)))
                checkpoint()
    make_figures(data, args.output)
    write_report(data, args.output)
    print(f"Read {args.output / 'report.md'}; raw data and figures are in the same directory.", flush=True)


if __name__ == "__main__":
    main()
