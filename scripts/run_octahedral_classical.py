#!/usr/bin/env python3
"""Exact equal-weight benchmark for the supplied octahedral covering notes."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.octahedral import LOCAL_STATES, OctahedralCoverings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[3, 4, 6, 8, 10, 12])
    parser.add_argument("--output", type=Path, default=ROOT / "results/octahedral-classical")
    args = parser.parse_args()
    if any(l < 3 for l in args.sizes):
        parser.error("Use sizes >=3 here; period-two clusters are validation cases, not bulk correlation benchmarks.")
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    records = []
    largest = None
    for length in sorted(set(args.sizes)):
        begin = perf_counter()
        model = OctahedralCoverings(length, length, length)
        setup_seconds = perf_counter()-begin
        cube_begin = perf_counter()
        cube_rank = model.cube_move_rank
        cube_seconds = perf_counter()-cube_begin
        assert model.dimension == model.ncells+3*length-1
        assert cube_rank == model.ncells-3*length+2
        records.append({"length": length, "octahedra": model.ncells,
                        "physical_sites": model.nsites, "constraint_rank": model.rank,
                        "log2_covering_count": model.dimension,
                        "cube_move_rank": cube_rank,
                        "log2_cube_move_sectors": model.dimension-cube_rank,
                        "constraint_setup_seconds": setup_seconds,
                        "cube_rank_seconds": cube_seconds})
        largest = model
        print(f"L={length}: {model.nsites} sites, 2^{model.dimension} coverings; "
              f"linear solve {setup_seconds:.4f}s, cube rank {cube_seconds:.4f}s", flush=True)
    model = largest
    distance = list(range(model.nx//2+1))
    origin = (0, 0, 0)
    signs = (1, 1, 1)
    start = perf_counter()
    local = [{"state": s.name, "support": s.incoming,
              "probability_exact": str(model.state_probability(origin, s)),
              "probability": float(model.state_probability(origin, s))}
             for s in LOCAL_STATES]
    means = {"Dx": model.dimer_probability(origin, 0),
             "Dy": model.dimer_probability(origin, 1),
             "Dz": model.dimer_probability(origin, 2),
             "any_trimer": model.trimer_probability(origin),
             "face_+++": model.face_probability(origin, signs)}
    series = {
        "Dx_along_x": [model.connected_dimer(origin, 0, (r, 0, 0), 0) for r in distance],
        "Dx_along_y": [model.connected_dimer(origin, 0, (0, r, 0), 0) for r in distance],
        "any_trimer_along_x": [model.connected_trimer(origin, (r, 0, 0)) for r in distance],
        "face_+++_along_x": [model.connected_face(origin, signs, (r, 0, 0), signs) for r in distance],
    }
    # Full catalog joint probabilities retain substantially more information
    # than testing only one dimer or summing the eight trimer faces.
    joint = {}
    for name, target in (("neighbor", (1, 0, 0)), ("separated", (1, 1, 0))):
        joint[name] = [[str(model.joint_state_probability(origin, s, target, t))
                        for t in LOCAL_STATES] for s in LOCAL_STATES]
    measurement_seconds = perf_counter()-start
    sample_model = OctahedralCoverings(4, 4, 4)
    sample = sample_model.sample(1, seed=20261020)[0]
    objects = sample_model.covering(sample)
    assert sum(len(obj.sites) for obj in objects) == sample_model.nsites
    sample_record = {"shape": sample_model.shape, "seed": 20261020,
                     "ownership_bits_hex": hex(sample),
                     "objects": [asdict(obj) for obj in objects],
                     "site_positions": sample_model.site_positions,
                     "local_states": [sample_model.local_state(sample, c).name
                                      for c in range(sample_model.ncells)]}
    (out / "sample-covering.json").write_text(json.dumps(sample_record, indent=2)+"\n")
    result = {
        "model": "Uniform classical owner-labelled octahedral dimer/trimer coverings",
        "assumption": "Every complete allowed covering has equal statistical weight; no fugacities or quantum overlaps.",
        "local_rule": "Any subset of three diameter dimers OR one of eight face trimers; two complementary trimers excluded.",
        "geometry": "Periodic cubic centres; physical covering sites at cubic-edge midpoints",
        "distance_unit": "Distance between neighboring octahedron centres",
        "sectors": "All allowed sectors included in uniform counting",
        "method": "Exact GF(2) linear algebra and rational event probabilities",
        "timings_scope": "Measured in this running environment; exclude tests, plotting and PDF generation",
        "sizes": records, "local_states": local,
        "means_exact": {k: str(v) for k, v in means.items()},
        "distance": distance,
        "correlations_exact": {k: [str(v) for v in vs] for k, vs in series.items()},
        "correlations": {k: [float(v) for v in vs] for k, vs in series.items()},
        "joint_catalog_probabilities_exact": joint,
        "measurement_seconds": measurement_seconds,
        "correlation_conclusion": "Exactly finite support for single-octahedron observables beyond nearest-neighbour centres, for all periodic dimensions >=3.",
        "exponential_correlation_length": None,
        "line_parity_constraint": "XOR of trimer indicators along a periodic coordinate line equals its length mod 2",
        "cube_sector_warning": "Classical move sectors are not by themselves proof of a quantum fracton phase.",
        "source_sha256": {"rvb/octahedral.py": hashlib.sha256((ROOT/"rvb/octahedral.py").read_bytes()).hexdigest(),
                          "scripts/run_octahedral_classical.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    (out / "data.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")

    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    labels = [r"$D_x$", r"$D_y$", r"$D_z$", "Any\ntrimer", "One specified\nface trimer"]
    colors = ["#276aa4"]*3+["#824c9d"]*2
    axes[0].bar(labels, [float(v) for v in means.values()], color=colors, width=.65)
    for index, value in enumerate(means.values()):
        axes[0].text(index, float(value)+.015, str(value), ha="center", fontsize=12)
    axes[0].set(title="Exact object probabilities", ylabel="Probability", ylim=(0, .6))
    keys = ["Dx_along_x", "any_trimer_along_x", "face_+++_along_x"]
    for key, label, color, marker in zip(keys, [r"Diameter dimer $D_x$", "Any trimer", "Specified face (+++)"],
                                       ["#276aa4", "#824c9d", "#c7832b"], ["o", "s", "^"]):
        axes[1].plot(distance, [float(v) for v in series[key]], marker=marker,
                     color=color, label=label, markersize=5, lw=1.4)
    axes[1].axhline(0, color="gray", lw=.7)
    axes[1].set(title="Connected occupation correlations along x", xlabel="Separation in cubic centre spacings",
                ylabel=r"$\langle O_0 O_r\rangle-\langle O_0\rangle\langle O_r\rangle$",
                xticks=distance, ylim=(-.085, .285))
    axes[1].legend(fontsize=9)
    axes[1].text(.42, .39, "Exact zeros beyond nearest neighbours\nNo exponential or power-law fit", fontsize=10,
                 transform=axes[1].transAxes)
    for ax in axes:
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Your 16-state octahedral model: equal-weight classical coverings", fontsize=14)
    fig.savefig(out/"exact-correlations.png", dpi=170, bbox_inches="tight")
    fig.savefig(out/"exact-correlations.pdf", bbox_inches="tight")

    overview = plt.figure(figsize=(8.27, 11.69))
    overview.text(.09, .93, "Octahedral dimer/trimer coverings", fontsize=21, color="#276aa4", weight="bold")
    overview.text(.09, .885, "First exact classical calculation", fontsize=17)
    paragraphs = [
        (.82, "Model used", "Corner-sharing octahedra with cubic centres. Covering sites lie at\n"
         "cubic-bond midpoints. Allow any subset of the three diameter dimers\n"
         "or one face trimer; forbid two complementary trimers.\n"
         "Every complete allowed covering is assigned equal statistical weight."),
        (.66, "Feasibility demonstrated", f"Computed periodic sizes {', '.join(map(str, sorted(set(args.sizes))))} cubed.\n"
         f"Largest: {model.ncells:,} octahedra and {model.nsites:,} covering sites.\n"
         f"The binary constraint solve took {records[-1]['constraint_setup_seconds']:.3f} seconds.\n"
         f"It represents 2^{model.dimension} coverings without enumerating them."),
        (.50, "Exact local results", "Each diameter dimer has probability 1/4. Any trimer: 1/2.\n"
         "Each specified triangular-face trimer: 1/16.\n"
         "Connected correlations of single-octahedron observables vanish\n"
         "beyond neighboring centres. There is no nonzero decay tail to fit."),
        (.34, "Nonlocal structure remains", "For N=Lx Ly Lz and S=Lx+Ly+Lz, the covering count is 2^(N+S-1).\n"
         "Cube-only moves leave 2^(2S-3) disconnected classical sectors.\n"
         "The trimer parity on every periodic coordinate line is fixed by\n"
         "the line length. Zero two-point tails do not remove these constraints."),
        (.18, "Scope", "This is a classical covering calculation. No quantum dimer/trimer\n"
         "state or Hamiltonian was assumed. Nonuniform weights would\n"
         "require new calculations. The cited ten-state Ising model has\n"
         "different constraints and its correlation results do not transfer.")]
    for y, title, body in paragraphs:
        overview.text(.09, y, title, fontsize=13, weight="bold")
        overview.text(.09, y-.028, body, fontsize=10.6, va="top", linespacing=1.55)
    overview.text(.09, .035, "Reproduce: python scripts/run_octahedral_classical.py\n"
                  "Raw results and a checked sample covering accompany this PDF.", fontsize=9, color="#68747f")
    pdf_path = out/"Octahedral_Classical_First_Results.pdf"
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(overview)
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(overview); plt.close(fig)

    rows = ["| L | Octahedra | Covering sites | log2(number of coverings) | Binary solve (s) |",
            "|---:|---:|---:|---:|---:|"]
    rows += [f"| {v['length']} | {v['octahedra']} | {v['physical_sites']} | {v['log2_covering_count']} | {v['constraint_setup_seconds']:.4f} |"
             for v in records]
    report = """# First exact calculation for the supplied octahedral coverings

The sixteen-state rule is computationally tractable here. This benchmark uses **equal statistical weight for every complete allowed covering**, summed over all sectors. No quantum spin state, singlet overlap, object fugacity, or Hamiltonian is assumed.

The geometry and local objects follow the supplied two-page notes: cubic octahedron centres, shared covering sites at cubic bond midpoints, any subset of three opposite-vertex dimers or one of eight face trimers. Two complementary face trimers on one octahedron are excluded. Crossings of dimer chords at an octahedron centre do not introduce additional sites.

## Actual computations

""" + "\n".join(rows) + """

Times measure construction and solution of the binary constraint system in this running environment. Computing the cube-move rank, observables, validation tests and plots has additional cost recorded separately in `data.json`. The enormous covering counts are handled as linear spaces; these configurations were not enumerated on the large tori.

## Exact probabilities and correlations

Every local catalog state has probability 1/16. Thus each diameter dimer has probability 1/4, any trimer has probability 1/2, and each specified face trimer has probability 1/16. Mean object counts per octahedron are 3/4 dimers and 1/2 trimers. The mean occupied-site count is therefore `2*(3/4)+3*(1/2)=3`, consistent with 3N globally covered sites.

![Exact local probabilities and connected correlations](exact-correlations.png)

For periodic lengths at least three, **any two single-octahedron catalog variables at non-neighboring centres are independent**. All connected local-object correlations therefore vanish beyond nearest-neighbour centres. In particular, the x-dimer covariance at the x neighbour is -1/16, the same specified face-trimer covariance there is +1/256, and the any-trimer covariance is zero for every distinct centre. The exact zeros are not a numerical floor and must not be fitted to a positive exponential length.

Neighboring octahedra have 128 compatible catalog pairs, each of probability 1/128; compatibility requires opposite ownership of the shared site. Non-neighbors have all 256 catalog pairs with probability 1/256. The full matrices are saved in `data.json`. Period-two clusters exhibit additional finite-size correlations and are used only for validation.

## Counting, sectors and nonlocal structure

Let N=Lx Ly Lz and S=Lx+Ly+Lz. There are 3N ownership bits and 2N parity equations, with S-1 linear dependencies. Hence the rank is 2N-S+1 and the number of coverings is **2^(N+S-1)**. This formula counts distinct owner-labelled objects; small periodic tori can identify physical endpoint sets while retaining different geometric object locations.

The elementary cube-flip span has rank N-S+2, leaving **2^(2S-3)** classical sectors inaccessible to cube-only updates. Independent free-bit sampling in this code includes all sectors and has no equilibration or local-update trapping.

Even though two-point local tails vanish, nonlocal constraints remain. If T_v indicates a trimer, its XOR along every periodic coordinate line equals the line length modulo two. The statistical ensemble is therefore not a product of independent octahedron states. Neither finite-range local correlations nor counting classical cube-move sectors by itself establishes a quantum fracton phase.

The counting formulas follow by writing a dependency of the xy and xz checks as lambda_xy=h(x)+f(z), lambda_xz=h(x)+g(y), with one redundant common constant. Cube-flip relations have coefficients u(x)+v(y)+w(z), with two redundant constants. Local independence can be seen by surrounding-cube flips: away from nearest neighbours they span the full four-bit catalog at either octahedron without changing the other.

## Relation to the cited papers

[Benton and Moessner, arXiv:2103.10817](https://arxiv.org/abs/2103.10817), Figure 3, supports the corner-sharing octahedral geometry. [Stern et al., arXiv:2603.12313v2](https://arxiv.org/abs/2603.12313v2), Section V, instead requires equal incoming counts along the three axes as integers: empty, eight faces and full, giving ten states. The present rule equates these counts only modulo two and includes six additional partial-dimer supports. The cited U(1) model's entropy, algebraic correlations and pinch-line features cannot be assigned to this sixteen-state model.

## Reproduction and next calculations

The independent validation in `tests/test_octahedral.py` enumerates geometric
coverings directly on a 2x2x2 periodic cluster, without using the binary
equations. It finds 8,192 distinct owner-labelled coverings, in one-to-one
agreement with the arrow solutions, and checks their joint probabilities.
The 19 tests also check exact site coverage, local rules, cube moves, rank
formulas, finite-range correlations and the nonlocal line-parity constraint.

```bash
/workspace/peps-environment/.venv/bin/python scripts/run_octahedral_classical.py
/workspace/peps-environment/.venv/bin/python -m pytest -q tests/test_octahedral.py
```

The implementation is `rvb/octahedral.py`; `data.json` contains exact fractions, timings and source hashes, and `sample-covering.json` contains a checked 4x4x4 covering with every site covered once. The PDF is [Octahedral_Classical_First_Results.pdf](Octahedral_Classical_First_Results.pdf).

Possible next steps are weighted dimers/trimers, restricted local catalogs, nonlocal cage observables and sector comparisons. The parity solver gives exact unweighted samples, but nonuniform weights require a separate weighted algorithm or controlled reweighting. The quantum interpretation remains postponed as requested.
"""
    (out/"report.md").write_text(report)
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
