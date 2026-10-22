#!/usr/bin/env python3
"""Exact spin-1 quantum validation on the exceptional 2x2x2 periodic cluster."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.octahedral_spin1 import Spin1OctahedralTwoCube


def number(value):
    return {"exact": str(value), "value": float(value)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results/octahedral-spin1")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    begin = perf_counter()
    model = Spin1OctahedralTwoCube()
    norm = model.norm_components()
    means = {
        "diameter_singlet_projector": model.dimer_projector_expectation((0,)),
        "face_singlet_projector": model.trimer_projector_expectation((0,)),
    }
    # This comparison uses the SAME physical projectors in an incoherent
    # mixture of normalized spin covering states. It is not an occupation
    # probability in the earlier orthogonal covering model.
    mixture_means = {
        "diameter_singlet_projector": Fraction(5, 9),
        "face_singlet_projector": Fraction(1, 9),
    }
    pair_cases = [
        ("meeting_at_vertex", 4, "Meet at\none vertex"),
        ("parallel_face", 2, "Parallel,\none face"),
        ("parallel_opposite", 3, "Parallel,\nopposite"),
        ("skew", 5, "Skew\nedges"),
    ]
    pair_records = []
    for name, edge, label in pair_cases:
        pair_records.append({
            "geometry": name, "plot_label": label, "pair_edges": [0, edge],
            "connected": number(model.dimer_projector_covariance(0, edge)),
            "incoherent_connected": number(Fraction(16, 81)),
        })
    face_records = []
    for vertex, label in ((1, "Edge\nseparation"), (3, "Face\ndiagonal"), (7, "Body\ndiagonal")):
        face_records.append({
            "centres": [[0, 0, 0], list(model.cells[vertex])],
            "plot_label": label,
            "connected": number(model.trimer_projector_covariance(0, vertex)),
        })
    seconds = perf_counter()-begin
    result = {
        "model": "Equal-coefficient spin-1 singlet dimer/trimer covering superposition",
        "shape": [2, 2, 2], "octahedra": 8, "physical_spin1_sites": 24,
        "owner_labelled_coverings": 8192,
        "dimer_only_coverings": 4096, "trimer_only_coverings": 4096,
        "singlet_convention": "Cartesian delta/sqrt(3) and epsilon/sqrt(6); trimer sites ordered x,y,z for every face.",
        "weights": "Equal coefficients for products of normalized singlets, including all covering sectors.",
        "method": "Exact integer tensor contractions and rational arithmetic",
        "norm_components": {key: number(value) for key, value in norm.items()},
        "means": {key: number(value) for key, value in means.items()},
        "incoherent_mixture_means": {key: number(value) for key, value in mixture_means.items()},
        "incoherent_mixture_definition": "rho=(1/8192) sum_C |C><C|; every spin covering |C> normalized.",
        "diameter_projector_covariances": pair_records,
        "reference_face_projector_covariances": face_records,
        "cube_pair_edges": model.edges,
        "reference_faces": "All-plus faces, with one outgoing physical site along each axis.",
        "period_two_artifact": "All 4096 dimer coverings represent one physical spin product; amplitudes add coherently.",
        "bulk_correlation_length": None,
        "bulk_decay_classification": "Undetermined: this exceptional small torus cannot diagnose a 3D bulk tail.",
        "calculation_seconds": seconds,
        "timing_scope": "Norm, means and saved connected correlations; excludes tests and plotting.",
        "source_sha256": {
            path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in ("rvb/octahedral_spin1.py", "scripts/run_octahedral_spin1.py")
        },
    }
    (out/"data.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")

    os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.7), layout="constrained")
    x = np.arange(2)
    coherent = [float(v) for v in means.values()]
    incoherent = [float(v) for v in mixture_means.values()]
    axes[0].bar(x-.18, coherent, width=.36, color="#276aa4", label="Coherent singlet sum")
    axes[0].bar(x+.18, incoherent, width=.36, color="#a4b2bc", label="Incoherent spin mixture")
    for i, (q, c) in enumerate(zip(coherent, incoherent)):
        axes[0].text(i-.18, q+.022, f"{q:.4f}", ha="center", fontsize=9)
        axes[0].text(i+.18, c+.022, f"{c:.4f}", ha="center", fontsize=9)
    axes[0].set(xticks=x, xticklabels=["Diameter-pair\nsinglet projector", "Face-trimer\nsinglet projector"],
                ylabel="Expectation value", title="Same physical operators", ylim=(0, 1.24))
    axes[0].legend(loc="upper right", fontsize=8)
    for ax, records, title, color in (
        (axes[1], pair_records, "Quantum pair-projector correlations", "#276aa4"),
        (axes[2], face_records, "Quantum trimer-projector correlations", "#824c9d"),
    ):
        values = [r["connected"]["value"] for r in records]
        ax.bar(range(len(records)), values, color=color, width=.62)
        ax.set(xticks=range(len(records)), xticklabels=[r["plot_label"] for r in records],
               title=title, ylabel=r"$\langle P_A P_B\rangle-\langle P_A\rangle\langle P_B\rangle$")
        ax.tick_params(axis="x", labelsize=8)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.set_ylim(0, max(values)*1.27)
        for i, value in enumerate(values):
            ax.text(i, value*1.025, f"{value:.3g}", ha="center", fontsize=8)
    for ax in axes:
        ax.grid(axis="y", alpha=.18)
    fig.suptitle("Spin-1 octahedral coverings: exact 24-spin validation only\n"
                 "2 × 2 × 2 periodic cluster — not a bulk decay or correlation-length calculation",
                 fontsize=13)
    fig.savefig(out/"quantum-validation.png", dpi=170, bbox_inches="tight")
    fig.savefig(out/"quantum-validation.pdf", bbox_inches="tight")

    overview = plt.figure(figsize=(8.27, 11.69))
    overview.text(.085, .94, "Spin-1 dimer and trimer singlets", fontsize=21, weight="bold", color="#276aa4")
    overview.text(.085, .895, "Exact quantum construction and first validation", fontsize=15)
    sections = [
        (.84, "Why spin 1 is the direct choice",
         "Both two and three spin-1 particles have a unique total-spin-zero state.\n"
         "Use normalized dimers and trimers with equal covering coefficients.\n"
         "The trimer is antisymmetric; its sites are ordered x, y, z on every face."),
        (.70, "The classical result does not transfer",
         "Different singlet coverings overlap. A single cube move between\n"
         "eight trimers and twelve dimers gives overlap 1/39,366.\n"
         "An exact tensor network with virtual dimension four includes these overlaps."),
        (.56, "What was actually calculated",
         "All 8,192 coverings on a periodic 2 × 2 × 2 cluster of 24 spin-1 sites.\n"
         "Norm and physical singlet-projector correlations use exact arithmetic.\n"
         f"Pair singlet expectation: {coherent[0]:.8f}. Face singlet: {coherent[1]:.8f}.\n"
         f"The reported quantum contractions took {seconds:.2f} seconds."),
        (.40, "Essential finite-size limitation",
         "On this torus all 4,096 dimer coverings become the same physical state.\n"
         "Their coherent enhancement strongly affects the measured values.\n"
         "These are construction checks. The bulk quantum correlation length\n"
         "and exponential-versus-power-law decay remain undetermined."),
        (.23, "What spin one-half would require",
         "Three spin-one-half particles have two doublets and one quartet,\n"
         "but no singlet. Choosing trimer spin 1/2 or 3/2 leaves an effective spin.\n"
         "Its internal state and coupling to other trimers would need specifying.\n"
         "That defines another quantum model, not an automatic simplification."),
    ]
    for y, title, body in sections:
        overview.text(.085, y, title, fontsize=12.5, weight="bold")
        overview.text(.085, y-.025, body, fontsize=10.2, va="top", linespacing=1.5)
    overview.text(.085, .065, "Reproduce: python scripts/run_octahedral_spin1.py\n"
                  "Exact fractions, source hashes and method details accompany this PDF.",
                  fontsize=9, color="#68747f")
    pdf_path = out/"Octahedral_Spin1_Quantum_Checks.pdf"
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(overview)
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(overview)
    plt.close(fig)

    report = f"""# Spin-1 singlet coverings: exact construction and a small quantum benchmark

Spin 1 gives a unique singlet for both dimers and trimers. We have implemented
the equal-coefficient covering superposition, including nonorthogonal overlaps.
The [model definition](../../docs/octahedral-spin1-model.md) specifies the
normalized tensors, phases and exact tensor-network representation.

**This is a 24-spin validation, not a calculation of the bulk quantum
correlation length.** The periodic 2x2x2 geometry has a severe identification
of dimer states, described below. Its correlations must not be fitted to a
bulk exponential or power law.

## State and observables

Each covering is a product of normalized spin-1 singlets: a symmetric
Cartesian delta/sqrt(3) for each dimer and epsilon/sqrt(6) for each trimer.
Trimer sites are ordered by the x, y, z axes for all eight faces.
All complete owner-labelled coverings have equal wavefunction coefficients.
There are no added orthogonality flags, fugacities or Hamiltonian assumptions.

The rank-six octahedron tensor and shared-site ownership tensor have virtual
dimension four (one vacuum plus three spin states). They define the state
on larger lattices too; the exact contraction below exploits a special small
torus and does not provide a scalable bulk contraction algorithm.

We measure physical singlet projectors, not classical object occupations:
`P_D = ((S_i · S_j)^2 - I)/3` for a diameter pair, and
`P_T = |t><t|` for a face trimer. Connected correlations are
`<P_A P_B> - <P_A><P_B>`.

## Exact overlap and complete covering sum

For the eight-trimer and twelve-dimer configurations related by one cube flip,
the overlap is **1/39,366**. It follows from the 24 positive three-edge-colourings
of the cube and the normalization factor `3^6 * 6^4`. This local overlap check
also applies on a larger torus.

On the 2x2x2 torus, the 8,192 coverings divide into 4,096 all-dimer coverings
and 4,096 all-trimer coverings. Let D be the normalized dimer product and T
the all-plus reference trimer product. The full unnormalized state is
`4096 D + G T`, where `G` is the product of `(I + Swap)` on twelve disjoint
physical spin pairs.

| Contribution to the unnormalized norm | Exact value |
|---|---:|
| Dimer part | {norm["dimer_part"]} |
| Trimer part | {norm["trimer_part"]} |
| Dimer–trimer interference | {norm["interference"]} |
| Total norm squared | {norm["total"]} |

These are integer tensor contractions with normalization restored as rational
fractions. No sampled covering average, tensor truncation or full `3^24`
state vector is used. Computing the reported norm, means and correlations
took {seconds:.3f} seconds in this environment, excluding plotting and tests.

## Physical projector results

| Observable | Coherent singlet superposition | Incoherent spin-covering mixture |
|---|---:|---:|
| Diameter-pair singlet projector | {coherent[0]:.10f} | 5/9 = {incoherent[0]:.10f} |
| Specified face-trimer singlet projector | {coherent[1]:.10f} | 1/9 = {incoherent[1]:.10f} |

The comparison mixture is `rho = (1/8192) sum_C |C><C|` of normalized
**physical spin states**, with the same projectors on both sides.
It is different from the earlier orthogonal-label model. For example,
the classical occupation probability 1/4 for a particular owner-labelled
dimer is not a physical singlet-projector expectation.

In the incoherent mixture, half the coverings are all dimers. Their physical
diameter pairs are singlets; the other half give 1/9 for a chosen pair.
For a face projector, the dimer sector gives 1/27. In the trimer sector,
the eight ownership patterns give one case with expectation 1, three with
1/9, and four with 1/27, averaging to 5/27. Hence the full mixture's face
expectation is 1/9.

![Exact quantum validation results](quantum-validation.png)

The diameter pairs are the twelve disjoint spin pairs represented by edges
of a cube. The four plotted classes describe their relative cube-edge
geometry, not four distances along a bulk lattice direction. All-plus
reference faces at different centres also have disjoint physical supports;
their plotted projector products commute. Exact fractions for every plotted
value are saved in `data.json`.

## What the finite-size artifact means

Opposite geometric chords have the same physical endpoints on a period-two
line. Therefore all 4,096 dimer coverings give the **same** normalized spin
state, so their amplitudes add before the norm is taken. This strongly enhances
the dimer part and explains why its pair-projector mean is near one.
That number is not evidence for bulk dimer order. The nonzero connected
correlations shown here likewise do not establish a bulk decay law.

The earlier exact finite-range result remains a result about uniform
classical covering occupations. It cannot be carried over to the spin-1
superposition. The **bulk quantum correlation length remains undetermined**;
it requires larger geometry with controlled contraction or a separately
validated sampling method.

## Spin one-half alternative

Three spin-one-half particles decompose into two spin-1/2 doublets and one
spin-3/2 quartet. Neither choice is a singlet. A doublet also has a
multiplicity degree of freedom, so choosing its total spin alone does not
fix the internal trimer state. The remaining effective spins on trimers
need an additional state or coupling prescription. Spin 1 is the direct
choice when every local dimer and trimer should itself be a singlet.

## Reproduce and inspect

The independent tests check total spin zero and singlet projectors using
spin matrices, enumerate all `3^12` cube-edge colourings for the overlap,
and compare tensor contractions with dense twelve-spin calculations.
They also verify all 8,192 covering labels, the special dimer-state
identification, and SU(2) invariance of the general dimension-four tensors.

```bash
/workspace/peps-environment/.venv/bin/python scripts/run_octahedral_spin1.py
/workspace/peps-environment/.venv/bin/python -m pytest -q tests/test_octahedral_spin1.py
```

- [PDF with the quantum checks and plots](Octahedral_Spin1_Quantum_Checks.pdf)
- [Exact fractions and calculation metadata](data.json)
- [State definition and tensor construction](../../docs/octahedral-spin1-model.md)
- [Quantum implementation](../../rvb/octahedral_spin1.py)
- [Classical results, a different ensemble](../octahedral-classical/report.md)
"""
    (out/"report.md").write_text(report)
    print(f"Calculation: {seconds:.3f} seconds")
    print(f"Pair projector: {means['diameter_singlet_projector']} = {coherent[0]:.10f}")
    print(f"Face projector: {means['face_singlet_projector']} = {coherent[1]:.10f}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
