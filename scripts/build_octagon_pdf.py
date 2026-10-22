#!/usr/bin/env python3
"""Build a self-contained, bookmarked PDF from the calculated vector figures."""
from __future__ import annotations

import hashlib
from io import BytesIO
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/workspace/peps-environment/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter, PageObject, Transformation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/octagon-square"
DEST = OUT / "Square_Octagon_RVB_Coverings_and_Correlations.pdf"
BLUE = "#245e9d"

FIGURES = [
    ("spin-lattice.pdf", "Physical lattice and spin positions",
     "Every colored dot is one spin 1/2. Blue: one per small-square edge; orange: two per dummy link.",
     "Green dashed triangles show allowed singlet bonds. Hollow gray vertices are not extra spins."),
    ("spin-unit-cell.pdf", "Enlarged cell: eight physical spins",
     "One simulation cell contains four blue spins and four orange endpoint spins.",
     "Each orange endpoint has a partner on the neighboring cell; both belong to the same dummy link."),
    ("dummy-link-states.pdf", "The two allowed dummy-link arrangements",
     "Either the two orange spins pair with one another, or both pair inside their respective triangles.",
     "There are no empty physical spin sites. All spins are paired in a complete covering."),
    ("mapping-singlets.pdf", "One explicit allowed singlet covering",
     "Green segments are occupied singlet bonds; each colored dot belongs to exactly one singlet.",
     "This is a 3 x 3 periodic example. Outgoing boundary stubs continue across the periodic boundary."),
    ("dummy-contraction.pdf", "Corrected mapping: contract the dummy links",
     "Covering -> reference-relative square-octagon loops -> shrink dummy links into square-lattice vertices.",
     "This map is one-to-one with inherited routing. It does not collapse the small squares."),
    ("correlations.pdf", "Calculated connected bond correlations",
     "Left: physical singlet RVB bond covariance. Right: exact orthogonal-covering occupation covariance.",
     "The RVB decay scale is about 0.8-0.9 simulation cells; the orthogonal correlations have exact finite range."),
    ("decay-diagnostics.pdf", "Diagnosis of exponential versus power-law decay",
     "A positive inverse local length and an increasing effective power exponent favor exponential-like decay.",
     "This conclusion concerns the resolved distances. Finite numerical truncation does not prove the bulk asymptote."),
    ("cylinder-lengths.pdf", "Independent finite-cylinder correlation lengths",
     "Exact transfer calculations at circumferences 2, 3 and 4 keep even and odd parity sectors separate.",
     "Width and sector dependence remain. The shaded CTMRG range is descriptive, not a confidence interval."),
]


def figure_pdf(fig):
    buffer = BytesIO()
    fig.savefig(buffer, format="pdf", facecolor="white")
    plt.close(fig)
    buffer.seek(0)
    return PdfReader(buffer).pages[0]


def add_figure(writer, filename, title, caption1, caption2, page_number, total):
    source = PdfReader(OUT / filename).pages[0]
    sw, sh = float(source.mediabox.width), float(source.mediabox.height)
    width, height = (595.276, 841.89) if sh > sw else (841.89, 595.276)
    fig = plt.figure(figsize=(width / 72, height / 72))
    fig.text(36/width, 1-35/height, title, fontsize=17, color=BLUE, va="top")
    fontsize = 9 if width < 700 else 10
    fig.text(36/width, 59/height, caption1, fontsize=fontsize, va="top")
    fig.text(36/width, 43/height, caption2, fontsize=fontsize, va="top")
    fig.text(36/width, 18/height, "Square-octagon RVB | New calculations for the specified covering rule",
             fontsize=8, color="#68747f")
    fig.text(1-36/width, 18/height, f"{page_number} / {total}",
             fontsize=8, color="#68747f", ha="right")
    frame = figure_pdf(fig)
    scale = min((width-72)/sw, (height-145)/sh)
    tx = (width-sw*scale)/2
    ty = 80 + (height-145-sh*scale)/2
    page = PageObject.create_blank_page(width=width, height=height)
    page.merge_page(frame)
    page.merge_transformed_page(source, Transformation().scale(scale).translate(tx, ty))
    writer.add_page(page)
    writer.add_outline_item(title, page_number-1)


def main():
    summary = json.loads((OUT / "figure-summary.json").read_text())
    models = summary["fit_window"]["models"]
    total = len(FIGURES) + 1
    writer = PdfWriter()
    fig = plt.figure(figsize=(8.2677, 11.6929))
    fig.text(.09, .92, "SQUARE-OCTAGON RVB", color=BLUE, fontsize=23, weight="bold")
    fig.text(.09, .865, "Lattice, singlet coverings\nand correlation lengths", fontsize=24, va="top")
    paragraphs = [
        (.75, "The state calculated", 14, "bold"),
        (.714, "Equal coefficients for all allowed complete singlet coverings.\n"
               "Triangle singlets point counterclockwise; dummy singlets point east/north.\n"
               "An unused dummy link contains a singlet between its two endpoint spins.\n"
               "No extra semionic loop-count sign is included.", 10.8, "normal"),
        (.605, "Main results", 14, "bold"),
        (.568, "Orthogonal coverings: occupation correlations are exactly zero\n"
               "beyond neighboring simulation cells.\n\n"
               "Physical singlet RVB: the calculated bond correlations support\n"
               "a short exponential-like decay scale of about 0.8-0.9 simulation cells.",
         11.5, "normal"),
        (.425, "How precise is that length?", 14, "bold"),
        (.389, f"For distances r = 12 to 18, a pure exponential fit gives xi = "
               f"{models['exponential']['xi']:.3f}.\n"
               f"Allowing a power prefactor gives xi = {models['power_exponential']['xi']:.3f}.\n"
               "These are finite-window fits, not precision bulk estimates. Tightening\n"
               "the numerical cutoff changes late correlations by approximately 3-8%.\n"
               "The exact cylinders also retain circumference and sector dependence.",
         10.8, "normal"),
        (.26, "Units and the corrected loop mapping", 14, "bold"),
        (.224, "One simulation cell contains one small square and eight physical spins.\n"
               "The measured translation is diagonal on the square graph obtained\n"
               "by contracting dummy links. The contraction is one-to-one relative\n"
               "to a reference covering, with inherited connections at four-leg vertices.\n"
               "Earlier small-square-collapse diagrams described a different projection.",
         10.8, "normal"),
    ]
    for y, content, size, weight in paragraphs:
        fig.text(.09, y, content, fontsize=size, weight=weight, va="top", linespacing=1.5)
    fig.text(.09, .075, "All curves are new calculations; no published correlation curves were copied.\n"
             "135 validation tests passed. Figures and source data are retained in the repository.",
             fontsize=9.2, color="#52616d", linespacing=1.5)
    fig.text(.91, .03, f"1 / {total}", ha="right", fontsize=8, color="#68747f")
    writer.add_page(figure_pdf(fig))
    writer.add_outline_item("Results, definitions and limitations", 0)
    for index, (filename, title, caption1, caption2) in enumerate(FIGURES, 2):
        add_figure(writer, filename, title, caption1, caption2, index, total)
    writer.add_metadata({
        "/Title": "Square-Octagon RVB: Coverings and Correlation Lengths",
        "/Subject": "Specified spin lattice, singlet coverings, corrected dummy-link contraction and calculated correlations",
        "/Author": "PEPS-CalculationsRVB",
        "/Keywords": "RVB, PEPS, square-octagon, singlet, dimers, correlation length",
    })
    writer.write(DEST)
    check = PdfReader(DEST)
    assert len(check.pages) == total
    for i, page in enumerate(check.pages):
        assert len(page.extract_text()) > 100, f"Unexpectedly empty page {i+1}"
    assert "0.8-0.9" in check.pages[0].extract_text()
    manifest = {
        "pdf": DEST.name, "pages": total, "bytes": DEST.stat().st_size,
        "sha256": hashlib.sha256(DEST.read_bytes()).hexdigest(),
        "source_figures": {name: hashlib.sha256((OUT/name).read_bytes()).hexdigest()
                           for name, *_ in FIGURES},
        "numerical_provenance": "figure-summary.json",
        "new_simulations_run_for_this_pdf": False,
    }
    (OUT / "pdf-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
