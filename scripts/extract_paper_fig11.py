#!/usr/bin/env python3
"""Recover displayed Fig.11 points from arXiv:1203.4816v3 PDF vector paths.

This is digitization of the published drawing, not original simulation data.
The script is deliberately specific to the uploaded v3 PDF and validates the
path counts and calibrated axes before assigning the ordered theta labels.
Requires pypdf; it is an optional paper-audit dependency, not a simulation one.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from pypdf import PdfReader


def extract(pdf: Path) -> dict:
    operations = PdfReader(pdf).pages[6].get_contents().operations
    # Raw drawing coordinates: r=1 and 2; log10(C)=-6 and -1.
    # These are the axis tick locations in this PDF's vector content stream.
    x1, x2 = 3630.52, 4191.34
    y_bottom, y_top = 1381.58, 2598.48
    curves = {60: [], 40: []}
    color = None
    path = []
    axis_frame_found = False
    for args, op in operations:
        if op == b"RG":
            color = tuple(map(float, args))
        if op == b"m":
            path = [tuple(map(float, args))]
        elif op == b"l":
            path.append(tuple(map(float, args)))
        elif op == b"c":
            # Circular markers are Bezier paths, not the correlation polyline.
            path = []
        elif op == b"S":
            if (len(path) == 5 and abs(path[0][0] - x1) < .02
                    and abs(path[0][1] - y_bottom) < .02
                    and abs(path[1][0] - x1) < .02
                    and abs(path[1][1] - y_top) < .02
                    and abs(path[2][1] - y_top) < .02
                    and abs(path[3][1] - y_bottom) < .02):
                axis_frame_found = True
            if (len(path) >= 2 and abs(path[0][0] - x1) < .02
                    and 2500 < path[0][1] < y_top
                    and path[1][0] > x1 + 100):
                dimension = 60 if color == (1., 0., 0.) else 40 if color == (0., 0., 0.) else None
                if dimension:
                    curves[dimension].append(path)
            path = []
    if not axis_frame_found:
        raise ValueError("Expected calibrated Figure 11 axis frame is absent")
    if [len(curves[d]) for d in (60, 40)] != [10, 10]:
        raise ValueError("Expected ten red and ten black correlation paths; PDF layout differs")
    output = []
    for dimension in (40, 60):
        paths = curves[dimension]
        # Figure labels list theta=.9,.8,...,0 from steepest to shallowest.
        # At the first common x the plotted polylines must have that ordering.
        slopes = [(p[1][1] - p[0][1]) / (p[1][0] - p[0][0]) for p in paths]
        if any(b < a for a, b in zip(slopes, slopes[1:])):
            raise ValueError("Curve drawing order does not match labeled slope order")
        for index, path in enumerate(paths):
            points = []
            for x, y in path:
                distance = 1 + (x - x1) / (x2 - x1)
                # Clipping intersections at the lower y boundary and beyond
                # r=4 are line endpoints, not displayed integer-distance data.
                if abs(distance - round(distance)) < .001 and 1 <= round(distance) <= 4:
                    log_value = -6 + 5 * (y - y_bottom) / (y_top - y_bottom)
                    points.append({"r": round(distance), "correlation": 10 ** log_value,
                                   "raw_pdf_xy": [x, y]})
            output.append({"theta": round(.9 - index / 10, 1),
                           "boundary_dimension": dimension, "points": points})
    return {
        "source": "arXiv:1203.4816v3, PDF page 7, Figure 11",
        "source_pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "method": "Extraction of plotted polyline vertices from the PDF vector content stream",
        "data_status": "Published plotted coordinates, NOT original raw numerical data or a new simulation",
        "observable_as_printed": "Connected correlation of two z tensor z operators, z=diag(1,-1,1)",
        "geometry_as_printed": "Infinite system; dimers on two blocked Fig.10b tensors in the same row",
        "calibration": {"r1_raw_x": x1, "r2_raw_x": x2,
                        "log10_minus6_raw_y": y_bottom, "log10_minus1_raw_y": y_top},
        "precision_note": "Coordinate rounding and plotting precision limit accuracy; report about 3 significant figures. No statistical error bars are available. Values clipped below 1e-6 are omitted, not zero.",
        "theta_assignment": "Paper labels order curves theta=0.9 down to 0 from steepest to shallowest; script checks this order from first line-segment slopes.",
        "curves": output,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = extract(args.pdf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")
    print(f"Extracted {sum(len(c['points']) for c in data['curves'])} plotted points in {len(data['curves'])} curves to {args.output}")
