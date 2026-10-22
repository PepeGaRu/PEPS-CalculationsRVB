# Revised collaborator notes

Open **Collaborator_Notes.pdf** (four A4 pages, revision 3).

This revision draws all ownership arrows in Figure 5, keeps only panels (a)
and (b) in Figure 6, removes the section headings, gathers all references at
the end, and explicitly states closure of the sixteen local patterns under XOR.

1. Upright gray square plaquettes, physical spin positions, and all eight local
   dummy-link coverings. Table triangles and gray edge corners follow the lattice
   orientation and the contracted square-loop ports.
2. Covering-to-loop mapping with green dashed triangle frames, plus the calculated
   2D bond correlations and decay diagnostics.
3. Complete geometric definition of the 3D octochlore covering model, its
   tetrahedral origin, all sixteen local arrow configurations, and their closure.
4. Cage-net bijection, normalized spin-1 singlet extension, and the corrected
   layer projection: the two vertical apex sites coincide in the top view.
   The correspondence to the two dummy-link slots is explained in the text.

Pages 3–4 incorporate the scientific content of José Garre Rubio's supplied
notes of 19 October 2026. This document can replace those notes.
The new diagram establishes the proposed geometric site-incidence connection.
It does not yet establish a quantum foliation: the full 3D covering rules differ
and neighbouring layers share vertical sites.

The 2D decay scale remains a finite-window estimate. The 3D bulk correlation
length is unresolved. This revision reuses the saved correlation data; no new
large tensor contractions were run.

## Edit and compile

Keep `Collaborator_Notes.tex` beside the seven figure PDFs. From this directory:

```bash
pdflatex -interaction=nonstopmode -halt-on-error Collaborator_Notes.tex
pdflatex -interaction=nonstopmode -halt-on-error Collaborator_Notes.tex
```

The source ZIP is self-contained for PDF compilation. It also contains the
plotting scripts and audit records. Regenerating figures requires the full
project and its configured Python environment:

```bash
python scripts/plot_collaborator_mapping.py
python scripts/plot_collaborator_correlations.py
python scripts/plot_collaborator_octahedral_catalogue.py
python scripts/plot_collaborator_layer.py
```

The correlation plot uses saved data in `results/octagon-square/ctmrg/`;
the files and settings are recorded in `correlations-2d-audit.json`.
The global loop example reuses `results/octagon-square/mapping-example.json`
and verifies its inverse dummy-link reconstruction. Local table supports
are recorded in `dummy-covering-table.json`; projected and opened coordinates
are specified in `layer-projection-coordinates.json`.

The bundle includes these audit documents:

- `docs/collaborator-2d-mapping-audit.md`
- `docs/octahedral-source-replacement.md`
- `docs/octahedral-layer-audit-v2.md`

The earlier three-page PDF is retained in the repository under
`archive/Collaborator_Notes_v1.pdf`.
