"""Plot finite-box quantum spin estimates; never fit a bulk decay length.

The errors are delete-one-block jackknife errors after discarding the first
quarter of each run. They do not account for incomplete mixing or finite size.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/octahedral_bulk_overlap"
OUTPUT = ROOT / "results/octahedral-spin1-bulk"


def estimate(records, name):
    weight = np.array([row.get("weight", 1.0) for row in records])
    sign = np.array([row["sign"] for row in records])
    values = np.array([row[name] for row in records])
    numerator, denominator = weight * values, weight * sign
    total_num, total_den = numerator.sum(), denominator.sum()
    if not total_den:
        return {"mean": None, "block_error": None}
    mean = total_num / total_den
    blocks = np.array_split(np.arange(len(records)), min(20, len(records)))
    jackknife = np.array([
        (total_num - numerator[block].sum()) /
        (total_den - denominator[block].sum()) for block in blocks
    ])
    error = np.sqrt((len(blocks) - 1) * np.mean((jackknife - jackknife.mean()) ** 2))
    return {"mean": float(mean), "block_error": float(error),
            "nonzero_samples": int(np.count_nonzero(values))}


def main():
    runs = []
    for path in sorted(SOURCE.glob("mcmc_*.json")):
        raw = json.loads(path.read_text())
        if raw["L"] < 4 or raw["steps"] < 20_000:
            continue
        records = raw["records"][len(raw["records"]) // 4:]
        run = {key: raw.get(key, 1.0 if key == "bias" else None)
               for key in ("L", "steps", "seed", "initial", "bias", "seconds")}
        run.update(source=str(path.relative_to(ROOT)),
                   measurements=len(records), burn_in_fraction=.25,
                   moves=raw["moves"], negative_proposals=raw["negative_proposals"],
                   correlations={name: estimate(records, name)
                                 for name in records[0] if name.startswith("spin_")},
                   raw_trimer_mismatch_fraction=float(np.mean(
                       [row["trimer_difference"] > 0 for row in records])))
        runs.append((run, raw))
    if not runs:
        raise RuntimeError("No completed even-box run with at least 20,000 steps")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary = {
        "status": "Exploratory finite boxes; incomplete mixing and finite-size errors remain.",
        "observable": "Physical S_i dot S_j, averaged over origins and cubic axes",
        "error_definition": "20-block jackknife after 25% burn-in; statistical errors only",
        "correlation_length": None,
        "runs": [run for run, _ in runs],
    }
    (OUTPUT / "spin-sampling-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), layout="constrained")
    colors = plt.get_cmap("tab10")
    zero_sample_points = False
    for index, (run, _) in enumerate(runs):
        method = ("enhanced sampling" if run["bias"] != 1
                  else "run B" if run["initial"] == "independent" else "run A")
        label = f'{run["L"]}³ ({3 * run["L"] ** 3} spins), {method}'
        distances = np.arange(1, run["L"] // 2 + 1)
        offset = (index - (len(runs) - 1) / 2) * .035
        for axis, direction in zip(axes, ("parallel", "perpendicular")):
            values = [run["correlations"][f"spin_{direction}_{r}"] for r in distances]
            axis.errorbar(distances + offset, [v["mean"] for v in values],
                          yerr=[v["block_error"] for v in values],
                          fmt="o", capsize=3, color=colors(index), label=label)
            # Zero estimates from never seeing a nonzero contribution cannot
            # establish a zero correlator; mark them separately, without a bar.
            for distance, value in zip(distances, values):
                if value["nonzero_samples"] == 0:
                    axis.plot(distance + offset, 0, "x", color="black", ms=8)
                    zero_sample_points = True
    for axis, title in zip(axes, ("Separation along the spin's coordinate line",
                                  "Separation perpendicular to that line")):
        axis.axhline(0, color=".55", lw=.8)
        axis.set(title=title, xlabel="Distance r (cubic-centre spacings)",
                 ylabel=r"$\langle\mathbf{S}_i\cdot\mathbf{S}_{i+r}\rangle$")
        axis.set_xticks(np.arange(1, max(r["L"] for r, _ in runs) // 2 + 1))
        axis.grid(alpha=.18)
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("Spin-1 correlations on periodic 3D boxes — no bulk ξ extracted",
                 fontsize=14, weight="bold")
    caption = ("Both sampling methods target the same quantum state. "
               "Bars: statistical errors only; size and mixing errors remain.")
    if zero_sample_points:
        caption += " Black ×: no nonzero contribution sampled, not proof of zero."
    fig.supxlabel(caption, fontsize=9)
    fig.savefig(OUTPUT / "spin-sampling.png", dpi=170)
    fig.savefig(OUTPUT / "spin-sampling.pdf")
    plt.close(fig)
    # The optional paper dependency is used only to assemble the two
    # already generated vector PDF pages into the user's requested file.
    from pypdf import PdfWriter
    combined = PdfWriter()
    for name in ("geometry.pdf", "spin-sampling.pdf"):
        combined.append(OUTPUT / name)
    combined.write(OUTPUT / "Three_Dimensional_Spin1_Status.pdf")


if __name__ == "__main__":
    main()
