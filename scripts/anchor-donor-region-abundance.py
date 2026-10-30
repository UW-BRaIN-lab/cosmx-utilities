"""Cell type abundance histograms for one donor, one panel per tissue region.

Reuses the label loading and the de novo / named split from
anchor-slide-abundance.py. Cells are assigned to a donor and region by
joining the curated per-FOV annotation sheet (via the slide crosswalk),
because a slide carries two tissue sections that can belong to different donors.
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SIBLING_SCRIPT = Path(__file__).with_name("anchor-slide-abundance.py")
OTHER_NAMED_BAR = "other named"
MIN_NAMED_FRACTION = 0.01
REGION_ORDER = ["Contralateral uninvolved", "Infiltrating edge", "Tumor bulk"]
PANEL_WIDTH_IN, PANEL_HEIGHT_IN = 14.0, 3.2
PNG_DPI = 200
LETTER_COLOR, NAMED_COLOR, OTHER_COLOR = "#4a6fa5", "#2e8b57", "#9a9a9a"
HIGHLIGHT_LETTER, HIGHLIGHT_COLOR = "b", "#c0392b"

spec = importlib.util.spec_from_file_location("anchor_slide_abundance", SIBLING_SCRIPT)
sibling = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sibling)


def donor_cells(cells, annotations_csv, crosswalk_csv, donor):
    annotations = pd.read_csv(annotations_csv)
    crosswalk = pd.read_csv(crosswalk_csv, usecols=["Canonical_slide_name", "AtoMx_flatfile_folder"])
    annotations = annotations.merge(crosswalk, left_on="Slide", right_on="Canonical_slide_name",
                                    validate="many_to_one")
    annotations = annotations[annotations.Case.astype(str) == str(donor)]
    parsed = cells.cell_id.str.extract(r"^(?P<folder>.+)_F(?P<fov>\d+)_C\d+$")
    cells = cells.assign(folder=parsed.folder, fov=parsed.fov.astype(int))
    return cells.merge(annotations[["AtoMx_flatfile_folder", "FOVs", "Block", "Region"]],
                       left_on=["folder", "fov"], right_on=["AtoMx_flatfile_folder", "FOVs"],
                       validate="many_to_one")


def label_table(cells):
    counts = cells.groupby(["Region", "cell_type"]).size().unstack(fill_value=0)
    fractions = counts.div(counts.sum(axis=1), axis=0)
    named = [c for c in counts if not sibling.is_denovo(c)]
    letters = sorted((c for c in counts if sibling.is_denovo(c)), key=lambda c: (len(c), c))
    kept_named = [c for c in named if (fractions[c] >= MIN_NAMED_FRACTION).any()]
    other = [c for c in named if c not in kept_named]
    shown = fractions[letters + sorted(kept_named)].copy()
    shown[OTHER_NAMED_BAR] = fractions[other].sum(axis=1)
    return counts, shown, letters


def plot_regions(counts, shown, letters, donor, path):
    regions = [r for r in REGION_ORDER if r in shown.index]
    fig, axes = plt.subplots(len(regions), 1, sharex=True, sharey=True,
                             figsize=(PANEL_WIDTH_IN, len(regions) * PANEL_HEIGHT_IN))
    colors = [HIGHLIGHT_COLOR if c == HIGHLIGHT_LETTER else LETTER_COLOR if c in letters else
              OTHER_COLOR if c == OTHER_NAMED_BAR else NAMED_COLOR for c in shown.columns]
    for ax, region in zip(list(axes), regions):
        ax.bar(range(shown.shape[1]), shown.loc[region], color=colors)
        ax.set_title(f"Donor {donor} - {region}  (n={counts.loc[region].sum():,})", fontsize=11)
        ax.set_ylabel("fraction of region's cells")
    axes[-1].set_xticks(range(shown.shape[1]))
    axes[-1].set_xticklabels(shown.columns, rotation=90, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=PNG_DPI)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typing-h5", required=True, type=Path)
    parser.add_argument("--anchor-cells-csv", required=True, type=Path)
    parser.add_argument("--fov-annotations-csv", required=True, type=Path)
    parser.add_argument("--crosswalk-csv", required=True, type=Path)
    parser.add_argument("--donor", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    labelled = sibling.load_labels(args.typing_h5, args.anchor_cells_csv)
    cells = donor_cells(labelled, args.fov_annotations_csv, args.crosswalk_csv, args.donor)
    print(cells.groupby(["Region", "Block", "slide_id"]).size())
    counts, shown, letters = label_table(cells)
    stem = f"donor_{args.donor}_by_region"
    counts.to_csv(args.output_dir / f"{stem}_counts.csv")
    plot_regions(counts, shown, letters, args.donor, args.output_dir / stem)


if __name__ == "__main__":
    main()
