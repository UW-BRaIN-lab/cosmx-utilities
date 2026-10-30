"""Per-slide cell-type abundance histograms for the k=27 anchor cohort.

De novo letters (lowercase-only labels) get one bar each; every named
reference type is lumped into one "named" bar. Slides are ordered by the
fraction of cells labelled `b`, so b-dominated slides cluster together.
"""
import argparse
import math
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

DOMINANT_LETTER = "b"
NAMED_BAR = "named"
SLIDES_PER_PAGE = 10
PANELS_PER_ROW = 2
PANEL_WIDTH_IN, PANEL_HEIGHT_IN = 7.0, 3.0
PNG_DPI = 200
HIGHLIGHT_COLOR, LETTER_COLOR, NAMED_COLOR = "#c0392b", "#4a6fa5", "#9a9a9a"


def is_denovo(label):
    return label.isalpha() and label.islower()


def load_labels(typing_h5, anchor_cells_csv):
    with h5py.File(typing_h5) as handle:
        labels = pd.DataFrame({"cell_id": handle["cell_id"][:].astype(str),
                               "cell_type": handle["cell_type"][:].astype(str)})
    slides = pd.read_csv(anchor_cells_csv, usecols=["cell_id", "slide_id"])
    return labels.merge(slides, on="cell_id", validate="one_to_one")


def fraction_table(cells):
    cells["bar"] = cells.cell_type.where(cells.cell_type.map(is_denovo), NAMED_BAR)
    counts = cells.groupby(["slide_id", "bar"]).size().unstack(fill_value=0)
    letters = sorted((c for c in counts if c != NAMED_BAR), key=lambda c: (len(c), c))
    counts = counts[letters + [NAMED_BAR]]
    return counts, counts.div(counts.sum(axis=1), axis=0)


def plot_page(counts, fractions, slides, path):
    rows = math.ceil(len(slides) / PANELS_PER_ROW)
    fig, axes = plt.subplots(rows, PANELS_PER_ROW, sharey=True, squeeze=False,
                             figsize=(PANELS_PER_ROW * PANEL_WIDTH_IN, rows * PANEL_HEIGHT_IN))
    colors = [HIGHLIGHT_COLOR if bar == DOMINANT_LETTER else
              NAMED_COLOR if bar == NAMED_BAR else LETTER_COLOR
              for bar in fractions.columns]
    for ax, slide in zip(axes.flat, slides):
        ax.bar(range(fractions.shape[1]), fractions.loc[slide], color=colors)
        ax.set_xticks(range(fractions.shape[1]))
        ax.set_xticklabels(fractions.columns, fontsize=8, rotation=90)
        ax.set_title(f"{slide}  (n={counts.loc[slide].sum():,}, "
                     f"b={fractions.loc[slide, DOMINANT_LETTER]:.0%})", fontsize=11)
        ax.set_ylim(0, 1)
    for ax in list(axes.flat)[len(slides):]:
        ax.axis("off")
    fig.supylabel("fraction of slide's anchor cells")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=PNG_DPI)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_pages(counts, fractions, order, output_dir):
    pages = [order[i:i + SLIDES_PER_PAGE] for i in range(0, len(order), SLIDES_PER_PAGE)]
    for number, slides in enumerate(pages, start=1):
        plot_page(counts, fractions, slides, output_dir / f"per_slide_histograms_p{number:02d}")


def plot_stacked(fractions, order, path):
    without_b = fractions.drop(columns=DOMINANT_LETTER)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
    top.bar(range(len(order)), fractions.loc[order, DOMINANT_LETTER], color=HIGHLIGHT_COLOR)
    top.set_ylabel(f"fraction `{DOMINANT_LETTER}`")
    renormalised = without_b.div(without_b.sum(axis=1), axis=0).loc[order]
    cmap = plt.get_cmap("tab20")
    bottom_edge = pd.Series(0.0, index=order)
    for i, bar in enumerate(renormalised.columns):
        color = NAMED_COLOR if bar == NAMED_BAR else cmap(i % 20)
        bottom.bar(range(len(order)), renormalised[bar], bottom=bottom_edge, color=color, label=bar)
        bottom_edge += renormalised[bar]
    bottom.set_ylabel(f"composition excluding `{DOMINANT_LETTER}`")
    bottom.set_xticks(range(len(order)))
    bottom.set_xticklabels(order, rotation=90, fontsize=6)
    bottom.legend(ncol=2, fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=PNG_DPI)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typing-h5", required=True, type=Path)
    parser.add_argument("--anchor-cells-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cells = load_labels(args.typing_h5, args.anchor_cells_csv)
    counts, fractions = fraction_table(cells)
    order = fractions[DOMINANT_LETTER].sort_values(ascending=False).index.tolist()
    counts.loc[order].to_csv(args.output_dir / "slide_celltype_counts.csv")
    fractions.loc[order].to_csv(args.output_dir / "slide_celltype_fractions.csv")
    plot_pages(counts, fractions, order, args.output_dir)
    plot_stacked(fractions, order, args.output_dir / "composition_overview.png")
    print(fractions[DOMINANT_LETTER].describe())


if __name__ == "__main__":
    main()
