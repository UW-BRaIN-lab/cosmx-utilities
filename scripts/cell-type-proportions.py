#!/usr/bin/env python3
"""
Generate a grouped bar chart comparing cell type proportions between two
CosMx slides. Uses the UMAP CSVs extracted by umap-cell-types.py.
"""

import argparse

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ── Configuration ──────────────────────────────────────────────────────────────

SLIDE_1_CSV = "/tmp/seurat_slide1_umap.csv"
SLIDE_1_LABEL = "GBM"
SLIDE_1_COLOR = "#E64B35"

SLIDE_2_CSV = "/tmp/seurat_slide2_umap.csv"
SLIDE_2_LABEL = "Normal"
SLIDE_2_COLOR = "#4DBBD5"
SLIDE_2_FOVS = set(range(1, 49)) | set(range(199, 201)) | set(range(148, 198))

NEURONAL_TYPES = [
    "Inhibitory.neuron.A",
    "Inhibitory.neuron.B",
    "Inhibitory.neuron.C",
    "L2.3.neuron",
    "L4.neuron",
    "L6.neuron",
]

MICROGLIA_TYPES = [
    "Microglia.A",
    "Microglia.B_1",
    "Microglia.B_2",
]

CELL_TYPE_DISPLAY_NAMES = {
    "Inhibitory.neuron.A": "Inh A (PVALB+)",
    "Inhibitory.neuron.B": "Inh B (CCK+)",
    "Inhibitory.neuron.C": "Inh C (SST+)",
    "L2.3.neuron": "L2/3",
    "L4.neuron": "L4",
    "L6.neuron": "L5/6",
}

OUTPUT_PATH = "cell-type-proportions.png"


def parse_args():
    parser = argparse.ArgumentParser(
        description="CosMx cell type density bar chart"
    )
    parser.add_argument(
        "--neurons-microglia", action="store_true",
        help="Restrict x-axis to neuronal and microglial cell types only",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    df1 = pd.read_csv(SLIDE_1_CSV)
    df1 = df1[df1["cell_type"].notna() & (df1["cell_type"] != "")]

    df2 = pd.read_csv(SLIDE_2_CSV)
    df2 = df2[df2["cell_type"].notna() & (df2["cell_type"] != "")]
    df2 = df2[df2["fov"].isin(SLIDE_2_FOVS)]

    # Compute mean cells per FOV for each cell type
    n_fovs_1 = df1["fov"].nunique()
    n_fovs_2 = df2["fov"].nunique()

    counts1 = df1["cell_type"].value_counts()
    counts2 = df2["cell_type"].value_counts()
    density1 = counts1 / n_fovs_1
    density2 = counts2 / n_fovs_2

    # Union of cell types, sorted alphabetically
    all_types = sorted(set(density1.index) | set(density2.index))

    if args.neurons_microglia:
        keep = set(NEURONAL_TYPES) | set(MICROGLIA_TYPES)
        all_types = [ct for ct in all_types if ct in keep]
        print(f"Filtering to neurons and microglia: {all_types}")

    vals1 = [density1.get(ct, 0) for ct in all_types]
    vals2 = [density2.get(ct, 0) for ct in all_types]

    display_labels = [CELL_TYPE_DISPLAY_NAMES.get(ct, ct) for ct in all_types]

    # ── Plot ───────────────────────────────────────────────────────────────
    x = np.arange(len(all_types))
    bar_width = 0.35

    fig, ax = plt.subplots(figsize=(14, 6))

    bars2 = ax.bar(x - bar_width / 2, vals2, bar_width,
                   label=f"{SLIDE_2_LABEL} ({n_fovs_2} FOVs, {len(df2):,} cells)",
                   color=SLIDE_2_COLOR, edgecolor="white", linewidth=0.5)
    bars1 = ax.bar(x + bar_width / 2, vals1, bar_width,
                   label=f"{SLIDE_1_LABEL} ({n_fovs_1} FOVs, {len(df1):,} cells)",
                   color=SLIDE_1_COLOR, edgecolor="white", linewidth=0.5)

    # Add value labels above bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height >= 1:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, height + 0.5,
                    f"{height:.1f}",
                    ha="center", va="bottom", fontsize=7,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(display_labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Mean cells per FOV", fontsize=11)
    ax.set_title("Cell type density by slide", fontsize=13)
    ax.legend(fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    output = OUTPUT_PATH
    if args.neurons_microglia:
        output = output.replace(".png", "-neurons-microglia.png")
    plt.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved to {output}")


if __name__ == "__main__":
    main()
