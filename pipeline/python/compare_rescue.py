#!/usr/bin/env python3
"""Stage 5h (Phase 3): did one more rescue iteration collapse the Low_signal sink? (48 -> X)

Reads the per-cell rescue labels (rescue_lowsignal.R) and reports how much of the Low_signal
pool leaves the flat/sink state when profiles are re-derived from the pool and cells reassigned:

  - RESCUE RATE — fraction of Low_signal cells that snapped to an existing NAMED type vs formed
    a de-novo cluster vs landed back on the Low_signal-sink profile. A big drop => the pool was
    still reference/transfer-limited (keep iterating before any hybrid interpretation); barely
    moves => the irreducible core (the hybrid / continuum / admixture story is what is left).
  - BY COMPARTMENT — where the rescued cells go (malignant / neuronal / immune / ...), from the
    hierarchy leaf->compartment map. If they resolve into real lineage programs, the drop is
    meaningful; if they only tile into new de-novo clusters, it may just be finer manifold tiling
    (the Phase-2 top-two structure is the arbiter — flagged, not decided, here).
  - CNV CONCORDANCE — cross-tab the rescued label against the Phase-1 CNV malignant call: do the
    cells rescued into malignant types coincide with the CNV-malignant Low_signal cells?

Inputs:
  --rescue       rescue_lowsignal.csv (cell_id, rescue_label, rescue_prob, is_denovo).
  --hierarchy    insitutree_hierarchy.json (leaf -> compartment; de-novo/unknown -> "denovo").
  --cell-table   cell_cnv_table.csv.gz from Phase 1 (optional; for CNV concordance).
Writes (--output-dir): rescue_summary.txt, rescue_by_compartment.csv, rescue_rate.png,
  rescue_vs_cnv.csv.

Usage:
    python pipeline/python/compare_rescue.py --rescue rescue/rescue_lowsignal.csv \\
        --hierarchy pipeline/reference/insitutree_hierarchy.json \\
        --cell-table diagnostics/cell_cnv_table.csv.gz --output-dir rescue
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compare_insitucnv_groups import DEFAULT_MALIGNANT


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rescue", type=Path, required=True)
    p.add_argument("--hierarchy", type=Path, required=True)
    p.add_argument("--cell-table", type=Path, default=None)
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--prob-threshold", type=float, default=0.5,
                   help="Min rescue_prob to count an assignment as confident.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def leaf_compartments(hierarchy_path: Path) -> dict[str, str]:
    tree = json.loads(hierarchy_path.read_text())
    out: dict[str, str] = {}

    def collect(node) -> list[str]:
        leaves: list[str] = []
        if isinstance(node, dict):
            for k, v in node.items():
                child = collect(v)
                leaves += child or [k]
        elif isinstance(node, list):
            for v in node:
                leaves += collect(v)
        elif isinstance(node, str):
            leaves += [node]
        return leaves

    for compartment, sub in tree.items():
        for leaf in collect(sub) or [compartment]:
            out[leaf] = compartment
    return out


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comp = leaf_compartments(args.hierarchy)
    lowsignal_denovo = f"{args.lowsignal_label}_denovo"

    rescue = pd.read_csv(args.rescue)
    n = len(rescue)
    conf = rescue["rescue_prob"] >= args.prob_threshold
    is_denovo = rescue["is_denovo"].astype(bool)
    # a cell "leaves the sink" if it snapped, confidently, to a NON-sink named type OR a
    # (non-sink) de-novo cluster; still-sink = the Low_signal-sink profile or an unconfident call
    to_sink = rescue["rescue_label"].astype(str).isin([args.lowsignal_label, lowsignal_denovo])
    named_rescue = conf & ~is_denovo & ~to_sink
    denovo_rescue = conf & is_denovo & ~to_sink
    still_sink = ~(named_rescue | denovo_rescue)

    lines = [f"Phase 3 — one more rescue on the Low_signal pool ({n:,} cells)\n"]
    lines.append(f"  rescued to a NAMED type:   {int(named_rescue.sum()):>9,}  "
                 f"({named_rescue.mean():.1%})")
    lines.append(f"  rescued to a NEW de-novo:  {int(denovo_rescue.sum()):>9,}  "
                 f"({denovo_rescue.mean():.1%})")
    lines.append(f"  still flat / sink / low-conf:{int(still_sink.sum()):>8,}  "
                 f"({still_sink.mean():.1%})")
    left = named_rescue.mean() + denovo_rescue.mean()
    lines.append(f"\n  48 -> {(1 - left) * 100:.0f}%-equivalent Low_signal remaining "
                 f"({left:.1%} of the pool left the sink).")
    lines.append("  => " + ("SUBSTANTIAL collapse — the pool was still reference/transfer-"
                            "limited; iterate before any hybrid interpretation (but confirm the "
                            "rescued cells carry real programs, not just new manifold tiles — see "
                            "the by-compartment split + Phase-2 top-two structure)."
                            if left > 0.33 else
                            "BARELY moves — this looks like the irreducible core; the hybrid / "
                            "continuum / admixture story is what genuinely remains."))

    # --- by compartment ---------------------------------------------------------------
    rescued = rescue[named_rescue].copy()
    rescued["compartment"] = rescued["rescue_label"].map(comp).fillna("other")
    by_comp = (rescued["compartment"].value_counts()
               .rename_axis("compartment").reset_index(name="n_cells"))
    by_comp["frac_of_pool"] = by_comp["n_cells"] / n
    by_comp.to_csv(args.output_dir / "rescue_by_compartment.csv", index=False)
    if len(by_comp):
        lines.append("\n  named rescues by compartment (fraction of the whole pool):")
        for _, r in by_comp.iterrows():
            lines.append(f"    {r['compartment']:<14} {int(r['n_cells']):>9,}  "
                         f"({r['frac_of_pool']:.1%})")

    # --- CNV concordance (optional) ---------------------------------------------------
    if args.cell_table is not None and args.cell_table.exists():
        malignant_types = set(DEFAULT_MALIGNANT)
        cells = pd.read_csv(args.cell_table, index_col=0)
        j = rescue.set_index("cell_id").join(
            cells[["is_malignant_call"]], how="inner")
        j["rescued_malignant"] = (j["rescue_prob"] >= args.prob_threshold) & \
            j["rescue_label"].astype(str).isin(malignant_types)
        j["cnv_malignant"] = j["is_malignant_call"].astype(bool)
        ct = pd.crosstab(j["rescued_malignant"], j["cnv_malignant"])
        ct.to_csv(args.output_dir / "rescue_vs_cnv.csv")
        both = int(((j["rescued_malignant"]) & (j["cnv_malignant"])).sum())
        n_resc_mal = int(j["rescued_malignant"].sum())
        n_cnv_mal = int(j["cnv_malignant"].sum())
        lines.append("\n  CNV concordance (of Low_signal cells also in the CNV table):")
        lines.append(f"    rescued->malignant AND CNV-malignant: {both:,}")
        lines.append(f"    rescued->malignant: {n_resc_mal:,}; CNV-malignant: {n_cnv_mal:,}")
        if n_resc_mal:
            lines.append(f"    => {both / n_resc_mal:.1%} of expression-rescued-malignant cells "
                         "are independently CNV-malignant.")

    # --- plot: rescue-rate stacked bar ------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.5, 5))
    parts = [("named rescue", named_rescue.mean(), "#1a9850"),
             ("de-novo rescue", denovo_rescue.mean(), "#66bd63"),
             ("still flat/sink", still_sink.mean(), "#bdbdbd")]
    bottom = 0.0
    for label, frac, color in parts:
        ax.bar(0, frac, bottom=bottom, color=color, label=f"{label} ({frac:.0%})")
        bottom += frac
    ax.set_xlim(-1, 1.6); ax.set_xticks([])
    ax.set_ylabel("fraction of the Low_signal pool")
    ax.set_title("Phase 3: does one more rescue collapse the sink?")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output_dir / "rescue_rate.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    (args.output_dir / "rescue_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nWrote summary + tables + plot to {args.output_dir}")


if __name__ == "__main__":
    main()
