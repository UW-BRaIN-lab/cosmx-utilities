#!/usr/bin/env python3
"""Where down the tree did InSituTree stop being confident?

InSituTree descends a hierarchy, choosing among a node's children at each step. It has no
de novo escape hatch, so when a node's children are indistinguishable at CosMx depth the
cells do not drop out -- they pile into whichever child wins a near-tie, and the result looks
like a confident call on an implausibly large population.

The number that separates a real call from a coin flip is the posterior AT THE STEP WHERE THE
CHOICE WAS MADE, read against chance for that node. A node with 5 children has chance 0.20;
a leaf holding 47,000 cells at a median posterior of 0.22 did not win, it was drawn out of a
hat. The same leaf at 0.85 is a finding.

So this reports, per final label: how many cells, how many siblings it was competing against,
chance for that contest, and the posterior distribution of the cells that landed there -- plus
the posterior at level 1, to separate "the branch call was wrong" from "the branch call was
right and only the leaf choice collapsed".

Reads the summaryAnnotation CSV written by R/insitutree_typing.R (cell_id, annotLevel_1..N,
probs_annotLevel_1..N) and the hierarchy JSON the run used.

Emits:
  resolution_by_label.csv    one row per final label, sorted by cells
  resolution_by_node.csv     one row per competition node: how decisive its contest was,
                             and how lopsided the outcome

Usage:
    uv run python pipeline/python/diagnose_insitutree_resolution.py \\
        --summary-csv insitutree_summary.csv \\
        --hierarchy pipeline/reference/retina_hierarchy.json \\
        --out-dir insitutree_resolution
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

LEVEL_RE = re.compile(r"^annotLevel_(\d+)$")
PROB_PREFIX = "probs_"
CELL_ID_COL = "cell_id"
# A contest is "at chance" when the winner's posterior is no better than this multiple of
# 1/n_siblings. 1.5x is deliberately generous: anything under it is not a call.
CHANCE_FACTOR = 1.5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary-csv", type=Path, required=True,
                   help="insitutree_summary.csv (summaryAnnotation, one row per cell).")
    p.add_argument("--hierarchy", type=Path, required=True,
                   help="The hierarchy JSON the run used.")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for the CSVs.")
    p.add_argument("--top", type=int, default=30,
                   help="Labels to print (all are written to CSV).")
    return p.parse_args()


def terminal_nodes(hierarchy: dict) -> dict[str, list[str]]:
    """{node path -> its leaf children}. These are the contests that pick a final label."""
    nodes: dict[str, list[str]] = {}

    def walk(node, path: str) -> None:
        if isinstance(node, dict):
            for name, child in node.items():
                walk(child, f"{path}/{name}" if path else name)
        elif isinstance(node, list):
            nodes[path] = list(node)

    walk(hierarchy, "")
    return nodes


def main() -> None:
    args = parse_args()

    summary = pd.read_csv(args.summary_csv)
    if CELL_ID_COL not in summary.columns:
        sys.exit(f"ERROR: {args.summary_csv} has no '{CELL_ID_COL}' column")

    levels = sorted(int(m.group(1)) for c in summary.columns
                    if (m := LEVEL_RE.match(c)))
    if not levels:
        sys.exit(f"ERROR: no annotLevel_* columns in {args.summary_csv}; "
                 f"found {sorted(summary.columns)[:20]}")
    finest = max(levels)
    label_col = f"annotLevel_{finest}"
    prob_cols = {lvl: f"{PROB_PREFIX}annotLevel_{lvl}" for lvl in levels}
    missing = [c for c in prob_cols.values() if c not in summary.columns]
    if missing:
        sys.exit(f"ERROR: {args.summary_csv} is missing paired probability column(s) {missing}")
    print(f"{len(summary):,} cells, levels {levels}, finest = {label_col}")

    nodes = terminal_nodes(json.loads(args.hierarchy.read_text()))
    siblings_of: dict[str, tuple[str, int]] = {}
    for path, leaves in nodes.items():
        for leaf in leaves:
            siblings_of[leaf] = (path, len(leaves))
    print(f"hierarchy: {len(nodes)} terminal nodes over {len(siblings_of)} leaves")

    labels = summary[label_col].astype(str)
    final_prob = summary[prob_cols[finest]].astype(float)
    top_prob = summary[prob_cols[min(levels)]].astype(float)

    rows = []
    for label, index in labels.groupby(labels).groups.items():
        node, n_siblings = siblings_of.get(label, ("(not a leaf)", 0))
        probs = final_prob.loc[index].dropna()
        chance = 1.0 / n_siblings if n_siblings else np.nan
        median = float(probs.median()) if len(probs) else np.nan
        rows.append({
            "label": label,
            "cells": len(index),
            "pct_of_cohort": round(len(index) / len(summary) * 100, 2),
            "node": node,
            "n_siblings": n_siblings,
            "chance": round(chance, 3) if n_siblings else np.nan,
            "median_prob": round(median, 3) if len(probs) else np.nan,
            "q25_prob": round(float(probs.quantile(0.25)), 3) if len(probs) else np.nan,
            "q75_prob": round(float(probs.quantile(0.75)), 3) if len(probs) else np.nan,
            "median_prob_level1": round(float(top_prob.loc[index].dropna().median()), 3)
                                  if top_prob.loc[index].notna().any() else np.nan,
            # The headline: how far above chance the typical cell actually sat.
            "x_chance": round(median / chance, 2) if n_siblings and median == median else np.nan,
            "at_chance": bool(n_siblings and median == median
                              and median <= chance * CHANCE_FACTOR),
        })

    by_label = pd.DataFrame(rows).sort_values("cells", ascending=False)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    by_label.to_csv(args.out_dir / "resolution_by_label.csv", index=False)

    # Per contest: how lopsided was it, and was the winner above chance?
    node_rows = []
    for path, leaves in nodes.items():
        sub = by_label[by_label["node"] == path]
        total = int(sub["cells"].sum())
        if not total:
            continue
        winner = sub.iloc[0]
        node_rows.append({
            "node": path,
            "n_leaves": len(leaves),
            "cells": total,
            "chance": round(1.0 / len(leaves), 3),
            "winner": winner["label"],
            "winner_share_pct": round(winner["cells"] / total * 100, 1),
            "winner_median_prob": winner["median_prob"],
            "winner_x_chance": winner["x_chance"],
            "leaves_under_1pct": int((sub["cells"] / total < 0.01).sum()),
            "at_chance": bool(winner["at_chance"]),
        })
    by_node = (pd.DataFrame(node_rows)
               .sort_values("winner_share_pct", ascending=False))
    by_node.to_csv(args.out_dir / "resolution_by_node.csv", index=False)

    pd.set_option("display.width", 250)
    print(f"\n=== terminal nodes, most lopsided first ===")
    print(by_node.to_string(index=False))
    print(f"\n=== top {args.top} labels by size ===")
    print(by_label.head(args.top).to_string(index=False))

    flagged = by_node[by_node["at_chance"]]
    print(f"\n{len(flagged)} of {len(by_node)} contests were AT CHANCE "
          f"(winner's median posterior <= {CHANCE_FACTOR}x 1/n_leaves):")
    for _, r in flagged.iterrows():
        print(f"  {r['node']:42s} winner {r['winner']} took {r['winner_share_pct']}% "
              f"at p={r['winner_median_prob']} vs chance {r['chance']}")
    print("\nA contest at chance did not resolve. Report those cells at the NODE, not the "
          "leaf — or raise QUANTILE_ABS so InSituTree stops there itself.")


if __name__ == "__main__":
    main()
