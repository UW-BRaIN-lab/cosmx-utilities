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

With --flat-typing-h5 it also compares each leaf against the FLAT InSituType run, which is
the test that actually works here. The posterior does not: a multinomial likelihood over a
few hundred counts and thousands of genes saturates, so every call comes back at ~1.000
whether it is right or wrong. What does discriminate is how a leaf fared when competition
was open. In the retina run the five cortical excitatory profiles won 244 cells between them
against 83 competitors; the tree then handed one of them 46,975. A profile that cannot win a
fair contest should not be trusted with a captive one.

Emits:
  resolution_by_label.csv    one row per final label, sorted by cells
  resolution_by_node.csv     one row per competition node: how decisive its contest was,
                             and how lopsided the outcome
  flat_vs_tree.csv           (--flat-typing-h5) per leaf and per node: cells won in the flat
                             run against cells assigned by the tree, and the ratio

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

import h5py
import numpy as np
import pandas as pd

LEVEL_RE = re.compile(r"^annotLevel_(\d+)$")
PROB_PREFIX = "probs_"
CELL_ID_COL = "cell_id"
# A contest is "at chance" when the winner's posterior is no better than this multiple of
# 1/n_siblings. 1.5x is deliberately generous: anything under it is not a call.
CHANCE_FACTOR = 1.5
# A node whose tree total dwarfs what its leaves won in open competition inherited those
# cells from dissolved de novo clusters; the finer the split it then reports, the less of it
# rests on any profile that ever demonstrated it could fit this data.
INHERITED_FOLD_FLAG = 5.0
DENOVO_KEEP_COL, DENOVO_LABEL_COL, DENOVO_NAME_COL = "keep", "denovo_label", "new_name"


def load_flat_counts(typing_h5: Path) -> pd.Series:
    """Cells per label in the flat InSituType run (named types AND de novo letters)."""
    with h5py.File(typing_h5, "r") as f:
        if "cell_type" not in f:
            sys.exit(f"ERROR: {typing_h5} has no /cell_type dataset")
        values = np.asarray(f["cell_type"][()])
    if values.dtype.kind == "S" or (values.size and isinstance(values.flat[0], bytes)):
        values = np.array([v.decode("utf-8") for v in values.ravel()])
    return pd.Series(values.astype(str)).value_counts()


def denovo_rename(annotations_csv: Path) -> dict[str, str]:
    """{tree leaf name -> the de novo letter it was built from}, for the kept clusters.

    Without this a kept de novo leaf looks like it won nothing in the flat run, because the
    flat run called those cells 'm', not 'denovo_astrocyte_1' -- and it would be flagged as
    a captive-audience leaf when it is the opposite: a profile fit to this very data.
    """
    ann = pd.read_csv(annotations_csv, dtype=str).fillna("")
    for col in (DENOVO_KEEP_COL, DENOVO_LABEL_COL, DENOVO_NAME_COL):
        if col not in ann.columns:
            sys.exit(f"ERROR: {annotations_csv} has no '{col}' column")
    kept = ann[ann[DENOVO_KEEP_COL].str.strip().isin({"1", "true", "True", "yes"})]
    return {r[DENOVO_NAME_COL].strip(): r[DENOVO_LABEL_COL].strip()
            for _, r in kept.iterrows() if r[DENOVO_NAME_COL].strip()}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary-csv", type=Path, required=True,
                   help="insitutree_summary.csv (summaryAnnotation, one row per cell).")
    p.add_argument("--hierarchy", type=Path, required=True,
                   help="The hierarchy JSON the run used.")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for the CSVs.")
    p.add_argument("--flat-typing-h5", type=Path, default=None,
                   help="The FLAT InSituType result (insitutype_result.h5). Adds the "
                        "flat-run count per leaf -- the diagnostic that works here, since "
                        "the posterior saturates at 1.0 and carries no information.")
    p.add_argument("--denovo-annotations", type=Path, default=None,
                   help="The keep/new_name CSV the profile matrix was built from, so kept "
                        "de novo leaves are matched to the letter they came from.")
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

    if args.flat_typing_h5 is not None:
        flat_counts = load_flat_counts(args.flat_typing_h5)
        source_label = (denovo_rename(args.denovo_annotations)
                        if args.denovo_annotations else {})
        print(f"\nFlat run: {int(flat_counts.sum()):,} cells, {len(flat_counts)} labels"
              + (f"; {len(source_label)} kept de novo leaves mapped to their letters"
                 if source_label else ""))

        tree_by_label = by_label.set_index("label")["cells"]
        rows = []
        for leaf, (node, _) in siblings_of.items():
            flat_name = source_label.get(leaf, leaf)
            flat = int(flat_counts.get(flat_name, 0))
            tree = int(tree_by_label.get(leaf, 0))
            rows.append({"leaf": leaf, "node": node, "flat_label": flat_name,
                         "flat_cells": flat, "tree_cells": tree,
                         "fold": round(tree / flat, 1) if flat else np.inf})
        flat_vs_tree = pd.DataFrame(rows).sort_values("tree_cells", ascending=False)

        # Per node: what its leaves won on their own, against what the tree gave them.
        per_node = (flat_vs_tree.groupby("node")[["flat_cells", "tree_cells"]].sum()
                    .rename(columns={"flat_cells": "node_flat_cells",
                                     "tree_cells": "node_tree_cells"}))
        # A cell too ambiguous to pick a leaf stops AT the node and carries the node's own
        # name, so it belongs to no leaf. Counting only leaves made single-child nodes read
        # as empty -- Mueller_glia showed 0 while holding 19,717 cells.
        stopped_at_node = {node: int(tree_by_label.get(node.split("/")[-1], 0))
                           for node in per_node.index}
        per_node["stopped_at_node"] = pd.Series(stopped_at_node)
        per_node["node_tree_cells"] += per_node["stopped_at_node"]
        per_node["node_fold"] = (per_node["node_tree_cells"]
                                 / per_node["node_flat_cells"].replace(0, np.nan)).round(1)
        per_node["inherited"] = per_node["node_fold"] >= INHERITED_FOLD_FLAG
        flat_vs_tree = flat_vs_tree.merge(per_node, on="node", how="left")
        flat_vs_tree.to_csv(args.out_dir / "flat_vs_tree.csv", index=False)

        print("\n=== what each contest's leaves won on their own, vs what the tree gave them ===")
        print(per_node.sort_values("node_fold", ascending=False).to_string())

        inherited = per_node[per_node["inherited"]].sort_values("node_fold", ascending=False)
        print(f"\n{len(inherited)} of {len(per_node)} contests inherited most of their cells "
              f"(tree total >= {INHERITED_FOLD_FLAG}x what their leaves won openly):")
        for node, r in inherited.iterrows():
            leaves = flat_vs_tree[flat_vs_tree["node"] == node].nlargest(1, "tree_cells")
            winner = leaves.iloc[0]
            print(f"  {node:42s} {int(r['node_flat_cells']):>7,} -> "
                  f"{int(r['node_tree_cells']):>7,} ({r['node_fold']}x); "
                  f"largest leaf {winner['leaf']} won {winner['flat_cells']:,} openly")
        print("\nThose contests are decided among profiles that never showed they fit this "
              "data. Report them at the NODE, or give the de novo cluster its own leaf.")

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
