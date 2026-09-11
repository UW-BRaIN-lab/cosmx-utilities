#!/usr/bin/env python3
"""Head-to-head: the de-novo-removed COUNTERFACTUAL against an independent SUPERVISED run.

Two different questions get asked about the same cells, and the PI wants both:

  A. COUNTERFACTUAL — "if the 27 de-novo letters had not been available to choose, which GBmap
     type would each cell have picked?" Answered by holding the k=27 semi-supervised fit fixed
     and taking the argmax over the 54 named columns of its own stored log-likelihood table
     (75c's forced_named_posteriors.csv). The profiles are the ones that fit converged on.

  B. FRESH SUPERVISED — "what would a supervised InSituType run against GBmap call this cell?"
     Answered by insitutype(n_clusts = 0) on the same cells from the published GBmap reference,
     with no de-novo clusters at any point. Because update_reference_profiles and rescale still
     run, the named profiles adapt to ALL 2.54M cells rather than only the ~6.7% a semi-supervised
     fit leaves named — so B's profiles genuinely differ from A's, and the two can disagree.

Neither is more correct; they answer different things. Where they AGREE, the letter's cells have
a GBmap identity robust to how you ask. Where they DIVERGE, the answer depends on whether the
de-novo clusters were there while the profiles were being fitted — which is itself worth knowing
before those letters are named.

Reports overall agreement, agreement per semi-supervised label, and the cross-tab of A against B.

Inputs:
  --typing-h5   anchor_typing.h5 from the k=27 fit — the semi-supervised label per cell.
  --source-a    A's per-cell call: CSV (cell_id, top1_type) or an InSituType result h5.
  --source-b    B's per-cell call, same two formats.
  --label-a / --label-b   names for the two columns in the output (default counterfactual/supervised).
Outputs:
  --output-dir  agreement_by_label.csv, a_vs_b_crosstab.csv; summaries print to stdout.

Usage:
    uv run python pipeline/python/compare_forced_call_sources.py \\
        --typing-h5 anchor_typing.h5 \\
        --source-a forced_named_posteriors.csv \\
        --source-b supervised_anchor_typing.h5 \\
        --annotations pipeline/reference/denovo_annotations/fullcohort_pruned_k27.csv \\
        --output-dir forced_vs_supervised/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from anchor_profiles import is_denovo, read_cell_calls

TOP_DESTINATIONS = 2


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typing-h5", type=Path, required=True)
    p.add_argument("--source-a", type=Path, required=True,
                   help="CSV (cell_id, top1_type) or InSituType result h5.")
    p.add_argument("--source-b", type=Path, required=True)
    p.add_argument("--label-a", default="counterfactual")
    p.add_argument("--label-b", default="supervised")
    p.add_argument("--annotations", type=Path,
                   help="De-novo annotation CSV, for readable letter labels.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_calls(path: Path) -> pd.Series:
    """Per-cell call from either format, as a Series indexed by cell_id."""
    if path.suffix.lower() in (".h5", ".hdf5"):
        df = read_cell_calls(path)
        return df.set_index("cell_id")["cell_type"].astype(str)
    df = pd.read_csv(path, usecols=["cell_id", "top1_type"])
    return df.set_index("cell_id")["top1_type"].astype(str)


def main() -> None:
    args = parse_args()

    semisup = read_cell_calls(args.typing_h5).set_index("cell_id")["cell_type"].astype(str)
    a = read_calls(args.source_a)
    b = read_calls(args.source_b)
    print(f"{args.label_a}: {len(a):,} cells, {a.nunique()} types")
    print(f"{args.label_b}: {len(b):,} cells, {b.nunique()} types")

    df = pd.DataFrame({"semisup": semisup, "a": a, "b": b}).dropna()
    if df.empty:
        sys.exit("ERROR: no cell_id overlap across the three inputs.")
    dropped = len(semisup) - len(df)
    if dropped:
        print(f"WARNING: {dropped:,} cells missing from one of the sources; excluded.")
    print(f"Compared on {len(df):,} cells\n")

    agree = df["a"] == df["b"]
    print(f"=== {args.label_a} vs {args.label_b} ===")
    print(f"  overall agreement: {100 * agree.mean():.1f}%")
    denovo_mask = df["semisup"].map(is_denovo)
    for name, mask in (("cells the fit left NAMED", ~denovo_mask),
                       ("cells the fit put in a LETTER", denovo_mask)):
        if mask.any():
            print(f"  {name:32s} {100 * agree[mask].mean():5.1f}%  ({int(mask.sum()):,} cells)")
    print("\n  Where they agree, the GBmap identity is robust to how the question is asked.")
    print("  Where they diverge, it depends on whether the de-novo clusters were present")
    print("  while the profiles were being fitted.\n")

    display = {}
    if args.annotations:
        ann = pd.read_csv(args.annotations).dropna(subset=["annotation"])
        display = dict(zip(ann["denovo_label"].astype(str).str.strip(),
                           ann["annotation"].astype(str).str.strip()))

    rows = []
    for label, grp in df.groupby("semisup", observed=True):
        top_a = grp["a"].value_counts(normalize=True)
        top_b = grp["b"].value_counts(normalize=True)
        row = {"semisup_label": label, "display_label": display.get(label, label),
               "is_denovo": is_denovo(label), "n_cells": len(grp),
               "agreement_pct": round(100 * (grp["a"] == grp["b"]).mean(), 2)}
        for rank in range(TOP_DESTINATIONS):
            for tag, top in ((args.label_a, top_a), (args.label_b, top_b)):
                row[f"{tag}_{rank + 1}"] = top.index[rank] if rank < len(top) else ""
                row[f"{tag}_{rank + 1}_pct"] = (round(100 * top.iloc[rank], 2)
                                                if rank < len(top) else float("nan"))
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values(["is_denovo", "n_cells"],
                                             ascending=[False, False])

    print("Per de-novo letter — the top call under each question:")
    dn = summary[summary["is_denovo"]]
    cols = ["display_label", "n_cells", f"{args.label_a}_1", f"{args.label_a}_1_pct",
            f"{args.label_b}_1", f"{args.label_b}_1_pct", "agreement_pct"]
    print(dn[cols].to_string(index=False))

    crosstab = pd.crosstab(df["a"], df["b"])
    crosstab.index.name = args.label_a
    crosstab.columns.name = args.label_b

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "agreement_by_label.csv", index=False)
    crosstab.to_csv(args.output_dir / "a_vs_b_crosstab.csv")
    print(f"\nWrote {args.output_dir}/agreement_by_label.csv ({len(summary)} labels)")
    print(f"Wrote {args.output_dir}/a_vs_b_crosstab.csv "
          f"({crosstab.shape[0]} x {crosstab.shape[1]})")


if __name__ == "__main__":
    main()
