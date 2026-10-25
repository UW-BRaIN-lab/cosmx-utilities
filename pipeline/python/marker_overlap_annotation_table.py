#!/usr/bin/env python3
"""Turn marker_overlap.csv into one evidence row per cluster, next to the current annotation.

Reads the long (cluster x reference type) table from marker_overlap_vs_reference.py and writes one row
per result type: its best reference match by concordance, how strong the support is, the reference
types that are nearly as good, and the best match by marker overlap. With --annotations-csv (e.g.
reference/denovo_annotations/stage4_k21.csv) the existing annotation and keep/new_name are joined on
for a side-by-side read. This is evidence for a human to check (against Napari), not a rename: it
never writes keep/new_name, so it cannot change what InSituTree is given.

support:
  strong   concordance >= 0.95 and the overlap p-value < 1e-3: the cluster's markers are that type's markers
  partial  the best-concordance type's overlap p-value < 0.05, but not strong
  none     no significant marker overlap with the best-concordance type: no SEA-AD name is proposed
Named result types (the Brain_6k reference types) are kept as controls: if they do not map to the right
reference type, the de novo rows mean little.

Usage:
    marker_overlap_annotation_table.py --marker-overlap-csv marker_overlap.csv \\
        --annotations-csv pipeline/reference/denovo_annotations/stage4_k21.csv --output table.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

STRONG_CONCORDANCE = 0.95
STRONG_P = 1e-3
PARTIAL_P = 0.05
ALSO_CLOSE_CONCORDANCE = 0.06    # other types this close to the best, and significant, are listed

STRONG, PARTIAL, NONE = "strong", "partial", "none"
REQUIRED = {"cluster", "kind", "n_cells", "n_markers_used", "reference_type", "concordance",
            "overlap", "overlap_p", "rank"}


def support_level(concordance: float, p: float) -> str:
    if concordance >= STRONG_CONCORDANCE and p < STRONG_P:
        return STRONG
    return PARTIAL if p < PARTIAL_P else NONE


def evidence_row(cluster: str, grp: pd.DataFrame, total_cells: int) -> dict:
    grp = grp.sort_values("rank")
    best = grp.iloc[0]
    close = grp.iloc[1:][(grp.iloc[1:]["concordance"] >= best["concordance"] - ALSO_CLOSE_CONCORDANCE)
                         & (grp.iloc[1:]["overlap_p"] < STRONG_P)]
    by_overlap = grp.sort_values(["overlap_p", "rank"]).iloc[0]
    level = support_level(float(best["concordance"]), float(best["overlap_p"]))
    return {
        "cluster": cluster, "kind": best["kind"], "n_cells": int(best["n_cells"]),
        "pct_cohort": round(100 * int(best["n_cells"]) / total_cells, 2),
        "support": level,
        "proposed_seaad_name": best["reference_type"] if level != NONE else "",
        "best_concordance": best["concordance"],
        "best_overlap": f"{int(best['overlap'])}/{int(best['n_markers_used'])}",
        "best_overlap_p": f"{best['overlap_p']:.1e}",
        "also_close": "; ".join(close["reference_type"]),
        "top_overlap_type": by_overlap["reference_type"],
        "top_overlap": f"{int(by_overlap['overlap'])}/{int(by_overlap['n_markers_used'])}",
        "top_overlap_p": f"{by_overlap['overlap_p']:.1e}",
    }


def build_table(long: pd.DataFrame, annotations: pd.DataFrame | None) -> pd.DataFrame:
    missing = REQUIRED - set(long.columns)
    if missing:
        sys.exit(f"ERROR: marker overlap table lacks columns {sorted(missing)}")
    sizes = long.drop_duplicates("cluster").set_index("cluster")["n_cells"]
    total = int(sizes.sum())
    table = pd.DataFrame([evidence_row(c, g, total) for c, g in long.groupby("cluster", sort=False)])
    table = table.sort_values(["kind", "n_cells"], ascending=[True, False]).reset_index(drop=True)
    if annotations is not None:
        ann = annotations.set_index("denovo_label")
        table["current_annotation"] = table["cluster"].map(ann["annotation"]).fillna("")
        table["current_new_name"] = table["cluster"].map(ann["new_name"]).fillna("")
        table["current_keep"] = table["cluster"].map(ann["keep"]).fillna("")
    return table


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--marker-overlap-csv", type=Path, required=True)
    p.add_argument("--annotations-csv", type=Path, default=None,
                   help="denovo_annotations CSV (denovo_label,keep,new_name,annotation,...).")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    long = pd.read_csv(args.marker_overlap_csv)
    annotations = pd.read_csv(args.annotations_csv, dtype=str) if args.annotations_csv else None
    table = build_table(long, annotations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    denovo = table[table["kind"] == "de novo"]
    total = int(table["n_cells"].sum())
    for level in (STRONG, PARTIAL, NONE):
        sel = denovo[denovo["support"] == level]
        print(f"de novo {level:<8} {len(sel):>3} clusters  {100 * sel['n_cells'].sum() / total:5.1f}% of cells")
    print(f"Wrote {args.output} ({len(table)} rows)")


if __name__ == "__main__":
    main()
