#!/usr/bin/env python3
"""Tests for marker_overlap_annotation_table.py: support levels, also_close, the overlap-best column,
the annotation join, and that it never emits keep/new_name for the clusters it summarises.

Run:  uv run --with pandas --with numpy python pipeline/python/tests/test_marker_overlap_annotation_table.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "marker_overlap_annotation_table.py"

# (cluster, kind, n_cells, [(reference_type, concordance, overlap, p), ...] best first by concordance)
CLUSTERS = [
    ("o", "de novo", 600, [("R_oligo", 0.99, 11, 2e-21), ("R_opc", 0.73, 1, 0.11)]),
    ("p", "de novo", 100, [("R_L6b", 0.99, 11, 2e-21), ("R_L6CT", 0.93, 12, 4e-24), ("R_NP", 0.51, 3, 2e-4)]),
    ("j", "de novo", 150, [("R_L23", 0.85, 4, 4e-6), ("R_L4", 0.84, 5, 7e-8)]),
    ("b", "de novo", 80, [("R_L4", 0.64, 1, 0.11), ("R_L23", 0.64, 3, 2e-4)]),
    ("i", "de novo", 500, [("R_L4", 0.72, 0, 1.0), ("R_L23", 0.71, 0, 1.0)]),
    ("Endothelial", "named", 70, [("R_endo", 0.95, 9, 2e-16), ("R_vlmc", 0.84, 3, 2e-4)]),
]
N_MARKERS = 20


def long_table() -> pd.DataFrame:
    rows = []
    for cluster, kind, n_cells, matches in CLUSTERS:
        for rank, (ref, conc, hits, p) in enumerate(matches, start=1):
            rows.append({"cluster": cluster, "kind": kind, "n_cells": n_cells, "n_markers_used": N_MARKERS,
                         "reference_type": ref, "concordance": conc, "overlap": hits, "overlap_p": p, "rank": rank})
    return pd.DataFrame(rows)


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        long_csv, ann_csv, out = tmp / "long.csv", tmp / "ann.csv", tmp / "table.csv"
        long_table().to_csv(long_csv, index=False)
        pd.DataFrame({"denovo_label": ["o", "i"], "keep": ["TRUE", "TRUE"],
                      "new_name": ["Oligodendrocyte_denovo", "Low_signal"],
                      "annotation": ["o - Oligodendrocyte", "i - Unresolved"]}).to_csv(ann_csv, index=False)

        cmd = [sys.executable, str(SCRIPT), "--marker-overlap-csv", str(long_csv),
               "--annotations-csv", str(ann_csv), "--output", str(out)]
        print("  $", " ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        print(res.stdout, res.stderr)
        assert res.returncode == 0
        t = pd.read_csv(out, keep_default_na=False).set_index("cluster")

        print("test_support_levels")
        assert t.loc["o", "support"] == "strong" and t.loc["o", "proposed_seaad_name"] == "R_oligo"
        assert t.loc["j", "support"] == "partial" and t.loc["j", "proposed_seaad_name"] == "R_L23"
        assert t.loc["i", "support"] == "none" and t.loc["i", "proposed_seaad_name"] == ""
        assert t.loc["b", "support"] == "none"          # best-by-concordance has p=0.11
        assert t.loc["Endothelial", "support"] == "strong" and t.loc["Endothelial", "kind"] == "named"

        print("test_also_close_lists_only_significant_near_ties")
        assert t.loc["p", "also_close"] == "R_L6CT", t.loc["p"]      # within 0.06 and p < 1e-3
        assert t.loc["o", "also_close"] == "" and t.loc["j", "also_close"] == "R_L4"

        print("test_top_overlap_type_catches_a_better_overlap_than_the_best_concordance")
        assert t.loc["b", "top_overlap_type"] == "R_L23" and t.loc["b", "top_overlap"] == "3/20"
        assert t.loc["p", "top_overlap_type"] == "R_L6CT"      # 12/20 beats the 11/20 of the concordance winner

        print("test_percentages_annotations_and_ordering")
        assert abs(t["pct_cohort"].sum() - 100) < 0.1
        assert t.loc["o", "current_new_name"] == "Oligodendrocyte_denovo" and t.loc["j", "current_annotation"] == ""
        assert list(t.index[:5]) == ["o", "i", "j", "p", "b"], list(t.index)   # de novo first, by size
        assert "keep" not in t.columns and "new_name" not in t.columns

        print("test_missing_columns_are_refused")
        long_table().drop(columns=["overlap_p"]).to_csv(tmp / "bad.csv", index=False)
        bad = subprocess.run([sys.executable, str(SCRIPT), "--marker-overlap-csv", str(tmp / "bad.csv"),
                              "--output", str(tmp / "x.csv")], capture_output=True, text=True)
        assert bad.returncode != 0 and "overlap_p" in (bad.stdout + bad.stderr)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
