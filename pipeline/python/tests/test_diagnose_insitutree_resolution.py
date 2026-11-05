#!/usr/bin/env python3
"""Tests for the InSituTree resolution diagnostic.

The fixture builds the two cases the diagnostic exists to tell apart:

  Resolved    a 2-leaf node whose winner sits at p=0.95 against chance 0.50
  At chance   a 5-leaf node whose winner takes ~everything at p=0.22 against chance 0.20

Both look identical in a count table -- one leaf holding nearly all the cells. Only the
posterior read against chance separates "a real population" from "a coin flip with a big
number attached".

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_diagnose_insitutree_resolution.py
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
_SCRIPT = _PY_DIR / "diagnose_insitutree_resolution.py"
_spec = importlib.util.spec_from_file_location("diagnose_insitutree_resolution", _SCRIPT)
dir_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dir_mod)

HIERARCHY = {
    "Glia": {"Astro": ["astro_a", "astro_b"]},
    "Cortical": {"Excitatory": [f"exc_{i}" for i in range(5)]},
}
RESOLVED_P, CHANCE_P = 0.95, 0.22
N_BIG, N_SMALL = 400, 5


def build(tmp: Path) -> tuple[Path, Path]:
    rows = []
    # Astro: a decisive 2-way contest.
    for i in range(N_BIG):
        rows.append({"cell_id": f"res{i}", "annotLevel_1": "Glia",
                     "annotLevel_2": "astro_a",
                     "probs_annotLevel_1": 0.99, "probs_annotLevel_2": RESOLVED_P})
    for i in range(N_BIG // 2):
        rows.append({"cell_id": f"res_b{i}", "annotLevel_1": "Glia",
                     "annotLevel_2": "astro_b",
                     "probs_annotLevel_1": 0.99, "probs_annotLevel_2": RESOLVED_P})
    # Excitatory: confident at level 1, a coin flip at the leaf, winner takes nearly all.
    for i in range(N_BIG):
        rows.append({"cell_id": f"ch{i}", "annotLevel_1": "Cortical",
                     "annotLevel_2": "exc_0",
                     "probs_annotLevel_1": 0.97, "probs_annotLevel_2": CHANCE_P})
    for leaf in range(1, 5):
        for i in range(N_SMALL):
            rows.append({"cell_id": f"ch_{leaf}_{i}", "annotLevel_1": "Cortical",
                         "annotLevel_2": f"exc_{leaf}",
                         "probs_annotLevel_1": 0.97, "probs_annotLevel_2": CHANCE_P})

    summary = tmp / "summary.csv"
    pd.DataFrame(rows).to_csv(summary, index=False)
    hierarchy = tmp / "hierarchy.json"
    hierarchy.write_text(json.dumps(HIERARCHY))
    return summary, hierarchy


def test_separates_resolved_from_at_chance(tmp: Path) -> None:
    summary, hierarchy = build(tmp)
    out = tmp / "out"
    result = subprocess.run(
        [sys.executable, str(_SCRIPT), "--summary-csv", str(summary),
         "--hierarchy", str(hierarchy), "--out-dir", str(out)],
        check=True, capture_output=True, text=True)

    by_label = pd.read_csv(out / "resolution_by_label.csv").set_index("label")
    assert by_label.loc["exc_0", "n_siblings"] == 5
    assert by_label.loc["exc_0", "chance"] == 0.2
    assert by_label.loc["exc_0", "at_chance"], "a 0.22 winner against chance 0.20 is a coin flip"
    assert not by_label.loc["astro_a", "at_chance"], "0.95 against chance 0.50 is a real call"
    # The branch call was fine in BOTH cases -- that is the distinction being drawn.
    assert by_label.loc["exc_0", "median_prob_level1"] == 0.97

    by_node = pd.read_csv(out / "resolution_by_node.csv").set_index("node")
    excitatory = by_node.loc["Cortical/Excitatory"]
    assert excitatory["winner"] == "exc_0"
    assert excitatory["winner_share_pct"] > 90
    assert bool(excitatory["at_chance"])
    assert not bool(by_node.loc["Glia/Astro", "at_chance"])
    assert "1 of 2 contests were AT CHANCE" in result.stdout
    print("  a lopsided contest at chance is flagged; a lopsided one above chance is not")


def test_size_alone_does_not_flag(tmp: Path) -> None:
    """A leaf can hold most of its branch and still be a real call -- size is not the test."""
    summary, hierarchy = build(tmp)
    frame = pd.read_csv(summary)
    # Make astro_a dominate as hard as exc_0 does, but keep its decisive posterior.
    frame = frame[~frame.annotLevel_2.eq("astro_b") | (frame.index % 40 == 0)]
    lopsided = tmp / "lopsided.csv"
    frame.to_csv(lopsided, index=False)

    out = tmp / "out_lopsided"
    subprocess.run(
        [sys.executable, str(_SCRIPT), "--summary-csv", str(lopsided),
         "--hierarchy", str(hierarchy), "--out-dir", str(out)],
        check=True, capture_output=True, text=True)
    by_node = pd.read_csv(out / "resolution_by_node.csv").set_index("node")
    astro = by_node.loc["Glia/Astro"]
    assert astro["winner_share_pct"] > 90, astro["winner_share_pct"]
    assert not bool(astro["at_chance"]), "dominance at a high posterior must NOT be flagged"
    print("  dominance alone is not flagged; only dominance at chance is")


def build_flat_run(tmp: Path) -> tuple[Path, Path]:
    """A flat run where the Excitatory profiles won almost nothing and a de novo letter won.

    This is the retina situation in miniature: `exc_0` took 3 cells in open competition and
    the tree then handed it 400, while `astro_a` is a kept de novo cluster that won openly
    under its letter and must NOT be flagged.
    """
    labels = (["exc_0"] * 3 + ["exc_1"] * 2            # the named profiles, barely winning
              + ["m"] * 600                             # the de novo letter that took them
              + ["astro_b"] * 200)
    typing = tmp / "flat.h5"
    with h5py.File(typing, "w") as f:
        dt = h5py.string_dtype(encoding="utf-8")
        f.create_dataset("cell_type", data=np.array(labels, dtype=object), dtype=dt)

    annotations = tmp / "annotations.csv"
    pd.DataFrame([
        {"denovo_label": "m", "keep": "1", "new_name": "astro_a"},
        {"denovo_label": "o", "keep": "0", "new_name": ""},
    ]).to_csv(annotations, index=False)
    return typing, annotations


def test_flat_vs_tree_flags_the_captive_contest(tmp: Path) -> None:
    summary, hierarchy = build(tmp)
    typing, annotations = build_flat_run(tmp)
    out = tmp / "out_flat"
    result = subprocess.run(
        [sys.executable, str(_SCRIPT), "--summary-csv", str(summary),
         "--hierarchy", str(hierarchy), "--flat-typing-h5", str(typing),
         "--denovo-annotations", str(annotations), "--out-dir", str(out)],
        check=True, capture_output=True, text=True)

    frame = pd.read_csv(out / "flat_vs_tree.csv").set_index("leaf")
    # exc_0: 3 cells won openly, 400 handed to it by the tree.
    assert frame.loc["exc_0", "flat_cells"] == 3
    assert frame.loc["exc_0", "tree_cells"] == N_BIG
    assert bool(frame.loc["exc_0", "inherited"]), "a contest decided among losing profiles"

    # astro_a is a KEPT de novo leaf: the flat run called those cells 'm'. Without the
    # rename it would look like it won nothing and be flagged -- the opposite of the truth.
    assert frame.loc["astro_a", "flat_label"] == "m"
    assert frame.loc["astro_a", "flat_cells"] == 600
    assert not bool(frame.loc["astro_a", "inherited"]), "a de novo leaf fit to this data"

    assert "contests inherited most of their cells" in result.stdout
    print("  flags the contest whose profiles never won openly, not the de novo leaf")


def test_missing_columns_error(tmp: Path) -> None:
    summary, hierarchy = build(tmp)
    frame = pd.read_csv(summary).drop(columns=["probs_annotLevel_2"])
    broken = tmp / "broken.csv"
    frame.to_csv(broken, index=False)
    result = subprocess.run(
        [sys.executable, str(_SCRIPT), "--summary-csv", str(broken),
         "--hierarchy", str(hierarchy), "--out-dir", str(tmp / "out_broken")],
        capture_output=True, text=True)
    assert result.returncode != 0
    assert "paired probability column" in result.stdout + result.stderr
    print("  a summary without paired probabilities errors")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("test_diagnose_insitutree_resolution")
        test_separates_resolved_from_at_chance(tmp)
        test_size_alone_does_not_flag(tmp)
        test_flat_vs_tree_flags_the_captive_contest(tmp)
        test_missing_columns_error(tmp)
    print("All tests passed.")


if __name__ == "__main__":
    main()
