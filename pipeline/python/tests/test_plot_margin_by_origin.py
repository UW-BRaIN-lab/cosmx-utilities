#!/usr/bin/env python3
"""Tests for plot_margin_by_origin.py — the Mac-side renderer for 75i Part 1.

The figure is only readable if the two groups are told apart correctly and each destination's
files are gathered into the right run, so those are what is tested here. The statistics
themselves belong to the R script and are covered by pipeline/R/tests/.

Runnable either under pytest or directly:
    PYTHONPATH=pipeline/python python pipeline/python/tests/test_plot_margin_by_origin.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location(
    "plot_margin_by_origin", _PY_DIR / "plot_margin_by_origin.py")
pm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pm)

DENSITY = pd.DataFrame({"scale": ["margin"] * 4, "group": ["A: l->Endo_capilar"] * 2
                        + ["B: Endo_capilar native"] * 2, "n": [100, 100, 50, 50],
                        "x": [-1.0, 1.0, -1.0, 1.0], "density": [0.1, 0.9, 0.9, 0.1]})


def test_group_labels_are_told_apart():
    """75i writes `A: l-><D>` and `B: <D> native`; only the second is the native group."""
    assert not pm.is_native("A: l->Endo_capilar")
    assert pm.is_native("B: Endo_capilar native")
    # The prefix alone must be enough — a destination whose NAME contains 'native' would
    # otherwise flip the A group's colour.
    assert not pm.is_native("A: l->Some_native_type")


def test_discover_groups_prefixed_files_by_destination():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        for dest in ("Endo_capilar", "Endo_arterial"):
            DENSITY.to_csv(root / f"l_vs_{dest}_margin_density_by_group.csv", index=False)
            (root / f"l_vs_{dest}_margin_summary.csv").write_text("group,n\nA,1\n")
            (root / f"l_vs_{dest}_argmax_consistency.csv").write_text("group,n\nA,1\n")
        runs = pm.discover(root)
    assert set(runs) == {"l_vs_Endo_capilar", "l_vs_Endo_arterial"}, runs
    assert runs["l_vs_Endo_capilar"]["summary"].name == "l_vs_Endo_capilar_margin_summary.csv"
    assert "argmax" in runs["l_vs_Endo_capilar"], runs


def test_discover_handles_a_single_unprefixed_run():
    """A directory fetched for one destination has no `<letter>_vs_<dest>_` prefix."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d) / "l_vs_Endo_capilar"
        root.mkdir()
        DENSITY.to_csv(root / "margin_density_by_group.csv", index=False)
        (root / "margin_summary.csv").write_text("group,n\nA,1\n")
        runs = pm.discover(root)
    assert list(runs) == ["l_vs_Endo_capilar"], runs
    assert runs["l_vs_Endo_capilar"]["summary"].name == "margin_summary.csv", runs


def test_discover_says_what_to_do_when_the_per_group_file_is_absent():
    """Directories fetched before the per-group patch have no curves — say so, don't crash."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "l_vs_Endo_capilar_margin_summary.csv").write_text("group,n\nA,1\n")
        try:
            pm.discover(root)
        except SystemExit as e:
            assert "re-run 75i" in str(e), e
        else:
            raise AssertionError("a directory with no density file must exit with guidance")


def test_summary_joins_argmax_consistency():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        DENSITY.to_csv(root / "l_vs_Endo_capilar_margin_density_by_group.csv", index=False)
        pd.DataFrame({"group": ["A: l->Endo_capilar", "B: Endo_capilar native"],
                      "n": [61721, 3456], "median_margin": [35.0, -1.2],
                      "pct_above_zero": [80.7, 47.9]}).to_csv(
            root / "l_vs_Endo_capilar_margin_summary.csv", index=False)
        pd.DataFrame({"group": ["A: l->Endo_capilar", "B: Endo_capilar native"],
                      "n": [61721, 3456], "pct_clust_is_argmax": [100.0, 63.9]}).to_csv(
            root / "l_vs_Endo_capilar_argmax_consistency.csv", index=False)
        summary = pm.summary_table(pm.discover(root))
    assert list(summary["run"].unique()) == ["l_vs_Endo_capilar"], summary
    native = summary[summary["group"].map(pm.is_native)]
    assert float(native["pct_clust_is_argmax"].iloc[0]) == 63.9, summary


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
