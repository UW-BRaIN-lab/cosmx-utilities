#!/usr/bin/env python3
"""Tests for plot_crosstab_sankey.py node ordering (no S3 / network).

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_plot_crosstab_sankey.py
"""
import importlib.util
from pathlib import Path

import pandas as pd

_SCRIPT = Path(__file__).resolve().parents[1] / "plot_crosstab_sankey.py"
_spec = importlib.util.spec_from_file_location("plot_crosstab_sankey", _SCRIPT)
pcs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pcs)


def test_natural_key_orders_numeric_labels_numerically():
    """Leiden clusters are numeric strings: plain string sort would give 1,10,2."""
    labels = ["10", "2", "1", "22", "3"]
    assert sorted(labels, key=pcs.natural_key) == ["1", "2", "3", "10", "22"]


def test_natural_key_handles_non_numeric_labels():
    labels = ["Oligodendrocyte", "astrocyte", "Low_signal"]
    assert sorted(labels, key=pcs.natural_key) == [
        "astrocyte", "Low_signal", "Oligodendrocyte",
    ]


def test_order_axis_size_sorts_descending_by_total():
    totals = pd.Series({"a": 5, "b": 100, "c": 50})
    assert pcs.order_axis(totals, "size") == ["b", "c", "a"]


def test_order_axis_natural_ignores_size():
    totals = pd.Series({"10": 900, "2": 1, "1": 500})
    assert pcs.order_axis(totals, "natural") == ["1", "2", "10"]


def test_order_axis_returns_every_node_exactly_once():
    """A dropped or duplicated node would silently corrupt the ribbon layout."""
    totals = pd.Series({str(i): i * 10 for i in range(15)})
    for how in ("size", "natural"):
        order = pcs.order_axis(totals, how)
        assert sorted(order) == sorted(totals.index), how
        assert len(order) == len(set(order)) == len(totals), how


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
