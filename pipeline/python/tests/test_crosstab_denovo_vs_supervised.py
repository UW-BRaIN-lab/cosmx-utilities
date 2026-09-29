#!/usr/bin/env python3
"""Tests for the de-novo vs supervised-GBmap cross-tab (no S3 / network).

Runnable either under pytest or directly:
    uv run --with pandas --with numpy --with h5py \\
        python pipeline/python/tests/test_crosstab_denovo_vs_supervised.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))  # crosstab_* imports its sibling anchor_profiles


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _PY_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ap = _load("anchor_profiles")
xt = _load("crosstab_denovo_vs_supervised")


def test_is_denovo_matches_only_letter_labels():
    """InSituType's cluster_name_pool is 1-2 lowercase letters; K=27 overflows to 'aa'."""
    assert ap.is_denovo("b")
    assert ap.is_denovo("aa")
    assert not ap.is_denovo("AC-like")
    assert not ap.is_denovo("Low_signal")
    assert not ap.is_denovo("TAM-MG")
    assert not ap.is_denovo("abc")


def test_split_named_denovo_preserves_order():
    named, denovo = ap.split_named_denovo(["AC-like", "b", "TAM-MG", "aa"])
    assert named == ["AC-like", "TAM-MG"]
    assert denovo == ["b", "aa"]


def test_summarize_row_ranks_destinations_by_share():
    counts = pd.Series({"AC-like": 70, "MES-like": 20, "Oligodendrocyte": 10})
    got = xt.summarize_row(counts, source="o")
    assert got["n_cells"] == 100
    assert got["gbmap_1"] == "AC-like" and got["gbmap_1_pct"] == 70.0
    assert got["gbmap_2"] == "MES-like" and got["gbmap_2_pct"] == 20.0
    assert got["gbmap_3"] == "Oligodendrocyte"


def test_summarize_row_pads_when_fewer_destinations_than_ranks():
    """A letter that maps onto a single GBmap type must not IndexError."""
    got = xt.summarize_row(pd.Series({"AC-like": 5, "MES-like": 0}), source="t")
    assert got["gbmap_1"] == "AC-like"
    assert got["gbmap_2"] == "" and np.isnan(got["gbmap_2_pct"])
    assert got["n_dest_90pct"] == 1


def test_n_dest_90pct_measures_dispersion():
    """1 = the letter is a rename of one GBmap type; large = the forced call is arbitrary."""
    clean = xt.summarize_row(pd.Series({"AC-like": 95, "MES-like": 5}), source="t")
    assert clean["n_dest_90pct"] == 1

    sprayed = xt.summarize_row(pd.Series({f"T{i}": 10 for i in range(10)}), source="b")
    assert sprayed["n_dest_90pct"] == 9  # nine types to reach 90%


def test_self_pct_is_named_only():
    """Named sources have a 'correct' forced answer; de-novo letters do not."""
    counts = pd.Series({"AC-like": 80, "MES-like": 20})
    assert xt.summarize_row(counts, source="AC-like")["self_pct"] == 80.0
    assert np.isnan(xt.summarize_row(counts, source="b")["self_pct"])


def test_load_display_labels_reads_annotation_column():
    labels = xt.load_display_labels(
        _PY_DIR.parents[0] / "reference/denovo_annotations/fullcohort_pruned_k27.csv")
    assert labels["b"] == "b - Low_signal sink"
    assert labels["o"] == "o - Hypoxia"


def test_load_display_labels_without_annotations_is_empty():
    assert xt.load_display_labels(None) == {}


def _write_map(tmp, rows):
    path = Path(tmp) / "map.csv"
    path.write_text("gbmap_type,compartment\n" + "\n".join(f"{a},{b}" for a, b in rows) + "\n")
    return path


def test_collapse_columns_conserves_cells_and_keeps_map_order():
    """Column order follows the map, not alphabetical — the Sankey axis must be stable."""
    ct = pd.DataFrame({"AC-like": [10, 0], "Neuron": [5, 20], "MES-like": [5, 0]},
                      index=["t", "n"])
    cmap = {"AC-like": "Malignant", "Neuron": "Normal_CNS", "MES-like": "Malignant"}
    got = xt.collapse_columns(ct, cmap)
    assert list(got.columns) == ["Malignant", "Normal_CNS"]
    assert got.values.sum() == ct.values.sum()
    assert got.loc["t", "Malignant"] == 15
    assert got.loc["n", "Normal_CNS"] == 20


def test_collapse_columns_refuses_to_silently_drop_unmapped_types():
    """An unmapped column would vanish from the roll-up and quietly change every percentage."""
    ct = pd.DataFrame({"AC-like": [10], "Mystery": [90]}, index=["t"])
    try:
        xt.collapse_columns(ct, {"AC-like": "Malignant"})
    except SystemExit as exc:
        assert "Mystery" in str(exc)
    else:
        raise AssertionError("expected SystemExit on an unmapped GBmap type")


def test_load_collapse_map_rejects_duplicate_types(tmp_path=None):
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        path = _write_map(tmp, [("AC-like", "Malignant"), ("AC-like", "Normal_CNS")])
        try:
            xt.load_collapse_map(path)
        except SystemExit as exc:
            assert "AC-like" in str(exc)
        else:
            raise AssertionError("expected SystemExit on a duplicated gbmap_type")


def test_shipped_compartment_map_is_complete_and_unambiguous():
    """The committed map must cover every GBmap type exactly once."""
    cmap = xt.load_collapse_map(
        _PY_DIR.parents[0] / "reference/gbmap_compartments.csv")
    assert len(cmap) == 54
    assert set(cmap.values()) == {"Malignant", "Normal_CNS", "Myeloid", "Lymphoid",
                                  "Vascular_mural", "Stress_sig"}


def test_summarize_row_self_key_scores_compartment_agreement():
    """A named source's 'unchanged' answer is its COMPARTMENT once the table is collapsed."""
    counts = pd.Series({"Malignant": 90, "Normal_CNS": 10})
    got = xt.summarize_row(counts, source="AC-like", prefix="compartment",
                           self_key="Malignant", n_top=2)
    assert got["compartment_1"] == "Malignant"
    assert got["self_pct"] == 90.0
    assert "compartment_3" not in got


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
