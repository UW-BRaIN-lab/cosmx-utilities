#!/usr/bin/env python3
"""Tests for the anchor-cohort UMAP: alignment, colouring and an end-to-end render (no S3).

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_plot_anchor_cohort_umap.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location("plot_anchor_cohort_umap",
                                               _PY_DIR / "plot_anchor_cohort_umap.py")
pacu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pacu)

CELL_IDS = [f"slide_F1_C{i}" for i in range(1, 11)]
COHORT_LABELS = {"slide_F1_C2": "b", "slide_F1_C4": "z", "slide_F1_C5": "aa", "slide_F1_C9": "OPC-like"}


def test_is_denovo_letter():
    assert pacu.is_denovo_letter("b") and pacu.is_denovo_letter("aa")
    assert not pacu.is_denovo_letter("OPC-like") and not pacu.is_denovo_letter("TAM-MG_aging_sig")


def test_alignment_masks_cohort_cells_in_embedding_order():
    labels = pd.Series(COHORT_LABELS)
    in_cohort, aligned = pacu.align_to_embedding(pd.Index(CELL_IDS), labels)
    assert in_cohort.sum() == 4
    assert list(aligned.index) == [c for c in CELL_IDS if c in COHORT_LABELS]
    assert list(aligned) == ["b", "z", "aa", "OPC-like"]


def test_missing_cohort_cell_raises():
    labels = pd.Series({**COHORT_LABELS, "slide_F9_C9": "b"})
    try:
        pacu.align_to_embedding(pd.Index(CELL_IDS), labels)
    except KeyError as error:
        assert "absent" in str(error)
    else:
        raise AssertionError("a cohort cell without a UMAP position must raise")


def test_denovo_only_colours_named_types_neutral():
    colors = pacu.label_colors(["b", "z", "OPC-like"], denovo_only=True)
    assert colors["OPC-like"] == pacu.NAMED_TYPES_COLOR
    assert colors["b"] != pacu.NAMED_TYPES_COLOR != colors["z"] and colors["b"] != colors["z"]
    all_colors = pacu.label_colors(["b", "z", "OPC-like"], denovo_only=False)
    assert all_colors["OPC-like"] != pacu.NAMED_TYPES_COLOR


def test_fixed_colors_used_and_missing_label_raises():
    fixed = {"b": "#111111", "z": "#222222", "OPC-like": "#333333"}
    colors = pacu.label_colors(["b", "OPC-like"], denovo_only=False, fixed_colors=fixed)
    assert colors == {"b": "#111111", "OPC-like": "#333333"}
    try:
        pacu.label_colors(["b", "q"], denovo_only=False, fixed_colors=fixed)
    except KeyError as error:
        assert "q" in str(error)
    else:
        raise AssertionError("a label with no Napari colour must raise")


def test_shipped_colour_table_covers_every_anchor_label():
    table = pacu.read_label_colors(_PY_DIR.parent / "reference" / "anchor_label_colors.csv")
    assert len(table) == 81 and len(set(table.values())) == 81


def test_select_panels_keeps_big_letters_and_reports_skipped():
    labels = pd.Series(["b"] * 5 + ["z"] * 3 + ["h"] * 1 + ["OPC-like"] * 9)
    panels, skipped = pacu.select_panels(labels, min_cells=3)
    assert panels == ["b", "z"]
    assert skipped == ["h"]
    named, named_skipped = pacu.select_panels(labels, min_cells=10, denovo=False)
    assert named == [] and named_skipped == ["OPC-like"]
    named, _ = pacu.select_panels(labels, min_cells=3, denovo=False)
    assert named == ["OPC-like"]


def test_end_to_end_render():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        rng = np.random.RandomState(0)
        adata = ad.AnnData(X=np.zeros((len(CELL_IDS), 2), dtype=np.float32),
                           obs=pd.DataFrame(index=CELL_IDS))
        adata.obsm["X_umap"] = rng.normal(size=(len(CELL_IDS), 2))
        adata.write_h5ad(tmp / "typed.h5ad")
        with h5py.File(tmp / "anchor_typing.h5", "w") as handle:
            handle["cell_id"] = np.array(list(COHORT_LABELS), dtype=h5py.string_dtype())
            handle["cell_type"] = np.array(list(COHORT_LABELS.values()), dtype=h5py.string_dtype())
        sys.argv = ["x", "--typed-h5ad", str(tmp / "typed.h5ad"),
                    "--anchor-typing-h5", str(tmp / "anchor_typing.h5"),
                    "--output-stem", str(tmp / "out"), "--min-panel-cells", "1"]
        pacu.main()
        for suffix in ("all_labels", "denovo_letters", "small_multiples", "small_multiples_named"):
            assert (tmp / f"out_{suffix}.png").stat().st_size > 0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            print(f"{name} ...", end=" ")
            fn()
            print("ok")
