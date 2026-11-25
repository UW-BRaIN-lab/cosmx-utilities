#!/usr/bin/env python3
"""Tests for pca_batch.py: the surgical batch split and Case/Block-as-text coercion used by stage 3a.

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_pca_batch.py
"""
import contextlib
import io
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
import pca_batch as pb  # noqa: E402


def _slide(obs: dict, prefix: str) -> ad.AnnData:
    n = len(next(iter(obs.values())))
    index = [f"{prefix}_F{i + 1}_C1" for i in range(n)]
    return ad.AnnData(np.ones((n, 1), dtype=np.float32), obs=pd.DataFrame(obs, index=index))


def test_mixed_case_types_crash_the_write_without_coercion():
    """The failure this module exists for: a numeric-Case slide concatenated with a Nerve slide."""
    numeric = _slide({"Case": [7488, 7488], "Region": ["Surgical"] * 2}, "S1")
    nerve = _slide({"Case": ["7706", "Nerve"], "Region": ["Contralateral uninvolved", "Nerve"]}, "S2")
    combined = ad.concat([numeric, nerve], join="outer", merge="same")
    print("uncoerced Case types:", combined.obs["Case"].map(type).tolist())
    with tempfile.TemporaryDirectory() as tmp:
        try:
            combined.write_h5ad(Path(tmp) / "combined.h5ad")
        except TypeError as err:
            print(f"write failed as expected: {err}")
        else:
            raise AssertionError("mixed int/str Case wrote cleanly; the coercion may be moot")


def test_coerced_slides_concatenate_and_write():
    numeric = _slide({"Case": [7488, 7488], "Block": ["A3", "A3"]}, "S1")
    float_case = _slide({"Case": [7706.0, np.nan], "Block": ["A7", np.nan]}, "S2")
    nerve = _slide({"Case": ["7706", "Nerve"], "Block": ["A7", "Nerve"]}, "S3")
    slides = []
    for slide in (numeric, float_case, nerve):
        slide.obs = pb.identifiers_as_text(slide.obs)
        slides.append(slide)
    combined = ad.concat(slides, join="outer", merge="same")
    print("coerced Case:", combined.obs["Case"].tolist())
    case = combined.obs["Case"].astype(object)
    assert case.iloc[:3].tolist() == ["7488", "7488", "7706"], "7706.0 must become '7706'"
    assert pd.isna(case.iloc[3]), "a missing Case stays missing, not the text 'nan'"
    assert case.iloc[4:].tolist() == ["7706", "Nerve"]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "combined.h5ad"
        combined.write_h5ad(path)
        back = ad.read_h5ad(path)
    print("round-tripped Case:", back.obs["Case"].astype(object).tolist())
    assert back.obs["Case"].astype(object).iloc[-1] == "Nerve"


def test_surgical_cells_get_their_own_batch():
    obs = pd.DataFrame({
        "Case": ["7488", "7488", "7302", "7302", "7302"],
        "Region": ["Infiltrating edge", "Surgical", "Tumor bulk", "Surgical", "Surgical"],
    })
    batch = pb.pca_batch_labels(obs, "Case", ["Surgical"])
    print(batch.tolist())
    assert batch.tolist() == ["7488", "7488_Surgical", "7302", "7302_Surgical", "7302_Surgical"]
    assert batch[obs.Region != "Surgical"].isin(batch[obs.Region == "Surgical"]).sum() == 0


def test_spaces_in_region_become_underscores():
    obs = pd.DataFrame({"Case": ["1"], "Region": ["Tumor bulk"]})
    assert pb.pca_batch_labels(obs, "Case", ["Tumor bulk"]).tolist() == ["1_Tumor_bulk"]


def test_absent_split_region_leaves_batches_alone():
    obs = pd.DataFrame({"Case": ["1", "2"], "Region": ["retina", "brain"]})
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        batch = pb.pca_batch_labels(obs, "Case", ["Surgical"])
    out = captured.getvalue()
    print(out)
    assert batch.tolist() == ["1", "2"]
    assert "not present" in out


def test_no_region_column_or_no_split_is_plain_batch():
    obs = pd.DataFrame({"slide_id": ["a", "b"]})
    assert pb.pca_batch_labels(obs, "slide_id", ["Surgical"]).tolist() == ["a", "b"]
    obs = pd.DataFrame({"Case": ["1"], "Region": ["Surgical"]})
    assert pb.pca_batch_labels(obs, "Case", []).tolist() == ["1"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
