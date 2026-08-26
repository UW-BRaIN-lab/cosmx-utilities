#!/usr/bin/env python3
"""Tests for scripts/atomx-leiden-typing-crosstab.py core logic (no S3 / network).

Runnable either under pytest or directly:
    uv run python scripts/tests/test_atomx_leiden_typing_crosstab.py
"""
import datetime as dt
import importlib.util
from pathlib import Path

import pandas as pd

_SCRIPT = Path(__file__).resolve().parents[1] / "atomx-leiden-typing-crosstab.py"
_spec = importlib.util.spec_from_file_location("atomx_leiden_typing_crosstab", _SCRIPT)
alt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(alt)

LEIDEN_HDR = "RNA_Neighbor.network.expression.space.1_1_cluster_RNA_Leiden.Clustering.1_1"


def _typing(version):
    return f"RNA_RNA_Cell.Typing.InSituType.{version}_clusters"


class _FakeS3:
    """Stands in for boto3's S3 client for the listing paginator only."""

    def __init__(self, objects):
        self.objects = objects

    def get_paginator(self, _name):
        outer = self

        class _P:
            def paginate(self, Bucket, Prefix):
                yield {"Contents": [
                    o for o in outer.objects if o["Key"].startswith(Prefix)
                ]}
        return _P()


def _obj(key, day):
    return {"Key": key, "LastModified": dt.datetime(2026, 8, day)}


# ---- authoritative-file selection -----------------------------------------

STUDY = "CosMx-Maddie/3D_study"
SLIDE = "20260708_UWA_599"
OUTER = f"{STUDY}/flatFiles/{SLIDE}/{SLIDE}_metadata_file.csv.gz"
NESTED = f"{STUDY}/flatFiles/rerun_20_08/flatFiles/{SLIDE}/{SLIDE}_metadata_file.csv.gz"


def test_nested_reexport_supersedes_the_outer_flat_file():
    """The real 3D study holds a nested re-export; the stitched metadata came
    from it, so the deeper file must win even though both exist."""
    s3 = _FakeS3([_obj(OUTER, 18), _obj(NESTED, 20)])
    keys = alt.newest_metadata_per_slide(s3, "b", STUDY)
    assert keys == {SLIDE: NESTED}, keys


def test_nested_wins_even_when_the_outer_file_is_newer():
    """Depth beats timestamp: a re-uploaded outer file is still superseded."""
    s3 = _FakeS3([_obj(OUTER, 25), _obj(NESTED, 20)])
    assert alt.newest_metadata_per_slide(s3, "b", STUDY) == {SLIDE: NESTED}


def test_newest_wins_among_equal_depth():
    a = f"{STUDY}/flatFiles/{SLIDE}/{SLIDE}_metadata_file.csv.gz"
    s3 = _FakeS3([_obj(a, 12)])
    assert alt.newest_metadata_per_slide(s3, "b", STUDY) == {SLIDE: a}


def test_each_slide_gets_its_own_file():
    other = "20260708_UWA_787"
    k1 = f"{STUDY}/flatFiles/{SLIDE}/{SLIDE}_metadata_file.csv.gz"
    k2 = f"{STUDY}/flatFiles/{other}/{other}_metadata_file.csv.gz"
    s3 = _FakeS3([_obj(k1, 12), _obj(k2, 12)])
    assert alt.newest_metadata_per_slide(s3, "b", STUDY) == {SLIDE: k1, other: k2}


def test_non_metadata_objects_are_ignored():
    s3 = _FakeS3([_obj(f"{STUDY}/flatFiles/{SLIDE}/{SLIDE}_exprMat_file.csv.gz", 12)])
    assert alt.newest_metadata_per_slide(s3, "b", STUDY) == {}


# ---- typing-run selection --------------------------------------------------

def test_highest_run_index_wins_not_header_order():
    """The real nested export lists 2_1 before 1_1, but position must not decide:
    reversed here so a position-based rule would pick the stale run."""
    headers = ["cell_id", _typing("1_1"), _typing("2_1"), LEIDEN_HDR]
    leiden, typing = alt.detect_columns(headers)
    assert typing == _typing("2_1"), typing
    assert leiden == LEIDEN_HDR


def test_minor_version_compared_numerically():
    """The 2D study's only run is 1_5; a string sort would rank 1_5 below 1_10."""
    headers = ["cell_id", _typing("1_5"), _typing("1_10"), LEIDEN_HDR]
    _, typing = alt.detect_columns(headers)
    assert typing == _typing("1_10"), typing


def test_single_run_is_selected():
    headers = ["cell_id", _typing("1_5"), LEIDEN_HDR]
    _, typing = alt.detect_columns(headers)
    assert typing == _typing("1_5")


def test_pinned_version_overrides_the_newest_rule():
    headers = ["cell_id", _typing("1_1"), _typing("2_1"), LEIDEN_HDR]
    _, typing = alt.detect_columns(headers, typing_version="1_1")
    assert typing == _typing("1_1"), typing


def test_pinned_version_absent_raises_and_lists_what_is_present():
    headers = ["cell_id", _typing("1_1"), _typing("2_1"), LEIDEN_HDR]
    try:
        alt.detect_columns(headers, typing_version="3_1")
    except alt.ColumnDetectionError as e:
        assert "3_1" in str(e)
        assert "1_1" in str(e) and "2_1" in str(e), str(e)
    else:
        raise AssertionError("expected ColumnDetectionError for an absent version")


def test_missing_leiden_column_raises():
    try:
        alt.detect_columns(["cell_id", _typing("1_1")])
    except alt.ColumnDetectionError as e:
        assert "Leiden" in str(e)
    else:
        raise AssertionError("expected ColumnDetectionError with no Leiden column")


def test_missing_typing_column_raises():
    try:
        alt.detect_columns(["cell_id", LEIDEN_HDR])
    except alt.ColumnDetectionError as e:
        assert "InSituType" in str(e)
    else:
        raise AssertionError("expected ColumnDetectionError with no typing column")


def test_explicit_leiden_column_must_exist():
    headers = ["cell_id", _typing("1_1"), LEIDEN_HDR]
    try:
        alt.detect_columns(headers, leiden_column="not_a_column")
    except alt.ColumnDetectionError as e:
        assert "not_a_column" in str(e)
    else:
        raise AssertionError("expected ColumnDetectionError for a bogus Leiden column")


# ---- purity ----------------------------------------------------------------

def test_leiden_purity_reports_dominant_type_and_share():
    counts = pd.DataFrame(
        {"g": [90, 10], "Astrocyte": [10, 40], "Unassigned": [0, 50]},
        index=["1", "2"],
    )
    purity = alt.leiden_purity(counts).set_index("leiden")
    assert purity.loc["1", "dominant_type"] == "g"
    assert purity.loc["1", "dominant_frac"] == 0.9
    assert purity.loc["1", "n_cells"] == 100
    assert purity.loc["1", "n_types_present"] == 2      # Unassigned is zero here
    assert purity.loc["2", "dominant_type"] == "Unassigned"
    assert purity.loc["2", "n_types_present"] == 3


def test_leiden_purity_is_ordered_by_cluster_size():
    counts = pd.DataFrame({"g": [5, 500, 50]}, index=["1", "2", "3"])
    assert alt.leiden_purity(counts)["leiden"].tolist() == ["2", "3", "1"]


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
