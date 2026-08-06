#!/usr/bin/env python3
"""Tests for scripts/generate-slide-metadata.py core logic (no S3 / network).

Runnable either under pytest or directly:  uv run python scripts/tests/test_generate_slide_metadata.py
"""
import csv
import gzip
import importlib.util
import os
import tempfile
from pathlib import Path

# Import the hyphenated script by file path.
_SCRIPT = Path(__file__).resolve().parents[1] / "generate-slide-metadata.py"
_spec = importlib.util.spec_from_file_location("gen_slide_metadata", _SCRIPT)
gsm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gsm)


def _write_gz(path, fieldnames, rows):
    with gzip.open(path, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _read_out(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


# ---- Fixtures as plain dicts ------------------------------------------------

FIELDS = ["cell_id", "cellSegmentationSetId", "TypeA", "TypeB", "Region"]
ROWS = [
    {"cell_id": "c_1_1_1", "cellSegmentationSetId": "seg-A", "TypeA": "astrocyte",
     "TypeB": "astrocyte", "Region": "retina"},
    {"cell_id": "c_1_1_2", "cellSegmentationSetId": "seg-A", "TypeA": "microglia",
     "TypeB": "microglia", "Region": "retina"},
    {"cell_id": "c_1_1_3", "cellSegmentationSetId": "seg-B", "TypeA": "rod",
     "TypeB": "rod", "Region": "brain"},
]


def test_multi_column_output_shape_and_colors():
    """Multiple --column specs produce cell_ID + (value,color) per column, and a
    label shared across two columns gets the SAME deterministic color."""
    specs = [
        gsm.ColumnSpec("celltype_norefit", "TypeA", "celltype_norefit_color"),
        gsm.ColumnSpec("celltype_refit", "TypeB", "celltype_refit_color"),
        gsm.ColumnSpec("Region", "Region", "Region_color"),
    ]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        out = os.path.join(d, "_metadata.csv")
        _write_gz(gz, FIELDS, ROWS)
        rows, seg_seen = gsm.read_rows(gz, None, specs)
        assert seg_seen == {"seg-A", "seg-B"}
        assert len(rows) == 3
        stats = gsm.write_output(out, specs, rows)
        assert stats["total_written"] == 3

        table = _read_out(out)
        header = table[0]
        assert header == [
            "cell_ID",
            "celltype_norefit", "celltype_norefit_color",
            "celltype_refit", "celltype_refit_color",
            "Region", "Region_color",
        ], header
        # cross-column color consistency: "astrocyte" appears in both type columns.
        idx = {h: i for i, h in enumerate(header)}
        astro = table[1]
        assert astro[idx["celltype_norefit"]] == "astrocyte"
        assert astro[idx["celltype_refit"]] == "astrocyte"
        assert astro[idx["celltype_norefit_color"]] == astro[idx["celltype_refit_color"]], \
            "same label in two columns must map to the same color"
        # colors are valid hex
        assert astro[idx["Region_color"]].startswith("#") and len(astro[idx["Region_color"]]) == 7


def test_seg_id_filter():
    """--seg-id keeps only rows whose cellSegmentationSetId matches."""
    specs = [gsm.ColumnSpec("Region", "Region", "Region_color")]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        _write_gz(gz, FIELDS, ROWS)
        rows, _ = gsm.read_rows(gz, {"seg-B"}, specs)
        assert len(rows) == 1 and rows[0][0] == "c_1_1_3"


def test_missing_source_column_is_empty_not_crash():
    """A requested source header absent from the file yields empty values."""
    specs = [gsm.ColumnSpec("niche", "DoesNotExist", "niche_color")]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        out = os.path.join(d, "_metadata.csv")
        _write_gz(gz, FIELDS, ROWS)
        rows, _ = gsm.read_rows(gz, None, specs)
        stats = gsm.write_output(out, specs, rows)
        assert stats["value_counts"] == {"niche": 0}
        table = _read_out(out)
        assert table[0] == ["cell_ID", "niche", "niche_color"]
        assert all(r[1] == "" and r[2] == "" for r in table[1:])


def test_legacy_mode_emits_cell_type_and_hex_color():
    """No --column falls back to the legacy cell_type + hex_color output."""
    class Args:
        column = None
        cell_type_column = "TypeA"
    specs = gsm.build_specs(Args())
    assert len(specs) == 1
    assert specs[0].out == "cell_type" and specs[0].color_col == "hex_color"
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        out = os.path.join(d, "_metadata.csv")
        _write_gz(gz, FIELDS, ROWS)
        rows, _ = gsm.read_rows(gz, None, specs)
        gsm.write_output(out, specs, rows)
        assert _read_out(out)[0] == ["cell_ID", "cell_type", "hex_color"]


def test_legacy_autodetect_source():
    """Legacy mode with no source auto-detects the InSituType clusters column."""
    class Args:
        column = None
        cell_type_column = None
    specs = gsm.build_specs(Args())
    assert specs[0].src is None
    fields = ["cell_id", "cellSegmentationSetId",
              "RNA_RNA_Cell.Typing.InSituType.No.Refit_1_clusters"]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        _write_gz(gz, fields, [
            {"cell_id": "c_1_1_1", "cellSegmentationSetId": "seg-A",
             "RNA_RNA_Cell.Typing.InSituType.No.Refit_1_clusters": "astrocyte"}])
        resolved = gsm.resolve_legacy_src(specs, gz)
        assert resolved[0].src == "RNA_RNA_Cell.Typing.InSituType.No.Refit_1_clusters"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
