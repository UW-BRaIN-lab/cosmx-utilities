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


def test_missing_source_column_becomes_unassigned_not_crash():
    """A requested source header absent from the file yields Unassigned, not blanks.

    An all-blank column is invisible in the viewer and, because pandas reads
    blanks as NaN, used to crash the color-by widget outright. Naming the
    category keeps those cells visible and countable.
    """
    specs = [gsm.ColumnSpec("niche", "DoesNotExist", "niche_color")]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        out = os.path.join(d, "_metadata.csv")
        _write_gz(gz, FIELDS, ROWS)
        rows, _ = gsm.read_rows(gz, None, specs)
        stats = gsm.write_output(out, specs, rows)
        assert stats["value_counts"] == {"niche": 1}, stats
        table = _read_out(out)
        assert table[0] == ["cell_ID", "niche", "niche_color"]
        assert all(r[1] == gsm.UNASSIGNED_LABEL for r in table[1:])
        # Every value carries a color, so nothing renders as an unlabelled gap.
        assert all(r[2].startswith("#") for r in table[1:])


def test_blank_typing_values_become_unassigned():
    """Cells that failed QC are left blank by cell typing; they must come out as
    an explicit category alongside the real types."""
    specs = [gsm.ColumnSpec("Cell Type", "TypeA", "Cell Type_color")]
    fields = ["cell_id", "cellSegmentationSetId", "TypeA"]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "in.csv.gz")
        out = os.path.join(d, "_metadata.csv")
        _write_gz(gz, fields, [
            {"cell_id": "c_1_1_1", "cellSegmentationSetId": "s", "TypeA": "Microglia.A"},
            {"cell_id": "c_1_1_2", "cellSegmentationSetId": "s", "TypeA": ""},
            {"cell_id": "c_1_1_3", "cellSegmentationSetId": "s", "TypeA": "Astrocyte.B"},
        ])
        rows, _ = gsm.read_rows(gz, None, specs)
        gsm.write_output(out, specs, rows)
        table = _read_out(out)
        values = [r[1] for r in table[1:]]
        assert values == ["Microglia.A", gsm.UNASSIGNED_LABEL, "Astrocyte.B"], values
        colors = {r[1]: r[2] for r in table[1:]}
        assert colors[gsm.UNASSIGNED_LABEL].startswith("#")
        # Unassigned must not collide with a real type's color.
        assert len(set(colors.values())) == 3, colors


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


# ---- Nested exports, source ranking, and cross-study annotation fill --------

class FakeS3:
    """Minimal S3 stand-in backed by {key: bytes}, for discovery/ranking tests."""

    def __init__(self, objects):
        self.objects = objects

    def list_objects_v2(self, Bucket, Prefix, Delimiter=None):
        children = set()
        for key in self.objects:
            if not key.startswith(Prefix):
                continue
            rest = key[len(Prefix):]
            if "/" in rest:
                children.add(Prefix + rest.split("/", 1)[0] + "/")
        return {"CommonPrefixes": [{"Prefix": p} for p in sorted(children)]}

    def head_object(self, Bucket, Key):
        if Key not in self.objects:
            from botocore.exceptions import ClientError
            raise ClientError({"Error": {"Code": "404"}}, "HeadObject")
        return {"ContentLength": len(self.objects[Key])}

    def get_object(self, Bucket, Key, Range=None):
        import io
        data = self.objects[Key]
        if Range:
            end = int(Range.split("-")[1])
            data = data[:end + 1]
        return {"Body": io.BytesIO(data)}


def _gz_bytes(fieldnames, rows):
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        text = io.StringIO()
        w = csv.DictWriter(text, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
        gz.write(text.getvalue().encode())
    return buf.getvalue()


BUCKET = "test-bucket"
STUDY = "CosMx-Maddie/RUN_3D"
SLIDE = "UWA_599_657"
NESTED_RUN = "RERUN_20260820"
OUTER_KEY = f"{STUDY}/flatFiles/{SLIDE}/{SLIDE}_metadata_file.csv.gz"
NESTED_KEY = (f"{STUDY}/flatFiles/{NESTED_RUN}/flatFiles/{SLIDE}/"
              f"{SLIDE}_metadata_file.csv.gz")


def test_find_metadata_file_discovers_nested_export():
    """A re-run AtoMx export nests under the original run's flatFiles and must
    still be found — and listed ahead of its parent."""
    s3 = FakeS3({
        OUTER_KEY: _gz_bytes(["cell_id"], [{"cell_id": "c_1_1_1"}]),
        NESTED_KEY: _gz_bytes(["cell_id"], [{"cell_id": "c_1_1_1"}]),
    })
    found = gsm.find_metadata_file(s3, BUCKET, "CosMx-Maddie", SLIDE)
    labels = [label for label, _key in found]
    assert labels == [f"RUN_3D/{NESTED_RUN}", "RUN_3D"], labels


def test_read_gz_header_from_ranged_get():
    """Headers come back without pulling the whole object down."""
    s3 = FakeS3({OUTER_KEY: _gz_bytes(
        ["cell_id", "fov", "TypeA"],
        [{"cell_id": "c_1_1_1", "fov": "1", "TypeA": "astrocyte"}])})
    assert gsm.read_gz_header(s3, BUCKET, OUTER_KEY) == ["cell_id", "fov", "TypeA"]


def test_rank_sources_prefers_file_with_requested_columns():
    """Both exports share a segmentation ID, so the one holding the requested
    typing column must win regardless of discovery order."""
    wanted = "RNA_RNA_Cell.Typing.InSituType.2_1_clusters"
    s3 = FakeS3({
        OUTER_KEY: _gz_bytes(["cell_id", "fov"], [{"cell_id": "c", "fov": "1"}]),
        NESTED_KEY: _gz_bytes(["cell_id", "fov", wanted],
                              [{"cell_id": "c", "fov": "1", wanted: "Microglia.A"}]),
    })
    specs = [gsm.ColumnSpec("Cell Type", wanted, "Cell Type_color")]
    # Deliberately hand them over in the "wrong" order.
    ranked = gsm.rank_sources(s3, BUCKET, [("outer", OUTER_KEY), ("nested", NESTED_KEY)], specs)
    assert [label for label, _ in ranked] == ["nested", "outer"]


def test_column_source_fallback_picks_first_present_header():
    """`A|B` takes whichever header the file actually has, so one command can
    serve studies that renamed a column between runs."""
    spec = gsm.ColumnSpec("Case Specific", "Case_specific_SORL1|Case_specific",
                          "Case Specific_color")
    assert spec.resolve(["cell_id", "Case_specific"]) == "Case_specific"
    assert spec.resolve(["cell_id", "Case_specific_SORL1"]) == "Case_specific_SORL1"
    # First listed wins when both are present.
    assert spec.resolve(["Case_specific", "Case_specific_SORL1"]) == "Case_specific_SORL1"
    assert spec.resolve(["cell_id"]) == ""


FILL_FIELDS = ["cell_id", "fov", "cellSegmentationSetId", "UWA"]


def test_fill_missing_annotations_joins_on_fov():
    """Blank annotations are filled per-FOV from a donor study of the same slide;
    values already present are never overwritten."""
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color")]
    with tempfile.TemporaryDirectory() as d:
        primary_gz = os.path.join(d, "primary.csv.gz")
        donor_gz = os.path.join(d, "donor.csv.gz")
        # Primary: FOV 1 blank, FOV 2 blank, FOV 3 already annotated.
        _write_gz(primary_gz, FILL_FIELDS, [
            {"cell_id": "c_1_1_1", "fov": "1", "cellSegmentationSetId": "s", "UWA": ""},
            {"cell_id": "c_1_1_2", "fov": "1", "cellSegmentationSetId": "s", "UWA": ""},
            {"cell_id": "c_1_2_1", "fov": "2", "cellSegmentationSetId": "s", "UWA": ""},
            {"cell_id": "c_1_3_1", "fov": "3", "cellSegmentationSetId": "s", "UWA": "keep"},
        ])
        # Donor has different cell IDs (different segmentation) but the same FOVs.
        _write_gz(donor_gz, FILL_FIELDS, [
            {"cell_id": "d_9_1_7", "fov": "1", "cellSegmentationSetId": "t", "UWA": "7796"},
            {"cell_id": "d_9_2_4", "fov": "2", "cellSegmentationSetId": "t", "UWA": "7665"},
            {"cell_id": "d_9_3_2", "fov": "3", "cellSegmentationSetId": "t", "UWA": "other"},
        ])
        rows, _ = gsm.read_rows(primary_gz, None, specs)
        donor_rows, _ = gsm.read_rows(donor_gz, None, specs)
        filled, skipped = gsm.fill_missing_annotations(rows, donor_rows, specs)

        assert filled == {"UWA": 3}, filled
        assert skipped == [], skipped
        by_cell = {cell_id: values[0] for cell_id, _fov, values in rows}
        assert by_cell["c_1_1_1"] == "7796"
        assert by_cell["c_1_1_2"] == "7796"
        assert by_cell["c_1_2_1"] == "7665"
        assert by_cell["c_1_3_1"] == "keep", "existing values must not be overwritten"


def test_fill_leaves_fovs_the_donor_lacks_blank():
    """An FOV absent from the donor stays blank rather than borrowing a neighbour."""
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color")]
    with tempfile.TemporaryDirectory() as d:
        primary_gz = os.path.join(d, "primary.csv.gz")
        donor_gz = os.path.join(d, "donor.csv.gz")
        _write_gz(primary_gz, FILL_FIELDS, [
            {"cell_id": "c_1_9_1", "fov": "9", "cellSegmentationSetId": "s", "UWA": ""},
        ])
        _write_gz(donor_gz, FILL_FIELDS, [
            {"cell_id": "d_1_1_1", "fov": "1", "cellSegmentationSetId": "t", "UWA": "7796"},
        ])
        rows, _ = gsm.read_rows(primary_gz, None, specs)
        donor_rows, _ = gsm.read_rows(donor_gz, None, specs)
        filled, _skipped = gsm.fill_missing_annotations(rows, donor_rows, specs)
        assert filled == {}
        assert rows[0][2][0] == ""


def test_fill_refuses_per_cell_columns():
    """A column that varies within an FOV is per-cell (e.g. cell typing, tied to
    its own segmentation) and must never be transferred by FOV join."""
    fields = ["cell_id", "fov", "cellSegmentationSetId", "UWA", "CellType"]
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color"),
             gsm.ColumnSpec("Cell Type", "CellType", "Cell Type_color")]
    with tempfile.TemporaryDirectory() as d:
        primary_gz = os.path.join(d, "primary.csv.gz")
        donor_gz = os.path.join(d, "donor.csv.gz")
        _write_gz(primary_gz, fields, [
            {"cell_id": "c_1_1_1", "fov": "1", "cellSegmentationSetId": "s",
             "UWA": "", "CellType": ""},
        ])
        # Donor FOV 1 holds one case but several cell types.
        _write_gz(donor_gz, fields, [
            {"cell_id": "d_1_1_1", "fov": "1", "cellSegmentationSetId": "t",
             "UWA": "7796", "CellType": "Microglia.A"},
            {"cell_id": "d_1_1_2", "fov": "1", "cellSegmentationSetId": "t",
             "UWA": "7796", "CellType": "Astrocyte.B"},
        ])
        rows, _ = gsm.read_rows(primary_gz, None, specs)
        donor_rows, _ = gsm.read_rows(donor_gz, None, specs)
        filled, skipped = gsm.fill_missing_annotations(rows, donor_rows, specs)

        assert filled == {"UWA": 1}, filled
        assert skipped == ["Cell Type"], skipped
        assert rows[0][2][0] == "7796", "FOV-constant case annotation fills"
        assert rows[0][2][1] == "", "per-cell typing must stay blank"


# ---- Per-FOV annotation sheets ---------------------------------------------

SHEET_FIELDS = ["Flow Cells", "FOVs", "Case_broad", "Case_specific_SORL1",
                "SORL1_mutation", "UWA"]
SHEET_SLIDE = "20260708_UWA_599_657_710_741"


def _write_sheet(path, rows, fieldnames=SHEET_FIELDS):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _sheet_rows(n=3, slide=SHEET_SLIDE):
    return [{"Flow Cells": slide, "FOVs": str(i), "Case_broad": "AD+LATE",
             "Case_specific_SORL1": "AD+LATE SORL1",
             "SORL1_mutation": "AD+LATE SORL1 R953C", "UWA": "741"}
            for i in range(1, n + 1)]


ANNOTATION_SPECS = [
    ("Case Broad", "Case_broad"),
    ("Case Specific", "Case_specific_SORL1|Case_specific"),
    ("SORL1 Mutation", "SORL1_mutation"),
    ("UWA", "UWA"),
]


def _annotation_specs():
    return [gsm.ColumnSpec(out, src, f"{out}_color") for out, src in ANNOTATION_SPECS]


def test_annotation_sheet_fills_blank_metadata():
    """A per-FOV sheet fills the annotations AtoMx never captured."""
    specs = _annotation_specs()
    fields = ["cell_id", "fov", "cellSegmentationSetId",
              "Case_broad", "Case_specific_SORL1", "SORL1_mutation", "UWA"]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "meta.csv.gz")
        sheet = os.path.join(d, "sheet.csv")
        _write_gz(gz, fields, [
            {"cell_id": "c_1_1_1", "fov": "1", "cellSegmentationSetId": "s",
             "Case_broad": "", "Case_specific_SORL1": "", "SORL1_mutation": "", "UWA": ""},
            {"cell_id": "c_1_2_1", "fov": "2", "cellSegmentationSetId": "s",
             "Case_broad": "", "Case_specific_SORL1": "", "SORL1_mutation": "", "UWA": ""},
        ])
        _write_sheet(sheet, _sheet_rows())

        rows, _ = gsm.read_rows(gz, None, specs)
        sheet_rows = gsm.read_annotation_csv(sheet, SHEET_SLIDE, specs)
        filled, skipped = gsm.fill_missing_annotations(rows, sheet_rows, specs)

        assert skipped == [], skipped
        assert filled == {"Case Broad": 2, "Case Specific": 2,
                          "SORL1 Mutation": 2, "UWA": 2}, filled
        assert rows[0][2] == ["AD+LATE", "AD+LATE SORL1", "AD+LATE SORL1 R953C", "741"]


def test_annotation_sheet_rejects_wrong_slide():
    """Pointing at another slide's sheet is the likely wiring mistake, so it
    must fail loudly rather than annotate cells with the wrong case."""
    specs = _annotation_specs()
    with tempfile.TemporaryDirectory() as d:
        sheet = os.path.join(d, "sheet.csv")
        _write_sheet(sheet, _sheet_rows(slide="20260708_UWA_787_795_6589_6745"))
        try:
            gsm.read_annotation_csv(sheet, SHEET_SLIDE, specs)
        except ValueError as e:
            assert "787_795_6589_6745" in str(e), e
        else:
            raise AssertionError("expected ValueError for a mismatched slide")


def test_annotation_sheet_accepts_alternate_fov_header():
    """Sheets label the FOV column FOVs/FOV/fov depending on who exported them."""
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color")]
    with tempfile.TemporaryDirectory() as d:
        sheet = os.path.join(d, "sheet.csv")
        _write_sheet(sheet,
                     [{"fov": "1", "UWA": "741"}, {"fov": "2", "UWA": "710"}],
                     fieldnames=["fov", "UWA"])
        rows = gsm.read_annotation_csv(sheet, SHEET_SLIDE, specs)
        assert [(r[1], r[2][0]) for r in rows] == [("1", "741"), ("2", "710")]


def test_annotation_sheet_without_fov_column_raises():
    """No FOV column means there is nothing to join on."""
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color")]
    with tempfile.TemporaryDirectory() as d:
        sheet = os.path.join(d, "sheet.csv")
        _write_sheet(sheet, [{"UWA": "741"}], fieldnames=["UWA"])
        try:
            gsm.read_annotation_csv(sheet, SHEET_SLIDE, specs)
        except ValueError as e:
            assert "FOV" in str(e), e
        else:
            raise AssertionError("expected ValueError when no FOV column is present")


def test_annotation_sheet_does_not_overwrite_atomx_values():
    """Where AtoMx already annotated a cell, the sheet leaves it alone."""
    specs = [gsm.ColumnSpec("UWA", "UWA", "UWA_color")]
    fields = ["cell_id", "fov", "cellSegmentationSetId", "UWA"]
    with tempfile.TemporaryDirectory() as d:
        gz = os.path.join(d, "meta.csv.gz")
        sheet = os.path.join(d, "sheet.csv")
        _write_gz(gz, fields, [
            {"cell_id": "c_1_1_1", "fov": "1", "cellSegmentationSetId": "s", "UWA": "existing"},
        ])
        _write_sheet(sheet, [{"fov": "1", "UWA": "741"}], fieldnames=["fov", "UWA"])
        rows, _ = gsm.read_rows(gz, None, specs)
        filled, _skipped = gsm.fill_missing_annotations(
            rows, gsm.read_annotation_csv(sheet, SHEET_SLIDE, specs), specs)
        assert filled == {}, filled
        assert rows[0][2][0] == "existing"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
