#!/usr/bin/env python3
"""Tests for scripts/cell-type-histograms.py core logic (no S3 / network).

Runnable either under pytest or directly:  uv run python scripts/tests/test_cell_type_histograms.py
"""
import importlib.util
import tempfile
from collections import Counter
from pathlib import Path

# Import the hyphenated script by file path.
_SCRIPT = Path(__file__).resolve().parents[1] / "cell-type-histograms.py"
_spec = importlib.util.spec_from_file_location("cell_type_histograms", _SCRIPT)
cth = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cth)


class _FakeBody:
    def __init__(self, text):
        self._text = text

    def read(self):
        return self._text.encode("utf-8")


class _FakeS3:
    """Stands in for boto3's S3 client, serving one metadata body per key."""

    def __init__(self, bodies):
        self.bodies = bodies
        self.requested_keys = []

    def get_object(self, Bucket, Key):
        self.requested_keys.append(Key)
        if Key not in self.bodies:
            raise KeyError(Key)
        return {"Body": _FakeBody(self.bodies[Key])}


BUCKET = "test-bucket"
PREFIX = "napari-stitched/Study/Run"
SLIDE = "slide_1"
KEY = f"{PREFIX}/{SLIDE}/_metadata.csv"

# Multi-column metadata as generate-slide-metadata.py writes it: a space in the
# column name, and each annotation paired with a <name>_color.
MULTI_COLUMN_CSV = (
    "cell_ID,Cell Type,Cell Type_color,Case Broad,Case Broad_color\n"
    "c_1_1_1,astrocyte,#111111,AD,#AAAAAA\n"
    "c_1_1_2,astrocyte,#111111,AD,#AAAAAA\n"
    "c_1_1_3,microglia,#222222,Control,#BBBBBB\n"
)

LEGACY_CSV = (
    "cell_ID,cell_type,hex_color\n"
    "c_1_1_1,astrocyte,#111111\n"
    "c_1_1_2,microglia,#222222\n"
)


def _load(csv_text, columns):
    s3 = _FakeS3({KEY: csv_text})
    return cth.load_metadata(s3, BUCKET, PREFIX, SLIDE, columns)


def test_slugify_makes_spaced_column_filename_safe():
    assert cth.slugify("Cell Type") == "Cell_Type"
    assert cth.slugify("SORL1 Mutation") == "SORL1_Mutation"
    assert cth.slugify("cell_type") == "cell_type"


def test_slugify_never_returns_empty():
    assert cth.slugify("///") == "column"


def test_deterministic_color_is_stable_and_value_specific():
    assert cth.deterministic_color("astrocyte") == cth.deterministic_color("astrocyte")
    assert cth.deterministic_color("astrocyte") != cth.deterministic_color("microglia")
    assert cth.deterministic_color("astrocyte").startswith("#")
    assert len(cth.deterministic_color("astrocyte")) == 7


def test_annotation_columns_finds_only_color_paired_columns():
    headers = ["cell_ID", "Cell Type", "Cell Type_color", "Case Broad", "Case Broad_color"]
    assert cth.annotation_columns(headers) == ["Cell Type", "Case Broad"]


def test_counts_a_spaced_multi_column_annotation():
    results = _load(MULTI_COLUMN_CSV, ["Cell Type"])
    counts, colors = results["Cell Type"]
    assert counts == Counter({"astrocyte": 2, "microglia": 1})
    assert colors == {"astrocyte": "#111111", "microglia": "#222222"}


def test_counts_several_columns_in_one_pass():
    """Two columns must cost one GET, not two: the files run to tens of MB."""
    s3 = _FakeS3({KEY: MULTI_COLUMN_CSV})
    results = cth.load_metadata(s3, BUCKET, PREFIX, SLIDE, ["Cell Type", "Case Broad"])
    assert s3.requested_keys == [KEY], s3.requested_keys
    assert results["Cell Type"][0] == Counter({"astrocyte": 2, "microglia": 1})
    assert results["Case Broad"][0] == Counter({"AD": 2, "Control": 1})
    assert results["Case Broad"][1] == {"AD": "#AAAAAA", "Control": "#BBBBBB"}


def test_legacy_single_column_metadata_still_works():
    """The pre-multi-column layout: cell_type paired with a generic hex_color."""
    results = _load(LEGACY_CSV, [cth.LEGACY_COLUMN])
    counts, colors = results[cth.LEGACY_COLUMN]
    assert counts == Counter({"astrocyte": 1, "microglia": 1})
    assert colors == {"astrocyte": "#111111", "microglia": "#222222"}


def test_missing_column_raises_and_names_available_columns():
    """The failure that motivated --column: a silent empty plot is worse than a stop."""
    try:
        _load(MULTI_COLUMN_CSV, ["cell_type"])
    except cth.MissingColumnError as e:
        message = str(e)
        assert "cell_type" in message, message
        assert "Cell Type" in message, message
        assert "Case Broad" in message, message
    else:
        raise AssertionError("expected MissingColumnError for an absent column")


def test_absent_color_column_falls_back_to_deterministic_color():
    csv_text = "cell_ID,Cell Type\nc_1_1_1,astrocyte\n"
    counts, colors = _load(csv_text, ["Cell Type"])["Cell Type"]
    assert counts == Counter({"astrocyte": 1})
    assert colors == {"astrocyte": cth.deterministic_color("astrocyte")}


def test_blank_values_are_not_counted():
    """Cells the annotation never reached must not become an empty-string bar."""
    csv_text = (
        "cell_ID,Cell Type,Cell Type_color\n"
        "c_1_1_1,astrocyte,#111111\n"
        "c_1_1_2,,\n"
        "c_1_1_3,   ,\n"
    )
    counts, colors = _load(csv_text, ["Cell Type"])["Cell Type"]
    assert counts == Counter({"astrocyte": 1}), counts
    assert "" not in colors


def test_values_are_stripped_of_surrounding_whitespace():
    csv_text = (
        "cell_ID,Cell Type,Cell Type_color\n"
        "c_1_1_1,astrocyte,#111111\n"
        "c_1_1_2,  astrocyte  ,#111111\n"
    )
    counts, _ = _load(csv_text, ["Cell Type"])["Cell Type"]
    assert counts == Counter({"astrocyte": 2}), counts


def test_plot_histogram_writes_a_png():
    counts = Counter({"astrocyte": 2, "microglia": 1})
    colors = {"astrocyte": "#111111", "microglia": "#222222"}
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "out.png"
        cth.plot_histogram(SLIDE, "Cell Type", counts, colors, str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
