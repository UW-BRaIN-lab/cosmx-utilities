"""Tests for napari_cosmx._colors (pandas-only; no napari/Qt needed).

Runnable under pytest or directly:
    uv run python napari-cosmx-fork/src/napari_cosmx/_tests/test_colors.py
"""
import io
import numpy as np
import pandas as pd

from napari_cosmx._colors import categorical_color_map, colorable_columns


def _multi_col_df():
    return pd.DataFrame({
        "cell_ID": ["c_1_1_1", "c_1_1_2", "c_1_1_3"],
        "fov": [1, 1, 1],
        "celltype_norefit": ["astrocyte", "microglia", "rod"],
        "celltype_norefit_color": ["#111111", "#222222", "#333333"],
        "Region": ["retina", "retina", "brain"],
        "Region_color": ["#aaaaaa", "#aaaaaa", "#bbbbbb"],
    })


def _legacy_df():
    return pd.DataFrame({
        "cell_ID": ["c_1_1_1", "c_1_1_2"],
        "cell_type": ["astrocyte", "microglia"],
        "hex_color": ["#111111", "#222222"],
    })


def test_colorable_columns_hides_ids_and_color_helpers():
    cols = colorable_columns(_multi_col_df().columns)
    assert cols == ["celltype_norefit", "Region"], cols
    # legacy hex_color is also hidden (ends with _color)
    assert colorable_columns(_legacy_df().columns) == ["cell_type"]


def test_per_column_color_preferred():
    df = _multi_col_df()
    assert categorical_color_map(df, "celltype_norefit") == {
        "astrocyte": "#111111", "microglia": "#222222", "rod": "#333333"}
    assert categorical_color_map(df, "Region") == {
        "retina": "#aaaaaa", "brain": "#bbbbbb"}


def test_legacy_hex_color_fallback():
    df = _legacy_df()
    assert categorical_color_map(df, "cell_type") == {
        "astrocyte": "#111111", "microglia": "#222222"}


def test_returns_none_when_no_color_column():
    df = pd.DataFrame({"cell_ID": ["c_1_1_1"], "Region": ["retina"]})
    assert categorical_color_map(df, "Region") is None


def test_nan_category_skipped():
    df = pd.DataFrame({
        "cell_ID": ["a", "b"],
        "Region": ["retina", np.nan],
        "Region_color": ["#aaaaaa", "#cccccc"],
    })
    result = categorical_color_map(df, "Region")
    assert result == {"retina": "#aaaaaa"}


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")


def test_blank_annotation_values_do_not_break_sorting():
    """A column with unassigned cells must still produce a sorted category list.

    AtoMx types only QC-passing cells, so a cell-type column legitimately has
    blanks. pandas reads a blank field as NaN, NaN is a float, and
    `sorted(np.unique(column))` then sorts str against float and raises
    TypeError. That fired while building the widget, so the viewer died before
    opening rather than degrading.
    """
    import numpy as np
    import pandas as pd

    # Two columns, so a blank is an empty *field* -- a blank line would be
    # skipped by read_csv and would not reproduce anything.
    csv_text = (
        "cell_ID,Cell Type\n"
        "c_1_1_1,Microglia.A\n"
        "c_1_1_2,\n"
        "c_1_1_3,Astrocyte.B\n"
        "c_1_1_4,\n"
        "c_1_1_5,Endothelial\n"
    )
    values = pd.read_csv(io.StringIO(csv_text))["Cell Type"]
    assert values.isna().sum() == 2, "fixture must contain blanks to be meaningful"

    try:
        sorted(np.unique(values))
    except TypeError as e:
        assert "not supported between instances" in str(e), e
    else:
        raise AssertionError("fixture no longer reproduces the original crash")

    # What the widget does now.
    assert sorted(pd.unique(values.dropna())) == [
        "Astrocyte.B", "Endothelial", "Microglia.A"]
