"""Tests for napari_cosmx._colors (pandas-only; no napari/Qt needed).

Runnable under pytest or directly:
    uv run python napari-cosmx-fork/src/napari_cosmx/_tests/test_colors.py
"""
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
