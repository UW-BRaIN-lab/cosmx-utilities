"""Shared helpers for coloring cells by a metadata column.

Kept dependency-light (pandas only) and free of napari imports so both
``gemini.py`` and ``_dock_widget.py`` can use it without a circular import.
"""
import pandas as pd

# Metadata columns that are identifiers/coordinates, never something to color by.
_NON_COLORABLE = ("cell_ID", "fov", "CellId")

# Suffix marking a companion color column (e.g. ``celltype_norefit_color``);
# ``hex_color`` also ends with this suffix, so both are hidden from color-by menus.
COLOR_COLUMN_SUFFIX = "_color"


def colorable_columns(columns):
    """Metadata columns a user may color cells by: drops identifiers and the
    ``*_color`` helper columns."""
    return [
        c for c in columns
        if c not in _NON_COLORABLE and not c.endswith(COLOR_COLUMN_SUFFIX)
    ]


def categorical_color_map(metadata, col_name):
    """Return ``{category value: hex color string}`` for a categorical column.

    Prefers a per-column ``<col_name>_color`` column (multi-annotation
    ``_metadata.csv``), then a generic ``hex_color`` column (legacy single
    cell-type files). Returns ``None`` when neither is present, signalling the
    caller to auto-generate colors. NaN category values are skipped.
    """
    for color_col in (f"{col_name}{COLOR_COLUMN_SUFFIX}", "hex_color"):
        if color_col in metadata.columns:
            pairs = metadata[[col_name, color_col]].drop_duplicates(subset=[col_name])
            return {
                row[col_name]: row[color_col]
                for _, row in pairs.iterrows()
                if pd.notna(row[col_name])
            }
    return None
