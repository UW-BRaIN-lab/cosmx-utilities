#!/usr/bin/env python3
"""Make every colour in existing Napari `_metadata.csv` files legible on a black background.

Applies napari_colors.legible_color to each `<name>_color` column (and the generic `hex_color`): a
colour that is too dark is lightened with its hue kept, every other colour is left alone. Labels and
all other columns are untouched, so the files stay the same shape and the change is only colour.
Every colour column is corrected, not just the ones that look dark, so a label keeps one colour in
every column it appears in. Safe to run twice.

Writes the corrected files to --out-dir (same file names) and a CSV listing each changed colour.

Usage:
    uv run --with pandas python scripts/recolor-napari-metadata.py \\
        --in-dir live --out-dir recolored --changes-csv recolor_changes.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from napari_colors import contrast_against_black, legible_color  # noqa: E402

COLOR_SUFFIX = "_color"
GENERIC_COLOR_COLUMN = "hex_color"
METADATA_SUFFIX = "_metadata.csv"


def color_columns(columns) -> list[str]:
    return [c for c in columns if c.endswith(COLOR_SUFFIX) or c == GENERIC_COLOR_COLUMN]


def label_column_for(color_column: str, columns) -> str:
    """The label column a colour column belongs to (`InSituTree_v5_color` -> `InSituTree_v5`)."""
    if color_column == GENERIC_COLOR_COLUMN:
        return "cell_type" if "cell_type" in columns else ""
    stem = color_column[: -len(COLOR_SUFFIX)]
    return stem if stem in columns else ""


def recolor_frame(frame: pd.DataFrame, source: str) -> tuple[pd.DataFrame, list[dict]]:
    """A corrected copy of `frame` and one record per (column, old colour) that changed."""
    out = frame.copy()
    changes: list[dict] = []
    for column in color_columns(frame.columns):
        mapping = {old: legible_color(old) for old in frame[column].dropna().unique() if isinstance(old, str)}
        label_col = label_column_for(column, frame.columns)
        for old, new in mapping.items():
            if new == old:
                continue
            rows = frame[column] == old
            labels = sorted(frame.loc[rows, label_col].dropna().astype(str).unique()) if label_col else []
            changes.append({"file": source, "column": column, "labels": " | ".join(labels), "old": old, "new": new,
                            "old_contrast": round(contrast_against_black(old), 2),
                            "new_contrast": round(contrast_against_black(new), 2), "cells": int(rows.sum())})
        out[column] = frame[column].map(lambda v: mapping.get(v, v))
    return out, changes


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", type=Path, required=True, help=f"Folder of *{METADATA_SUFFIX} files.")
    p.add_argument("--out-dir", type=Path, required=True, help="Where the corrected files go (same names).")
    p.add_argument("--changes-csv", type=Path, required=True, help="Listing of every colour that changed.")
    args = p.parse_args()

    files = sorted(args.in_dir.glob(f"*{METADATA_SUFFIX}"))
    if not files:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} files in {args.in_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_changes: list[dict] = []
    for path in files:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])
        frame = frame.replace({"": None})
        corrected, changes = recolor_frame(frame, path.name)
        corrected.to_csv(args.out_dir / path.name, index=False, na_rep="")
        all_changes += changes
        print(f"{path.name}: {len(changes)} colours changed in {len({c['column'] for c in changes})} columns")
    pd.DataFrame(all_changes).to_csv(args.changes_csv, index=False)
    print(f"Wrote {len(files)} files to {args.out_dir} and {len(all_changes)} changes to {args.changes_csv}")


if __name__ == "__main__":
    main()
