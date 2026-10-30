#!/usr/bin/env python3
"""Rename one label in chosen columns of existing Napari `_metadata.csv` files.

Every cell carrying --old in a named column gets --new, and that column's `<column>_color` for those
cells becomes the colour the shared rule (napari_colors.deterministic_color) gives --new, so the file
still passes validate-napari-label-columns.py. Nothing else changes. A column missing from a file, or
--new already present in the column, is an error: the rename would silently merge two labels.

Usage:
    uv run --with pandas --with duckdb python scripts/rename-napari-label.py \\
        --in-dir live --out-dir renamed --column InSituTree_v5 --column InSituTree_v5_collapsed \\
        --old Excitatory_upper --new Excitatory_L2_3_mixed
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from napari_colors import deterministic_color  # noqa: E402

COLOR_SUFFIX = "_color"
METADATA_SUFFIX = "_metadata.csv"


def rename_label(frame: pd.DataFrame, columns: list[str], old: str, new: str) -> tuple[pd.DataFrame, dict[str, int]]:
    """A copy of `frame` with `old` -> `new` in each of `columns` (colour column updated); cells renamed per column."""
    out = frame.copy()
    renamed: dict[str, int] = {}
    new_color = deterministic_color(new)
    for column in columns:
        color_column = column + COLOR_SUFFIX
        missing = [c for c in (column, color_column) if c not in frame.columns]
        if missing:
            raise KeyError(f"columns missing: {missing}")
        if (frame[column] == new).any():
            raise ValueError(f"{column!r} already holds {new!r}; renaming {old!r} would merge two labels")
        rows = frame[column] == old
        out.loc[rows, column] = new
        out.loc[rows, color_column] = new_color
        renamed[column] = int(rows.sum())
    return out, renamed


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", type=Path, required=True, help=f"Folder of *{METADATA_SUFFIX} files.")
    p.add_argument("--out-dir", type=Path, required=True, help="Where the renamed files go (same names).")
    p.add_argument("--column", action="append", required=True, help="A label column to rename in. Repeatable.")
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    args = p.parse_args()

    files = sorted(args.in_dir.glob(f"*{METADATA_SUFFIX}"))
    if not files:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} files in {args.in_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for path in files:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])
        try:
            out, renamed = rename_label(frame, args.column, args.old, args.new)
        except (KeyError, ValueError) as err:
            sys.exit(f"ERROR in {path.name}: {err}")
        out.to_csv(args.out_dir / path.name, index=False)
        print(f"{path.name}: " + ", ".join(f"{c}: {n:,} cells" for c, n in renamed.items()))
    print(f"Wrote {len(files)} files to {args.out_dir}")


if __name__ == "__main__":
    main()
