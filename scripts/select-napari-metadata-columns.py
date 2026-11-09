#!/usr/bin/env python3
"""Keep only the chosen label columns in stitched Napari `_metadata.csv` files.

Napari-CosMx puts EVERY value column in the hover tooltip. Each typing pass adds one, so a
file that has accumulated the AtoMx types, the Leiden partition, the niches and several of our
own fits produces a tooltip too long to read -- and the reviewer then cannot tell which colour
belongs to the fit they were asked to look at. Dropping the columns a given review does not
need is what makes the tooltip usable again.

This runs on the FINISHED metadata files, not on a typing result: the expensive per-cell joins
have already happened, and a review set is just a projection of them. Keep the full file --
it is the one that can be re-projected for the next question.

Each value column is paired with `<name>_color`, which napari-cosmx expects, so the pair is
kept or dropped together. A kept column whose `_color` partner is missing is an ERROR, not a
warning: napari-cosmx falls back to the file-wide `hex_color` for any column without its own,
which renders that column as a single colour over the whole slide -- a plausible-looking
picture carrying no information.

`cell_ID` is always kept and always first; it is the key napari-cosmx joins on.

Usage:
    uv run --with pandas python scripts/select-napari-metadata-columns.py \\
        --metadata-dir ./napari_metadata_out --out-dir ./napari_metadata_review \\
        --columns insitutree_75,insitutree_66,denovo_split,Region
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

CELL_ID_COLUMN = "cell_ID"
COLOR_SUFFIX = "_color"
METADATA_SUFFIX = "_metadata.csv"


def color_column(column: str) -> str:
    return f"{column}{COLOR_SUFFIX}"


def columns_to_keep(available: list[str], requested: list[str], source: Path) -> list[str]:
    """`cell_ID`, then each requested column followed by its `_color`, validated.

    Both failure modes exit rather than degrade: a requested column that is absent means the
    caller is thinking of a different file, and a kept column without its `_color` would be
    drawn in one flat colour.
    """
    if CELL_ID_COLUMN not in available:
        sys.exit(f"ERROR: {source} has no '{CELL_ID_COLUMN}' column")

    absent = [c for c in requested if c not in available]
    if absent:
        value_columns = [c for c in available
                         if c != CELL_ID_COLUMN and not c.endswith(COLOR_SUFFIX)]
        sys.exit(f"ERROR: {source} has no column(s) {absent}.\n"
                 f"  Columns it does have: {value_columns}")

    uncoloured = [c for c in requested if color_column(c) not in available]
    if uncoloured:
        sys.exit(f"ERROR: {source} has no {COLOR_SUFFIX} partner for {uncoloured}.\n"
                 f"  napari-cosmx would fall back to the file-wide hex_color and draw each "
                 f"of them as one colour over the whole slide.")

    kept = [CELL_ID_COLUMN]
    for column in requested:
        kept += [column, color_column(column)]
    return kept


def select_columns(source: Path, destination: Path, requested: list[str]) -> tuple[int, list[str]]:
    """Write `source` to `destination` with only `requested` (plus colours). Returns (rows, dropped)."""
    # dtype=str / keep_default_na=False: labels and hex colours pass through untouched, and a
    # blank stays blank instead of becoming the string "nan" on the way out.
    available = list(pd.read_csv(source, nrows=0).columns)
    kept = columns_to_keep(available, requested, source)
    dropped = [c for c in available if c not in kept and not c.endswith(COLOR_SUFFIX)]

    table = pd.read_csv(source, usecols=kept, dtype=str, keep_default_na=False)
    table[kept].to_csv(destination, index=False)
    return len(table), dropped


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata-dir", type=Path, required=True,
                   help=f"Directory of stitched <slide>{METADATA_SUFFIX} files.")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Where the trimmed copies are written, under the same names.")
    p.add_argument("--columns", required=True,
                   help="Comma-separated value columns to KEEP, in tooltip order. Their "
                        f"'{COLOR_SUFFIX}' partners come along automatically; "
                        f"'{CELL_ID_COLUMN}' is always kept.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    requested = [c.strip() for c in args.columns.split(",") if c.strip()]
    if not requested:
        sys.exit("ERROR: --columns named no columns")

    sources = sorted(args.metadata_dir.glob(f"*{METADATA_SUFFIX}"))
    if not sources:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} in {args.metadata_dir}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Keeping {CELL_ID_COLUMN} + {requested} in {len(sources)} file(s)")

    reported_dropped = False
    for source in sources:
        rows, dropped = select_columns(source, args.out_dir / source.name, requested)
        if not reported_dropped:
            print(f"Dropping {len(dropped)} column(s): {dropped}")
            reported_dropped = True
        print(f"  {source.name}: {rows:,} cells")

    print(f"Wrote {len(sources)} file(s) to {args.out_dir}")
    print(f"The full files stay in {args.metadata_dir} -- re-run this for a different review.")


if __name__ == "__main__":
    main()
