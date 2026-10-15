#!/usr/bin/env python3
"""Check merged Napari `_metadata.csv` files against the labels they were built from, by an independent join.

celltypes-to-napari-metadata.py joins the labels to Napari through the flat files, on (fov, cell_ID). This
checks that result a different way, so a mistake in one route cannot hide behind the same mistake in the other:
the pipeline index `<slide>_F<fov>_C<cell>` is rebuilt straight from Napari's own id `c_<n>_<fov>_<cell>` and the
label is looked up directly in the source (cell_id, label) CSV. It fails (exit 1) on:

  - rows or `cell_ID` order that differ from the original stitched file, or placeholder columns still present;
  - any column missing, or lacking its `<name>_color` pair;
  - a label that disagrees with the direct lookup (a cell with no source label must be blank);
  - a colour that is not the shared hash of its label, or a label with more than one colour.

A mismatch count EQUAL to the number of cells checked means the check is wrong, not the data (one earlier check
forgot a label prefix and flagged every cell), so that case prints a hint to suspect the check.

Usage:
    uv run python scripts/validate-napari-label-columns.py --merged-dir merged --original-dir stitched \\
        --column Donor=donor.csv --column SORL1_status=sorl1_status.csv --labels-dir napari_labels
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import duckdb
import pandas as pd

CELL_ID_COLUMN = "cell_ID"
COLOR_SUFFIX = "_color"
NAPARI_ID_RE = re.compile(r"^c_(?P<slide>\d+)_(?P<fov>\d+)_(?P<cell>\d+)$")
PLACEHOLDERS = ("cell_type", "hex_color")


def deterministic_color(value: str) -> str:
    """Same hash as generate-slide-metadata.py / celltypes-to-napari-metadata.py."""
    return duckdb.sql("SELECT printf('#%06X', abs(hash($1)) % 16777216)", params=[value]).fetchone()[0]


def pipeline_index(napari_ids: pd.Series, slide: str) -> pd.Series:
    """`c_<n>_<fov>_<cell>` -> `<slide>_F<fov>_C<cell>`; unparseable ids become NaN."""
    parts = napari_ids.astype(str).str.extract(NAPARI_ID_RE)
    return slide + "_F" + parts["fov"] + "_C" + parts["cell"]


def check_slide(slide: str, merged: pd.DataFrame, original: pd.DataFrame,
                sources: dict[str, dict[str, str]]) -> list[str]:
    """Return failure messages for one slide (empty = pass)."""
    fails = []
    if len(merged) != len(original) or not merged[CELL_ID_COLUMN].equals(original[CELL_ID_COLUMN]):
        fails.append(f"rows or {CELL_ID_COLUMN} order differ from the original ({len(original):,} -> {len(merged):,})")
    for col in PLACEHOLDERS:
        if col in merged.columns:
            fails.append(f"placeholder column '{col}' is still present")
    index = pipeline_index(merged[CELL_ID_COLUMN], slide)
    if index.isna().any():
        fails.append(f"{int(index.isna().sum()):,} cell ids do not look like c_<n>_<fov>_<cell>")
    for col, lookup in sources.items():
        if col not in merged.columns or col + COLOR_SUFFIX not in merged.columns:
            fails.append(f"column '{col}' or '{col}{COLOR_SUFFIX}' is missing")
            continue
        expected = index.map(lookup).fillna("")
        got = merged[col].fillna("").astype(str)
        wrong = int((got != expected).sum())
        if wrong:
            note = "  (equals the cells checked: suspect the CHECK, not the data)" if wrong == len(merged) else ""
            fails.append(f"{col}: {wrong:,} of {len(merged):,} labels disagree with the direct lookup{note}")
        colors = merged[col + COLOR_SUFFIX].fillna("").astype(str)
        per_label = pd.DataFrame({"label": got, "color": colors}).drop_duplicates()
        multi = per_label.groupby("label")["color"].nunique()
        if (multi > 1).any():
            fails.append(f"{col}: a label has more than one colour: {list(multi[multi > 1].index[:3])}")
        bad_color = [(label, color) for label, color in per_label.itertuples(index=False)
                     if label and color != deterministic_color(label)]
        if bad_color:
            fails.append(f"{col}: colour is not the shared hash for {bad_color[:3]}")
        if (colors[got == ""] != "").any():
            fails.append(f"{col}: blank labels carry a colour")
    return fails


def load_sources(labels_dir: Path, column_specs: list[str]) -> dict[str, dict[str, str]]:
    sources = {}
    for spec in column_specs:
        name, _, filename = spec.partition("=")
        if not name or not filename:
            sys.exit(f"ERROR: --column needs NAME=FILE, got {spec!r}")
        df = pd.read_csv(labels_dir / filename, dtype=str)
        if df["cell_id"].duplicated().any():
            sys.exit(f"ERROR: {filename} has duplicated cell ids.")
        sources[name] = dict(zip(df["cell_id"], df["label"]))
    return sources


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged-dir", type=Path, required=True, help="Directory of the finished <slide>_metadata.csv.")
    p.add_argument("--original-dir", type=Path, required=True, help="Directory of the unmodified stitched files.")
    p.add_argument("--labels-dir", type=Path, required=True, help="Directory of the (cell_id, label) CSVs.")
    p.add_argument("--column", action="append", default=[], required=True, metavar="NAME=FILE",
                   help="A Napari column and the label CSV it was built from. Repeatable.")
    args = p.parse_args()

    sources = load_sources(args.labels_dir, args.column)
    slides = sorted(f.name[: -len("_metadata.csv")] for f in args.merged_dir.glob("*_metadata.csv"))
    if not slides:
        sys.exit(f"ERROR: no *_metadata.csv in {args.merged_dir}")
    failed = 0
    for slide in slides:
        original_path = args.original_dir / f"{slide}_metadata.csv"
        if not original_path.exists():
            print(f"FAIL {slide}: no original at {original_path}")
            failed += 1
            continue
        merged = pd.read_csv(args.merged_dir / f"{slide}_metadata.csv", dtype=str, keep_default_na=False)
        original = pd.read_csv(original_path, dtype=str, keep_default_na=False)
        fails = check_slide(slide, merged, original, sources)
        labelled = {c: int((merged[c] != "").sum()) for c in sources if c in merged.columns}
        print(f"{'FAIL' if fails else 'ok  '} {slide:34s} {len(merged):7,d} cells  labelled: "
              + ", ".join(f"{c} {n:,}" for c, n in labelled.items()))
        for message in fails:
            print(f"       - {message}")
        failed += bool(fails)
    print(f"\n{len(slides) - failed} of {len(slides)} slides pass")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
