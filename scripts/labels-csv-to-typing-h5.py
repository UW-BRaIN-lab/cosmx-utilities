#!/usr/bin/env python3
"""Turn a per-cell label CSV into the small h5 that celltypes-to-napari-metadata.py reads.

That script takes its labels from a stage-4 InSituType result (`/cell_id`, `/cell_type`,
optionally `/prob`), which is the right input for a cell-typing column. Some of the labels
worth seeing on a slide are not cell-typing results at all -- the anchor cohort's Leiden
cluster is a column in a sidecar CSV (`anchor_cells.csv`: cell_id, slide_id, cluster) -- so
this converts any (id, label) pair of columns into the same three-dataset shape.

Doing it this way rather than hand-merging CSVs keeps every column on the one tested join.
The cell-id join between the pipeline index and Napari's own key is the part that is easy to
get quietly wrong, and it lives in celltypes-to-napari-metadata.py; this script deliberately
does nothing but reshape, so it cannot introduce a second way of getting it wrong.

Labels are written as strings. Integer cluster ids become "0", "1", ... which Napari will
colour but which read poorly in a legend, so --label-prefix turns them into e.g. "leiden_12".

Usage:
    uv run python scripts/labels-csv-to-typing-h5.py \\
        --csv anchor_cells.csv --label-column cluster --label-prefix leiden_ \\
        --output anchor_leiden.h5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", type=Path, required=True,
                   help="Per-cell CSV; one row per cell.")
    p.add_argument("--id-column", default="cell_id",
                   help="Column holding the pipeline cell index "
                        "'<slide>_F<fov>_C<cell_ID>' (default cell_id).")
    p.add_argument("--label-column", required=True,
                   help="Column holding the label to show in Napari.")
    p.add_argument("--label-prefix", default="",
                   help="Prepended to every label, so numeric cluster ids read as "
                        "'leiden_12' rather than '12'.")
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.csv)
    for col in (args.id_column, args.label_column):
        if col not in df.columns:
            sys.exit(f"ERROR: {args.csv} has no column {col!r}. Present: "
                     f"{', '.join(map(str, df.columns))}")

    df = df[[args.id_column, args.label_column]].dropna()
    if df.empty:
        sys.exit(f"ERROR: no rows left after dropping missing {args.id_column}/"
                 f"{args.label_column}.")
    dupes = int(df[args.id_column].duplicated().sum())
    if dupes:
        sys.exit(f"ERROR: {dupes:,} duplicated cell ids. One row per cell is required, or the "
                 f"Napari join would silently pick one arbitrarily.")

    ids = df[args.id_column].astype(str).to_numpy()
    labels = (args.label_prefix + df[args.label_column].astype(str)).to_numpy()

    # Fixed-length bytes truncate silently, which would corrupt cell ids and break the join
    # with no error, so size both datasets to what is actually present.
    id_width = max(len(s) for s in ids)
    label_width = max(len(s) for s in labels)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as f:
        f.create_dataset("cell_id", data=np.array(ids, dtype=f"S{id_width}"))
        f.create_dataset("cell_type", data=np.array(labels, dtype=f"S{label_width}"))

    print(f"{len(df):,} cells, {df[args.label_column].nunique()} distinct labels")
    print(f"  id width {id_width}, label width {label_width}")
    print(f"  example: {ids[0]} -> {labels[0]}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
