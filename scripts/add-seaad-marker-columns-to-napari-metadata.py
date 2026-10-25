#!/usr/bin/env python3
"""Add SEA-AD marker evidence to stitched Napari `_metadata.csv` files, so it can be checked on the slide.

`pipeline/reference/denovo_annotations/stage4_k21_seaad_markers.csv` says which SEA-AD MTG cell type each
K=21 cluster's marker genes match, and how strongly. This writes that as two columns next to the K=21
column already in each file, so Maddie can toggle between "what the cluster was called" and "what its
markers say" and judge whether the layer position fits:

  SEAAD_markers_K21          the best-matching SEA-AD type ("Oligodendrocyte", "L6b", ...), or
                             "No SEA-AD support"; coloured per SEA-AD type, so a de novo cluster and a
                             named Brain_6k type that map to the same SEA-AD type share a colour.
  SEAAD_markers_K21_support  strong / partial / none, the evidence table's `support` level.

This is evidence, not a rename: it is derived from marker genes, the Brain_6k named types are included
as controls (Brain_6k `L6 neuron` maps to Oligodendrocyte, `L4 neuron` to L5 IT), and nothing in the
existing columns is changed.

THE JOIN needs no cell ids. Each file already carries the K=21 call as text in --source-column
(`InSituType_K21`): the annotation string for a de novo cluster ("o - Oligodendrocyte") or the Brain_6k
type name for a named one. The evidence table lists exactly those strings (`current_annotation` for de novo
clusters, `cluster` for named types), so the mapping is a lookup on the column's own values. A value with
no entry in the table is an error, not a blank, because that would hide a stale or mismatched table; cells
with an empty source value stay empty and Napari draws them uncoloured.

Every value column gets a `<name>_color` pair. The colour function is the SAME hash
generate-slide-metadata.py and add-annotations-to-napari-metadata.py use, so a label keeps its colour in
every layer.

Expects `<slide>_metadata.csv` files in --metadata-dir (the stitched `_metadata.csv` of each slide
directory, copied out under the slide name). Existing files are never modified; the result goes to
--out-dir, and a file that already has either new column is refused.

Usage:
    uv run python scripts/add-seaad-marker-columns-to-napari-metadata.py \\
        --metadata-dir ./stitched_metadata --out-dir ./napari_metadata \\
        --evidence-table pipeline/reference/denovo_annotations/stage4_k21_seaad_markers.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd

CELL_ID_COLUMN = "cell_ID"
COLOR_SUFFIX = "_color"
METADATA_SUFFIX = "_metadata.csv"
DEFAULT_SOURCE_COLUMN = "InSituType_K21"
NAME_COLUMN = "SEAAD_markers_K21"
SUPPORT_COLUMN = "SEAAD_markers_K21_support"
NO_SUPPORT_LABEL = "No SEA-AD support"
REFERENCE_PREFIX = "SEAAD_"
EVIDENCE_COLUMNS = {"cluster", "kind", "support", "proposed_seaad_name", "current_annotation"}
DE_NOVO, NAMED = "de novo", "named"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata-dir", type=Path, required=True,
                   help="Directory of <slide>_metadata.csv.")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--evidence-table", type=Path, required=True,
                   help="stage4_k21_seaad_markers.csv from marker_overlap_annotation_table.py.")
    p.add_argument("--source-column", default=DEFAULT_SOURCE_COLUMN,
                   help=f"Column holding the K=21 call (default {DEFAULT_SOURCE_COLUMN}).")
    return p.parse_args()


def deterministic_color(value: str) -> str:
    """Deterministic hex colour -- identical to generate-slide-metadata.py."""
    return duckdb.sql("SELECT printf('#%06X', abs(hash($1)) % 16777216)",
                      params=[value]).fetchone()[0]


def load_mapping(evidence_table: Path) -> pd.DataFrame:
    """source value -> (SEA-AD name, support), keyed the way the K=21 column spells its values."""
    table = pd.read_csv(evidence_table, dtype=str, keep_default_na=False)
    missing = EVIDENCE_COLUMNS - set(table.columns)
    if missing:
        sys.exit(f"ERROR: {evidence_table} lacks columns {sorted(missing)}")
    table["source_value"] = table["current_annotation"].where(table["kind"] == DE_NOVO, table["cluster"])
    if (table["source_value"] == "").any():
        blank = table.loc[table["source_value"] == "", "cluster"].tolist()
        sys.exit(f"ERROR: de novo cluster(s) {blank} have no current_annotation to match the metadata on.")
    dup = table["source_value"][table["source_value"].duplicated()].tolist()
    if dup:
        sys.exit(f"ERROR: the evidence table maps these source values more than once: {dup}")
    table[NAME_COLUMN] = (table["proposed_seaad_name"].str.removeprefix(REFERENCE_PREFIX)
                          .where(table["proposed_seaad_name"] != "", NO_SUPPORT_LABEL))
    table[SUPPORT_COLUMN] = table["support"]
    return table.set_index("source_value")[[NAME_COLUMN, SUPPORT_COLUMN]]


def add_columns(meta: pd.DataFrame, source_column: str, mapping: pd.DataFrame, slide: str) -> pd.DataFrame:
    source = meta[source_column].astype(str)
    unknown = sorted(set(source[source != ""]) - set(mapping.index))
    if unknown:
        sys.exit(f"ERROR: {slide}: {source_column} has value(s) the evidence table does not cover: "
                 f"{unknown}. The table and the metadata are from different runs.")
    out = meta.copy()
    for column in (NAME_COLUMN, SUPPORT_COLUMN):
        out[column] = source.map(mapping[column]).fillna("").to_numpy()
        colors = {v: deterministic_color(v) for v in sorted(set(out[column])) if v}
        out[column + COLOR_SUFFIX] = out[column].map(colors).fillna("").to_numpy()
    return out


def main() -> None:
    args = parse_args()
    mapping = load_mapping(args.evidence_table)
    files = sorted(args.metadata_dir.glob(f"*{METADATA_SUFFIX}"))
    if not files:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} files in {args.metadata_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"{len(files)} slide file(s); adding {NAME_COLUMN} and {SUPPORT_COLUMN} from {args.source_column}\n")
    print(f"  {'slide':36s} {'cells':>9s} {'empty source':>13s}")
    for path in files:
        slide = path.name[: -len(METADATA_SUFFIX)]
        # keep_default_na=False: text such as "NA" or "None" in an existing column must come back as written.
        meta = pd.read_csv(path, dtype=str, keep_default_na=False)
        for required in (CELL_ID_COLUMN, args.source_column):
            if required not in meta.columns:
                sys.exit(f"ERROR: {path} has no '{required}' column; columns are {list(meta.columns)}")
        clash = [c for c in (NAME_COLUMN, SUPPORT_COLUMN) if c in meta.columns]
        if clash:
            sys.exit(f"ERROR: {path.name} already has {clash}. Refusing to overwrite.")
        out = add_columns(meta, args.source_column, mapping, slide)
        if len(out) != len(meta) or not out[CELL_ID_COLUMN].equals(meta[CELL_ID_COLUMN]):
            sys.exit(f"ERROR: rows or cell_ID changed for {slide}.")
        out.to_csv(args.out_dir / path.name, index=False)
        n_empty = int((out[NAME_COLUMN] == "").sum())
        print(f"  {slide:36s} {len(out):>9,d} {n_empty:>13,d}")
    print(f"\nWrote {len(files)} file(s) to {args.out_dir}")


if __name__ == "__main__":
    main()
