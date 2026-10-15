#!/usr/bin/env python3
"""Turn the SORL1 per-cell labels table into one (cell_id, label) CSV per Napari dropdown column.

Input is `cell_labels.csv.gz` written by pipeline/slurm/54_export_cell_labels.sh: one row per cell, keyed
by the pipeline index `<slide>_F<fov>_C<cell_ID>`. Each output CSV feeds labels-csv-to-typing-h5.py, and from
there the tested join in celltypes-to-napari-metadata.py, so every column rides the same join.

Columns written (file -> Napari column it is meant to become):
  donor.csv          Donor             "UWA_741"
  sorl1_status.csv   SORL1_status      "SORL1 R953C" / "SORL1 C467Y" / "No SORL1 variant"
  insitutype_k21.csv InSituType_K21    reference types by name; de novo letters with their reading
                                       ("o - Oligodendrocyte"), taken from the K=21 annotation table
  insitutree.csv     InSituTree        the tree's own leaf / group label for the run given
  opc_confident.csv  OPC_confidence    only for cells the tree called OPC: "OPC (confident)" or
                                       "OPC (low confidence)"; every other cell is absent (blank in Napari)

The SORL1 status comes from `sorl1_mutation`, and is checked against `case_group`: a disagreement
between the two columns stops the run, because a wrong carrier label on a slide is worse than none.
The OPC flag is read from the table (column `opc_confident_<run>`, made by the export job) and never
recomputed here, so there is one definition of it.

Usage:
    uv run python scripts/prepare-sorl1-napari-labels.py --labels cell_labels.csv.gz \\
        --k21-annotations stage4_k21.csv --tree-run tree_v2 --out-dir napari_labels

Run order for a stitched study (all local; `stitched/` holds each slide's bare `_metadata.csv` copied to
`<slide>_metadata.csv`; `flat/` holds the export's `<slide>_metadata_file.csv.gz`):
  1. this script                                 -> napari_labels/*.csv
  2. per column: labels-csv-to-typing-h5.py --csv napari_labels/X.csv --label-column label --output X.h5
  3. per column, chained with --merge-into the previous step's directory:
       celltypes-to-napari-metadata.py --typing-h5 X.h5 --flatfiles flat --column NAME --merge-into <prev> --out-dir stepN
     with `--drop-columns cell_type,hex_color` on the FIRST step (the stitcher's placeholders)
  4. validate-napari-label-columns.py --merged-dir <last step> --original-dir stitched --labels-dir napari_labels --column ...
  5. back up each `_metadata.csv` server-side to `_metadata.csv.bak` (check it is absent first), upload, download back and
     compare md5. Never re-stitch a slide afterwards: the stitcher replaces its whole output directory.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

DONOR_PREFIX = "UWA_"
NO_VARIANT = "No SORL1 variant"
VARIANT_RE = re.compile(r"SORL1\s+([A-Z]\d+[A-Z])")
DENOVO_RE = re.compile(r"^[a-z]{1,2}$")
OPC_CONFIDENT = "OPC (confident)"
OPC_LOW = "OPC (low confidence)"
OUTPUTS = {"donor": "donor.csv", "sorl1_status": "sorl1_status.csv", "insitutype_k21": "insitutype_k21.csv",
           "insitutree": "insitutree.csv", "opc_confidence": "opc_confident.csv"}


def sorl1_status(mutation: pd.Series, case_group: pd.Series) -> pd.Series:
    """Variant-level SORL1 status, checked against the independent case_group column."""
    variant = mutation.astype(str).str.extract(VARIANT_RE)[0]
    status = variant.map(lambda v: f"SORL1 {v}" if isinstance(v, str) else NO_VARIANT)
    carrier_by_mutation = variant.notna()
    carrier_by_group = case_group.astype(str).str.contains("SORL1", na=False)
    disagree = carrier_by_mutation != carrier_by_group
    if disagree.any():
        raise ValueError(f"{int(disagree.sum()):,} cells where sorl1_mutation and case_group disagree on "
                         f"carrier status, e.g. {mutation[disagree].iloc[0]!r} vs {case_group[disagree].iloc[0]!r}")
    return status


def k21_names(labels: pd.Series, annotation_csv: Path) -> pd.Series:
    """Readable K=21 labels: de novo letters become the annotation text ('o - Oligodendrocyte')."""
    ann = pd.read_csv(annotation_csv, dtype=str)
    reading = dict(zip(ann["denovo_label"].str.strip(), ann["annotation"].str.strip()))
    letters = {x for x in labels.dropna().unique() if DENOVO_RE.match(str(x))}
    missing = sorted(letters - set(reading))
    if missing:
        raise ValueError(f"de novo labels with no row in {annotation_csv.name}: {missing}")
    return labels.map(lambda x: reading.get(x, x) if DENOVO_RE.match(str(x)) else x)


def opc_confidence(table: pd.DataFrame, run: str, opc_labels: list[str]) -> pd.Series:
    flag_col = f"opc_confident_{run}"
    if flag_col not in table.columns:
        raise ValueError(f"column {flag_col!r} is missing from the labels table. Rerun the export job with "
                         f"OPC_FLAG=\"{run}=OPC\" (pipeline/slurm/54_export_cell_labels.sh) and pull the new table.")
    called = table[f"celltype_{run}"].isin(opc_labels)
    out = pd.Series(pd.NA, index=table.index, dtype="object")
    out[called & table[flag_col].astype(bool)] = OPC_CONFIDENT
    out[called & ~table[flag_col].astype(bool)] = OPC_LOW
    return out


def build_label_frames(table: pd.DataFrame, k21_annotation_csv: Path, tree_run: str,
                       opc_labels: list[str]) -> dict[str, pd.Series]:
    return {
        "donor": DONOR_PREFIX + table["donor"].astype(str),
        "sorl1_status": sorl1_status(table["sorl1_mutation"], table["case_group"]),
        "insitutype_k21": k21_names(table["celltype_k21"], k21_annotation_csv),
        "insitutree": table[f"celltype_{tree_run}"],
        "opc_confidence": opc_confidence(table, tree_run, opc_labels),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--labels", type=Path, required=True, help="cell_labels.csv.gz from the export job.")
    p.add_argument("--k21-annotations", type=Path, required=True, help="denovo_annotations/stage4_k21.csv")
    p.add_argument("--tree-run", default="tree_v2", help="InSituTree run name in the labels table.")
    p.add_argument("--opc-labels", default="OPC", help="Comma-separated tree labels that mean OPC.")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    table = pd.read_csv(args.labels, index_col="cell_id", dtype={"donor": str})
    if table.index.duplicated().any():
        sys.exit("ERROR: duplicated cell ids in the labels table.")
    try:
        frames = build_label_frames(table, args.k21_annotations, args.tree_run, args.opc_labels.split(","))
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for key, series in frames.items():
        out = series.dropna().rename("label").rename_axis("cell_id").reset_index()
        out.to_csv(args.out_dir / OUTPUTS[key], index=False)
        print(f"{OUTPUTS[key]:20s} {len(out):8,d} cells, {out['label'].nunique():3d} labels")
    print(f"Wrote {len(frames)} label files to {args.out_dir}")


if __name__ == "__main__":
    main()
