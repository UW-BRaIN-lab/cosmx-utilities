#!/usr/bin/env python3
"""Add Region / Case / Block to stitched Napari `_metadata.csv` files, from the curated FOV sheet.

The stitcher writes `_metadata.csv` from whatever cell typing AtoMx supplied. With none, the
file holds only a placeholder `cell_type` of "Unassigned" and a generic `hex_color`, so there is
nothing to colour by. Region -- Tumor bulk / Infiltrating edge / Contralateral uninvolved -- is
the context that makes a scatter of typed cells readable on a slide, and it is already curated
per FOV in the AtoMx annotation reference, so it can be added without touching the flat files.

THE JOIN. Napari's `cell_ID` is `c_<slide#>_<FOV>_<cell>`, so the FOV is recoverable from the id
itself. The annotation reference keys on canonical slide names while the stitched directories
use AtoMx export folder names, hence the crosswalk. A slide carries two tissue pieces from two
different cases, so Case is genuinely informative: it is what says which half is which donor.

Every annotation column is paired with a `<name>_color`, which napari-cosmx expects for each
value column. The colour function is the SAME hash generate-slide-metadata.py and
celltypes-to-napari-metadata.py use, so a category keeps its colour whichever script wrote it.

Cells whose FOV has no annotation get an empty value, so Napari still draws them uncoloured,
and the count is reported per slide -- silently dropping them would hide a broken join.

Expects `<slide>_metadata.csv` files in --metadata-dir, named by the AtoMx flat-file folder
(e.g. 7134A77439A6_metadata.csv). Stitched output stores a bare `_metadata.csv` inside each
slide directory, so copy those out under the slide name first.

Usage:
    uv run python scripts/add-annotations-to-napari-metadata.py \\
        --metadata-dir ./stitched_metadata --out-dir ./napari_metadata \\
        --annotations pipeline/reference/gbm_fov_annotations.csv \\
        --crosswalk pipeline/reference/gbm_slide_name_crosswalk.csv
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
# Napari's cell_ID is c_<slide number>_<FOV>_<cell within FOV>.
NAPARI_ID_RE = re.compile(r"^c_(?P<slide>\d+)_(?P<fov>\d+)_(?P<cell>\d+)$")
CROSSWALK_FOLDER = "AtoMx_flatfile_folder"
CROSSWALK_CANONICAL = "Canonical_slide_name"
ALLOWED_COLUMNS = ("Region", "Case", "Block")
METADATA_SUFFIX = "_metadata.csv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata-dir", type=Path, required=True,
                   help="Directory of <slide>_metadata.csv, named by AtoMx flat-file folder.")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--annotations", type=Path, required=True,
                   help="AtoMx annotation reference: Slide, Case, Block, Region, FOVs.")
    p.add_argument("--crosswalk", type=Path, required=True,
                   help="Crosswalk from canonical slide name to flat-file folder name.")
    p.add_argument("--columns", default="Region",
                   help=f"Comma-separated subset of {', '.join(ALLOWED_COLUMNS)} "
                        f"(default Region).")
    return p.parse_args()


def deterministic_color(value: str) -> str:
    """Deterministic hex colour -- identical to generate-slide-metadata.py.

    Sharing the implementation matters: a label must keep the same colour whichever script
    wrote the column, or the same category changes colour between layers.
    """
    return duckdb.sql("SELECT printf('#%06X', abs(hash($1)) % 16777216)",
                      params=[value]).fetchone()[0]


def load_fov_table(annotations: Path, crosswalk: Path, columns: list[str]) -> pd.DataFrame:
    """Per-FOV annotation indexed by (flat-file folder, FOV)."""
    # The upstream AtoMx file carries a UTF-8 BOM that would corrupt the first header.
    ann = pd.read_csv(annotations, encoding="utf-8-sig")
    need = {"Slide", "FOVs", *columns}
    if not need.issubset(ann.columns):
        sys.exit(f"ERROR: {annotations} is missing {sorted(need - set(ann.columns))}.")
    xw = pd.read_csv(crosswalk)
    for col in (CROSSWALK_FOLDER, CROSSWALK_CANONICAL):
        if col not in xw.columns:
            sys.exit(f"ERROR: {crosswalk} is missing {col}.")

    ann = ann.merge(xw[[CROSSWALK_CANONICAL, CROSSWALK_FOLDER]],
                    left_on="Slide", right_on=CROSSWALK_CANONICAL, how="left")
    unmapped = ann[CROSSWALK_FOLDER].isna()
    if unmapped.any():
        print(f"WARNING: {ann.loc[unmapped, 'Slide'].nunique()} annotated slide(s) are not in "
              f"the crosswalk and are ignored.", file=sys.stderr)
        ann = ann[~unmapped]
    ann["FOVs"] = ann["FOVs"].astype(int)
    if ann.duplicated([CROSSWALK_FOLDER, "FOVs"]).any():
        sys.exit("ERROR: the annotation reference has more than one row per (slide, FOV); a FOV "
                 "must belong to exactly one case.")
    return ann.set_index([CROSSWALK_FOLDER, "FOVs"])[columns]


def annotate_slide(meta: pd.DataFrame, slide: str, fovs: pd.DataFrame,
                   columns: list[str]) -> tuple[pd.DataFrame, int, int]:
    """Add the columns (and their _color pairs) to one slide. Returns (frame, n_unmatched_ids,
    n_unannotated_cells)."""
    ids = meta[CELL_ID_COLUMN].astype(str).str.extract(NAPARI_ID_RE)
    bad_id = ids["fov"].isna()
    fov = ids["fov"].astype("Int64")

    # A slide absent from the reference gets an empty lookup, so every cell comes out blank and
    # the per-slide report flags it, rather than raising partway through a directory of slides.
    known = slide in fovs.index.get_level_values(0)
    slide_fovs = fovs.loc[slide] if known else None

    out = meta.copy()
    n_unannotated = 0
    for col in columns:
        lookup = slide_fovs[col] if known else pd.Series(dtype=object)
        values = fov.map(lookup).astype(object)
        values = values.where(values.notna(), "")
        values[bad_id.to_numpy()] = ""
        out[col] = values.astype(str).to_numpy()
        colors = {v: deterministic_color(v) for v in sorted(set(out[col])) if v}
        out[col + COLOR_SUFFIX] = out[col].map(colors).fillna("").to_numpy()
        n_unannotated = max(n_unannotated, int((out[col] == "").sum()))
    return out, int(bad_id.sum()), n_unannotated


def main() -> None:
    args = parse_args()
    columns = [c.strip() for c in args.columns.split(",") if c.strip()]
    unknown = [c for c in columns if c not in ALLOWED_COLUMNS]
    if unknown or not columns:
        sys.exit(f"ERROR: --columns must be a subset of {list(ALLOWED_COLUMNS)}; got {columns}")

    fovs = load_fov_table(args.annotations, args.crosswalk, columns)
    files = sorted(args.metadata_dir.glob(f"*{METADATA_SUFFIX}"))
    if not files:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} files in {args.metadata_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"{len(files)} slide file(s); adding {', '.join(columns)}\n")
    print(f"  {'slide':22s} {'cells':>9s} {'unannotated':>12s}")
    total_bad = 0
    for path in files:
        slide = path.name[: -len(METADATA_SUFFIX)]
        meta = pd.read_csv(path, dtype=str)
        if CELL_ID_COLUMN not in meta.columns:
            sys.exit(f"ERROR: {path} has no '{CELL_ID_COLUMN}' column; columns are "
                     f"{list(meta.columns)}")
        clash = [c for c in columns if c in meta.columns]
        if clash:
            sys.exit(f"ERROR: {path.name} already has {clash}. Refusing to overwrite.")
        out, bad_ids, n_un = annotate_slide(meta, slide, fovs, columns)
        if len(out) != len(meta):
            sys.exit(f"ERROR: row count changed for {slide} ({len(meta):,} -> {len(out):,}).")
        out.to_csv(args.out_dir / path.name, index=False)
        total_bad += bad_ids
        flag = "   <-- CHECK: none matched" if n_un == len(out) else ""
        print(f"  {slide:22s} {len(out):>9,d} {n_un:>12,d}{flag}")
    if total_bad:
        print(f"\nWARNING: {total_bad:,} cell ids did not look like c_<slide>_<fov>_<cell> and "
              f"were left blank.", file=sys.stderr)
    print(f"\nWrote {len(files)} file(s) to {args.out_dir}")


if __name__ == "__main__":
    main()
