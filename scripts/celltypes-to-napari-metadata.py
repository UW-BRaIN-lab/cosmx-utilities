#!/usr/bin/env python3
"""Write Napari `_metadata.csv` files coloured by our own InSituType cell typing.

`generate-slide-metadata.py` builds Napari metadata from the AtoMx flat files, so it can
only show AtoMx's own cell-typing column. This writes the same shape of file from a
pipeline stage-4 result instead, so a de novo cluster can be viewed in situ on the slide --
which is what lets a neuropathologist say what `o` or `m` actually is.

THE CELL-ID JOIN IS THE WHOLE PROBLEM. Stage 1 replaces the per-FOV `cell_ID` with a
globally unique index (`<slide_id>_F<fov>_C<cell_ID>`) so per-slide AnnDatas can be
concatenated, but Napari keys on the flat file's own `cell_id` (e.g. `c_1_2_345`). Neither
is derivable from the other, so this parses `fov` and `cell_ID` back out of the pipeline
index and joins them against the slide's `_metadata_file.csv.gz` on that pair -- the same
(fov, cell_ID) key `flatfiles_to_anndata.py` used to build the index in the first place.
Any cell that fails to join is reported, never silently dropped.

Reads:
  --typing-h5   stage-4 insitutype_result.h5 (/cell_id, /cell_type, /prob)
  --flatfiles   directory of <slide>_metadata_file.csv.gz, one per slide

Writes one `<slide>_metadata.csv` per slide with `cell_id`, the cell-type column, and its
posterior, ready to drop next to the Napari slide directory.

Usage:
    uv run python scripts/celltypes-to-napari-metadata.py \\
        --typing-h5 retina_k19.h5 --flatfiles ./flat --out-dir ./napari_metadata
    # only the de novo clusters, everything else blanked out, to make them pop:
    uv run python scripts/celltypes-to-napari-metadata.py ... --denovo-only
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Stage 1 builds the obs index as "<slide_id>_F<fov>_C<cell_ID>"; both fields are needed to
# get back to the flat file's own cell_id. slide_id itself may contain underscores.
INDEX_RE = re.compile(r"^(?P<slide>.+)_F(?P<fov>\d+)_C(?P<cell>.+)$")
# InSituType names de novo clusters with bare letters; everything else is a reference type.
DENOVO_RE = re.compile(r"^[a-z]{1,2}$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="Stage-4 insitutype_result.h5.")
    p.add_argument("--flatfiles", type=Path, required=True,
                   help="Directory holding <slide>_metadata_file.csv.gz.")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Directory to write <slide>_metadata.csv into.")
    p.add_argument("--column", default="cell_type",
                   help="Name of the cell-type column in the output (default: cell_type).")
    p.add_argument("--denovo-only", action="store_true",
                   help="Blank every named reference type, keeping only de novo clusters, "
                        "so they stand out against an uncoloured background.")
    p.add_argument("--min-prob", type=float, default=0.0,
                   help="Blank calls below this posterior (default 0 = keep all).")
    return p.parse_args()


def _decode(arr) -> np.ndarray:
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def load_typing(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        df = pd.DataFrame({"cell": _decode(f["cell_id"]),
                           "cell_type": _decode(f["cell_type"]),
                           "prob": np.asarray(f["prob"][()], dtype=float)})
    parts = df["cell"].str.extract(INDEX_RE)
    unparsed = int(parts["slide"].isna().sum())
    if unparsed:
        print(f"ERROR: {unparsed:,} cell ids do not match "
              f"'<slide>_F<fov>_C<cell_ID>', e.g. "
              f"{df.loc[parts['slide'].isna(), 'cell'].head(3).tolist()}", file=sys.stderr)
        sys.exit(1)
    df["slide_id"] = parts["slide"]
    df["fov"] = parts["fov"].astype(int)
    df["cell_ID"] = parts["cell"].astype(str)
    return df


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    typed = load_typing(args.typing_h5)
    print(f"{len(typed):,} typed cells across {typed.slide_id.nunique()} slides, "
          f"{typed.cell_type.nunique()} types")

    if args.min_prob > 0:
        low = typed["prob"] < args.min_prob
        print(f"blanking {int(low.sum()):,} calls below posterior {args.min_prob}")
        typed.loc[low, "cell_type"] = ""
    if args.denovo_only:
        named = ~typed["cell_type"].str.fullmatch(DENOVO_RE, na=False)
        print(f"--denovo-only: blanking {int(named.sum()):,} named calls, keeping "
              f"{int((~named).sum()):,} de novo")
        typed.loc[named, "cell_type"] = ""

    total_written = 0
    for slide_id, grp in typed.groupby("slide_id", observed=True):
        meta_path = args.flatfiles / f"{slide_id}_metadata_file.csv.gz"
        if not meta_path.exists():
            print(f"WARN: no flat-file metadata for {slide_id} at {meta_path}; skipping",
                  file=sys.stderr)
            continue
        # Only the join keys and the Napari key are needed; the rest of the file is wide.
        flat = pd.read_csv(meta_path, usecols=["fov", "cell_ID", "cell_id"],
                           dtype={"cell_ID": str, "cell_id": str})
        flat["fov"] = flat["fov"].astype(int)

        merged = flat.merge(grp[["fov", "cell_ID", "cell_type", "prob"]],
                            on=["fov", "cell_ID"], how="left")
        n_typed = int(merged["cell_type"].notna().sum())
        # Cells present on the slide but absent from the typing were dropped by stage-3a QC.
        # They stay in the file with an empty call so Napari still draws them uncoloured.
        n_missing = len(merged) - n_typed
        out = pd.DataFrame({
            "cell_id": merged["cell_id"],
            args.column: merged["cell_type"].fillna(""),
            f"{args.column}_prob": merged["prob"].round(4).fillna(""),
        })
        out_path = args.out_dir / f"{slide_id}_metadata.csv"
        out.to_csv(out_path, index=False)
        total_written += len(out)
        print(f"  {slide_id:14s} {len(out):7,d} cells  {n_typed:7,d} typed  "
              f"{n_missing:6,d} not in typing (QC-dropped)  -> {out_path.name}")

    print(f"\nWrote {total_written:,} rows to {args.out_dir}")


if __name__ == "__main__":
    main()
