#!/usr/bin/env python3
"""Write Napari `_metadata.csv` files coloured by our own InSituType cell typing.

`generate-slide-metadata.py` builds Napari metadata from the AtoMx flat files, so it can
only show AtoMx's own cell-typing column. This writes the same shape of file from a
pipeline stage-4 result instead, so a de novo cluster can be viewed in situ on the slide --
which is what lets a neuropathologist say what `o` or `m` actually is.

THE CELL-ID JOIN IS THE WHOLE PROBLEM. Stage 1 replaces the per-FOV `cell_ID` with a
globally unique index (`<slide_id>_F<fov>_C<cell_ID>`) so per-slide AnnDatas can be
concatenated, while Napari keys on the flat file's own `cell_id` (e.g. `c_2_1_96`). Neither
is derivable from the other. Two ways across:

  --combined-h5ad  PREFERRED. Stage 1 now keeps the flat file's key as
                   obs['flatfile_cell_id'], so the join is a straight lookup on the
                   pipeline index -- nothing is reconstructed and the flat files are not
                   needed at all.
  --flatfiles      FALLBACK, for results built before stage 1 kept that column. Parses fov
                   and cell_ID back out of the pipeline index and re-joins the slide's
                   `_metadata_file.csv.gz` on that pair -- the same (fov, cell_ID) key
                   flatfiles_to_anndata.py used to build the index in the first place.

Either way, cells on the slide with no typing call (stage-3a QC dropped them) are written
with an empty value so Napari still draws them uncoloured; the count is always reported.

Reads:
  --typing-h5   stage-4 insitutype_result.h5 (/cell_id, /cell_type, /prob)

Writes one `<slide>_metadata.csv` per slide with `cell_id`, the cell-type column, and its
posterior, ready to drop next to the Napari slide directory. With `--merge-into` it instead
adds those columns to an existing Napari metadata file, so one slide can carry AtoMx's
cell-typing column and ours side by side and be toggled between them.

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

import duckdb
import h5py
import numpy as np
import pandas as pd

# Stage 1 builds the obs index as "<slide_id>_F<fov>_C<cell_ID>"; both fields are needed to
# get back to the flat file's own cell_id. slide_id itself may contain underscores.
INDEX_RE = re.compile(r"^(?P<slide>.+)_F(?P<fov>\d+)_C(?P<cell>.+)$")
# InSituType names de novo clusters with bare letters; everything else is a reference type.
DENOVO_RE = re.compile(r"^[a-z]{1,2}$")
# napari-cosmx conventions, matching scripts/generate-slide-metadata.py: the key column is
# spelled `cell_ID` even though it holds the flat file's `cell_id` VALUES, and every value
# column must be paired with a `<name>_color` or the layer renders uncoloured.
CELL_ID_COLUMN = "cell_ID"
COLOR_SUFFIX = "_color"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="Stage-4 insitutype_result.h5.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--combined-h5ad", type=Path,
                     help="Stage-3a combined_qc.h5ad carrying obs['flatfile_cell_id'] — "
                          "the direct route, no re-joining. Preferred.")
    src.add_argument("--flatfiles", type=Path,
                     help="Directory of <slide>_metadata_file.csv.gz. Fallback for data "
                          "built before stage 1 kept flatfile_cell_id: rebuilds the key by "
                          "re-joining on (fov, cell_ID).")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Directory to write <slide>_metadata.csv into.")
    p.add_argument("--column", default="cell_type",
                   help="Name of the cell-type column in the output (default: cell_type).")
    p.add_argument("--denovo-only", action="store_true",
                   help="Blank every named reference type, keeping only de novo clusters, "
                        "so they stand out against an uncoloured background.")
    p.add_argument("--min-prob", type=float, default=0.0,
                   help="Blank calls below this posterior (default 0 = keep all).")
    p.add_argument("--with-prob", action="store_true",
                   help="Also write a <column>_prob column. Off by default: napari-cosmx "
                        "expects every value column to have a paired _color, and a bare "
                        "numeric column has none.")
    p.add_argument("--split-spec", type=Path, default=None,
                   help="CSV (label, column, values, new_name) breaking part of a call out "
                        "as its own Napari value — e.g. the cells of one unsupervised "
                        "cluster inside a de novo letter. Needs --combined-h5ad pointed at "
                        "an AnnData carrying that column (cosmx_clustered.h5ad for leiden).")
    p.add_argument("--split-obs-h5ad", type=Path, default=None,
                   help="A second AnnData supplying the obs columns --split-spec keys on, "
                        "joined by cell index. Needed because the Napari join key "
                        "(flatfile_cell_id) lives in the stage-3a combined AnnData while "
                        "obs['leiden'] lives in the stage-3c clustered one.")
    p.add_argument("--merge-into", type=Path,
                   help="Directory of existing <slide>_metadata.csv to add our column(s) "
                        "to, rather than writing a fresh file. Keeps every existing column "
                        "-- so a slide can carry AtoMx's call and ours side by side and be "
                        "toggled between them in Napari.")
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


def deterministic_color(value: str) -> str:
    """Deterministic hex colour for a label — identical to generate-slide-metadata.py.

    Sharing the implementation matters: a label must keep the same colour whichever script
    wrote the column, or the same cell type changes colour between layers.
    """
    return duckdb.sql("SELECT printf('#%06X', abs(hash($1)) % 16777216)",
                      params=[value]).fetchall()[0][0]


def napari_frame(args, cell_ids: pd.Series, calls: pd.Series,
                 probs: pd.Series) -> pd.DataFrame:
    """Assemble one slide's columns in napari-cosmx's expected shape."""
    calls = calls.fillna("").astype(str)
    colors = {v: deterministic_color(v) for v in sorted(set(calls)) if v}
    frame = pd.DataFrame({
        CELL_ID_COLUMN: cell_ids.astype(str).to_numpy(),
        args.column: calls.to_numpy(),
        args.column + COLOR_SUFFIX: calls.map(colors).fillna("").to_numpy(),
    })
    if args.with_prob:
        frame[f"{args.column}_prob"] = probs.round(4).fillna("").to_numpy()
    return frame


def write_slide(args, slide_id: str, out: pd.DataFrame) -> int:
    """Write one slide's metadata, optionally merged into an existing Napari file.

    Merging is a LEFT join from the existing file, so its rows, their order, and every
    column it already had are preserved exactly -- Napari keeps working with the old
    colouring and simply gains another column to switch to. A name collision is a hard
    error rather than an overwrite: silently replacing AtoMx's own cell-typing column is
    the one outcome that would destroy the comparison this flag exists to enable.
    """
    out_path = args.out_dir / f"{slide_id}_metadata.csv"
    if not args.merge_into:
        out.to_csv(out_path, index=False)
        return len(out)

    existing_path = args.merge_into / f"{slide_id}_metadata.csv"
    if not existing_path.exists():
        print(f"WARN: nothing to merge into at {existing_path}; writing our columns alone",
              file=sys.stderr)
        out.to_csv(out_path, index=False)
        return len(out)

    existing = pd.read_csv(existing_path, dtype=str)
    if CELL_ID_COLUMN not in existing.columns:
        print(f"ERROR: {existing_path} has no '{CELL_ID_COLUMN}' column to join on; "
              f"columns are {list(existing.columns)}", file=sys.stderr)
        sys.exit(1)
    clash = [c for c in out.columns if c != CELL_ID_COLUMN and c in existing.columns]
    if clash:
        print(f"ERROR: {existing_path} already has column(s) {clash}. Refusing to "
              f"overwrite — pass a different --column.", file=sys.stderr)
        sys.exit(1)

    merged = existing.merge(out, on=CELL_ID_COLUMN, how="left")
    if len(merged) != len(existing):
        # Duplicate cell_ids on either side would silently multiply rows.
        print(f"ERROR: merge changed the row count for {slide_id} "
              f"({len(existing):,} -> {len(merged):,}); duplicate cell_id?", file=sys.stderr)
        sys.exit(1)
    # Count non-EMPTY, not non-null: untyped cells are written as "" so Napari still draws
    # them, so notna() would report every row as matched and hide a broken join.
    filled = int((merged[args.column].fillna("").astype(str) != "").sum())
    print(f"    merged into {existing_path.name}: {filled:,} of {len(existing):,} "
          f"existing rows got a call, {len(existing.columns)} existing column(s) kept")
    merged.fillna("").to_csv(out_path, index=False)
    return len(merged)


SPLIT_COLUMNS = ("label", "column", "values", "new_name")
SPLIT_VALUE_SEP = ";"


def load_split_spec(spec_path: Path) -> pd.DataFrame:
    """Read and validate the split spec. Mirrors merge_reference_types.py's spec style."""
    spec = pd.read_csv(spec_path, dtype=str).fillna("")
    for col in SPLIT_COLUMNS:
        if col not in spec.columns:
            sys.exit(f"ERROR: split spec {spec_path} missing '{col}' column")
    for col in SPLIT_COLUMNS:
        spec[col] = spec[col].str.strip()
    blank = spec[(spec == "").any(axis=1)]
    if not blank.empty:
        sys.exit(f"ERROR: split spec {spec_path} has blank fields in row(s) "
                 f"{blank.index.tolist()}")
    clash = sorted(set(spec["new_name"]) & set(spec["label"]))
    if clash:
        sys.exit(f"ERROR: split new_name(s) collide with a label being split: {clash}")
    return spec


def apply_splits(calls: pd.Series, obs: pd.DataFrame, spec: pd.DataFrame) -> pd.Series:
    """Rename the subset of `label` whose `column` is in `values` to `new_name`.

    Every clause is validated and must move at least one cell. A split that silently moves
    nothing would leave the population it was meant to surface invisible in Napari, which
    is the one failure the caller cannot see by looking at the output.
    """
    out = calls.copy()
    for _, row in spec.iterrows():
        label, column, new_name = row["label"], row["column"], row["new_name"]
        values = {v.strip() for v in row["values"].split(SPLIT_VALUE_SEP) if v.strip()}
        if column not in obs.columns:
            sys.exit(f"ERROR: split spec names obs column '{column}', which this AnnData "
                     f"does not have. Point --combined-h5ad at one that does "
                     f"(cosmx_clustered.h5ad carries 'leiden'). Present: "
                     f"{sorted(obs.columns)[:20]}")
        present = set(obs[column].astype(str).unique())
        unknown = sorted(values - present)
        if unknown:
            sys.exit(f"ERROR: split spec value(s) {unknown} absent from obs['{column}']")
        if not (out == label).any():
            sys.exit(f"ERROR: split spec label '{label}' matches no call")

        moving = (out == label) & obs[column].astype(str).isin(values)
        if not moving.any():
            sys.exit(f"ERROR: split '{label}' x {column} in {sorted(values)} -> "
                     f"'{new_name}' moves no cell; the population would be invisible")
        out[moving] = new_name
        print(f"  split {int(moving.sum()):,} of {int((out == label).sum()) + int(moving.sum()):,} "
              f"'{label}' cells into '{new_name}' ({column} in {sorted(values)})")
    return out


def load_split_columns(h5ad_path: Path, columns: list[str]) -> pd.DataFrame:
    """The obs columns a split keys on, indexed by the pipeline cell index.

    Read separately from the Napari key so the two can come from different files: stage 3c
    carries obs['leiden'] but not obs['flatfile_cell_id'], and stage 3a the reverse.
    """
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        absent = [c for c in columns if c not in adata.obs.columns]
        if absent:
            print(f"ERROR: {h5ad_path} obs has no {absent}. Point --split-obs-h5ad at an "
                  f"AnnData that has them (cosmx_clustered.h5ad carries 'leiden'); present: "
                  f"{sorted(adata.obs.columns)[:20]}", file=sys.stderr)
            sys.exit(1)
        return adata.obs[columns].copy()
    finally:
        adata.file.close()


def write_from_obs(args, typed: pd.DataFrame) -> None:
    """Direct path: obs already carries the Napari key, so just look it up.

    Requires stage 1 to have kept `flatfile_cell_id`. Errors rather than falling back
    silently, since a silent fallback would hide that the AnnData predates that change.
    """
    import anndata as ad

    print(f"Reading {args.combined_h5ad} (backed)")
    adata = ad.read_h5ad(args.combined_h5ad, backed="r")
    missing = [c for c in ("flatfile_cell_id", "slide_id") if c not in adata.obs]
    if missing:
        print(f"ERROR: obs is missing {missing}. This AnnData predates stage 1 keeping the "
              f"flat-file key — re-run stage 1, or use --flatfiles instead.", file=sys.stderr)
        sys.exit(1)

    obs = adata.obs[["flatfile_cell_id", "slide_id"]].copy()
    joined = obs.join(typed.set_index("cell")[["cell_type", "prob"]], how="left")
    n_typed = int(joined["cell_type"].notna().sum())
    print(f"{n_typed:,} of {len(joined):,} cells have a typing call "
          f"({len(joined) - n_typed:,} without)")

    total = 0
    for slide_id, grp in joined.groupby("slide_id", observed=True):
        out = napari_frame(args, grp["flatfile_cell_id"], grp["cell_type"], grp["prob"])
        print(f"  {str(slide_id):14s} {len(out):7,d} cells  "
              f"{int(grp['cell_type'].notna().sum()):7,d} typed")
        total += write_slide(args, str(slide_id), out)
    print(f"\nWrote {total:,} rows to {args.out_dir}")


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

    if args.split_spec:
        spec = load_split_spec(args.split_spec)
        source = args.split_obs_h5ad or args.combined_h5ad
        if source is None:
            print("ERROR: --split-spec needs --split-obs-h5ad (or --combined-h5ad) to read "
                  "the obs column it keys on.", file=sys.stderr)
            sys.exit(1)
        print(f"Applying {len(spec)} split(s) from {args.split_spec}, "
              f"keyed on obs from {source}")
        extra = load_split_columns(source, sorted(spec["column"].unique()))
        aligned = extra.reindex(typed["cell"]).reset_index(drop=True)
        unmatched = int(aligned.isna().all(axis=1).sum())
        if unmatched:
            print(f"ERROR: {unmatched:,} of {len(typed):,} typed cells are absent from "
                  f"{source}; the split would apply to some cells and silently not others.",
                  file=sys.stderr)
            sys.exit(1)
        typed["cell_type"] = apply_splits(typed["cell_type"], aligned, spec)

    if args.combined_h5ad:
        write_from_obs(args, typed)
        return

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
        out = napari_frame(args, merged["cell_id"], merged["cell_type"], merged["prob"])
        print(f"  {slide_id:14s} {len(out):7,d} cells  {n_typed:7,d} typed  "
              f"{n_missing:6,d} not in typing (QC-dropped)")
        total_written += write_slide(args, str(slide_id), out)

    print(f"\nWrote {total_written:,} rows to {args.out_dir}")


if __name__ == "__main__":
    main()
