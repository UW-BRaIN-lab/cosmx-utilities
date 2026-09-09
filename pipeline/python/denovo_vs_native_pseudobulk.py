#!/usr/bin/env python3
"""Marker heatmap inputs: a de-novo letter's cells, split by the GBmap class they are FORCED
into, against the cells that got that same GBmap class NATIVELY.

The PI's question, in her words: for the cells that were "different enough" to be put into a
de-novo cluster when the semi-supervised fit allowed one, how do they compare to the cells that
received a GBmap call without needing a de-novo bin? Take de-novo `t`. The forced supervised run
scatters it across OPC-like (32%), AC-like (30%) and MES-like_hypoxia_MHC (19%). So build, for
each destination D:

    t->D        cells with semisup == t   AND forced == D
    D [native]  cells with semisup == D   (the fit named them; no de-novo bin needed)

and put the pair side by side in one heatmap. Two things are then readable at once:

  1. DOWN each pair — how t->D differs from native D. Note this difference is guaranteed by
     construction: a cell landed in `t` precisely because it fit t's profile better than any
     GBmap profile. So "they differ" is definitional, NOT a finding. What is informative is
     WHICH genes differ, and whether they form a coherent programme (hypoxia, proliferation,
     stress) or look like scatter.
  2. ACROSS the t->D columns — whether the letter's cells that got sent to different GBmap
     classes actually differ from EACH OTHER. If t->OPC-like and t->AC-like are
     transcriptionally the same, the supervised split of `t` is arbitrary and `t` is one
     population. If they differ, the forced call is recovering real substructure inside `t`.
     This is usually the more interesting axis, and it is why every destination is drawn in one
     figure rather than one figure per destination.

A `<letter> [all]` column is included by default so each subset can also be read against the
letter as a whole.

Inputs:
  --counts-h5     anchor_input.h5 (stage 4a): /counts CSC genes x cells, /genes, /cell_id.
  --typing-h5     anchor_typing.h5 — semi-supervised label per cell (/cell_id, /cell_type).
  --forced-csv    75c's forced_named_posteriors.csv (cell_id, top1_type) — the forced GBmap
                  class, read off the fit's own stored logliks. Do NOT use 75's re-scored
                  supervised_gbmap_posteriors.csv here; it disagrees with the fit on 45% of
                  cells (see 75c).
  --letter        the de-novo letter to dissect, e.g. t
  --destinations  comma-separated GBmap classes to compare against, e.g.
                  OPC-like,AC-like,MES-like_hypoxia_MHC

Writes exactly what R/marker_heatmap.R reads in its no-region mode:
  <out>/marker_heatmap_zmatrix.csv    genes x groups, z-scored
  <out>/top_markers_per_cluster.csv   gene, cluster (row split)
  <out>/group_sizes.csv               group, n_cells, dropped

Usage:
    uv run python pipeline/python/denovo_vs_native_pseudobulk.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv \\
        --letter t --destinations OPC-like,AC-like,MES-like_hypoxia_MHC \\
        --output-dir t_vs_native
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from anchor_profiles import decode, read_cell_calls
from pseudobulk_core import (DEFAULT_SCALE_FACTOR, group_means, log_normalize, onehot,
                             select_markers, zscore_rows)

NATIVE_SUFFIX = " [native]"
ALL_SUFFIX = " [all]"
FORCED_SEP = "->"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True,
                   help="anchor_input.h5: /counts (CSC genes x cells), /genes, /cell_id.")
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="anchor_typing.h5: semi-supervised /cell_type per /cell_id.")
    p.add_argument("--forced-csv", type=Path, required=True,
                   help="forced_named_posteriors.csv from 75c (cell_id, top1_type).")
    p.add_argument("--letter", required=True, help="De-novo letter to dissect, e.g. t")
    p.add_argument("--destinations", required=True,
                   help="Comma-separated GBmap classes to compare against.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=8,
                   help="Top markers selected per group (default 8).")
    p.add_argument("--min-group-n", type=int, default=50,
                   help="Drop groups with fewer cells than this (default 50).")
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    p.add_argument("--no-all-column", action="store_true",
                   help="Omit the '<letter> [all]' whole-letter reference column.")
    return p.parse_args()


def read_counts(path: Path) -> tuple[sp.csc_matrix, np.ndarray, np.ndarray]:
    """Read the stage-4a CSC genes x cells counts, returning (matrix, genes, cell_ids)."""
    with h5py.File(path, "r") as f:
        for key in ("counts/shape", "counts/data", "counts/indices", "counts/indptr",
                    "genes", "cell_id"):
            if key not in f:
                sys.exit(f"ERROR: /{key} missing in {path}.")
        n_genes, n_cells = (int(x) for x in f["counts/shape"][()])
        mat = sp.csc_matrix(
            (np.asarray(f["counts/data"][()], dtype=np.float64),
             np.asarray(f["counts/indices"][()], dtype=np.int64),
             np.asarray(f["counts/indptr"][()], dtype=np.int64)),
            shape=(n_genes, n_cells))
        genes = np.asarray(decode(f["genes"][()]))
        cell_id = np.asarray(decode(f["cell_id"][()]))
    return mat, genes, cell_id


def build_groups(semisup: pd.Series, forced: pd.Series, letter: str,
                 destinations: list[str], include_all: bool) -> pd.Series:
    """Per-cell group label, or NaN for cells in none of the compared groups.

    A cell belongs to at most one group: `<letter>-><D>` when the fit called it `letter` and the
    forced run sends it to D, or `<D> [native]` when the fit named it D outright. Cells of the
    letter that went to some other destination are unlabelled (they are not part of this
    comparison) unless the whole-letter reference column is on, which additionally tags every
    `letter` cell.
    """
    group = pd.Series(pd.NA, index=semisup.index, dtype="object")
    is_letter = semisup == letter
    for dest in destinations:
        group[is_letter & (forced == dest)] = f"{letter}{FORCED_SEP}{dest}"
        group[semisup == dest] = f"{dest}{NATIVE_SUFFIX}"
    # The whole-letter column reuses the letter's cells, so it is built as a separate pass
    # (a cell can contribute to both its <letter>-><D> column and the <letter> [all] column).
    return group, is_letter if include_all else None


def main() -> None:
    args = parse_args()
    destinations = [d.strip() for d in args.destinations.split(",") if d.strip()]
    if not destinations:
        sys.exit("ERROR: --destinations is empty.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]]
    forced = pd.read_csv(args.forced_csv, usecols=["cell_id", "top1_type"])
    labels = calls.merge(forced, on="cell_id", how="inner").set_index("cell_id")
    print(f"Labels joined for {len(labels):,} cells")

    missing = [d for d in destinations if d not in set(labels["top1_type"])]
    if missing:
        sys.exit(f"ERROR: destination(s) never appear as a forced call: {missing}")
    if args.letter not in set(labels["cell_type"]):
        sys.exit(f"ERROR: --letter {args.letter!r} is not a semi-supervised label.")

    group, letter_mask = build_groups(labels["cell_type"], labels["top1_type"],
                                      args.letter, destinations, not args.no_all_column)

    counts, genes, cell_id = read_counts(args.counts_h5)
    print(f"Counts: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")

    # Align the labels to the counts' cell order, then keep only the cells we compare.
    group = group.reindex(cell_id)
    letter_mask = letter_mask.reindex(cell_id).fillna(False).to_numpy() \
        if letter_mask is not None else None
    keep = group.notna().to_numpy()
    if letter_mask is not None:
        keep = keep | letter_mask
    if not keep.any():
        sys.exit("ERROR: no cells matched the requested letter/destinations.")
    print(f"Comparing {int(keep.sum()):,} cells")

    # Subset columns on the CSC (cells are columns -> cheap) BEFORE transposing to cells x genes.
    norm = log_normalize(counts[:, keep].T.tocsr(), args.scale_factor)
    kept_group = group[keep]
    kept_letter = letter_mask[keep] if letter_mask is not None else None

    # Mean log-norm profile per group (genes x groups) over the paired groups only; the
    # whole-letter column is appended after, since its cells overlap the <letter>-><D> ones.
    paired = kept_group.notna().to_numpy()
    oh, glabels = onehot(kept_group[paired].to_numpy())
    profile = pd.DataFrame(group_means(norm[paired], oh).T, index=genes, columns=glabels)
    sizes = pd.Series(np.asarray(oh.sum(axis=0)).ravel(), index=glabels, dtype=int)

    all_col = f"{args.letter}{ALL_SUFFIX}"
    if kept_letter is not None and kept_letter.any():
        profile[all_col] = np.asarray(norm[kept_letter].mean(axis=0)).ravel()
        sizes[all_col] = int(kept_letter.sum())

    # Interleave each forced subset with its native counterpart so the pair reads together,
    # then the whole-letter reference last.
    order = []
    for dest in destinations:
        for col in (f"{args.letter}{FORCED_SEP}{dest}", f"{dest}{NATIVE_SUFFIX}"):
            if col in profile.columns:
                order.append(col)
    if all_col in profile.columns:
        order.append(all_col)
    profile = profile[order]

    small = sizes[sizes < args.min_group_n].index.tolist()
    if small:
        print(f"Dropping {len(small)} group(s) under {args.min_group_n} cells: {small}")
        profile = profile.drop(columns=[c for c in small if c in profile.columns])
    if profile.shape[1] < 2:
        sys.exit("ERROR: fewer than 2 groups survive; nothing to compare.")

    print("Group sizes:")
    for col in profile.columns:
        print(f"  {col:44s} {sizes[col]:>9,}")

    markers, gene_to_group = select_markers(profile, list(profile.columns), args.top_n)
    pb_z = zscore_rows(profile.loc[markers])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pb_z.to_csv(args.output_dir / "marker_heatmap_zmatrix.csv")
    pd.DataFrame({"cluster": [gene_to_group[g] for g in markers],
                  "gene": markers}).to_csv(
        args.output_dir / "top_markers_per_cluster.csv", index=False)
    sizes.rename("n_cells").rename_axis("group").reset_index().assign(
        dropped=lambda d: d["group"].isin(small)).to_csv(
        args.output_dir / "group_sizes.csv", index=False)

    print(f"\nWrote {args.output_dir} — {pb_z.shape[0]} markers x {pb_z.shape[1]} groups.")
    print("Render on the Mac (the container has no ComplexHeatmap):")
    print(f"  Rscript pipeline/R/marker_heatmap.R {args.output_dir} <out_dir> \"\" "
          f"\"{args.letter} forced into GBmap classes vs cells natively called those classes\"")


if __name__ == "__main__":
    main()
