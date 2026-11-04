#!/usr/bin/env python3
"""Decompose a de novo InSituType cluster against an unsupervised partition, with QC depth.

A de novo letter that spans several unsupervised clusters is either (a) one real population
the typing over-merged with others, or (b) a pool of low-signal cells that share no biology
beyond being faint. Those two call for opposite actions -- (a) wants the hierarchy to pull
the components apart, (b) wants the cells excluded -- and a composition table alone cannot
tell them apart, because both look like "this letter spans four Leiden clusters".

What separates them is depth. A genuinely novel population sits at cohort-typical counts and
simply has no reference profile; a low-signal pool sits well below cohort median on counts
and genes detected, and often above it on negative-probe proportion. So every component row
here carries its QC distribution beside its composition, and a `__cohort__` baseline row is
emitted first so the comparison is on the page rather than in the reader's head.

Scope note: the input AnnData is the QC-PASSING cohort, so cells dropped at stage 1 are not
represented here at all. This answers "are the cells that survived QC faint?", not "how many
failed it" -- that second question is answered from the stitched Napari metadata, which
carries every segmented cell.

Reads the stage-3c clustered AnnData (obs: the partition column, the QC columns written by
concat_qc_anndata.py, Region) BACKED, so only obs is materialized, and the stage-4 InSituType
result HDF5 (/cell_id, /cell_type, /prob) written by insitutype_typing.R.

With --reference-csv each component also gets its nearest reference types by log1p cosine.
That is the column that decides what to DO about an over-merged cluster: components whose
nearest types sit under different hierarchy branches will be separated by the tree on their
own, while only the components nearest the same branch belong in a new leaf built for it.
Building that leaf from the whole cluster instead would bake the other components' signal
into its profile, and a leaf fit to our own data then out-competes the atlas profiles it
should never have been competing with.

Emits, per --cluster:
  <cluster>_components.csv   one row per partition component (plus the cohort baseline):
                             cells, % of the cluster, % of that component, median of each
                             --depth-cols, median negprobe proportion, median area, median
                             posterior, dominant Region and its share, and (with
                             --reference-csv) the nearest reference types and cosines.

Usage:
    uv run python pipeline/python/diagnose_denovo_components.py \\
        --clustered-h5ad cosmx_clustered.h5ad \\
        --typing-h5 insitutype_result.h5 \\
        --cluster o --cluster a \\
        --reference-csv pipeline/reference/retina_combined_panel.csv \\
        --out-dir denovo_components
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# The log-normalize / one-hot / group-mean trio already exists for the marker heatmap and is
# the same arithmetic here; importing keeps one definition of "mean log-norm profile".
_MP_PATH = Path(__file__).resolve().parent / "marker_pseudobulk.py"
_mp_spec = importlib.util.spec_from_file_location("marker_pseudobulk", _MP_PATH)
mp = importlib.util.module_from_spec(_mp_spec)
_mp_spec.loader.exec_module(mp)

DEFAULT_PARTITION_COL = "leiden"
# Same name and default as cluster_batch_entropy.py, so the two read alike.
DEFAULT_DEPTH_COLS = "total_counts,qc_genes_detected"
NEGPROBE_COL = "qc_negprobe_prop"
AREA_COL = "qc_area"
REGION_COL = "Region"
POSTERIOR_COL = "prob"
BASELINE_LABEL = "__cohort__"
# Components below this share of the cluster are rolled up rather than listed: a long tail of
# 50-cell components buries the structure the table exists to show.
DEFAULT_MIN_COMPONENT_PCT = 1.0
TAIL_LABEL = "__other__"
DEFAULT_TOP_REFERENCE = 2
GENE_PROBE_VALUE = "gene"
PROBE_TYPE_COL = "probe_type"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="Stage-3c cosmx_clustered.h5ad (obs: partition col + QC cols).")
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="Stage-4 insitutype_result.h5 (/cell_id, /cell_type, /prob).")
    p.add_argument("--cluster", action="append", required=True, metavar="LABEL",
                   help="De novo cluster to decompose (repeat for several, e.g. -c o -c a).")
    p.add_argument("--partition-col", default=DEFAULT_PARTITION_COL,
                   help=f"obs column holding the unsupervised partition "
                        f"(default: {DEFAULT_PARTITION_COL}).")
    p.add_argument("--depth-cols", default=DEFAULT_DEPTH_COLS,
                   help=f"Comma-separated obs columns to report medians for "
                        f"(default: {DEFAULT_DEPTH_COLS}).")
    p.add_argument("--min-component-pct", type=float, default=DEFAULT_MIN_COMPONENT_PCT,
                   help=f"Roll components below this %% of the cluster into "
                        f"'{TAIL_LABEL}' (default: {DEFAULT_MIN_COMPONENT_PCT}).")
    p.add_argument("--reference-csv", type=Path, default=None,
                   help="InSituType reference profiles (genes x types). Adds each component's "
                        "nearest reference types by log1p cosine -- which is what says WHICH "
                        "branch a component belongs under, not merely that it is distinct.")
    p.add_argument("--top-reference", type=int, default=DEFAULT_TOP_REFERENCE,
                   help=f"Nearest reference types to report per component "
                        f"(default: {DEFAULT_TOP_REFERENCE}).")
    p.add_argument("--scale-factor", type=float, default=mp.DEFAULT_SCALE_FACTOR,
                   help="Log-normalize scale factor (matches the marker heatmap).")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for the CSVs.")
    return p.parse_args()


def _decode(arr) -> np.ndarray:
    """hdf5r writes variable-length UTF-8; h5py hands them back as bytes or str."""
    values = np.asarray(arr[()])
    if values.dtype.kind == "S" or (values.size and isinstance(values.flat[0], bytes)):
        return np.array([v.decode("utf-8") for v in values.ravel()])
    return values.astype(str).ravel()


def load_typing(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        missing = [k for k in ("cell_id", "cell_type") if k not in f]
        if missing:
            sys.exit(f"ERROR: {path} is missing dataset(s) {missing}; expected an "
                     f"insitutype_typing.R result with /cell_id and /cell_type")
        typed = pd.DataFrame(
            {"cell_type": _decode(f["cell_type"])},
            index=pd.Index(_decode(f["cell_id"]), name="cell"),
        )
        # Older results predate the posterior being written; the table degrades to NaN
        # rather than failing, since depth is the load-bearing column here.
        typed[POSTERIOR_COL] = (np.asarray(f[POSTERIOR_COL][()], dtype=np.float64)
                                if POSTERIOR_COL in f else np.nan)
    return typed


def load_obs(path: Path, required: list[str], optional: list[str]) -> pd.DataFrame:
    """Read obs only, from a backed handle closed before returning (X stays on disk)."""
    adata = ad.read_h5ad(path, backed="r")
    try:
        absent = [c for c in required if c not in adata.obs.columns]
        if absent:
            sys.exit(f"ERROR: {path} obs is missing column(s) {absent}; "
                     f"present: {sorted(adata.obs.columns)[:40]}")
        for col in optional:
            if col not in adata.obs.columns:
                print(f"WARN: obs has no '{col}'; that column is omitted", file=sys.stderr)
        present = required + [c for c in optional if c in adata.obs.columns]
        return adata.obs[present].copy()
    finally:
        adata.file.close()


def component_reference_similarity(h5ad: Path, cell_ids: pd.Index,
                                   component_of: pd.Series, reference_csv: Path,
                                   scale_factor: float, top_n: int) -> pd.DataFrame:
    """Nearest reference types per component, by cosine of log-space mean profiles.

    Cosine is taken in log space on both sides -- log1p of the reference, mean log-norm of
    the cells -- because a linear-space cosine between a 6k panel's profiles is dominated by
    a handful of very high-count genes and reads ~0.99 for everything.

    Only the cells of the components in question are materialized; the rest of the matrix
    stays on disk.
    """
    reference = pd.read_csv(reference_csv, index_col=0)
    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        wanted = adata.obs_names.isin(set(cell_ids))
        subset = adata[wanted].to_memory()
    finally:
        adata.file.close()

    if PROBE_TYPE_COL in subset.var.columns:
        gene_mask = (subset.var[PROBE_TYPE_COL] == GENE_PROBE_VALUE).to_numpy()
    else:
        print(f"WARN: var has no '{PROBE_TYPE_COL}'; treating every column as a gene probe",
              file=sys.stderr)
        gene_mask = np.ones(subset.n_vars, dtype=bool)

    counts = sp.csr_matrix(subset.X)
    norm = mp.log_normalize(counts[:, gene_mask], scale_factor)
    gene_names = pd.Index(subset.var_names[gene_mask].to_numpy())

    shared = gene_names.intersection(reference.index)
    if len(shared) < 50:
        sys.exit(f"ERROR: only {len(shared)} genes shared between the AnnData and "
                 f"{reference_csv}; check they describe the same panel")
    print(f"  reference: {reference.shape[1]} types, {len(shared):,} genes shared")

    labels = component_of.reindex(subset.obs_names).to_numpy()
    onehot, order = mp.onehot(labels)
    # group_means returns groups x genes; the cosine below wants genes down the rows.
    profiles = pd.DataFrame(mp.group_means(norm, onehot).T, index=gene_names, columns=order)

    def unit(frame: pd.DataFrame) -> np.ndarray:
        matrix = frame.loc[shared].to_numpy(dtype=np.float64)
        norms = np.linalg.norm(matrix, axis=0)
        norms[norms == 0] = 1.0
        return matrix / norms

    cosines = pd.DataFrame(unit(profiles).T @ unit(np.log1p(reference)),
                           index=order, columns=reference.columns)
    rows = []
    for component in order:
        ranked = cosines.loc[component].sort_values(ascending=False)
        row: dict[str, object] = {"component": component}
        for rank in range(min(top_n, len(ranked))):
            suffix = "" if rank == 0 else f"_{rank + 1}"
            row[f"nearest{suffix}"] = ranked.index[rank]
            row[f"cosine{suffix}"] = round(float(ranked.iloc[rank]), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, depth_cols: list[str], label: str,
              n_parent: int, n_component: int | None) -> dict:
    """One output row: composition on the left, the QC distribution on the right."""
    row: dict[str, object] = {
        "component": label,
        "cells": len(df),
        "pct_of_cluster": round(len(df) / n_parent * 100, 1) if n_parent else np.nan,
        # What share of the partition component this cluster accounts for -- the reverse
        # direction, which says whether the component is mostly this letter or mostly
        # something else. NaN for the cohort baseline, where it is not meaningful.
        "pct_of_component": (round(len(df) / n_component * 100, 1)
                             if n_component else np.nan),
    }
    for col in depth_cols:
        row[f"median_{col}"] = round(float(df[col].median()), 1)
    if NEGPROBE_COL in df:
        row[f"median_{NEGPROBE_COL}"] = round(float(df[NEGPROBE_COL].median()), 4)
    if AREA_COL in df:
        row[f"median_{AREA_COL}"] = round(float(df[AREA_COL].median()), 1)
    if POSTERIOR_COL in df:
        row["median_posterior"] = round(float(df[POSTERIOR_COL].median()), 3)
    if REGION_COL in df and len(df):
        shares = df[REGION_COL].astype(str).value_counts(normalize=True) * 100
        row["dominant_region"] = shares.index[0]
        row["dominant_region_pct"] = round(float(shares.iloc[0]), 1)
    return row


def main() -> None:
    args = parse_args()
    depth_cols = [c.strip() for c in args.depth_cols.split(",") if c.strip()]

    typed = load_typing(args.typing_h5)
    print(f"Typing: {len(typed):,} cells, {typed['cell_type'].nunique()} labels")

    obs = load_obs(args.clustered_h5ad,
                   required=[args.partition_col, *depth_cols],
                   optional=[NEGPROBE_COL, AREA_COL, REGION_COL])
    print(f"Clustered: {len(obs):,} cells, "
          f"{obs[args.partition_col].nunique()} {args.partition_col} components")

    df = obs.join(typed, how="inner")
    if df.empty:
        sys.exit("ERROR: no cell ids in common between the AnnData and the typing result; "
                 "check that both come from the same run")
    if len(df) < len(typed):
        print(f"NOTE: {len(typed) - len(df):,} typed cells absent from the AnnData "
              f"(joined {len(df):,})")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    component_sizes = df[args.partition_col].astype(str).value_counts()

    for cluster in args.cluster:
        sub = df[df["cell_type"].astype(str) == cluster]
        if sub.empty:
            print(f"WARN: no cells labelled '{cluster}'; skipping", file=sys.stderr)
            continue

        rows = [summarize(df, depth_cols, BASELINE_LABEL, n_parent=0, n_component=None)]
        rows[0]["pct_of_cluster"] = np.nan

        parts = sub[args.partition_col].astype(str)
        shares = parts.value_counts()
        keep = shares[shares / len(sub) * 100 >= args.min_component_pct].index
        for name in keep:
            rows.append(summarize(sub[parts == name], depth_cols, name,
                                  n_parent=len(sub),
                                  n_component=int(component_sizes.get(name, 0))))
        tail = sub[~parts.isin(keep)]
        if not tail.empty:
            rows.append(summarize(tail, depth_cols, TAIL_LABEL,
                                  n_parent=len(sub), n_component=None))

        table = pd.DataFrame(rows)

        if args.reference_csv is not None:
            # Label every cell of this cluster with the component it was counted under, so
            # the rolled-up tail gets a profile too and the rows line up with the table.
            component_of = parts.where(parts.isin(keep), TAIL_LABEL)
            similarity = component_reference_similarity(
                args.clustered_h5ad, sub.index, component_of,
                args.reference_csv, args.scale_factor, args.top_reference)
            table = table.merge(similarity, on="component", how="left")

        out = args.out_dir / f"{cluster}_components.csv"
        table.to_csv(out, index=False)

        print(f"\n=== cluster '{cluster}': {len(sub):,} cells, "
              f"{parts.nunique()} {args.partition_col} components "
              f"({len(keep)} at >={args.min_component_pct}% of the cluster) ===")
        print(table.to_string(index=False))
        print(f"-> {out}")

    print("\nRead each component's depth against the __cohort__ row: a component at or above "
          "it is a population the reference lacks; one well below it is low-signal.")
    if args.reference_csv is not None:
        print("Then read 'nearest': components whose nearest reference types sit under "
              "different hierarchy branches are the ones the typing over-merged, and only "
              "those under the SAME branch belong in a new leaf for it.")


if __name__ == "__main__":
    main()
