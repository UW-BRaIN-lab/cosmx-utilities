#!/usr/bin/env python3
"""Per-cluster marker selection + pseudobulk z-score matrix for a marker heatmap.

The compute half of the marker-gene heatmap (the render half is pipeline/R/
marker_heatmap.R). Adapted from the lab's InSituType marker_heatmap.R, but driven
by the Stage 3c Leiden clusters instead of InSituType cell types, and reading the
cohort AnnData directly rather than a small local counts.mtx.

Reads cosmx_clustered.h5ad (gene + neg + falsecode probes, obs has the cluster key,
Region, Case) and, on gene probes only:

  1. Log-normalizes: log1p(counts / per-cell gene total * scale_factor)  (= Seurat
     LogNormalize / scanpy normalize_total+log1p), done sparsely.
  2. Builds a per-cluster mean log-norm profile (genes x clusters) via a one-hot
     cluster matrix — no dense cells x genes matrix is ever formed.
  3. Picks the top-N markers per cluster by differential score
     (cluster_mean - mean(other clusters)), deduped in cluster order -- OR, with
     --genes-csv, skips selection entirely and uses a curated gene list.
  4. Pseudobulks the marker genes by cluster x Region (mean log-norm per group),
     drops groups with < min_group_n cells, and z-scores each gene across groups.
     --no-region-split pseudobulks by cluster alone (columns are bare clusters).

Writes small CSVs for the R renderer:
  <out>/marker_heatmap_zmatrix.csv   rows = markers, cols = "<cluster> | <Region>"
                                     (or bare "<cluster>" with --no-region-split)
  <out>/top_markers_per_cluster.csv  columns: cluster, gene  (row-split order; with
                                     --gene-group-column the "cluster" column holds
                                     the curated module instead of a marked cluster)
  <out>/group_sizes.csv              columns: group, n_cells (incl. dropped flag)
  <out>/gene_detection_rate.csv      rows = genes, cols = groups; fraction of cells
                                     with a nonzero count. Z-scores of a gene seen in
                                     ~0% of cells are noise -- read them together.
  <out>/genes_missing.csv            curated genes absent from the panel (--genes-csv)

Swap --group-key from leiden to an InSituType cell_type column once Stage 4 lands.

Usage:
    uv run python pipeline/python/marker_pseudobulk.py \\
        --clustered-h5ad cosmx_clustered.h5ad --output-dir marker_heatmap \\
        --group-key leiden --top-n 5 --min-group-n 10

    # curated gene panel, rows split by module, columns = Leiden clusters
    uv run python pipeline/python/marker_pseudobulk.py \\
        --clustered-h5ad cosmx_typed.h5ad --output-dir synaptic_heatmap \\
        --group-key leiden --no-region-split \\
        --genes-csv cosmx_channel_synaptic_gene_audit.csv --gene-group-column module
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Tumor -> edge -> normal; column order in the heatmap and the default region set.
REGION_ORDER = ["Tumor bulk", "Infiltrating edge", "Contralateral uninvolved"]
DEFAULT_TOP_N = 5
DEFAULT_MIN_GROUP_N = 10
DEFAULT_SCALE_FACTOR = 1e4
# Below this cohort-wide per-cell detection rate a gene's pseudobulk z-score is driven
# by a handful of counts; flagged in the log rather than silently plotted.
MIN_INFORMATIVE_DETECTION = 0.01


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="cosmx_clustered.h5ad from stage 3c.")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Directory for the heatmap CSVs.")
    p.add_argument("--group-key", default="leiden",
                   help="obs column to group cells by (default: leiden clusters).")
    p.add_argument("--region-key", default="Region",
                   help="obs column with the tissue region.")
    p.add_argument("--regions", nargs="+", default=None,
                   help="Regions to include (default: all of REGION_ORDER present).")
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                   help="Top markers per cluster.")
    p.add_argument("--min-group-n", type=int, default=DEFAULT_MIN_GROUP_N,
                   help="Drop cluster x Region groups with fewer than this many cells.")
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR,
                   help="LogNormalize scale factor.")
    p.add_argument("--clusters", default=None,
                   help="Comma-separated subset of group-key values to SELECT markers "
                        "for and show as heatmap columns (e.g. the de novo InSituType "
                        "letters 'a,b,c,...'). Differential scores are still computed "
                        "against ALL clusters, so the markers are genuinely "
                        "distinguishing; this only restricts which clusters get their "
                        "top-N picked (and dedup priority) and which appear in the "
                        "heatmap. Default: all clusters.")
    p.add_argument("--genes-csv", type=Path, default=None,
                   help="CSV of curated genes to plot INSTEAD of data-driven top-N "
                        "markers. Genes absent from the panel are reported and skipped; "
                        "row order follows the CSV.")
    p.add_argument("--gene-column", default="gene",
                   help="Column of --genes-csv holding the gene symbol.")
    p.add_argument("--gene-group-column", default=None,
                   help="Column of --genes-csv to split heatmap rows by (e.g. 'module'). "
                        "Default: one unsplit block.")
    p.add_argument("--no-region-split", action="store_true",
                   help="Pseudobulk by --group-key alone; columns are bare cluster "
                        "labels with no Region sub-split (and --region-key is ignored).")
    return p.parse_args()


def load_curated_genes(csv_path: Path, gene_col: str,
                       group_col: str | None) -> tuple[list[str], dict[str, str]]:
    """Read a curated gene list -> (genes in file order, gene -> row-split group)."""
    table = pd.read_csv(csv_path)
    for col in (gene_col, *( [group_col] if group_col else [] )):
        if col not in table.columns:
            print(f"ERROR: --genes-csv has no column '{col}'. "
                  f"Available: {list(table.columns)}", file=sys.stderr)
            sys.exit(1)
    table = table[table[gene_col].notna()]
    genes, gene_to_group = [], {}
    for _, row in table.iterrows():
        gene = str(row[gene_col]).strip()
        if not gene or gene in gene_to_group:
            continue
        genes.append(gene)
        gene_to_group[gene] = str(row[group_col]).strip() if group_col else "curated"
    return genes, gene_to_group


def log_normalize(counts: sp.csr_matrix, scale_factor: float) -> sp.csr_matrix:
    """log1p(counts / per-cell total * scale_factor), sparse (Seurat LogNormalize)."""
    counts = counts.tocsr().astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    inv = np.where(totals > 0, scale_factor / totals, 0.0)
    norm = sp.diags(inv) @ counts          # row-scale to scale_factor
    norm = norm.tocsr()
    norm.data = np.log1p(norm.data)        # log1p(0)=0, so only stored entries change
    return norm


def onehot(labels: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    """Cells x categories one-hot (sparse) + the category labels (first-seen order)."""
    cats = pd.Categorical(labels)
    codes = cats.codes
    n, k = len(codes), len(cats.categories)
    oh = sp.csr_matrix((np.ones(n), (np.arange(n), codes)), shape=(n, k))
    return oh, np.asarray(cats.categories)


def group_means(norm: sp.csr_matrix, oh: sp.csr_matrix) -> np.ndarray:
    """Mean of `norm` rows within each one-hot group -> (n_groups x n_genes) dense."""
    sums = np.asarray((oh.T @ norm).todense())   # n_groups x n_genes
    sizes = np.asarray(oh.sum(axis=0)).ravel()
    return sums / sizes[:, None]


def main() -> None:
    args = parse_args()

    print(f"Reading {args.clustered_h5ad}")
    adata = ad.read_h5ad(args.clustered_h5ad)
    split_region = not args.no_region_split
    required_keys = (args.group_key, args.region_key) if split_region else (args.group_key,)
    for key in required_keys:
        if key not in adata.obs:
            print(f"ERROR: obs is missing '{key}'. Available: {list(adata.obs.columns)}",
                  file=sys.stderr)
            sys.exit(1)
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing.", file=sys.stderr)
        sys.exit(1)

    # Region filter (default: all of REGION_ORDER that are present). Skipped entirely
    # for --no-region-split: the subset .copy() is the peak-memory step at cohort scale.
    if split_region:
        region = adata.obs[args.region_key].astype(str).to_numpy()
        wanted = args.regions or [r for r in REGION_ORDER if r in set(region)]
        keep = np.isin(region, wanted)
        adata = adata[keep].copy()
        region = adata.obs[args.region_key].astype(str).to_numpy()
        print(f"{adata.n_obs:,} cells across regions {wanted}")
    else:
        region, wanted = None, []
        print(f"{adata.n_obs:,} cells (no Region split)")
    cluster = adata.obs[args.group_key].astype(str).to_numpy()

    gene_mask = (adata.var["probe_type"] == "gene").to_numpy()
    genes = adata.var_names[gene_mask].to_numpy()
    print(f"Log-normalizing {int(gene_mask.sum())} gene probes")
    norm = log_normalize(adata.X.tocsr()[:, gene_mask], args.scale_factor)

    # Per-cluster mean log-norm profile (genes x clusters).
    cl_oh, cl_labels = onehot(cluster)
    profile = group_means(norm, cl_oh).T          # genes x clusters
    profile = pd.DataFrame(profile, index=genes, columns=cl_labels)

    # Order clusters numerically when they look like integers (leiden), else lexically.
    try:
        cl_order = [str(c) for c in sorted(cl_labels, key=lambda x: int(x))]
    except ValueError:
        cl_order = sorted(cl_labels)
    profile = profile[cl_order]

    # Which clusters to select markers for / show. The differential is always computed
    # against ALL clusters (profile keeps every column); --clusters only narrows which
    # clusters get their top-N picked (and the dedup priority) and which appear in the
    # heatmap, so a subset — e.g. the de novo InSituType letters — isn't starved of
    # markers by the named types claiming shared genes first.
    if args.clusters:
        requested = {c.strip() for c in args.clusters.split(",") if c.strip()}
        missing = sorted(requested - set(cl_order))
        if missing:
            print(f"ERROR: --clusters values not in '{args.group_key}': {missing}. "
                  f"Available: {cl_order}", file=sys.stderr)
            sys.exit(1)
        select_order = [c for c in cl_order if c in requested]
        print(f"Selecting markers for {len(select_order)} requested clusters "
              f"(differential still vs all {len(cl_order)})")
    else:
        select_order = cl_order

    args.output_dir.mkdir(parents=True, exist_ok=True)
    missing: list[str] = []
    if args.genes_csv:
        requested_genes, gene_to_group = load_curated_genes(
            args.genes_csv, args.gene_column, args.gene_group_column)
        on_panel = set(genes)
        ordered_markers = [g for g in requested_genes if g in on_panel]
        missing = [g for g in requested_genes if g not in on_panel]
        gene_to_cluster = {g: gene_to_group[g] for g in ordered_markers}
        print(f"Curated list: {len(ordered_markers)}/{len(requested_genes)} genes on panel")
        if missing:
            print(f"  {len(missing)} NOT on the panel (skipped): {', '.join(missing)}")
        if not ordered_markers:
            print("ERROR: no curated gene is on the panel.", file=sys.stderr)
            sys.exit(1)
        pd.DataFrame({"gene": missing,
                      "group": [gene_to_group[g] for g in missing]}).to_csv(
            args.output_dir / "genes_missing.csv", index=False)
    else:
        print(f"Selecting top {args.top_n} markers per cluster")
        ordered_markers = []
        gene_to_cluster = {}
        for c in select_order:
            others_mean = profile.drop(columns=c).mean(axis=1)
            diff = (profile[c] - others_mean).sort_values(ascending=False)
            for g in diff.index[: args.top_n]:
                if g not in gene_to_cluster:
                    ordered_markers.append(g)
                    gene_to_cluster[g] = c

    # Pseudobulk the marker genes by cluster x Region, over the selected clusters' cells.
    sel_mask = np.isin(cluster, select_order)
    marker_idx = pd.Index(genes).get_indexer(ordered_markers)
    if split_region:
        group = np.char.add(np.char.add(cluster[sel_mask].astype(str), " | "),
                            region[sel_mask].astype(str))
    else:
        group = cluster[sel_mask].astype(str)
    grp_oh, grp_labels = onehot(group)
    sub = norm[sel_mask][:, marker_idx]                         # cells x markers
    pb = group_means(sub, grp_oh).T                             # markers x groups
    pb = pd.DataFrame(pb, index=ordered_markers, columns=grp_labels)
    grp_sizes = pd.Series(np.asarray(grp_oh.sum(axis=0)).ravel(), index=grp_labels)

    # Per-group detection rate (fraction of cells with a nonzero count). log1p keeps the
    # sparsity pattern, so the binarised log-norm matrix is the raw detection indicator.
    detected = sub.copy()
    detected.data = np.ones_like(detected.data)
    det = pd.DataFrame(group_means(detected, grp_oh).T,
                       index=ordered_markers, columns=grp_labels)

    small = grp_sizes[grp_sizes < args.min_group_n].index.tolist()
    if small:
        print(f"Dropping {len(small)} groups with < {args.min_group_n} cells")
        pb = pb.drop(columns=small)
        det = det.drop(columns=small)

    # Order columns by (cluster, region); z-score each marker (row) across groups.
    region_rank = {r: i for i, r in enumerate(wanted)}
    cl_rank = {c: i for i, c in enumerate(select_order)}
    def _col_key(g: str):
        c, _, r = g.partition(" | ")
        return (cl_rank.get(c, len(cl_rank)), region_rank.get(r, len(region_rank)))
    pb = pb[sorted(pb.columns, key=_col_key)]
    det = det[pb.columns]

    mean = pb.mean(axis=1)
    sd = pb.std(axis=1, ddof=1)
    pb_z = pb.sub(mean, axis=0).div(sd.where(sd > 0, 1.0), axis=0)
    pb_z[sd <= 0] = 0.0  # constant markers -> flat (avoid NaN)

    zpath = args.output_dir / "marker_heatmap_zmatrix.csv"
    mpath = args.output_dir / "top_markers_per_cluster.csv"
    spath = args.output_dir / "group_sizes.csv"
    print(f"Writing {zpath}")
    pb_z.to_csv(zpath)
    pd.DataFrame({"cluster": [gene_to_cluster[g] for g in ordered_markers],
                  "gene": ordered_markers}).to_csv(mpath, index=False)
    grp_sizes.rename("n_cells").rename_axis("group").reset_index().assign(
        dropped=lambda d: d["group"].isin(small)).to_csv(spath, index=False)
    dpath = args.output_dir / "gene_detection_rate.csv"
    print(f"Writing {dpath}")
    det.to_csv(dpath)

    # A z-score is only as meaningful as the counts under it: call out genes the panel
    # technically carries but essentially never detects, so they aren't over-read.
    overall_det = det.mul(grp_sizes[det.columns], axis=1).sum(axis=1) / grp_sizes[det.columns].sum()
    sparse = overall_det[overall_det < MIN_INFORMATIVE_DETECTION].sort_values()
    if len(sparse):
        print(f"NOTE: {len(sparse)} gene(s) detected in < "
              f"{MIN_INFORMATIVE_DETECTION:.1%} of cells -- their z-scores are noise: "
              f"{', '.join(sparse.index)}")

    region_note = f" x {len(wanted)} regions" if split_region else ""
    print(f"Done. {pb_z.shape[0]} genes x {pb_z.shape[1]} groups "
          f"({len(select_order)} clusters{region_note}).")


if __name__ == "__main__":
    main()
