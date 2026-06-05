#!/usr/bin/env python3
"""Concatenate per-slide CosMx AnnDatas, QC, and emit the scPearsonPCA inputs.

Stage 3a of the analysis pipeline. Reads the per-slide .h5ad files produced by
stage 1 (flatfiles_to_anndata.py), applies the per-cell quality filters from the
Bruker CosMx Scratch Space preprocessing vignette, concatenates the survivors, and
writes the two artifacts the rest of stage 3 needs:

  - <output>.h5ad        combined raw gene counts (sparse, QC-passed cells x ALL
                         gene probes), obs carrying the batch column + QC metrics +
                         all-gene total counts. Canonical record; input to the GPU
                         clustering stage (3c) and to downstream cell typing.

  - <output>.pca_input.h5   compact, R-friendly input to the Pearson-PCA stage (3b),
                         following the scPearsonPCA patient-batch vignette, which
                         passes an HVG subset PLUS precomputed all-gene statistics:
                           /counts        HVG counts as CSC (genes_hvg x cells) slots
                           /genes         HVG gene names
                           /cell_id       cell ids (column order of /counts)
                           /batch         per-cell batch label (e.g. patient "Case")
                           /total_counts  per-cell total over ALL gene probes (tc)
                           /genefreq/matrix  genes_hvg x batches gene-frequency
                                             (numerator = HVG counts per batch,
                                              denominator = ALL-gene counts per batch)
                           /genefreq/batch   batch labels (columns of the matrix)

Why this shape: scPearsonPCA::sparse_quasipoisson_pca_seurat_batch does the SVD on a
genes x genes cross-product, but it still materialises the genes x cells counts (and
a working copy) as a single R dgCMatrix, whose 32-bit index slots cap total nonzeros
at 2**31. Subsetting to ~2000 HVGs (the vignette's choice for the 6k panel) keeps
that matrix small enough, while tc and the per-batch gene frequency are computed here
over ALL genes (scipy uses 64-bit indices, so the full cohort is fine) and passed in
precomputed, exactly as the vignette does. The PCA never touches a dense matrix.

QC metrics are computed from X + var["probe_type"] (not named count columns) so the
thresholds apply regardless of CosMx export column naming. Only gene probes enter the
counts / tc; negprobes feed only the negprobe-proportion QC.

Memory: concatenating the cohort holds it in RAM — run on a big-memory CPU node.

Usage:
    uv run python pipeline/python/concat_qc_anndata.py \\
        --anndata-dir /path/to/per_slide_h5ads \\
        --output /path/to/combined_qc \\
        --batch-col Case --n-hvg 2000 \\
        --min-gene-counts 50 --max-area 30000 --max-negprobe-prop 0.1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# An R dgCMatrix stores its @i / @p index slots as 32-bit signed integers.
DGC_NNZ_LIMIT = 2**31 - 1

# Vignette defaults for the 6k panel.
DEFAULT_BATCH_COL = "Case"        # patient-level batch (Danaher: slide effects are minor)
DEFAULT_N_HVG = 2000              # FindVariableFeatures(nfeatures = 2000) for 6k
DEFAULT_MIN_GENE_COUNTS = 50      # per-cell gene-count floor (~50-100 for 6k)
DEFAULT_MAX_AREA_PX = 30000
DEFAULT_MAX_NEGPROBE_PROP = 0.1


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--anndata-dir", type=Path, required=True,
                   help="Directory of per-slide <slide>.h5ad files from stage 1.")
    p.add_argument("--output", type=Path, required=True,
                   help="Output prefix; writes <output>.h5ad and <output>.pca_input.h5.")
    p.add_argument("--batch-col", default=DEFAULT_BATCH_COL,
                   help="obs column used as the PCA batch variable (patient/Case).")
    p.add_argument("--cohort", type=Path, default=None,
                   help="Cohort CSV with Donor + Block columns (e.g. "
                        "pipeline/cohort_wenyu.csv). If given, keep only cells whose "
                        "(batch-col, Block) matches a cohort pair, dropping the "
                        "extraneous tissues co-mounted on each slide.")
    p.add_argument("--n-hvg", type=int, default=DEFAULT_N_HVG,
                   help="Number of highly variable genes (seurat_v3) for the PCA.")
    p.add_argument("--min-gene-counts", type=int, default=DEFAULT_MIN_GENE_COUNTS,
                   help="Drop cells whose summed gene-probe counts fall below this.")
    p.add_argument("--max-area", type=float, default=DEFAULT_MAX_AREA_PX,
                   help="Drop cells whose segmented area (px) exceeds this. "
                        "Skipped if no area column is present.")
    p.add_argument("--max-negprobe-prop", type=float,
                   default=DEFAULT_MAX_NEGPROBE_PROP,
                   help="Drop cells whose negprobe fraction of total counts exceeds this.")
    p.add_argument("--min-genes", type=int, default=0,
                   help="Optional floor on distinct genes detected per cell (0 disables).")
    p.add_argument("--area-col", default="Area",
                   help="obs column holding segmented cell area in pixels.")
    return p.parse_args()


def load_slides(anndata_dir: Path) -> ad.AnnData:
    paths = sorted(anndata_dir.glob("*.h5ad"))
    if not paths:
        print(f"ERROR: no .h5ad files under {anndata_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Reading {len(paths)} per-slide AnnDatas from {anndata_dir}")
    slides = []
    for path in paths:
        a = ad.read_h5ad(path)
        print(f"  {path.name}: {a.shape[0]} cells x {a.shape[1]} probes")
        slides.append(a)
    adata = ad.concat(slides, join="outer", merge="same", index_unique=None)
    if not adata.obs.index.is_unique:
        print("ERROR: concatenated obs index is not unique across slides",
              file=sys.stderr)
        sys.exit(1)
    print(f"Concatenated: {adata.shape[0]} cells x {adata.shape[1]} probes")
    return adata


def compute_qc(adata: ad.AnnData, area_col: str) -> pd.DataFrame:
    if "probe_type" not in adata.var:
        print("ERROR: var['probe_type'] missing; was this built by stage 1?",
              file=sys.stderr)
        sys.exit(1)
    gene_mask = (adata.var["probe_type"] == "gene").to_numpy()
    neg_mask = (adata.var["probe_type"] == "negprobe").to_numpy()
    X = adata.X.tocsr()

    total_all = np.asarray(X.sum(axis=1)).ravel()
    gene_counts = np.asarray(X[:, gene_mask].sum(axis=1)).ravel()
    neg_counts = np.asarray(X[:, neg_mask].sum(axis=1)).ravel()
    genes_detected = np.asarray((X[:, gene_mask] > 0).sum(axis=1)).ravel()

    with np.errstate(divide="ignore", invalid="ignore"):
        neg_prop = np.where(total_all > 0, neg_counts / total_all, 0.0)

    qc = pd.DataFrame(
        {
            "qc_gene_counts": gene_counts.astype(np.int64),
            "qc_genes_detected": genes_detected.astype(np.int64),
            "qc_negprobe_prop": neg_prop.astype(np.float32),
        },
        index=adata.obs.index,
    )
    if area_col in adata.obs:
        qc["qc_area"] = adata.obs[area_col].to_numpy()
    return qc


def qc_keep_mask(qc: pd.DataFrame, args: argparse.Namespace) -> np.ndarray:
    keep = qc["qc_gene_counts"].to_numpy() >= args.min_gene_counts
    keep &= qc["qc_negprobe_prop"].to_numpy() <= args.max_negprobe_prop
    if args.min_genes > 0:
        keep &= qc["qc_genes_detected"].to_numpy() >= args.min_genes
    if "qc_area" in qc:
        area = qc["qc_area"].to_numpy()
        keep &= ~(area > args.max_area)  # NaN areas are kept (filter in-evaluable)
    else:
        print(f"WARN: no '{args.area_col}' column; skipping the area filter",
              file=sys.stderr)
    return keep


def apply_cohort_filter(adata: ad.AnnData, cohort_path: Path,
                        batch_col: str) -> ad.AnnData:
    """Keep only cells whose (batch_col, Block) is a cohort (Donor, Block) pair.

    Two tissues are mounted per CosMx slide, so each slide carries extraneous
    donors/blocks outside the study cohort. The metadata Block is authoritative
    (e.g. donor 7464 is correctly 'A4' even though its slide name says 'A5'), so we
    match (Case, Block) directly against the cohort table — no slide-name parsing.
    """
    cohort = pd.read_csv(cohort_path, dtype=str)
    for col in ("Donor", "Block"):
        if col not in cohort.columns:
            print(f"ERROR: cohort CSV {cohort_path} missing '{col}' column",
                  file=sys.stderr)
            sys.exit(1)
    if "Block" not in adata.obs or batch_col not in adata.obs:
        print(f"ERROR: obs needs '{batch_col}' and 'Block' for cohort filtering; "
              f"have {list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)

    cohort_pairs = pd.MultiIndex.from_arrays(
        [cohort["Donor"].str.strip(), cohort["Block"].str.strip()])
    obs_pairs = pd.MultiIndex.from_arrays([
        adata.obs[batch_col].astype(str).str.strip(),
        adata.obs["Block"].astype(str).str.strip(),
    ])
    keep = obs_pairs.isin(cohort_pairs)

    n_before = adata.shape[0]
    kept_donors = set(adata.obs[batch_col].astype(str).str.strip()[keep].unique())
    cohort_donors = set(cohort["Donor"].str.strip().unique())
    absent = sorted(cohort_donors - kept_donors)
    print(f"Cohort filter: kept {int(keep.sum()):,} / {n_before:,} cells across "
          f"{len(kept_donors)} / {len(cohort_donors)} cohort donors")
    if absent:
        print(f"WARN: {len(absent)} cohort donor(s) absent from the data "
              f"(slides not yet processed?): {absent}", file=sys.stderr)
    return adata[keep].copy()


def select_hvgs(gene_adata: ad.AnnData, n_hvg: int) -> np.ndarray:
    """seurat_v3 HVGs on raw counts — the scanpy analog of FindVariableFeatures(vst)."""
    n_hvg = min(n_hvg, gene_adata.n_vars)
    sc.pp.highly_variable_genes(gene_adata, flavor="seurat_v3", n_top_genes=n_hvg)
    return gene_adata.var["highly_variable"].to_numpy()


def per_batch_gene_frequency(
    hvg_counts: sp.csr_matrix, tc: np.ndarray, batch: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Gene frequency per batch, vignette-style: HVG counts per batch over ALL-gene
    counts per batch. Returns (genes_hvg x batches matrix, batch labels)."""
    batches = pd.unique(batch)  # first-appearance order
    mat = np.zeros((hvg_counts.shape[1], len(batches)), dtype=np.float64)
    hvg_csc = hvg_counts.tocsc()
    for j, b in enumerate(batches):
        in_b = batch == b
        denom = float(tc[in_b].sum())  # all-gene total counts in this batch
        if denom <= 0:
            continue
        numer = np.asarray(hvg_csc[in_b].sum(axis=0)).ravel()  # per-HVG counts in batch
        mat[:, j] = numer / denom
    return mat, np.asarray(batches, dtype=object)


def _write_strings(grp, name: str, values: np.ndarray) -> None:
    # h5py cannot write numpy fixed-width unicode ('<U…') to a variable-length
    # string dtype ("No conversion path"). Pass an object array of Python str.
    grp.create_dataset(name, data=np.asarray(values, dtype=object),
                       dtype=h5py.string_dtype())


def write_pca_input(path: Path, hvg_counts: sp.csr_matrix, genes: np.ndarray,
                    cell_id: np.ndarray, batch: np.ndarray, tc: np.ndarray,
                    genefreq: np.ndarray, batch_labels: np.ndarray) -> None:
    # hvg_counts is cells x genes CSR; its (indptr, indices, data) ARE the
    # (@p, @i, @x) of a genes x cells CSC dgCMatrix — no transpose needed.
    hvg_counts.sort_indices()
    hvg_counts.eliminate_zeros()
    with h5py.File(path, "w") as f:
        g = f.create_group("counts")
        g.create_dataset("data", data=hvg_counts.data.astype(np.int32))
        g.create_dataset("indices", data=hvg_counts.indices.astype(np.int32))
        g.create_dataset("indptr", data=hvg_counts.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(
            [hvg_counts.shape[1], hvg_counts.shape[0]], dtype=np.int64))  # genes x cells
        _write_strings(f, "genes", genes)
        _write_strings(f, "cell_id", cell_id)
        _write_strings(f, "batch", batch)
        f.create_dataset("total_counts", data=tc.astype(np.float64))
        gf = f.create_group("genefreq")
        gf.create_dataset("matrix", data=genefreq)
        _write_strings(gf, "batch", batch_labels)


def main() -> None:
    args = parse_args()

    adata = load_slides(args.anndata_dir)

    if args.batch_col not in adata.obs:
        print(f"ERROR: batch column '{args.batch_col}' not in obs. Available: "
              f"{list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)

    qc = compute_qc(adata, args.area_col)
    adata.obs = adata.obs.join(qc)

    keep = qc_keep_mask(qc, args)
    n_before = adata.shape[0]
    adata = adata[keep].copy()
    print(f"QC: kept {adata.shape[0]} / {n_before} cells "
          f"({n_before - adata.shape[0]} dropped)")

    # Restrict to the study cohort before computing tc / HVGs / gene frequency, so
    # those statistics reflect only the cells that enter the PCA.
    if args.cohort:
        adata = apply_cohort_filter(adata, args.cohort, args.batch_col)

    # All-gene matrix: PCA, tc, and gene frequency all use gene probes only.
    gene_mask = (adata.var["probe_type"] == "gene").to_numpy()
    gene_adata = adata[:, gene_mask].copy()
    gene_X = gene_adata.X.tocsr()
    tc = np.asarray(gene_X.sum(axis=1)).ravel().astype(np.float64)  # all-gene totals
    adata.obs["total_counts"] = tc  # gene_adata shares adata's obs order
    print(f"{gene_adata.n_vars} gene probes; median tc = {np.median(tc):.0f}")

    hvg_mask = select_hvgs(gene_adata, args.n_hvg)
    print(f"Selected {int(hvg_mask.sum())} HVGs (seurat_v3)")
    hvg_counts = gene_X[:, hvg_mask]
    hvg_genes = gene_adata.var_names[hvg_mask].to_numpy()

    batch = adata.obs[args.batch_col].astype(str).to_numpy()
    genefreq, batch_labels = per_batch_gene_frequency(hvg_counts, tc, batch)
    print(f"Gene frequency over {len(batch_labels)} '{args.batch_col}' batches")

    hvg_counts = hvg_counts.tocsr()
    hvg_counts.sort_indices()
    hvg_counts.eliminate_zeros()
    if hvg_counts.nnz > DGC_NNZ_LIMIT:
        print(
            f"ERROR: HVG counts have {hvg_counts.nnz:,} nonzeros, over the R "
            f"dgCMatrix 32-bit limit ({DGC_NNZ_LIMIT:,}). Reduce --n-hvg or split "
            f"the cohort.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_h5ad = args.output.with_suffix(".h5ad")
    out_pca = args.output.with_suffix(".pca_input.h5")

    # Full all-probe matrix is the canonical record: 3c clusters on it and stage 4
    # InSituType needs the negprobes. var['probe_type'] lets consumers subset.
    print(f"Writing {out_h5ad} (all {adata.n_vars} probes, gene + neg + falsecode)")
    adata.write_h5ad(out_h5ad, compression="gzip")

    print(f"Writing {out_pca} ({hvg_counts.shape[0]:,} cells x "
          f"{hvg_counts.shape[1]} HVGs, {hvg_counts.nnz:,} nonzeros)")
    write_pca_input(
        out_pca, hvg_counts, hvg_genes,
        adata.obs.index.to_numpy(), batch, tc, genefreq, batch_labels,
    )

    print(f"Done. {adata.shape[0]:,} cells across {len(batch_labels)} "
          f"'{args.batch_col}' batches.")


if __name__ == "__main__":
    main()
