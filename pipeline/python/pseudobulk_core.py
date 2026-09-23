#!/usr/bin/env python3
"""Shared pseudobulk + marker-selection machinery for the marker heatmaps.

marker_pseudobulk.py (Leiden x Region, from the cohort AnnData) and
denovo_vs_native_pseudobulk.py (de-novo letter vs its forced GBmap class, from the anchor
counts) run the same four steps on different inputs and groupings, so those steps live here:

  log_normalize   log1p(counts / per-cell total * scale_factor)  — Seurat LogNormalize
  onehot          cells x groups indicator, for mean-by-group without a dense matrix
  group_means     mean log-norm per group
  select_markers  top-N per group by (group mean - mean of the others), deduped
  zscore_rows     z-score each gene across the displayed groups

Everything here is pure: no argparse, no I/O.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

DEFAULT_SCALE_FACTOR = 10_000.0


def log_normalize(counts: sp.csr_matrix, scale_factor: float = DEFAULT_SCALE_FACTOR
                  ) -> sp.csr_matrix:
    """log1p(counts / per-cell total * scale_factor), sparse (Seurat LogNormalize)."""
    counts = counts.tocsr().astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    # Guard the division itself rather than np.where'ing after it: np.where evaluates both
    # branches, so an empty cell would warn about dividing by zero even though 0.0 wins.
    inv = np.divide(scale_factor, totals, out=np.zeros_like(totals, dtype=np.float64),
                    where=totals > 0)
    norm = sp.diags(inv) @ counts
    norm.data = np.log1p(norm.data)
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
    sums = np.asarray((oh.T @ norm).todense())
    sizes = np.asarray(oh.sum(axis=0)).ravel()
    return sums / np.where(sizes > 0, sizes, 1)[:, None]


def select_markers(profile: pd.DataFrame, select_order: list[str], top_n: int
                   ) -> tuple[list[str], dict[str, str]]:
    """Top-N markers per group by (group mean - mean of the other groups), deduped.

    profile is genes x groups of mean log-norm expression. The differential is computed
    against every column of `profile`, while `select_order` chooses which groups get their
    top-N picked and in what priority — so a subset of interest is not starved of markers by
    other groups claiming shared genes first. Returns the marker genes in selection order and
    {gene -> the group it was picked for} (for the heatmap's row split).
    """
    ordered: list[str] = []
    gene_to_group: dict[str, str] = {}
    for g in select_order:
        others = profile.drop(columns=g).mean(axis=1)
        diff = (profile[g] - others).sort_values(ascending=False)
        for gene in diff.index[:top_n]:
            if gene not in gene_to_group:
                ordered.append(gene)
                gene_to_group[gene] = g
    return ordered, gene_to_group


def zscore_rows(pb: pd.DataFrame) -> pd.DataFrame:
    """Z-score each row across columns; constant rows go flat rather than NaN."""
    mean = pb.mean(axis=1)
    sd = pb.std(axis=1, ddof=1)
    out = pb.sub(mean, axis=0).div(sd.where(sd > 0, 1.0), axis=0)
    out[sd <= 0] = 0.0
    return out
