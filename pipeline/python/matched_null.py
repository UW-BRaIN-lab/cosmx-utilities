#!/usr/bin/env python3
"""Matched-null gene sets: the control that separates a real deficit from shallow cells.

A gene set's mean expression is compared between two groups of cells all over this pipeline, and
a negative difference is read as "this group lacks that programme". That reading is only safe
once DEPTH has been ruled out. Scores here are library-size normalised, which removes the
first-order effect, but normalisation cannot recover a transcript that was never detected: a
shallower group drops lightly expressed genes to zero more often, a zero contributes 0 to the
mean, and the difference goes negative for genes of ANY biology. The de-novo letters run as low
as 0.37x their comparators' depth, so this is the live confound rather than a hypothetical one.

The control is to re-run the same measurement on gene sets matched to the real one on cohort
abundance AND detection rate. Those sets carry the dropout behaviour and none of the biology, so
where the observed value falls among them is the answer: inside the null means depth, outside it
means the programme.

A THIRD AXIS. 75q measured, per gene, how strongly its expression tracks the per-cell negprobe
load. Genes at the bottom of that spread deflate when background inflates a cell's library size,
because the normalisation divides by a total that includes the background. CDK4, MDM2 and OS9 --
three of the five 12q amplicon genes -- sit in the panel's bottom 0.25%, so a population carrying
more background than its comparator shows a deficit in them for reasons that have nothing to do
with copy number. Abundance and detection do not capture that, and matching cannot fix it: only ~15 genes on a
6,175-gene panel are that far down, so any workable candidate pool reaches back into mid-pack.
A correction for this is UNSOLVED -- see test_background_dilution_fools_the_two_axis_null, which
documents the failure mode rather than a fix. Before reaching for one, check whether the
population even carries more background than its comparator (75q block 1); if it does not, the
mechanism never fires.

Used by diagnose_program_by_section.py (75e) and negprobe_background_test.py (75q).
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from anchor_profiles import decode

# Guards log10 of a gene that is never detected in the cohort.
MEAN_COUNT_EPS = 1e-6
DEFAULT_NULL_NEIGHBOURS = 100


def load_counts(counts_h5: Path) -> tuple[sp.csc_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Load the genes x cells CSC once, with per-cell library sizes.

    The matched-null control scores hundreds of gene sets against the same matrix, so the file is
    read once here instead of once per set.
    """
    with h5py.File(counts_h5, "r") as f:
        n_genes, n_cells = (int(x) for x in f["counts/shape"][()])
        mat = sp.csc_matrix(
            (np.asarray(f["counts/data"][()], dtype=np.float64),
             np.asarray(f["counts/indices"][()], dtype=np.int64),
             np.asarray(f["counts/indptr"][()], dtype=np.int64)),
            shape=(n_genes, n_cells))
        genes = np.asarray(decode(f["genes"][()]))
        cell_id = np.asarray(decode(f["cell_id"][()]))
    totals = np.asarray(mat.sum(axis=0)).ravel()
    return mat, genes, cell_id, totals


def resolve_genes(genes: np.ndarray, program: list[str]) -> np.ndarray:
    """Panel positions of the programme genes, failing loudly on anything off-panel."""
    missing = [g for g in program if g not in set(genes)]
    if missing:
        sys.exit(f"ERROR: programme gene(s) not on the panel: {missing}")
    return pd.Index(genes).get_indexer(program)


def score_gene_set(mat: sp.csc_matrix, totals: np.ndarray, idx: np.ndarray,
                   scale_factor: float) -> np.ndarray:
    """Mean log-norm expression of the genes at `idx`, one value per cell, in column order.

    Only those rows are ever materialised, so this stays light at cohort scale and is cheap
    enough to call once per null set.
    """
    sub = np.asarray(mat[idx, :].todense())
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.log1p(np.divide(sub * scale_factor, totals,
                                  out=np.zeros_like(sub), where=totals > 0))
    return norm.mean(axis=0)


def gene_statistics(mat: sp.csc_matrix) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene mean raw count and detection rate, read off the CSC's nonzeros alone.

    `mat` is genes x cells and column-compressed, so `indices` already holds the GENE of every
    nonzero: one bincount gives the detected-cell count and a weighted one gives the total, with
    no transpose and nothing densified.
    """
    n_genes, n_cells = mat.shape
    detected = np.bincount(mat.indices, minlength=n_genes).astype(np.float64)
    summed = np.bincount(mat.indices, weights=mat.data, minlength=n_genes)
    return summed / n_cells, detected / n_cells


def matched_pools(mean_counts: np.ndarray, detection: np.ndarray, program_idx: np.ndarray,
                  n_neighbours: int) -> list[np.ndarray]:
    """For each programme gene, the `n_neighbours` genes most like it in abundance + detection.

    DEPTH is the confound this exists for. A population sequenced shallower drops a lightly
    expressed gene to zero more often, and a zero contributes 0 to the mean score — so ANY such
    gene shows a deficit in that population whatever its biology. Matching the null on detection
    rate as well as abundance is what lets a null set reproduce that artefact: if it does, the
    programme's own deficit is a statement about depth, not about the programme.
    """
    axes = np.vstack([np.log10(mean_counts + MEAN_COUNT_EPS), detection])
    axes = (axes - axes.mean(axis=1, keepdims=True)) / axes.std(axis=1, keepdims=True)
    pools = []
    for gene in program_idx:
        distance = np.linalg.norm(axes - axes[:, [gene]], axis=0)
        distance[program_idx] = np.inf          # a programme gene can never stand in for another
        pools.append(np.argsort(distance)[:n_neighbours])
    return pools


def draw_null_sets(pools: list[np.ndarray], n_sets: int,
                   rng: np.random.Generator) -> list[list[int]]:
    """One matched stand-in per programme gene per set, never the same gene twice in a set."""
    sets = []
    for _ in range(n_sets):
        chosen: list[int] = []
        for pool in pools:
            free = [g for g in pool if g not in chosen]
            if not free:
                sys.exit("ERROR: a matched pool was exhausted; raise --null-neighbours.")
            chosen.append(int(rng.choice(free)))
        sets.append(chosen)
    return sets


def empirical_p(observed: float, null: np.ndarray) -> float:
    """One-sided p that a deficit this large comes from matched genes alone.

    (hits + 1) / (n + 1): the observed set is itself one draw from the null hypothesis, so a p of
    exactly 0 is not available however many sets are drawn.
    """
    null = null[~np.isnan(null)]
    if not len(null) or np.isnan(observed):
        return float("nan")
    return float((np.sum(null <= observed) + 1) / (len(null) + 1))


def read_neg(counts_h5: Path, cell_id: np.ndarray) -> np.ndarray:
    """Per-cell mean negprobe count, in the matrix's column order."""
    with h5py.File(counts_h5, "r") as f:
        if "neg" not in f:
            sys.exit(f"ERROR: /neg missing in {counts_h5}; it is written by "
                     f"split_insitutype_inputs.py and is what this test measures against.")
        neg = np.asarray(f["neg"][()], dtype=np.float64)
    if neg.size != cell_id.size:
        sys.exit(f"ERROR: /neg has {neg.size:,} entries for {cell_id.size:,} cells.")
    return neg


def background_affinity(mat, totals: np.ndarray, neg: np.ndarray, cells: np.ndarray,
                        n_strata: int, scale_factor: float,
                        min_stratum: int) -> tuple[np.ndarray, int]:
    """Per-gene correlation of log-norm expression with `/neg`, pooled over depth strata.

    Accumulates sums off each stratum's nonzeros, so the only dense objects are per-gene
    accumulators. Zeros are carried correctly: they contribute nothing to the sums and are
    accounted for by dividing by the stratum's full cell count.
    """
    n_genes = mat.shape[0]
    sub_totals = totals[cells]
    edges = np.quantile(sub_totals, np.linspace(0, 1, n_strata + 1))
    stratum = np.clip(np.searchsorted(edges, sub_totals, side="right") - 1, 0, n_strata - 1)

    weighted_r = np.zeros(n_genes)
    weight = 0
    used = 0
    for s in range(n_strata):
        take = cells[stratum == s]
        if take.size < min_stratum:
            continue
        block = mat[:, take].tocsc()
        n = take.size
        col = np.repeat(np.arange(n), np.diff(block.indptr))
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.log1p(block.data * scale_factor / totals[take][col])
        neg_s = neg[take]
        sum_x = np.bincount(block.indices, weights=x, minlength=n_genes)
        sum_x2 = np.bincount(block.indices, weights=x * x, minlength=n_genes)
        sum_xn = np.bincount(block.indices, weights=x * neg_s[col], minlength=n_genes)

        mean_x, mean_n = sum_x / n, neg_s.mean()
        var_x = np.maximum(sum_x2 / n - mean_x ** 2, 0.0)
        var_n = max(float(neg_s.var()), 0.0)
        denom = np.sqrt(var_x * var_n)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(denom > 0, (sum_xn / n - mean_x * mean_n) / denom, 0.0)
        weighted_r += r * n
        weight += n
        used += 1
    if not weight:
        sys.exit("ERROR: no depth stratum had enough cells; lower --depth-strata.")
    return weighted_r / weight, used
