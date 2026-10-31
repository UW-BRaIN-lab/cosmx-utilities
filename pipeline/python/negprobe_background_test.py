#!/usr/bin/env python3
"""Is a de-novo letter's "carries instead" signature a cell state, or probe background?

75g found that letter `b` (778k cells, 31% of the anchor cohort, 62% of the Low_signal sink) is
NOT flat — at matched depth its cells are 1.44-1.72x MORE concentrated in their top genes than
the types they would be forced into. So they have a profile. The question this answers is what
that profile IS, because the genes carrying it look wrong for biology: SNRPA1, A1BG, SNRPD1,
VPS29, NDUFA3, PTEN, GPX1, F2RL3, B9D2 recur in `l`'s and `e`'s lists too, i.e. wherever real
signal is weak, while `b` is MISSING the panel's most abundant transcripts (TMSB4X, B2M, UBB).

TWO DIFFERENT "BACKGROUND" HYPOTHESES, and they predict opposite things:

  AMBIENT RNA — free transcript from lysed cells, assigned to whatever cell it lands in. Its
  composition mirrors the ambient pool, which is dominated by the ABUNDANT transcripts. An
  ambient-dominated cell is therefore ENRICHED in TMSB4X/B2M/UBB, not depleted. `b` is depleted,
  so ambient already looks like the wrong story — but it is worth measuring rather than asserting.

  PROBE BACKGROUND — non-specific binding and optical false calls. Its composition mirrors probe
  chemistry, not the transcriptome, and that is exactly what the negative-control probes measure.
  `/neg` is the per-cell mean negprobe count, so a gene whose signal is probe background should
  track `/neg` ACROSS CELLS once depth is held fixed.

So this reports three things:

  1. BACKGROUND BURDEN BY CELL TYPE. Per type: median depth, median `/neg`, and the ratio of a
     negprobe's count to an average gene's count. That last one is unitless and depth-normalised
     by construction, and it is the "how much of this population is background" number.
  2. PER-GENE BACKGROUND AFFINITY. Within one group of cells, correlation of each gene's
     log-norm expression with that cell's `/neg`, computed WITHIN DEPTH STRATA and pooled. Depth
     drives both expression and negprobe counts, so an unstratified correlation would be
     confounded for every gene at once.
  3. THE SIGNATURE AGAINST A MATCHED NULL. Mean affinity of the signature genes versus gene sets
     matched on abundance and detection rate. Affinity for every gene is computed once, so the
     null costs nothing. A contrast set (the abundant genes `b` is missing) is scored the same
     way, since the two hypotheses above separate on it.

READ: signature affinity well ABOVE its matched null = probe background, and the signature is an
artefact of weak signal rather than a cell state. Affinity INSIDE the null = the genes are not
tracking background, and the signature has to be explained some other way.

Inputs:
  --counts-h5  anchor_input.h5: /counts (CSC genes x cells), /genes, /cell_id, /neg
  --typing-h5  anchor_typing.h5: /cell_type per /cell_id
Outputs (--output-dir):
  background_by_celltype.csv, gene_background_affinity.csv, signature_vs_null.csv

Usage:
    uv run python pipeline/python/negprobe_background_test.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --letter b --output-dir out/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from anchor_profiles import read_cell_calls
from matched_null import (DEFAULT_NULL_NEIGHBOURS, draw_null_sets, empirical_p, gene_statistics,
                          load_counts, matched_pools)
from pseudobulk_core import DEFAULT_SCALE_FACTOR

# `b`'s recurring "carries instead" set from 75g, which also appears in `l`'s and `e`'s lists.
DEFAULT_SIGNATURE = ("SNRPA1,A1BG,SNRPD1,VPS29,NDUFA3,PTEN,GPX1,F2RL3,B9D2,H2AX,DUSP2,RBCK1,CES3")
# The panel's most abundant transcripts, which `b` is MISSING. Ambient RNA predicts the opposite,
# so these separate the two hypotheses.
DEFAULT_CONTRAST = "TMSB4X,B2M,UBB,COX1,COX2,HLA-A,HLA-B,ACTB"
DEFAULT_STRATA = 10
MIN_CELLS_PER_TYPE = 200
# Below this many cells a stratum's correlation is too noisy to pool.
MIN_STRATUM_CELLS = 50


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True)
    p.add_argument("--typing-h5", type=Path, required=True)
    p.add_argument("--letter", default="b",
                   help="Cell type whose cells the per-gene affinity is computed WITHIN "
                        "(default b). The background-burden table covers every type regardless.")
    p.add_argument("--signature-genes", default=DEFAULT_SIGNATURE,
                   help="The set under test (default: b's recurring 'carries instead' genes).")
    p.add_argument("--contrast-genes", default=DEFAULT_CONTRAST,
                   help="A second set scored the same way (default: the abundant transcripts b "
                        "is missing, which ambient RNA predicts should be ENRICHED instead).")
    p.add_argument("--depth-strata", type=int, default=DEFAULT_STRATA,
                   help=f"Depth quantile bins the correlation is computed within (default "
                        f"{DEFAULT_STRATA}). Depth drives expression and negprobe counts alike, "
                        f"so an unstratified correlation is confounded for every gene at once.")
    p.add_argument("--null-sets", type=int, default=200)
    p.add_argument("--null-neighbours", type=int, default=DEFAULT_NULL_NEIGHBOURS)
    p.add_argument("--null-seed", type=int, default=0)
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


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


def resolve_present(genes: np.ndarray, wanted: str, label: str) -> list[int]:
    """Panel positions of `wanted`, warning about absences rather than exiting.

    Unlike a programme under test, these sets are defaults carried from an older run and a gene
    that was renamed or pruned should not stop the job.
    """
    names = [g.strip() for g in wanted.split(",") if g.strip()]
    present = pd.Index(genes)
    idx = present.get_indexer(names)
    missing = [n for n, i in zip(names, idx) if i < 0]
    if missing:
        print(f"  WARNING: {label} gene(s) not on the panel, dropped: {', '.join(missing)}")
    kept = [int(i) for i in idx if i >= 0]
    if not kept:
        sys.exit(f"ERROR: none of the {label} genes is on the panel.")
    return kept


def background_by_celltype(totals: np.ndarray, neg: np.ndarray, cell_type: np.ndarray,
                           n_genes: int, min_cells: int) -> pd.DataFrame:
    """Per cell type: depth, negprobe level, and a negprobe's count over an average gene's.

    The ratio is the comparable number. `/neg` alone rises with depth like everything else, so a
    deep population can carry more background in absolute terms while being cleaner.
    """
    frame = pd.DataFrame({"cell_type": cell_type, "depth": totals, "neg": neg})
    frame = frame[frame["depth"] > 0]
    rows = []
    for name, grp in frame.groupby("cell_type", observed=True):
        if len(grp) < min_cells:
            continue
        # MEAN, not median. /neg is a mean count over ~20-30 negprobes and most cells score 0 on
        # all of them, so its median is 0 for nearly every population and ranks on it are ties.
        rows.append({
            "cell_type": name, "n_cells": len(grp),
            "median_depth": float(grp["depth"].median()),
            "mean_neg": float(grp["neg"].mean()),
            "pct_cells_with_neg": 100 * float((grp["neg"] > 0).mean()),
            "background_ratio": float(grp["neg"].mean() / (grp["depth"].mean() / n_genes)),
        })
    return pd.DataFrame(rows).sort_values("background_ratio", ascending=False)


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


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mat, genes, cell_id, totals = load_counts(args.counts_h5)
    neg = read_neg(args.counts_h5, cell_id)
    calls = read_cell_calls(args.typing_h5).set_index("cell_id")["cell_type"]
    cell_type = calls.reindex(cell_id).to_numpy()
    n_genes = mat.shape[0]
    print(f"{cell_id.size:,} cells x {n_genes:,} genes; "
          f"{int(pd.isna(cell_type).sum()):,} cells without a call")

    print("\n=== 1. background burden by cell type ===")
    burden = background_by_celltype(totals, neg, cell_type, n_genes, MIN_CELLS_PER_TYPE)
    print(f"  {'cell_type':30s} {'cells':>10s} {'depth':>8s} {'mean_neg':>9s} {'%>0':>6s} "
          f"{'bg_ratio':>9s}")
    shown = pd.concat([burden.head(8), burden[burden["cell_type"] == args.letter],
                       burden.tail(4)]).drop_duplicates("cell_type")
    for _, row in shown.iterrows():
        mark = "  <== letter" if row["cell_type"] == args.letter else ""
        print(f"  {str(row['cell_type'])[:28]:30s} {row['n_cells']:>10,} "
              f"{row['median_depth']:>8.0f} {row['mean_neg']:>9.4f} "
              f"{row['pct_cells_with_neg']:>6.1f} {row['background_ratio']:>9.3f}{mark}")
    if args.letter in set(burden["cell_type"]):
        rank = int((burden["cell_type"] == args.letter).argmax()) + 1
        print(f"  {args.letter} ranks {rank} of {len(burden)} types by background ratio "
              f"(1 = most background-dominated)")

    cells = np.flatnonzero(cell_type == args.letter)
    if cells.size < MIN_STRATUM_CELLS * args.depth_strata:
        sys.exit(f"ERROR: only {cells.size:,} cells labelled {args.letter!r}.")
    print(f"\n=== 2. per-gene background affinity within {cells.size:,} {args.letter} cells ===")
    affinity, used = background_affinity(mat, totals, neg, cells, args.depth_strata,
                                         args.scale_factor, MIN_STRATUM_CELLS)
    print(f"  pooled over {used} depth strata; "
          f"all-gene affinity median {np.median(affinity):+.4f}, "
          f"10th-90th {np.percentile(affinity, 10):+.4f} to {np.percentile(affinity, 90):+.4f}")

    mean_counts, detection = gene_statistics(mat)
    rank_pct = pd.Series(affinity).rank(pct=True).to_numpy()
    gene_table = pd.DataFrame({"gene": genes, "background_affinity": affinity,
                               "affinity_pct": 100 * rank_pct,
                               "mean_count": mean_counts, "detection_rate": detection})

    print("\n=== 3. signature and contrast against matched nulls ===")
    rng = np.random.default_rng(args.null_seed)
    rows = []
    for label, spec in (("signature", args.signature_genes), ("contrast", args.contrast_genes)):
        idx = np.asarray(resolve_present(genes, spec, label))
        pools = matched_pools(mean_counts, detection, idx, args.null_neighbours)
        draws = draw_null_sets(pools, args.null_sets, rng)
        observed = float(affinity[idx].mean())
        null = np.array([affinity[np.asarray(d)].mean() for d in draws])
        # One-sided the other way from a deficit: background predicts affinity ABOVE the null.
        p = empirical_p(-observed, -null)
        # Where the set sits in the WHOLE panel's affinity spread. The matched-null p says the
        # direction is reliable; it does not say the effect is large. Probe background should put
        # a set near the top of the panel, not mid-pack, and with millions of cells a 0.02
        # correlation difference is reliable and meaningless at the same time.
        panel_pct = 100 * float((affinity < observed).mean())
        rows.append({"set": label, "n_genes": idx.size, "observed_affinity": observed,
                     "null_median": float(np.median(null)),
                     "null_p95": float(np.percentile(null, 95)),
                     "p_value_above": p, "panel_percentile": panel_pct,
                     "genes": ",".join(genes[idx])})
        print(f"  {label:10s} {idx.size:>3} genes   affinity {observed:+.4f}   "
              f"null median {np.median(null):+.4f}   p(above) {p:.3f}   "
              f"panel pct {panel_pct:>5.1f}")
    print("  READ: probe background needs BOTH a matched-null p <= 0.05 AND a high panel")
    print("  percentile. The p alone only says the direction is reliable -- across millions of")
    print("  cells a 0.02 difference in correlation is reliable and trivial at once. A set that")
    print("  is above its matched peers but mid-pack on the panel is NOT background, and the")
    print("  signature then needs another explanation.")
    print("  The contrast set is the AMBIENT-RNA check: ambient predicts the abundant transcripts")
    print(f"  are ENRICHED in {args.letter}, so see the burden table and their affinity together.")

    burden.to_csv(args.output_dir / "background_by_celltype.csv", index=False)
    gene_table.sort_values("background_affinity", ascending=False).to_csv(
        args.output_dir / "gene_background_affinity.csv", index=False)
    pd.DataFrame(rows).to_csv(args.output_dir / "signature_vs_null.csv", index=False)
    print(f"\nWrote 3 tables to {args.output_dir}")


if __name__ == "__main__":
    main()
