#!/usr/bin/env python3
"""Is a de-novo letter's gene programme a cell state, or a slide-level handling artefact?

75d showed that de-novo `t` is separated from every natively-called GBmap class by a heat-shock
programme (HSPA1A/HSPA1B/HSPB1/DNAJB1/HSP90AA1/HSPH1). Heat-shock is exactly the signature warm
ischaemia and slow fixation produce, so before calling it a tumour state we have to rule out
pre-analytical handling. Same machinery works for any letter and any gene set.

The discriminating comparison is WITHIN a slide, not across slides:

  * If the programme is a CELL STATE, the letter's cells carry it wherever they occur — the
    letter-vs-other gap is positive on slide after slide, and the two groups' slide-level scores
    are only loosely coupled.
  * If it is HANDLING, the whole section is elevated — the letter's cells and everyone else's
    rise and fall together, so the slide-level scores are tightly correlated and the within-slide
    gap collapses toward zero.

So this reports, per slide: the programme score in the letter's cells, the score in every other
cell on that same slide, and the gap. Plus three summaries — the share of slides with a positive
gap, the correlation between the two groups' per-slide scores, and how much of the total score
variance sits BETWEEN slides versus between the groups within a slide.

Cheap: only the programme genes and per-cell totals are needed, so the counts matrix is never
transposed or densified.

Inputs:
  --counts-h5    anchor_input.h5: /counts (CSC genes x cells), /genes, /cell_id.
  --typing-h5    anchor_typing.h5: semi-supervised /cell_type per /cell_id.
  --cells-csv    anchor_cells.csv from 71 (cell_id, slide_id, cluster).
  --letter       the de-novo letter under test, e.g. t
  --genes        comma-separated programme genes (default: the heat-shock set 75d surfaced).
Outputs:
  --output-csv   per-slide table; the summaries print to stdout.

Usage:
    uv run python pipeline/python/diagnose_program_by_slide.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --cells-csv anchor_cells.csv --letter t --output-csv t_heatshock_by_slide.csv
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
from pseudobulk_core import DEFAULT_SCALE_FACTOR

# The programme 75d surfaced as separating `t` from every natively-called class.
DEFAULT_PROGRAM = "HSPA1A,HSPA1B,HSPB1,DNAJB1,HSP90AA1,HSPH1"
MIN_SLIDE_CELLS = 50


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True)
    p.add_argument("--typing-h5", type=Path, required=True)
    p.add_argument("--cells-csv", type=Path, required=True,
                   help="anchor_cells.csv: cell_id + slide_id.")
    p.add_argument("--letter", required=True)
    p.add_argument("--genes", default=DEFAULT_PROGRAM,
                   help=f"Comma-separated programme genes (default: {DEFAULT_PROGRAM}).")
    p.add_argument("--slide-key", default="slide_id")
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--min-slide-cells", type=int, default=MIN_SLIDE_CELLS,
                   help="Skip slides with fewer than this many cells in either group.")
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    return p.parse_args()


def program_score(counts_h5: Path, program: list[str], scale_factor: float
                  ) -> pd.Series:
    """Mean log-norm expression of `program` per cell, without transposing the matrix.

    Per-cell totals come from column sums of the CSC (genes x cells), and only the programme
    genes' rows are ever materialised — so this stays light at cohort scale.
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

    missing = [g for g in program if g not in set(genes)]
    if missing:
        sys.exit(f"ERROR: programme gene(s) not on the panel: {missing}")
    idx = pd.Index(genes).get_indexer(program)

    totals = np.asarray(mat.sum(axis=0)).ravel()          # per-cell library size
    sub = np.asarray(mat[idx, :].todense())               # len(program) x cells
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.log1p(np.divide(sub * scale_factor, totals,
                                  out=np.zeros_like(sub), where=totals > 0))
    return pd.Series(norm.mean(axis=0), index=cell_id, name="score")


def main() -> None:
    args = parse_args()
    program = [g.strip() for g in args.genes.split(",") if g.strip()]
    if not program:
        sys.exit("ERROR: --genes is empty.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]].set_index("cell_id")
    cells = pd.read_csv(args.cells_csv)
    id_col = "cell_id" if "cell_id" in cells.columns else cells.columns[0]
    if args.slide_key not in cells.columns:
        sys.exit(f"ERROR: {args.cells_csv} has no '{args.slide_key}'. "
                 f"Columns: {list(cells.columns)}")
    slide = cells.set_index(id_col)[args.slide_key].astype(str)

    print(f"Scoring the {len(program)}-gene programme: {', '.join(program)}")
    score = program_score(args.counts_h5, program, args.scale_factor)

    df = pd.DataFrame({"score": score})
    df["cell_type"] = calls["cell_type"].reindex(df.index)
    df["slide"] = slide.reindex(df.index)
    df = df.dropna(subset=["cell_type", "slide"])
    if df.empty:
        sys.exit("ERROR: no cells left after joining labels and slides — check the ids.")
    df["is_letter"] = df["cell_type"] == args.letter
    if not df["is_letter"].any():
        sys.exit(f"ERROR: no cells labelled {args.letter!r}.")
    print(f"{len(df):,} cells, {df['slide'].nunique()} slides, "
          f"{int(df['is_letter'].sum()):,} labelled {args.letter}")

    grouped = df.groupby(["slide", "is_letter"], observed=True)["score"].agg(["mean", "size"])
    tbl = grouped.unstack("is_letter")
    tbl.columns = [f"{a}_{'letter' if b else 'other'}" for a, b in tbl.columns]
    tbl = tbl.rename(columns={"mean_letter": "score_letter", "mean_other": "score_other",
                              "size_letter": "n_letter", "size_other": "n_other"})
    for col in ("score_letter", "score_other", "n_letter", "n_other"):
        if col not in tbl:
            tbl[col] = np.nan
    tbl["n_letter"] = tbl["n_letter"].fillna(0).astype(int)
    tbl["n_other"] = tbl["n_other"].fillna(0).astype(int)
    tbl["gap"] = tbl["score_letter"] - tbl["score_other"]
    tbl["letter_pct"] = (100 * tbl["n_letter"] / (tbl["n_letter"] + tbl["n_other"])).round(1)

    usable = tbl[(tbl["n_letter"] >= args.min_slide_cells)
                 & (tbl["n_other"] >= args.min_slide_cells)].copy()
    if usable.empty:
        sys.exit(f"ERROR: no slide has >= {args.min_slide_cells} cells in both groups.")

    pos = float((usable["gap"] > 0).mean())
    r = float(usable["score_letter"].corr(usable["score_other"]))
    # Variance sitting BETWEEN slides (both groups moving together = handling) versus the
    # within-slide letter-vs-other gap (a cell state that travels with the cells).
    between = float(usable[["score_letter", "score_other"]].mean(axis=1).var(ddof=1))
    within = float((usable["gap"] / 2).pow(2).mean())

    print(f"\n=== {args.letter}: programme score, {len(usable)} usable slides "
          f"(>= {args.min_slide_cells} cells in both groups) ===")
    print(f"  overall  {args.letter}: {df.loc[df.is_letter,'score'].mean():.3f}   "
          f"other: {df.loc[~df.is_letter,'score'].mean():.3f}")
    print(f"  slides where {args.letter} scores ABOVE the rest of its own slide: "
          f"{100*pos:.1f}%  ({int((usable['gap']>0).sum())}/{len(usable)})")
    print(f"  median within-slide gap: {usable['gap'].median():+.3f}  "
          f"(10th-90th {usable['gap'].quantile(.1):+.3f} to {usable['gap'].quantile(.9):+.3f})")
    print(f"  corr(slide score in {args.letter}, slide score in others): r = {r:+.2f}")
    print(f"  variance between slides: {between:.4f}   within-slide gap: {within:.4f}")
    print("\n  READ: a gap positive on nearly every slide, with a modest r, says the programme\n"
          "  travels with the cells = a CELL STATE. A high r with the gap collapsing toward\n"
          "  zero, and a few slides carrying the signal, says the SECTION is elevated = handling.")

    print(f"\n  Top 8 slides by {args.letter} score:")
    top = usable.sort_values("score_letter", ascending=False).head(8)
    print(f"    {'slide':34s} {'score_' + args.letter:>10s} {'other':>8s} {'gap':>8s} "
          f"{'%' + args.letter:>7s} {'n':>8s}")
    for s, row in top.iterrows():
        print(f"    {str(s)[:32]:34s} {row['score_letter']:>10.3f} {row['score_other']:>8.3f} "
              f"{row['gap']:>+8.3f} {row['letter_pct']:>7.1f} {row['n_letter']:>8,}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    tbl.sort_values("score_letter", ascending=False).to_csv(args.output_csv)
    print(f"\nWrote {args.output_csv} ({len(tbl)} slides)")


if __name__ == "__main__":
    main()
