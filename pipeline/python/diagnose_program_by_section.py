#!/usr/bin/env python3
"""Is a de-novo letter's gene programme a cell state, or a slide-level handling artefact?

75d showed that de-novo `t` is separated from every natively-called GBmap class by a heat-shock
programme (HSPA1A/HSPA1B/HSPB1/DNAJB1/HSP90AA1/HSPH1). Heat-shock is exactly the signature warm
ischaemia and slow fixation produce, so before calling it a tumour state we have to rule out
pre-analytical handling. Same machinery works for any letter and any gene set.

SECTION, NOT SLIDE. Two tissue sections are mounted per CosMx slide and they need not even
come from the same donor, so slide_id is the WRONG unit for a handling question: it pools two
pieces with different blocks, and possibly different collection and fixation histories. Worse, if
the letter's cells sit mostly in one of the two pieces, a "within-slide" letter-vs-rest gap is
partly a between-piece comparison and would read as a cell state when it is not.

Getting that unit is the hard part. Splitting a slide's FOVs into contiguous runs (what
tissue_section_gap.py reasons about, e.g. 1-135 then 181-205) does NOT work on this cohort: run
39972319 found 56 of 57 slides with a single run, because the FOV numbering is continuous 1-200
across both pieces. So --unit fov-run is kept but is not the default, and its
sections-per-slide line is the check that it applies at all.

The default unit is instead the FOV itself. An FOV lies entirely within one piece by
construction, so it controls for the section AND for position inside it — strictly more
conservative than a section, and available straight from the stage-1 cell id
("<slide_id>_F<fov>_C<cell_ID>") with no extra join. The cost is smaller groups, so
--min-section-cells is lower by default and FOVs without enough of both groups are skipped.

When an authoritative per-Region or per-section annotation exists, pass it with --sections-csv
and it overrides everything above.

The discriminating comparison is then WITHIN a section:

  * If the programme is a CELL STATE, the letter's cells carry it wherever they occur — the
    letter-vs-other gap is positive section after section, and the two groups' section-level
    scores are only loosely coupled.
  * If it is HANDLING, the whole piece is elevated — the letter's cells and everyone else's rise
    and fall together, so the section-level scores are tightly correlated and the within-section
    gap collapses toward zero.

So this reports, per section: the programme score in the letter's cells, the score in every other
cell of that same piece, and the gap. Plus three summaries — the share of sections with a positive
gap, the correlation between the two groups' per-section scores, and how much of the total score
variance sits BETWEEN sections versus between the groups within one.

Cheap: only the programme genes and per-cell totals are needed, so the counts matrix is never
transposed or densified.

Inputs:
  --counts-h5    anchor_input.h5: /counts (CSC genes x cells), /genes, /cell_id.
  --typing-h5    anchor_typing.h5: semi-supervised /cell_type per /cell_id.
  --letter       the de-novo letter under test, e.g. t
  --sections-csv OPTIONAL cell_id -> section override, when an authoritative per-section or
                 per-Region annotation exists. Without it, sections are derived from the ids.
  --genes        comma-separated programme genes (default: the heat-shock set 75d surfaced).
Outputs:
  --output-csv   per-slide table; the summaries print to stdout.

Usage:
    uv run python pipeline/python/diagnose_program_by_slide.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --letter t --output-csv t_heatshock_by_section.csv
"""

from __future__ import annotations

import argparse
import re
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
MIN_SECTION_CELLS = 25
# A gap larger than this in a slide's FOV numbering starts a new tissue piece. The observed
# runs are far apart (1-135 then 181-205), so this is not a delicate threshold.
DEFAULT_FOV_GAP = 10
CELL_ID_RE = re.compile(r"^(?P<slide>.+)_F(?P<fov>\d+)_C\d+$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True)
    p.add_argument("--typing-h5", type=Path, required=True)
    p.add_argument("--letter", required=True)
    p.add_argument("--unit", choices=("fov", "fov-run"), default="fov",
                   help="Grouping unit derived from the cell ids. 'fov' (default) is one FOV, "
                        "which lies entirely within one tissue piece. 'fov-run' splits a slide "
                        "into contiguous FOV runs — verify its sections-per-slide line says 2 "
                        "before trusting it; on this cohort the numbering is continuous and it "
                        "collapses to the slide.")
    p.add_argument("--sections-csv", type=Path,
                   help="Optional cell_id -> section override (columns cell_id, section), for "
                        "when an authoritative per-Region annotation exists. Overrides --unit.")
    p.add_argument("--fov-gap", type=int, default=DEFAULT_FOV_GAP,
                   help=f"With --unit fov-run, the FOV numbering gap that starts a new piece "
                        f"(default {DEFAULT_FOV_GAP}).")
    p.add_argument("--genes", default=DEFAULT_PROGRAM,
                   help=f"Comma-separated programme genes (default: {DEFAULT_PROGRAM}).")
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--min-section-cells", type=int, default=MIN_SECTION_CELLS,
                   help="Skip sections with fewer than this many cells in either group.")
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


def parse_cell_ids(cell_ids: pd.Index) -> pd.DataFrame:
    """Split "<slide_id>_F<fov>_C<cell_ID>" into slide + fov, dropping ids that do not match."""
    parsed = pd.Series(cell_ids, index=cell_ids).str.extract(CELL_ID_RE)
    bad = parsed["slide"].isna()
    if bad.all():
        sys.exit("ERROR: no cell id matched '<slide_id>_F<fov>_C<cell_ID>'. Pass --sections-csv "
                 f"instead. First id: {cell_ids[0]!r}")
    if bad.any():
        print(f"WARNING: {int(bad.sum()):,} cell ids did not parse; they are dropped.")
    parsed = parsed[~bad].copy()
    parsed["fov"] = parsed["fov"].astype(int)
    return parsed


def fovs_from_cell_ids(cell_ids: pd.Index) -> pd.Series:
    """One group per FOV — the conservative unit, always inside a single tissue piece."""
    parsed = parse_cell_ids(cell_ids)
    return parsed["slide"] + ":F" + parsed["fov"].astype(str)


def sections_from_cell_ids(cell_ids: pd.Index, fov_gap: int) -> pd.Series:
    """Derive a tissue-section label per cell from "<slide_id>_F<fov>_C<cell_ID>".

    Two pieces are mounted per slide and show up as two contiguous FOV runs, so within each
    slide the sorted FOVs are split wherever the numbering jumps by more than `fov_gap`. The
    label is "<slide>:F<first>-<last>" so it stays readable in the output table.
    """
    parsed = parse_cell_ids(cell_ids)
    section = pd.Series(index=parsed.index, dtype="object")
    for slide, grp in parsed.groupby("slide", sort=False):
        fovs = np.sort(grp["fov"].unique())
        # Start a new run wherever the FOV numbering jumps.
        starts = np.concatenate([[0], np.where(np.diff(fovs) > fov_gap)[0] + 1])
        ends = np.concatenate([starts[1:] - 1, [len(fovs) - 1]])
        run_of = {}
        for lo_i, hi_i in zip(starts, ends):
            lo, hi = fovs[lo_i], fovs[hi_i]
            for f in fovs[lo_i:hi_i + 1]:
                run_of[f] = f"{slide}:F{lo}-{hi}"
        section.loc[grp.index] = grp["fov"].map(run_of)
    return section


def main() -> None:
    args = parse_args()
    program = [g.strip() for g in args.genes.split(",") if g.strip()]
    if not program:
        sys.exit("ERROR: --genes is empty.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]].set_index("cell_id")

    print(f"Scoring the {len(program)}-gene programme: {', '.join(program)}")
    score = program_score(args.counts_h5, program, args.scale_factor)

    if args.sections_csv:
        override = pd.read_csv(args.sections_csv)
        if not {"cell_id", "section"}.issubset(override.columns):
            sys.exit(f"ERROR: {args.sections_csv} needs cell_id + section columns.")
        section = override.set_index("cell_id")["section"].astype(str)
        print(f"Grouping by section from {args.sections_csv.name}: "
              f"{section.nunique()} distinct")
    elif args.unit == "fov":
        section = fovs_from_cell_ids(pd.Index(score.index))
        slides = section.str.rsplit(":", n=1).str[0]
        print(f"Grouping by FOV: {section.nunique():,} FOVs across {slides.nunique()} slides "
              f"(an FOV lies within one tissue piece, so this controls for the section)")
    else:
        section = sections_from_cell_ids(pd.Index(score.index), args.fov_gap)
        per_slide = section.drop_duplicates().str.split(":").str[0].value_counts()
        print(f"Grouping by FOV run: {section.nunique()} across {per_slide.size} slides")
        print("  runs per slide: "
              + ", ".join(f"{k} slide(s) with {v}" for v, k in per_slide.value_counts().items())
              + "   (2 is the expected mounting)")
        if (per_slide == 1).mean() > 0.5:
            print("  WARNING: most slides yielded ONE run, so this has collapsed to the slide "
                  "and does NOT separate the two mounted pieces. The FOV numbering is likely "
                  "continuous across them. Use --unit fov, or supply --sections-csv.")

    df = pd.DataFrame({"score": score})
    df["cell_type"] = calls["cell_type"].reindex(df.index)
    df["slide"] = section.reindex(df.index)
    df = df.dropna(subset=["cell_type", "slide"])
    if df.empty:
        sys.exit("ERROR: no cells left after joining labels and slides — check the ids.")
    df["is_letter"] = df["cell_type"] == args.letter
    if not df["is_letter"].any():
        sys.exit(f"ERROR: no cells labelled {args.letter!r}.")
    print(f"{len(df):,} cells, {df['slide'].nunique()} sections, "
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

    usable = tbl[(tbl["n_letter"] >= args.min_section_cells)
                 & (tbl["n_other"] >= args.min_section_cells)].copy()
    if usable.empty:
        sys.exit(f"ERROR: no group has >= {args.min_section_cells} cells of BOTH the letter and "
                 f"the rest. Lower --min-section-cells, or use a coarser --unit.")

    # How much of the letter survives the size filter. A verdict drawn from a small, biased
    # slice of the letter would not generalise, so this is reported rather than assumed.
    kept_letter = usable["n_letter"].sum() / max(tbl["n_letter"].sum(), 1)
    kept_other = usable["n_other"].sum() / max(tbl["n_other"].sum(), 1)

    pos = float((usable["gap"] > 0).mean())
    r = float(usable["score_letter"].corr(usable["score_other"]))
    # Variance sitting BETWEEN sections (both groups moving together = handling) versus the
    # within-section letter-vs-other gap (a cell state that travels with the cells).
    between = float(usable[["score_letter", "score_other"]].mean(axis=1).var(ddof=1))
    within = float((usable["gap"] / 2).pow(2).mean())

    print(f"\n=== {args.letter}: programme score, {len(usable)} usable sections "
          f"(>= {args.min_section_cells} cells in both groups) ===")
    print(f"  overall  {args.letter}: {df.loc[df.is_letter,'score'].mean():.3f}   "
          f"other: {df.loc[~df.is_letter,'score'].mean():.3f}")
    print(f"  coverage: usable groups hold {100*kept_letter:.1f}% of the {args.letter} cells "
          f"and {100*kept_other:.1f}% of the rest"
          + ("   <-- THIN, treat the verdict as provisional" if kept_letter < 0.5 else ""))
    print(f"  sections where {args.letter} scores ABOVE the rest of its own section: "
          f"{100*pos:.1f}%  ({int((usable['gap']>0).sum())}/{len(usable)})")
    print(f"  median within-section gap: {usable['gap'].median():+.3f}  "
          f"(10th-90th {usable['gap'].quantile(.1):+.3f} to {usable['gap'].quantile(.9):+.3f})")
    print(f"  corr(section score in {args.letter}, section score in others): r = {r:+.2f}")
    print(f"  variance between sections: {between:.4f}   within-section gap: {within:.4f}")
    print("\n  READ: a gap positive on nearly every section, with a modest r, says the\n"
          "  programme travels with the cells = a CELL STATE. A high r with the gap collapsing\n"
          "  toward zero, and a few pieces carrying the signal, = HANDLING.")

    print(f"\n  Top 8 sections by {args.letter} score:")
    top = usable.sort_values("score_letter", ascending=False).head(8)
    print(f"    {'section':34s} {'score_' + args.letter:>10s} {'other':>8s} {'gap':>8s} "
          f"{'%' + args.letter:>7s} {'n':>8s}")
    for s, row in top.iterrows():
        print(f"    {str(s)[:32]:34s} {row['score_letter']:>10.3f} {row['score_other']:>8.3f} "
              f"{row['gap']:>+8.3f} {row['letter_pct']:>7.1f} {row['n_letter']:>8,}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    tbl.sort_values("score_letter", ascending=False).to_csv(args.output_csv)
    print(f"\nWrote {args.output_csv} ({len(tbl)} sections)")


if __name__ == "__main__":
    main()
