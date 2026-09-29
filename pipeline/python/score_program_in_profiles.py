#!/usr/bin/env python3
"""Score a gene programme across every column of a reference profile matrix.

This answers a question that needs no cohort and no cluster time: does the REFERENCE itself
stratify on the programme? If it does, then a diagnostic that compares a de-novo letter against
"every other cell" is measuring the units' cell-type composition as much as the letter, and the
comparison group has to be chosen explicitly (diagnose_program_by_section.py --compare-to).

That is not hypothetical. On gbmap_level4_panel.csv with the amplicon set
EGFR,CDK4,MDM2,NUP107,OS9,SLC35E3, all nine malignant Core-L4 columns take ranks 1-9 of 54
(mean 1.93 against 0.72 for the other 45), so an all-other-cells pool — mostly non-malignant —
sits about 1.2 log units below any malignant column and dilutes a malignant-vs-malignant
difference to nothing.

Columns are put on a common scale before scoring (each is expected counts, so the column sums
differ), and the score is the mean log1p CP10K over the programme genes.

Usage:
    uv run python pipeline/python/score_program_in_profiles.py \\
        --profiles pipeline/reference/gbmap_level4_panel.csv \\
        --genes EGFR,CDK4,MDM2,NUP107,OS9,SLC35E3 --highlight malignant
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pseudobulk_core import DEFAULT_SCALE_FACTOR

# The nine malignant Core-L4 columns. RG and Stress_sig are excluded: both are borderline and
# neither is a Neftel state. Kept in step with diagnose_program_by_section.MALIGNANT_CORE_L4.
MALIGNANT_CORE_L4 = ("AC-like", "AC-like_Prolif", "MES-like_hypoxia_independent",
                     "MES-like_hypoxia_MHC", "NPC-like_OPC", "NPC-like_Prolif",
                     "NPC-like_neural", "OPC-like", "OPC-like_Prolif")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles", type=Path, required=True,
                   help="Reference matrix, genes x types, first column the gene name.")
    p.add_argument("--genes", required=True, help="Comma-separated programme genes.")
    p.add_argument("--highlight", default="",
                   help="Comma-separated types to flag and summarise as a group, or 'malignant' "
                        "for the nine malignant Core-L4 columns.")
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    p.add_argument("--output-csv", type=Path)
    return p.parse_args()


def score_profiles(profiles: pd.DataFrame, program: list[str], scale_factor: float) -> pd.Series:
    """Mean log1p normalised expression of `program` in each column of `profiles`."""
    missing = [g for g in program if g not in profiles.index]
    if missing:
        sys.exit(f"ERROR: programme gene(s) not in the profile matrix: {missing}")
    normed = np.log1p(profiles / profiles.sum(axis=0) * scale_factor)
    return normed.loc[program].mean(axis=0)


def main() -> None:
    args = parse_args()
    program = [g.strip() for g in args.genes.split(",") if g.strip()]
    if not program:
        sys.exit("ERROR: --genes is empty.")

    profiles = pd.read_csv(args.profiles, index_col=0)
    print(f"{args.profiles.name}: {profiles.shape[0]:,} genes x {profiles.shape[1]} types")
    print(f"Scoring the {len(program)}-gene programme: {', '.join(program)}")

    score = score_profiles(profiles, program, args.scale_factor)

    if args.highlight.strip().lower() == "malignant":
        flagged = [c for c in MALIGNANT_CORE_L4 if c in score.index]
    else:
        flagged = [c.strip() for c in args.highlight.split(",")
                   if c.strip() and c.strip() in score.index]

    out = pd.DataFrame({"score": score.round(3),
                        "rank": score.rank(ascending=False).astype(int),
                        "highlighted": [c in flagged for c in score.index]})
    out = out.sort_values("score", ascending=False)
    print(f"\n=== {len(out)} types ranked by the programme score ===")
    print(out.to_string())

    if flagged:
        rest = [c for c in score.index if c not in flagged]
        worst = int(out.loc[flagged, "rank"].max())
        print(f"\nhighlighted ({len(flagged)}): mean {score[flagged].mean():.3f}, "
              f"worst rank {worst} of {len(out)}")
        print(f"others     ({len(rest)}): mean {score[rest].mean():.3f}")
        if worst <= len(flagged):
            print("  => the highlighted types SWEEP the top of the ranking, so this programme "
                  "separates\n     them from everything else. A comparison group of 'every other "
                  "cell' is dominated by\n     the others and will dilute the contrast — choose "
                  "the comparator explicitly.")

    print("\n=== per gene, highlighted columns ===")
    normed = np.log1p(profiles / profiles.sum(axis=0) * args.scale_factor)
    print(normed.loc[program, flagged if flagged else list(score.index)[:10]].round(2).to_string())

    if args.output_csv:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output_csv)
        print(f"\nWrote {args.output_csv}")


if __name__ == "__main__":
    main()
