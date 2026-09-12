#!/usr/bin/env python3
"""Is a letter's programme deficit carried by the whole cohort, or by one slide?

The 75e per-type summary pools over units, so a verdict can rest on a handful of FOVs from a
single slide without saying so. In the amplicon runs the usable FOVs for b, t and e were all
topped by 7495G37302G3 — which was also the amplicon-high outlier — so "the letters lack the
amplicon" had to be separated from "one amplicon-high patient dominates the comparison".

Takes the LONG per-unit x comparator tables that 75e writes with --per-unit-csv, subtracts the
matched-control run from the amplicon run FOV by FOV, and reports the excess per slide plus a
leave-one-slide-out sweep. Excess is the quantity that matters: the plain amplicon gap is
inflated by the letters simply being dim, and only the control subtraction removes that.

NOTE ON UNITS: a CosMx slide carries two tissue pieces that need not share a donor, and this
cohort's FOV numbering is continuous across both, so FOVs cannot be assigned to one of them.
Slide is therefore the finest honest grouping here — report it as slide, not donor.

Usage:
    uv run python pipeline/python/program_gap_by_slide.py \\
        --amplicon b_amplicon_program_by_unit.csv \\
        --control  b_control_program_by_unit.csv \\
        --comparator OPC-like
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MIN_UNITS_PER_SLIDE = 3


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--amplicon", type=Path, required=True,
                   help="Per-unit table from the run scoring the programme under test.")
    p.add_argument("--control", type=Path, required=True,
                   help="Per-unit table from the matched-control run (same letter, same "
                        "--compare-to, control gene set).")
    p.add_argument("--comparator", default="OPC-like",
                   help="Which comparator row to analyse (default OPC-like).")
    p.add_argument("--min-units", type=int, default=MIN_UNITS_PER_SLIDE,
                   help=f"Slides with fewer matched FOVs are listed but not ranked "
                        f"(default {MIN_UNITS_PER_SLIDE}).")
    p.add_argument("--output-csv", type=Path)
    return p.parse_args()


def load(path: Path, comparator: str, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"unit", "slide", "comparator", "gap"}
    if not need.issubset(df.columns):
        sys.exit(f"ERROR: {path} is missing {sorted(need - set(df.columns))}. It must be a "
                 f"--per-unit-csv table from 75e.")
    sub = df[df["comparator"] == comparator]
    if sub.empty:
        sys.exit(f"ERROR: no rows for comparator {comparator!r} in {path}. Present: "
                 f"{', '.join(sorted(df['comparator'].unique()))}")
    return sub.set_index("unit")[["slide", "gap"]].rename(columns={"gap": label})


def main() -> None:
    args = parse_args()
    amp = load(args.amplicon, args.comparator, "gap_amplicon")
    ctl = load(args.control, args.comparator, "gap_control")

    # Inner join on the FOV: only units measured in BOTH runs support a matched subtraction.
    both = amp.join(ctl.drop(columns="slide"), how="inner")
    dropped = len(amp) - len(both)
    if both.empty:
        sys.exit("ERROR: the two runs share no unit for this comparator.")
    both["excess"] = both["gap_amplicon"] - both["gap_control"]
    if dropped:
        print(f"NOTE: {dropped} unit(s) present in the amplicon run only, dropped from the "
              f"matched comparison.")

    print(f"=== {args.comparator}: {len(both)} matched FOVs across "
          f"{both['slide'].nunique()} slides ===")
    print(f"  overall median excess: {both['excess'].median():+.3f}   "
          f"(amplicon {both['gap_amplicon'].median():+.3f}, "
          f"control {both['gap_control'].median():+.3f})")
    print(f"  FOVs with a NEGATIVE excess: {100 * (both['excess'] < 0).mean():.1f}%")

    per_slide = (both.groupby("slide")
                     .agg(units=("excess", "size"), median_excess=("excess", "median"),
                          median_amplicon=("gap_amplicon", "median"),
                          median_control=("gap_control", "median"))
                     .sort_values("median_excess"))
    print(f"\n  === per slide (>= {args.min_units} FOVs ranked first) ===")
    print(f"    {'slide':28s} {'units':>6s} {'excess':>8s} {'amplicon':>9s} {'control':>8s}")
    ranked = per_slide[per_slide["units"] >= args.min_units]
    thin = per_slide[per_slide["units"] < args.min_units]
    for name, row in pd.concat([ranked, thin]).iterrows():
        flag = "" if row["units"] >= args.min_units else "   (thin)"
        print(f"    {str(name)[:26]:28s} {int(row['units']):>6d} {row['median_excess']:>+8.3f} "
              f"{row['median_amplicon']:>+9.3f} {row['median_control']:>+8.3f}{flag}")

    if len(ranked) >= 2:
        share = 100 * (ranked["median_excess"] < 0).mean()
        print(f"\n  slides with >= {args.min_units} FOVs and a negative median excess: "
              f"{share:.0f}% ({int((ranked['median_excess'] < 0).sum())}/{len(ranked)})")

    print("\n  === leave one slide out ===")
    print("  If the effect is carried by one slide, dropping it collapses the overall excess.")
    print(f"    {'dropped slide':28s} {'units left':>11s} {'median excess':>14s} {'change':>9s}")
    base = both["excess"].median()
    loo = []
    for slide in per_slide.index:
        rest = both[both["slide"] != slide]
        if rest.empty:
            continue
        med = rest["excess"].median()
        loo.append({"slide": slide, "n_left": len(rest), "median_excess": med,
                    "change": med - base})
    loo_df = pd.DataFrame(loo).sort_values("change", ascending=False)
    for _, row in loo_df.iterrows():
        print(f"    {str(row['slide'])[:26]:28s} {int(row['n_left']):>11d} "
              f"{row['median_excess']:>+14.3f} {row['change']:>+9.3f}")
    worst = loo_df.iloc[0]
    print(f"\n  Worst case: dropping {worst['slide']} moves the median excess to "
          f"{worst['median_excess']:+.3f} (from {base:+.3f}).")
    print("  A verdict that survives every single-slide deletion is not one slide's artefact.")

    if args.output_csv:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        per_slide.to_csv(args.output_csv)
        print(f"\nWrote {args.output_csv}")


if __name__ == "__main__":
    main()
