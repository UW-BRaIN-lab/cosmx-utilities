#!/usr/bin/env python3
"""Is a letter's programme deficit cohort-wide, or carried by one donor, slide or region?

The 75e per-type summary pools over FOVs, so a verdict can rest on a handful of FOVs from a
single slide without saying so. In the amplicon runs the usable FOVs for b, t and e were all
topped by 7495G37302G3 — which was also the amplicon-high outlier — so "the letters lack the
amplicon" has to be separated from "one amplicon-high patient dominates the comparison".

Takes the LONG per-unit x comparator tables that 75e writes with --per-unit-csv, subtracts the
matched-control run from the amplicon run FOV by FOV, and reports the excess per group plus a
leave-one-group-out sweep. Excess is the quantity that matters: the plain amplicon gap is
inflated by the letters simply being dim, and only the control subtraction removes that.

DONOR, NOT SLIDE. A CosMx slide carries two tissue pieces from two different cases, and this
cohort's FOV numbering is continuous 1-200 across both, so the pieces cannot be separated from
the cell ids alone. They CAN be separated from the AtoMx annotation reference, which is one row
per (slide, FOV) with Case, Block and Region — pass --annotations and --crosswalk and group by
case or region. The crosswalk is needed because the annotation file keys on canonical slide
names while the flat files (and therefore our cell ids) use AtoMx export folder names.

Region matters as much as donor here: the amplicon question is about tumour genetics, so a
deficit that only appears in contralateral uninvolved tissue would mean something quite
different from one that holds in tumour bulk.

Usage:
    uv run python pipeline/python/program_gap_by_group.py \\
        --amplicon b_amplicon_program_by_unit.csv \\
        --control  b_control_program_by_unit.csv \\
        --comparator OPC-like --group-by case
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from fov_annotations import (DEFAULT_ANNOTATIONS, DEFAULT_CROSSWALK,
                             load_fov_annotations)

MIN_UNITS_PER_GROUP = 3
# A sign test alone over-claims: a group at -0.03 counts as "negative" while being substantively
# null. Groups are also counted against this magnitude, so a near-zero group reads as null.
MATERIAL_EXCESS = 0.1
GROUP_COLUMNS = {"case": "Case", "region": "Region", "block": "Block", "slide": "slide"}


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
    p.add_argument("--group-by", choices=sorted(GROUP_COLUMNS), default="case",
                   help="Grouping for the breakdown and the leave-one-out (default case). "
                        "Anything but 'slide' needs --annotations and --crosswalk.")
    p.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS,
                   help=f"AtoMx annotation reference: one row per (Slide, FOV) with Case, Block, "
                        f"Region (default the committed {DEFAULT_ANNOTATIONS.name}).")
    p.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK,
                   help=f"Slide-name crosswalk mapping the flat-file folder name in our cell ids "
                        f"to the annotation file's canonical slide name (default the committed "
                        f"{DEFAULT_CROSSWALK.name}).")
    p.add_argument("--regions", default="",
                   help="Comma-separated Region values to keep (e.g. 'Tumor bulk,Infiltrating "
                        "edge'). Default keeps all. Needs --annotations.")
    p.add_argument("--material", type=float, default=MATERIAL_EXCESS,
                   help=f"Magnitude below which a group's median excess is reported as null "
                        f"rather than as supporting (default {MATERIAL_EXCESS}). A bare sign test "
                        f"counts a group at -0.03 as negative, which over-claims.")
    p.add_argument("--min-units", type=int, default=MIN_UNITS_PER_GROUP,
                   help=f"Groups with fewer matched FOVs are listed but not ranked "
                        f"(default {MIN_UNITS_PER_GROUP}).")
    p.add_argument("--output-csv", type=Path)
    return p.parse_args()


def load_per_unit(path: Path, comparator: str, label: str) -> pd.DataFrame:
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


def load_annotations(annotations: Path, crosswalk: Path) -> pd.DataFrame:
    """Per-FOV Case/Block/Region. `slide` is dropped: the per-unit tables already carry it."""
    return load_fov_annotations(annotations, crosswalk).drop(columns="slide")


def main() -> None:
    args = parse_args()
    group_col = GROUP_COLUMNS[args.group_by]
    if args.group_by != "slide":
        missing = [str(f) for f in (args.annotations, args.crosswalk) if not f.is_file()]
        if missing:
            sys.exit(f"ERROR: --group-by {args.group_by} needs --annotations and --crosswalk; "
                     f"not found: {', '.join(missing)}")

    amp = load_per_unit(args.amplicon, args.comparator, "gap_amplicon")
    ctl = load_per_unit(args.control, args.comparator, "gap_control")

    # Inner join on the FOV: only units measured in BOTH runs support a matched subtraction.
    both = amp.join(ctl.drop(columns="slide"), how="inner")
    dropped = len(amp) - len(both)
    if both.empty:
        sys.exit("ERROR: the two runs share no unit for this comparator.")
    both["excess"] = both["gap_amplicon"] - both["gap_control"]
    if dropped:
        print(f"NOTE: {dropped} unit(s) present in the amplicon run only, dropped from the "
              f"matched comparison.")

    if args.annotations.is_file() and args.crosswalk.is_file():
        ann = load_annotations(args.annotations, args.crosswalk)
        before = len(both)
        joined = both.join(ann, how="left")
        unannotated = int(joined["Case"].isna().sum())
        if unannotated == before and args.group_by == "slide":
            # A cohort with no matching annotation is fine when grouping by slide — that is the
            # one grouping the cell ids already support on their own.
            print(f"NOTE: none of the {before} FOVs is in the annotation reference; continuing "
                  f"without it, which --group-by slide does not need.")
        else:
            both = joined
            if unannotated:
                print(f"WARNING: {unannotated}/{before} matched FOVs have no annotation and are "
                      f"dropped from the grouped view.")
                both = both[both["Case"].notna()]
        if args.regions.strip():
            keep = [r.strip() for r in args.regions.split(",") if r.strip()]
            unknown = set(keep) - set(ann["Region"].unique())
            if unknown:
                sys.exit(f"ERROR: unknown region(s) {sorted(unknown)}. Present: "
                         f"{', '.join(sorted(ann['Region'].unique()))}")
            both = both[both["Region"].isin(keep)]
            print(f"Kept {len(both)} FOVs in region(s): {', '.join(keep)}")
        if both.empty:
            sys.exit("ERROR: no annotated FOVs left after filtering.")

    print(f"\n=== {args.comparator}: {len(both)} matched FOVs, "
          f"{both[group_col].nunique()} {args.group_by}(s) ===")
    print(f"  overall median excess: {both['excess'].median():+.3f}   "
          f"(amplicon {both['gap_amplicon'].median():+.3f}, "
          f"control {both['gap_control'].median():+.3f})")
    print(f"  FOVs with a NEGATIVE excess: {100 * (both['excess'] < 0).mean():.1f}%")

    if "Region" in both.columns and args.group_by != "region":
        print("\n  by region (the amplicon question is about tumour, so this must not be "
              "contralateral-only):")
        for region, g in both.groupby("Region"):
            print(f"    {str(region)[:30]:32s} {len(g):>5d} FOVs   median excess "
                  f"{g['excess'].median():+.3f}")

    per_group = (both.groupby(group_col)
                     .agg(units=("excess", "size"), median_excess=("excess", "median"),
                          median_amplicon=("gap_amplicon", "median"),
                          median_control=("gap_control", "median"))
                     .sort_values("median_excess"))
    print(f"\n  === per {args.group_by} (>= {args.min_units} FOVs ranked first) ===")
    print(f"    {args.group_by:28s} {'units':>6s} {'excess':>8s} {'amplicon':>9s} {'control':>8s}")
    ranked = per_group[per_group["units"] >= args.min_units]
    thin = per_group[per_group["units"] < args.min_units]
    for name, row in pd.concat([ranked, thin]).iterrows():
        flag = "" if row["units"] >= args.min_units else "   (thin)"
        print(f"    {str(name)[:26]:28s} {int(row['units']):>6d} {row['median_excess']:>+8.3f} "
              f"{row['median_amplicon']:>+9.3f} {row['median_control']:>+8.3f}{flag}")

    if len(ranked) >= 2:
        n_neg = int((ranked["median_excess"] <= -args.material).sum())
        n_pos = int((ranked["median_excess"] >= args.material).sum())
        n_null = len(ranked) - n_neg - n_pos
        print(f"\n  {args.group_by}s with >= {args.min_units} FOVs "
              f"(|excess| >= {args.material} to count either way):")
        print(f"    materially NEGATIVE {n_neg}/{len(ranked)}   "
              f"null {n_null}/{len(ranked)}   materially POSITIVE {n_pos}/{len(ranked)}")
        if n_null:
            nulls = ranked.index[ranked["median_excess"].abs() < args.material].tolist()
            print(f"    null {args.group_by}(s): "
                  f"{', '.join(str(x) for x in nulls)} — these do NOT support the effect")

    # With a single group there is nothing to leave out, and dropping it empties the frame.
    if len(per_group) < 2:
        print(f"\n  === leave one {args.group_by} out: SKIPPED ===")
        print(f"  Only one {args.group_by} is present, so there is no deletion to make and no "
              f"robustness claim\n  to be had from this run. Widen --regions, or group by "
              f"something with more levels.")
        if args.output_csv:
            args.output_csv.parent.mkdir(parents=True, exist_ok=True)
            per_group.to_csv(args.output_csv)
            print(f"\nWrote {args.output_csv}")
        return

    print(f"\n  === leave one {args.group_by} out ===")
    print(f"  If the effect is carried by one {args.group_by}, dropping it collapses the excess.")
    print(f"    {'dropped':28s} {'units left':>11s} {'median excess':>14s} {'change':>9s}")
    base = both["excess"].median()
    loo = []
    for name in per_group.index:
        rest = both[both[group_col] != name]
        if rest.empty:
            continue
        med = rest["excess"].median()
        loo.append({"group": name, "n_left": len(rest), "median_excess": med,
                    "change": med - base})
    loo_df = pd.DataFrame(loo).sort_values("change", ascending=False)
    for _, row in loo_df.iterrows():
        print(f"    {str(row['group'])[:26]:28s} {int(row['n_left']):>11d} "
              f"{row['median_excess']:>+14.3f} {row['change']:>+9.3f}")
    worst = loo_df.iloc[0]
    print(f"\n  Worst case: dropping {worst['group']} moves the median excess to "
          f"{worst['median_excess']:+.3f} (from {base:+.3f}).")
    print(f"  A verdict that survives every single-{args.group_by} deletion is not one "
          f"{args.group_by}'s artefact.")

    if args.output_csv:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        per_group.to_csv(args.output_csv)
        print(f"\nWrote {args.output_csv}")


if __name__ == "__main__":
    main()
