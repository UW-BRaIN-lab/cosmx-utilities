#!/usr/bin/env python3
"""Do a de-novo letter's cells and the natively-called cells of its destination come from the
SAME donors, or different ones?

THE QUESTION. In the semi-supervised fit every cell chose from the same 81 profiles — 54 named
GBmap types plus the 27 de-novo clusters. Some cells picked de-novo `l`; others picked
`Endo_arterial` outright. Same menu, same competition. So what is different about the cells that
needed a de-novo bin?

    A   semisup == <letter>  AND forced == D     the letter's cells that GBmap would call D
    B   semisup == D                             cells that got D without needing a de-novo bin

Comparing A and B on markers is what 75d does, and its own header records the catch: a cell
landed in the letter precisely BECAUSE it fit the letter's profile better, so "A and B differ" is
definitional, not a finding. This script asks a question whose answer is not built in.

WHY DONOR IS THE DISCRIMINATING AXIS. The 54 named profiles are FIXED — rescaled reference
columns that cannot move. The 27 de-novo profiles are FITTED, estimated by EM from this cohort's
own cells. A fitted profile goes to where the cells are; a fixed one cannot. So a de-novo cluster
wins wherever there is a systematic offset between this cohort's cells and the reference —
platform, batch, patient, fixation, segmentation — and the letter would then be the part of that
offset large enough to beat a fixed profile, not a cell state at all.

That reading makes a prediction the marker comparison cannot: if the letter is absorbing a
cohort or patient offset, A and B SEPARATE BY DONOR — some donors send their D-like cells to the
letter, others to D. If instead the letter is a real cell state, every donor should contribute
both, in similar proportion.

THE NULL. "Overdispersed across donors" needs a yardstick, because cell-type composition varies
between patients for ordinary biological reasons. The yardstick used here is the cohort's own:
every pair of NAMED types inside D's compartment. Both members of such a pair are fixed
profiles, so their donor-to-donor ratio carries real composition variation and sampling noise but
NOT the fitted-vs-fixed asymmetry. If A-vs-B is overdispersed far beyond that null, the letter is
tracking donors in a way ordinary composition does not explain.

Dispersion is Pearson's phi: mean squared standardised residual of each group's share, where
phi ~ 1 means donors are homogeneous and phi >> 1 means the split is donor-structured.

Grouping by `export` instead of `case` turns the same test into a batch test — Export_source
names which of the five flat-file exports a slide came from.

Usage:
    uv run python pipeline/python/denovo_vs_native_by_group.py \\
        --typing-h5 anchor_typing.h5 --forced-csv forced_named_posteriors.csv \\
        --letter l --destinations Endo_arterial,Endo_capilar,Pericyte \\
        --group-by case --output-dir l_by_donor
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from anchor_profiles import read_cell_calls
from fov_annotations import DEFAULT_ANNOTATIONS, DEFAULT_CROSSWALK, annotate_cells

_REFERENCE = Path(__file__).resolve().parents[1] / "reference"
DEFAULT_COMPARTMENTS = _REFERENCE / "gbmap_compartments.csv"
GROUP_COLUMNS = {"case": "Case", "slide": "slide", "region": "Region",
                 "block": "Block", "export": "Export_source"}
MIN_CELLS_PER_GROUP = 25
MIN_GROUPS = 3
# A group whose share sits outside this band is "polarised" — it sent essentially all of its
# D-like cells one way. A handful of those is the readable signature of a donor-driven split.
POLARISED = 0.10


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="anchor_typing.h5 — the semi-supervised label per cell.")
    p.add_argument("--forced-csv", type=Path,
                   help="75c's forced_named_posteriors.csv (cell_id, top1_type). Without it the "
                        "letter is compared as a whole rather than split by destination.")
    p.add_argument("--letter", required=True, help="De-novo letter to dissect, e.g. l")
    p.add_argument("--destinations", default="",
                   help="Comma-separated GBmap types to compare against. Default: the letter's "
                        "own largest forced destinations (needs --forced-csv).")
    p.add_argument("--top-dest", type=int, default=3,
                   help="How many destinations to take when --destinations is omitted.")
    p.add_argument("--group-by", choices=sorted(GROUP_COLUMNS), default="case",
                   help="case (default), slide, region, block, or export for a batch test.")
    p.add_argument("--compartments", type=Path, default=DEFAULT_COMPARTMENTS,
                   help="gbmap_compartments.csv, for the named-pair null.")
    p.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    p.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    p.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_GROUP,
                   help=f"Groups with fewer than this many cells across A+B are dropped "
                        f"(default {MIN_CELLS_PER_GROUP}).")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def dispersion(table: pd.DataFrame, min_cells: int) -> dict:
    """Pearson dispersion of the A-share across groups, plus the readable descriptives.

    `table` has one row per group with columns n_a and n_b. phi ~ 1 means every group splits in
    the same proportion (homogeneous); phi >> 1 means the split tracks the grouping.
    """
    t = table[(table["n_a"] + table["n_b"]) >= min_cells].copy()
    if len(t) < MIN_GROUPS:
        return {"groups": len(t), "phi": float("nan"), "share": float("nan"),
                "polarised": 0, "lo": float("nan"), "hi": float("nan")}
    t["n"] = t["n_a"] + t["n_b"]
    share = t["n_a"].sum() / t["n"].sum()
    if share in (0.0, 1.0):
        return {"groups": len(t), "phi": float("nan"), "share": share,
                "polarised": 0, "lo": float("nan"), "hi": float("nan")}
    t["share"] = t["n_a"] / t["n"]
    # Pearson chi-square over (groups - 1) degrees of freedom: the standard dispersion estimate.
    chi2 = (((t["share"] - share) ** 2) * t["n"] / (share * (1 - share))).sum()
    polarised = int(((t["share"] <= POLARISED) | (t["share"] >= 1 - POLARISED)).sum())
    return {"groups": len(t), "phi": chi2 / (len(t) - 1), "share": share,
            "polarised": polarised, "lo": t["share"].min(), "hi": t["share"].max()}


def counts_by_group(cells: pd.DataFrame, group_col: str,
                    mask_a: pd.Series, mask_b: pd.Series) -> pd.DataFrame:
    """One row per group: how many cells fell in A and in B."""
    a = cells.loc[mask_a].groupby(group_col, observed=True).size().rename("n_a")
    b = cells.loc[mask_b].groupby(group_col, observed=True).size().rename("n_b")
    return pd.concat([a, b], axis=1).fillna(0).astype(int)


def named_pair_null(cells: pd.DataFrame, group_col: str, siblings: list[str],
                    min_cells: int) -> pd.DataFrame:
    """Dispersion for every pair of NAMED types in the compartment — the matched yardstick."""
    rows = []
    for x, y in itertools.combinations(siblings, 2):
        tbl = counts_by_group(cells, group_col,
                              cells["cell_type"] == x, cells["cell_type"] == y)
        d = dispersion(tbl, min_cells)
        if not np.isnan(d["phi"]):
            rows.append({"pair": f"{x} vs {y}", **d})
    return pd.DataFrame(rows).sort_values("phi") if rows else pd.DataFrame()


def main() -> None:
    args = parse_args()
    group_col = GROUP_COLUMNS[args.group_by]
    extra = ("Export_source",) if args.group_by == "export" else ()

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]].set_index("cell_id")
    print(f"{len(calls):,} cells, {calls['cell_type'].nunique()} types")

    ann = annotate_cells(pd.Index(calls.index), annotations=args.annotations,
                         crosswalk=args.crosswalk, extra_columns=extra)
    cells = calls.join(ann, how="inner")
    if group_col not in cells.columns:
        sys.exit(f"ERROR: grouping column {group_col} not available.")
    print(f"{len(cells):,} annotated cells across {cells[group_col].nunique()} "
          f"{args.group_by}(s)")

    if args.letter not in set(cells["cell_type"]):
        sys.exit(f"ERROR: no cells labelled {args.letter!r}. Present de-novo labels: "
                 f"{sorted(t for t in cells['cell_type'].unique() if len(t) == 1)}")

    forced = None
    if args.forced_csv:
        f = pd.read_csv(args.forced_csv)
        col = "top1_type" if "top1_type" in f.columns else f.columns[1]
        forced = f.set_index("cell_id")[col]
        cells = cells.join(forced.rename("forced"))
        print(f"forced call joined for {int(cells['forced'].notna().sum()):,} cells")

    is_letter = cells["cell_type"] == args.letter
    if args.destinations.strip():
        destinations = [d.strip() for d in args.destinations.split(",") if d.strip()]
    elif forced is not None:
        destinations = (cells.loc[is_letter, "forced"].value_counts()
                             .head(args.top_dest).index.tolist())
        print(f"destinations (the letter's own largest): {', '.join(destinations)}")
    else:
        sys.exit("ERROR: pass --destinations, or --forced-csv to take them automatically.")

    comp = pd.read_csv(args.compartments).set_index("gbmap_type")["compartment"]
    present = set(cells["cell_type"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries, per_group = [], []
    for dest in destinations:
        if dest not in present:
            print(f"  (skipped {dest}: no cells natively called it)")
            continue
        mask_a = is_letter & (cells["forced"] == dest) if forced is not None else is_letter
        mask_b = cells["cell_type"] == dest
        tbl = counts_by_group(cells, group_col, mask_a, mask_b)
        d = dispersion(tbl, args.min_cells)
        label = f"{args.letter}->{dest}" if forced is not None else args.letter

        siblings = [t for t in comp[comp == comp.get(dest)].index
                    if t in present and t != dest]
        null = named_pair_null(cells, group_col, [dest] + siblings, args.min_cells)
        pct = (float((null["phi"] < d["phi"]).mean() * 100)
               if not null.empty and not np.isnan(d["phi"]) else float("nan"))

        summaries.append({"comparison": f"{label} vs {dest} [native]", "destination": dest,
                          "n_a": int(tbl["n_a"].sum()), "n_b": int(tbl["n_b"].sum()),
                          **d, "null_pairs": len(null),
                          "null_median_phi": float(null["phi"].median()) if not null.empty
                          else float("nan"), "percentile_vs_null": pct})
        tbl = tbl.assign(comparison=label,
                         share_a=(tbl["n_a"] / (tbl["n_a"] + tbl["n_b"])).round(4))
        per_group.append(tbl)

        print(f"\n=== {label}  vs  {dest} [native] — by {args.group_by} ===")
        print(f"  cells: A {int(tbl['n_a'].sum()):,}   B {int(tbl['n_b'].sum()):,}")
        print(f"  usable {args.group_by}s (>= {args.min_cells} cells): {d['groups']}")
        print(f"  overall share going to the letter: {100*d['share']:.1f}%")
        print(f"  per-{args.group_by} share: {100*d['lo']:.1f}% to {100*d['hi']:.1f}%")
        print(f"  polarised {args.group_by}s (<= {100*POLARISED:.0f}% or >= "
              f"{100*(1-POLARISED):.0f}%): {d['polarised']}/{d['groups']}")
        print(f"  dispersion phi = {d['phi']:.1f}   (1.0 = donors split identically)")
        if not null.empty:
            print(f"  matched null — {len(null)} named pairs in the {comp.get(dest)} "
                  f"compartment: median phi {null['phi'].median():.1f}, "
                  f"range {null['phi'].min():.1f}-{null['phi'].max():.1f}")
            print(f"  => this split is more donor-structured than {pct:.0f}% of them")

    if not summaries:
        sys.exit("ERROR: no destination produced a comparison.")
    pd.DataFrame(summaries).to_csv(args.output_dir / "dispersion_summary.csv", index=False)
    pd.concat(per_group).rename_axis(args.group_by).to_csv(
        args.output_dir / f"counts_by_{args.group_by}.csv")
    print(f"\nWrote {args.output_dir}/dispersion_summary.csv and "
          f"counts_by_{args.group_by}.csv")
    print("\nREAD: phi near the null's range = the letter is NOT donor-structured beyond ordinary\n"
          "composition, so a cell-state reading survives. phi far above it, with polarised\n"
          "donors, = the letter is absorbing a donor/batch offset the fixed profiles cannot move to.")


if __name__ == "__main__":
    main()
