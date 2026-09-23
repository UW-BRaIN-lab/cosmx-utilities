#!/usr/bin/env python3
"""If we drop a de-novo letter from the hierarchy, where do its cells actually go?

THE QUESTION BEHIND THE DROP DECISIONS. Dropping a letter means not giving it a leaf and letting
GBmap's named profiles absorb its cells. The triage justified each drop by showing the letter's
markers match some GBmap type -- but a matching leaf is necessary, not sufficient. It does not
show the FIXED profile can capture those cells, and the sharded InSituTree run is the experiment
for exactly that: typed against fixed profiles only, it sent 52.9% of the cohort to Low_signal.

So this joins the two runs on cell id and reads the answer off directly:

    anchor fit          semi-supervised, 27 fitted de-novo letters + 54 fixed named profiles
    InSituTree run      the SAME cells typed against FIXED profiles only, Low_signal available

For every anchor label, what did the fixed-profile run call those same cells? A dropped letter
whose cells land on its matching leaf was a safe drop. One whose cells land in Low_signal was
not: the leaf existed and still could not hold them, and that letter needs its profile kept.

This needs no rebuild and no counts matrix -- both inputs are label tables that already exist.

THE ONE THING THAT CAN SINK IT is the cell-id join. The two runs must share the stage-1
"<slide>_F<fov>_C<cell>" ids (per-cell InSituCNV, for instance, does NOT join to the full cohort
because its id format differs). The overlap is therefore reported and checked before anything is
computed, rather than assumed.

Inputs:
  --anchor-h5   stage4_anchor_pruned/anchor/anchor_typing.h5 (cell_id, cell_type = letters + named)
  --typed-h5ad  stage4_insitutree/cosmx_typed.h5ad (obs.index = cell_id, obs[cell_type] fixed-profile)
Outputs (--output-dir):
  fate_by_anchor_label.csv   per anchor label: n, % Low_signal, top fixed-profile destinations
  fate_matrix.csv            full anchor label x fixed-profile label contingency
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from anchor_profiles import read_cell_calls

LOW_LABEL = "Low_signal"
MIN_OVERLAP = 0.05          # below this the join is broken, not merely partial
DENOVO_MAX_LEN = 2          # InSituType's cluster_name_pool: 1-2 lowercase letters


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--anchor-h5", type=Path, required=True)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--celltype-key", default="cell_type",
                   help="obs column holding the fixed-profile label (default cell_type).")
    p.add_argument("--low-label", default=LOW_LABEL)
    p.add_argument("--top-dest", type=int, default=3)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def is_denovo(label: str) -> bool:
    return len(label) <= DENOVO_MAX_LEN and label.isalpha() and label.islower()


def main() -> None:
    args = parse_args()

    anchor = read_cell_calls(args.anchor_h5)[["cell_id", "cell_type"]]
    anchor = anchor.set_index("cell_id")["cell_type"].rename("anchor")
    print(f"anchor: {len(anchor):,} cells, {anchor.nunique()} labels")

    # Only obs is needed; backed mode keeps the 7.5M-cell matrix off the heap.
    typed = ad.read_h5ad(args.typed_h5ad, backed="r")
    if args.celltype_key not in typed.obs.columns:
        sys.exit(f"ERROR: obs has no column {args.celltype_key!r}. Present: "
                 f"{', '.join(map(str, typed.obs.columns))}")
    fixed = pd.Series(typed.obs[args.celltype_key].astype(str).to_numpy(),
                      index=typed.obs.index.to_numpy(), name="fixed")
    print(f"fixed-profile run: {len(fixed):,} cells, {fixed.nunique()} labels")
    if args.low_label not in set(fixed):
        print(f"WARNING: {args.low_label!r} does not appear in the fixed-profile labels. "
              f"Largest are: {', '.join(fixed.value_counts().head(5).index)}")

    # The join is the load-bearing assumption, so verify it before computing anything.
    common = anchor.index.intersection(fixed.index)
    frac = len(common) / max(len(anchor), 1)
    print(f"\ncell-id overlap: {len(common):,} of {len(anchor):,} anchor cells ({100*frac:.1f}%)")
    if frac < MIN_OVERLAP:
        print(f"  anchor id example: {anchor.index[0]!r}")
        print(f"  typed  id example: {fixed.index[0]!r}")
        sys.exit(f"ERROR: only {100*frac:.1f}% of anchor cells are in the fixed-profile run. The "
                 f"two runs do not share a cell-id format, so this comparison cannot be made. "
                 f"Compare the id examples above.")
    if frac < 0.9:
        print(f"  NOTE: {100*(1-frac):.1f}% of anchor cells are absent from the typed run; the "
              f"anchor is a stratified subset, so partial overlap is expected, but read the "
              f"per-label counts rather than assuming every label is equally covered.")

    df = pd.DataFrame({"anchor": anchor.reindex(common), "fixed": fixed.reindex(common)})
    matrix = pd.crosstab(df["anchor"], df["fixed"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output_dir / "fate_matrix.csv")

    baseline = 100 * (df["fixed"] == args.low_label).mean()
    print(f"\nBASELINE: {baseline:.1f}% of all joined cells are {args.low_label} under fixed "
          f"profiles.\n  Read every letter against this, not against zero.\n")

    rows = []
    for label, row in matrix.iterrows():
        n = int(row.sum())
        low = int(row.get(args.low_label, 0))
        dest = row.drop(labels=[args.low_label], errors="ignore").sort_values(ascending=False)
        top = dest.head(args.top_dest)
        rows.append({"anchor_label": label, "denovo": is_denovo(str(label)), "n": n,
                     "pct_low_signal": round(100 * low / n, 1) if n else float("nan"),
                     "top_destination": top.index[0] if len(top) else "",
                     "pct_top_destination": round(100 * top.iloc[0] / n, 1) if len(top) else 0.0,
                     "next_destinations": "; ".join(f"{k} {100*v/n:.0f}%"
                                                    for k, v in top.iloc[1:].items())})
    out = pd.DataFrame(rows).sort_values(["denovo", "pct_low_signal"], ascending=[False, False])
    out.to_csv(args.output_dir / "fate_by_anchor_label.csv", index=False)

    letters = out[out["denovo"]]
    print(f"=== DE-NOVO LETTERS under fixed profiles ({len(letters)}) ===")
    print(f"{'letter':>7s} {'n':>10s} {'% Low_signal':>13s} {'top destination':>26s} {'%':>6s}")
    for _, r in letters.iterrows():
        flag = "  <-- worse than baseline" if r.pct_low_signal > baseline else ""
        print(f"{r.anchor_label:>7s} {r.n:>10,} {r.pct_low_signal:>12.1f}% "
              f"{r.top_destination[:26]:>26s} {r.pct_top_destination:>5.1f}%{flag}")

    named = out[~out["denovo"]]
    if len(named):
        print(f"\n=== CONTROL: cells the anchor fit NAMED, under fixed profiles ===")
        print("  (these already matched a fixed profile once, so they set the ceiling)")
        print(f"  median % Low_signal: {named['pct_low_signal'].median():.1f}%   "
              f"n labels: {len(named)}")
        agree = named.apply(lambda r: r.top_destination == r.anchor_label, axis=1).mean()
        print(f"  fixed-profile run re-calls the SAME label for {100*agree:.0f}% of them")

    print(f"\nWrote {args.output_dir}/fate_by_anchor_label.csv and fate_matrix.csv")
    print("\nREAD: a dropped letter whose cells land on its matching leaf was a safe drop. One whose\n"
          "cells go to Low_signal at well above baseline was not -- the leaf existed and still could\n"
          "not hold them, so that letter needs its profile kept.")


if __name__ == "__main__":
    main()
