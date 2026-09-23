#!/usr/bin/env python3
"""The A-vs-B "origin" split that the 75k/75l/75m diagnostics all run on, built once.

75d established the two groups, and the column names its heatmap uses:

    A   <letter>-><D>   semisup == letter AND the forced supervised run calls it D
    B   <D> [native]    semisup == D — the fit named it outright, no de-novo bin needed

75d compares them on markers, 75h on donors, 75i on the likelihood margin. The three
diagnostics added after those — spatial (75k), QC (75l), co-embedding (75m) — need the same
split and nothing else from 75d, so the construction lives here instead of being re-derived
three times. The labels then stay identical across every figure in the set: a group called
`l->Endo_capilar` means the same cells wherever the PI reads it.

The frame is LONG — one row per (cell, group) membership — because the optional whole-letter
`<letter> [all]` group deliberately re-uses cells that are already in their `<letter>-><D>`
group. Group statistics are then a plain groupby; anything that must run once per cell (a
spatial neighbourhood, an embedding) should compute on the unique cell ids first and merge.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from anchor_profiles import read_cell_calls
from denovo_vs_native_pseudobulk import (ALL_SUFFIX, FORCED_SEP, NATIVE_SUFFIX, build_groups,
                                         destinations_from_crosstab)

FORCED_ORIGIN = "forced"
NATIVE_ORIGIN = "native"
ALL_ORIGIN = "all"


def add_origin_arguments(p) -> None:
    """The flags every origin diagnostic shares — spelled as 75d spells them."""
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="anchor_typing.h5 — the semi-supervised label per cell.")
    p.add_argument("--forced-csv", type=Path, required=True,
                   help="75c's forced_named_posteriors.csv (cell_id, top1_type). NOT 75's "
                        "re-scored file, which disagrees with the fit on 45% of cells.")
    p.add_argument("--letter", default="l", help="De-novo letter under test (default l).")
    p.add_argument("--destinations", default="",
                   help="Comma-separated GBmap types to compare against. Default: the letter's "
                        "own largest, taken from --crosstab.")
    p.add_argument("--crosstab", type=Path, default=None,
                   help="denovo_vs_gbmap_crosstab.csv, to derive --destinations automatically.")
    p.add_argument("--top-destinations", type=int, default=3)
    p.add_argument("--min-destination-pct", type=float, default=1.0)
    p.add_argument("--no-all-column", action="store_true",
                   help="Omit the '<letter> [all]' whole-letter reference group.")


def group_order(letter: str, destinations: list[str], include_all: bool) -> list[str]:
    """Forced subset interleaved with its native counterpart, so each pair reads together."""
    order = []
    for dest in destinations:
        order += [f"{letter}{FORCED_SEP}{dest}", f"{dest}{NATIVE_SUFFIX}"]
    if include_all:
        order.append(f"{letter}{ALL_SUFFIX}")
    return order


def load_origin(args) -> tuple[pd.DataFrame, list[str]]:
    """Long frame (cell_id, group, origin, destination) plus the groups in reading order."""
    if args.destinations.strip():
        destinations = [d.strip() for d in args.destinations.split(",") if d.strip()]
    elif args.crosstab is not None:
        destinations = destinations_from_crosstab(
            args.crosstab, args.letter, args.top_destinations, args.min_destination_pct)
    else:
        sys.exit("ERROR: pass --destinations, or --crosstab to derive them.")
    if not destinations:
        sys.exit("ERROR: no destinations selected.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]]
    forced = pd.read_csv(args.forced_csv, usecols=["cell_id", "top1_type"])
    labels = calls.merge(forced, on="cell_id", how="inner").set_index("cell_id")
    print(f"Labels joined for {len(labels):,} cells")

    if args.letter not in set(labels["cell_type"]):
        present = sorted(t for t in labels["cell_type"].unique() if len(str(t)) <= 2)
        sys.exit(f"ERROR: --letter {args.letter!r} is not a semi-supervised label. "
                 f"De-novo labels present: {present}")
    missing = [d for d in destinations if d not in set(labels["top1_type"])]
    if missing:
        sys.exit(f"ERROR: destination(s) never appear as a forced call: {missing}")

    include_all = not args.no_all_column
    group, letter_mask = build_groups(labels["cell_type"], labels["top1_type"],
                                      args.letter, destinations, include_all)

    paired = group.dropna()
    rows = [pd.DataFrame({"group": paired.to_numpy()}, index=paired.index)]
    if letter_mask is not None and letter_mask.any():
        all_ids = labels.index[letter_mask.to_numpy()]
        rows.append(pd.DataFrame({"group": f"{args.letter}{ALL_SUFFIX}"}, index=all_ids))

    tidy = pd.concat(rows)
    tidy.index.name = "cell_id"
    is_forced = tidy["group"].str.contains(FORCED_SEP, regex=False)
    is_native = tidy["group"].str.endswith(NATIVE_SUFFIX)
    tidy["origin"] = ALL_ORIGIN
    tidy.loc[is_forced, "origin"] = FORCED_ORIGIN
    tidy.loc[is_native, "origin"] = NATIVE_ORIGIN
    tidy["destination"] = pd.NA
    tidy.loc[is_forced, "destination"] = (
        tidy.loc[is_forced, "group"].str.split(FORCED_SEP, n=1).str[1])
    tidy.loc[is_native, "destination"] = (
        tidy.loc[is_native, "group"].str.slice(0, -len(NATIVE_SUFFIX)))

    order = [g for g in group_order(args.letter, destinations, include_all)
             if g in set(tidy["group"])]
    sizes = tidy["group"].value_counts()
    print(f"{tidy.index.nunique():,} distinct cells in {len(order)} groups:")
    for g in order:
        print(f"  {g:<32} {sizes[g]:>9,}")
    return tidy, order
