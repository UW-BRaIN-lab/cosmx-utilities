#!/usr/bin/env python3
"""Is each de novo cluster an unanchored known type, or something distinct / flat?

For every de novo cluster (a, b, c, ...) in an InSituType result h5 this reports its size, its
cosine similarity to every named reference profile, the best and second-best named match, and its
cosine to the average named profile. Reading the table:

  * best_cosine high (>= --known-threshold) with a margin over the runner-up: the cluster looks like
    a known type the fit did not anchor. Merge it into that name (InSituTree does this) rather than
    hunting for a better reference. Retina de novo clusters of this kind sat at cosine 0.87-0.94.
  * best_cosine high but margin tiny: it matches several near-identical named types (the excitatory
    neuron profiles are 0.90-0.98 to each other), so the name is ambiguous.
  * cosine_to_mean_named close to best_cosine: the cluster is flat, an average of everything rather
    than a population, which is what low-signal cells look like.
  * best_cosine low: a distinct population the reference does not cover.

Named profiles come from --reference-csv when given (the reference as built, restricted to the
genes in the result) or else from the result's own fitted named columns. With
update_reference=TRUE those fitted columns have been pulled toward the data, which inflates the
cosines, so prefer --reference-csv. Cosine is on raw profile values, the same quantity as the
quoted retina/Brain_6k figures; the centred_best_cosine column re-does it after subtracting the mean
profile across all types, since raw cosines between expression profiles are high for any pair.

Writes <output-dir>/denovo_vs_reference.csv (one row per de novo cluster) and
<output-dir>/denovo_vs_reference_matrix.csv (de novo x named cosines), and prints a summary of how
much of the cohort sits in known-like, ambiguous, flat and distinct de novo clusters.

Usage:
    denovo_vs_reference_cosine.py --result-h5 insitutype_result.h5 \\
        --reference-csv seaad_mtg_panel.csv --output-dir denovo_vs_reference/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from anchor_profiles import read_cell_calls, read_profiles, split_named_denovo

DEFAULT_KNOWN_THRESHOLD = 0.85
DEFAULT_AMBIGUOUS_MARGIN = 0.02
# "Flat": close to the average named profile (>= known threshold) and about as close to it as to
# its best match. Closeness to the average must be high in absolute terms, or a cluster unlike
# everything would count as flat too.
FLAT_GAP = 0.02

KNOWN_LIKE = "known-like"
AMBIGUOUS = "known-like, ambiguous"
FLAT = "flat"
DISTINCT = "distinct"


def cosine_matrix(rows: pd.DataFrame, cols: pd.DataFrame) -> pd.DataFrame:
    """Cosine between every column of `rows` and every column of `cols` (both genes x types)."""
    a = rows.to_numpy(dtype=float)
    b = cols.to_numpy(dtype=float)
    a = a / np.maximum(np.linalg.norm(a, axis=0, keepdims=True), 1e-12)
    b = b / np.maximum(np.linalg.norm(b, axis=0, keepdims=True), 1e-12)
    return pd.DataFrame(a.T @ b, index=rows.columns, columns=cols.columns)


def named_reference(profiles: pd.DataFrame, named: list[str],
                    reference_csv: Path | None) -> pd.DataFrame:
    """Named profiles on the result's genes: the reference as built, else the fitted columns."""
    if reference_csv is None:
        return profiles[named]
    reference = pd.read_csv(reference_csv, index_col=0)
    missing = [t for t in named if t not in reference.columns]
    if missing:
        sys.exit(f"ERROR: named types in the result are absent from {reference_csv}: {missing}")
    shared = profiles.index.intersection(reference.index)
    if len(shared) < 0.5 * len(profiles.index):
        sys.exit(f"ERROR: only {len(shared)} of {len(profiles.index)} result genes are in "
                 f"{reference_csv}; gene naming differs.")
    print(f"Using {len(shared):,} genes shared with {reference_csv.name}")
    return reference.loc[shared, named]


def classify(best: float, margin: float, to_mean: float,
             known_threshold: float, ambiguous_margin: float) -> str:
    if to_mean >= known_threshold and to_mean >= best - FLAT_GAP:
        return FLAT
    if best < known_threshold:
        return DISTINCT
    return AMBIGUOUS if margin < ambiguous_margin else KNOWN_LIKE


def denovo_table(profiles: pd.DataFrame, calls: pd.DataFrame, reference_csv: Path | None,
                 known_threshold: float, ambiguous_margin: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    named, denovo = split_named_denovo(list(profiles.columns))
    if not named or not denovo:
        sys.exit(f"ERROR: need named and de novo types; found {len(named)} named, "
                 f"{len(denovo)} de novo.")
    reference = named_reference(profiles, named, reference_csv)
    fitted = profiles.loc[reference.index, denovo]

    raw = cosine_matrix(fitted, reference)
    # Gene-wise: shift every profile by the mean over all named+de novo types, then cosine again.
    centre = pd.concat([reference, fitted], axis=1).mean(axis=1)
    centred = cosine_matrix(fitted.sub(centre, axis=0), reference.sub(centre, axis=0))
    mean_named = reference.mean(axis=1).to_frame("mean_named")
    to_mean = cosine_matrix(fitted, mean_named)["mean_named"]

    sizes = calls["cell_type"].value_counts()
    total = int(sizes.sum())
    rows = []
    for t in denovo:
        order = raw.loc[t].sort_values(ascending=False)
        best, second = float(order.iloc[0]), float(order.iloc[1]) if len(order) > 1 else float("nan")
        margin = best - second
        n = int(sizes.get(t, 0))
        rows.append({
            "cluster": t, "n_cells": n, "pct_cohort": round(100 * n / total, 2),
            "best_named": order.index[0], "best_cosine": round(best, 3),
            "second_named": order.index[1] if len(order) > 1 else "",
            "second_cosine": round(second, 3), "margin": round(margin, 3),
            "centred_best_cosine": round(float(centred.loc[t, order.index[0]]), 3),
            "cosine_to_mean_named": round(float(to_mean[t]), 3),
            "verdict": classify(best, margin, float(to_mean[t]), known_threshold, ambiguous_margin),
        })
    table = pd.DataFrame(rows).sort_values("n_cells", ascending=False).reset_index(drop=True)
    return table, raw.round(3)


def print_summary(table: pd.DataFrame, calls: pd.DataFrame) -> None:
    total = len(calls)
    named_cells = total - int(table["n_cells"].sum())
    print(f"\n{total:,} cells; named {100 * named_cells / total:.1f}%; "
          f"de novo {100 * table['n_cells'].sum() / total:.1f}% in {len(table)} clusters")
    for verdict, grp in table.groupby("verdict"):
        print(f"  {verdict:<24} {len(grp):>3} clusters  {100 * grp['n_cells'].sum() / total:5.1f}% of cells")
    print("\n" + table[["cluster", "n_cells", "pct_cohort", "best_named", "best_cosine",
                        "margin", "cosine_to_mean_named", "verdict"]].to_string(index=False))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--result-h5", type=Path, required=True, help="InSituType insitutype_result.h5.")
    p.add_argument("--reference-csv", type=Path, default=None,
                   help="genes x types reference as built (default: the result's fitted named columns).")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--known-threshold", type=float, default=DEFAULT_KNOWN_THRESHOLD,
                   help="best cosine at or above this = looks like a known type.")
    p.add_argument("--ambiguous-margin", type=float, default=DEFAULT_AMBIGUOUS_MARGIN,
                   help="best minus second-best cosine below this = ambiguous among named types.")
    args = p.parse_args()

    profiles = read_profiles(args.result_h5)
    calls = read_cell_calls(args.result_h5)
    table, matrix = denovo_table(profiles, calls, args.reference_csv,
                                 args.known_threshold, args.ambiguous_margin)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "denovo_vs_reference.csv", index=False)
    matrix.to_csv(args.output_dir / "denovo_vs_reference_matrix.csv")
    print_summary(table, calls)
    print(f"\nWrote {args.output_dir}/denovo_vs_reference.csv and denovo_vs_reference_matrix.csv")


if __name__ == "__main__":
    main()
