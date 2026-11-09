#!/usr/bin/env python3
"""How far apart are the reference cell types, on the genes the fit actually sees?

InSituType assigns a cell by likelihood against every reference profile at once. Two profiles
that are nearly identical split the evidence between them: neither wins convincingly, the cell
falls to a de novo cluster or to a low posterior, and BOTH types lose. The fix is to merge the
pair into one type -- but only where the distinction was never recoverable, which is what this
measures.

The number that matters is the cosine on the PRUNED gene set (73_select_genes.sh), not on the
full panel: the fit never sees the dropped genes, so a pair that is separable only through them
is not separable at all. Pass --keep-genes and the report adds a `delta` column -- how much
each pair moved when the panel was pruned. A large negative delta is a pair the gene selection
made WORSE, which is worth knowing before the selection is re-run on a new cohort.

Cosine is on log1p profiles. Raw cosine between expression profiles is dominated by a handful
of highly expressed genes and reads 0.9+ for almost any pair; log1p spreads the weight across
the panel without centring, so the ranking reflects the whole profile. Read it as a ranking,
not as an absolute: the question "is this pair distinguishable" is answered by where the pair
sits relative to pairs we know the fit does separate.

Merge candidates are the connected components of the graph "cosine >= --threshold". A
component can CHAIN -- a-b and b-c both above the bar with a-c below it -- so components above
--max-component-size are reported as needing a manual cut rather than as a proposed merge.

Writes, under --output-dir:
    reference_collinearity_matrix.csv    types x types, the evaluated gene set
    reference_collinearity_nearest.csv   one row per type: its nearest neighbour
    reference_collinearity_pairs.csv     every pair at or above --report-floor
    reference_merge_candidates.csv       connected components at --threshold

Usage:
    reference_collinearity.py --reference retina_combined_panel.csv \\
        --keep-genes kept_genes.txt --output-dir reference_collinearity/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Pairs at or above this are reported individually; below it only the summary counts apply.
DEFAULT_REPORT_FLOOR = 0.95
# The bar for "these two were never separable". Set from the measured retina reference, where
# every pair the fit demonstrably confuses sits above 0.99.
DEFAULT_THRESHOLD = 0.99
# A merge group larger than this has chained through intermediates and needs a human cut.
DEFAULT_MAX_COMPONENT = 6
SUMMARY_BANDS = (0.99, 0.98, 0.96, 0.95)


def load_reference(path: Path) -> pd.DataFrame:
    """genes x types, first column the gene index."""
    frame = pd.read_csv(path, index_col=0)
    if frame.empty or frame.shape[1] < 2:
        sys.exit(f"ERROR: {path} needs at least two cell-type columns")
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        bad = numeric.columns[numeric.isna().any()].tolist()[:5]
        sys.exit(f"ERROR: non-numeric values in {path} column(s) {bad}")
    return numeric


def restrict(reference: pd.DataFrame, keep_genes: Path | None) -> tuple[pd.DataFrame, int]:
    """The reference on the kept genes, plus how many kept genes it could not supply.

    A kept gene absent from the reference is not an error -- 73 selects on the panel, and a
    panel gene with no reference profile is exactly the case criterion 1 exists to catch.
    """
    if keep_genes is None:
        return reference, 0
    wanted = [line.strip() for line in keep_genes.read_text().splitlines() if line.strip()]
    if not wanted:
        sys.exit(f"ERROR: {keep_genes} is empty")
    shared = [gene for gene in wanted if gene in reference.index]
    if not shared:
        sys.exit(f"ERROR: none of the {len(wanted)} genes in {keep_genes} are in the reference")
    return reference.loc[shared], len(wanted) - len(shared)


def cosine_matrix(reference: pd.DataFrame) -> pd.DataFrame:
    """log1p cosine between every pair of type columns."""
    values = np.log1p(reference.to_numpy(dtype=float))
    norms = np.linalg.norm(values, axis=0)
    if (norms == 0).any():
        dead = reference.columns[norms == 0].tolist()
        sys.exit(f"ERROR: all-zero profile(s) in the reference: {dead}")
    unit = values / norms
    return pd.DataFrame(unit.T @ unit, index=reference.columns, columns=reference.columns)


def pair_table(matrix: pd.DataFrame, floor: float) -> pd.DataFrame:
    """Every unordered pair at or above `floor`, strongest first."""
    values = matrix.to_numpy()
    rows, cols = np.triu_indices(len(matrix), k=1)
    cosines = values[rows, cols]
    keep = cosines >= floor
    table = pd.DataFrame({
        "type_a": matrix.index.to_numpy()[rows[keep]],
        "type_b": matrix.columns.to_numpy()[cols[keep]],
        "cosine": cosines[keep],
    })
    return table.sort_values("cosine", ascending=False, ignore_index=True)


def nearest_table(matrix: pd.DataFrame) -> pd.DataFrame:
    """Each type's closest other type. Ascending: the most distinguishable types first."""
    values = matrix.to_numpy().copy()
    np.fill_diagonal(values, -np.inf)
    best = values.argmax(axis=1)
    table = pd.DataFrame({
        "type": matrix.index,
        "nearest": matrix.columns.to_numpy()[best],
        "cosine": values[np.arange(len(matrix)), best],
    })
    return table.sort_values("cosine", ignore_index=True)


def components(matrix: pd.DataFrame, threshold: float) -> list[list[str]]:
    """Connected components of `cosine >= threshold`, largest first."""
    types = list(matrix.index)
    adjacency = matrix.to_numpy() >= threshold
    np.fill_diagonal(adjacency, False)
    seen: set[int] = set()
    groups: list[list[str]] = []
    for start in range(len(types)):
        if start in seen:
            continue
        stack, group = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            group.append(node)
            for other in np.flatnonzero(adjacency[node]):
                if other not in seen:
                    seen.add(int(other))
                    stack.append(int(other))
        if len(group) > 1:
            groups.append(sorted(types[i] for i in group))
    return sorted(groups, key=len, reverse=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", type=Path, required=True,
                        help="genes x cell types CSV, first column the gene index")
    parser.add_argument("--keep-genes", type=Path,
                        help="one gene per line (73_select_genes.sh kept_genes.txt)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                        help=f"merge-candidate bar (default {DEFAULT_THRESHOLD})")
    parser.add_argument("--report-floor", type=float, default=DEFAULT_REPORT_FLOOR,
                        help=f"list pairs at or above this (default {DEFAULT_REPORT_FLOOR})")
    parser.add_argument("--max-component-size", type=int, default=DEFAULT_MAX_COMPONENT,
                        help="larger groups are flagged as chained, not proposed")
    args = parser.parse_args(argv)

    full = load_reference(args.reference)
    print(f"reference: {full.shape[0]:,} genes x {full.shape[1]} cell types")

    evaluated, unmatched = restrict(full, args.keep_genes)
    if args.keep_genes is not None:
        print(f"pruned panel: {evaluated.shape[0]:,} genes from {args.keep_genes} "
              f"({unmatched:,} kept genes have no reference profile)")

    matrix = cosine_matrix(evaluated)
    nearest = nearest_table(matrix)
    pairs = pair_table(matrix, args.report_floor)

    if args.keep_genes is not None:
        full_matrix = cosine_matrix(full)
        full_pairs = pair_table(full_matrix, 0.0).set_index(["type_a", "type_b"])["cosine"]
        pairs["cosine_full_panel"] = (
            pairs.set_index(["type_a", "type_b"]).index.map(full_pairs).to_numpy())
        pairs["delta"] = pairs["cosine"] - pairs["cosine_full_panel"]
        nearest["cosine_full_panel"] = [
            full_matrix.loc[t, n] for t, n in zip(nearest["type"], nearest["nearest"])]
        nearest["delta"] = nearest["cosine"] - nearest["cosine_full_panel"]

    total_pairs = len(matrix) * (len(matrix) - 1) // 2
    print(f"\n{total_pairs:,} pairs. Above:")
    for band in SUMMARY_BANDS:
        n = int((pair_table(matrix, band)).shape[0])
        print(f"  >= {band:.2f}   {n:5d}  ({100 * n / total_pairs:4.1f}%)")

    print("\n================ most distinguishable types first ================")
    with pd.option_context("display.max_rows", None, "display.width", 160):
        print(nearest.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    groups = components(matrix, args.threshold)
    proposed = [g for g in groups if len(g) <= args.max_component_size]
    chained = [g for g in groups if len(g) > args.max_component_size]
    print(f"\n================ merge candidates at cosine >= {args.threshold} ================")
    if not groups:
        print("  none -- every pair is below the bar")
    for group in proposed:
        print(f"  [{len(group)}] " + "  |  ".join(group))
    for group in chained:
        print(f"  CHAINED [{len(group)}] -- cut this one by hand:")
        for member in group:
            print(f"      {member}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output_dir / "reference_collinearity_matrix.csv")
    nearest.to_csv(args.output_dir / "reference_collinearity_nearest.csv", index=False)
    pairs.to_csv(args.output_dir / "reference_collinearity_pairs.csv", index=False)
    pd.DataFrame([{"component": i, "size": len(group), "type": member,
                   "chained": len(group) > args.max_component_size}
                  for i, group in enumerate(groups) for member in group]
                 ).to_csv(args.output_dir / "reference_merge_candidates.csv", index=False)
    print(f"\nWrote 4 tables to {args.output_dir}")

    print("\nREAD: a merge group is a claim that the fit could never have told those types "
          "apart.\nConfirm each one against reference_calibration.csv -- a type that DID "
          "attract pinned\ncells was separable in practice, whatever its cosine says.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
