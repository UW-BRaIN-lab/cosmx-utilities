#!/usr/bin/env python3
"""Merge or drop cell-type columns of an InSituType reference profile matrix.

A combined reference built from several source atlases (see pipeline/reference_builder/)
prefixes type names per source, so one biological cell type recurs once per source. Some
of those copies carry real signal and should stay as separate, competing types; others are
indistinguishable at CosMx depth and only split a population by dataset provenance. This
applies a reviewed decision about which is which.

Whether a duplicate group is distinguishable is a signal-to-noise question, not a cosine
question: two profiles at cosine 0.99 can still differ on a handful of high-count genes,
and two at 0.95 can differ only on genes that are near-zero in a real cell. The test that
decided the committed spec was, per group, how many expected counts of a 400-count cell
sit on genes that differ >=2x between copies -- e.g. astrocyte 117/400 (kept separate)
versus oligodendrocyte 10/400 (merged).

MERGE takes the unweighted MEAN of the source columns. The reference builder already
rescales each source independently (profiles / quantile(profiles, 0.99) * 1000), so the
columns are on a comparable scale and a plain mean is well defined. The mean is also
mildly desirable on its own terms: source-level assay offsets (e.g. NEAT1, which is
nuclear-retained and so runs high in single-nucleus atlases) partially cancel, and CosMx
is neither assay, so an intermediate profile is no worse a target than either extreme.

DROP removes a column outright. Use it for taxonomy QC categories rather than cell types
-- the Allen taxonomy's `Splatter` (low-quality/doublet dump) and `Miscellaneous` sit at
cosine 0.98-0.99 against every cortical neuron type, which makes them generic-neuron
attractors that would soak up real cells from other tissues.

Spec CSV columns:
  action        merge | drop
  source_type   an existing column of --reference
  merged_name   the output column name (merge only; blank for drop)

Usage:
    uv run python pipeline/python/merge_reference_types.py \\
        --reference pipeline/reference/retina_combined_panel.csv \\
        --spec pipeline/reference/retina_type_merges.csv \\
        --output pipeline/reference/retina_combined_panel.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

VALID_ACTIONS = ("merge", "drop")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", type=Path, required=True,
                   help="Reference profile CSV (genes x cell types, gene index first).")
    p.add_argument("--spec", type=Path, required=True,
                   help="Merge/drop spec CSV (action, source_type, merged_name).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output CSV. May be the same path as --reference.")
    return p.parse_args()


def load_spec(spec_path: Path, columns: pd.Index) -> tuple[dict[str, list[str]], list[str]]:
    """Return ({merged_name: [source columns]}, [columns to drop]), fully validated."""
    spec = pd.read_csv(spec_path, dtype=str).fillna("")
    for col in ("action", "source_type", "merged_name"):
        if col not in spec.columns:
            sys.exit(f"ERROR: spec {spec_path} missing '{col}' column")
    spec["action"] = spec["action"].str.strip().str.lower()
    spec["source_type"] = spec["source_type"].str.strip()
    spec["merged_name"] = spec["merged_name"].str.strip()

    bad = sorted(set(spec["action"]) - set(VALID_ACTIONS))
    if bad:
        sys.exit(f"ERROR: unknown action(s) {bad}; expected one of {list(VALID_ACTIONS)}")

    absent = [t for t in spec["source_type"] if t not in columns]
    if absent:
        sys.exit(f"ERROR: spec source_type(s) absent from the reference: {absent}")

    # Each source type may be spoken for once, or the intent is ambiguous.
    dupes = spec["source_type"][spec["source_type"].duplicated()].tolist()
    if dupes:
        sys.exit(f"ERROR: spec lists source_type(s) more than once: {sorted(set(dupes))}")

    merges: dict[str, list[str]] = {}
    for _, row in spec[spec["action"] == "merge"].iterrows():
        if not row["merged_name"]:
            sys.exit(f"ERROR: merge row for '{row['source_type']}' has no merged_name")
        merges.setdefault(row["merged_name"], []).append(row["source_type"])

    drops = spec.loc[spec["action"] == "drop", "source_type"].tolist()
    for _, row in spec[spec["action"] == "drop"].iterrows():
        if row["merged_name"]:
            sys.exit(f"ERROR: drop row for '{row['source_type']}' must not set merged_name")

    # A single-source "merge" is a rename, which is legal but almost always a typo in a
    # duplicate-collapsing spec, so say so rather than silently renaming.
    for name, sources in merges.items():
        if len(sources) < 2:
            print(f"WARN: merge '{name}' has only one source ({sources[0]}) — "
                  f"this is a rename, not a merge", file=sys.stderr)

    # A merged name must not collide with a column that survives untouched.
    consumed = set(spec["source_type"])
    survivors = set(columns) - consumed
    clash = sorted(set(merges) & survivors)
    if clash:
        sys.exit(f"ERROR: merged_name(s) collide with untouched reference columns: {clash}")

    return merges, drops


def main() -> None:
    args = parse_args()

    print(f"Reading {args.reference}")
    ref = pd.read_csv(args.reference, index_col=0)
    print(f"  {ref.shape[0]:,} genes x {ref.shape[1]} cell types")

    merges, drops = load_spec(args.spec, ref.columns)

    out = ref.copy()
    for merged_name, sources in merges.items():
        # Unweighted mean of already per-source-rescaled profiles; see module docstring.
        out[merged_name] = ref[sources].mean(axis=1)
        out = out.drop(columns=sources)
        print(f"  merged {len(sources)} -> '{merged_name}': {', '.join(sources)}")
    if drops:
        out = out.drop(columns=drops)
        for d in drops:
            print(f"  dropped '{d}'")

    # Keep a stable, readable column order rather than "merged columns last".
    out = out[sorted(out.columns)]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.index.name = ref.index.name or "gene"
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} cell types; "
          f"{ref.shape[1]} -> {out.shape[1]})")


if __name__ == "__main__":
    main()
