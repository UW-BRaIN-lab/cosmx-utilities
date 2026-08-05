#!/usr/bin/env python3
"""Concatenate per-slide sharded InSituTree results into one cohort result h5.

The sharded typing array (85b_insitutree_sharded.sh) writes one compact result h5 per slide,
each with the schema pipeline/R/insitutree_typing.R emits and write_celltypes.py consumes:

  /cell_id    cell ids
  /cell_type  finest resolved leaf label ("" for cells unresolved at level 1)
  /prob       posterior probability of the assigned label (NaN when unresolved)

This stacks them (in sorted slide order) into a single insitutree_result.h5 with the SAME
three datasets, so the existing 90_write_celltypes.sh reads it unchanged
(RESULT_BASENAME=insitutree_result.h5). It also prints the cohort-wide label distribution and
the Low_signal fraction — the headline "how many Low_signal cells do we truly have" number —
and writes it alongside as a small CSV.

Usage:
    uv run python pipeline/python/concat_typing_results.py \\
        --results-dir per_slide_results \\
        --output insitutree_result.h5 \\
        --summary-csv insitutree_label_counts.csv
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

# The InSituTree hierarchy parks untypeable cells in this top-level sink (see
# pipeline/reference/insitutree_hierarchy.json). Single-child hierarchy nodes resolve to the
# NODE name, so the sink label is the node "Low_signal", not the leaf "Low_signal_denovo".
LOW_SIGNAL_LABEL = "Low_signal"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, required=True,
                   help="Directory of per-slide result h5s (each /cell_id,/cell_type,/prob).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output cohort insitutree_result.h5.")
    p.add_argument("--summary-csv", type=Path, default=None,
                   help="Optional CSV of label -> count, cohort-wide.")
    return p.parse_args()


def _read_strings(dset) -> np.ndarray:
    """Read an h5 string dataset as a numpy array of Python str (decode bytes if needed)."""
    vals = dset[:]
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, (bytes, bytearray)) else str(v)
                         for v in vals], dtype=object)
    return vals.astype(object)


def _write_strings(f, name: str, values: np.ndarray) -> None:
    f.create_dataset(name, data=np.asarray(values, dtype=object),
                     dtype=h5py.string_dtype())


def main() -> None:
    args = parse_args()

    result_files = sorted(args.results_dir.glob("*.h5"))
    if not result_files:
        print(f"ERROR: no per-slide result h5s in {args.results_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Concatenating {len(result_files)} per-slide results from {args.results_dir}")

    cell_ids: list[np.ndarray] = []
    cell_types: list[np.ndarray] = []
    probs: list[np.ndarray] = []
    for i, rf in enumerate(result_files, 1):
        with h5py.File(rf, "r") as f:
            for key in ("cell_id", "cell_type", "prob"):
                if key not in f:
                    print(f"ERROR: {rf.name} missing /{key}", file=sys.stderr)
                    sys.exit(1)
            cid = _read_strings(f["cell_id"])
            cty = _read_strings(f["cell_type"])
            prb = f["prob"][:].astype(np.float64)
        if not (len(cid) == len(cty) == len(prb)):
            print(f"ERROR: {rf.name} dataset lengths disagree "
                  f"({len(cid)}/{len(cty)}/{len(prb)})", file=sys.stderr)
            sys.exit(1)
        cell_ids.append(cid)
        cell_types.append(cty)
        probs.append(prb)
        print(f"  [{i:>2}/{len(result_files)}] {rf.stem}: {len(cid):,} cells")

    cell_id = np.concatenate(cell_ids)
    cell_type = np.concatenate(cell_types)
    prob = np.concatenate(probs)

    n_dup = len(cell_id) - len(set(cell_id.tolist()))
    if n_dup:
        # cell_id is globally unique in the pipeline; a collision means a slide was typed twice.
        print(f"ERROR: {n_dup:,} duplicate cell_ids across shards — overlapping/duplicated "
              f"per-slide results?", file=sys.stderr)
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output} ({len(cell_id):,} cells)")
    with h5py.File(args.output, "w") as f:
        _write_strings(f, "cell_id", cell_id)
        _write_strings(f, "cell_type", cell_type)
        f.create_dataset("prob", data=prob)

    # Cohort-wide label distribution + the Low_signal headline.
    counts = Counter(("" if t == "" else str(t)) for t in cell_type.tolist())
    n = len(cell_type)
    n_low = counts.get(LOW_SIGNAL_LABEL, 0)
    n_unresolved = counts.get("", 0)
    print("\nCohort-wide label distribution (top 25):")
    for label, c in sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:25]:
        shown = label if label else "(unresolved level-1)"
        print(f"  {shown:<28} {c:>10,}  {100.0 * c / n:5.1f}%")
    print(f"\nLow_signal sink: {n_low:,} / {n:,} = {100.0 * n_low / n:.1f}%")
    print(f"Unresolved at level 1 (empty label): {n_unresolved:,} "
          f"= {100.0 * n_unresolved / n:.2f}%")

    if args.summary_csv is not None:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        import csv
        with open(args.summary_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cell_type", "count", "fraction"])
            for label, c in sorted(counts.items(), key=lambda kv: kv[1], reverse=True):
                w.writerow([label if label else "(unresolved)", c, c / n])
        print(f"Wrote {args.summary_csv}")

    print(f"\nDone. {len(cell_id):,} cells over {len(result_files)} slides.")


if __name__ == "__main__":
    main()
