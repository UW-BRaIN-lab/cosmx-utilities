#!/usr/bin/env python3
"""Extract the NAMED-only (GBmap) profile matrix from an InSituType run, for a fully
supervised re-score of the same cells.

The question this serves: our k=27 pruned anchor fit is SEMI-supervised — every cell picks
either a named GBmap Core-L4 type or one of the 27 de-novo letters. The PI wants to know
what each de-novo letter would be CALLED if it were not allowed to be de-novo: force every
cell onto its best GBmap type and see where the letters land.

To make that a fair question, the supervised profiles must be the SAME profiles the
semi-supervised fit converged on — not the raw scRNA-seq GBmap panel. InSituType's
`insitutypeML` (the supervised scorer) does NO platform rescaling, so handing it the raw
scRNA-seq reference would mix a scale correction into the answer; the de-novo letters would
then be re-called partly because of the scRNA->CosMx shift rather than because of biology.
The anchor run's `$profiles` are already CosMx-rescaled (it ran rescale = TRUE) and already
restricted to the pruned gene panel, so dropping its de-novo columns and keeping the named
ones gives exactly "the GBmap reference, as this cohort sees it".

So: read /profiles, keep every named column, drop every de-novo letter column.

Reads (stage the h5 down from Kopah first, e.g.:
    s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage4_anchor_pruned/anchor/anchor_typing.h5 .):
  --profiles-h5  InSituType result h5 with /profiles, /profile_genes, /profile_types.
Writes:
  --output       genes x named-cell-types CSV (gene = row index), linear scale, in the
                 format R/flat_posteriors.R --profiles expects.

Usage:
    uv run python pipeline/python/prep_supervised_profiles.py \\
        --profiles-h5 anchor_typing.h5 \\
        --output gbmap_named_profiles.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from anchor_profiles import read_profiles, split_named_denovo


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True,
                   help="InSituType run result h5 with /profiles, /profile_genes, "
                        "/profile_types (from R/insitutype_typing.R).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output named-only profile CSV (genes x cell types).")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    profiles = read_profiles(args.profiles_h5)
    named, denovo = split_named_denovo(list(profiles.columns))
    print(f"Profiles: {profiles.shape[0]:,} genes x {profiles.shape[1]} types "
          f"({len(named)} named + {len(denovo)} de-novo)")
    if not named:
        sys.exit("ERROR: no named types in this run's profiles — nothing to score against.")
    print(f"  dropping {len(denovo)} de-novo columns: {denovo}")

    out = profiles[named]

    # An all-zero column would make insitutypeML's log-likelihood for that type degenerate,
    # and would silently never win — better to fail loudly than to ship a dead reference.
    dead = [t for t in named if not (out[t] > 0).any()]
    if dead:
        sys.exit(f"ERROR: named profile column(s) are all zero: {dead}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} named cell types)")


if __name__ == "__main__":
    main()
