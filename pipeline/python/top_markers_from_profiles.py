#!/usr/bin/env python3
"""Top marker genes for EVERY type in an InSituType result (named + de novo), from the profiles.

Reads /profiles (types x genes), /profile_genes and /profile_types from an InSituType result h5
and, for each type, ranks genes by log2((profile + eps) / (mean profile over the OTHER types + eps)).
An expression floor (--min-profile) keeps near-zero genes from topping the list: with no floor, a
gene at 0.01 counts/cell in one type and 0.0005 elsewhere is a "100x marker" of nothing. The
profile value and the other-types mean are written next to each fold-change so a reader can see
how strongly the gene is actually expressed.

Caveat that belongs with any list from this script: profiles are InSituType's fitted mean counts per
cell. With update_reference=TRUE the named types are the reference profiles after rescaling and the
de novo types are fitted to the data, so a named-versus-de-novo comparison is approximate. These are
profile-based markers, not a differential-expression test on single cells.

Usage:
    top_markers_from_profiles.py --profiles-h5 insitutype_result.h5 --output top_markers.csv \\
        [--names names.csv] [--top-n 20] [--min-profile 0.2]
`--names` is a CSV with columns `type,name` that labels types (e.g. a de novo letter -> a reading).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DEFAULT_TOP_N = 20
DEFAULT_MIN_PROFILE = 0.2   # about the median non-zero profile value in the SORL1 K=21 fit
PSEUDOCOUNT = 0.05


def _decode(arr) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]


def read_profiles(path: Path) -> pd.DataFrame:
    """genes x types, whichever orientation the file stores."""
    with h5py.File(path, "r") as f:
        mat = f["profiles"][()]
        genes = _decode(f["profile_genes"][()])
        types = _decode(f["profile_types"][()])
    if mat.shape == (len(types), len(genes)):
        mat = mat.T
    elif mat.shape != (len(genes), len(types)):
        raise ValueError(f"/profiles shape {mat.shape} fits neither orientation of "
                         f"{len(genes)} genes x {len(types)} types")
    return pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types).astype(float)


def top_markers(profiles: pd.DataFrame, top_n: int, min_profile: float) -> pd.DataFrame:
    """Long table: type, rank, gene, log2fc, profile, other_mean."""
    n_types = profiles.shape[1]
    if n_types < 2:
        raise ValueError("need at least two types to rank markers against the others")
    totals = profiles.sum(axis=1)
    rows = []
    for t in profiles.columns:
        own = profiles[t]
        other_mean = (totals - own) / (n_types - 1)
        log2fc = np.log2((own + PSEUDOCOUNT) / (other_mean + PSEUDOCOUNT))
        eligible = own >= min_profile
        ranked = log2fc[eligible].sort_values(ascending=False).head(top_n)
        for rank, (gene, score) in enumerate(ranked.items(), start=1):
            rows.append({"type": t, "rank": rank, "gene": gene, "log2fc": round(float(score), 2),
                         "profile": round(float(own[gene]), 3),
                         "other_mean": round(float(other_mean[gene]), 3)})
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True, help="Long-format CSV.")
    p.add_argument("--names", type=Path, default=None, help="CSV with columns type,name.")
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    p.add_argument("--min-profile", type=float, default=DEFAULT_MIN_PROFILE)
    args = p.parse_args()

    profiles = read_profiles(args.profiles_h5)
    table = top_markers(profiles, args.top_n, args.min_profile)
    if args.names:
        names = pd.read_csv(args.names, dtype=str).set_index("type")["name"]
        table.insert(1, "name", table["type"].map(names).fillna(""))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    short = table.groupby("type")["gene"].count()
    print(f"{profiles.shape[1]} types x {profiles.shape[0]:,} genes; floor {args.min_profile}; "
          f"wrote {args.output} ({len(table)} rows)")
    thin = short[short < args.top_n]
    if len(thin):
        print(f"NOTE: fewer than {args.top_n} genes pass the floor for: "
              f"{', '.join(f'{t} ({n})' for t, n in thin.items())}", file=sys.stderr)


if __name__ == "__main__":
    main()
