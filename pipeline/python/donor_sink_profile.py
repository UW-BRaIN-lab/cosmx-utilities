#!/usr/bin/env python3
"""Per-donor de-novo sink rate, by region, across every donor in a typed AnnData.

Repeats the Wenyu-cohort Pattern A/B split at full-cohort scope. Reads obs ONLY -- the
counts matrix is never touched -- so a 12 GB .h5ad costs a few seconds and a few hundred MB.

The sink label is NOT hardcoded, and no label convention is assumed. Runs differ: bare
de-novo letters awaiting SME annotation, renamed "<name>_denovo" clusters, or a named
Low_signal type. So this ALWAYS prints the real cell_type inventory first, then measures
whichever label you nominate -- and refuses, listing the candidates, if that label is not
in the object. Which cluster is the low-signal sink is a belief, so it must stay a
parameter rather than an assumption.

    python3 donor_sink_profile.py <typed.h5ad> --sink '<label>' [--top-labels 6]
"""
import argparse
import collections
import sys

import h5py
import numpy as np

# InSituType names de-novo clusters from cluster_name_pool = c(letters, "aa", "ab", ...),
# i.e. one or two lowercase letters. Named reference types always carry an uppercase
# letter, digit or separator, so this never collides with a GBmap type.
DENOVO = lambda s: 1 <= len(s) <= 2 and s.isalpha() and s.islower()

CONTRA_HINT = "contralateral"


def decode(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def read_obs_column(obs, name):
    """Return a numpy array of python strings for an obs column, categorical or not."""
    if name not in obs:
        return None
    node = obs[name]
    if isinstance(node, h5py.Group) and "categories" in node:          # modern anndata
        cats = np.array([decode(c) for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
        return out.astype(object)
    if "__categories" in obs and name in obs["__categories"]:           # legacy anndata
        cats = np.array([decode(c) for c in obs["__categories"][name][:]], dtype=object)
        codes = node[:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "").astype(object)
    return np.array([decode(v) for v in node[:]], dtype=object)         # plain column


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad")
    ap.add_argument("--sink", default="b", help="cell_type label believed to be the sink")
    ap.add_argument("--top-labels", type=int, default=6,
                    help="how many cell_type labels to list in the inventory")
    ap.add_argument("--celltype-key", default="cell_type")
    ap.add_argument("--donor-key", default="Case")
    ap.add_argument("--region-key", default="Region")
    args = ap.parse_args()

    with h5py.File(args.h5ad, "r") as f:
        obs = f["obs"]
        cell_type = read_obs_column(obs, args.celltype_key)
        donor = read_obs_column(obs, args.donor_key)
        region = read_obs_column(obs, args.region_key)
    if cell_type is None or donor is None:
        sys.exit(f"obs lacks {args.celltype_key!r} or {args.donor_key!r}")
    if region is None:
        region = np.full(len(donor), "", dtype=object)
    n = len(donor)
    print(f"{n:,} cells, {len(set(donor))} donors\n")

    # ALWAYS show the real label inventory. Label conventions differ between runs -- bare
    # de-novo letters, renamed "<name>_denovo", a named Low_signal type -- and a heuristic
    # that silently matches nothing would otherwise report 0.0% for every donor.
    label_total = collections.Counter(cell_type)
    print(f"cell_type inventory ({len(label_total)} distinct labels), top {args.top_labels}:")
    for lab, cnt in label_total.most_common(args.top_labels):
        mark = "  <-- nominated sink" if lab == args.sink else ""
        print(f"  {lab:28} {cnt:10,}  {cnt / n:6.1%}{mark}")

    denovo_total = collections.Counter(t for t in cell_type if DENOVO(t))
    if denovo_total:
        print(f"\n  bare-letter de-novo labels: {sum(denovo_total.values()):,} "
              f"({sum(denovo_total.values()) / n:.1%}) across {len(denovo_total)}: "
              f"{', '.join(sorted(denovo_total))}")
    else:
        print("\n  no bare-letter de-novo labels -- this run's clusters are already renamed")

    if args.sink not in label_total:
        print(f"\nERROR: --sink {args.sink!r} is not a cell_type in this object.")
        print("Pick one from the inventory above and re-run with --sink '<label>'.")
        sys.exit(2)
    print()

    # per donor x region: cells, sink cells, any-de-novo cells
    agg = collections.defaultdict(lambda: [0, 0, 0])
    for t, d, r in zip(cell_type, donor, region):
        for key in ((d, r), (d, "*")):
            a = agg[key]
            a[0] += 1
            if t == args.sink:
                a[1] += 1
            if DENOVO(t):
                a[2] += 1

    regions = sorted({r for (_, r) in agg if r not in ("*",) and r})
    contra = next((r for r in regions if CONTRA_HINT in r.lower()), None)
    print(f"regions: {regions}")
    print(f"treating {contra!r} as the uninvolved comparator\n"
          if contra else "no contralateral region found -- pattern split skipped\n")

    rows = []
    for d in {d for (d, _) in agg}:
        overall = agg[(d, "*")]
        con = agg.get((d, contra)) if contra else None
        con_rate = con[1] / con[0] if con and con[0] else None
        tumour = [agg[(d, r)] for r in regions if r != contra and (d, r) in agg]
        tum_rate = (sum(a[1] for a in tumour) / sum(a[0] for a in tumour)) if tumour else None
        rows.append((con_rate if con_rate is not None else -1, d, overall,
                     con_rate, tum_rate))

    print(f"{'donor':7} {'cells':>10} {'sink':>8} {'de-novo':>8} "
          f"{'contra':>8} {'tumour':>8} {'lift':>8}  pattern")
    for _, d, overall, con_rate, tum_rate in sorted(rows, reverse=True):
        pct = lambda v: f"{v:7.1%}" if v is not None else "      -"
        lift = (tum_rate - con_rate) if (con_rate is not None and tum_rate is not None) else None
        if con_rate is None:
            pat = "no uninvolved section"
        elif con_rate > 0.60:
            pat = "B: all regions high -> tissue quality"
        elif lift is not None and lift > 0.25:
            pat = "A: tumour-specific"
        elif con_rate < 0.20 and (lift is None or lift < 0.25):
            pat = "clean"
        else:
            pat = "intermediate"
        print(f"{d:7} {overall[0]:10,} {overall[1] / overall[0]:7.1%} "
              f"{overall[2] / overall[0]:7.1%} {pct(con_rate)} {pct(tum_rate)} "
              f"{pct(lift)}  {pat}")

    b_donors = [r[1] for r in rows if r[3] is not None and r[3] > 0.60]
    if b_donors:
        cells = sum(agg[(d, "*")][0] for d in b_donors)
        sink = sum(agg[(d, "*")][1] for d in b_donors)
        all_sink = sum(a[1] for k, a in agg.items() if k[1] == "*")
        print(f"\nPattern B donors ({len(b_donors)}): {', '.join(sorted(b_donors))}")
        print(f"  {cells:,} cells = {cells / n:.1%} of cohort")
        print(f"  {sink:,} sink cells = {sink / all_sink:.1%} of ALL sink cells")


if __name__ == "__main__":
    main()
