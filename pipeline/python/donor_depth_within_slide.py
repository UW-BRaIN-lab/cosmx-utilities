#!/usr/bin/env python3
"""Within-slide comparison of SORL1 carriers vs non-carriers on RNA depth.

Carriers sit only on the earlier-run pilot slides, so a pooled carrier-vs-non-carrier depth gap
mixes run with carrier status. Here each slide is its own control: only slides holding BOTH
carriers and non-carriers are used, and every donor's depth metric is centred on its slide's
non-carrier mean. The statistic is the mean centred value over carriers.

The p-value is an exact permutation test: carrier labels are shuffled WITHIN slides (the number
of carriers per slide is kept), over every arrangement, so it needs no distributional
assumption. With six carriers this is a handful of hundred arrangements. It tests only whether
the within-slide gap is larger than label noise; it cannot say WHY (tissue quality vs biology).

Usage:
    donor_depth_within_slide.py --donor-depth donor_depth.csv [--metric shallow_frac]
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ("shallow_frac", "median_total_counts")


def usable_slides(table: pd.DataFrame) -> list[str]:
    """Slides that hold at least one carrier AND one non-carrier."""
    per_slide = table.groupby("slide_id")["carrier"].agg(["sum", "size"])
    return per_slide.index[(per_slide["sum"] > 0) & (per_slide["sum"] < per_slide["size"])].tolist()


def centred_gap(table: pd.DataFrame, metric: str, carrier: pd.Series) -> float:
    """Mean over carriers of (value - mean of that slide's non-carriers)."""
    baseline = table[~carrier].groupby("slide_id")[metric].mean()
    carriers = table[carrier]
    return float((carriers[metric] - carriers["slide_id"].map(baseline)).mean())


def exact_permutation_p(table: pd.DataFrame, metric: str) -> tuple[float, float, int]:
    """(observed gap, two-sided exact p, number of arrangements), shuffling within slides."""
    observed_carrier = table["carrier"].to_numpy()
    observed = centred_gap(table, metric, table["carrier"])
    slide_arrangements = []
    for _, group in table.groupby("slide_id"):
        k = int(group["carrier"].sum())
        slide_arrangements.append([(group.index[list(combo)], group.index)
                                   for combo in itertools.combinations(range(len(group)), k)])
    gaps = []
    for choice in itertools.product(*slide_arrangements):
        carrier = pd.Series(False, index=table.index)
        for chosen, _ in choice:
            carrier[chosen] = True
        gaps.append(centred_gap(table, metric, carrier))
    gaps = np.array(gaps)
    p = float((np.abs(gaps) >= abs(observed) - 1e-12).mean())
    return observed, p, len(gaps)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--donor-depth", type=Path, required=True, help="donor_depth.csv from 53.")
    args = p.parse_args()

    table = pd.read_csv(args.donor_depth, dtype={"donor": str})
    table["carrier"] = table["carrier"].astype(bool)
    slides = usable_slides(table)
    used = table[table["slide_id"].isin(slides)].reset_index(drop=True)
    print(f"{len(slides)} slides hold both carriers and non-carriers: "
          f"{used['carrier'].sum()} carriers vs {(~used['carrier']).sum()} non-carriers")

    for metric in METRICS:
        print(f"\n=== {metric} ===")
        rows = []
        for slide, group in used.groupby("slide_id"):
            rows.append({"slide": slide,
                         "carriers": ",".join(group.loc[group.carrier, "donor"]),
                         "carrier_mean": group.loc[group.carrier, metric].mean(),
                         "noncarrier_mean": group.loc[~group.carrier, metric].mean()})
        per_slide = pd.DataFrame(rows)
        per_slide["gap"] = per_slide["carrier_mean"] - per_slide["noncarrier_mean"]
        print(per_slide.round(3).to_string(index=False))
        observed, pval, n_arr = exact_permutation_p(used, metric)
        print(f"within-slide centred gap (mean over carriers) = {observed:+.3f}; "
              f"exact two-sided permutation p = {pval:.3f} over {n_arr} arrangements")


if __name__ == "__main__":
    main()
