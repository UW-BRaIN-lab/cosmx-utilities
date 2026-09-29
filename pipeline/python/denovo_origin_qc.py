#!/usr/bin/env python3
"""75l: is the boring explanation the right one — are the forced cells just shallower?

A cell needs a de-novo bin when no fixed profile fits it well. Low RNA depth does that on its
own: a cell with 80 counts matches everything equally badly, and the fitted de-novo profile,
free to move, collects it. If `l->Endo_capilar` were simply thin versions of native
Endo_capilar, every other comparison in the 75 series would be measuring depth, and the marker
attenuation 75d reports would have a purely technical cause.

So this checks the four QC axes the cohort carries, per origin group:

    qc_gene_counts     panel RNA depth — the axis the Stage-3a min-50 floor was applied to
    qc_genes_detected  distinct genes seen; separates "thin" from "flat but deep"
    qc_negprobe_prop   negative-probe fraction — background, i.e. how noisy the cell is
    qc_area            segmented area; a depth deficit at equal area is biology, not slicing

The comparison that matters is each forced group against ITS OWN native counterpart, so the
verdict block prints those pairs rather than a cohort-wide ranking (51_qc_by_celltype.sh already
does the cohort-wide version, by cell type).

Reads obs in backed mode — the expression matrix is never loaded.

Inputs:
  --typed-h5ad  cosmx_typed.h5ad (or any AnnData carrying the qc_* columns)
  the shared origin flags (--typing-h5, --forced-csv, --letter, --destinations/--crosstab)
Writes (--output-dir):
  qc_by_origin.csv   per group x metric: n, median, q10/q25/q75/q90
  qc_by_origin.png   box plots, one panel per metric, forced next to native

Usage:
    python pipeline/python/denovo_origin_qc.py \\
        --typed-h5ad cosmx_typed.h5ad --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv --letter l \\
        --destinations Endo_capilar,Endo_arterial,Pericyte --output-dir l_qc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from denovo_origin import (ALL_ORIGIN, FORCED_ORIGIN, NATIVE_ORIGIN, NATIVE_SUFFIX,
                           add_origin_arguments, load_origin)

# (label, obs-column candidates in priority order, log-scale the axis)
METRICS = [
    ("panel counts", ["qc_gene_counts", "nCount_RNA", "total_counts", "nCount"], True),
    ("genes detected", ["qc_genes_detected", "nFeature_RNA", "n_genes_by_counts"], True),
    ("negprobe fraction", ["qc_negprobe_prop", "negprobe_prop"], False),
    ("segmented area", ["qc_area", "Area", "area"], True),
]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    add_origin_arguments(p)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def resolve_metrics(cols) -> list[tuple[str, str, bool]]:
    """(label, column, log) for every metric the file actually carries."""
    found = []
    for label, candidates, log in METRICS:
        for c in candidates:
            if c in cols:
                found.append((label, c, log))
                break
        else:
            print(f"  (no column for '{label}'; tried {candidates})")
    if not found:
        sys.exit(f"ERROR: none of the QC columns are in obs. Present: {list(cols)[:30]}")
    return found


def plot_panels(values: dict, order: list[str], metrics, out: Path) -> None:
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.1 * len(metrics), 0.42 * len(order) + 3))
    axes = np.atleast_1d(axes)
    y = np.arange(len(order))[::-1]
    for ax, (label, col, log) in zip(axes, metrics):
        data = [values[col][g] for g in order]
        colours = ["firebrick" if NATIVE_SUFFIX not in g else "steelblue" for g in order]
        bp = ax.boxplot(data, positions=y, vert=False, widths=0.62, showfliers=False,
                        patch_artist=True, medianprops={"color": "black"})
        for patch, colour in zip(bp["boxes"], colours):
            patch.set_facecolor(colour)
            patch.set_alpha(0.55)
        if log:
            ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels(order if ax is axes[0] else [])
        ax.set_xlabel(f"{label}  ({col})", fontsize=8)
        ax.grid(axis="x", alpha=0.3)
    axes[0].set_title("forced (red) vs native (blue), per QC axis", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tidy, order = load_origin(args)

    print(f"Reading obs from {args.typed_h5ad} (backed)")
    obs = ad.read_h5ad(args.typed_h5ad, backed="r").obs
    metrics = resolve_metrics(obs.columns)
    print(f"QC axes: {', '.join(c for _, c, _ in metrics)}")

    overlap = tidy.index.unique().intersection(obs.index)
    print(f"cell-id overlap: {len(overlap):,} of {tidy.index.nunique():,} compared cells "
          f"({len(overlap) / max(tidy.index.nunique(), 1):.1%})")
    if len(overlap) == 0:
        sys.exit("ERROR: no compared cell ids are present in the typed run — check the join.")

    numeric = pd.DataFrame(
        {c: pd.to_numeric(obs[c], errors="coerce") for _, c, _ in metrics}, index=obs.index)

    values, rows = {c: {} for _, c, _ in metrics}, []
    for grp in order:
        ids = tidy.index[tidy["group"] == grp].unique().intersection(numeric.index)
        block = numeric.loc[ids]
        for _, col, _ in metrics:
            v = block[col].dropna().to_numpy()
            values[col][grp] = v
            q = np.quantile(v, QUANTILES) if len(v) else np.full(len(QUANTILES), np.nan)
            rows.append({"group": grp, "metric": col, "n": len(v),
                         **{f"q{int(p * 100)}": val for p, val in zip(QUANTILES, q)}})
    table = pd.DataFrame(rows)
    table.to_csv(args.output_dir / "qc_by_origin.csv", index=False)
    plot_panels(values, order, metrics, args.output_dir / "qc_by_origin.png")

    print("\n=== each forced group against its own native counterpart (medians) ===")
    pairs = (tidy[tidy["origin"] == FORCED_ORIGIN]
             .groupby("destination", observed=True)["group"].first())
    for dest, forced_grp in pairs.items():
        native_grp = f"{dest}{NATIVE_SUFFIX}"
        if native_grp not in values[metrics[0][1]]:
            continue
        print(f"\n  {forced_grp}  vs  {native_grp}")
        for label, col, _ in metrics:
            a, b = values[col][forced_grp], values[col][native_grp]
            if not len(a) or not len(b):
                continue
            ma, mb = float(np.median(a)), float(np.median(b))
            ratio = ma / mb if mb else np.nan
            # A tolerance, not an exact test: medians of integer counts tie routinely, and
            # "forced is lower" on a 1.00x ratio reads as a finding when it is a tie.
            if not np.isfinite(ratio) or abs(ratio - 1) < 0.02:
                arrow = "the same"
            else:
                arrow = "higher" if ratio > 1 else "lower"
            print(f"    {label:<20} forced {ma:>10.4g}   native {mb:>10.4g}   "
                  f"{ratio:>5.2f}x  forced is {arrow}")

    print("\nREAD: forced cells at LOWER depth and FEWER genes than their native counterpart "
          "means\ndepth is the parsimonious explanation, and every marker comparison in the set "
          "is\nconfounded by it. Forced cells at EQUAL OR HIGHER depth rule that out — the "
          "de-novo\nbin is then collecting something other than thin cells.")
    print(f"\nWrote {args.output_dir}/qc_by_origin.csv and .png")


if __name__ == "__main__":
    main()
