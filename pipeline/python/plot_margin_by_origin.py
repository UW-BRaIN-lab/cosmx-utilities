#!/usr/bin/env python3
"""Render 75i Part 1: the margin distribution per origin group — the PI's lead figure.

75i computes everything this draws and writes it as CSV; nothing in the pipeline ever rendered
it, so the figure the question was actually asked for has never been seen. This is the Mac-side
renderer (the InSituType container has no plotting stack — same split as 75d's heatmaps).

WHAT THE FIGURE SHOWS. The assignment is argmax over the 81 stored log-likelihoods, so

    margin = loglik(<letter>) - loglik(<destination>)

is literally why each cell went where it did. Group A (the letter's cells) sits above zero and
group B (the natively-called cells) below it BY CONSTRUCTION, so the figure is NOT asking
whether they separate. Two things on it are worth reading:

  1. WHERE ZERO FALLS in the pooled density. Through the body of one peak = one continuum cut
     at an arbitrary threshold. In a trough between two modes = two populations.
  2. HOW WIDE each group is, and how much of B sits on the wrong side of zero. A native group
     straddling zero is not a confident population that the letter is stealing from; it is the
     tail of the same distribution.

Posterior probability is deliberately NOT plotted alongside: it saturates at 1.000 for every
group in this cohort and carries no information. The margin is the usable quantity.

Reads a directory of 75i outputs, prefixed as Kopah stores them
(`<letter>_vs_<destination>_margin_density_by_group.csv`) or unprefixed for a single run.

Writes (--output-dir):
  margin_by_origin.png   per-group density, one panel per destination, zero marked
  margin_boundary.png    pooled density per destination, with the trough verdict
  margin_by_origin.csv   the summary table behind them, argmax consistency joined on

Usage:
    python pipeline/python/plot_margin_by_origin.py \\
        --input-dir margin_decomposition --output-dir margin_figures
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DENSITY_SUFFIX = "margin_density_by_group.csv"
RUN_PATTERN = re.compile(rf"^(?P<run>.+?)_{DENSITY_SUFFIX}$")
FORCED_COLOUR, NATIVE_COLOUR = "firebrick", "steelblue"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", type=Path, required=True,
                   help="Directory of 75i outputs (one or more destinations).")
    p.add_argument("--scale", default="margin", choices=["margin", "margin_per_count"],
                   help="Raw margin (default), or divided by the cell's depth.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def discover(input_dir: Path) -> dict[str, dict[str, Path]]:
    """Map run name -> its 75i CSVs. A run is one letter-vs-destination comparison."""
    runs: dict[str, dict[str, Path]] = {}
    for path in sorted(input_dir.glob(f"*{DENSITY_SUFFIX}")):
        m = RUN_PATTERN.match(path.name)
        run = m.group("run") if m else input_dir.name
        runs[run] = {"density": path}
    if not runs:
        sys.exit(f"ERROR: no *{DENSITY_SUFFIX} in {input_dir}.\n"
                 f"       75i writes it only since the per-group patch — re-run 75i with "
                 f"SKIP_GENES=1 if this directory predates that.")
    for run, files in runs.items():
        stem = f"{run}_" if (input_dir / f"{run}_{DENSITY_SUFFIX}").exists() else ""
        for key, name in (("summary", "margin_summary.csv"),
                          ("argmax", "argmax_consistency.csv"),
                          ("verdict", "density_at_zero.csv"),
                          ("pooled", "margin_histogram_margin.csv")):
            candidate = input_dir / f"{stem}{name}"
            if candidate.exists():
                files[key] = candidate
    return runs


def is_native(group: str) -> bool:
    """75i labels its two groups `A: <letter>-><D>` and `B: <D> native`.

    The `A:`/`B:` prefix is built unconditionally by the R script, so match on it alone.
    Falling back to a "native" substring looks harmless and is not: a destination whose own
    NAME contains the word would flip the forced group's colour and its reading.
    """
    return group.strip().lower().startswith("b:")


def summary_table(runs: dict) -> pd.DataFrame:
    rows = []
    for run, files in runs.items():
        if "summary" not in files:
            continue
        summ = pd.read_csv(files["summary"])
        if "argmax" in files:
            argmax = pd.read_csv(files["argmax"])[["group", "pct_clust_is_argmax"]]
            summ = summ.merge(argmax, on="group", how="left")
        summ.insert(0, "run", run)
        rows.append(summ)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_ridges(runs: dict, scale: str, summary: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(len(runs), 1, figsize=(9.0, 2.9 * len(runs)), squeeze=False,
                             sharex=False)
    for ax, (run, files) in zip(axes[:, 0], runs.items()):
        curves = pd.read_csv(files["density"])
        curves = curves[curves["scale"] == scale]
        if curves.empty:
            ax.text(0.5, 0.5, f"no '{scale}' curves for {run}", ha="center", va="center")
            continue
        labels = sorted(curves["group"].unique())
        if sum(is_native(g) for g in labels) != 1:
            print(f"  WARNING: {run} has {labels} — expected exactly one `B: ... native` "
                  f"group. Colours and the annotation below may be wrong.", file=sys.stderr)
        for group, block in curves.groupby("group", sort=True):
            colour = NATIVE_COLOUR if is_native(group) else FORCED_COLOUR
            n = int(block["n"].iloc[0])
            ax.fill_between(block["x"], block["density"], color=colour, alpha=0.40)
            ax.plot(block["x"], block["density"], color=colour, lw=1.4,
                    label=f"{group}  (n={n:,})")
        ax.axvline(0, color="black", lw=1.2, ls="--")
        ax.set_yticks([])
        ax.set_ylabel("density", fontsize=8)
        ax.set_title(run.replace("_", " "), fontsize=10, loc="left")
        ax.legend(fontsize=7, loc="upper right", frameon=False)
        # The share of the native group on the wrong side of zero is the number that decides
        # whether B is a confident population at all — annotate it rather than making the
        # reader estimate it off the curve.
        native = summary[(summary["run"] == run)
                         & summary["group"].map(is_native)] if len(summary) else summary
        if len(native) and "pct_above_zero" in native:
            pct = float(native["pct_above_zero"].iloc[0])
            ax.annotate(f"{pct:.0f}% of the native group scores HIGHER on the letter",
                        xy=(0.01, 0.92), xycoords="axes fraction", fontsize=7.5, color="0.3")
    axes[-1, 0].set_xlabel("loglik(letter) - loglik(destination)"
                           + ("  per count" if scale.endswith("count") else ""))
    fig.suptitle("Why each cell went where it did — margin by origin", y=0.997, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_boundary(runs: dict, out: Path) -> None:
    usable = {r: f for r, f in runs.items() if "pooled" in f}
    if not usable:
        print("  (no margin_histogram_margin.csv; skipping the boundary figure)")
        return
    fig, axes = plt.subplots(len(usable), 1, figsize=(9.0, 2.5 * len(usable)), squeeze=False)
    for ax, (run, files) in zip(axes[:, 0], usable.items()):
        pooled = pd.read_csv(files["pooled"])
        ax.fill_between(pooled["x"], pooled["density"], color="0.6", alpha=0.5)
        ax.plot(pooled["x"], pooled["density"], color="0.25", lw=1.4)
        ax.axvline(0, color="firebrick", lw=1.4, ls="--")
        note = ""
        if "verdict" in files:
            v = pd.read_csv(files["verdict"])
            v = v[v["scale"] == "margin"]
            if len(v):
                ratio, trough = float(v["ratio_at_zero"].iloc[0]), bool(v["trough_at_zero"].iloc[0])
                note = (f"ratio_at_zero {ratio:.3f}  |  trough at zero: {trough}\n"
                        f"calibration: 1.000 one continuum, 0.423 lopsided single population, "
                        f"0.000 + trough two populations")
        ax.annotate(note, xy=(0.01, 0.72), xycoords="axes fraction", fontsize=7.5, color="0.25")
        ax.set_yticks([])
        ax.set_title(run.replace("_", " "), fontsize=10, loc="left")
    axes[-1, 0].set_xlabel("loglik(letter) - loglik(destination), both groups pooled")
    fig.suptitle("Does the boundary fall in a trough, or through a peak?", y=0.997, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    runs = discover(args.input_dir)
    summary = summary_table(runs)
    # Largest destination first, so the panels read in the order the Sankey does rather than
    # alphabetically -- each destination is its own 75i job, so filename order is meaningless.
    if len(summary):
        size = (summary[~summary["group"].map(is_native)]
                .groupby("run")["n"].max().sort_values(ascending=False))
        runs = {r: runs[r] for r in list(size.index) + [r for r in runs if r not in size.index]}
    print(f"{len(runs)} run(s): {', '.join(runs)}")
    if len(summary):
        summary.to_csv(args.output_dir / "margin_by_origin.csv", index=False)
        cols = [c for c in ("run", "group", "n", "median_margin", "pct_above_zero",
                            "median_depth", "pct_clust_is_argmax") if c in summary]
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print("\n" + summary[cols].to_string(index=False))
        if "pct_clust_is_argmax" in summary:
            print("\nA group well below 100% argmax-consistent was NOT assigned by the stored "
                  "logliks,\nso for those cells the margin is not the reason they went where "
                  "they did — read that\ngroup's curve with the caveat.")

    plot_ridges(runs, args.scale, summary, args.output_dir / "margin_by_origin.png")
    plot_boundary(runs, args.output_dir / "margin_boundary.png")
    print(f"\nWrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
