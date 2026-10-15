#!/usr/bin/env python3
"""How each cell type splits across disease groups, for one or more labelings of the same cells.

Built to answer: for a given cell type, what share of its cells is Control vs AD Pure vs AD+LATE vs AD+LATE SORL1, and
does a coarser labeling (InSituTree) hide a group difference the finer one (InSituType K=21) shows? Three views per
labeling, because the first is easy to misread:

  pct_of_type    row percentages: of this type's cells, the % in each group. Mostly mirrors group SIZE (the groups
                 hold 28/29/25/17% of cells), so on its own it says little.
  enrichment     the same percentages divided by each group's share of all cells (1.0 = no skew; 2.0 = twice as many
                 cells as the group's size predicts). Cell-level, so one donor with many cells can dominate.
  donor_level    for each donor, the % of ITS cells that are this type; averaged per group, with a Kruskal-Wallis test
                 across the groups on the per-donor values and Benjamini-Hochberg q across types. The donor is the unit,
                 which is what a group comparison needs. Exploratory: groups are small (6 to 12 donors), and the
                 test says a group differs, not which one.

`crosswalk` follows each label of the FIRST labeling to the leaves of the SECOND that its cells went to, with that leaf's
own donor-level result, so a gradient that does not survive the move is visible.

Confounds that are NOT corrected here: the SORL1 donors all sit on the July-run slides, and depth differs by group.
Usage:
    celltype_by_group.py --labels cell_labels.csv.gz --label celltype_k21 --label celltype_tree_v2 \\
        --names-csv stage4_k21.csv --out-dir composition_by_group
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal

GROUP_COLUMN = "case_group"
DONOR_COLUMN = "donor"
GROUPS = ["Control", "AD Pure", "AD+LATE", "AD+LATE SORL1"]
FDR_THRESHOLD = 0.10


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(p) + 1)
    q = p * len(p) / ranks
    q[order[::-1]] = np.minimum.accumulate(q[order[::-1]])
    return np.minimum(q, 1.0)


def group_shares(table: pd.DataFrame, groups: list[str]) -> pd.Series:
    return table[GROUP_COLUMN].value_counts(normalize=True).reindex(groups)


def pct_and_enrichment(table: pd.DataFrame, label: str, groups: list[str]):
    counts = pd.crosstab(table[label], table[GROUP_COLUMN]).reindex(columns=groups).fillna(0).astype(int)
    row_share = counts.div(counts.sum(axis=1), axis=0)
    enrichment = row_share.div(group_shares(table, groups), axis=1)
    pct = (100 * row_share).round(2)
    pct.insert(0, "cells", counts.sum(axis=1))
    return counts, pct, enrichment.round(2)


def donor_level(table: pd.DataFrame, label: str, groups: list[str]) -> pd.DataFrame:
    """Per-group mean of the per-donor % of cells in each type, Kruskal-Wallis across groups, BH q."""
    per_donor = pd.crosstab(table[DONOR_COLUMN], table[label], normalize="index") * 100
    donor_group = table.drop_duplicates(DONOR_COLUMN).set_index(DONOR_COLUMN)[GROUP_COLUMN].reindex(per_donor.index)
    rows = []
    for t in per_donor.columns:
        values = [per_donor.loc[donor_group == g, t].to_numpy() for g in groups]
        usable = [v for v in values if len(v)]
        try:
            p = float(kruskal(*usable)[1])
        except ValueError:                      # every value identical across all donors
            p = 1.0
        rows.append([t] + [round(float(v.mean()), 2) if len(v) else np.nan for v in values] + [p])
    out = pd.DataFrame(rows, columns=["type"] + [f"mean_pct_{g}" for g in groups] + ["p"]).set_index("type")
    out["q"] = benjamini_hochberg(out["p"].to_numpy())
    out["donors_per_group"] = ", ".join(f"{g}: {int((donor_group == g).sum())}" for g in groups)
    return out


def crosswalk(table: pd.DataFrame, first: str, second: str, first_dl: pd.DataFrame, second_dl: pd.DataFrame,
              groups: list[str]) -> pd.DataFrame:
    """For each label of `first`: its donor-level result, and the top destination leaf of `second` with its own result."""
    x = pd.crosstab(table[first], table[second])
    rows = []
    for k in x.index:
        dest = (x.loc[k] / x.loc[k].sum()).sort_values(ascending=False)
        top, share = dest.index[0], float(dest.iloc[0])
        second_best = dest.index[1] if len(dest) > 1 else ""
        rows.append({
            "type": k, "cells": int(x.loc[k].sum()),
            **{f"first_mean_pct_{g}": first_dl.loc[k, f"mean_pct_{g}"] for g in groups},
            "first_q": round(float(first_dl.loc[k, "q"]), 4),
            "top_destination": top, "top_destination_share": round(share, 3),
            "second_destination": second_best, "second_destination_share": round(float(dest.get(second_best, 0)), 3),
            **{f"dest_mean_pct_{g}": second_dl.loc[top, f"mean_pct_{g}"] for g in groups},
            "dest_q": round(float(second_dl.loc[top, "q"]), 4),
        })
    return pd.DataFrame(rows).set_index("type")


def apply_names(table: pd.DataFrame, label: str, names_csv: Path | None) -> pd.Series:
    """Readable labels: de novo letters become 'o - Oligodendrocyte' when a names table is given."""
    if names_csv is None:
        return table[label]
    names = pd.read_csv(names_csv, dtype=str)
    reading = dict(zip(names["denovo_label"].str.strip(), names["annotation"].str.strip()))
    return table[label].map(lambda x: reading.get(x, x))


def write_workbook(path: Path, results: dict, sheet_names: list[str], cw: pd.DataFrame | None, notes: list[str]) -> None:
    """One workbook: a README, then per labeling the three views, then the crosswalk."""
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame({"README": notes}).to_excel(xw, sheet_name="README", index=False)
        for (label, (_, _, pct, enr, dl)), short in zip(results.items(), sheet_names):
            pct.to_excel(xw, sheet_name=f"{short}_pct_of_type")
            enr.to_excel(xw, sheet_name=f"{short}_enrichment")
            dl.to_excel(xw, sheet_name=f"{short}_donor_level")
        if cw is not None:
            cw.to_excel(xw, sheet_name="crosswalk")
        for ws in xw.book.worksheets:
            ws.column_dimensions["A"].width = 46 if ws.title != "README" else 150
            ws.freeze_panes = "B2"


def plot_enrichment(results: dict, titles: list[str], out_png: Path, not_in_second: set[str]) -> None:
    """Heatmaps of log2(observed / expected) cells per group, one panel per labeling; * marks donor-level q < FDR."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = list(results.items())
    heights = [max(len(r[3]), 4) for _, r in panels]
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2 * len(panels), 0.3 * max(heights) + 2), squeeze=False,
                             gridspec_kw={"wspace": 1.1})
    for ax, (label, (_, counts, _, enr, dl)), title in zip(axes[0], panels, titles):
        log2 = np.log2(enr.clip(lower=0.05))
        order = (log2[GROUPS[-1]] - log2[GROUPS[0]]).sort_values(ascending=False).index
        log2 = log2.loc[order]
        image = ax.imshow(log2.to_numpy(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
        ax.set_xticks(range(len(GROUPS)), GROUPS, rotation=30, ha="right", fontsize=8)
        names = [f"{'* ' if dl.loc[t, 'q'] < FDR_THRESHOLD else ''}{t}{' (not in tree)' if t in not_in_second else ''}"
                 for t in log2.index]
        ax.set_yticks(range(len(names)), names, fontsize=7.5)
        for i, t in enumerate(log2.index):
            for j, g in enumerate(GROUPS):
                ax.text(j, i, f"{enr.loc[t, g]:.1f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(log2.iloc[i, j]) > 1.5 else "black")
        ax.set_title(title, fontsize=10)
    fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.4, label="log2(observed / expected cells)")
    fig.suptitle("Cells per type by group, relative to each group's share of all cells (numbers = observed/expected; "
                 "* = donor-level q < 0.10)", fontsize=9)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--labels", type=Path, required=True, help="cell_labels.csv.gz (donor, case_group, labelings).")
    p.add_argument("--label", action="append", required=True, help="A label column. Repeatable; the first is the finer one.")
    p.add_argument("--names-csv", type=Path, default=None, help="denovo_annotations CSV to name de novo letters.")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--sheet-name", action="append", default=[], help="Short sheet prefix per --label, in order (e.g. K21).")
    p.add_argument("--title", action="append", default=[], help="Figure panel title per --label, in order.")
    p.add_argument("--xlsx", type=Path, default=None, help="Also write one workbook.")
    p.add_argument("--figure", type=Path, default=None, help="Also write the enrichment heatmap PNG.")
    args = p.parse_args()

    table = pd.read_csv(args.labels, dtype={DONOR_COLUMN: str}, usecols=["cell_id", DONOR_COLUMN, GROUP_COLUMN, *args.label])
    unknown = sorted(set(table[GROUP_COLUMN].dropna()) - set(GROUPS))
    if unknown or table[GROUP_COLUMN].isna().any():
        sys.exit(f"ERROR: groups outside {GROUPS}: {unknown}, or missing values")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for label in args.label:
        shown = table.assign(**{label: apply_names(table, label, args.names_csv if label == args.label[0] else None)})
        counts, pct, enrichment = pct_and_enrichment(shown, label, GROUPS)
        dl = donor_level(shown, label, GROUPS)
        results[label] = (shown, counts, pct, enrichment, dl)
        pct.to_csv(args.out_dir / f"{label}_pct_of_type.csv")
        enrichment.to_csv(args.out_dir / f"{label}_enrichment.csv")
        dl.to_csv(args.out_dir / f"{label}_donor_level.csv")
        print(f"{label}: {len(dl)} types, donor-level q<{FDR_THRESHOLD}: {int((dl['q'] < FDR_THRESHOLD).sum())}")
    if len(args.label) >= 2:
        first, second = args.label[0], args.label[1]
        # the two labelings must be compared on the SAME cells, with the finer one's readable names as the key
        both = results[first][0].assign(**{second: table[second]})
        cw = crosswalk(both, first, second, results[first][4], results[second][4], GROUPS)
        cw.to_csv(args.out_dir / "crosswalk.csv")
        print(f"crosswalk: {len(cw)} {first} types -> {second}")
    else:
        cw = None
    names_kept = None
    if args.names_csv is not None:
        ann = pd.read_csv(args.names_csv, dtype=str)
        if "keep" in ann.columns:
            names_kept = {a.strip() for a, k in zip(ann["annotation"], ann["keep"]) if k.strip().lower() not in ("true", "t", "1", "yes", "y")}
    if args.xlsx:
        sheets = args.sheet_name or [lab.replace("celltype_", "")[:12] for lab in args.label]
        notes = [__doc__.split("Usage:")[0].strip(), "",
                 "Sheets per labeling: pct_of_type (row %), enrichment (observed/expected cells), donor_level (per-donor % averaged by group, Kruskal-Wallis p, BH q).",
                 "crosswalk: each label of the first labeling and the leaf of the second that most of its cells went to, with that leaf's donor-level result."]
        write_workbook(args.xlsx, results, sheets, cw, notes)
        print(f"Wrote {args.xlsx}")
    if args.figure:
        titles = args.title or args.label
        plot_enrichment(results, titles, args.figure, names_kept or set())
        print(f"Wrote {args.figure}")


if __name__ == "__main__":
    main()
