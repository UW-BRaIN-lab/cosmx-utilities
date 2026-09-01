#!/usr/bin/env python3
"""Cross-tab the semi-supervised anchor calls against a fully supervised GBmap re-score.

Answers the PI's question directly: our k=27 pruned anchor fit let each cell choose either a
named GBmap Core-L4 type or one of 27 de-novo letters. If we take the SAME cells and force
every one onto its best GBmap type (75_supervised_gbmap.sh -> insitutypeML, no de-novo
option), what does each letter get called?

Joins per-cell on cell_id:
  LEFT   semi-supervised label from the anchor typing h5 (/cell_type): letters + named types.
  RIGHT  forced GBmap label from the supervised posteriors CSV (top1_type).

Emits three tables:
  denovo_vs_gbmap_crosstab.csv  counts, rows = de-novo letters only. This is the Sankey the
                                PI asked for; feed it straight to plot_crosstab_sankey.py.
  all_vs_gbmap_crosstab.csv     counts, rows = every semi-supervised label. The named rows are
                                the control: cells the semi-supervised fit already named
                                should mostly keep that name under the forced run, which is
                                what makes the letters' destinations interpretable.
  supervised_mapping_summary.csv  one row per semi-supervised label:
                                n_cells, is_denovo, the top three GBmap destinations with
                                their shares, self_pct (named rows only: fraction re-called as
                                themselves), median forced-call confidence, and n_dest_90pct
                                — how many GBmap types it takes to cover 90% of the row.
                                n_dest_90pct is the dispersion readout: 1 means the letter is
                                essentially a rename of one GBmap type; a large value means it
                                is genuinely off-reference and the forced call is arbitrary.

With --collapse-map, two more tables are emitted at COMPARTMENT resolution
(denovo_vs_compartment_crosstab.csv, all_vs_compartment_crosstab.csv) and the summary gains
compartment_1 / compartment_1_pct / n_comp_90pct. This is usually the resolution to read:
GBmap's fine leaves are mutually collinear on a panel this size, so the named control
re-calls only ~79% of cells onto the same LEAF but ~97% onto the same COMPARTMENT. A letter
needing many leaves but few compartments is scattering WITHIN one lineage — a reference
artefact, not a finding. Both control numbers are printed so the gap stays visible.

CAVEAT worth carrying into the figure: insitutypeML has no "unassigned" option, so EVERY cell
gets a GBmap name whether or not it fits. The confidence columns (and --min-prob) are what
separate "this letter really is that type" from "this letter had to be called something".

Usage:
    uv run python pipeline/python/crosstab_denovo_vs_supervised.py \\
        --typing-h5 anchor_typing.h5 \\
        --posteriors supervised_gbmap_posteriors.csv \\
        --annotations pipeline/reference/denovo_annotations/fullcohort_pruned_k27.csv \\
        --output-dir supervised_gbmap/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from anchor_profiles import is_denovo, read_cell_calls

N_TOP_DESTINATIONS = 3
COVERAGE_FRAC = 0.90
COVERAGE_EPS = 1e-9


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="Semi-supervised anchor typing h5 (/cell_id, /cell_type).")
    p.add_argument("--posteriors", type=Path, required=True,
                   help="Supervised flat_posteriors.R CSV (cell_id, top1_type, top1_prob, ...).")
    p.add_argument("--annotations", type=Path,
                   help="De-novo annotation CSV (denovo_label, annotation, ...) used to give "
                        "the letters readable Sankey labels. Optional.")
    p.add_argument("--collapse-map", type=Path,
                   help="CSV mapping each GBmap type to a coarser compartment (columns "
                        "gbmap_type, compartment; see reference/gbmap_compartments.csv). "
                        "Adds compartment-resolution cross-tabs and summary columns.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--min-prob", type=float, default=0.0,
                   help="Drop forced calls whose top1 posterior is below this before "
                        "cross-tabbing (default 0.0 = keep every cell, since insitutypeML "
                        "always assigns one). Use e.g. 0.8 for a confidence-gated view.")
    return p.parse_args()


def load_display_labels(annotations_path: Path | None) -> dict[str, str]:
    """{letter -> readable label} from an annotations CSV's `annotation` column.

    The annotation strings already carry their letter ("o - Hypoxia"), so they are used
    verbatim; letters missing from the CSV keep their bare letter.
    """
    if annotations_path is None:
        return {}
    ann = pd.read_csv(annotations_path)
    if not {"denovo_label", "annotation"}.issubset(ann.columns):
        sys.exit(f"ERROR: {annotations_path} needs denovo_label + annotation columns.")
    ann = ann.dropna(subset=["annotation"])
    return dict(zip(ann["denovo_label"].astype(str).str.strip(),
                    ann["annotation"].astype(str).str.strip()))


def summarize_row(counts: pd.Series, source: str, prefix: str = "gbmap",
                  self_key: str | None = None, n_top: int = N_TOP_DESTINATIONS) -> dict:
    """Top destinations, dispersion and self-agreement for one semi-supervised label.

    prefix names the emitted destination columns (gbmap_1 / compartment_1 / ...). self_key is
    the destination that counts as "unchanged" for a named source — the label itself at leaf
    resolution, its compartment at compartment resolution; None means use `source`.
    """
    total = counts.sum()
    ranked = counts[counts > 0].sort_values(ascending=False)
    out = {"n_cells": int(total)}
    for rank in range(n_top):
        dest = ranked.index[rank] if rank < len(ranked) else ""
        share = ranked.iloc[rank] / total if rank < len(ranked) else np.nan
        out[f"{prefix}_{rank + 1}"] = dest
        out[f"{prefix}_{rank + 1}_pct"] = round(100 * share, 2) if rank < len(ranked) else np.nan
    # How many destinations to cover COVERAGE_FRAC of the row: 1 = a clean rename, many =
    # the forced call is spread thin and therefore arbitrary. The epsilon keeps a cumulative
    # sum that lands exactly on the threshold (10 equal shares -> 0.8999...) from counting
    # one destination too many.
    cumulative = (ranked / total).cumsum()
    out["n_dest_90pct"] = (int((cumulative < COVERAGE_FRAC - COVERAGE_EPS).sum() + 1)
                           if len(ranked) else 0)
    # Named sources have a "correct" answer under the forced run; de-novo letters do not.
    key = source if self_key is None else self_key
    out["self_pct"] = (round(100 * counts.get(key, 0) / total, 2)
                       if not is_denovo(source) and key else np.nan)
    return out


def load_collapse_map(path: Path) -> dict[str, str]:
    """{gbmap_type -> compartment} from the reference CSV, in first-appearance order."""
    m = pd.read_csv(path)
    if not {"gbmap_type", "compartment"}.issubset(m.columns):
        sys.exit(f"ERROR: {path} needs gbmap_type + compartment columns.")
    m["gbmap_type"] = m["gbmap_type"].astype(str).str.strip()
    m["compartment"] = m["compartment"].astype(str).str.strip()
    dupes = m["gbmap_type"][m["gbmap_type"].duplicated()].tolist()
    if dupes:
        sys.exit(f"ERROR: {path} maps these types more than once: {sorted(set(dupes))}")
    return dict(zip(m["gbmap_type"], m["compartment"]))


def collapse_columns(crosstab: pd.DataFrame, cmap: dict[str, str]) -> pd.DataFrame:
    """Sum a leaf-resolution cross-tab into compartments, keeping the map's column order."""
    unmapped = sorted(set(crosstab.columns) - set(cmap))
    if unmapped:
        sys.exit(f"ERROR: --collapse-map does not cover {len(unmapped)} GBmap type(s): "
                 f"{unmapped}. Add them to the map rather than letting cells vanish.")
    collapsed = crosstab.T.groupby(crosstab.columns.map(cmap)).sum().T
    order = list(dict.fromkeys(cmap[t] for t in crosstab.columns if cmap[t] in collapsed.columns))
    return collapsed[order]


def main() -> None:
    args = parse_args()

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]]
    calls = calls.rename(columns={"cell_type": "semisup"})
    print(f"Semi-supervised anchor calls: {len(calls):,} cells, "
          f"{calls['semisup'].nunique()} labels")

    forced = pd.read_csv(args.posteriors, usecols=["cell_id", "top1_type", "top1_prob"])
    print(f"Supervised GBmap posteriors: {len(forced):,} cells, "
          f"{forced['top1_type'].nunique()} GBmap types")

    joined = calls.merge(forced, on="cell_id", how="inner")
    if joined.empty:
        sys.exit("ERROR: no cell_id overlap — the two runs were not scored on the same cells.")
    unmatched = len(calls) - len(joined)
    if unmatched:
        print(f"WARNING: {unmatched:,} anchor cells ({100 * unmatched / len(calls):.2f}%) had "
              f"no supervised posterior; they are excluded from the cross-tabs.")
    print(f"Joined: {len(joined):,} cells")

    if args.min_prob > 0:
        before = len(joined)
        joined = joined[joined["top1_prob"] >= args.min_prob]
        print(f"--min-prob {args.min_prob}: kept {len(joined):,} / {before:,} cells "
              f"({100 * len(joined) / before:.1f}%)")
        if joined.empty:
            sys.exit("ERROR: --min-prob dropped every cell.")

    confidence = joined.groupby("semisup", observed=True)["top1_prob"].median()
    crosstab = pd.crosstab(joined["semisup"], joined["top1_type"])

    summary = pd.DataFrame(
        [{"semisup_label": src, **summarize_row(crosstab.loc[src], src)}
         for src in crosstab.index]
    )

    # Optional compartment roll-up. GBmap's fine leaves are mutually collinear, so this is
    # usually the resolution the forced call is actually trustworthy at.
    cmap = load_collapse_map(args.collapse_map) if args.collapse_map else None
    comp_ct = collapse_columns(crosstab, cmap) if cmap else None
    if cmap:
        comp_summary = pd.DataFrame(
            [{"semisup_label": src,
              **{k: v for k, v in
                 summarize_row(comp_ct.loc[src], src, prefix="compartment",
                               self_key=cmap.get(src), n_top=2).items()
                 if k != "n_cells"}}
             for src in comp_ct.index]
        ).rename(columns={"n_dest_90pct": "n_comp_90pct", "self_pct": "self_pct_compartment"})
        summary = summary.merge(comp_summary, on="semisup_label")
    summary["is_denovo"] = summary["semisup_label"].map(is_denovo)
    summary["median_forced_prob"] = summary["semisup_label"].map(confidence).round(3)

    display = load_display_labels(args.annotations)
    summary["display_label"] = [display.get(lbl, lbl) for lbl in summary["semisup_label"]]
    summary = summary.sort_values(["is_denovo", "n_cells"], ascending=[False, False])
    columns = ["semisup_label", "display_label", "is_denovo", "n_cells"]
    if cmap:
        columns += ["compartment_1", "compartment_1_pct", "compartment_2", "compartment_2_pct",
                    "n_comp_90pct", "self_pct_compartment"]
    columns += ["gbmap_1", "gbmap_1_pct", "gbmap_2", "gbmap_2_pct", "gbmap_3", "gbmap_3_pct",
                "n_dest_90pct", "median_forced_prob", "self_pct"]
    summary = summary[columns]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_path = args.output_dir / "all_vs_gbmap_crosstab.csv"
    crosstab.rename(index=lambda lbl: display.get(lbl, lbl)).to_csv(all_path)

    denovo_rows = [lbl for lbl in crosstab.index if is_denovo(lbl)]
    if not denovo_rows:
        sys.exit("ERROR: no de-novo letters in the anchor calls — wrong typing h5?")
    denovo_ct = crosstab.loc[denovo_rows]
    # Drop GBmap types no letter ever lands on, so the Sankey's right axis is not padded
    # with dozens of empty nodes.
    denovo_ct = denovo_ct.loc[:, denovo_ct.sum(0) > 0]
    denovo_path = args.output_dir / "denovo_vs_gbmap_crosstab.csv"
    denovo_ct.rename(index=lambda lbl: display.get(lbl, lbl)).to_csv(denovo_path)

    comp_paths = []
    if cmap:
        comp_all = args.output_dir / "all_vs_compartment_crosstab.csv"
        comp_ct.rename(index=lambda lbl: display.get(lbl, lbl)).to_csv(comp_all)
        comp_denovo = args.output_dir / "denovo_vs_compartment_crosstab.csv"
        comp_ct.loc[denovo_rows].rename(index=lambda lbl: display.get(lbl, lbl)).to_csv(comp_denovo)
        comp_paths = [(comp_denovo, comp_ct.loc[denovo_rows].shape), (comp_all, comp_ct.shape)]

    summary_path = args.output_dir / "supervised_mapping_summary.csv"
    summary.to_csv(summary_path, index=False)

    named_mask = ~summary["is_denovo"]
    if named_mask.any():
        named_cells = summary.loc[named_mask, "n_cells"].sum()

        def agreement(col: str) -> float:
            hits = (summary.loc[named_mask, col].fillna(0) / 100
                    * summary.loc[named_mask, "n_cells"]).sum()
            return 100 * hits / named_cells

        print(f"\nCONTROL: of {int(named_cells):,} cells the semi-supervised fit already named,")
        print(f"  {agreement('self_pct'):.1f}% keep the same LEAF type under the forced run.")
        if cmap:
            print(f"  {agreement('self_pct_compartment'):.1f}% keep the same COMPARTMENT.")
            print("  Read compartment-level destinations; treat leaf-level ones as a hint.")

    denovo_cells = summary.loc[summary["is_denovo"], "n_cells"].sum()
    print(f"De-novo letters: {len(denovo_rows)} labels, {int(denovo_cells):,} cells "
          f"({100 * denovo_cells / len(joined):.1f}% of joined)\n")
    shown = ["display_label", "n_cells"]
    if cmap:
        shown += ["compartment_1", "compartment_1_pct", "n_comp_90pct"]
    shown += ["gbmap_1", "gbmap_1_pct", "n_dest_90pct"]
    print(summary[summary["is_denovo"]][shown].to_string(index=False))

    print(f"\nWrote {denovo_path} ({denovo_ct.shape[0]} letters x "
          f"{denovo_ct.shape[1]} GBmap types)")
    print(f"Wrote {all_path} ({crosstab.shape[0]} x {crosstab.shape[1]})")
    for path, shape in comp_paths:
        print(f"Wrote {path} ({shape[0]} x {shape[1]} compartments)")
    print(f"Wrote {summary_path} ({len(summary)} rows)")


if __name__ == "__main__":
    main()
