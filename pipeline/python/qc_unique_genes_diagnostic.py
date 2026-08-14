#!/usr/bin/env python3
"""Unique-genes-per-cell by cell type — and what a min-genes QC filter would cost.

Two questions this answers:
  1. Does the Low_signal pool have fewer UNIQUE genes per cell than other types? (descriptive)
  2. Should we add a min-genes QC filter?  -> the collateral test: of the cells a candidate
     threshold would DROP, how many are CNV-malignant tumour we'd be deleting?

Uses the QC metric already computed at stage 1 (qc_genes_detected = # gene probes detected per
cell, joined into obs by concat_qc_anndata.py) — no recompute. cell_type comes from obs or a
typing result h5 (/cell_id, /cell_type). An optional per-cell CNV table splits Low_signal into
malignant vs normal, which is the crux for the filter decision (the InSituCNV finding is that
much of Low_signal is real low-transcript CNV-malignant tumour, [[project_insitucnv_lowsignal]]).

The CNV table is InSituCNV's cell_cnv_table.csv.gz. Its malignant flag is `is_malignant_call`
(= mal_sig > sig_thr, sig_thr = 95th pct of diploid-reference + contralateral-Low_signal cells);
the continuous score is `mal_sig` (or `cnv_score`), usable via --cnv-score-col + --cnv-threshold.

Reads obs in BACKED mode, so a 7.5M-cell h5ad only pulls its obs table, not X.

Outputs (into --output-dir):
  unique_genes_by_celltype.csv   per cell_type: n, median/quartiles of genes + counts
  category_rollup.csv            Low_signal vs malignant vs other-named medians (answers Q1)
  filter_collateral.csv          per candidate threshold: cells dropped, % Low_signal, and
                                  (with --cnv) % CNV-malignant among the dropped (answers Q2)
  lowsignal_by_cnv.csv           (with --cnv) unique-genes distribution of Low_signal split by CNV
  unique_genes_by_celltype.png   box plot of unique genes per cell type (Low_signal highlighted)
  filter_collateral.png          (with --cnv) dropped-cell composition vs threshold

Usage:
    uv run --no-project --with anndata --with pandas --with numpy --with matplotlib python \\
        pipeline/python/qc_unique_genes_diagnostic.py \\
            --h5ad cosmx_clustered.h5ad --celltype-col cell_type \\
            [--typing insitutree_result.h5] \\
            [--cnv cell_cnv_table.csv.gz --cnv-malignant-col is_malignant_call] \\
            [ or  --cnv cell_cnv_table.csv.gz --cnv-score-col mal_sig --cnv-threshold <sig_thr> ] \\
            [--malignant-labels AC-like,MES-like,NPC-like,OPC-like,Hypoxia_denovo,...] \\
            --output-dir qc_unique_genes
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
DEFAULT_THRESHOLDS = "10,20,30,50,75,100,150,200"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", type=Path, required=True, help="QC'd/clustered h5ad (obs metrics).")
    p.add_argument("--genes-col", default="qc_genes_detected", help="obs col: unique genes/cell.")
    p.add_argument("--counts-col", default="qc_gene_counts", help="obs col: gene counts/cell.")
    p.add_argument("--celltype-col", default="cell_type", help="obs col: cell type (or via --typing).")
    p.add_argument("--typing", type=Path, default=None,
                   help="Optional typing result h5 (/cell_id,/cell_type) or CSV to supply cell_type.")
    p.add_argument("--lowsignal-label", default="Low_signal",
                   help="Substring identifying the Low_signal / sink label(s).")
    p.add_argument("--malignant-labels", default=None,
                   help="Comma list of tumour cell types for the category rollup (else inferred "
                        "as any label containing 'like'/'MES'/'NPC'/'OPC'/'denovo').")
    p.add_argument("--cnv", type=Path, default=None, help="Per-cell CNV CSV (cell_id + flag/score).")
    p.add_argument("--cnv-malignant-col", default="malignant",
                   help="Boolean CNV column; or use --cnv-score-col + --cnv-threshold.")
    p.add_argument("--cnv-score-col", default=None)
    p.add_argument("--cnv-threshold", type=float, default=None)
    p.add_argument("--thresholds", default=DEFAULT_THRESHOLDS, help="Candidate min-genes cutoffs.")
    p.add_argument("--output-dir", type=Path, default=Path("qc_unique_genes"))
    return p.parse_args()


def _decode(a) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in a])


def load_cell_type(typing: Path) -> pd.Series:
    if typing.suffix in (".h5", ".hdf5"):
        import h5py
        with h5py.File(typing, "r") as f:
            return pd.Series(_decode(f["cell_type"][()]), index=_decode(f["cell_id"][()]))
    df = pd.read_csv(typing)
    idc = "cell_id" if "cell_id" in df.columns else df.columns[0]
    ctc = "cell_type" if "cell_type" in df.columns else df.columns[1]
    return df.set_index(idc)[ctc]


def is_malignant(labels: pd.Series, explicit: str | None) -> np.ndarray:
    if explicit:
        wanted = {s.strip() for s in explicit.split(",") if s.strip()}
        return np.asarray(labels.isin(wanted))
    pat = r"like|MES|NPC|OPC|_denovo|Hypoxia|Astrocytic|Interferon"
    return np.asarray(labels.str.contains(pat, case=False, regex=True, na=False))


def stats_frame(gb) -> pd.DataFrame:
    """Per-group n / mean / quantiles for a SeriesGroupBy (robust across pandas versions)."""
    out = pd.DataFrame({"n": gb.size().astype(int), "mean": gb.mean().round(1)})
    qs = gb.quantile(QUANTILES).unstack()
    qs.columns = [f"p{int(round(c * 100))}" for c in qs.columns]
    return out.join(qs)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    obs = ad.read_h5ad(args.h5ad, backed="r").obs.copy()
    print(f"obs: {len(obs):,} cells, {obs.shape[1]} columns")
    for col in (args.genes_col, args.counts_col):
        if col not in obs.columns:
            sys.exit(f"ERROR: obs has no '{col}'. Available: {list(obs.columns)}")

    # cell_type from obs or the typing result
    if args.typing is not None:
        ct = load_cell_type(args.typing).reindex(obs.index)
    elif args.celltype_col in obs.columns:
        ct = obs[args.celltype_col].astype(str)
    else:
        sys.exit(f"ERROR: no '{args.celltype_col}' in obs and no --typing given.")
    df = pd.DataFrame({
        "cell_type": ct.astype(str).to_numpy(),
        "genes": pd.to_numeric(obs[args.genes_col], errors="coerce").to_numpy(),
        "counts": pd.to_numeric(obs[args.counts_col], errors="coerce").to_numpy(),
    }, index=obs.index).dropna(subset=["genes"])
    df["is_lowsignal"] = df["cell_type"].str.contains(args.lowsignal_label, case=False, na=False)
    df["is_malignant_type"] = is_malignant(df["cell_type"], args.malignant_labels)

    # --- Q1: unique genes by cell type -------------------------------------------------
    by_type = (stats_frame(df.groupby("cell_type")["genes"])
               .join(df.groupby("cell_type")["counts"].median().rename("counts_median"))
               .sort_values("p50"))
    by_type.to_csv(args.output_dir / "unique_genes_by_celltype.csv")

    df["category"] = np.where(df["is_lowsignal"], "Low_signal",
                              np.where(df["is_malignant_type"], "malignant", "other_named"))
    rollup = stats_frame(df.groupby("category")["genes"])
    rollup["counts_median"] = df.groupby("category")["counts"].median()
    rollup.to_csv(args.output_dir / "category_rollup.csv")
    print("\n=== unique genes/cell by category (p50 = median) ===")
    print(rollup[["n", "p10", "p50", "p90", "counts_median"]].to_string())

    # --- optional CNV join -------------------------------------------------------------
    have_cnv = False
    if args.cnv is not None:
        cnv = pd.read_csv(args.cnv)
        idc = "cell_id" if "cell_id" in cnv.columns else cnv.columns[0]
        cnv = cnv.set_index(idc)
        if args.cnv_score_col and args.cnv_threshold is not None:
            mal = (pd.to_numeric(cnv[args.cnv_score_col], errors="coerce") >= args.cnv_threshold)
        elif args.cnv_malignant_col in cnv.columns:
            mal = cnv[args.cnv_malignant_col].astype(str).str.lower().isin({"true", "1", "yes"})
        else:
            sys.exit(f"ERROR: CNV table needs '{args.cnv_malignant_col}' or --cnv-score-col.")
        df["cnv_malignant"] = mal.reindex(df.index).fillna(False).to_numpy()
        have_cnv = True
        ls = df[df["is_lowsignal"]].copy()
        ls["cnv_status"] = np.where(ls["cnv_malignant"], "CNV_malignant", "CNV_normal")
        ls_split = stats_frame(ls.groupby("cnv_status")["genes"])
        ls_split.to_csv(args.output_dir / "lowsignal_by_cnv.csv")
        print("\n=== Low_signal unique genes split by CNV status ===")
        print(ls_split[["n", "p10", "p50", "p90"]].to_string())

    # --- Q2: filter collateral ---------------------------------------------------------
    thresholds = [int(t) for t in args.thresholds.split(",") if t.strip()]
    rows = []
    n_total = len(df)
    for t in thresholds:
        dropped = df[df["genes"] < t]
        row = {"min_genes": t, "cells_dropped": len(dropped),
               "pct_of_all": round(100 * len(dropped) / n_total, 2),
               "pct_dropped_lowsignal": round(100 * dropped["is_lowsignal"].mean(), 1) if len(dropped) else 0.0}
        if have_cnv:
            row["dropped_cnv_malignant"] = int(dropped["cnv_malignant"].sum())
            row["pct_dropped_cnv_malignant"] = round(100 * dropped["cnv_malignant"].mean(), 1) if len(dropped) else 0.0
        rows.append(row)
    collateral = pd.DataFrame(rows)
    collateral.to_csv(args.output_dir / "filter_collateral.csv", index=False)
    print("\n=== filter collateral: what each min-genes cutoff would drop ===")
    print(collateral.to_string(index=False))

    _plot(df, collateral, have_cnv, args)
    print(f"\nWrote diagnostics to {args.output_dir}/")


def _plot(df: pd.DataFrame, collateral: pd.DataFrame, have_cnv: bool, args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # box plot of unique genes per cell type (smallest-median first; Low_signal in red)
    order = df.groupby("cell_type")["genes"].median().sort_values().index.tolist()
    data = [df.loc[df["cell_type"] == t, "genes"].to_numpy() for t in order]
    fig, ax = plt.subplots(figsize=(max(7, len(order) * 0.34), 6))
    bp = ax.boxplot(data, vert=True, showfliers=False, patch_artist=True, widths=0.6)
    for t, box in zip(order, bp["boxes"]):
        ls = args.lowsignal_label.lower() in t.lower()
        box.set(facecolor="#c0392b" if ls else "#4c78a8", alpha=0.85)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel("unique genes per cell")
    ax.set_title("Unique genes per cell by cell type (Low_signal in red)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "unique_genes_by_celltype.png", dpi=130)
    plt.close(fig)

    if have_cnv:
        fig, ax = plt.subplots(figsize=(7, 5))
        x = collateral["min_genes"]
        ax.plot(x, collateral["cells_dropped"], "-o", color="#666", label="cells dropped")
        ax.plot(x, collateral["dropped_cnv_malignant"], "-o", color="#c0392b",
                label="of those, CNV-malignant")
        ax.set_xlabel("min-genes filter threshold")
        ax.set_ylabel("cells")
        ax.set_title("A min-genes filter deletes real (CNV-malignant) tumour")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.output_dir / "filter_collateral.png", dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
