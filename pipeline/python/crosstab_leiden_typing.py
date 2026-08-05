#!/usr/bin/env python3
"""Cross-tabulate InSituTree cell types against the Stage-3c Leiden clustering.

Step-3 Low_signal diagnosis (the cheap, already-computed route): every cohort cell already
carries an unsupervised Stage-3c Leiden cluster (obs['leiden']) AND an InSituTree label
(obs['cell_type']). Crossing them answers "which unsupervised neighbourhood does each
Low_signal cell sit in, and who are its well-typed neighbours?" without any new clustering.

A Low_signal cell whose Leiden cluster is otherwise, say, 80% Oligodendrocyte is most likely
an oligodendrocyte that fell below the lineage-signal floor on the 6k panel; one sitting in a
cluster dominated by malignant states is most likely tumour. That co-residency is the
"neighbours in expression space" signal.

Reads cosmx_typed.h5ad (obs: leiden, cell_type[, insitutype_prob, Region]; obsm['X_umap'] if
present), in backed mode (obs/obsm load to memory; the big X stays on disk). Emits:

  leiden_celltype_counts.csv       leiden x cell_type contingency (counts)
  leiden_celltype_fractions.csv    row-normalised (each leiden cluster's composition)
  lowsignal_by_leiden.csv          per-cluster Low_signal diagnosis (n, %, dominant co-resident
                                   non-Low_signal type = the neighbour hypothesis)
  leiden_celltype_heatmap.png      composition heatmap (leiden x top cell types)
  lowsignal_neighbours.png         Low_signal count per cluster, annotated with the dominant
                                   co-resident typed identity
  lowsignal_umap.png               (if obsm['X_umap']) UMAP, Low_signal vs typed (InSituTree call)
  leiden_umap.png                  (if obsm['X_umap']) same UMAP coloured by Leiden cluster, labelled
  flatcore_umap.png                (if obsm['X_umap']) largest-Low_signal Leiden cluster highlighted —
                                   shows the flat core is one contiguous cluster

Usage:
    uv run python pipeline/python/crosstab_leiden_typing.py \\
        --typed-h5ad cosmx_typed.h5ad \\
        --output-dir leiden_crosstab \\
        [--clustered-h5ad cosmx_clustered.h5ad]   # only if leiden not carried into typed
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402

LOW_SIGNAL_LABEL = "Low_signal"
UMAP_PLOT_MAX = 400_000   # subsample cap for the scatter PNG (7.5M points don't render)
UMAP_SEED = 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True,
                   help="cosmx_typed.h5ad (obs has cell_type; leiden if carried from stage 3c).")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Directory for the crosstab CSVs + figures.")
    p.add_argument("--clustered-h5ad", type=Path, default=None,
                   help="Optional cosmx_clustered.h5ad to source leiden if it is not in the "
                        "typed obs (joined by cell_id).")
    p.add_argument("--leiden-key", default="leiden", help="obs column with the Leiden cluster.")
    p.add_argument("--celltype-key", default="cell_type", help="obs column with the type label.")
    p.add_argument("--low-signal-label", default=LOW_SIGNAL_LABEL,
                   help=f"Label of the InSituTree sink (default {LOW_SIGNAL_LABEL}).")
    p.add_argument("--top-types", type=int, default=20,
                   help="Cell types (by cohort count) to show in the heatmap.")
    return p.parse_args()


def _load_obs(args) -> tuple[pd.DataFrame, np.ndarray | None]:
    print(f"Reading {args.typed_h5ad} (backed)")
    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    obs = adata.obs
    if args.celltype_key not in obs:
        print(f"ERROR: obs has no '{args.celltype_key}'; present: {list(obs.columns)}",
              file=sys.stderr)
        sys.exit(1)

    leiden = None
    if args.leiden_key in obs:
        leiden = obs[args.leiden_key].astype(str).to_numpy()
    elif args.clustered_h5ad is not None:
        print(f"leiden not in typed obs; joining from {args.clustered_h5ad}")
        clu = ad.read_h5ad(args.clustered_h5ad, backed="r")
        if args.leiden_key not in clu.obs:
            print(f"ERROR: '{args.leiden_key}' absent from clustered obs too.", file=sys.stderr)
            sys.exit(1)
        leiden_by_id = clu.obs[args.leiden_key].astype(str)
        leiden = leiden_by_id.reindex(obs.index).to_numpy()
    else:
        print(f"ERROR: '{args.leiden_key}' not in typed obs and no --clustered-h5ad given.",
              file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame({
        "cell_id": obs.index.to_numpy(),
        "leiden": leiden,
        "cell_type": obs[args.celltype_key].astype(str).to_numpy(),
    })
    if "Region" in obs:
        df["Region"] = obs["Region"].astype(str).to_numpy()

    umap = None
    if "X_umap" in adata.obsm:
        umap = np.asarray(adata.obsm["X_umap"])
        print(f"UMAP embedding present: {umap.shape}")

    n_missing = int(pd.isna(df["leiden"]).sum())
    if n_missing:
        print(f"WARN: {n_missing:,} cells have no leiden (dropped from crosstab).")
        keep = ~pd.isna(df["leiden"]).to_numpy()
        df = df[keep].reset_index(drop=True)
        if umap is not None:
            umap = umap[keep]
    print(f"{len(df):,} cells; {df['leiden'].nunique()} leiden clusters; "
          f"{df['cell_type'].nunique()} cell types")
    return df, umap


def _lowsignal_diagnosis(df: pd.DataFrame, low_label: str) -> pd.DataFrame:
    """Per Leiden cluster: Low_signal load + the dominant co-resident NON-Low_signal identity."""
    rows = []
    for cl, g in df.groupby("leiden", sort=False):
        n = len(g)
        n_low = int((g["cell_type"] == low_label).sum())
        non_low = g.loc[g["cell_type"] != low_label, "cell_type"]
        if len(non_low):
            vc = non_low.value_counts()
            dom_type, dom_n = vc.index[0], int(vc.iloc[0])
            dom_frac_of_nonlow = dom_n / len(non_low)
        else:
            dom_type, dom_frac_of_nonlow = "(cluster is all Low_signal)", np.nan
        rows.append({
            "leiden": cl,
            "n_cells": n,
            "n_low_signal": n_low,
            "frac_cluster_low_signal": n_low / n,
            "dominant_coresident_type": dom_type,
            "dominant_coresident_frac_of_nonlow": dom_frac_of_nonlow,
        })
    out = pd.DataFrame(rows).sort_values("n_low_signal", ascending=False).reset_index(drop=True)
    total_low = out["n_low_signal"].sum()
    out["frac_of_all_low_signal"] = out["n_low_signal"] / total_low if total_low else np.nan
    return out


def _heatmap(frac: pd.DataFrame, top_types: list[str], out_path: Path) -> None:
    m = frac.reindex(columns=top_types).fillna(0.0)
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(top_types)),
                                    max(4, 0.35 * len(m))))
    im = ax.imshow(m.to_numpy(), aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(top_types)))
    ax.set_xticklabels(top_types, rotation=90, fontsize=7)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(m.index, fontsize=7)
    ax.set_xlabel("InSituTree cell type")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Per-Leiden-cluster composition (row-normalised)")
    fig.colorbar(im, ax=ax, fraction=0.025, label="fraction of cluster")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _lowsignal_bar(diag: pd.DataFrame, out_path: Path, low_label: str) -> None:
    d = diag[diag["n_low_signal"] > 0].head(30)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.3 * len(d))))
    y = range(len(d))
    ax.barh(list(y), d["n_low_signal"].to_numpy(), color="#c0392b")
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"leiden {c}" for c in d["leiden"]], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(f"# {low_label} cells in cluster")
    ax.set_title(f"{low_label} load per Leiden cluster (annotated: dominant co-resident type)")
    for i, (_, r) in enumerate(d.iterrows()):
        ax.text(r["n_low_signal"], i,
                f"  {r['dominant_coresident_type']} "
                f"({100 * r['frac_cluster_low_signal']:.0f}% LS)",
                va="center", fontsize=6)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _umap_highlight(umap: np.ndarray, df: pd.DataFrame, out_path: Path, low_label: str) -> None:
    n = len(df)
    rng = np.random.default_rng(UMAP_SEED)
    idx = (rng.choice(n, size=UMAP_PLOT_MAX, replace=False) if n > UMAP_PLOT_MAX
           else np.arange(n))
    sub_umap = umap[idx]
    is_low = (df["cell_type"].to_numpy()[idx] == low_label)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(sub_umap[~is_low, 0], sub_umap[~is_low, 1], s=1, c="#cccccc",
               linewidths=0, label="typed")
    ax.scatter(sub_umap[is_low, 0], sub_umap[is_low, 1], s=1, c="#c0392b",
               linewidths=0, label=low_label)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{low_label} cells on the cohort UMAP "
                 f"(n={len(idx):,} of {n:,} shown)")
    ax.legend(markerscale=6, loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_subsample(umap: np.ndarray, n: int) -> np.ndarray:
    rng = np.random.default_rng(UMAP_SEED)
    return (rng.choice(n, size=UMAP_PLOT_MAX, replace=False) if n > UMAP_PLOT_MAX
            else np.arange(n))


def _umap_leiden(umap: np.ndarray, df: pd.DataFrame, out_path: Path) -> None:
    """Cohort UMAP coloured by Leiden cluster, each cluster labelled at its centroid — shows the
    unsupervised cluster structure the crosstab is built on (so 'the flat core is cluster N' is
    visible, not just tabulated)."""
    n = len(df)
    idx = _plot_subsample(umap, n)
    su = umap[idx]
    le = df["leiden"].to_numpy()[idx]
    # Order clusters numerically when possible so colours/labels are stable across runs.
    def _key(c):
        try:
            return (0, int(c))
        except ValueError:
            return (1, c)
    cats = sorted(pd.unique(df["leiden"]), key=_key)
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    cmap = {c: base[i % len(base)] for i, c in enumerate(cats)}
    fig, ax = plt.subplots(figsize=(8, 8))
    for c in cats:
        m = le == c
        if m.any():
            ax.scatter(su[m, 0], su[m, 1], s=1, c=[cmap[c]], linewidths=0)
    # Centroid labels (computed on the full data, not the subsample, so they sit true).
    le_full = df["leiden"].to_numpy()
    for c in cats:
        m = le_full == c
        if m.any():
            ax.text(umap[m, 0].mean(), umap[m, 1].mean(), str(c),
                    fontsize=8, fontweight="bold", ha="center", va="center",
                    color="black",
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Cohort UMAP by Leiden cluster ({len(cats)} clusters)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _umap_flatcore(umap: np.ndarray, df: pd.DataFrame, cluster: str, out_path: Path,
                   low_label: str, ls_frac: float) -> None:
    """Highlight the single largest-Low_signal Leiden cluster (the 'flat core') in red, the rest of
    the Low_signal pool in salmon, and typed cells in grey — the direct visual answer to 'is the
    flat core one contiguous Leiden cluster?'."""
    n = len(df)
    idx = _plot_subsample(umap, n)
    su = umap[idx]
    le = df["leiden"].to_numpy()[idx]
    is_low = (df["cell_type"].to_numpy()[idx] == low_label)
    core = le == cluster
    other_low = is_low & ~core
    typed = ~is_low & ~core
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(su[typed, 0], su[typed, 1], s=1, c="#d2d5da", linewidths=0, label="typed")
    ax.scatter(su[other_low, 0], su[other_low, 1], s=1, c="#f0a58f", linewidths=0,
               label=f"other {low_label}")
    ax.scatter(su[core, 0], su[core, 1], s=1, c="#c0392b", linewidths=0,
               label=f"Leiden {cluster} — flat core")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Flat core = Leiden cluster {cluster} "
                 f"({100 * ls_frac:.0f}% {low_label}, one contiguous region)")
    ax.legend(markerscale=6, loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df, umap = _load_obs(args)

    counts = pd.crosstab(df["leiden"], df["cell_type"])
    fractions = counts.div(counts.sum(axis=1), axis=0)
    counts.to_csv(args.output_dir / "leiden_celltype_counts.csv")
    fractions.to_csv(args.output_dir / "leiden_celltype_fractions.csv")
    print(f"Wrote contingency ({counts.shape[0]} clusters x {counts.shape[1]} types)")

    diag = _lowsignal_diagnosis(df, args.low_signal_label)
    diag.to_csv(args.output_dir / "lowsignal_by_leiden.csv", index=False)
    print("\nLow_signal by Leiden cluster (top 15 by Low_signal count):")
    with pd.option_context("display.max_rows", 15, "display.width", 140):
        print(diag.head(15).to_string(index=False))

    top_types = counts.sum(axis=0).sort_values(ascending=False).head(args.top_types).index.tolist()
    _heatmap(fractions, top_types, args.output_dir / "leiden_celltype_heatmap.png")
    _lowsignal_bar(diag, args.output_dir / "lowsignal_neighbours.png", args.low_signal_label)
    if umap is not None:
        _umap_highlight(umap, df, args.output_dir / "lowsignal_umap.png", args.low_signal_label)
        _umap_leiden(umap, df, args.output_dir / "leiden_umap.png")
        # diag is sorted by Low_signal count desc, so row 0 is the largest-Low_signal cluster.
        core = diag.iloc[0]
        _umap_flatcore(umap, df, core["leiden"], args.output_dir / "flatcore_umap.png",
                       args.low_signal_label, core["frac_cluster_low_signal"])
        print(f"Wrote lowsignal_umap.png, leiden_umap.png, flatcore_umap.png "
              f"(flat core = leiden {core['leiden']}, "
              f"{100 * core['frac_cluster_low_signal']:.0f}% {args.low_signal_label})")

    print(f"\nDone. Crosstab + diagnosis in {args.output_dir}")


if __name__ == "__main__":
    main()
