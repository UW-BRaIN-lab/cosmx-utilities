#!/usr/bin/env python3
"""Post-typing diagnostics: spatial concordance, and whether sibling types are real.

Two questions, both answerable straight after stage 4b with no ground truth:

1. SPATIAL CONCORDANCE. Cross-tabulate the call against an anatomical obs column
   (`Region` for the retina study: Retina / Optic nerve / Gray matter / White matter /
   Adjacent soft tissue). A type confined to one region is behaving; a type spread evenly
   across every region is the "cells inappropriately everywhere" failure. Reported
   data-driven -- modal region and its share per type, plus the same aggregated by source
   atlas -- so no expected-region mapping has to be asserted up front.

2. SIBLING ADJUDICATION. A combined reference built from several source atlases carries
   the same biological cell type once per source. Where those copies retain discriminating
   signal we keep them competing (see pipeline/reference/README.md), but the signal is
   confounded: it mixes real regional biology with a source-level ASSAY offset. NEAT1 is
   the clearest example -- nuclear-retained, so it runs high in single-nucleus atlases, and
   its per-source means differ 6-fold across ALL types, not just astrocytes.

   So for every terminal competition node in the hierarchy (a node whose children are all
   leaves), this asks what actually explains the choice among siblings:

     - if the choice tracks the anatomical column  -> plausibly real regional biology,
       keep the siblings split;
     - if it tracks cell area / total counts / genes detected -> it is the assay offset
       sorting cells by segmentation quality and depth, so collapse the node to its parent.

   Association is normalized mutual information, and the technical covariates are
   discretized into as many quantile bins as there are anatomical categories so the values
   are comparable. ONLY THE SECOND FAILURE MODE LOOKS LIKE SUCCESS -- a node whose siblings
   sort by area will still produce a confident-looking, spatially-structured call, because
   segmentation quality itself varies by tissue. Hence both are always reported together.

Groups come from the hierarchy JSON rather than a separate list, so the thing tested is
exactly the competition the hierarchy defines.

Inputs are the stage-3a cohort AnnData (obs + expression) and the stage-4b result h5. It
does NOT need stage 3b/3c or the stage-4c writeback, so it can run as soon as 80 finishes.

Usage:
    uv run python pipeline/python/diagnose_typing_concordance.py \\
        --combined-h5ad combined_qc.h5ad \\
        --typing-h5 insitutype_result.h5 \\
        --hierarchy pipeline/reference/retina_hierarchy.json \\
        --region-col Region \\
        --out-dir typing_diagnostics
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

# Nuclear-retained lncRNA: the single clearest single-nucleus-vs-whole-cell assay marker,
# and the largest discriminator between the astrocyte copies in the combined reference.
# If sibling assignment is being driven by the assay offset, it shows up here.
DEFAULT_MARKERS = ("NEAT1",)

# Technical covariates written by stage 3a. If sibling choice tracks these rather than
# anatomy, the split is segmentation/depth, not biology.
TECHNICAL_COLS = ("qc_area", "total_counts", "qc_genes_detected", "qc_gene_counts")

# Below this modal share, a type is spatially diffuse -- present everywhere in similar
# proportion, which is what a low-signal catch-all looks like.
DIFFUSE_MODAL_SHARE = 0.40


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-h5ad", type=Path, required=True,
                   help="Stage-3a combined_qc.h5ad (obs + expression).")
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="Stage-4b insitutype_result.h5 (/cell_id, /cell_type, /prob).")
    p.add_argument("--hierarchy", type=Path,
                   help="Hierarchy JSON; terminal competition nodes become the sibling "
                        "groups to adjudicate. Omit to run concordance only.")
    p.add_argument("--region-col", default="Region",
                   help="obs column holding the anatomical label (default: Region).")
    p.add_argument("--markers", default=",".join(DEFAULT_MARKERS),
                   help="Comma-separated genes to report per sibling; '' to skip.")
    p.add_argument("--min-group-cells", type=int, default=200,
                   help="Skip sibling groups with fewer assigned cells than this.")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for CSV output.")
    return p.parse_args()


def _decode(arr) -> np.ndarray:
    """hdf5r writes variable-length UTF-8; h5py hands them back as bytes or str."""
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def terminal_groups(node, path: str = "") -> dict[str, list[str]]:
    """Nodes whose children are ALL leaves -> {node path: sibling leaf names}.

    A node that mixes leaf arrays and sub-dicts still recurses into the sub-dicts, so a
    deep hierarchy yields one group per competition that actually decides between leaves.
    """
    out: dict[str, list[str]] = {}
    if isinstance(node, list):
        return out
    leaf_children = {k: v for k, v in node.items() if isinstance(v, list)}
    for name, leaves in leaf_children.items():
        if len(leaves) >= 2:
            out[f"{path}/{name}" if path else name] = list(leaves)
    for name, child in node.items():
        if not isinstance(child, list):
            out.update(terminal_groups(child, f"{path}/{name}" if path else name))
    return out


def normalized_mutual_info(a: pd.Series, b: pd.Series) -> float:
    """NMI(a, b) = MI / sqrt(H(a) H(b)); 0 when either side is constant."""
    tab = pd.crosstab(a, b).to_numpy(dtype=np.float64)
    n = tab.sum()
    if n == 0:
        return float("nan")
    pij = tab / n
    pi = pij.sum(axis=1, keepdims=True)
    pj = pij.sum(axis=0, keepdims=True)
    nz = pij > 0
    mi = float((pij[nz] * np.log(pij[nz] / (pi @ pj)[nz])).sum())
    hi = float(-(pi[pi > 0] * np.log(pi[pi > 0])).sum())
    hj = float(-(pj[pj > 0] * np.log(pj[pj > 0])).sum())
    if hi <= 0 or hj <= 0:
        return 0.0
    return mi / np.sqrt(hi * hj)


def quantile_bin(values: pd.Series, n_bins: int) -> pd.Series:
    """Quantile-bin a covariate to n_bins so its NMI is comparable to the region NMI.

    NMI is sensitive to the number of categories, so comparing a 5-region association
    against a continuous covariate is meaningless unless the covariate is coarsened to the
    same cardinality. Ties can collapse bins; that only weakens the technical side, which
    is the conservative direction here.
    """
    try:
        return pd.qcut(values, n_bins, labels=False, duplicates="drop").astype("float")
    except ValueError:
        return pd.Series(np.zeros(len(values)), index=values.index)


def load_marker_matrix(adata, genes: list[str]) -> pd.DataFrame:
    """Per-cell counts for a few genes, read column-wise from a backed AnnData."""
    present = [g for g in genes if g in adata.var_names]
    absent = [g for g in genes if g not in adata.var_names]
    if absent:
        print(f"Note: marker(s) not in the panel, skipped: {absent}")
    if not present:
        return pd.DataFrame(index=adata.obs.index)
    sub = adata[:, present].to_memory()
    X = sub.X
    dense = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return pd.DataFrame(dense, index=adata.obs.index, columns=present)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {args.typing_h5}")
    with h5py.File(args.typing_h5, "r") as f:
        typed = pd.DataFrame(
            {"cell_type": _decode(f["cell_type"]),
             "prob": np.asarray(f["prob"][()], dtype=np.float64)},
            index=pd.Index(_decode(f["cell_id"]), name="cell"))
    print(f"  {len(typed):,} cells, {typed['cell_type'].nunique()} types, "
          f"median posterior {typed['prob'].median():.3f}")

    # Backed: this only needs obs plus a couple of gene columns, never the full matrix.
    print(f"Reading {args.combined_h5ad} (backed)")
    adata = ad.read_h5ad(args.combined_h5ad, backed="r")
    if args.region_col not in adata.obs:
        sys.exit(f"ERROR: obs has no '{args.region_col}'. Available: "
                 f"{sorted(adata.obs.columns)[:40]} ...")

    obs = adata.obs.copy()
    markers = [g.strip() for g in args.markers.split(",") if g.strip()]
    if markers:
        obs = obs.join(load_marker_matrix(adata, markers))

    df = obs.join(typed, how="inner")
    n_unmatched = len(obs) - len(df)
    if n_unmatched:
        print(f"WARN: {n_unmatched:,} of {len(obs):,} cells have no typing call "
              f"(excluded)", file=sys.stderr)
    df[args.region_col] = df[args.region_col].astype(str)
    regions = sorted(df[args.region_col].unique())
    n_bins = len(regions)
    print(f"  {len(df):,} typed cells over {n_bins} '{args.region_col}' values: {regions}")

    # --- 1. spatial concordance -------------------------------------------------
    tab = pd.crosstab(df["cell_type"], df[args.region_col])
    share = tab.div(tab.sum(axis=1), axis=0)
    conc = pd.DataFrame({
        "n_cells": tab.sum(axis=1),
        "modal_region": share.idxmax(axis=1),
        "modal_share": share.max(axis=1),
        "median_posterior": df.groupby("cell_type", observed=True)["prob"].median(),
    }).join(share.add_prefix("frac_"))
    conc = conc.sort_values("n_cells", ascending=False)
    conc.to_csv(args.out_dir / "region_concordance.csv")
    tab.to_csv(args.out_dir / "region_crosstab.csv")

    print(f"\n=== 1. spatial concordance ({args.region_col}) ===")
    print(f"{'cell_type':44s} {'n':>9s} {'modal region':>20s} {'share':>7s} {'post':>6s}")
    for t, r in conc.head(25).iterrows():
        print(f"{str(t)[:44]:44s} {int(r['n_cells']):9,d} "
              f"{str(r['modal_region'])[:20]:>20s} {r['modal_share']:7.2f} "
              f"{r['median_posterior']:6.3f}")
    if len(conc) > 25:
        print(f"  ... {len(conc)-25} more types in region_concordance.csv")

    diffuse = conc[(conc["modal_share"] < DIFFUSE_MODAL_SHARE) & (conc["n_cells"] >= 1000)]
    print(f"\nSpatially DIFFUSE types (modal share < {DIFFUSE_MODAL_SHARE:.0%}, "
          f">=1000 cells) — the 'everywhere' failure mode:")
    if diffuse.empty:
        print("  none")
    for t, r in diffuse.iterrows():
        print(f"  {str(t)[:44]:44s} {int(r['n_cells']):9,d} cells, "
              f"modal {r['modal_share']:.2f} ({r['modal_region']})")

    # Same view aggregated by source atlas: a source's types should concentrate in the
    # tissue that atlas came from.
    # Merged duplicate types carry no source prefix, so they appear under their own name
    # -- which is correct: they are pan-tissue by construction and have no home atlas.
    src = df["cell_type"].astype(str).str.split("_").str[0].rename("source_atlas")
    by_src = pd.crosstab(src, df[args.region_col])
    by_src_share = by_src.div(by_src.sum(axis=1), axis=0)
    by_src_share.to_csv(args.out_dir / "region_by_source_atlas.csv")
    print(f"\nRegion share by source atlas (row-normalized):")
    print(by_src_share.round(3).to_string())

    # --- 2. sibling adjudication ------------------------------------------------
    if not args.hierarchy:
        print("\n(no --hierarchy given; skipping sibling adjudication)")
        return

    groups = terminal_groups(json.loads(args.hierarchy.read_text()))
    print(f"\n=== 2. sibling adjudication ({len(groups)} terminal competition nodes) ===")
    print("NMI of the sibling choice against anatomy vs against technical covariates.")
    print("region > technical -> real biology, keep split.  technical > region -> assay "
          "offset, collapse to parent.\n")

    rows = []
    for node, leaves in sorted(groups.items()):
        sub = df[df["cell_type"].isin(leaves)]
        if len(sub) < args.min_group_cells:
            print(f"  {node:34s} SKIPPED ({len(sub):,} cells < {args.min_group_cells})")
            continue
        copy = sub["cell_type"].astype(str)
        if copy.nunique() < 2:
            print(f"  {node:34s} SKIPPED (only 1 sibling used: {copy.iloc[0]})")
            continue
        rec = {"node": node, "n_siblings_used": copy.nunique(),
               "n_siblings_defined": len(leaves), "n_cells": len(sub),
               "median_posterior": sub["prob"].median(),
               f"nmi_{args.region_col}": normalized_mutual_info(
                   copy, sub[args.region_col])}
        for col in TECHNICAL_COLS:
            if col in sub:
                rec[f"nmi_{col}"] = normalized_mutual_info(
                    copy, quantile_bin(sub[col], n_bins))
        tech = {k: v for k, v in rec.items()
                if k.startswith("nmi_") and k != f"nmi_{args.region_col}"}
        rec["max_technical_nmi"] = max(tech.values()) if tech else float("nan")
        rec["worst_technical"] = max(tech, key=tech.get)[4:] if tech else ""
        rec["verdict"] = ("keep split (region-driven)"
                          if rec[f"nmi_{args.region_col}"] >= rec["max_technical_nmi"]
                          else "COLLAPSE? (technical-driven)")
        for g in markers:
            if g in sub:
                per = sub.groupby(copy, observed=True)[g].mean()
                rec[f"{g}_min"], rec[f"{g}_max"] = per.min(), per.max()
                rec[f"{g}_fold"] = per.max() / per.min() if per.min() > 0 else np.inf
                rec[f"{g}_high_in"] = per.idxmax()
        rows.append(rec)

        print(f"  {node:34s} n={len(sub):>8,d}  {args.region_col} NMI="
              f"{rec[f'nmi_{args.region_col}']:.3f}  "
              f"technical NMI={rec['max_technical_nmi']:.3f}"
              f" ({rec['worst_technical']})  post={rec['median_posterior']:.3f}")
        print(f"  {'':34s} -> {rec['verdict']}")
        for g in markers:
            if f"{g}_fold" in rec:
                print(f"  {'':34s}    {g} {rec[f'{g}_fold']:.2f}x across siblings, "
                      f"highest in {rec[f'{g}_high_in']}")
        counts = copy.value_counts()
        print(f"  {'':34s}    split: " +
              ", ".join(f"{k[:30]} {v/len(sub):.0%}" for k, v in counts.items()))

    if rows:
        out = pd.DataFrame(rows).set_index("node")
        out.to_csv(args.out_dir / "sibling_adjudication.csv")
        flagged = out[out["verdict"].str.startswith("COLLAPSE")]
        print(f"\n{len(flagged)} of {len(out)} nodes look technical-driven"
              + (f": {list(flagged.index)}" if len(flagged) else ""))
    print(f"\nWrote CSVs to {args.out_dir}")


if __name__ == "__main__":
    main()
