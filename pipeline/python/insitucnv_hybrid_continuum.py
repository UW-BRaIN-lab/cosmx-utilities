#!/usr/bin/env python3
"""Stage 5g (Phase 2): hybrid state vs continuum among the CNV-high Low_signal cells.

CNV says malignant-or-not; it CANNOT separate a genuine HYBRID cell (one cell co-running
programs from different compartments — the notable case for the lab's NLGN3 neuron-to-glioma
thread) from a CONTINUUM cell (one cell intermediate between adjacent malignant states). This
script does that discrimination on the CNV-high Low_signal working set, using three readouts:

  1. TOP-TWO STRUCTURE (primary) — each cell's top-two profile assignments (from the flat
     posteriors, flat_posteriors.R). Both top-two malignant => CONTINUUM (intermediate on the
     malignant manifold). Malignant + Neuron (cross-compartment) => candidate HYBRID. A weak /
     ambiguous top-1 => reference-limited, iterate (Phase 3).
  2. CO-EXPRESSION CONFIRMATION — score each cell on the malignant modules + the Neuronal_NLGN3
     program (sc.tl.score_genes). A true hybrid co-elevates BOTH, not one-or-the-other.
  3. SPATIAL NULL — a genuine hybrid's neuronal co-expression is cell-intrinsic and spatially
     dispersed; a spillover pseudo-hybrid sits at boundaries with real neurons. Compare each
     candidate's LOCAL NEURON DENSITY (fraction of spatial neighbours typed Neuron) against the
     malignant-cell background: candidates enriched for neuron neighbours are spillover-suspect.

Inputs:
  --typed-h5ad   cosmx_typed.h5ad (raw counts; module scoring).
  --cell-table   cell_cnv_table.csv.gz from insitucnv_lowsignal_diagnostics.py (Phase 1).
  --posteriors   flat_posteriors.csv from flat_posteriors.R (top-K per cell).
  --signatures   gene_signatures.csv (module,gene).
  --hierarchy    insitutree_hierarchy.json (leaf -> compartment map).
Writes (--output-dir): hybrid_continuum_cells.csv.gz (per working-set cell, all readouts),
  toptwo_structure.csv, coexpression.png, spatial_null.png, HYBRID_SUMMARY.txt.

Usage:
    python pipeline/python/insitucnv_hybrid_continuum.py \\
        --typed-h5ad cosmx_typed.h5ad --cell-table diagnostics/cell_cnv_table.csv.gz \\
        --posteriors posteriors/flat_posteriors.csv \\
        --signatures pipeline/reference/gene_signatures.csv \\
        --hierarchy pipeline/reference/insitutree_hierarchy.json --output-dir hybrid
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors

from compare_insitucnv_groups import DEFAULT_MALIGNANT

MALIGNANT_COMPARTMENT = "Malignant"
NEURON_COMPARTMENT = "Neuron"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--cell-table", type=Path, required=True)
    p.add_argument("--posteriors", type=Path, required=True)
    p.add_argument("--signatures", type=Path, required=True)
    p.add_argument("--hierarchy", type=Path, required=True)
    p.add_argument("--malignant-groups", default=",".join(DEFAULT_MALIGNANT))
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--neuronal-module", default="Neuronal_NLGN3")
    p.add_argument("--neuron-type", default="Neuron",
                   help="Profile/cell_type name of mature neurons (for the spatial null).")
    p.add_argument("--k-neighbors", type=int, default=30)
    p.add_argument("--score-percentile", type=float, default=75.0,
                   help="A module score is 'elevated' above this cohort percentile.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_signatures(path: Path) -> dict[str, list[str]]:
    df = pd.read_csv(path, comment="#")
    return {m: g["gene"].astype(str).tolist() for m, g in df.groupby("module", sort=False)}


def leaf_compartments(hierarchy_path: Path) -> dict[str, str]:
    """Map every leaf profile name to its top-level compartment from the hierarchy JSON."""
    tree = json.loads(hierarchy_path.read_text())
    out: dict[str, str] = {}

    def collect(node) -> list[str]:
        leaves: list[str] = []
        if isinstance(node, dict):
            for k, v in node.items():
                child = collect(v)
                leaves += child or [k]
        elif isinstance(node, list):
            for v in node:
                leaves += collect(v)
        elif isinstance(node, str):
            leaves += [node]
        return leaves

    for compartment, sub in tree.items():
        for leaf in collect(sub) or [compartment]:
            out[leaf] = compartment
    return out


def score_modules(typed_h5ad: Path, signatures: dict[str, list[str]]) -> pd.DataFrame:
    """Log-normalise the raw-count typed cohort and score each module (sc.tl.score_genes)."""
    adata = ad.read_h5ad(typed_h5ad)
    if "probe_type" in adata.var:
        adata = adata[:, (adata.var["probe_type"] == "gene").to_numpy()].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    scores = {}
    for module, genes in signatures.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 3:
            print(f"  WARN: module {module} has only {len(present)} genes on panel; skipped",
                  file=sys.stderr)
            continue
        sc.tl.score_genes(adata, present, score_name=f"score_{module}", ctrl_size=50)
        scores[f"score_{module}"] = adata.obs[f"score_{module}"].to_numpy()
        print(f"  scored {module}: {len(present)}/{len(genes)} genes on panel")
    return pd.DataFrame(scores, index=adata.obs_names)


def local_neuron_density(df: pd.DataFrame, neuron_type: str, k: int) -> pd.Series:
    """Per tissue section, fraction of each cell's k spatial neighbours typed as `neuron_type`."""
    out = pd.Series(np.nan, index=df.index)
    if not {"spatial_x", "spatial_y", "tissue_section", "cell_type"}.issubset(df.columns):
        return out
    is_neuron = (df["cell_type"].astype(str) == neuron_type).to_numpy()
    for _, sub in df.groupby("tissue_section", observed=True):
        rows = sub.index
        xy = sub[["spatial_x", "spatial_y"]].to_numpy()
        finite = np.isfinite(xy).all(axis=1)
        if finite.sum() < 2:
            continue
        xy_f = xy[finite]
        kk = min(k, len(xy_f) - 1)
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(xy_f)
        _, idx = nn.kneighbors(xy_f)
        neigh = idx[:, 1:]
        neur = is_neuron[df.index.get_indexer(rows[finite])]
        out.loc[rows[finite]] = neur[neigh].mean(axis=1)
    return out


def classify_toptwo(row, comp, malignant_types, lowsignal_label, prob_gap=0.10) -> str:
    """Top-two structure -> continuum / hybrid_neuronal / cross_other / ambiguous_iterate."""
    t1, t2 = row.get("top1_type"), row.get("top2_type")
    p1, p2 = row.get("top1_prob", np.nan), row.get("top2_prob", np.nan)
    if not isinstance(t1, str):
        return "unknown"
    c1, c2 = comp.get(t1, "other"), comp.get(t2, "other")
    # a weak / flat top-1 (dominated by the Low_signal sink, or no clear winner) -> iterate
    if t1 == f"{lowsignal_label}_denovo" or (np.isfinite(p1) and p1 < 0.5) \
            or (np.isfinite(p1) and np.isfinite(p2) and (p1 - p2) < prob_gap):
        return "ambiguous_iterate"
    comps = {c1, c2}
    if comps == {MALIGNANT_COMPARTMENT}:
        return "continuum"
    if comps == {MALIGNANT_COMPARTMENT, NEURON_COMPARTMENT}:
        return "hybrid_neuronal"
    if MALIGNANT_COMPARTMENT in comps:
        return "cross_other"
    return "non_malignant_toptwo"


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    malignant_types = {s.strip() for s in args.malignant_groups.split(",") if s.strip()}
    signatures = read_signatures(args.signatures)
    comp = leaf_compartments(args.hierarchy)
    mal_modules = [f"score_{m}" for m in ("AClike", "MESlike", "OPClike", "NPClike")
                   if m in signatures]
    neuro_col = f"score_{args.neuronal_module}"

    print("Scoring gene-signature modules on the typed cohort...")
    scores = score_modules(args.typed_h5ad, signatures)

    print("Loading Phase-1 cell table + flat posteriors...")
    cells = pd.read_csv(args.cell_table, index_col=0)
    post = pd.read_csv(args.posteriors, index_col="cell_id")
    df = cells.join(scores, how="left").join(
        post[[c for c in post.columns if c.startswith("top")]], how="left")

    # working set: CNV-high Low_signal (the flat-but-malignant cells Stage 3 must resolve)
    ct = df["cell_type"].astype(str)
    ws = df[(ct == args.lowsignal_label) & df.get("is_malignant_call", False).astype(bool)].copy()
    print(f"Working set (CNV-high Low_signal): {len(ws):,} cells.")
    if not len(ws):
        sys.exit("No CNV-high Low_signal cells in the table; nothing to classify.")

    # (1) top-two structure
    ws["toptwo_class"] = ws.apply(
        lambda r: classify_toptwo(r, comp, malignant_types, args.lowsignal_label), axis=1)
    struct = ws["toptwo_class"].value_counts()
    struct.to_csv(args.output_dir / "toptwo_structure.csv")

    # (2) co-expression: elevated malignant module AND elevated neuronal module
    lines = [f"Hybrid vs continuum — CNV-high Low_signal working set: {len(ws):,} cells\n"]
    lines.append("TOP-TWO STRUCTURE (primary discriminator):")
    for k, v in struct.items():
        lines.append(f"  {k:<22} {v:>8,}  ({v / len(ws):.1%})")

    have_scores = neuro_col in ws.columns and any(m in ws.columns for m in mal_modules)
    if have_scores:
        mal_score = ws[[m for m in mal_modules if m in ws.columns]].max(axis=1)
        neuro_score = ws[neuro_col]
        mal_hi_thr = float(np.nanpercentile(df[[m for m in mal_modules if m in df.columns]]
                                            .max(axis=1), args.score_percentile))
        neuro_hi_thr = float(np.nanpercentile(df[neuro_col], args.score_percentile))
        ws["malignant_module_max"] = mal_score
        ws["neuronal_module"] = neuro_score
        ws["coexpressing"] = (mal_score > mal_hi_thr) & (neuro_score > neuro_hi_thr)
        n_co = int(ws["coexpressing"].sum())
        lines.append(f"\nCO-EXPRESSION (malignant module AND {args.neuronal_module} both above "
                     f"cohort {args.score_percentile:g}th pct):")
        lines.append(f"  co-expressing cells: {n_co:,} ({n_co / len(ws):.1%} of working set)")
        hybrid_top = ws["toptwo_class"] == "hybrid_neuronal"
        if hybrid_top.any():
            conf = float((ws.loc[hybrid_top, "coexpressing"]).mean())
            lines.append(f"  of the top-two 'hybrid_neuronal' candidates, {conf:.1%} also "
                         "co-express by module score (confirmation rate).")

        fig, ax = plt.subplots(figsize=(6.5, 6))
        sub = ws.dropna(subset=["malignant_module_max", "neuronal_module"])
        col = np.where(sub["toptwo_class"] == "hybrid_neuronal", "#d73027",
              np.where(sub["toptwo_class"] == "continuum", "#4575b4", "#bbbbbb"))
        ax.scatter(sub["malignant_module_max"], sub["neuronal_module"], s=3, alpha=0.3, c=col)
        ax.axvline(mal_hi_thr, color="k", ls="--", lw=1); ax.axhline(neuro_hi_thr, color="k", ls="--", lw=1)
        ax.set_xlabel("malignant module score (max of AC/MES/OPC/NPC-like)")
        ax.set_ylabel(f"{args.neuronal_module} score")
        ax.set_title("Co-expression: malignant vs neuronal program\n"
                     "(red=top-two hybrid, blue=continuum; top-right = coincident elevation)")
        fig.tight_layout()
        fig.savefig(args.output_dir / "coexpression.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
    else:
        lines.append("\nCO-EXPRESSION: module scores unavailable; skipped.")

    # (3) spatial null: local neuron density of candidates vs malignant background
    ws["local_neuron_density"] = local_neuron_density(
        df, args.neuron_type, args.k_neighbors).reindex(ws.index)
    mal_bg = local_neuron_density(
        df[df["cell_type"].astype(str).isin(malignant_types)], args.neuron_type,
        args.k_neighbors)
    cand = ws[ws["toptwo_class"] == "hybrid_neuronal"] if have_scores else ws.iloc[:0]
    if have_scores:
        cand = ws[(ws["toptwo_class"] == "hybrid_neuronal") | ws["coexpressing"]]
    lines.append("\nSPATIAL NULL (local neuron-neighbour density; spillover sits next to neurons):")
    if len(cand) and cand["local_neuron_density"].notna().any():
        cand_med = float(cand["local_neuron_density"].median())
        bg_med = float(mal_bg.median()) if len(mal_bg) else float("nan")
        lines.append(f"  hybrid candidates median neuron-neighbour frac: {cand_med:.3f}")
        lines.append(f"  malignant-cell background median:               {bg_med:.3f}")
        lines.append("  => " + ("candidates are ENRICHED for neuron neighbours vs background — "
                                "spillover cannot be excluded for much of the signal."
                                if np.isfinite(bg_med) and cand_med > 2 * max(bg_med, 1e-3)
                                else "candidates are NOT preferentially adjacent to neurons — "
                                "the neuronal co-expression is spatially dispersed / "
                                "cell-intrinsic, consistent with a real hybrid state."))
        fig, ax = plt.subplots(figsize=(7, 5))
        bins = np.linspace(0, max(0.05, float(np.nanmax(
            np.r_[cand["local_neuron_density"].to_numpy(), mal_bg.to_numpy()]))), 40)
        ax.hist(mal_bg.dropna(), bins=bins, density=True, histtype="step", lw=2,
                color="#999999", label=f"malignant background (n={mal_bg.notna().sum():,})")
        ax.hist(cand["local_neuron_density"].dropna(), bins=bins, density=True, histtype="step",
                lw=2, color="#d73027", label=f"hybrid candidates (n={cand['local_neuron_density'].notna().sum():,})")
        ax.set_xlabel(f"local {args.neuron_type}-neighbour fraction")
        ax.set_ylabel("density"); ax.legend(fontsize=8)
        ax.set_title("Spatial null: are hybrid candidates just sitting next to neurons?")
        fig.tight_layout()
        fig.savefig(args.output_dir / "spatial_null.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
    else:
        lines.append("  no hybrid candidates / neuron-density unavailable; skipped.")

    ws.to_csv(args.output_dir / "hybrid_continuum_cells.csv.gz", index=True)
    (args.output_dir / "HYBRID_SUMMARY.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nWrote tables + plots + HYBRID_SUMMARY.txt to {args.output_dir}")


if __name__ == "__main__":
    main()
