#!/usr/bin/env python3
"""Stage 5d (diagnostic): edge vs bulk — why does infiltrating-edge Low_signal show a
stronger malignant CNV signature than tumor-bulk Low_signal? Biology or a depth artifact?

Tumor bulk is more necrotic/hypoxic, so its cells can have degraded RNA and thus weaker
CNV inference regardless of true clonality. This rebuilds the per-cell malignant-signature
(cosine to the malignant CNV consensus, exactly as compare_insitucnv_groups.py) from the
per-section CNV outputs, joins per-cell RNA depth from the typed cohort, and for Low_signal
cells split by Region reports:
  - signature + depth distributions per region;
  - DEPTH-STRATIFIED signature: bulk vs edge WITHIN matched depth bins. If edge stays above
    bulk at equal depth, the edge>bulk gap is biological, not a depth/quality artifact;
  - per-donor consistency (is edge>bulk in most donors, or driven by a few?).

Reads (--cnv-dir): the per-section *_cnv.h5ad from run_insitucnv.py.
      (--typed-h5ad): cosmx_typed.h5ad, for per-cell total_counts / nFeature (obs-only).
Writes (--output-dir): edge_vs_bulk_summary.txt, sig_hist_by_region.png, sig_vs_depth.png.

Usage:
    python pipeline/python/insitucnv_edge_vs_bulk.py \\
        --cnv-dir persection --typed-h5ad cosmx_typed.h5ad \\
        --reference-file pipeline/reference/insitucnv_reference_types.txt --output-dir out
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
import scipy.sparse as sp

from compare_insitucnv_groups import CONTRALATERAL, DEFAULT_MALIGNANT, load_concat

BULK, EDGE = "Tumor bulk", "Infiltrating edge"
DEPTH_CANDIDATES = ["total_counts", "nCount_RNA", "nCount", "qc_gene_counts"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cnv-dir", type=Path, required=True)
    p.add_argument("--typed-h5ad", type=Path, required=True,
                   help="cosmx_typed.h5ad (obs-only read for per-cell RNA depth).")
    p.add_argument("--reference-file", type=Path, required=True,
                   help="(unused for the signature but kept for interface symmetry).")
    p.add_argument("--malignant-groups", default=",".join(DEFAULT_MALIGNANT))
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--n-depth-bins", type=int, default=4)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def malignant_signature(X, obs, malignant_types, celltype_key):
    """Per-cell cosine to the mean CNV profile of the confidently-malignant cells."""
    mal = np.flatnonzero(obs[celltype_key].astype(str).isin(malignant_types).to_numpy())
    if not mal.size:
        sys.exit("ERROR: no malignant-class cells to build the consensus.")
    centroid = np.asarray(X[mal].mean(axis=0)).ravel()
    cnorm = float(np.linalg.norm(centroid))
    if cnorm == 0:
        sys.exit("ERROR: malignant centroid is all-zero.")
    unit = centroid / cnorm
    dots = np.asarray(X @ unit).ravel()
    rn = np.sqrt(np.asarray(X.multiply(X).sum(axis=1)).ravel())
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(rn > 0, dots / rn, 0.0)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    malignant_types = {s.strip() for s in args.malignant_groups.split(",") if s.strip()}

    X, obs, _ = load_concat(args.cnv_dir, args.celltype_key, args.region_key)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    obs["mal_sig"] = malignant_signature(X, obs, malignant_types, args.celltype_key)

    # join per-cell RNA depth from the typed cohort (obs-only, backed)
    tobs = ad.read_h5ad(args.typed_h5ad, backed="r").obs
    depth_col = next((c for c in DEPTH_CANDIDATES if c in tobs.columns), None)
    if depth_col is None:
        sys.exit(f"ERROR: no depth column in typed obs (tried {DEPTH_CANDIDATES}).")
    obs["depth"] = pd.to_numeric(tobs[depth_col].reindex(obs.index), errors="coerce").to_numpy()
    print(f"Depth column: '{depth_col}'; joined for {int(obs['depth'].notna().sum()):,} cells.")

    # Low_signal cells by region
    ls = obs[obs[args.celltype_key].astype(str) == args.lowsignal_label].copy()
    ls["region"] = ls[args.region_key].astype(str)
    ls["donor"] = ls["tissue_section"].astype(str).str.split("__").str[1]

    lines = ["Edge vs bulk — Low_signal malignant-signature by region\n"]
    lines.append(f"{'region':<26} {'n':>10} {'sig_median':>11} {'sig_q75':>9} "
                 f"{'depth_median':>13}")
    for r in (BULK, EDGE, CONTRALATERAL):
        sub = ls[ls["region"] == r]
        if not len(sub):
            continue
        lines.append(f"{r:<26} {len(sub):>10,} {sub['mal_sig'].median():>11.3f} "
                     f"{sub['mal_sig'].quantile(0.75):>9.3f} {sub['depth'].median():>13.0f}")

    # ---- DEPTH-STRATIFIED: bulk vs edge within matched depth bins --------------------
    be = ls[ls["region"].isin([BULK, EDGE])].dropna(subset=["depth"])
    qs = np.quantile(be["depth"], np.linspace(0, 1, args.n_depth_bins + 1))
    qs[0], qs[-1] = -np.inf, np.inf
    be["depth_bin"] = pd.cut(be["depth"], bins=np.unique(qs), include_lowest=True)
    lines.append("\nDEPTH-STRATIFIED median signature (bulk vs edge at matched depth):")
    lines.append(f"  {'depth bin (counts)':<26} {'bulk_med':>9} {'n_bulk':>9} "
                 f"{'edge_med':>9} {'n_edge':>9}  edge>bulk?")
    edge_wins = 0
    bins = [b for b in be["depth_bin"].cat.categories]
    for b in bins:
        bb = be[(be["depth_bin"] == b) & (be["region"] == BULK)]["mal_sig"]
        ee = be[(be["depth_bin"] == b) & (be["region"] == EDGE)]["mal_sig"]
        if not len(bb) or not len(ee):
            continue
        win = ee.median() > bb.median()
        edge_wins += int(win)
        lines.append(f"  {str(b):<26} {bb.median():>9.3f} {len(bb):>9,} "
                     f"{ee.median():>9.3f} {len(ee):>9,}  {'YES' if win else 'no'}")
    verdict = ("edge > bulk at EVERY matched depth bin -> the gap is BIOLOGICAL, not a "
               "depth artifact." if edge_wins == len([b for b in bins]) else
               f"edge > bulk in {edge_wins}/{len(bins)} depth bins -> partly depth-driven; "
               "inspect the bins where bulk catches up.")
    lines.append(f"  => {verdict}")

    # ---- per-donor consistency -------------------------------------------------------
    lines.append("\nPer-donor Low_signal median signature (bulk vs edge):")
    lines.append(f"  {'donor':<8} {'bulk_med':>9} {'edge_med':>9}  edge>bulk?")
    consistent = 0
    donors = sorted(ls["donor"].dropna().unique())
    for d in donors:
        bb = ls[(ls["donor"] == d) & (ls["region"] == BULK)]["mal_sig"]
        ee = ls[(ls["donor"] == d) & (ls["region"] == EDGE)]["mal_sig"]
        if len(bb) < 50 or len(ee) < 50:
            continue
        win = ee.median() > bb.median()
        consistent += int(win)
        lines.append(f"  {d:<8} {bb.median():>9.3f} {ee.median():>9.3f}  "
                     f"{'YES' if win else 'no'}")
    lines.append(f"  => edge > bulk in {consistent} donors with both regions "
                 f"(≥50 cells each).")

    (args.output_dir / "edge_vs_bulk_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))

    # ---- plots -----------------------------------------------------------------------
    colors = {BULK: "#d95f0e", EDGE: "#7570b3", CONTRALATERAL: "#1b9e77"}
    fig, ax = plt.subplots(figsize=(8, 5))
    bins_h = np.linspace(-0.4, 1.0, 60)
    for r in (CONTRALATERAL, BULK, EDGE):
        sub = ls[ls["region"] == r]["mal_sig"].dropna()
        if len(sub):
            ax.hist(sub, bins=bins_h, density=True, histtype="step", lw=2,
                    color=colors[r], label=f"{r} (n={len(sub):,})")
    ax.set_xlabel("malignant-signature (cosine to malignant CNV consensus)")
    ax.set_ylabel("density"); ax.legend(fontsize=8)
    ax.set_title("Low_signal malignant-signature by region")
    fig.tight_layout()
    fig.savefig(args.output_dir / "sig_hist_by_region.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    centers = [be[be["depth_bin"] == b]["depth"].median() for b in bins]
    for r in (BULK, EDGE):
        meds = [be[(be["depth_bin"] == b) & (be["region"] == r)]["mal_sig"].median()
                for b in bins]
        ax.plot(centers, meds, "o-", color=colors[r], label=r)
    ax.set_xlabel(f"RNA depth ({depth_col}, bin median)")
    ax.set_ylabel("median malignant-signature")
    ax.legend(fontsize=9); ax.set_title("Signature vs depth: edge vs bulk (matched depth)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "sig_vs_depth.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote summary + plots to {args.output_dir}")


if __name__ == "__main__":
    main()
