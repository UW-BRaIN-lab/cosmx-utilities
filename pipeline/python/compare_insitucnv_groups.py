#!/usr/bin/env python3
"""Stage 5c: aggregate per-section InSituCNV results and test whether Low_signal is tumor.

Concatenates the per-tissue-section CNV h5ads (obsm['X_cnv']), forms per-GROUP mean CNV
profiles, and asks the core question: do the Low_signal cells carry the GBM malignant
copy-number signature (chr7 gain / chr10 loss / chr9 loss), clustering with the known
tumor states rather than the diploid reference?

Groups: obs['cell_type'], except Low_signal is split by Region ("Low_signal | Tumor bulk"
etc.) so its tumour-region vs contralateral contrast is visible. Each group is labelled a
CLASS: reference (diploid baseline), malignant (positive control), low_signal (the test),
or other. Per-group CNV is noisy on a targeted panel, so the readout is aggregate (Moldia's
own guidance): mean profiles, their cosine similarity, chromosome-arm summaries, and a
per-cell cnv_score compared against control-derived expectations.

Full methods — design choices, the malignant-signature metric, thresholds, field-effect
controls, and limitations — are documented in pipeline/INSITUCNV_METHODS.md.

Reads (--cnv-dir): all *_cnv.h5ad from run_insitucnv.py (shared X_cnv windows + uns['cnv']).
Writes (--output-dir):
  group_mean_cnv.csv          groups x genomic windows (mean CNV)
  cosine_similarity.csv       group x group cosine similarity of mean profiles
  chr_arm_summary.csv         per-group mean CNV per chromosome (+ class)
  cnv_score_by_group.csv      per-group cnv_score quantiles + % above malignant threshold
  chromosome_heatmap.png      groups x windows heatmap, chromosomes delimited, class colour-bar
  cosine_similarity.png       clustered similarity heatmap
  chr7_chr10.png              chr7 (gain) vs chr10 (loss) mean CNV per group
  cnv_score.png               cnv_score distribution per class
  SUMMARY.txt                 the plain-language verdict + calibration numbers

Usage:
    python pipeline/python/compare_insitucnv_groups.py \\
        --cnv-dir persection --reference-file pipeline/reference/insitucnv_reference_types.txt \\
        --output-dir out --min-cells 200
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
from sklearn.metrics.pairwise import cosine_similarity

# The 14 confidently-malignant InSituTree states (GBmap Neftel + our de-novo tumour).
# These are the POSITIVE CONTROL — they must show chr7 gain + chr10 loss.
DEFAULT_MALIGNANT = [
    "AC-like", "AC-like_Prolif", "MES-like_hypoxia_MHC", "MES-like_hypoxia_independent",
    "NPC-like_OPC", "NPC-like_Prolif", "NPC-like_neural", "OPC-like", "OPC-like_Prolif",
    "Hypoxia_denovo", "MES_AClike_denovo", "MESlike_denovo", "OPClike_denovo", "Stress_denovo",
]
CONTRALATERAL = "Contralateral uninvolved"  # obs['Region'] value marking uninvolved brain
TUMOR_BULK = "Tumor bulk"
INFILTRATING_EDGE = "Infiltrating edge"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cnv-dir", type=Path, required=True,
                   help="Directory of *_cnv.h5ad from run_insitucnv.py.")
    p.add_argument("--reference-file", type=Path, required=True)
    p.add_argument("--malignant-groups", default=",".join(DEFAULT_MALIGNANT),
                   help="Comma-separated malignant cell types (positive control).")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--min-cells", type=int, default=200,
                   help="Drop groups smaller than this from the profile comparison.")
    p.add_argument("--donor-threshold", type=float, default=None,
                   help="Write an ADDITIONAL per-donor heatmap at this malignant cutoff "
                        "(e.g. the bimodal trough ~0.45) to *_thr<val> files, alongside the "
                        "strict sig_thr version. Omit for strict-only.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_reference_types(path: Path) -> list[str]:
    return [ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()
            if ln.split("#", 1)[0].strip()]


def malignant_signature(X, obs, malignant_types, celltype_key="cell_type"):
    """Per-cell cosine similarity to the mean CNV profile of the confidently-malignant
    cells (the directional chr-loss discriminator; expression-based CNV on a targeted panel
    detects losses far better than gains, and the non-directional L2 burden barely separates
    tumour from normal). Returns a float array aligned to X's rows, or ``None`` when there is
    no usable malignant consensus (no malignant cells, or an all-zero centroid)."""
    mal = np.flatnonzero(obs[celltype_key].astype(str).isin(malignant_types).to_numpy())
    if not mal.size:
        return None
    centroid = np.asarray(X[mal].mean(axis=0)).ravel()
    cnorm = float(np.linalg.norm(centroid))
    if cnorm == 0:
        return None
    unit = centroid / cnorm
    dots = np.asarray(X @ unit).ravel()
    xsq = X.multiply(X) if hasattr(X, "multiply") else np.square(X)
    rn = np.sqrt(np.asarray(xsq.sum(axis=1)).ravel())
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(rn > 0, dots / rn, 0.0)


def load_concat(cnv_dir: Path, celltype_key: str, region_key: str, with_spatial: bool = False):
    """Concatenate per-section CNV results into arrays (no full AnnData needed downstream).

    Carries the small set of obs columns downstream code needs (typing/region/donor/block/
    slide + cnv_score + tissue_section). When ``with_spatial`` and the section h5ads have
    ``obsm['spatial']``, the per-cell centroid is added to obs as ``spatial_x``/``spatial_y``
    (per-slide global-px, so only meaningful WITHIN a tissue_section)."""
    files = sorted(cnv_dir.glob("*_cnv.h5ad"))
    if not files:
        sys.exit(f"ERROR: no *_cnv.h5ad in {cnv_dir}")
    xs, obs_parts, chr_pos, width = [], [], None, None
    for f in files:
        a = ad.read_h5ad(f)
        if "X_cnv" not in a.obsm:
            print(f"  WARN: {f.name} has no X_cnv; skipped", file=sys.stderr)
            continue
        x = a.obsm["X_cnv"]
        x = sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr()
        if width is None:
            width = x.shape[1]
            chr_pos = a.uns.get("cnv", {}).get("chr_pos")
        elif x.shape[1] != width:
            print(f"  WARN: {f.name} X_cnv width {x.shape[1]} != {width}; skipped",
                  file=sys.stderr)
            continue
        xs.append(x)
        cols = {c: a.obs[c].astype(str).to_numpy() for c in
                (celltype_key, region_key, "cnv_score", "tissue_section",
                 "Case", "Block", "slide_id") if c in a.obs}
        part = pd.DataFrame(cols, index=a.obs.index)
        if with_spatial and "spatial" in a.obsm:
            xy = np.asarray(a.obsm["spatial"], dtype=np.float64)
            part["spatial_x"], part["spatial_y"] = xy[:, 0], xy[:, 1]
        obs_parts.append(part)
        print(f"  {f.name}: {a.n_obs:,} cells")
    if not xs:
        sys.exit("ERROR: no usable CNV sections.")
    X = sp.vstack(xs, format="csr")
    obs = pd.concat(obs_parts, axis=0)
    obs["cnv_score"] = pd.to_numeric(obs["cnv_score"], errors="coerce")
    print(f"Concatenated: {X.shape[0]:,} cells x {X.shape[1]} windows from {len(xs)} sections.")
    return X, obs, chr_pos


def window_chromosomes(chr_pos: dict, n_windows: int) -> np.ndarray:
    """Map each X_cnv column to its chromosome using uns['cnv']['chr_pos'] (chrom->start col)."""
    if not chr_pos:
        return np.array(["?"] * n_windows)
    items = sorted(chr_pos.items(), key=lambda kv: kv[1])
    labels = np.empty(n_windows, dtype=object)
    for i, (chrom, start) in enumerate(items):
        end = items[i + 1][1] if i + 1 < len(items) else n_windows
        labels[start:end] = chrom
    return labels


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference_types = set(read_reference_types(args.reference_file))
    malignant_types = set(s.strip() for s in args.malignant_groups.split(",") if s.strip())

    X, obs, chr_pos = load_concat(args.cnv_dir, args.celltype_key, args.region_key)
    win_chrom = window_chromosomes(chr_pos, X.shape[1])

    # --- group + class labels -------------------------------------------------------
    ct = obs[args.celltype_key].to_numpy()
    region = obs[args.region_key].to_numpy() if args.region_key in obs else np.array([""] * len(obs))
    is_ls = ct == args.lowsignal_label
    group = np.where(is_ls, [f"{args.lowsignal_label} | {r}" for r in region], ct)
    obs["group"] = group

    def classify(cell_type: str) -> str:
        if cell_type == args.lowsignal_label:
            return "low_signal"
        if cell_type in malignant_types:
            return "malignant"
        if cell_type in reference_types:
            return "reference"
        return "other"
    obs["class"] = [classify(c) for c in ct]
    group_class = obs.groupby("group", observed=True)["class"].first()

    # --- per-group mean CNV profiles (drop small groups) ----------------------------
    counts = obs["group"].value_counts()
    keep_groups = counts[counts >= args.min_cells].index.tolist()
    grp_arr = obs["group"].to_numpy()
    means, order = [], []
    for g in keep_groups:
        rows = np.flatnonzero(grp_arr == g)
        means.append(np.asarray(X[rows].mean(axis=0)).ravel())
        order.append(g)
    M = np.vstack(means)  # groups x windows
    prof = pd.DataFrame(M, index=order)
    prof.to_csv(args.output_dir / "group_mean_cnv.csv")

    # --- cosine similarity between group-mean profiles ------------------------------
    sim = pd.DataFrame(cosine_similarity(M), index=order, columns=order)
    sim.to_csv(args.output_dir / "cosine_similarity.csv")

    # --- per-chromosome arm summary -------------------------------------------------
    chrom_order = sorted(set(win_chrom), key=lambda c: int(c.replace("chr", "")) if
                         c.replace("chr", "").isdigit() else 99)
    arm = pd.DataFrame({c: prof.loc[:, win_chrom == c].mean(axis=1) for c in chrom_order})
    arm.insert(0, "class", group_class.reindex(order).to_numpy())
    arm.insert(1, "n_cells", counts.reindex(order).to_numpy())
    arm.to_csv(args.output_dir / "chr_arm_summary.csv")

    # --- panel gene coverage per chromosome vs malignant signal ---------------------
    # Reviewer/PI point: losses can concentrate in gene-poor regions, so a targeted panel
    # might under-sample them. Cross the panel gene count per chromosome (from a section's
    # var; all sections share the gene set) against the malignant-consensus mean CNV per
    # chromosome, to show the signal-carrying arms (chr7 gain, chr10/9/14 loss) are well
    # covered — i.e. no arm's signal (or its absence) is a gene-density/coverage artifact.
    try:
        first_var = ad.read_h5ad(sorted(Path(args.cnv_dir).glob("*_cnv.h5ad"))[0],
                                 backed="r").var
    except Exception:
        first_var = None
    mal_rows = np.flatnonzero(obs[args.celltype_key].astype(str).isin(malignant_types).to_numpy())
    if first_var is not None and "chromosome" in first_var and mal_rows.size:
        genes_per_chr = first_var["chromosome"].astype(str).value_counts()
        centroid = np.asarray(X[mal_rows].mean(axis=0)).ravel()
        cov = pd.DataFrame([
            {"chromosome": c, "n_panel_genes": int(genes_per_chr.get(c, 0)),
             "n_windows": int((win_chrom == c).sum()),
             "malignant_mean_cnv": (float(centroid[win_chrom == c].mean())
                                    if (win_chrom == c).any() else np.nan)}
            for c in chrom_order])
        cov.to_csv(args.output_dir / "chr_coverage_vs_signal.csv", index=False)
        fig, ax = plt.subplots(figsize=(9, 6))
        yy = np.arange(len(cov)); vals = cov["malignant_mean_cnv"].to_numpy()
        ax.barh(yy, vals, color=["#d62728" if v > 0 else "#1f77b4" for v in vals])
        for i, r in cov.iterrows():
            v = r["malignant_mean_cnv"]
            ax.text((0.0006 if v >= 0 else -0.0006), i, f"{int(r['n_panel_genes'])} genes",
                    va="center", ha="left" if v >= 0 else "right", fontsize=7)
        ax.set_yticks(yy); ax.set_yticklabels([c.replace("chr", "chr ") for c in cov["chromosome"]],
                                              fontsize=8)
        ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("malignant-consensus mean CNV  (red = gain, blue = loss)")
        ax.set_title("Per-chromosome malignant signal vs panel gene coverage\n"
                     "(panel gene count labelled on each bar)", fontsize=10)
        fig.tight_layout()
        fig.savefig(args.output_dir / "chr_coverage_vs_signal.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # --- cnv_score: threshold from controls, then score every group -----------------
    neg_mask = (obs["class"] == "reference").to_numpy() | \
               ((obs[args.celltype_key] == args.lowsignal_label).to_numpy() &
                (region == CONTRALATERAL))
    pos_mask = (obs["class"] == "malignant").to_numpy()
    neg_scores = obs.loc[neg_mask, "cnv_score"].dropna()
    threshold = float(np.percentile(neg_scores, 95)) if len(neg_scores) else float("nan")
    sc_summary = obs.groupby("group", observed=True)["cnv_score"].agg(
        n="size", median="median", q75=lambda s: s.quantile(0.75))
    sc_summary["frac_above_malignant_thr"] = obs.groupby("group", observed=True)["cnv_score"] \
        .apply(lambda s: float((s > threshold).mean()))
    sc_summary["class"] = group_class.reindex(sc_summary.index).to_numpy()
    sc_summary = sc_summary.sort_values("median", ascending=False)
    sc_summary.to_csv(args.output_dir / "cnv_score_by_group.csv")

    # --- malignant CNV-SIGNATURE score (directional — the real discriminator) -------
    # Expression-based CNV on a targeted panel detects LOSSES far better than gains, and the
    # L2 burden above is non-directional and noise-dominated, so it barely separates tumour
    # from normal. Score each cell instead by how well its CNV profile matches the malignant
    # CONSENSUS direction (the chr10/arm-loss pattern): cosine similarity to the mean profile
    # of the confidently-malignant cells. Calibrated on the negative controls.
    ls_contra = ((obs[args.celltype_key] == args.lowsignal_label).to_numpy()
                 & (region == CONTRALATERAL))
    mal_rows = np.flatnonzero((obs["class"] == "malignant").to_numpy())
    sig_thr, sig_summary = float("nan"), None
    if mal_rows.size:
        sig = malignant_signature(X, obs, malignant_types, args.celltype_key)
        if sig is not None:
            obs["mal_sig"] = sig
            neg_sig = obs.loc[neg_mask, "mal_sig"].dropna()
            sig_thr = float(np.percentile(neg_sig, 95)) if len(neg_sig) else float("nan")
            sig_summary = obs.groupby("group", observed=True)["mal_sig"].agg(
                n="size", median="median")
            sig_summary["frac_above"] = obs.groupby("group", observed=True)["mal_sig"] \
                .apply(lambda s: float((s > sig_thr).mean()))
            sig_summary["class"] = group_class.reindex(sig_summary.index).to_numpy()
            sig_summary = sig_summary.sort_values("median", ascending=False)
            sig_summary.to_csv(args.output_dir / "mal_signature_by_group.csv")

    # ================================ PLOTS ========================================
    class_rank = {"reference": 0, "low_signal": 1, "malignant": 2, "other": 3}
    class_color = {"reference": "#2c7fb8", "low_signal": "#d95f0e",
                   "malignant": "#c51b8a", "other": "#999999"}
    plot_order = sorted(order, key=lambda g: (class_rank.get(group_class.get(g, "other"), 3), g))
    pidx = [order.index(g) for g in plot_order]

    # chromosome heatmap (groups x windows)
    fig, ax = plt.subplots(figsize=(13, max(4, 0.28 * len(plot_order) + 2)))
    im = ax.imshow(M[pidx], aspect="auto", cmap="bwr", vmin=-0.4, vmax=0.4,
                   interpolation="nearest")
    ax.set_yticks(range(len(plot_order)))
    ax.set_yticklabels(plot_order, fontsize=6)
    # chromosome boundaries + centered labels
    bounds = [chr_pos[c] for c in chrom_order if c in (chr_pos or {})]
    for b in bounds:
        ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.4)
    centers = [(chr_pos[c] + (chr_pos[chrom_order[i + 1]] if i + 1 < len(chrom_order) else M.shape[1])) / 2
               for i, c in enumerate(chrom_order) if c in (chr_pos or {})]
    ax.set_xticks(centers)
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=6)
    ax.set_xlabel("chromosome (genomic windows)")
    for tick, g in zip(ax.get_yticklabels(), plot_order):
        tick.set_color(class_color.get(group_class.get(g, "other"), "#333"))
    ax.set_title("InSituCNV mean profile per group (blue=loss, red=gain; "
                 "GBM tumour = chr7 gain + chr10 loss)")
    fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01, label="mean CNV")
    fig.tight_layout()
    fig.savefig(args.output_dir / "chromosome_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # cosine similarity heatmap
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(plot_order) + 2),
                                    max(5, 0.3 * len(plot_order) + 2)))
    S = sim.loc[plot_order, plot_order].to_numpy()
    im = ax.imshow(S, cmap="viridis", vmin=-1, vmax=1)
    ax.set_xticks(range(len(plot_order))); ax.set_yticks(range(len(plot_order)))
    ax.set_xticklabels(plot_order, fontsize=6, rotation=90)
    ax.set_yticklabels(plot_order, fontsize=6)
    for ticks in (ax.get_xticklabels(), ax.get_yticklabels()):
        for tick, g in zip(ticks, plot_order):
            tick.set_color(class_color.get(group_class.get(g, "other"), "#333"))
    ax.set_title("Cosine similarity of mean CNV profiles")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.tight_layout()
    fig.savefig(args.output_dir / "cosine_similarity.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # chr7 vs chr10
    if "chr7" in arm.columns and "chr10" in arm.columns:
        fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * len(plot_order) + 2)))
        y = np.arange(len(plot_order))
        ax.barh(y - 0.2, arm.loc[plot_order, "chr7"], height=0.4, color="#d62728", label="chr7 (gain)")
        ax.barh(y + 0.2, arm.loc[plot_order, "chr10"], height=0.4, color="#1f77b4", label="chr10 (loss)")
        ax.set_yticks(y); ax.set_yticklabels(plot_order, fontsize=6)
        for tick, g in zip(ax.get_yticklabels(), plot_order):
            tick.set_color(class_color.get(group_class.get(g, "other"), "#333"))
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("mean CNV"); ax.legend(fontsize=8)
        ax.set_title("GBM signature: chr7 gain (+) & chr10 loss (-) per group")
        fig.tight_layout()
        fig.savefig(args.output_dir / "chr7_chr10.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # cnv_score by class
    fig, ax = plt.subplots(figsize=(7, 5))
    classes = ["reference", "low_signal", "malignant", "other"]
    data = [obs.loc[obs["class"] == c, "cnv_score"].dropna().to_numpy() for c in classes]
    data = [(d if len(d) else np.array([np.nan])) for d in data]
    bp = ax.boxplot(data, vert=True, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black"))
    for c, box in zip(classes, bp["boxes"]):
        box.set_facecolor(class_color[c])
    ax.set_xticklabels(classes, fontsize=9)
    if np.isfinite(threshold):
        ax.axhline(threshold, color="red", ls="--", lw=1,
                   label=f"malignant thr (ref 95th pct = {threshold:.3g})")
        ax.legend(fontsize=8)
    ax.set_ylabel("per-cell cnv_score (L2 CNV burden)")
    ax.set_title("CNV burden by class")
    fig.tight_layout()
    fig.savefig(args.output_dir / "cnv_score.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # Low_signal rate + magnitude by region (typing-level; where the flat fraction comes
    # from). Bar height = total cells in the region, split typed vs Low_signal.
    region_order = [r for r in (CONTRALATERAL, INFILTRATING_EDGE, TUMOR_BULK) if r in set(region)]
    if region_order and args.lowsignal_label in set(ct):
        is_ls_all = ct == args.lowsignal_label
        rc = pd.DataFrame({
            "region": region_order,
            "total_cells": [int((region == r).sum()) for r in region_order],
            "low_signal": [int((is_ls_all & (region == r)).sum()) for r in region_order],
        })
        rc["typed"] = rc["total_cells"] - rc["low_signal"]
        rc["ls_rate"] = rc["low_signal"] / rc["total_cells"]
        tot_ls = int(rc["low_signal"].sum())
        rc["pct_of_all_lowsignal"] = rc["low_signal"] / tot_ls if tot_ls else 0.0
        rc.to_csv(args.output_dir / "lowsignal_by_region.csv", index=False)

        fig, ax = plt.subplots(figsize=(7, 5))
        x = np.arange(len(region_order))
        ax.bar(x, rc["typed"], color="#bdbdbd", label="typed")
        ax.bar(x, rc["low_signal"], bottom=rc["typed"], color="#d95f0e", label="Low_signal")
        for i in range(len(region_order)):
            ax.text(i, rc["total_cells"][i], f"{rc['ls_rate'][i]:.0%} LS\n"
                    f"({int(rc['low_signal'][i]):,} / {int(rc['pct_of_all_lowsignal'][i] * 100)}% of LS)",
                    ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(region_order, rotation=15, fontsize=9)
        ax.set_ylabel("cells"); ax.legend(loc="upper left")
        ax.set_ylim(top=rc["total_cells"].max() * 1.18)
        ax.set_title("Low_signal rate & magnitude by region")
        fig.tight_layout()
        fig.savefig(args.output_dir / "lowsignal_by_region.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # within-region signature: reference vs Low_signal vs malignant, per region. Both
    # reference and Low_signal cells here are equally exposed to the tissue field effect, so
    # Low_signal sitting ABOVE reference within a tumour region is tumour-specific signal.
    if sig_summary is not None:
        cls_arr = obs["class"].to_numpy()
        is_ls = ct == args.lowsignal_label
        sig_all = obs["mal_sig"].to_numpy()
        panel = [(CONTRALATERAL, "contralateral"), (INFILTRATING_EDGE, "infiltrating edge"),
                 (TUMOR_BULK, "tumor bulk")]
        fig, axes = plt.subplots(1, len(panel), figsize=(4 * len(panel), 5), sharey=True)
        for ax, (r, title) in zip(np.atleast_1d(axes), panel):
            in_r = region == r
            groups_pc = [("reference", cls_arr == "reference"), ("Low_signal", is_ls),
                         ("malignant", cls_arr == "malignant")]
            data = []
            for _, m in groups_pc:
                v = sig_all[in_r & m]
                v = v[np.isfinite(v)]
                data.append(v if v.size else np.array([np.nan]))
            bp = ax.boxplot(data, showfliers=False, patch_artist=True,
                            medianprops=dict(color="black"))
            for box, col in zip(bp["boxes"], ["#2c7fb8", "#d95f0e", "#c51b8a"]):
                box.set_facecolor(col)
            ax.set_xticklabels([g for g, _ in groups_pc], fontsize=8, rotation=20)
            ax.set_title(title, fontsize=10)
            # field-effect floor for THIS region = median reference (diploid) signature;
            # a Low_signal cell only reads as tumour if it clears this floor, not just >0.
            ref_v = sig_all[in_r & (cls_arr == "reference")]
            ref_v = ref_v[np.isfinite(ref_v)]
            if ref_v.size:
                ax.axhline(float(np.median(ref_v)), color="#2c7fb8", ls=":", lw=1.6,
                           label="field-effect floor (ref median)")
            # malignant call threshold(s): strict (sig_thr) + optional sensitive (trough)
            ax.axhline(sig_thr, color="#c51b8a", ls="--", lw=1.1,
                       label=f"malignant thr {sig_thr:.2f}")
            if args.donor_threshold is not None:
                ax.axhline(args.donor_threshold, color="#d95f0e", ls="--", lw=1.1,
                           label=f"sensitive thr {args.donor_threshold:.2f}")
        np.atleast_1d(axes)[0].set_ylabel("malignant-signature (cosine)")
        np.atleast_1d(axes)[-1].legend(fontsize=7, loc="upper right")
        fig.suptitle("Within-region signature vs field-effect floor & malignant threshold\n"
                     "(Low_signal is tumour only where it clears the region's reference floor)",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(args.output_dir / "within_region_signature.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # Low_signal CNV resolution by region: splits each region into typed / Low_signal-CNV-
    # normal / Low_signal-CNV-malignant — layers the malignant call onto the region typing so
    # you see how much of each region is transcriptionally-flat-but-CNV-malignant tumor.
    if sig_summary is not None:
        r_order = [r for r in (CONTRALATERAL, INFILTRATING_EDGE, TUMOR_BULK) if r in set(region)]
        is_ls = ct == args.lowsignal_label
        sig_all = obs["mal_sig"].to_numpy()

        def region_cnv_breakdown(thr, suffix):
            rows = []
            for r in r_order:
                in_r = region == r
                n_tot = int(in_r.sum())
                ls_r = in_r & is_ls
                sig_r = sig_all[ls_r]
                sig_r = sig_r[np.isfinite(sig_r)]
                n_ls, n_mal = int(ls_r.sum()), int((sig_r > thr).sum())
                rows.append(dict(region=r, total=n_tot, typed=n_tot - n_ls,
                                 ls_normal=n_ls - n_mal, ls_malignant=n_mal,
                                 malignant_pct_of_region=(n_mal / n_tot if n_tot else 0.0),
                                 malignant_pct_of_lowsignal=(n_mal / n_ls if n_ls else 0.0)))
            rcn = pd.DataFrame(rows)
            rcn.to_csv(args.output_dir / f"lowsignal_cnv_by_region{suffix}.csv", index=False)
            fig, ax = plt.subplots(figsize=(7, 5))
            x = np.arange(len(r_order))
            typed, lsn, lsm = (rcn[c].to_numpy() for c in ("typed", "ls_normal", "ls_malignant"))
            ax.bar(x, typed, color="#bdbdbd", label="typed")
            ax.bar(x, lsn, bottom=typed, color="#4575b4", label="Low_signal · CNV-normal")
            ax.bar(x, lsm, bottom=typed + lsn, color="#d73027",
                   label="Low_signal · CNV-malignant")
            for i in range(len(r_order)):
                ax.text(i, rcn["total"][i], f"{rcn['malignant_pct_of_region'][i]:.0%} malig\n"
                        f"({int(lsm[i]):,})", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(x); ax.set_xticklabels(r_order, rotation=15, fontsize=9)
            ax.set_ylabel("cells"); ax.legend(loc="upper left")
            ax.set_ylim(top=rcn["total"].max() * 1.18)
            ax.set_title(f"Low_signal CNV resolution by region (malig thr {thr:.3g})")
            fig.tight_layout()
            fig.savefig(args.output_dir / f"lowsignal_cnv_by_region{suffix}.png", dpi=180,
                        bbox_inches="tight")
            plt.close(fig)

        if r_order:
            region_cnv_breakdown(sig_thr, "")           # strict / conservative call
            if args.donor_threshold is not None:        # sensitive (e.g. bimodal trough)
                region_cnv_breakdown(args.donor_threshold, f"_thr{args.donor_threshold:g}")

    # Field-effect check: median malignant-signature of each REFERENCE (non-malignant) cell
    # type, split by region. Real CNV is a genotype property and can't depend on location; a
    # field effect (ambient tumour RNA pulled in by spatial smoothing) does. So a reference
    # type elevated in tumour regions but FLAT in contralateral proves its apparent CNV is
    # neighbourhood contamination, not copy number — the direct rebuttal to "normal cells
    # show CNV". Quiescent, dispersed normals stay flat everywhere; tumour-resident TME
    # subsets (SMC_prolif, Mono_hypoxia, TAM-hypoxia, Tip-like) light up only in tumour.
    if sig_summary is not None:
        sig_all = obs["mal_sig"].to_numpy()
        rcols = [r for r in (CONTRALATERAL, INFILTRATING_EDGE, TUMOR_BULK) if r in set(region)]
        rows = []
        for t in sorted(reference_types):
            mt = ct == t
            if int(mt.sum()) < args.min_cells:
                continue
            rec = {"cell_type": t}
            for r in rcols:
                v = sig_all[mt & (region == r)]
                v = v[np.isfinite(v)]
                rec[r] = float(np.median(v)) if v.size >= 50 else np.nan
                rec[f"n_{r}"] = int(v.size)
            rows.append(rec)
        if rows:
            nb = pd.DataFrame(rows)
            nb["_srt"] = nb[[c for c in (INFILTRATING_EDGE, TUMOR_BULK) if c in nb]].max(axis=1)
            nb = nb.sort_values("_srt", ascending=False).drop(columns="_srt")
            nb.to_csv(args.output_dir / "normal_types_by_region.csv", index=False)
            H = nb[rcols].to_numpy(dtype=float)
            fig, ax = plt.subplots(figsize=(6.5, max(4, 0.34 * len(nb) + 1.5)))
            im = ax.imshow(H, aspect="auto", cmap="Reds", vmin=0, vmax=0.6)
            ax.set_xticks(range(len(rcols))); ax.set_xticklabels(rcols, rotation=20, fontsize=8)
            ax.set_yticks(range(len(nb))); ax.set_yticklabels(nb["cell_type"], fontsize=7)
            for i in range(H.shape[0]):
                for j in range(H.shape[1]):
                    if np.isfinite(H[i, j]):
                        ax.text(j, i, f"{H[i, j]:+.2f}", ha="center", va="center", fontsize=6,
                                color="white" if H[i, j] > 0.35 else "black")
            ax.set_title("Reference (non-malignant) cell types: malignant-signature by region\n"
                         "(flat in contralateral, elevated only in tumor = field effect, not CNV)",
                         fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="median malignant-signature")
            fig.tight_layout()
            fig.savefig(args.output_dir / "normal_types_by_region.png", dpi=180,
                        bbox_inches="tight")
            plt.close(fig)

    # ============================== VERDICT ========================================
    lines = []
    lines.append(f"InSituCNV group comparison — {X.shape[0]:,} cells, {len(order)} groups "
                 f">= {args.min_cells} cells.\n")
    lsg = [x for x in order if x.startswith(args.lowsignal_label)]

    # PRIMARY: directional malignant-signature (cosine to the malignant CNV consensus).
    if sig_summary is not None:
        def fs(mask):
            s = obs.loc[mask, "mal_sig"].dropna()
            return float((s > sig_thr).mean()) if len(s) else float("nan")
        lines.append("MALIGNANT-SIGNATURE = per-cell cosine to the malignant CNV consensus "
                     "(the chr-loss pattern; the right discriminator on a targeted panel).")
        lines.append(f"  threshold (95th pct of diploid+contralateral controls) = {sig_thr:.4g}")
        lines.append("  CALIBRATION (should hold if the run is trustworthy):")
        lines.append(f"    malignant (positive control): {fs(pos_mask):.1%} above  (expect HIGH)")
        lines.append(f"    diploid reference:            "
                     f"{fs((obs['class'] == 'reference').to_numpy()):.1%}  (expect ~5%)")
        lines.append(f"    Low_signal | {CONTRALATERAL}:  {fs(ls_contra):.1%}  (expect LOW)")
        lines.append("  THE TEST — Low_signal signature-high fraction by region:")
        for g in lsg:
            fa = sig_summary.loc[g, "frac_above"] if g in sig_summary.index else float("nan")
            md = sig_summary.loc[g, "median"] if g in sig_summary.index else float("nan")
            lines.append(f"    {g:<34} sig-high={fa:.1%}  median={md:+.3f}  (n={int(counts[g]):,})")
        lines.append("")

        # WITHIN-REGION contrast: Low_signal vs REFERENCE in the SAME region. Both are equally
        # exposed to the tissue field effect (ambient tumour RNA bleeding into reference cells
        # inside the tumour), so the fraction of Low_signal above the SAME-region reference
        # 95th pct is a FIELD-CORRECTED malignant estimate; malignant-in-region is the ceiling.
        cls_arr = obs["class"].to_numpy()
        sig_arr = obs["mal_sig"].to_numpy()
        is_ls = ct == args.lowsignal_label
        wr_rows = []
        for r in (TUMOR_BULK, INFILTRATING_EDGE, CONTRALATERAL):
            in_r = region == r
            ref_r = sig_arr[in_r & (cls_arr == "reference")]
            ls_r = sig_arr[in_r & is_ls]
            mal_r = sig_arr[in_r & (cls_arr == "malignant")]
            ref_r, ls_r, mal_r = (v[np.isfinite(v)] for v in (ref_r, ls_r, mal_r))
            if not ref_r.size or not ls_r.size:
                continue
            thr_r = float(np.percentile(ref_r, 95))
            wr_rows.append(dict(
                region=r, n_ref=int(ref_r.size), n_lowsignal=int(ls_r.size),
                n_malignant=int(mal_r.size),
                ref_median=float(np.median(ref_r)), lowsignal_median=float(np.median(ls_r)),
                malignant_median=float(np.median(mal_r)) if mal_r.size else float("nan"),
                same_region_thr=thr_r,
                lowsignal_above_ref=float(np.mean(ls_r > thr_r)),
                malignant_above_ref=float(np.mean(mal_r > thr_r)) if mal_r.size else float("nan")))
        if wr_rows:
            pd.DataFrame(wr_rows).to_csv(args.output_dir / "within_region_contrast.csv",
                                         index=False)
            lines.append("  WITHIN-REGION contrast — Low_signal vs REFERENCE in the SAME region "
                         "(field-effect control; above-ref = fraction over the same-region "
                         "reference 95th pct):")
            for x in wr_rows:
                lines.append(
                    f"    {x['region']:<26} ref_med={x['ref_median']:+.3f}  "
                    f"LS_med={x['lowsignal_median']:+.3f}  mal_med={x['malignant_median']:+.3f}"
                    f"  | LS above-ref={x['lowsignal_above_ref']:.1%}  "
                    f"mal above-ref={x['malignant_above_ref']:.1%}")
            lines.append("    (LS above-ref >> 5% => tumour-region Low_signal is malignant "
                         "BEYOND the field effect; mal above-ref = positive-control ceiling.)")
        lines.append("")

    # SECONDARY: chromosome-arm means (losses are the detectable GBM signal on this panel).
    loss_arms = [c for c in ("chr10", "chr14", "chr15", "chr22") if c in arm.columns]
    lines.append(f"CHROMOSOME-ARM means (secondary; chr7 gain is weak by expression, loss "
                 f"arms {loss_arms} carry the signal):")
    for g in lsg:
        vals = "  ".join(f"{c}={arm.loc[g, c]:+.4f}"
                         for c in (["chr7"] + loss_arms) if c in arm.columns)
        lines.append(f"  {g:<34} {vals}")

    # nearest class by cosine of the full mean profiles.
    lines.append("\nNearest class by mean cosine similarity of profiles:")
    for g in lsg:
        sims = {cl: sim.loc[g, [o for o in order if group_class.get(o) == cl and o != g]].mean()
                for cl in ("malignant", "reference")}
        verdict = max(sims, key=lambda k: (sims[k] if np.isfinite(sims[k]) else -9))
        lines.append(f"  {g:<34} malignant={sims['malignant']:.3f}  "
                     f"reference={sims['reference']:.3f}  -> closer to {verdict.upper()}")

    # ---- PER-DONOR Low_signal malignancy (patient-specificity) ----------------------
    # The malignant fraction of Low_signal is patient-specific, so a single cohort number
    # hides the biology. Report, per donor, the fraction of Low_signal cells above the
    # malignant threshold (overall + by region) and a donor x region heatmap.
    if sig_summary is not None and "tissue_section" in obs:
        sig_arr = obs["mal_sig"].to_numpy()
        donor_arr = obs["tissue_section"].astype(str).str.split("__").str[1].to_numpy()
        ls_mask = (ct == args.lowsignal_label) & np.isfinite(sig_arr)
        dfd = pd.DataFrame({"donor": donor_arr[ls_mask], "region": region[ls_mask],
                            "mal_sig": sig_arr[ls_mask]})

        def donor_breakdown(thr, suffix, summary_lines=None):
            fmal = lambda s: float((s > thr).mean())
            bd = (dfd.groupby("donor")["mal_sig"]
                  .agg(n="size", median="median", malignant_frac=fmal)
                  .sort_values("malignant_frac", ascending=False))
            bd.to_csv(args.output_dir / f"lowsignal_by_donor{suffix}.csv")
            pv = (dfd.groupby(["donor", "region"])["mal_sig"].apply(fmal)
                  .unstack().reindex(index=bd.index))
            pv.to_csv(args.output_dir / f"lowsignal_by_donor_region{suffix}.csv")
            if summary_lines is not None:
                summary_lines.append(f"\nPER-DONOR Low_signal malignant fraction (> {thr:.3g}; "
                                     "patient-specificity — varies hugely, no single number):")
                for d, r in bd.iterrows():
                    summary_lines.append(f"  donor {d:<8} malignant_frac={r['malignant_frac']:.1%}"
                                         f"  median={r['median']:+.3f}  (n={int(r['n']):,})")
            cols = [c for c in (TUMOR_BULK, INFILTRATING_EDGE, CONTRALATERAL) if c in pv.columns]
            Hm = pv[cols].to_numpy(dtype=float)
            fig, ax = plt.subplots(figsize=(6, max(4, 0.36 * len(bd) + 1.5)))
            im = ax.imshow(Hm, aspect="auto", cmap="magma", vmin=0, vmax=1)
            ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=20, fontsize=8)
            ax.set_yticks(range(len(bd))); ax.set_yticklabels(bd.index, fontsize=8)
            for i in range(Hm.shape[0]):
                for j in range(Hm.shape[1]):
                    if np.isfinite(Hm[i, j]):
                        ax.text(j, i, f"{Hm[i, j]:.0%}", ha="center", va="center", fontsize=7,
                                color="white" if Hm[i, j] < 0.6 else "black")
            ax.set_ylabel("donor (sorted by overall Low_signal malignant fraction)")
            title = "Low_signal malignant fraction by donor x region"
            ax.set_title(title + (f" (thr {thr:.3g})" if suffix else ""))
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="malignant fraction")
            fig.tight_layout()
            fig.savefig(args.output_dir / f"lowsignal_by_donor{suffix}.png", dpi=180,
                        bbox_inches="tight")
            plt.close(fig)

        # strict (sig_thr) version — unchanged filenames + SUMMARY block
        donor_breakdown(sig_thr, "", summary_lines=lines)
        (args.output_dir / "SUMMARY.txt").write_text("\n".join(lines) + "\n")
        # optional sensitive version at a chosen cutoff (e.g. the bimodal trough), NEW files
        if args.donor_threshold is not None:
            donor_breakdown(args.donor_threshold, f"_thr{args.donor_threshold:g}")
    else:
        (args.output_dir / "SUMMARY.txt").write_text("\n".join(lines) + "\n")

    print("\n" + "\n".join(lines))
    print(f"\nWrote tables + plots + SUMMARY.txt to {args.output_dir}")


if __name__ == "__main__":
    main()
