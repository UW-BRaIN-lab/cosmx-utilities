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
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_reference_types(path: Path) -> list[str]:
    return [ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()
            if ln.split("#", 1)[0].strip()]


def load_concat(cnv_dir: Path, celltype_key: str, region_key: str):
    """Concatenate per-section CNV results into arrays (no full AnnData needed downstream)."""
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
                (celltype_key, region_key, "cnv_score", "tissue_section") if c in a.obs}
        obs_parts.append(pd.DataFrame(cols, index=a.obs.index))
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
        centroid = np.asarray(X[mal_rows].mean(axis=0)).ravel()
        cnorm = float(np.linalg.norm(centroid))
        if cnorm > 0:
            unit = centroid / cnorm
            dots = np.asarray(X @ unit).ravel()
            rn = np.sqrt(np.asarray(X.multiply(X).sum(axis=1)).ravel())
            with np.errstate(divide="ignore", invalid="ignore"):
                obs["mal_sig"] = np.where(rn > 0, dots / rn, 0.0)
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

    (args.output_dir / "SUMMARY.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nWrote tables + plots + SUMMARY.txt to {args.output_dir}")


if __name__ == "__main__":
    main()
