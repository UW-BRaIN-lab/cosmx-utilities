#!/usr/bin/env python3
"""Stage 5e (diagnostic): characterise the Low_signal population beyond the malignant/normal
CNV split — Phase 1 of the PI's "No Dominant Identity" workflow.

The malignant-vs-normal CNV cut (compare_insitucnv_groups.py) is the FIRST cut, not the whole
case. This script layers the cheap, CNV-reuse diagnostics on top of it — no new typing run,
no new CNV array — by rebuilding the per-cell malignant-signature from the per-section CNV
outputs (exactly as compare/edge_vs_bulk) and joining per-cell morphology + depth from the
typed cohort. It answers, in order:

  STAGE 0  reference / negative-control integrity
    - are the diploid-reference cells and the negative-control (contralateral Low_signal)
      cells non-overlapping populations? (must be, or the calibration is circular);
    - is the "contralateral uninvolved" tissue truly uninvolved? report the per-donor CNV
      floor and test whether the few contralateral CNV-high cells are spatially SCATTERED
      (noise / false-positive floor) or CLUSTERED (a real micro-focus of infiltration that
      would compress dynamic range).

  STAGE 1  edge dilution (the key CNV failure mode at the margin)
    - neighbourhood smoothing pulls an infiltrating tumour cell's CNV toward diploid because
      its neighbours are diploid. Plot malignant-signature vs LOCAL MALIGNANT DENSITY *within
      expression state*: if known-malignant cells in normal-dominated neighbourhoods show a
      depressed signature, that is a false-negative gradient (dilution), not biology — and the
      sparse Low_signal infiltrators are under-called by the same factor;
    - report the informative-gene count per chromosome so a sparse arm is not over-read.

  STAGE 2  doublet / segmentation-artifact screen (area x count only; the spatial-spillover
    half needs per-cell programs and is deferred to Phase 2)
    - doublets run large-and-high-count. Compare the joint (cell area x total counts)
      distribution of CNV-high Low_signal against confidently-typed singlets and quantify the
      large-high tail. Prior: Low_signal median depth is LOWER than typed, so expect a modest
      doublet fraction — this measures it.

  SPATIAL MAP  do the Low_signal CNV-malignant cells trace the infiltrating margin? (If so,
    that is the headline, not an artifact.) Faceted spatial scatter per tissue section.

Reads (--cnv-dir): per-section *_cnv.h5ad from run_insitucnv.py (X_cnv + spatial + obs).
      (--typed-h5ad): cosmx_typed.h5ad, for per-cell Area / total_counts / typing prob.
Writes (--output-dir): cell_cnv_table.csv.gz (per-cell master, reused by later phases) +
      stage0_contralateral_floor.{csv,png}, edge_dilution.{csv,png}, chr_informative_genes.csv,
      doublet_screen.{csv,png}, spatial_margin_map.png, and DIAGNOSTICS_SUMMARY.txt.

Usage:
    python pipeline/python/insitucnv_lowsignal_diagnostics.py \\
        --cnv-dir persection --typed-h5ad cosmx_typed.h5ad \\
        --reference-file pipeline/reference/insitucnv_reference_types.txt \\
        --output-dir diagnostics --sensitive-threshold 0.45
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
from sklearn.neighbors import NearestNeighbors

from compare_insitucnv_groups import (
    CONTRALATERAL, INFILTRATING_EDGE, TUMOR_BULK, DEFAULT_MALIGNANT,
    load_concat, malignant_signature, read_reference_types)

AREA_CANDIDATES = ["Area", "Area.um2", "area", "cell_area", "Area_um2"]
DEPTH_CANDIDATES = ["total_counts", "nCount_RNA", "nCount", "qc_gene_counts"]
PROB_CANDIDATES = ["insitutype_prob", "prob", "typing_prob"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cnv-dir", type=Path, required=True,
                   help="Directory of *_cnv.h5ad from run_insitucnv.py.")
    p.add_argument("--typed-h5ad", type=Path, required=True,
                   help="cosmx_typed.h5ad (obs-only read for Area / depth / typing prob).")
    p.add_argument("--reference-file", type=Path, required=True,
                   help="Diploid reference cell_type list (defines the calibration controls).")
    p.add_argument("--malignant-groups", default=",".join(DEFAULT_MALIGNANT),
                   help="Comma-separated malignant cell types (positive control).")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--sig-threshold", type=float, default=None,
                   help="Malignant-call cutoff on the signature. Default: 95th pct of the "
                        "negative controls (diploid reference + contralateral Low_signal), "
                        "matching compare_insitucnv_groups.py.")
    p.add_argument("--sensitive-threshold", type=float, default=None,
                   help="Optional 2nd, more sensitive cutoff (e.g. the bimodal trough ~0.45) "
                        "reported alongside the strict one.")
    p.add_argument("--k-neighbors", type=int, default=30,
                   help="Spatial neighbours for local-density / dilution (per tissue section).")
    p.add_argument("--map-sections", type=int, default=6,
                   help="How many tissue sections to render in the spatial margin map "
                        "(largest edge-containing sections first).")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


# ------------------------------------------------------------------ joins / setup ----------
def join_typed_obs(obs: pd.DataFrame, typed_h5ad: Path) -> tuple[str | None, str | None]:
    """Join per-cell Area / total_counts / typing prob from the typed cohort (obs-only,
    backed) onto ``obs`` in place. Returns the resolved (area_col, depth_col) names."""
    tobs = ad.read_h5ad(typed_h5ad, backed="r").obs
    area_col = next((c for c in AREA_CANDIDATES if c in tobs.columns), None)
    depth_col = next((c for c in DEPTH_CANDIDATES if c in tobs.columns), None)
    prob_col = next((c for c in PROB_CANDIDATES if c in tobs.columns), None)
    if area_col:
        obs["area"] = pd.to_numeric(tobs[area_col].reindex(obs.index), errors="coerce").to_numpy()
    if depth_col:
        obs["depth"] = pd.to_numeric(tobs[depth_col].reindex(obs.index), errors="coerce").to_numpy()
    if prob_col:
        obs["typing_prob"] = pd.to_numeric(tobs[prob_col].reindex(obs.index),
                                           errors="coerce").to_numpy()
    print(f"Joined typed obs: area='{area_col}', depth='{depth_col}', prob='{prob_col}' "
          f"for {obs.index.size:,} cells.")
    return area_col, depth_col


def per_section_neighbor_stats(obs: pd.DataFrame, k: int) -> pd.DataFrame:
    """Per tissue section (global-px coords only collide across slides), compute for every
    cell: local malignant-class density (fraction of k spatial NN that are malignant-class),
    local malignant-CALL density (fraction of k NN called malignant by CNV), and distance to
    the nearest malignant-class cell. Sections without spatial coords are left as NaN."""
    n = obs.index.size
    local_mal_density = np.full(n, np.nan)
    local_call_density = np.full(n, np.nan)
    dist_to_malignant = np.full(n, np.nan)
    if not {"spatial_x", "spatial_y"}.issubset(obs.columns):
        print("  WARN: no spatial coords in obs; skipping neighbour stats.", file=sys.stderr)
        return pd.DataFrame({"local_mal_density": local_mal_density,
                             "local_call_density": local_call_density,
                             "dist_to_malignant": dist_to_malignant}, index=obs.index)
    pos = {cell: i for i, cell in enumerate(obs.index)}
    is_mal_class = obs["is_malignant_class"].to_numpy()
    is_mal_call = obs["is_malignant_call"].to_numpy()
    for section, sub in obs.groupby("tissue_section", observed=True):
        rows = np.array([pos[c] for c in sub.index])
        xy = sub[["spatial_x", "spatial_y"]].to_numpy()
        finite = np.isfinite(xy).all(axis=1)
        if finite.sum() < 2:
            continue
        rows_f, xy_f = rows[finite], xy[finite]
        kk = min(k, len(xy_f) - 1)
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(xy_f)  # +1: self is nearest
        _, idx = nn.kneighbors(xy_f)
        neigh = idx[:, 1:]  # drop self
        mc = is_mal_class[rows_f]
        cc = is_mal_call[rows_f]
        local_mal_density[rows_f] = mc[neigh].mean(axis=1)
        local_call_density[rows_f] = cc[neigh].mean(axis=1)
        mal_idx = np.flatnonzero(mc)
        if mal_idx.size:
            nn_mal = NearestNeighbors(n_neighbors=1).fit(xy_f[mal_idx])
            d, _ = nn_mal.kneighbors(xy_f)
            dist_to_malignant[rows_f] = d.ravel()
    return pd.DataFrame({"local_mal_density": local_mal_density,
                         "local_call_density": local_call_density,
                         "dist_to_malignant": dist_to_malignant}, index=obs.index)


# ------------------------------------------------------------------ stage 0 ----------------
def stage0_reference_integrity(obs, reference_types, lowsignal_label, lines) -> None:
    """Confirm the diploid-reference cells and the negative-control (contralateral
    Low_signal) cells are non-overlapping populations."""
    ct = obs["cell_type"].astype(str)
    ref_ids = set(obs.index[ct.isin(reference_types).to_numpy()])
    negctrl_ids = set(obs.index[((ct == lowsignal_label).to_numpy())
                                & (obs["Region"].to_numpy() == CONTRALATERAL)])
    overlap = ref_ids & negctrl_ids
    lines.append("STAGE 0 — reference / negative-control integrity")
    lines.append(f"  diploid-reference cells:            {len(ref_ids):,}")
    lines.append(f"  negative-control (contra Low_signal): {len(negctrl_ids):,}")
    lines.append(f"  overlap (must be 0):                {len(overlap):,}  "
                 f"=> {'OK, disjoint' if not overlap else 'PROBLEM: circular calibration'}")


def stage0_contralateral_floor(obs, lowsignal_label, sig_thr, out_dir, lines) -> None:
    """Per-donor contralateral CNV floor + scattered-vs-clustered test for the few
    contralateral CNV-high Low_signal cells (real micro-focus vs noise floor)."""
    ct = obs["cell_type"].astype(str).to_numpy()
    contra_ls = (ct == lowsignal_label) & (obs["Region"].to_numpy() == CONTRALATERAL)
    sub = obs.loc[contra_ls].copy()
    lines.append("\nSTAGE 0 — is contralateral truly uninvolved? (empirical CNV floor)")
    if not len(sub):
        lines.append("  no contralateral Low_signal cells found.")
        return
    call = (sub["mal_sig"] > sig_thr).to_numpy()
    rows = []
    for donor, g in sub.groupby("Case", observed=True):
        c = (g["mal_sig"] > sig_thr).to_numpy()
        rows.append(dict(donor=str(donor), n=int(len(g)), floor=float(c.mean()),
                         median_sig=float(g["mal_sig"].median())))
    fdf = pd.DataFrame(rows).sort_values("floor", ascending=False)
    # clustered vs scattered: for contralateral CNV-high cells, how enriched are their
    # neighbours for other CNV-high cells vs the global contralateral CNV-high rate? ratio
    # ~1 => scattered (noise floor); >>1 => clustered (a real infiltrating micro-focus).
    base = float(call.mean())
    hi = sub.loc[call]
    clus = float(hi["local_call_density"].mean()) if len(hi) and \
        "local_call_density" in hi and hi["local_call_density"].notna().any() else float("nan")
    ratio = (clus / base) if base > 0 and np.isfinite(clus) else float("nan")
    fdf.to_csv(out_dir / "stage0_contralateral_floor.csv", index=False)
    lines.append(f"  overall contralateral CNV-high floor: {base:.2%} "
                 f"(n={len(sub):,}, thr={sig_thr:.3g})")
    lines.append(f"  per-donor floor range: {fdf['floor'].min():.2%} .. {fdf['floor'].max():.2%}")
    lines.append(f"  scattered-vs-clustered ratio (neighbour CNV-high enrichment): {ratio:.2f} "
                 f"=> {'SCATTERED (noise floor)' if np.isfinite(ratio) and ratio < 2 else 'CLUSTERED — inspect for a real micro-focus' if np.isfinite(ratio) else 'n/a'}")

    fig, ax = plt.subplots(figsize=(7, max(3, 0.32 * len(fdf) + 1.5)))
    ax.barh(np.arange(len(fdf)), fdf["floor"], color="#1b9e77")
    ax.set_yticks(np.arange(len(fdf))); ax.set_yticklabels(fdf["donor"], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(base, color="k", ls="--", lw=1, label=f"cohort floor {base:.1%}")
    ax.set_xlabel("contralateral Low_signal CNV-high fraction (false-positive floor)")
    ax.set_title("Stage 0: contralateral CNV floor per donor")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "stage0_contralateral_floor.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ stage 1 ----------------
def stage1_edge_dilution(obs, lowsignal_label, out_dir, lines) -> None:
    """malignant-signature vs local malignant density, WITHIN expression state. A drop for
    the known-malignant positive control at low density = neighbourhood-smoothing dilution."""
    lines.append("\nSTAGE 1 — edge dilution (signature vs local malignant density)")
    if "local_mal_density" not in obs or obs["local_mal_density"].isna().all():
        lines.append("  no neighbour stats available; skipped.")
        return
    ct = obs["cell_type"].astype(str).to_numpy()
    strata = [("malignant", obs["is_malignant_class"].to_numpy()),
              (lowsignal_label, ct == lowsignal_label),
              ("reference", obs["class"].to_numpy() == "reference")]
    bins = np.linspace(0, 1, 6)  # local malignant-density bins 0-0.2-...-1.0
    centers = 0.5 * (bins[:-1] + bins[1:])
    rows = []
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {"malignant": "#c51b8a", lowsignal_label: "#d95f0e", "reference": "#2c7fb8"}
    for name, mask in strata:
        d = obs.loc[mask, ["local_mal_density", "mal_sig"]].dropna()
        if not len(d):
            continue
        which = np.digitize(d["local_mal_density"].to_numpy(), bins[1:-1])
        meds = [d["mal_sig"].to_numpy()[which == i]
                for i in range(len(centers))]
        med_vals = [float(np.median(v)) if len(v) else np.nan for v in meds]
        n_vals = [int(len(v)) for v in meds]
        ax.plot(centers, med_vals, "o-", color=colors.get(name, "#555"),
                label=f"{name} (n={len(d):,})")
        for c, m, nb in zip(centers, med_vals, n_vals):
            rows.append(dict(stratum=name, density_bin_center=float(c),
                             median_sig=m, n=nb))
    pd.DataFrame(rows).to_csv(out_dir / "edge_dilution.csv", index=False)
    ax.set_xlabel("local malignant density (fraction of spatial neighbours that are malignant)")
    ax.set_ylabel("median malignant-signature")
    ax.legend(fontsize=8)
    ax.set_title("Stage 1: does the signature drop where malignant neighbours are sparse?\n"
                 "(a drop for the MALIGNANT control = edge-dilution false-negative gradient)")
    fig.tight_layout()
    fig.savefig(out_dir / "edge_dilution.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # verdict: compare malignant-control signature in the lowest vs highest density bin
    mal = obs.loc[obs["is_malignant_class"].to_numpy(), ["local_mal_density", "mal_sig"]].dropna()
    if len(mal):
        lo = mal.loc[mal["local_mal_density"] < 0.2, "mal_sig"].median()
        hi = mal.loc[mal["local_mal_density"] >= 0.8, "mal_sig"].median()
        if not (np.isfinite(lo) and np.isfinite(hi)):
            lines.append(f"  malignant-control signature: low-density(<0.2) median={lo:+.3f}  "
                         f"high-density(>=0.8) median={hi:+.3f}")
            lines.append("  => insufficient malignant cells spanning the density range to "
                         "assess dilution (expect a full range only where a bulk/margin exists).")
        else:
            drop = hi - lo
            lines.append(f"  malignant-control signature: low-density(<0.2) median={lo:+.3f}  "
                         f"high-density(>=0.8) median={hi:+.3f}  drop={drop:+.3f}")
            lines.append("  => " + ("DILUTION present — known-malignant cells are under-scored in "
                                    "normal-dominated neighbourhoods, so sparse Low_signal "
                                    "infiltrators at the margin are under-called by the same factor."
                                    if drop > 0.05 else
                                    "little dilution — signature roughly flat across neighbourhood "
                                    "composition."))


def chr_informative_genes(cnv_dir: Path, out_dir: Path, lines) -> None:
    """How many panel genes inform each chromosome (arm-level CNV on a targeted panel is only
    as trustworthy as its gene count). Read one section's var (chromosome column)."""
    files = sorted(cnv_dir.glob("*_cnv.h5ad"))
    if not files:
        return
    var = ad.read_h5ad(files[0], backed="r").var
    chrom_col = next((c for c in ("chromosome", "chrom", "chr") if c in var.columns), None)
    lines.append("\nSTAGE 1 — informative genes per chromosome (sparse arms = over-read risk)")
    if chrom_col is None:
        lines.append("  no chromosome column in section var; skipped.")
        return
    counts = (var[chrom_col].astype(str).value_counts()
              .rename_axis("chromosome").reset_index(name="n_genes"))
    counts["chromosome"] = counts["chromosome"].str.replace("chr", "", regex=False)
    counts = counts[counts["chromosome"].str.isdigit()]
    counts["chr_num"] = counts["chromosome"].astype(int)
    counts = counts.sort_values("chr_num")
    counts.to_csv(out_dir / "chr_informative_genes.csv", index=False)
    sparse = counts.loc[counts["n_genes"] < 20, "chromosome"].tolist()
    key = {c: int(counts.loc[counts["chromosome"] == c, "n_genes"].iloc[0])
           for c in ("7", "9", "10") if c in set(counts["chromosome"])}
    lines.append(f"  GBM-relevant arms: " + "  ".join(f"chr{c}={n}" for c, n in key.items()))
    lines.append(f"  sparse arms (<20 genes; interpret trends cautiously): "
                 f"chr{', chr'.join(sparse) if sparse else '(none)'}")


# ------------------------------------------------------------------ stage 2 ----------------
def stage2_doublet_screen(obs, lowsignal_label, sig_thr, out_dir, lines, suffix="") -> None:
    """Joint (cell area x total counts) of CNV-high Low_signal vs typed singlets. Doublets run
    large-and-high-count; the large-high tail is the suspect fraction."""
    lines.append(f"\nSTAGE 2 — doublet / segmentation screen (area x count; thr={sig_thr:.3g})")
    if "area" not in obs or "depth" not in obs:
        lines.append("  area or depth column unavailable in typed obs; skipped.")
        return
    ct = obs["cell_type"].astype(str).to_numpy()
    typed = obs.loc[(ct != lowsignal_label)].dropna(subset=["area", "depth"])
    ls_hi = obs.loc[(ct == lowsignal_label) & (obs["mal_sig"] > sig_thr)].dropna(
        subset=["area", "depth"])
    if not len(typed) or not len(ls_hi):
        lines.append("  not enough cells with area+depth; skipped.")
        return
    # "large-high" doublet quadrant defined by the typed-singlet 95th percentiles
    area_hi = float(np.percentile(typed["area"], 95))
    depth_hi = float(np.percentile(typed["depth"], 95))
    def frac_quadrant(df):
        return float(((df["area"] > area_hi) & (df["depth"] > depth_hi)).mean())
    f_typed, f_lshi = frac_quadrant(typed), frac_quadrant(ls_hi)
    pd.DataFrame([
        dict(group="typed_singlets", n=int(len(typed)),
             area_median=float(typed["area"].median()),
             depth_median=float(typed["depth"].median()), doublet_quadrant_frac=f_typed),
        dict(group="lowsignal_CNVhigh", n=int(len(ls_hi)),
             area_median=float(ls_hi["area"].median()),
             depth_median=float(ls_hi["depth"].median()), doublet_quadrant_frac=f_lshi),
    ]).to_csv(out_dir / f"doublet_screen{suffix}.csv", index=False)
    lines.append(f"  large-high quadrant = area>{area_hi:.0f} AND depth>{depth_hi:.0f} "
                 f"(typed-singlet 95th pcts)")
    depth_ratio = float(ls_hi["depth"].median() / typed["depth"].median()) \
        if typed["depth"].median() else float("nan")
    enriched = f_lshi > 2 * max(f_typed, 0.01)  # doublet quadrant over the typed-singlet baseline
    lines.append(f"  typed singlets in quadrant:      {f_typed:.2%}  "
                 f"(area_med={typed['area'].median():.0f}, depth_med={typed['depth'].median():.0f})")
    lines.append(f"  CNV-high Low_signal in quadrant:  {f_lshi:.2%}  "
                 f"(area_med={ls_hi['area'].median():.0f}, depth_med={ls_hi['depth'].median():.0f}; "
                 f"depth ratio vs typed = {depth_ratio:.2f})")
    lines.append("  => " + (
        f"elevated large-high tail ({f_lshi:.1%} vs {f_typed:.1%} baseline) — a non-trivial "
        "doublet fraction; flag these for tighter segmentation / exclusion."
        if enriched else
        "doublet-suspect fraction is NOT enriched over typed singlets"
        + (" and depth runs lower" if depth_ratio < 1 else "")
        + " — segmentation doublets are not the bulk of CNV-high Low_signal."))

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.hexbin(np.log1p(typed["depth"]), typed["area"], gridsize=45, cmap="Greys",
              bins="log", mincnt=1)
    ax.scatter(np.log1p(ls_hi["depth"]), ls_hi["area"], s=4, alpha=0.35, color="#d73027",
               label=f"CNV-high Low_signal (n={len(ls_hi):,})")
    ax.axvline(np.log1p(depth_hi), color="k", ls="--", lw=1)
    ax.axhline(area_hi, color="k", ls="--", lw=1)
    ax.set_xlabel("log1p(total counts)"); ax.set_ylabel("cell area")
    ax.set_title("Stage 2: area x count — typed singlets (grey) vs CNV-high Low_signal (red)\n"
                 "top-right quadrant = large-high doublet suspects")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / f"doublet_screen{suffix}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ stage 2b ---------------
def stage2b_density_screen(obs, lowsignal_label, out_dir, lines) -> None:
    """Bloated / over-segmented mask screen — a failure mode ORTHOGONAL to doublets. Doublets sit
    top-RIGHT of area x count (large + high count); a ballooned or over-segmented mask sits top-
    LEFT (large area, low count), so it passes the doublet filter while carrying few transcripts
    per unit area (background / debris / a sparse cell whose mask over-grew). Single-axis QC
    (min counts, max area, max negprobe) cannot catch it — it is a JOINT area/count anomaly, i.e.
    low transcript DENSITY = counts / area. Flags the suspects, checks whether they carry the same
    malignant signature (a sensitivity on the CNV-high fraction), and emits their cell IDs for
    mask-level (Napari) inspection."""
    lines.append("\nSTAGE 2b — bloated-mask / low-transcript-density screen (orthogonal to doublets)")
    if "density" not in obs or "area" not in obs:
        lines.append("  density/area unavailable; skipped.")
        return
    ct = obs["cell_type"].astype(str).to_numpy()
    typed = obs.loc[ct != lowsignal_label].dropna(subset=["area", "depth", "density"])
    ws = obs.loc[(ct == lowsignal_label) & obs["is_malignant_call"].to_numpy()].dropna(
        subset=["area", "depth", "density"])
    if not len(typed) or not len(ws):
        lines.append("  not enough cells with area+depth; skipped.")
        return
    area_hi = float(np.percentile(typed["area"], 95))
    depth_hi = float(np.percentile(typed["depth"], 95))
    dens_lo = float(np.percentile(typed["density"], 5))      # singlet low-density floor
    large_lowcount = (ws["area"] > area_hi) & (ws["depth"] < depth_hi)   # the top-left cloud
    low_density = ws["density"] < dens_lo
    bloated = (ws["area"] > area_hi) & low_density           # strict large-and-low-density corner
    flag = large_lowcount | bloated

    def pct(m):
        return f"{int(m.sum()):>7,} ({m.mean():.2%})"
    lines.append(f"  typed-singlet density: median={typed['density'].median():.4f}, "
                 f"5th pct={dens_lo:.4f}; CNV-high Low_signal median={ws['density'].median():.4f}")
    lines.append(f"  large area (>{area_hi:.0f}) + low count (<{depth_hi:.0f}):  {pct(large_lowcount)}")
    lines.append(f"  low transcript density (< singlet 5th pct):     {pct(low_density)}")
    lines.append(f"  bloated-mask corner (large AND low-density):    {pct(bloated)}")
    lines.append(f"  flagged mal_sig median={ws.loc[flag, 'mal_sig'].median():.3f} vs "
                 f"rest={ws.loc[~flag, 'mal_sig'].median():.3f}  "
                 f"(=> excluding them moves the CNV-high count by {flag.mean():.2%})")
    if "Region" in ws:
        rc = ws.loc[flag, "Region"].value_counts()
        lines.append("  flagged region: " + ", ".join(f"{k} {int(v):,}" for k, v in rc.items()))
    lines.append("  => a distinct QC flag (bloated/over-segmented masks); single-axis QC "
                 "(min-count/max-area/negprobe) misses it because it is an area/count RATIO "
                 "anomaly. Same malignant signature => set aside for mask-level review, not "
                 "auto-excluded.")

    cols = [c for c in ("cell_type", "Region", "Case", "slide_id", "tissue_section",
                        "spatial_x", "spatial_y", "area", "depth", "density", "mal_sig")
            if c in ws.columns]
    flagged = ws.loc[flag, cols].copy()
    flagged.index.name = "cell_id"
    flagged.to_csv(out_dir / "flagged_lowdensity_cells.csv")
    lines.append(f"  wrote {len(flagged):,} flagged cell IDs -> flagged_lowdensity_cells.csv "
                 "(for mask-level / Napari inspection).")

    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(0, float(np.percentile(typed["density"], 99)), 60)
    ax.hist(typed["density"], bins=bins, density=True, histtype="step", lw=2, color="#888888",
            label=f"typed singlets (n={len(typed):,})")
    ax.hist(ws["density"], bins=bins, density=True, histtype="step", lw=2, color="#d73027",
            label=f"CNV-high Low_signal (n={len(ws):,})")
    ax.axvline(dens_lo, color="k", ls="--", lw=1, label=f"singlet 5th pct ({dens_lo:.3f})")
    ax.set_xlabel("transcript density (total counts / cell area)")
    ax.set_ylabel("density"); ax.legend(fontsize=8)
    ax.set_title("Stage 2b: transcript density — bloated / over-segmented masks (low = suspect)")
    fig.tight_layout()
    fig.savefig(out_dir / "density_screen.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ spatial map ------------
def spatial_margin_map(obs, lowsignal_label, sig_thr, n_sections, out_dir, lines,
                       suffix="") -> None:
    """Faceted per-section spatial scatter: do the Low_signal CNV-malignant cells trace the
    infiltrating margin? Prioritise the largest sections that contain an infiltrating edge."""
    lines.append(f"\nSPATIAL MAP — do CNV-malignant Low_signal cells trace the margin? "
                 f"(thr={sig_thr:.3g})")
    if not {"spatial_x", "spatial_y"}.issubset(obs.columns):
        lines.append("  no spatial coords; skipped.")
        return
    ct = obs["cell_type"].astype(str).to_numpy()
    region = obs["Region"].to_numpy()
    edge_sections = (obs.loc[region == INFILTRATING_EDGE, "tissue_section"]
                     .value_counts())
    if not len(edge_sections):
        edge_sections = obs["tissue_section"].value_counts()
    chosen = edge_sections.index[:n_sections].tolist()
    dropped = len(edge_sections) - len(chosen)
    lines.append(f"  rendered {len(chosen)} of {len(edge_sections)} edge-containing sections "
                 f"({dropped} not shown — largest-first).")

    ncol = min(3, len(chosen))
    nrow = int(np.ceil(len(chosen) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 5 * nrow), squeeze=False)
    is_ls = ct == lowsignal_label
    is_mal_class = obs["is_malignant_class"].to_numpy()
    call = (obs["mal_sig"] > sig_thr).to_numpy()
    for ax, section in zip(axes.ravel(), chosen):
        m = (obs["tissue_section"] == section).to_numpy()
        sx = obs.loc[m, "spatial_x"].to_numpy(); sy = obs.loc[m, "spatial_y"].to_numpy()
        ls_m, mc_m, cl_m = is_ls[m], is_mal_class[m], call[m]
        other = ~ls_m & ~mc_m
        ax.scatter(sx[other], sy[other], s=1, color="#dddddd", label="other typed")
        ax.scatter(sx[mc_m], sy[mc_m], s=1.5, color="#fbb4b9", label="malignant (typed)")
        ax.scatter(sx[ls_m & ~cl_m], sy[ls_m & ~cl_m], s=2, color="#4575b4",
                   label="Low_signal · CNV-normal")
        ax.scatter(sx[ls_m & cl_m], sy[ls_m & cl_m], s=3, color="#d73027",
                   label="Low_signal · CNV-malignant")
        ax.set_title(section, fontsize=7)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        ax.invert_yaxis()
    for ax in axes.ravel()[len(chosen):]:
        ax.axis("off")
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8, markerscale=4)
    fig.suptitle("Spatial map: Low_signal CNV-malignant vs the typed tumour bulk")
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    fig.savefig(out_dir / f"spatial_margin_map{suffix}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ main -------------------
def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference_types = set(read_reference_types(args.reference_file))
    malignant_types = {s.strip() for s in args.malignant_groups.split(",") if s.strip()}

    X, obs, _ = load_concat(args.cnv_dir, args.celltype_key, args.region_key, with_spatial=True)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    obs = obs.rename(columns={args.celltype_key: "cell_type", args.region_key: "Region"})

    sig = malignant_signature(X, obs, malignant_types, "cell_type")
    if sig is None:
        sys.exit("ERROR: no usable malignant consensus (no malignant cells or zero centroid).")
    obs["mal_sig"] = sig

    # classes + boolean helpers
    ct = obs["cell_type"].astype(str)
    obs["is_malignant_class"] = ct.isin(malignant_types).to_numpy()
    obs["class"] = np.where(ct.to_numpy() == args.lowsignal_label, "low_signal",
                    np.where(obs["is_malignant_class"], "malignant",
                     np.where(ct.isin(reference_types).to_numpy(), "reference", "other")))

    # threshold: 95th pct of negative controls (reference + contralateral Low_signal)
    neg_mask = (obs["class"] == "reference").to_numpy() | \
               ((ct.to_numpy() == args.lowsignal_label) & (obs["Region"].to_numpy() == CONTRALATERAL))
    neg_sig = obs.loc[neg_mask, "mal_sig"].dropna()
    sig_thr = args.sig_threshold if args.sig_threshold is not None else (
        float(np.percentile(neg_sig, 95)) if len(neg_sig) else float("nan"))
    obs["is_malignant_call"] = (obs["mal_sig"] > sig_thr).to_numpy()
    print(f"Malignant-call threshold = {sig_thr:.4g} "
          f"({'user-supplied' if args.sig_threshold is not None else '95th pct of neg controls'}).")

    join_typed_obs(obs, args.typed_h5ad)
    if {"area", "depth"}.issubset(obs.columns):  # transcript density (for the bloated-mask screen)
        obs["density"] = (obs["depth"] / obs["area"]).replace([np.inf, -np.inf], np.nan)
    neigh = per_section_neighbor_stats(obs, args.k_neighbors)
    obs = obs.join(neigh)

    # ---- per-cell master table (reused by later phases) ----
    keep = ["cell_type", "class", "Region", "Case", "Block", "slide_id", "tissue_section",
            "spatial_x", "spatial_y", "cnv_score", "mal_sig", "is_malignant_call",
            "is_malignant_class", "area", "depth", "density", "typing_prob",
            "local_mal_density", "local_call_density", "dist_to_malignant"]
    table = obs[[c for c in keep if c in obs.columns]].copy()
    table.index.name = "cell_id"
    table.to_csv(args.output_dir / "cell_cnv_table.csv.gz", index=True)
    print(f"Wrote per-cell master table: {len(table):,} cells x {table.shape[1]} cols.")

    # ---- diagnostics ----
    lines = [f"InSituCNV Low_signal diagnostics (Phase 1) — {len(obs):,} cells, "
             f"threshold {sig_thr:.4g}\n"]
    stage0_reference_integrity(obs, reference_types, args.lowsignal_label, lines)
    stage0_contralateral_floor(obs, args.lowsignal_label, sig_thr, args.output_dir, lines)
    stage1_edge_dilution(obs, args.lowsignal_label, args.output_dir, lines)
    chr_informative_genes(args.cnv_dir, args.output_dir, lines)
    stage2_doublet_screen(obs, args.lowsignal_label, sig_thr, args.output_dir, lines)
    stage2b_density_screen(obs, args.lowsignal_label, args.output_dir, lines)
    spatial_margin_map(obs, args.lowsignal_label, sig_thr, args.map_sections,
                       args.output_dir, lines)
    if args.sensitive_threshold is not None:
        st = args.sensitive_threshold
        suffix = f"_thr{st:g}"
        stage2_doublet_screen(obs, args.lowsignal_label, st, args.output_dir, lines, suffix)
        spatial_margin_map(obs, args.lowsignal_label, st, args.map_sections,
                           args.output_dir, lines, suffix)

    (args.output_dir / "DIAGNOSTICS_SUMMARY.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nWrote tables + plots + DIAGNOSTICS_SUMMARY.txt to {args.output_dir}")


if __name__ == "__main__":
    main()
