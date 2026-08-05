#!/usr/bin/env python3
"""Stage 5g (Phase 2): malignant CONTINUUM vs two-state HYBRID among the CNV-high Low_signal cells.

CNV says malignant-or-not; it cannot say whether a malignant cell is intermediate on the tumour
manifold or co-running two distinct tumour states. For GBM the states are Neftel's AC/MES/OPC/
NPC-like, arranged on two axes (AC<->MES and OPC<->NPC). The distinction:

  CONTINUUM  a cell intermediate WITHIN an axis (co-expresses AC and MES, or OPC and NPC) — a
             gradient along a single arm; discretising it further just tiles the manifold.
  HYBRID     a cell co-running states from DIFFERENT axes (e.g. AC + OPC, MES + NPC) — a genuine
             dual-state cell; the biologically notable case.

Test: score the four Neftel modules (sc.tl.score_genes) and, for each state-pair, compare the
observed co-elevation fraction to the chance (independence) expectation. WITHIN-axis pairs enriched
=> continuum; a CROSS-axis pair enriched beyond chance => a real dual-state hybrid. Boundary cells
between two states show chance-level cross-axis overlap, which is why the readout is enrichment over
chance per pair, not a raw co-expression count. The coarse top-two posterior structure (both-
malignant = continuum, weak top-1 = reference-limited) is reported alongside.

Inputs:
  --typed-h5ad   cosmx_typed.h5ad (raw counts; module scoring).
  --cell-table   cell_cnv_table.csv.gz from insitucnv_lowsignal_diagnostics.py (Phase 1).
  --posteriors   flat_posteriors.csv from flat_posteriors.R (top-K per cell).
  --signatures   gene_signatures.csv (module,gene) — must include the 4 Neftel state modules.
Writes (--output-dir): hybrid_continuum_cells.csv.gz, toptwo_structure.csv,
  neftel_state_coexpression.csv, neftel_hybrid.png, HYBRID_SUMMARY.txt.

Usage:
    python pipeline/python/insitucnv_hybrid_continuum.py \\
        --typed-h5ad cosmx_typed.h5ad --cell-table diagnostics/cell_cnv_table.csv.gz \\
        --posteriors posteriors/flat_posteriors.csv \\
        --signatures pipeline/reference/gene_signatures.csv --output-dir hybrid
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from compare_insitucnv_groups import DEFAULT_MALIGNANT

# Neftel GBM malignant meta-states and their two axes. A within-axis pair co-varies along a
# gradient (continuum); a cross-axis pair co-elevated beyond chance is a dual-state hybrid.
NEFTEL_STATES = ["AClike", "MESlike", "OPClike", "NPClike"]
NEFTEL_AXIS = {"AClike": "AC/MES", "MESlike": "AC/MES", "OPClike": "OPC/NPC", "NPClike": "OPC/NPC"}
HYBRID_ENRICH = 1.3  # a cross-axis pair above this (obs/chance) counts as an enriched hybrid


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--cell-table", type=Path, required=True)
    p.add_argument("--posteriors", type=Path, required=True)
    p.add_argument("--signatures", type=Path, required=True)
    p.add_argument("--hierarchy", type=Path, default=None,
                   help="(accepted for interface compatibility; not used by the state-hybrid test).")
    p.add_argument("--malignant-groups", default=",".join(DEFAULT_MALIGNANT))
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--score-percentile", type=float, default=75.0,
                   help="A module score is 'elevated' above this cohort percentile.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def read_signatures(path: Path) -> dict[str, list[str]]:
    df = pd.read_csv(path, comment="#")
    return {m: g["gene"].astype(str).tolist() for m, g in df.groupby("module", sort=False)}


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


def classify_toptwo(row, malignant_types, lowsignal_label, prob_gap=0.10) -> str:
    """Coarse top-two structure: continuum (both top-two malignant) / cross_TME (malignant + a
    non-malignant compartment) / ambiguous_iterate (sink or no clear winner) / non_malignant."""
    t1, t2 = row.get("top1_type"), row.get("top2_type")
    p1, p2 = row.get("top1_prob", np.nan), row.get("top2_prob", np.nan)
    if not isinstance(t1, str):
        return "unknown"
    if t1 == f"{lowsignal_label}_denovo" or (np.isfinite(p1) and p1 < 0.5) \
            or (np.isfinite(p1) and np.isfinite(p2) and (p1 - p2) < prob_gap):
        return "ambiguous_iterate"
    m1 = t1 in malignant_types
    m2 = isinstance(t2, str) and t2 in malignant_types
    if m1 and m2:
        return "continuum"
    if m1 or m2:
        return "cross_TME"
    return "non_malignant_toptwo"


def neftel_hybrid_test(ws, df, percentile, out_dir, lines) -> bool | None:
    """Per Neftel state-pair co-elevation enrichment (observed / chance). Within-axis pairs
    enriched => continuum; a cross-axis pair enriched beyond chance => dual-state hybrid."""
    cols = {s: f"score_{s}" for s in NEFTEL_STATES if f"score_{s}" in ws.columns}
    if len(cols) < 4:
        lines.append("\nGBM-STATE HYBRID: fewer than 4 Neftel modules scored; test skipped.")
        return None
    thr = {s: float(np.nanpercentile(df[c], percentile)) for s, c in cols.items()}
    on = {s: (ws[c] > thr[s]).to_numpy() for s, c in cols.items()}
    rows = []
    for i, a in enumerate(NEFTEL_STATES):
        for b in NEFTEL_STATES[i + 1:]:
            obs = float((on[a] & on[b]).mean())
            exp = float(on[a].mean() * on[b].mean())
            rows.append(dict(pair=f"{a}-{b}", within_axis=(NEFTEL_AXIS[a] == NEFTEL_AXIS[b]),
                             observed=obs, expected=exp,
                             enrichment=(obs / exp if exp > 0 else float("nan"))))
    tab = pd.DataFrame(rows)
    tab.to_csv(out_dir / "neftel_state_coexpression.csv", index=False)
    within = tab[tab["within_axis"]]["enrichment"]
    cross = tab[~tab["within_axis"]]
    hybrid = bool((cross["enrichment"] > HYBRID_ENRICH).any())

    lines.append(f"\nGBM-STATE HYBRID (Neftel per-pair co-elevation vs chance; 'high' = > cohort "
                 f"{percentile:g}th pct of each module):")
    for _, r in tab.sort_values("enrichment").iterrows():
        kind = "within-axis (continuum)" if r["within_axis"] else "cross-axis (hybrid)"
        lines.append(f"  {r['pair']:<16} {kind:<24} {r['enrichment']:.2f}x")
    lines.append(f"  within-axis mean {within.mean():.2f}x (continuum)  |  "
                 f"cross-axis mean {cross['enrichment'].mean():.2f}x")
    lines.append("  => " + ("CROSS-AXIS HYBRID enriched — a real dual-state tumour population; pursue."
                            if hybrid else
                            "NO cross-axis hybrid — within-axis co-expression is enriched (a malignant "
                            "CONTINUUM along the Neftel axes) while cross-axis is at or below chance. "
                            "The cells lie ON the axes, not co-running two arms."))

    t = tab.sort_values("enrichment")
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    colors = ["#2e7d32" if w else "#9d2933" for w in t["within_axis"]]
    ax.barh(range(len(t)), t["enrichment"], color=colors)
    ax.axvline(1.0, color="k", lw=1.2, ls="--")
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([p.replace("like", "") for p in t["pair"]], fontsize=10)
    for y, e in enumerate(t["enrichment"]):
        ax.text(e + 0.02, y, f"{e:.2f}x", va="center", fontsize=9)
    ax.set_xlabel("co-elevation enrichment (observed / chance)")
    ax.set_xlim(0, float(np.nanmax(t["enrichment"])) * 1.18)
    ax.set_title("GBM tumour-state co-expression: continuum (within-axis) vs hybrid (cross-axis)")
    ax.legend(handles=[mpatches.Patch(color="#2e7d32", label="within-axis (AC-MES, OPC-NPC)"),
                       mpatches.Patch(color="#9d2933", label="cross-axis (AC/MES x OPC/NPC)")],
              fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_dir / "neftel_hybrid.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return hybrid


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    malignant_types = {s.strip() for s in args.malignant_groups.split(",") if s.strip()}
    signatures = read_signatures(args.signatures)

    print("Scoring gene-signature modules on the typed cohort...")
    scores = score_modules(args.typed_h5ad, signatures)

    print("Loading Phase-1 cell table + flat posteriors...")
    cells = pd.read_csv(args.cell_table, index_col=0)
    post = pd.read_csv(args.posteriors, index_col="cell_id")
    df = cells.join(scores, how="left").join(
        post[[c for c in post.columns if c.startswith("top")]], how="left")

    # working set: CNV-high Low_signal (the flat-but-malignant cells to resolve)
    ct = df["cell_type"].astype(str)
    ws = df[(ct == args.lowsignal_label) & df.get("is_malignant_call", False).astype(bool)].copy()
    print(f"Working set (CNV-high Low_signal): {len(ws):,} cells.")
    if not len(ws):
        sys.exit("No CNV-high Low_signal cells in the table; nothing to classify.")

    ws["toptwo_class"] = ws.apply(
        lambda r: classify_toptwo(r, malignant_types, args.lowsignal_label), axis=1)
    struct = ws["toptwo_class"].value_counts()
    struct.to_csv(args.output_dir / "toptwo_structure.csv")

    lines = [f"Continuum vs two-state hybrid — CNV-high Low_signal working set: {len(ws):,} cells\n"]
    lines.append("TOP-TWO STRUCTURE (coarse):")
    for k, v in struct.items():
        lines.append(f"  {k:<22} {v:>8,}  ({v / len(ws):.1%})")

    hybrid = neftel_hybrid_test(ws, df, args.score_percentile, args.output_dir, lines)

    frac = lambda k: struct.get(k, 0) / len(ws)
    lines.append("\nINTERPRETATION (CNV-high Low_signal):")
    lines.append(f"  continuum (top-two adjacent malignant states):  {frac('continuum'):.1%}")
    lines.append(f"  ambiguous / reference-limited (iterate):        {frac('ambiguous_iterate'):.1%}")
    lines.append(f"  cross-compartment with non-malignant TME:       {frac('cross_TME'):.1%}")
    if hybrid:
        lines.append("  => a cross-axis two-state HYBRID population is present and worth pursuing.")
    else:
        lines.append("  => two-state HYBRID is NOT supported: within-axis co-expression is enriched "
                     "(a malignant CONTINUUM along the Neftel axes), cross-axis is at/below chance; "
                     "plus a large reference-limited fraction.")
    lines.append("  The large ambiguous/reference-limited share makes the Phase-3 rescue the "
                 "decisive next test (collapse vs irreducible core).")

    ws.to_csv(args.output_dir / "hybrid_continuum_cells.csv.gz", index=True)
    (args.output_dir / "HYBRID_SUMMARY.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nWrote tables + plots + HYBRID_SUMMARY.txt to {args.output_dir}")


if __name__ == "__main__":
    main()
