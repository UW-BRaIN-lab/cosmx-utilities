#!/usr/bin/env python3
"""Stage 5h (Phase 3): did one more rescue iteration collapse the Low_signal sink? (48 -> X)

Reads the per-cell rescue labels (rescue_lowsignal.R) and asks how much of the Low_signal pool
genuinely RESOLVES when profiles are re-derived from the pool and cells reassigned.

The trap (the PI flagged it): a de-novo rescue assigns EVERY cell to some data-derived cluster,
so a high "left the sink" number is guaranteed and does NOT mean identities were resolved — a
rescue iteration can just TILE a flat/continuous pool into arbitrary bins. So we do not treat
de-novo assignment as resolution. Instead:
  - RESOLUTION = the fraction that snapped to an existing NAMED type (real reference identity),
    broken down by compartment.
  - DE-NOVO clusters are scored for MARKER COHERENCE from the rescued profile matrix: a cluster
    whose top markers are dominated by housekeeping / heat-shock / hypoxia genes is continuum
    "tiling" (not a new identity); one with a distinct lineage program is a candidate real
    off-reference population. Only named + candidate-de-novo count as genuinely resolved.
  - CNV CONCORDANCE: cross-tab the expression rescue-to-malignant against the Phase-1 CNV call.

A big genuine-resolution fraction => the pool was reference/transfer-limited (iterate). A small
one with the bulk in tiling/sink => an irreducible flat core, whatever the raw collapse number.

Inputs:
  --rescue       rescue_lowsignal.csv (cell_id, rescue_label, rescue_prob, is_denovo).
  --hierarchy    insitutree_hierarchy.json (leaf -> compartment; also defines the NAMED types).
  --profiles     rescued_profiles.csv (genes x profiles) from rescue_lowsignal.R — for the
                 de-novo marker-coherence check. Optional; without it de-novo stays unclassified.
  --cell-table   cell_cnv_table.csv.gz from Phase 1 (optional; for CNV concordance).
Writes (--output-dir): rescue_summary.txt, rescue_by_compartment.csv, denovo_marker_coherence.csv,
  rescue_rate.png, rescue_vs_cnv.csv.

Usage:
    python pipeline/python/compare_rescue.py --rescue rescue/rescue_lowsignal.csv \\
        --hierarchy pipeline/reference/insitutree_hierarchy.json \\
        --profiles rescue/rescued_profiles.csv \\
        --cell-table diagnostics/cell_cnv_table.csv.gz --output-dir rescue
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compare_insitucnv_groups import DEFAULT_MALIGNANT

# Genes that mark the transcriptionally-flat sink itself, not a lineage: ribosomal / core
# housekeeping + heat-shock (stress) + hypoxia. A de-novo cluster whose top markers are mostly
# these is re-tiling the flat/stressed pool, not resolving a new cell identity.
HK_STRESS = {
    "NDUFA3", "SNRPD1", "SNRPA1", "SNRPG", "SNRPE", "A1BG", "UBB", "UBC", "VPS29", "YWHAZ",
    "TMSB4X", "TMSB10", "ACTB", "ACTG1", "GAPDH", "B2M", "TPT1", "FTL", "FTH1", "EEF1A1",
    "EEF2", "NACA", "PPIA", "NDUFA4", "COX4I1", "ATP5F1E",
    "HSPB1", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA8", "HSPA5", "HSPE1", "HSPD1",
    "DNAJB1", "DNAJA1", "BAG3", "UBC", "VEGFA", "MT3", "NEAT1", "MALAT1",
    "FOS", "JUN", "JUNB", "EGR1", "DUSP1", "DUSP2",
}
_HK_PREFIX = ("RPL", "RPS", "MRPL", "MRPS", "MT-", "MT1", "MT2")


def is_hk_stress(gene: str) -> bool:
    return gene in HK_STRESS or gene.startswith(_HK_PREFIX)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rescue", type=Path, required=True)
    p.add_argument("--hierarchy", type=Path, required=True)
    p.add_argument("--profiles", type=Path, default=None,
                   help="rescued_profiles.csv (genes x profiles) for the de-novo marker check.")
    p.add_argument("--cell-table", type=Path, default=None)
    p.add_argument("--lowsignal-label", default="Low_signal")
    p.add_argument("--prob-threshold", type=float, default=0.5,
                   help="Min rescue_prob to count an assignment as confident.")
    p.add_argument("--top-markers", type=int, default=8,
                   help="Top genes per de-novo profile inspected for the housekeeping/stress "
                        "fraction.")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def leaf_compartments(hierarchy_path: Path) -> dict[str, str]:
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


def classify_denovo_profiles(profiles_path: Path, named_types: set[str], top_k: int):
    """Score each de-novo profile (a column NOT in the named InSituTree set) for marker
    coherence. Returns (per-cluster DataFrame, {cluster: 'tiling'|'candidate'}). A cluster whose
    top markers are dominated (>=50%) by housekeeping/heat-shock/hypoxia genes is 'tiling'."""
    prof = pd.read_csv(profiles_path, index_col=0)
    denovo_cols = [c for c in prof.columns if c not in named_types]
    rows, cls = [], {}
    for c in denovo_cols:
        v = prof[c].astype(float)
        top = list(v.sort_values(ascending=False).head(top_k).index)
        hk_frac = float(np.mean([is_hk_stress(g) for g in top]))
        peak = float(v.max() / (v.mean() + 1e-9))
        label = "tiling" if hk_frac >= 0.5 else "candidate"
        cls[c] = label
        rows.append(dict(cluster=c, peakiness=peak, housekeeping_frac=hk_frac,
                         classification=label, top_markers=", ".join(top)))
    return pd.DataFrame(rows).sort_values("housekeeping_frac"), cls


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comp = leaf_compartments(args.hierarchy)
    named_types = set(comp)
    lowsignal_denovo = f"{args.lowsignal_label}_denovo"

    rescue = pd.read_csv(args.rescue)
    n = len(rescue)
    conf = rescue["rescue_prob"] >= args.prob_threshold
    is_denovo = rescue["is_denovo"].astype(bool)
    to_sink = rescue["rescue_label"].astype(str).isin([args.lowsignal_label, lowsignal_denovo])
    named_rescue = conf & ~is_denovo & ~to_sink
    denovo_rescue = conf & is_denovo & ~to_sink

    # ---- de-novo marker coherence: tiling vs candidate-real ---------------------------
    denovo_cls, coherence_df = {}, None
    if args.profiles is not None and args.profiles.exists():
        coherence_df, denovo_cls = classify_denovo_profiles(
            args.profiles, named_types, args.top_markers)
        coherence_df.to_csv(args.output_dir / "denovo_marker_coherence.csv", index=False)

    lbl = rescue["rescue_label"].astype(str)
    denovo_tiling = denovo_rescue & lbl.map(lambda x: denovo_cls.get(x) == "tiling").fillna(False)
    denovo_candidate = denovo_rescue & lbl.map(
        lambda x: denovo_cls.get(x) == "candidate").fillna(False)
    denovo_unclassified = denovo_rescue & ~denovo_tiling & ~denovo_candidate
    still_sink = ~(named_rescue | denovo_rescue)
    # RESOLUTION = a real reference identity OR a coherent new program (NOT stress/HK tiling)
    resolved = named_rescue | denovo_candidate

    lines = [f"Phase 3 — one more rescue on the Low_signal pool ({n:,} cells)\n"]
    lines.append("OUTCOME (de-novo assignment is NOT resolution — every cell gets a cluster; the")
    lines.append("meaningful metric is snapping to a NAMED type or a marker-COHERENT new program):")
    lines.append(f"  named type (real reference identity):     {int(named_rescue.sum()):>9,}  "
                 f"({named_rescue.mean():.1%})")
    if denovo_cls:
        lines.append(f"  de-novo, coherent candidate program:      "
                     f"{int(denovo_candidate.sum()):>9,}  ({denovo_candidate.mean():.1%})")
        lines.append(f"  de-novo, housekeeping/stress TILING:      "
                     f"{int(denovo_tiling.sum()):>9,}  ({denovo_tiling.mean():.1%})")
        if denovo_unclassified.any():
            lines.append(f"  de-novo, unclassified:                    "
                         f"{int(denovo_unclassified.sum()):>9,}  ({denovo_unclassified.mean():.1%})")
    else:
        lines.append(f"  de-novo cluster (UNCLASSIFIED — no profiles given): "
                     f"{int(denovo_rescue.sum()):>9,}  ({denovo_rescue.mean():.1%})")
    lines.append(f"  still flat / sink / low-confidence:       {int(still_sink.sum()):>9,}  "
                 f"({still_sink.mean():.1%})")

    raw_left = (named_rescue | denovo_rescue).mean()
    res_frac = resolved.mean()
    lines.append(f"\n  RAW 'left the sink' (incl. tiling): {raw_left:.1%}   |   GENUINE resolution "
                 f"(named + coherent de-novo): {res_frac:.1%}")
    if denovo_cls:
        n_tiling = sum(1 for v in denovo_cls.values() if v == "tiling")
        lines.append(f"  de-novo clusters: {n_tiling}/{len(denovo_cls)} are housekeeping/stress "
                     f"TILING (not new identities).")
    if res_frac > 0.33:
        lines.append("  => SUBSTANTIAL genuine resolution — the pool was still reference/off-"
                     "reference-limited; another native rescue iteration is worthwhile.")
    else:
        lines.append("  => the raw collapse is dominated by de-novo TILING of a flat/stressed "
                     "continuum; only a small fraction resolves to a real identity (see below). "
                     "The residue is an IRREDUCIBLE flat malignant continuum, not clean missable "
                     "cell types — iterating will not dissolve it.")

    # ---- named rescues by compartment ------------------------------------------------
    rescued_named = rescue[named_rescue].copy()
    rescued_named["compartment"] = rescued_named["rescue_label"].map(comp).fillna("other")
    by_comp = (rescued_named["compartment"].value_counts()
               .rename_axis("compartment").reset_index(name="n_cells"))
    by_comp["frac_of_pool"] = by_comp["n_cells"] / n
    by_comp.to_csv(args.output_dir / "rescue_by_compartment.csv", index=False)
    if len(by_comp):
        lines.append("\n  named rescues by compartment (fraction of the whole pool):")
        for _, r in by_comp.iterrows():
            lines.append(f"    {r['compartment']:<14} {int(r['n_cells']):>9,}  "
                         f"({r['frac_of_pool']:.1%})")

    if coherence_df is not None and len(coherence_df):
        lines.append("\n  de-novo cluster marker coherence (sorted; low housekeeping_frac = more "
                     "lineage-like):")
        for _, r in coherence_df.iterrows():
            lines.append(f"    {r['cluster']:<6} {r['classification']:<10} hk_frac="
                         f"{r['housekeeping_frac']:.2f}  peak={r['peakiness']:.0f}  "
                         f"[{r['top_markers']}]")

    # ---- CNV concordance (optional) --------------------------------------------------
    if args.cell_table is not None and args.cell_table.exists():
        malignant_types = set(DEFAULT_MALIGNANT)
        cells = pd.read_csv(args.cell_table, index_col=0)
        j = rescue.set_index("cell_id").join(cells[["is_malignant_call"]], how="inner")
        j["rescued_malignant"] = (j["rescue_prob"] >= args.prob_threshold) & \
            j["rescue_label"].astype(str).isin(malignant_types)
        j["cnv_malignant"] = j["is_malignant_call"].astype(bool)
        ct = pd.crosstab(j["rescued_malignant"], j["cnv_malignant"])
        ct.to_csv(args.output_dir / "rescue_vs_cnv.csv")
        both = int(((j["rescued_malignant"]) & (j["cnv_malignant"])).sum())
        n_resc_mal = int(j["rescued_malignant"].sum())
        n_cnv_mal = int(j["cnv_malignant"].sum())
        union = n_resc_mal + n_cnv_mal - both
        lines.append("\n  CNV concordance (Low_signal cells in the CNV table):")
        lines.append(f"    rescued->malignant AND CNV-malignant: {both:,}")
        lines.append(f"    rescued->malignant: {n_resc_mal:,}; CNV-malignant: {n_cnv_mal:,}; "
                     f"Jaccard {both / union:.0%}" if union else "")
        if n_resc_mal:
            lines.append(f"    => {both / n_resc_mal:.1%} of expression-rescued-malignant cells "
                         "are independently CNV-malignant (orthogonal axes, partial overlap).")

    # ---- plot: outcome stacked bar ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.8, 5))
    parts = [("named type", named_rescue.mean(), "#1a9850")]
    if denovo_cls:
        parts += [("de-novo: coherent", denovo_candidate.mean(), "#66bd63"),
                  ("de-novo: HK/stress tiling", denovo_tiling.mean(), "#fdae61")]
        if denovo_unclassified.any():
            parts += [("de-novo: unclassified", denovo_unclassified.mean(), "#fee08b")]
    else:
        parts += [("de-novo (unclassified)", denovo_rescue.mean(), "#fdae61")]
    parts += [("still flat/sink", still_sink.mean(), "#bdbdbd")]
    bottom = 0.0
    for label, frac, color in parts:
        ax.bar(0, frac, bottom=bottom, color=color, label=f"{label} ({frac:.0%})")
        bottom += frac
    ax.set_xlim(-1, 1.8); ax.set_xticks([])
    ax.set_ylabel("fraction of the Low_signal pool")
    ax.set_title("Phase 3: rescue outcome\n(de-novo tiling is NOT resolution)")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output_dir / "rescue_rate.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    (args.output_dir / "rescue_summary.txt").write_text("\n".join(l for l in lines if l) + "\n")
    print("\n" + "\n".join(l for l in lines if l))
    print(f"\nWrote summary + tables + plot to {args.output_dir}")


if __name__ == "__main__":
    main()
