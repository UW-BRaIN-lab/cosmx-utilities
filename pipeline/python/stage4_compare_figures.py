#!/usr/bin/env python3
"""Stage 4 lab-meeting comparison figures (PNG) for the four Wenyu InSituType runs.

Generates three slide-ready PNGs:
  composition_stacked.png  - % named vs de-novo per run, with the refit rare-type
                             artifact slab broken out.
  per_run_breakdown.png    - top-12 cell types per run, de-novo clusters highlighted.
  lineage_sankey.png       - coarse-lineage flow, Core-L4-rescale -> Extended-L3-rescale.

Summary numbers are embedded with provenance: each run's stage-4b assignment-count log
(stage4 / stage4_refit / stage4_ext_l3 / stage4_extl3_rescale) and the
comparisons/stage4_vs_stage4_extl3_rescale lineage cross-tab. De-novo clusters render in
pink (#D4537E) so they read distinctly from the red rare-type-artifact bars.

Run anywhere with matplotlib:
    uv run python pipeline/python/stage4_compare_figures.py --output-dir stage4_qc/figures
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

PINK, TEAL, RED, GRAY = "#D4537E", "#1D9E75", "#E24B4A", "#888780"
CAT_COL = {"denovo": PINK, "named": TEAL, "artifact": RED}
LINEAGE_COL = {"Tumor": "#D85A30", "Glia": "#1D9E75", "Neuronal": "#7F77DD",
               "Myeloid": "#BA7517", "Lymphoid": "#D4537E", "Vascular": "#378ADD",
               "Unresolved/stress": "#888780"}

RUN_LABELS = ["Core L4\nrescale", "Core L4\nrefit", "Ext L3\nrefit",
              "Ext L3 rescale\n(keeper)"]

# 4-segment composition, % of 2.33M cells (named-mapped, named-artifact, denovo-interp, denovo-lowsig)
COMPOSITION = {
    "Named — mapped to GBmap":      ([10.4, 40.7, 27.7, 13.9], TEAL),
    "Named — rare-type artifact":   ([0.0, 19.6, 27.7, 0.0], RED),
    "De novo — interpretable type": ([38.6, 36.0, 16.6, 33.1], PINK),
    "De novo — low-signal / stress":([51.0, 3.7, 27.9, 53.0], GRAY),
}

# top-12 cell types per run: (label, count, category)
BREAKDOWN = {
    "Core L4 · rescale  (66 types)": [
        ("j — Low-signal/generic", 931427, "denovo"), ("l — Oligodendrocyte", 227055, "denovo"),
        ("f — Hypoxia/angiogenic", 170275, "denovo"), ("b — Heat-shock/stress", 152244, "denovo"),
        ("k — MES-like", 147772, "denovo"), ("h — OPC-like tumor", 138557, "denovo"),
        ("i — Neuronal", 109685, "denovo"), ("e — Astrocyte", 90657, "denovo"),
        ("d — MES/AC-like tumor", 44874, "denovo"), ("a — Myeloid/TAM", 35424, "denovo"),
        ("g — Neuronal", 26987, "denovo"), ("c — Vascular/stroma", 14304, "denovo")],
    "Core L4 · refit  (63 types)": [
        ("Reg_T", 326992, "artifact"), ("k — Heat-shock/stress", 284986, "denovo"),
        ("l — Oligodendrocyte", 142209, "denovo"), ("c — OPC-like tumor", 133356, "denovo"),
        ("AC-like_Prolif", 132450, "named"), ("DC3", 129582, "artifact"),
        ("OPC-like_Prolif", 119220, "named"), ("d — Hypoxia/MES-like", 107318, "denovo"),
        ("AC-like", 78191, "named"), ("CD8_EM", 78185, "named"),
        ("b — Neuronal", 70969, "denovo"), ("f — Mixed (myeloid/MES)", 56050, "denovo")],
    "Ext L3 · refit  (33 types)": [
        ("B_cell", 362572, "artifact"), ("k — MES/AC-like tumor", 289617, "denovo"),
        ("RG", 242807, "artifact"), ("AC-like", 159360, "named"),
        ("d — Astrocyte+heat-shock", 155141, "denovo"), ("a — OPC-like tumor", 141101, "denovo"),
        ("Oligodendrocyte", 131967, "named"), ("i — Hypoxia/collagen", 98324, "denovo"),
        ("b — Neuronal", 87014, "denovo"), ("l — OPC-like/stem tumor", 69751, "denovo"),
        ("c — Mixed/low-signal", 69110, "denovo"), ("OPC-like", 57160, "named")],
    "Ext L3 · rescale  (33 types) — keeper": [
        ("k — Stressed/hypoxic tumor", 831946, "denovo"), ("g — Low-signal/generic", 188159, "denovo"),
        ("h — Metallothionein", 180861, "denovo"), ("f — OPC-like tumor", 144076, "denovo"),
        ("l — Astrocyte", 131689, "denovo"), ("e — Hypoxia/angiogenic", 130029, "denovo"),
        ("Oligodendrocyte", 127680, "named"), ("d — Neuronal", 106759, "denovo"),
        ("b — Myeloid/TAM", 92190, "denovo"), ("c — OPC-like/stem tumor", 87781, "denovo"),
        ("OPC-like", 66568, "named"), ("i — Vascular/fibroblast", 51584, "denovo")],
}

LIN_ORDER = ["Tumor", "Glia", "Neuronal", "Myeloid", "Lymphoid", "Vascular", "Unresolved/stress"]
LIN_LEFT = {"Tumor": 605, "Glia": 331, "Neuronal": 139, "Myeloid": 108, "Lymphoid": 7,
            "Vascular": 58, "Unresolved/stress": 1084}   # Core-L4-rescale, thousands
LIN_RIGHT = {"Tumor": 478, "Glia": 267, "Neuronal": 143, "Myeloid": 144, "Lymphoid": 8,
             "Vascular": 91, "Unresolved/stress": 1201}  # Ext-L3-rescale, thousands
LIN_FLOW = {  # src -> dst, thousands
    "Tumor": {"Tumor": 405, "Glia": 6, "Myeloid": 1, "Vascular": 3, "Unresolved/stress": 190},
    "Glia": {"Tumor": 7, "Glia": 235, "Neuronal": 5, "Myeloid": 2, "Unresolved/stress": 81},
    "Neuronal": {"Tumor": 3, "Neuronal": 135, "Unresolved/stress": 1},
    "Myeloid": {"Tumor": 3, "Myeloid": 103, "Unresolved/stress": 2},
    "Lymphoid": {"Lymphoid": 6, "Unresolved/stress": 1},
    "Vascular": {"Tumor": 1, "Glia": 1, "Myeloid": 1, "Vascular": 53, "Unresolved/stress": 2},
    "Unresolved/stress": {"Tumor": 59, "Glia": 24, "Neuronal": 2, "Myeloid": 38,
                          "Lymphoid": 2, "Vascular": 34, "Unresolved/stress": 925},
}


def fig_composition(out: Path) -> None:
    x = np.arange(4)
    fig, ax = plt.subplots(figsize=(8, 5.2))
    bottom = np.zeros(4)
    for label, (vals, col) in COMPOSITION.items():
        vals = np.array(vals)
        ax.bar(x, vals, bottom=bottom, color=col, label=label, width=0.62,
               edgecolor="white", linewidth=0.7)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 4:
                ax.text(xi, b + v / 2, f"{round(v)}%", ha="center", va="center",
                        color="white", fontsize=9)
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(RUN_LABELS, fontsize=10)
    ax.set_ylim(0, 100); ax.set_ylabel("% of 2.33M cells")
    ax.set_title("Cell-type composition across InSituType runs", fontsize=13)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=9, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)


def fig_breakdown(out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    for ax, (name, rows) in zip(axes.ravel(), BREAKDOWN.items()):
        rows = rows[::-1]  # largest at top after barh
        counts = [r[1] for r in rows]
        cols = [CAT_COL[r[2]] for r in rows]
        ax.barh(range(len(rows)), counts, color=cols)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("cells")
        ax.xaxis.set_major_formatter(lambda v, _: f"{v/1000:.0f}k" if v else "0")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, color=c) for c in (PINK, TEAL, RED)]
    fig.legend(handles, ["De novo cluster", "Named (GBmap)", "Named — rare-type artifact"],
               loc="upper center", ncol=3, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Top 12 cell types per run", fontsize=14, y=1.06)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)


def fig_sankey(out: Path) -> None:
    total = sum(LIN_LEFT.values())
    gap = total * 0.02

    def layout(tot):
        y, pos = 0.0, {}
        for n in LIN_ORDER:
            pos[n] = (y, tot[n]); y += tot[n] + gap
        return pos, y

    Lp, yL = layout(LIN_LEFT); Rp, yR = layout(LIN_RIGHT)
    xL0, xL1, xR0, xR1 = 0.0, 0.03, 0.97, 1.0
    xm = (xL1 + xR0) / 2
    fig, ax = plt.subplots(figsize=(10, 6.2))
    Loff = {n: 0.0 for n in LIN_ORDER}; Roff = {n: 0.0 for n in LIN_ORDER}
    for s in LIN_ORDER:
        for d in LIN_ORDER:
            f = LIN_FLOW.get(s, {}).get(d, 0)
            if not f:
                continue
            y1 = Lp[s][0] + Loff[s]; y2 = Rp[d][0] + Roff[d]; Loff[s] += f; Roff[d] += f
            verts = [(xL1, y1), (xm, y1), (xm, y2), (xR0, y2),
                     (xR0, y2 + f), (xm, y2 + f), (xm, y1 + f), (xL1, y1 + f), (xL1, y1)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                     MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), facecolor=LINEAGE_COL[s],
                                   edgecolor="none", alpha=0.32))
    for n in LIN_ORDER:
        ax.add_patch(Rectangle((xL0, Lp[n][0]), xL1 - xL0, Lp[n][1], color=LINEAGE_COL[n]))
        ax.add_patch(Rectangle((xR0, Rp[n][0]), xR1 - xR0, Rp[n][1], color=LINEAGE_COL[n]))
        ax.text(xL0 - 0.012, Lp[n][0] + Lp[n][1] / 2, n, ha="right", va="center", fontsize=9)
        ax.text(xR1 + 0.012, Rp[n][0] + Rp[n][1] / 2, n, ha="left", va="center", fontsize=9)
    ax.set_xlim(-0.18, 1.18); ax.set_ylim(0, max(yL, yR)); ax.invert_yaxis(); ax.axis("off")
    ax.text(xL1, -gap * 1.5, "Core L4 · rescale", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text(xR0, -gap * 1.5, "Ext L3 · rescale (keeper)", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.set_title("Lineage flow between the two rescale runs", fontsize=13, pad=24)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", type=Path, default=Path("stage4_qc/figures"))
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_composition(args.output_dir / "composition_stacked.png")
    fig_breakdown(args.output_dir / "per_run_breakdown.png")
    fig_sankey(args.output_dir / "lineage_sankey.png")
    print(f"Wrote 3 figures to {args.output_dir}")


if __name__ == "__main__":
    main()
