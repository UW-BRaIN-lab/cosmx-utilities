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
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Patch, Rectangle
from matplotlib.path import Path as MPath

DBLUE, TEAL, RED, GRAY = "#185FA5", "#1D9E75", "#E24B4A", "#888780"
CAT_COL = {"denovo": DBLUE, "named": TEAL, "artifact": RED}
LINEAGE_COL = {"Tumor": "#D85A30", "Glia": "#1D9E75", "Neuronal": "#7F77DD",
               "Myeloid": "#BA7517", "Lymphoid": "#D4537E", "Vascular": "#378ADD",
               "Unresolved/stress": "#888780"}

RUN_LABELS = ["Core L4\nrescale", "Core L4\nrefit", "Ext L3\nrefit",
              "Ext L3 rescale\n(keeper)"]

# 4-segment composition, % of 2.33M cells (named-mapped, named-artifact, denovo-interp, denovo-lowsig)
COMPOSITION = {
    "Named — mapped to GBmap":      ([10.4, 40.7, 27.7, 13.9], TEAL),
    "Named — rare-type artifact":   ([0.0, 19.6, 27.7, 0.0], RED),
    "De novo — interpretable type": ([38.6, 36.0, 16.6, 33.1], DBLUE),
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


# run display name -> (Kopah dir / counts-CSV stem, annotation table, artifact named-types)
RUNS_META = [
    ("Core L4 · rescale", "stage4", "stage4.csv", set()),
    ("Core L4 · refit", "stage4_refit", "stage4_refit.csv", {"Reg_T", "DC3"}),
    ("Ext L3 · refit", "stage4_ext_l3", "stage4_ext_l3.csv", {"B_cell", "RG", "Mast"}),
    ("Ext L3 · rescale — keeper", "stage4_extl3_rescale", "stage4_extl3_rescale.csv", set()),
]


def _categorize(label: str, artifacts: set[str]) -> str:
    if len(label) == 1 and label.islower():
        return "denovo"
    if label in artifacts:
        return "artifact"
    return "named"


def build_rows_from_counts(counts_dir: Path, anno_dir: Path | None):
    """Per run: read counts_<dir>.csv (cell_type,count) -> sorted (label, count, cat) rows,
    relabeling de-novo letters via the run's annotation table when anno_dir is given."""
    rows_by_run = {}
    for name, stem, anno_file, artifacts in RUNS_META:
        cpath = counts_dir / f"counts_{stem}.csv"
        if not cpath.exists():
            print(f"WARN: missing {cpath}; skipping {name}")
            continue
        df = pd.read_csv(cpath, index_col=0)
        counts = df.iloc[:, 0]
        anno = {}
        if anno_dir is not None and (anno_dir / anno_file).exists():
            a = pd.read_csv(anno_dir / anno_file)
            anno = dict(zip(a["denovo_label"].astype(str), a["annotation"].astype(str)))
        rows = []
        for lab, cnt in counts.sort_values(ascending=False).items():
            lab = str(lab)
            cat = _categorize(lab, artifacts)
            disp = anno.get(lab, lab) if cat == "denovo" else lab
            rows.append((disp, int(cnt), cat))
        rows_by_run[name] = rows
    return rows_by_run


def fig_breakdown(out: Path, rows_by_run: dict, top_n: int | None, suptitle: str,
                  show_counts: bool = True) -> None:
    maxbars = max(len(r if top_n is None else r[:top_n]) for r in rows_by_run.values())
    fig, axes = plt.subplots(2, 2, figsize=(16, max(9, maxbars * 0.22 * 2)))
    fontsize = 8 if maxbars <= 14 else 6
    for ax, (name, rows) in zip(axes.ravel(), rows_by_run.items()):
        show = (rows if top_n is None else rows[:top_n])[::-1]  # largest at top after barh
        cols = [CAT_COL[r[2]] for r in show]
        ax.barh(range(len(show)), [r[1] for r in show], color=cols)
        ax.set_yticks(range(len(show))); ax.set_yticklabels([r[0] for r in show], fontsize=fontsize)
        ax.set_ylim(-0.6, len(show) - 0.4)
        if show_counts:
            ndn = sum(r[2] == "denovo" for r in rows)
            ax.set_title(f"{name}  ({len(rows)} types, {ndn} de novo)", fontsize=11)
        else:
            ax.set_title(name, fontsize=11)
        ax.set_xlabel("cells")
        ax.xaxis.set_major_formatter(lambda v, _: f"{v/1000:.0f}k" if v else "0")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.ravel()[len(rows_by_run):]:
        ax.set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, color=c) for c in (DBLUE, TEAL, RED)]
    fig.legend(handles, ["De novo cluster", "Named (GBmap)", "Named — rare-type artifact"],
               loc="upper center", ncol=3, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(suptitle, fontsize=14, y=1.02)
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


# ---------------------------------------------------------------------------
# De-novo-resolved lineage Sankey: named GBmap types collapse into coarse
# lineages, but every de-novo cluster (a-l) gets its own labelled node instead
# of vanishing into a single Unresolved/stress bar. Flows are computed from the
# real Core-L4-rescale vs Ext-L3-rescale cluster cross-tab, so nothing is
# hand-transcribed. De-novo nodes are coloured by their interpreted lineage and
# hatched so they read distinctly from the named-lineage nodes.

LINEAGE_RULES = [  # (regex on named label, lineage); first match wins
    (r"AC-like|MES-like|NPC-like|OPC-like", "Tumor"),
    (r"^(Astrocyte|Oligodendrocyte|RG|OPC)$", "Glia"),
    (r"^Neuron$", "Neuronal"),
    (r"^(TAM-|Mono|DC\d|cDC|pDC|Mast|Neutrophil|DC$)", "Myeloid"),
    (r"^(B_cell|Plasma_B|CD4|CD8|NK|Reg_T|Prolif_T)", "Lymphoid"),
    (r"^(Endo|Pericyte|Perivascular|SMC|VLMC|Mural|Scavenging|Tip-like)", "Vascular"),
]


def named_lineage(label: str) -> str:
    for pat, lin in LINEAGE_RULES:
        if re.search(pat, label):
            return lin
    return "Unresolved/stress"


def annotation_lineage(annotation: str) -> str:
    """Coarse lineage for a de-novo cluster from its annotation text (tumour wins
    over a co-occurring stress/hypoxia term, e.g. 'Stressed/hypoxic tumor-bulk')."""
    a = annotation.lower()
    if re.search(r"tumor|tumour|mes|ac-like|opc-like|stem", a):
        return "Tumor"
    if re.search(r"astrocyte|oligodendrocyte|\brg\b|\bopc\b", a):
        return "Glia"
    if "neuron" in a:
        return "Neuronal"
    if re.search(r"myeloid|tam", a):
        return "Myeloid"
    if re.search(r"vascular|fibroblast|angiogenic|stroma|endothel", a):
        return "Vascular"
    return "Unresolved/stress"


def is_denovo(label: str) -> bool:
    return len(label) == 1 and label.islower()


def load_anno(path: Path) -> dict[str, str]:
    """letter -> display label ('k · Stressed/hypoxic tumor-bulk')."""
    df = pd.read_csv(path)
    out = {}
    for letter, anno in zip(df["denovo_label"].astype(str), df["annotation"].astype(str)):
        out[letter] = anno.replace(" - ", " · ")
    return out


def _denovo_node(letter: str, anno: dict[str, str]) -> tuple[str, str]:
    """Return (display label, interpreted lineage) for a de-novo letter."""
    disp = anno.get(letter, letter)
    lin = annotation_lineage(disp)
    return disp, lin


def fig_sankey_denovo(crosstab: Path, left_anno: Path, right_anno: Path,
                      out: Path, min_frac: float = 0.0012) -> None:
    ct = pd.read_csv(crosstab, index_col=0)
    lanno, ranno = load_anno(left_anno), load_anno(right_anno)

    def side_node(label, anno):  # -> (display name, lineage, is_denovo)
        if is_denovo(label):
            disp, lin = _denovo_node(label, anno)
            return disp, lin, True
        return named_lineage(label), named_lineage(label), False

    # Aggregate the fine cross-tab onto (left node, right node).
    flows: dict[tuple[str, str], float] = {}
    lmeta: dict[str, tuple[str, bool]] = {}  # node -> (lineage, is_denovo)
    rmeta: dict[str, tuple[str, bool]] = {}
    for r in ct.index:
        ln, llin, ldn = side_node(r, lanno)
        lmeta[ln] = (llin, ldn)
        for c in ct.columns:
            v = float(ct.loc[r, c])
            if v <= 0:
                continue
            rn, rlin, rdn = side_node(c, ranno)
            rmeta[rn] = (rlin, rdn)
            flows[(ln, rn)] = flows.get((ln, rn), 0.0) + v

    ltot = {n: sum(v for (s, _), v in flows.items() if s == n) for n in lmeta}
    rtot = {n: sum(v for (_, d), v in flows.items() if d == n) for n in rmeta}

    def order(meta, tot):  # named lineage first, then its de-novo nodes (size desc)
        nodes = []
        for lin in LIN_ORDER:
            named = [n for n, (l, dn) in meta.items() if l == lin and not dn and tot[n] > 0]
            denovo = sorted((n for n, (l, dn) in meta.items()
                             if l == lin and dn and tot[n] > 0), key=lambda n: -tot[n])
            nodes += sorted(named, key=lambda n: -tot[n]) + denovo
        return nodes

    L, R = order(lmeta, ltot), order(rmeta, rtot)
    total = sum(ltot.values())
    gap = total * 0.012
    thresh = min_frac * total

    def layout(nodes, tot):
        y, pos = 0.0, {}
        for n in nodes:
            pos[n] = (y, tot[n]); y += tot[n] + gap
        return pos, y

    Lp, yL = layout(L, ltot); Rp, yR = layout(R, rtot)
    xL0, xL1, xR0, xR1 = 0.0, 0.022, 0.978, 1.0
    xm = (xL1 + xR0) / 2
    fig, ax = plt.subplots(figsize=(13, max(10, 0.34 * max(len(L), len(R)) + 2)))
    Loff = {n: 0.0 for n in L}; Roff = {n: 0.0 for n in R}
    for s in L:                                   # ribbons big-first per source
        for d in sorted(R, key=lambda d: -flows.get((s, d), 0)):
            f = flows.get((s, d), 0)
            if f <= 0:
                continue
            y1, y2 = Lp[s][0] + Loff[s], Rp[d][0] + Roff[d]
            Loff[s] += f; Roff[d] += f
            if f < thresh:
                continue
            verts = [(xL1, y1), (xm, y1), (xm, y2), (xR0, y2),
                     (xR0, y2 + f), (xm, y2 + f), (xm, y1 + f), (xL1, y1 + f), (xL1, y1)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                     MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), facecolor=LINEAGE_COL[lmeta[s][0]],
                                   edgecolor="none", alpha=0.40))

    def draw_nodes(nodes, pos, meta, x0, x1, txt_x, ha):
        for n in nodes:
            lin, dn = meta[n]
            ax.add_patch(Rectangle((x0, pos[n][0]), x1 - x0, pos[n][1],
                                   facecolor=LINEAGE_COL[lin], edgecolor="white",
                                   linewidth=0.4, hatch="////" if dn else None))
            ax.text(txt_x, pos[n][0] + pos[n][1] / 2, n, ha=ha, va="center", fontsize=7.5)

    draw_nodes(L, Lp, lmeta, xL0, xL1, xL0 - 0.01, "right")
    draw_nodes(R, Rp, rmeta, xR0, xR1, xR1 + 0.01, "left")
    ax.set_xlim(-0.34, 1.34); ax.set_ylim(0, max(yL, yR)); ax.invert_yaxis(); ax.axis("off")
    ax.text(xL1, -gap * 1.5, "Core L4 · rescale", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text(xR0, -gap * 1.5, "Ext L3 · rescale (keeper)", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.set_title(f"Lineage flow between the two rescale runs — de novo clusters broken out  "
                 f"({int(total):,} cells)", fontsize=12.5, pad=22)
    leg = [Patch(facecolor=LINEAGE_COL[l], label=l) for l in LIN_ORDER] + \
          [Patch(facecolor="white", edgecolor="#444", hatch="////", label="de novo cluster")]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=4, fontsize=8, frameon=False)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}  ({len(L)} left x {len(R)} right nodes)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", type=Path, default=Path("stage4_qc/figures"))
    p.add_argument("--crosstab", type=Path,
                   default=Path("stage4_qc/core_vs_ext/run_comparison_cluster_xtab.csv"),
                   help="Core-L4-rescale (rows) vs Ext-L3-rescale (cols) cluster cross-tab; "
                        "drives the de-novo-resolved lineage Sankey.")
    p.add_argument("--counts-dir", type=Path, default=None,
                   help="Dir with per-run counts_<run>.csv (cell_type,count) — dump via "
                        "obs value_counts. Enables the all-types breakdown + count-accurate "
                        "top-12. Without it, the top-12 uses the embedded summary.")
    p.add_argument("--anno-dir", type=Path,
                   default=Path("pipeline/reference/denovo_annotations"),
                   help="denovo_annotations dir, to relabel de-novo letters in the "
                        "count-CSV breakdowns.")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_composition(args.output_dir / "composition_stacked.png")
    fig_sankey(args.output_dir / "lineage_sankey.png")
    if args.crosstab.exists():
        fig_sankey_denovo(
            args.crosstab,
            args.anno_dir / "stage4.csv",
            args.anno_dir / "stage4_extl3_rescale.csv",
            args.output_dir / "lineage_denovo_sankey.png",
        )
    else:
        print(f"WARN: {args.crosstab} missing; skipping de-novo-resolved Sankey")

    if args.counts_dir is not None:
        rows = build_rows_from_counts(args.counts_dir, args.anno_dir)
        fig_breakdown(args.output_dir / "per_run_breakdown.png", rows, top_n=12,
                      suptitle="Top 12 cell types per run", show_counts=True)
        fig_breakdown(args.output_dir / "per_run_breakdown_all.png", rows, top_n=None,
                      suptitle="All cell types per run", show_counts=True)
        print(f"Wrote 4 figures (incl. all-types) to {args.output_dir}")
    else:
        fig_breakdown(args.output_dir / "per_run_breakdown.png", BREAKDOWN, top_n=12,
                      suptitle="Top 12 cell types per run", show_counts=False)
        print(f"Wrote 3 figures to {args.output_dir} "
              f"(pass --counts-dir for the all-types breakdown)")


if __name__ == "__main__":
    main()
