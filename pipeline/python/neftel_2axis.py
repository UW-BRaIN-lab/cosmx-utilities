#!/usr/bin/env python3
"""Neftel 2-axis (butterfly) scoring of the GBM cohort — a CONTINUOUS malignant-state map.

Complements the discrete top-two hybrid readout (R/flat_posteriors.R): those peaked posteriors
flag only cells the classifier itself cannot separate (a conservative lower bound on hybrids).
This places every cell on the Neftel et al. (2019, Cell) malignant-state plane, so gradient cells
between corners are visible as continuous position rather than a hard call.

Per cell, score the four Neftel-like modules (scanpy sc.tl.score_genes) then, following Neftel:
  D  = max(SC_OPC, SC_NPC) - max(SC_AC, SC_MES)
  y  = sign(D) * log2(|D| + 1)                  # OPC/NPC-like (top) vs AC/MES-like (bottom)
  dx = (SC_OPC - SC_NPC) if D > 0 else (SC_AC - SC_MES)
  x  = sign(dx) * log2(|dx| + 1)                # within the winning pair
Corners: top-right OPC, top-left NPC, bottom-right AC, bottom-left MES. Cells near the origin are
intermediate / between-state.

Modules come from pipeline/reference/gene_signatures.csv (AClike/MESlike/OPClike/NPClike), already
panel-restricted. Malignant + Low_signal cell sets are read from insitutree_hierarchy.json so this
tracks the typing.

Inputs:
  --typed-h5ad   cosmx_typed.h5ad (obs['cell_type']; raw counts in layers['counts'] or X)
  --signatures   gene_signatures.csv (module,gene)
  --hierarchy    insitutree_hierarchy.json (to identify the malignant states + the sink)
Outputs (--output-dir):
  neftel_coords.csv.gz     per-cell cell_id, cell_type, SC_AC/MES/OPC/NPC, neftel_x, neftel_y
  butterfly_state.png      malignant cells, coloured by InSituTree meta-state (corners occupied?)
  butterfly_lowsignal.png  malignant density (grey) + Low_signal overlay (does the sink sit central?)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MODULES = {"AClike": "AC", "MESlike": "MES", "OPClike": "OPC", "NPClike": "NPC"}
# malignant leaf -> Neftel meta-state (for colouring the state butterfly)
META = {
    "AC-like": "AC", "AC-like_Prolif": "AC",
    "MES-like_hypoxia_independent": "MES", "MES-like_hypoxia_MHC": "MES",
    "MESlike_denovo": "MES", "MES_AClike_denovo": "MES", "Hypoxia_denovo": "MES",
    "NPC-like_OPC": "NPC", "NPC-like_Prolif": "NPC", "NPC-like_neural": "NPC",
    "OPC-like": "OPC", "OPC-like_Prolif": "OPC", "OPClike_denovo": "OPC",
    "Stress_sig": "Stress", "Stress_denovo": "Stress",
}
STATE_COLOR = {"AC": "#d1495b", "MES": "#e3873c", "NPC": "#3a7ca5", "OPC": "#2a9d8f",
               "Stress": "#8a8d93"}
INTERMEDIATE_R = 1.0   # |x|,|y| within this radius = "intermediate / between-state"
PLOT_MAX = 200_000


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    p.add_argument("--signatures", type=Path, required=True)
    p.add_argument("--hierarchy", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--low-signal-label", default="Low_signal")
    return p.parse_args()


def load_modules(path: Path) -> dict[str, list[str]]:
    rows = [l for l in path.read_text().splitlines()
            if l.strip() and not l.lstrip().startswith("#")]
    import csv, io
    r = csv.DictReader(io.StringIO("\n".join(rows)))
    mods: dict[str, list[str]] = {}
    for row in r:
        mods.setdefault(row["module"], []).append(row["gene"])
    return mods


def ensure_lognorm(adata: ad.AnnData) -> None:
    """Put log-normalised expression in adata.X, from raw counts. Prefer layers['counts'];
    otherwise decide from X: integer-valued => raw (normalise); else assume already log."""
    if "counts" in adata.layers:
        print("Using layers['counts'] as the raw source.")
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return
    X = adata.X
    sample = X[:2000]
    data = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
    data = data[np.isfinite(data)]
    is_int = data.size > 0 and np.allclose(data, np.round(data))
    xmax = float(data.max()) if data.size else 0.0
    print(f"X sample: max={xmax:.3f}, integer-valued={is_int}")
    if is_int:
        print("X looks like raw counts -> normalize_total + log1p.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("X is non-integer -> assuming already log-normalised; using as-is. "
              "(If corners look collapsed, re-run from a raw-counts source.)")


def neftel_coords(sc_ac, sc_mes, sc_opc, sc_npc):
    D = np.maximum(sc_opc, sc_npc) - np.maximum(sc_ac, sc_mes)
    y = np.sign(D) * np.log2(np.abs(D) + 1)
    dx = np.where(D > 0, sc_opc - sc_npc, sc_ac - sc_mes)
    x = np.sign(dx) * np.log2(np.abs(dx) + 1)
    return x, y


def _subsample(n: int) -> np.ndarray:
    rng = np.random.default_rng(0)
    return rng.choice(n, size=PLOT_MAX, replace=False) if n > PLOT_MAX else np.arange(n)


def _corners(ax):
    for (xx, yy, lab, col) in [(1, 1, "OPC", "#2a9d8f"), (-1, 1, "NPC", "#3a7ca5"),
                               (1, -1, "AC", "#d1495b"), (-1, -1, "MES", "#e3873c")]:
        ax.text(xx, yy, lab, ha="center", va="center", fontsize=13, fontweight="bold",
                color=col, transform=ax.transData)
    ax.axhline(0, color="#999", lw=.6, zorder=0)
    ax.axvline(0, color="#999", lw=.6, zorder=0)


def butterfly_state(x, y, meta, out_path):
    idx = _subsample(len(x))
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    for state, col in STATE_COLOR.items():
        m = meta[idx] == state
        if m.any():
            ax.scatter(x[idx][m], y[idx][m], s=2, c=col, linewidths=0, alpha=.5, label=state)
    _corners(ax)
    lim = np.nanpercentile(np.abs(np.concatenate([x, y])), 99.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("within-pair  (NPC ← → OPC   |   MES ← → AC)")
    ax.set_ylabel("AC/MES-like  ← → OPC/NPC-like")
    ax.set_title(f"Neftel 2-axis map — malignant cells by state (n={len(x):,})")
    ax.legend(markerscale=4, loc="upper left", fontsize=8, framealpha=.9)
    fig.tight_layout(); fig.savefig(out_path, dpi=150); plt.close(fig)


def butterfly_lowsignal(xm, ym, xl, yl, out_path, low_label):
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    lim = np.nanpercentile(np.abs(np.concatenate([xm, ym])), 99.5)
    ax.hexbin(xm, ym, gridsize=80, cmap="Greys", bins="log", extent=(-lim, lim, -lim, lim))
    li = _subsample(len(xl))
    ax.scatter(xl[li], yl[li], s=2, c="#c0392b", linewidths=0, alpha=.35, label=low_label)
    _corners(ax)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("within-pair  (NPC ← → OPC   |   MES ← → AC)")
    ax.set_ylabel("AC/MES-like  ← → OPC/NPC-like")
    ax.set_title(f"Malignant density (grey) + {low_label} overlay (red)")
    ax.legend(markerscale=4, loc="upper left", fontsize=8, framealpha=.9)
    fig.tight_layout(); fig.savefig(out_path, dpi=150); plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    hier = json.loads(args.hierarchy.read_text())
    malignant = set(hier.get("Malignant", []))
    print(f"{len(malignant)} malignant states; sink label = {args.low_signal_label}")

    print(f"Reading {args.typed_h5ad}")
    adata = ad.read_h5ad(args.typed_h5ad)
    if args.celltype_key not in adata.obs:
        print(f"ERROR: obs has no '{args.celltype_key}'", file=sys.stderr); sys.exit(1)
    ensure_lognorm(adata)

    mods = load_modules(args.signatures)
    scores = {}
    for module, state in MODULES.items():
        genes = [g for g in mods.get(module, []) if g in adata.var_names]
        if len(genes) < 3:
            print(f"ERROR: module {module} has <3 genes on the panel", file=sys.stderr); sys.exit(1)
        sc.tl.score_genes(adata, genes, score_name=f"SC_{state}", use_raw=False)
        scores[state] = adata.obs[f"SC_{state}"].to_numpy()
        print(f"  scored {module} -> SC_{state} ({len(genes)} panel genes)")

    x, y = neftel_coords(scores["AC"], scores["MES"], scores["OPC"], scores["NPC"])
    ct = adata.obs[args.celltype_key].astype(str).to_numpy()

    out = pd.DataFrame({
        "cell_id": adata.obs_names.to_numpy(), "cell_type": ct,
        "SC_AC": scores["AC"], "SC_MES": scores["MES"],
        "SC_OPC": scores["OPC"], "SC_NPC": scores["NPC"],
        "neftel_x": x, "neftel_y": y,
    })
    coords_path = args.output_dir / "neftel_coords.csv.gz"
    out.to_csv(coords_path, index=False)
    print(f"Wrote {coords_path} ({len(out):,} cells)")

    is_mal = np.isin(ct, list(malignant))
    is_low = ct == args.low_signal_label
    print(f"malignant cells: {int(is_mal.sum()):,} | {args.low_signal_label}: {int(is_low.sum()):,}")

    meta = pd.Series(ct).map(META).to_numpy()
    if is_mal.any():
        butterfly_state(x[is_mal], y[is_mal], meta[is_mal],
                        args.output_dir / "butterfly_state.png")
        # continuous analog of the hybrid fraction: malignant cells near the origin
        inter = is_mal & (np.abs(x) < INTERMEDIATE_R) & (np.abs(y) < INTERMEDIATE_R)
        print(f"intermediate (|x|,|y|<{INTERMEDIATE_R}) among malignant: "
              f"{int(inter.sum()):,} ({100*inter.sum()/is_mal.sum():.1f}%)")
    if is_mal.any() and is_low.any():
        butterfly_lowsignal(x[is_mal], y[is_mal], x[is_low], y[is_low],
                            args.output_dir / "butterfly_lowsignal.png", args.low_signal_label)
    print(f"Done. Neftel 2-axis map in {args.output_dir}")


if __name__ == "__main__":
    main()
