#!/usr/bin/env python3
"""Compare an external cell-typing against ours, cell-by-cell.

Built to compare a collaborator's InSituType/typing of the same CosMx cohort against our
keeper run, but general to any two per-cell labelings that share the underlying cells.
Aligns cells by a canonical (slide, fov, cell) key parsed from each side's IDs (CosMx
IDs differ only in formatting — ours `<slide>_F<fov>_C<cell>`, others often
`<slide>_<fov>_<cell>`), then on the shared cells reports:

  - adjusted Rand index (whole-labeling agreement)
  - per external-type -> top of-ours mapping (where each of their types lands in ours)
  - a row-normalized cross-tab heatmap (their types x ours), ordered by abundance
  - the raw cross-tab + the agreement table as CSVs

Inputs:
  --external-h5ad   their typed AnnData (obs has --external-key); read backed (obs only).
  --ours-csv        our per-cell labels CSV: index = our cell id, one column of labels
                    (e.g. `ad.read_h5ad(keeper).obs[['cell_type']].to_csv(...)`).
Our annotated labels of the form "x - identity" are reduced to "identity" for display.

Run anywhere with anndata + pandas + scikit-learn + matplotlib (Mac venv or APPTAINER_RSC).
Memory: external obs loads fully (backed) — fine for cohort scale; ours is a light CSV.

Usage:
    uv run python pipeline/python/compare_external_typing.py \\
        --external-h5ad wenyu.h5ad --external-key cell_type --label-external Wenyu \\
        --ours-csv stage4_qc/keeper_celltypes.csv --label-ours keeper \\
        --output-dir stage4_qc/figures
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import adjusted_rand_score

OUR_ID_RE = r"(.+)_F(\d+)_C(\d+)$"        # <slide>_F<fov>_C<cell>
EXTERNAL_ID_RE = r"(.+)_(\d+)_(\d+)$"      # <slide>_<fov>_<cell>


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--external-h5ad", type=Path, required=True)
    p.add_argument("--external-key", default="cell_type")
    p.add_argument("--ours-csv", type=Path, required=True,
                   help="Per-cell labels CSV: index = our cell id, one label column.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--label-external", default="external")
    p.add_argument("--label-ours", default="ours")
    p.add_argument("--our-id-regex", default=OUR_ID_RE,
                   help=r"Regex with 3 groups (slide, fov, cell) for OUR ids.")
    p.add_argument("--external-id-regex", default=EXTERNAL_ID_RE,
                   help=r"Regex with 3 groups (slide, fov, cell) for EXTERNAL ids.")
    return p.parse_args()


def keyer(regex: str):
    rx = re.compile(regex)
    def key(idx: str):
        m = rx.match(str(idx))
        return f"{m.group(1)}|{m.group(2)}|{m.group(3)}" if m else None
    return key


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ours = pd.read_csv(args.ours_csv, index_col=0)
    ours = ours.iloc[:, [0]]; ours.columns = ["ours"]
    ours["ours"] = ours["ours"].astype(str).map(lambda s: s.split(" - ", 1)[1] if " - " in s else s)
    ours["key"] = [keyer(args.our_id_regex)(i) for i in ours.index]
    nbad = ours["key"].isna().sum()
    print(f"ours: {len(ours):,} cells ({nbad} ids unparsed)")

    a = ad.read_h5ad(args.external_h5ad, backed="r")
    if args.external_key not in a.obs:
        raise SystemExit(f"external obs missing '{args.external_key}'")
    ext = pd.DataFrame({"ext": a.obs[args.external_key].astype(str).values},
                       index=a.obs.index.astype(str))
    ext["key"] = [keyer(args.external_id_regex)(i) for i in ext.index]
    print(f"{args.label_external}: {len(ext):,} cells ({ext['key'].isna().sum()} ids unparsed)")

    m = ours.dropna(subset=["key"]).merge(ext.dropna(subset=["key"]), on="key", how="inner")
    print(f"shared cells: {len(m):,}")
    if not len(m):
        raise SystemExit("no shared cells — check the id regexes")
    print(f"adjusted Rand index: {adjusted_rand_score(m['ours'], m['ext']):.4f}")

    # per external-type -> top of-ours mapping
    rows = []
    for et, sub in m.groupby("ext", observed=True):
        top = sub["ours"].value_counts(normalize=True)
        rows.append({"external_type": et, "n_cells": len(sub),
                     "top_ours": top.index[0], "top_frac": round(float(top.iloc[0]), 3)})
    agree = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    apath = args.output_dir / "external_vs_ours_agreement.csv"
    agree.to_csv(apath, index=False)
    print(f"wrote {apath}\n", agree.to_string(index=False))

    ct = pd.crosstab(m["ext"], m["ours"])
    ct.to_csv(args.output_dir / "external_vs_ours_crosstab.csv")
    rn = ct.div(ct.sum(1), axis=0)
    rorder = [r for r in m["ext"].value_counts().index if r in rn.index]
    corder = [c for c in m["ours"].value_counts().index if c in rn.columns]
    rn = rn.loc[rorder, corder]

    fig, ax = plt.subplots(figsize=(max(10, 0.4 * rn.shape[1] + 4),
                                    max(9, 0.32 * rn.shape[0] + 3)))
    im = ax.imshow(rn.values, aspect="auto", cmap="magma_r", vmin=0, vmax=1)
    ax.set_xticks(range(rn.shape[1])); ax.set_xticklabels(rn.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(rn.shape[0])); ax.set_yticklabels(rn.index, fontsize=7)
    ax.set_xlabel(args.label_ours); ax.set_ylabel(args.label_external)
    ax.set_title(f"{args.label_external} cell type → {args.label_ours} "
                 f"(row-normalized, {len(m):,} shared cells)", fontsize=12)
    fig.colorbar(im, ax=ax, fraction=0.025, label="fraction of external type")
    fig.tight_layout()
    hpath = args.output_dir / "external_vs_ours_crosstab.png"
    fig.savefig(hpath, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {hpath}")


if __name__ == "__main__":
    main()
