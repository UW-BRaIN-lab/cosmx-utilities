#!/usr/bin/env python3
"""Export the top-N marker genes per cluster from an InSituType result to a tidy CSV.

The marker heatmap caps at a handful of genes per cluster for legibility; this writes the
full top-N (default 20) to a shareable table. Two ranking columns per cluster:
  - specificity (fold-enrichment) = profile / panel-wide mean  -> what's DISTINCTIVE (the
    identity markers), the same metric the heatmap and inspect_denovo_profiles.py use.
  - expression = raw profile value                             -> what the cluster is MADE of.

Reads an InSituType result h5 (/profiles, /profile_genes, /profile_types, optional /cell_type).
Optionally joins a denovo_annotations CSV (denovo_label, annotation, new_name, keep) so each
cluster carries its human label. De-novo clusters (1-2 lowercase letters) by default; add
--include-named for every cell type.

Output (long/tidy): one row per (cluster, rank) with cluster, annotation, new_name, keep, size,
rank, gene, fold_enrichment, mean_expression. Long format sorts/filters cleanly in Excel.

Usage:
    uv run --no-project --with h5py --with numpy --with pandas python \\
        pipeline/python/export_top_markers.py \\
            --profiles-h5 anchor_typing.h5 \\
            --annotations pipeline/reference/denovo_annotations/fullcohort_pruned_k27.csv \\
            --top-n 20 --output top20_markers_per_cluster.csv [--include-named]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DENOVO_RE = re.compile(r"^[a-z]{1,2}$")  # InSituType de-novo labels = 1-2 lowercase letters


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True)
    p.add_argument("--annotations", type=Path, default=None,
                   help="denovo_annotations CSV to join labels (denovo_label + annotation/...).")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--include-named", action="store_true",
                   help="Also export named cell types, not just the de-novo clusters.")
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


def _decode(a) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in a]


def read_profiles(h5: Path) -> tuple[pd.DataFrame, "pd.Series | None"]:
    with h5py.File(h5, "r") as f:
        for k in ("profiles", "profile_genes", "profile_types"):
            if k not in f:
                sys.exit(f"ERROR: /{k} missing in {h5}")
        mat = f["profiles"][()]
        genes = _decode(f["profile_genes"][()])
        types = _decode(f["profile_types"][()])
        cell_type = _decode(f["cell_type"][()]) if "cell_type" in f else None
    ng, nt = len(genes), len(types)
    if mat.shape == (nt, ng):
        mat = mat.T
    elif mat.shape != (ng, nt):
        sys.exit(f"ERROR: /profiles shape {mat.shape} matches neither ({ng},{nt}) nor its transpose.")
    df = pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types).astype(float)
    sizes = pd.Series(cell_type).value_counts() if cell_type is not None else None
    return df, sizes


def main() -> None:
    args = parse_args()
    profiles, sizes = read_profiles(args.profiles_h5)
    types = list(profiles.columns)
    cols = types if args.include_named else [t for t in types if DENOVO_RE.match(t)]
    if not cols:
        sys.exit("No clusters to export (no de-novo labels found; try --include-named).")

    # annotation join (optional)
    ann = {}
    if args.annotations is not None:
        a = pd.read_csv(args.annotations)
        a["denovo_label"] = a["denovo_label"].astype(str).str.strip()
        a = a.set_index("denovo_label")
        for col in ("annotation", "new_name", "keep"):
            if col in a.columns:
                ann[col] = a[col].to_dict()

    panel_mean = profiles.mean(axis=1).replace(0, np.nan)
    rows = []
    for c in cols:
        fold = (profiles[c] / panel_mean).replace([np.inf, -np.inf], np.nan).dropna()
        top = fold.sort_values(ascending=False).head(args.top_n)
        for rank, (gene, fe) in enumerate(top.items(), start=1):
            rows.append({
                "cluster": c,
                "annotation": ann.get("annotation", {}).get(c, ""),
                "new_name": ann.get("new_name", {}).get(c, ""),
                "keep": ann.get("keep", {}).get(c, ""),
                "size": int(sizes.get(c, 0)) if sizes is not None else np.nan,
                "rank": rank,
                "gene": gene,
                "fold_enrichment": round(float(fe), 2),
                "mean_expression": round(float(profiles.loc[gene, c]), 3),
            })
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote {args.output}: {out['cluster'].nunique()} clusters x top-{args.top_n} "
          f"= {len(out)} rows.")


if __name__ == "__main__":
    main()
