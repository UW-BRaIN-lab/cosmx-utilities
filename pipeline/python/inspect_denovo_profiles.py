#!/usr/bin/env python3
"""Inspect the de-novo clusters in an InSituType result, to triage them for the InSituTree reference.

Phase B of the reference rebuild: a semi-supervised InSituType rescale (72_anchor_typing.sh) produces
$profiles = the named GBmap types re-scaled to this cohort + de-novo clusters (single lowercase
letters). We keep the de-novo clusters that are REAL tumor programs (distinct from any conserved
GBmap type) and drop those that merely duplicate a named type. This script surfaces, per de-novo
cluster, the evidence for that decision:

  - top SPECIFIC markers  = genes most enriched in this cluster vs the panel-wide mean (identity)
  - top EXPRESSED markers = highest raw profile value (what it's made of)
  - nearest NAMED type    = max profile cosine to a named GBmap type + that cosine (a high cosine to
                            a conserved type => likely a duplicate to DROP; low => a distinct program)
  - size                  = cells assigned (from /cell_type, if present)

Reads an InSituType result h5 (R/insitutype_typing.R): /profiles (genes x types), /profile_genes,
/profile_types, and optionally /cell_type for sizes. Emits a CSV + a readable printout to annotate.

Usage:
    uv run --no-project --with h5py --with numpy --with pandas python \\
        pipeline/python/inspect_denovo_profiles.py \\
            --profiles-h5 anchor_typing.h5 --output denovo_inspection.csv --top-k 12
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DENOVO_RE = re.compile(r"^[a-z]$")   # de-novo clusters are single lowercase letters


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True,
                   help="InSituType result h5 with /profiles, /profile_genes, /profile_types.")
    p.add_argument("--output", type=Path, required=True, help="Output inspection CSV.")
    p.add_argument("--top-k", type=int, default=12, help="Markers to list per de-novo cluster.")
    return p.parse_args()


def _decode(arr) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]


def read_profiles(h5: Path) -> tuple[pd.DataFrame, np.ndarray | None]:
    """genes x types DataFrame (orientation-robust), + optional per-cell cell_type for sizes."""
    with h5py.File(h5, "r") as f:
        for k in ("profiles", "profile_genes", "profile_types"):
            if k not in f:
                print(f"ERROR: /{k} missing in {h5}", file=sys.stderr); sys.exit(1)
        mat = f["profiles"][()]
        genes = _decode(f["profile_genes"][()])
        types = _decode(f["profile_types"][()])
        cell_type = _decode(f["cell_type"][()]) if "cell_type" in f else None
    ng, nt = len(genes), len(types)
    if mat.shape == (nt, ng):
        mat = mat.T
    elif mat.shape != (ng, nt):
        print(f"ERROR: /profiles shape {mat.shape} matches neither ({ng},{nt}) nor transpose.",
              file=sys.stderr); sys.exit(1)
    df = pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types).astype(float)
    return df, (np.array(cell_type) if cell_type is not None else None)


def main() -> None:
    args = parse_args()
    profiles, cell_type = read_profiles(args.profiles_h5)
    types = list(profiles.columns)
    denovo = [t for t in types if DENOVO_RE.match(t)]
    named = [t for t in types if not DENOVO_RE.match(t)]
    print(f"{profiles.shape[0]:,} genes x {len(types)} types: {len(named)} named + "
          f"{len(denovo)} de-novo {denovo}")
    if not denovo:
        print("No de-novo clusters found.", file=sys.stderr); sys.exit(1)

    sizes = (pd.Series(cell_type).value_counts() if cell_type is not None else None)

    # Specificity = profile / panel-wide mean expression per gene (fold-enrichment). Cosine to
    # named types uses L2-normalised profile vectors.
    panel_mean = profiles.mean(axis=1).replace(0, np.nan)
    named_mat = profiles[named].to_numpy()
    named_norm = named_mat / (np.linalg.norm(named_mat, axis=0, keepdims=True) + 1e-12)

    rows = []
    for d in denovo:
        col = profiles[d]
        fold = (col / panel_mean).replace([np.inf, -np.inf], np.nan).dropna()
        top_spec = fold.sort_values(ascending=False).head(args.top_k).index.tolist()
        top_expr = col.sort_values(ascending=False).head(args.top_k).index.tolist()
        v = col.to_numpy(); v = v / (np.linalg.norm(v) + 1e-12)
        cos = named_norm.T @ v
        j = int(np.argmax(cos))
        rows.append({
            "denovo_label": d,
            "size": int(sizes.get(d, 0)) if sizes is not None else np.nan,
            "nearest_named": named[j],
            "cosine_to_nearest": round(float(cos[j]), 3),
            "top_specific_markers": ",".join(top_spec),
            "top_expressed_markers": ",".join(top_expr),
        })
    out = pd.DataFrame(rows)
    if sizes is not None:
        out = out.sort_values("size", ascending=False)
    out.to_csv(args.output, index=False)

    print(f"\n{'lbl':<4}{'size':>9}  {'nearest named (cos)':<34}{'top specific markers'}")
    print("-" * 100)
    for _, r in out.iterrows():
        sz = f"{r['size']:,}" if pd.notna(r["size"]) else "  n/a"
        near = f"{r['nearest_named']} ({r['cosine_to_nearest']})"
        print(f"{r['denovo_label']:<4}{sz:>9}  {near:<34}{r['top_specific_markers']}")
    print(f"\nWrote {args.output}")
    print("HINT: high cosine_to_nearest (>~0.9) => likely a duplicate of that named type (DROP); "
          "low cosine + coherent tumor markers => a distinct program to KEEP + name.")


if __name__ == "__main__":
    main()
