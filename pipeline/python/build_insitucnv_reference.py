#!/usr/bin/env python3
"""Stage 5b-ref: build a clean diploid reference for InSituCNV — per donor, plus a global pool.

Tumor-bulk tissue sections have scarce, tumor-surrounded normal cells, so baselining
inferCNV against a section's OWN reference cells contaminates the baseline (the normal
cells' spatially-smoothed profiles pick up tumor RNA) and suppresses the CNV signal. This
builds a diploid reference from reference-type cells in CONTRALATERAL-UNINVOLVED tissue —
genuinely diploid, in clean (non-tumor) neighborhoods — processed through the SAME
normalize -> spatial-smooth -> log recipe as the test cells, then averaged per gene.

PER-DONOR MATCHED reference is the CNV gold standard: comparing a donor's tumor to THAT
donor's own contralateral normal cancels the donor's germline/expression baseline and
batch (a cross-donor pool leaves those in, and they can masquerade as CNV). So we emit one
reference column per donor that has enough contralateral reference cells, PLUS a GLOBAL
column pooling all donors — the fallback for donors with no (or too little) contralateral
tissue. run_insitucnv.py picks its section's donor column when present, else GLOBAL.

Reads (--sections-dir): per-section h5ads from prep_insitucnv_input.py (obs cell_type,
Region, Case; obsm spatial; X raw gene counts; var chr positions).
Writes (--output): CSV indexed by gene, one column per qualifying donor + a 'GLOBAL' column
(each an M_log1p reference profile).

Usage:
    python pipeline/python/build_insitucnv_reference.py \\
        --sections-dir sections --reference-file pipeline/reference/insitucnv_reference_types.txt \\
        --output reference_vector.csv --n-neighbors 20
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

from insitucnv.tl.moments import smooth_data_for_cnv


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sections-dir", type=Path, required=True,
                   help="Directory of per-section h5ads from prep_insitucnv_input.py.")
    p.add_argument("--reference-file", type=Path, required=True,
                   help="Diploid reference cell_type list (one per line; '#' comments ok).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output reference CSV (gene x [donors + GLOBAL]).")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--region-key", default="Region")
    p.add_argument("--donor-key", default="Case")
    p.add_argument("--contralateral-label", default="Contralateral uninvolved",
                   help="obs['Region'] value marking clean, uninvolved tissue.")
    p.add_argument("--n-neighbors", type=int, default=20,
                   help="Spatial smoothing neighbors — MUST match the run_insitucnv setting.")
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--min-ref-cells", type=int, default=20,
                   help="Skip a section contributing fewer clean reference cells than this.")
    p.add_argument("--min-donor-cells", type=int, default=50,
                   help="Emit a per-donor column only if the donor has at least this many "
                        "contralateral reference cells; otherwise it falls back to GLOBAL.")
    return p.parse_args()


def read_reference_types(path: Path) -> set[str]:
    return {ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()
            if ln.split("#", 1)[0].strip()}


def main() -> None:
    args = parse_args()
    ref_types = read_reference_types(args.reference_file)
    files = sorted(args.sections_dir.glob("*.h5ad"))
    if not files:
        sys.exit(f"ERROR: no *.h5ad in {args.sections_dir}")
    print(f"{len(files)} sections; diploid reference types: {len(ref_types)}")

    genes = None
    donor_sum: dict[str, np.ndarray] = {}
    donor_n: dict[str, int] = {}
    global_sum, global_n = None, 0

    for f in files:
        adata = ad.read_h5ad(f)
        if not {args.region_key, args.celltype_key, args.donor_key} <= set(adata.obs):
            continue
        mask = (adata.obs[args.celltype_key].astype(str).isin(ref_types).to_numpy()
                & (adata.obs[args.region_key].astype(str).to_numpy() == args.contralateral_label))
        n = int(mask.sum())
        if n < args.min_ref_cells:
            continue

        # same recipe as run_insitucnv: normalize -> spatial-smooth -> renorm + log1p
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=args.target_sum)
        scv.pp.neighbors(adata, use_rep="spatial", n_neighbors=args.n_neighbors)
        smooth_data_for_cnv(adata, n_neighbors=args.n_neighbors)
        adata.layers["M_log1p"] = adata.layers["M"].copy()
        sc.pp.normalize_total(adata, target_sum=args.target_sum, layer="M_log1p")
        sc.pp.log1p(adata, layer="M_log1p")

        if genes is None:
            genes = adata.var_names.to_numpy()
        elif not np.array_equal(genes, adata.var_names.to_numpy()):
            sys.exit(f"ERROR: gene set of {f.name} differs from earlier sections.")

        M = adata.layers["M_log1p"][mask]            # (n_ref, n_genes)
        donors = adata.obs[args.donor_key].astype(str).to_numpy()[mask]
        section_sum = np.asarray(M.sum(axis=0)).ravel()
        global_sum = section_sum if global_sum is None else global_sum + section_sum
        global_n += n
        for d in np.unique(donors):
            rows = donors == d
            s = np.asarray(M[rows].sum(axis=0)).ravel()
            donor_sum[d] = s if d not in donor_sum else donor_sum[d] + s
            donor_n[d] = int(rows.sum()) + donor_n.get(d, 0)
        print(f"  {f.name}: +{n:,} contralateral reference cells "
              f"(donors {sorted(np.unique(donors))})")

    if global_sum is None or global_n == 0:
        sys.exit("ERROR: no contralateral reference cells found across sections — cannot "
                 "build a reference. (Is the --contralateral-label correct?)")

    cols = {"GLOBAL": global_sum / global_n}
    kept, dropped = [], []
    for d, n in sorted(donor_n.items()):
        if n >= args.min_donor_cells:
            cols[d] = donor_sum[d] / n
            kept.append(f"{d}({n})")
        else:
            dropped.append(f"{d}({n})")
    out = pd.DataFrame(cols, index=pd.Index(genes, name="gene"))
    out.to_csv(args.output)
    print(f"\nReference: GLOBAL from {global_n:,} cells; per-donor columns kept "
          f"({args.min_donor_cells}+ cells): {kept}")
    if dropped:
        print(f"  donors below threshold -> GLOBAL fallback: {dropped}")
    print(f"  -> {args.output}  ({out.shape[1]} columns: {out.shape[1] - 1} donors + GLOBAL)")


if __name__ == "__main__":
    main()
