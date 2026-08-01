#!/usr/bin/env python3
"""Stage 5b: InSituCNV copy-number inference for ONE tissue section.

Runs the Moldia InSituCNV recipe (spatial neighbor smoothing + infercnvpy) on a single
per-section AnnData produced by prep_insitucnv_input.py. Processing one tissue section at
a time is the whole point: the spatial neighbor graph must never bridge two donors (see
prep_insitucnv_input.py), so each section is smoothed and baselined in isolation.

Recipe (Moldia manuscript, Colorectal_cancer_CosMxWTx/01_insituCNV):
  1. raw counts -> normalize_total (NOT logged); keep raw in layer 'counts'.
  2. spatial KNN graph (scvelo) -> icv.tl.smooth_data_for_cnv: M = connectivities @ Xnorm
     (layer 'M'), pooling counts across spatial neighbors to recover CNV-usable depth.
  3. re-normalize + log1p the smoothed layer -> layer 'M_log1p'.
  4. infercnvpy.tl.infercnv against the diploid reference cell types -> obsm['X_cnv'].
  5. per-cell cnv_score = mean(X_cnv**2) (an L2 CNV burden) -> obs['cnv_score'].

Gene positions (var['chromosome','start','end'], autosomes) are already annotated by prep,
so every section shares identical infercnv windows and the per-section X_cnv matrices
concatenate cleanly in the compare step.

Reads:
  --input          one sections/<id>.h5ad (X raw gene counts, obsm['spatial'], var positions,
                   obs['cell_type']).
  --reference-file diploid reference cell_type list (superset; types absent from THIS section
                   are skipped with a warning).
Writes:
  --output         <id>_cnv.h5ad with obsm['X_cnv'], uns['cnv'], obs[..., 'cnv_score']; the
                   heavy smoothed layers are dropped to keep the file small.

Usage:
    python pipeline/python/run_insitucnv.py \\
        --input sections/SLIDE__DONOR.h5ad \\
        --reference-file pipeline/reference/insitucnv_reference_types.txt \\
        --output SLIDE__DONOR_cnv.h5ad \\
        --n-neighbors 20 --window-size 100 --step 10 --dynamic-threshold 1.5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

# Moldia package (on PYTHONPATH via the SIF). Its __init__ is lazy — it declares
# __all__=["pp","tl"] but does NOT import the submodules, so `import insitucnv; insitucnv.tl`
# raises AttributeError. Import the smoothing function directly (as the container build-check
# does) instead of going through the package namespace.
from insitucnv.tl.moments import smooth_data_for_cnv


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--reference-file", type=Path, required=True)
    p.add_argument("--reference-vector", type=Path, default=None,
                   help="CSV of diploid references in the M_log1p space (gene x [donors + "
                        "GLOBAL], from build_insitucnv_reference.py). If given, inferCNV "
                        "baselines against this section's DONOR column (matched normal) when "
                        "present, else the GLOBAL pool — instead of the section's own cells, "
                        "avoiding reference contamination in tumor-bulk sections.")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--donor-key", default="Case",
                   help="obs column identifying the donor, to pick the matched reference column.")
    p.add_argument("--spatial-key", default="spatial")
    p.add_argument("--n-neighbors", type=int, default=20,
                   help="Spatial neighbors for the smoothing graph (MUST match the reference "
                        "build in build_insitucnv_reference.py / 86b).")
    p.add_argument("--window-size", type=int, default=100, help="infercnv window size.")
    p.add_argument("--step", type=int, default=10, help="infercnv step.")
    p.add_argument("--dynamic-threshold", type=float, default=1.5,
                   help="infercnv denoising threshold (SDs).")
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--chunksize", type=int, default=5000, help="infercnv chunk size.")
    p.add_argument("--n-jobs", type=int, default=None, help="infercnv threads.")
    return p.parse_args()


def read_reference_types(path: Path) -> list[str]:
    out = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def main() -> None:
    args = parse_args()

    print(f"Reading {args.input}")
    adata = ad.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in ("chromosome", "start", "end"):
        if col not in adata.var:
            sys.exit(f"ERROR: var['{col}'] missing; run prep_insitucnv_input.py first.")
    if args.spatial_key not in adata.obsm:
        sys.exit(f"ERROR: obsm['{args.spatial_key}'] missing; run prep first.")
    if args.celltype_key not in adata.obs:
        sys.exit(f"ERROR: obs['{args.celltype_key}'] missing.")

    # reference cell types present IN THIS SECTION
    wanted = read_reference_types(args.reference_file)
    adata.obs[args.celltype_key] = adata.obs[args.celltype_key].astype("category")
    present = set(adata.obs[args.celltype_key].cat.categories)
    ref = [r for r in wanted if r in present]
    missing = [r for r in wanted if r not in present]
    if missing:
        print(f"  reference types absent from this section (skipped): {len(missing)}")
    n_ref = int(adata.obs[args.celltype_key].isin(ref).sum())
    if args.reference_vector is None:
        if not ref:
            sys.exit("ERROR: none of the diploid reference cell types are present in this "
                     "section — cannot baseline. (prep should have skipped it, or pass "
                     "--reference-vector for a global baseline.)")
        print(f"  diploid reference (per-section): {len(ref)} types, {n_ref:,} cells")
    else:
        print(f"  diploid reference: GLOBAL vector '{args.reference_vector.name}' "
              f"(section also has {n_ref:,} local reference cells)")

    # 1. raw -> normalized (not logged); stash raw
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=args.target_sum)

    # 2. spatial neighbor graph + neighborhood smoothing (M = connectivities @ Xnorm)
    print(f"Spatial neighbor graph (n_neighbors={args.n_neighbors}) + smoothing ...")
    scv.pp.neighbors(adata, use_rep=args.spatial_key, n_neighbors=args.n_neighbors)
    smooth_data_for_cnv(adata, n_neighbors=args.n_neighbors)  # -> layers['M']
    if "M" not in adata.layers:
        sys.exit("ERROR: smooth_data_for_cnv did not create layer 'M'.")

    # 3. re-normalize + log1p the smoothed layer
    adata.layers["M_log1p"] = adata.layers["M"].copy()
    sc.pp.normalize_total(adata, target_sum=args.target_sum, layer="M_log1p")
    sc.pp.log1p(adata, layer="M_log1p")

    # 4. inferCNV against the diploid reference
    print(f"inferCNV (window={args.window_size}, step={args.step}, "
          f"dynamic_threshold={args.dynamic_threshold}) ...")
    infercnv_kwargs = dict(
        window_size=args.window_size,
        step=args.step,
        dynamic_threshold=args.dynamic_threshold,
        layer="M_log1p",
        chunksize=args.chunksize,
        n_jobs=args.n_jobs,
    )
    if args.reference_vector is not None:
        # Diploid reference (M_log1p space): this section's DONOR column (matched normal) if
        # present, else the GLOBAL pool. Aligned to this section's gene order.
        rv = pd.read_csv(args.reference_vector, index_col=0)
        rv.columns = rv.columns.astype(str)
        donor = None
        if args.donor_key in adata.obs:
            vc = adata.obs[args.donor_key].astype(str).value_counts()
            donor = str(vc.index[0]) if len(vc) else None
        if donor is not None and donor in rv.columns:
            col = donor
            print(f"  reference: matched donor column '{donor}'")
        elif "GLOBAL" in rv.columns:
            col = "GLOBAL"
            print(f"  reference: GLOBAL pool (no matched column for donor '{donor}')")
        else:
            col = rv.columns[0]
            print(f"  reference: column '{col}' (no donor/GLOBAL column found)")
        ref_aligned = rv[col].reindex(adata.var_names).to_numpy(dtype=float)
        n_missing = int(np.isnan(ref_aligned).sum())
        if n_missing:
            print(f"  WARN: {n_missing} genes absent from reference vector; filled with mean")
            ref_aligned = np.nan_to_num(ref_aligned, nan=float(np.nanmean(ref_aligned)))
        cnv.tl.infercnv(adata, reference=ref_aligned.reshape(1, -1), **infercnv_kwargs)
    else:
        cnv.tl.infercnv(adata, reference_key=args.celltype_key, reference_cat=ref,
                        **infercnv_kwargs)
    if "X_cnv" not in adata.obsm:
        sys.exit("ERROR: infercnv did not populate obsm['X_cnv'].")

    # 5. per-cell CNV burden (L2 over genomic windows)
    xcnv = adata.obsm["X_cnv"]
    sq = xcnv.multiply(xcnv) if hasattr(xcnv, "multiply") else np.square(xcnv)
    adata.obs["cnv_score"] = np.asarray(sq.mean(axis=1)).ravel().astype(np.float32)

    # drop the heavy intermediate layers before writing (keep raw counts + X_cnv)
    for lyr in ("M", "M_log1p"):
        adata.layers.pop(lyr, None)
    adata.X = adata.layers.pop("counts")  # restore X to raw counts

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}  "
          f"(X_cnv: {adata.obsm['X_cnv'].shape}, cnv_score median "
          f"{np.median(adata.obs['cnv_score']):.4g})")
    adata.write_h5ad(args.output, compression="gzip")
    print("Done: InSituCNV for one tissue section.")


if __name__ == "__main__":
    main()
