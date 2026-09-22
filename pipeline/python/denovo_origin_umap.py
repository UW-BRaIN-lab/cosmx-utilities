#!/usr/bin/env python3
"""75m: co-embed a de-novo letter's cells with the natively-called cells of its destinations.

The PI's third leg: compute one embedding over the vascular cells alone and ask whether the
forced cells intermingle with the native cluster they were assigned to, or sit as a satellite
island beside it. A co-embedding is the right instrument because the question is about
NEIGHBOURHOODS in expression space, and a cohort-wide UMAP has no resolution inside a
compartment that is 4% of the cells.

The figure answers it by eye. The number that goes with it is read off the k-nearest-neighbour
graph the embedding already built, in PCA space rather than in the distorted 2-D projection:

    for every cell of group G, what fraction of its graph neighbours belong to group H?

divided by H's share of the co-embedded cells, so 1.0 means "no more often than chance". A
forced group that intermingles has a near-1 entry against its native counterpart and no strong
self-enrichment; a satellite island has a large diagonal and a small off-diagonal, whatever its
markers look like.

Note what this can and cannot settle. Intermingling is evidence the two groups are one
population. Separation is NOT evidence of shoehorning by itself — the assignment put these cells
in different bins BECAUSE they differ, so some separation is built in (75d's header makes the
same point about markers). The informative comparison is against the native-vs-native entries in
the same matrix: two native vascular types are two real, distinct populations, so their mutual
enrichment sets the scale for what "genuinely different" looks like here.

Inputs:
  --typed-h5ad  cosmx_typed.h5ad — read backed, then subset to the compared cells.
  the shared origin flags (--typing-h5, --forced-csv, --letter, --destinations/--crosstab)
Writes (--output-dir):
  umap_group.png / umap_origin.png / umap_donor.png / umap_depth.png
  neighbour_enrichment.csv   group x group, observed/expected on the kNN graph
  (--output-h5ad) the embedded AnnData, obsm['X_umap']

Usage:
    python pipeline/python/denovo_origin_umap.py \\
        --typed-h5ad cosmx_typed.h5ad --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv --letter l \\
        --destinations Endo_capilar,Endo_arterial,Pericyte --output-dir l_umap
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# keep numba off the GPU JIT path that trips Hyak's driver (see cluster_embedding.py)
os.environ.setdefault("NUMBA_CUDA_ENABLE_PYNVJITLINK", "0")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from denovo_origin import add_origin_arguments, load_origin
from umap_lowsignal import embed, scatter

DEPTH_CANDIDATES = ["qc_gene_counts", "nCount_RNA", "total_counts", "nCount"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True)
    add_origin_arguments(p)
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--output-h5ad", type=Path, default=None)
    return p.parse_args()


def neighbour_enrichment(graph: sp.spmatrix, groups: pd.Series,
                         order: list[str]) -> pd.DataFrame:
    """Row G, column H: share of G's graph neighbours in H, over H's share of all cells.

    The graph is binarised first — scanpy's connectivities are UMAP-fuzzed weights, and
    weighting a neighbour count by them would mix "how many" with "how close".
    """
    adj = (graph > 0).astype(np.float32)
    codes = pd.Categorical(groups, categories=order)
    onehot = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32),
         (np.arange(len(codes)), codes.codes)), shape=(len(codes), len(order)))
    counts = pd.DataFrame(np.asarray((adj @ onehot).todense()),
                          index=groups.index, columns=order)
    observed = counts.groupby(groups.to_numpy(), observed=True).sum()
    observed = observed.div(observed.sum(axis=1), axis=0)
    expected = pd.Series(codes.value_counts(), index=order) / len(codes)
    return observed.reindex(order).div(expected, axis=1)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tidy, order = load_origin(args)
    # A cell in the whole-letter reference group is already in its own pair group; the
    # embedding needs each cell once, so the pair label wins and the [all] group is dropped.
    paired = tidy[tidy["group"].isin([g for g in order if not g.endswith("[all]")])]
    per_cell = paired.groupby(level=0)["group"].first()
    embed_order = [g for g in order if g in set(per_cell)]

    print(f"Reading {args.typed_h5ad} (backed) and subsetting to the compared cells")
    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    keep = adata.obs_names.isin(per_cell.index)
    print(f"cell-id overlap: {int(keep.sum()):,} of {len(per_cell):,} compared cells "
          f"({keep.sum() / max(len(per_cell), 1):.1%})")
    if not keep.any():
        sys.exit("ERROR: no compared cell ids are present in the typed run — check the join.")
    adata = adata[keep].to_memory()
    if "probe_type" in adata.var:
        adata = adata[:, (adata.var["probe_type"] == "gene").to_numpy()].copy()
    adata.obs["group"] = per_cell.reindex(adata.obs_names).astype("object")
    adata.obs["origin"] = tidy.groupby(level=0)["origin"].first().reindex(
        adata.obs_names).astype("object")
    print(f"Co-embedding {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    depth_col = next((c for c in DEPTH_CANDIDATES if c in adata.obs), None)
    adata = embed(adata, args)
    xy = np.asarray(adata.obsm["X_umap"])

    scatter(xy, adata.obs["group"], f"{args.letter} and its destinations — group",
            args.output_dir / "umap_group.png", order=embed_order)
    scatter(xy, adata.obs["origin"], f"{args.letter} and its destinations — forced vs native",
            args.output_dir / "umap_origin.png")
    if "Case" in adata.obs:
        scatter(xy, adata.obs["Case"].astype(str), "donor",
                args.output_dir / "umap_donor.png")
    if depth_col is not None:
        scatter(xy, adata.obs[depth_col], f"panel depth ({depth_col})",
                args.output_dir / "umap_depth.png", categorical=False)

    graph = adata.obsp.get("connectivities", adata.obsp.get("distances"))
    if graph is None:
        print("WARNING: no kNN graph in obsp; skipping the enrichment matrix.", file=sys.stderr)
    else:
        enr = neighbour_enrichment(sp.csr_matrix(graph), adata.obs["group"], embed_order)
        enr.to_csv(args.output_dir / "neighbour_enrichment.csv")
        print("\n=== kNN neighbour enrichment (observed / expected) ===")
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(enr.round(2))
        print("\nREAD down each forced row: the entry against its own native column is how far "
              "the two\nintermingle (1.0 = at chance). Calibrate against the native-vs-native "
              "entries -- those are\ntwo genuinely distinct populations measured the same way. "
              "Expected is a share of ALL\nco-embedded cells, so one group sitting off on an "
              "island lifts every other entry above 1;\nthe native-vs-native baseline moves "
              "with them, which is why it is the thing to read against.")

    if args.output_h5ad is not None:
        keep_obs = [c for c in ("group", "origin", "cell_type", "Region", "Case", depth_col)
                    if c and c in adata.obs]
        # Carry the kNN graph, not just the projection. Every figure here is a scatter of
        # X_umap, but the enrichment matrix is computed in PCA space off this graph -- without
        # it a re-render would need the GPU job run again just to recover a number.
        out = ad.AnnData(X=None, obs=adata.obs[keep_obs].copy(), obsm={"X_umap": xy},
                         uns={"denovo_origin_umap": args.letter})
        if graph is not None:
            out.obsp["connectivities"] = sp.csr_matrix(graph)
        out.write_h5ad(args.output_h5ad)
    print(f"\nWrote UMAP figures to {args.output_dir}")


if __name__ == "__main__":
    main()
