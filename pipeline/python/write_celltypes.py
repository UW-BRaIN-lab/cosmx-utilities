#!/usr/bin/env python3
"""Write InSituType cell types back into the clustered cohort AnnData.

Stage 4c of the analysis pipeline. Joins the stage-4b InSituType result onto the
stage-3c clustered AnnData by cell id, producing cosmx_typed.h5ad: the same cells and
UMAP as cosmx_clustered.h5ad, with two new obs columns

  cell_type        InSituType assignment (named GBmap types + de novo letter clusters)
  insitutype_prob  per-cell posterior confidence in that assignment

and leiden retained. Downstream marker/compare tooling then re-renders by cell type with
the existing one-flag flip (marker_pseudobulk.py --group-key cell_type). A UMAP coloured
by cell_type is emitted for review via the shared plot_qc.make_qc_plots.

Optional, for a typed AnnData that can be handed to someone else:

  --label-map NAME=CSV   add obs cell_type_<NAME> by mapping cell_type through a (leaf, collapsed_label)
                         CSV (the same maps celltype_by_group.py takes). Repeatable. cell_type itself is
                         left as the run wrote it.
  --annotations CSV      add the donor columns (case_broad, case_group, sorl1_mutation) from the FOV table,
                         joined on obs[Region] (the donor id)
  --opc-labels A,B       add obs opc_confident: called one of these labels AND (Leiden 0 OR total_counts >= 600)

Usage:
    uv run python pipeline/python/write_celltypes.py \\
        --clustered-h5ad cosmx_clustered.h5ad \\
        --insitutype-h5 insitutype_result.h5 \\
        --output cosmx_typed.h5ad \\
        --qc-plots-dir qc_plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from celltype_by_group import collapse_labels
from export_cell_labels import DONOR_COLUMNS, opc_confident
from plot_qc import make_qc_plots

DONOR_OBS_KEY = "Region"
OPC_RUN_NAME = "typed"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True,
                   help="cosmx_clustered.h5ad from stage 3c (obs has leiden + UMAP).")
    p.add_argument("--insitutype-h5", type=Path, required=True,
                   help="insitutype_result.h5 from stage 4b.")
    p.add_argument("--output", type=Path, required=True,
                   help="Output cosmx_typed.h5ad.")
    p.add_argument("--qc-plots-dir", type=Path, default=None,
                   help="If set, write a UMAP-by-cell_type PNG here.")
    p.add_argument("--label-map", action="append", default=[], metavar="NAME=CSV",
                   help="Add obs cell_type_<NAME> from a (leaf, collapsed_label) CSV. Repeatable.")
    p.add_argument("--annotations", type=Path, default=None,
                   help="FOV table (donor, case_broad, case_group, sorl1_mutation) joined on obs[Region].")
    p.add_argument("--opc-labels", default=None, metavar="LABEL[,LABEL]",
                   help="Add obs opc_confident for these OPC labels.")
    return p.parse_args()


def add_label_maps(obs: pd.DataFrame, label_maps: dict[str, Path]) -> None:
    """cell_type_<name> = cell_type through each map; unmapped labels and NA are kept."""
    for name, path in label_maps.items():
        obs[f"cell_type_{name}"] = collapse_labels(obs["cell_type"].astype(object), path).astype("category")


def add_donor_annotations(obs: pd.DataFrame, annotations: pd.DataFrame) -> None:
    donors = annotations.astype({"donor": str}).drop_duplicates("donor").set_index("donor")
    joined = donors[DONOR_COLUMNS].reindex(obs[DONOR_OBS_KEY].astype(str))
    for column in DONOR_COLUMNS:
        obs[column] = pd.Categorical(joined[column].to_numpy())
    absent = sorted(set(obs[DONOR_OBS_KEY].astype(str)) - set(donors.index))
    if absent:
        print(f"WARN: donors with no annotation row (left NA): {absent}", file=sys.stderr)


def add_opc_confident(obs: pd.DataFrame, labels: list[str]) -> None:
    frame = pd.DataFrame({f"celltype_{OPC_RUN_NAME}": obs["cell_type"].astype(object),
                          "leiden": obs["leiden"], "total_counts": obs["total_counts"]})
    obs["opc_confident"] = opc_confident(frame, OPC_RUN_NAME, labels).to_numpy()


def _decode(arr) -> np.ndarray:
    """hdf5r writes variable-length UTF-8 strings; h5py reads them as bytes or str."""
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def main() -> None:
    args = parse_args()

    print(f"Reading {args.insitutype_h5}")
    with h5py.File(args.insitutype_h5, "r") as f:
        cell_id = _decode(f["cell_id"])
        cell_type = _decode(f["cell_type"])
        prob = np.asarray(f["prob"][()], dtype=np.float64)
    typed = pd.DataFrame({"cell_type": cell_type, "insitutype_prob": prob},
                         index=pd.Index(cell_id, name="cell_id"))
    print(f"  {len(typed):,} typed cells across {typed['cell_type'].nunique()} cell types")

    print(f"Reading {args.clustered_h5ad}")
    adata = ad.read_h5ad(args.clustered_h5ad)

    aligned = typed.reindex(adata.obs.index)
    n_missing = int(aligned["cell_type"].isna().sum())
    if n_missing:
        print(f"WARN: {n_missing:,} / {adata.n_obs:,} clustered cells have no "
              f"InSituType assignment (left as NA)", file=sys.stderr)
    n_extra = int(len(typed) - aligned["cell_type"].notna().sum())
    if n_extra:
        print(f"WARN: {n_extra:,} typed cells not present in the clustered AnnData "
              f"(dropped)", file=sys.stderr)

    adata.obs["cell_type"] = aligned["cell_type"].astype("category")
    adata.obs["insitutype_prob"] = aligned["insitutype_prob"].astype(np.float32)

    label_maps = {}
    for spec in args.label_map:
        name, _, path = spec.partition("=")
        if not name or not path:
            sys.exit(f"--label-map needs NAME=CSV, got {spec!r}")
        label_maps[name] = Path(path)
    add_label_maps(adata.obs, label_maps)
    if args.annotations is not None:
        add_donor_annotations(adata.obs, pd.read_csv(args.annotations, dtype=str))
    if args.opc_labels:
        add_opc_confident(adata.obs, args.opc_labels.split(","))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}")
    adata.write_h5ad(args.output, compression="gzip")

    if args.qc_plots_dir is not None:
        print(f"Rendering UMAP-by-cell_type to {args.qc_plots_dir}")
        make_qc_plots(adata, ["cell_type"], args.qc_plots_dir)

    print(f"Done. {adata.n_obs:,} cells; cell_type + insitutype_prob in obs.")


if __name__ == "__main__":
    main()
