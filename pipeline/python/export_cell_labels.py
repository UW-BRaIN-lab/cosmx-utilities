#!/usr/bin/env python3
"""One per-cell table: metadata, UMAP and the cell type from every InSituType run.

Built so the typing can be used offline (the per-cell labels, not the 6-10 GB AnnData): one
row per cell with slide, donor, diagnosis, Leiden, UMAP, depth, and a `celltype_<run>` +
`prob_<run>` column pair for each run given. Reads obs + obsm only (backed mode) from the
stage-3c clustered AnnData and the compact insitutype_result.h5 of each run, so no
expression matrix is loaded.

Cells are joined on cell id; a run that typed fewer cells than the AnnData leaves NA and
says so. Writes gzip-compressed CSV.

Usage:
    export_cell_labels.py --clustered-h5ad cosmx_clustered.h5ad \\
        --annotations fov_annotations_sorl1.csv --output cell_labels.csv.gz \\
        --run k21=result_k21.h5 --run k30=result_k30.h5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

OBS_COLUMNS = {
    "slide_id": "slide_id", "Region": "donor", "fov": "fov", "leiden": "leiden",
    "total_counts": "total_counts", "qc_genes_detected": "genes_detected", "qc_area": "area",
}
DONOR_COLUMNS = ["case_broad", "case_group", "sorl1_mutation"]


def _decode(arr) -> np.ndarray:
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def read_run(path: Path) -> pd.DataFrame:
    """cell_type + prob per cell id from a compact InSituType result h5."""
    with h5py.File(path, "r") as f:
        return pd.DataFrame({"cell_type": _decode(f["cell_type"]),
                             "prob": np.asarray(f["prob"][()], dtype=np.float64)},
                            index=pd.Index(_decode(f["cell_id"]), name="cell_id"))


def build_table(obs: pd.DataFrame, umap: np.ndarray, annotations: pd.DataFrame,
                runs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    table = obs[list(OBS_COLUMNS)].rename(columns=OBS_COLUMNS).copy()
    table.index.name = "cell_id"
    table["donor"] = table["donor"].astype(str)
    table["umap1"], table["umap2"] = umap[:, 0], umap[:, 1]

    donors = annotations.astype({"donor": str}).drop_duplicates("donor").set_index("donor")
    table = table.join(donors[DONOR_COLUMNS], on="donor")

    for name, typed in runs.items():
        aligned = typed.reindex(table.index)
        n_missing = int(aligned["cell_type"].isna().sum())
        if n_missing:
            print(f"WARN: run '{name}' has no label for {n_missing:,} of {len(table):,} cells",
                  file=sys.stderr)
        table[f"celltype_{name}"] = aligned["cell_type"]
        table[f"prob_{name}"] = aligned["prob"].round(4)
    return table


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustered-h5ad", type=Path, required=True)
    p.add_argument("--annotations", type=Path, required=True,
                   help="fov_annotations_sorl1.csv (donor -> case_broad, case_group, sorl1_mutation).")
    p.add_argument("--output", type=Path, required=True, help="Output .csv.gz")
    p.add_argument("--run", action="append", default=[], metavar="NAME=RESULT_H5",
                   help="A typing run: name=path/to/insitutype_result.h5. Repeatable.")
    args = p.parse_args()

    import anndata as ad

    runs = {}
    for spec in args.run:
        name, _, path = spec.partition("=")
        if not name or not path:
            p.error(f"--run needs NAME=PATH, got {spec!r}")
        runs[name] = read_run(Path(path))
        print(f"run '{name}': {len(runs[name]):,} typed cells, "
              f"{runs[name]['cell_type'].nunique()} cell types")

    adata = ad.read_h5ad(args.clustered_h5ad, backed="r")
    absent = [c for c in OBS_COLUMNS if c not in adata.obs]
    if absent or "X_umap" not in adata.obsm:
        print(f"ERROR: obs missing {absent}; X_umap present: {'X_umap' in adata.obsm}",
              file=sys.stderr)
        sys.exit(1)
    table = build_table(adata.obs[list(OBS_COLUMNS)].copy(), np.asarray(adata.obsm["X_umap"]),
                        pd.read_csv(args.annotations, dtype=str), runs)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, compression="gzip")
    print(f"Wrote {args.output}: {len(table):,} cells x {table.shape[1]} columns")


if __name__ == "__main__":
    main()
