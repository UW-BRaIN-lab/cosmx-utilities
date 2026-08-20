#!/usr/bin/env python

from napari_cosmx import DEFAULT_NDIM, DEFAULT_Z_STEP_UM
import numpy as np
import pandas as pd
import anndata as ad
import argparse
import os
import sys
import numpy as np
import pandas as pd
import os.path

def main(args_list=None):
    parser = argparse.ArgumentParser(description='Create AnnData object from counts matrix and metadata')
    parser.add_argument("-X", "--counts",
        help="cell x gene counts (raw) matrix in MatrixMarket format")
    parser.add_argument("--obs",
        help="cell metadata file in csv format, with first column as index")
    parser.add_argument("--var",
        help="feature/gene metadata in csv format, with first column as index")
    parser.add_argument("--coords",
        help="Spatial coords CSV. The first columns must be in axis order: y,x for --ndim 2 and z,y,x for --ndim 3.")
    parser.add_argument("--umap",
        help="umap dims file in csv format")
    parser.add_argument("-o", "--outputdir",
        help="Where to write h5ad file",
        default=".")
    parser.add_argument("--filename",
        help="Name for h5ad file",
        default="adata.h5ad")
    parser.add_argument("-n", "--name",
        help="Name of anndata object",
        default="CosMx study")
    parser.add_argument("--ndim",
        help="Optional: Dimensionality of spatial coordinates. Axis order is y,x for 2D and z,y,x for 3D.",
        choices=[2, 3],
        default=DEFAULT_NDIM,
        type=int)
    parser.add_argument("--z-step-um",
        help="Optional: Spacing between z planes in microns.",
        default=DEFAULT_Z_STEP_UM,
        type=float)
    parser.add_argument("--colors",
        help="csv files with colors to import",
        nargs='*')
    args = parser.parse_args(args=args_list)

    if not any([args.counts, args.obs]):
        sys.exit("Need counts or obs to create AnnData object")
    X = obs = var = None
    if args.counts is not None:
        X = ad.read_mtx(args.counts, dtype=np.int32).X
    if args.obs is not None:
        obs = pd.read_csv(args.obs, index_col=0)
    if args.var is not None:
        var = pd.read_csv(args.var, index_col=0)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    if args.coords is not None:
        coords = pd.read_csv(args.coords)
        if coords.shape[1] < args.ndim:
            axis_order = "y,x" if args.ndim == 2 else "z,y,x"
            sys.exit(f"Expected at least {args.ndim} coordinate columns in {args.coords} for axis order {axis_order}, found {coords.shape[1]}")
        adata.obsm['spatial'] = coords.iloc[:, :args.ndim].to_numpy()
    if args.umap is not None:
        adata.obsm['umap'] = pd.read_csv(args.umap).to_numpy()
    adata.uns['name'] = args.name
    adata.uns['CosMx'] = {
        'ndim': args.ndim,
        'z_step_um': args.z_step_um,
    }
    adata.strings_to_categoricals()
    if args.colors is not None:
        for i in args.colors:
            file = os.path.basename(i)
            cat = file.rpartition("_colors.csv")[0]
            if cat in adata.obs:
                if adata.obs[cat].dtype.name == 'category':
                    cols = pd.read_csv(i, header=None)
                    colors_dict = dict(zip(cols[0], cols[1]))
                    adata.uns[cat + "_colors"] = [colors_dict[k] for k in adata.obs[cat].cat.categories]
                else:
                    print(f"{cat} is not a categorical column in AnnData object")
            else:
                print(f"{cat} not found in AnnData obs")

    adata.write(os.path.join(args.outputdir, args.filename), compression="gzip")

if __name__ == '__main__':
    sys.exit(main())