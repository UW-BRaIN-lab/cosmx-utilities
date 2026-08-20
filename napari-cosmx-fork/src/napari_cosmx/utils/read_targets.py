#!/usr/bin/env python

import numpy as np
import pandas as pd
import anndata as ad
import argparse
import os
import sys
import numpy as np
import pandas as pd
import vaex
import glob

def main(args_list=None):
    parser = argparse.ArgumentParser(description='Read decoded targets and write to hdf5')
    parser.add_argument("folder",
        help="Voting folder")
    parser.add_argument("-o", "--outputdir",
        help="Where to write hdf5 file",
        default=".")
    parser.add_argument("--filename",
        help="Name for hdf5 file",
        default="targets.hdf5")
    args = parser.parse_args(args=args_list)

    def read_targets_file(filename):
        print(f"Reading targets from {filename}...")
        header = pd.read_csv(filename, nrows=0)
        usecols = ["fov", "CellId", "x", "y", "target", "CellComp"]
        dtype = {
            "fov": int,
            "CellId": int,
            "target": "category",
            "x": float,
            "y": float,
            "CellComp": "category"
        }
        if "z" in header.columns:
            usecols.insert(4, "z")
            dtype["z"] = int
        return pd.read_csv(
            filename,
            usecols=usecols,
            dtype=dtype)

    def find_target_call_files(path, fov_folder="", summary_folder=""):
        return glob.glob(os.path.join(
            path,
            fov_folder,
            summary_folder,
            "[a-zA-Z0-9_]*__complete_code_cell_target_call_coord.csv"
        ))

    res = find_target_call_files(args.folder, fov_folder="FOV[0-9]*", summary_folder="FOV_Analysis_Summary")
    if len(res) == 0:
        res = find_target_call_files(args.folder, fov_folder="FOV[0-9]*")
    if len(res) == 0:
        res = find_target_call_files(args.folder)
    if len(res) == 0:
        sys.exit(f"No target call files found at {args.folder}")

    targets = pd.concat([read_targets_file(i) for i in res])
    if 'z' not in targets:
        targets['z'] = 0

    df = vaex.from_pandas(targets)
    # using v5 preview features for categorical encoding
    df = df._future()
    df.ordinal_encode("target", inplace=True)
    df.ordinal_encode("CellComp", inplace=True)
    output_path = os.path.join(args.outputdir, args.filename)
    print(f"Writing targets to {output_path}...")
    df.export_hdf5(
        output_path,
        mode='w'
    )

if __name__ == '__main__':
    sys.exit(main())