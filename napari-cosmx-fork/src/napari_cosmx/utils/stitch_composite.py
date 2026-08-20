#!/usr/bin/env python

from napari_cosmx import DEFAULT_NDIM, DEFAULT_Z_STEP_UM
from napari_cosmx.utils import _stitch as stitch
from napari_cosmx.utils._patterns import get_fov_number, get_zslice_number
import argparse
import os
import sys
import re
import numpy as np
import pandas as pd
import zarr
import dask.array as da
from skimage import io

COMPOSITE_2D_PATTERN = re.compile(r"CELLCOMPOSITE_F[0-9]+\.JPG")
COMPOSITE_3D_PATTERN = re.compile(r"CELLCOMPOSITE_F[0-9]+_Z[0-9]+\.JPG")

def _collect_composite_tiles(root_dir):
    tile_map = {}
    has_z = False
    for root, dirs, files in os.walk(root_dir):
        for filename in files:
            upper = filename.upper()
            if COMPOSITE_3D_PATTERN.match(upper):
                is_match = True
            elif COMPOSITE_2D_PATTERN.match(upper):
                is_match = True
            else:
                is_match = False
            if not is_match:
                continue
            path = os.path.join(root, filename)
            fov = get_fov_number(path)
            zslice = get_zslice_number(path)
            has_z = has_z or (zslice is not None)
            key = (fov, 0 if zslice is None else zslice)
            tile_map.setdefault(key, []).append(path)
    return tile_map, has_z

def main(args_list=None):
    parser = argparse.ArgumentParser(description='Tile CellComposite images.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inputdir",
        help="Required: Path to CellComposite images.",
        default=".")
    parser.add_argument("-o", "--outputdir",
        help="Required: Path to existing stitched output.",
        default=".")
    parser.add_argument("-u", "--umperpx",
        help="Optional: Override image scale in um per pixel.\n"+
        "Instrument-specific values to use:\n-> beta04 = 0.1228",
        default=None,
        type=float)
    args = parser.parse_args(args=args_list)

    # Check output directory
    if not os.path.exists(os.path.join(args.outputdir, "images")):
        sys.exit(f"{args.outputdir}/images path does not exist.\nRun stitch-images first to create it.")
    store = os.path.join(args.outputdir, "images")
    grp = zarr.open(store, mode = 'a')

    fov_offsets = pd.DataFrame.from_dict(grp.attrs['CosMx']['fov_offsets'])
    output_ndim = int(grp.attrs['CosMx'].get('ndim', DEFAULT_NDIM))
    z_slices = [int(z) for z in grp.attrs['CosMx'].get('z_slices', [0])]
    if 'scale_um' in grp.attrs['CosMx']:
        scale_dict = stitch.get_scales(um_per_px=grp.attrs['CosMx']['scale_um'])
    elif args.umperpx is not None:
        scale_dict = stitch.get_scales(um_per_px=args.umperpx)
    else:
        sys.exit("No um_per_px in metadata or provided as argument")
    scale_dict['z_step_um'] = float(grp.attrs['CosMx'].get('z_step_um', scale_dict.get('z_step_um', DEFAULT_Z_STEP_UM)))

    composite_tiles, composite_have_z = _collect_composite_tiles(args.inputdir)
    composite_res = [path for paths in composite_tiles.values() for path in paths]

    fov_height = grp.attrs['CosMx']['fov_height']
    fov_width = grp.attrs['CosMx']['fov_width']
    dash = (fov_height/fov_width) != 1
    
    top_origin_px, left_origin_px, height, width = stitch.base(
        fov_offsets, fov_height, fov_width, scale_dict, dash)
    
    if len(composite_res) != 0:
        if output_ndim == 3:
            im = da.zeros((len(z_slices), height, width, 3), dtype=np.uint8, chunks=(1,) + stitch.CHUNKS + (3,))
            for z_index, zslice in enumerate(z_slices):
                for fov in fov_offsets['FOV']:
                    tile_path = composite_tiles.get((int(fov), int(zslice)), [])
                    if len(tile_path) == 0:
                        print(f"Could not find CellComposite image for FOV {fov}, z {zslice}")
                        continue
                    if len(tile_path) > 1:
                        print(f"Multiple CellComposite images for FOV {fov}, z {zslice}; using {tile_path[0]}")
                    tile = io.imread(tile_path[0])
                    y, x = stitch.fov_origin(fov_offsets, fov, top_origin_px, left_origin_px, fov_height, scale_dict, dash)
                    im[z_index, y:y+tile.shape[0], x:x+tile.shape[1], :] = tile
                    print(f"Added composite for FOV {fov}, z {zslice}")
        else:
            if composite_have_z:
                target_z = z_slices[0] if len(z_slices) > 0 else 0
            else:
                target_z = 0
            im = da.zeros((height, width, 3), dtype=np.uint8, chunks=stitch.CHUNKS + (3,))
            for fov in fov_offsets['FOV']:
                tile_path = composite_tiles.get((int(fov), int(target_z)), [])
                if len(tile_path) == 0:
                    print(f"Could not find CellComposite image for FOV {fov}")
                    continue
                if len(tile_path) > 1:
                    print(f"Multiple CellComposite images for FOV {fov}; using {tile_path[0]}")
                tile = io.imread(tile_path[0])
                y, x = stitch.fov_origin(fov_offsets, fov, top_origin_px, left_origin_px, fov_height, scale_dict, dash)
                im[y:y+tile.shape[0], x:x+tile.shape[1], :] = tile
                print(f"Added composite for FOV {fov}")
        
        stitch.write_pyramid(im, scale_dict, store=store, path="composite")
    else:
        print(f"No CellComposite images found at {args.inputdir}")

if __name__ == '__main__':
    sys.exit(main())