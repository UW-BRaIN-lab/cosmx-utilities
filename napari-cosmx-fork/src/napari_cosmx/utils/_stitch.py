from napari_cosmx import DASH_UM_PER_PX, ALPHA_UM_PER_PX, BETA_UM_PER_PX, DEFAULT_COLORMAPS, DEFAULT_Z_STEP_UM

import pandas as pd
import dask.array as da
import os
import zarr
from skimage.transform import resize
import tifffile
import json
import math
from numcodecs import Zlib
from functools import partial
from tqdm.auto import tqdm as std_tqdm

zarr.storage.default_compressor = Zlib()
CHUNKS = (8192, 8192)  # 'auto' or tuple

fov_tqdm = partial(
        std_tqdm, desc='Added FOV', unit=" FOVs", ncols=40, mininterval=1.2,
        bar_format="{desc} {n_fmt}/{total_fmt}|{bar}|{percentage:3.0f}%")

def read_image_description_metadata(tiff_path):
    try:
        with tifffile.TiffFile(tiff_path) as im:
            image_description = im.pages[0].tags.get('ImageDescription')
            if image_description is None:
                return None
            return json.loads(image_description.value)
    except Exception:
        return None

def get_z_step_um(tiff_path):
    metadata = read_image_description_metadata(tiff_path)
    if not metadata:
        return None
    z_step_um = metadata.get('ZStackStepSize_um')
    return float(z_step_um) if z_step_um is not None else None

def offsets(offsetsdir: str):
    """Reads FOV coordinates data.

    Args:
        offsetsdir (str): Directory location containing a file
        ending in "fov_positions_file.csv.gz" (AtoMx SIP exported format), 
        or "FOV_Locations.csv" or legacy format "latest.fovs.csv".

    Returns:
        DataFrame: coordinates of each FOV
    """
    legacy_format = True
    atomx_format = False
    for filename in os.listdir(offsetsdir):
        if filename.endswith("FOV_Locations.csv"):
            print(f"Using FOV locations from {filename}")
            legacy_format = False
            atomx_format = False
            df = pd.read_csv(os.path.join(offsetsdir, filename))
            # Slide	X_mm	Y_mm	Z_um    FOV	Order
            if 'Z_mm' not in df.columns and 'Z_um' in df.columns:
                df['Z_um'] = df['Z_um'] / 1e3 # convert to millimeters
                df = df.rename(columns={'Z_um':'Z_mm'})
            z_mm_index = df.columns.get_loc('Z_mm')
            df.insert(z_mm_index + 1, 'ZOffset_mm', -2.0) # hardcoded
            df.insert(z_mm_index + 2, 'ROI', 0) # hardcoded
            break
        if filename.endswith("fov_positions_file.csv.gz"):
            print(f"Using AtoMx 2.x format of FOV locations from {filename}")
            legacy_format = False
            atomx_format = True
            break
            
    if legacy_format:
        print(f"Using legacy format to read FOVs")
        df = pd.read_csv(os.path.join(offsetsdir, "latest.fovs.csv"), header=None)
        cols = {k: v for k, v in enumerate(
            ["Slide", "X_mm", "Y_mm", "Z_mm", "ZOffset_mm", "ROI", "FOV", "Order"]
            )}
        df = df.rename(columns=cols)
    if atomx_format:
        df = pd.read_csv(os.path.join(offsetsdir, filename))
        # FOV  x_global_px  y_global_px  x_global_mm  y_global_mm
        df = df.rename(columns={'x_global_mm':'X_mm', 'y_global_mm':'Y_mm'})
        df.drop(['x_global_px', 'y_global_px'], axis=1, inplace=True, errors='ignore')
        df.insert(0, 'Slide', 1) # hardcoded 
        df['Z_mm'] = 0.0 # hardcoded
        z_mm_index = df.columns.get_loc('Z_mm')
        df.insert(z_mm_index + 1, 'ZOffset_mm', -2.0) # hardcoded
        df.insert(z_mm_index + 2, 'ROI', 0) # hardcoded
        df.insert(z_mm_index + 3, 'Order', range(len(df))) # hardcoded
        # align columns to legacy format to put "FOV" before "Order" for consistency
        if 'FOV' in df.columns:
            fov_column = df.pop('FOV')
            order_index = df.columns.get_loc('Order')
            df.insert(order_index - 1, 'FOV', fov_column)

    return df

def _resize(image):
    if image.ndim == 2:
        output_shape = (
            max(1, image.shape[0]//2),
            max(1, image.shape[1]//2),
        )
    elif image.ndim == 3:
        output_shape = (
            image.shape[0],
            max(1, image.shape[1]//2),
            max(1, image.shape[2]//2),
        )
    else:
        raise ValueError(f"Unsupported image ndim for pyramid resize: {image.ndim}")
    return resize(
        image,
        output_shape=output_shape,
        order=0,
        preserve_range=True,
        anti_aliasing=False
    )

def _resize_rgb(image):
    return resize(
        image,
        output_shape=(
            max(1, image.shape[0]//2),
            max(1, image.shape[1]//2),
            image.shape[2],
        ),
        order=0,
        preserve_range=True,
        anti_aliasing=False
    )

def _resize_rgb_3d(image):
    return resize(
        image,
        output_shape=(
            image.shape[0],
            max(1, image.shape[1]//2),
            max(1, image.shape[2]//2),
            image.shape[3],
        ),
        order=0,
        preserve_range=True,
        anti_aliasing=False
    )

def _downsample_chunks(chunks, preserve_axes=None, preserve_last=False):
    if preserve_axes is None:
        preserve_axes = set()
    new_chunks = []
    for axis, axis_chunks in enumerate(chunks):
        if axis in preserve_axes or (preserve_last and axis == len(chunks) - 1):
            new_chunks.append(tuple(axis_chunks))
            continue
        new_chunks.append(tuple(max(1, chunk//2) for chunk in axis_chunks))
    return tuple(new_chunks)

def _multiscale_datasets(levels, dimensions, um_per_px, z_step_um):
    """OME-NGFF dataset entries, one per pyramid level, at doubling scale."""
    datasets = []
    pyramid_scale = 1
    for i in range(levels):
        if dimensions == ['z', 'y', 'x']:
            scale = [z_step_um, um_per_px*pyramid_scale, um_per_px*pyramid_scale]
        else:
            scale = [um_per_px*pyramid_scale]*len(dimensions)
        datasets.append({'path': str(i),
                         'coordinateTransformations': [{'type': 'scale',
                          'scale': scale}]})
        pyramid_scale *= 2
    return datasets


def write_pyramid_by_plane(plane_builder, shape, chunks, dtype, scale_dict,
                           store, path):
    """Write a 3D multiscale pyramid one z plane at a time.

    ``plane_builder(z_index)`` returns the ``(1, y, x)`` dask array for a plane.

    Assembling a whole volume as a single dask array is what the obvious
    implementation does, but dask adds a graph layer per ``__setitem__``, so the
    graph grows with tiles x chunks across the *entire* volume. For a 200-FOV
    slide at 8 z planes that is ~3.1M graph keys against ~49K for the same slide
    in 2D, and it exhausts 60 GB while still assembling -- before a single chunk
    is written. Building one plane at a time holds the graph at exactly the 2D
    size no matter how many planes there are, and hands each finished plane
    straight to zarr.

    Pyramid levels are likewise derived by reading the previous level back from
    zarr rather than chaining ``map_blocks`` onto the in-memory graph, so each
    level starts from a shallow graph instead of re-deriving every level below it.
    """
    n_planes, height, width = shape
    levels = max(1, math.floor(math.log2(max(height, width)/256)))
    um_per_px = scale_dict["um_per_px"]
    z_step_um = scale_dict.get("z_step_um", DEFAULT_Z_STEP_UM)
    dimensions = ['z', 'y', 'x']

    print(f"Writing {path} multiscale output to zarr.")
    grp = zarr.open(store, mode='a')
    base = grp.require_dataset(f"{path}/0", shape=shape, chunks=chunks,
                               dtype=dtype, dimension_separator="/",
                               exact=False, overwrite=True)
    print(f"Writing level 1 of {levels}, shape: {shape}, chunksize: {chunks}")
    for z_index in range(n_planes):
        plane = plane_builder(z_index)
        da.to_zarr(plane, base,
                   region=(slice(z_index, z_index + 1), slice(None), slice(None)))

    for i in range(1, levels):
        previous = da.from_zarr(store, component=f"{path}/{i-1}")
        new_chunks = _downsample_chunks(previous.chunks, preserve_axes={0},
                                        preserve_last=False)
        level = previous.map_blocks(_resize, dtype=previous.dtype,
                                    chunks=new_chunks)
        print(f"Writing level {i+1} of {levels}, shape: {level.shape}, "
              f"chunksize: {level.chunksize}")
        level.to_zarr(store, component=f"{path}/{i}", overwrite=True,
                      write_empty_chunks=False, dimension_separator="/")

    grp = zarr.open(store, mode='r+')
    grp[path].attrs['multiscales'] = [{
        'axes': [{'name': dim, 'type': 'space', 'unit': 'micrometer'}
                 for dim in dimensions],
        'datasets': _multiscale_datasets(levels, dimensions, um_per_px, z_step_um),
        'type': 'resize'
        }]


def write_pyramid(image, scale_dict, store, path):
    spatial_shape = image.shape[-3:-1] if path == "composite" else image.shape[-2:]
    PYRAMID_LEVELS = max(1, math.floor(math.log2(max(spatial_shape)/256)))
    um_per_px = scale_dict["um_per_px"]
    z_step_um = scale_dict.get("z_step_um", DEFAULT_Z_STEP_UM)
    pyramid_scale = 1
    if path == "composite" and image.ndim == 3:
        dimensions = ['y', 'x']
        resize_fn = _resize_rgb
        preserve_axes = set()
        preserve_last = True
    elif path == "composite" and image.ndim == 4:
        dimensions = ['z', 'y', 'x']
        resize_fn = _resize_rgb_3d
        preserve_axes = {0}
        preserve_last = True
    elif image.ndim == 2:
        dimensions = ['y', 'x']
        resize_fn = _resize
        preserve_axes = set()
        preserve_last = False
    elif image.ndim == 3:
        dimensions = ['z', 'y', 'x']
        resize_fn = _resize
        preserve_axes = {0}
        preserve_last = False
    else:
        raise ValueError(f"Unsupported image ndim for pyramid write: {image.ndim}")
    datasets = [{}]*PYRAMID_LEVELS
    print(f"Writing {path} multiscale output to zarr.")
    for i in range(PYRAMID_LEVELS):
        print(f"Writing level {i+1} of {PYRAMID_LEVELS}, shape: {image.shape}, chunksize: {image.chunksize}")
        image.to_zarr(store, component=path+f"/{i}", overwrite=True, write_empty_chunks=False, dimension_separator="/")
        new_chunks = _downsample_chunks(image.chunks, preserve_axes=preserve_axes, preserve_last=preserve_last)
        image = image.map_blocks(resize_fn, dtype=image.dtype, chunks=new_chunks)
        if dimensions == ['z', 'y', 'x']:
            scale = [z_step_um, um_per_px*pyramid_scale, um_per_px*pyramid_scale]
        else:
            scale = [um_per_px*pyramid_scale]*len(dimensions)
        datasets[i] = {'path': str(i), 
                       'coordinateTransformations':[{'type':'scale', 
                       'scale': scale}]} 
        pyramid_scale *= 2
    grp = zarr.open(store, mode = 'r+')
    grp[path].attrs['multiscales'] = [{
        'axes':[{'name': dim, 'type': 'space', 'unit': 'micrometer'} for dim in dimensions],
        'datasets': datasets,
        'type': 'resize'
        }]
    channel_name = os.path.splitext(path)[0]
    # write image intensity stats as omero metadata
    if channel_name not in ['labels', 'composite']:
        window = {}
        print("Calculating contrast limits")
        window['min'], window['max'] = int(da.min(image)), int(da.max(image))
        window['start'],window['end'] = [int(x) for x in da.percentile(image.ravel()[image.ravel()!=0], (0.1, 99.9))]
        if window['start'] - window['end'] == 0:
            if window['end'] == 0:
                print(f"\nWARNING: {channel_name} image is empty!")
                window['end'] = 1000
            else:
                window['start'] = 0
        print(f"Writing omero metadata...\n{str(window)}")
        color = DEFAULT_COLORMAPS[channel_name] if channel_name in DEFAULT_COLORMAPS else DEFAULT_COLORMAPS[None]
        grp[path].attrs['omero'] = {'name':channel_name, 'channels': [{
            'label':channel_name,
            'window': window,
            'color': color
            }]}

def base(fov_offsets, fov_height, fov_width, scale_dict, dash):
    px_per_mm = scale_dict["px_per_mm"]
    if dash:
        top_origin_px = max(fov_offsets['X_mm'])*px_per_mm + fov_height
        left_origin_px = min(fov_offsets['Y_mm'])*px_per_mm
        height = round(top_origin_px - min(fov_offsets['X_mm'])*px_per_mm)
        width = round((max(fov_offsets['Y_mm'])*px_per_mm + fov_width) - left_origin_px)
    else:
        top_origin_px = min(fov_offsets['Y_mm'])*px_per_mm - fov_height
        left_origin_px = max(fov_offsets['X_mm'])*px_per_mm
        height = round(max(fov_offsets['Y_mm'])*px_per_mm - top_origin_px)
        width = round((left_origin_px + fov_width) - min(fov_offsets['X_mm'])*px_per_mm)
    return top_origin_px, left_origin_px, height, width
 
def fov_origin(fov_offsets, fov, top_origin_px, left_origin_px, fov_height, scale_dict, dash):
    px_per_mm = scale_dict["px_per_mm"]
    if dash:
        y = round(top_origin_px - (fov_offsets[fov_offsets['FOV'] == fov].iloc[0, ]["X_mm"]*px_per_mm + fov_height))
        x = round(fov_offsets[fov_offsets['FOV'] == fov].iloc[0, ]["Y_mm"]*px_per_mm - left_origin_px)
    else:
        y = round((fov_offsets[fov_offsets['FOV'] == fov].iloc[0, ]["Y_mm"]*px_per_mm - fov_height) - top_origin_px)
        x = round(left_origin_px - fov_offsets[fov_offsets['FOV'] == fov].iloc[0, ]["X_mm"]*px_per_mm)
    return y, x

def get_scales(tiff_path=None, um_per_px=None, scale=1):
    if um_per_px is None:
        with tifffile.TiffFile(tiff_path) as im:
            try:
                j = read_image_description_metadata(tiff_path)
                assert j is not None
                Magnification, PixelSize_um = j['Magnification'], j['PixelSize_um']
                um_per_px = round(PixelSize_um/Magnification, 4)
                print(f"Reading pixel size and magnification from metadata... scale = {um_per_px:.4f} um/px")
            except:
                im_shape = im.pages[0].shape
                fov_height,fov_width = im_shape[0],im_shape[1]
                dash = (fov_height/fov_width) != 1
                if dash:
                    instrument = 'DASH'
                    um_per_px = DASH_UM_PER_PX
                else:
                    beta = fov_height%133 == fov_width%133 == 0
                    if beta:
                        instrument = 'BETA'
                        um_per_px = BETA_UM_PER_PX
                    else:
                        instrument = 'ALPHA'
                        um_per_px = ALPHA_UM_PER_PX
                print(f"Pixel size and magnification not found in metadata, reverting to {instrument} default: {um_per_px:.4f} um/px.")
    um_per_px = round(um_per_px/scale, 4)
    mm_per_px = um_per_px/1000
    px_per_mm = 1/mm_per_px
    if scale != 1:
        print(f"Scaling by {scale} based on user input...")
        print(f"New scale = {um_per_px:.4f} um/px")
    return {"um_per_px":um_per_px, "mm_per_px":mm_per_px, "px_per_mm":px_per_mm}  
