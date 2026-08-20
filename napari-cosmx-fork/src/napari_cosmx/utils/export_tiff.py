#!/usr/bin/env python

import numpy as np
from tifffile import TiffWriter, imread
import zarr
import dask.array as da
import argparse
import importlib
import os
import sys
from scipy import ndimage
from tqdm import tqdm
from skimage.transform import resize
from pathlib import Path
from dask.diagnostics import ProgressBar
import tempfile
from napari_cosmx.pairing import pair
import pandas as pd
import re


def _build_parser():
    parser = argparse.ArgumentParser(
        description=(
            'Export stitched Zarr to OME-TIFF. '
            'Choose what to export using --channels, --proteins, or --all.'
        )
    )
    parser.add_argument("-i", "--inputdir",
        help="Required: Path to existing stitched output.",
        default=".")
    parser.add_argument("-o", "--outputdir",
        help="Required: Path to write OME-TIFF file.",
        default=".")
    parser.add_argument("--filename",
        help="Name for OME-TIFF file, use ome.tif extension.",
        default="cosmx-wsi.ome.tif")
    parser.add_argument("--compression",
        help="Passed to TiffWriter, default is 'zlib'. "+
        "Other options include 'lzma' (smallest), 'lzw', and 'none'",
        default='zlib')
    parser.add_argument("-b", "--batchsize",
        help="Maximum number of channels/proteins per OME-TIFF batch (labels from --segmentation do not count). Recommended = 5 or fewer.\n",
        default=5,
        type=int)
    parser.add_argument("-s", "--segmentation",
        help="\nOptional: Include a segmentation border layer in each export (3D uses per-Z borders).",
        action='store_true')
    parser.add_argument("-c", "--channels",
        help="Optional: Export only specific morphology channels (space- or comma-separated).",
        nargs="*",
        default=None)
    parser.add_argument("-p", "--proteins",
        help="Optional: Export only specific proteins (space- or comma-separated).",
        nargs="*",
        default=None)
    parser.add_argument("--all",
        help="Export all available morphology channels and proteins (large file). Required when not using --channels/--proteins.",
        action='store_true')
    parser.add_argument("--levels",
        help="Optional: Specify number of pyramid levels. Capped at available zarr levels.\n",
        default=8,
        type=int)
    parser.add_argument("--volumetrics",
        help="For 3D inputs, write a single volumetric CZYX OME-TIFF instead of the default per-Z CYX split output.",
        action='store_true')
    parser.add_argument("--z-planes",
        help="Optional: Export only selected source z planes for 3D inputs. Example: --z-planes 1 3 or --z-planes 1,3",
        nargs="*",
        default=None)
    parser.add_argument("-v", "--verbose",
        help="Print verbose output?",
        action='store_true')
    parser.add_argument("--dry-run",
        help="Preview export plan (batches, levels, estimated output file count) and exit without writing TIFF files.",
        action='store_true')
    parser.add_argument("--libvips",
        help="\nOptional: Experimental path using libvips for 2D exports and 3D per-Z split exports. Depends on compatible tifffile/libvips builds and is disabled when --volumetrics is used.",
        action='store_true')
    parser.add_argument("--vipshome",
        help="Optional: Path to libvips binaries. Windows-only helper when vips and DLLs are not already in PATH.",
        default=None)
    parser.add_argument("--vipsconcurrency",
        help="Optional: Specify number of threads for vips.\n",
        default=8,
        type=int)
    return parser


def _normalize_compression(compression):
    if compression is None:
        return None
    if str(compression).lower() == 'none':
        return None
    return compression


def _output_axes(data):
    if data.ndim == 3:
        return 'CYX'
    if data.ndim == 4:
        return 'CZYX'
    raise ValueError(f"Unsupported stacked data shape for OME-TIFF export: {data.shape}")


def _build_item_metadata(names, window, pixelsize, axes, z_step_um=None):
    metadata = {
        'axes': axes,
        'Channel': {'Name': names},
        'PhysicalSizeX': pixelsize,
        'PhysicalSizeXUnit': 'µm',
        'PhysicalSizeY': pixelsize,
        'PhysicalSizeYUnit': 'µm',
        'ContrastLimits': [window['min'], window['max']],
        'Window': {'Start': window['start'], 'End': window['end']}
    }
    if axes == 'CZYX' and z_step_um is not None:
        metadata['PhysicalSizeZ'] = z_step_um
        metadata['PhysicalSizeZUnit'] = 'µm'
    return metadata


def _tiff_resolution_from_pixelsize_um(pixelsize_um):
    if pixelsize_um is None:
        return None
    pixelsize_um = float(pixelsize_um)
    if pixelsize_um <= 0:
        return None
    # TIFF resolution is in pixels per unit; use centimeters for broad tool support.
    px_per_cm = (1.0 / pixelsize_um) * 10000.0
    return (px_per_cm, px_per_cm)


def _libvips_available(vipshome=None):
    if os.name == "nt" and vipshome is not None:
        os.environ['PATH'] = vipshome + ';' + os.environ['PATH']

    try:
        importlib.import_module('pyvips')
    except Exception as exc:
        return False, str(exc)

    return True, None


def _insert_suffix_before_tiff(path, suffix):
    path = Path(path)
    filename = path.name
    lower_name = filename.lower()
    if lower_name.endswith('.ome.tif'):
        base = filename[:-8]
        ext = '.ome.tif'
    else:
        base = path.stem
        ext = ''.join(path.suffixes) if path.suffixes else '.tif'
    return str(path.with_name(f"{base}{suffix}{ext}"))


def _selected_z_indices(z_values, z_planes_arg):
    if z_planes_arg is None:
        return list(range(len(z_values)))

    requested = []
    for token in z_planes_arg:
        for part in str(token).split(','):
            part = part.strip()
            if not part:
                continue
            requested.append(int(part))

    if len(requested) == 0:
        return list(range(len(z_values)))

    z_to_index = {int(z): idx for idx, z in enumerate(z_values)}
    missing = [z for z in requested if z not in z_to_index]
    if missing:
        raise ValueError(
            f"Requested z planes not found: {missing}. Available z planes: {list(map(int, z_values))}"
        )

    selected = []
    seen = set()
    for z in requested:
        if z in seen:
            continue
        seen.add(z)
        selected.append(z_to_index[z])
    return selected


def _infer_pixelsize_from_groups(zarr_array):
    candidate_keys = [k for k in zarr_array.keys() if k not in ['labels', 'protein']]
    if len(candidate_keys) == 0 and 'protein' in zarr_array:
        candidate_keys = [f"protein/{k}" for k in zarr_array['protein'].group_keys()]
    for key in candidate_keys:
        attrs = zarr_array[key].attrs
        if 'multiscales' not in attrs:
            continue
        datasets = attrs['multiscales'][0].get('datasets', [])
        if len(datasets) == 0:
            continue
        scales = datasets[0].get('coordinateTransformations', [{}])[0].get('scale', [])
        if len(scales) >= 2:
            # scale is [z, y, x] for 3D or [y, x] for 2D; use Y spacing.
            return float(scales[-2])
    return None


def _logical_shape_from_ome_sizes(pixels):
    size_c = pixels.get('SizeC')
    size_y = pixels.get('SizeY')
    size_x = pixels.get('SizeX')
    if size_c is None or size_y is None or size_x is None:
        return None

    size_c = int(size_c)
    size_y = int(size_y)
    size_x = int(size_x)
    size_z = int(pixels.get('SizeZ', 1))

    if size_z <= 1:
        return (size_c, size_y, size_x)
    return (size_c, size_z, size_y, size_x)


def verify_ome_tiff(path, *, expected_axes=None, expected_channels=None,
                    expected_shape=None, expected_pyramid_levels=None,
                    expected_pixel_size_um=None, verbose=False):
    """Verify key structural and metadata properties of an exported OME-TIFF.

    Parameters
    ----------
    path : str or Path
        Path to the OME-TIFF file to verify.
    expected_axes : str, optional
        Expected OME axes string, e.g. ``'CYX'`` or ``'CZYX'``.
    expected_channels : list of str, optional
        Expected channel names in order.
    expected_shape : tuple, optional
        Expected full-resolution data shape (C[Z]YX).
    expected_pyramid_levels : int, optional
        Minimum number of pyramid levels expected (IFDs or SubIFDs).
    expected_pixel_size_um : float, optional
        Expected PhysicalSizeX/Y in micrometres (checked within 1 % tolerance).
    verbose : bool, optional
        If True, print a summary of what was found.

    Returns
    -------
    dict
        A report dict with keys ``axes``, ``shape``, ``n_levels``,
        ``channels``, ``pixel_size_um``, ``z_step_um``, ``has_ome_xml``,
        and ``issues`` (list of str describing any discrepancies found).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    from tifffile import TiffFile

    path = str(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    issues = []

    with TiffFile(path) as tif:
        has_ome_xml = tif.is_ome
        if not has_ome_xml:
            issues.append("File does not contain OME-XML metadata.")

        # Axes from OME-XML
        axes = None
        channels = []
        pixel_size_um = None
        z_step_um = None
        shape = None

        if has_ome_xml and tif.ome_metadata:
            import xml.etree.ElementTree as ET
            ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
            try:
                root = ET.fromstring(tif.ome_metadata)
                image = root.find('ome:Image', ns) or root.find('Image')
                if image is not None:
                    pixels = image.find('ome:Pixels', ns) or image.find('Pixels')
                    if pixels is not None:
                        axes = pixels.get('DimensionOrder')
                        shape = _logical_shape_from_ome_sizes(pixels)
                        px = pixels.get('PhysicalSizeX')
                        if px is not None:
                            pixel_size_um = float(px)
                        pz = pixels.get('PhysicalSizeZ')
                        if pz is not None:
                            z_step_um = float(pz)
                        channels = [
                            (ch.get('Name') or ch.get('ID', ''))
                            for ch in (pixels.findall('ome:Channel', ns) or pixels.findall('Channel'))
                        ]
            except Exception as exc:
                issues.append(f"Could not parse OME-XML: {exc}")

        # Shape and pyramid level count from page series
        n_levels = 0
        if tif.series:
            series = tif.series[0]
            if shape is None:
                shape = series.shape
            # Count pyramid levels: base level counts as 1, each SubIFD adds more
            try:
                n_levels = len(series.levels)
            except AttributeError:
                n_levels = 1  # tifffile version without .levels

        # Validate against expectations
        if expected_axes is not None and axes is not None:
            # OME DimensionOrder is written as reversed traversal order (e.g. XYCZT).
            # We compare our logical read order (CYX/CZYX) against the OME convention.
            ome_to_logical = axes[::-1] if axes else axes
            if ome_to_logical != expected_axes:
                issues.append(
                    f"Axes mismatch: expected '{expected_axes}', "
                    f"OME DimensionOrder='{axes}' (logical='{ome_to_logical}')."
                )

        if expected_shape is not None and shape is not None:
            if tuple(shape) != tuple(expected_shape):
                issues.append(
                    f"Shape mismatch: expected {tuple(expected_shape)}, got {tuple(shape)}."
                )

        if expected_channels is not None and channels:
            if list(channels) != list(expected_channels):
                issues.append(
                    f"Channel names mismatch: expected {expected_channels}, got {list(channels)}."
                )

        if expected_pyramid_levels is not None:
            if n_levels < expected_pyramid_levels:
                issues.append(
                    f"Pyramid levels: expected at least {expected_pyramid_levels}, got {n_levels}."
                )

        if expected_pixel_size_um is not None and pixel_size_um is not None:
            rel_err = abs(pixel_size_um - expected_pixel_size_um) / expected_pixel_size_um
            if rel_err > 0.01:
                issues.append(
                    f"Pixel size mismatch: expected {expected_pixel_size_um} µm, "
                    f"got {pixel_size_um} µm (>{rel_err*100:.1f}% error)."
                )

    report = {
        'path': path,
        'has_ome_xml': has_ome_xml,
        'axes': axes,
        'shape': shape,
        'n_levels': n_levels,
        'channels': channels,
        'pixel_size_um': pixel_size_um,
        'z_step_um': z_step_um,
        'issues': issues,
    }

    if verbose:
        print(f"\n--- OME-TIFF verification report ---")
        print(f"  path       : {path}")
        print(f"  has OME-XML: {has_ome_xml}")
        print(f"  axes       : {axes}")
        print(f"  shape      : {shape}")
        print(f"  n_levels   : {n_levels}")
        print(f"  channels   : {channels}")
        print(f"  px size µm : {pixel_size_um}")
        print(f"  z step µm  : {z_step_um}")
        if issues:
            print(f"  ISSUES ({len(issues)}):")
            for issue in issues:
                print(f"    - {issue}")
        else:
            print(f"  OK — no issues found.")
        print()

    return report


def _write_ome_pyramid(path, levels, pyramid_levels, item_metadata, compression,
                       resolution_px_per_cm=None, verbose=False):
    compression = _normalize_compression(compression)
    with TiffWriter(path, bigtiff=True, ome=True) as tif:
        options = dict(
            tile=(1024, 1024),
            metadata=item_metadata,
            compression=compression
        )
        if resolution_px_per_cm is not None:
            options['resolutionunit'] = 'CENTIMETER'
        if pyramid_levels > 1:
            options['subifds'] = pyramid_levels - 1
        base_y = int(levels[0].shape[-2])
        base_x = int(levels[0].shape[-1])
        for i in tqdm(range(pyramid_levels), ncols=60, smoothing=1):
            write_options = dict(options)
            if resolution_px_per_cm is not None:
                cur_y = int(levels[i].shape[-2])
                cur_x = int(levels[i].shape[-1])
                downsample_x = max(1.0, base_x / max(1, cur_x))
                downsample_y = max(1.0, base_y / max(1, cur_y))
                write_options['resolution'] = (
                    resolution_px_per_cm[0] / downsample_x,
                    resolution_px_per_cm[1] / downsample_y,
                )
            if i > 0:
                write_options['subfiletype'] = 1
            tif.write(
                data=levels[i],
                **write_options
            )
            if i == 0:
                del options['metadata']
                if 'subifds' in options:
                    del options['subifds']
    if verbose:
        print(f"Wrote {path}")


def _write_ome_with_libvips(path, data, item_metadata, pyramid_levels,
                            outputdir, vipshome=None, vipsconcurrency=8,
                            resolution_px_per_cm=None, verbose=False):
    available, error_message = _libvips_available(vipshome=vipshome)
    if not available:
        raise RuntimeError(
            "pyvips/libvips is unavailable: "
            f"{error_message}"
        )

    os.environ['VIPS_PROGRESS'] = "1"
    os.environ['VIPS_CONCURRENCY'] = str(vipsconcurrency)
    import pyvips # must import after updating path

    tmpdirname = tempfile.mkdtemp(dir=outputdir)
    tmptif = os.path.join(tmpdirname, 'tmp.ome.tif')

    def _align16(value):
        value = max(16, int(value))
        return ((value + 15) // 16) * 16

    # Keep temporary tile dimensions valid for small images.
    tmp_tile_h = _align16(min(1024, int(data.shape[-2])))
    tmp_tile_w = _align16(min(1024, int(data.shape[-1])))

    with TiffWriter(tmptif, bigtiff=True, ome=True) as tif:
        print(f"Writing uncompressed empty tiff, shape: {data.shape}.\nThis could take a while.")
        options = dict(
            tile=(tmp_tile_h, tmp_tile_w),
            metadata=item_metadata,
            compression=None,
        )
        tif.write(
            data=None,
            dtype=data.dtype,
            shape=data.shape,
            **options
        )

    # Open TIFF as Zarr and stream pixels in.
    store = imread(tmptif, mode='r+', aszarr=True)
    z = zarr.open(store, mode='r+')
    if data.shape[0] == 1 and len(data.shape) == 3:
        print("Only 1 layer for this ome-tiff file. Reshaping.")
        data = data.squeeze(axis=0)
    print("Writing data to tiff as zarr")
    with ProgressBar():
        da.to_zarr(arr=data, url=z)
    store.close()

    print(f"Creating pyramidal OME-TIFF at {path}")
    image = pyvips.Image.new_from_file(
        tmptif,
        n=int(data.shape[0]) if len(data.shape) == 3 else 1,
    )
    tile_height = _align16(data.shape[-2] // 2**(pyramid_levels - 1))
    tile_width = _align16(data.shape[-1] // 2**(pyramid_levels - 1))
    save_kwargs = dict(
        compression="deflate", bigtiff=True,
        tile=True, tile_width=tile_width, tile_height=tile_height,
        pyramid=True, subifd=True
    )
    if resolution_px_per_cm is not None:
        # pyvips xres/yres are in pixels/mm; convert from pixels/cm.
        px_per_mm = resolution_px_per_cm[0] / 10.0
        save_kwargs['xres'] = px_per_mm
        save_kwargs['yres'] = px_per_mm
        save_kwargs['resunit'] = 'cm'
    image.tiffsave(path, **save_kwargs)
    if verbose:
        print(f"Wrote {path} via libvips")

def split_list_element(input_list):
    """
    Splits a list element by spaces or commas.

    Args:
        input_list: A list of length 1 containing a string.

    Returns:
        A list of strings, or the original list if the input is invalid.
        Returns an empty list if the input list is empty.
    """

    if not isinstance(input_list, list):
        return "Input must be a list"  # Or raise an exception

    if not input_list:
        return [] #Handle empty list

    if len(input_list) != 1:
        return "Input list must have length 1"  # Or raise an exception

    element = input_list[0]

    if "Cathepsin B" in element:
        element = element.replace("Cathepsin B", "CATHEPSIN_B_PLACEHOLDER") # Use a placeholder
        split_elements = re.split(r'[ ,]+', element)
        cleaned_elements = [item.replace("CATHEPSIN_B_PLACEHOLDER", "Cathepsin B") if item == "CATHEPSIN_B_PLACEHOLDER" else item for item in split_elements if item]
    else:
        split_elements = re.split(r'[ ,]+', element)
        cleaned_elements = [item for item in split_elements if item]

    return cleaned_elements


def _parse_cli_values(values):
    if values is None:
        return None
    if not isinstance(values, list):
        return split_list_element(values)
    if not values:
        return []

    cleaned_values = []
    for value in values:
        split_value = split_list_element([str(value)])
        if isinstance(split_value, list):
            cleaned_values.extend(split_value)
    return cleaned_values

class BatchStorage:
    def __init__(self, batch_size):
        self.batch_size = batch_size
        self.storage = {}
        self.current_key = None
        self.item_count = 0
        self.labels = None  # Store labels separately

    def set_labels(self, labels):
        """Sets the labels to be added to each batch."""
        self.labels = labels

    def add_item(self, item, new_batch=False):
        """Adds an item to the storage, including handling labels."""

        if not self.storage:
            self.current_key = 0
            self.storage[self.current_key] = []
            self.item_count = 0
            if self.labels is not None:  # Add labels to the first batch if available
                self.storage[self.current_key].append(self.labels)

        if new_batch or self.item_count >= self.batch_size:
            self.current_key = max(self.storage.keys()) + 1 if self.storage else 0
            self.storage[self.current_key] = []
            self.item_count = 0
            if self.labels is not None: # Add labels to the new batch if available
                self.storage[self.current_key].append(self.labels)

        self.storage[self.current_key].append(item)
        self.item_count += 1

    def get_batch(self, key):
        return self.storage.get(key)

    def __str__(self):
        return str(self.storage)

def _edges(x):
    x = np.asarray(x)
    if x.ndim < 2:
        return np.zeros_like(x, dtype=np.uint8)

    kernel_2d = np.ones((3, 3), dtype=np.int8)
    kernel_2d[1, 1] = -8
    # Apply 2D boundary detection over the last two axes; keep earlier axes independent.
    kernel = kernel_2d.reshape((1,) * (x.ndim - 2) + kernel_2d.shape)
    arr = ndimage.convolve(x, kernel, output=np.int32)
    arr[arr != 0] = 1
    return arr.astype(np.uint8)

def _scale_edges(x):
    # Scale the binary [0, 1] to [100, 10000]
    x_scaled = x * 9900 + 100  # Linear scaling
    return x_scaled.astype('uint16') #return as uint16, as it is larger than uint8

def Parse():
    parser = _build_parser()
    args = parser.parse_args()

    return(args)

def main(args_list=None):
    if args_list is None or isinstance(args_list, (list, tuple)):
        parser = _build_parser()
        args = parser.parse_args(args=args_list)
    else:
        # Backward compatible support for tests or callers that pass a namespace-like object.
        args = args_list
        if not hasattr(args, 'volumetrics'):
            args.volumetrics = False
        if not hasattr(args, 'libvips'):
            args.libvips = False
        if not hasattr(args, 'vipshome'):
            args.vipshome = None

    # Check output directory
    if not os.path.exists(args.outputdir):
        print(f"Output path does not exist, creating {args.outputdir}")
        os.mkdir(args.outputdir)
    store = os.path.join(args.inputdir, "images")
    if not os.path.exists(store):
        sys.exit(f"Could not find images directory at {args.inputdir}")

    libvips_enabled = bool(getattr(args, 'libvips', False))
    if libvips_enabled:
        available, error_message = _libvips_available(vipshome=getattr(args, 'vipshome', None))
        if not available:
            print(
                "--libvips requested but pyvips/libvips is unavailable; "
                f"using tifffile writer instead. Details: {error_message}"
            )
            libvips_enabled = False
            args.libvips = False

    # Top-level input, checking, initializing
    zarr_array = zarr.open(store, mode='r')
    cosmx_meta = zarr_array.attrs.get('CosMx', {})
    if 'scale_um' in cosmx_meta:
        pixelsize = cosmx_meta['scale_um']
    else:
        pixelsize = _infer_pixelsize_from_groups(zarr_array)
        if pixelsize is None:
            sys.exit("Could not find scaling information from top-level zarr or group multiscales. Error 1.")
        if args.verbose:
            print(f"Using inferred pixel size from group multiscales: {pixelsize} um/px")
    

    batches = BatchStorage(args.batchsize) # intialize

    ### Segmentation
    has_labels = False
    idx = 0
    if args.segmentation:
        for item in zarr_array.items():
            if item[0] == "labels":
                has_labels = True
                idx = 1
        if not has_labels:
            sys.exit(f"Error. Segmentation labels were requested but the directory 'labels' could not be found.")
        else:
            if args.verbose:
                print(f"Adding segmentation to each batch.")
            batches.set_labels('labels')

    ### Channels
    if args.channels is not None and not args.all:
        valid_channels = [key for key in zarr_array.keys() if key not in ["labels", "protein"]]
        if len(args.channels) == 0:
            channels_to_process = valid_channels
        else:
            cleaned_list = _parse_cli_values(args.channels)
            channels_to_process = []
            for x in cleaned_list:
                if x in valid_channels:
                    channels_to_process.append(x)
                else:
                    print(f"Warning! {x} is not a valid channel and will be ignored.")
        if len(channels_to_process) == 0:
            print(f"--channels were requested but no valid channels were found in the zarr store.") 
        else:
            [batches.add_item(x) for x in channels_to_process]
    
    ### Proteins
    has_proteins = False
    for item in zarr_array.items():
        if item[0] == "protein":
            has_proteins = True

    if args.proteins is not None and has_proteins and not args.all:
        valid_proteins = ['protein/' + x for x in zarr_array['protein'].group_keys()]
        if len(args.proteins) == 0:
            # requests all proteins
            proteins_to_process = valid_proteins
        else:
            cleaned_list = ['protein/' + x for x in _parse_cli_values(args.proteins)]
            proteins_to_process = []
            for x in cleaned_list:
                if x in valid_proteins:
                    proteins_to_process.append(x)
                else:
                    print(f"Warning! {x} is not a valid protein and will be ignored (check spelling).")
        if len(proteins_to_process) == 0:
            print(f"--proteins were specified but no valid protein names were given.")
        else:
            [batches.add_item(x) for x in proteins_to_process]
    elif args.proteins and not args.all and args.verbose:
        print(f"Requested protein export but no protein zarr found. Ignoring.")

    # Explicit export-all mode.
    if args.all:
        if (args.channels is not None or args.proteins is not None) and args.verbose:
            print("--all was provided; ignoring --channels/--proteins filters and exporting all available items.")
        default_items = [key for key in zarr_array.keys() if key not in ["labels", "protein"]]
        if has_proteins:
            default_items.extend([f"protein/{x}" for x in zarr_array['protein'].group_keys()])

        if default_items:
            for item_name in default_items:
                batches.add_item(item_name)
            if args.verbose:
                print(
                    f"Exporting all available items "
                    f"({len(default_items)} total)."
                )

    if len(batches.storage) == 0:
        print(
            "No exportable items were selected, so nothing was written. "
            "Use --channels/--proteins to choose content, or pass --all to export all available items."
        )
        return 0

    if args.dry_run:
        print("Dry run: no files will be written.")
        total_estimated_files = 0
        for key, items in batches.storage.items():
            if len(items) <= idx:
                print(f"  batch {key}: skipped (no data items found)")
                continue

            first_item = items[idx]
            levels_available = 0
            pyramid_levels = 0
            estimated_files = 1
            try:
                attrs = zarr_array[first_item].attrs
                datasets = attrs["multiscales"][0]["datasets"]
                levels_available = len(datasets)
                requested_levels = args.levels if args.levels is not None else levels_available
                pyramid_levels = min(max(1, requested_levels), levels_available)

                base = da.from_zarr(store + f"/{first_item}", component=datasets[0]["path"])
                if base.ndim == 3:
                    z_values = cosmx_meta.get('z_slices', None)
                    if z_values is None or len(z_values) != int(base.shape[0]):
                        z_values = list(range(int(base.shape[0])))
                    selected_indices = _selected_z_indices(z_values, args.z_planes)
                    if not args.volumetrics:
                        estimated_files = max(1, len(selected_indices))
            except Exception:
                pass

            total_estimated_files += estimated_files
            content_items = items[idx:]
            has_seg_in_batch = idx == 1 and len(items) > 0 and items[0] == 'labels'
            print(
                f"  batch {key}: items={len(content_items)}, segmentation={has_seg_in_batch}, "
                f"levels={pyramid_levels if pyramid_levels else 'unknown'}, "
                f"estimated_files={estimated_files}"
            )
            if args.verbose:
                print(f"    contents: {content_items}")

        print(f"Estimated output files: {total_estimated_files}")
        return 0

    ### Processing
    for key, items in batches.storage.items():
        if args.verbose:
            print(f"Processing batch number {str(key)}")
        # get attributes from first item (skip labels if applicable)
        attrs = zarr_array[items[idx]].attrs
        datasets = attrs["multiscales"][0]["datasets"] # dimensions of each level
        omero = attrs["omero"]
        window = omero['channels'][0]['window']
        names = [x.replace(".zarr", "").replace("protein/", "") for x in items]
        
        available_levels = len(datasets)
        requested_levels = args.levels if args.levels is not None else available_levels
        pyramid_levels = min(max(1, requested_levels), available_levels)
        if args.verbose and requested_levels > available_levels:
            print(
                f"Requested {requested_levels} pyramid levels but only {available_levels} are available. "
                f"Using {pyramid_levels}."
            )

        levels = []
        for d in datasets:
            arrays = []
            for i in items:
                i_array = da.from_zarr(store + f"/{i}", component=d["path"])
                if has_labels and i == "labels":
                    i_array = i_array.map_blocks(_edges, dtype=np.uint8).map_blocks(_scale_edges, dtype=np.uint16) # binary mask, bounded
                arrays.append(i_array)
            stacked_array = da.stack(arrays)
            levels.append(stacked_array)

        data = levels[0]
        axes = _output_axes(data)
        z_step_um = cosmx_meta.get('z_step_um', None)
        resolution_px_per_cm = _tiff_resolution_from_pixelsize_um(pixelsize)
        item_metadata = _build_item_metadata(
            names=names,
            window=window,
            pixelsize=pixelsize,
            axes=axes,
            z_step_um=z_step_um
        )

        z_values = None
        if data.ndim == 4:
            z_values = cosmx_meta.get('z_slices', None)
            if z_values is None or len(z_values) != int(data.shape[1]):
                z_values = list(range(int(data.shape[1])))
            try:
                selected_indices = _selected_z_indices(z_values, args.z_planes)
            except ValueError as exc:
                sys.exit(str(exc))

            if len(selected_indices) != int(data.shape[1]):
                levels = [level[:, selected_indices, :, :] for level in levels]
                z_values = [z_values[i] for i in selected_indices]
                data = levels[0]
                axes = _output_axes(data)
                item_metadata = _build_item_metadata(
                    names=names,
                    window=window,
                    pixelsize=pixelsize,
                    axes=axes,
                    z_step_um=z_step_um
                )

        path = os.path.join(args.outputdir, f"batch_{key}_{args.filename}")
        if args.verbose:
            print(f"Writing {path}.")

        split_3d_for_export = (data.ndim == 4) and (not args.volumetrics)
        if data.ndim == 4 and split_3d_for_export:
            if args.verbose:
                print("3D stack detected; exporting one CYX OME-TIFF per z plane.")
            use_libvips_split = libvips_enabled
            for z_index, z_val in enumerate(z_values):
                z_levels = [level[:, z_index, :, :] for level in levels[:pyramid_levels]]
                z_metadata = _build_item_metadata(
                    names=[f"{n} [z={z_val}]" for n in names],
                    window=window,
                    pixelsize=pixelsize,
                    axes='CYX',
                    z_step_um=None
                )
                z_path = _insert_suffix_before_tiff(path, f"_z{int(z_val):03d}")
                if use_libvips_split:
                    _write_ome_with_libvips(
                        path=z_path,
                        data=z_levels[0],
                        item_metadata=z_metadata,
                        pyramid_levels=pyramid_levels,
                        outputdir=args.outputdir,
                        vipshome=args.vipshome,
                        vipsconcurrency=args.vipsconcurrency,
                        resolution_px_per_cm=resolution_px_per_cm,
                        verbose=args.verbose,
                    )
                    validation = verify_ome_tiff(
                        z_path,
                        expected_shape=tuple(int(x) for x in z_levels[0].shape),
                        expected_channels=z_metadata.get('Channel', {}).get('Name'),
                        expected_pixel_size_um=pixelsize,
                        verbose=args.verbose,
                    )
                    if validation['issues']:
                        raise SystemExit(
                            "libvips split-z export validation failed for "
                            f"{z_path}: {validation['issues']}. "
                            "Please rerun without --libvips to use the tifffile writer."
                        )
                else:
                    _write_ome_pyramid(
                        path=z_path,
                        levels=z_levels,
                        pyramid_levels=pyramid_levels,
                        item_metadata=z_metadata,
                        compression=args.compression,
                        resolution_px_per_cm=resolution_px_per_cm,
                        verbose=args.verbose
                    )
            continue

        use_libvips = libvips_enabled
        if args.volumetrics and use_libvips:
            # Keep volumetric output on tifffile path for metadata consistency.
            if args.verbose:
                print("--volumetrics requested; disabling libvips and using tifffile writer.")
            use_libvips = False
        if data.ndim == 4 and use_libvips:
            # libvips path currently assumes 2D channel stacks; keep 3D output on tifffile path.
            if args.verbose:
                print("3D stack detected; using tifffile writer (libvips path is 2D-only).")
            use_libvips = False

        if not use_libvips:
            _write_ome_pyramid(
                path=path,
                levels=levels,
                pyramid_levels=pyramid_levels,
                item_metadata=item_metadata,
                compression=args.compression,
                resolution_px_per_cm=resolution_px_per_cm,
                verbose=args.verbose
            )
        else:
            _write_ome_with_libvips(
                path=path,
                data=data,
                item_metadata=item_metadata,
                pyramid_levels=pyramid_levels,
                outputdir=args.outputdir,
                vipshome=args.vipshome,
                vipsconcurrency=args.vipsconcurrency,
                resolution_px_per_cm=resolution_px_per_cm,
                verbose=args.verbose,
            )

if __name__ == '__main__':
    sys.exit(main())