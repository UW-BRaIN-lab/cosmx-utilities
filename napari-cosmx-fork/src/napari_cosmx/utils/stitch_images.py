#!/usr/bin/env python

from napari_cosmx import DEFAULT_NDIM, DEFAULT_Z_STEP_UM
from napari_cosmx.pairing import pair_np
from napari_cosmx.utils import _stitch as stitch
from napari_cosmx.utils._patterns import get_fov_number, get_zslice_number
from napari_cosmx.utils._stitch import fov_tqdm
from tqdm.auto import tqdm
from importlib.metadata import PackageNotFoundError, version
import argparse
import os
import sys
import re
import numpy as np
import tifffile
import zarr
import dask.array as da
import json

LABELS_2D_PATTERN = re.compile(r"CELLLABELS_F[0-9]+\.TIF")
LABELS_3D_PATTERN = re.compile(r"CELLLABELS_F[0-9]+_Z[0-9]+\.TIF")
MORPHOLOGY_2D_PATTERN = re.compile(r".*C902_P99_N99_F[0-9]+\.TIF")
MORPHOLOGY_3D_PATTERN = re.compile(r".*C902_P99_N99_F[0-9]+_Z[0-9]+\.TIF")


# This fork ships under its own distribution name, so upstream's
# version('napari_cosmx') lookup misses. Try the fork first, then upstream.
PACKAGE_NAMES = ("napari-cosmx-fork", "napari_cosmx")


def _package_version():
    """Installed version string, or "unknown" when no distribution is found
    (e.g. running straight from a source checkout)."""
    for name in PACKAGE_NAMES:
        try:
            return version(name)
        except PackageNotFoundError:
            continue
    return "unknown"


def _integer_resize_nearest(tile, target_shape):
    """Resize 2D tile with nearest-neighbor semantics.

    Uses fast repeat when scale factors are integral; otherwise falls back to
    index remapping that preserves dtype and avoids interpolation artifacts.
    """
    src_h, src_w = tile.shape
    dst_h, dst_w = int(target_shape[0]), int(target_shape[1])
    if (src_h, src_w) == (dst_h, dst_w):
        return tile

    if dst_h % src_h == 0 and dst_w % src_w == 0:
        rep_h = dst_h // src_h
        rep_w = dst_w // src_w
        return np.repeat(np.repeat(tile, rep_h, axis=0), rep_w, axis=1)

    row_idx = np.floor(np.arange(dst_h, dtype=np.float64) * (src_h / dst_h)).astype(np.int64)
    col_idx = np.floor(np.arange(dst_w, dtype=np.float64) * (src_w / dst_w)).astype(np.int64)
    row_idx = np.clip(row_idx, 0, src_h - 1)
    col_idx = np.clip(col_idx, 0, src_w - 1)
    return tile[row_idx[:, None], col_idx[None, :]]


def _metadata_binning_from_first_page(first_page):
    """Extract integer morphology binning from ImageDescription if available."""
    tags = {}
    try:
        for tag in first_page.tags.values():
            tags[tag.name] = tag.value
        if 'ImageDescription' not in tags:
            return None
        payload = json.loads(tags['ImageDescription'])
        channels = payload.get('Channels', [])
        if not channels:
            return None
        values = []
        for channel in channels:
            if 'Binning' in channel:
                values.append(int(channel['Binning']))
        if not values:
            return None
        # if mixed values exist, use largest to avoid undersizing target canvas
        return max(values)
    except Exception:
        return None


def _infer_morphology_binning(first_morphology_tif, label_shape=None):
    """Infer morphology binning from metadata first, then shape ratio fallback."""
    with tifffile.TiffFile(first_morphology_tif) as tif:
        morph_shape = tif.pages[0].shape
        binning = _metadata_binning_from_first_page(tif.pages[0])
        if binning is not None and binning > 1:
            return int(binning), morph_shape

    if label_shape is None:
        return 1, morph_shape

    lh, lw = int(label_shape[0]), int(label_shape[1])
    mh, mw = int(morph_shape[0]), int(morph_shape[1])
    if mh == 0 or mw == 0:
        return 1, morph_shape
    if lh % mh == 0 and lw % mw == 0 and (lh // mh) == (lw // mw):
        factor = lh // mh
        if factor > 1:
            return int(factor), morph_shape
    return 1, morph_shape

def _collect_tiles(root_dir, pattern):
    tile_map = {}
    has_z_tokens = False
    for root, dirs, files in os.walk(root_dir):
        for filename in files:
            if not pattern.match(filename.upper()):
                continue
            path = os.path.join(root, filename)
            fov = get_fov_number(path)
            zslice = get_zslice_number(path)
            has_z_tokens = has_z_tokens or (zslice is not None)
            key = (fov, 0 if zslice is None else zslice)
            tile_map.setdefault(key, []).append(path)
    return tile_map, has_z_tokens

def _collect_label_tiles(inputdir, pattern):
    tile_map = {}
    has_z_tokens = False
    for root, dirs, files in os.walk(inputdir, topdown=True):
        rel = os.path.relpath(root, inputdir)
        if rel == '.':
            dirs[:] = [d for d in dirs if d.upper().startswith('FOV')]
            allowed_here = True
        else:
            allowed_here = rel.split(os.sep)[0].upper().startswith('FOV')
        if not allowed_here:
            continue
        for filename in files:
            if not pattern.match(filename.upper()):
                continue
            path = os.path.join(root, filename)
            fov = get_fov_number(path)
            zslice = get_zslice_number(path)
            has_z_tokens = has_z_tokens or (zslice is not None)
            key = (fov, 0 if zslice is None else zslice)
            tile_map.setdefault(key, []).append(path)
    return tile_map, has_z_tokens

def _collect_segmentation_label_tiles(inputdir, subdirs, pattern):
    """Collect label tiles from explicit ``Segmentation_*`` subdirectories.

    AtoMx writes one ``CellStatsDir/Segmentation_<uuid>_00N`` directory per
    segmentation run, and a run's analysis is only valid against its own
    segmentation. Subdirs are searched in the given order and claim whole FOVs
    first-come-first-served, so the primary segmentation always wins and later
    ones only gap-fill FOVs it does not cover. Any FOV still uncovered falls
    back to the base ``FOV*/`` labels.

    Claiming is per-FOV rather than per (FOV, z) so a 3D FOV always comes from
    a single segmentation instead of being mixed plane-by-plane across versions.
    """
    tile_map = {}
    has_z_tokens = False
    covered_fovs = set()

    for subdir in subdirs:
        search_dir = os.path.join(inputdir, subdir)
        if not os.path.isdir(search_dir):
            sys.exit(f"ERROR: --celllabels-subdir not found: {search_dir}")
        sub_tiles, sub_has_z = _collect_label_tiles(search_dir, pattern)
        new_fovs = {fov for (fov, _z) in sub_tiles} - covered_fovs
        for (fov, zslice), paths in sub_tiles.items():
            if fov in new_fovs:
                tile_map.setdefault((fov, zslice), []).extend(paths)
        if new_fovs:
            has_z_tokens = has_z_tokens or sub_has_z
        covered_fovs |= new_fovs
        print(f"  {subdir}: {len(new_fovs)} FOVs")

    base_tiles, base_has_z = _collect_label_tiles(inputdir, pattern)
    fallback_fovs = {fov for (fov, _z) in base_tiles} - covered_fovs
    if fallback_fovs:
        has_z_tokens = has_z_tokens or base_has_z
        for (fov, zslice), paths in base_tiles.items():
            if fov in fallback_fovs:
                tile_map.setdefault((fov, zslice), []).extend(paths)
        covered_fovs |= fallback_fovs
        print(f"  Base CellLabels (fallback): {len(fallback_fovs)} FOVs")

    print(f"Using CellLabels from: {', '.join(subdirs)} ({len(covered_fovs)} total FOVs)")
    return tile_map, has_z_tokens

def _resolve_tile(tile_map, fov, zslice, label):
    tile_path = tile_map.get((int(fov), int(zslice)), [])
    if len(tile_path) == 0:
        tqdm.write(f"Could not find {label} for FOV {fov}, z {zslice}")
        return None
    if len(tile_path) > 1:
        tqdm.write(f"Multiple {label} files found for FOV {fov}, z {zslice}\nUsing {tile_path[0]}")
    return tile_path[0]

def _detect_z_slices(*tile_maps):
    z_slices = sorted({key[1] for tile_map in tile_maps for key in tile_map})
    return z_slices or [0]

def _map_z_to_morphology(z_slices, morphology_z_slices):
    """Map each written morphology plane to the source z plane to draw it from.

    When both labels and morphology are 3D the planes match up directly, and
    anything unmatched falls back to the nearest available plane — z sets are
    not always contiguous or aligned.

    A single 2D morphology acquisition under 3D labels is written once rather
    than repeated per plane (see morphology_is_2d in main), so this is called
    with just that one plane; napari draws it beneath the volumetric labels.
    """
    available = sorted(morphology_z_slices)
    if not available:
        return {z: None for z in z_slices}
    if len(available) == 1:
        return {z: available[0] for z in z_slices}
    return {
        z: (z if z in available else min(available, key=lambda a: abs(a - z)))
        for z in z_slices
    }

def _pick_imagesdir(inputdir, input_ndim, imagesdir):
    """Locate the morphology directory, tolerating a labels/morphology ndim mismatch.

    Prefers ``Morphology<input_ndim>D`` but falls back to whichever
    ``Morphology*D`` directory exists, since a 3D resegmentation is exported
    alongside the original 2D morphology acquisition.
    """
    if imagesdir is not None:
        return imagesdir
    preferred = os.path.join(inputdir, f"Morphology{input_ndim}D")
    if os.path.isdir(preferred):
        return preferred
    for ndim in (2, 3):
        candidate = os.path.join(inputdir, f"Morphology{ndim}D")
        if os.path.isdir(candidate):
            print(f"No Morphology{input_ndim}D directory; using "
                  f"{os.path.basename(candidate)} instead.")
            return candidate
    sys.exit(
        f"Expected morphology images under {preferred}. "
        "Pass --imagesdir to use a custom morphology location."
    )

def main(args_list=None):
    parser = argparse.ArgumentParser(description='Tile CellLabels and morphology TIFFs.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inputdir",
        help="Required: Path to CellLabels and morphology images.",
        default=".")
    parser.add_argument("--imagesdir",
        help="Optional: Path to morphology images, if different than inputdir.",
        default=None)
    parser.add_argument("--celllabels-subdir",
        help="Optional: Subdirectory within inputdir containing CellLabels.\n"
             "Use to select a specific segmentation version, e.g.\n"
             "Segmentation_uuid_003. Repeat as a comma-separated list to let\n"
             "older versions gap-fill FOVs the first one does not cover.\n"
             "When omitted, searches the base FOV*/ directories of inputdir.",
        default=None)
    parser.add_argument("-o", "--outputdir",
        help="Required: Where to create zarr output.",
        default=".")
    parser.add_argument("-f", "--offsetsdir",
        help="Required: Path to directory location containing a file ending in fov_positions_file.csv.gz (AtoMx SIP exported format), FOV_Locations.csv or legacy format latest.fovs.csv.",
        default=".")
    parser.add_argument("-l", "--labels",
        help="\nOptional: Only stitch labels.",
        action='store_true')
    parser.add_argument("-u", "--umperpx",
        help="Optional: Override image scale in um per pixel.\n"+
        "Instrument-specific values to use:\n-> beta04 = 0.1228",
        default=None,
        type=float)
    parser.add_argument("--input-ndim",
        help="Optional: Dimensionality to detect in CellLabels filenames\n"
             "(CellLabels_F###.tif vs CellLabels_F###_Z###.tif).",
        choices=[2, 3],
        default=DEFAULT_NDIM,
        type=int)
    parser.add_argument("--morphology-ndim",
        help="Optional: Dimensionality to detect in morphology filenames.\n"
             "Defaults to auto-detect, since a 3D resegmentation is exported\n"
             "against the original 2D morphology acquisition. A single 2D\n"
             "morphology plane is stored once, not repeated per z plane.",
        choices=[2, 3],
        default=None,
        type=int)
    parser.add_argument("--output-ndim",
        help="Optional: Dimensionality of stitched output.",
        choices=[2, 3],
        default=DEFAULT_NDIM,
        type=int)
    parser.add_argument("--ndim",
        help=argparse.SUPPRESS,
        choices=[2, 3],
        default=None,
        type=int)
    parser.add_argument("-z", "--zslice",
        help="Optional: Z plane to use for 3D-input to 2D-output mode.",
        default=None,
        type=int)
    parser.add_argument("--z-step-um",
        help="Optional: Spacing between z planes in microns.",
        default=None,
        type=float)
    parser.add_argument("--dotzarr",
        help="\nOptional: Add .zarr extension on multiscale pyramids.",
        action='store_true')
    parser.add_argument("--channel-name",
        help="Override the kit's channel-to-marker mapping. Repeatable.\n"
             "Format: CH=MARKER, where CH is one of B, G, Y, R, U.\n"
             "Example: --channel-name B=AT8 --channel-name G=6E10\n"
             "Used when the MorphologyKit metadata in the TIFFs is wrong or\n"
             "inconsistent across studies for the same actual panel.",
        action='append',
        default=[])
    args = parser.parse_args(args=args_list)

    if args.ndim is not None:
        args.output_ndim = args.ndim

    if args.input_ndim == 2 and args.output_ndim == 3:
        sys.exit("Cannot produce 3D output from 2D input filename pattern.")
    if args.input_ndim == 2 and args.zslice is not None:
        sys.exit("--zslice only applies when --input-ndim 3 and --output-ndim 2.")
    if args.input_ndim == 3 and args.output_ndim == 3 and args.zslice is not None:
        sys.exit("--zslice is only for 3D-input to 2D-output mode.")

    # Check output directory
    if not os.path.exists(args.outputdir):
        print(f"Output path does not exist, creating {args.outputdir}")
        os.mkdir(args.outputdir)
    store = os.path.join(args.outputdir, "images")
    if not os.path.exists(store):
        os.mkdir(store)

    args.imagesdir = _pick_imagesdir(args.inputdir, args.input_ndim, args.imagesdir)

    # Read FOV locations file
    fov_offsets = stitch.offsets(args.offsetsdir)

    labels_pattern = LABELS_3D_PATTERN if args.input_ndim == 3 else LABELS_2D_PATTERN

    if args.celllabels_subdir:
        subdirs = [s.strip() for s in args.celllabels_subdir.split(",") if s.strip()]
        labels_tiles, labels_have_z = _collect_segmentation_label_tiles(
            args.inputdir, subdirs, labels_pattern)
    else:
        labels_tiles, labels_have_z = _collect_label_tiles(args.inputdir, labels_pattern)
    labels_res = [path for paths in labels_tiles.values() for path in paths]

    # Check input directory for images and get image dimensions.
    # Labels define target stitched geometry whenever available.
    label_shape = None
    scale_ref_tif = None
    zstep_ref_tif = None
    if len(labels_res) == 0:
        suffix = "_Z###" if args.input_ndim == 3 else ""
        print(f"No CellLabels_FXXX{suffix}.tif files found at {args.inputdir}")
    else:
        label_ref_tif = labels_res[0]
        label_shape = tifffile.TiffFile(label_ref_tif).pages[0].shape
        scale_ref_tif = label_ref_tif
        zstep_ref_tif = label_ref_tif

    ihc_tiles = {}
    ihc_res = []
    ihc_have_z = False
    if not args.labels:
        if args.morphology_ndim is not None:
            morphology_ndim = args.morphology_ndim
            morphology_pattern = (MORPHOLOGY_3D_PATTERN if morphology_ndim == 3
                                  else MORPHOLOGY_2D_PATTERN)
            ihc_tiles, ihc_have_z = _collect_tiles(args.imagesdir, morphology_pattern)
        else:
            # Auto-detect: labels and morphology can disagree, so try the
            # dimensionality implied by the labels first and fall back.
            morphology_ndim = args.input_ndim
            for candidate in (args.input_ndim, 2 if args.input_ndim == 3 else 3):
                morphology_pattern = (MORPHOLOGY_3D_PATTERN if candidate == 3
                                      else MORPHOLOGY_2D_PATTERN)
                ihc_tiles, ihc_have_z = _collect_tiles(args.imagesdir, morphology_pattern)
                if ihc_tiles:
                    morphology_ndim = candidate
                    break
            if morphology_ndim != args.input_ndim:
                print(f"Detected {morphology_ndim}D morphology images alongside "
                      f"{args.input_ndim}D CellLabels.")
        ihc_res = [path for paths in ihc_tiles.values() for path in paths]
        morphology_binning = 1
        morphology_native_shape = None

        if args.zslice is not None:
            z_string = f"_Z{args.zslice:03}"
        else:
            z_string = ""

        if len(ihc_res) == 0:
            print(f"No _FXXX{z_string}.TIF images found at {args.imagesdir}")
        else:
            morph_ref_tif = ihc_res[0]
            morphology_binning, morphology_native_shape = _infer_morphology_binning(
                morph_ref_tif,
                label_shape=label_shape,
            )
            if morphology_binning > 1:
                print(
                    "Detected morphology binning "
                    f"x{morphology_binning}; resizing morphology tiles "
                    "to label geometry before stitching."
                )
            with tifffile.TiffFile(morph_ref_tif) as im:
                n = len(im.pages)
                if n <= 1:
                    sys.exit("Expecting multi-channel TIFFs")
                if label_shape is None:
                    label_shape = im.pages[0].shape
                if label_shape is None:
                    sys.exit("No images found, exiting.")
            # get morphology kit metadata
                channels = ['B','G','Y','R','U']
                # Default channel names (B: Histone, G: Empty, Y: rRNA, R: GFAP, U: DAPI)
                markers = ['Histone','Empty','rRNA','GFAP','DAPI']
                tif_tags = {}
                try:
                    for tag in im.pages[0].tags.values():
                        tif_tags[tag.name] = tag.value
                    j = json.loads(tif_tags['ImageDescription'])
                    reagents = j['MorphologyKit']['MorphologyReagents']
                    mkit = {}
                    for r in reagents:
                        channel = r['Fluorophore']['ChannelId']
                        target = r['BiologicalTarget'].replace("/", "_")
                        mkit[channel] = target
                    markers = [mkit[c] for c in channels]
                except:
                    pass # channel names left as default

                # Apply --channel-name overrides. The kit's MorphologyReagents
                # metadata is sometimes wrong (e.g. mislabeling 6E10 as CD68/CD3
                # or AT8 as Histone) and even inconsistent across studies, so
                # callers can pin the true panel.
                if args.channel_name:
                    overrides = {}
                    for pair in args.channel_name:
                        ch, sep, name = pair.partition("=")
                        if not sep or ch not in channels:
                            sys.exit(f"Invalid --channel-name {pair!r}: expected CH=MARKER with CH in {channels}")
                        overrides[ch] = name
                    markers = [overrides.get(c, m) for c, m in zip(channels, markers)]
            if scale_ref_tif is None:
                scale_ref_tif = morph_ref_tif
            # z-step comes from morphology metadata when available.
            zstep_ref_tif = morph_ref_tif

    if label_shape is None:
        sys.exit("No labels or morphology images found, exiting.")

    if args.input_ndim == 3 and not (labels_have_z or ihc_have_z):
        sys.exit("Requested --input-ndim 3, but no filenames with _Z### were found.")

    if args.input_ndim == 3:
        # Labels define the output z planes whenever present; morphology may be
        # a single 2D acquisition whose z=0 key would otherwise add a phantom plane.
        detected_z_slices = (_detect_z_slices(labels_tiles) if labels_tiles
                             else _detect_z_slices(ihc_tiles))
        if args.output_ndim == 2:
            if args.zslice is None:
                middle_index = (len(detected_z_slices) - 1) // 2
                z_slices = [detected_z_slices[middle_index]]
                print(f"No --zslice specified; using middle z plane: {z_slices[0]}")
            else:
                z_slices = [args.zslice]
        else:
            z_slices = detected_z_slices
    else:
        z_slices = [0]
        
    fov_height = int(label_shape[0])
    fov_width = int(label_shape[1])
    dash = (fov_height/fov_width) != 1
    if args.umperpx == None:
        scale_dict = stitch.get_scales(tiff_path=scale_ref_tif)
    else:
        scale_dict = stitch.get_scales(um_per_px=args.umperpx)

    if args.z_step_um is None:
        z_step_um = stitch.get_z_step_um(zstep_ref_tif)
        if z_step_um is not None:
            print(f"Reading z-step from image metadata: {z_step_um:.4f} um")
        else:
            z_step_um = DEFAULT_Z_STEP_UM
            print(f"No z-step found in image metadata, using default: {z_step_um:.4f} um")
    else:
        z_step_um = args.z_step_um
    scale_dict['z_step_um'] = z_step_um
    
    top_origin_px, left_origin_px, height, width = stitch.base(
        fov_offsets, fov_height, fov_width, scale_dict, dash)

    array_shape = (height, width) if args.output_ndim == 2 else (len(z_slices), height, width)
    chunks = stitch.CHUNKS if args.output_ndim == 2 else (1,) + stitch.CHUNKS

    # A 2D morphology acquisition under a 3D labels volume is stored once rather
    # than copied to every plane. The copies would be byte-identical and carry no
    # extra information, but morphology dominates a mosaic's size — for a typical
    # slide that is ~27 GiB duplicated 8 times, which alone exceeds Fargate's
    # 200 GiB ephemeral storage ceiling. napari renders the single plane beneath
    # the 3D labels just as well.
    morphology_z_slices = sorted({zslice for (_fov, zslice) in ihc_tiles})
    morphology_is_2d = (args.output_ndim == 3 and len(morphology_z_slices) == 1)
    if morphology_is_2d:
        print(f"Morphology is a single 2D plane (z={morphology_z_slices[0]}); "
              "storing it once instead of repeating it across all "
              f"{len(z_slices)} label planes.")
    morphology_shape = (height, width) if (args.output_ndim == 2 or morphology_is_2d) else array_shape
    morphology_chunks = stitch.CHUNKS if (args.output_ndim == 2 or morphology_is_2d) else chunks

    if len(labels_res) != 0:
        print("Stitching cell segmentation labels.")

        def _build_labels_plane(z_index, plane_shape, plane_chunks):
            """Assemble one z plane's FOV tiles into a dask array."""
            plane = da.zeros(plane_shape, dtype=np.uint32, chunks=plane_chunks)
            zslice = z_slices[z_index]
            for fov in fov_tqdm(fov_offsets['FOV']):
                tile_path = _resolve_tile(labels_tiles, fov, zslice, "CellLabels image")
                if tile_path is None:
                    continue
                tile = tifffile.imread(tile_path).astype(np.uint32)
                pair_np(fov, tile)
                y, x = stitch.fov_origin(fov_offsets, fov, top_origin_px, left_origin_px, fov_height, scale_dict, dash)
                if len(plane_shape) == 2:
                    plane[y:y+tile.shape[0], x:x+tile.shape[1]] = tile
                else:
                    plane[0, y:y+tile.shape[0], x:x+tile.shape[1]] = tile
            return plane

        if args.output_ndim == 2:
            im = _build_labels_plane(0, array_shape, chunks)
            stitch.write_pyramid(im, scale_dict, store=store, path="labels")
        else:
            # One plane at a time: a whole-volume dask graph grows with
            # tiles x chunks across every plane and exhausts memory while
            # still being assembled. See write_pyramid_by_plane.
            stitch.write_pyramid_by_plane(
                lambda z_index: _build_labels_plane(z_index, (1, height, width), chunks),
                shape=array_shape, chunks=chunks, dtype=np.uint32,
                scale_dict=scale_dict, store=store, path="labels")
        #TODO: Add .zarr extension to labels if --dotzarr is used. May not be recognized by previous reader versions.
             # Needs more work before readable by napari-ome-zarr anyway

    print("Saving metadata")
    grp = zarr.open(store, mode = 'a')
    grp.attrs['CosMx'] = {
        'fov_height': fov_height,
        'fov_width': fov_width,
        'fov_offsets': fov_offsets.to_dict(),
        'input_ndim': args.input_ndim,
        'ndim': args.output_ndim,
        # Morphology can be flat while labels are volumetric; the reader uses
        # this to pick 2D vs 3D scale/translate per image layer.
        'morphology_ndim': 2 if (args.output_ndim == 2 or morphology_is_2d) else 3,
        'z_slices': z_slices,
        'z_step_um': z_step_um,
        'scale_um': scale_dict['um_per_px'],
        'version': _package_version()
    }

    if len(ihc_res) != 0:
        write_planes = [morphology_z_slices[0]] if morphology_is_2d else z_slices
        morphology_z_for = _map_z_to_morphology(write_planes, set(morphology_z_slices))
        flat_output = (args.output_ndim == 2 or morphology_is_2d)
        for i in range(n):
            im = da.zeros(morphology_shape, dtype=np.uint16, chunks=morphology_chunks)
            print(f"Stitching images for {markers[i]}.")
            for z_index, zslice in enumerate(write_planes):
                morphology_z = morphology_z_for[zslice]
                for fov in fov_tqdm(fov_offsets['FOV']):
                    tile_path = _resolve_tile(ihc_tiles, fov, morphology_z, "image")
                    if tile_path is None:
                        continue
                    with tifffile.TiffFile(tile_path) as my_tiff:
                        tile = my_tiff.pages[i].asarray()
                    if tile.shape != (fov_height, fov_width):
                        tile = _integer_resize_nearest(tile, (fov_height, fov_width))
                    y, x = stitch.fov_origin(fov_offsets, fov, top_origin_px, left_origin_px, fov_height, scale_dict, dash)
                    if flat_output:
                        im[y:y+tile.shape[0], x:x+tile.shape[1]] = tile
                    else:
                        im[z_index, y:y+tile.shape[0], x:x+tile.shape[1]] = tile
            if args.dotzarr:
                markers[i] += ".zarr"
            stitch.write_pyramid(im, scale_dict, store=store, path=f"{markers[i]}")

if __name__ == '__main__':
    sys.exit(main())