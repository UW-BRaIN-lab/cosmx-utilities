import napari
from napari.experimental import link_layers
from napari_cosmx import DASH_MM_PER_PX, ALPHA_MM_PER_PX, BETA_MM_PER_PX, DEFAULT_COLORMAPS, DEFAULT_NDIM, DEFAULT_PROTEIN_Z_STEP_UM, DEFAULT_Z_STEP_UM, OTHER_KEYS
from superqt.utils import qdebounced
from skimage import io
from skimage.segmentation import find_boundaries
from scipy import ndimage
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import os
import warnings
import re
import glob
import anndata as ad
import zarr
import vaex
import dask.array as da
import pickle
from napari_cosmx.pairing import unpair, pair
from napari_cosmx._colors import categorical_color_map
from napari_cosmx._dock_widget import GeminiQWidget
from sklearn import preprocessing
from skimage.draw import polygon, ellipse
from scipy import ndimage as ndi
from napari.utils.colormaps import AVAILABLE_COLORMAPS, label_colormap, color_dict_to_colormap
from napari.utils.colormaps.standardize_color import transform_color
from vispy.color.colormap import Colormap
from inspect import signature
from os import listdir
from importlib.metadata import version

def _find_boundaries_by_slice(labels):
    if labels.ndim == 2:
        return find_boundaries(labels)
    if labels.ndim == 3:
        return np.stack([find_boundaries(labels[index]) for index in range(labels.shape[0])])
    raise ValueError(f"Unsupported labels ndim for boundary extraction: {labels.ndim}")


def _normalize_window(window):
    """Return a napari-safe contrast window derived from Omero-style values.

    This helper preserves Omero as the source of truth but repairs edge cases
    (for example min==max on constant images) so contrast limit ranges are
    strictly increasing and accepted by newer napari versions.
    """
    normalized = {
        'min': int(window['min']),
        'max': int(window['max']),
        'start': int(window['start']),
        'end': int(window['end']),
    }

    if normalized['max'] <= normalized['min']:
        if normalized['max'] <= 0:
            normalized['min'], normalized['max'] = 0, 1
        else:
            normalized['min'] = 0
        if normalized['max'] <= normalized['min']:
            normalized['max'] = normalized['min'] + 1

    if normalized['end'] <= normalized['start']:
        if normalized['end'] <= 0:
            normalized['start'], normalized['end'] = 0, 1
        else:
            normalized['start'] = 0
        if normalized['end'] <= normalized['start']:
            normalized['end'] = normalized['start'] + 1

    normalized['min'] = min(normalized['min'], normalized['start'])
    normalized['max'] = max(normalized['max'], normalized['end'])
    if normalized['max'] <= normalized['min']:
        normalized['max'] = normalized['min'] + 1

    return normalized

class Gemini:
    """Initialize new instance and launch viewer

    Args:
        path (str): path to adata, transcripts, and/or images or .h5ad path
        viewer (napari.viewer.Viewer): If None, napari will be launched.
    """
    def __init__(self, path, viewer=None):
        assert os.path.exists(path), f"Could not find {path}"
        self._rotate = 0
        self.folder = path
        self.adata = None
        self.fov_labels = None
        self.segmentation_layer = None
        self.cells_layer = None
        self.color_cells_layer = None
        if path.endswith(".h5ad"):
            self.adata = ad.read_h5ad(path)
            self.folder = os.path.dirname(path)
        files = next(os.walk(self.folder))[2]
        if self.adata is None:
            res = [f for f in files if f.endswith('.h5ad')]
            if len(res) != 1:
                print("Could not find AnnData .h5ad file")
            else:
                self.adata = ad.read_h5ad(os.path.join(self.folder, res[0]))
        assert os.path.exists(os.path.join(self.folder, "images")), f"No images directory found at {self.folder}"
        self.grp = zarr.open(os.path.join(self.folder, "images"), mode = 'r+',)
        cosmx_attrs = self.grp.attrs['CosMx']
        self.ndim = int(cosmx_attrs.get('ndim', DEFAULT_NDIM))
        self.is_3d = self.ndim == 3
        self.z_slices = list(cosmx_attrs.get('z_slices', [0]))
        self.z_slice_to_index = {int(zslice): index for index, zslice in enumerate(self.z_slices)}
        # Main CosMx images/transcripts use the unified z_step_um key; protein can
        # still carry its own z spacing in the protein group metadata and override it below.
        self.z_step_um = float(cosmx_attrs.get('z_step_um', cosmx_attrs.get('rna_z_step_um', DEFAULT_Z_STEP_UM)))
        self.rna_z_step_um = self.z_step_um
        self.protein_z_step_um = float(cosmx_attrs.get('protein_z_step_um', self.z_step_um))
        self.rna_z_step_mm = self.rna_z_step_um / 1000
        self.protein_z_step_mm = self.protein_z_step_um / 1000
        self.fov_height = cosmx_attrs['fov_height']
        self.fov_width = cosmx_attrs['fov_width']
        self.fov_offsets = pd.DataFrame.from_dict(cosmx_attrs['fov_offsets'])
        self.dash = (self.fov_height/self.fov_width) != 1
        self.beta = self.fov_height%133 == self.fov_width%133 == 0
        # try to get scale from CosMx metadata, otherwise revert to instrument defaults
        if 'scale_um' in self.grp.attrs['CosMx']:
            um_per_px = self.grp.attrs['CosMx']['scale_um']
            print(f"Reading CosMx scale from .attrs: {um_per_px:.4f} um/px")
            self.mm_per_px = um_per_px/1000
        elif 'scale_mm' in self.grp.attrs['CosMx']:
            self.mm_per_px = self.grp.attrs['CosMx']['scale_mm']
            print(f"Reading CosMx scale from .attrs: {self.mm_per_px} mm/px")
        else:
            if self.dash:
                self.mm_per_px = DASH_MM_PER_PX
                self.instrument = 'DASH'
            elif self.beta:
                self.mm_per_px = BETA_MM_PER_PX
                self.instrument = 'BETA'
            else:
                self.mm_per_px = ALPHA_MM_PER_PX
                self.instrument = 'ALPHA'
            print(f"No CosMx scale found in .attrs, using {self.instrument} default: {(self.mm_per_px*1000):.4f} um/px")
        
        self._channels = sorted(set([s.replace(r'.zarr', '') for s in self.grp.group_keys()]) - set(OTHER_KEYS))
        self._proteins = []
        self.expr_scale = 1
        self.targets = self._load_targets()
        self._load_protein_state()
        if (self.adata is not None) and ('name' in self.adata.uns):
            self.name = self.adata.uns['name']
        else:
            folder = os.path.normpath(self.folder)
            slide = os.path.basename(folder)
            parent = os.path.basename(os.path.dirname(folder))
            # Strip trailing date/number suffixes (e.g. "22426_06_04_2026_12_30_21_218")
            short_parent = re.sub(r'\d[\d_]*$', '', parent).rstrip('_')
            self.name = f"{short_parent}/{slide}" if short_parent else slide

        if self.adata is None:
            # Accept either legacy *_metdata.csv or gzipped *_metadata_file.csv.gz.
            # Prefer plain CSV when both naming styles are present.
            metadata_csv = sorted(glob.glob(os.path.join(self.folder, "*_metadata.csv")))
            metadata_gz = sorted(glob.glob(os.path.join(self.folder, "*_metadata_file.csv.gz")))
            res = metadata_csv if metadata_csv else metadata_gz
            if len(res) == 1:
                self.read_metadata(res[0])
            elif len(res) > 1:
                self.read_metadata(res[0])
                print(
                    f"Found multiple metadata files matching expected patterns; "
                    f"using {os.path.basename(res[0])}"
                )
            else:
                self.metadata = None
                print(
                    f"No AnnData .h5ad or matching metadata file found at {self.folder}. "
                    "Expected *_metdata.csv or *_metadata_file.csv.gz"
                )
        else:
            if self.adata.obs.index.name == "cell_ID":
                self.adata.obs['cell_ID'] = self.adata.obs.index
                res = [re.search("c_.*_(.*)_(.*)", i) for i in self.adata.obs.index]
                self.adata.obs['UID'] = [pair(int(i.group(1)), int(i.group(2))) for i in res]
                self.adata.obs.set_index('UID', inplace=True)
                self.metadata = self.adata.obs
            elif self.adata.obs.index.name != "UID":
                print(f"Expected index cell_ID or UID in AnnData obs, not {self.adata.obs.index}, \
                    unable to read metadata")
            else:
                self.metadata = self.adata.obs

        # launch viewer
        if viewer:
            self.viewer = viewer
        else:
            self.viewer = napari.Viewer()
        self.viewer.layers.events.removed.connect(self._layer_removed)
        self.update_console()
        self.viewer.title = self.name
        self.viewer.scale_bar.visible = True
        self.viewer.scale_bar.unit = "mm"
        self._configure_z_plane_overlay()
        if any("labels" in x for x in self.grp.group_keys()):
            self.add_cell_labels()
            self.add_segmentation()
        self.add_fov_labels()

    def _load_targets(self):
        target_path = os.path.join(self.folder, "targets.hdf5")
        if not os.path.exists(target_path):
            print(f"No targets.hdf5 file found at {self.folder}")
            return None

        df = vaex.open(target_path)
        # global_z must be added to the underlying DataFrameLocal (df) before the
        # join, because after join() the result is a DataFrameJoin and add_column /
        # numpy array __setitem__ require a DataFrameLocal.
        if self.is_3d and 'z' in df.get_column_names():
            z_values = df.z.evaluate().astype(np.int64)
            global_z = z_values.copy()
            for k, v in self.z_slice_to_index.items():
                global_z[z_values == int(k)] = int(v)
            df.add_column('global_z', global_z)
        offsets = vaex.from_pandas(self.fov_offsets.loc[:, ['X_mm', 'Y_mm', 'FOV']])
        targets = df.join(offsets, left_on='fov', right_on='FOV')
        if self.dash:
            targets['global_x'] = targets.x + targets.Y_mm*(1/self.mm_per_px) - self._top_left_mm()[1]*(1/self.mm_per_px)
            targets['global_y'] = targets.y - targets.X_mm*(1/self.mm_per_px) - self._top_left_mm()[0]*(1/self.mm_per_px)
        else:
            targets['global_x'] = targets.x - targets.X_mm*(1/self.mm_per_px) - self._top_left_mm()[1]*(1/self.mm_per_px)
            targets['global_y'] = targets.y + targets.Y_mm*(1/self.mm_per_px) - self._top_left_mm()[0]*(1/self.mm_per_px)
        return targets

    def _load_protein_state(self):
        if 'protein' not in self.grp.group_keys():
            return

        self._proteins = [s.replace(r'.zarr', '') for s in self.grp['protein'].group_keys()]
        if 'CosMx' in self.grp['protein'].attrs:
            self.expr_scale = self.grp['protein'].attrs['CosMx']['scale']
            self.protein_z_step_um = float(self.grp['protein'].attrs['CosMx'].get('z_step_um', self.protein_z_step_um if self.protein_z_step_um is not None else DEFAULT_PROTEIN_Z_STEP_UM))
            self.protein_z_step_mm = self.protein_z_step_um / 1000
    
    def __repr__(self):
        return f'napari_cosmx.gemini.Gemini(r"{self.folder}")'

    def show_widget(self):
        if self.viewer is not None:
            widgets = [w for w in self.viewer.window._dock_widgets.values()
                if isinstance(w.widget(), GeminiQWidget)]
            if len(widgets) != 0:
                for w in widgets:
                    w.show()
            else:
                self.viewer.window.add_dock_widget(GeminiQWidget(self.viewer, self),
                    area='right',
                    name=self.name)

    def show_group_manager(self):
        """Open the Layer Group Manager as a dock widget in the viewer."""
        from napari_cosmx._dock_widget import LayerGroupManagerWidget
        widgets = [w for w in self.viewer.window._dock_widgets.values()
            if isinstance(w.widget(), LayerGroupManagerWidget)]
        if widgets:
            for w in widgets:
                w.show()
        else:
            self.viewer.window.add_dock_widget(
                LayerGroupManagerWidget(self.viewer),
                area='right',
                name='Layer Groups'
            )

    def _top_left_mm(self):
        """Return (y, x) tuple of mm translation for image layers
        """
        if self.dash:
            return (-max(self.fov_offsets['X_mm']), min(self.fov_offsets['Y_mm']))
        else:
            return (min(self.fov_offsets['Y_mm']), -max(self.fov_offsets['X_mm']))

    def _image_translate(self, layer_ndim=None):
        """Translation for an image layer.

        `layer_ndim` is the dimensionality of the array actually being added,
        which can be 2 even in a 3D dataset when morphology is a single plane.
        """
        ndim = self.ndim if layer_ndim is None else layer_ndim
        return self._top_left_mm() if ndim != 3 else (0.0, *self._top_left_mm())

    def _image_scale(self, protein=False, layer_ndim=None):
        """Scale for an image layer; see _image_translate for `layer_ndim`."""
        ndim = self.ndim if layer_ndim is None else layer_ndim
        if protein:
            xy_scale = self.mm_per_px*(1/self.expr_scale)
            if ndim == 3:
                return (self.protein_z_step_mm, xy_scale, xy_scale)
            return (xy_scale, xy_scale)
        if ndim == 3:
            return (self.rna_z_step_mm, self.mm_per_px, self.mm_per_px)
        return (self.mm_per_px, self.mm_per_px)

    def _fov_shape_scale(self):
        if self.is_3d:
            return (self.rna_z_step_mm, 1.0, 1.0)
        return (1.0, 1.0)

    def _resolve_single_z_index(self, z_plane=None):
        z_indices = self._resolve_z_indices(
            max_slices=len(self.z_slices),
            z_plane=z_plane,
            use_current_if_none=None,
        )
        if len(z_indices) != 1:
            raise ValueError("Expected a single z plane for this operation; pass one slice (e.g. 2 or 'current').")
        return int(z_indices[0])

    def _slice_multiscale_at_z(self, image_pyramid, z_index):
        sliced = []
        for image in image_pyramid:
            if image.ndim == 2:
                sliced.append(image)
            elif image.ndim == 3:
                sliced.append(image[z_index])
            elif image.ndim == 4:
                sliced.append(image[z_index, :, :, :])
            else:
                raise ValueError(f"Unsupported image ndim for z slicing: {image.ndim}")
        return sliced

    def _xy_scale(self, protein=False):
        if protein:
            xy_scale = self.mm_per_px * (1 / self.expr_scale)
            return (xy_scale, xy_scale)
        return (self.mm_per_px, self.mm_per_px)

    def z_plane_mapping(self):
        """Return mapping from array z index to source z plane value."""
        return {index: int(z_value) for index, z_value in enumerate(self.z_slices)}

    def z_plane_to_index(self, z_plane):
        """Map source z plane value to array z index."""
        z_plane = int(z_plane)
        if z_plane not in self.z_slice_to_index:
            raise ValueError(f"z plane {z_plane} not in available z planes: {self.z_slices}")
        return int(self.z_slice_to_index[z_plane])

    def _z_index_from_value_or_index(self, z_plane):
        z_plane = int(z_plane)
        if z_plane in self.z_slice_to_index:
            return int(self.z_slice_to_index[z_plane])
        if 0 <= z_plane < len(self.z_slices):
            return z_plane
        raise ValueError(
            f"z plane {z_plane} is neither a known z value nor a valid z index. "
            f"Available z planes: {self.z_slices}"
        )

    def _current_z_index(self):
        current_step = getattr(self.viewer.dims, 'current_step', ())
        z_index = int(current_step[0]) if len(current_step) > 0 else 0
        return max(0, min(z_index, len(self.z_slices) - 1))

    def _format_z_planes(self):
        if len(self.z_slices) <= 8:
            return ",".join(str(int(z)) for z in self.z_slices)
        head = ",".join(str(int(z)) for z in self.z_slices[:4])
        tail = ",".join(str(int(z)) for z in self.z_slices[-2:])
        return f"{head},...,{tail}"

    def _update_z_plane_overlay(self, event=None):
        if not self.is_3d:
            return
        overlay = getattr(self.viewer, 'text_overlay', None)
        if overlay is None:
            return
        z_index = self._current_z_index()
        z_value = int(self.z_slices[z_index])
        ndisplay = int(getattr(self.viewer.dims, 'ndisplay', 2))
        overlay.text = (
            f"Z display: {'2D' if ndisplay == 2 else '3D'} | "
            f"index {z_index}/{len(self.z_slices)-1} -> z={z_value} | "
            f"available: [{self._format_z_planes()}]"
        )

    def _configure_z_plane_overlay(self):
        if not self.is_3d:
            return
        overlay = getattr(self.viewer, 'text_overlay', None)
        if overlay is None:
            return
        overlay.visible = True
        # Keep z-readout away from the scale bar for readability.
        for attr in ('position', 'anchor'):
            if hasattr(overlay, attr):
                try:
                    setattr(overlay, attr, 'top_left')
                except Exception:
                    pass
        try:
            self.viewer.dims.events.current_step.connect(self._update_z_plane_overlay)
            self.viewer.dims.events.ndisplay.connect(self._update_z_plane_overlay)
        except Exception:
            pass
        self._update_z_plane_overlay()

    def _resolve_layer_name(self, base_name, modality, z_value=None):
        """Return a display name with modality tags only when needed.

        Disambiguate real collisions across loaded naming domains so napari does
        not fall back to auto-generated suffixes like ``[1]``.
        """
        name = base_name
        if modality == 'protein':
            if base_name in set(self.channels):
                name = f"{base_name} [protein]"
            elif self.has_rna and base_name in set(self.genes):
                name = f"{base_name} [protein]"
        elif modality == 'rna':
            if base_name in set(self.proteins):
                name = f"{base_name} [rna]"
        if z_value is not None:
            name = f"{name} [z={int(z_value)}]"
        return name

    def _tag_layer(self, layer, modality, channel=None, z_value=None):
        """Attach CosMx provenance metadata to a napari layer."""
        layer.metadata['cosmx_dataset'] = self.folder
        layer.metadata['cosmx_modality'] = modality
        if channel is not None:
            layer.metadata['cosmx_channel'] = channel
        if z_value is not None:
            layer.metadata['cosmx_z_plane'] = int(z_value)

    def omero(self, channel_name, protein=False, auto=True):
        """Return omero metadata
        https://ngff.openmicroscopy.org/latest/#omero-md
            - contrast limits for image as tuple
            - colormap name

        Args:
            channel_name (str): Name of channel (e.g. "DAPI" or "4-1BB")
            protein (bool): In protein group
            auto (bool): Calculate contrast limits from data if not in metadata
        """
        dirname = channel_name
        if protein:
            dirname = f"protein/{dirname}"
        if dirname not in self.grp:
            dirname = f"{dirname}.zarr"
        assert dirname in self.grp, f"{'protein' if protein else 'channel'} {channel_name} not found"

        window = {}
        color = DEFAULT_COLORMAPS[channel_name] if channel_name in DEFAULT_COLORMAPS else DEFAULT_COLORMAPS[None]
        metadata = self.grp[f"{dirname}"].attrs

        if "omero" in metadata:
            window = metadata["omero"]["channels"][0]["window"]
            if "color" in metadata["omero"]["channels"][0]:
                color = metadata["omero"]["channels"][0]["color"]
        elif auto: # calculate limits on the fly if no ome metadata
            print(f"Could not find omero metadata for {channel_name}, calculating contrast limits now.")
            datasets = metadata["multiscales"][0]["datasets"]
            component = datasets[-1]["path"] # use top-level of pyramid for speed
            image = da.from_zarr(os.path.join(self.folder, "images", f"{dirname}"), component=component)
            window = {}
            window['min'], window['max'] = int(da.min(image)), int(da.max(image))
            window['start'],window['end'] = [int(x) for x in da.percentile(image.ravel()[image.ravel()!=0], (0.1, 99.9))]
        else:
            window['min'], window['max'] = 0, 2**16 - 1
            window['start'],window['end'] = 0, 2**16 - 1

        if window['start'] - window['end'] == 0:
            if window['end'] == 0:
                print(f"WARNING: {channel_name} image is empty.")
                window['end'] = 1000
            else:
                window['start'] = 0
        window = _normalize_window(window)
        return {
            'window': window,
            'color': color
        }

    @qdebounced(timeout=400)
    def _update_omero_metadata(self, layer):
        self.grp[layer.metadata['path']].attrs['omero'] = layer.metadata['omero']
    
    def _layer_removed(self, event):
        layer = event.value
        print(f"{layer.name} removed")

    def _on_contrast_limits_changed(self, event):
        layer = event.source
        if 'omero' in layer.metadata:
            # update contrast limits
            window = {}
            window['min'], window['max'] = tuple(int(i) for i in layer.contrast_limits_range)
            window['start'], window['end'] = tuple(int(i) for i in layer.contrast_limits)
            layer.metadata['omero']['channels'][0]['window'] = window
            self._update_omero_metadata(layer)

    def _on_colormap_changed(self, event):
        layer = event.source
        if 'omero' in layer.metadata:
            # update colormap
            layer.metadata['omero']['channels'][0]['color'] = layer.colormap.name
            self._update_omero_metadata(layer)
    
    def add_channel(self, name, colormap=None, z_plane=None):
        """Add the requested morphology channel layer

        Args:
            name (str): Name of channel
            z_plane (int or str, optional): For 3D datasets, add a single z plane
                as a 2D layer. Supports explicit index, mapped z value, or
                'current'. If None, add the full 3D multiscale layer.
        """
        assert name in self.channels, f"{name} not one of available channels: {self.channels}"
        dirname = name if name in self.grp.group_keys() else f"{name}.zarr" 
        metadata = self.grp[f"{dirname}"].attrs
        # track updates to contrast limits and colormap
        sync_metadata = 'omero' in metadata
        try:
            metadata['path'] = dirname
        except PermissionError:
            print("Error updating metadata, local changes will not sync")
            sync_metadata = False
        datasets = metadata["multiscales"][0]["datasets"]
        im = [da.from_zarr(os.path.join(self.folder, "images", f"{dirname}"), component=d["path"]) for d in datasets]
        omero = self.omero(name)
        window = omero['window']
        if colormap is None:
            # use default for channel
            colormap = omero['color']
        elif 'omero' in metadata:
            # set as new default for channel
            metadata['omero']['channels'][0]['color'] = colormap
        z_value = None
        layer_name = name
        # Morphology may be stored as one 2D plane even in a 3D dataset (a 3D
        # resegmentation of a 2D-imaged slide), so place it by the array's own
        # dimensionality rather than the dataset's.
        stored_ndim = im[0].ndim
        layer_scale = self._image_scale(layer_ndim=stored_ndim)
        layer_translate = self._image_translate(layer_ndim=stored_ndim)
        if self.is_3d and stored_ndim == 3 and z_plane is not None:
            z_index = self._resolve_single_z_index(z_plane)
            im = self._slice_multiscale_at_z(im, z_index)
            z_value = int(self.z_slices[z_index])
            layer_name = f"{name} [z={z_value}]"
            layer_scale = self._xy_scale()
            layer_translate = self._top_left_mm()
        layer = self.viewer.add_image(im, colormap=colormap, blending="additive", name=layer_name,
            contrast_limits = (window['start'],window['end']),
            scale=layer_scale,
            translate=layer_translate,
            rotate=self.rotate,
            multiscale=True,
            rgb=False,
            metadata=metadata)
        layer.contrast_limits_range = window['min'],window['max']
        if sync_metadata:
            self._update_omero_metadata(layer)
            layer.events.contrast_limits.connect(self._on_contrast_limits_changed)
            layer.events.colormap.connect(self._on_colormap_changed)
        self._tag_layer(layer, 'morphology', name, z_value)
        self.order_layers()
        return layer

    def add_composite(self, z_plane=None):
        """Add RGB composite layer

        Args:
            z_plane (int or str, optional): For 3D datasets, add a single z plane
                as a 2D RGB layer. Supports explicit index, mapped z value, or
                'current'. If None, add the full multiscale layer.
        """
        assert self.grp, "No zarr images directory found"
        assert 'composite' in self.grp.group_keys(), "No composite layer found, use stitch-composite to create it"
        metadata = self.grp["composite"].attrs
        datasets = metadata["multiscales"][0]["datasets"]
        im = [da.from_zarr(os.path.join(self.folder, "images", "composite"), component=d["path"]) for d in datasets]
        z_value = None
        layer_name = "Composite"
        stored_ndim = im[0].ndim
        layer_scale = self._image_scale(layer_ndim=stored_ndim)
        layer_translate = self._image_translate(layer_ndim=stored_ndim)
        if self.is_3d and stored_ndim == 3 and z_plane is not None:
            z_index = self._resolve_single_z_index(z_plane)
            im = self._slice_multiscale_at_z(im, z_index)
            z_value = int(self.z_slices[z_index])
            layer_name = f"Composite [z={z_value}]"
            layer_scale = self._xy_scale()
            layer_translate = self._top_left_mm()
        layer = self.viewer.add_image(im, name=layer_name,
            scale=layer_scale,
            translate=layer_translate,
            rotate=self.rotate,
            rgb=True)
        self._tag_layer(layer, 'morphology', 'composite', z_value)
        return layer

    def add_protein(self, name, colormap=None, visible=True, z_plane=None):
        """Add the requested protein expression image

        Args:
            name (str): Name of protein
            z_plane (int or str, optional): For 3D datasets, add a single z plane
                as a 2D layer. Supports explicit index, mapped z value, or
                'current'. If None, add the full 3D multiscale layer.
        """        
        assert self.proteins, "No proteins found"
        assert name in self.proteins, f"{name} not found in proteins"
        dirname = name if name in self.grp[f"protein"].group_keys() else f"{name}.zarr"
        metadata = self.grp[f"protein/{dirname}"].attrs
        # track updates to contrast limits and colormap
        sync_metadata = 'omero' in metadata
        try:
            metadata['path'] = f"protein/{dirname}"
        except PermissionError:
            print("Error updating metadata, local changes will not sync")
            sync_metadata = False
        datasets = metadata["multiscales"][0]["datasets"]
        im = [da.from_zarr(os.path.join(self.folder, "images", "protein", dirname), component=d["path"]) for d in datasets]
        omero = self.omero(name, protein=True)
        window = omero['window']
        if colormap is None:
            # use default for channel
            colormap = omero['color']
        elif 'omero' in metadata:
            # set as new default for channel
            metadata['omero']['channels'][0]['color'] = colormap
        z_value = None
        layer_name = self._resolve_layer_name(name, 'protein')
        layer_scale = self._image_scale(protein=True)
        layer_translate = self._image_translate()
        if self.is_3d and z_plane is not None:
            z_index = self._resolve_single_z_index(z_plane)
            im = self._slice_multiscale_at_z(im, z_index)
            z_value = int(self.z_slices[z_index])
            layer_name = self._resolve_layer_name(name, 'protein', z_value)
            layer_scale = self._xy_scale(protein=True)
            layer_translate = self._top_left_mm()
        layer = self.viewer.add_image(im, name=layer_name, multiscale=True,
            colormap=colormap, blending="additive",
            contrast_limits = (window['start'],window['end']),
            scale=layer_scale,
            translate=layer_translate, visible=visible,
            rotate=self.rotate,
            rgb=False,
            metadata=metadata)
        layer.contrast_limits_range = window['min'],window['max']
        if sync_metadata:
            self._update_omero_metadata(layer)
            layer.events.contrast_limits.connect(self._on_contrast_limits_changed)
            layer.events.colormap.connect(self._on_colormap_changed)
        self._tag_layer(layer, 'protein', name, z_value)
        self.order_layers()
        return layer

    def add_segmentation(self, z_plane=None):
        """Add the cell segmentation image layer.

        Args:
            z_plane (int or str, optional): For 3D datasets, add a single z plane
                as a 2D layer. Supports explicit index, mapped z value, or
                'current'. If None, add the full 3D multiscale layer.
        """
        assert any("labels" in x for x in self.grp.group_keys()),\
            f"labels not found in zarr keys: {self.grp.group_keys()}"
        labelsdir = 'labels' if 'labels' in self.grp.group_keys() else 'labels.zarr'
        datasets = self.grp[labelsdir].attrs["multiscales"][0]["datasets"]
        labels = []
        for d in datasets:
            arr = da.from_zarr(os.path.join(self.folder, "images", labelsdir), component=d["path"])
            if self.is_3d and z_plane is not None:
                z_index = self._resolve_single_z_index(z_plane)
                arr = arr[z_index]
            labels.append(arr.map_blocks(_find_boundaries_by_slice))
        cm = Colormap(['transparent', 'cyan'], controls = [0.0, 1.0])
        z_value = None
        layer_name = 'Segmentation'
        layer_scale = self._image_scale()
        layer_translate = self._image_translate()
        if self.is_3d and z_plane is not None:
            z_index = self._resolve_single_z_index(z_plane)
            z_value = int(self.z_slices[z_index])
            layer_name = f"Segmentation [z={z_value}]"
            layer_scale = self._xy_scale()
            layer_translate = self._top_left_mm()
        layer = self.viewer.add_image(labels, contrast_limits=(0, 1), colormap=cm,
            scale=layer_scale, translate=layer_translate, blending="translucent",
            rotate=self.rotate,
            rgb=False)
        layer.opacity = 0.75
        self.segmentation_layer = layer
        layer.name = layer_name
        self._tag_layer(layer, 'segmentation', z_value=z_value)
        return layer
        
    def add_cell_labels(self, z_plane=None):
        """Add the cell labels layer.

        Args:
            z_plane (int or str, optional): For 3D datasets, add a single z plane
                as a 2D labels layer. Supports explicit index, mapped z value, or
                'current'. If None, add the full 3D multiscale labels layer.
        """
        assert any("labels" in x for x in self.grp.group_keys()),\
            f"labels not found in zarr keys: {self.grp.group_keys()}"
        labelsdir = 'labels' if 'labels' in self.grp.group_keys() else 'labels.zarr'
        datasets = self.grp[labelsdir].attrs["multiscales"][0]["datasets"]
        labels = []
        for d in datasets:
            arr = da.from_zarr(os.path.join(self.folder, "images", labelsdir), component=d["path"])
            if self.is_3d and z_plane is not None:
                z_index = self._resolve_single_z_index(z_plane)
                arr = arr[z_index]
            labels.append(arr)
        # TODO: need to scale for protein expression images
        z_value = None
        layer_scale = self._image_scale()
        layer_translate = self._image_translate()
        layer_name = 'Metadata'
        if self.is_3d and z_plane is not None:
            z_index = self._resolve_single_z_index(z_plane)
            z_value = int(self.z_slices[z_index])
            layer_name = f"Metadata [z={z_value}]"
            layer_scale = self._xy_scale()
            layer_translate = self._top_left_mm()
        layer = self.viewer.add_labels(labels, scale=layer_scale, translate=layer_translate,
            rotate=self.rotate)
        if z_plane is None:
            self.cells_layer = layer
            self.cells_layer.opacity = 1.0
            self.cells_layer.visible = False
        if self.metadata is not None:
            df = self.metadata.copy()
            df['index'] = self.metadata.index
            layer.features = pd.concat([pd.DataFrame.from_dict({'index': [0]}), df])
        layer.name = layer_name
        layer.metadata['label_info'] = {
            'contour': 0,
            'label_color_index': {None: 0}
        }
        layer.editable = False
        self._tag_layer(layer, 'segmentation', z_value=z_value)
        return layer

    def _resolve_labels_layer(self, labels_layer=None):
        if labels_layer is None:
            assert self.cells_layer is not None, "No labels layer selected"
            return self.cells_layer
        if isinstance(labels_layer, str):
            assert labels_layer in self.viewer.layers, f"Layer {labels_layer} not found"
            layer = self.viewer.layers[labels_layer]
        else:
            layer = labels_layer
        assert isinstance(layer, napari.layers.Labels), "labels_layer must be a Labels layer"
        return layer

    def _labels_for_shape_selection(self):
        labels_data = self.cells_layer.data
        if hasattr(labels_data, 'shapes'):
            level = 1 if len(labels_data) > 1 else 0
            labels_arr = labels_data[level]
            dims = np.asarray(labels_data.shapes, dtype=float)
            spatial_scale = dims[0, -2:] / dims[level, -2:]
            return labels_arr, spatial_scale
        return labels_data, np.ones(2, dtype=float)

    def _resolve_z_indices(self, max_slices, z_plane=None, default_index=None, use_current_if_none=None):
        if z_plane == 'current':
            z_plane = None
            default_index = None
            use_current_if_none = True
        elif z_plane == 'auto':
            z_plane = None
            use_current_if_none = None
        if z_plane in {'any', 'all'}:
            return list(range(max_slices))
        if isinstance(z_plane, range):
            z_indices = [self._z_index_from_value_or_index(z_value) for z_value in z_plane]
            if not z_indices:
                raise ValueError("z_plane range must contain at least one slice")
            if min(z_indices) < 0 or max(z_indices) >= max_slices:
                raise IndexError(f"z plane range {list(z_plane)} is out of bounds for labels with {max_slices} slices")
            return sorted(set(z_indices))
        if isinstance(z_plane, (tuple, list)):
            if len(z_plane) == 2:
                start, stop = [int(z_value) for z_value in z_plane]
                if stop < start:
                    raise ValueError("z_plane range stop must be greater than or equal to start")
                z_indices = [
                    index for index, z_value in enumerate(self.z_slices)
                    if start <= int(z_value) <= stop
                ]
                if not z_indices:
                    raise ValueError(
                        f"No available z planes in requested range [{start}, {stop}]. "
                        f"Available z planes: {self.z_slices}"
                    )
            else:
                z_indices = [self._z_index_from_value_or_index(z_value) for z_value in z_plane]
            if not z_indices:
                raise ValueError("z_plane list must contain at least one slice")
            if min(z_indices) < 0 or max(z_indices) >= max_slices:
                raise IndexError(f"z plane selection {z_plane} is out of bounds for labels with {max_slices} slices")
            return sorted(set(z_indices))
        if z_plane is None:
            if default_index is not None:
                z_index = int(default_index)
            else:
                if use_current_if_none is None:
                    ndisplay = getattr(self.viewer.dims, 'ndisplay', 2)
                    use_current_if_none = ndisplay == 2
                if use_current_if_none:
                    current_step = getattr(self.viewer.dims, 'current_step', ())
                    z_index = int(current_step[0]) if len(current_step) > 0 else 0
                else:
                    return list(range(max_slices))
        else:
            z_value = int(z_plane)
            z_index = self._z_index_from_value_or_index(z_value)
        if z_index < 0 or z_index >= max_slices:
            raise IndexError(f"z plane {z_plane} is out of bounds for labels with {max_slices} slices")
        return [z_index]

    def _format_rna_z_selection_suffix(self, z_indices):
        """Format a concise RNA z-selection suffix for layer names."""
        if not self.is_3d:
            return None
        z_indices = sorted({int(i) for i in z_indices})
        if len(z_indices) == 0:
            return "3D no-z"
        if z_indices == list(range(len(self.z_slices))):
            return "3D all-z"
        z_values = [int(self.z_slices[index]) for index in z_indices]
        if len(z_values) == 1:
            return f"3D z={z_values[0]}"
        if all((z_values[i + 1] - z_values[i]) == 1 for i in range(len(z_values) - 1)):
            return f"3D z={z_values[0]}-{z_values[-1]}"
        return "3D z=" + ",".join(str(value) for value in z_values)

    def _shape_z_indices(self, labels_arr, z_plane=None, inferred_z_index=None):
        if not self.is_3d:
            return [None]
        return self._resolve_z_indices(
            max_slices=labels_arr.shape[0],
            z_plane=z_plane,
            default_index=inferred_z_index,
            use_current_if_none=None,
        )

    def _mask_pixels(self, mask_shape, shape_type, crop_shape):
        if shape_type == 'ellipse':
            r_radius = (mask_shape[-1, 0] - mask_shape[0, 0]) / 2
            r = mask_shape[0, 0] + r_radius
            c_radius = (mask_shape[1, 1] - mask_shape[0, 1]) / 2
            c = mask_shape[0, 1] + c_radius
            return ellipse(r, c, r_radius, c_radius, crop_shape)
        r, c = mask_shape.T
        return polygon(r, c, crop_shape)
  
    def cells_in_shape(self, shape_layer, idx=None, z_plane=None):
        """Find unique cell IDs from closed shape(s)

        Cells returned are from segmentation, some may be missing in metadata.

        Args:
            shape_layer (Shapes layer): Shapes layer drawn in napari
            idx (list or int): Which shapes from layer to use, if None use all
            z_plane (int, str, tuple, list, or range, optional): For 3D datasets,
                use a single z plane, 'any', 'current', 'auto', or a z-plane range.
                Two-element tuple or list inputs are treated as an inclusive range.
                If None or 'auto', use current z in 2D display mode and all z slices
                in 3D display mode, unless the shape stores a z index.

        Returns:
            array: cell UIDs
        """
        cells = set()
        assert isinstance(shape_layer, napari.layers.Shapes), "Not a shapes layer"
        if idx is None:
            idx = range(len(shape_layer.data))
        elif not isinstance(idx, list):
            idx = [idx]
        if self.cells_layer is not None:
            labels_arr, spatial_scale = self._labels_for_shape_selection()
            for i in idx:
                if shape_layer.shape_type[i] not in ['ellipse', 'polygon', 'rectangle']:
                    continue
                shape_data = np.asarray(shape_layer.data[i], dtype=float)
                inferred_z_index = None
                if self.is_3d and shape_data.shape[1] == 3:
                    unique_z = np.unique(np.rint(shape_data[:, 0]).astype(int))
                    if len(unique_z) == 1:
                        inferred_z_index = int(unique_z[0])
                    shape_data = shape_data[:, 1:]
                elif shape_data.shape[1] != 2:
                    raise ValueError(f"Unsupported shape dimensionality: {shape_data.shape[1]}")
                mask_shape = (shape_data * self.mm_per_px - self._top_left_mm()) * (1/self.mm_per_px)
                mask_shape = mask_shape / spatial_scale
                mask_shape = mask_shape.astype(int)
                z_indices = self._shape_z_indices(labels_arr, z_plane=z_plane, inferred_z_index=inferred_z_index)
                for z_index in z_indices:
                    labels_plane = labels_arr if z_index is None else labels_arr[z_index]
                    ymin = max(np.min(mask_shape[:, 0]), 0)
                    ymax = min(np.max(mask_shape[:, 0]), labels_plane.shape[0] - 1)
                    xmin = max(np.min(mask_shape[:, 1]), 0)
                    xmax = min(np.max(mask_shape[:, 1]), labels_plane.shape[1] - 1)
                    crop = labels_plane[ymin:ymax, xmin:xmax]
                    crop_shape = mask_shape.copy()
                    crop_shape[:, 0] -= ymin
                    crop_shape[:, 1] -= xmin
                    rr, cc = self._mask_pixels(crop_shape, shape_layer.shape_type[i], crop.shape)
                    uids = np.unique(crop.compute()[rr, cc])
                    cells.update(uids[uids != 0])
        return np.array([int(i) for i in cells], dtype='uint32')

    def cell_target_count(self, cells):
        """Return frequency table for targets in cells

        For the given cells create frequency table
        for the targets assigned to those cells.

        Args:
            array: cell UIDs

        Returns:
            pandas.Series: series with counts and index of target
        """
        cells_df = pd.DataFrame([unpair(i) for i in cells],
                                columns=['fov', 'CellId'])
        self.targets['join_id'] = self.targets['fov'].astype(str) + '_' + self.targets['CellId'].astype(str)
        cells_df['join_id'] = cells_df['fov'].astype(str) + '_' + cells_df['CellId'].astype(str)
        df = self.targets.join(vaex.from_pandas(cells_df), on='join_id', how='inner', lsuffix='left')
        del self.targets['join_id']
        return df['target'].to_pandas_series().value_counts()

    def cells_at_points(self, points_layer, z_plane=None):
        """Find unique cell IDs from points

        Cells returned are from segmentation, some may be missing in metadata.

        Args:
            shape_layer (Points layer): Points layer drawn in napari
            z_plane (int, str, tuple, list, or range, optional): For 3D datasets,
                use a single z plane, 'any', 'current', 'auto', or a z-plane range.
                If None or 'auto', use current z in 2D display mode and all z slices
                in 3D display mode.

        Returns:
            array: cell UIDs
        """
        assert isinstance(points_layer, napari.layers.Points), "Not a points layer"
        assert self.cells_layer is not None, "No labels layers"
        layer_coords = points_layer.data
        scale = np.array(self._image_scale())
        translate = np.array(self._image_translate())
        layer_image_coords = (layer_coords*scale - translate) / scale
        layer_image_coords = layer_image_coords.astype('int')
        if self.is_3d:
            if layer_image_coords.shape[1] == 2:
                current_step = getattr(self.viewer.dims, 'current_step', ())
                z_coords = np.full(layer_image_coords.shape[0], int(current_step[0]) if len(current_step) > 0 else 0, dtype=int)
                y_coords = layer_image_coords[:, 0]
                x_coords = layer_image_coords[:, 1]
            else:
                z_coords = layer_image_coords[:, 0]
                y_coords = layer_image_coords[:, 1]
                x_coords = layer_image_coords[:, 2]
            selected_z = self._resolve_z_indices(
                max_slices=self.cells_layer.data[0].shape[0],
                z_plane=z_plane,
                use_current_if_none=None,
            )
            keep = np.isin(z_coords, selected_z)
            cells = self.cells_layer.data[0].vindex[z_coords[keep], y_coords[keep], x_coords[keep]].compute()
        else:
            cells = self.cells_layer.data[0].vindex[layer_image_coords[:, 0], layer_image_coords[:, 1]].compute()
        return cells[cells != 0]

    def get_offsets(self, fov):
        """Get offsets for given FOV

        Args:
            fov (int): FOV number

        Returns:
            tuple: x and y offsets in mm
        """
        offset = self.fov_offsets[self.fov_offsets['FOV'] == fov]
        if self.dash:
            x_offset = -offset.iloc[0, ]["X_mm"]
            y_offset = offset.iloc[0, ]["Y_mm"]
        else:
            x_offset = offset.iloc[0, ]["Y_mm"]
            y_offset = -offset.iloc[0, ]["X_mm"]
        return (x_offset, y_offset)

    @property
    def channels(self):
        """List: Available morphology channels to display."""
        return self._channels

    @property
    def has_rna(self):
        """bool: Transcript modality is available for this dataset."""
        return self.targets is not None

    @property
    def has_protein(self):
        """bool: Protein modality is available for this dataset."""
        return len(self._proteins) > 0

    @property
    def modalities(self):
        """tuple: Available omics modalities for this dataset."""
        modalities = []
        if self.has_rna:
            modalities.append('rna')
        if self.has_protein:
            modalities.append('protein')
        return tuple(modalities)

    @property
    def is_multiomics(self):
        """bool: Both RNA and protein modalities are available."""
        return self.has_rna and self.has_protein

    @property
    def is_protein(self):
        """bool: Backward-compatible alias for datasets that include protein data."""
        return self.has_protein

    @property
    def proteins(self):
        """List: Available protein channels to display."""
        return self._proteins

    @property
    def rotate(self):
        """float: Degrees to rotate all layers CCW for display."""
        return self._rotate

    @rotate.setter
    def rotate(self, angle):
        for i in self.viewer.layers:
            i.rotate = angle
        self._rotate = angle


    @property
    def genes(self):
        """List: If transcripts are loaded, return gene names."""
        if self.has_rna:
            return sorted([i for i in self.targets.category_labels('target')
                if not i.startswith('FalseCode') and
                not i.startswith('NegPrb') and
                not i.startswith('SystemControl')])
        return []

    def update_console(self):
        if self.viewer.window._qt_viewer.dockConsole.isVisible():
            self.viewer.update_console({'gem': self})
        else:
            self.viewer.window._qt_viewer.dockConsole.visibilityChanged.connect(
                lambda: self.viewer.update_console({'gem': self})
            )

    def rect_for_fov(self, fov):
        fov_height = self.fov_height*self.mm_per_px
        fov_width = self.fov_width*self.mm_per_px
        rect = np.array([
            list(self.get_offsets(fov)),
            list(map(sum, zip(self.get_offsets(fov), (fov_height, 0)))),
            list(map(sum, zip(self.get_offsets(fov), (fov_height, fov_width)))),
            list(map(sum, zip(self.get_offsets(fov), (0, fov_width))))
        ])
        y_offset = self._top_left_mm()[0]
        x_offset = self._top_left_mm()[1]
        return [[i[0] - y_offset, i[1] - x_offset] for i in rect]

    def rect_for_fov_z(self, fov, z_index):
        rect = self.rect_for_fov(fov)
        return [[z_index, y, x] for y, x in rect]

    # keep reference to layer and keep on top
    # keep selected FOVs in object
    def add_fov_labels(self):
        if self.is_3d:
            rects = [
                self.rect_for_fov_z(fov, z_index)
                for z_index in range(len(self.z_slices))
                for fov in self.fov_offsets['FOV']
            ]
            labels = np.array([
                int(fov)
                for z_index in range(len(self.z_slices))
                for fov in self.fov_offsets['FOV']
            ])
            z_planes = np.array([
                int(self.z_slices[z_index])
                for z_index in range(len(self.z_slices))
                for fov in self.fov_offsets['FOV']
            ])
            shape_properties = {
                'label': labels,
                'z_plane': z_planes,
            }
        else:
            rects = [self.rect_for_fov(i) for i in self.fov_offsets['FOV']]
            shape_properties = {
                'label': self.fov_offsets['FOV'].to_numpy()
            }
        text_parameters = {
            'text': 'label',
            'size': 12,
            'color': 'white'
        }
        shapes_layer = self.viewer.add_shapes(rects,
            face_color='#90ee90',
            edge_color='white',
            edge_width=0.02,
            properties=shape_properties,
            text = text_parameters,
            name = 'FOV labels',
            scale=self._fov_shape_scale(),
            translate=self._image_translate(),
            rotate=self.rotate)
        shapes_layer.opacity = 0.5
        shapes_layer.editable = False
        self.fov_labels = shapes_layer

    def center_fov(self, fov:int, buffer:float=1.0):
        """Center FOV in canvas and zoom to fill

        Args:
            fov (int): FOV number
            buffer (float): Buffer size for zoom. < 1 equals zoom out.
        """        
        extent = [np.min(self.rect_for_fov(fov), axis=0) + self._top_left_mm(),
            np.max(self.rect_for_fov(fov), axis=0) + self._top_left_mm()]
        size = extent[1] - extent[0]
        self.viewer.camera.center = np.add(extent[0], np.divide(size, 2))
        self.viewer.camera.zoom = np.min(np.array(self.viewer._canvas_size) / size) * buffer
        
    def move_to_top(self, layer):
        """Move layer to top of layer list. Layer must be in list.

        Args:
            layer (Layer): Layer to move
        """        
        idx = self.viewer.layers.index(layer)
        self.viewer.layers.move(idx, -1)

    def order_layers(self):
        if self.cells_layer is not None:
            self.move_to_top(self.cells_layer)
        if self.color_cells_layer is not None:
            self.move_to_top(self.color_cells_layer)
        points_layers = [l for l in self.viewer.layers if isinstance(l, napari.layers.Points)]
        for l in points_layers:
            self.move_to_top(l) 
        if self.segmentation_layer is not None:
            self.move_to_top(self.segmentation_layer)
        if self.fov_labels is not None:
            self.move_to_top(self.fov_labels)
        if self.cells_layer is not None:
            self.viewer.layers.selection.active = self.cells_layer

    def read_metadata(self, path):
        df = pd.read_csv(path)
        if not 'UID' in df:
            uid_built = False

            # Prefer parsing from string cell labels when available, e.g.
            # "c_1_130_1010" -> fov=130, cell=1010.
            for cell_col in ('cell', 'cell_ID', 'cell_id'):
                if cell_col not in df.columns:
                    continue
                parsed = df[cell_col].astype(str).str.extract(r"c_(?:.*_)?(\d+)_(\d+)$")
                if parsed.notna().all(axis=1).all():
                    fov_vals = parsed[0].astype(np.int64).to_numpy()
                    cell_vals = parsed[1].astype(np.int64).to_numpy()
                    df['UID'] = [pair(int(fov), int(cell)) for fov, cell in zip(fov_vals, cell_vals)]
                    uid_built = True
                    break

            # Fallback for tables where cell_ID are numeric and fov is explicit.
            if not uid_built and {'fov', 'cell_ID'}.issubset(df.columns):
                fov_vals = pd.to_numeric(df['fov'], errors='coerce')
                cell_vals = pd.to_numeric(df['cell_ID'], errors='coerce')
                if fov_vals.notna().all() and cell_vals.notna().all():
                    fov_vals = fov_vals.astype(np.int64).to_numpy()
                    cell_vals = cell_vals.astype(np.int64).to_numpy()
                    df['UID'] = [pair(int(fov), int(cell)) for fov, cell in zip(fov_vals, cell_vals)]
                    uid_built = True

            assert uid_built, (
                "Need UID column or parseable cell identifiers. Supported inputs are "
                "(1) string labels in 'cell'/'cell_ID' like c_*_<fov>_<cell>, or "
                "(2) numeric 'fov' + numeric 'cell_ID' columns."
            )
        df.set_index('UID', inplace=True)   
        self.metadata = df
        if self.cells_layer is not None:
            # update features
            df = self.metadata.copy()
            df['index'] = self.metadata.index
            self.cells_layer.features = pd.concat([pd.DataFrame.from_dict({'index': [0]}), df])
    
    def export_omero(self, path=None):
        """Export image metadata including colormap and contrast limits for all available channels.

        Args:
            path (str, optional): If provided write results to CSV. Defaults to None.

        Returns:
            DataFrame: Image metadata in pandas DataFrame.
        """
        omero = [self.omero(i, auto=False) for i in self.channels]
        df = pd.DataFrame(data={
            'name': self.channels,
            'type': 'channel',
            'color': [i['color'] for i in omero],
            'start': [i['window']['start'] for i in omero],
            'end': [i['window']['end'] for i in omero]
        })
        if self.proteins:
            omero = [self.omero(i, protein=True, auto=False) for i in self.proteins]
            protein_df = pd.DataFrame(data={
                'name': self.proteins,
                'type': 'protein',
                'color': [i['color'] for i in omero],
                'start': [i['window']['start'] for i in omero],
                'end': [i['window']['end'] for i in omero]
            })
            df = pd.concat([df, protein_df])
        if path:
            df.to_csv(path, index=False)
        return df
    
    def is_categorical_metadata(self, col_name):
        return not is_numeric_dtype(self.metadata[col_name]) or len(pd.unique(self.metadata[col_name])) < 30

    def color_cells(self, col_name, color=None, contour=0, subset=None, labels_layer=None, keep_previous=False):
        """Change cell labels layer based on metadata

        Args:
            col_name (str): Column name in metadata, "all" colors cells the same color.
            color (str|dict): (1) Color name if col_name is "all", or
                (2) dictionary with keys being metadata values and value being color name, or
                (3) colormap for continuous metadata, https://matplotlib.org/stable/tutorials/colors/colormaps.html
            contour (int): Labels layer contour, 0 is filled, otherwise thickness of lines.
            subset (list of int): List of cell UIDs, if given only color these cells.
            labels_layer (str or Labels layer, optional): Labels layer to color.
                If None, use the default Metadata labels layer.
            keep_previous (bool, optional): If True, keep previous colored layers instead of replacing.
                Defaults to False.
        """
        base_labels_layer = self._resolve_labels_layer(labels_layer)
        label_colors = {}
        if self.metadata is None:        
            assert col_name == "all", "No metadata loaded"
            color = transform_color(color)[0]
            label_colors = {1: color, None: transform_color('transparent')[0]}
        else:
            cells = subset if (subset is not None) else self.metadata.index
            if len(cells) == 0:
                return
            if col_name == "all":
                color = transform_color(color)[0]
                label_colors = {k:color for k in cells}
            else:
                assert col_name in self.metadata, f"{col_name} not in metadata"
                if self.is_categorical_metadata(col_name):
                    if color is None:
                        if self.adata is not None and col_name + "_colors" in self.adata.uns:
                            # get colors from AnnData object
                            categories = self.adata.obs[col_name].cat.categories
                            cat_colors = [transform_color(i)[0] for i in self.adata.uns[col_name + "_colors"]]
                            color = dict(zip(categories, cat_colors))
                        else:
                            # Prefer a per-column '<col>_color' (multi-annotation
                            # _metadata.csv), then legacy 'hex_color', for consistent
                            # colors across slides; otherwise auto-generate.
                            hex_map = categorical_color_map(self.metadata, col_name)
                            if hex_map is not None:
                                color = {k: transform_color(v)[0] for k, v in hex_map.items()}
                            else:
                                vals = np.unique(self.metadata[col_name])
                                cm = label_colormap(len(vals)+1)
                                color = dict(zip(vals, cm.colors[1:]))
                    else:
                        assert isinstance(color, dict), "color needs to be dict for categorical metadata"
                        color = {k:transform_color(v)[0] for k,v in color.items()}
                    # NaN != NaN in Python, so dict lookup with `v in color` fails
                    # for NaN keys. Use pd.isna to handle NaN values explicitly.
                    nan_color = next((c for k, c in color.items() if pd.isna(k)), None)
                    def _lookup(v):
                        if pd.isna(v):
                            return nan_color if nan_color is not None else transform_color('transparent')[0]
                        return color[v] if v in color else transform_color('transparent')[0]
                    label_colors = {k: _lookup(v)
                        for k, v in zip(cells, self.metadata.loc[cells][col_name])}
                else:
                    if color is None:
                        color = 'gray'
                    assert color in AVAILABLE_COLORMAPS, f"{color} not in {AVAILABLE_COLORMAPS.keys()}"
                    cm = AVAILABLE_COLORMAPS[color]
                    min_max_scaler = preprocessing.MinMaxScaler()
                    # normalize to 0-1 with full range present in metadata
                    x = pd.Series(
                        min_max_scaler.fit_transform(self.metadata[col_name].values.reshape(-1, 1))[:, 0],
                        self.metadata.index,
                        dtype=np.float32
                    )
                    label_colors = {k:v for k,v in zip(cells, cm.map(x.loc[cells]))}
            label_colors[None] = transform_color('transparent')[0]
        name = 'Cells' if col_name == 'all' else col_name
        # when labels are fast in napari can set labels layer colors directly
        self._color_cells(colors=label_colors, name=name, contour=contour, labels_layer=base_labels_layer, keep_previous=keep_previous)
        self.order_layers()

    def _map_labels_to_colors(self, arr):
        label_info = getattr(self, '_active_label_info', None)
        if label_info is None:
            label_info = self.cells_layer.metadata['label_info']
        contour = label_info['contour']
        label_color_index = label_info['label_color_index']
        # following logic in napari Labels layer
        arr_modified = arr
        if contour > 0:
            arr_modified = np.zeros_like(arr)
            struct_elem = ndi.generate_binary_structure(arr.ndim, 1)
            thick_struct_elem = ndi.iterate_structure(
                struct_elem, contour
            ).astype(bool)
            boundaries = ndi.grey_dilation(
                arr, footprint=struct_elem
            ) != ndi.grey_erosion(arr, footprint=thick_struct_elem)
            arr_modified[boundaries] = arr[boundaries]
        u, inv = np.unique(arr_modified, return_inverse=True)
        image = np.array(
            [
                label_color_index[x]
                if x in label_color_index
                else label_color_index[None]
                for x in u
            ]
        )[inv].reshape(arr_modified.shape)
        return image

    def _split_label_levels(self, labels_data):
        """Normalize labels data into per-level arrays and multiscale flag.

        napari may expose labels data as a plain list, a MultiScaleData
        container, a dask array, or a numpy array depending on version and
        data source. This helper provides a uniform `(levels, is_multiscale)`
        representation for downstream color mapping.
        """
        if isinstance(labels_data, list):
            return labels_data, True
        if hasattr(labels_data, 'shapes'):
            return list(labels_data), True
        return [labels_data], False

    def _map_label_level_to_colors(self, arr):
        """Map one labels level to RGBA colors for dask or numpy inputs."""
        if hasattr(arr, 'map_blocks'):
            return arr.map_blocks(self._map_labels_to_colors, dtype=float)
        return self._map_labels_to_colors(np.asarray(arr)).astype(float, copy=False)

    def _normalize_label_colors(self, colors):
        """Normalize label-color dict values to RGBA arrays.

        Older napari releases can fail when `color_dict_to_colormap` receives
        string values directly in a label-color dict. Converting values to RGBA
        here keeps behavior consistent across napari versions.
        """
        normalized = {}
        for label, value in colors.items():
            rgba = transform_color('transparent')[0] if value is None else transform_color(value)[0]
            normalized[label] = rgba
        if None not in normalized:
            normalized[None] = transform_color('transparent')[0]
        return normalized

    def _color_cells(self, colors, name='Cells', contour=0, labels_layer=None, keep_previous=False):
        if not keep_previous and self.color_cells_layer in self.viewer.layers:
            self.viewer.layers.remove(self.color_cells_layer)
        base_labels_layer = self._resolve_labels_layer(labels_layer)
        color_layer_name = name
        z_value = base_labels_layer.metadata.get('cosmx_z_plane')
        if z_value is None:
            z_match = re.search(r'\[z=\s*(-?\d+)\s*\]', getattr(base_labels_layer, 'name', ''))
            if z_match is not None:
                z_value = int(z_match.group(1))
        if z_value is not None:
            color_layer_name = f"{name} [z={int(z_value)}]"
        labels_data, is_multiscale = self._split_label_levels(base_labels_layer.data)
        normalized_colors = self._normalize_label_colors(colors)
        # from napari Labels layer
        with warnings.catch_warnings():
            # suppress known warning
            warnings.filterwarnings(action='ignore', message='.*distinct colors.*')
            custom_colormap, label_color_index = color_dict_to_colormap(normalized_colors)
        if 'label_info' not in base_labels_layer.metadata:
            base_labels_layer.metadata['label_info'] = {'contour': 0, 'label_color_index': {None: 0}}
        base_labels_layer.metadata['label_info']['contour'] = contour
        base_labels_layer.metadata['label_info']['label_color_index'] = label_color_index
        self._active_label_info = base_labels_layer.metadata['label_info']
        im_levels = [self._map_label_level_to_colors(arr) for arr in labels_data]
        self._active_label_info = None
        add_image_sig = signature(self.viewer.add_image)
        image_kwargs = {
            'name': color_layer_name,
            'scale': tuple(base_labels_layer.scale),
            'translate': tuple(base_labels_layer.translate),
            'blending': 'translucent',
            'rotate': self.rotate,
            'contrast_limits': (0, 1),
            'colormap': custom_colormap,
            'cache': True,
            'rgb': False,
        }
        if 'interpolation' in add_image_sig.parameters:
            image_kwargs['interpolation'] = 'nearest'
        image_data = im_levels if is_multiscale else im_levels[0]
        layer = self.viewer.add_image(image_data, **image_kwargs)
        self.color_cells_layer = layer

    def plot_transcripts(self, gene, color, point_size=5, z_plane=None):
        """Plot targets as dots

        Args:
            gene (str): Target to plot
            color (str): Color for points.
            point_size (int, optional): Point size. Defaults to 5.
            z_plane (int, str, tuple, list, or range, optional): For 3D datasets,
                use a single z plane, 'all', 'current', 'auto', or a z-plane range.
                ('any' is also accepted as an alias of 'all').
                If None or 'auto', plot all z slices.

        Returns:
            napari.layers.Points: A points layer for transcripts.
        """
        assert self.targets, "No targets found, use read_targets.py first to create targets.hdf5 file."
        targets = self.targets[self.targets.target == gene]
        selected_z = None
        if self.is_3d and 'global_z' in targets.get_column_names():
            selected_z = self._resolve_z_indices(
                max_slices=len(self.z_slices),
                z_plane=z_plane,
                use_current_if_none=False,
            )
            targets = targets[targets.global_z.isin(selected_z)]
        y = targets.global_y.evaluate()
        x = targets.global_x.evaluate()
        if self.is_3d and 'global_z' in targets.get_column_names():
            z = targets.global_z.evaluate()
            points = np.column_stack((z, y, x))
        else:
            points = np.column_stack((y, x))
        add_points_sig = signature(self.viewer.add_points)
        edge_kw = {}
        if 'edge_width' in add_points_sig.parameters:
            edge_kw['edge_width'] = 0
        elif 'border_width' in add_points_sig.parameters:
            edge_kw['border_width'] = 0
        points_layer = self.viewer.add_points(points,
            size = point_size,
            face_color=color,
            scale=self._image_scale(),
            translate=self._image_translate(),
            rotate=self.rotate,
            **edge_kw)
        points_layer_name = self._resolve_layer_name(gene, 'rna')
        if selected_z is not None:
            points_layer_name = f"{points_layer_name} [{self._format_rna_z_selection_suffix(selected_z)}]"
        points_layer.name = points_layer_name
        points_layer.opacity = 1.0
        selected_z_values = None
        if selected_z is not None:
            selected_z_values = [int(self.z_slices[index]) for index in sorted(selected_z)]
        single_z_value = selected_z_values[0] if selected_z_values is not None and len(selected_z_values) == 1 else None
        self._tag_layer(points_layer, 'rna', gene, single_z_value)
        if selected_z_values is not None:
            points_layer.metadata['cosmx_z_selection'] = selected_z_values
        self.order_layers()
        return points_layer

    def add_points(self, fov=None, color=None, gray_no_cell=None, point_size=1, z_plane=None):
        """Add all targets

        Args:
            fov (int or list, optional): Only add points for specified fov(s).
            color (str, optional): Color for points. If None will color by gene.
            gray_no_cell (bool, optional): Whether to color targets with no CellID gray. Defaults to True if color specified.
            point_size (int, optional): Point size. Defaults to 1.
            z_plane (int, str, tuple, list, or range, optional): For 3D datasets,
                use a single z plane, 'all', 'current', 'auto', or a z-plane range.
                ('any' is also accepted as an alias of 'all').
                If None or 'auto', plot all z slices.

        Returns:
            napari.layers.Points: A points layer for transcripts.
        """
        assert self.targets, "No targets found, use read_targets.py first to create targets.hdf5 file."
        if gray_no_cell is None:
            gray_no_cell = color is not None
        targets = self.targets[~self.targets.target.str.startswith("NegPrb") &
            ~self.targets.target.str.startswith("FalseCode") &
            ~self.targets.target.str.startswith("SystemControl")]
        if fov is not None:
            targets = targets[targets.fov.isin(fov if isinstance(fov, list) else [fov])]
        z = None
        selected_z = None
        if self.is_3d and 'global_z' in targets.get_column_names():
            selected_z = self._resolve_z_indices(
                max_slices=len(self.z_slices),
                z_plane=z_plane,
                use_current_if_none=False,
            )
            targets = targets[targets.global_z.isin(selected_z)]
        # Materialize only after all filters are applied. Convert only fields that
        # require NumPy-specific behavior to reduce copies and keep memory lower.
        y = targets.global_y.evaluate()
        x = targets.global_x.evaluate()
        id = targets.CellId.evaluate()
        genes = targets.target.evaluate()
        if self.is_3d and 'global_z' in targets.get_column_names():
            z = targets.global_z.evaluate()
        u, inv = np.unique(genes, return_inverse=True)
        cm = label_colormap(len(self.targets.category_labels('target'))+1)
        face_color = cm.colors[1:][inv] if color is None else np.full(shape=(len(genes), 4), fill_value=transform_color(color))
        if gray_no_cell:
            face_color[id == 0] = transform_color('gray')
        if z is not None:
            points = np.column_stack((z, y, x))
        else:
            points = np.column_stack((y, x))
        print(f"Plotting {len(genes):,} transcripts")
        add_points_sig = signature(self.viewer.add_points)
        edge_kw = {}
        if 'edge_width' in add_points_sig.parameters:
            edge_kw['edge_width'] = 0
        elif 'border_width' in add_points_sig.parameters:
            edge_kw['border_width'] = 0
        points_layer = self.viewer.add_points(points,
            size = point_size,
            properties={
                'CellId': id,
                'Gene': genes
            },
            scale=self._image_scale(),
            face_color=face_color,
            translate=self._image_translate(),
            rotate=self.rotate,
            **edge_kw)
        points_layer.opacity = 0.7
        points_layer_name = "Targets"
        if selected_z is not None:
            points_layer_name = f"Targets [{self._format_rna_z_selection_suffix(selected_z)}]"
        points_layer.name = points_layer_name
        selected_z_values = None
        if selected_z is not None:
            selected_z_values = [int(self.z_slices[index]) for index in sorted(selected_z)]
        single_z_value = selected_z_values[0] if selected_z_values is not None and len(selected_z_values) == 1 else None
        self._tag_layer(points_layer, 'rna', 'Targets', single_z_value)
        if selected_z_values is not None:
            points_layer.metadata['cosmx_z_selection'] = selected_z_values
        return points_layer

    def layers_to_metadata(self:'napari_cosmx.gemini.Gemini', layer_names:'list', shape_z_plane=None, points_z_plane=None)-> 'pandas.core.frame.DataFrame':
        """Create metadata columns to indicate if cells overlap with the shapes or points in layers.

        Args:
            layer_names (list): Names of shapes or points layers, layer names will be used as column names in metadata.
            shape_z_plane (int, str, tuple, list, or range, optional): For 3D
                shape layers, use a single z plane, 'any', 'current', 'auto', or
                a z-plane range.
            points_z_plane (int, str, tuple, list, or range, optional): For 3D
                points layers, use a single z plane, 'any', 'current', 'auto', or
                a z-plane range.

        Returns:
            pandas.core.frame.DataFrame: Metadata object with additional columns appended, True if cell overlaps with shape, False otherwise.
        """
        assert self.metadata is not None, "No metadata found"
        meta = self.metadata.copy()
        invalid = [x for x in layer_names if x not in self.viewer.layers]
        assert not invalid, f"No layers found named {invalid}."
        layers = [self.viewer.layers[x] for x in layer_names]
        invalid = [x for x in layers if not isinstance(x, napari.layers.Shapes) and not isinstance(x, napari.layers.Points)]
        assert not invalid, f"These layers are not Shapes or Points layers: {invalid}."
        for layer in layer_names:
            meta[layer] = False
            if isinstance(self.viewer.layers[layer], napari.layers.Shapes):
                cells_in_layer = self.cells_in_shape(self.viewer.layers[layer], z_plane=shape_z_plane)
            elif isinstance(self.viewer.layers[layer], napari.layers.Points):
                cells_in_layer = self.cells_at_points(self.viewer.layers[layer], z_plane=points_z_plane)
            meta.loc[meta.index.isin(cells_in_layer), layer] = True
        return meta

    def save_layers(self:'napari_cosmx.gemini.Gemini', layer_names:'list', out_dir:'str'=None):
        """Save Shapes and Points layers as pickle objects. Re-open with load_layer.

        Args:
            layer_names (list): The names of the layers to save.
            out_dir (str, optional): Directory to save pickle file. Defaults to slide folder.
        """
        if out_dir is None:
            out_dir = self.folder
        invalid = [x for x in layer_names if x not in self.viewer.layers]
        assert not invalid, f"No layers found named {invalid}."
        layers = [self.viewer.layers[x] for x in layer_names]
        invalid = [x for x in layers if not isinstance(x, napari.layers.Shapes) and not isinstance(x, napari.layers.Points)]
        assert not invalid, f"These layers are not Shapes or Points layers: {invalid}."
        for layer_name in layer_names:
            print('Saving ' + layer_name + ' to ' + out_dir)
            layer_tuple = self.viewer.layers[layer_name].as_layer_data_tuple()
            # save a pickle
            with open(os.path.join(out_dir, layer_name+'.pickle'), 'wb') as f:
                pickle.dump(layer_tuple, f)

    def load_layers(self:'napari_cosmx.gemini.Gemini', layer_names:'list'=[], input_dir:'str'=None):
        """Loads previously saved layer (pickle file) into napari.

        Args:
            layer_names (list, optional): If not empty, load only these layers as long as they end with ".pickle". Defaults to [].
            input_dir (str, optional): directory of the pickle files. Defaults to slide folder.
        """
        if input_dir is None:
            input_dir = self.folder
        if len(layer_names) == 0:
            layer_names = [i for i in listdir(input_dir) if i.endswith(".pickle")]
        else:
            layer_names = [i+'.pickle' for i in layer_names]
        
        for layer_name in layer_names:
            with open(os.path.join(input_dir, layer_name), 'rb') as f:
                layer_tuple = pickle.load(f)
            layer_metadata = layer_tuple[1]
            print('Adding ' + layer_tuple[2] + ' for ' + layer_name)

            if layer_tuple[2] == 'shapes':
                valid_kwargs = signature(self.viewer.add_shapes).parameters
            elif layer_tuple[2] == 'points':
                valid_kwargs = signature(self.viewer.add_points).parameters

            layer_kwargs = {k:v for k,v in layer_metadata.items() if k in valid_kwargs}
            if version('napari') <= '0.4.16':
                layer_kwargs['text']['color'] = None

            if layer_tuple[2] == 'shapes':
                self.viewer.add_shapes(data=layer_tuple[0], **layer_kwargs)
            elif layer_tuple[2] == 'points':
                self.viewer.add_points(data=layer_tuple[0], **layer_kwargs)

