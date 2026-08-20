import os
import numpy as np
import pandas as pd
import zarr
import pytest

from napari_cosmx.gemini import Gemini
from napari_cosmx import DEFAULT_Z_STEP_UM, DEFAULT_NDIM

_REAL_DATA = 'data/Liver-S2'


# --- Synthetic zarr helper ---

def _make_gemini_zarr(
    tmp_path,
    *,
    ndim=2,
    z_slices=None,
    channels=('B', 'DNA'),
    with_labels=True,
    with_protein=False,
    protein_channels=('CD3',),
    scale_um=0.12,
):
    """Build a minimal zarr images tree that Gemini can open without running
    the stitch pipeline.

    Parameters
    ----------
    scale_um:
        Pass ``None`` to omit the optional scale key and exercise the
        instrument-default fallback path.
    z_slices:
        For 3D stores, the list of physical z-plane values.  Also sets
        ``z_step_um=0.8`` and records ``ndim=3`` in the CosMx attrs.
    """
    images_path = tmp_path / 'images'
    grp = zarr.open(images_path.as_posix(), mode='w')

    fov_offsets = pd.DataFrame(
        {'FOV': [1], 'X_mm': [0.0], 'Y_mm': [0.0]}
    ).to_dict()
    cosmx = {
        'fov_height': 4,
        'fov_width': 4,
        'fov_offsets': fov_offsets,
        'ndim': ndim,
    }
    if scale_um is not None:
        cosmx['scale_um'] = scale_um
    if z_slices is not None:
        cosmx['z_slices'] = list(z_slices)
        cosmx['z_step_um'] = 0.8
    grp.attrs['CosMx'] = cosmx

    use_3d = ndim == 3 and z_slices is not None
    nz = len(z_slices) if z_slices is not None else 1
    shape = (nz, 4, 4) if use_3d else (4, 4)

    def _add_image_group(parent, name):
        g = parent.create_group(name)
        # Non-zero fill so omero() percentile does not operate on an empty array
        g.create_dataset('0', data=np.full(shape, 100, dtype=np.uint16))
        g.attrs['multiscales'] = [{'datasets': [{'path': '0'}]}]
        return g

    for ch in channels:
        _add_image_group(grp, ch)

    if with_labels:
        lab = grp.create_group('labels')
        lab.create_dataset('0', data=np.zeros(shape, dtype=np.uint16))
        lab.attrs['multiscales'] = [{'datasets': [{'path': '0'}]}]

    if with_protein:
        prot_root = grp.create_group('protein')
        prot_root.attrs['CosMx'] = {'scale': 1.0, 'z_step_um': 0.8}
        for pch in protein_channels:
            _add_image_group(prot_root, pch)

    return tmp_path


# --- Real-data tests (guarded by data availability) ---

@pytest.mark.skipif(not os.path.exists(_REAL_DATA), reason='Liver-S2 data not found')
def test_launch_viewer(make_napari_viewer):
    viewer = make_napari_viewer()
    gem = Gemini(_REAL_DATA, viewer=viewer)
    assert isinstance(gem, Gemini)


@pytest.mark.skipif(not os.path.exists(_REAL_DATA), reason='Liver-S2 data not found')
def test_add_channel(make_napari_viewer):
    viewer = make_napari_viewer()
    gem = Gemini(_REAL_DATA, viewer=viewer)
    gem.add_channel('DNA', colormap='blue')
    assert 'DNA' in viewer.layers, 'Morphology layer not added to viewer'
    cm = viewer.layers['DNA'].colormap.name
    assert cm == 'blue', f'Colormap is {cm}, not blue as requested'


# --- Synthetic 2-D property tests ---

def test_gemini_2d_synthetic_properties(make_napari_viewer, tmp_path):
    """2-D synthetic store without targets or protein: verify key properties."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B', 'DNA'), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.is_3d is False
    assert gem.z_slices == [0]           # default when z_slices not in attrs
    assert gem.has_rna is False          # no targets.hdf5
    assert gem.has_protein is False
    assert gem.is_multiomics is False
    assert gem.is_protein is False
    assert gem.modalities == ()
    assert set(gem.channels) == {'B', 'DNA'}
    assert gem.proteins == []


def test_gemini_backward_compat_no_scale_um(make_napari_viewer, tmp_path):
    """Old stores without scale_um/scale_mm load gracefully with instrument defaults."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), scale_um=None)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.is_3d is False
    assert gem.z_slices == [0]
    assert gem.z_step_um == DEFAULT_Z_STEP_UM
    assert gem.mm_per_px > 0


def test_gemini_backward_compat_no_ndim(make_napari_viewer, tmp_path):
    """Stores without ndim/z_slices/z_step_um keys default to 2-D behaviour."""
    images_path = tmp_path / 'images'
    grp = zarr.open(images_path.as_posix(), mode='w')
    fov_offsets = pd.DataFrame({'FOV': [1], 'X_mm': [0.0], 'Y_mm': [0.0]}).to_dict()
    grp.attrs['CosMx'] = {
        'fov_height': 4,
        'fov_width': 4,
        'fov_offsets': fov_offsets,
        'scale_um': 0.12,
        # intentionally missing: 'ndim', 'z_slices', 'z_step_um'
    }
    ch_grp = grp.create_group('B')
    ch_grp.create_dataset('0', data=np.full((4, 4), 100, dtype=np.uint16))
    ch_grp.attrs['multiscales'] = [{'datasets': [{'path': '0'}]}]

    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.ndim == DEFAULT_NDIM
    assert gem.is_3d is False
    assert gem.z_slices == [0]
    assert gem.z_step_um == DEFAULT_Z_STEP_UM


# --- Layer metadata tag tests ---

def test_gemini_layer_metadata_tags_2d(make_napari_viewer, tmp_path):
    """_tag_layer attaches cosmx provenance metadata to napari layers (2-D case)."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_labels=True)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    ch_layer = gem.add_channel('B')
    assert ch_layer.metadata['cosmx_modality'] == 'morphology'
    assert ch_layer.metadata['cosmx_channel'] == 'B'
    assert ch_layer.metadata['cosmx_dataset'] == tmp_path.as_posix()
    assert 'cosmx_z_plane' not in ch_layer.metadata  # 2-D: no z tag

    seg_layer = gem.add_segmentation()
    assert seg_layer.metadata['cosmx_modality'] == 'segmentation'
    assert 'cosmx_channel' not in seg_layer.metadata
    assert 'cosmx_z_plane' not in seg_layer.metadata


def test_gemini_layer_metadata_tags_3d_z_plane(make_napari_viewer, tmp_path):
    """3-D z_plane slices receive a cosmx_z_plane tag; full 3-D layers do not."""
    _make_gemini_zarr(tmp_path, ndim=3, z_slices=[1, 2], channels=('B',), with_labels=True)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    ch_layer = gem.add_channel('B', z_plane=1)
    assert ch_layer.metadata['cosmx_z_plane'] == 1
    assert ch_layer.metadata['cosmx_channel'] == 'B'
    assert ch_layer.name == 'B [z=1]'

    seg_layer = gem.add_segmentation(z_plane=2)
    assert seg_layer.metadata['cosmx_z_plane'] == 2
    assert seg_layer.name == 'Segmentation [z=2]'

    full_layer = gem.add_channel('B')
    assert 'cosmx_z_plane' not in full_layer.metadata
    assert full_layer.name == 'B'


# --- 3-D capability tests ---

def test_gemini_3d_z_plane_mapping(make_napari_viewer, tmp_path):
    """z_plane_mapping() and z_plane_to_index() reflect non-contiguous z slices."""
    _make_gemini_zarr(tmp_path, ndim=3, z_slices=[1, 3, 5], channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.is_3d is True
    assert gem.z_slices == [1, 3, 5]
    assert gem.z_plane_mapping() == {0: 1, 1: 3, 2: 5}
    assert gem.z_plane_to_index(1) == 0
    assert gem.z_plane_to_index(3) == 1
    assert gem.z_plane_to_index(5) == 2
    with pytest.raises(ValueError):
        gem.z_plane_to_index(2)  # 2 not in z_slices


def test_gemini_3d_add_channel_full_and_z_sliced(make_napari_viewer, tmp_path):
    """add_channel returns a 3-D layer by default and a 2-D slice when z_plane is given."""
    _make_gemini_zarr(tmp_path, ndim=3, z_slices=[1, 2], channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    full_layer = gem.add_channel('B')
    assert full_layer.data[0].shape == (2, 4, 4)  # ZYX
    assert full_layer.name == 'B'
    assert 'cosmx_z_plane' not in full_layer.metadata

    z1_layer = gem.add_channel('B', z_plane=1)
    assert z1_layer.data[0].shape == (4, 4)        # YX slice
    assert z1_layer.name == 'B [z=1]'
    assert z1_layer.metadata['cosmx_z_plane'] == 1

    z2_layer = gem.add_channel('B', z_plane=2)
    assert z2_layer.name == 'B [z=2]'
    assert z2_layer.metadata['cosmx_z_plane'] == 2


# --- Multi-omics property tests ---

def test_gemini_protein_only_properties(make_napari_viewer, tmp_path):
    """has_protein, proteins list, and add_protein layer tags work without RNA."""
    _make_gemini_zarr(
        tmp_path, ndim=2, channels=('B',), with_labels=False,
        with_protein=True, protein_channels=('CD3', 'CD8'),
    )
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.has_protein is True
    assert gem.is_protein is True          # backward-compat alias
    assert sorted(gem.proteins) == ['CD3', 'CD8']
    assert gem.has_rna is False
    assert gem.is_multiomics is False
    assert gem.modalities == ('protein',)

    prot_layer = gem.add_protein('CD3')
    assert prot_layer.metadata['cosmx_modality'] == 'protein'
    assert prot_layer.metadata['cosmx_channel'] == 'CD3'
    assert prot_layer.metadata['cosmx_dataset'] == tmp_path.as_posix()


def test_gemini_resolve_layer_name_without_collision(make_napari_viewer, tmp_path):
    """_resolve_layer_name returns the base name when there is no name collision."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_protein=True, protein_channels=('CD3',))
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    # Not multiomics (no targets.hdf5) - no modality tag regardless
    assert gem._resolve_layer_name('B', 'rna') == 'B'
    assert gem._resolve_layer_name('CD3', 'protein') == 'CD3'
    # z_value suffix is always appended when provided
    assert gem._resolve_layer_name('B', 'morphology', z_value=3) == 'B [z=3]'


def test_gemini_resolve_layer_name_with_collision(make_napari_viewer, tmp_path):
    """When gene and protein share a name in multiomics, _resolve_layer_name disambiguates."""
    import vaex

    # 'B' present in both morphology channels and protein group -> name collision
    _make_gemini_zarr(
        tmp_path, ndim=2, channels=('B', 'DNA'),
        with_labels=False, with_protein=True, protein_channels=('B',),
    )
    # Minimal targets.hdf5 so that has_rna=True and gene 'B' appears in targets
    df = pd.DataFrame({
        'fov': [1],
        'CellId': [1],
        'x': [1.0],
        'y': [1.0],
        'target': pd.Series(['B'], dtype='category'),
        'CellComp': pd.Series(['Nucleus'], dtype='category'),
    })
    t = vaex.from_pandas(df.assign())
    t = t._future()
    t.ordinal_encode('target', inplace=True)
    t.ordinal_encode('CellComp', inplace=True)
    t.export_hdf5((tmp_path / 'targets.hdf5').as_posix(), mode='w')

    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.is_multiomics is True
    assert 'B' in gem.genes
    assert 'B' in gem.proteins

    # Collision - modality tag added
    assert gem._resolve_layer_name('B', 'rna') == 'B [rna]'
    assert gem._resolve_layer_name('B', 'protein') == 'B [protein]'
    # No collision for 'DNA' (only in channels, not proteins)
    assert gem._resolve_layer_name('DNA', 'rna') == 'DNA'

    # add_protein uses _resolve_layer_name internally
    prot_layer = gem.add_protein('B')
    assert prot_layer.name == 'B [protein]'


def test_gemini_resolve_layer_name_with_morphology_protein_collision(make_napari_viewer, tmp_path):
    """Protein names colliding with morphology channels are disambiguated."""
    _make_gemini_zarr(
        tmp_path, ndim=2, channels=('CD45',),
        with_labels=False, with_protein=True, protein_channels=('CD45',),
    )
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    assert gem.has_rna is False
    assert gem._resolve_layer_name('CD45', 'protein') == 'CD45 [protein]'

    prot_layer = gem.add_protein('CD45')
    assert prot_layer.name == 'CD45 [protein]'


def test_plot_transcripts_3d_default_all_and_naming(make_napari_viewer, tmp_path):
    """3-D transcript plotting defaults to all z slices and names layers by z selection."""
    import vaex

    _make_gemini_zarr(tmp_path, ndim=3, z_slices=[1, 2], channels=('B',), with_labels=False)
    df = pd.DataFrame({
        'fov': [1, 1],
        'CellId': [1, 2],
        'x': [1.0, 2.0],
        'y': [1.0, 2.0],
        'z': [1, 2],
        'target': pd.Series(['GeneA', 'GeneA'], dtype='category'),
        'CellComp': pd.Series(['Nucleus', 'Nucleus'], dtype='category'),
    })
    t = vaex.from_pandas(df.assign())
    t = t._future()
    t.ordinal_encode('target', inplace=True)
    t.ordinal_encode('CellComp', inplace=True)
    t.export_hdf5((tmp_path / 'targets.hdf5').as_posix(), mode='w')

    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)
    viewer.dims.ndisplay = 2
    if hasattr(viewer.dims, 'set_current_step'):
        viewer.dims.set_current_step(0, 0)

    default_layer = gem.plot_transcripts('GeneA', color='white', point_size=3)
    assert default_layer.data.shape[0] == 2
    assert default_layer.name == 'GeneA [3D all-z]'
    assert default_layer.metadata['cosmx_z_selection'] == [1, 2]

    single_layer = gem.plot_transcripts('GeneA', color='white', point_size=3, z_plane=1)
    assert single_layer.data.shape[0] == 1
    assert single_layer.name == 'GeneA [3D z=1]'
    assert single_layer.metadata['cosmx_z_selection'] == [1]
    assert single_layer.metadata['cosmx_z_plane'] == 1


def test_color_cells_all_accepts_string_color_name(make_napari_viewer, tmp_path):
    """Public API accepts a string color name for col_name='all'."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    labels = [
        np.array([[0, 1], [1, 0]], dtype=np.uint32),
        np.array([[0]], dtype=np.uint32),
    ]
    legacy_labels_layer = viewer.add_labels(labels, name='LegacyLabels')

    gem.color_cells('all', color='red', labels_layer=legacy_labels_layer)

    assert gem.color_cells_layer is not None
    assert gem.color_cells_layer in viewer.layers
    assert gem.color_cells_layer.name == 'Cells'


def test_color_cells_categorical_accepts_dict_color_names(make_napari_viewer, tmp_path):
    """Public API accepts categorical color dict with string color names."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    labels = [
        np.array([[0, 1], [1, 0]], dtype=np.uint32),
        np.array([[0]], dtype=np.uint32),
    ]
    legacy_labels_layer = viewer.add_labels(labels, name='LegacyLabels')
    gem.metadata = pd.DataFrame({'major': ['A']}, index=pd.Index([1], name='UID'))

    gem.color_cells('major', color={'A': 'red'}, labels_layer=legacy_labels_layer)

    assert gem.color_cells_layer is not None
    assert gem.color_cells_layer in viewer.layers
    assert gem.color_cells_layer.name == 'major'


def test_color_cells_continuous_accepts_colormap_name(make_napari_viewer, tmp_path):
    """Public API accepts colormap name for continuous metadata columns."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    # 30 unique numeric values forces the continuous metadata path.
    labels_0 = np.arange(1, 31, dtype=np.uint32).reshape(5, 6)
    labels_1 = np.arange(1, 13, dtype=np.uint32).reshape(3, 4)
    legacy_labels_layer = viewer.add_labels([labels_0, labels_1], name='LegacyLabels')
    gem.metadata = pd.DataFrame(
        {'score': np.linspace(0.0, 1.0, 30, dtype=np.float32)},
        index=pd.Index(np.arange(1, 31), name='UID')
    )

    gem.color_cells('score', color='gray', labels_layer=legacy_labels_layer)

    assert gem.color_cells_layer is not None
    assert gem.color_cells_layer in viewer.layers
    assert gem.color_cells_layer.name == 'score'


def test_color_cells_preserves_z_suffix_in_layer_name(make_napari_viewer, tmp_path):
    """When coloring a z-specific labels layer, include that z in output name."""
    _make_gemini_zarr(tmp_path, ndim=2, channels=('B',), with_labels=False)
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)

    labels = [
        np.array([[0, 1], [1, 0]], dtype=np.uint32),
        np.array([[0]], dtype=np.uint32),
    ]
    legacy_labels_layer = viewer.add_labels(labels, name='LegacyLabels [z=7]')
    gem.metadata = pd.DataFrame({'major': ['A']}, index=pd.Index([1], name='UID'))

    gem.color_cells('major', color={'A': 'red'}, labels_layer=legacy_labels_layer)

    assert gem.color_cells_layer is not None
    assert gem.color_cells_layer in viewer.layers
    assert gem.color_cells_layer.name == 'major [z=7]'