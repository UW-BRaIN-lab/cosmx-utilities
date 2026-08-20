from napari_cosmx.utils.stitch_images import main as stitch_images
from napari_cosmx.utils.stitch_composite import main as stitch_composite
from napari_cosmx.utils.create_anndata import main as create_anndata
from napari_cosmx import DEFAULT_Z_STEP_UM
from napari_cosmx._dock_widget import GeminiQWidget
from napari_cosmx.gemini import Gemini
from napari_cosmx.pairing import pair
import anndata as ad
import json
import numpy as np
import pandas as pd
import tifffile
import vaex
import zarr
from skimage import io
from qtpy.QtWidgets import QGroupBox
import pytest


def _write_morphology_tiff(path, channels):
     with tifffile.TiffWriter(path) as tif:
          for channel in channels:
               tif.write(channel)


def _write_morphology_tiff_with_binning(path, channels, binning):
     description = {
          'Channels': [
               {'UID': uid, 'Binning': int(binning)}
               for uid in ['B', 'G', 'Y', 'R', 'U']
          ]
     }
     with tifffile.TiffWriter(path) as tif:
          for index, channel in enumerate(channels):
               if index == 0:
                    tif.write(channel, description=json.dumps(description))
               else:
                    tif.write(channel)

def test_stitch_images(make_napari_viewer, tmp_path):
    stitch_images(args_list=['-i', 'data/Liver-S2-source/CellStatsDir',
         '-f', 'data/Liver-S2-source/RunSummary',
         '-o', tmp_path.as_posix()])
    viewer = make_napari_viewer()
    gem = Gemini(tmp_path.as_posix(), viewer=viewer)
    assert isinstance(gem, Gemini)

def test_stitch_images_3d_ingestion(make_napari_viewer, tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.ones((4, 4), dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z002.tif', np.full((4, 4), 2, dtype=np.uint16))
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.full((4, 4), 10, dtype=np.uint16),
          np.full((4, 4), 20, dtype=np.uint16),
     ])
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z002.TIF', [
          np.full((4, 4), 30, dtype=np.uint16),
          np.full((4, 4), 40, dtype=np.uint16),
     ])

     targets = vaex.from_pandas(pd.DataFrame({
          'fov': [1, 1],
          'CellId': [1, 2],
          'x': [1.0, 2.0],
          'y': [1.0, 2.0],
          'z': [1, 2],
          'target': pd.Series(['GeneA', 'GeneB'], dtype='category'),
          'CellComp': pd.Series(['Nucleus', 'Cytoplasm'], dtype='category'),
     }).assign())
     targets = targets._future()
     targets.ordinal_encode('target', inplace=True)
     targets.ordinal_encode('CellComp', inplace=True)
     targets.export_hdf5((tmp_path / 'targets.hdf5').as_posix(), mode='w')

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', tmp_path.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '3',
          '--z-step-um', '0.8',
     ])

     grp = zarr.open((tmp_path / 'images').as_posix(), mode='r')
     assert grp.attrs['CosMx']['ndim'] == 3
     assert grp.attrs['CosMx']['z_slices'] == [1, 2]
     assert grp['labels']['0'].shape == (2, 4, 4)
     assert grp['B']['0'].shape == (2, 4, 4)

     viewer = make_napari_viewer()
     gem = Gemini(tmp_path.as_posix(), viewer=viewer)
     assert gem.is_3d is True
     assert 'global_z' in gem.targets.get_column_names()
     assert gem.targets.global_z.evaluate().tolist() == [0, 1]
     assert gem.fov_labels is not None
     assert len(gem.fov_labels.data) == len(gem.z_slices) * len(gem.fov_offsets)
     assert gem.fov_labels.data[0].shape[1] == 3
     shape = viewer.add_shapes([
          np.array([[0, 0], [0, 3], [3, 3], [3, 0]])
     ], shape_type='polygon', name='ROI')
     assert set(gem.cells_in_shape(shape, z_plane=1).tolist()) == {pair(1, 1)}
     assert set(gem.cells_in_shape(shape, z_plane=2).tolist()) == {pair(1, 2)}
     assert set(gem.cells_in_shape(shape, z_plane=(1, 2)).tolist()) == {pair(1, 1), pair(1, 2)}
     assert set(gem.cells_in_shape(shape, z_plane='any').tolist()) == {pair(1, 1), pair(1, 2)}
     all_points = gem.add_points(point_size=1, z_plane='any')
     assert all_points.data.shape[0] == 2
     z2_points = gem.add_points(point_size=1, z_plane=(2, 2))
     assert z2_points.data.shape[0] == 1
     gene_a_points = gem.plot_transcripts('GeneA', color='white', point_size=1, z_plane='any')
     assert gene_a_points.data.shape[0] == 1
     gene_b_z2_points = gem.plot_transcripts('GeneB', color='white', point_size=1, z_plane=(2, 2))
     assert gene_b_z2_points.data.shape[0] == 1
     viewer.dims.ndisplay = 2
     if hasattr(viewer.dims, 'set_current_step'):
          viewer.dims.set_current_step(0, 0)
     else:
          viewer.dims.current_step = (0, 0, 0)
     default_points_2d = gem.add_points(point_size=1)
     # 3D data defaults to all z-slices regardless of display mode
     assert default_points_2d.data.shape[0] == 2
     viewer.dims.ndisplay = 3
     default_points_3d = gem.add_points(point_size=1)
     assert default_points_3d.data.shape[0] == 2
     z1_layer = gem.add_channel('B', z_plane=1)
     z2_layer = gem.add_channel('B', z_plane=2)
     assert z1_layer.name == 'B [z=1]'
     assert z2_layer.name == 'B [z=2]'
     z1_data = z1_layer.data[0] if isinstance(z1_layer.data, list) else z1_layer.data
     z2_data = z2_layer.data[0] if isinstance(z2_layer.data, list) else z2_layer.data
     assert z1_data.shape == (4, 4)
     assert z2_data.shape == (4, 4)
     z1_labels = gem.add_cell_labels(z_plane=1)
     z2_labels = gem.add_cell_labels(z_plane=2)
     z1_seg = gem.add_segmentation(z_plane=1)
     assert z1_labels.name == 'Metadata [z=1]'
     assert z2_labels.name == 'Metadata [z=2]'
     assert z1_seg.name == 'Segmentation [z=1]'
     z1_labels_data = z1_labels.data[0] if isinstance(z1_labels.data, list) else z1_labels.data
     z1_seg_data = z1_seg.data[0] if isinstance(z1_seg.data, list) else z1_seg.data
     assert z1_labels_data.shape == (4, 4)
     assert z1_seg_data.shape == (4, 4)
     gem.color_cells('all', color='white', labels_layer='Metadata [z=1]')
     assert len(gem.color_cells_layer.scale) == 2
     assert gem.color_cells_layer.name == 'Cells [z=1]'
     gem.color_cells('all', color='white')
     assert len(gem.color_cells_layer.scale) == 3
     assert gem.color_cells_layer.name == 'Cells'

def test_gemini_multiomics_capabilities(make_napari_viewer, tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     offsetsdir = tmp_path / 'RunSummary'
     inputdir.mkdir()
     offsetsdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001.tif', np.ones((4, 4), dtype=np.uint16))
     _write_morphology_tiff(inputdir / 'Morphology_C902_P99_N99_F001.TIF', [
          np.full((4, 4), 10, dtype=np.uint16),
          np.full((4, 4), 20, dtype=np.uint16),
     ])

     targets = vaex.from_pandas(pd.DataFrame({
          'fov': [1],
          'CellId': [1],
          'x': [1.0],
          'y': [1.0],
          'target': pd.Series(['GeneA'], dtype='category'),
          'CellComp': pd.Series(['Nucleus'], dtype='category'),
     }).assign())
     targets = targets._future()
     targets.ordinal_encode('target', inplace=True)
     targets.ordinal_encode('CellComp', inplace=True)
     targets.export_hdf5((tmp_path / 'targets.hdf5').as_posix(), mode='w')

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '--imagesdir', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', tmp_path.as_posix(),
     ])

     grp = zarr.open((tmp_path / 'images').as_posix(), mode='r+')
     protein_root = grp.create_group('protein')
     protein_root.attrs['CosMx'] = {
          'scale': 2.0,
          'z_step_um': DEFAULT_Z_STEP_UM,
     }
     protein_channel = protein_root.create_group('CD3')
     protein_channel.create_dataset('0', shape=(4, 4), dtype=np.uint16, data=np.full((4, 4), 25, dtype=np.uint16))
     protein_channel.attrs['multiscales'] = [{
          'datasets': [{'path': '0'}]
     }]

     viewer = make_napari_viewer()
     gem = Gemini(tmp_path.as_posix(), viewer=viewer)

     assert gem.has_rna is True
     assert gem.has_protein is True
     assert gem.is_multiomics is True
     assert gem.modalities == ('rna', 'protein')
     assert gem.genes == ['GeneA']
     assert gem.proteins == ['CD3']
     rna_layer = gem.plot_transcripts('GeneA', color='white', point_size=1)
     assert rna_layer.data.shape[1] == 2
     protein_layer = gem.add_protein('CD3')
     assert protein_layer.data[0].shape == (4, 4)
     widget = GeminiQWidget(viewer, gem)
     titles = {group_box.title() for group_box in widget.findChildren(QGroupBox)}
     assert 'Dataset' in titles
     assert 'RNA Transcripts' in titles
     assert 'Protein Expression' in titles

def test_stitch_images_3d_input_2d_output_single_z(tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     outdir = tmp_path / 'out2d'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     outdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.ones((4, 4), dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z002.tif', np.full((4, 4), 2, dtype=np.uint16))
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.full((4, 4), 10, dtype=np.uint16),
          np.full((4, 4), 20, dtype=np.uint16),
     ])
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z002.TIF', [
          np.full((4, 4), 30, dtype=np.uint16),
          np.full((4, 4), 40, dtype=np.uint16),
     ])

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', outdir.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '2',
          '--zslice', '2',
     ])

     grp = zarr.open((outdir / 'images').as_posix(), mode='r')
     assert grp.attrs['CosMx']['input_ndim'] == 3
     assert grp.attrs['CosMx']['ndim'] == 2
     assert grp.attrs['CosMx']['z_slices'] == [2]
     assert grp['labels']['0'].shape == (4, 4)

def test_stitch_images_3d_noncontiguous_z_mapping(make_napari_viewer, tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.full((4, 4), 1, dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z003.tif', np.full((4, 4), 3, dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z004.tif', np.full((4, 4), 4, dtype=np.uint16))
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.full((4, 4), 10, dtype=np.uint16),
          np.full((4, 4), 20, dtype=np.uint16),
     ])
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z003.TIF', [
          np.full((4, 4), 30, dtype=np.uint16),
          np.full((4, 4), 40, dtype=np.uint16),
     ])
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z004.TIF', [
          np.full((4, 4), 50, dtype=np.uint16),
          np.full((4, 4), 60, dtype=np.uint16),
     ])

     targets = vaex.from_pandas(pd.DataFrame({
          'fov': [1, 1, 1],
          'CellId': [1, 3, 4],
          'x': [1.0, 2.0, 3.0],
          'y': [1.0, 2.0, 3.0],
          'z': [1, 3, 4],
          'target': pd.Series(['GeneA', 'GeneB', 'GeneC'], dtype='category'),
          'CellComp': pd.Series(['Nucleus', 'Cytoplasm', 'Membrane'], dtype='category'),
     }).assign())
     targets = targets._future()
     targets.ordinal_encode('target', inplace=True)
     targets.ordinal_encode('CellComp', inplace=True)
     targets.export_hdf5((tmp_path / 'targets.hdf5').as_posix(), mode='w')

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', tmp_path.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '3',
     ])

     viewer = make_napari_viewer()
     gem = Gemini(tmp_path.as_posix(), viewer=viewer)
     assert gem.z_slices == [1, 3, 4]
     assert gem.z_plane_mapping() == {0: 1, 1: 3, 2: 4}
     assert gem.z_plane_to_index(3) == 1
     shape = viewer.add_shapes([
          np.array([[0, 0], [0, 3], [3, 3], [3, 0]])
     ], shape_type='polygon', name='ROI-noncontig')
     assert set(gem.cells_in_shape(shape, z_plane=(1, 4)).tolist()) == {pair(1, 1), pair(1, 3), pair(1, 4)}
     with pytest.raises(ValueError):
          gem.cells_in_shape(shape, z_plane=(2, 2))


def test_stitch_images_3d_binned_morphology_metadata_detection(tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.ones((4, 4), dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z002.tif', np.full((4, 4), 2, dtype=np.uint16))
     _write_morphology_tiff_with_binning(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.array([[10, 11], [12, 13]], dtype=np.uint16),
          np.array([[20, 21], [22, 23]], dtype=np.uint16),
     ], binning=2)
     _write_morphology_tiff_with_binning(morphology_dir / 'Morphology_C902_P99_N99_F001_Z002.TIF', [
          np.array([[30, 31], [32, 33]], dtype=np.uint16),
          np.array([[40, 41], [42, 43]], dtype=np.uint16),
     ], binning=2)

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', tmp_path.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '3',
     ])

     grp = zarr.open((tmp_path / 'images').as_posix(), mode='r')
     assert grp['labels']['0'].shape == (2, 4, 4)
     assert grp['B']['0'].shape == (2, 4, 4)
     # nearest-neighbor upscale preserves 2x2 blocks from binned source
     b_z1 = grp['B']['0'][0]
     assert b_z1[0, 0] == 10
     assert b_z1[0, 1] == 10
     assert b_z1[1, 0] == 10
     assert b_z1[2, 2] == 13


def test_stitch_images_3d_binned_morphology_shape_fallback(tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     outdir = tmp_path / 'out_shape_fallback'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     outdir.mkdir()
     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.ones((4, 4), dtype=np.uint16))
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.array([[7, 8], [9, 10]], dtype=np.uint16),
          np.array([[17, 18], [19, 20]], dtype=np.uint16),
     ])

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', outdir.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '3',
     ])

     grp = zarr.open((outdir / 'images').as_posix(), mode='r')
     assert grp['B']['0'].shape == (1, 4, 4)
     b = grp['B']['0'][0]
     assert b[0, 0] == 7
     assert b[1, 1] == 7
     assert b[3, 3] == 10

def test_create_anndata_3d_metadata(tmp_path):
     obs_path = tmp_path / 'obs.csv'
     coords_path = tmp_path / 'coords.csv'
     outdir = tmp_path / 'out'
     outdir.mkdir()

     pd.DataFrame({'cell_type': ['A', 'B']}, index=['c_1_1', 'c_1_2']).to_csv(obs_path)
     pd.DataFrame({'z': [0, 1], 'y': [10, 20], 'x': [30, 40]}).to_csv(coords_path, index=False)

     create_anndata(args_list=[
          '--obs', obs_path.as_posix(),
          '--coords', coords_path.as_posix(),
          '--ndim', '3',
          '--z-step-um', '0.8',
          '--outputdir', outdir.as_posix(),
     ])

     adata = ad.read_h5ad((outdir / 'adata.h5ad').as_posix())
     assert adata.obsm['spatial'].shape == (2, 3)
     assert adata.uns['CosMx']['ndim'] == 3
     assert adata.uns['CosMx']['z_step_um'] == 0.8

def test_stitch_composite_2d_defaults(tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     offsetsdir = tmp_path / 'RunSummary'
     outdir = tmp_path / 'out'
     inputdir.mkdir()
     offsetsdir.mkdir()
     outdir.mkdir()

     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)
     tifffile.imwrite(inputdir / 'CellLabels_F001.tif', np.ones((4, 4), dtype=np.uint16))

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '--imagesdir', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', outdir.as_posix(),
          '--labels',
     ])

     grp = zarr.open((outdir / 'images').as_posix(), mode='r+')
     meta = dict(grp.attrs['CosMx'])
     meta.pop('ndim', None)
     meta.pop('z_step_um', None)
     grp.attrs['CosMx'] = meta

     tile = np.zeros((4, 4, 3), dtype=np.uint8)
     tile[:, :, 0] = 255
     io.imsave((inputdir / 'CellComposite_F001.jpg').as_posix(), tile, check_contrast=False)

     stitch_composite(args_list=['-i', inputdir.as_posix(), '-o', outdir.as_posix()])

     out = zarr.open((outdir / 'images').as_posix(), mode='r')
     assert out['composite']['0'].shape == (4, 4, 3)

def test_stitch_composite_3d_default_z_step(tmp_path):
     inputdir = tmp_path / 'CellStatsDir'
     morphology_dir = inputdir / 'Morphology3D'
     offsetsdir = tmp_path / 'RunSummary'
     outdir = tmp_path / 'out3d'
     inputdir.mkdir()
     morphology_dir.mkdir()
     offsetsdir.mkdir()
     outdir.mkdir()

     pd.DataFrame([[0, 0.0, 0.0, 0.0, 0.0, 0, 1, 0]]).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)

     tifffile.imwrite(inputdir / 'CellLabels_F001_Z001.tif', np.ones((4, 4), dtype=np.uint16))
     tifffile.imwrite(inputdir / 'CellLabels_F001_Z002.tif', np.full((4, 4), 2, dtype=np.uint16))
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z001.TIF', [
          np.full((4, 4), 10, dtype=np.uint16),
          np.full((4, 4), 20, dtype=np.uint16),
     ])
     _write_morphology_tiff(morphology_dir / 'Morphology_C902_P99_N99_F001_Z002.TIF', [
          np.full((4, 4), 30, dtype=np.uint16),
          np.full((4, 4), 40, dtype=np.uint16),
     ])

     stitch_images(args_list=[
          '-i', inputdir.as_posix(),
          '-f', offsetsdir.as_posix(),
          '-o', outdir.as_posix(),
          '--input-ndim', '3',
          '--output-ndim', '3',
          '--z-step-um', '0.8',
     ])

     grp = zarr.open((outdir / 'images').as_posix(), mode='r+')
     meta = dict(grp.attrs['CosMx'])
     meta.pop('z_step_um', None)
     grp.attrs['CosMx'] = meta

     tile_z1 = np.zeros((4, 4, 3), dtype=np.uint8)
     tile_z2 = np.zeros((4, 4, 3), dtype=np.uint8)
     tile_z1[:, :, 1] = 127
     tile_z2[:, :, 2] = 255
     io.imsave((inputdir / 'CellComposite_F001_Z001.jpg').as_posix(), tile_z1, check_contrast=False)
     io.imsave((inputdir / 'CellComposite_F001_Z002.jpg').as_posix(), tile_z2, check_contrast=False)

     stitch_composite(args_list=['-i', inputdir.as_posix(), '-o', outdir.as_posix()])

     out = zarr.open((outdir / 'images').as_posix(), mode='r')
     assert out['composite']['0'].shape == (2, 4, 4, 3)
     multiscales = out['composite'].attrs['multiscales'][0]['datasets']
     assert multiscales[0]['coordinateTransformations'][0]['scale'][0] == DEFAULT_Z_STEP_UM
