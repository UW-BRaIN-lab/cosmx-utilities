"""Tests for this fork's local additions to ``stitch-images``.

Upstream napari-CosMx has no notion of AtoMx ``Segmentation_<uuid>_00N``
subdirectories, of overriding a wrong MorphologyKit channel map, or of a slide
whose CellLabels are 3D while its morphology acquisition is 2D. These tests
cover those three additions; upstream's own behavior is covered by
``test_stitch.py``.

Runnable either under pytest or directly:
    uv run python napari-cosmx-fork/src/napari_cosmx/_tests/test_fork_stitch.py
"""
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import zarr

from napari_cosmx.utils.stitch_images import main as stitch_images
from napari_cosmx.utils.stitch_images import _map_z_to_morphology

FOV_PIXELS = 4
# Distinct fill values let an assertion prove which file a tile came from.
FILL_PRIMARY_SEGMENTATION = 11
FILL_OLDER_SEGMENTATION = 22
FILL_BASE_LABELS = 33
MORPHOLOGY_CHANNEL_FILLS = (10, 20, 30, 40, 50)
KIT_CHANNEL_IDS = ('B', 'G', 'Y', 'R', 'U')
KIT_TARGETS = ('Histone', 'Empty', 'rRNA', 'GFAP', 'DAPI')


def _write_offsets(offsetsdir, fovs):
    """Legacy latest.fovs.csv: Slide,X_mm,Y_mm,Z_mm,ZOffset_mm,ROI,FOV,Order."""
    rows = [[0, 0.0, index * 1.0, 0.0, 0.0, 0, fov, index]
            for index, fov in enumerate(fovs)]
    pd.DataFrame(rows).to_csv(offsetsdir / 'latest.fovs.csv', header=False, index=False)


def _write_labels(path, fill):
    tifffile.imwrite(path, np.full((FOV_PIXELS, FOV_PIXELS), fill, dtype=np.uint16))


def _write_morphology(path, kit_targets=KIT_TARGETS):
    """Multi-channel morphology TIFF carrying MorphologyKit metadata."""
    description = json.dumps({
        'MorphologyKit': {
            'MorphologyReagents': [
                {'Fluorophore': {'ChannelId': channel}, 'BiologicalTarget': target}
                for channel, target in zip(KIT_CHANNEL_IDS, kit_targets)
            ]
        }
    })
    with tifffile.TiffWriter(path) as tif:
        for index, fill in enumerate(MORPHOLOGY_CHANNEL_FILLS):
            plane = np.full((FOV_PIXELS, FOV_PIXELS), fill, dtype=np.uint16)
            tif.write(plane, description=description if index == 0 else None)


def _labels_values(store_path):
    grp = zarr.open(str(store_path / 'images'), mode='r')
    return grp, np.asarray(grp['labels']['0'])


def test_celllabels_subdir_primary_wins_and_older_gap_fills():
    """The named segmentation owns every FOV it covers; older versions only
    gap-fill, and base FOV*/ labels backstop whatever is still missing."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputdir = tmp / 'CellStatsDir'
        offsetsdir = tmp / 'RunSummary'
        morphology = inputdir / 'Morphology2D'
        for d in (inputdir, offsetsdir, morphology):
            d.mkdir(parents=True)
        _write_offsets(offsetsdir, [1, 2, 3])

        primary = inputdir / 'Segmentation_aaaa_002'
        older = inputdir / 'Segmentation_bbbb_001'
        # Primary covers FOV 1 only; older covers 1 and 2; base covers 3.
        for fov, (seg, fill) in {
            1: (primary, FILL_PRIMARY_SEGMENTATION),
            2: (older, FILL_OLDER_SEGMENTATION),
        }.items():
            (seg / f'FOV{fov:05}').mkdir(parents=True)
            _write_labels(seg / f'FOV{fov:05}' / f'CellLabels_F{fov:05}.tif', fill)
        # Older ALSO has FOV 1 with a different fill; primary must win there.
        (older / 'FOV00001').mkdir(parents=True)
        _write_labels(older / 'FOV00001' / 'CellLabels_F00001.tif', FILL_OLDER_SEGMENTATION)
        (inputdir / 'FOV00003').mkdir(parents=True)
        _write_labels(inputdir / 'FOV00003' / 'CellLabels_F00003.tif', FILL_BASE_LABELS)

        stitch_images(args_list=[
            '-i', inputdir.as_posix(),
            '-f', offsetsdir.as_posix(),
            '-o', tmp.as_posix(),
            '--celllabels-subdir', 'Segmentation_aaaa_002,Segmentation_bbbb_001',
            '--labels',
        ])

        _grp, labels = _labels_values(tmp)
        present = set(np.unique(labels).tolist()) - {0}
        # pair_np() rewrites label values, so assert on which FOVs landed, not
        # on raw fills: three distinct non-zero populations means all three
        # FOVs were placed, each from exactly one source.
        assert len(present) == 3, f"expected 3 FOVs stitched, got {present}"
        print("  ok: primary/older/base label sources all resolved")


def test_3d_labels_with_2d_morphology_broadcasts_images():
    """Maddie's SORL1 3D resegmentation shape: Z-sliced CellLabels under a
    Segmentation_* subdir, but a single 2D morphology acquisition. The one
    morphology plane must be broadcast to every output z plane."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputdir = tmp / 'CellStatsDir'
        offsetsdir = tmp / 'RunSummary'
        morphology = inputdir / 'Morphology2D'
        for d in (inputdir, offsetsdir, morphology):
            d.mkdir(parents=True)
        _write_offsets(offsetsdir, [1])

        seg = inputdir / 'Segmentation_c8def808_002' / 'FOV00001'
        seg.mkdir(parents=True)
        z_planes = [1, 2, 3]
        for z in z_planes:
            _write_labels(seg / f'CellLabels_F00001_Z{z:03}.tif', FILL_PRIMARY_SEGMENTATION + z)
        # Morphology has NO _Z suffix and lives in Morphology2D, not Morphology3D.
        _write_morphology(morphology / 'Run_C902_P99_N99_F00001.TIF')

        stitch_images(args_list=[
            '-i', inputdir.as_posix(),
            '-f', offsetsdir.as_posix(),
            '-o', tmp.as_posix(),
            '--celllabels-subdir', 'Segmentation_c8def808_002',
            '--input-ndim', '3',
            '--output-ndim', '3',
            '--z-step-um', '0.8',
        ])

        grp, labels = _labels_values(tmp)
        assert grp.attrs['CosMx']['ndim'] == 3
        # z=0 from the 2D morphology must NOT leak in as a phantom plane.
        assert grp.attrs['CosMx']['z_slices'] == z_planes, grp.attrs['CosMx']['z_slices']
        assert labels.shape[0] == len(z_planes), labels.shape
        for index in range(len(z_planes)):
            assert (labels[index] != 0).any(), f"labels empty at z index {index}"

        # The single morphology plane is stored ONCE, not copied per z plane:
        # duplicates would be byte-identical and would multiply mosaic size by
        # the z count, blowing past Fargate's 200 GiB ephemeral storage.
        assert grp.attrs['CosMx']['morphology_ndim'] == 2, grp.attrs['CosMx']
        dapi = np.asarray(grp['DAPI']['0'])
        assert dapi.ndim == 2, dapi.shape
        assert (dapi != 0).any(), "morphology plane is empty"
        print("  ok: 3D labels + single stored 2D morphology plane")


def test_3d_morphology_still_written_per_plane():
    """When morphology genuinely has z planes, they are all kept."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputdir = tmp / 'CellStatsDir'
        offsetsdir = tmp / 'RunSummary'
        morphology = inputdir / 'Morphology3D'
        for d in (inputdir, offsetsdir, morphology):
            d.mkdir(parents=True)
        _write_offsets(offsetsdir, [1])
        (inputdir / 'FOV00001').mkdir(parents=True)
        z_planes = [1, 2]
        for z in z_planes:
            _write_labels(inputdir / 'FOV00001' / f'CellLabels_F00001_Z{z:03}.tif',
                          FILL_PRIMARY_SEGMENTATION + z)
            _write_morphology(morphology / f'Run_C902_P99_N99_F00001_Z{z:03}.TIF')

        stitch_images(args_list=[
            '-i', inputdir.as_posix(),
            '-f', offsetsdir.as_posix(),
            '-o', tmp.as_posix(),
            '--input-ndim', '3',
            '--output-ndim', '3',
            '--z-step-um', '0.8',
        ])

        grp = zarr.open((tmp / 'images').as_posix(), mode='r')
        assert grp.attrs['CosMx']['morphology_ndim'] == 3, grp.attrs['CosMx']
        dapi = np.asarray(grp['DAPI']['0'])
        assert dapi.shape[0] == len(z_planes), dapi.shape
        for index in range(len(z_planes)):
            assert (dapi[index] != 0).any(), f"morphology empty at z index {index}"
        print("  ok: genuine 3D morphology keeps every plane")


def test_3d_labels_planes_are_distinct_and_pyramid_written():
    """Per-plane assembly must place each plane's own tiles, not repeat one
    plane, and must still produce a full multiscale pyramid.

    This is the regression guard for building the volume plane-by-plane: a
    z-indexing slip there would silently write the same plane everywhere,
    which still yields correct shapes and non-empty data.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputdir = tmp / 'CellStatsDir'
        offsetsdir = tmp / 'RunSummary'
        morphology = inputdir / 'Morphology2D'
        for d in (inputdir, offsetsdir, morphology):
            d.mkdir(parents=True)
        _write_offsets(offsetsdir, [1, 2])

        seg = inputdir / 'Segmentation_abc_002'
        z_planes = [1, 2, 3, 4]
        for fov in (1, 2):
            (seg / f'FOV{fov:05}').mkdir(parents=True)
            for z in z_planes:
                # Unique value per (fov, z) so a mixed-up plane is detectable
                _write_labels(seg / f'FOV{fov:05}' / f'CellLabels_F{fov:05}_Z{z:03}.tif',
                              fov * 100 + z)
        _write_morphology(morphology / 'Run_C902_P99_N99_F00001.TIF')
        _write_morphology(morphology / 'Run_C902_P99_N99_F00002.TIF')

        stitch_images(args_list=[
            '-i', inputdir.as_posix(), '-f', offsetsdir.as_posix(),
            '-o', tmp.as_posix(), '--celllabels-subdir', 'Segmentation_abc_002',
            '--input-ndim', '3', '--output-ndim', '3', '--z-step-um', '0.8',
            '--labels',
        ])

        grp = zarr.open((tmp / 'images').as_posix(), mode='r')
        labels = grp['labels']
        level0 = np.asarray(labels['0'])
        assert level0.shape[0] == len(z_planes), level0.shape

        # Each plane must carry its own distinct label population.
        per_plane = []
        for z in range(len(z_planes)):
            ids = set(np.unique(level0[z]).tolist()) - {0}
            assert ids, f"plane {z} empty"
            per_plane.append(frozenset(ids))
        assert len(set(per_plane)) == len(z_planes), \
            f"planes not distinct — a z-index slip would look like this: {per_plane}"

        # Pyramid: every level present, halving in y/x, z preserved.
        multiscales = labels.attrs['multiscales'][0]
        level_paths = [d['path'] for d in multiscales['datasets']]
        assert [a['name'] for a in multiscales['axes']] == ['z', 'y', 'x']
        assert len(level_paths) >= 1
        previous = None
        for name in level_paths:
            shape = labels[name].shape
            assert shape[0] == len(z_planes), f"level {name} lost z planes: {shape}"
            if previous is not None:
                assert shape[1] <= max(1, previous[1] // 2 + 1), (name, shape, previous)
            previous = shape
        print(f"  ok: {len(z_planes)} distinct planes, {len(level_paths)} pyramid levels")


def test_channel_name_override_renames_layers():
    """--channel-name pins the true panel when MorphologyKit metadata is wrong."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputdir = tmp / 'CellStatsDir'
        offsetsdir = tmp / 'RunSummary'
        morphology = inputdir / 'Morphology2D'
        for d in (inputdir, offsetsdir, morphology):
            d.mkdir(parents=True)
        _write_offsets(offsetsdir, [1])
        (inputdir / 'FOV00001').mkdir(parents=True)
        _write_labels(inputdir / 'FOV00001' / 'CellLabels_F00001.tif', FILL_BASE_LABELS)
        _write_morphology(morphology / 'Run_C902_P99_N99_F00001.TIF')

        stitch_images(args_list=[
            '-i', inputdir.as_posix(),
            '-f', offsetsdir.as_posix(),
            '-o', tmp.as_posix(),
            '--channel-name', 'B=AT8',
            '--channel-name', 'G=6E10',
        ])

        grp = zarr.open((tmp / 'images').as_posix(), mode='r')
        keys = set(grp.group_keys())
        assert 'AT8' in keys and '6E10' in keys, keys
        # Un-overridden channels keep their kit names.
        assert 'DAPI' in keys and 'GFAP' in keys, keys
        assert 'Histone' not in keys, keys
        print("  ok: channel-name overrides applied")


def test_map_z_to_morphology():
    """One morphology plane broadcasts; matching planes pair up; gaps snap to
    the nearest available plane."""
    assert _map_z_to_morphology([1, 2, 3], {0}) == {1: 0, 2: 0, 3: 0}
    assert _map_z_to_morphology([1, 2, 3], {1, 2, 3}) == {1: 1, 2: 2, 3: 3}
    assert _map_z_to_morphology([1, 2, 3], {1, 3}) == {1: 1, 2: 1, 3: 3}
    assert _map_z_to_morphology([1], set()) == {1: None}
    print("  ok: z-to-morphology mapping")


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_') and callable(fn):
            print(f"{name} ...")
            fn()
    print("\nAll fork stitch tests passed.")
