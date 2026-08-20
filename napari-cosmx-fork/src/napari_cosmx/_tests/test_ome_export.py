import os

import numpy as np
import pytest
import zarr
import napari_cosmx.utils.export_tiff as export_tiff_mod

from napari_cosmx.utils.export_tiff import (
    BatchStorage,
    _edges,
    _infer_pixelsize_from_groups,
    _normalize_compression,
    _output_axes,
    _scale_edges,
    _selected_z_indices,
    main,
    split_list_element,
    verify_ome_tiff,
)

test_data_location = os.path.join('data', 'Liver-S2')
LIBVIPS_AVAILABLE, _ = export_tiff_mod._libvips_available()
RUN_LIBVIPS_SMOKE = os.getenv("NAPARI_COSMX_RUN_LIBVIPS_SMOKE", "1").lower() not in {
    "0",
    "false",
    "no",
}

def test_split_list_element():
    assert split_list_element(["a,b,c"]) == ["a", "b", "c"]
    assert split_list_element(["a b c"]) == ["a", "b", "c"]
    assert split_list_element(["a, b,c"]) == ["a", "b", "c"]
    assert split_list_element(["a"]) == ["a"]
    with pytest.raises(TypeError):  # Correctly test for TypeError
        split_list_element()
    assert split_list_element([""]) == []
    assert split_list_element([]) == []
    assert split_list_element("a") == "Input must be a list"
    assert split_list_element(["a", "b"]) == "Input list must have length 1"
    assert split_list_element(0) == "Input must be a list"


def test_batch_storage():
    storage = BatchStorage(2)
    storage.add_item("item1")
    storage.add_item("item2")
    storage.add_item("item3")
    assert len(storage.storage) == 2
    assert storage.get_batch(0) == ['item1', 'item2']
    assert storage.get_batch(1) == ['item3']

    assert str(print(storage)) == "None"

    storage_with_labels = BatchStorage(2)
    labels = "labels_data"
    storage_with_labels.set_labels(labels)
    storage_with_labels.add_item("item1")
    storage_with_labels.add_item("item2")
    storage_with_labels.add_item("item3")
    assert storage_with_labels.get_batch(0) == ['labels_data', 'item1', 'item2']
    assert storage_with_labels.get_batch(1) == ['labels_data', 'item3']


def test_edges():
    x = np.array([[1, 1, 2], [1, 1, 1], [1,1,1]])
    edges = _edges(x)
    expected = np.array([[0, 1, 1], [0, 1, 1], [0,0,0]], dtype=np.uint8)
    np.testing.assert_array_equal(edges, expected)

def test_scaling():
    x = np.array([[1, 1, 2], [1, 1, 1], [1,1,1]])
    edges = _edges(x)
    scaled = _scale_edges(edges)
    expected = np.array(
        [[100, 10000, 10000], [100, 10000, 10000], [100, 100, 100]],
        dtype=np.uint16,
    )
    np.testing.assert_array_equal(scaled, expected)


def test_normalize_compression():
    assert _normalize_compression("zlib") == "zlib"
    assert _normalize_compression("none") is None
    assert _normalize_compression(None) is None


def test_output_axes():
    assert _output_axes(np.zeros((2, 16, 16), dtype=np.uint8)) == "CYX"
    assert _output_axes(np.zeros((2, 3, 16, 16), dtype=np.uint8)) == "CZYX"
    with pytest.raises(ValueError):
        _output_axes(np.zeros((16, 16), dtype=np.uint8))


def test_selected_z_indices_default_and_subset():
    z_values = [1, 3, 4]
    assert _selected_z_indices(z_values, None) == [0, 1, 2]
    assert _selected_z_indices(z_values, ["3", "1"]) == [1, 0]
    assert _selected_z_indices(z_values, ["3,4"]) == [1, 2]
    assert _selected_z_indices(z_values, ["3", "3", "4"]) == [1, 2]


def test_selected_z_indices_missing_plane_raises():
    with pytest.raises(ValueError, match="Requested z planes not found"):
        _selected_z_indices([1, 3, 4], ["2"])


def test_infer_pixelsize_from_groups_non_protein(tmp_path):
    root = zarr.open((tmp_path / "images").as_posix(), mode="w")
    grp = root.create_group("DNA")
    grp.attrs["multiscales"] = [{
        "datasets": [{
            "path": "0",
            "coordinateTransformations": [{"type": "scale", "scale": [0.8, 0.12, 0.12]}],
        }]
    }]
    assert _infer_pixelsize_from_groups(root) == pytest.approx(0.12)


def test_infer_pixelsize_from_groups_protein_only(tmp_path):
    root = zarr.open((tmp_path / "images").as_posix(), mode="w")
    protein = root.create_group("protein")
    pgrp = protein.create_group("CD3")
    pgrp.attrs["multiscales"] = [{
        "datasets": [{
            "path": "0",
            "coordinateTransformations": [{"type": "scale", "scale": [0.2, 0.2]}],
        }]
    }]
    assert _infer_pixelsize_from_groups(root) == pytest.approx(0.2)


@pytest.mark.skipif(not os.path.exists(test_data_location), reason="Test data not found. Skipping.") #Skip test if data is not found
def test_main_cli_and_basic_export(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["-i", "not_a_valid_zarr_path", "-o", tmp_path.as_posix(), "--filename", "test.ome.tif"])

    assert "Could not find images directory at not_a_valid_zarr_path" in str(excinfo.value)

    filename = "test.ome.tif"
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    main([
        "-i", test_data_location,
        "-o", outdir.as_posix(),
        "--filename", filename,
        "-s",
        "-c", "foo", "DNA",
    ])
    assert (outdir / f"batch_0_{filename}").exists()

    outdir_invalid = tmp_path / "out_invalid"
    outdir_invalid.mkdir(parents=True, exist_ok=True)
    main([
        "-i", test_data_location,
        "-o", outdir_invalid.as_posix(),
        "--filename", filename,
        "-s",
        "-c", "no valid IFs",
        "-v",
    ])
    assert outdir_invalid.exists()
    assert not (outdir_invalid / f"batch_0_{filename}").exists()


# ---------------------------------------------------------------------------
# Helpers for building a minimal synthetic zarr fixture
# ---------------------------------------------------------------------------

def _make_synthetic_zarr(root_path, *, n_channels=2, height=64, width=64,
                          z_slices=None, pixel_size_um=0.12, include_labels=False):
    """Write a minimal zarr store that export_tiff can consume.

    Returns the path to the store directory (images/).
    """
    import numpy as np
    import zarr

    store = zarr.open(str(root_path / "images"), mode="w")
    is_3d = z_slices is not None and len(z_slices) > 1

    channel_names = [f"CH{i}" for i in range(n_channels)]
    for ch in channel_names:
        if is_3d:
            data = (np.random.randint(0, 1000, (len(z_slices), height, width), dtype=np.uint16))
        else:
            data = (np.random.randint(0, 1000, (height, width), dtype=np.uint16))

        grp = store.create_group(ch)
        grp.create_dataset("0", data=data)
        scale = ([pixel_size_um * 0.001, pixel_size_um * 0.001] if not is_3d
                 else [0.0008, pixel_size_um * 0.001, pixel_size_um * 0.001])
        dim_names = ["y", "x"] if not is_3d else ["z", "y", "x"]
        grp.attrs["multiscales"] = [{
            "axes": [{"name": d, "type": "space", "unit": "micrometer"} for d in dim_names],
            "datasets": [{
                "path": "0",
                "coordinateTransformations": [{"type": "scale", "scale": scale}],
            }],
        }]
        lo, hi = int(data.min()), int(data.max())
        grp.attrs["omero"] = {
            "name": ch,
            "channels": [{
                "label": ch,
                "window": {"min": lo, "max": hi, "start": lo, "end": hi},
                "color": "FFFFFF",
            }],
        }

    if include_labels:
        if is_3d:
            labels = np.random.randint(0, 16, (len(z_slices), height, width), dtype=np.uint16)
        else:
            labels = np.random.randint(0, 16, (height, width), dtype=np.uint16)
        labels_grp = store.create_group("labels")
        labels_grp.create_dataset("0", data=labels)
        scale = ([pixel_size_um * 0.001, pixel_size_um * 0.001] if not is_3d
                 else [0.0008, pixel_size_um * 0.001, pixel_size_um * 0.001])
        dim_names = ["y", "x"] if not is_3d else ["z", "y", "x"]
        labels_grp.attrs["multiscales"] = [{
            "axes": [{"name": d, "type": "space", "unit": "micrometer"} for d in dim_names],
            "datasets": [{
                "path": "0",
                "coordinateTransformations": [{"type": "scale", "scale": scale}],
            }],
        }]

    fov_offsets = {"FOV": [1], "X_mm": [0.0], "Y_mm": [0.0]}
    store.attrs["CosMx"] = {
        "scale_um": pixel_size_um,
        "fov_height": height,
        "fov_width": width,
        "fov_offsets": fov_offsets,
        "ndim": 3 if is_3d else 2,
        "z_slices": list(z_slices) if z_slices else [0],
    }
    return root_path


def test_verify_ome_tiff_file_not_found():
    with pytest.raises(FileNotFoundError):
        verify_ome_tiff("/nonexistent/path/output.ome.tif")


def test_verify_ome_tiff_2d_roundtrip(tmp_path):
    """Full 2D export -> verify: axes CYX, shape, channel names, pixel size."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=2, height=64, width=64,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "test2d.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0", "CH1",
        "--levels", "1",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "OME-TIFF was not created"

    report = verify_ome_tiff(out_file, verbose=True)
    assert report["has_ome_xml"], "No OME-XML in exported file"
    assert report["issues"] == [], f"Unexpected issues: {report['issues']}"

    # Channel count should match
    assert len(report["channels"]) == 2

    # Pixel size written to OME-XML
    assert report["pixel_size_um"] == pytest.approx(0.12, rel=0.01)


def test_main_no_selectors_writes_nothing(tmp_path):
    """Without --channels/--proteins/--all, exporter should write nothing."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=2, height=32, width=32,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "default_selectors.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "--levels", "1",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert not out_file.exists(), "Expected no output when no selectors are provided"


def test_main_all_flag_exports_all_items(tmp_path):
    """--all exports every available channel/protein selection."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=2, height=32, width=32,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "all_items.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "--levels", "1",
        "--all",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "Expected --all to create an OME-TIFF"


def test_main_all_overrides_explicit_selectors_without_duplicates(tmp_path):
    """--all should not duplicate channels/proteins even if selectors are also provided."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=2, height=32, width=32,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "all_no_dupes.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "--levels", "1",
        "--all",
        "-c", "CH0",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "Expected --all export output"
    report = verify_ome_tiff(out_file, verbose=True)
    assert len(report["channels"]) == 2, (
        f"Expected each channel once with --all, got {report['channels']}"
    )


def test_main_dry_run_prints_plan_and_writes_nothing(tmp_path, capsys):
    """--dry-run should report an export plan and not write files."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=2, height=32, width=32,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "dry_run.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "--levels", "1",
        "--all",
        "--dry-run",
    ])

    captured = capsys.readouterr()
    assert "Dry run: no files will be written." in captured.out
    assert "Estimated output files:" in captured.out
    out_file = outdir / f"batch_0_{filename}"
    assert not out_file.exists(), "Dry-run should not create output files"


def test_verify_ome_tiff_3d_roundtrip(tmp_path):
    """Full 3D volumetric export (--volumetrics): axes CZYX, z step metadata present."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 2, 3]
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=z_slices, pixel_size_um=0.18)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "test3d.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0",
        "--levels", "1",
        "--volumetrics",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "OME-TIFF was not created"

    report = verify_ome_tiff(out_file, verbose=True)
    assert report["has_ome_xml"], "No OME-XML in exported file"
    assert report["issues"] == [], f"Unexpected issues: {report['issues']}"
    assert report["pixel_size_um"] == pytest.approx(0.18, rel=0.01)


def test_verify_ome_tiff_3d_volumetric_with_segmentation(tmp_path):
    """3D volumetric export with --segmentation succeeds and preserves OME metadata."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 2, 3]
    _make_synthetic_zarr(
        store_root,
        n_channels=1,
        height=32,
        width=32,
        z_slices=z_slices,
        pixel_size_um=0.18,
        include_labels=True,
    )

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "test3d_seg.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0",
        "--segmentation",
        "--levels", "1",
        "--volumetrics",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "OME-TIFF was not created"

    report = verify_ome_tiff(out_file, verbose=True)
    assert report["has_ome_xml"], "No OME-XML in exported file"
    assert report["issues"] == [], f"Unexpected issues: {report['issues']}"
    assert report["pixel_size_um"] == pytest.approx(0.18, rel=0.01)


def test_verify_ome_tiff_split_z_planes(tmp_path):
    """Default 3D export: one file per z plane, each 2D CYX, named with suffix."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 3]
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=z_slices, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "test_split.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0",
        "--levels", "1",
        # Split-z is the default behavior for 3D; no flag needed
    ])

    for z_val in z_slices:
        z_file = outdir / f"batch_0_test_split_z{int(z_val):03d}.ome.tif"
        assert z_file.exists(), f"Expected per-z file not found: {z_file.name}"
        report = verify_ome_tiff(z_file, verbose=True)
        assert report["has_ome_xml"], f"No OME-XML in {z_file.name}"
        assert report["issues"] == [], f"Issues in {z_file.name}: {report['issues']}"


def test_always_writes_tiff_resolution_tags(tmp_path):
    """TIFF resolution tags (pixels/cm) are always written for 2D exports."""
    from tifffile import TiffFile

    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=None, pixel_size_um=0.2)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "tags_test.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0",
        "--levels", "1",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "Expected OME-TIFF output file"

    with TiffFile(str(out_file)) as tif:
        page0 = tif.pages[0]
        assert "XResolution" in page0.tags
        assert "YResolution" in page0.tags
        # 3 is TIFF RESUNIT.CENTIMETER
        assert int(page0.tags["ResolutionUnit"].value) == 3


def test_verify_ome_tiff_z_plane_selection(tmp_path):
    """--z-planes with --volumetrics exports selected z planes into a CZYX file."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 2, 3]
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=z_slices, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()
    filename = "test_zsel.ome.tif"
    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", filename,
        "-c", "CH0",
        "--levels", "1",
        "--z-planes", "1", "3",
        "--volumetrics",
    ])

    out_file = outdir / f"batch_0_{filename}"
    assert out_file.exists(), "OME-TIFF was not created for z-plane selection"
    report = verify_ome_tiff(out_file, verbose=True)
    assert report["has_ome_xml"], "No OME-XML in exported file"
    assert report["issues"] == [], f"Unexpected issues: {report['issues']}"
    # Shape should reflect only the 2 selected z planes (CZYX: 1 channel, 2 z, H, W)
    if report["shape"] is not None:
        assert report["shape"][1] == 2, (
            f"Expected 2 z planes after selection, got shape {report['shape']}"
        )


@pytest.mark.skipif(not LIBVIPS_AVAILABLE, reason="pyvips/libvips is unavailable")
def test_main_3d_split_with_libvips_runs_validation(tmp_path, monkeypatch):
    """For default 3D split-z with --libvips, run per-file validation and continue when clean."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 3]
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=z_slices, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()

    calls = {"write": 0, "validate": 0}

    def _fake_write_ome_with_libvips(path, data, item_metadata, pyramid_levels,
                                     outputdir, vipshome=None, vipsconcurrency=8,
                                     resolution_px_per_cm=None, verbose=False):
        calls["write"] += 1
        # Create a placeholder so later checks on path existence are possible if needed.
        with open(path, "wb") as _:
            pass

    def _fake_verify_ome_tiff(path, **kwargs):
        calls["validate"] += 1
        return {"issues": []}

    monkeypatch.setattr(export_tiff_mod, "_write_ome_with_libvips", _fake_write_ome_with_libvips)
    monkeypatch.setattr(export_tiff_mod, "verify_ome_tiff", _fake_verify_ome_tiff)

    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", "libvips_split.ome.tif",
        "-c", "CH0",
        "--levels", "1",
        "--libvips",
    ])

    assert calls["write"] == len(z_slices)
    assert calls["validate"] == len(z_slices)


@pytest.mark.skipif(
    (not LIBVIPS_AVAILABLE) or (not RUN_LIBVIPS_SMOKE),
    reason="pyvips/libvips is unavailable or libvips smoke test disabled",
)
def test_libvips_per_z_3d_synthetic_roundtrip(tmp_path):
    """Synthetic smoke test for real libvips split-z export on 3D input.

    Dev note: this path is experimental and can be sensitive to tifffile/libvips
    version combinations across environments. 
    Set environment variable NAPARI_COSMX_RUN_LIBVIPS_SMOKE=0 to skip if needed.
    """
    store_root = tmp_path / "store"
    store_root.mkdir()
    z_slices = [1, 3]
    _make_synthetic_zarr(
        store_root,
        n_channels=1,
        height=32,
        width=32,
        z_slices=z_slices,
        pixel_size_um=0.12,
    )

    outdir = tmp_path / "out"
    outdir.mkdir()

    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", "libvips_smoke.ome.tif",
        "-c", "CH0",
        "--levels", "1",
        "--libvips",
    ])

    for z_val in z_slices:
        z_file = outdir / f"batch_0_libvips_smoke_z{int(z_val):03d}.ome.tif"
        assert z_file.exists(), f"Expected per-z output not found: {z_file.name}"
        report = verify_ome_tiff(z_file)
        assert report["has_ome_xml"], f"No OME-XML in {z_file.name}"
        assert report["issues"] == [], f"Unexpected issues in {z_file.name}: {report['issues']}"


@pytest.mark.skipif(not LIBVIPS_AVAILABLE, reason="pyvips/libvips is unavailable")
def test_main_3d_split_with_libvips_validation_failure_exits(tmp_path, monkeypatch):
    """If split-z libvips validation fails, exporter exits with rerun guidance."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=[1, 3], pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()

    def _fake_write_ome_with_libvips(path, data, item_metadata, pyramid_levels,
                                     outputdir, vipshome=None, vipsconcurrency=8,
                                     resolution_px_per_cm=None, verbose=False):
        with open(path, "wb") as _:
            pass

    def _fake_verify_ome_tiff(path, **kwargs):
        return {"issues": ["axes mismatch"]}

    monkeypatch.setattr(export_tiff_mod, "_write_ome_with_libvips", _fake_write_ome_with_libvips)
    monkeypatch.setattr(export_tiff_mod, "verify_ome_tiff", _fake_verify_ome_tiff)

    with pytest.raises(SystemExit, match="Please rerun without --libvips"):
        main([
            "-i", str(store_root),
            "-o", str(outdir),
            "--filename", "libvips_split_fail.ome.tif",
            "-c", "CH0",
            "--levels", "1",
            "--libvips",
        ])


@pytest.mark.skipif(not LIBVIPS_AVAILABLE, reason="pyvips/libvips is unavailable")
def test_main_volumetrics_disables_libvips(tmp_path, monkeypatch):
    """With --volumetrics on 3D input, libvips is disabled and tifffile path is used."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=[1, 2, 3], pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()

    calls = {"vips": 0, "pyr": 0}

    def _fake_write_ome_with_libvips(*args, **kwargs):
        calls["vips"] += 1

    def _fake_write_ome_pyramid(path, levels, pyramid_levels, item_metadata,
                                compression, resolution_px_per_cm=None, verbose=False):
        calls["pyr"] += 1
        with open(path, "wb") as _:
            pass

    monkeypatch.setattr(export_tiff_mod, "_write_ome_with_libvips", _fake_write_ome_with_libvips)
    monkeypatch.setattr(export_tiff_mod, "_write_ome_pyramid", _fake_write_ome_pyramid)

    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", "volumetric_no_vips.ome.tif",
        "-c", "CH0",
        "--levels", "1",
        "--volumetrics",
        "--libvips",
    ])

    assert calls["vips"] == 0
    assert calls["pyr"] == 1


def test_main_disables_libvips_when_unavailable(tmp_path, monkeypatch, capsys):
    """When --libvips is requested but unavailable, exporter falls back to tifffile."""
    store_root = tmp_path / "store"
    store_root.mkdir()
    _make_synthetic_zarr(store_root, n_channels=1, height=32, width=32,
                          z_slices=None, pixel_size_um=0.12)

    outdir = tmp_path / "out"
    outdir.mkdir()

    calls = {"vips": 0, "pyr": 0}

    def _fake_libvips_available(vipshome=None):
        return False, "No module named 'pyvips'"

    def _fake_write_ome_with_libvips(*args, **kwargs):
        calls["vips"] += 1

    def _fake_write_ome_pyramid(path, levels, pyramid_levels, item_metadata,
                                compression, resolution_px_per_cm=None, verbose=False):
        calls["pyr"] += 1
        with open(path, "wb") as _:
            pass

    monkeypatch.setattr(export_tiff_mod, "_libvips_available", _fake_libvips_available)
    monkeypatch.setattr(export_tiff_mod, "_write_ome_with_libvips", _fake_write_ome_with_libvips)
    monkeypatch.setattr(export_tiff_mod, "_write_ome_pyramid", _fake_write_ome_pyramid)

    main([
        "-i", str(store_root),
        "-o", str(outdir),
        "--filename", "fallback_no_vips.ome.tif",
        "-c", "CH0",
        "--levels", "1",
        "--libvips",
    ])

    assert calls["vips"] == 0
    assert calls["pyr"] == 1
    assert "--libvips requested but pyvips/libvips is unavailable" in capsys.readouterr().out



