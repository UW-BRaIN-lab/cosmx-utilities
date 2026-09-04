#!/usr/bin/env python3
"""Tests for morphology_localization.py (no S3 / network).

Synthetic FOVs are built with known ground truth -- a crisp nuclear stain, a
uniform wash, a cytoplasm-seeking stain, and a clipped one -- and each test asserts
the metric separates the case it exists to separate.

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_morphology_localization.py
"""
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import tifffile

_SCRIPT = Path(__file__).resolve().parents[1] / "morphology_localization.py"
_spec = importlib.util.spec_from_file_location("morphology_localization", _SCRIPT)
ml = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ml)

SIZE = 200
BACKGROUND_LEVEL = 100
NUCLEAR_LEVEL = 2000
CEILING = 65535


def synthetic_labels() -> np.ndarray:
    """0 = background, 1 = nuclear, 2 = cytoplasm, in concentric bands."""
    labels = np.zeros((SIZE, SIZE), dtype=np.uint16)
    labels[40:160, 40:160] = 2          # cytoplasm
    labels[70:130, 70:130] = 1          # nuclear inside it
    return labels


def synthetic_planes(labels: np.ndarray) -> np.ndarray:
    """Five channels, in CHANNEL_ORDER, with deliberately different behaviour."""
    nuclear, cytoplasm, background = labels == 1, labels == 2, labels == 0

    crisp = np.full(labels.shape, BACKGROUND_LEVEL, dtype=np.uint16)
    crisp[nuclear] = NUCLEAR_LEVEL
    crisp[cytoplasm] = BACKGROUND_LEVEL * 2

    wash = np.full(labels.shape, NUCLEAR_LEVEL, dtype=np.uint16)  # uniform everywhere

    cytoplasmic = np.full(labels.shape, BACKGROUND_LEVEL, dtype=np.uint16)
    cytoplasmic[cytoplasm] = NUCLEAR_LEVEL
    cytoplasmic[nuclear] = BACKGROUND_LEVEL * 3

    clipped = np.full(labels.shape, BACKGROUND_LEVEL, dtype=np.uint16)
    clipped[nuclear | cytoplasm] = CEILING

    dim = np.full(labels.shape, BACKGROUND_LEVEL, dtype=np.uint16)
    return np.stack([crisp, wash, cytoplasmic, clipped, dim])


def write_fov(root: Path, fov: int = 1) -> None:
    labels = synthetic_labels()
    tifffile.imwrite(root / f"CompartmentLabels_F{fov:05d}.tif", labels)
    tifffile.imwrite(root / f"20250101_000000_S1_C902_P99_N99_F{fov:05d}.TIF",
                     synthetic_planes(labels))


def rows_by_marker() -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_fov(root)
        pairs = ml.fov_pairs(root)
        assert len(pairs) == 1, f"expected one FOV pair, got {pairs}"
        fov, morphology, compartment = pairs[0]
        overrides = dict(zip(ml.CHANNEL_ORDER,
                             ["crisp", "wash", "cytoplasmic", "clipped", "dim"]))
        rows = ml.analyse_fov("SlideX", fov, morphology, compartment, overrides)
    return {r["marker"]: r for r in rows}


def test_crisp_nuclear_stain_scores_high_localization():
    row = rows_by_marker()["crisp"]
    print(f"crisp: localization={row['localization_index']:.1f} "
          f"nuc/cyto={row['nuclear_over_cytoplasm']:.1f}")
    assert row["localization_index"] > 10
    assert row["nuclear_over_cytoplasm"] > 5


def test_uniform_wash_collapses_localization_to_one():
    """The watercolor case: signal everywhere, marking nothing."""
    row = rows_by_marker()["wash"]
    print(f"wash: localization={row['localization_index']:.2f} "
          f"background_share={row['background_share']:.2f}")
    assert abs(row["localization_index"] - 1.0) < 0.01
    assert abs(row["background_share"] - 1.0) < 0.01


def test_cytoplasmic_stain_is_caught_by_nuclear_over_cytoplasm():
    """Localization alone cannot tell 'in cells' from 'in nuclei'; this can."""
    row = rows_by_marker()["cytoplasmic"]
    print(f"cytoplasmic: localization={row['localization_index']:.1f} "
          f"nuc/cyto={row['nuclear_over_cytoplasm']:.2f}")
    assert row["localization_index"] > 1, "it is still enriched inside cells"
    assert row["nuclear_over_cytoplasm"] < 0.5, "but not in nuclei"


def test_clipping_is_reported_so_it_is_not_mistaken_for_chemistry():
    rows = rows_by_marker()
    print(f"clipped saturated_frac={rows['clipped']['saturated_frac']:.3f} "
          f"crisp={rows['crisp']['saturated_frac']:.3f}")
    assert rows["clipped"]["saturated_frac"] > 0.2
    assert rows["crisp"]["saturated_frac"] == 0.0


def test_compartment_masks_name_the_conventional_codes():
    masks = ml.compartment_masks(synthetic_labels())
    assert set(masks) == {"background", "nuclear", "cytoplasm"}


def test_compartment_masks_keep_unknown_codes_rather_than_guessing():
    labels = synthetic_labels()
    labels[:60, :60] = 7
    masks = ml.compartment_masks(labels)
    assert "code_7" in masks, f"unexpected encoding dropped: {sorted(masks)}"


def test_compartment_masks_ignore_specks_below_the_pixel_floor():
    labels = synthetic_labels()
    labels[0, 0] = 9  # one pixel
    assert "code_9" not in ml.compartment_masks(labels)


def test_channel_override_wins_over_kit_metadata():
    """Kit metadata is known to mislabel targets, so the override must win."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_fov(root)
        _, morphology, _ = ml.fov_pairs(root)[0]
        with tifffile.TiffFile(morphology) as tif:
            markers = ml.channel_markers(tif, {"U": "PinnedDAPI"})
    assert markers[ml.CHANNEL_ORDER.index("U")] == "PinnedDAPI"
    assert markers[ml.CHANNEL_ORDER.index("B")] == ml.DEFAULT_MARKERS[0]


def test_fov_without_compartment_labels_is_skipped_not_fatal():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_fov(root, fov=1)
        labels = synthetic_labels()
        tifffile.imwrite(root / "20250101_000000_S1_C902_P99_N99_F00002.TIF",
                         synthetic_planes(labels))
        pairs = ml.fov_pairs(root)
    assert [p[0] for p in pairs] == ["1"]


def test_shape_mismatch_raises_rather_than_comparing_garbage():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_fov(root)
        # Rewrite the labels at half size, as a binning mismatch would.
        tifffile.imwrite(root / "CompartmentLabels_F00001.tif",
                         synthetic_labels()[:100, :100])
        fov, morphology, compartment = ml.fov_pairs(root)[0]
        try:
            ml.analyse_fov("SlideX", fov, morphology, compartment, {})
        except ml.ChannelCountError:
            return
    raise AssertionError("expected ChannelCountError on a shape mismatch")


if __name__ == "__main__":
    failures = 0
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            try:
                func()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    print(f"\n{failures} failure(s)")
    raise SystemExit(1 if failures else 0)
