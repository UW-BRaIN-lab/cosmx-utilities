#!/usr/bin/env python3
"""Tests for segmentation_homogeneity.py (no S3 / network).

Synthetic cohorts with a KNOWN batch effect: one group of slides segmented with
larger cell boundaries than the other. Each test asserts the tool separates that
from an evenly-segmented cohort, since a false "all fine" here would let a real
batch effect through.

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_segmentation_homogeneity.py
"""
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT = Path(__file__).resolve().parents[1] / "segmentation_homogeneity.py"
_spec = importlib.util.spec_from_file_location("segmentation_homogeneity", _SCRIPT)
sh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sh)

N_CELLS = 3000
BASE_AREA = 300.0


def write_slide(root: Path, slide_id: str, area_scale: float = 1.0,
                version: str = "1.0", seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    area = rng.lognormal(np.log(BASE_AREA * area_scale), 0.35, N_CELLS)
    pd.DataFrame({
        "fov": rng.integers(1, 21, N_CELLS),
        "cell_ID": np.arange(N_CELLS),
        "Area": area,
        "AspectRatio": rng.uniform(0.6, 1.8, N_CELLS),
        "Width": np.sqrt(area) * rng.uniform(0.9, 1.1, N_CELLS),
        "Height": np.sqrt(area) * rng.uniform(0.9, 1.1, N_CELLS),
        "version": version,
        "cellSegmentationSetId": f"uuid-{slide_id}",
        "Run_name": "run",
    }).to_csv(root / f"{slide_id}_metadata_file.csv.gz", index=False)


def cohort(tmp: Path, split: bool) -> pd.DataFrame:
    """split=True gives 4 slides at 1.0x area and 4 at 1.35x, tagged by version."""
    for i in range(4):
        write_slide(tmp, f"early{i}", 1.0, "1.0", seed=i)
    for i in range(4):
        write_slide(tmp, f"late{i}", 1.35 if split else 1.0,
                    "2.0" if split else "1.0", seed=10 + i)
    rows = [sh.slide_summary(s, p) for s, p in sh.metadata_paths(tmp, [])]
    return pd.DataFrame([r for r in rows if r])


def test_a_boundary_shift_shows_up_as_a_group_difference():
    """The batch effect we care about: same config, different cell sizes."""
    with tempfile.TemporaryDirectory() as tmp:
        frame = cohort(Path(tmp), split=True)
    grouping = sh.choose_grouping(frame, None)
    assert grouping == "version", f"expected version to be chosen, got {grouping}"
    summary = sh.report_groups(frame, grouping)
    early, late = summary.loc["1.0", "Area_p50"], summary.loc["2.0", "Area_p50"]
    print(f"Area_p50 early={early:.1f} late={late:.1f} ratio={late / early:.2f}")
    assert late / early > 1.25


def test_an_evenly_segmented_cohort_is_not_flagged():
    """A false alarm here would send someone chasing a batch effect that isn't there."""
    with tempfile.TemporaryDirectory() as tmp:
        frame = cohort(Path(tmp), split=False)
    assert sh.choose_grouping(frame, None) is None, "nothing should separate them"
    assert len(sh.report_outliers(frame)) == 0
    print(f"even cohort: {len(frame)} slides, no grouping, no outliers")


def test_outliers_are_named_when_one_slide_drifts():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for i in range(6):
            write_slide(root, f"normal{i}", 1.0, "1.0", seed=i)
        write_slide(root, "drifted", 1.6, "1.0", seed=99)
        rows = [sh.slide_summary(s, p) for s, p in sh.metadata_paths(root, [])]
        frame = pd.DataFrame([r for r in rows if r])
    flagged = sh.report_outliers(frame)
    assert "drifted" in set(flagged["slide_id"]), flagged
    print(f"flagged: {sorted(set(flagged['slide_id']))}")


def test_area_threshold_shows_a_fixed_cut_is_not_neutral():
    """Large-celled slides lose several percent to a fixed cut; narrow ones lose none."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_slide(root, "narrow", 1.0, seed=1)
        write_slide(root, "wide", 2.5, seed=2)
        cut = BASE_AREA * 3
        rows = {s: sh.slide_summary(s, p, area_threshold=cut)
                for s, p in sh.metadata_paths(root, [])}
    narrow = rows["narrow"]["pct_above_threshold"]
    wide = rows["wide"]["pct_above_threshold"]
    print(f"cut={cut:.0f}: narrow loses {narrow}%, wide loses {wide}%")
    assert wide > narrow * 3, "the cut must bite the large-celled slide much harder"
    assert narrow < 5.0


def test_area_threshold_is_absent_unless_requested():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_slide(root, "s1", 1.0)
        row = sh.slide_summary("s1", root / "s1_metadata_file.csv.gz")
    assert "pct_above_threshold" not in row


def test_choose_grouping_ignores_a_column_unique_to_every_slide():
    """cellSegmentationSetId is per-slide, so it groups nothing and must be skipped."""
    frame = pd.DataFrame({
        "slide_id": ["a", "b", "c"],
        "version": ["1.0", "1.0", "1.0"],
        "cellSegmentationSetId": ["u1", "u2", "u3"],
    })
    assert sh.choose_grouping(frame, None) is None


def test_choose_grouping_honours_an_external_request():
    frame = pd.DataFrame({"slide_id": ["a", "b"], "seg_date": ["2026-04-16", "2026-05-08"]})
    assert sh.choose_grouping(frame, "seg_date") == "seg_date"
    try:
        sh.choose_grouping(frame, "nope")
    except sh.MissingColumnError:
        return
    raise AssertionError("expected MissingColumnError for an unknown --group-by")


def test_summary_reports_provenance_and_density():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        write_slide(root, "s1", 1.0, "1.0")
        row = sh.slide_summary("s1", root / "s1_metadata_file.csv.gz")
    assert row["version"] == "1.0"
    assert row["cellSegmentationSetId"] == "uuid-s1"
    assert row["n_cells"] == N_CELLS
    assert 0 < row["cells_per_fov"] <= N_CELLS
    print(f"cells_per_fov={row['cells_per_fov']} Area_p50={row['Area_p50']}")


def test_missing_area_column_raises_rather_than_returning_nothing():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "x_metadata_file.csv.gz"
        pd.DataFrame({"fov": [1, 2, 3]}).to_csv(p, index=False)
        try:
            sh.slide_summary("x", p)
        except sh.MissingColumnError:
            return
    raise AssertionError("expected MissingColumnError when Area is absent")


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
