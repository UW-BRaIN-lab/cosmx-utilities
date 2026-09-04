#!/usr/bin/env python3
"""Tests for morphology_channel_contrast.py (no S3 / network).

The point of these tests is not that the arithmetic runs, but that each metric
actually separates the hypothesis it is supposed to separate. Synthetic channels
are built with known ground truth -- a crisp stain, the same stain buried under a
diffuse haze, the same stain with granules in a minority of cells, and a stain that
tracks cytoplasm -- and each test asserts the metric moves the predicted way.

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_morphology_channel_contrast.py
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT = Path(__file__).resolve().parents[1] / "morphology_channel_contrast.py"
_spec = importlib.util.spec_from_file_location("morphology_channel_contrast", _SCRIPT)
mcc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mcc)

N_CELLS = 4000
DIFFUSE_HAZE = 600.0
GRANULE_FRACTION = 0.25


def synthetic_frame(seed: int = 0) -> pd.DataFrame:
    """Four channels with known, deliberately different failure modes."""
    rng = np.random.default_rng(seed)
    area = rng.lognormal(np.log(300), 0.4, N_CELLS)
    rrna = rng.lognormal(np.log(500), 0.5, N_CELLS)

    # Crisp: bright nuclear signal on a low floor, independent of cell size.
    crisp_mean = rng.lognormal(np.log(400), 0.55, N_CELLS) + 5.0
    crisp_max = crisp_mean * rng.uniform(1.5, 2.5, N_CELLS)

    # Hazy: identical signal, plus a large diffuse offset in both Mean and Max.
    hazy_mean = crisp_mean + DIFFUSE_HAZE
    hazy_max = crisp_max + DIFFUSE_HAZE

    # Granular: crisp signal, but a minority of cells carry bright puncta that lift
    # Max only. This is the lipofuscin shape.
    granular_mean = crisp_mean.copy()
    granular_max = crisp_max.copy()
    granules = rng.random(N_CELLS) < GRANULE_FRACTION
    granular_max[granules] = granular_mean[granules] * rng.uniform(
        8.0, 20.0, int(granules.sum()))

    # Cytoplasmic: intensity accumulates with cell area and tracks the rRNA channel.
    cyto_mean = 20.0 + 0.5 * area + 0.3 * rrna + rng.normal(0, 10.0, N_CELLS)
    cyto_max = cyto_mean * rng.uniform(1.4, 2.0, N_CELLS)

    return pd.DataFrame({
        "fov": rng.integers(1, 6, N_CELLS),
        "Area": area,
        "Mean.rRNA": rrna,
        "Max.rRNA": rrna * rng.uniform(1.4, 2.0, N_CELLS),
        "Mean.Crisp": crisp_mean, "Max.Crisp": crisp_max,
        "Mean.Hazy": hazy_mean, "Max.Hazy": hazy_max,
        "Mean.Granular": granular_mean, "Max.Granular": granular_max,
        "Mean.Cyto": cyto_mean, "Max.Cyto": cyto_max,
    })


def metrics_for(channel: str, frame: pd.DataFrame | None = None) -> dict:
    frame = synthetic_frame() if frame is None else frame
    return mcc.channel_metrics(frame, channel, "rRNA")


def test_diffuse_haze_collapses_the_contrast_index():
    """A diffuse offset lifts the background floor far more than the signal."""
    crisp, hazy = metrics_for("Crisp"), metrics_for("Hazy")
    print(f"contrast_index crisp={crisp['contrast_index']:.2f} "
          f"hazy={hazy['contrast_index']:.2f}")
    assert hazy["contrast_index"] < crisp["contrast_index"] / 2
    assert hazy["background_p05"] > crisp["background_p05"] * 4


def test_diffuse_haze_pushes_peakedness_down_not_up():
    """The metric that separates haze from granularity: haze drives Max/Mean to 1."""
    crisp, hazy = metrics_for("Crisp"), metrics_for("Hazy")
    print(f"peakedness crisp={crisp['peakedness']:.2f} hazy={hazy['peakedness']:.2f}")
    assert hazy["peakedness"] < crisp["peakedness"]
    assert hazy["peakedness_tail"] <= crisp["peakedness_tail"] * 1.05


def test_granules_in_a_minority_show_in_the_tail_not_the_median():
    """Why peakedness_tail exists: the median is blind to a 25% granular subset."""
    crisp, granular = metrics_for("Crisp"), metrics_for("Granular")
    print(f"peakedness median crisp={crisp['peakedness']:.2f} "
          f"granular={granular['peakedness']:.2f}; "
          f"tail crisp={crisp['peakedness_tail']:.2f} "
          f"granular={granular['peakedness_tail']:.2f}")
    assert granular["peakedness"] < crisp["peakedness"] * 1.2, \
        "median should NOT flag a minority subset"
    assert granular["peakedness_tail"] > crisp["peakedness_tail"] * 2


def test_cytoplasmic_signal_tracks_cell_area():
    """A true nuclear stain's per-cell mean should not grow with cell size."""
    crisp, cyto = metrics_for("Crisp"), metrics_for("Cyto")
    print(f"area_rho crisp={crisp['area_rho']:.3f} cyto={cyto['area_rho']:.3f}")
    assert abs(crisp["area_rho"]) < 0.1
    assert cyto["area_rho"] > 0.4


def test_rrna_coupling_survives_partialling_out_area():
    """The DAPI-binds-RNA test: coupling to rRNA beyond the shared area effect."""
    crisp, cyto = metrics_for("Crisp"), metrics_for("Cyto")
    print(f"rrna_rho_given_area crisp={crisp['rrna_rho_given_area']:.3f} "
          f"cyto={cyto['rrna_rho_given_area']:.3f}")
    assert abs(crisp["rrna_rho_given_area"]) < 0.1
    assert cyto["rrna_rho_given_area"] > 0.4


def test_partial_spearman_removes_a_pure_common_cause():
    """Two channels driven by area plus independent noise share nothing else."""
    rng = np.random.default_rng(1)
    area = pd.Series(rng.lognormal(np.log(300), 0.4, N_CELLS))
    left = area * 2.0 + rng.normal(0, 150.0, N_CELLS)
    right = area * 0.5 + rng.normal(0, 40.0, N_CELLS)
    raw = mcc.spearman(left, right)
    partial = mcc.partial_spearman(left, right, area)
    print(f"common cause: raw={raw:.3f} partial={partial:.3f}")
    assert raw > 0.5
    assert abs(partial) < 0.1


def test_partial_spearman_is_undefined_under_perfect_collinearity():
    """Unguarded, this returns a spurious 0.5 that reads as a real coupling."""
    rng = np.random.default_rng(2)
    area = pd.Series(rng.lognormal(np.log(300), 0.4, N_CELLS))
    left, right = area * 2.0, area * 0.5
    partial = mcc.partial_spearman(left, right, area)
    print(f"collinear: partial={partial}")
    assert np.isnan(partial), f"expected NaN, got {partial}"


def test_spearman_is_one_for_a_monotonic_nonlinear_relationship():
    values = pd.Series(np.arange(1.0, 501.0))
    assert mcc.spearman(values, values ** 3) > 0.999


def test_discover_channels_requires_both_mean_and_max():
    columns = ["fov", "Area", "Mean.DAPI", "Max.DAPI", "Mean.Orphan", "Max.Ghost"]
    assert mcc.discover_channels(columns) == ["DAPI"]


def test_resolve_role_detects_by_substring_and_honours_overrides():
    channels = ["DAPI", "Histone", "rRNA", "GFAP"]
    assert mcc.resolve_role(channels, mcc.NUCLEAR_ALIASES, None, "nuclear") == "DAPI"
    assert mcc.resolve_role(channels, mcc.REFERENCE_ALIASES, None, "reference") == "Histone"
    assert mcc.resolve_role(channels, mcc.RNA_ALIASES, None, "rrna") == "rRNA"
    # Kit metadata is unreliable, so an override must win over detection.
    assert mcc.resolve_role(channels, mcc.NUCLEAR_ALIASES, "GFAP", "nuclear") == "GFAP"


def test_resolve_role_rejects_an_override_that_is_not_a_channel():
    try:
        mcc.resolve_role(["DAPI"], mcc.NUCLEAR_ALIASES, "Nope", "nuclear")
    except mcc.MissingChannelError:
        return
    raise AssertionError("expected MissingChannelError for an unknown override")


def test_resolve_role_returns_none_when_nothing_matches():
    assert mcc.resolve_role(["B", "G", "Y"], mcc.NUCLEAR_ALIASES, None, "nuclear") is None


def test_contrast_ratio_table_is_reference_over_nuclear():
    frame = pd.DataFrame([
        {"slide_id": "s1", "channel": "DAPI", "contrast_index": 2.0},
        {"slide_id": "s1", "channel": "Histone", "contrast_index": 8.0},
        {"slide_id": "s2", "channel": "DAPI", "contrast_index": 5.0},
        {"slide_id": "s2", "channel": "Histone", "contrast_index": 5.0},
    ])
    table = mcc.contrast_ratio_table(frame, "DAPI", "Histone")
    assert table.loc["s1", "reference_over_nuclear"] == 4.0
    assert table.loc["s2", "reference_over_nuclear"] == 1.0
    # Sorted worst-DAPI first, so the slides driving the effect are at the top.
    assert list(table.index) == ["s1", "s2"]


def test_zero_background_yields_nan_not_infinity():
    """An empty channel must not produce an infinite contrast index."""
    frame = synthetic_frame()
    frame["Mean.Empty"] = 0.0
    frame["Max.Empty"] = 0.0
    assert np.isnan(mcc.channel_metrics(frame, "Empty", "rRNA")["contrast_index"])


def test_metadata_paths_extracts_slide_id_and_skips_strangers(tmp_path=None):
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "SlideA_metadata_file.csv.gz").touch()
        (root / "SlideB_metadata_file.csv").touch()
        (root / "SlideA_exprMat_file.csv.gz").touch()
        found = dict(mcc.metadata_paths(root, []))
        assert set(found) == {"SlideA", "SlideB"}


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
