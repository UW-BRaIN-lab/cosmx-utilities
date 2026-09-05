#!/usr/bin/env python3
"""Tests for segmentation_versions.py (no S3 / network).

The fingerprint is subtle -- a build is identified by which KEYS exist, not by
their values -- so these fix that behaviour against realistic payloads modelled on
the two AtoMx profiles we actually saw (2026-04-16 and 2026-05-08).

Runnable either under pytest or directly:
    uv run python pipeline/python/tests/test_segmentation_versions.py
"""
import importlib.util
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "segmentation_versions.py"
_spec = importlib.util.spec_from_file_location("segmentation_versions", _SCRIPT)
sv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sv)

BASE_PARAMS = {
    "NuclearDiameterUm": 10.8, "CellDiameterUm": 14.4, "CellDilationUm": 0.72,
    "MinCellSizeUm": 3.5, "MaxCellSizeUm": 42.8,
    "NucleiModel": "bsbNuc", "CytoplasmModel": "bsbNeuro3Ch",
    "ForegroundThreshold": -1,
}


def april_payload():
    """The 2026-04-16 profile: no Run3DSegmentation, no assignRNAtoNearest."""
    return [{"Cellsegmentationsetconfig": "B - Brain Tissue (RNA)",
             "Parameters": dict(BASE_PARAMS)}]


def may_payload():
    """The 2026-05-08 profile: identical parameters, two extra keys."""
    params = dict(BASE_PARAMS, Run3DSegmentation=False, assignRNAtoNearest=False)
    return [{"Cellsegmentationsetconfig": "B - Brain Tissue (RNA)",
             "Parameters": params}]


def test_the_newer_build_is_identified_by_key_presence():
    facts = sv.extract_facts(may_payload())
    print(f"may -> {facts['atomx_build']} via {facts['fingerprint_keys']}")
    assert facts["atomx_build"] == sv.BUILD_WITH
    assert "Run3DSegmentation" in facts["fingerprint_keys"]


def test_the_older_build_is_identified_by_their_absence():
    facts = sv.extract_facts(april_payload())
    print(f"april -> {facts['atomx_build']}")
    assert facts["atomx_build"] == sv.BUILD_WITHOUT
    assert facts["fingerprint_keys"] == ""


def test_a_false_value_still_counts_as_present():
    """Both real profiles carry these as false; the VALUE says nothing, presence does."""
    assert sv.extract_facts(may_payload())["atomx_build"] == sv.BUILD_WITH


def test_identical_parameters_are_reported_for_both_builds():
    """The whole point: same config either side of the update."""
    april, may = sv.extract_facts(april_payload()), sv.extract_facts(may_payload())
    for key in sv.PARAMETER_KEYS:
        assert april[key] == may[key], f"{key} differs: {april[key]} vs {may[key]}"
    print(f"parameters identical across builds: {sv.PARAMETER_KEYS[:3]}...")


def test_parameters_are_found_at_any_nesting_depth():
    nested = {"a": {"b": [{"c": {"CellDiameterUm": 14.4}}]}}
    assert sv.extract_facts(nested)["CellDiameterUm"] == 14.4


def test_a_date_field_is_captured_when_present():
    payload = {"DateCreated": "2026-05-08T10:59:05", "Parameters": dict(BASE_PARAMS)}
    facts = sv.extract_facts(payload)
    assert "2026-05-08" in facts["date_field"], facts["date_field"]
    print(f"date_field={facts['date_field']}")


def test_missing_date_is_empty_not_absent():
    """Downstream joins on a stable schema, so the column must always exist."""
    assert sv.extract_facts(april_payload())["date_field"] == ""


def test_cellstats_prefix_appends_to_the_slide_scan_path():
    """CellStatsDir is a CHILD of decoded_prefix, not a sibling of DecodedFiles."""
    assert sv.cellstats_prefix("a/b/DecodedFiles/S1/20250403_213640_S3/") == \
        "a/b/DecodedFiles/S1/20250403_213640_S3/CellStatsDir"


def test_walk_scalars_keeps_the_first_value_not_the_last():
    payload = {"x": {"CellDiameterUm": 14.4}, "y": {"CellDiameterUm": 99.9}}
    assert sv.walk_scalars(payload, {})["CellDiameterUm"] == 14.4


def test_manifest_reader_strips_the_crlf_from_the_last_column():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "m.csv"
        p.write_bytes(b"slide_id,decoded_prefix\r\nS1,a/b/S1/scan/\r\n")
        frame = sv.read_manifest(p)
    assert frame["decoded_prefix"].iloc[0] == "a/b/S1/scan/", repr(frame["decoded_prefix"].iloc[0])
    print("CRLF stripped from decoded_prefix")


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
