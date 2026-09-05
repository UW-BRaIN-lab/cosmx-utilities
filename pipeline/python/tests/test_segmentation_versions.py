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


def test_the_key_fingerprint_is_recorded_but_does_not_classify():
    """It was tried and failed: the JSON stores neither key, so it cannot group."""
    april, may = sv.extract_facts(april_payload()), sv.extract_facts(may_payload())
    assert april["fingerprint_keys"] == ""
    assert "Run3DSegmentation" in may["fingerprint_keys"]
    assert "atomx_build" not in april, "build must come from the date, not the keys"
    print("fingerprint kept as a diagnostic only")


def test_the_stored_date_is_parsed_into_a_groupable_column():
    facts = sv.extract_facts({"Datecreated": "2026-05-08T17:59:05.381",
                              "Parameters": dict(BASE_PARAMS)})
    print(f"segmentation_date={facts['segmentation_date']}")
    assert facts["segmentation_date"] == "2026-05-08"


def test_a_missing_date_yields_an_empty_string_not_a_wrong_group():
    facts = sv.extract_facts(april_payload())
    assert facts["segmentation_date"] == ""


def test_date_strings_sort_correctly_against_the_cutoff():
    """Grouping compares ISO strings, so it must order as dates do."""
    assert "2026-04-17" < sv.DEFAULT_BUILD_CUTOFF <= "2026-05-08"


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


def test_active_uuid_is_read_from_a_truncated_gzip_stream():
    """The metadata is range-requested, so decompression hits a truncated stream."""
    import gzip, io

    class FakeBody:
        def __init__(self, data): self._data = data
        def read(self): return self._data

    class FakeClient:
        def __init__(self, blob): self._blob = blob
        def get_object(self, **kwargs):
            start, end = kwargs["Range"].removeprefix("bytes=").split("-")
            return {"Body": FakeBody(self._blob[int(start):int(end) + 1])}

    rows = ["fov,cell_ID,cellSegmentationSetId,Area"]
    rows += [f"{i},{i}, \"uuid-abc-123\",100" for i in range(4000)]
    blob = gzip.compress("\n".join(rows).encode())
    assert len(blob) > 200, "fixture should be a real gzip stream"

    uuid = sv.active_segmentation_uuid(FakeClient(blob), "b", "flat/", "S1")
    print(f"recovered uuid={uuid!r} from {len(blob)} gzip bytes")
    assert uuid == "uuid-abc-123", uuid


def test_active_uuid_is_empty_when_the_column_is_absent():
    import gzip

    class FakeBody:
        def __init__(self, d): self._d = d
        def read(self): return self._d

    class FakeClient:
        def __init__(self, b): self._b = b
        def get_object(self, **kwargs): return {"Body": FakeBody(self._b)}

    blob = gzip.compress(b"fov,cell_ID,Area\n1,1,100\n")
    assert sv.active_segmentation_uuid(FakeClient(blob), "b", "flat/", "S1") == ""


def test_active_uuid_survives_a_read_failure():
    class Boom:
        def get_object(self, **kwargs): raise RuntimeError("403")
    assert sv.active_segmentation_uuid(Boom(), "b", "flat/", "S1") == ""


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
