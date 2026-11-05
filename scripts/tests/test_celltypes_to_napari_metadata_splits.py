#!/usr/bin/env python3
"""Tests for --split-spec: breaking part of a call out as its own Napari value.

The load-bearing property is that every failure is LOUD. A split exists to make a
population visible in Napari; one that quietly moves no cell leaves that population
indistinguishable from the cluster around it, and the output looks perfectly normal —
the reviewer simply never sees the thing they were asked to look at.

Runnable under pytest or directly:
    uv run python scripts/tests/test_celltypes_to_napari_metadata_splits.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "celltypes_to_napari_metadata", _SCRIPTS / "celltypes-to-napari-metadata.py")
ctn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ctn)


def make_obs() -> pd.DataFrame:
    return pd.DataFrame(
        {"leiden": ["22", "22", "5", "5", "22"]},
        index=[f"c{i}" for i in range(5)])


def write_spec(path: Path, **row) -> Path:
    pd.DataFrame([row]).to_csv(path, index=False)
    return path


def test_moves_only_the_intersection(tmp: Path) -> None:
    obs = make_obs()
    calls = pd.Series(["a", "a", "a", "b", "b"], index=obs.index)
    spec = write_spec(tmp / "ok.csv", label="a", column="leiden", values="22",
                      new_name="a_deep_retinal")
    out = ctn.apply_splits(calls, obs, ctn.load_split_spec(spec))
    # c0, c1 are 'a' AND leiden 22. c2 is 'a' but leiden 5; c4 is leiden 22 but 'b'.
    assert list(out) == ["a_deep_retinal", "a_deep_retinal", "a", "b", "b"], list(out)
    print("  splits the (label x component) intersection and nothing else")


def test_noop_split_is_fatal(tmp: Path) -> None:
    obs = make_obs()
    calls = pd.Series(["a", "a", "a", "b", "a"], index=obs.index)   # no 'b' in leiden 22
    spec = write_spec(tmp / "noop.csv", label="b", column="leiden", values="22",
                      new_name="zz")
    try:
        ctn.apply_splits(calls, obs, ctn.load_split_spec(spec))
    except SystemExit as exc:
        assert "moves no cell" in str(exc), exc
    else:
        raise AssertionError("a split that moves nothing must exit, not pass silently")
    print("  a split that moves no cell exits")


def test_unknown_inputs_are_fatal(tmp: Path) -> None:
    obs = make_obs()
    calls = pd.Series(["a"] * 5, index=obs.index)

    for name, row, expected in [
        ("value", dict(label="a", column="leiden", values="999", new_name="zz"),
         "absent from obs"),
        ("column", dict(label="a", column="not_a_column", values="22", new_name="zz"),
         "does not have"),
        ("label", dict(label="zzz", column="leiden", values="22", new_name="zz"),
         "matches no call"),
    ]:
        spec = write_spec(tmp / f"bad_{name}.csv", **row)
        try:
            ctn.apply_splits(calls, obs, ctn.load_split_spec(spec))
        except SystemExit as exc:
            assert expected in str(exc), (name, str(exc))
        else:
            raise AssertionError(f"an unknown {name} must exit")
    print("  unknown value / column / label each exit with a specific message")


def test_spec_validation(tmp: Path) -> None:
    # A new_name equal to a label being split would make the split unreadable.
    spec = write_spec(tmp / "clash.csv", label="a", column="leiden", values="22",
                      new_name="a")
    try:
        ctn.load_split_spec(spec)
    except SystemExit as exc:
        assert "collide" in str(exc), exc
    else:
        raise AssertionError("new_name colliding with a split label must exit")

    pd.DataFrame([{"label": "a", "column": "leiden", "values": "22"}]).to_csv(
        tmp / "short.csv", index=False)
    try:
        ctn.load_split_spec(tmp / "short.csv")
    except SystemExit as exc:
        assert "new_name" in str(exc), exc
    else:
        raise AssertionError("a spec missing a column must exit")
    print("  spec validation catches name collisions and missing columns")


def test_committed_retina_spec_is_valid() -> None:
    """The spec we actually ship must parse — a typo in it fails a 12-slide run."""
    committed = _SCRIPTS.parent / "pipeline" / "reference" / "napari_splits_retina.csv"
    if not committed.exists():
        print("  (no committed retina spec; skipped)")
        return
    spec = ctn.load_split_spec(committed)
    assert not spec.empty
    print(f"  committed retina spec parses: "
          f"{', '.join(spec['label'] + ' -> ' + spec['new_name'])}")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("test_celltypes_to_napari_metadata_splits")
        test_moves_only_the_intersection(tmp)
        test_noop_split_is_fatal(tmp)
        test_unknown_inputs_are_fatal(tmp)
        test_spec_validation(tmp)
        test_committed_retina_spec_is_valid()
    print("All tests passed.")


if __name__ == "__main__":
    main()
