#!/usr/bin/env python3
"""Tests for select-napari-metadata-columns.py (no S3 / network).

What matters: the tooltip ends up holding exactly the columns asked for, in that order; every
kept column keeps its `_color` partner (a column without one is drawn as one flat colour, which
looks like a result); row order and cell ids survive, since napari-cosmx joins positionally on
nothing but trusts `cell_ID`; and a blank label stays blank rather than becoming "nan", which
would appear in the tooltip as a cell type.

Run:  uv run --with pandas python scripts/tests/test_select_napari_metadata_columns.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "select-napari-metadata-columns.py"

VALUE_COLUMNS = ("celltype_refit", "Region", "leiden", "insitutree_75", "denovo_split")


def _fixture(tmp: Path, uncoloured: tuple[str, ...] = ()) -> Path:
    """Two slides with the full column set; `uncoloured` names columns written WITHOUT a _color."""
    metadata_dir = tmp / "full"
    metadata_dir.mkdir()
    for slide, cell_count in (("SLIDEA", 4), ("SLIDEB", 3)):
        table = {"cell_ID": [f"c_1_{i}_{i}" for i in range(cell_count)]}
        for position, column in enumerate(VALUE_COLUMNS):
            # One deliberate blank per column: an unlabelled cell must stay unlabelled.
            table[column] = [f"{column}_v{i}" if i else "" for i in range(cell_count)]
            if column not in uncoloured:
                table[f"{column}_color"] = [f"#00000{position}"] * cell_count
        pd.DataFrame(table).to_csv(metadata_dir / f"{slide}_metadata.csv", index=False)
    return metadata_dir


def _run(metadata_dir: Path, out_dir: Path, columns: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--metadata-dir", str(metadata_dir),
         "--out-dir", str(out_dir), "--columns", columns],
        capture_output=True, text=True)


def test_keeps_requested_columns_in_order() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        out_dir = tmp / "review"
        result = _run(_fixture(tmp), out_dir, "insitutree_75,denovo_split")
        assert result.returncode == 0, result.stderr

        written = sorted(p.name for p in out_dir.glob("*_metadata.csv"))
        assert written == ["SLIDEA_metadata.csv", "SLIDEB_metadata.csv"], written

        out = pd.read_csv(out_dir / "SLIDEA_metadata.csv", dtype=str, keep_default_na=False)
        assert list(out.columns) == [
            "cell_ID",
            "insitutree_75", "insitutree_75_color",
            "denovo_split", "denovo_split_color",
        ], list(out.columns)
        print(f"  columns: {list(out.columns)}")


def test_drops_the_rest_and_keeps_rows() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        metadata_dir, out_dir = _fixture(tmp), tmp / "review"
        assert _run(metadata_dir, out_dir, "Region").returncode == 0

        source = pd.read_csv(metadata_dir / "SLIDEB_metadata.csv", dtype=str, keep_default_na=False)
        out = pd.read_csv(out_dir / "SLIDEB_metadata.csv", dtype=str, keep_default_na=False)

        assert len(out) == len(source) == 3, (len(out), len(source))
        assert out["cell_ID"].tolist() == source["cell_ID"].tolist()
        assert out["Region"].tolist() == source["Region"].tolist()
        for gone in ("celltype_refit", "leiden", "insitutree_75", "denovo_split"):
            assert gone not in out.columns, gone
            assert f"{gone}_color" not in out.columns, gone
        print(f"  dropped 4 value columns and their colours; {len(out)} rows intact")


def test_blank_label_stays_blank() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        out_dir = tmp / "review"
        assert _run(_fixture(tmp), out_dir, "insitutree_75").returncode == 0

        raw = (out_dir / "SLIDEA_metadata.csv").read_text().splitlines()
        first_cell = raw[1].split(",")
        assert first_cell[1] == "", f"blank label became {first_cell[1]!r}"
        print(f"  blank row written as {raw[1]!r}")


def test_unknown_column_is_fatal() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        result = _run(_fixture(tmp), tmp / "review", "insitutree_66")
        assert result.returncode != 0, "a column that is not in the file must not pass silently"
        assert "insitutree_66" in result.stderr, result.stderr
        # The message has to say what IS available, or the caller cannot fix the typo.
        assert "insitutree_75" in result.stderr, result.stderr
        print(f"  {result.stderr.strip().splitlines()[0]}")


def test_column_without_color_partner_is_fatal() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        metadata_dir = _fixture(tmp, uncoloured=("leiden",))
        result = _run(metadata_dir, tmp / "review", "leiden")
        assert result.returncode != 0, "a column with no _color would render as one flat colour"
        assert "leiden" in result.stderr, result.stderr
        print(f"  {result.stderr.strip().splitlines()[0]}")


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        print(f"{test.__name__}:")
        test()
    print(f"\nAll {len(tests)} tests passed")
