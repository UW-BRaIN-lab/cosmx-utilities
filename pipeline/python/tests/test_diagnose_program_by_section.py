#!/usr/bin/env python3
"""Tests for section derivation in diagnose_program_by_section.py (no S3 / network).

The one that matters is test_section_grouping_is_not_fooled_by_a_hot_section: two tissue
sections are mounted per CosMx slide, so grouping by slide pools two pieces. If one piece is
hot (ischaemia / slow fixation) and the letter's cells are concentrated in it, slide-level
grouping reports a confident CELL STATE that is not there.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy --with scipy \\
        python pipeline/python/tests/test_diagnose_program_by_section.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location(
    "diagnose_program_by_section", _PY_DIR / "diagnose_program_by_section.py")
dp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dp)


def test_sections_split_a_slide_at_the_fov_gap():
    """Two pieces per slide show up as two contiguous FOV runs (e.g. 1-135 and 181-205)."""
    ids = pd.Index([f"SL01A_F{f}_C{i}" for i, f in enumerate([1, 2, 3, 181, 182, 183])])
    got = dp.sections_from_cell_ids(ids, fov_gap=10)
    assert got.nunique() == 2
    assert got.iloc[0] == "SL01A:F1-3"
    assert got.iloc[-1] == "SL01A:F181-183"


def test_sections_keep_a_contiguous_run_together():
    """A slide with one continuous piece must not be chopped up by small FOV steps."""
    ids = pd.Index([f"SL02B_F{f}_C{i}" for i, f in enumerate([1, 2, 3, 4, 5])])
    assert dp.sections_from_cell_ids(ids, fov_gap=10).nunique() == 1


def test_sections_never_merge_across_slides():
    """Two slides sharing FOV numbers are still two different pieces of tissue."""
    ids = pd.Index(["SL01A_F1_C1", "SL02B_F1_C2"])
    got = dp.sections_from_cell_ids(ids, fov_gap=10)
    assert got.nunique() == 2


def _fixture_hot_section(rng):
    """One section per slide is hot; the letter is 80% concentrated in that hot section."""
    rows = []
    for s in range(10):
        slide = f"SL{s:02d}A{s:02d}B"
        for fov_lo, hot in ((1, True), (181, False)):
            n_letter = 320 if hot else 80
            for is_letter, n in ((True, n_letter), (False, 600)):
                score = rng.normal(8.0 if hot else 1.0, 0.5, size=n)
                rows.append(pd.DataFrame({
                    "score": score,
                    "t": is_letter,
                    "cell_id": [f"{slide}_F{fov_lo + (i % 100)}_C{rng.integers(1e9)}"
                                for i in range(n)],
                }))
    return pd.concat(rows, ignore_index=True)


def _verdict(df, unit):
    g = df.groupby([unit, "t"])["score"].mean().unstack()
    gap = (g[True] - g[False]).dropna()
    return float((gap > 0).mean()), float(gap.median())


def test_section_grouping_is_not_fooled_by_a_hot_section():
    """The regression this file exists for. Ground truth is HANDLING, not a cell state."""
    df = _fixture_hot_section(np.random.default_rng(3))
    df["section"] = dp.sections_from_cell_ids(pd.Index(df["cell_id"]), 10).to_numpy()
    df["slide"] = pd.Series(df["section"]).str.split(":").str[0]

    slide_pos, slide_gap = _verdict(df, "slide")
    sec_pos, sec_gap = _verdict(df, "section")

    # Grouping by slide pools the hot and cold pieces and manufactures a cell-state signal.
    assert slide_pos == 1.0, "slide-level should (wrongly) look unanimous"
    assert slide_gap > 1.0, "slide-level should (wrongly) show a large gap"
    # Grouping by section compares like with like, and the signal disappears.
    assert abs(sec_gap) < 0.2, f"section-level gap should collapse, got {sec_gap:+.3f}"
    assert 0.2 < sec_pos < 0.8, f"section-level should be a coin toss, got {sec_pos:.2f}"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
