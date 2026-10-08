#!/usr/bin/env python3
"""apply_fov_annotations must accept a donor-id (numeric) region.

SORL1 sets Region = donor id (741, 6803). Read as int, the composition report crashed
(`ValueError: Unknown format code 's'`) on all 20 stage-1 tasks, and a numeric Region would
have reached the batch column. Pins: no crash, Region is text, every cell is reassigned.

Run:  uv run --with pandas --with numpy --with anndata --with scipy python pipeline/python/tests/test_fov_annotations_numeric_region.py
"""
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flatfiles_to_anndata import apply_fov_annotations  # noqa: E402

SLIDE = "SLIDE_A"


def test_numeric_region_stays_text_and_reassigns_every_cell():
    obs = pd.DataFrame({"fov": [1, 1, 2, 3, 3, 3], "Region": [""] * 6})
    table = pd.DataFrame({
        "slide_id": [SLIDE] * 3, "fov": [1, 2, 3], "region": [741, 741, 6803],
        "mixed_adjacent": 0, "exclude": 0,
    })
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "ann.csv"
        table.to_csv(path, index=False)
        out = apply_fov_annotations(obs, path, SLIDE)

    print("Region values:", out["Region"].tolist(), "dtype:", out["Region"].dtype)
    assert out["Region"].tolist() == ["741", "741", "741", "6803", "6803", "6803"]
    assert all(isinstance(r, str) for r in out["Region"])
    assert out["Region_excluded"].sum() == 0


if __name__ == "__main__":
    test_numeric_region_stays_text_and_reassigns_every_cell()
    print("OK")
