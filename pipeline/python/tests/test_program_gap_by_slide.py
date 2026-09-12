#!/usr/bin/env python3
"""Tests for program_gap_by_slide.py — the leave-one-slide-out check (no S3 / network).

The pair that matters is the two discriminating cases: a deficit carried by ONE slide must show
up as a large move when that slide is dropped, and a deficit carried by the whole cohort must
survive every deletion. Without both, a leave-one-out that always looks reassuring would be
worthless — and the amplicon result it is used on has exactly this failure mode available, since
7495G37302G3 dominates the usable FOVs for b, t and e.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy \\
        python pipeline/python/tests/test_program_gap_by_slide.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "program_gap_by_slide.py"
# One big slide and three small ones, the shape of the real usable-FOV distribution.
SLIDES = {"HOT": 40, "S2": 10, "S3": 10, "S4": 10}
CONTROL_GAP = -0.3


def _write_table(path: Path, excess_of) -> None:
    """A per-unit table shaped like 75e's --per-unit-csv output."""
    rng = np.random.default_rng(0)
    rows = [{"unit": f"{slide}:F{i}", "slide": slide, "comparator": "OPC-like",
             "gap": rng.normal(CONTROL_GAP, 0.05) + excess_of(slide)}
            for slide, n in SLIDES.items() for i in range(n)]
    pd.DataFrame(rows).to_csv(path, index=False)


def _leave_one_out_changes(tmp: Path, excess_of) -> dict[str, float]:
    """Run the script and pull the per-slide 'change' column out of its leave-one-out block."""
    _write_table(tmp / "amp.csv", excess_of)
    _write_table(tmp / "ctl.csv", lambda slide: 0.0)
    out = subprocess.run([sys.executable, str(SCRIPT), "--amplicon", str(tmp / "amp.csv"),
                          "--control", str(tmp / "ctl.csv")], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    body = out.stdout.split("=== leave one slide out ===")[1]
    changes = {}
    for line in body.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in SLIDES:
            changes[parts[0]] = float(parts[3])
    assert set(changes) == set(SLIDES), f"did not parse every slide: {changes}"
    return changes


def test_a_one_slide_artefact_is_caught():
    """Only HOT carries the deficit, so dropping HOT must move the median excess sharply."""
    with tempfile.TemporaryDirectory() as d:
        changes = _leave_one_out_changes(Path(d), lambda s: -0.8 if s == "HOT" else 0.0)
    assert changes["HOT"] > 0.3, f"dropping the one hot slide should move the excess up: {changes}"


def test_a_cohort_wide_effect_survives_every_deletion():
    """Every slide carries it, so no single deletion may meaningfully move the excess."""
    with tempfile.TemporaryDirectory() as d:
        changes = _leave_one_out_changes(Path(d), lambda s: -0.8)
    worst = max(abs(c) for c in changes.values())
    assert worst < 0.15, f"a cohort-wide effect must be robust to deletion, got {changes}"


def test_missing_comparator_fails_loudly():
    """A typo'd comparator must not silently analyse an empty frame."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_table(tmp / "amp.csv", lambda s: -0.8)
        _write_table(tmp / "ctl.csv", lambda s: 0.0)
        out = subprocess.run([sys.executable, str(SCRIPT), "--amplicon", str(tmp / "amp.csv"),
                              "--control", str(tmp / "ctl.csv"), "--comparator", "AC-like"],
                             capture_output=True, text=True)
    assert out.returncode != 0
    assert "no rows for comparator" in (out.stdout + out.stderr)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
