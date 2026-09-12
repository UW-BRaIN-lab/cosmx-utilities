#!/usr/bin/env python3
"""Tests for program_gap_by_group.py — the leave-one-slide-out check (no S3 / network).

The pair that matters is the two discriminating cases: a deficit carried by ONE slide must show
up as a large move when that slide is dropped, and a deficit carried by the whole cohort must
survive every deletion. Without both, a leave-one-out that always looks reassuring would be
worthless — and the amplicon result it is used on has exactly this failure mode available, since
7495G37302G3 dominates the usable FOVs for b, t and e.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy \\
        python pipeline/python/tests/test_program_gap_by_group.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "program_gap_by_group.py"
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
                          "--control", str(tmp / "ctl.csv"), "--group-by", "slide"],
                         capture_output=True, text=True)
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
                              "--control", str(tmp / "ctl.csv"), "--comparator", "AC-like",
                              "--group-by", "slide"], capture_output=True, text=True)
    assert out.returncode != 0
    assert "no rows for comparator" in (out.stdout + out.stderr)


# --- the annotation join: one slide, two donors -------------------------------------------
# This is the case the whole rename is for. A slide carries two cases with FOV numbering running
# continuously across both, so slide-level grouping pools them; only the annotation reference can
# split them. If one case carries the deficit and the other does not, grouping by slide hides it.
TWO_DONOR_SLIDE = "7495G37302G3"
CANONICAL = "7495 G3 7302 G3"


def _write_annotation_fixture(tmp: Path) -> tuple[Path, Path]:
    """FOVs 1-20 are case 7495 (Tumor bulk), 21-40 are case 7302 (Contralateral)."""
    rows = []
    for fov in range(1, 41):
        first = fov <= 20
        rows.append({"Slide": CANONICAL, "Case": 7495 if first else 7302,
                     "Block": "G3", "Region": "Tumor bulk" if first else
                     "Contralateral uninvolved", "FOVs": fov})
    ann = tmp / "annotations.csv"
    # The real file carries a UTF-8 BOM; the loader must survive it.
    pd.DataFrame(rows).to_csv(ann, index=False, encoding="utf-8-sig")
    xw = tmp / "crosswalk.csv"
    pd.DataFrame([{"Canonical_slide_name": CANONICAL,
                   "AtoMx_flatfile_folder": TWO_DONOR_SLIDE}]).to_csv(xw, index=False)
    return ann, xw


def _write_two_donor_table(path: Path, excess_of) -> None:
    rng = np.random.default_rng(0)
    rows = [{"unit": f"{TWO_DONOR_SLIDE}:F{fov}", "slide": TWO_DONOR_SLIDE,
             "comparator": "OPC-like", "gap": rng.normal(CONTROL_GAP, 0.05) + excess_of(fov)}
            for fov in range(1, 41)]
    pd.DataFrame(rows).to_csv(path, index=False)


def _run_two_donor(tmp: Path, excess_of, *extra: str) -> str:
    ann, xw = _write_annotation_fixture(tmp)
    _write_two_donor_table(tmp / "amp.csv", excess_of)
    _write_two_donor_table(tmp / "ctl.csv", lambda fov: 0.0)
    out = subprocess.run([sys.executable, str(SCRIPT), "--amplicon", str(tmp / "amp.csv"),
                          "--control", str(tmp / "ctl.csv"), "--annotations", str(ann),
                          "--crosswalk", str(xw), *extra], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout


def test_annotations_split_one_slide_into_its_two_donors():
    """Only case 7495 carries the deficit; grouping by case must separate the two."""
    with tempfile.TemporaryDirectory() as d:
        txt = _run_two_donor(Path(d), lambda fov: -0.8 if fov <= 20 else 0.0)
    block = txt.split("=== per case")[1].split("=== leave one")[0]
    rows = {p[0]: float(p[2]) for p in (l.split() for l in block.splitlines()) if len(p) >= 5
            and p[0] in {"7495", "7302"}}
    assert set(rows) == {"7495", "7302"}, f"both cases must appear separately: {block}"
    assert rows["7495"] < -0.5, rows
    assert abs(rows["7302"]) < 0.15, rows


def test_region_breakdown_is_reported():
    """A deficit confined to contralateral tissue would mean something else entirely."""
    with tempfile.TemporaryDirectory() as d:
        txt = _run_two_donor(Path(d), lambda fov: -0.8 if fov <= 20 else 0.0)
    assert "by region" in txt
    assert "Tumor bulk" in txt and "Contralateral uninvolved" in txt


def test_region_filter_drops_the_other_region():
    with tempfile.TemporaryDirectory() as d:
        txt = _run_two_donor(Path(d), lambda fov: -0.8 if fov <= 20 else 0.0,
                             "--regions", "Tumor bulk")
    assert "Kept 20 FOVs in region(s): Tumor bulk" in txt
    # Only case 7495 survives that filter, so there is nothing to leave out and the script must
    # say so rather than crashing or implying robustness it cannot show.
    assert "leave one case out: SKIPPED" in txt


def test_unknown_region_fails_loudly():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        ann, xw = _write_annotation_fixture(tmp)
        _write_two_donor_table(tmp / "amp.csv", lambda fov: -0.8)
        _write_two_donor_table(tmp / "ctl.csv", lambda fov: 0.0)
        out = subprocess.run([sys.executable, str(SCRIPT), "--amplicon", str(tmp / "amp.csv"),
                              "--control", str(tmp / "ctl.csv"), "--annotations", str(ann),
                              "--crosswalk", str(xw), "--regions", "Necrosis"],
                             capture_output=True, text=True)
    assert out.returncode != 0
    assert "unknown region" in (out.stdout + out.stderr)


def test_grouping_by_case_requires_the_annotation_files():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_two_donor_table(tmp / "amp.csv", lambda fov: -0.8)
        _write_two_donor_table(tmp / "ctl.csv", lambda fov: 0.0)
        out = subprocess.run([sys.executable, str(SCRIPT), "--amplicon", str(tmp / "amp.csv"),
                              "--control", str(tmp / "ctl.csv")], capture_output=True, text=True)
    assert out.returncode != 0
    assert "needs --annotations" in (out.stdout + out.stderr)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
