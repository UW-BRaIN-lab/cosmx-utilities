#!/usr/bin/env python3
"""Tests for denovo_vs_native_by_group.py — the donor-structure test (no S3 / network).

The pair that matters is the two cases that must come out opposite: a letter whose cells are
interleaved with the native ones inside every donor must NOT look donor-structured, and a letter
that is a per-donor offset must. Without both, a dispersion statistic that always looked alarming
would be useless — and the whole point of the test is to distinguish those two worlds.

Runnable either under pytest or directly:
    uv run --with pandas --with numpy \\
        python pipeline/python/tests/test_denovo_vs_native_by_group.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))

_spec = importlib.util.spec_from_file_location(
    "denovo_vs_native_by_group", _PY_DIR / "denovo_vs_native_by_group.py")
dn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dn)

DONORS = [f"D{i}" for i in range(10)]
PER_DONOR = 200


def _cells(share_of: dict[str, float]) -> pd.DataFrame:
    """One row per cell: donor, and whether it went to the letter (A) or native D (B)."""
    rng = np.random.default_rng(7)
    rows = []
    for donor in DONORS:
        n_a = rng.binomial(PER_DONOR, share_of[donor])
        rows += [{"Case": donor, "cell_type": "l"}] * n_a
        rows += [{"Case": donor, "cell_type": "Endo_arterial"}] * (PER_DONOR - n_a)
    return pd.DataFrame(rows)


def _phi(cells: pd.DataFrame) -> dict:
    tbl = dn.counts_by_group(cells, "Case", cells["cell_type"] == "l",
                             cells["cell_type"] == "Endo_arterial")
    return dn.dispersion(tbl, min_cells=25)


def test_interleaved_letter_is_not_donor_structured():
    """Every donor sends the same 40% to the letter: phi must sit near 1."""
    got = _phi(_cells({d: 0.4 for d in DONORS}))
    assert got["groups"] == len(DONORS), got
    assert got["phi"] < 3, f"a homogeneous split must not look donor-structured: {got}"
    assert got["polarised"] == 0, got
    assert 0.35 < got["share"] < 0.45, got


def test_a_per_donor_offset_is_caught():
    """Half the donors send nearly everything to the letter, half nearly nothing."""
    share = {d: (0.97 if i % 2 else 0.03) for i, d in enumerate(DONORS)}
    got = _phi(_cells(share))
    assert got["phi"] > 50, f"a donor-split letter must be flagged: {got}"
    assert got["polarised"] == len(DONORS), got


def test_the_two_worlds_are_orders_of_magnitude_apart():
    """The statistic has to separate them by a lot, not marginally."""
    even = _phi(_cells({d: 0.4 for d in DONORS}))
    split = _phi(_cells({d: (0.97 if i % 2 else 0.03) for i, d in enumerate(DONORS)}))
    assert split["phi"] > 20 * even["phi"], (even, split)


def test_thin_groups_are_dropped_not_averaged():
    """A donor with a handful of cells must not contribute a wild share to the dispersion."""
    cells = _cells({d: 0.4 for d in DONORS})
    cells = pd.concat([cells, pd.DataFrame([{"Case": "TINY", "cell_type": "l"}] * 3)])
    tbl = dn.counts_by_group(cells, "Case", cells["cell_type"] == "l",
                             cells["cell_type"] == "Endo_arterial")
    assert "TINY" in tbl.index
    assert dn.dispersion(tbl, min_cells=25)["groups"] == len(DONORS)


def test_too_few_groups_returns_nan_rather_than_a_number():
    """Two donors cannot support a dispersion claim; it must refuse rather than mislead."""
    cells = pd.DataFrame([{"Case": "D0", "cell_type": "l"}] * 50
                         + [{"Case": "D0", "cell_type": "Endo_arterial"}] * 50
                         + [{"Case": "D1", "cell_type": "l"}] * 50
                         + [{"Case": "D1", "cell_type": "Endo_arterial"}] * 50)
    tbl = dn.counts_by_group(cells, "Case", cells["cell_type"] == "l",
                             cells["cell_type"] == "Endo_arterial")
    assert np.isnan(dn.dispersion(tbl, min_cells=25)["phi"])


def test_a_one_sided_split_does_not_divide_by_zero():
    """If every cell went one way the share is degenerate; phi must be nan, not an exception."""
    cells = pd.DataFrame([{"Case": d, "cell_type": "l"} for d in DONORS for _ in range(50)])
    tbl = dn.counts_by_group(cells, "Case", cells["cell_type"] == "l",
                             cells["cell_type"] == "Endo_arterial")
    got = dn.dispersion(tbl, min_cells=25)
    assert np.isnan(got["phi"]) and got["share"] == 1.0, got


def test_named_pair_null_uses_only_named_siblings():
    """The yardstick must be built from fixed-profile pairs, never from the letter itself."""
    rng = np.random.default_rng(3)
    rows = []
    for donor in DONORS:
        for cell_type, n in (("l", 40), ("Endo_arterial", 80), ("Endo_capilar", 80),
                             ("Pericyte", 80)):
            rows += [{"Case": donor, "cell_type": cell_type}] * int(rng.normal(n, 4))
    cells = pd.DataFrame(rows)
    null = dn.named_pair_null(cells, "Case",
                              ["Endo_arterial", "Endo_capilar", "Pericyte"], min_cells=25)
    assert len(null) == 3, null           # 3 choose 2
    # Substring checks false-positive here ("arteriaL VS Pericyte"), so compare the members.
    members = {m for pair in null["pair"] for m in pair.split(" vs ")}
    assert "l" not in members, members
    assert null["phi"].max() < 5, f"evenly-composed named pairs must give a quiet null: {null}"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
