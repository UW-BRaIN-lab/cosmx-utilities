#!/usr/bin/env python3
"""Tests for the pinned-cell split/exclude in denovo_vs_native_pseudobulk.py.

The scenario is the hypothesis being tested: natives selected as InSituType's best exemplars
carry their type's markers MORE strongly than ordinary cells of that type. Planted here as
pinned natives at 2x the marker rate of unpinned natives and of the letter's forced cells, so
the pooled native group looks stronger than the forced cells (the "half amplitude" pattern)
while the UNPINNED natives match them. The real script must reproduce exactly that, on a
synthetic counts/typing/forced/pinned fileset run through the CLI.

Run: python pipeline/python/tests/test_denovo_vs_native_pinned.py   (or under pytest)
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PY_DIR))
import denovo_vs_native_pseudobulk as dvn  # noqa: E402

SCRIPT = PY_DIR / "denovo_vs_native_pseudobulk.py"
DEST_A, DEST_B = "Endo_arterial", "Pericyte"
LETTER = "l"
N_GENES = 40
MARKERS = {DEST_A: range(0, 6), DEST_B: range(6, 12)}
LOW_RATE, HIGH_RATE, BASE_RATE = 6.0, 12.0, 1.0
N = 150   # cells per planted group


def _cells(prefix: str, n: int) -> list[str]:
    return [f"{prefix}_F1_C{i}" for i in range(n)]


def _plant(tmp: Path):
    rng = np.random.default_rng(3)
    groups = {  # name -> (semisup label, forced label, pinned?, marker dest, marker rate)
        "forced_A": (LETTER, DEST_A, False, DEST_A, LOW_RATE),
        "forced_B": (LETTER, DEST_B, False, DEST_B, LOW_RATE),
        "nat_A_pinned": (DEST_A, DEST_A, True, DEST_A, HIGH_RATE),
        "nat_A_unpinned": (DEST_A, DEST_A, False, DEST_A, LOW_RATE),
        "nat_B_pinned": (DEST_B, DEST_B, True, DEST_B, HIGH_RATE),
        "nat_B_unpinned": (DEST_B, DEST_B, False, DEST_B, LOW_RATE),
    }
    ids, semisup, forced, pinned, cols = [], [], [], [], []
    for gname, (sem, frc, pin, mdest, rate) in groups.items():
        cid = _cells(gname, N)
        ids += cid; semisup += [sem] * N; forced += [frc] * N
        if pin:
            pinned += cid
        lam = np.full((N_GENES, N), BASE_RATE)
        lam[list(MARKERS[mdest]), :] = rate
        cols.append(rng.poisson(lam))
    counts = sp.csc_matrix(np.hstack(cols).astype(float))
    genes = [f"G{i}" for i in range(N_GENES)]

    with h5py.File(tmp / "counts.h5", "w") as f:
        f["counts/shape"] = np.array(counts.shape)
        f["counts/data"] = counts.data
        f["counts/indices"] = counts.indices
        f["counts/indptr"] = counts.indptr
        f["genes"] = np.array(genes, dtype="S")
        f["cell_id"] = np.array(ids, dtype="S")
    with h5py.File(tmp / "typing.h5", "w") as f:
        f["cell_id"] = np.array(ids, dtype="S")
        f["cell_type"] = np.array(semisup, dtype="S")
    pd.DataFrame({"cell_id": ids, "top1_type": forced}).to_csv(tmp / "forced.csv", index=False)
    pd.DataFrame({"cell_id": pinned, "pinned_type": "x", "label": "x"}).to_csv(
        tmp / "pinned.csv", index=False)
    return ids, pinned


def _run(tmp: Path, *extra: str) -> tuple[Path, str]:
    out = tmp / "out"
    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--counts-h5", str(tmp / "counts.h5"),
         "--typing-h5", str(tmp / "typing.h5"), "--forced-csv", str(tmp / "forced.csv"),
         "--letter", LETTER, "--destinations", f"{DEST_A},{DEST_B}", "--top-n", "6",
         "--min-group-n", "20", "--output-dir", str(out), *extra],
        capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout, res.stderr)
        raise AssertionError("script failed")
    return out, res.stdout


def test_split_relabels_only_native_groups():
    group = pd.Series(["l->A", "A [native]", "A [native]", "B [native]", pd.NA],
                      index=["c1", "c2", "c3", "c4", "c5"], dtype="object")
    got = dvn.split_natives_by_pinned(group, {"c2", "c4", "c1"}, "split")
    assert list(got[:4]) == ["l->A", "A [native, pinned]", "A [native, unpinned]",
                             "B [native, pinned]"], list(got)
    assert pd.isna(got["c5"])
    ex = dvn.split_natives_by_pinned(group, {"c2", "c1"}, "exclude")
    assert pd.isna(ex["c2"]) and ex["c3"] == "A [native]" and ex["c1"] == "l->A", list(ex)


def test_bad_mode_is_rejected():
    try:
        dvn.split_natives_by_pinned(pd.Series(["A [native]"], index=["c"]), set(), "keep")
    except ValueError:
        return
    raise AssertionError("an unknown mode must raise")


def test_split_shows_the_selection_effect_in_amplitude():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _plant(tmp)
        out, stdout = _run(tmp, "--pinned-csv", str(tmp / "pinned.csv"),
                           "--pinned-mode", "split")
        amp = pd.read_csv(out / "marker_amplitude.csv")
    print(amp.to_string(index=False))
    a = amp[amp.destination == DEST_A].set_index("group")["amplitude"]
    forced, pooled = a[f"{LETTER}->{DEST_A}"], a[f"{DEST_A} [native, all]"]
    pinned, unpinned = a[f"{DEST_A} [native, pinned]"], a[f"{DEST_A} [native, unpinned]"]
    # The pooled native group outranks the forced cells only because of the pinned exemplars ...
    # (log-normalised, so compare differences rather than ratios)
    assert pinned - forced > 0.2, (pinned, forced)
    assert pooled > forced, (pooled, forced)
    # ... and the natives that were NOT selected look like the forced cells: their gap is noise
    # next to the pinned gap.
    assert abs(unpinned - forced) < 0.05, (unpinned, forced)
    assert (pinned - forced) > 5 * abs(unpinned - forced), (pinned, unpinned, forced)
    assert "expected 0" in stdout


def test_exclude_drops_pinned_cells_and_keeps_the_pooled_baseline():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _plant(tmp)
        out, _ = _run(tmp, "--pinned-csv", str(tmp / "pinned.csv"), "--pinned-mode", "exclude")
        sizes = pd.read_csv(out / "group_sizes.csv").set_index("group")["n_cells"]
        amp = pd.read_csv(out / "marker_amplitude.csv")
    assert sizes[f"{DEST_A} [native]"] == N, sizes.to_dict()      # only the unpinned half
    assert not any("pinned" in g for g in sizes.index), list(sizes.index)
    # the pooled row still reflects ALL natives, so it stays above the unpinned-only column
    a = amp[amp.destination == DEST_A].set_index("group")["amplitude"]
    assert a[f"{DEST_A} [native, all]"] > a[f"{DEST_A} [native]"], a.to_dict()


def test_without_a_pinned_file_the_original_columns_are_unchanged():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _plant(tmp)
        out, _ = _run(tmp)
        sizes = pd.read_csv(out / "group_sizes.csv").set_index("group")["n_cells"]
    assert sizes[f"{DEST_A} [native]"] == 2 * N and sizes[f"{LETTER}->{DEST_A}"] == N
    assert not any("pinned" in g for g in sizes.index)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
