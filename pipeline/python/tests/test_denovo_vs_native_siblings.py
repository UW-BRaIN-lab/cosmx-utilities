#!/usr/bin/env python3
"""Tests for sibling_vs_destination and the excluded-gene amplitude in denovo_vs_native_pseudobulk.py.

Two planted worlds, run through the real CLI on synthetic counts/typing/forced files:

  (top_n is kept at 6 for the sibling tests so each group's markers are its real signal; with a
   larger top_n every group also picks noise genes and the correlations measure noise.)

  STRESS   every forced cell of the letter carries one shared stress programme and almost none of
           its destination's markers (the shape of de-novo `t`), while B's native cells carry
           their own markers PLUS some of the same stress genes. Siblings must then resemble each
           other far more than they resemble their destinations, and B's amplitude must be
           inflated by the stress genes it shares -- so excluding them must shrink it.
  DEST     each forced subset looks like its own destination and unlike its sibling (real
           substructure). The sibling-minus-destination sign must flip.

Without both, a statistic that always said "siblings win" would be worthless.

Run: python pipeline/python/tests/test_denovo_vs_native_siblings.py   (or under pytest)
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
DEST_A, DEST_B, LETTER = "Endo_arterial", "Pericyte", "l"
N_GENES, N = 40, 150
MARK = {DEST_A: range(0, 6), DEST_B: range(6, 12)}
STRESS = list(range(20, 26))
STRESS_NAMES = [f"G{i}" for i in STRESS]
BASE = 1.0


def _write(tmp: Path, groups: dict) -> None:
    """groups: name -> (semisup, forced, {gene index: rate})."""
    rng = np.random.default_rng(11)
    ids, semisup, forced, cols = [], [], [], []
    for name, (sem, frc, rates) in groups.items():
        cid = [f"{name}_F1_C{i}" for i in range(N)]
        ids += cid; semisup += [sem] * N; forced += [frc] * N
        lam = np.full((N_GENES, N), BASE)
        for g, r in rates.items():
            lam[g, :] = r
        cols.append(rng.poisson(lam))
    counts = sp.csc_matrix(np.hstack(cols).astype(float))
    with h5py.File(tmp / "counts.h5", "w") as f:
        f["counts/shape"] = np.array(counts.shape)
        f["counts/data"] = counts.data
        f["counts/indices"] = counts.indices
        f["counts/indptr"] = counts.indptr
        f["genes"] = np.array([f"G{i}" for i in range(N_GENES)], dtype="S")
        f["cell_id"] = np.array(ids, dtype="S")
    with h5py.File(tmp / "typing.h5", "w") as f:
        f["cell_id"] = np.array(ids, dtype="S")
        f["cell_type"] = np.array(semisup, dtype="S")
    pd.DataFrame({"cell_id": ids, "top1_type": forced}).to_csv(tmp / "forced.csv", index=False)
    pinned = [c for c in ids if c.startswith("nat_") and "_pinned_" in c]
    pd.DataFrame({"cell_id": pinned}).to_csv(tmp / "pinned.csv", index=False)


def _rates(marker_dest=None, marker_rate=0.0, stress_rate=0.0):
    r = {}
    if marker_dest:
        r.update({g: marker_rate for g in MARK[marker_dest]})
    r.update({g: stress_rate for g in STRESS} if stress_rate else {})
    return r


def _stress_world() -> dict:
    return {
        "forced_A": (LETTER, DEST_A, _rates(DEST_A, 1.5, 8.0)),
        "forced_B": (LETTER, DEST_B, _rates(DEST_B, 1.5, 8.0)),
        "nat_A_unpinned": (DEST_A, DEST_A, _rates(DEST_A, 8.0)),
        "nat_A_pinned": (DEST_A, DEST_A, _rates(DEST_A, 9.0)),
        "nat_B_unpinned": (DEST_B, DEST_B, _rates(DEST_B, 8.0, 5.0)),
        "nat_B_pinned": (DEST_B, DEST_B, _rates(DEST_B, 9.0, 5.0)),
    }


def _dest_world() -> dict:
    return {
        "forced_A": (LETTER, DEST_A, _rates(DEST_A, 7.0)),
        "forced_B": (LETTER, DEST_B, _rates(DEST_B, 7.0)),
        "nat_A_unpinned": (DEST_A, DEST_A, _rates(DEST_A, 8.0)),
        "nat_A_pinned": (DEST_A, DEST_A, _rates(DEST_A, 9.0)),
        "nat_B_unpinned": (DEST_B, DEST_B, _rates(DEST_B, 8.0)),
        "nat_B_pinned": (DEST_B, DEST_B, _rates(DEST_B, 9.0)),
    }


def _run(groups: dict, *extra: str, top_n: int = 6) -> Path:
    tmp = Path(tempfile.mkdtemp())
    _write(tmp, groups)
    out = tmp / "out"
    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--counts-h5", str(tmp / "counts.h5"),
         "--typing-h5", str(tmp / "typing.h5"), "--forced-csv", str(tmp / "forced.csv"),
         "--letter", LETTER, "--destinations", f"{DEST_A},{DEST_B}", "--top-n", str(top_n),
         "--min-group-n", "20", "--pinned-csv", str(tmp / "pinned.csv"),
         "--output-dir", str(out), *extra],
        capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout, res.stderr)
        raise AssertionError("script failed")
    return out


def test_keyword_and_gene_list_resolution():
    assert dvn.resolve_excluded_genes(None) == set()
    assert dvn.resolve_excluded_genes("heat-shock") == set(dvn.HEAT_SHOCK_GENES)
    got = dvn.resolve_excluded_genes("Heat-Shock, FOS")
    assert got == set(dvn.HEAT_SHOCK_GENES) | {"FOS"}, got


def test_stress_world_siblings_beat_destinations():
    sib = pd.read_csv(_run(_stress_world()) / "sibling_vs_destination.csv")
    print(sib[["forced", "native_ref", "r_siblings_mean", "r_native", "siblings_minus_native"]]
          .to_string(index=False))
    unpinned = sib[sib.native_ref.str.endswith("unpinned]")]
    assert len(unpinned) == 2, unpinned
    assert (unpinned["siblings_minus_native"] > 0.3).all(), unpinned["siblings_minus_native"]
    # the z-score across k columns has a negative null r of -1/(k-1)
    assert np.allclose(sib["null_r"], -1.0 / (sib["n_columns"] - 1))


def test_destination_world_flips_the_sign():
    sib = pd.read_csv(_run(_dest_world()) / "sibling_vs_destination.csv")
    unpinned = sib[sib.native_ref.str.endswith("unpinned]")]
    print(unpinned[["forced", "r_siblings_mean", "r_native", "siblings_minus_native"]]
          .to_string(index=False))
    assert (unpinned["siblings_minus_native"] < 0).all(), unpinned["siblings_minus_native"]


def test_excluding_shared_stress_genes_shrinks_the_inflated_amplitude():
    # top_n=12 so B's marker set holds its 6 real markers AND the 6 stress genes it shares.
    out = _run(_stress_world(), "--exclude-genes", ",".join(STRESS_NAMES), top_n=12)
    amp = pd.read_csv(out / "marker_amplitude.csv")
    forced_b = amp[(amp.destination == DEST_B) & (amp.group == f"{LETTER}->{DEST_B}")
                   ].set_index("marker_set")["amplitude"]
    print(forced_b.to_dict())
    assert set(forced_b.index) == {"all", "excluding_listed"}
    assert forced_b["all"] - forced_b["excluding_listed"] > 0.3, forced_b.to_dict()
    # a destination that shares no stress genes is unaffected
    forced_a = amp[(amp.destination == DEST_A) & (amp.group == f"{LETTER}->{DEST_A}")
                   ].set_index("marker_set")["amplitude"]
    assert abs(forced_a["all"] - forced_a["excluding_listed"]) < 0.05, forced_a.to_dict()
    genes = pd.read_csv(out / "marker_amplitude_genes.csv")
    flagged = set(genes.loc[genes.is_excluded, "gene"])
    assert flagged <= set(STRESS_NAMES) and flagged, flagged
    # the per-gene contributions of the flagged genes explain the inflation
    fb = genes[(genes.destination == DEST_B) & (genes.group == f"{LETTER}->{DEST_B}")]
    assert fb.loc[fb.is_excluded, "contribution"].mean() > fb.loc[~fb.is_excluded,
                                                                 "contribution"].mean()


def test_without_exclude_genes_only_one_marker_set_is_reported():
    amp = pd.read_csv(_run(_stress_world()) / "marker_amplitude.csv")
    assert set(amp["marker_set"]) == {"all"}


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
