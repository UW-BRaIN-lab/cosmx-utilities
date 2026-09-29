#!/usr/bin/env python3
"""Tests for denovo_under_fixed_profiles.py (no S3 / network).

The three that matter:
  * a letter whose cells land on its matching leaf reads as a SAFE drop;
  * a letter whose cells go to Low_signal reads as an UNSAFE one;
  * a broken cell-id join FAILS LOUDLY rather than reporting an empty comparison as a result.
The third is the point of the whole design -- per-cell InSituCNV already turned out not to join
to the full cohort, so a silent zero-overlap answer is a real failure mode here, not a
hypothetical one.

Runnable under pytest or directly:
    uv run --with pandas --with anndata python pipeline/python/tests/test_denovo_under_fixed_profiles.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
SCRIPT = _PY_DIR / "denovo_under_fixed_profiles.py"

_spec = importlib.util.spec_from_file_location("dufp", SCRIPT)
dufp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dufp)


def test_denovo_letters_are_recognised_by_shape():
    """InSituType's cluster_name_pool is 1-2 lowercase letters; K=27 overflows to 'aa'."""
    for label in ("l", "b", "aa", "z"):
        assert dufp.is_denovo(label), label
    for label in ("Neuron", "Endo_capilar", "AC-like", "TAM-MG_prolif", "OPC"):
        assert not dufp.is_denovo(label), label


def _write_inputs(tmp: Path, fixed_of, ids=None, typed_ids=None):
    """An anchor label table and a fixed-profile typed h5ad over the same cells."""
    import anndata as ad
    import h5py
    ids = ids or [f"SL01_F{i % 50}_C{i}" for i in range(600)]
    anchor = ["l"] * 200 + ["aa"] * 200 + ["Oligodendrocyte"] * 200
    with h5py.File(tmp / "anchor.h5", "w") as f:
        f.create_dataset("cell_id", data=np.array(ids, dtype="S64"))
        f.create_dataset("cell_type", data=np.array(anchor, dtype="S64"))
    tids = typed_ids if typed_ids is not None else ids
    obs = pd.DataFrame({"cell_type": [fixed_of(a) for a in anchor]}, index=pd.Index(tids))
    ad.AnnData(X=np.zeros((len(tids), 1), dtype="float32"), obs=obs).write_h5ad(tmp / "typed.h5ad")
    return tmp / "anchor.h5", tmp / "typed.h5ad"


def _run(tmp: Path, *extra) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), "--anchor-h5", str(tmp / "anchor.h5"),
                           "--typed-h5ad", str(tmp / "typed.h5ad"),
                           "--output-dir", str(tmp / "out"), *extra],
                          capture_output=True, text=True)


def test_a_safe_drop_shows_its_cells_landing_on_the_leaf():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_inputs(tmp, lambda a: {"l": "Endo_capilar", "aa": "Oligodendrocyte"}.get(a, a))
        out = _run(tmp)
        assert out.returncode == 0, out.stderr
        got = pd.read_csv(tmp / "out" / "fate_by_anchor_label.csv").set_index("anchor_label")
        assert got.loc["l", "top_destination"] == "Endo_capilar"
        assert got.loc["l", "pct_low_signal"] == 0.0
        assert got.loc["aa", "top_destination"] == "Oligodendrocyte"


def test_an_unsafe_drop_shows_its_cells_in_the_sink():
    """The leaf exists, the fixed profile still cannot hold the cells."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_inputs(tmp, lambda a: "Low_signal" if a == "l" else a)
        out = _run(tmp)
        assert out.returncode == 0, out.stderr
        got = pd.read_csv(tmp / "out" / "fate_by_anchor_label.csv").set_index("anchor_label")
        assert got.loc["l", "pct_low_signal"] == 100.0
        assert "worse than baseline" in out.stdout


def test_a_broken_cell_id_join_fails_loudly():
    """The failure mode that actually happened once: two runs, incompatible id formats."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_inputs(tmp, lambda a: a,
                      typed_ids=[f"c{i}" for i in range(600)])   # a different id scheme
        out = _run(tmp)
        assert out.returncode != 0
        combined = out.stdout + out.stderr
        assert "do not share a cell-id format" in combined
        assert "id example" in combined


def test_a_missing_celltype_column_is_named_not_guessed():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _write_inputs(tmp, lambda a: a)
        out = _run(tmp, "--celltype-key", "not_a_column")
        assert out.returncode != 0
        assert "no column" in (out.stdout + out.stderr)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
