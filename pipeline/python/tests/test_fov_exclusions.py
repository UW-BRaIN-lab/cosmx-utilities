#!/usr/bin/env python3
"""Tests for the shared FOV exclusion, including the InSituType/InSituTree input path.

The load-bearing property is that a mismatch is FATAL. An exclusion that silently matches
nothing hands the caller contaminated data labelled as filtered, which is worse than not
offering the option at all.

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_fov_exclusions.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location("fov_exclusions", _PY_DIR / "fov_exclusions.py")
fx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fx)

SLIDE, OTHER_SLIDE = "uwa0001eyes", "uwa0002eyes"
EXCLUDED_FOVS = (2, 3)
N_PER_FOV = 10


def build_obs() -> pd.DataFrame:
    rows = []
    for slide in (SLIDE, OTHER_SLIDE):
        for fov in (1, 2, 3, 4):
            for cell in range(N_PER_FOV):
                rows.append({"slide_id": slide, "fov": fov,
                             "cell": f"{slide}_F{fov}_C{cell}"})
    return pd.DataFrame(rows).set_index("cell")


def write_annotations(path: Path, slide: str = SLIDE) -> Path:
    pd.DataFrame([{"slide_id": slide, "fov": fov, "region": "Retina",
                   "mixed_adjacent": 0, "exclude": int(fov in EXCLUDED_FOVS)}
                  for fov in (1, 2, 3, 4)]).to_csv(path, index=False)
    return path


def test_drops_only_the_named_slide(tmp_path: Path) -> None:
    obs = build_obs()
    kept = fx.drop_excluded_fovs(obs, write_annotations(tmp_path / "fov.csv"))
    assert len(kept) == len(obs) - len(EXCLUDED_FOVS) * N_PER_FOV
    # The other slide shares FOV numbers 2 and 3 and must be untouched -- the exclusion is
    # keyed on (slide, FOV), and dropping by FOV alone would silently take 20 more cells.
    assert (kept.slide_id == OTHER_SLIDE).sum() == 4 * N_PER_FOV
    assert not ((kept.slide_id == SLIDE) & kept.fov.isin(EXCLUDED_FOVS)).any()
    print("  drops the right (slide, FOV) pairs and nothing else")


def test_mismatch_is_fatal(tmp_path: Path) -> None:
    obs = build_obs()
    annotations = write_annotations(tmp_path / "bad.csv", slide="not_a_slide")
    try:
        fx.drop_excluded_fovs(obs, annotations)
    except SystemExit as exc:
        assert "matched any cell" in str(exc)
    else:
        raise AssertionError("a slide-id mismatch must exit, not silently drop nothing")

    # Nothing marked excluded is a legitimate no-op, not an error.
    pd.DataFrame([{"slide_id": SLIDE, "fov": 1, "region": "Retina",
                   "mixed_adjacent": 0, "exclude": 0}]).to_csv(tmp_path / "none.csv", index=False)
    assert len(fx.drop_excluded_fovs(obs, tmp_path / "none.csv")) == len(obs)
    print("  a mismatch exits; an empty exclusion list is a no-op")


def test_insitutype_input_is_filtered(tmp_path: Path) -> None:
    """End to end: the h5 that InSituType and InSituTree read must be short those cells."""
    obs = build_obs()
    n_genes = 6
    var = pd.DataFrame(
        {"probe_type": ["gene"] * (n_genes - 2) + ["negprobe"] * 2},
        index=[f"GENE{g}" for g in range(n_genes - 2)] + ["NegPrb1", "NegPrb2"])
    adata = ad.AnnData(X=sp.csr_matrix(np.ones((len(obs), n_genes), dtype=np.float32)),
                       obs=obs, var=var)
    h5ad = tmp_path / "combined_qc.h5ad"
    adata.write_h5ad(h5ad)

    out = tmp_path / "insitutype_input.h5"
    subprocess.run(
        [sys.executable, str(_PY_DIR / "prep_insitutype_inputs.py"),
         "--combined-h5ad", str(h5ad),
         "--fov-annotations", str(write_annotations(tmp_path / "fov2.csv")),
         "--output", str(out)],
        check=True, capture_output=True, text=True,
    )
    expected = len(obs) - len(EXCLUDED_FOVS) * N_PER_FOV
    with h5py.File(out, "r") as f:
        assert f["counts/shape"][()][1] == expected, f["counts/shape"][()]
        assert len(f["cell_id"]) == expected
        assert len(f["neg"]) == expected
    print(f"  insitutype_input.h5 carries {expected} cells, not {len(obs)}")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("test_fov_exclusions")
        test_drops_only_the_named_slide(tmp)
        test_mismatch_is_fatal(tmp)
        test_insitutype_input_is_filtered(tmp)
    print("All tests passed.")


if __name__ == "__main__":
    main()
