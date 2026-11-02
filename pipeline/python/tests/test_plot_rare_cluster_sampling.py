#!/usr/bin/env python3
"""Tests for the rare-cluster sampling figure: the joins, the guards and an end-to-end render (no S3).

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_plot_rare_cluster_sampling.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location("plot_rare_cluster_sampling",
                                               _PY_DIR / "plot_rare_cluster_sampling.py")
prcs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(prcs)

CELL_IDS = pd.Index([f"slideA_F1_C{i}" for i in range(1, 9)] + [f"slideB_F1_C{i}" for i in range(1, 5)])
LEIDEN = pd.Series(["0"] * 6 + ["21"] * 3 + ["0"] * 3)
ANCHOR_IDS = pd.Series(["slideA_F1_C1", "slideA_F1_C2", "slideA_F1_C7", "slideA_F1_C8", "slideB_F1_C2"])


def test_membership_marks_anchor_cells_in_adata_order():
    mask = prcs.anchor_membership(CELL_IDS, ANCHOR_IDS)
    print("membership:", mask.astype(int))
    assert mask.sum() == len(ANCHOR_IDS)
    assert list(CELL_IDS[mask]) == ["slideA_F1_C1", "slideA_F1_C2", "slideA_F1_C7", "slideA_F1_C8",
                                    "slideB_F1_C2"]


def test_membership_refuses_mismatched_runs_and_duplicates():
    for bad, needle in ((pd.Series(["other_F1_C1", "other_F1_C2"]), "only"),
                        (pd.Series(["slideA_F1_C1", "slideA_F1_C1"]), "duplicated")):
        try:
            prcs.anchor_membership(CELL_IDS, bad)
        except ValueError as error:
            print("raised:", error)
            assert needle in str(error)
        else:
            raise AssertionError("a bad anchor list must raise")


def test_per_slide_labels_are_refused():
    try:
        prcs.reject_per_slide_labels(pd.Series(["slideA|3", "slideA|4"]))
    except ValueError as error:
        assert "per-slide" in str(error)
    else:
        raise AssertionError("slide|n labels must be refused")
    prcs.reject_per_slide_labels(LEIDEN)


def test_sampling_table_counts_and_orders_by_size():
    mask = prcs.anchor_membership(CELL_IDS, ANCHOR_IDS)
    table = prcs.cluster_sampling_table(LEIDEN, mask)
    print(table)
    assert list(table.index) == ["0", "21"]
    assert table.loc["21", "n_cells"] == 3 and table.loc["21", "n_sampled"] == 2
    assert abs(table.loc["21", "share_sampled"] - 2 / 3) < 1e-9
    assert table["n_cells"].sum() == len(CELL_IDS)


def test_pick_slide_is_the_one_with_most_target_cells():
    slides = pd.Series(["slideA"] * 8 + ["slideB"] * 4)
    assert prcs.pick_slide(slides, (LEIDEN == "21").to_numpy()) == "slideA"


STRATA = pd.DataFrame({
    "slide_cluster": ["slideA|0"] * 5 + ["slideA|1"] * 3 + ["slideA|2"] * 2 + ["slideB|0"] * 5,
    "is_target": [False] * 5 + [True, True, False] + [True, True] + [False] * 5})


def test_sampling_cap_is_the_largest_stratum():
    assert prcs.sampling_cap(STRATA["slide_cluster"]) == 5


def test_slide_strata_table_splits_target_and_rest_per_slide():
    table = prcs.slide_strata_table(STRATA, "slideA")
    print(table)
    assert list(table.index) == ["0", "1", "2"]
    assert list(table["n_sampled"]) == [5, 3, 2]
    assert list(table["n_target"]) == [0, 2, 2]
    assert prcs.slide_strata_table(STRATA, "slideB")["n_sampled"].sum() == 5


def test_slide_strata_table_unknown_slide_raises():
    try:
        prcs.slide_strata_table(STRATA, "nope")
    except ValueError as error:
        assert "no anchor cells" in str(error)
    else:
        raise AssertionError("an unknown slide must raise")


def test_spatial_columns_priority_and_missing():
    assert prcs.spatial_columns(pd.Index(["CenterX_global_px", "CenterY_global_px", "x"])) == (
        "CenterX_global_px", "CenterY_global_px")
    try:
        prcs.spatial_columns(pd.Index(["leiden"]))
    except KeyError:
        return
    raise AssertionError("missing centroid columns must raise")


def test_end_to_end_render():
    rng = np.random.default_rng(1)
    n = 400
    ids = [f"s{i % 3}_F1_C{i}" for i in range(n)]
    obs = pd.DataFrame({"leiden": rng.choice(["0", "1", "21"], n, p=[.6, .35, .05]),
                        "slide_id": [i.split("_")[0] for i in ids],
                        "CenterX_global_px": rng.uniform(0, 1000, n),
                        "CenterY_global_px": rng.uniform(0, 1000, n)}, index=ids)
    adata = ad.AnnData(np.zeros((n, 2), dtype=np.float32), obs=obs)
    adata.obsm["X_umap"] = rng.normal(size=(n, 2))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        adata.write_h5ad(tmp / "typed.h5ad")
        kept = ids[::2]
        pd.DataFrame({"cell_id": kept,
                      "slide_cluster": [f"{cell.split('_')[0]}|{k % 4}" for k, cell in enumerate(kept)]}
                     ).to_csv(tmp / "anchor_cells.csv", index=False)
        result = subprocess.run(
            [sys.executable, str(_PY_DIR / "plot_rare_cluster_sampling.py"), "--typed-h5ad",
             str(tmp / "typed.h5ad"), "--anchor-cells-csv", str(tmp / "anchor_cells.csv"),
             "--cluster", "21", "--output-stem", str(tmp / "fig")],
            capture_output=True, text=True)
        print(result.stdout, result.stderr)
        assert result.returncode == 0
        assert (tmp / "fig.png").stat().st_size > 10_000
        summary = pd.read_csv(tmp / "fig_sampling_by_cluster.csv", dtype={"leiden": str})
        assert set(summary["leiden"]) == {"0", "1", "21"}
        assert summary["n_cells"].sum() == n


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            print(f"== {name}")
            fn()
    print("all passed")
