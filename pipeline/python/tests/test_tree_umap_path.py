#!/usr/bin/env python3
"""The InSituTree typing must flow through 90 -> 95 with no code change.

concat_typing_results.py writes the same /cell_id, /cell_type, /prob schema that
write_celltypes.py reads from an InSituType result, so stage 90 can write the tree calls
into the clustered AnnData and stage 95 can draw them. The only thing stopping it was
names: both stages hard-coded cosmx_typed.h5ad, so a tree run would have overwritten the
InSituType run's file. This test pins the contract that makes the reuse legal.

Two things it guards that are easy to lose:
  - X_umap must survive stage 90, or stage 95 has nothing to draw.
  - stage 95 must work with NO annotation table, since a tree typing has no de-novo
    letters to rename.

Runnable under pytest or directly:
    uv run python pipeline/python/tests/test_tree_umap_path.py
"""
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
N_CELLS, N_GENES, CELLS_PER_FOV = 600, 6, 20
# A mix of reference leaves and de-novo leaves, as insitutree_75 carries.
LEAVES = ("denovo_oligodendrocyte_1", "OPC", "HRCA_retinal rod cell",
          "denovo_astrocyte_1", "Mueller_glia")


def _build(tmp_path: Path, n_typed: int = N_CELLS) -> tuple[Path, Path]:
    rng = np.random.default_rng(0)
    cells = [f"c_1_{i // CELLS_PER_FOV + 1}_{i % CELLS_PER_FOV}" for i in range(N_CELLS)]
    obs = pd.DataFrame({"slide_id": "s1",
                        "fov": [i // CELLS_PER_FOV + 1 for i in range(N_CELLS)],
                        "leiden": pd.Categorical([str(i % 4) for i in range(N_CELLS)])},
                       index=pd.Index(cells, name="cell_id"))
    adata = ad.AnnData(X=sp.csr_matrix(rng.poisson(3, (N_CELLS, N_GENES)).astype(np.float32)),
                       obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(N_GENES)]))
    adata.obsm["X_umap"] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    clustered = tmp_path / "cosmx_clustered.h5ad"
    adata.write_h5ad(clustered)

    result = tmp_path / "insitutree_result.h5"
    with h5py.File(result, "w") as f:
        for name, values in (("cell_id", cells[:n_typed]),
                             ("cell_type", [LEAVES[i % len(LEAVES)] for i in range(n_typed)])):
            f.create_dataset(name, data=np.asarray(values, dtype=object),
                             dtype=h5py.string_dtype("utf-8"))
        f.create_dataset("prob", data=rng.uniform(0.8, 1.0, n_typed))
    return clustered, result


def _run(script: str, *args: str) -> None:
    done = subprocess.run([sys.executable, str(_PY_DIR / script), *args],
                          capture_output=True, text=True, cwd=_PY_DIR)
    if done.returncode:
        raise AssertionError(f"{script} failed:\n{done.stdout[-2000:]}\n{done.stderr[-2000:]}")


def test_stage90_writes_tree_calls_and_keeps_the_umap(tmp_path: Path) -> None:
    clustered, result = _build(tmp_path)
    typed_path = tmp_path / "cosmx_typed_tree.h5ad"
    # No --qc-plots-dir: that path needs scanpy, which the Mac venv does not carry.
    _run("write_celltypes.py", "--clustered-h5ad", str(clustered),
         "--insitutype-h5", str(result), "--output", str(typed_path))

    typed = ad.read_h5ad(typed_path)
    assert typed.n_obs == N_CELLS
    assert set(typed.obs["cell_type"]) == set(LEAVES)
    assert "X_umap" in typed.obsm, "stage 95 draws X_umap; losing it breaks the figure"
    print(f"  90: {typed.n_obs} cells, {typed.obs['cell_type'].nunique()} leaves, X_umap kept")


def test_stage95_draws_without_an_annotation_table(tmp_path: Path) -> None:
    clustered, result = _build(tmp_path)
    typed_path = tmp_path / "cosmx_typed_tree.h5ad"
    _run("write_celltypes.py", "--clustered-h5ad", str(clustered),
         "--insitutype-h5", str(result), "--output", str(typed_path))

    umap, panels = tmp_path / "umap.png", tmp_path / "panels.png"
    _run("plot_celltype_umap.py", "--h5ad", str(typed_path), "--color-key", "cell_type",
         "--napari-colors", "--output", str(umap), "--title", "InSituTree, 75 leaves",
         "--panels-output", str(panels))
    for png in (umap, panels):
        assert png.stat().st_size > 5000, f"{png.name} is too small to be a real figure"
    print("  95: drew the UMAP and the per-label panels with Napari colours, no mapping")


def test_cells_the_tree_never_typed_do_not_break_the_figure(tmp_path: Path) -> None:
    """The tree typed 835,077 of the cohort's 849,510 cells, so NaNs are the normal case.

    The 14,433 difference is the ciliary/cornea block, excluded before typing. A label
    column with NaNs is the shape stage 95 will actually be handed, and an NA category
    rendering as a 15th "cell type" would be worse than it being dropped.
    """
    n_typed = N_CELLS - 100
    clustered, result = _build(tmp_path, n_typed=n_typed)
    typed_path = tmp_path / "cosmx_typed_tree.h5ad"
    _run("write_celltypes.py", "--clustered-h5ad", str(clustered),
         "--insitutype-h5", str(result), "--output", str(typed_path))

    typed = ad.read_h5ad(typed_path)
    assert int(pd.isna(typed.obs["cell_type"]).sum()) == N_CELLS - n_typed
    umap = tmp_path / "umap.png"
    _run("plot_celltype_umap.py", "--h5ad", str(typed_path), "--color-key", "cell_type",
         "--napari-colors", "--output", str(umap))
    assert umap.stat().st_size > 5000
    print(f"  untyped: {N_CELLS - n_typed} cells carry no call; the figure still draws")


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        print("test_tree_umap_path")
        for name, case in (("a", test_stage90_writes_tree_calls_and_keeps_the_umap),
                           ("b", test_stage95_draws_without_an_annotation_table),
                           ("c", test_cells_the_tree_never_typed_do_not_break_the_figure)):
            # pytest's tmp_path already exists; make one per case so direct runs match.
            case_dir = Path(d) / name
            case_dir.mkdir()
            case(case_dir)
    print("All tests passed.")


if __name__ == "__main__":
    main()
