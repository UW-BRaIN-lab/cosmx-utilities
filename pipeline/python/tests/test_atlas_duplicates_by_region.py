#!/usr/bin/env python3
"""Tests for atlas_duplicates_by_region.py.

The fixture is built so the naive reading is wrong. Gray matter holds 70% of the cells, so
raw percentages put nearly every type's modal region there -- including one whose cells are
genuinely enriched in retina. Only the cohort-normalised statistic recovers that, and in a
region that large enrichment is capped at 1/0.70 = 1.43, so a fixed enrichment bar could not
be met at all. obs['Region'] is deliberately stale and disagrees with the sheet everywhere.

    uv run python pipeline/python/tests/test_atlas_duplicates_by_region.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location(
    "atlas_duplicates_by_region", _PY_DIR / "atlas_duplicates_by_region.py")
ad_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ad_mod)

SLIDE = "uwa0001eyes"
STALE_REGION = "Cortex"          # what the AnnData claims; the sheet must win
# fov -> (region, excluded)
FOVS = {1: ("Gray matter", 0), 2: ("Gray matter", 0), 3: ("Gray matter", 0),
        4: ("Retina", 0), 5: ("Optic nerve", 0), 6: ("Optic nerve", 1)}

# type -> {region: n cells}. Gray matter ends up with 70% of the cohort.
POPULATION = {
    # REGIONAL: each atlas copy owns its own tissue.
    "HRCA_astrocyte": {"Retina": 750, "Gray matter": 200, "Optic nerve": 50},
    "Allen_Astrocyte": {"Gray matter": 1900, "Retina": 50, "Optic nerve": 50},
    "Mona_astrocyte": {"Optic nerve": 800, "Gray matter": 150, "Retina": 50},
    # PARTIAL, and the reason the raw table misleads: more of this type's cells sit in gray
    # matter than anywhere else, but it is retina where it is over-represented.
    "HRCA_microglial cell": {"Gray matter": 500, "Retina": 400, "Optic nerve": 100},
    "Allen_Microglia": {"Gray matter": 1900, "Retina": 50, "Optic nerve": 50},
    # MIXED: both copies live in the same tissue, so one of them is not earning its column.
    "Mona_fibroblast": {"Gray matter": 700, "Retina": 100, "Optic nerve": 100},
    "Allen_Fibroblast": {"Gray matter": 800, "Retina": 100, "Optic nerve": 100},
    # TOO FEW CELLS to judge either way.
    "Mona_pericyte": {"Gray matter": 60, "Retina": 20, "Optic nerve": 20},
    "Allen_pericyte": {"Gray matter": 60, "Retina": 20, "Optic nerve": 20},
    # Filler so gray matter dominates the cohort the way the real one does.
    "Oligodendrocyte": {"Gray matter": 4000, "Retina": 200, "Optic nerve": 200},
}
# Cells in the EXCLUDED optic-nerve FOV. On their own they would move Allen_Fibroblast's home
# to optic nerve and turn the fibroblast verdict REGIONAL, so the default must drop them.
EXCLUDED_CELLS = {"Allen_Fibroblast": 2000, "Mona_fibroblast": 20}

REGION_FOV = {"Gray matter": [1, 2, 3], "Retina": [4], "Optic nerve": [5]}
CLASSES = [("astrocyte", ["HRCA_astrocyte", "Allen_Astrocyte", "Mona_astrocyte"]),
           ("microglia", ["HRCA_microglial cell", "Allen_Microglia"]),
           ("fibroblast", ["Mona_fibroblast", "Allen_Fibroblast"]),
           ("pericyte", ["Mona_pericyte", "Allen_pericyte"])]


def write_fixture(tmp: Path) -> tuple[Path, Path, Path]:
    slides, fovs, types = [], [], []
    for cell_type, by_region in POPULATION.items():
        for region, n in by_region.items():
            choices = REGION_FOV[region]
            for i in range(n):
                slides.append(SLIDE)
                fovs.append(choices[i % len(choices)])
                types.append(cell_type)
    for cell_type, n in EXCLUDED_CELLS.items():
        slides.extend([SLIDE] * n)
        fovs.extend([6] * n)
        types.extend([cell_type] * n)

    obs = pd.DataFrame({
        "slide_id": slides,
        "fov": np.array(fovs, dtype=int),
        "cell_type": pd.Categorical(types),
        "Region": STALE_REGION,
    }, index=[f"c{i:06d}" for i in range(len(slides))])
    adata = ad.AnnData(X=sp.csr_matrix((len(obs), 2), dtype=np.float32), obs=obs)
    h5ad = tmp / "cosmx_typed.h5ad"
    adata.write_h5ad(h5ad)

    sheet = tmp / "fov_annotations.csv"
    pd.DataFrame([{"slide_id": SLIDE, "fov": fov, "region": region,
                   "mixed_adjacent": 0, "exclude": excluded}
                  for fov, (region, excluded) in FOVS.items()]).to_csv(sheet, index=False)

    classes = tmp / "classes.csv"
    pd.DataFrame([{"class": name, "type": t, "note": ""}
                  for name, members in CLASSES for t in members]).to_csv(classes, index=False)
    return h5ad, sheet, classes


def run(tmp: Path, *extra: str) -> pd.DataFrame:
    h5ad, sheet, classes = write_fixture(tmp)
    out = tmp / ("out" + "".join(extra).replace("-", ""))
    assert ad_mod.main(["--typed-h5ad", str(h5ad), "--fov-annotations", str(sheet),
                        "--classes", str(classes), "--output-dir", str(out), *extra]) == 0
    return pd.read_csv(out / "atlas_duplicates_summary.csv"), out


def test_verdicts():
    with tempfile.TemporaryDirectory() as td:
        summary, _ = run(Path(td))
        verdicts = summary.groupby("class")["verdict"].first().to_dict()
        print("\nverdicts:", verdicts)
        assert verdicts["astrocyte"] == "REGIONAL"
        assert verdicts["fibroblast"] == "MIXED"
        assert verdicts["microglia"] == "PARTIAL"
        assert verdicts["pericyte"] == "TOO FEW CELLS"


def test_home_region_is_not_the_modal_region():
    with tempfile.TemporaryDirectory() as td:
        summary, _ = run(Path(td))
        row = summary.set_index("type").loc["HRCA_microglial cell"]
        print(f"\nHRCA_microglial cell: modal={row['modal_region']} home={row['home_region']} "
              f"concentration={row['home_concentration']:.3f} "
              f"enrichment={row['home_enrichment']:.2f}")
        assert row["modal_region"] == "Gray matter", "most of its cells ARE in gray matter"
        assert row["home_region"] == "Retina", "but retina is where it is over-represented"
        assert row["home_concentration"] < ad_mod.STRONG_CONCENTRATION

        # And the cap the fixed-enrichment bar would have hit: nothing in gray matter can
        # exceed 1/0.70, so the largest region would be unjudgeable on enrichment alone.
        allen = summary.set_index("type").loc["Allen_Astrocyte"]
        assert allen["home_region"] == "Gray matter"
        assert allen["home_enrichment"] < 1.5
        assert allen["home_concentration"] >= ad_mod.STRONG_CONCENTRATION


def test_the_sheet_beats_the_stale_obs_region():
    with tempfile.TemporaryDirectory() as td:
        _, out = run(Path(td))
        table = pd.read_csv(out / "atlas_duplicates_astrocyte.csv", header=[0, 1], index_col=0)
        assert STALE_REGION not in table.index, "obs['Region'] must not reach the output"
        assert set(table.index) == {"Gray matter", "Retina", "Optic nerve"}


def test_excluded_fovs_are_dropped_by_default():
    with tempfile.TemporaryDirectory() as td:
        default, _ = run(Path(td))
        assert default.groupby("class")["verdict"].first()["fibroblast"] == "MIXED"
        assert int(default.set_index("type").loc["Allen_Fibroblast", "n_cells"]) == 1000

    with tempfile.TemporaryDirectory() as td:
        included, _ = run(Path(td), "--include-excluded")
        row = included.set_index("type").loc["Allen_Fibroblast"]
        print(f"\nwith excluded FOVs: Allen_Fibroblast n={row['n_cells']} "
              f"home={row['home_region']} verdict={row['verdict']}")
        assert int(row["n_cells"]) == 3000
        assert row["home_region"] == "Optic nerve"
        assert row["verdict"] != "MIXED", "the excluded FOV alone would flip the call"


if __name__ == "__main__":
    test_verdicts()
    test_home_region_is_not_the_modal_region()
    test_the_sheet_beats_the_stale_obs_region()
    test_excluded_fovs_are_dropped_by_default()
    print("\nall atlas_duplicates_by_region tests passed")
