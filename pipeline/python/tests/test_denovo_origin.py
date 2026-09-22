#!/usr/bin/env python3
"""Tests for the origin diagnostics — denovo_origin, 75k's spatial test, 75m's mixing metric.

Each statistic is tested on the two worlds it has to tell apart, because a statistic that looked
alarming (or quiet) in both would be useless for the question these jobs exist to answer:

  spatial (75k)   a group planted ON vessels must come out enriched; one scattered through
                  parenchyma must not — with the SAME FOVs and the same mural density, so only
                  localisation can explain the difference
  mixing  (75m)   two groups sharing a blob must show no self-preference; a satellite island
                  must show a large diagonal against near-zero off-diagonals

Runnable either under pytest or directly:
    PYTHONPATH=pipeline/python python pipeline/python/tests/test_denovo_origin.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _PY_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


origin = _load("denovo_origin")
spatial = _load("denovo_origin_spatial")
umap = _load("denovo_origin_umap")

FOVS = [f"S1_F{i}" for i in range(6)]
PER_FOV = 400
FOV_SIZE = 4000.0


# --------------------------------------------------------------------------- 75k, spatial


def _fov_layout(rng, n_on_vessel: int, n_scattered: int, n_mural: int):
    """One FOV: mural cells and an on-vessel group along a line, a scattered group off it."""
    def on_line(n):
        t = rng.uniform(0, FOV_SIZE, n)
        return np.column_stack([t, 0.6 * FOV_SIZE - 0.35 * t]) + rng.normal(0, 55.0, (n, 2))

    def anywhere(n):
        return rng.uniform(0, FOV_SIZE, (n, 2))

    xy = np.vstack([on_line(n_mural), on_line(n_on_vessel), anywhere(n_scattered),
                    anywhere(PER_FOV - n_mural - n_on_vessel - n_scattered)])
    kind = (["mural"] * n_mural + ["on_vessel"] * n_on_vessel
            + ["scattered"] * n_scattered + ["filler"] * (len(xy) - n_mural - n_on_vessel
                                                          - n_scattered))
    return xy, kind


def _spatial_cells() -> pd.DataFrame:
    rng = np.random.default_rng(3)
    blocks, kinds, fov_of = [], [], []
    for fov in FOVS:
        xy, kind = _fov_layout(rng, n_on_vessel=60, n_scattered=60, n_mural=70)
        blocks.append(xy)
        kinds += kind
        fov_of += [fov] * len(xy)
    xy = np.vstack(blocks)
    return pd.DataFrame({"slide": "S1", "fov": fov_of, "x": xy[:, 0], "y": xy[:, 1],
                         "kind": kinds, "is_mural": [k == "mural" for k in kinds]})


def _ratio(cells, pools, fov_key, kind, rng):
    ids = cells.index[cells["kind"] == kind]
    observed = float(np.nanmean(cells.loc[ids, "mural_frac"].to_numpy()))
    null = spatial.null_distribution(pools, fov_key.loc[ids].value_counts(), 200, rng)
    return observed / float(np.nanmean(null))


def _spatial_ratios():
    cells = _spatial_cells()
    cells["mural_frac"] = spatial.mural_fraction(cells, k=15)
    fov_key = cells["fov"]
    pools = spatial.fov_pools(cells["mural_frac"], fov_key)
    rng = np.random.default_rng(11)
    return (_ratio(cells, pools, fov_key, "on_vessel", rng),
            _ratio(cells, pools, fov_key, "scattered", rng))


def test_on_vessel_group_is_called_perivascular():
    on_vessel, _ = _spatial_ratios()
    assert on_vessel > 1.4, f"cells planted on vessels must be enriched: ratio {on_vessel:.2f}"


def test_scattered_group_is_not_called_perivascular():
    _, scattered = _spatial_ratios()
    assert scattered < 1.1, f"cells scattered off vessels must not be: ratio {scattered:.2f}"


def test_the_two_spatial_worlds_are_far_apart():
    """Marginal separation would be useless — the whole point is a decisive call."""
    on_vessel, scattered = _spatial_ratios()
    assert on_vessel > 2 * scattered, (on_vessel, scattered)


def test_null_draws_are_size_matched_and_without_replacement():
    """A draw the size of its whole pool must reproduce the pool mean exactly, every time."""
    pool = np.arange(50, dtype=float)
    pools = {"S1_F0": pool}
    got = spatial.null_distribution(pools, pd.Series({"S1_F0": 50}), 25,
                                    np.random.default_rng(0))
    assert np.allclose(got, pool.mean()), got


def test_null_ignores_fovs_the_group_never_occupies():
    """The null is drawn from the group's OWN FOVs; a vessel-free FOV must not dilute it."""
    pools = {"busy": np.full(100, 0.5), "empty": np.zeros(100)}
    got = spatial.null_distribution(pools, pd.Series({"busy": 20}), 20,
                                    np.random.default_rng(0))
    assert np.allclose(got, 0.5), got


def test_slide_fov_falls_back_to_the_cell_id():
    """obs without slide_id/fov still resolves, via `<slide>_F<fov>_C<cell>`."""
    obs = pd.DataFrame(index=["S1_F3_C1", "S1_F3_C2", "S2_F11_C9"])
    got = spatial.resolve_slide_fov(obs)
    assert list(got["slide"]) == ["S1", "S1", "S2"], got
    assert list(got["fov"]) == [3, 3, 11], got


def _clique_cells() -> pd.DataFrame:
    """A group that is tightly clustered AND self-labelled mural, but NOT on the vessels.

    This is the `c` situation: 72% of c's cells are Pericyte under the fixed-profile run, so
    they sit in the anchor set and count toward each other's neighbourhoods. A group like that
    scores as perivascular from its own members alone.
    """
    rng = np.random.default_rng(5)
    blocks, kinds, fov_of = [], [], []
    for fov in FOVS:
        xy, kind = _fov_layout(rng, n_on_vessel=0, n_scattered=0, n_mural=70)
        # the clique: a tight blob in a corner, far from the vessel line
        clique = rng.normal(0, 60.0, (60, 2)) + np.array([3400.0, 3400.0])
        blocks += [xy, clique]
        kinds += kind + ["clique"] * len(clique)
        fov_of += [fov] * (len(xy) + len(clique))
    xy = np.vstack(blocks)
    cells = pd.DataFrame({"slide": "S1", "fov": fov_of, "x": xy[:, 0], "y": xy[:, 1],
                          "kind": kinds})
    # BOTH the real mural cells and the clique carry an anchor type.
    cells["is_mural"] = cells["kind"].isin(["mural", "clique"])
    return cells


def _clique_ratio(exclude: bool) -> float:
    cells = _clique_cells()
    reference = None
    if exclude:
        reference = pd.Series(cells["kind"] != "clique", index=cells.index)
    cells["mural_frac"] = spatial.mural_fraction(cells, k=15, reference=reference)
    fov_key = cells["fov"]
    pools = spatial.fov_pools(cells["mural_frac"], fov_key)
    return _ratio(cells, pools, fov_key, "clique", np.random.default_rng(13))


def test_a_self_mural_clique_looks_perivascular_without_the_control():
    """The circularity is real — without holding it out, the clique scores as vessel-bound."""
    got = _clique_ratio(exclude=False)
    assert got > 2.0, f"expected the uncontrolled statistic to be inflated, got {got:.2f}"


def test_exclude_compared_removes_the_self_supplied_evidence():
    """Held out of the reference, the clique is judged only by cells that are not under test."""
    got = _clique_ratio(exclude=True)
    assert got < 1.0, f"a clique away from the vessels must not look perivascular: {got:.2f}"


def test_the_control_changes_the_verdict_not_just_the_number():
    uncontrolled, controlled = _clique_ratio(False), _clique_ratio(True)
    assert uncontrolled > 3 * controlled, (uncontrolled, controlled)


def test_reference_none_matches_the_original_self_excluding_behaviour():
    """With no reference every cell is its own query AND a neighbour; the self-hit must go."""
    cells = pd.DataFrame({"slide": "S1", "fov": "S1_F0",
                          "x": np.arange(60.0), "y": np.zeros(60),
                          "is_mural": [True] * 30 + [False] * 30})
    frac = spatial.mural_fraction(cells, k=4)
    # Cell 0 sits at the mural end: all four nearest others are mural, and it must not count
    # ITSELF among them (which would be indistinguishable here, so check the far end too).
    assert frac.iloc[0] == 1.0, frac.head()
    assert frac.iloc[-1] == 0.0, frac.tail()


class _FovArgs:
    n_example_fovs = 3
    min_native_cells = 10
    mural_quantiles = "0.25,0.75"


def _fov_selection_inputs():
    """Five FOVs spanning the mural range, plus one that is legible but has no native cells."""
    spec = [("wall_to_wall", 0.60, 40), ("dense", 0.40, 40), ("mid", 0.20, 40),
            ("sparse", 0.10, 40), ("vessel_free", 0.02, 40), ("no_natives", 0.20, 0)]
    rows, tidy_rows = [], []
    for fov, mural_share, n_native in spec:
        n = 500
        for i in range(n):
            rows.append({"cell": f"{fov}_C{i}", "fov": fov,
                         "is_mural": i < int(n * mural_share)})
        for i in range(n_native):
            tidy_rows.append({"cell": f"{fov}_C{i + 200}", "group": "D [native]",
                              "origin": origin.NATIVE_ORIGIN})
        for i in range(60):
            tidy_rows.append({"cell": f"{fov}_C{i + 300}", "group": "l->D",
                              "origin": origin.FORCED_ORIGIN})
    cells = pd.DataFrame(rows).set_index("cell")
    tidy = pd.DataFrame(tidy_rows).set_index("cell")
    return cells, tidy, cells["fov"]


def test_fov_selection_requires_native_cells():
    cells, tidy, fov_key = _fov_selection_inputs()
    chosen = spatial.pick_example_fovs(cells, tidy, fov_key, _FovArgs())
    assert "no_natives" not in chosen, chosen


def test_fov_selection_skips_the_extremes_of_mural_density():
    """The old rule picked the densest FOVs, where mural is wall-to-wall and nothing shows."""
    cells, tidy, fov_key = _fov_selection_inputs()
    chosen = spatial.pick_example_fovs(cells, tidy, fov_key, _FovArgs())
    assert "wall_to_wall" not in chosen, chosen
    assert "vessel_free" not in chosen, chosen
    assert "mid" in chosen, chosen


def test_fov_selection_falls_back_loudly_rather_than_drawing_nothing():
    cells, tidy, fov_key = _fov_selection_inputs()

    class Strict(_FovArgs):
        min_native_cells = 10_000

    chosen = spatial.pick_example_fovs(cells, tidy, fov_key, Strict())
    assert len(chosen) == Strict.n_example_fovs, chosen


# --------------------------------------------------------------------------- 75m, mixing


def _enrichment(centres: dict, n_per_group: int = 300):
    rng = np.random.default_rng(5)
    order = list(centres)
    groups = np.repeat(order, n_per_group)
    xy = np.vstack([np.array(centres[g]) + rng.normal(0, 1.0, (n_per_group, 2))
                    for g in order])
    _, idx = cKDTree(xy).query(xy, k=16)
    rows = np.repeat(np.arange(len(xy)), 15)
    graph = sp.csr_matrix((np.ones(rows.size, dtype=np.float32), (rows, idx[:, 1:].ravel())),
                          shape=(len(xy), len(xy)))
    return umap.neighbour_enrichment(graph, pd.Series(groups), order)


def test_intermingled_groups_show_no_self_preference():
    """Two groups in one blob: each sees the other as often as itself."""
    enr = _enrichment({"A": (0, 0), "B": (0, 0)})
    assert abs(enr.loc["A", "B"] - enr.loc["A", "A"]) < 0.2, enr
    assert 0.8 < enr.loc["A", "B"] < 1.2, enr


def test_a_satellite_island_is_caught():
    """A group on its own island: large diagonal, near-zero against everything else."""
    enr = _enrichment({"A": (0, 0), "B": (0, 0), "island": (60, 60)})
    assert enr.loc["island", "island"] > 2.5, enr
    assert enr.loc["island", "A"] < 0.05, enr
    assert enr.loc["A", "island"] < 0.05, enr


def test_enrichment_rows_are_normalised_against_group_size():
    """A big group must not look enriched just for being big."""
    rng = np.random.default_rng(2)
    groups = np.array(["big"] * 900 + ["small"] * 100)
    xy = rng.normal(0, 1.0, (1000, 2))          # one blob: everything is at chance
    _, idx = cKDTree(xy).query(xy, k=16)
    rows = np.repeat(np.arange(1000), 15)
    graph = sp.csr_matrix((np.ones(rows.size, dtype=np.float32), (rows, idx[:, 1:].ravel())),
                          shape=(1000, 1000))
    enr = umap.neighbour_enrichment(graph, pd.Series(groups), ["big", "small"])
    assert np.allclose(enr.to_numpy(), 1.0, atol=0.25), enr


# --------------------------------------------------------------------------- the shared split


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _origin_frame(tmp: Path, no_all: bool = False):
    ids = [f"c{i}" for i in range(10)]
    semisup = ["l", "l", "l", "Endo_capilar", "Endo_capilar", "Pericyte", "l", "q", "q", "q"]
    forced = ["Endo_capilar", "Endo_capilar", "Pericyte", "Endo_capilar", "Endo_capilar",
              "Pericyte", "MES-like_hypoxia_MHC", "Endo_capilar", "Pericyte", "Pericyte"]
    import h5py
    with h5py.File(tmp / "typing.h5", "w") as f:
        dt = h5py.special_dtype(vlen=str)
        f.create_dataset("cell_id", data=np.array(ids, dtype=object), dtype=dt)
        f.create_dataset("cell_type", data=np.array(semisup, dtype=object), dtype=dt)
    pd.DataFrame({"cell_id": ids, "top1_type": forced}).to_csv(tmp / "forced.csv", index=False)
    return origin.load_origin(_Args(
        typing_h5=tmp / "typing.h5", forced_csv=tmp / "forced.csv", letter="l",
        destinations="Endo_capilar,Pericyte", crosstab=None, top_destinations=3,
        min_destination_pct=1.0, no_all_column=no_all))


def test_origin_split_assigns_forced_native_and_all():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        tidy, order = _origin_frame(Path(d))
    counts = tidy["group"].value_counts()
    assert counts["l->Endo_capilar"] == 2, counts        # c0, c1
    assert counts["Endo_capilar [native]"] == 2, counts  # c3, c4
    assert counts["l->Pericyte"] == 1, counts            # c2
    assert counts["Pericyte [native]"] == 1, counts      # c5
    # c6 is an `l` cell whose forced call is neither destination: [all] only.
    assert counts["l [all]"] == 4, counts
    assert order[0] == "l->Endo_capilar" and order[1] == "Endo_capilar [native]", order


def test_origin_marks_destination_and_origin():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        tidy, _ = _origin_frame(Path(d))
    forced = tidy[tidy["origin"] == origin.FORCED_ORIGIN]
    assert set(forced["destination"]) == {"Endo_capilar", "Pericyte"}, forced
    native = tidy[tidy["origin"] == origin.NATIVE_ORIGIN]
    assert set(native["destination"]) == {"Endo_capilar", "Pericyte"}, native
    assert set(tidy[tidy["origin"] == origin.ALL_ORIGIN]["group"]) == {"l [all]"}


def test_no_all_column_leaves_each_cell_once():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        tidy, order = _origin_frame(Path(d), no_all=True)
    assert tidy.index.is_unique, tidy
    assert not any(g.endswith("[all]") for g in order), order


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
