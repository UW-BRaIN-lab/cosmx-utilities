#!/usr/bin/env python3
"""Tests for celltype_by_group.py (no S3 / network).

Plants a known structure: type X is 10% of cells in Control donors and 30% in SORL1 donors, type Y is the reverse, type Z is
flat. Pins: row percentages sum to 100, enrichment is observed/expected (a type that mirrors group size is ~1), the donor level
uses the donor as the unit (a donor with 10x the cells must not dominate), Kruskal-Wallis finds X and Y but not Z, BH q is
monotone and correct on a known vector, and the crosswalk follows cells to the other labeling.

Run:  uv run --with pandas --with numpy --with scipy python pipeline/python/tests/test_celltype_by_group.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import celltype_by_group as cg  # noqa: E402


def _table(seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for g, n_donors, px, py in (("Control", 6, .10, .30), ("AD Pure", 6, .15, .25), ("AD+LATE", 6, .20, .20),
                                ("AD+LATE SORL1", 6, .30, .10)):
        for d in range(n_donors):
            n = 3000 if d else 30000                      # the first donor of every group has 10x the cells
            z = 0.2
            probs = [px, py, z, 1 - px - py - z]
            labels = rng.choice(["X", "Y", "Z", "W"], size=n, p=probs)
            rows += [(f"{g}_{d}_{i}", f"{g}_{d}", g, lab, "coarse_" + ("A" if lab in ("X", "Y") else "B"))
                     for i, lab in enumerate(labels)]
    return pd.DataFrame(rows, columns=["cell_id", "donor", "case_group", "fine", "coarse"])


def test_bh_matches_a_hand_computed_case():
    q = cg.benjamini_hochberg(np.array([0.01, 0.04, 0.03, 0.20]))
    print(q)
    assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.20], atol=1e-6)


def test_percentages_enrichment_and_donor_unit():
    t = _table()
    counts, pct, enr = cg.pct_and_enrichment(t, "fine", cg.GROUPS)
    assert np.allclose(pct[cg.GROUPS].sum(axis=1), 100, atol=0.05)
    # the big donor makes cell-level and donor-level answers differ: cell level is pulled to each group's first donor
    dl = cg.donor_level(t, "fine", cg.GROUPS)
    print(dl[[f"mean_pct_{g}" for g in cg.GROUPS] + ["q"]].round(3))
    assert dl.loc["X", "mean_pct_Control"] < 15 and dl.loc["X", "mean_pct_AD+LATE SORL1"] > 25    # planted 10% -> 30%
    assert dl.loc["Y", "mean_pct_Control"] > 25 and dl.loc["Y", "mean_pct_AD+LATE SORL1"] < 15    # planted 30% -> 10%
    assert dl.loc["X", "q"] < 0.1 and dl.loc["Y", "q"] < 0.1 and dl.loc["Z", "q"] > 0.1           # Z is flat
    # trend across the ordered Control < AD Pure < AD+LATE axis (planted X 10->15->20%, Y 30->25->20%, Z flat)
    assert dl.loc["X", "trend_rho"] > 0.5 and dl.loc["Y", "trend_rho"] < -0.5, dl[["trend_rho", "trend_q"]]
    assert dl.loc["X", "trend_q"] < 0.05 and dl.loc["Z", "trend_q"] > 0.1


def test_crosswalk_follows_cells_and_reports_the_destination_result():
    t = _table()
    f_dl = cg.donor_level(t, "fine", cg.GROUPS); c_dl = cg.donor_level(t, "coarse", cg.GROUPS)
    cw = cg.crosswalk(t, "fine", "coarse", f_dl, c_dl, cg.GROUPS)
    print(cw[["cells", "top_destination", "top_destination_share", "first_q", "dest_q"]])
    assert cw.loc["X", "top_destination"] == "coarse_A" and cw.loc["Z", "top_destination"] == "coarse_B"
    assert cw.loc["X", "top_destination_share"] == 1.0
    # merging X (rises with disease) with Y (falls) cancels: the coarse leaf shows no group difference
    assert cw.loc["X", "first_q"] < 0.1 and cw.loc["X", "dest_q"] > 0.1


def test_workbook_and_figure_are_written_and_readable():
    import tempfile
    t = _table(); results = {}
    for lab in ("fine", "coarse"):
        counts, pct, enr = cg.pct_and_enrichment(t, lab, cg.GROUPS)
        results[lab] = (t, counts, pct, enr, cg.donor_level(t, lab, cg.GROUPS))
    cw = cg.crosswalk(t, "fine", "coarse", results["fine"][4], results["coarse"][4], cg.GROUPS)
    with tempfile.TemporaryDirectory() as d:
        xlsx, png = Path(d) / "o.xlsx", Path(d) / "o.png"
        cg.write_workbook(xlsx, results, ["F", "C"], cw, ["note"])
        cg.plot_enrichment(results, ["fine", "coarse"], png, {"W"})
        sheets = pd.ExcelFile(xlsx).sheet_names
        print(sheets)
        assert sheets == ["README", "F_pct_of_type", "F_enrichment", "F_donor_level", "C_pct_of_type", "C_enrichment",
                          "C_donor_level", "crosswalk"]
        assert png.stat().st_size > 5000


def test_a_group_with_no_donors_left_is_blank_not_an_error():
    """Leaving out the slides that hold a group's donors (the September-only comparison) must not crash."""
    t = _table(); t = t[t.case_group != "AD+LATE SORL1"]
    dl = cg.donor_level(t, "fine", cg.GROUPS)
    print(dl[["mean_pct_AD+LATE SORL1", "p", "trend_rho"]])
    assert dl["mean_pct_AD+LATE SORL1"].isna().all()
    assert dl.loc["X", "q"] < 0.1 and dl.loc["X", "trend_rho"] > 0.5        # the other three groups still compare


def test_named_letters_and_unknown_groups():
    names = Path("/tmp/_names_test.csv"); pd.DataFrame({"denovo_label": ["o"], "annotation": ["o - Oligodendrocyte"]}).to_csv(names, index=False)
    t = pd.DataFrame({"x": ["o", "Endothelial"]})
    assert cg.apply_names(t, "x", names).tolist() == ["o - Oligodendrocyte", "Endothelial"]
    assert cg.apply_names(t, "x", None).tolist() == ["o", "Endothelial"]


if __name__ == "__main__":
    test_bh_matches_a_hand_computed_case()
    test_percentages_enrichment_and_donor_unit()
    test_crosswalk_follows_cells_and_reports_the_destination_result()
    test_workbook_and_figure_are_written_and_readable()
    test_a_group_with_no_donors_left_is_blank_not_an_error()
    test_named_letters_and_unknown_groups()
    print("OK")
