#!/usr/bin/env python3
"""Tests for donor_depth_within_slide.py.

Plants (a) a clear within-slide carrier gap that a pooled comparison would also see, and
(b) a pure SLIDE effect: carriers live on the shallow slide, with no gap within it. The pooled
difference is large in (b); the within-slide statistic must be ~0 and non-significant, which is
the whole point of the script.

Run:  uv run --with pandas --with numpy python pipeline/python/tests/test_donor_depth_within_slide.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from donor_depth_within_slide import exact_permutation_p, usable_slides  # noqa: E402


def _table(rows):
    return pd.DataFrame(rows, columns=["donor", "slide_id", "carrier", "shallow_frac"])


def test_only_slides_with_both_groups_are_used():
    t = _table([("a", "S1", True, .9), ("b", "S1", False, .5), ("c", "S2", False, .5),
                ("d", "S2", False, .6), ("e", "S3", True, .8), ("f", "S3", True, .7)])
    print(usable_slides(t))
    assert usable_slides(t) == ["S1"]


def test_clear_within_slide_gap_is_detected():
    rows = []
    for s in ("S1", "S2", "S3"):
        rows += [(f"{s}c", s, True, 0.9), (f"{s}n1", s, False, 0.5), (f"{s}n2", s, False, 0.5)]
    observed, p, n = exact_permutation_p(_table(rows), "shallow_frac")
    print(f"gap={observed:+.3f} p={p:.4f} arrangements={n}")
    assert abs(observed - 0.4) < 1e-9 and n == 27
    assert p <= 1 / 27 * 2 + 1e-9     # the observed arrangement and its mirror at most


def test_pure_slide_effect_gives_no_within_slide_gap():
    # carriers all sit on the shallow slide S1 (0.9) but are no different from S1's non-carriers;
    # S2 is deep with no carriers (dropped as unusable by the caller). Pooled gap is huge.
    rows = [("c1", "S1", True, 0.9), ("c2", "S1", True, 0.9), ("n1", "S1", False, 0.9),
            ("n2", "S1", False, 0.9), ("n3", "S2", False, 0.3), ("n4", "S2", False, 0.3)]
    t = _table(rows)
    pooled = t[t.carrier].shallow_frac.mean() - t[~t.carrier].shallow_frac.mean()
    used = t[t.slide_id.isin(usable_slides(t))].reset_index(drop=True)
    observed, p, _ = exact_permutation_p(used, "shallow_frac")
    print(f"pooled gap={pooled:+.3f}  within-slide gap={observed:+.3f} p={p:.3f}")
    assert pooled > 0.25 and abs(observed) < 1e-9 and p == 1.0


if __name__ == "__main__":
    test_only_slides_with_both_groups_are_used()
    test_clear_within_slide_gap_is_detected()
    test_pure_slide_effect_gives_no_within_slide_gap()
    print("OK")
