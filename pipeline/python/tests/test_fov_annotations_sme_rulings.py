#!/usr/bin/env python3
"""The committed FOV sheet must still say what the neuropathologist ruled.

Her answers are the only source for these FOVs -- there is nothing in the data to re-derive
them from -- so an edit that quietly reverts one is unrecoverable without asking her again.
This pins the rulings of 2026-11-16 and the invariants of the mixture column.

    uv run python pipeline/python/tests/test_fov_annotations_sme_rulings.py
"""
import sys
from pathlib import Path

import pandas as pd

SHEET = Path(__file__).resolve().parents[2] / "reference" / "fov_annotations.csv"

SOFT, RETINA, ON, WM, GM = ("Adjacent soft tissue", "Retina", "Optic nerve",
                            "White matter", "Gray matter")
VOCABULARY = {SOFT, RETINA, ON, WM, GM, "Choroid", "EXCLUDE_not_in_atlas", "Unlocated"}

# (slide, fov) -> (primary region, mixed flag, mixture or "")
RULINGS = {
    ("uwa7634eyes", 135): (SOFT, 0, ""),
    ("uwa7634eyes", 178): (WM, 0, ""),
    ("uwa7634eyes", 179): (WM, 0, ""),
    ("uwa7509eyes", 94): (ON, 1, f"{ON} + {SOFT}"),
    ("uwa7509eyes", 104): (ON, 1, f"{SOFT} + {RETINA} + {ON}"),
    ("uwa7509eyes", 217): (WM, 1, f"{WM} + {GM}"),
    ("eyes7597", 67): (ON, 0, ""),
    ("UWA7689eyes", 86): (SOFT, 0, ""),
    ("UWA7689eyes", 94): (SOFT, 1, f"{SOFT} + {RETINA}"),
    ("UWA7689eyes", 102): (SOFT, 0, ""),
    ("uwa7753eyes", 113): (SOFT, 1, f"{SOFT} + {RETINA}"),
    ("uwa7753eyes", 117): (SOFT, 1, f"{SOFT} + {RETINA} + {ON}"),
    ("uwa7753eyes", 132): (ON, 1, f"{SOFT} + {RETINA} + {ON}"),
    ("uwa7753eyes", 196): (WM, 0, ""),
    ("UWA7740eyes", 37): (ON, 0, ""),
    ("UWA7740eyes", 46): (ON, 1, f"{ON} + {SOFT}"),
    ("UWA7697eyes", 45): (RETINA, 0, ""),
    ("UWA7697eyes", 46): (RETINA, 1, f"{RETINA} + {SOFT} + {ON}"),
    ("UWA7697eyes", 134): ("Unlocated", 0, ""),
    ("eyes7597", 57): (ON, 0, ""),
    ("eyes7597", 206): (WM, 0, ""),
    ("uwa7611eyes", 82): (ON, 1, f"{ON} + {SOFT}"),
    ("uwa7509eyes", 73): (SOFT, 1, f"{SOFT} + {ON}"),
    ("UWA7575eyes", 49): (ON, 1, f"{ON} + {SOFT}"),
    ("eyes7517", 81): (ON, 1, f"{ON} + {SOFT}"),
}
# The 46 FOVs that were in no table of the deck, by slide.
PREVIOUSLY_UNRULED = {
    "UWA7689eyes": list(range(86, 103)),
    "uwa7753eyes": [110, 111, 112, 113, 116, 117, 193, 194, 196, 197],
    "UWA7740eyes": [36, 37, 38, 39, 40, 43, 44, 45, 46],
    "UWA7697eyes": [45, 46, 69, 134],
    "eyes7597": [57, 206],
    "uwa7611eyes": [82], "uwa7509eyes": [73], "UWA7575eyes": [49], "eyes7517": [81],
}


def load() -> pd.DataFrame:
    sheet = pd.read_csv(SHEET, dtype={"slide_id": str}).fillna({"region_mixture": ""})
    for col in ("slide_id", "fov", "region", "mixed_adjacent", "exclude", "region_mixture"):
        assert col in sheet.columns, f"{SHEET} has no {col!r} column"
    return sheet


def test_every_ruling_is_in_the_sheet():
    sheet = load().set_index(["slide_id", "fov"])
    for key, (region, mixed, mixture) in RULINGS.items():
        assert key in sheet.index, f"{key} is missing from the sheet"
        row = sheet.loc[key]
        got = (row["region"], int(row["mixed_adjacent"]), row["region_mixture"])
        assert got == (region, mixed, mixture), f"{key}: {got} != {(region, mixed, mixture)}"
    print(f"  {len(RULINGS)} rulings intact")


def test_the_46_unruled_fovs_all_have_a_row_now():
    sheet = load()
    keyed = set(zip(sheet.slide_id, sheet.fov))
    missing = [(s, f) for s, fovs in PREVIOUSLY_UNRULED.items()
               for f in fovs if (s, f) not in keyed]
    assert not missing, f"still unruled: {missing}"
    total = sum(len(v) for v in PREVIOUSLY_UNRULED.values())
    print(f"  all {total} formerly unruled FOVs are annotated")


def test_the_primary_region_is_always_a_member_of_its_mixture():
    """The mixture is extra information about the FOV, never a contradiction of its label."""
    sheet = load()
    mixed = sheet[sheet["region_mixture"] != ""]
    bad = [(r.slide_id, r.fov, r.region, r.region_mixture) for r in mixed.itertuples()
           if r.region not in [t.strip() for t in r.region_mixture.split("+")]]
    assert not bad, f"primary region absent from its own mixture: {bad}"
    assert (mixed["mixed_adjacent"] == 1).all(), "a recorded mixture implies mixed_adjacent=1"
    print(f"  {len(mixed)} FOVs carry a mixture, every primary a member of it")


def test_vocabulary_and_uniqueness():
    sheet = load()
    unknown = set(sheet["region"]) - VOCABULARY
    assert not unknown, f"unknown region value(s): {unknown}"
    dupes = sheet[sheet.duplicated(["slide_id", "fov"], keep=False)]
    assert dupes.empty, dupes[["slide_id", "fov"]].to_dict("records")
    # The one FOV she could not find must stay excluded, or its cells rejoin the cohort
    # carrying a region nobody verified.
    unlocated = sheet[sheet["region"] == "Unlocated"]
    assert len(unlocated) == 1 and int(unlocated.iloc[0]["exclude"]) == 1
    print(f"  {len(sheet)} rows, {len(set(sheet.region))} region values, "
          f"{int((sheet.exclude == 1).sum())} excluded")


if __name__ == "__main__":
    for test in (test_every_ruling_is_in_the_sheet,
                 test_the_46_unruled_fovs_all_have_a_row_now,
                 test_the_primary_region_is_always_a_member_of_its_mixture,
                 test_vocabulary_and_uniqueness):
        print(test.__name__)
        test()
    print("\nall SME-ruling tests passed")
