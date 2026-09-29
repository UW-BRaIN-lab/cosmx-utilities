#!/usr/bin/env python3
"""Validate merged Napari `_metadata.csv` files for the anchor cohort, with an INDEPENDENT join.

The merge scripts join through the AtoMx flat files, so checking their output with the same join
would only confirm the join agrees with itself. This rebuilds each cell's pipeline id straight
from the Napari id and looks the labels up directly:

    Napari id   c_<slide#>_<FOV>_<cell>          e.g.  c_8_1_12
    pipeline id <slide folder>_F<FOV>_C<cell>    e.g.  7347G17347G3_F1_C12

so the flat files, the crosswalk and the merge script are all bypassed for the labels. Per slide
it checks that:

  * the expected columns exist, each label column with its `<name>_color` twin;
  * rows and `cell_ID` order equal the stitched original (--original-dir);
  * `gbmap_type` and `leiden` equal a direct lookup, for EVERY cell (untyped cells must be blank);
  * `gbmap_denovo` equals `gbmap_type` where that is a de novo letter and is blank elsewhere;
  * no cell carries a Leiden or letter label without a GBmap call (a stray label);
  * the number of typed cells equals the anchor cells on that slide (--anchor-cells-csv);
  * Region / Case equal a direct lookup in the curated FOV sheet (--annotations, --crosswalk);
  * every label has exactly one colour, and it is the shared hash colour.

WHEN A MISMATCH COUNT EQUALS THE NUMBER OF TYPED CELLS EXACTLY, SUSPECT THE CHECK, NOT THE DATA.
Every typed cell disagreeing is what a missing label prefix looks like -- it happened here, when
Napari's `leiden_9` was compared against the raw `9`. The report says so when it sees that pattern.

Exit status is non-zero if any slide fails, so it can gate an upload.

Usage:
    uv run python scripts/validate-napari-anchor-metadata.py \\
        --merged-dir build/step4 --original-dir build/stitched \\
        --typing-h5 anchor_typing.h5 --leiden-csv anchor_cohort_leiden.csv \\
        --anchor-cells-csv anchor_cells.csv \\
        --annotations pipeline/reference/gbm_fov_annotations.csv \\
        --crosswalk pipeline/reference/gbm_slide_name_crosswalk.csv
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

METADATA_SUFFIX = "_metadata.csv"
COLOR_SUFFIX = "_color"
CELL_ID_COLUMN = "cell_ID"
NAPARI_ID_RE = re.compile(r"^c_\d+_(?P<fov>\d+)_(?P<cell>\d+)$")
DENOVO_RE = re.compile(r"^[a-z]{1,2}$")     # InSituType's letters; K=27 overflows to "aa"
ANNOTATION_COLUMNS = ("Region", "Case")


def _load_writer():
    """The annotation script owns the colour hash and the FOV-sheet loader; import them rather
    than re-implement, so the validator cannot drift from the writer."""
    path = Path(__file__).with_name("add-annotations-to-napari-metadata.py")
    spec = importlib.util.spec_from_file_location("add_annotations_writer", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged-dir", type=Path, required=True,
                   help="Directory of the final <slide>_metadata.csv to validate.")
    p.add_argument("--original-dir", type=Path,
                   help="Directory of the stitched originals (same file names): rows and cell_ID "
                        "order must be unchanged.")
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="anchor_typing.h5 (/cell_id, /cell_type): the k=27 fit.")
    p.add_argument("--leiden-csv", type=Path,
                   help="Cohort-wide Leiden per anchor cell (cell_id, leiden).")
    p.add_argument("--leiden-prefix", default="leiden_",
                   help="Prefix the merge put in front of each Leiden value (default leiden_).")
    p.add_argument("--anchor-cells-csv", type=Path,
                   help="anchor_cells.csv (cell_id, slide_id): typed cells per slide must equal it.")
    p.add_argument("--annotations", type=Path, help="AtoMx annotation reference (Region/Case check).")
    p.add_argument("--crosswalk", type=Path, help="Slide-name crosswalk, needed with --annotations.")
    p.add_argument("--gbmap-column", default="gbmap_type")
    p.add_argument("--denovo-column", default="gbmap_denovo")
    p.add_argument("--leiden-column", default="leiden")
    return p.parse_args()


def _decode(arr) -> np.ndarray:
    vals = np.asarray(arr[()])
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def load_typing(path: Path) -> pd.Series:
    with h5py.File(path, "r") as f:
        return pd.Series(_decode(f["cell_type"]), index=_decode(f["cell_id"]), name="cell_type")


def load_leiden(path: Path, prefix: str) -> pd.Series:
    df = pd.read_csv(path, dtype=str)
    if not {"cell_id", "leiden"} <= set(df.columns):
        sys.exit(f"ERROR: {path} needs cell_id and leiden columns; has {list(df.columns)}")
    return (prefix + df["leiden"]).set_axis(df["cell_id"]).rename("leiden")


def hint_if_all_wrong(mismatches: int, typed: int, what: str) -> str:
    """A mismatch count equal to the typed count is the fingerprint of a broken check."""
    if typed > 0 and mismatches == typed:
        return (f"{what}: {mismatches:,} mismatches == the typed-cell count exactly. That is the "
                f"signature of a bug in the CHECK (e.g. a missing label prefix), not of bad data.")
    return ""


def validate_slide(slide: str, d: pd.DataFrame, original: pd.DataFrame | None, typing: pd.Series,
                   leiden: pd.Series | None, anchor_n: int | None, fovs: pd.DataFrame | None,
                   args, hash_color) -> tuple[list[str], dict]:
    """Return (problems, stats) for one slide. `d` is read with keep_default_na=False."""
    problems: list[str] = []
    gb, dn, ld = args.gbmap_column, args.denovo_column, args.leiden_column
    required = [CELL_ID_COLUMN, gb, gb + COLOR_SUFFIX, dn, dn + COLOR_SUFFIX]
    if leiden is not None:
        required += [ld, ld + COLOR_SUFFIX]
    missing = [c for c in required if c not in d.columns]
    if missing:
        return [f"missing columns {missing}"], {}

    if original is not None:
        if len(d) != len(original):
            problems.append(f"rows {len(d):,} != original {len(original):,}")
        elif not (d[CELL_ID_COLUMN].values == original[CELL_ID_COLUMN].values).all():
            problems.append("cell_ID order differs from the stitched original")

    # THE INDEPENDENT JOIN: pipeline id from the Napari id, labels looked up directly.
    ids = d[CELL_ID_COLUMN].str.extract(NAPARI_ID_RE)
    bad_ids = int(ids["fov"].isna().sum())
    if bad_ids:
        problems.append(f"{bad_ids:,} cell ids do not look like c_<slide>_<fov>_<cell>")
    pid = slide + "_F" + ids["fov"] + "_C" + ids["cell"]
    typed = int((d[gb] != "").sum())

    gb_mism = int((d[gb] != pid.map(typing).fillna("")).sum())
    if gb_mism:
        problems.append(f"{gb}: {gb_mism:,} cells differ from a direct lookup")
        h = hint_if_all_wrong(gb_mism, typed, gb)
        if h:
            problems.append(h)
    ld_mism = 0
    if leiden is not None:
        ld_mism = int((d[ld] != pid.map(leiden).fillna("")).sum())
        if ld_mism:
            problems.append(f"{ld}: {ld_mism:,} cells differ from a direct lookup")
            h = hint_if_all_wrong(ld_mism, typed, ld)
            if h:
                problems.append(h)

    is_letter = d[gb].map(lambda v: bool(DENOVO_RE.fullmatch(v)) if v else False)
    if not (d[dn] == d[gb].where(is_letter, "")).all():
        problems.append(f"{dn} is not '{gb} where a de novo letter, blank elsewhere'")
    if leiden is not None:
        stray = int(((d[gb] == "") & ((d[ld] != "") | (d[dn] != ""))).sum())
        if stray:
            problems.append(f"{stray:,} cells carry a {ld}/{dn} label with no {gb} call")

    if anchor_n is not None and typed != anchor_n:
        problems.append(f"typed cells {typed:,} != anchor cells on this slide {anchor_n:,}")

    reg_mism = 0
    if fovs is not None:
        known = slide in fovs.index.get_level_values(0)
        for col in ANNOTATION_COLUMNS:
            if col not in d.columns:
                continue
            if not known:
                problems.append(f"slide {slide} is not in the annotation reference")
                break
            want = ids["fov"].astype("Int64").map(fovs.loc[slide][col]).astype(object)
            want = want.where(want.notna(), "").astype(str)
            n = int((d[col] != want).sum())
            reg_mism += n
            if n:
                problems.append(f"{col}: {n:,} cells differ from a direct lookup in the FOV sheet")

    # Colours: one per label, and it must be the shared hash.
    for col in [c for c in (gb, dn, ld, *ANNOTATION_COLUMNS) if c in d.columns
                and c + COLOR_SUFFIX in d.columns]:
        sub = d[d[col] != ""]
        for label, colours in sub.groupby(col)[col + COLOR_SUFFIX]:
            if colours.nunique() != 1:
                problems.append(f"{col}={label!r} has {colours.nunique()} different colours")
            elif colours.iloc[0] != hash_color(label):
                problems.append(f"{col}={label!r} colour {colours.iloc[0]} != shared hash "
                                f"{hash_color(label)}")

    return problems, {"rows": len(d), "typed": typed, "gb": gb_mism, "leiden": ld_mism,
                      "region": reg_mism}


def main() -> None:
    args = parse_args()
    if args.annotations and not args.crosswalk:
        sys.exit("ERROR: --annotations needs --crosswalk (canonical slide name -> folder).")
    writer = _load_writer()
    hash_cache: dict[str, str] = {}

    def hash_color(label: str) -> str:
        if label not in hash_cache:
            hash_cache[label] = writer.deterministic_color(label)
        return hash_cache[label]

    typing = load_typing(args.typing_h5)
    leiden = load_leiden(args.leiden_csv, args.leiden_prefix) if args.leiden_csv else None
    anchor_counts = None
    if args.anchor_cells_csv:
        ac = pd.read_csv(args.anchor_cells_csv, usecols=["slide_id"])
        anchor_counts = ac.groupby("slide_id").size()
    fovs = (writer.load_fov_table(args.annotations, args.crosswalk, list(ANNOTATION_COLUMNS))
            if args.annotations else None)

    files = sorted(args.merged_dir.glob(f"*{METADATA_SUFFIX}"))
    if not files:
        sys.exit(f"ERROR: no *{METADATA_SUFFIX} files in {args.merged_dir}")
    print(f"{len(files)} slide file(s); typing {len(typing):,} cells"
          + (f", Leiden {len(leiden):,}" if leiden is not None else "") + "\n")
    print(f"  {'slide':22s} {'rows':>9s} {'typed':>8s} {'gbmap✗':>7s} {'leiden✗':>8s} {'region✗':>8s}")

    failed: dict[str, list[str]] = {}
    for path in files:
        slide = path.name[: -len(METADATA_SUFFIX)]
        d = pd.read_csv(path, dtype=str, keep_default_na=False)
        original = None
        if args.original_dir:
            opath = args.original_dir / path.name
            if not opath.exists():
                failed[slide] = [f"no stitched original at {opath}"]
                print(f"  {slide:22s} FAIL (no original)")
                continue
            original = pd.read_csv(opath, dtype=str, keep_default_na=False)
        anchor_n = None
        if anchor_counts is not None:
            anchor_n = int(anchor_counts.get(slide, 0))
        problems, st = validate_slide(slide, d, original, typing, leiden, anchor_n, fovs, args,
                                      hash_color)
        if problems:
            failed[slide] = problems
        if st:
            print(f"  {slide:22s} {st['rows']:>9,d} {st['typed']:>8,d} {st['gb']:>7d} "
                  f"{st['leiden']:>8d} {st['region']:>8d}   {'FAIL' if problems else 'ok'}")
        else:
            print(f"  {slide:22s} FAIL")

    if failed:
        print(f"\n{len(failed)} of {len(files)} slide(s) FAILED:")
        for slide, problems in failed.items():
            print(f"  {slide}")
            for prob in problems:
                print(f"      - {prob}")
        sys.exit(1)
    print(f"\nALL {len(files)} SLIDES PASS EVERY CHECK")


if __name__ == "__main__":
    main()
