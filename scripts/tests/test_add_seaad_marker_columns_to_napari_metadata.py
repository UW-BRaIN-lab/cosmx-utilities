#!/usr/bin/env python3
"""Tests for add-seaad-marker-columns-to-napari-metadata.py on synthetic stitched metadata.

Pins the mapping (de novo cluster by annotation text, named type by its own name), that a de novo
cluster and a named type mapping to one SEA-AD type share a colour, that no-support clusters get
"No SEA-AD support", that empty source cells stay empty, and the refusals: an unmapped K=21 value,
a file that already has the columns, a duplicate source value in the table.

Run:  uv run --with pandas --with duckdb python scripts/tests/test_add_seaad_marker_columns_to_napari_metadata.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "add-seaad-marker-columns-to-napari-metadata.py"

EVIDENCE = pd.DataFrame([
    # cluster, kind, support, proposed_seaad_name, current_annotation
    ("o", "de novo", "strong", "SEAAD_Oligodendrocyte", "o - Oligodendrocyte"),
    ("p", "de novo", "strong", "SEAAD_L6b", "p - Neuron (SEMA3E/HS3ST4)"),
    ("s", "de novo", "partial", "SEAAD_L5 IT", "s - Excitatory neuron, L4-like (RORB)"),
    ("i", "de novo", "none", "", "i - Unresolved: low-information (LTBP2/FOXJ1 set)"),
    ("Oligodendrocyte", "named", "strong", "SEAAD_Oligodendrocyte", ""),
    ("L6 neuron", "named", "strong", "SEAAD_Oligodendrocyte", ""),
], columns=["cluster", "kind", "support", "proposed_seaad_name", "current_annotation"])

SOURCE = ["o - Oligodendrocyte", "Oligodendrocyte", "L6 neuron", "i - Unresolved: low-information (LTBP2/FOXJ1 set)",
          "s - Excitatory neuron, L4-like (RORB)", "p - Neuron (SEMA3E/HS3ST4)", ""]


def metadata(source=SOURCE) -> pd.DataFrame:
    return pd.DataFrame({"cell_ID": [f"c_5_1_{i}" for i in range(len(source))],
                         "Donor": ["UWA_747"] * (len(source) - 1) + ["NA"], "Donor_color": ["#6CD090"] * len(source),
                         "InSituType_K21": source})


def run(*args: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPT), *args]
    print("  $", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-1200:], res.stderr[-800:])
    return res


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        meta_dir, out_dir, ev = tmp / "meta", tmp / "out", tmp / "evidence.csv"
        meta_dir.mkdir()
        EVIDENCE.to_csv(ev, index=False)
        metadata().to_csv(meta_dir / "SLIDE_A_metadata.csv", index=False)
        base = ["--metadata-dir", str(meta_dir), "--out-dir", str(out_dir), "--evidence-table", str(ev)]

        print("test_mapping_colours_and_untouched_columns")
        assert run(*base).returncode == 0
        out = pd.read_csv(out_dir / "SLIDE_A_metadata.csv", dtype=str, keep_default_na=False)
        assert list(out["SEAAD_markers_K21"]) == ["Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte",
                                                  "No SEA-AD support", "L5 IT", "L6b", ""], list(out["SEAAD_markers_K21"])
        assert list(out["SEAAD_markers_K21_support"]) == ["strong", "strong", "strong", "none", "partial", "strong", ""]
        colour = out["SEAAD_markers_K21_color"]
        assert colour[0] == colour[1] == colour[2] != "", "same SEA-AD type must share one colour"
        assert colour[3] != colour[0] and colour[6] == "", list(colour)
        assert out["SEAAD_markers_K21_support_color"][0] == out["SEAAD_markers_K21_support_color"][5]
        for untouched in ("cell_ID", "Donor", "Donor_color", "InSituType_K21"):
            assert list(out[untouched]) == list(metadata()[untouched]), untouched
        assert out["Donor"].iloc[-1] == "NA", "text 'NA' in an existing column must survive the round trip"
        assert (meta_dir / "SLIDE_A_metadata.csv").read_text() == metadata().to_csv(index=False), "input modified"

        print("test_unmapped_value_is_refused_not_blanked")
        metadata(SOURCE + ["z - Something new"]).to_csv(meta_dir / "SLIDE_B_metadata.csv", index=False)
        res = run(*base)
        assert res.returncode != 0 and "z - Something new" in (res.stdout + res.stderr)
        (meta_dir / "SLIDE_B_metadata.csv").unlink()

        print("test_file_that_already_has_the_columns_is_refused")
        out.to_csv(meta_dir / "SLIDE_C_metadata.csv", index=False)
        res = run(*base)
        assert res.returncode != 0 and "Refusing to overwrite" in (res.stdout + res.stderr)
        (meta_dir / "SLIDE_C_metadata.csv").unlink()

        print("test_duplicate_source_value_in_table_is_refused")
        dup = pd.concat([EVIDENCE, EVIDENCE.iloc[[4]]])
        dup.to_csv(tmp / "dup.csv", index=False)
        res = run("--metadata-dir", str(meta_dir), "--out-dir", str(tmp / "o2"), "--evidence-table", str(tmp / "dup.csv"))
        assert res.returncode != 0 and "more than once" in (res.stdout + res.stderr)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
