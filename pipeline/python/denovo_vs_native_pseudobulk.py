#!/usr/bin/env python3
"""Marker heatmap inputs: a de-novo letter's cells, split by the GBmap class they are FORCED
into, against the cells that got that same GBmap class NATIVELY.

The PI's question, in her words: for the cells that were "different enough" to be put into a
de-novo cluster when the semi-supervised fit allowed one, how do they compare to the cells that
received a GBmap call without needing a de-novo bin? Take de-novo `t`. The forced supervised run
scatters it across OPC-like (32%), AC-like (30%) and MES-like_hypoxia_MHC (19%). So build, for
each destination D:

    t->D        cells with semisup == t   AND forced == D
    D [native]  cells with semisup == D   (the fit named them; no de-novo bin needed)

and put the pair side by side in one heatmap. Two things are then readable at once:

  1. DOWN each pair — how t->D differs from native D. Note this difference is guaranteed by
     construction: a cell landed in `t` precisely because it fit t's profile better than any
     GBmap profile. So "they differ" is definitional, NOT a finding. What is informative is
     WHICH genes differ, and whether they form a coherent programme (hypoxia, proliferation,
     stress) or look like scatter.
  2. ACROSS the t->D columns — whether the letter's cells that got sent to different GBmap
     classes actually differ from EACH OTHER. If t->OPC-like and t->AC-like are
     transcriptionally the same, the supervised split of `t` is arbitrary and `t` is one
     population. If they differ, the forced call is recovering real substructure inside `t`.
     This is usually the more interesting axis, and it is why every destination is drawn in one
     figure rather than one figure per destination.

PINNED CELLS. The `D [native]` group is not a clean "the likelihood chose D" population:
insitutype() derives its own best-exemplar cells and overwrites their labels (75o measured
30.6% of named cells pinned, up to 97.8% of a native group). `--pinned-csv` takes 75o's
pinned_cells.csv and either SPLITS each native group into `D [native, pinned]` and
`D [native, unpinned]` (--pinned-mode split, the default) or DROPS the pinned cells
(--pinned-mode exclude). Either way marker_amplitude.csv is written: for each destination,
how strongly each group carries D's own markers, in one unit across groups, so "the l->D cells
carry D's markers at half the strength" can be checked against natives that were not selected
as exemplars. Letter cells are never pinned (pinned cells always carry a named label).

Two further readouts answer "is the resemblance to the destination real?":

  * sibling_vs_destination.csv — for each forced subset, its Pearson r (over the marker genes, on
    the z-scored matrix) to its SIBLINGS (the letter's other forced subsets) against its r to the
    destination's native cells, split pinned/unpinned. Siblings, not the `[all]` column: the whole
    letter contains the subset, so r to it is inflated. The z-scoring across k columns gives a
    null r of -1/(k-1), reported as `null_r` so a small negative r is not read as anti-correlation.
  * marker_amplitude_genes.csv — the per-gene contribution behind each amplitude, and with
    --exclude-genes (e.g. heat-shock) the amplitude recomputed without those genes.

A `<letter> [all]` column is included by default so each subset can also be read against the
letter as a whole.

Inputs:
  --counts-h5     anchor_input.h5 (stage 4a): /counts CSC genes x cells, /genes, /cell_id.
  --typing-h5     anchor_typing.h5 — semi-supervised label per cell (/cell_id, /cell_type).
  --forced-csv    75c's forced_named_posteriors.csv (cell_id, top1_type) — the forced GBmap
                  class, read off the fit's own stored logliks (see 75c).
  --letter        the de-novo letter to dissect, e.g. t
  --destinations  comma-separated GBmap classes to compare against, e.g.
                  OPC-like,AC-like,MES-like_hypoxia_MHC. Or omit it and pass --crosstab to
                  take the letter's own largest destinations automatically, which is what you
                  want when sweeping several letters.

Writes exactly what R/marker_heatmap.R reads in its no-region mode:
  <out>/marker_heatmap_zmatrix.csv    genes x groups, z-scored
  <out>/top_markers_per_cluster.csv   gene, cluster (row split)
  <out>/group_sizes.csv               group, n_cells, dropped
  <out>/marker_amplitude.csv          destination, group, n_cells, amplitude, marker_set
  <out>/marker_amplitude_genes.csv    destination, group, gene, contribution, is_excluded
  <out>/sibling_vs_destination.csv    per forced subset: r to siblings vs r to native destination
  <out>/column_correlations.csv       full column x column Pearson r on the marker z-matrix

Usage:
    uv run python pipeline/python/denovo_vs_native_pseudobulk.py \\
        --counts-h5 anchor_input.h5 --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv \\
        --letter t --destinations OPC-like,AC-like,MES-like_hypoxia_MHC \\
        --output-dir t_vs_native
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from anchor_profiles import decode, read_cell_calls
from pseudobulk_core import (DEFAULT_SCALE_FACTOR, group_means, log_normalize, onehot,
                             select_markers, zscore_rows)

NATIVE_SUFFIX = " [native]"
ALL_SUFFIX = " [all]"
FORCED_SEP = "->"
PINNED_MODES = ("split", "exclude")
PINNED_SUFFIX = " [native, pinned]"
UNPINNED_SUFFIX = " [native, unpinned]"
# The heat-shock block the PI page names as t's defining programme (HSPA1A, HSPA1B, HSPB1, DNAJB1,
# HSP90AA1, HSPH1). --exclude-genes heat-shock expands to this.
HEAT_SHOCK_GENES = ("HSPA1A", "HSPA1B", "HSPB1", "DNAJB1", "HSP90AA1", "HSPH1")
HEAT_SHOCK_KEYWORD = "heat-shock"
GENE_DETAIL_N = 30


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-h5", type=Path, required=True,
                   help="anchor_input.h5: /counts (CSC genes x cells), /genes, /cell_id.")
    p.add_argument("--typing-h5", type=Path, required=True,
                   help="anchor_typing.h5: semi-supervised /cell_type per /cell_id.")
    p.add_argument("--forced-csv", type=Path, required=True,
                   help="forced_named_posteriors.csv from 75c (cell_id, top1_type).")
    p.add_argument("--letter", required=True, help="De-novo letter to dissect, e.g. t")
    p.add_argument("--destinations",
                   help="Comma-separated GBmap classes to compare against. Omit to derive "
                        "them from --crosstab.")
    p.add_argument("--crosstab", type=Path,
                   help="denovo_vs_gbmap_crosstab.csv — derive the destinations from this "
                        "letter's own largest ones instead of listing them by hand.")
    p.add_argument("--top-destinations", type=int, default=3,
                   help="With --crosstab, how many destinations to take (default 3).")
    p.add_argument("--min-destination-pct", type=float, default=5.0,
                   help="With --crosstab, ignore destinations below this %% of the letter "
                        "(default 5); keeps thin columns out of a swept figure.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=8,
                   help="Top markers selected per group (default 8).")
    p.add_argument("--min-group-n", type=int, default=50,
                   help="Drop groups with fewer cells than this (default 50).")
    p.add_argument("--scale-factor", type=float, default=DEFAULT_SCALE_FACTOR)
    p.add_argument("--pinned-csv", type=Path,
                   help="75o's pinned_cells.csv (cell_id, ...): InSituType's pinned cells.")
    p.add_argument("--pinned-mode", choices=PINNED_MODES, default="split",
                   help="With --pinned-csv: split each native group into pinned/unpinned "
                        "(default) or exclude the pinned cells.")
    p.add_argument("--exclude-genes",
                   help="Comma-separated genes (or the keyword 'heat-shock') to leave out of the "
                        "marker set used for marker_amplitude.csv. Amplitude is then reported for "
                        "both marker sets (marker_set = all / excluding_listed), so a gap that "
                        "shrinks when these genes go is one they were carrying.")
    p.add_argument("--gene-detail-n", type=int, default=GENE_DETAIL_N,
                   help="Genes per destination in marker_amplitude_genes.csv (default 30).")
    p.add_argument("--no-all-column", action="store_true",
                   help="Omit the '<letter> [all]' whole-letter reference column.")
    return p.parse_args()


def read_counts(path: Path) -> tuple[sp.csc_matrix, np.ndarray, np.ndarray]:
    """Read the stage-4a CSC genes x cells counts, returning (matrix, genes, cell_ids)."""
    with h5py.File(path, "r") as f:
        for key in ("counts/shape", "counts/data", "counts/indices", "counts/indptr",
                    "genes", "cell_id"):
            if key not in f:
                sys.exit(f"ERROR: /{key} missing in {path}.")
        n_genes, n_cells = (int(x) for x in f["counts/shape"][()])
        mat = sp.csc_matrix(
            (np.asarray(f["counts/data"][()], dtype=np.float64),
             np.asarray(f["counts/indices"][()], dtype=np.int64),
             np.asarray(f["counts/indptr"][()], dtype=np.int64)),
            shape=(n_genes, n_cells))
        genes = np.asarray(decode(f["genes"][()]))
        cell_id = np.asarray(decode(f["cell_id"][()]))
    return mat, genes, cell_id


def destinations_from_crosstab(path: Path, letter: str, top_n: int, min_pct: float) -> list[str]:
    """The letter's largest forced destinations, as a fallback for --destinations.

    Rows of the cross-tab carry the readable display label ("t - AC-like"), so the letter is
    matched on the token before the first " - ".
    """
    ct = pd.read_csv(path, index_col=0)
    match = [i for i in ct.index if str(i).split(" - ")[0].strip() == letter]
    if not match:
        sys.exit(f"ERROR: letter {letter!r} not among the cross-tab rows: {list(ct.index)[:5]}...")
    row = ct.loc[match[0]]
    pct = row / row.sum() * 100
    keep = pct[pct >= min_pct].sort_values(ascending=False).head(top_n)
    if keep.empty:
        sys.exit(f"ERROR: no destination of {letter!r} reaches {min_pct}%.")
    print(f"Destinations for {letter} from {path.name}: "
          + ", ".join(f"{d} ({v:.1f}%)" for d, v in keep.items()))
    return list(keep.index)


def build_groups(semisup: pd.Series, forced: pd.Series, letter: str,
                 destinations: list[str], include_all: bool) -> pd.Series:
    """Per-cell group label, or NaN for cells in none of the compared groups.

    A cell belongs to at most one group: `<letter>-><D>` when the fit called it `letter` and the
    forced run sends it to D, or `<D> [native]` when the fit named it D outright. Cells of the
    letter that went to some other destination are unlabelled (they are not part of this
    comparison) unless the whole-letter reference column is on, which additionally tags every
    `letter` cell.
    """
    group = pd.Series(pd.NA, index=semisup.index, dtype="object")
    is_letter = semisup == letter
    for dest in destinations:
        group[is_letter & (forced == dest)] = f"{letter}{FORCED_SEP}{dest}"
        group[semisup == dest] = f"{dest}{NATIVE_SUFFIX}"
    # The whole-letter column reuses the letter's cells, so it is built as a separate pass
    # (a cell can contribute to both its <letter>-><D> column and the <letter> [all] column).
    return group, is_letter if include_all else None


def split_natives_by_pinned(group: pd.Series, pinned: set[str], mode: str) -> pd.Series:
    """Re-label native groups by whether InSituType pinned the cell.

    Only `<D> [native]` groups are touched. `split` turns each into `<D> [native, pinned]` and
    `<D> [native, unpinned]`; `exclude` drops the pinned cells from the comparison (NaN).
    """
    if mode not in PINNED_MODES:
        raise ValueError(f"mode must be one of {PINNED_MODES}, got {mode!r}")
    is_native = group.str.endswith(NATIVE_SUFFIX, na=False).to_numpy()
    is_pinned = group.index.isin(pinned)
    out = group.copy()
    if mode == "exclude":
        out[is_native & is_pinned] = pd.NA
        return out
    base = group.str[:-len(NATIVE_SUFFIX)]
    out[is_native & is_pinned] = base[is_native & is_pinned] + PINNED_SUFFIX
    out[is_native & ~is_pinned] = base[is_native & ~is_pinned] + UNPINNED_SUFFIX
    return out


def resolve_excluded_genes(spec: str | None) -> set[str]:
    """--exclude-genes value -> gene set; 'heat-shock' expands to the page's six genes."""
    if not spec:
        return set()
    genes: set[str] = set()
    for token in (t.strip() for t in spec.split(",") if t.strip()):
        genes.update(HEAT_SHOCK_GENES if token.lower() == HEAT_SHOCK_KEYWORD else [token])
    return genes


def marker_amplitude(profile: pd.DataFrame, pooled_native: pd.DataFrame,
                     groups_by_dest: dict[str, list[str]], sizes: pd.Series,
                     top_n: int, min_group_n: int, exclude_genes: frozenset[str] = frozenset(),
                     detail_n: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """How strongly each group carries destination D's own markers, in one unit.

    Markers of D = the top_n genes by (pooled native D - mean of the OTHER destinations' pooled
    natives), where "pooled" is every cell the fit named D, pinned or not, so the marker set is
    the same however the natives are split. A group's amplitude = mean over those genes of
    (its mean log-norm expression - the same baseline). All groups of one destination share the
    gene set and the baseline, so their amplitudes are directly comparable.

    `exclude_genes` are removed BEFORE the top_n are taken, so the set refills with the next-best
    markers rather than shrinking. Returns (amplitude table, per-gene table); the per-gene table
    covers the top `detail_n` markers of the UNfiltered set, with `is_excluded` marking the genes
    that would be removed, and is empty when detail_n is 0.
    """
    rows, gene_rows = [], []
    for dest, groups in groups_by_dest.items():
        others = [o for o in pooled_native.columns if o != dest]
        if not others:
            continue
        baseline = pooled_native[others].mean(axis=1)
        ranked = (pooled_native[dest] - baseline).sort_values(ascending=False)
        markers = ranked[~ranked.index.isin(exclude_genes)].index[:top_n]
        candidates = {f"{dest} [native, all]": pooled_native[dest]}
        candidates.update({g: profile[g] for g in groups if g in profile.columns})
        for name, col in candidates.items():
            n = int(sizes.get(name, 0)) if name in sizes.index else np.nan
            rows.append({"destination": dest, "group": name, "n_cells": n,
                         "amplitude": float((col[markers] - baseline[markers]).mean()),
                         "below_min_group_n": bool(n < min_group_n) if n == n else False,
                         "marker_set": "excluding_listed" if exclude_genes else "all"})
            if detail_n:
                for rank, gene in enumerate(ranked.index[:detail_n], start=1):
                    gene_rows.append({"destination": dest, "group": name, "gene": gene,
                                      "marker_rank": rank,
                                      "contribution": float(col[gene] - baseline[gene]),
                                      "is_excluded": gene in exclude_genes})
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def sibling_vs_destination(pb_z: pd.DataFrame, letter: str, destinations: list[str]
                           ) -> pd.DataFrame:
    """Each forced subset's r to its siblings vs its r to the destination's native cells.

    r is Pearson over the marker genes on the z-scored matrix. Siblings are the letter's OTHER
    forced subsets; the `[all]` column is deliberately not used because it contains the subset.
    A z-score across k columns gives a null r of -1/(k-1), reported as `null_r`.
    A positive `siblings_minus_native` means the subset resembles its siblings more than the
    cells natively called its destination.
    """
    k = pb_z.shape[1]
    forced_cols = [f"{letter}{FORCED_SEP}{d}" for d in destinations
                   if f"{letter}{FORCED_SEP}{d}" in pb_z.columns]
    all_col = f"{letter}{ALL_SUFFIX}"
    rows = []
    for dest in destinations:
        col = f"{letter}{FORCED_SEP}{dest}"
        if col not in pb_z.columns:
            continue
        siblings = [c for c in forced_cols if c != col]
        r_sib = [float(pb_z[col].corr(pb_z[c])) for c in siblings]
        for variant in (NATIVE_SUFFIX, UNPINNED_SUFFIX, PINNED_SUFFIX):
            ref = f"{dest}{variant}"
            if ref not in pb_z.columns:
                continue
            r_native = float(pb_z[col].corr(pb_z[ref]))
            rows.append({
                "forced": col, "native_ref": ref, "n_siblings": len(siblings),
                "r_siblings_mean": float(np.mean(r_sib)) if r_sib else np.nan,
                "r_siblings_max": float(np.max(r_sib)) if r_sib else np.nan,
                "r_native": r_native,
                "siblings_minus_native": (float(np.mean(r_sib)) - r_native) if r_sib else np.nan,
                "r_letter_all_contains_subset": (float(pb_z[col].corr(pb_z[all_col]))
                                                 if all_col in pb_z.columns else np.nan),
                "null_r": -1.0 / (k - 1), "n_columns": k})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if args.destinations:
        destinations = [d.strip() for d in args.destinations.split(",") if d.strip()]
    elif args.crosstab:
        destinations = destinations_from_crosstab(
            args.crosstab, args.letter, args.top_destinations, args.min_destination_pct)
    else:
        sys.exit("ERROR: pass --destinations, or --crosstab to derive them.")
    if not destinations:
        sys.exit("ERROR: no destinations selected.")

    calls = read_cell_calls(args.typing_h5)[["cell_id", "cell_type"]]
    forced = pd.read_csv(args.forced_csv, usecols=["cell_id", "top1_type"])
    labels = calls.merge(forced, on="cell_id", how="inner").set_index("cell_id")
    print(f"Labels joined for {len(labels):,} cells")

    missing = [d for d in destinations if d not in set(labels["top1_type"])]
    if missing:
        sys.exit(f"ERROR: destination(s) never appear as a forced call: {missing}")
    if args.letter not in set(labels["cell_type"]):
        sys.exit(f"ERROR: --letter {args.letter!r} is not a semi-supervised label.")

    group, letter_mask = build_groups(labels["cell_type"], labels["top1_type"],
                                      args.letter, destinations, not args.no_all_column)
    if args.pinned_csv:
        pinned = set(pd.read_csv(args.pinned_csv, usecols=["cell_id"])["cell_id"])
        n_pinned_letter = int(group.index.isin(pinned)[(labels["cell_type"] == args.letter)
                                                       .to_numpy()].sum())
        print(f"Pinned cells: {len(pinned):,} in the file; {n_pinned_letter:,} of them carry "
              f"the letter's label (expected 0)")
        group = split_natives_by_pinned(group, pinned, args.pinned_mode)

    counts, genes, cell_id = read_counts(args.counts_h5)
    print(f"Counts: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")

    # Align the labels to the counts' cell order, then keep only the cells we compare.
    group = group.reindex(cell_id)
    letter_mask = letter_mask.reindex(cell_id).fillna(False).to_numpy() \
        if letter_mask is not None else None
    keep = group.notna().to_numpy()
    if letter_mask is not None:
        keep = keep | letter_mask
    # Every cell the fit named a compared destination is needed for the pooled-native marker
    # baseline, even when --pinned-mode exclude drops it from the groups.
    semisup_by_count = labels["cell_type"].reindex(cell_id)
    native_any = semisup_by_count.isin(destinations).to_numpy()
    keep = keep | native_any
    if not keep.any():
        sys.exit("ERROR: no cells matched the requested letter/destinations.")
    print(f"Comparing {int(keep.sum()):,} cells")

    # Subset columns on the CSC (cells are columns -> cheap) BEFORE transposing to cells x genes.
    norm = log_normalize(counts[:, keep].T.tocsr(), args.scale_factor)
    kept_group = group[keep]
    kept_letter = letter_mask[keep] if letter_mask is not None else None
    kept_semisup = semisup_by_count[keep]

    # Mean log-norm profile per group (genes x groups) over the paired groups only; the
    # whole-letter column is appended after, since its cells overlap the <letter>-><D> ones.
    paired = kept_group.notna().to_numpy()
    oh, glabels = onehot(kept_group[paired].to_numpy())
    profile = pd.DataFrame(group_means(norm[paired], oh).T, index=genes, columns=glabels)
    sizes = pd.Series(np.asarray(oh.sum(axis=0)).ravel(), index=glabels, dtype=int)

    all_col = f"{args.letter}{ALL_SUFFIX}"
    if kept_letter is not None and kept_letter.any():
        profile[all_col] = np.asarray(norm[kept_letter].mean(axis=0)).ravel()
        sizes[all_col] = int(kept_letter.sum())

    pooled_native = pd.DataFrame(
        {d: np.asarray(norm[(kept_semisup == d).to_numpy()].mean(axis=0)).ravel()
         for d in destinations if (kept_semisup == d).any()}, index=genes)
    native_variants = (NATIVE_SUFFIX, PINNED_SUFFIX, UNPINNED_SUFFIX)
    groups_by_dest = {d: [f"{args.letter}{FORCED_SEP}{d}"] + [f"{d}{v}" for v in native_variants]
                      for d in pooled_native.columns}
    exclude_genes = frozenset(resolve_excluded_genes(args.exclude_genes))
    absent = sorted(g for g in exclude_genes if g not in set(genes))
    if absent:
        print(f"NOTE: --exclude-genes not on this panel, ignored: {absent}")
    amplitude, amplitude_genes = marker_amplitude(
        profile, pooled_native, groups_by_dest, sizes, args.top_n, args.min_group_n,
        detail_n=args.gene_detail_n)
    if exclude_genes:
        amp_ex, _ = marker_amplitude(profile, pooled_native, groups_by_dest, sizes,
                                     args.top_n, args.min_group_n, exclude_genes=exclude_genes)
        amplitude = pd.concat([amplitude, amp_ex], ignore_index=True)
        amplitude_genes["is_excluded"] = amplitude_genes["gene"].isin(exclude_genes)

    # Interleave each forced subset with its native counterpart so the pair reads together,
    # then the whole-letter reference last.
    order = []
    for dest in destinations:
        for col in (f"{args.letter}{FORCED_SEP}{dest}", *(f"{dest}{v}" for v in native_variants)):
            if col in profile.columns:
                order.append(col)
    if all_col in profile.columns:
        order.append(all_col)
    profile = profile[order]

    small = sizes[sizes < args.min_group_n].index.tolist()
    if small:
        print(f"Dropping {len(small)} group(s) under {args.min_group_n} cells: {small}")
        profile = profile.drop(columns=[c for c in small if c in profile.columns])
    if profile.shape[1] < 2:
        sys.exit("ERROR: fewer than 2 groups survive; nothing to compare.")

    print("Group sizes:")
    for col in profile.columns:
        print(f"  {col:44s} {sizes[col]:>9,}")

    markers, gene_to_group = select_markers(profile, list(profile.columns), args.top_n)
    pb_z = zscore_rows(profile.loc[markers])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pb_z.to_csv(args.output_dir / "marker_heatmap_zmatrix.csv")
    pd.DataFrame({"cluster": [gene_to_group[g] for g in markers],
                  "gene": markers}).to_csv(
        args.output_dir / "top_markers_per_cluster.csv", index=False)
    sizes.rename("n_cells").rename_axis("group").reset_index().assign(
        dropped=lambda d: d["group"].isin(small)).to_csv(
        args.output_dir / "group_sizes.csv", index=False)

    sib = sibling_vs_destination(pb_z, args.letter, destinations)
    if not sib.empty:
        sib.to_csv(args.output_dir / "sibling_vs_destination.csv", index=False)
        print("\nSiblings vs destination (r over the marker genes; null r = "
              f"{sib['null_r'].iloc[0]:.3f}):")
        print(sib.drop(columns=["n_columns"]).to_string(index=False,
                                                        float_format=lambda x: f"{x:.3f}"))
    pb_z.corr().to_csv(args.output_dir / "column_correlations.csv")
    if not amplitude_genes.empty:
        amplitude_genes.to_csv(args.output_dir / "marker_amplitude_genes.csv", index=False)
    if not amplitude.empty:
        amplitude.to_csv(args.output_dir / "marker_amplitude.csv", index=False)
        print("\nMarker amplitude (each destination's own markers; comparable within a "
              "destination):")
        print(amplitude.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print(f"\nWrote {args.output_dir} — {pb_z.shape[0]} markers x {pb_z.shape[1]} groups.")
    print("Render on the Mac (the container has no ComplexHeatmap):")
    print(f"  Rscript pipeline/R/marker_heatmap.R {args.output_dir} <out_dir> \"\" "
          f"\"{args.letter} forced into GBmap classes vs cells natively called those classes\"")


if __name__ == "__main__":
    main()
