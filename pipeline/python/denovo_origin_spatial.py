#!/usr/bin/env python3
"""75k: are a de-novo letter's cells where their assigned cell type actually lives?

THE POINT OF THIS ONE. Every other diagnostic in the 75 series reads a quantity the fit itself
produced — markers (75d), donors (75h), the likelihood margin (75i). This one reads tissue
coordinates, which the typing never saw. If `l->Endo_capilar` cells are real capillary
endothelium they lie ON VESSELS, and vessels are multicellular tubes: endothelium is wrapped by
pericytes and smooth muscle. A tumour cell shoehorned into an endothelial profile has no reason
to sit next to a pericyte. So the test is a neighbourhood one:

    for each cell, what fraction of its k nearest spatial neighbours are MURAL?

Mural, not endothelial, is the deliberate choice — asking whether cells called endothelial sit
next to cells called endothelial would just re-read the typing back out. Pericytes are a
different population, and adjacency to them is not implied by anything in the assignment.

The neighbourhood is built from the FIXED-PROFILE run's 7.5M cells (cosmx_typed.h5ad), not the
2.54M anchor: the anchor is a Leiden-STRATIFIED subsample with a cap per stratum, so its local
density is a function of how abundant each cluster was, and neighbourhood composition measured
on it would be an artefact of the sampling. The full typed run is the real tissue.

THE NULL. "40% of its neighbours are mural" means nothing without a yardstick, because a FOV
that happens to be vessel-rich raises that number for every cell in it. So each group is
compared against size-matched random cells drawn FROM THE SAME FOVs, which holds both tissue
density and FOV composition fixed and leaves only the group's own localisation.

Inputs:
  --typed-h5ad  cosmx_typed.h5ad — obs only (backed): centroids, fov, cell_type.
  the shared origin flags (--typing-h5, --forced-csv, --letter, --destinations/--crosstab)
Writes (--output-dir):
  spatial_neighbourhood.csv   per group: observed mural fraction, null mean/sd, z, ratio
  spatial_neighbourhood.png   observed vs null per group
  fov_<slide>_F<fov>.png      example FOVs, cells coloured by origin over the mural layer

Usage:
    python pipeline/python/denovo_origin_spatial.py \\
        --typed-h5ad cosmx_typed.h5ad --typing-h5 anchor_typing.h5 \\
        --forced-csv forced_named_posteriors.csv --letter l \\
        --destinations Endo_capilar,Endo_arterial,Pericyte --output-dir l_spatial
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from denovo_origin import (ALL_ORIGIN, NATIVE_ORIGIN, NATIVE_SUFFIX, add_origin_arguments,
                           load_origin)
from prep_insitucnv_input import pick_spatial_cols

# Vessel-wall cells. Endothelial types are deliberately absent — see the module docstring.
DEFAULT_ANCHOR_TYPES = "Pericyte,SMC,SMC_COL,Perivascular_fibroblast,Scavenging_pericyte"
CELL_ID_PATTERN = r"^(?P<slide>.+)_F(?P<fov>\d+)_C\d+$"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--typed-h5ad", type=Path, required=True,
                   help="cosmx_typed.h5ad — read backed, obs only.")
    add_origin_arguments(p)
    p.add_argument("--celltype-key", default="cell_type",
                   help="obs column holding the fixed-profile label (default cell_type).")
    p.add_argument("--anchor-types", default=DEFAULT_ANCHOR_TYPES,
                   help=f"Comma-separated types defining the vessel wall (default "
                        f"{DEFAULT_ANCHOR_TYPES}).")
    p.add_argument("--k", type=int, default=15,
                   help="Spatial neighbours per cell (default 15).")
    p.add_argument("--exclude-compared", action="store_true",
                   help="Hold every compared cell out of the neighbourhood REFERENCE. Required "
                        "whenever the letter's own cells carry an --anchor-type under the "
                        "fixed-profile run (c is ~72%% Pericyte), or each group partly supplies "
                        "its own evidence. The overlap is reported either way.")
    p.add_argument("--n-permutations", type=int, default=200,
                   help="Size-matched within-FOV draws for the null (default 200).")
    p.add_argument("--n-example-fovs", type=int, default=6,
                   help="FOVs to draw (default 6). See --min-native-cells for how they are "
                        "chosen — NOT simply the ones with the most forced cells.")
    p.add_argument("--min-native-cells", type=int, default=10,
                   help="A FOV must hold at least this many NATIVE cells to be drawn (default "
                        "10). Native cells are the scarce side (a few thousand across 57 "
                        "slides), so a panel without them shows nothing to compare.")
    p.add_argument("--mural-quantiles", default="0.25,0.75",
                   help="Keep only FOVs whose mural share falls between these quantiles OF THE "
                        "ELIGIBLE FOVs (default 0.25,0.75). Vessel architecture is legible in "
                        "the middle of the range: a FOV that is mural wall-to-wall has no "
                        "contrast, and one with no vessels has nothing to sit on.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def resolve_slide_fov(obs: pd.DataFrame) -> pd.DataFrame:
    """slide + fov per cell, from obs when it carries them and from the cell id otherwise.

    Cell ids are `<slide>_F<fov>_C<cell>`. A slide holds TWO tissue sections, but they sit far
    apart in global-px, so a k-nearest-neighbour query inside a slide does not cross the gap.
    """
    out = pd.DataFrame(index=obs.index)
    parsed = obs.index.to_series().str.extract(CELL_ID_PATTERN)
    if "slide_id" in obs:
        out["slide"] = obs["slide_id"].astype(str)
    else:
        out["slide"] = parsed["slide"]
    if "fov" in obs:
        out["fov"] = pd.to_numeric(obs["fov"], errors="coerce")
    else:
        out["fov"] = pd.to_numeric(parsed["fov"], errors="coerce")
    bad = out["slide"].isna() | out["fov"].isna()
    if bad.all():
        sys.exit("ERROR: could not resolve slide/fov from obs or from the cell-id pattern "
                 f"{CELL_ID_PATTERN!r}. obs columns: {list(obs.columns)[:20]}")
    if bad.any():
        print(f"  ({int(bad.sum()):,} cells dropped: no resolvable slide/fov)")
    return out[~bad]


def mural_fraction(cells: pd.DataFrame, k: int, reference: pd.Series | None = None) -> pd.Series:
    """Fraction of each cell's k nearest neighbours (same slide) that are mural.

    `reference` selects which cells may SERVE AS neighbours. Every cell is still a query point —
    the null needs the same statistic for the background cells — but with the compared cells
    held out of the reference the answer cannot be built from the cells under test.

    One KD-tree per slide, queried for every cell on that slide, then the index array is thrown
    away — only the float32 fraction is kept, so peak memory stays at one slide's worth.
    """
    frac = pd.Series(np.nan, index=cells.index, dtype="float32")
    radius = []
    for slide, block in cells.groupby("slide", observed=True, sort=False):
        in_ref = (np.ones(len(block), dtype=bool) if reference is None
                  else reference.reindex(block.index).fillna(False).to_numpy())
        ref = block[in_ref]
        if len(ref) <= k + 1:
            continue
        tree = cKDTree(ref[["x", "y"]].to_numpy(dtype=np.float64))
        # k + 1 neighbours, then drop one: for a query point that is ITSELF in the reference the
        # first hit is the point itself, and for one that is not there is no self-hit to drop —
        # taking the same column from both would either keep a self-match or discard a real
        # neighbour, and the two cases coexist whenever the compared cells are held out.
        dist, idx = tree.query(block[["x", "y"]].to_numpy(dtype=np.float64), k=k + 1, workers=-1)
        take = np.where(in_ref[:, None], idx[:, 1:], idx[:, :-1])
        edge = np.where(in_ref, dist[:, -1], dist[:, -2])
        is_mural = ref["is_mural"].to_numpy()
        frac.loc[block.index] = is_mural[take].mean(axis=1).astype(np.float32)
        radius.append(np.median(edge))
    if radius:
        print(f"Median radius of the {k}-neighbour disc: {np.median(radius):,.0f} px "
              f"(a CosMx pixel is ~0.12 um, a cell ~10 um across)")
    return frac


def fov_pools(frac: pd.Series, fov_key: pd.Series) -> dict:
    """Per-FOV array of mural fractions — the draw pool for the null, built ONCE.

    Built by a single groupby rather than per group: a boolean scan per FOV would be
    O(FOVs x cells), which at 7.5M cells across a few thousand FOVs is the whole job.
    """
    valid = frac.notna().to_numpy()
    return {fov: block.to_numpy()
            for fov, block in frac[valid].groupby(fov_key[valid], observed=True)}


def null_distribution(pools: dict, member_fovs: pd.Series, n_perm: int,
                      rng: np.random.Generator) -> np.ndarray:
    """Mean mural fraction of size-matched random cells drawn from the group's own FOVs.

    `member_fovs` counts the group's cells per FOV; each draw takes that many cells at random
    from every cell in the same FOV, so FOV composition is held fixed and only the group's
    localisation within those FOVs is left to explain a difference.

    Drawn WITHOUT replacement, one independent permutation per row: a group can be a large
    fraction of its FOV, and sampling with replacement would then overstate the null's spread
    by dropping the finite-population correction.
    """
    totals, weights = np.zeros(n_perm), 0
    for fov, n_members in member_fovs.items():
        pool = pools.get(fov)
        if pool is None or len(pool) == 0:
            continue
        n = int(min(n_members, len(pool)))
        idx = rng.permuted(np.tile(np.arange(len(pool)), (n_perm, 1)), axis=1)[:, :n]
        totals += pool[idx].sum(axis=1)
        weights += n
    return totals / weights if weights else np.full(n_perm, np.nan)


def plot_summary(table: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 0.55 * len(table) + 2.0))
    y = np.arange(len(table))[::-1]
    ax.barh(y, table["null_mean"], color="0.85", height=0.62, label="size-matched null")
    ax.errorbar(table["null_mean"], y, xerr=table["null_sd"], fmt="none",
                ecolor="0.45", elinewidth=1.2, capsize=3)
    ax.scatter(table["observed"], y, s=60, color="firebrick", zorder=3, label="observed")
    ax.set_yticks(y)
    ax.set_yticklabels(table.index)
    ax.set_xlabel("fraction of k nearest spatial neighbours that are mural")
    ax.set_title("Perivascular localisation, by origin")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_fov(block: pd.DataFrame, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(block["x"], block["y"], s=5, color="0.90", rasterized=True, label="all cells")
    mural = block[block["is_mural"]]
    # A light tint, not the near-black this used to be: in a vessel-rich FOV the mural layer is
    # most of the field, and drawn dark it buries the groups it is supposed to be context for.
    ax.scatter(mural["x"], mural["y"], s=9, color="#9fb3c8", rasterized=True, label="mural")
    palette = plt.get_cmap("tab10")
    for i, (grp, sub) in enumerate(block.dropna(subset=["group"]).groupby("group",
                                                                         observed=True)):
        # Native groups are a handful of cells per FOV, so they get a bigger marker and an
        # outline; the forced groups are abundant and stay small.
        native = str(grp).endswith(NATIVE_SUFFIX)
        ax.scatter(sub["x"], sub["y"], s=44 if native else 15, color=palette(i % 10),
                   label=str(grp), rasterized=True,
                   edgecolors="black" if native else "none",
                   linewidths=0.6 if native else 0.0, zorder=3 if native else 2)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10)
    ax.legend(markerscale=1.8, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)


def pick_example_fovs(cells: pd.DataFrame, tidy: pd.DataFrame, fov_key: pd.Series,
                      args) -> list[str]:
    """FOVs chosen for LEGIBILITY, not for holding the most forced cells.

    Ranking by forced-cell count picks the most vessel-dense FOVs in the cohort — precisely the
    ones where mural cells are wall-to-wall and no vessel architecture is visible — and says
    nothing about whether any NATIVE cells are present to compare against. Since the native
    groups are a few thousand cells across 57 slides, such a panel routinely held four or five
    of them, which is not a comparison.

    So: require native cells, keep the middle of the mural-density range, then rank by how many
    native cells there are, because they are the limiting side.
    """
    per_fov = pd.DataFrame({
        "mural_share": cells.groupby(fov_key, observed=True)["is_mural"].mean(),
        "n_cells": cells.groupby(fov_key, observed=True).size()})
    native_ids = tidy.index[tidy["origin"] == NATIVE_ORIGIN].unique()
    forced_ids = tidy.index[(tidy["origin"] != NATIVE_ORIGIN)
                            & (tidy["origin"] != ALL_ORIGIN)].unique()
    per_fov["n_native"] = fov_key.loc[cells.index.intersection(native_ids)].value_counts()
    per_fov["n_forced"] = fov_key.loc[cells.index.intersection(forced_ids)].value_counts()
    per_fov = per_fov.fillna({"n_native": 0, "n_forced": 0})

    eligible = per_fov[(per_fov["n_native"] >= args.min_native_cells)
                       & (per_fov["n_forced"] > 0)]
    if eligible.empty:
        print(f"\nWARNING: no FOV holds {args.min_native_cells} native cells; falling back to "
              f"the FOVs with the most forced cells, which are the least legible ones.",
              file=sys.stderr)
        return list(per_fov.sort_values("n_forced", ascending=False)
                    .head(args.n_example_fovs).index)

    lo_q, hi_q = (float(q) for q in args.mural_quantiles.split(","))
    lo, hi = eligible["mural_share"].quantile([lo_q, hi_q])
    banded = eligible[eligible["mural_share"].between(lo, hi)]
    if banded.empty:
        banded = eligible
    chosen = banded.sort_values("n_native", ascending=False).head(args.n_example_fovs)

    print(f"\nExample FOVs — {len(eligible)} of {len(per_fov)} FOVs hold "
          f">={args.min_native_cells} native cells; keeping mural share in "
          f"[{lo:.1%}, {hi:.1%}] (quantiles {lo_q:g}-{hi_q:g} of those):")
    print(f"  {'fov':<26}{'cells':>7}{'mural':>8}{'forced':>8}{'native':>8}")
    for fov, r in chosen.iterrows():
        print(f"  {fov:<26}{int(r['n_cells']):>7,}{r['mural_share']:>8.1%}"
              f"{int(r['n_forced']):>8,}{int(r['n_native']):>8,}")
    return list(chosen.index)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    tidy, order = load_origin(args)

    print(f"Reading obs from {args.typed_h5ad} (backed)")
    adata = ad.read_h5ad(args.typed_h5ad, backed="r")
    obs = adata.obs
    if args.celltype_key not in obs:
        sys.exit(f"ERROR: obs['{args.celltype_key}'] missing. Present: {list(obs.columns)[:20]}")
    xcol, ycol = pick_spatial_cols(obs.columns)
    print(f"Centroids from ('{xcol}', '{ycol}'); {len(obs):,} cells in the fixed-profile run")

    anchor_types = {t.strip() for t in args.anchor_types.split(",") if t.strip()}
    cells = resolve_slide_fov(obs)
    cells["x"] = pd.to_numeric(obs[xcol], errors="coerce").reindex(cells.index)
    cells["y"] = pd.to_numeric(obs[ycol], errors="coerce").reindex(cells.index)
    cells["is_mural"] = obs[args.celltype_key].astype(str).reindex(cells.index).isin(
        anchor_types).to_numpy()
    cells = cells.dropna(subset=["x", "y"])
    present = anchor_types & set(obs[args.celltype_key].astype(str).unique())
    if not present:
        sys.exit(f"ERROR: none of --anchor-types {sorted(anchor_types)} appear in "
                 f"obs['{args.celltype_key}'].")
    print(f"Mural anchor: {int(cells['is_mural'].sum()):,} cells "
          f"({cells['is_mural'].mean():.2%} of the run) across types {sorted(present)}")

    overlap = tidy.index.unique().intersection(cells.index)
    print(f"cell-id overlap: {len(overlap):,} of {tidy.index.nunique():,} compared cells "
          f"({len(overlap) / max(tidy.index.nunique(), 1):.1%})")
    if len(overlap) == 0:
        sys.exit("ERROR: no compared cell ids are present in the typed run — check the join.")

    # How circular is this? A group whose own cells carry an anchor type under the fixed-profile
    # run partly supplies its own evidence, and the ratio is inflated for BOTH sides of that
    # pair. Report it always, so the caveat is visible even when the flag is off.
    print("\nShare of each group that is ITSELF mural under the fixed-profile run:")
    for grp in order:
        ids = tidy.index[tidy["group"] == grp].unique().intersection(cells.index)
        if len(ids):
            print(f"  {grp:<32} {cells.loc[ids, 'is_mural'].mean():>6.1%}")
    print("  (a large share means the group is part of its own neighbourhood — "
          "re-run with --exclude-compared)")

    reference = None
    if args.exclude_compared:
        reference = pd.Series(True, index=cells.index)
        reference[cells.index.isin(tidy.index.unique())] = False
        print(f"\nHolding {int((~reference).sum()):,} compared cells out of the reference; "
              f"{int(reference.sum()):,} cells remain available as neighbours.")

    print(f"\nComputing the {args.k}-neighbour mural fraction, one KD-tree per slide...")
    cells["mural_frac"] = mural_fraction(cells, args.k, reference)
    fov_key = cells["slide"].astype(str) + "_F" + cells["fov"].astype(int).astype(str)
    pools = fov_pools(cells["mural_frac"], fov_key)

    rows = []
    for grp in order:
        ids = tidy.index[tidy["group"] == grp].unique().intersection(cells.index)
        sub = cells.loc[ids]
        observed = float(np.nanmean(sub["mural_frac"].to_numpy()))
        member_fovs = fov_key.loc[ids].value_counts()
        null = null_distribution(pools, member_fovs, args.n_permutations, rng)
        mu, sd = float(np.nanmean(null)), float(np.nanstd(null))
        rows.append({"group": grp, "n_cells": len(ids), "observed": observed,
                     "null_mean": mu, "null_sd": sd,
                     "z": (observed - mu) / sd if sd > 0 else np.nan,
                     "ratio": observed / mu if mu > 0 else np.nan,
                     "n_fovs": len(member_fovs)})
    table = pd.DataFrame(rows).set_index("group").loc[[r["group"] for r in rows]]
    table.to_csv(args.output_dir / "spatial_neighbourhood.csv")
    plot_summary(table, args.output_dir / "spatial_neighbourhood.png")

    print("\n=== perivascular localisation (mural fraction among k nearest neighbours) ===")
    print(f"{'group':<32}{'n':>9}{'observed':>10}{'null':>9}{'ratio':>8}{'z':>9}")
    for grp, r in table.iterrows():
        print(f"{grp:<32}{int(r['n_cells']):>9,}{r['observed']:>10.3f}{r['null_mean']:>9.3f}"
              f"{r['ratio']:>8.2f}{r['z']:>9.1f}")
    print("\nREAD: ratio > 1 with a large z = the group sits on vessels beyond what its own "
          "FOVs\nexplain. A forced group whose ratio matches its native counterpart's is "
          "localised like it.\nA forced group at ratio ~1 is scattered through parenchyma, "
          "whatever its markers say.")

    chosen = pick_example_fovs(cells, tidy, fov_key, args)
    grouping = tidy[tidy["origin"] != ALL_ORIGIN].groupby(level=0)["group"].first()
    for fov in chosen:
        block = cells[(fov_key == fov).to_numpy()].copy()
        block["group"] = grouping.reindex(block.index)
        plot_fov(block, f"{fov} — {args.letter} and its native counterparts",
                 args.output_dir / f"fov_{fov}.png")
    print(f"\nWrote {len(chosen)} example FOV panels to {args.output_dir}")


if __name__ == "__main__":
    main()
