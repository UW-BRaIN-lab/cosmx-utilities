#!/usr/bin/env python3
"""Compare per-cell contrast between CosMx morphology channels.

DAPI reads hazy and unreliable in aged human brain while the Histone marker reads
crisp on the same sections. Three explanations predict different signatures in the
per-cell morphology intensities AtoMx already exports (`Mean.<channel>` and
`Max.<channel>` in `<slide>_metadata_file.csv.gz`), so the existing flat files can
separate them without a new acquisition:

  lipofuscin autofluorescence     Granular, ~0.5-2um puncta, and only in the cells
                                  that carry granules (large neurons). That subset
                                  shows up as a heavy right tail in Max/Mean, so the
                                  signature is `peakedness_tail`, NOT the median
                                  peakedness, which the unaffected majority hides.
                                  Affects every channel, worst in violet, and tracks
                                  donor age.
  fixation autofluorescence       Diffuse. Raises the background floor (`p05`) and
                                  collapses `contrast_index` in EVERY channel, while
                                  PUSHING peakedness DOWN (a constant offset added to
                                  both Max and Mean drives their ratio toward 1).
  DAPI bound to cytoplasmic RNA   DAPI-specific. `area_rho` > 0 (a nuclear stain's
                                  per-cell mean should not grow with cell area) and
                                  `rrna_rho_given_area` high. Histone is an antibody,
                                  so it shares no chemistry with the rRNA probe and
                                  a DAPI-only coupling cannot be explained away.

The headline number is `contrast_index` (p95/p05 of the per-cell mean) and its
Histone-over-DAPI ratio, which turns "Histone segments better" into a number.

Usage:
    uv run python pipeline/python/morphology_channel_contrast.py \
        --flatfiles-dir staged_flatfiles \
        --manifest pipeline/manifest.csv \
        --output morphology_contrast_by_slide.csv \
        --fov-output morphology_contrast_by_fov.csv \
        --figure morphology_contrast.png
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MEAN_PREFIX = "Mean."
MAX_PREFIX = "Max."
AREA_COLUMN = "Area"
FOV_COLUMN = "fov"
METADATA_SUFFIX_PATTERN = re.compile(r"^(?P<slide>.+)_metadata_file\.csv(?:\.gz)?$")

# Role detection is by case-insensitive substring because the MorphologyKit
# `BiologicalTarget` metadata is inconsistent across studies (see the
# --channel-name note in napari-cosmx-fork/utils/stitch_images.py). Every role can
# be pinned explicitly on the command line when detection picks wrong.
NUCLEAR_ALIASES = ("dapi",)
REFERENCE_ALIASES = ("histone", "h3")
RNA_ALIASES = ("rrna", "18s", "ribosom")

# Percentiles defining the contrast index. p05 stands in for the per-FOV background
# floor and p95 for real signal; both are far enough from the tails to survive the
# handful of saturated cells every CosMx FOV carries.
BACKGROUND_PERCENTILE = 5
SIGNAL_PERCENTILE = 95
# Lipofuscin sits in a minority of cells, so its granularity lives in the tail of
# the Max/Mean distribution; the median is deliberately blind to it.
PEAKEDNESS_TAIL_PERCENTILE = 90
# Below this residual variance the control explains essentially all of a variable
# and the partial correlation stops being identified. Left unguarded, three
# near-perfect correlations collapse to r/(1+r) -> 0.5, which reads as a real
# moderate coupling instead of "not answerable from this data".
MIN_RESIDUAL_VARIANCE = 1e-6
MIN_CELLS_PER_GROUP = 50


class MissingChannelError(Exception):
    """A requested morphology channel is absent from a slide's metadata file."""


class MissingColumnError(Exception):
    """A structural column (Area, fov) is absent from a slide's metadata file."""


def discover_channels(columns: list[str]) -> list[str]:
    """Morphology channels carrying BOTH a Mean. and a Max. column, in file order."""
    means = [c[len(MEAN_PREFIX):] for c in columns if c.startswith(MEAN_PREFIX)]
    maxes = {c[len(MAX_PREFIX):] for c in columns if c.startswith(MAX_PREFIX)}
    return [name for name in means if name in maxes]


def resolve_role(channels: list[str], aliases: tuple[str, ...],
                 override: str | None, role: str) -> str | None:
    """Pick the channel for a role, preferring an explicit override."""
    if override:
        if override not in channels:
            raise MissingChannelError(
                f"--{role}-channel {override!r} not among discovered channels: {channels}")
        return override
    for channel in channels:
        lowered = channel.lower()
        if any(alias in lowered for alias in aliases):
            return channel
    return None


def spearman(left: pd.Series, right: pd.Series) -> float:
    """Rank correlation. Ranking then taking Pearson avoids a scipy dependency."""
    if len(left) < MIN_CELLS_PER_GROUP:
        return float("nan")
    ranked_left, ranked_right = left.rank(), right.rank()
    if ranked_left.std() == 0 or ranked_right.std() == 0:
        return float("nan")
    return float(ranked_left.corr(ranked_right))


def partial_spearman(target: pd.Series, other: pd.Series,
                     control: pd.Series) -> float:
    """Rank correlation of target vs other with `control` partialled out.

    Cell area drives both intensities mechanically (a bigger cell integrates more
    of anything), so the raw rRNA coupling is uninterpretable until area is removed.
    """
    r_to = spearman(target, other)
    r_tc = spearman(target, control)
    r_oc = spearman(other, control)
    if any(np.isnan(v) for v in (r_to, r_tc, r_oc)):
        return float("nan")
    residual_target = 1 - r_tc ** 2
    residual_other = 1 - r_oc ** 2
    if residual_target < MIN_RESIDUAL_VARIANCE or residual_other < MIN_RESIDUAL_VARIANCE:
        return float("nan")
    return float((r_to - r_tc * r_oc)
                 / np.sqrt(residual_target * residual_other))


def channel_metrics(frame: pd.DataFrame, channel: str,
                    rna_channel: str | None) -> dict[str, float]:
    """Contrast, peakedness and coupling metrics for one channel over one group."""
    mean_values = frame[MEAN_PREFIX + channel]
    max_values = frame[MAX_PREFIX + channel]
    area = frame[AREA_COLUMN]

    background = float(np.percentile(mean_values, BACKGROUND_PERCENTILE))
    signal = float(np.percentile(mean_values, SIGNAL_PERCENTILE))
    q25, median, q75 = (float(np.percentile(mean_values, q)) for q in (25, 50, 75))

    # A zero background floor means the channel is empty (or the export clipped it);
    # the ratio is undefined rather than infinite, and NaN keeps it out of summaries.
    contrast_index = signal / background if background > 0 else float("nan")
    robust_cv = (q75 - q25) / median if median > 0 else float("nan")

    positive = mean_values > 0
    if positive.any():
        ratio = max_values[positive] / mean_values[positive]
        peakedness = float(ratio.median())
        peakedness_p90 = float(np.percentile(ratio, PEAKEDNESS_TAIL_PERCENTILE))
        # Normalising the tail by the median separates "a minority of cells are
        # granular" from "this whole channel is punctate".
        peakedness_tail = peakedness_p90 / peakedness if peakedness > 0 else float("nan")
    else:
        peakedness = peakedness_p90 = peakedness_tail = float("nan")

    metrics = {
        "n_cells": int(len(frame)),
        "background_p05": background,
        "median": median,
        "signal_p95": signal,
        "contrast_index": contrast_index,
        "robust_cv": robust_cv,
        "peakedness": peakedness,
        "peakedness_p90": peakedness_p90,
        "peakedness_tail": peakedness_tail,
        "area_rho": spearman(mean_values, area),
    }
    if rna_channel and rna_channel != channel:
        rna_values = frame[MEAN_PREFIX + rna_channel]
        metrics["rrna_rho"] = spearman(mean_values, rna_values)
        metrics["rrna_rho_given_area"] = partial_spearman(mean_values, rna_values, area)
    else:
        metrics["rrna_rho"] = float("nan")
        metrics["rrna_rho_given_area"] = float("nan")
    return metrics


def metadata_paths(flatfiles_dir: Path | None,
                   explicit: list[Path]) -> list[tuple[str, Path]]:
    """(slide_id, path) pairs, slide_id taken from the flat-file naming convention."""
    paths = list(explicit)
    if flatfiles_dir:
        paths.extend(sorted(flatfiles_dir.rglob("*_metadata_file.csv*")))

    resolved: list[tuple[str, Path]] = []
    seen: set[Path] = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        match = METADATA_SUFFIX_PATTERN.match(path.name)
        if not match:
            print(f"WARN: skipping {path.name}, not a *_metadata_file.csv[.gz]",
                  file=sys.stderr)
            continue
        resolved.append((match.group("slide"), path))
    return resolved


def analyse_slide(slide_id: str, path: Path, args) -> tuple[list[dict], list[dict]]:
    """Per-slide and per-FOV metric rows for one metadata file."""
    frame = pd.read_csv(path)

    for column in (AREA_COLUMN, FOV_COLUMN):
        if column not in frame.columns:
            raise MissingColumnError(f"{path.name} has no {column!r} column")

    channels = discover_channels(list(frame.columns))
    if not channels:
        raise MissingChannelError(f"{path.name} has no Mean./Max. channel pairs")

    rna_channel = resolve_role(channels, RNA_ALIASES, args.rrna_channel, "rrna")

    # Zero-area cells are segmentation artefacts and make every ratio meaningless.
    # --max-area drops the merged-cell blobs that dominate the upper tail.
    keep = frame[AREA_COLUMN] > 0
    if args.max_area is not None:
        keep &= frame[AREA_COLUMN] <= args.max_area
    dropped = int((~keep).sum())
    if dropped:
        print(f"  {slide_id}: dropped {dropped:,} of {len(frame):,} cells on area")
    frame = frame.loc[keep]

    if len(frame) < MIN_CELLS_PER_GROUP:
        print(f"WARN: {slide_id} has only {len(frame)} usable cells, skipping",
              file=sys.stderr)
        return [], []

    slide_rows = []
    for channel in channels:
        row = {"slide_id": slide_id, "channel": channel}
        row.update(channel_metrics(frame, channel, rna_channel))
        slide_rows.append(row)

    fov_rows = []
    if args.fov_output:
        for fov, group in frame.groupby(FOV_COLUMN, sort=True):
            if len(group) < MIN_CELLS_PER_GROUP:
                continue
            for channel in channels:
                row = {"slide_id": slide_id, "fov": int(fov), "channel": channel}
                row.update(channel_metrics(group, channel, rna_channel))
                fov_rows.append(row)

    return slide_rows, fov_rows


def contrast_ratio_table(slide_frame: pd.DataFrame, nuclear: str,
                         reference: str) -> pd.DataFrame:
    """Per-slide reference-over-nuclear contrast ratio, the headline comparison."""
    wide = slide_frame.pivot(index="slide_id", columns="channel", values="contrast_index")
    missing = [c for c in (nuclear, reference) if c not in wide.columns]
    if missing:
        raise MissingChannelError(f"no contrast_index for channel(s): {missing}")
    ratio = (wide[reference] / wide[nuclear]).rename("reference_over_nuclear")
    return pd.concat([wide[[nuclear, reference]], ratio], axis=1).sort_values(
        "reference_over_nuclear", ascending=False)


def write_figure(slide_frame: pd.DataFrame, ratios: pd.DataFrame,
                 nuclear: str, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    channels = list(slide_frame["channel"].unique())
    panels = [
        ("contrast_index", "Contrast index (p95/p05)\nhigher = crisper"),
        ("peakedness_tail", "Peakedness tail (p90/p50 of Max/Mean)\nhigher = granular subset"),
        ("area_rho", "Intensity vs cell area (rho)\n>0 = cytoplasmic"),
        ("rrna_rho_given_area", "Coupling to rRNA | area (rho)\n>0 = tracks cytoplasmic RNA"),
    ]

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(panels) + 1, figsize=(4.6 * (len(panels) + 1), 4.4))
    for axis, (metric, title) in zip(axes, panels):
        values = [slide_frame.loc[slide_frame["channel"] == c, metric].dropna()
                  for c in channels]
        # Set tick labels separately: boxplot's own keyword was renamed in
        # Matplotlib 3.9, and the container and the Mac are on different versions.
        axis.boxplot(values, showfliers=False)
        axis.set_xticks(range(1, len(channels) + 1))
        axis.set_xticklabels(channels)
        for position, series in enumerate(values, start=1):
            jitter = rng.normal(0, 0.04, len(series))
            axis.plot(position + jitter, series, ".", alpha=0.45, markersize=5)
        axis.set_title(title, fontsize=10)
        axis.tick_params(axis="x", rotation=45)
        if metric in ("area_rho", "rrna_rho_given_area"):
            axis.axhline(0, color="0.4", linewidth=0.8, linestyle="--")

    axis = axes[-1]
    axis.barh(range(len(ratios)), ratios["reference_over_nuclear"], color="#4c72b0")
    axis.axvline(1, color="0.4", linewidth=0.8, linestyle="--")
    axis.set_yticks(range(len(ratios)))
    axis.set_yticklabels(ratios.index, fontsize=6)
    axis.invert_yaxis()
    axis.set_title(f"Contrast ratio over {nuclear}\nper slide", fontsize=10)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--flatfiles-dir", type=Path,
                        help="Directory searched recursively for *_metadata_file.csv[.gz]")
    parser.add_argument("--metadata", type=Path, action="append", default=[],
                        help="Explicit metadata file. Repeatable.")
    parser.add_argument("--manifest", type=Path,
                        help="Pipeline manifest CSV; joins run_date and instrument_id "
                             "onto the per-slide table so batch and instrument effects "
                             "can be separated from tissue effects.")
    parser.add_argument("--nuclear-channel",
                        help="Channel treated as the nuclear stain (default: detect DAPI)")
    parser.add_argument("--reference-channel",
                        help="Channel compared against it (default: detect Histone)")
    parser.add_argument("--rrna-channel",
                        help="rRNA channel used for the coupling test (default: detect)")
    parser.add_argument("--max-area", type=float,
                        help="Drop cells above this area (merged-cell segmentation blobs)")
    parser.add_argument("--output", type=Path, required=True,
                        help="Per-slide, per-channel metrics CSV")
    parser.add_argument("--fov-output", type=Path,
                        help="Per-FOV, per-channel metrics CSV (optional)")
    parser.add_argument("--figure", type=Path, help="Summary figure PNG (optional)")
    args = parser.parse_args()

    if not args.flatfiles_dir and not args.metadata:
        parser.error("give --flatfiles-dir or at least one --metadata")

    targets = metadata_paths(args.flatfiles_dir, args.metadata)
    if not targets:
        print("ERROR: no metadata files found", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(targets)} metadata file(s)")

    slide_rows: list[dict] = []
    fov_rows: list[dict] = []
    for slide_id, path in targets:
        print(f"Reading {path.name}")
        slide, fov = analyse_slide(slide_id, path, args)
        slide_rows.extend(slide)
        fov_rows.extend(fov)

    if not slide_rows:
        print("ERROR: no slide produced metrics", file=sys.stderr)
        sys.exit(1)

    slide_frame = pd.DataFrame(slide_rows)
    channels = list(slide_frame["channel"].unique())
    nuclear = resolve_role(channels, NUCLEAR_ALIASES, args.nuclear_channel, "nuclear")
    reference = resolve_role(channels, REFERENCE_ALIASES, args.reference_channel,
                             "reference")
    if not nuclear or not reference:
        raise MissingChannelError(
            f"could not identify nuclear ({nuclear}) and reference ({reference}) "
            f"channels among {channels}; pin them with --nuclear-channel / "
            f"--reference-channel")

    if args.manifest and args.manifest.exists():
        manifest = pd.read_csv(args.manifest, dtype={"slide_id": str})
        columns = [c for c in ("slide_id", "run_date", "instrument_id", "export_batch")
                   if c in manifest.columns]
        slide_frame = slide_frame.merge(manifest[columns], on="slide_id", how="left")

    slide_frame.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

    if fov_rows:
        pd.DataFrame(fov_rows).to_csv(args.fov_output, index=False)
        print(f"Wrote {args.fov_output}")

    ratios = contrast_ratio_table(slide_frame, nuclear, reference)

    print(f"\nPer-channel medians across {slide_frame['slide_id'].nunique()} slide(s):")
    summary = slide_frame.groupby("channel")[
        ["contrast_index", "peakedness", "peakedness_tail", "area_rho",
         "rrna_rho_given_area"]].median()
    print(summary.round(3).to_string())

    print(f"\n{reference} over {nuclear} contrast ratio: "
          f"median {ratios['reference_over_nuclear'].median():.2f}, "
          f"range {ratios['reference_over_nuclear'].min():.2f}-"
          f"{ratios['reference_over_nuclear'].max():.2f}")
    print(f"Slides where {nuclear} matches or beats {reference}: "
          f"{int((ratios['reference_over_nuclear'] <= 1).sum())} of {len(ratios)}")

    if args.figure:
        write_figure(slide_frame, ratios, nuclear, args.figure)


if __name__ == "__main__":
    main()
