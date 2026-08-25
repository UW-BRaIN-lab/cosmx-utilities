#!/usr/bin/env python3
"""Measure the physical gap between two FOV runs on one slide, from cell centroids.

Answers: are these two FOV runs one continuous tissue mass, or two separate pieces?
That decides whether they may share one InSituCNV spatial neighbour graph -- the graph
is kNN (scv.pp.neighbors, n_neighbors=20), so it bridges any gap smaller than the local
20-NN radius no matter how far apart the FOV grids look.

Stdlib only, and Python 3.6+ -- runs on a Hyak login node, which has neither
pandas/numpy nor a recent interpreter.

    python3 tissue_section_gap.py <workdir> <slide_id> <run_a> <run_b>
    e.g.  python3 tissue_section_gap.py /gscratch/scrubbed/$USER 7353A6A18A17 1-135 181-205
"""
import csv
import gzip
import math
import sys
from collections import defaultdict


def distance(first, second):
    """Euclidean distance; math.hypot rather than math.dist, which needs 3.8+."""
    return math.hypot(first[0] - second[0], first[1] - second[1])

GRID_PX = 1000          # spatial hash cell size for the nearest-neighbour search
MAX_SEARCH_PX = 20000   # give up beyond this; anything this far apart is clearly separate
DENSITY_SAMPLE = 20     # k for the local-radius estimate, matching --n-neighbors


def parse_run(text):
    low, high = text.split("-")
    return range(int(low), int(high) + 1)


def load_scale(positions_path):
    """px per mm, derived from the two coordinate systems in fov_positions."""
    with gzip.open(positions_path, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for earlier in rows:
        for later in rows:
            dx_px = float(later["x_global_px"]) - float(earlier["x_global_px"])
            dx_mm = float(later["x_global_mm"]) - float(earlier["x_global_mm"])
            if abs(dx_mm) > 0.5:
                return abs(dx_px / dx_mm)
    return None


def load_centroids(metadata_path, run_a, run_b):
    want = {}
    for fov in run_a:
        want[fov] = "a"
    for fov in run_b:
        want[fov] = "b"
    points = {"a": [], "b": []}
    with gzip.open(metadata_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            side = want.get(int(row["fov"]))
            if side is None:
                continue
            points[side].append((float(row["CenterX_global_px"]),
                                 float(row["CenterY_global_px"])))
    return points["a"], points["b"]


def build_grid(points):
    grid = defaultdict(list)
    for point in points:
        grid[(int(point[0]) // GRID_PX, int(point[1]) // GRID_PX)].append(point)
    return grid


def nearest_distance(point, grid):
    """Distance from point to the closest grid member, or None past MAX_SEARCH_PX."""
    cx, cy = int(point[0]) // GRID_PX, int(point[1]) // GRID_PX
    best = None
    for ring in range(0, MAX_SEARCH_PX // GRID_PX + 1):
        if best is not None and best <= ring * GRID_PX:
            break                      # no closer point can hide in a further ring
        for gx in range(cx - ring, cx + ring + 1):
            for gy in range(cy - ring, cy + ring + 1):
                if ring and max(abs(gx - cx), abs(gy - cy)) != ring:
                    continue           # interior already scanned
                for other in grid.get((gx, gy), ()):
                    d = distance(point, other)
                    if best is None or d < best:
                        best = d
    return best


def local_knn_radius(points, grid, sample_every):
    """Median distance to the DENSITY_SAMPLE-th neighbour within the same tissue."""
    radii = []
    for i in range(0, len(points), sample_every):
        point = points[i]
        cx, cy = int(point[0]) // GRID_PX, int(point[1]) // GRID_PX
        near = [other for gx in range(cx - 1, cx + 2)
                for gy in range(cy - 1, cy + 2)
                for other in grid.get((gx, gy), ())]
        if len(near) <= DENSITY_SAMPLE:
            continue
        dists = sorted(distance(point, other) for other in near)
        radii.append(dists[DENSITY_SAMPLE])   # skip self at index 0
    radii.sort()
    return radii[len(radii) // 2] if radii else None


def main():
    workdir, slide_id, run_a_text, run_b_text = sys.argv[1:5]
    run_a, run_b = parse_run(run_a_text), parse_run(run_b_text)

    try:
        px_per_mm = load_scale(f"{workdir}/{slide_id}_fov_positions_file.csv.gz")
    except FileNotFoundError:
        px_per_mm = None   # distances still comparable, just reported in px
    def as_um(px):
        return f"{px / px_per_mm * 1000:8.1f} um" if px_per_mm else f"{px:8.1f} px"

    points_a, points_b = load_centroids(
        f"{workdir}/{slide_id}_metadata_file.csv.gz", run_a, run_b)
    print(f"=== {slide_id}: FOVs {run_a_text} vs {run_b_text} ===")
    print(f"  scale: {px_per_mm:.1f} px/mm" if px_per_mm else "  scale: unknown")
    print(f"  run A ({run_a_text}): {len(points_a):,} cells")
    print(f"  run B ({run_b_text}): {len(points_b):,} cells")
    if not points_a or not points_b:
        sys.exit("one run has no cells -- nothing to compare")

    grid_a = build_grid(points_a)
    grid_b = build_grid(points_b)

    print(f"\n  Local {DENSITY_SAMPLE}-NN radius WITHIN each run "
          f"(the distance the CNV graph actually reaches):")
    for name, points, grid in (("A", points_a, grid_a), ("B", points_b, grid_b)):
        radius = local_knn_radius(points, grid, max(1, len(points) // 2000))
        print(f"    run {name}: {as_um(radius)}" if radius else f"    run {name}: too sparse")

    print(f"\n  Cross-run nearest-neighbour distances (B -> A), "
          f"sampled over {min(len(points_b), 4000):,} cells:")
    step = max(1, len(points_b) // 4000)
    dists = sorted(d for d in (nearest_distance(p, grid_a) for p in points_b[::step])
                   if d is not None)
    if not dists:
        print(f"    no run-A cell within {MAX_SEARCH_PX:,} px of any run-B cell "
              f"-- unambiguously separate tissue")
        return
    for label, value in (("minimum", dists[0]),
                         ("1st pct", dists[len(dists) // 100]),
                         ("median ", dists[len(dists) // 2])):
        print(f"    {label}: {as_um(value)}")
    print(f"\n  VERDICT: compare the minimum above against the within-run "
          f"{DENSITY_SAMPLE}-NN radius.\n"
          f"  minimum >> radius  -> separate pieces, the kNN graph cannot bridge them\n"
          f"  minimum <= radius  -> continuous tissue, one section is correct")


if __name__ == "__main__":
    main()
