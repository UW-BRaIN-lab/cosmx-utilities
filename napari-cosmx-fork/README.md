This is a fork of Napari CosMX to speed development of our napari-cosmx-docker project.

## Upstream baseline

Currently tracking **napari-CosMx 0.5.0.0**, taken from the wheel published at
`assets/napari-cosmx releases/napari_CosMx-0.5.0.0-py3-none-any.whl` in
[Nanostring-Biostats/CosMx-Analysis-Scratch-Space](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space)
(announced in the [3D + multiomics post](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/napari-3d-multiomics/)).
Upstream ships wheels rather than source, so a version bump means unpacking the
new wheel over `src/napari_cosmx/` and re-applying the patches below.

What 0.5.0.0 brought us: 3D cell segmentation with Z-plane navigation, multi-omics
layers, a layer group manager, and support for napari up to 0.6.6.

**Reading older stitched output is unaffected.** `gemini.py` defaults `ndim` to 2
and `z_slices` to `[0]` when a zarr has no such attributes, so mosaics written by
0.4.17.4 — which is most of our studies, acquired before the instrument could do
3D segmentation — keep loading as plain 2D. No re-stitching required.

## Local patches

Everything this fork changes relative to stock 0.5.0.0. Keep this list current;
it is the checklist for the next upstream bump.

| Area | Change | Why |
| --- | --- | --- |
| `__init__.py` | Guard the `_reader` / `_function` imports | The headless Fargate image has no Qt, and the plugin hooks would otherwise fail at import |
| `napari.yaml` | Plugin renamed `napari-cosmx-fork` | Lets the fork coexist with an installed stock plugin instead of colliding on the manifest name |
| `gemini.py`, `_reader.py` | Window title is `<study>/<slide>` | The leaf directory alone is not enough to tell slides apart when several viewers are open |
| `_colors.py` (new), `gemini.py`, `_dock_widget.py` | Prefer a `<column>_color` column, then legacy `hex_color` | Keeps a category the same color across every slide; also hides `*_color` helper columns from the color-by menus and fixes NaN category lookup |
| `utils/stitch_images.py` | `--celllabels-subdir` | AtoMx writes one `Segmentation_<uuid>_00N` directory per segmentation run, and a run's analysis is only valid against its own. Upstream's collector prunes every non-`FOV*` directory, so these are otherwise unreachable |
| `utils/stitch_images.py` | `--channel-name CH=MARKER` | MorphologyKit metadata is sometimes wrong (6E10 labeled CD68/CD3, AT8 labeled Histone) and inconsistent across studies for the same panel. Also restores our named channel defaults, which 0.5.0.0 replaced with raw `B/G/Y/R/U` |
| `utils/stitch_images.py` | `--morphology-ndim` plus auto-detection | A 3D resegmentation is exported against the *original 2D* morphology acquisition — `Morphology2D/`, no `_Z###` in filenames. Upstream's single `--input-ndim` drives labels and morphology together, so 3D labels would yield zero morphology images. The single 2D plane is broadcast across all output z planes |
| `utils/stitch_images.py` | Version lookup tries the fork's distribution name | This fork installs as `napari-cosmx-fork`, so upstream's bare `version('napari_cosmx')` raises `PackageNotFoundError` |

## Tests

`src/napari_cosmx/_tests/test_fork_stitch.py` and `test_colors.py` cover the
additions above; upstream's own `test_stitch.py` is kept as shipped.

```bash
uv run python napari-cosmx-fork/src/napari_cosmx/_tests/test_fork_stitch.py
```

Tests that build a real napari viewer need a display; they segfault under
offscreen Qt on macOS (stock 0.5.0.0 does too, so it is not something this fork
introduced). Run those on a machine with a display server.
