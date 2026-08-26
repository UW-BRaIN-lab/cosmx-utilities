#!/usr/bin/env bash
# Open one Napari window per stitched slide, for side-by-side spatial comparison.
#
# Run this inside the DCV desktop session on a napari EC2 instance.
#
#   ./ec2/open-napari-slides.sh                         # discover & open (small sets)
#   ./ec2/open-napari-slides.sh <slide-dir> [<slide-dir> ...]   # open specific slides
#
# Slide discovery finds every directory containing an `images/` subdir under
# STITCHED_DIR, at any depth — so it works for both a flat layout
# (/mnt/local/stitched/<slide>/) and the nested experiment/run layout
# (/mnt/local/stitched/<run>/<slide>/).
#
# Safety: with no arguments, if more than MAX_AUTO slides are found it lists them
# and stops instead of opening all at once (many viewers at once will swamp a
# software-rendered desktop). Pass the specific slide dirs you want to compare.
#
# Why not `napari --plugin napari-cosmx-fork <slide>`? The napari-cosmx reader
# adds its layers as a side effect and returns no layer data, so the headless
# `--plugin` path opens nothing and exits. We launch the viewer programmatically
# (Gemini + GeminiQWidget + napari.run()), which also skips the reader dialog.
set -uo pipefail

REPO_DIR="${REPO_DIR:-/opt/cosmx-utilities}"
STITCHED_DIR="${STITCHED_DIR:-/mnt/local/stitched}"
MAX_AUTO="${MAX_AUTO:-6}"   # max slides to auto-open when no args are given

open_slide() {
    local slide="${1%/}"
    local name
    name="$(basename "$slide")"
    echo "opening $name"
    (
        cd "$REPO_DIR" && uv run python - "$slide" <<'PY'
import sys
import napari
from napari.settings import get_settings
from napari_cosmx.gemini import Gemini
from napari_cosmx._dock_widget import GeminiQWidget

# napari 0.6 defaults this off, which drops the per-cell readout on hover --
# cell_ID, cell type and the rest live in the Metadata layer's features and are
# surfaced through this tooltip.
get_settings().appearance.layer_tooltip_visibility = True

# Mirror napari_cosmx's reader_function: create the viewer, load the slide, and
# add the Gemini control panel (Color Cells + channels) on the right.
viewer = napari.Viewer()
gem = Gemini(sys.argv[1], viewer=viewer)
# show_stitching_widget=False as the reader does. That panel builds mosaics
# rather than viewing them, so it is dead weight here -- and its height pushes
# the window's minimum past the desktop's, leaving a window taller than the
# screen that xfwm will only resize sideways, with the status bar and console
# unreachable below the bottom edge.
viewer.window.add_dock_widget(
    GeminiQWidget(viewer, gem, show_stitching_widget=False),
    area="right", name=gem.name)
# napari opens at 640x480; size to the desktop so the whole UI is reachable.
window = viewer.window._qt_window
available = window.screen().availableGeometry()
window.move(available.left(), available.top())
window.resize(available.width(), available.height())
napari.run()
PY
    ) >"/tmp/napari-${name}.log" 2>&1 &
}

if [ "$#" -gt 0 ]; then
    # Explicit slide dirs: open exactly what was asked.
    n=0
    for slide in "$@"; do
        slide="${slide%/}"
        [ -d "$slide/images" ] || { echo "skip (no images/): $slide"; continue; }
        open_slide "$slide"; n=$((n + 1))
    done
    echo "launched $n viewer(s)."
    exit 0
fi

# No args: discover slide dirs (any depth) by locating their images/ subdirs.
found=()
while IFS= read -r d; do found+=("$d"); done < <(
    find "$STITCHED_DIR" -maxdepth 3 -type d -name images -printf '%h\n' 2>/dev/null | sort
)

if [ "${#found[@]}" -eq 0 ]; then
    echo "No slides found under $STITCHED_DIR (looked for */images and */*/images)."
    exit 1
fi

if [ "${#found[@]}" -gt "$MAX_AUTO" ]; then
    echo "Found ${#found[@]} slides under $STITCHED_DIR — too many to open at once"
    echo "(software rendering; opening all would swamp the desktop)."
    echo "Re-run with the specific slides you want to compare, e.g.:"
    echo "  $0 ${found[0]} ${found[1]}"
    echo ""
    echo "Available slides:"
    printf '  %s\n' "${found[@]}"
    exit 2
fi

for slide in "${found[@]}"; do open_slide "$slide"; done
echo "launched ${#found[@]} viewer(s); each opens its own window."
