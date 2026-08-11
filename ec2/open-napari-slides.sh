#!/usr/bin/env bash
# Open one Napari window per stitched slide, for side-by-side spatial comparison.
#
# Run this inside the DCV desktop session on a napari EC2 instance. With no
# arguments it opens every slide under /mnt/local/stitched; otherwise it opens
# the slide directories you pass.
#
#   ./ec2/open-napari-slides.sh                       # all synced slides
#   ./ec2/open-napari-slides.sh /mnt/local/stitched/slide2 /mnt/local/stitched/eyes7517
#
# Why not `napari --plugin napari-cosmx-fork <slide>`? The napari-cosmx reader
# adds its layers as a side effect and returns no layer data, so the headless
# `--plugin` path opens nothing and exits. We launch the viewer programmatically
# (Gemini(...) + napari.run()) instead, which also skips the reader-choice dialog.
set -uo pipefail

REPO_DIR="${REPO_DIR:-/opt/cosmx-utilities}"
STITCHED_DIR="${STITCHED_DIR:-/mnt/local/stitched}"

slides=("$@")
if [ ${#slides[@]} -eq 0 ]; then
    slides=("$STITCHED_DIR"/*/)
fi

for slide in "${slides[@]}"; do
    slide="${slide%/}"
    [ -d "$slide/images" ] || { echo "skip (no images/): $slide"; continue; }
    name="$(basename "$slide")"
    echo "opening $name"
    (
        cd "$REPO_DIR" && uv run python - "$slide" <<'PY'
import sys
import napari
from napari_cosmx.gemini import Gemini
from napari_cosmx._dock_widget import GeminiQWidget

# Mirror napari_cosmx's reader_function: create the viewer, load the slide, and
# add the Gemini control panel (Color Cells + channels) on the right. Launching
# Gemini() alone loads the data but omits the dock widget.
viewer = napari.Viewer()
gem = Gemini(sys.argv[1], viewer=viewer)
viewer.window.add_dock_widget(GeminiQWidget(viewer, gem), area="right", name=gem.name)
napari.run()
PY
    ) >"/tmp/napari-${name}.log" 2>&1 &
done

echo "launched ${#slides[@]} viewer(s); each opens its own window."
