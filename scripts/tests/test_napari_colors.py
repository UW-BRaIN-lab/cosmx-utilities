"""Tests for napari_colors.py: the one hex-colour rule for Napari metadata.

Run:  uv run --with duckdb --with pytest python -m pytest scripts/tests/test_napari_colors.py -s
"""
import colorsys
import importlib.util
import random
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("napari_colors", SCRIPTS / "napari_colors.py")
nc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nc)

PREVIOUSLY_INVISIBLE = "#0A0404"      # Excitatory_L4_5_IT_like, as Maddie saw it
DARK_BLUE = "#190B93"                 # Inhibitory neuron C
ALREADY_VISIBLE = "#F5C242"
RANDOM_COLOR_COUNT = 5000


def _hue(hex_color):
    r, g, b = (c / 255 for c in nc._channels(hex_color))
    return colorsys.rgb_to_hls(r, g, b)[0]


def test_contrast_endpoints():
    print(nc.contrast_against_black("#000000"), nc.contrast_against_black("#FFFFFF"))
    assert nc.contrast_against_black("#000000") == 1.0
    assert abs(nc.contrast_against_black("#FFFFFF") - 21.0) < 1e-9


def test_a_visible_colour_is_returned_unchanged():
    assert nc.legible_color(ALREADY_VISIBLE) == ALREADY_VISIBLE
    assert nc.legible_color(ALREADY_VISIBLE.lower()) == ALREADY_VISIBLE


def test_a_near_black_colour_is_lightened_to_the_minimum_contrast():
    fixed = nc.legible_color(PREVIOUSLY_INVISIBLE)
    print(PREVIOUSLY_INVISIBLE, "->", fixed, round(nc.contrast_against_black(fixed), 2))
    assert fixed != PREVIOUSLY_INVISIBLE
    assert nc.contrast_against_black(fixed) >= nc.MIN_CONTRAST_AGAINST_BLACK


def test_hue_is_kept_for_a_dark_saturated_colour():
    fixed = nc.legible_color(DARK_BLUE)
    print(DARK_BLUE, "->", fixed)
    assert abs(_hue(fixed) - _hue(DARK_BLUE)) < 0.02
    assert nc.contrast_against_black(fixed) >= nc.MIN_CONTRAST_AGAINST_BLACK


def test_legible_color_is_idempotent_and_never_darkens():
    rng = random.Random(0)
    changed = 0
    for _ in range(RANDOM_COLOR_COUNT):
        color = "#%06X" % rng.randrange(nc.COLOR_SPACE)
        once = nc.legible_color(color)
        assert nc.legible_color(once) == once
        assert nc.contrast_against_black(once) >= nc.MIN_CONTRAST_AGAINST_BLACK - 1e-9
        assert nc.relative_luminance(once) >= nc.relative_luminance(color) - 1e-9
        changed += once != color
    share = changed / RANDOM_COLOR_COUNT
    print(f"share of random colours changed: {share:.3f}")
    assert 0.10 < share < 0.25     # only the dark tail moves


def test_deterministic_color_is_the_legible_hash_and_stable():
    label = "Excitatory_L4_5_IT_like"
    raw = nc.hash_color(label)
    print(label, raw, "->", nc.deterministic_color(label))
    assert raw == PREVIOUSLY_INVISIBLE
    assert nc.deterministic_color(label) == nc.legible_color(raw) != raw
    assert nc.deterministic_color(label) == nc.deterministic_color(label)
    assert nc.deterministic_color("Microglia A") != nc.deterministic_color("Microglia B")
