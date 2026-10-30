"""The one hex-colour rule for every Napari `_metadata.csv` column.

A label's colour is a hash of its name, so a label keeps its colour whichever script wrote the
column. The hash alone picks near-black colours for some labels (`Excitatory_L4_5_IT_like` was
#0A0404), which cannot be seen on Napari's black background. `legible_color` lightens such a colour
and keeps its hue, so only the colours that were too dark change; every other label keeps the colour
the hash gave it.

Import this instead of redefining `deterministic_color`: the writers and the validators must agree.
"""

from __future__ import annotations

import colorsys

import duckdb

HEX_DIGITS_PER_CHANNEL = 2
MAX_CHANNEL = 255
COLOR_SPACE = 16777216  # 2**24 hex colours
# WCAG contrast ratio of a colour against black: (luminance + 0.05) / 0.05. 3:1 is the WCAG minimum for
# graphical objects. Colours below it (about 18% of the hash space) are lightened up to it.
MIN_CONTRAST_AGAINST_BLACK = 3.0
BLACK_LUMINANCE_OFFSET = 0.05
LIGHTNESS_SEARCH_STEPS = 24


def _channels(hex_color: str) -> tuple[int, int, int]:
    digits = hex_color.lstrip("#")
    step = HEX_DIGITS_PER_CHANNEL
    return tuple(int(digits[i:i + step], 16) for i in range(0, 3 * step, step))


def _to_hex(red: float, green: float, blue: float) -> str:
    r, g, b = (round(c * MAX_CHANNEL) for c in (red, green, blue))
    return f"#{r:02X}{g:02X}{b:02X}"


def _linear(channel: float) -> float:
    return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4


def relative_luminance(hex_color: str) -> float:
    red, green, blue = (_linear(c / MAX_CHANNEL) for c in _channels(hex_color))
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def contrast_against_black(hex_color: str) -> float:
    return (relative_luminance(hex_color) + BLACK_LUMINANCE_OFFSET) / BLACK_LUMINANCE_OFFSET


def legible_color(hex_color: str, min_contrast: float = MIN_CONTRAST_AGAINST_BLACK) -> str:
    """`hex_color` itself if it is visible on black, else the same hue lightened until it is.

    Idempotent, so it can be applied to a file that was already corrected.
    """
    if contrast_against_black(hex_color) >= min_contrast:
        return hex_color.upper()
    red, green, blue = (c / MAX_CHANNEL for c in _channels(hex_color))
    hue, lightness, saturation = colorsys.rgb_to_hls(red, green, blue)
    low, high = lightness, 1.0
    for _ in range(LIGHTNESS_SEARCH_STEPS):
        mid = (low + high) / 2
        if contrast_against_black(_to_hex(*colorsys.hls_to_rgb(hue, mid, saturation))) >= min_contrast:
            high = mid
        else:
            low = mid
    return _to_hex(*colorsys.hls_to_rgb(hue, high, saturation))


def hash_color(value: str) -> str:
    """The raw hash colour, before the legibility rule (DuckDB's hash, as every writer used)."""
    return duckdb.sql(f"SELECT printf('#%06X', abs(hash($1)) % {COLOR_SPACE})",
                      params=[value]).fetchone()[0]


def deterministic_color(value: str) -> str:
    """The colour a label gets in every Napari metadata column: its hash colour, made legible on black."""
    return legible_color(hash_color(value))
