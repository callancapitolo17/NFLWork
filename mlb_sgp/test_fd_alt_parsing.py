"""Unit tests for FD alt-market name parsing.

These guard the regex + unicode-minus normalizer in scraper_fanduel_singles.
They are intentionally tiny — the function under test is just regex matching
plus a string replace — but they catch the silent-fail mode where FD swaps
ASCII '-' for U+2212 (true minus) or an em/en-dash, and the alt-spread
regex (`[+-]`) refuses to match. Without these tests, a single character
swap on FD's side could drop every F5/F3 alt-spread row.
"""
import sys
from pathlib import Path

# Allow `from scraper_fanduel_singles import ...` when pytest runs from the
# repo root (the package is conventionally imported via mlb_sgp/, but this
# file lives inside mlb_sgp so we splice in the dir.
sys.path.insert(0, str(Path(__file__).parent))

from scraper_fanduel_singles import (
    _FD_ALT_SPREAD_RE,
    _FD_ALT_TOTAL_RE,
    _normalize_alt_name,
)


def test_ascii_alt_spread():
    m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox +1.5"))
    assert m is not None
    assert m.group("line") == "+1.5"


def test_unicode_minus_alt_spread():
    m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox −1.5"))
    assert m is not None
    assert m.group("line") == "-1.5"


def test_em_dash_alt_spread():
    m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox —1.5"))
    assert m is not None
    assert m.group("line") == "-1.5"


def test_en_dash_alt_spread():
    m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox –1.5"))
    assert m is not None
    assert m.group("line") == "-1.5"


def test_extra_whitespace_alt_total():
    m = _FD_ALT_TOTAL_RE.match(_normalize_alt_name("Over   (8.5)"))
    assert m is not None
    assert m.group("line") == "8.5"


def test_no_parens_alt_total_over():
    # FD's full-game "Alternate Total Runs" runners are named WITHOUT parens
    # ("Over 7.5"), unlike the F5 ladder which uses "Over (7.5)". Both must
    # parse or FG alt-totals silently drop (every runner fails -> 0 rows).
    m = _FD_ALT_TOTAL_RE.match(_normalize_alt_name("Over 7.5"))
    assert m is not None
    assert m.group("line") == "7.5"


def test_no_parens_alt_total_under():
    m = _FD_ALT_TOTAL_RE.match(_normalize_alt_name("Under 10.5"))
    assert m is not None
    assert m.group("line") == "10.5"


def test_parens_alt_total_still_parses():
    # Regression guard: the F5 paren format must keep working after the fix.
    m = _FD_ALT_TOTAL_RE.match(_normalize_alt_name("Over (2.5)"))
    assert m is not None
    assert m.group("line") == "2.5"


def test_garbage_returns_none():
    m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name("garbage input"))
    assert m is None
