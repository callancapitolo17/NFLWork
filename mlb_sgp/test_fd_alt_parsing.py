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


# --- Aggregate parse-failure tripwire ("N seen, 0 parsed") -----------------
# The FG alt-total bug shipped because every runner individually failed the
# regex and the scraper still exited 0 with a normal-looking summary. The
# tripwire counts, per alt market type, how many runners needed a name-parse
# vs how many succeeded — "seen > 0, ok == 0" is the signature of FD format
# drift and must be loudly distinguishable from "FD posted no alts today".

from datetime import datetime, timezone

from fd_client import Event, Runner
from scraper_fanduel_singles import new_alt_parse_stats, parse_runners_to_wide_rows

_EVENT = Event(
    event_id="e1",
    home_team="New York Yankees",
    away_team="Boston Red Sox",
    start_time="2026-06-10T23:05:00.000Z",
)
_FETCH = datetime.now(timezone.utc)


def _alt_total_runner(name: str, rid: str) -> Runner:
    return Runner(runner_id=rid, market_id="m1", name=name,
                  line=None, american_odds=-110)


def test_stats_all_alt_totals_fail_to_parse():
    # Simulated format drift: FD renames runners so nothing matches the regex.
    stats = new_alt_parse_stats()
    runners = [_alt_total_runner("Total Runs Above 7.5", "r1"),
               _alt_total_runner("Total Runs Below 7.5", "r2")]
    rows = parse_runners_to_wide_rows(
        _EVENT, runners, {"m1": ("FG", "alternate_totals")}, _FETCH,
        alt_parse_stats=stats)
    assert rows == []
    assert stats["alternate_totals"]["seen"] == 2
    assert stats["alternate_totals"]["ok"] == 0


def test_stats_healthy_alt_totals_parse():
    stats = new_alt_parse_stats()
    runners = [_alt_total_runner("Over 7.5", "r1"),
               _alt_total_runner("Under 7.5", "r2")]
    rows = parse_runners_to_wide_rows(
        _EVENT, runners, {"m1": ("FG", "alternate_totals")}, _FETCH,
        alt_parse_stats=stats)
    assert len(rows) == 1  # Over + Under collapse into one wide row
    assert stats["alternate_totals"]["seen"] == 2
    assert stats["alternate_totals"]["ok"] == 2


def test_stats_param_is_optional():
    # Existing callers that don't pass stats must keep working unchanged.
    runners = [_alt_total_runner("Over 7.5", "r1"),
               _alt_total_runner("Under 7.5", "r2")]
    rows = parse_runners_to_wide_rows(
        _EVENT, runners, {"m1": ("FG", "alternate_totals")}, _FETCH)
    assert len(rows) == 1
