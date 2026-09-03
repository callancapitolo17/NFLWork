"""The shared start-time disambiguator behind every game resolution.

A Kalshi MLB event's team pair is a filter, not an identity — see
``leg_types.unique_game_by_start``. This file pins the four outcomes the
callers branch on, because each one has a different consequence: a unique
match prices, ``unmatched`` and ``ambiguous`` decline, and ``no_kalshi_start``
means the ticker itself was unreadable.
"""
from datetime import datetime, timedelta, timezone

from kalshi_common.leg_types import (SCHEDULE_START_TOLERANCE_MIN,
                                     as_naive_utc, parse_suffix_start_utc,
                                     unique_game_by_start)

# 2026-09-04, 14:10 and 19:15 ET -> 18:10 and 23:15 UTC.
DH1 = datetime(2026, 9, 4, 18, 10)
DH2 = datetime(2026, 9, 4, 23, 15)


def test_a_lone_candidate_matches():
    assert unique_game_by_start(DH1, [("g", DH1)]) == ("g", "")


def test_each_doubleheader_game_picks_its_own_row():
    rows = [("game-1", DH1), ("game-2", DH2)]
    assert unique_game_by_start(DH1, rows)[0] == "game-1"
    assert unique_game_by_start(DH2, rows)[0] == "game-2"


def test_two_rows_inside_the_tolerance_fail_closed():
    near = DH1 + timedelta(minutes=SCHEDULE_START_TOLERANCE_MIN - 1)
    assert unique_game_by_start(DH1, [("a", DH1), ("b", near)]) == (None, "ambiguous")


def test_nothing_in_range_is_unmatched_not_nearest_neighbour():
    far = DH1 + timedelta(hours=6)
    assert unique_game_by_start(DH1, [("a", far)]) == (None, "unmatched")


def test_an_unreadable_kalshi_start_declines():
    assert unique_game_by_start(None, [("a", DH1)]) == (None, "no_kalshi_start")


def test_empty_candidates_are_unmatched():
    assert unique_game_by_start(DH1, []) == (None, "unmatched")


def test_aware_and_naive_starts_compare_without_raising():
    """The Odds API sends tz-aware; mlb_target_lines stores naive UTC.

    Comparing them raw raises TypeError, which the callers' fail-safe
    ``except`` would swallow into a silent decline on EVERY game.
    """
    aware = DH1.replace(tzinfo=timezone.utc)
    assert unique_game_by_start(DH1, [("g", aware)]) == ("g", "")
    assert unique_game_by_start(aware, [("g", DH1)]) == ("g", "")


def test_a_row_with_no_start_time_is_skipped_not_matched():
    assert unique_game_by_start(DH1, [("a", None)]) == (None, "unmatched")


def test_the_real_skew_between_the_two_feeds_is_inside_the_tolerance():
    """Live 2026-09-01: Kalshi's suffix minute and the Odds API differ by 1.

    Suffix 26SEP012210STLLAD is 22:10 ET; the Odds API said 2026-09-02T02:11Z.
    An exact-equality match would decline every game in the league.
    """
    kalshi = parse_suffix_start_utc("26SEP012210STLLAD")
    odds_api = datetime(2026, 9, 2, 2, 11, tzinfo=timezone.utc)
    assert abs((kalshi - as_naive_utc(odds_api)).total_seconds()) == 60
    assert unique_game_by_start(kalshi, [("g", odds_api)]) == ("g", "")
