"""Issue #103 Phase 3 — deliberate start-time horizon on target-line
enumeration.

`warm_cycle` fans out one structure fetch per game per book per pass over
`SELECT DISTINCT game_id FROM mlb_target_lines` with no time filter. Today
that is ~5-7 games only because the Odds API `/events` window is short —
an accidental rate limiter. Phase 4 swaps in the full Kalshi board (35-48
open games, ~3 days out), so the window must become deliberate FIRST.

One filter, one place: `enumerate_kalshi_targets` drops any matched game
whose commence_time is after `now + horizon_hours`. `mlb_target_lines` is
a full DELETE+INSERT of that list, so every reader inherits the bound.
"""
from __future__ import annotations

import logging
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common import sgp_runner

NOW = datetime(2026, 9, 10, 15, 0, tzinfo=timezone.utc)
HORIZON_HOURS = 24.0
# Kalshi suffix grammar: YYMMMDDHHMM{Away}{Home}, HH:MM in US/Eastern.
# 11:00 ET on Sep 10 = 15:00 UTC = NOW.
SUFFIX_BY_GAME = {
    "today": "26SEP101100BOSNYY",       # 15:00Z, 0h ahead
    "tomorrow": "26SEP111100BOSNYY",    # 15:00Z +1d = exactly NOW + 24h
    "day_after": "26SEP121100BOSNYY",   # 15:00Z +2d
}
TEAMS = ("New York Yankees", "Boston Red Sox")   # (home, away)


def _schedule_row(game_id: str, commence_time: datetime) -> dict:
    return {"game_id": game_id, "home_team": TEAMS[0], "away_team": TEAMS[1],
            "commence_time": commence_time}


class _FrozenDatetime(datetime):
    """`datetime` stand-in whose `now()` is pinned to NOW — the module reads
    the wall clock itself (no injectable clock, by user decision)."""
    @classmethod
    def now(cls, tz=None):
        return NOW if tz is not None else NOW.replace(tzinfo=None)


def _install_fakes(monkeypatch, games: dict[str, datetime]):
    """Pin the clock to NOW; Kalshi lists every game in `games`; the Odds
    API schedules each at its commence_time; each game carries one spread x
    one total."""
    monkeypatch.setattr(sgp_runner, "datetime", _FrozenDatetime)
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events", lambda: [
        {"event_ticker": f"KXMLBGAME-{SUFFIX_BY_GAME[g]}"} for g in games])
    monkeypatch.setattr(sgp_runner, "_fetch_schedule_from_odds_api",
                        lambda: [_schedule_row(g, ct) for g, ct in games.items()])
    line_fetches: list[str] = []

    def fake_spreads(suffix, **kw):
        line_fetches.append(suffix)
        return [(-1.5, "home")]

    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_spread_lines", fake_spreads)
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_total_lines", lambda suffix: [8.5])
    return line_fetches


THREE_DAYS = {
    "today": NOW,
    "tomorrow": NOW + timedelta(hours=24),
    "day_after": NOW + timedelta(hours=48),
}


def test_inside_horizon_kept_outside_dropped(monkeypatch):
    _install_fakes(monkeypatch, THREE_DAYS)
    targets = sgp_runner.enumerate_kalshi_targets(
        horizon_hours=HORIZON_HOURS)
    assert sorted({t.game_id for t in targets}) == ["today", "tomorrow"]


def test_boundary_game_exactly_at_horizon_is_kept(monkeypatch):
    """`<=` at the deadline: a game starting exactly now + horizon is in."""
    _install_fakes(monkeypatch, {"tomorrow": NOW + timedelta(hours=24)})
    targets = sgp_runner.enumerate_kalshi_targets(
        horizon_hours=HORIZON_HOURS)
    assert [t.game_id for t in targets] == ["tomorrow"]


def test_one_second_past_horizon_is_dropped(monkeypatch):
    # The Kalshi suffix is minute-precise; only the Odds API side carries
    # seconds. Both feeds still match (30-min tolerance), so the horizon is
    # the only thing that can drop this game.
    _install_fakes(monkeypatch,
                   {"tomorrow": NOW + timedelta(hours=24, seconds=1)})
    targets = sgp_runner.enumerate_kalshi_targets(
        horizon_hours=HORIZON_HOURS)
    assert targets == []


def test_dropped_game_costs_zero_kalshi_line_fetches(monkeypatch):
    """The horizon check runs BEFORE the per-game spread/total fetches, so
    a far-out game is zero requests, not just zero rows."""
    line_fetches = _install_fakes(monkeypatch, THREE_DAYS)
    sgp_runner.enumerate_kalshi_targets(horizon_hours=HORIZON_HOURS)
    assert line_fetches == [SUFFIX_BY_GAME["today"], SUFFIX_BY_GAME["tomorrow"]]


def test_naive_commence_time_is_compared_as_utc(monkeypatch):
    """A naive Odds API commence_time (tests seed them; `write_target_lines`
    stores naive UTC) must compare as UTC, not raise against the aware clock."""
    _install_fakes(monkeypatch, {
        "today": NOW.replace(tzinfo=None),
        "day_after": (NOW + timedelta(hours=48)).replace(tzinfo=None),
    })
    targets = sgp_runner.enumerate_kalshi_targets(horizon_hours=HORIZON_HOURS)
    assert [t.game_id for t in targets] == ["today"]


def test_horizon_drop_count_logged_at_info(monkeypatch, caplog):
    """The knob must be visible every cycle in bot.log (INFO, not print)."""
    _install_fakes(monkeypatch, THREE_DAYS)
    with caplog.at_level(logging.INFO, logger="kalshi_common.sgp_runner"):
        sgp_runner.enumerate_kalshi_targets(
            horizon_hours=HORIZON_HOURS)
    lines = [r for r in caplog.records if "[horizon]" in r.getMessage()]
    assert len(lines) == 1
    assert lines[0].levelno == logging.INFO
    assert lines[0].getMessage() == (
        "[horizon] horizon_hours=24.0 kept_games=2 dropped_games=1")


def test_zero_disables_the_window_loudly(monkeypatch, caplog):
    """0 = no window (matches FLIGHT_HORIZON_HOURS' convention), but NEVER
    silently: a WARNING every cycle, because 0 turns the #103 load limiter
    off."""
    _install_fakes(monkeypatch, THREE_DAYS)
    with caplog.at_level(logging.INFO, logger="kalshi_common.sgp_runner"):
        targets = sgp_runner.enumerate_kalshi_targets(horizon_hours=0)
    assert sorted({t.game_id for t in targets}) == ["day_after", "today", "tomorrow"]
    warnings = [r for r in caplog.records
                if r.levelno == logging.WARNING and "DISABLED" in r.getMessage()]
    assert len(warnings) == 1
    assert "[horizon]" in warnings[0].getMessage()


def test_negative_horizon_raises_before_any_kalshi_call(monkeypatch):
    def boom():
        raise AssertionError("Kalshi must not be called on a bad horizon")

    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events", boom)
    with pytest.raises(ValueError, match="horizon_hours must be >= 0"):
        sgp_runner.enumerate_kalshi_targets(horizon_hours=-1)


def test_horizon_is_required_not_defaulted():
    """No hidden default: a caller that forgets the horizon must fail at the
    call, not fall through to an unbounded slate (the pre-#103 accident)."""
    with pytest.raises(TypeError, match="horizon_hours"):
        sgp_runner.enumerate_kalshi_targets()   # type: ignore[call-arg]
    with pytest.raises(TypeError, match="horizon_hours"):
        sgp_runner.target_line_cycle(bot_market_db="unused.duckdb")  # type: ignore[call-arg]


def test_target_line_cycle_threads_horizon_into_enumeration(tmp_path, monkeypatch):
    """The one-filter-one-place contract: `target_line_cycle` passes the
    horizon through and writes exactly what enumeration kept, so warming
    and the taker's parlay cache inherit the bound from the table."""
    import duckdb
    _install_fakes(monkeypatch, THREE_DAYS)
    db_path = str(tmp_path / "market.duckdb")
    out = sgp_runner.target_line_cycle(bot_market_db=db_path,
                                       horizon_hours=HORIZON_HOURS)
    assert sorted({t.game_id for t in out}) == ["today", "tomorrow"]
    con = duckdb.connect(db_path, read_only=True)
    try:
        games = [r[0] for r in con.execute(
            "SELECT DISTINCT game_id FROM mlb_target_lines ORDER BY game_id"
        ).fetchall()]
    finally:
        con.close()
    assert games == ["today", "tomorrow"]
