"""Tests for kalshi_mlb_rfq/sgp_runner.py - bot's SGP scrape orchestration."""
from datetime import datetime, timedelta, timezone

import duckdb
import pytest

from kalshi_mlb_rfq.sgp_runner import should_scrape, latest_sgp_fetch_time


def test_should_scrape_no_last_fetch():
    """Empty table / no prior scrape -> always scrape."""
    assert should_scrape(last_fetch_time=None, now=datetime.now(timezone.utc),
                          min_interval_sec=30) is True


def test_should_scrape_recent_fetch():
    """Scraped < min_interval ago -> skip."""
    now = datetime.now(timezone.utc)
    last = now - timedelta(seconds=10)
    assert should_scrape(last, now, min_interval_sec=30) is False


def test_should_scrape_stale_fetch():
    """Scraped > min_interval ago -> scrape."""
    now = datetime.now(timezone.utc)
    last = now - timedelta(seconds=120)
    assert should_scrape(last, now, min_interval_sec=30) is True


def test_should_scrape_normalizes_naive_datetime():
    """Both naive and tz-aware datetimes should work."""
    now = datetime.now(timezone.utc)
    last_naive = (now - timedelta(seconds=120)).replace(tzinfo=None)
    assert should_scrape(last_naive, now, min_interval_sec=30) is True


def test_latest_sgp_fetch_time_empty(tmp_path):
    db = str(tmp_path / "empty.duckdb")
    duckdb.connect(db).close()
    assert latest_sgp_fetch_time(db) is None


def test_latest_sgp_fetch_time_returns_max(tmp_path):
    db = str(tmp_path / "f.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,'2026-05-13 10:00','dk_direct',-1.5,8.5)")
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','B','FG','dk',2.0,100,'2026-05-13 10:05','dk_direct',-1.5,8.5)")
    con.close()
    result = latest_sgp_fetch_time(db)
    assert result is not None
    assert result.hour == 10 and result.minute == 5


def test_enumerate_kalshi_targets_returns_target_lines(monkeypatch, tmp_path):
    """When Kalshi API returns 1 open game with 2 spreads x 2 totals,
    enumerate_kalshi_targets yields 4 TargetLine entries (one per combo)."""
    from kalshi_mlb_rfq import sgp_runner

    # Kalshi event ticker format: YYMMMDDHHMM{AwayCode}{HomeCode} (11-char date prefix).
    # Schedule below has home=NYY, away=BOS, so ticker tail must be BOSNYY.
    fake_event = {
        "event_ticker": "KXMLBGAME-26MAY132300BOSNYY",
        "status": "open",
    }
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events",
                         lambda: [fake_event])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_spread_lines",
                         lambda suffix: [(-1.5, "home"), (-2.5, "home")])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_total_lines",
                         lambda suffix: [8.5, 9.5])
    schedule_db = str(tmp_path / "mlb.duckdb")
    con = duckdb.connect(schedule_db)
    con.execute("""
        CREATE TABLE mlb_odds_temp (
            date VARCHAR, id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            total_line DOUBLE, consensus_prob_home DOUBLE, commence_time VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_odds_temp VALUES ('2026-05-13','g1','New York Yankees','Boston Red Sox',8.5,0.55,'2026-05-13T23:00:00+00:00')")
    con.close()

    targets = sgp_runner.enumerate_kalshi_targets(schedule_db_path=schedule_db)
    assert len(targets) == 4
    spreads = sorted({t.spread for t in targets})
    totals = sorted({t.total for t in targets})
    assert spreads == [-2.5, -1.5]
    assert totals == [8.5, 9.5]
    assert all(t.period == "FG" for t in targets)
    assert all(t.game_id == "g1" for t in targets)


def test_enumerate_kalshi_targets_no_events_returns_empty(monkeypatch):
    from kalshi_mlb_rfq import sgp_runner
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events", lambda: [])
    assert sgp_runner.enumerate_kalshi_targets(schedule_db_path="/nonexistent.duckdb") == []


def test_enumerate_kalshi_targets_skips_unknown_team_codes(monkeypatch, tmp_path):
    """If event_ticker has unknown team codes, skip the event."""
    from kalshi_mlb_rfq import sgp_runner
    fake_event = {"event_ticker": "KXMLBGAME-26MAY132300ZZZXXX"}
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events", lambda: [fake_event])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_spread_lines",
                         lambda suffix: [(-1.5, "home")])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_total_lines", lambda suffix: [8.5])
    schedule_db = str(tmp_path / "mlb.duckdb")
    con = duckdb.connect(schedule_db)
    con.execute("CREATE TABLE mlb_odds_temp (id VARCHAR, home_team VARCHAR, away_team VARCHAR, commence_time VARCHAR)")
    con.close()
    assert sgp_runner.enumerate_kalshi_targets(schedule_db_path=schedule_db) == []
