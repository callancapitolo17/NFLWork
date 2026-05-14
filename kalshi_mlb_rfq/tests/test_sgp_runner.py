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
