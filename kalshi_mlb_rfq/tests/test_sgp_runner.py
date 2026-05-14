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


def test_enumerate_kalshi_targets_returns_target_lines(monkeypatch):
    """When Kalshi API returns 1 open game with 2 spreads x 2 totals,
    enumerate_kalshi_targets yields 4 TargetLine entries (one per combo).
    Schedule comes from the Odds API (mocked here)."""
    from datetime import datetime, timezone
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
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)
    fake_schedule = {
        ("New York Yankees", "Boston Red Sox"): {
            "game_id": "g1",
            "home_team": "New York Yankees",
            "away_team": "Boston Red Sox",
            "commence_time": ct,
        }
    }
    monkeypatch.setattr(sgp_runner, "_fetch_schedule_from_odds_api",
                         lambda: fake_schedule)

    targets = sgp_runner.enumerate_kalshi_targets()
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
    assert sgp_runner.enumerate_kalshi_targets() == []


def test_enumerate_kalshi_targets_skips_unknown_team_codes(monkeypatch):
    """If event_ticker has unknown team codes, skip the event."""
    from kalshi_mlb_rfq import sgp_runner
    fake_event = {"event_ticker": "KXMLBGAME-26MAY132300ZZZXXX"}
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events", lambda: [fake_event])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_spread_lines",
                         lambda suffix: [(-1.5, "home")])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_total_lines", lambda suffix: [8.5])
    monkeypatch.setattr(sgp_runner, "_fetch_schedule_from_odds_api", lambda: {})
    assert sgp_runner.enumerate_kalshi_targets() == []


def test_fetch_schedule_from_odds_api_handles_missing_key(monkeypatch, tmp_path):
    """When ODDS_API_KEY is unset (env + ~/.Renviron absent), returns {}."""
    from kalshi_mlb_rfq import sgp_runner
    monkeypatch.delenv("ODDS_API_KEY", raising=False)
    # Point HOME to a tmp dir with no .Renviron so the fallback finds nothing.
    monkeypatch.setenv("HOME", str(tmp_path))
    assert sgp_runner._fetch_schedule_from_odds_api() == {}


def test_write_target_lines_atomic_replace(tmp_path):
    """write_target_lines DELETEs all existing rows then INSERTs new ones,
    in a single transaction."""
    from kalshi_mlb_rfq.sgp_runner import write_target_lines
    from mlb_sgp._shared import TargetLine
    from datetime import datetime, timezone

    db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)
    first = [
        TargetLine("g1", "NYY", "BOS", ct, "FG", -1.5, 8.5),
        TargetLine("g1", "NYY", "BOS", ct, "FG", -1.5, 9.5),
    ]
    write_target_lines(first, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*) FROM mlb_target_lines").fetchone()[0]
    con.close()
    assert out == 2

    # Second write (different lines) → must replace, not append
    second = [TargetLine("g1", "NYY", "BOS", ct, "FG", -2.5, 8.5)]
    write_target_lines(second, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*), MIN(spread) FROM mlb_target_lines").fetchone()
    con.close()
    assert out == (1, -2.5)


def test_write_target_lines_empty_clears_existing(tmp_path):
    """Empty list still clears the table."""
    from kalshi_mlb_rfq.sgp_runner import write_target_lines
    from mlb_sgp._shared import TargetLine
    from datetime import datetime, timezone

    db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)
    write_target_lines([TargetLine("g1", "NYY", "BOS", ct, "FG", -1.5, 8.5)], db_path=db)
    write_target_lines([], db_path=db)
    con = duckdb.connect(db, read_only=True)
    assert con.execute("SELECT COUNT(*) FROM mlb_target_lines").fetchone() == (0,)
    con.close()


def test_run_scrapers_subprocess_plumbing(tmp_path):
    """Verify subprocess invocation: correct env, success path."""
    from kalshi_mlb_rfq.sgp_runner import run_scrapers
    import sys

    scraper_dir = tmp_path / "fake_scrapers"
    scraper_dir.mkdir()
    fake = scraper_dir / "fake_scraper.py"
    fake.write_text(
        "import os\n"
        "open(os.environ['TARGET'], 'w').write(os.environ.get('MLB_SGP_DB_PATH', 'unset'))\n"
    )
    marker = tmp_path / "marker.txt"

    result = run_scrapers(
        scraper_dir=str(scraper_dir),
        scraper_names=["fake_scraper.py"],
        venv_python=sys.executable,
        timeout_sec=10,
        env={"MLB_SGP_DB_PATH": "/tmp/test.duckdb", "TARGET": str(marker)},
    )
    assert result == {"fake_scraper.py": 0}
    assert marker.read_text() == "/tmp/test.duckdb"


def test_run_scrapers_handles_timeout(tmp_path):
    """If a scraper exceeds timeout_sec, kill it and return -1."""
    from kalshi_mlb_rfq.sgp_runner import run_scrapers
    import sys

    scraper_dir = tmp_path / "fake"
    scraper_dir.mkdir()
    fake = scraper_dir / "slow_scraper.py"
    fake.write_text("import time; time.sleep(60)\n")

    result = run_scrapers(
        scraper_dir=str(scraper_dir), scraper_names=["slow_scraper.py"],
        venv_python=sys.executable, timeout_sec=2, env={},
    )
    assert result["slow_scraper.py"] != 0


def test_read_priced_rows_filters_by_staleness(tmp_path):
    from kalshi_mlb_rfq.sgp_runner import read_priced_rows
    db = str(tmp_path / "r.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,NOW(),'dk_direct',-1.5,8.5)")
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,NOW() - INTERVAL 10 MINUTE,'dk_direct',-2.5,8.5)")
    con.close()
    df = read_priced_rows(db, max_age_sec=60)
    assert len(df) == 1
    assert df.iloc[0]["spread_line"] == -1.5


def test_read_priced_rows_missing_db(tmp_path):
    from kalshi_mlb_rfq.sgp_runner import read_priced_rows
    df = read_priced_rows(str(tmp_path / "nonexistent.duckdb"), max_age_sec=60)
    assert df.empty


def test_sgp_cycle_orchestrates_full_tick(monkeypatch, tmp_path):
    """sgp_cycle() composes enumerate → write → run_scrapers → return rc map.
    Stubs out the API + subprocess to verify orchestration ordering."""
    from kalshi_mlb_rfq import sgp_runner

    call_order = []
    def fake_enum():
        call_order.append("enum")
        return []
    def fake_write(targets, db_path):
        call_order.append("write")
    def fake_run(**kw):
        call_order.append("scrape")
        return {}

    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets", fake_enum)
    monkeypatch.setattr(sgp_runner, "write_target_lines", fake_write)
    monkeypatch.setattr(sgp_runner, "run_scrapers", fake_run)

    rcs = sgp_runner.sgp_cycle(
        bot_market_db=str(tmp_path / "bot.duckdb"),
        scraper_dir=str(tmp_path / "scrapers"),
        venv_python="python",
        timeout_sec=60,
    )
    assert call_order == ["enum", "write", "scrape"]
    assert isinstance(rcs, dict)


def test_sgp_cycle_passes_env_to_scrapers(monkeypatch, tmp_path):
    """sgp_cycle must set MLB_SGP_DB_PATH and MLB_SGP_PERIODS=FG in env."""
    from kalshi_mlb_rfq import sgp_runner

    captured_env = {}
    def fake_run(scraper_dir, scraper_names, venv_python, timeout_sec, env=None):
        captured_env.update(env or {})
        return {}

    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets", lambda: [])
    monkeypatch.setattr(sgp_runner, "write_target_lines", lambda targets, db_path: None)
    monkeypatch.setattr(sgp_runner, "run_scrapers", fake_run)

    sgp_runner.sgp_cycle(
        bot_market_db="/tmp/bot.duckdb",
        scraper_dir="/tmp/scrapers",
        venv_python="python",
        timeout_sec=60,
    )
    assert captured_env.get("MLB_SGP_DB_PATH") == "/tmp/bot.duckdb"
    assert captured_env.get("MLB_SGP_PERIODS") == "FG"
