"""Tests for kalshi_mlb_rfq/main.py cache shape after line-source pivot."""
from datetime import datetime, timezone
import duckdb
import pytest


def test_build_parlay_lines_cache_merges_schedule_and_targets(tmp_path):
    """_build_parlay_lines_cache reads schedule from mlb.duckdb::mlb_odds_temp
    and target lines from bot DB::mlb_target_lines, returns multi-line cache."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache

    sched_db = str(tmp_path / "mlb.duckdb")
    bot_db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)

    con = duckdb.connect(sched_db)
    con.execute("""
        CREATE TABLE mlb_odds_temp (
            date VARCHAR, id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            total_line DOUBLE, consensus_prob_home DOUBLE, commence_time VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_odds_temp VALUES ('2026-05-13','g1','New York Yankees','Boston Red Sox',8.5,0.55,'2026-05-13T23:00:00+00:00')")
    con.close()

    con = duckdb.connect(bot_db)
    con.execute("""
        CREATE TABLE mlb_target_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time TIMESTAMP, period VARCHAR,
            spread DOUBLE, total DOUBLE, written_at TIMESTAMP
        )
    """)
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','NYY','BOS', ?, 'FG', -1.5, 8.5, NOW())", [ct])
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','NYY','BOS', ?, 'FG', -2.5, 8.5, NOW())", [ct])
    con.close()

    cache = _build_parlay_lines_cache(schedule_db=sched_db, bot_db=bot_db)
    assert "g1" in cache
    entry = cache["g1"]
    assert entry["home_team"] == "New York Yankees"
    assert entry["away_team"] == "Boston Red Sox"
    lines_set = set(entry["fg_lines"])
    assert lines_set == {(-1.5, 8.5), (-2.5, 8.5)}


def test_build_parlay_lines_cache_drops_games_missing_from_schedule(tmp_path):
    """Game in target_lines but not in mlb_odds_temp → dropped."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache

    sched_db = str(tmp_path / "mlb.duckdb")
    bot_db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)

    con = duckdb.connect(sched_db)
    con.execute("CREATE TABLE mlb_odds_temp (date VARCHAR, id VARCHAR, home_team VARCHAR, away_team VARCHAR, total_line DOUBLE, consensus_prob_home DOUBLE, commence_time VARCHAR)")
    con.close()  # empty schedule

    con = duckdb.connect(bot_db)
    con.execute("CREATE TABLE mlb_target_lines (game_id VARCHAR, home_team VARCHAR, away_team VARCHAR, commence_time TIMESTAMP, period VARCHAR, spread DOUBLE, total DOUBLE, written_at TIMESTAMP)")
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','NYY','BOS', ?, 'FG', -1.5, 8.5, NOW())", [ct])
    con.close()

    cache = _build_parlay_lines_cache(schedule_db=sched_db, bot_db=bot_db)
    assert cache == {}, "game missing from schedule should be dropped"


def test_build_parlay_lines_cache_empty_dbs(tmp_path):
    """Both DBs missing → empty cache, no errors."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache
    cache = _build_parlay_lines_cache(
        schedule_db=str(tmp_path / "nonexistent_sched.duckdb"),
        bot_db=str(tmp_path / "nonexistent_bot.duckdb"),
    )
    assert cache == {}
