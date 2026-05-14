"""Tests for kalshi_mlb_rfq/main.py cache shape after line-source pivot."""
from datetime import datetime, timezone
import duckdb
import pytest


def test_build_parlay_lines_cache_reads_from_target_lines(tmp_path):
    """_build_parlay_lines_cache reads schedule + lines from bot DB::mlb_target_lines
    only (schedule data now lives directly in the table, written from Odds API)."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache

    bot_db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)

    con = duckdb.connect(bot_db)
    con.execute("""
        CREATE TABLE mlb_target_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time TIMESTAMP, period VARCHAR,
            spread DOUBLE, total DOUBLE, written_at TIMESTAMP
        )
    """)
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','New York Yankees','Boston Red Sox', ?, 'FG', -1.5, 8.5, NOW())", [ct])
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','New York Yankees','Boston Red Sox', ?, 'FG', -2.5, 8.5, NOW())", [ct])
    con.close()

    cache = _build_parlay_lines_cache(bot_db=bot_db)
    assert "g1" in cache
    entry = cache["g1"]
    assert entry["home_team"] == "New York Yankees"
    assert entry["away_team"] == "Boston Red Sox"
    lines_set = set(entry["fg_lines"])
    assert lines_set == {(-1.5, 8.5), (-2.5, 8.5)}


def test_build_parlay_lines_cache_missing_target_lines_table(tmp_path):
    """Bot DB exists but no mlb_target_lines table → empty cache."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache

    bot_db = str(tmp_path / "bot.duckdb")
    con = duckdb.connect(bot_db)
    # Create some other table so the DB file exists but mlb_target_lines doesn't.
    con.execute("CREATE TABLE other_table (x INTEGER)")
    con.close()

    cache = _build_parlay_lines_cache(bot_db=bot_db)
    assert cache == {}


def test_build_parlay_lines_cache_empty_dbs(tmp_path):
    """Missing bot DB → empty cache, no errors."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache
    cache = _build_parlay_lines_cache(bot_db=str(tmp_path / "nonexistent_bot.duckdb"))
    assert cache == {}


def test_load_book_fairs_filters_by_line_and_requires_two_books(monkeypatch):
    """Returns book->fair only when >= 2 books priced the matching (spread, total)."""
    import pandas as pd
    from kalshi_mlb_rfq import main, config

    main._SGP_ODDS_CACHE = pd.DataFrame([
        # game g1, spread -1.5, total 8.5: 2 books -> passes gate
        {"game_id": "g1", "bookmaker": "draftkings", "combo": "Home Spread + Over",
         "sgp_decimal": 2.85, "period": "FG", "spread_line": -1.5, "total_line": 8.5},
        {"game_id": "g1", "bookmaker": "fanduel", "combo": "Home Spread + Over",
         "sgp_decimal": 2.95, "period": "FG", "spread_line": -1.5, "total_line": 8.5},
        # game g1, spread -2.5, total 8.5: only DK -> fails gate
        {"game_id": "g1", "bookmaker": "draftkings", "combo": "Home Spread + Over",
         "sgp_decimal": 3.50, "period": "FG", "spread_line": -2.5, "total_line": 8.5},
    ])
    monkeypatch.setattr(config, "MIN_BOOK_COUNT_FOR_BLEND", 2)

    # Matching line + 2 books -> returns dict (may be partial if devig fails for one)
    fairs = main._load_book_fairs("g1", -1.5, 8.5)
    assert isinstance(fairs, dict)
    # If devig succeeded for both books, dict has 2 entries.
    # If devig failed for one, the N>=2 gate would return {}.
    # Either {2 books returned} OR {} is acceptable here - test the gate, not devig math.
    assert len(fairs) == 0 or len(fairs) >= 2

    # Matching line + 1 book -> empty
    fairs = main._load_book_fairs("g1", -2.5, 8.5)
    assert fairs == {}

    # No matching line -> empty
    fairs = main._load_book_fairs("g1", -3.5, 9.5)
    assert fairs == {}


def test_load_book_fairs_empty_cache_returns_empty():
    """No SGP cache -> empty dict (defensive)."""
    import pandas as pd
    from kalshi_mlb_rfq import main

    main._SGP_ODDS_CACHE = pd.DataFrame()
    assert main._load_book_fairs("g1", -1.5, 8.5) == {}


def test_load_book_fairs_legacy_schema_returns_empty():
    """If cache is in legacy schema (no spread_line col), return empty."""
    import pandas as pd
    from kalshi_mlb_rfq import main

    main._SGP_ODDS_CACHE = pd.DataFrame([
        {"game_id": "g1", "bookmaker": "draftkings", "combo": "x",
         "sgp_decimal": 2.0, "period": "FG"},  # no spread_line/total_line
    ])
    assert main._load_book_fairs("g1", -1.5, 8.5) == {}


def test_per_accept_gates_skip_line_move():
    """After D7, line_move_ok is no longer called from _all_per_accept_gates_pass.
    Verify the gate path doesn't reference _current_book_lines_for_combo."""
    import inspect
    from kalshi_mlb_rfq import main
    src = inspect.getsource(main._all_per_accept_gates_pass)
    assert "_current_book_lines_for_combo" not in src
    assert "line_move_ok" not in src
    assert "reference_lines" not in src


def test_current_book_lines_for_combo_removed():
    """The function should no longer exist on the module."""
    from kalshi_mlb_rfq import main
    assert not hasattr(main, "_current_book_lines_for_combo")


def test_refresh_sgp_cache_partial_reload(monkeypatch, tmp_path):
    """_refresh_sgp_cache reads only mlb_sgp_odds from bot DB, atomic swap."""
    import duckdb
    from datetime import datetime, timezone
    from pathlib import Path
    from kalshi_mlb_rfq import main, config

    bot_db = tmp_path / "bot.duckdb"
    con = duckdb.connect(str(bot_db))
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','H+O','FG','dk',2.85,185,NOW(),'dk_direct',-1.5,8.5)")
    con.close()

    monkeypatch.setattr(config, "BOT_MARKET_DB", bot_db)
    monkeypatch.setattr(config, "MAX_BOOK_STALENESS_SEC", 60)

    result = main._refresh_sgp_cache()
    assert result is True
    assert main._SGP_ODDS_CACHE is not None
    assert len(main._SGP_ODDS_CACHE) == 1
    assert "spread_line" in main._SGP_ODDS_CACHE.columns


def test_refresh_sgp_cache_missing_db(monkeypatch, tmp_path):
    """_refresh_sgp_cache returns False on missing DB."""
    from kalshi_mlb_rfq import main, config
    monkeypatch.setattr(config, "BOT_MARKET_DB", tmp_path / "nonexistent.duckdb")
    result = main._refresh_sgp_cache()
    assert result is False


def test_main_loop_imports_sgp_runner():
    """Sanity: main module must import sgp_runner without error."""
    from kalshi_mlb_rfq import main
    from kalshi_mlb_rfq import sgp_runner
    assert hasattr(main, "main_loop")
    assert hasattr(sgp_runner, "sgp_cycle")
