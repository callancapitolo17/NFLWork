import duckdb
import pytest
import tempfile
from pathlib import Path
from nfl_draft.lib import db as db_module


def test_write_connection_opens_and_closes(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    with db_module.write_connection() as con:
        con.execute("CREATE TABLE t (x INT)")
        con.execute("INSERT INTO t VALUES (1)")
    # Connection released — open a new read connection
    with db_module.read_connection() as con:
        result = con.execute("SELECT x FROM t").fetchone()
    assert result == (1,)


def test_init_schema_creates_all_tables(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    with db_module.read_connection() as con:
        tables = {row[0] for row in con.execute("SHOW TABLES").fetchall()}
    expected = {
        "draft_markets", "draft_odds", "kalshi_trades", "kalshi_poll_state",
        "draft_bets", "players", "player_aliases", "teams", "team_aliases",
        "market_map", "draft_odds_unmapped", "draft_odds_unmapped_players",
    }
    assert expected.issubset(tables), f"Missing: {expected - tables}"


def test_init_schema_is_idempotent(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    db_module.init_schema()  # second call must not error
