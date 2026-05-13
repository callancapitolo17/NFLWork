"""Tests for mlb_sgp/db.py — schema migration + writers."""
import importlib

import duckdb

import mlb_sgp.db as db_mod
from mlb_sgp.db import ensure_table


def test_ensure_table_fresh_db(tmp_path):
    db = str(tmp_path / "fresh.duckdb")
    ensure_table(db_path=db)
    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
    assert "game_id" in cols


def test_ensure_table_migrates_legacy_schema(tmp_path):
    """If the table exists with the old (no spread/total cols) schema, ALTER it."""
    db = str(tmp_path / "legacy.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR,
            bookmaker VARCHAR, sgp_decimal DOUBLE, sgp_american INTEGER,
            fetch_time TIMESTAMP, source VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','Home Spread + Over','FG','draftkings',2.85,185,NOW(),'draftkings_direct')")
    con.close()

    ensure_table(db_path=db)

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    legacy_row = con.execute("SELECT spread_line, total_line FROM mlb_sgp_odds").fetchone()
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
    assert legacy_row == (None, None), "pre-existing rows have NULL line cols"


def test_ensure_table_idempotent(tmp_path):
    db = str(tmp_path / "idem.duckdb")
    ensure_table(db_path=db)
    ensure_table(db_path=db)  # second call must not error
    ensure_table(db_path=db)  # third call must not error either


def test_mlb_db_default_when_env_unset(monkeypatch):
    monkeypatch.delenv("MLB_SGP_DB_PATH", raising=False)
    importlib.reload(db_mod)
    assert str(db_mod.MLB_DB).endswith("mlb_mm.duckdb")


def test_mlb_db_uses_env_override(monkeypatch, tmp_path):
    override = str(tmp_path / "override.duckdb")
    monkeypatch.setenv("MLB_SGP_DB_PATH", override)
    importlib.reload(db_mod)
    assert str(db_mod.MLB_DB) == override
