"""Tests for migration 003: composite PK on placed_bets."""
from __future__ import annotations
import importlib.util
from pathlib import Path

import duckdb
import pytest


def _load_migration():
    """Import the migration module by file path (folder name has a space)."""
    here = Path(__file__).parent
    spec = importlib.util.spec_from_file_location(
        "_m003", here / "003_placed_bets_composite_pk.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Schema mirrors the live placed_bets after migration 002.
OLD_SCHEMA_DDL = """
CREATE TABLE placed_bets (
    bet_hash         VARCHAR PRIMARY KEY,
    game_id          VARCHAR,
    home_team        VARCHAR,
    away_team        VARCHAR,
    market           VARCHAR,
    bet_on           VARCHAR,
    line             DOUBLE,
    odds             INTEGER,
    actual_size      DOUBLE,
    recommended_size DOUBLE,
    bookmaker        VARCHAR,
    model_prob       DOUBLE,
    model_ev         DOUBLE,
    account          VARCHAR,
    status           VARCHAR DEFAULT 'placed',
    ticket_number    VARCHAR,
    error_msg        VARCHAR,
    error_msg_key    VARCHAR,
    wz_odds_at_place INTEGER
)
"""


def _seed_old_db(path: Path) -> None:
    con = duckdb.connect(str(path))
    try:
        con.execute(OLD_SCHEMA_DDL)
    finally:
        con.close()


def _pk_cols(path: Path) -> set[str]:
    con = duckdb.connect(str(path))
    try:
        rows = con.execute("""
            SELECT UNNEST(constraint_column_names)
            FROM duckdb_constraints()
            WHERE table_name = 'placed_bets' AND constraint_type = 'PRIMARY KEY'
        """).fetchall()
    finally:
        con.close()
    return {r[0] for r in rows}


def test_backfills_null_and_empty_account_to_wagerzon(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", None, "placed"])
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h2", "", "placed"])
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h3", "WagerzonJ", "placed"])
    con.close()

    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    rows = sorted(con.execute("SELECT bet_hash, account FROM placed_bets").fetchall())
    con.close()
    assert rows == [("h1", "Wagerzon"), ("h2", "Wagerzon"), ("h3", "WagerzonJ")]


def test_promotes_pk_to_composite(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    _load_migration().run(str(db))
    assert _pk_cols(db) == {"bet_hash", "account"}


def test_same_hash_different_account_is_allowed(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    con.close()

    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "WagerzonJ", "placed"])
    count = con.execute("SELECT COUNT(*) FROM placed_bets WHERE bet_hash = 'h1'").fetchone()[0]
    con.close()
    assert count == 2


def test_duplicate_hash_and_account_is_rejected(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                    ["h1", "Wagerzon", "placed"])
    con.close()


def test_idempotent(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    con.close()

    migrate = _load_migration().run
    migrate(str(db))
    migrate(str(db))   # second run must be a no-op

    con = duckdb.connect(str(db))
    count = con.execute("SELECT COUNT(*) FROM placed_bets").fetchone()[0]
    con.close()
    assert count == 1
    assert _pk_cols(db) == {"bet_hash", "account"}
