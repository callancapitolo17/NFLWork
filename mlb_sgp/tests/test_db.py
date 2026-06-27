"""Tests for mlb_sgp/db.py — schema migration + writers."""
import importlib

import duckdb
import pytest

import mlb_sgp.db as db_mod
from mlb_sgp.db import ensure_table


@pytest.fixture(autouse=True)
def _reset_db_module():
    """Restore db_mod.MLB_DB to its default after each test.

    Several tests use importlib.reload(db_mod) to exercise env-var-driven
    module init. Without this fixture, the last reload leaves MLB_DB
    pointing at a tmp_path that no longer exists, leaking state to any
    later test that imports MLB_DB.
    """
    yield
    importlib.reload(db_mod)


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


from datetime import datetime, timezone
from mlb_sgp._shared import PricedRow
from mlb_sgp.db import upsert_priced_rows, ensure_table


def test_upsert_priced_rows_inserts(tmp_path):
    db = str(tmp_path / "u.duckdb")
    ensure_table(db_path=db)
    rows = [
        PricedRow(
            game_id="g1", combo="Home Spread + Over", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=2.85, sgp_american=185,
            fetch_time=datetime.now(timezone.utc),
        ),
        PricedRow(
            game_id="g1", combo="Home Spread + Under", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=3.20, sgp_american=220,
            fetch_time=datetime.now(timezone.utc),
        ),
    ]
    upsert_priced_rows(rows, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute(
        "SELECT game_id, combo, spread_line, total_line, sgp_american FROM mlb_sgp_odds ORDER BY combo"
    ).fetchall()
    con.close()
    assert len(out) == 2
    assert out[0] == ("g1", "Home Spread + Over", -1.5, 8.5, 185)
    assert out[1] == ("g1", "Home Spread + Under", -1.5, 8.5, 220)


def test_upsert_priced_rows_replaces_same_key(tmp_path):
    """Same (game, combo, period, spread, total, bookmaker, source) → replace."""
    db = str(tmp_path / "r.duckdb")
    ensure_table(db_path=db)
    def make(price):
        return PricedRow(
            game_id="g1", combo="Home Spread + Over", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=price, sgp_american=185,
            fetch_time=datetime.now(timezone.utc),
        )
    upsert_priced_rows([make(2.50)], db_path=db)
    upsert_priced_rows([make(2.75)], db_path=db)
    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*), MAX(sgp_decimal) FROM mlb_sgp_odds").fetchone()
    con.close()
    assert out == (1, 2.75)


def test_upsert_priced_rows_empty_is_noop(tmp_path):
    db = str(tmp_path / "e.duckdb")
    ensure_table(db_path=db)
    upsert_priced_rows([], db_path=db)  # must not error


def test_upsert_dedups_null_spread_line(tmp_path):
    """Moneyline×total rows carry spread_line=None. Re-upserting the same key
    must REPLACE (NULL-safe dedup), not accumulate a duplicate every cycle."""
    from datetime import datetime, timezone
    from mlb_sgp._shared import PricedRow
    from mlb_sgp.db import upsert_priced_rows

    db = str(tmp_path / "ml.duckdb")

    def row(dec):
        return PricedRow(
            game_id="g1", combo="Home ML + Over", period="FG",
            spread_line=None, total_line=8.5, bookmaker="draftkings",
            source="draftkings_direct", sgp_decimal=dec, sgp_american=100,
            fetch_time=datetime(2026, 6, 23, tzinfo=timezone.utc))

    upsert_priced_rows([row(2.4)], db_path=db)
    upsert_priced_rows([row(2.7)], db_path=db)   # same key, fresher price
    con = duckdb.connect(db, read_only=True)
    n = con.execute(
        "SELECT COUNT(*) FROM mlb_sgp_odds WHERE combo='Home ML + Over'").fetchone()[0]
    price = con.execute(
        "SELECT sgp_decimal FROM mlb_sgp_odds WHERE combo='Home ML + Over'").fetchone()[0]
    sl = con.execute(
        "SELECT spread_line FROM mlb_sgp_odds WHERE combo='Home ML + Over'").fetchone()[0]
    con.close()
    assert n == 1            # deduped, not duplicated
    assert price == 2.7      # replaced with the newer price
    assert sl is None        # stored as NULL, no sentinel leak


def test_upsert_mixed_null_and_nonnull_batch(tmp_path):
    """A single batch mixing spread+total rows (non-NULL lines) and ML×total rows
    (NULL spread) must dedup each independently."""
    from datetime import datetime, timezone
    from mlb_sgp._shared import PricedRow
    from mlb_sgp.db import upsert_priced_rows

    db = str(tmp_path / "mix.duckdb")
    ft = datetime(2026, 6, 23, tzinfo=timezone.utc)

    def batch(dec):
        return [
            PricedRow("g1", "Home Spread + Over", "FG", -1.5, 8.5,
                      "draftkings", "draftkings_direct", dec, 100, ft),
            PricedRow("g1", "Home ML + Over", "FG", None, 8.5,
                      "draftkings", "draftkings_direct", dec, 100, ft),
        ]

    upsert_priced_rows(batch(2.0), db_path=db)
    upsert_priced_rows(batch(2.5), db_path=db)
    con = duckdb.connect(db, read_only=True)
    total = con.execute("SELECT COUNT(*) FROM mlb_sgp_odds").fetchone()[0]
    con.close()
    assert total == 2        # both keys deduped, one row each
