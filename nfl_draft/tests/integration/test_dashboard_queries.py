"""Regression tests for legacy dashboard tabs after the kalshi_odds rename."""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module


@pytest.fixture
def seeded_db(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    # Schema for legacy kalshi_odds (matches migrated structure)
    with db_module.write_connection() as con:
        con.execute("""
            CREATE TABLE kalshi_odds (
                fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
                ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
                yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
                last_price INTEGER, volume BIGINT, volume_24h BIGINT,
                liquidity BIGINT, open_interest INTEGER
            )
        """)
        now = datetime.now()
        con.execute(
            "INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFT1', 'EVT1', 'TKR1', 'Title', 'Cam Ward', 50, 51, 49, 50, 50, 1000, 10000, 5000, 200)",
            [now],
        )
    return tmp_path / "test.duckdb"


def test_get_latest_odds_returns_one_row(seeded_db, monkeypatch):
    # Import after monkeypatching DB_PATH
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    df = legacy_db.get_latest_odds()
    assert df is not None and len(df) == 1
    assert df.iloc[0]["candidate"] == "Cam Ward"


def test_get_price_history_returns_one_row(seeded_db, monkeypatch):
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    df = legacy_db.get_price_history()
    assert df is not None and len(df) == 1


def test_get_snapshot_count_returns_one(seeded_db, monkeypatch):
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    snapshots, first, last = legacy_db.get_snapshot_count()
    assert snapshots == 1
