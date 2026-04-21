import duckdb
import pytest
from pathlib import Path
from nfl_draft.lib import db as db_module
from nfl_draft.lib import migrate_from_kalshi_draft as migrate


@pytest.fixture
def fake_legacy_db(tmp_path):
    """Build a kalshi_draft.duckdb with sample data matching the real schema."""
    legacy_path = tmp_path / "kalshi_draft.duckdb"
    con = duckdb.connect(str(legacy_path))
    con.execute("""
        CREATE TABLE draft_odds (
            fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
            ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
            yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
            last_price INTEGER, volume BIGINT, volume_24h BIGINT,
            liquidity BIGINT, open_interest INTEGER
        )
    """)
    con.execute("INSERT INTO draft_odds VALUES (CURRENT_TIMESTAMP, 'KXNFLDRAFT1', 'EVT1', 'TKR1', 'Title', 'Cam Ward', 50, 51, 49, 50, 50, 1000, 10000, 5000, 200)")
    con.execute("CREATE TABLE draft_series (series_ticker VARCHAR PRIMARY KEY, title VARCHAR, category VARCHAR, discovered_at TIMESTAMP)")
    con.execute("INSERT INTO draft_series VALUES ('KXNFLDRAFT1', 'Test Series', 'NFL', CURRENT_TIMESTAMP)")
    con.execute("CREATE TABLE market_info (ticker VARCHAR PRIMARY KEY, title VARCHAR, subtitle VARCHAR, series_ticker VARCHAR, expiration_time TIMESTAMP, close_time TIMESTAMP, updated_at TIMESTAMP)")
    con.execute("CREATE TABLE positions (fetch_time TIMESTAMP, ticker VARCHAR, position INTEGER, market_exposure DOUBLE, realized_pnl DOUBLE, resting_orders_count INTEGER, total_traded DOUBLE, fees_paid DOUBLE)")
    con.execute("CREATE TABLE resting_orders (fetch_time TIMESTAMP, order_id VARCHAR, ticker VARCHAR, side VARCHAR, type VARCHAR, yes_price INTEGER, no_price INTEGER, remaining_count INTEGER, created_time TIMESTAMP, expiration_time TIMESTAMP)")
    con.execute("CREATE TABLE consensus_board (fetch_time TIMESTAMP, rank INTEGER, player_name VARCHAR, position VARCHAR, school VARCHAR, source VARCHAR)")
    con.execute("CREATE TABLE detected_edges (fetch_time TIMESTAMP, edge_type VARCHAR, description VARCHAR, market_a VARCHAR, market_b VARCHAR, price_a DOUBLE, price_b DOUBLE, implied_edge DOUBLE, confidence VARCHAR)")
    con.close()
    return legacy_path


def test_migration_renames_draft_odds_to_kalshi_odds(fake_legacy_db, tmp_path, monkeypatch):
    new_db = tmp_path / "nfl_draft.duckdb"
    monkeypatch.setattr(db_module, "DB_PATH", new_db)
    db_module.init_schema()
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)
    with db_module.read_connection() as con:
        tables = {row[0] for row in con.execute("SHOW TABLES").fetchall()}
        kalshi_count = con.execute("SELECT COUNT(*) FROM kalshi_odds").fetchone()[0]
    assert "kalshi_odds" in tables
    assert kalshi_count == 1


def test_migration_is_idempotent(fake_legacy_db, tmp_path, monkeypatch):
    new_db = tmp_path / "nfl_draft.duckdb"
    monkeypatch.setattr(db_module, "DB_PATH", new_db)
    db_module.init_schema()
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)  # re-run
    with db_module.read_connection() as con:
        count = con.execute("SELECT COUNT(*) FROM kalshi_odds").fetchone()[0]
    assert count == 1  # not duplicated
