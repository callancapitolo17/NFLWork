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


# ---------------------------------------------------------------------------
# Portal dashboard query tests (Tasks 22-23)
# ---------------------------------------------------------------------------


def test_cross_book_grid_query_outlier_flags(monkeypatch, tmp_path):
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
            "VALUES ('first_qb_cam-ward', 'first_at_position', 'Cam Ward', 'QB')"
        )
        for book, prob in [
            ("kalshi", 0.50),
            ("fanduel", 0.51),
            ("bookmaker", 0.49),
            ("wagerzon", 0.52),
            ("draftkings", 0.30),
        ]:
            con.execute(
                "INSERT INTO draft_odds VALUES ('first_qb_cam-ward', ?, 100, ?, ?, ?)",
                [book, prob, prob, now],
            )
    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=10)
    assert len(rows) == 1
    assert rows[0]["outlier_count"] == 1
    assert rows[0]["flags"]["draftkings"] is True
    assert rows[0]["flags"]["kalshi"] is False


def test_ev_candidates_ranks_by_abs_delta(monkeypatch, tmp_path):
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type) VALUES ('m1', 'prop')")
        con.execute("INSERT INTO draft_markets (market_id, market_type) VALUES ('m2', 'prop')")
        for book, prob in [("kalshi", 0.50), ("dk", 0.50), ("fd", 0.50), ("bm", 0.25)]:
            con.execute("INSERT INTO draft_odds VALUES ('m1', ?, 100, ?, ?, ?)", [book, prob, prob, now])
        for book, prob in [("kalshi", 0.50), ("dk", 0.50), ("fd", 0.50), ("bm", 0.65)]:
            con.execute("INSERT INTO draft_odds VALUES ('m2', ?, 100, ?, ?, ?)", [book, prob, prob, now])
    from nfl_draft.lib.queries import ev_candidates
    rows = ev_candidates(threshold_pp=10)
    assert len(rows) == 2
    assert rows[0]["market_id"] == "m1"
    assert rows[0]["delta"] == pytest.approx(-0.25)
    assert rows[1]["market_id"] == "m2"
    assert rows[1]["delta"] == pytest.approx(0.15)


def test_trade_tape_marks_large_fills(monkeypatch, tmp_path):
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t1', 'TKR', 'yes', 50, 200, 1000.0, ?, ?)",
            [now, now],
        )
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t2', 'TKR', 'yes', 50, 20, 100.0, ?, ?)",
            [now, now],
        )
    from nfl_draft.lib.queries import trade_tape
    rows = trade_tape(limit=10, large_threshold_usd=500.0)
    assert len(rows) == 2
    by_id = {r["trade_id"]: r for r in rows}
    assert by_id["t1"]["is_large"] is True
    assert by_id["t2"]["is_large"] is False


def test_single_venue_market_no_crash(monkeypatch, tmp_path):
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type) VALUES ('solo', 'prop')")
        con.execute("INSERT INTO draft_odds VALUES ('solo', 'kalshi', 100, 0.6, 0.6, ?)", [now])
    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=10)
    assert len(rows) == 1
    assert rows[0]["outlier_count"] == 0
    assert rows[0]["flags"] == {}
