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


def test_get_latest_odds_returns_one_row(seeded_db):
    # The seeded_db fixture monkeypatches nfl_draft.lib.db.DB_PATH; legacy
    # kalshi_draft.db helpers now route through read_connection() so they
    # automatically pick up that patched location — no extra patch needed.
    from kalshi_draft import db as legacy_db
    df = legacy_db.get_latest_odds()
    assert df is not None and len(df) == 1
    assert df.iloc[0]["candidate"] == "Cam Ward"


def test_get_price_history_returns_one_row(seeded_db):
    from kalshi_draft import db as legacy_db
    df = legacy_db.get_price_history()
    assert df is not None and len(df) == 1


def test_get_snapshot_count_returns_one(seeded_db):
    from kalshi_draft import db as legacy_db
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
    from datetime import datetime, timedelta
    from nfl_draft.lib.db import write_connection
    # Two *different* logical trades — use distinct traded_at so the fill-dedup
    # GROUP BY in trade_tape() keeps them separate. Notional is recomputed at
    # read time (count * taker-side price * 0.01) so pick sizes that cross the
    # $500 large-fill threshold both ways: t1 = 2000*50¢ = $1000 (large),
    # t2 = 40*50¢ = $20 (not large).
    t1 = datetime.now()
    t2 = t1 - timedelta(seconds=30)
    with write_connection() as con:
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t1', 'TKR', 'yes', 50, 2000, 1000.0, ?, ?)",
            [t1, t1],
        )
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t2', 'TKR', 'yes', 50, 40, 20.0, ?, ?)",
            [t2, t2],
        )
    from nfl_draft.lib.queries import trade_tape
    rows = trade_tape(limit=10, large_threshold_usd=500.0)
    assert len(rows) == 2
    by_id = {r["trade_id"]: r for r in rows}
    assert by_id["t1"]["is_large"] is True
    assert by_id["t2"]["is_large"] is False


def test_trade_tape_joins_market_info_and_filters_by_size(monkeypatch, tmp_path):
    """Trade Tape should (a) LEFT JOIN to market_info so rows carry
    human-readable title/subtitle, and (b) hard-filter by ``min_size_usd``
    so small retail fills can be hidden. Tickers with no market_info row
    must still appear (LEFT JOIN) with title/subtitle = None.
    """
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "t.duckdb")
    db_module.init_schema()

    from datetime import datetime, timedelta
    from nfl_draft.lib.db import write_connection
    with write_connection() as con:
        # init_schema creates market_info; just seed a row for TKR-1.
        con.execute(
            "INSERT INTO market_info (ticker, title, subtitle) "
            "VALUES ('TKR-1', 'Who will be #1?', 'Fernando Mendoza')"
        )
        # Three *distinct* trades: two big on the known ticker, one small on
        # an unknown ticker (to verify LEFT JOIN returns nulls). Use distinct
        # traded_at timestamps so the fill-dedup GROUP BY keeps them separate.
        now = datetime.now()
        t1 = now
        t2 = now - timedelta(seconds=30)
        t3 = now - timedelta(seconds=60)
        # NB: t3 is side='no' with yes_price=1; at read time the query flips
        # it to taker-side (100-1 = 99). Taker notional = 1000 * 99¢ = $990.
        # That's intentional — we're verifying taker-side math + LEFT JOIN.
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t1', 'TKR-1', 'yes', 99, 100, 99.0, ?, ?)",
            [t1, t1],
        )
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t2', 'TKR-1', 'yes', 99, 500, 495.0, ?, ?)",
            [t2, t2],
        )
        con.execute(
            "INSERT INTO kalshi_trades VALUES ('t3', 'TKR-2-UNKNOWN', 'no', 1, 1000, 10.0, ?, ?)",
            [t3, t3],
        )
    from nfl_draft.lib.queries import trade_tape

    # No filter: all 3 trades returned, names populated where the JOIN hits.
    all_rows = trade_tape(limit=10, min_size_usd=0)
    assert len(all_rows) == 3
    for r in all_rows:
        if r["ticker"] == "TKR-1":
            assert r["market_title"] == "Who will be #1?"
            assert r["market_subtitle"] == "Fernando Mendoza"
        else:
            # LEFT JOIN miss — title/subtitle come back null.
            assert r["market_title"] is None
            assert r["market_subtitle"] is None

    # Filter at $100: t2 ($495 taker notional) and t3 (1000 * 99¢ taker
    # notional = $990 after the side-flip) both survive; t1 ($99) is below.
    big = trade_tape(limit=10, min_size_usd=100)
    assert {r["trade_id"] for r in big} == {"t2", "t3"}


def test_trade_tape_aggregates_multi_fill_trades(monkeypatch, tmp_path):
    """Multiple Kalshi trade_ids with same (ticker, traded_at, price, side) should collapse to 1 row.

    Regression guard for the Trade Tape fill-dedup fix: Kalshi emits each
    order-match as its own trade_id, so a 1000-contract taker order that
    matched 3 resting makers shows up in the raw table as 3 rows. The query
    must GROUP BY the natural trade key and SUM count / notional so the
    dashboard shows 1 row per logical trade.
    """
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "t.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    t = datetime(2026, 4, 20, 14, 0, 0)
    with write_connection() as con:
        # 3 fills of a single 1000-contract trade
        con.execute("INSERT INTO kalshi_trades VALUES ('t1', 'TKR', 'yes', 50, 400, 200.0, ?, ?)", [t, t])
        con.execute("INSERT INTO kalshi_trades VALUES ('t2', 'TKR', 'yes', 50, 300, 150.0, ?, ?)", [t, t])
        con.execute("INSERT INTO kalshi_trades VALUES ('t3', 'TKR', 'yes', 50, 300, 150.0, ?, ?)", [t, t])
    from nfl_draft.lib.queries import trade_tape
    rows = trade_tape(limit=10, min_size_usd=0)
    assert len(rows) == 1  # collapsed
    r = rows[0]
    assert r["count"] == 1000
    assert r["price_cents"] == 50  # yes side, same as yes_price
    assert r["notional_usd"] == 500.0  # 1000 * 50¢


def test_trade_tape_taker_side_pricing_for_no(monkeypatch, tmp_path):
    """When side=no, price should be 100 - yes_price and notional should match taker cost.

    Regression guard for the taker-side pricing fix: the DB stores YES-side
    price regardless of taker side, so a NO buyer paying 1¢ shows up in the
    raw row as ``price_cents=99``. The query must flip it at read time.
    """
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "t.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    t = datetime(2026, 4, 20, 14, 0, 0)
    with write_connection() as con:
        # NO buyer: yes_price=99 (NO cost=1), count=1000, notional=$10
        con.execute("INSERT INTO kalshi_trades VALUES ('tX', 'TKR', 'no', 99, 1000, 10.0, ?, ?)", [t, t])
    from nfl_draft.lib.queries import trade_tape
    rows = trade_tape(limit=10, min_size_usd=0)
    assert len(rows) == 1
    r = rows[0]
    assert r["side"] == "no"
    assert r["price_cents"] == 1  # taker-side (100 - 99)
    assert abs(r["notional_usd"] - 10.0) < 0.01  # 1000 × 1¢


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


# ---------------------------------------------------------------------------
# Lock-contention graceful-degradation tests (C-3)
# ---------------------------------------------------------------------------


def test_queries_return_locked_sentinel_on_lock_error(monkeypatch):
    """When a DuckDB lock error fires, query functions return QueryLocked
    rather than crashing. Dash callbacks interpret this as PreventUpdate.

    We simulate the lock error by patching the duckdb.connect used inside
    read_connection() to raise an IOException with the classic lock-error
    signature.
    """
    import duckdb
    from nfl_draft.lib import db as db_module
    from nfl_draft.lib import queries

    def fake_connect(*args, **kwargs):
        raise duckdb.IOException(
            "IO Error: Could not set lock on file "
            "\"/tmp/x.duckdb\": Conflicting lock is held in /tmp/x.duckdb"
        )

    # read_connection lives in db_module; it calls duckdb.connect via the
    # duckdb module imported at the top of that file.
    monkeypatch.setattr(db_module.duckdb, "connect", fake_connect)

    assert isinstance(queries.cross_book_grid(threshold_pp=10), queries.QueryLocked)
    assert isinstance(queries.ev_candidates(threshold_pp=10), queries.QueryLocked)
    assert isinstance(queries.trade_tape(), queries.QueryLocked)
    assert isinstance(queries.bet_log_rows(), queries.QueryLocked)
    assert isinstance(queries.all_market_ids(), queries.QueryLocked)
    assert isinstance(queries.latest_max_fetched_at("draft_odds"), queries.QueryLocked)


def test_kalshi_scrape_writes_to_legacy_kalshi_odds(monkeypatch, tmp_path):
    """fetch_draft_odds must write a fresh snapshot to kalshi_odds on every
    scrape — the legacy Price History / Edge Detection / Consensus tabs
    depend on that write path. Regression guard for the Feb 19 staleness
    bug: the adapter used to rely on a side effect of fetch_markets_for_series,
    but that helper only fetches. Writes now happen explicitly in the adapter
    via _write_legacy_kalshi_odds.
    """
    # Redirect DB to a temp file so we don't touch the real portal DB.
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "t.duckdb")
    db_module.init_schema()

    # Legacy kalshi_odds + market_info + draft_series tables aren't created
    # by init_schema — they're legacy shapes. Seed them with the same
    # shapes the production DB uses (no PKs, per DESCRIBE).
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
        # market_info is now created by init_schema (IF NOT EXISTS), so
        # this CREATE would collide. Skip — init_schema already made it.
        con.execute("""
            CREATE TABLE draft_series (
                series_ticker VARCHAR, title VARCHAR, category VARCHAR,
                discovered_at TIMESTAMP
            )
        """)

    # Mock the network layer + series discovery.
    from nfl_draft.scrapers import kalshi as kalshi_mod
    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher,
        "discover_draft_series",
        lambda: [{"series_ticker": "KXFAKE", "title": "Fake Series"}],
    )
    fake_response = {
        "markets": [
            {
                "ticker": "KXFAKE-TKR1",
                "event_ticker": "KXFAKE-EVT1",
                "title": "Mock Market",
                "yes_sub_title": "Mock Candidate",
                "yes_bid": 50, "yes_ask": 51,
                "no_bid": 49, "no_ask": 50,
                "last_price": 50, "volume": 100,
                "volume_24h": 500, "liquidity_dollars": "250.0",
                "open_interest": 75,
            },
        ],
        "cursor": None,
    }
    monkeypatch.setattr(
        kalshi_mod, "public_request", lambda path: fake_response
    )

    # Run the scrape.
    rows = kalshi_mod.fetch_draft_odds()

    # OddsRow output sanity.
    assert len(rows) == 1
    assert rows[0].book == "kalshi"

    # Legacy kalshi_odds write happened.
    from datetime import datetime, timedelta
    with db_module.write_connection() as con:
        count = con.execute("SELECT COUNT(*) FROM kalshi_odds").fetchone()[0]
        last_fetch = con.execute("SELECT MAX(fetch_time) FROM kalshi_odds").fetchone()[0]
    assert count == 1, f"expected 1 row in kalshi_odds, got {count}"
    assert last_fetch is not None
    # fetch_time should be ~now (within the last minute).
    assert datetime.now() - last_fetch < timedelta(minutes=1), \
        f"kalshi_odds fetch_time {last_fetch} is not recent"


def test_cross_book_grid_excludes_stale_rows(monkeypatch, tmp_path):
    """Rows older than MAX_AGE_HOURS must be filtered out of cross_book_grid.

    Regression guard for the pre-draft audit finding: FanDuel hadn't scraped
    in 22h but its rows still polluted the grid. After the fix, a row with
    fetched_at 3 days old should be excluded entirely, leaving only the
    fresh row in the per-market book map.
    """
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "t.duckdb")
    db_module.init_schema()

    from datetime import datetime, timedelta
    from nfl_draft.lib.db import write_connection

    now = datetime.now()
    stale = now - timedelta(days=3)
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type) VALUES ('m1', 'prop')"
        )
        # Fresh row: kalshi @ 0.50
        con.execute(
            "INSERT INTO draft_odds VALUES ('m1', 'kalshi', 100, 0.5, 0.5, ?)",
            [now],
        )
        # Stale row: draftkings @ 0.30 (would flag as outlier if included)
        con.execute(
            "INSERT INTO draft_odds VALUES ('m1', 'draftkings', 100, 0.3, 0.3, ?)",
            [stale],
        )

    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=10)

    assert len(rows) == 1
    # Only kalshi should be present; draftkings row is stale and excluded.
    assert list(rows[0]["books"].keys()) == ["kalshi"]


def test_queries_propagate_non_lock_errors(monkeypatch):
    """A non-lock DuckDB error (e.g. missing table) must still raise so
    genuine bugs surface loudly — only lock-family errors get swallowed."""
    import duckdb
    from nfl_draft.lib import db as db_module
    from nfl_draft.lib import queries

    def fake_connect(*args, **kwargs):
        # A schema error, not a lock error
        raise duckdb.Error("Catalog Error: Table does not exist")

    monkeypatch.setattr(db_module.duckdb, "connect", fake_connect)

    import pytest as _pytest
    with _pytest.raises(duckdb.Error):
        queries.cross_book_grid(threshold_pp=10)
