import datetime
from unabated_edge import storage, config


def test_snapshot_flag_and_research(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    storage.snapshot_lines("soccer", [{"ts": now, "event_id": 9, "market_source_id": 7, "bet_type": "bt1", "side": "1", "price": -150.0, "points": None}])
    storage.log_flagged({"ts": now, "sport": "soccer", "event_id": 9, "market_ticker": "X", "outcome": "home",
                         "fair_prob": 0.6, "yes_ask": 0.5, "ev_pct": 0.1, "kelly_contracts": 3, "dry_run": True})
    storage.emit("candidate_priced", "soccer", event_id=9, ev_pct=0.1)
    storage.flush()
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM flagged_edges").fetchone()[0] == 1
    with storage.connect(config.RESEARCH_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM research_events").fetchone()[0] == 1


def _book_row(now, ticker="T-O25"):
    return {"ts": now, "market_ticker": ticker, "floor_strike": 2.5,
            "yes_bid": 0.40, "yes_bid_qty": 120.0, "no_bid": 0.58, "no_bid_qty": 65.0,
            "yes_ask": 0.42, "no_ask": 0.60, "volume": 1000, "open_interest": 400,
            "depth": {"yes_bids": [[0.40, 120.0]], "no_bids": [[0.58, 65.0]]}}


def test_snapshot_books_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    storage.snapshot_books("soccer", [_book_row(now)])
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        row = c.execute("SELECT sport, market_ticker, yes_bid, no_bid_qty, "
                        "json_extract(depth, '$.yes_bids[0][0]') FROM book_snapshots").fetchone()
    assert row[0] == "soccer" and row[1] == "T-O25"
    assert row[2] == 0.40 and row[3] == 65.0
    assert float(row[4]) == 0.40          # full depth ladder survives as queryable JSON


def test_insert_trades_dedups_on_trade_id(tmp_path, monkeypatch):
    """Overlapping poll windows re-deliver the same trades; PK + INSERT OR
    IGNORE must keep exactly one row per trade_id."""
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    t = {"trade_id": "abc", "market_ticker": "T-O25", "created_time": "2026-07-10T17:00:00Z",
         "yes_price": 0.42, "count": 10.0, "taker_side": "yes"}
    storage.insert_trades("soccer", [t], now)
    storage.insert_trades("soccer", [t, {**t, "trade_id": "def"}], now)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM kalshi_trades").fetchone()[0] == 2
