"""Unit test: fetch_trades must open exactly ONE write_connection per call,
regardless of how many tickers or batches it processes.

Before this fix, fetch_trades opened a write_connection per ticker-batch
AND per series poll_state update — dozens to hundreds per 15-second cycle.
Batching into one connection reduces DuckDB connect/close churn and
within-process contention with the periodic scrape."""
from datetime import datetime
from unittest.mock import MagicMock


class _ConnCounter:
    """Stand-in for write_connection that counts context-manager entries."""
    def __init__(self, fake_conn):
        self.fake_conn = fake_conn
        self.enter_count = 0

    def __call__(self):
        return self

    def __enter__(self):
        self.enter_count += 1
        return self.fake_conn

    def __exit__(self, *args):
        return False


def test_fetch_trades_opens_single_write_connection(monkeypatch):
    """With two series and three tickers total (one batch each),
    fetch_trades must open exactly one write connection for the whole call."""
    from nfl_draft.scrapers import kalshi as kalshi_mod

    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher, "discover_draft_series",
        lambda: [{"series_ticker": "S1"}, {"series_ticker": "S2"}],
    )
    monkeypatch.setattr(
        kalshi_mod, "_tickers_for_series",
        lambda rcon, st: ["T1", "T2"] if st == "S1" else ["T3"],
    )

    class _ReadCtx:
        def __enter__(self):
            m = MagicMock()
            m.execute.return_value.fetchall.return_value = []
            return m
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "read_connection", lambda: _ReadCtx())

    def _fake_public_request(path):
        if "cursor=" in path:
            return None
        return {"trades": [{"fake": True}], "cursor": None}
    monkeypatch.setattr(kalshi_mod, "public_request", _fake_public_request)

    now = datetime.now()
    TradeRow = kalshi_mod.TradeRow

    def _fake_parse(raw):
        if not raw:
            return []
        return [TradeRow(
            trade_id=f"id-{id(raw)}", ticker="T?",
            side="yes", price_cents=50, count=1,
            traded_at=now, fetched_at=now,
        )]
    monkeypatch.setattr(kalshi_mod, "parse_trades_response", _fake_parse)

    fake_conn = MagicMock()
    counter = _ConnCounter(fake_conn)
    monkeypatch.setattr(kalshi_mod, "write_connection", counter)

    kalshi_mod.fetch_trades()

    assert counter.enter_count == 1, (
        f"fetch_trades opened write_connection {counter.enter_count} times; "
        "must be exactly 1 regardless of ticker or batch count."
    )


def test_fetch_trades_writes_trades_and_poll_state_on_shared_conn(monkeypatch):
    """The single write connection must receive both INSERT OR IGNORE INTO
    kalshi_trades and the kalshi_poll_state upsert."""
    from nfl_draft.scrapers import kalshi as kalshi_mod

    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher, "discover_draft_series",
        lambda: [{"series_ticker": "S1"}],
    )
    monkeypatch.setattr(
        kalshi_mod, "_tickers_for_series",
        lambda rcon, st: ["T1"],
    )

    class _ReadCtx:
        def __enter__(self):
            m = MagicMock()
            m.execute.return_value.fetchall.return_value = []
            return m
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "read_connection", lambda: _ReadCtx())

    def _fake_public_request(path):
        if "cursor=" in path:
            return None
        return {"trades": [{"fake": True}], "cursor": None}
    monkeypatch.setattr(kalshi_mod, "public_request", _fake_public_request)

    now = datetime.now()
    TradeRow = kalshi_mod.TradeRow
    monkeypatch.setattr(
        kalshi_mod, "parse_trades_response",
        lambda raw: [TradeRow(
            trade_id="id-1", ticker="T1",
            side="yes", price_cents=50, count=1,
            traded_at=now, fetched_at=now,
        )] if raw else [],
    )

    fake_conn = MagicMock()
    class _CtxWrap:
        def __enter__(self): return fake_conn
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "write_connection", lambda: _CtxWrap())

    kalshi_mod.fetch_trades()

    exec_sqls = [call.args[0] for call in fake_conn.execute.call_args_list]
    assert any("kalshi_trades" in sql for sql in exec_sqls), (
        "Expected INSERT OR IGNORE INTO kalshi_trades; got SQL: " + repr(exec_sqls)
    )
    assert any("kalshi_poll_state" in sql for sql in exec_sqls), (
        "Expected kalshi_poll_state upsert; got SQL: " + repr(exec_sqls)
    )
