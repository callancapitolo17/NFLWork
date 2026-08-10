"""Discovery-pass starvation regression tests (2026-08-10 incident).

One _discovery_tick over the WS mirror's exchange-wide open-RFQ snapshot
(4k+ tickers) ran 30+ minutes — one REST get_market per unseen ticker —
starving every other main-loop arm and making SIGTERM wedge until SIGKILL.

Pins the three fixes:
  1. SCOPE_FETCH_BUDGET_PER_TICK bounds REST scope fetches per tick;
     tickers past the budget wait for a later tick.
  2. A cleared _running event breaks the per-RFQ loop promptly.
  3. _SCOPE_CACHE overflow evicts the oldest quarter, never clear().

No network, no live DBs — fake source/gateway, tmp_path DuckDB.
"""
import importlib

import pytest


def _fresh_db(monkeypatch, tmp_path, name):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    return db


class CountingSource:
    """N out-of-scope RFQs; counts get_market calls."""

    def __init__(self, n_rfqs):
        self.fetches = 0
        self._rfqs = [
            {"id": f"r{i}", "market_ticker": f"TICK-{i}", "contracts": 1}
            for i in range(n_rfqs)
        ]

    def poll(self):
        return self._rfqs

    def get_market(self, ticker):
        self.fetches += 1
        return {"mve_selected_legs": [{"market_ticker": "NOT-MLB"}]}


class NeverQuoteGateway:
    def submit_quote(self, *a, **kw):
        raise AssertionError("out-of-scope RFQs must never reach submit_quote")


def test_scope_fetch_budget_bounds_rest_calls_per_tick(monkeypatch, tmp_path):
    _fresh_db(monkeypatch, tmp_path, "budget.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 7)
    main._SCOPE_CACHE.clear()
    src = CountingSource(n_rfqs=50)

    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    assert src.fetches == 7, (
        f"expected exactly the budget (7) of REST fetches in one tick, "
        f"got {src.fetches}")

    # Later ticks keep draining the backlog: 7 more, none re-fetched.
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    assert src.fetches == 14


def test_running_cleared_breaks_pass_promptly(monkeypatch, tmp_path):
    db = _fresh_db(monkeypatch, tmp_path, "sigterm.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 10_000)
    main._SCOPE_CACHE.clear()

    class StopAfterThree(CountingSource):
        # Simulates SIGTERM landing mid-pass: the handler clears _running.
        def get_market(self, ticker):
            if self.fetches == 3:
                main._running.clear()
            return super().get_market(ticker)

    src = StopAfterThree(n_rfqs=500)
    try:
        main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    finally:
        main._running.set()   # never leak a cleared event into other tests
    assert src.fetches <= 4, (
        f"pass must stop within one RFQ of _running.clear(), "
        f"processed {src.fetches} fetches")
    with db.connect(read_only=True) as con:
        n = con.execute("SELECT COUNT(*) FROM quote_decisions").fetchone()[0]
    assert n <= 4


def test_scope_cache_overflow_evicts_partially_not_clear(monkeypatch, tmp_path):
    _fresh_db(monkeypatch, tmp_path, "cache.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 10_000)
    monkeypatch.setattr(main, "_SCOPE_CACHE_MAX", 8)
    main._SCOPE_CACHE.clear()

    src = CountingSource(n_rfqs=9)   # one past the cap
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)

    # Cap 8, insert #9 triggers eviction of the oldest 2 (8 // 4): the cache
    # must retain the recent majority, not be wiped to a single survivor.
    assert len(main._SCOPE_CACHE) == 7
    assert "TICK-0" not in main._SCOPE_CACHE     # oldest evicted
    assert "TICK-8" in main._SCOPE_CACHE         # newest present
    main._SCOPE_CACHE.clear()
