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


def test_discovery_pass_uses_constant_db_connections(monkeypatch, tmp_path):
    """2026-08-10 round 2: on the ~1GB live state DB every connection close
    pays a full CHECKPOINT (~0.5s), so the per-RFQ seen-SELECT +
    per-decision INSERT pattern starved the loop even after scope checks
    went wire-free. A pass must open O(1) connections, not O(n_rfqs)."""
    import kalshi_mlb_mm.db as db_mod
    _fresh_db(monkeypatch, tmp_path, "connes.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 10_000)
    main._SCOPE_CACHE.clear()

    real_connect = db_mod.connect
    calls = {"n": 0}

    def counting_connect(*a, **kw):
        calls["n"] += 1
        return real_connect(*a, **kw)

    monkeypatch.setattr(main.db, "connect", counting_connect)
    src = CountingSource(n_rfqs=200)   # 200 brand-new out-of-scope RFQs
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)

    # Pre-pass seen lookup + hoisted open-count/fills reads + one flush —
    # a handful, regardless of RFQ count. 200 RFQs must NOT mean 400+.
    assert calls["n"] <= 8, (
        f"discovery pass opened {calls['n']} DB connections for 200 RFQs — "
        "per-RFQ connect pattern has regressed")
    with db_mod.connect(read_only=True) as con:
        n_dec = con.execute("SELECT COUNT(*) FROM quote_decisions").fetchone()[0]
        n_seen = con.execute("SELECT COUNT(*) FROM seen_rfqs").fetchone()[0]
    assert n_dec == 200 and n_seen == 200, (
        f"buffered writes must all flush: decisions={n_dec} seen={n_seen}")
    main._SCOPE_CACHE.clear()


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


def test_rfq_inline_legs_resolve_scope_with_zero_fetches(monkeypatch, tmp_path):
    """WS rfq_created frames and the REST gap-fill both carry
    mve_selected_legs on the RFQ itself (verified 2026-08-10), so scope must
    resolve without ANY get_market call — the budgeted fetch is only a
    fallback for payloads missing the field."""
    db = _fresh_db(monkeypatch, tmp_path, "inline.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 0)   # any fetch would defer
    main._SCOPE_CACHE.clear()

    class InlineLegsSource:
        def poll(self):
            return [{
                "id": f"r{i}", "market_ticker": f"INLINE-{i}", "contracts": 1,
                "mve_selected_legs": [
                    {"event_ticker": "KXNBAGAME-X", "market_ticker": "NOT-MLB",
                     "side": "yes"},
                ],
            } for i in range(20)]

        def get_market(self, ticker):
            raise AssertionError(
                "scope must come from the RFQ's own mve_selected_legs — "
                f"get_market({ticker}) should never be called")

    main._discovery_tick(InlineLegsSource(), NeverQuoteGateway(), dry_run=True)
    with db.connect(read_only=True) as con:
        n = con.execute(
            "SELECT COUNT(*) FROM quote_decisions WHERE decision='skipped'"
        ).fetchone()[0]
    assert n == 20, f"all 20 inline-legs RFQs must be decided, got {n}"
    main._SCOPE_CACHE.clear()


def test_scope_cache_overflow_evicts_partially_not_clear(monkeypatch, tmp_path):
    _fresh_db(monkeypatch, tmp_path, "cache.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 10_000)
    monkeypatch.setattr(main, "_SCOPE_CACHE_MAX", 8)
    monkeypatch.setattr(main, "_DISCOVERY_CURSOR", 0)   # rotation would reorder inserts
    main._SCOPE_CACHE.clear()

    src = CountingSource(n_rfqs=9)   # one past the cap
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)

    # Cap 8, insert #9 triggers eviction of the oldest 2 (8 // 4): the cache
    # must retain the recent majority, not be wiped to a single survivor.
    assert len(main._SCOPE_CACHE) == 7
    assert "TICK-0" not in main._SCOPE_CACHE     # oldest evicted
    assert "TICK-8" in main._SCOPE_CACHE         # newest present
    main._SCOPE_CACHE.clear()


def test_stale_rfqs_skipped_before_any_work(monkeypatch, tmp_path):
    """MAX_RFQ_AGE_SEC: an RFQ resting past the age cap is skipped before
    scope work, decision logging, or any fetch — quoting it late is adverse
    selection, and classifying it burns lap budget. Missing created_ts
    passes through (fail-open)."""
    db = _fresh_db(monkeypatch, tmp_path, "stale.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "MAX_RFQ_AGE_SEC", 600)
    monkeypatch.setattr(main, "_DISCOVERY_CURSOR", 0)
    main._SCOPE_CACHE.clear()

    class MixedAgeSource:
        def __init__(self):
            self.fetches = 0

        def poll(self):
            from datetime import datetime, timedelta, timezone
            old = (datetime.now(timezone.utc) - timedelta(seconds=3600)).isoformat()
            return (
                [{"id": f"old{i}", "market_ticker": f"OLD-{i}", "contracts": 1,
                  "created_ts": old} for i in range(30)]
                + [{"id": "fresh1", "market_ticker": "FRESH-1", "contracts": 1}]
            )

        def get_market(self, ticker):
            self.fetches += 1
            return {"mve_selected_legs": [{"market_ticker": "NOT-MLB"}]}

    src = MixedAgeSource()
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    with db.connect(read_only=True) as con:
        rows = con.execute("SELECT rfq_id FROM quote_decisions").fetchall()
    assert [r[0] for r in rows] == ["fresh1"], (
        f"only the fresh RFQ may be decided, got {rows}")
    assert src.fetches == 1, "stale RFQs must not trigger market fetches"
    main._SCOPE_CACHE.clear()


def test_pass_timebox_bounds_lap_and_rotates(monkeypatch, tmp_path):
    """DISCOVERY_PASS_BUDGET_SEC: a lap past its wall budget defers the rest
    of the snapshot but always makes progress (>=1 RFQ), and the next lap
    resumes at the rotation cursor instead of re-starting from the head."""
    db = _fresh_db(monkeypatch, tmp_path, "timebox.duckdb")
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "SCOPE_FETCH_BUDGET_PER_TICK", 10_000)
    monkeypatch.setattr(cfg, "DISCOVERY_PASS_BUDGET_SEC", 0.0)   # expire instantly
    monkeypatch.setattr(main, "_DISCOVERY_CURSOR", 0)
    main._SCOPE_CACHE.clear()

    src = CountingSource(n_rfqs=10)
    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    with db.connect(read_only=True) as con:
        n1 = con.execute("SELECT COUNT(*) FROM quote_decisions").fetchone()[0]
    assert n1 == 1, f"zero budget must still process exactly one RFQ, got {n1}"
    assert main._DISCOVERY_CURSOR == 1, (
        f"cursor must advance past the processed RFQ, got {main._DISCOVERY_CURSOR}")

    main._discovery_tick(src, NeverQuoteGateway(), dry_run=True)
    with db.connect(read_only=True) as con:
        ids = sorted(r[0] for r in con.execute(
            "SELECT rfq_id FROM quote_decisions").fetchall())
    assert ids == ["r0", "r1"], (
        f"second lap must resume at the cursor (r1), got {ids}")
    main._SCOPE_CACHE.clear()


def test_coverage_tick_emits_rfq_ingestion_summary(monkeypatch):
    """Non-MLB flow must stay measured after the ingestion filter: the
    coverage arm persists the drop counter + mirror size as a periodic
    `rfq_ingestion_summary` research event (user requirement 2026-08-11 —
    the per-RFQ records are gone, the aggregate must not be)."""
    from kalshi_mlb_mm import main

    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda et, **kw: emitted.append((et, kw)))

    class FakeCoverage:
        def snapshot_and_reset(self):
            return {}

    class FakeSGPService:
        coverage = FakeCoverage()

    class FakeSource:
        ingestion_dropped = 12345
        def poll(self):
            return [{"id": "r1"}, {"id": "r2"}]

    main._coverage_summary_tick(sgp_service=FakeSGPService(),
                                rfq_source=FakeSource())
    kinds = [et for et, _ in emitted]
    assert "rfq_ingestion_summary" in kinds
    payload = dict(emitted[kinds.index("rfq_ingestion_summary")][1])["payload"]
    assert payload["non_mlb_dropped_total"] == 12345
    assert payload["mirror_size"] == 2
