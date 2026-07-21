"""Issue #18 — correctness batch (audit B4, B5, O-1, O-2):
1. B4: accept-time fill size parses `contracts_fp` (fixed-point string), never
   defaults to 1; unparseable size records 0 conservatively.
2. B5: a transient `get_market` failure must NOT poison the scope cache; the
   cache is size-bounded.
3. O-1: `_commence_time` / `_resolve_game_for_legs` results are cached in
   memory and invalidated on SGP refresh (no DuckDB connect per call).
4. O-2: maker tables migrate to TIMESTAMPTZ; instants preserved; idempotent.
"""
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# 1. B4 — accept-time fill size
# ---------------------------------------------------------------------------
def test_accept_fill_contracts_parses_fp_string():
    from kalshi_mlb_mm.main import _accept_fill_contracts
    assert _accept_fill_contracts({"contracts_fp": "25.00"}) == 25.0
    # fp takes precedence over the legacy int field
    assert _accept_fill_contracts({"contracts_fp": "21.50", "contracts": 3}) == 21.5
    # fp zero/absent falls back to legacy int
    assert _accept_fill_contracts({"contracts_fp": "0.00", "contracts": 3}) == 3.0
    assert _accept_fill_contracts({"contracts": 5}) == 5.0
    # unknown size → None (caller records conservatively, never 1)
    assert _accept_fill_contracts({"contracts_fp": "garbage"}) is None
    assert _accept_fill_contracts({}) is None
    assert _accept_fill_contracts(None) is None


def _confirm_env(monkeypatch, tmp_path, db_name, quote_response):
    """Minimal accepted-quote confirm-path env (mirrors the N5 fast-fill test)."""
    import importlib, json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.router as router_mod
    from kalshi_mlb_mm import main
    from kalshi_common import auth_client

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    importlib.reload(db)
    db.init_database()

    _evt = "KXMLBGAME-25JUN271905TEXLAA"
    legs = [{"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
             "event_ticker": _evt, "side": "yes"},
            {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
             "event_ticker": _evt, "side": "yes"}]
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-b4", "r-b4", "COMBO-B4", "g4", 0.45, 0.45, None, 0.55, 0.55,
             "open", datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs "
            "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
            "first_seen_at, last_decision) VALUES (?,?,?,?,?,?,?)",
            ["r-b4", "COMBO-B4", True, "g4", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    monkeypatch.setattr(router_mod, "combo_fair", lambda *a, **k: 0.56)

    def fake_api(method, path, *a, **kw):
        if path.startswith("/communications/quotes/"):
            return 200, {"quote": quote_response}, None
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    class GW:
        def confirm(self, qid):
            return True

        def cancel(self, qid):
            return True

    return db, main, GW()


def test_confirm_records_contracts_fp_size(monkeypatch, tmp_path):
    """THE regression test for B4: quote-status response denominates the fill
    as contracts_fp:"25.00" (no legacy int field) → 25 must be recorded.
    Pre-fix code read only `contracts` and defaulted to 1."""
    db, main, gw = _confirm_env(
        monkeypatch, tmp_path, "b4_fp.duckdb",
        {"status": "accepted", "accepted_side": "no", "contracts_fp": "25.00"})
    main._confirm_tick(gw, dry_run=False)
    with db.connect(read_only=True) as con:
        ct = con.execute(
            "SELECT contracts FROM fills WHERE quote_id='qid-b4'").fetchone()[0]
    assert ct == 25.0, f"expected 25 contracts from contracts_fp, got {ct}"


def test_confirm_records_zero_when_size_unparseable(monkeypatch, tmp_path):
    """No parseable size in the response → record 0 (NOT 1) with
    reconciled=FALSE, so the N8 convention counts it at the per-fill cap and
    the reconcile sweep sets the true size from Kalshi positions."""
    db, main, gw = _confirm_env(
        monkeypatch, tmp_path, "b4_unparse.duckdb",
        {"status": "accepted", "accepted_side": "no"})
    main._confirm_tick(gw, dry_run=False)
    with db.connect(read_only=True) as con:
        ct, rec = con.execute(
            "SELECT contracts, reconciled FROM fills WHERE quote_id='qid-b4'").fetchone()
    assert ct == 0.0, f"unparseable size must record 0, got {ct}"
    assert rec is False


# ---------------------------------------------------------------------------
# 2. B5 — scope-cache poisoning + bound
# ---------------------------------------------------------------------------
def _scope_env(monkeypatch, tmp_path, db_name):
    import importlib
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    import pandas as pd

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(main, "_SCOPE_CACHE", {})
    return db, main


class _NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("should not submit in this test")


def test_transient_market_fetch_failure_not_cached(monkeypatch, tmp_path):
    """THE regression test for B5: get_market returns None (HTTP blip) → the
    verdict must NOT be cached (pre-fix: in_scope=False cached forever) and the
    RFQ must not be recorded out-of-scope. A later successful fetch decides."""
    db, main = _scope_env(monkeypatch, tmp_path, "b5_transient.duckdb")

    class FailingSrc:
        def poll(self):
            return [{"id": "r-b5", "market_ticker": "COMBO-B5", "contracts": 1}]

        def get_market(self, t):
            return None

    main._discovery_tick(FailingSrc(), _NoQuoteGW(), dry_run=True)
    assert "COMBO-B5" not in main._SCOPE_CACHE, \
        "a transient fetch failure must not poison the scope cache"
    with db.connect(read_only=True) as con:
        seen = con.execute(
            "SELECT in_scope FROM seen_rfqs WHERE rfq_id='r-b5'").fetchone()
    assert seen is None, "fetch failure must not record the RFQ as decided"

    # Later tick with a working fetch reaches a real verdict and caches it.
    class WorkingSrc(FailingSrc):
        def get_market(self, t):
            return {"mve_selected_legs": [{"market_ticker": "FOO"}]}

    main._discovery_tick(WorkingSrc(), _NoQuoteGW(), dry_run=True)
    assert "COMBO-B5" in main._SCOPE_CACHE, \
        "a successful fetch must cache its verdict"


def test_scope_cache_size_bounded(monkeypatch, tmp_path):
    db, main = _scope_env(monkeypatch, tmp_path, "b5_bound.duckdb")
    monkeypatch.setattr(main, "_SCOPE_CACHE_MAX", 3)
    main._SCOPE_CACHE.update({f"T{i}": (False, None, None) for i in range(3)})

    class Src:
        def poll(self):
            return [{"id": "r-bound", "market_ticker": "COMBO-BOUND", "contracts": 1}]

        def get_market(self, t):
            return {"mve_selected_legs": [{"market_ticker": "FOO"}]}

    main._discovery_tick(Src(), _NoQuoteGW(), dry_run=True)
    assert "COMBO-BOUND" in main._SCOPE_CACHE
    assert len(main._SCOPE_CACHE) <= 3, \
        f"scope cache must stay bounded, got {len(main._SCOPE_CACHE)}"


# ---------------------------------------------------------------------------
# 3. O-1 — game-metadata caches (no DuckDB connect per call)
# ---------------------------------------------------------------------------
def test_commence_time_cached_until_sgp_refresh(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main

    main._invalidate_game_caches()
    calls = []
    ct = datetime(2026, 7, 20, 23, 0, tzinfo=timezone.utc)
    monkeypatch.setattr(main, "_commence_time_uncached",
                        lambda gid: (calls.append(gid), ct)[1])

    assert main._commence_time("g1") == ct
    assert main._commence_time("g1") == ct
    assert calls == ["g1"], "second call must hit the cache"

    # SGP refresh invalidates (MARKET_DB absent → refresh early-returns AFTER
    # clearing the caches).
    monkeypatch.setattr(cfg, "MARKET_DB", tmp_path / "absent.duckdb")
    main._refresh_sgp()
    assert main._commence_time("g1") == ct
    assert calls == ["g1", "g1"], "SGP refresh must invalidate the cache"


def test_commence_time_failure_not_cached(monkeypatch):
    from kalshi_mlb_mm import main
    main._invalidate_game_caches()
    calls = []
    monkeypatch.setattr(main, "_commence_time_uncached",
                        lambda gid: (calls.append(gid), None)[1])
    assert main._commence_time("g1") is None
    assert main._commence_time("g1") is None
    assert calls == ["g1", "g1"], "a None (lookup failure) must not be cached"


def test_resolve_game_cached_by_event_ticker(monkeypatch):
    from kalshi_mlb_mm import main
    from kalshi_common.legset import CanonicalLeg
    main._invalidate_game_caches()
    calls = []
    monkeypatch.setattr(main, "_resolve_game_for_legs_uncached",
                        lambda gl: (calls.append(gl[0].game_id), "gid-1")[1])
    legs = [CanonicalLeg("KXMLBGAME-25JUL201905TEXLAA", "ml", None, "home")]
    assert main._resolve_game_for_legs(legs) == "gid-1"
    assert main._resolve_game_for_legs(legs) == "gid-1"
    assert calls == ["KXMLBGAME-25JUL201905TEXLAA"], "second call must hit the cache"


# ---------------------------------------------------------------------------
# 4. O-2 — TIMESTAMPTZ migration
# ---------------------------------------------------------------------------
def test_timestamptz_migration_idempotent_and_instant_preserving(monkeypatch, tmp_path):
    import importlib
    import duckdb
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db

    db_path = tmp_path / "migrate.duckdb"
    monkeypatch.setattr(cfg, "DB_PATH", db_path)
    importlib.reload(db)

    # Simulate a pre-migration DB: old naive-TIMESTAMP fills table with one row
    # written the old way (aware UTC bind → stored naive-LOCAL by DuckDB).
    aware_utc = datetime(2026, 7, 20, 18, 0, 0, tzinfo=timezone.utc)
    con = duckdb.connect(str(db_path))
    con.execute("""
        CREATE TABLE fills (
            fill_id VARCHAR PRIMARY KEY, quote_id VARCHAR NOT NULL,
            rfq_id VARCHAR NOT NULL, combo_market_ticker VARCHAR NOT NULL,
            game_id VARCHAR NOT NULL, side_held VARCHAR NOT NULL,
            contracts DOUBLE NOT NULL, price DOUBLE NOT NULL, fee DOUBLE NOT NULL,
            model_fair_at_quote DOUBLE, book_fair_at_quote DOUBLE,
            blended_fair_at_quote DOUBLE, fair_at_confirm DOUBLE,
            realized_pnl DOUBLE, filled_at TIMESTAMP NOT NULL,
            reconciled BOOLEAN NOT NULL DEFAULT FALSE)""")
    con.execute(
        "INSERT INTO fills VALUES ('f1','q1','r1','T','g','yes',1,0.5,0.01,"
        "NULL,NULL,NULL,NULL,NULL,?,FALSE)", [aware_utc])
    con.close()

    db.init_database()          # runs schema + migration
    db.init_database()          # idempotent — must not raise

    with db.connect(read_only=True) as con:
        col_type = con.execute(
            "SELECT data_type FROM duckdb_columns() "
            "WHERE table_name='fills' AND column_name='filled_at'").fetchone()[0]
        migrated = con.execute("SELECT filled_at FROM fills WHERE fill_id='f1'").fetchone()[0]
    assert col_type == "TIMESTAMP WITH TIME ZONE", col_type
    assert migrated == aware_utc, \
        f"migration must preserve the instant: {migrated!r} != {aware_utc!r}"

    # New writes round-trip UTC exactly.
    later = datetime(2026, 7, 21, 1, 30, 0, tzinfo=timezone.utc)
    with db.connect() as con:
        con.execute("UPDATE fills SET filled_at=? WHERE fill_id='f1'", [later])
    with db.connect(read_only=True) as con:
        rt = con.execute("SELECT filled_at FROM fills WHERE fill_id='f1'").fetchone()[0]
    assert rt == later


def test_all_maker_timestamp_columns_are_timestamptz(monkeypatch, tmp_path):
    import importlib
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "fresh_tz.duckdb")
    importlib.reload(db)
    db.init_database()
    with db.connect(read_only=True) as con:
        naive_cols = con.execute(
            "SELECT table_name, column_name FROM duckdb_columns() "
            "WHERE data_type = 'TIMESTAMP'").fetchall()
    assert naive_cols == [], f"naive TIMESTAMP columns remain: {naive_cols}"
