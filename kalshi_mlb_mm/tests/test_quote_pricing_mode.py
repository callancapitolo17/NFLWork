"""Issue #54 — QUOTE_PRICING_MODE routing: cached / shadow / live.

PIN TESTS first (written against pre-#54 code, must stay green forever):
`cached` mode is byte-identical to the pre-#54 quote path — a grid RFQ is
priced from the sweep cache with ZERO on-demand engine traffic.

Regression tests for the live/shadow modes follow; they FAIL on pre-#54
code by construction (grids never touched the engine before this ticket).

Scaffolding style follows test_main_on_demand.py; leg dicts go through
conftest.leg() so the Kalshi event-ticker shape is real (#71).
"""
import importlib

import pandas as pd
import pytest

from kalshi_common import legset
from mlb_sgp._shared import GameRef

from kalshi_mlb_mm.tests.conftest import leg

# 2-leg same-game spread+total -> grid_spread_total route. LAA is the home
# team of 25JUN271905TEXLAA, so LAA2 is home -1.5 (target "Home Spread + Over").
GRID_LEGS = [
    leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
    leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes"),
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)

# Balanced 4-cell grid: every cell devigs to ~0.25, so the CACHED consensus
# fair is ~0.2375 — far from the 0.30/0.32 LIVE fairs the FakeEngine serves,
# which lets every assertion prove WHICH source priced the quote.
ST_CELLS = {"Home Spread + Over": 4.2, "Home Spread + Under": 4.2,
            "Away Spread + Over": 3.8, "Away Spread + Under": 3.8}
CACHED_FAIR = pytest.approx(0.2375, abs=0.01)
LIVE_FAIRS = {"draftkings": 0.30, "fanduel": 0.32}
LIVE_FAIR = pytest.approx(0.31, abs=0.001)


def _grid_rows(book, game_id, spread_line, total_line, family_decimals):
    return [
        {"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
         "sgp_decimal": d, "fetch_time": None,
         "spread_line": spread_line, "total_line": total_line}
        for c, d in family_decimals.items()
    ]


class FakeEngine:
    """Duck-typed OnDemandEngine: serves canned live fairs, records traffic."""

    def __init__(self, fairs=None, empty_landed=False):
        self.fairs = fairs                    # {book: fair} | None
        self.empty_landed = empty_landed      # fresh landing with zero books
        self.ensure_calls = []
        self.refetch_calls = []

    def lookup(self, h):
        return self.fairs

    def lookup_results(self, h):
        if not self.fairs:
            return None
        from mlb_sgp._shared import OnDemandBookResult
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=4, latency_sec=1.2)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0 if (self.fairs or self.empty_landed) else None

    def landed_empty(self, h):
        return self.empty_landed

    def result_age_sec(self, h):
        return 2.0 if (self.fairs or self.empty_landed) else None

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return True

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None


class Src:
    def poll(self):
        return [{"id": "r-grid", "market_ticker": "COMBO-GRID", "contracts": 1}]

    def get_market(self, t):
        return {}


class NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not submit in these tests")


def _setup(monkeypatch, tmp_path, engine, db_name, mode):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # raising=False: the attribute does not exist on pre-#54 config, and the
    # pin tests must run (and pass) against that code.
    monkeypatch.setattr(cfg, "QUOTE_PRICING_MODE", mode, raising=False)
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame(
        _grid_rows("dk", "game1", -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", "game1", -1.5, 8.5, ST_CELLS)))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-GRID": (True, None, GRID_LEGS)})
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    return main, db, emitted


def _last_decision(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()


def _quote_priced_payload(emitted):
    payloads = [kw["payload"] for ev, kw in emitted if ev == "quote_priced"]
    return payloads[-1] if payloads else None


# --------------------------------------------------------------------------- #
# PIN — cached mode is the pre-#54 quote path, byte-identical                  #
# --------------------------------------------------------------------------- #

def test_cached_mode_grid_prices_from_cache_no_engine_traffic(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=LIVE_FAIRS)     # engine HAS fairs — must be ignored
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "pin1.duckdb",
                               mode="cached")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    assert eng.ensure_calls == [], "cached mode must never touch the engine"
    assert eng.refetch_calls == []
    payload = _quote_priced_payload(emitted)
    assert payload is not None
    assert payload["blended_fair"] == CACHED_FAIR, \
        "cached mode must price from the sweep grid, not the live fairs"


def test_cached_mode_grid_quotes_even_when_engine_has_nothing(monkeypatch, tmp_path):
    # The stronger pin: with an engine that has NO results at all, cached mode
    # still quotes (pre-#54 grids never depended on the engine).
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "pin2.duckdb",
                               mode="cached")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    assert eng.ensure_calls == []
    assert _quote_priced_payload(emitted)["blended_fair"] == CACHED_FAIR


# --------------------------------------------------------------------------- #
# LIVE mode — every in-scope shape prices from a fetch initiated at RFQ time   #
# --------------------------------------------------------------------------- #

def test_live_mode_grid_pending_then_priced_from_live_fairs(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "live1.duckdb",
                               mode="live")
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)

    # Tick 1: no live result yet — the GRID sub-combo must be fed to the
    # engine (pre-#54 code prices it from cache instead) and the RFQ skipped.
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert eng.ensure_calls and eng.ensure_calls[0][0] == h
    assert eng.ensure_calls[0][1] == GREF
    assert _last_decision(db) == ("skipped", "on_demand_pending")

    # Tick 2: live fairs landed — the quote must trace to THEM, not the cache.
    eng.fairs = dict(LIVE_FAIRS)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    payload = _quote_priced_payload(emitted)
    assert payload["blended_fair"] == LIVE_FAIR, \
        "live mode must price from the live fetch, never the sweep cache"
    assert payload["pricing_mode"] == "live"


def test_live_mode_timeout_declines_never_falls_back_to_cache(monkeypatch, tmp_path):
    # A fetch landed with ZERO books (every book over budget / declined). The
    # sweep cache holds a perfectly quotable grid — live mode must still
    # DECLINE with live_fetch_timeout, never silently quote the cache.
    eng = FakeEngine(fairs=None, empty_landed=True)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "live2.duckdb",
                               mode="live")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "live_fetch_timeout")
    assert _quote_priced_payload(emitted) is None
    # An empty landing is fresh for QUOTE_FRESH_SEC — do not hammer the books
    # with an immediate re-fetch; retry resumes once the landing ages out.
    assert eng.ensure_calls == []


def test_live_mode_one_book_declines_live_too_few_books(monkeypatch, tmp_path):
    # One book answered in time; MIN_AGREEING_BOOKS=2. The cache (2 books)
    # must not top it up — live mode declines with its own reason so the
    # monitor can tell "live fetch thin" from the sweep-era too_few_books.
    eng = FakeEngine(fairs={"draftkings": 0.30})
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "live3.duckdb",
                               mode="live")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "live_too_few_books")
    assert _quote_priced_payload(emitted) is None


def test_live_mode_confirm_refetches_grid_and_voids_on_failure(monkeypatch, tmp_path):
    # Confirm last look in live mode: a GRID quote's re-price must come from a
    # synchronous live re-fetch; a failed re-fetch voids — even though the
    # sweep cache could still price the combo.
    import json
    from datetime import datetime, timezone
    eng = FakeEngine(fairs=None)                       # refetch_now -> False
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "live4.duckdb",
                               mode="live")
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-lv", "r-grid", "COMBO-GRID", "game1", 0.5, 0.43,
             0.31, 0.31, 0.31, "open", now, None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-grid", "COMBO-GRID", True, "game1", json.dumps(GRID_LEGS), now,
             "quoted", ""])
        con.execute("UPDATE live_quotes SET leg_prices_json = "
                    "'{\"L\": {\"yes_bid\": 0.5, \"yes_ask\": 0.52}}'")
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    monkeypatch.setattr(main.auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "yes",
                                                         "contracts": 1}}, None))

    class GW:
        def confirm(self, qid):
            raise AssertionError("failed live re-fetch must void, not confirm")

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)
    assert eng.refetch_calls, "live mode confirm must live-refetch grid shapes"
    jobs, budget = eng.refetch_calls[0]
    assert jobs[0][0] == h and jobs[0][1] == GREF
    assert budget == main.CONFIRM_REFETCH_BUDGET_SEC
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-lv'").fetchone()
    assert st[0] == "voided"


# --------------------------------------------------------------------------- #
# SHADOW mode — quote from cache, ALSO live-fetch and log both fairs           #
# --------------------------------------------------------------------------- #

def test_shadow_mode_quotes_from_cache_and_fetches_live(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "sh1.duckdb",
                               mode="shadow")
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    # Quoting is NOT blocked on the live fetch...
    assert _last_decision(db) == ("dry_run_quote", None)
    assert _quote_priced_payload(emitted)["blended_fair"] == CACHED_FAIR
    # ...but the fetch was fired for the comparison dataset.
    assert eng.ensure_calls and eng.ensure_calls[0][0] == h


def test_shadow_mode_logs_cached_vs_live_pair_once_per_landing(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=dict(LIVE_FAIRS))
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "sh2.duckdb",
                               mode="shadow")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    assert _quote_priced_payload(emitted)["blended_fair"] == CACHED_FAIR
    comps = [kw["payload"] for ev, kw in emitted
             if ev == "shadow_fair_comparison"]
    assert len(comps) == 1, "shadow mode must log the cached-vs-live pair"
    assert comps[0]["cached_fair"] == CACHED_FAIR
    assert comps[0]["live_fair"] == LIVE_FAIR
    # Same landing consumed again next tick -> no duplicate comparison row.
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    comps = [1 for ev, _ in emitted if ev == "shadow_fair_comparison"]
    assert len(comps) == 1


# --------------------------------------------------------------------------- #
# Router + engine + config units                                               #
# --------------------------------------------------------------------------- #

def test_router_live_routing_prices_grids_from_lookup_not_sgp_df():
    from kalshi_mlb_mm import router
    canon = legset.parse_legs(GRID_LEGS)
    sgp_df = pd.DataFrame(_grid_rows("dk", "game1", -1.5, 8.5, ST_CELLS)
                          + _grid_rows("fd", "game1", -1.5, 8.5, ST_CELLS))
    lookup = lambda h: dict(LIVE_FAIRS)
    cons, reason = router.subcombo_consensus(
        "game1", canon, sgp_df, 2, 0.07,
        on_demand_fairs=lookup, live_routing=True)
    assert reason == "ok"
    assert cons.fair == LIVE_FAIR
    # Without a lookup, live routing must fail closed — never read sgp_df.
    cons, reason = router.subcombo_consensus(
        "game1", canon, sgp_df, 2, 0.07,
        on_demand_fairs=None, live_routing=True)
    assert cons is None and reason == "unpriceable"


def test_router_live_routing_param_is_keyword_only():
    import inspect
    from kalshi_mlb_mm import router
    for fn in (router.subcombo_consensus, router.combo_fair_detail,
               router.combo_fair):
        param = inspect.signature(fn).parameters["live_routing"]
        assert param.kind is inspect.Parameter.KEYWORD_ONLY, fn.__name__
        assert param.default is False, fn.__name__


def test_engine_landed_empty_and_result_age(monkeypatch):
    from kalshi_mlb_mm.on_demand import OnDemandEngine, QUOTE_FRESH_SEC

    clock = [100.0]

    class SvcEmpty:               # every book fails -> {} lands
        books = ("dk",)
        on_demand_deadline_sec = 10.0

        def price_on_demand(self, book, game, legs):
            return None

    eng = OnDemandEngine(SvcEmpty(), now_fn=lambda: clock[0], autostart=False)
    assert eng.landed_empty("h1") is False          # nothing stored yet
    assert eng.result_age_sec("h1") is None
    eng.ensure_fetch("h1", GREF, legset.parse_legs(GRID_LEGS))
    eng._drain_once()
    assert eng.landed_empty("h1") is True           # fresh empty landing
    assert eng.result_age_sec("h1") == pytest.approx(0.0)
    assert eng.lookup("h1") is None
    clock[0] += QUOTE_FRESH_SEC + 1.0               # aged out -> re-fetchable
    assert eng.landed_empty("h1") is False
    assert eng.result_age_sec("h1") == pytest.approx(QUOTE_FRESH_SEC + 1.0)


def test_config_mode_default_shadow_and_validation():
    import kalshi_mlb_mm.config as cfg
    # Worktrees have no kalshi_mlb_mm/.env (gitignored), so the default applies.
    assert cfg.QUOTE_PRICING_MODE == "shadow"
    with pytest.raises(ValueError):
        cfg.validate_quote_pricing_mode("bogus")
    for mode in ("live", "cached", "shadow"):
        assert cfg.validate_quote_pricing_mode(mode) == mode
