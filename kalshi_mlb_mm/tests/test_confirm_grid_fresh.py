"""Issue #17: confirm-window live re-fetch must cover GRID-priced quotes.

Phase 2 gave on-demand sub-combos a synchronous live re-fetch at confirm
time, but grid-priced quotes (the bulk of the priceable universe) still
last-looked against the ~150s `_SGP_ODDS` scrape cache — so for an accept
arriving within the same scrape cycle, cur_fair == prev_fair by construction
and the drift/EV gate was a no-op exactly in the fast-pickoff window.

These tests seed a HEALTHY stale cache that prices the combo fine (so the
pre-#17 code would confirm) and then assert the confirm decision follows the
FRESH engine fetch instead. Scaffolding follows test_main_on_demand.py.
"""
import importlib
import json
from datetime import datetime, timezone

import pandas as pd

from kalshi_common import legset
from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
from mlb_sgp._shared import GameRef

EVT = "KXMLBGAME-25JUN271905TEXLAA"
GRID_LEGS = [  # 2-leg same-game spread+total -> grid_spread_total route
    {"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
     "event_ticker": EVT, "side": "yes"},
    {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
     "event_ticker": EVT, "side": "yes"},
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)

# Canonical form: spread line -1.5 home, total 8.5 over -> "Home Spread + Over".
GRID_CANON = legset.parse_legs(GRID_LEGS)
GRID_HASH = legset.leg_set_hash(GRID_CANON)

# 4 equal decimals -> probit devig gives exactly 0.25 per cell, per book.
STALE_CACHE_FAIR = 0.25


def _stale_grid_df():
    rows = []
    for book in ("draftkings", "fanduel"):
        for cell in SPREAD_TOTAL_FAMILY:
            rows.append(dict(game_id="game1", combo=cell, period="FG",
                             bookmaker=book, sgp_decimal=3.8, fetch_time=None,
                             spread_line=-1.5, total_line=8.5))
    return pd.DataFrame(rows)


class FakeEngine:
    """Serves `fairs` from lookup after refetch_now was called; fairs=None
    models a failed/late fetch (store stays stale -> lookup None)."""

    def __init__(self, fairs=None):
        self.fairs = fairs
        self.refetch_calls = []

    def lookup(self, h):
        return self.fairs

    def lookup_results(self, h):
        if not self.fairs:
            return None
        from mlb_sgp._shared import OnDemandBookResult
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=4, latency_sec=1.0)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0 if self.fairs else None

    def ensure_fetch(self, h, game, legs):
        return False

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None


def _setup(monkeypatch, tmp_path, engine, db_name, *, yes_bid=0.28, no_bid=0.68):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", _stale_grid_df())
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(main, "_ENGINE", engine)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-grid", "r-grid", "COMBO-GRID", "game1", yes_bid, no_bid,
             None, STALE_CACHE_FAIR, STALE_CACHE_FAIR, "open", now, None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-grid", "COMBO-GRID", True, "game1", json.dumps(GRID_LEGS),
             now, "quoted", ""])
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "yes",
                                                         "contracts": 1}}, None))
    return main, db


class MustNotConfirmGW:
    def confirm(self, qid):
        raise AssertionError("stale-cache grid fair must not drive a confirm")

    def cancel(self, qid):
        return True


class ConfirmGW:
    def __init__(self):
        self.confirmed = []

    def confirm(self, qid):
        self.confirmed.append(qid)
        return True

    def cancel(self, qid):
        return True


def test_confirm_grid_fresh_fair_moved_voids(monkeypatch, tmp_path):
    """THE #17 regression: fresh fetch says fair collapsed 0.25 -> 0.10 while
    the stale cache still shows 0.25. Pre-#17 code read only the cache
    (drift == 0) and confirmed; the fresh check must void on drift."""
    eng = FakeEngine(fairs={"draftkings": 0.10, "fanduel": 0.10})
    main, db = _setup(monkeypatch, tmp_path, eng, "grid_moved.duckdb")
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    assert eng.refetch_calls, "confirm must live re-fetch the grid game"
    jobs, budget = eng.refetch_calls[0]
    assert jobs[0][0] == GRID_HASH and jobs[0][1] == GREF
    assert budget == main.CONFIRM_REFETCH_BUDGET_SEC
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-grid'").fetchone()[0]
        dec = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='q-grid' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()[0]
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
    assert st == "voided"
    assert dec == "voided_last_look"
    assert n_fills == 0


def test_confirm_grid_refetch_failure_voids(monkeypatch, tmp_path):
    """Fresh fetch timed out / landed too few books -> lookup stays None ->
    void, even though the stale cache would have priced the combo fine."""
    eng = FakeEngine(fairs=None)
    main, db = _setup(monkeypatch, tmp_path, eng, "grid_timeout.duckdb")
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    assert eng.refetch_calls, "confirm must attempt the live re-fetch"
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-grid'").fetchone()[0]
        dec = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='q-grid' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()[0]
    assert st == "voided"
    assert dec == "voided_no_fresh_books"


def test_confirm_grid_engine_absent_voids(monkeypatch, tmp_path):
    """No engine at confirm time -> can't re-price live -> void (fail-safe),
    never fall back to the stale cache."""
    main, db = _setup(monkeypatch, tmp_path, None, "grid_noengine.duckdb")
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-grid'").fetchone()[0]
        dec = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='q-grid' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()[0]
    assert st == "voided"
    assert dec == "voided_no_fresh_books"


def test_confirm_grid_fresh_within_tolerance_confirms_with_fresh_fair(
        monkeypatch, tmp_path):
    """Fresh fair 0.26 (drift 0.01 <= 0.02, still +EV) -> confirm, and
    fills.fair_at_confirm stores the FRESH 0.26, not the stale-cache 0.25."""
    eng = FakeEngine(fairs={"draftkings": 0.26, "fanduel": 0.26})
    main, db = _setup(monkeypatch, tmp_path, eng, "grid_confirm.duckdb")
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    gw = ConfirmGW()
    main._confirm_tick(gw, dry_run=False)
    assert gw.confirmed == ["q-grid"]
    with db.connect(read_only=True) as con:
        fair_at_confirm, blended_at_q = con.execute(
            "SELECT fair_at_confirm, blended_fair_at_quote FROM fills "
            "WHERE quote_id='q-grid'").fetchone()
    assert abs(fair_at_confirm - 0.26) < 1e-9, (
        f"fair_at_confirm must be the FRESH fair, got {fair_at_confirm}")
    assert abs(blended_at_q - STALE_CACHE_FAIR) < 1e-9
    # Firehose: the fresh-check event carries age/latency/delta measurement.
    checks = [kw for ev, kw in emitted if ev == "confirm_fresh_check"]
    assert len(checks) == 1
    payload = checks[0]["payload"]
    assert payload["quote_age_sec"] is not None and payload["quote_age_sec"] >= 0
    assert abs(payload["fresh_fair"] - 0.26) < 1e-9
    assert abs(payload["stale_fair"] - STALE_CACHE_FAIR) < 1e-9
    assert payload["per_book"][GRID_HASH]["draftkings"]["latency_sec"] == 1.0


def test_confirm_fresh_check_emitted_on_void_too(monkeypatch, tmp_path):
    """The measurement event fires on voids as well — that's the 'how often
    the fresh check saves us' numerator."""
    eng = FakeEngine(fairs={"draftkings": 0.10, "fanduel": 0.10})
    main, db = _setup(monkeypatch, tmp_path, eng, "grid_evt.duckdb")
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    checks = [kw for ev, kw in emitted if ev == "confirm_fresh_check"]
    assert len(checks) == 1
    payload = checks[0]["payload"]
    assert abs(payload["fresh_fair"] - 0.10) < 1e-9
    assert abs(payload["stale_fair"] - STALE_CACHE_FAIR) < 1e-9


# --------------------------------------------------------------------------- #
# Router unit level: fresh_grid_fairs plumbing                                #
# --------------------------------------------------------------------------- #

def test_router_grid_uses_fresh_lookup_not_cache():
    from kalshi_mlb_mm import router
    fresh = {GRID_HASH: {"draftkings": 0.11, "fanduel": 0.11}}
    fair = router.combo_fair(GRID_LEGS, _stale_grid_df(), lambda gl: "game1",
                             min_books=2, band=0.02,
                             fresh_grid_fairs=lambda h: fresh.get(h))
    assert fair is not None and abs(fair - 0.11) < 1e-9, (
        f"grid route must price from the fresh lookup, got {fair}")


def test_router_grid_fresh_lookup_miss_fails_closed():
    from kalshi_mlb_mm import router
    fair = router.combo_fair(GRID_LEGS, _stale_grid_df(), lambda gl: "game1",
                             min_books=2, band=0.02,
                             fresh_grid_fairs=lambda h: None)
    assert fair is None, "fresh-lookup miss must NOT fall back to the cache"


def test_router_grid_default_path_unchanged():
    from kalshi_mlb_mm import router
    fair = router.combo_fair(GRID_LEGS, _stale_grid_df(), lambda gl: "game1",
                             min_books=2, band=0.02)
    assert fair is not None and abs(fair - STALE_CACHE_FAIR) < 1e-9
