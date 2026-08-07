"""Issue #55 — quorum quoting: quote on 2 fresh books, stragglers refine.

The engine lands each book's result INCREMENTALLY (a partial set is servable
the moment >=1 book returned), so the 2s discovery tick quotes as soon as
MIN_AGREEING_BOOKS fresh books pass the dispersion gate — never waiting on
the slowest book. Stragglers then refine through the EXISTING tick machinery:
within QUOTE_HYSTERESIS -> hold; beyond -> replace (issue #16 orphan-safety);
dispersion bust -> the resting quote is PULLED (`quorum_dispersion_bust`).
Two-book quotes ride a noisy 2-sample sigma, so they carry one extra margin
term (QUORUM_MARGIN_ADDON) inside the #19 floor, logged as a component.

Scaffolding style follows test_live_pricing.py; leg dicts go through
conftest.leg() so the Kalshi event-ticker shape is real (#71).
"""
import importlib
import inspect
import threading
import time

import pandas as pd
import pytest

from kalshi_common import legset
from mlb_sgp._shared import GameRef, OnDemandBookResult

from kalshi_mlb_mm.tests.conftest import leg

GRID_LEGS = [
    leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
    leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes"),
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)
GREF_B = GameRef(game_id="gB", home_team="Los Angeles Dodgers",
                 away_team="San Diego Padres", commence_time=None)
SINGLE_LEG_B = [leg("KXMLBGAME-25JUN272005SDLAD-LAD", "yes")]

# Fast-book fairs chosen to pass the #20 gate as a pair (sigma_z ~ 0.040).
FAST_FAIRS = {"fast_a": 0.30, "fast_b": 0.32}


def _result(book, fair):
    return OnDemandBookResult(book=book, fair=fair, route="partition",
                              n_cells_priced=4, latency_sec=0.1)


class StaggeredService:
    """price_on_demand stub with one book blocked on an Event until the test
    releases it. `fairs[book] is None` means the book fails (returns None)."""

    def __init__(self, fairs, release, slow_books=("slow",),
                 deadline_sec=5.0):
        self.books = tuple(fairs)
        self.fairs = dict(fairs)
        self.release = release
        self.slow_books = set(slow_books)
        self.on_demand_deadline_sec = deadline_sec
        self.calls = []

    def price_on_demand(self, book, game, legs):
        self.calls.append(book)
        if book in self.slow_books:
            self.release.wait(timeout=8.0)
        fair = self.fairs.get(book)
        if fair is None:
            return None
        return _result(book, fair)


def _poll_until(predicate, timeout_sec=3.0):
    deadline = time.monotonic() + timeout_sec
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    return predicate()


# --------------------------------------------------------------------------- #
# Engine — incremental landing semantics                                       #
# --------------------------------------------------------------------------- #

def test_engine_serves_partial_results_mid_flight():
    """THE #55 engine regression: with one book still in flight, lookup must
    already serve the books that HAVE landed (on pre-quorum code the store
    only fills at completion, so this poll times out)."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    svc = StaggeredService({**FAST_FAIRS, "slow": 0.31}, release)
    eng = OnDemandEngine(svc, autostart=False)
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    eng.ensure_fetch(h, GREF, canon)
    worker = threading.Thread(target=eng._drain_once, daemon=True)
    worker.start()
    try:
        assert _poll_until(lambda: eng.lookup(h) is not None
                           and len(eng.lookup(h)) == 2), \
            "partial (2-book) result must be servable while the slow book " \
            "is still in flight"
        fairs = eng.lookup(h)
        assert fairs == FAST_FAIRS
        # Mid-flight is NOT an empty landing — live_fetch_timeout semantics
        # belong to COMPLETED zero-book flights only.
        assert eng.landed_empty(h) is False
    finally:
        release.set()
    worker.join(timeout=8.0)
    assert eng.lookup(h) == {**FAST_FAIRS, "slow": 0.31}


def test_engine_zero_book_mid_flight_stays_pending():
    """No partial landed yet -> no store entry: lookup None, landed_empty
    False (the poll keeps the RFQ on_demand_pending, exactly as today)."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    svc = StaggeredService({"slow_a": None, "slow_b": None}, release,
                           slow_books=("slow_a", "slow_b"))
    eng = OnDemandEngine(svc, autostart=False)
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    eng.ensure_fetch(h, GREF, canon)
    worker = threading.Thread(target=eng._drain_once, daemon=True)
    worker.start()
    try:
        assert _poll_until(lambda: len(svc.calls) == 2, timeout_sec=3.0)
        assert eng.lookup(h) is None
        assert eng.landed_empty(h) is False
        assert eng.landed_at(h) is None
    finally:
        release.set()
    worker.join(timeout=8.0)
    # Both books failed -> the COMPLETED flight is an empty landing.
    assert eng.lookup(h) is None
    assert eng.landed_empty(h) is True


def test_engine_late_book_discarded_after_completion():
    """A book that returns after its flight finished (dropped at the budget)
    must not mutate the landed entry — post-completion immutability keeps
    refetch_now's landed-after-entry guarantee airtight."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    svc = StaggeredService({"fast_a": 0.30, "fast_b": 0.32, "slow": 0.90},
                           release, deadline_sec=0.4)
    eng = OnDemandEngine(svc, autostart=False)
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    eng.ensure_fetch(h, GREF, canon)
    eng._drain_once()                     # synchronous; slow book dropped
    assert eng.lookup(h) == FAST_FAIRS
    landed = eng.landed_at(h)
    release.set()                         # slow book now finishes... too late
    time.sleep(0.3)
    assert eng.lookup(h) == FAST_FAIRS, "late book must be discarded"
    assert eng.landed_at(h) == landed


def test_engine_records_budget_dropped_books():
    """#55 item 4 observability: a book over its #50 budget is dropped from
    THIS flight and recorded as such; a book that merely FAILED (returned
    None inside the budget) is not 'dropped'."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    svc = StaggeredService({"fast_a": 0.30, "failer": None, "slow": 0.31},
                           release, deadline_sec=0.4)
    eng = OnDemandEngine(svc, autostart=False)
    canon = legset.parse_legs(GRID_LEGS)
    h = legset.leg_set_hash(canon)
    eng.ensure_fetch(h, GREF, canon)
    try:
        eng._drain_once()
        assert eng.dropped_books(h) == ("slow",)
        assert eng.lookup(h) == {"fast_a": 0.30}
    finally:
        release.set()


# --------------------------------------------------------------------------- #
# Tick-level harness (dry_run=False; style follows test_live_pricing._setup)   #
# --------------------------------------------------------------------------- #

class Src:
    def poll(self):
        return [{"id": "r-q", "market_ticker": "COMBO-Q", "contracts": 1}]

    def get_market(self, t):
        return {}


class RecordingGW:
    def __init__(self):
        self.submits = []
        self.cancels = []

    def submit_quote(self, rid, yb, nb):
        self.submits.append((rid, yb, nb))
        return f"qid-{len(self.submits)}"

    def cancel(self, qid):
        self.cancels.append(qid)
        return True


class MutableFakeEngine:
    """Duck-typed engine whose served fairs the test mutates between ticks
    to simulate staggered book arrivals."""

    def __init__(self, fairs):
        self.fairs = dict(fairs)
        self.ensure_calls = []

    def lookup(self, h):
        return dict(self.fairs) if self.fairs else None

    def lookup_results(self, h):
        if not self.fairs:
            return None
        return {b: _result(b, f) for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0 if self.fairs else None

    def landed_empty(self, h):
        return False

    def result_age_sec(self, h):
        return 1.0 if self.fairs else None

    def dropped_books(self, h):
        return ()

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append(h)
        return True

    def refetch_now(self, jobs, deadline_sec):
        return True


def _setup(monkeypatch, tmp_path, engine, db_name, legs=GRID_LEGS):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # One open quote's worst-case exposure saturates the default per-combo
    # cap and would mask the refine paths (same as test_quote_replace).
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 200.0)
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame(
        {"game_id": ["game1"], "combo": ["c"], "period": ["FG"],
         "bookmaker": ["dk"], "sgp_decimal": [2.0], "fetch_time": [None],
         "spread_line": [-1.5], "total_line": [8.5]}))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(
        main, "_resolve_game_for_legs",
        lambda gl: "game1" if gl[0].game_id == "25JUN271905TEXLAA" else "gB")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(
        main, "_game_ref", lambda gid: GREF if gid == "game1" else GREF_B)
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs_: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-Q": (True, None, legs)})
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


def _quote_rows(db):
    with db.connect(read_only=True) as con:
        return dict(con.execute(
            "SELECT quote_id, status FROM live_quotes").fetchall())


# --------------------------------------------------------------------------- #
# Quote at quorum — the tick submits while the slow book is still in flight    #
# --------------------------------------------------------------------------- #

def test_quote_fires_at_quorum_while_slow_book_in_flight(monkeypatch, tmp_path):
    """THE #55 tick regression: 2 fast books landed, DK-shaped straggler
    still in flight -> the very next tick QUOTES (pre-quorum code keeps the
    RFQ on_demand_pending until the whole pool completes)."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    svc = StaggeredService({**FAST_FAIRS, "slow": 0.31}, release)
    eng = OnDemandEngine(svc, autostart=False)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "quorum1.duckdb")
    gw = RecordingGW()

    main._discovery_tick(Src(), gw, dry_run=False)        # queue the fetch
    assert _last_decision(db) == ("skipped", "on_demand_pending")
    worker = threading.Thread(target=eng._drain_once, daemon=True)
    worker.start()
    try:
        canon = legset.parse_legs(GRID_LEGS)
        h = legset.leg_set_hash(canon)
        assert _poll_until(lambda: eng.lookup(h) is not None)
        main._discovery_tick(Src(), gw, dry_run=False)    # quorum tick
        assert not release.is_set() and worker.is_alive(), \
            "test invariant: the slow book must still be in flight"
        assert len(gw.submits) == 1, "quote must fire at quorum (2 books)"
        assert _last_decision(db) == ("quoted", None)
        payload = _quote_priced_payload(emitted)
        assert payload["quorum_size"] == 2
        assert payload["live_games"], \
            "a quorum quote must keep the #54 post-RFQ live trace"
    finally:
        release.set()
    worker.join(timeout=8.0)


def test_all_books_over_budget_declines_live_fetch_timeout(monkeypatch, tmp_path):
    """All books blow the budget -> the completed flight is an empty landing
    -> decline live_fetch_timeout (#54 vocabulary, unchanged by quorum)."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    release = threading.Event()
    release.set()                                   # books "run" but all fail
    svc = StaggeredService({"a": None, "b": None}, release, slow_books=())
    eng = OnDemandEngine(svc, autostart=False)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "quorum2.duckdb")
    gw = RecordingGW()
    main._discovery_tick(Src(), gw, dry_run=False)
    eng._drain_once()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert _last_decision(db) == ("skipped", "live_fetch_timeout")
    assert gw.submits == []


# --------------------------------------------------------------------------- #
# Straggler refinement — hold / replace / pull                                 #
# --------------------------------------------------------------------------- #

def test_straggler_within_hysteresis_no_action(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    eng = MutableFakeEngine({"fast_a": 0.31, "fast_b": 0.31})
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "quorum3.duckdb")
    monkeypatch.setattr(cfg, "QUORUM_MARGIN_ADDON", 0.0)  # isolate hysteresis
    gw = RecordingGW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.submits) == 1
    # Straggler agrees exactly -> consensus unchanged -> hold the resting quote.
    eng.fairs["slow"] = 0.31
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.submits) == 1, "within hysteresis the quote must rest"
    assert _quote_rows(db) == {"qid-1": "open"}
    # The hold must come from HYSTERESIS (no new decision row), not from an
    # exposure cap quietly blocking the re-quote.
    assert _last_decision(db) == ("quoted", None)
    assert _quote_priced_payload(emitted)["refine_action"] == "hold"


def test_straggler_beyond_hysteresis_replaces(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    eng = MutableFakeEngine(dict(FAST_FAIRS))          # median 0.31
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "quorum4.duckdb")
    monkeypatch.setattr(cfg, "QUORUM_MARGIN_ADDON", 0.0)  # isolate hysteresis
    gw = RecordingGW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.submits) == 1
    # Straggler moves the median 0.31 -> 0.32 (> QUOTE_HYSTERESIS=0.005).
    eng.fairs["slow"] = 0.33
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.submits) == 2, "beyond hysteresis the quote must replace"
    assert _quote_rows(db) == {"qid-1": "replaced", "qid-2": "open"}
    payload = _quote_priced_payload(emitted)
    assert payload["refine_action"] == "replace"
    assert payload["resting_quote_id"] == "qid-1"


def test_straggler_dispersion_bust_pulls_quote(monkeypatch, tmp_path):
    """THE #55 pull regression: a straggler that blows the #20 dispersion
    gate must PULL the resting quote (`quorum_dispersion_bust`) — pre-quorum
    code just skipped the re-quote and left the stale quote resting."""
    eng = MutableFakeEngine(dict(FAST_FAIRS))
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "quorum5.duckdb")
    gw = RecordingGW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.submits) == 1
    eng.fairs["slow"] = 0.55                 # sigma_z >> SIGMA_Z_MAX
    main._discovery_tick(Src(), gw, dry_run=False)
    assert gw.cancels == ["qid-1"], "dispersion bust must pull the quote"
    assert _quote_rows(db) == {"qid-1": "cancelled"}
    assert _last_decision(db) == ("pulled", "quorum_dispersion_bust")
    pulls = [kw for ev, kw in emitted if ev == "quote_pulled"]
    assert pulls and pulls[0]["payload"]["reason"] == "quorum_dispersion_bust"


# --------------------------------------------------------------------------- #
# Quorum margin add-on — one extra term inside the #19 floor                   #
# --------------------------------------------------------------------------- #

def test_quorum_addon_widens_quote_and_is_logged():
    from kalshi_mlb_mm import pricing
    base = pricing.quote(0.5, 0.03, sigma_pts=0.0, n_games=1,
                         min_margin_pts=0.01, k_sigma=1.0)
    quorum = pricing.quote(0.5, 0.03, sigma_pts=0.0, n_games=1,
                           min_margin_pts=0.01, k_sigma=1.0,
                           quorum_addon_pts=0.01)
    assert quorum.floor_pts == pytest.approx(base.floor_pts + 0.01)
    assert quorum.quorum_addon_pts == 0.01
    assert base.quorum_addon_pts == 0.0
    assert quorum.yes_bid < base.yes_bid and quorum.no_bid < base.no_bid, \
        "the add-on must WIDEN both sides"


def test_quorum_addon_only_at_exactly_two_books(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    # 2 books -> quorum quote carries the add-on as a logged margin component.
    eng2 = MutableFakeEngine({"fast_a": 0.31, "fast_b": 0.31})
    main, db, emitted = _setup(monkeypatch, tmp_path, eng2, "quorum6.duckdb")
    main._discovery_tick(Src(), RecordingGW(), dry_run=True)
    payload = _quote_priced_payload(emitted)
    assert payload["quorum_size"] == 2
    assert payload["quorum_addon_pts"] == pytest.approx(cfg.QUORUM_MARGIN_ADDON)
    assert cfg.QUORUM_MARGIN_ADDON > 0.0


def test_no_addon_at_three_books(monkeypatch, tmp_path):
    eng3 = MutableFakeEngine({"fast_a": 0.31, "fast_b": 0.31, "slow": 0.31})
    main, db, emitted = _setup(monkeypatch, tmp_path, eng3, "quorum7.duckdb")
    main._discovery_tick(Src(), RecordingGW(), dry_run=True)
    payload = _quote_priced_payload(emitted)
    assert payload["quorum_size"] == 3
    assert payload["quorum_addon_pts"] == 0.0


def test_combo_fair_detail_carries_min_n_books():
    """Cross-game combo: game1 answered with 2 books, gB with 3 -> the combo
    is a quorum (2-book) quote and min_n_books must say so."""
    from kalshi_mlb_mm import router
    canon = legset.parse_legs(GRID_LEGS + SINGLE_LEG_B)
    by_game = legset.partition_by_game(canon)
    hashes = {gl[0].game_id: legset.leg_set_hash(gl) for gl in by_game.values()}
    fairs_by_hash = {
        hashes["25JUN271905TEXLAA"]: {"fast_a": 0.30, "fast_b": 0.32},
        hashes["25JUN272005SDLAD"]: {"a": 0.50, "b": 0.51, "c": 0.52},
    }
    resolve = lambda gl: gl[0].game_id
    detail, reason = router.combo_fair_detail(
        GRID_LEGS + SINGLE_LEG_B, None, resolve, 2, 0.07,
        on_demand_fairs=lambda h: fairs_by_hash.get(h), live_routing=True)
    assert reason == "ok"
    assert detail.min_n_books == 2


# --------------------------------------------------------------------------- #
# Interface discipline — new params keyword-only, config key present           #
# --------------------------------------------------------------------------- #

def test_new_params_keyword_only_and_config_key():
    import dataclasses
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import pricing, router
    param = inspect.signature(pricing.quote).parameters["quorum_addon_pts"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default == 0.0
    assert isinstance(cfg.QUORUM_MARGIN_ADDON, float)
    assert cfg.QUORUM_MARGIN_ADDON >= 0.0
    fields = {f.name for f in dataclasses.fields(router.ComboFair)}
    assert "min_n_books" in fields
