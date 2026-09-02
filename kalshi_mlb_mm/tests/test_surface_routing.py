"""Issue #98 — the router reads the leg surface for CROSS-GAME combos.

The acceptance criteria are all about WHERE a number came from and what it
cost, so every test here asserts on the FakeEngine's `ensure_calls` counter
(outbound book traffic) rather than on the price alone:

  1. a cross-game RFQ prices with ZERO outbound book requests
  2. a same-game RFQ still takes the live path, with a populated surface
     sitting right there untouched (the mirror image of #54's "cache present
     but never consulted")
  3. a mixed combo fetches for its 2-leg game ONLY

Scaffolding follows test_live_pricing.py; leg dicts go through conftest.leg()
so the Kalshi event-ticker shape is real (#71).
"""
import importlib
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common import legset
from mlb_sgp._shared import GameRef

from kalshi_mlb_mm.leg_surface.store import LegSurface, SurfaceRow
from kalshi_mlb_mm.tests.conftest import leg

# Two DIFFERENT games. TEXLAA: away TEX, home LAA. NYYBOS: away NYY, home BOS.
GAME_A = "25JUN271905TEXLAA"
GAME_B = "25JUN272005NYYBOS"

LEG_A = leg(f"KXMLBTOTAL-{GAME_A}-9", "yes")            # Over 8.5, FG
LEG_B_SPREAD = leg(f"KXMLBSPREAD-{GAME_B}-BOS2", "yes")  # home -1.5, FG
LEG_B_TOTAL = leg(f"KXMLBTOTAL-{GAME_B}-9", "yes")       # Over 8.5, FG

CROSS_LEGS = [LEG_A, LEG_B_SPREAD]                 # 1 leg + 1 leg -> all surface
MIXED_LEGS = [LEG_A, LEG_B_SPREAD, LEG_B_TOTAL]    # 1 leg + 2-leg same-game grid
SAME_GAME_LEGS = [LEG_B_SPREAD, LEG_B_TOTAL]       # one game, 2 legs -> live only

GREF = GameRef(game_id="game1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)

# Surface fairs: two books, close enough to clear SIGMA_Z_MAX=0.07 in z-space
# (sigma_z ~ 0.018). Per-leg fair ~0.555, so a 2-game product lands ~0.308 —
# comfortably inside [MIN_FAIR_PROB, MAX_FAIR_PROB].
SURFACE_FAIRS = {"draftkings": 0.55, "fanduel": 0.56}
SURFACE_COMBO_FAIR = pytest.approx(0.555 * 0.555, abs=0.002)
# Deliberately far from the surface number, so any assertion can prove which
# source priced a quote.
LIVE_FAIRS = {"betmgm": 0.30, "novig": 0.32}


def _row(book, canon_leg, fair, *, route="structure", built_at=None):
    """One SurfaceRow for a CanonicalLeg. Raw prices are cosmetic here — the
    router only ever reads fair_prob — but they are filled in rather than
    zeroed so the row is a realistic object."""
    return SurfaceRow(
        book=book, game_id=canon_leg.game_id,
        game_start_time=datetime(2025, 6, 27, 23, 5),
        period=canon_leg.period, market_type=canon_leg.market_type,
        line=canon_leg.line, side=canon_leg.side, fair_prob=fair,
        raw_decimal=1.0 / fair, raw_decimal_opp=1.0 / (1.0 - fair),
        raw_overround=1.05, route=route,
        built_at=built_at or datetime.now(timezone.utc))


def _surface_for(leg_dicts, fairs=SURFACE_FAIRS, *, route="structure",
                 built_at=None) -> LegSurface:
    """A LegSurface holding every leg in `leg_dicts` at every book in `fairs`."""
    surface = LegSurface()
    canon = legset.parse_legs(leg_dicts)
    for book, fair in fairs.items():
        surface.publish(book, route,
                        [_row(book, c, fair, route=route, built_at=built_at)
                         for c in canon])
    return surface


class FakeEngine:
    """Duck-typed OnDemandEngine that COUNTS outbound fetch requests."""

    def __init__(self, fairs=None):
        self.fairs = fairs
        self.ensure_calls = []
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

    def landed_empty(self, h):
        return False

    def result_age_sec(self, h):
        return 1.0 if self.fairs else None

    def completed_fetch_age_sec(self, h):
        return 1.0

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return True

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None


class Src:
    def __init__(self, ticker="COMBO-X"):
        self.ticker = ticker

    def poll(self):
        return [{"id": "r-1", "market_ticker": self.ticker, "contracts": 1}]

    def get_market(self, t):
        return {}


class NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not submit in these tests")


def _setup(monkeypatch, tmp_path, engine, db_name, legs, surface):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_first_pitch_utc", lambda gl: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SURFACE", surface)
    monkeypatch.setattr(cfg, "SURFACE_ENABLED", surface is not None)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-X": (True, None, legs)})
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
# Acceptance 1 — a cross-game RFQ costs ZERO outbound book requests            #
# --------------------------------------------------------------------------- #

def test_cross_game_prices_from_surface_with_zero_book_requests(monkeypatch,
                                                                tmp_path):
    # The engine has NO fairs at all: pre-#98 this RFQ would have queued two
    # fetches and skipped on_demand_pending. It must now price in ONE tick.
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf1.duckdb",
                               CROSS_LEGS, _surface_for(CROSS_LEGS))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert eng.ensure_calls == [], \
        "a cross-game combo must generate zero outbound book requests"
    assert _last_decision(db) == ("dry_run_quote", None)
    payload = _quote_priced_payload(emitted)
    assert payload["blended_fair"] == SURFACE_COMBO_FAIR


def test_cross_game_quote_records_surface_row_ages(monkeypatch, tmp_path):
    # The trace records the age of the rows that actually backed the quote —
    # research_queries.sql query 18 reads it, and #99's age gate is set from
    # that distribution. 12s: comfortably inside SURFACE_MAX_AGE_SEC (30s), so
    # this test measures the trace, not the gate (that is test_surface_age_gate).
    built = datetime.now(timezone.utc) - timedelta(seconds=12)
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(
        monkeypatch, tmp_path, eng, "surf2.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=built))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    payload = _quote_priced_payload(emitted)
    surface_games = payload["surface_games"]
    assert len(surface_games) == 2, "both games are surface-priced"
    for game in surface_games.values():
        assert set(game["books"]) == set(SURFACE_FAIRS)
        for book in game["books"].values():
            assert book["age_sec"] == pytest.approx(12, abs=5)
        assert game["excluded_by_age"] is None, "nothing was stale here"
        assert game["oldest_used_age_sec"] == pytest.approx(12, abs=5)
    # live_games stays untouched — report.py reads it.
    assert payload["live_games"] is None


# --------------------------------------------------------------------------- #
# Acceptance 2 — a same-game RFQ still takes the live path                     #
# --------------------------------------------------------------------------- #

def test_same_game_grid_ignores_a_populated_surface(monkeypatch, tmp_path):
    # The surface holds BOTH legs of this game at quotable prices. Same-game
    # legs are correlated, so multiplying them would be wrong — the grid must
    # still fetch live, and price from the live number.
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf3.duckdb",
                               SAME_GAME_LEGS, _surface_for(SAME_GAME_LEGS))
    canon = legset.parse_legs(SAME_GAME_LEGS)
    h = legset.leg_set_hash(canon)

    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert [c[0] for c in eng.ensure_calls] == [h], \
        "a same-game grid must still fetch live"
    assert _last_decision(db) == ("skipped", "on_demand_pending")

    eng.fairs = dict(LIVE_FAIRS)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    payload = _quote_priced_payload(emitted)
    assert payload["blended_fair"] == pytest.approx(0.31, abs=0.001), \
        "same-game must price from the live fetch, never the leg surface"
    assert payload["surface_games"] is None


# --------------------------------------------------------------------------- #
# Acceptance 3 — a mixed combo splits                                          #
# --------------------------------------------------------------------------- #

def test_mixed_combo_fetches_only_the_same_game_half(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf4.duckdb",
                               MIXED_LEGS, _surface_for(MIXED_LEGS))
    canon = legset.parse_legs(MIXED_LEGS)
    by_game = legset.partition_by_game(canon)
    hash_b = legset.leg_set_hash(by_game[GAME_B])

    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert [c[0] for c in eng.ensure_calls] == [hash_b], \
        "only the 2-leg game may cost a book request"
    assert _last_decision(db) == ("skipped", "on_demand_pending")

    # Game B lands live; game A stays on the surface. The combo fair is the
    # product of the two, proving both halves contributed.
    eng.fairs = dict(LIVE_FAIRS)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    payload = _quote_priced_payload(emitted)
    assert payload["blended_fair"] == pytest.approx(0.555 * 0.31, abs=0.002)
    assert len(payload["surface_games"]) == 1
    assert len(payload["live_games"]) == 1


# --------------------------------------------------------------------------- #
# Gate semantics — same #20 gate, cached input, distinguishable declines       #
# --------------------------------------------------------------------------- #

def test_one_book_on_the_surface_declines(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(
        monkeypatch, tmp_path, eng, "surf5.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, fairs={"draftkings": 0.55}))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "surface_too_few_books")
    assert eng.ensure_calls == [], "a thin surface must not trigger a fetch"
    assert _quote_priced_payload(emitted) is None


def test_dispersed_surface_books_decline(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(
        monkeypatch, tmp_path, eng, "surf6.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, fairs={"draftkings": 0.40, "fanduel": 0.70}))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "surface_dispersion")
    assert _quote_priced_payload(emitted) is None


def test_missing_leg_on_the_surface_declines_without_falling_back(monkeypatch,
                                                                 tmp_path):
    # Only game B's leg is on the surface (e.g. an alt line no book posts for
    # game A). The missing half must DECLINE, never trigger a live fetch —
    # a fallback would reintroduce the per-RFQ traffic #98 removes.
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf7.duckdb",
                               CROSS_LEGS, _surface_for([LEG_B_SPREAD]))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "surface_too_few_books")
    assert eng.ensure_calls == []


# --------------------------------------------------------------------------- #
# Confirm last look — a surface combo has no flight to re-fetch                #
# --------------------------------------------------------------------------- #

def _seed_accepted_quote(main, db, legs, fair):
    import json
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            # The taker hits our NO bid, so we hold YES at (1 - no_bid) = 0.25
            # against a ~0.308 fair — comfortably +EV, so last_look_ok turns on
            # the re-price alone, which is this test's subject.
            ["q-s", "r-1", "COMBO-X", "game1", 0.20, 0.75,
             fair, fair, fair, "open", now, None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-1", "COMBO-X", True, "game1", json.dumps(legs), now,
             "quoted", ""])
        con.execute("UPDATE live_quotes SET leg_prices_json = "
                    "'{\"L\": {\"yes_bid\": 0.5, \"yes_ask\": 0.52}}'")


def test_confirm_reprices_a_cross_game_combo_off_the_surface(monkeypatch,
                                                             tmp_path):
    # Pre-#98 this path demanded a synchronous live re-fetch of every
    # sub-combo and voided without one. A surface combo has no flight — the
    # ingest loop keeps its rows current — so it must confirm on the surface's
    # current number instead of voiding.
    eng = FakeEngine(fairs=None)          # refetch_now would return False
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf9.duckdb",
                               CROSS_LEGS, _surface_for(CROSS_LEGS))
    _seed_accepted_quote(main, db, CROSS_LEGS, 0.555 * 0.555)
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    monkeypatch.setattr(main.auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "no",
                                                         "contracts": 1}}, None))
    confirmed = []

    class GW:
        def confirm(self, qid):
            confirmed.append(qid)
            return True

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)
    assert eng.refetch_calls == [], \
        "a surface-routed combo must not trigger a confirm-time book fetch"
    assert confirmed == ["q-s"]


# --------------------------------------------------------------------------- #
# Post-fill cooldown — "the books were re-asked" without a flight              #
# --------------------------------------------------------------------------- #

def test_post_fill_cooldown_waits_for_newer_surface_rows(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    filled_at = datetime.now(timezone.utc)
    stale = _surface_for(CROSS_LEGS,
                         built_at=filled_at - timedelta(seconds=5))
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf10.duckdb",
                               CROSS_LEGS, stale)
    by_game = legset.partition_by_game(legset.parse_legs(CROSS_LEGS))

    # Rows built BEFORE the fill are the very snapshot that got picked off.
    assert main._post_fill_live_refresh_landed(by_game, filled_at) is False

    fresh = _surface_for(CROSS_LEGS,
                         built_at=filled_at + timedelta(seconds=5))
    monkeypatch.setattr(main, "_SURFACE", fresh)
    assert main._post_fill_live_refresh_landed(by_game, filled_at) is True


def test_post_fill_cooldown_needs_min_agreeing_books_refreshed(monkeypatch,
                                                               tmp_path):
    # One refreshed book cannot re-price the combo, so one is not enough.
    eng = FakeEngine(fairs=None)
    filled_at = datetime.now(timezone.utc)
    surface = LegSurface()
    canon = legset.parse_legs(CROSS_LEGS)
    surface.publish("draftkings", "singles",
                    [_row("draftkings", c, 0.55, route="singles",
                          built_at=filled_at + timedelta(seconds=5))
                     for c in canon])
    surface.publish("fanduel", "structure",
                    [_row("fanduel", c, 0.56,
                          built_at=filled_at - timedelta(seconds=5))
                     for c in canon])
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf11.duckdb",
                               CROSS_LEGS, surface)
    by_game = legset.partition_by_game(canon)
    assert main._post_fill_live_refresh_landed(by_game, filled_at) is False


# --------------------------------------------------------------------------- #
# Router units                                                                 #
# --------------------------------------------------------------------------- #

def test_routes_to_surface_only_for_single_leg_groups():
    from kalshi_mlb_mm import router
    by_game = legset.partition_by_game(legset.parse_legs(MIXED_LEGS))
    assert router.routes_to_surface(by_game[GAME_A]) is True
    assert router.routes_to_surface(by_game[GAME_B]) is False


def test_f5_tie_leg_never_routes_to_surface():
    # classify_subcombo's F5-TIE guard runs BEFORE its n == 1 branch, so a
    # lone TIE leg stays unpriceable rather than becoming a surface lookup.
    from kalshi_mlb_mm import router
    tie = legset.parse_legs([leg(f"KXMLBF5-{GAME_A}-TIE", "yes")])
    assert router.routes_to_surface(tie) is False


def test_fanduel_two_routes_count_as_one_book(monkeypatch, tmp_path):
    # FD publishes ml/spread/I1 via structure and FG/F5 totals via singles.
    # LegSurface.book_fairs collapses them; if it did not, FD alone would
    # satisfy MIN_AGREEING_BOOKS=2 and this RFQ would quote.
    surface = LegSurface()
    canon = legset.parse_legs(CROSS_LEGS)
    surface.publish("fanduel", "structure",
                    [_row("fanduel", c, 0.55, route="structure") for c in canon])
    surface.publish("fanduel", "singles",
                    [_row("fanduel", c, 0.56, route="singles") for c in canon])
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "surf8.duckdb",
                               CROSS_LEGS, surface)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "surface_too_few_books")


def test_surface_fairs_none_reproduces_pre_98_routing():
    """The rollback oracle: without a surface lookup a single-leg group must
    route to the live engine exactly as it did before this ticket."""
    from kalshi_mlb_mm import router
    canon = legset.parse_legs(CROSS_LEGS)
    single = legset.partition_by_game(canon)[GAME_A]
    seen = []

    def od_lookup(h):
        seen.append(h)
        return dict(LIVE_FAIRS)

    cons, reason = router.subcombo_consensus(
        "game1", single, None, 2, 0.07,
        on_demand_fairs=od_lookup, live_routing=True)
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.31, abs=0.001)
    assert seen == [legset.leg_set_hash(single)]


def test_surface_route_never_consults_the_on_demand_lookup():
    from kalshi_mlb_mm import router
    canon = legset.parse_legs(CROSS_LEGS)
    single = legset.partition_by_game(canon)[GAME_A]

    def od_lookup(h):
        raise AssertionError("surface-routed groups must not read the engine")

    cons, reason = router.subcombo_consensus(
        "game1", single, None, 2, 0.07,
        on_demand_fairs=od_lookup, live_routing=True,
        surface_fairs=lambda gl: dict(SURFACE_FAIRS))
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.555, abs=0.001)


def test_combo_fair_detail_skips_game_resolution_for_surface_groups():
    """Surface rows key on CanonicalLeg.game_id, so mlb_target_lines is not
    consulted — an unresolvable game must not become a false decline."""
    from kalshi_mlb_mm import router

    def resolve(gl):
        raise AssertionError("surface-routed groups must not resolve a game_id")

    detail, reason = router.combo_fair_detail(
        CROSS_LEGS, None, resolve, 2, 0.07,
        on_demand_fairs=lambda h: None, live_routing=True,
        surface_fairs=lambda gl: dict(SURFACE_FAIRS))
    assert reason == "ok"
    assert detail.fair == SURFACE_COMBO_FAIR


def test_combo_fair_detail_still_resolves_live_games():
    """The mixed case: game B is live-routed, so it still resolves — an
    unresolved live game keeps declining exactly as before."""
    from kalshi_mlb_mm import router
    detail, reason = router.combo_fair_detail(
        MIXED_LEGS, None, lambda gl: None, 2, 0.07,
        on_demand_fairs=lambda h: dict(LIVE_FAIRS), live_routing=True,
        surface_fairs=lambda gl: dict(SURFACE_FAIRS))
    assert detail is None
    assert reason == "unresolved_game"
