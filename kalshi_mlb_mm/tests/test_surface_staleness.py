"""Issue #99 — the leg surface's two freshness guards.

Guard 1, the AGE GATE: a quote may only use surface rows younger than
SURFACE_MAX_AGE_SEC, and a stale row must produce a DISTINCT, countable
decline — `surface_stale`, not #98's `surface_too_few_books` (which means no
book published the leg at all). The two point at opposite fixes, so the tests
here assert on the reason string, not merely on "no quote".

Guard 2, the CONSTITUENT FRESHNESS VETO: Kalshi's own single-leg market for
the same leg trades in real time. If it moved AFTER the surface row was built,
the book's cached number is stale by construction however young its clock
says it is. The veto must cost ZERO extra Kalshi calls, so it reads a baseline
out of the tape fed by the reads the bot already makes.

Scaffolding is imported from test_surface_routing (#98) so both tickets
exercise the same fake engine, RFQ source and leg shapes.
"""
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common import legset
from kalshi_mlb_mm.constituent_tape import ConstituentTape
from kalshi_mlb_mm.tests.test_surface_routing import (
    CROSS_LEGS, FakeEngine, NoQuoteGW, SURFACE_FAIRS, Src, _last_decision,
    _quote_priced_payload, _setup, _surface_for)


def _fresh(seconds):
    return datetime.now(timezone.utc) - timedelta(seconds=seconds)


def _reasons(db) -> set:
    """Every decline reason the tick recorded. Used instead of "the last
    decision" for the veto's control cases: a tick that gets all the way to a
    failed submit records no row at all, so "last decision" would be None."""
    with db.connect(read_only=True) as con:
        return {r[0] for r in con.execute(
            "SELECT reason FROM quote_decisions").fetchall()}


# --------------------------------------------------------------------------- #
# Guard 1 — the age gate                                                       #
# --------------------------------------------------------------------------- #

def test_stale_rows_decline_as_surface_stale_not_too_few_books(monkeypatch,
                                                               tmp_path):
    """The ticket's first acceptance criterion: stale rows produce a distinct,
    countable reason. 42s > SURFACE_MAX_AGE_SEC(30) at BOTH books."""
    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "stale1.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(42)))
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert _last_decision(db) == ("skipped", "surface_stale")
    assert eng.ensure_calls == [], \
        "a stale surface must NOT fall back to live fetches — #98's rule"


def test_absent_rows_still_decline_as_surface_too_few_books(monkeypatch,
                                                            tmp_path):
    """The control for the test above: an EMPTY surface is a coverage problem
    and keeps #98's reason. If these two ever collapsed onto one string the
    gate would be unmeasurable."""
    from kalshi_mlb_mm.leg_surface.store import LegSurface
    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(monkeypatch, tmp_path, eng, "stale2.duckdb",
                                CROSS_LEGS, LegSurface())
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert _last_decision(db) == ("skipped", "surface_too_few_books")


def test_one_stale_book_leaves_the_fresh_ones_pricing(monkeypatch, tmp_path):
    """Per-BOOK, not per-combo: a book that goes dark drops out and the rest
    still quote. This is what makes DraftKings' structural exclusion
    survivable rather than fatal."""
    from kalshi_mlb_mm.leg_surface.store import LegSurface
    stale = _surface_for(CROSS_LEGS, {"draftkings": 0.55}, built_at=_fresh(90))
    fresh = _surface_for(CROSS_LEGS, {"fanduel": 0.56, "betmgm": 0.55},
                         built_at=_fresh(5))
    merged = LegSurface()
    for (book, route), slice_ in {**stale._by_slice, **fresh._by_slice}.items():
        merged.publish(book, route, list(slice_.values()))

    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "stale3.duckdb",
                               CROSS_LEGS, merged)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert _last_decision(db) == ("dry_run_quote", None)
    payload = _quote_priced_payload(emitted)
    for game in payload["surface_games"].values():
        assert set(game["books"]) == {"fanduel", "betmgm"}
        assert set(game["excluded_by_age"]) == {"draftkings"}


def test_zero_max_age_disables_the_gate(monkeypatch, tmp_path):
    """The documented rollback: SURFACE_MAX_AGE_SEC <= 0 is pre-#99 behaviour,
    config only. A 10-minute-old row prices again."""
    import kalshi_mlb_mm.config as cfg
    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "stale4.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(600)))
    monkeypatch.setattr(cfg, "SURFACE_MAX_AGE_SEC", 0.0)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert _last_decision(db) == ("dry_run_quote", None)


def test_age_tally_counts_used_and_excluded_per_book(monkeypatch, tmp_path):
    """A book ageing out must be a NUMBER, not an absence — this tally is what
    the periodic surface_age_summary drains."""
    from kalshi_mlb_mm.leg_surface.store import LegSurface
    stale = _surface_for(CROSS_LEGS, {"draftkings": 0.55}, built_at=_fresh(90))
    fresh = _surface_for(CROSS_LEGS, {"fanduel": 0.56, "betmgm": 0.55},
                         built_at=_fresh(5))
    merged = LegSurface()
    for (book, route), slice_ in {**stale._by_slice, **fresh._by_slice}.items():
        merged.publish(book, route, list(slice_.values()))

    eng = FakeEngine(fairs=None)
    main, _db, emitted = _setup(monkeypatch, tmp_path, eng, "stale5.duckdb",
                                CROSS_LEGS, merged)
    monkeypatch.setattr(main, "_SURFACE_AGE_TALLY", {})
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    main._surface_age_summary_tick()

    summaries = [kw["payload"] for ev, kw in emitted
                 if ev == "surface_age_summary"]
    assert len(summaries) == 1
    books = summaries[-1]["books"]
    assert books["draftkings"]["excluded_by_age"] == 2   # both games
    assert books["draftkings"]["used"] == 0
    assert books["fanduel"]["used"] == 2
    assert books["fanduel"]["excluded_by_age"] == 0
    assert main._SURFACE_AGE_TALLY == {}, "the tally drains on emit"


def test_startup_warns_per_route_not_per_book(monkeypatch, caplog):
    """FanDuel runs BOTH routes on different cadences; warning on the book
    alone would implicate its 20s structure route, which clears the gate."""
    import logging
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "SURFACE_ENABLED", True)
    monkeypatch.setattr(cfg, "SURFACE_MAX_AGE_SEC", 30.0)
    with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm"):
        main._warn_structurally_excluded_books()
    warned = {line.split("leg surface: ")[1].split()[0]
              for line in caplog.text.splitlines() if "cadence" in line}
    assert "draftkings/singles" in warned
    assert "fanduel/structure" not in warned


# --------------------------------------------------------------------------- #
# Guard 2 — the constituent freshness veto                                     #
# --------------------------------------------------------------------------- #

def _tape_with(ticker, yes_bid, yes_ask, at):
    tape = ConstituentTape(180.0, 64)
    tape.record({ticker: {"yes_bid": yes_bid, "yes_ask": yes_ask}}, at)
    return tape


def _snapshot(legs, yes_bid, yes_ask):
    return {str(l["market_ticker"]): {"yes_bid": yes_bid, "yes_ask": yes_ask}
            for l in legs}


def test_veto_fires_when_kalshi_moved_after_the_row_was_built(monkeypatch,
                                                              tmp_path):
    """The ticket's second acceptance criterion. Rows are 10s old (INSIDE the
    age gate — this is the case the age gate cannot see); Kalshi was ~0.41 when
    they were built and is ~0.61 now, a 0.20 move against a 0.03 threshold."""
    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "veto1.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    snapshot = _snapshot(CROSS_LEGS, 0.60, 0.62)
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    tape = ConstituentTape(180.0, 64)
    tape.record(_snapshot(CROSS_LEGS, 0.40, 0.42), _fresh(20))
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE", tape)

    main._discovery_tick(Src(), NoQuoteGW(), dry_run=False)

    assert _last_decision(db) == ("skipped", "surface_constituent_moved")


def test_quiet_constituent_does_not_veto(monkeypatch, tmp_path):
    """The control: same setup, Kalshi unchanged. The veto must not be a
    blanket decline of every surface combo."""
    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "veto2.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    snapshot = _snapshot(CROSS_LEGS, 0.40, 0.42)
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    tape = ConstituentTape(180.0, 64)
    tape.record(_snapshot(CROSS_LEGS, 0.40, 0.42), _fresh(20))
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE", tape)

    main._discovery_tick(Src(), _AcceptGW(), dry_run=False)

    assert "surface_constituent_moved" not in _reasons(db)


def test_no_baseline_fails_open_and_is_counted(monkeypatch, tmp_path):
    """An empty tape is 'no signal', not 'no move' — the same contract that
    makes an unreadable ticker safe in singles.jumped_tickers. Failing closed
    would decline every combo on a cold process."""
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(
        monkeypatch, tmp_path, eng, "veto3.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: _snapshot(CROSS_LEGS, 0.40, 0.42))
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE",
                        ConstituentTape(180.0, 64))

    main._discovery_tick(Src(), _AcceptGW(), dry_run=False)

    assert "surface_constituent_moved" not in _reasons(db)
    checks = [kw["payload"] for ev, kw in emitted
              if ev == "surface_constituent_check"]
    assert checks, "the miss must be recorded, not silent"
    assert {l["verdict"] for l in checks[-1]["legs"]} == {"no_baseline"}


def test_veto_disabled_by_config(monkeypatch, tmp_path):
    """Rollback: the check still runs and still records, but never declines."""
    import kalshi_mlb_mm.config as cfg
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(
        monkeypatch, tmp_path, eng, "veto4.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: _snapshot(CROSS_LEGS, 0.60, 0.62))
    tape = ConstituentTape(180.0, 64)
    tape.record(_snapshot(CROSS_LEGS, 0.40, 0.42), _fresh(20))
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE", tape)
    monkeypatch.setattr(cfg, "SURFACE_CONSTITUENT_VETO_ENABLED", False)

    main._discovery_tick(Src(), _AcceptGW(), dry_run=False)

    assert "surface_constituent_moved" not in _reasons(db)
    checks = [kw["payload"] for ev, kw in emitted
              if ev == "surface_constituent_check"]
    assert {l["verdict"] for l in checks[-1]["legs"]} == {"moved"}


def test_veto_costs_zero_extra_kalshi_calls(monkeypatch, tmp_path):
    """The ticket's third acceptance criterion. #17's snapshot is the ONLY
    Kalshi read in the quote path, and the veto reuses it."""
    eng = FakeEngine(fairs=None)
    main, _db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "veto5.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    calls = []

    def _counting(legs):
        calls.append(legs)
        return _snapshot(CROSS_LEGS, 0.40, 0.42)

    monkeypatch.setattr(main, "_leg_market_prices", _counting)
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE", ConstituentTape(180.0, 64))
    main._discovery_tick(Src(), _AcceptGW(), dry_run=False)

    assert len(calls) == 1, "the veto must add no Kalshi reads of its own"


def test_quote_time_snapshot_feeds_the_tape(monkeypatch, tmp_path):
    """The tape is only affordable because it is fed by reads already paid
    for. Without this the veto would have no baseline for the NEXT RFQ."""
    eng = FakeEngine(fairs=None)
    main, _db, _emitted = _setup(
        monkeypatch, tmp_path, eng, "veto6.duckdb", CROSS_LEGS,
        _surface_for(CROSS_LEGS, built_at=_fresh(10)))
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: _snapshot(CROSS_LEGS, 0.40, 0.42))
    tape = ConstituentTape(180.0, 64)
    monkeypatch.setattr(main, "_CONSTITUENT_TAPE", tape)

    main._discovery_tick(Src(), _AcceptGW(), dry_run=False)

    assert tape.ticker_count() == len(CROSS_LEGS)


class _AcceptGW:
    """Quote gateway that accepts submissions (the veto tests need to get PAST
    the quote step, unlike #98's NoQuoteGW)."""

    def submit_quote(self, rfq_id, yes_bid, no_bid):
        return None      # a falsy quote id: recorded as a failed submit, no DB row


# --------------------------------------------------------------------------- #
# The tape itself                                                              #
# --------------------------------------------------------------------------- #

def test_tape_returns_the_newest_point_at_or_before():
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(180.0, 64)
    tape.record({"T": {"yes_bid": 0.40, "yes_ask": 0.42}},
                now - timedelta(seconds=40))
    tape.record({"T": {"yes_bid": 0.50, "yes_ask": 0.52}},
                now - timedelta(seconds=10))

    older, at = tape.price_at_or_before("T", now - timedelta(seconds=20))
    assert older == pytest.approx(0.41, abs=0.01)
    assert at == now - timedelta(seconds=40)
    newer, _ = tape.price_at_or_before("T", now)
    assert newer == pytest.approx(0.51, abs=0.01)


def test_tape_has_no_answer_before_its_first_point():
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(180.0, 64)
    tape.record({"T": {"yes_bid": 0.40, "yes_ask": 0.42}}, now)
    assert tape.price_at_or_before("T", now - timedelta(seconds=1)) is None


def test_tape_skips_degenerate_books():
    """yes_ask = 1.00 is the known divide-by-zero devig gotcha; recording it
    would let the veto compare against a number devigged_yes itself refuses."""
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(180.0, 64)
    assert tape.record({"T": {"yes_bid": 0.0, "yes_ask": 1.0}}, now) == 0
    assert tape.price_at_or_before("T", now) is None


def test_tape_prunes_outside_the_retention_window():
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(30.0, 64)
    tape.record({"T": {"yes_bid": 0.40, "yes_ask": 0.42}},
                now - timedelta(seconds=120))
    tape.record({"T": {"yes_bid": 0.50, "yes_ask": 0.52}}, now)
    assert tape.point_count() == 1


def test_tape_is_bounded_per_ticker():
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(3600.0, 4)
    for i in range(20):
        tape.record({"T": {"yes_bid": 0.40, "yes_ask": 0.42}},
                    now + timedelta(seconds=i))
    assert tape.point_count() == 4


def test_tape_drops_out_of_order_arrivals():
    """An out-of-order append would break the bisect in price_at_or_before,
    and a deque has no useful way to reorder."""
    now = datetime.now(timezone.utc)
    tape = ConstituentTape(180.0, 64)
    tape.record({"T": {"yes_bid": 0.50, "yes_ask": 0.52}}, now)
    assert tape.record({"T": {"yes_bid": 0.40, "yes_ask": 0.42}},
                       now - timedelta(seconds=5)) == 0
    assert tape.point_count() == 1


def test_leg_ticker_map_round_trips_canonical_legs():
    """CanonicalLeg carries no ticker (and #86 re-encodes F5 winners), so the
    veto rebuilds the map by parsing each raw leg on its own."""
    from kalshi_mlb_mm import main
    mapping = main._leg_ticker_map(CROSS_LEGS)
    canon = legset.parse_legs(CROSS_LEGS)
    assert set(mapping) == set(canon)
    assert set(mapping.values()) == {l["market_ticker"] for l in CROSS_LEGS}


def test_a_stale_book_elsewhere_does_not_mislabel_a_thin_game(monkeypatch,
                                                              tmp_path):
    """Precision of the reason, not just its existence. Game A loses a stale
    book and still prices; game B is thin for want of COVERAGE. The decline
    must stay `surface_too_few_books` — pointing an operator at the cadence
    when the fix is the ingest matrix would waste the whole gate."""
    from kalshi_mlb_mm.leg_surface.store import LegSurface
    from kalshi_mlb_mm.tests.test_surface_routing import (LEG_A, LEG_B_SPREAD,
                                                          _row)

    canon_a = legset.parse_legs([LEG_A])[0]
    canon_b = legset.parse_legs([LEG_B_SPREAD])[0]
    # One publish per book (publish REPLACES a (book, route) slice), so each
    # book's rows for both games have to go in together.
    plan = {"fanduel": ([(canon_a, 0.56, 5), (canon_b, 0.56, 5)]),
            "betmgm": ([(canon_a, 0.55, 5)]),          # game B: no coverage
            "draftkings": ([(canon_a, 0.55, 90)])}     # game A: stale
    merged = LegSurface()
    for book, rows in plan.items():
        merged.publish(book, "structure",
                       [_row(book, c, fair, built_at=_fresh(age))
                        for c, fair, age in rows])

    eng = FakeEngine(fairs=None)
    main, db, _emitted = _setup(monkeypatch, tmp_path, eng, "stale6.duckdb",
                                CROSS_LEGS, merged)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)

    assert _last_decision(db) == ("skipped", "surface_too_few_books")


# --------------------------------------------------------------------------- #
# The confirm last look — the fill moment, the strictest place                 #
# --------------------------------------------------------------------------- #

def _confirm_setup(monkeypatch, tmp_path, db_name, surface):
    """An OPEN, ACCEPTED quote on a CROSS-GAME combo whose legs the surface
    prices. Mirrors test_confirm_singles_veto._setup; the fresh Kalshi read is
    identical to the snapshot so #17's own veto passes and the surface gate is
    what the test actually exercises."""
    import importlib
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main
    from kalshi_mlb_mm.tests.test_surface_routing import GREF

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(main, "_ENGINE", None)
    monkeypatch.setattr(main, "_SURFACE", surface)
    monkeypatch.setattr(cfg, "SURFACE_ENABLED", True)
    snapshot = _snapshot(CROSS_LEGS, 0.40, 0.42)
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, leg_prices_json) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-x", "r-x", "COMBO-X", "game1", 0.28, 0.68, None, 0.308,
             0.308, "open", now, None, json.dumps(snapshot)])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-x", "COMBO-X", True, "game1", json.dumps(CROSS_LEGS), now,
             "quoted", ""])
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "yes",
                                                         "contracts": 1}}, None))
    return main, db


class _MustNotConfirmGW:
    def confirm(self, qid):
        raise AssertionError("a stale surface must void, never confirm")

    def cancel(self, qid):
        return True


def _confirm_outcome(db):
    with db.connect(read_only=True) as con:
        status = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-x'").fetchone()[0]
        decision = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='q-x' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()[0]
        fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
    return status, decision, fills


def test_confirm_voids_on_a_stale_surface(monkeypatch, tmp_path):
    """A surface combo has no flight to re-fetch, so row age is its ENTIRE
    freshness proof at the fill moment. Confirming a fill off a 90s row would
    leave the last look weaker than the quote gate that preceded it."""
    main, db = _confirm_setup(monkeypatch, tmp_path, "cf1.duckdb",
                              _surface_for(CROSS_LEGS, built_at=_fresh(90)))
    main._confirm_tick(_MustNotConfirmGW(), dry_run=False)

    assert _confirm_outcome(db) == ("voided", "voided_surface_stale", 0)


def test_confirm_proceeds_on_a_fresh_surface(monkeypatch, tmp_path):
    """The control: identical setup, rows inside the gate. Without this the
    test above would pass on a confirm path that never works at all."""
    main, db = _confirm_setup(monkeypatch, tmp_path, "cf2.duckdb",
                              _surface_for(CROSS_LEGS, built_at=_fresh(5)))
    gw = _ConfirmGW()
    main._confirm_tick(gw, dry_run=False)

    _status, decision, _fills = _confirm_outcome(db)
    assert decision != "voided_surface_stale"


class _ConfirmGW:
    def __init__(self):
        self.confirmed = []

    def confirm(self, qid):
        self.confirmed.append(qid)
        return True

    def cancel(self, qid):
        return True
