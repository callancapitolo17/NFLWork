"""Regression tests for issue #33 — a book's transport failure must be
accounted as a FAILURE on both the sweep path and the on-demand path.

Contract:
  * ``SGPService.refresh``      transport error -> result ``None``, strike++,
                                client torn down for reinit on the 3rd strike.
  * ``sgp_cycle``               result ``None`` -> ``clear_source`` is NOT run,
                                the book's previously-written rows survive.
                                result ``[]``   -> rows cleared (valid empty
                                slate; unchanged correct behavior).
  * ``price_on_demand``         transport error -> ``None`` + the SAME strike
                                path, so sweep and on-demand share one
                                recovery. A clean "book won't price this"
                                (event miss) must NOT count a strike.

No HTTP anywhere: the ``runners`` / ``on_demand_hooks`` test seams inject
callables, and sgp_cycle gets a fake service.
"""
from datetime import datetime, timezone

import duckdb
import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import (BookTransportError, GameRef, PricedRow,
                             ResolvedLeg, TargetLine)

CT = datetime(2026, 7, 25, 23, 0, tzinfo=timezone.utc)
TARGETS = [TargetLine(game_id="g1", home_team="H", away_team="A",
                      commence_time=CT, period="FG", spread=-1.5, total=8.5)]


def _row(book):
    return PricedRow(game_id="g1", combo="Home Spread + Over", period="FG",
                     spread_line=-1.5, total_line=8.5, bookmaker=book,
                     source=f"{book}_direct", sgp_decimal=3.5,
                     sgp_american=250, fetch_time=CT)


def _dead_book(book="draftkings", stage="events", status=403):
    def runner(_targets):
        raise BookTransportError(book, stage, status_code=status)
    return runner


# --------------------------------------------------------------------------- #
# refresh(): transport error is a failure, not an empty slate                  #
# --------------------------------------------------------------------------- #

def test_transport_error_is_reported_as_failure_not_empty_rows():
    svc = SGPService(books=("draftkings",),
                     runners={"draftkings": _dead_book()})
    out = svc.refresh(TARGETS)
    assert out["draftkings"] is None, \
        "a 403 must surface as None (failure), never as [] (empty slate)"


def test_transport_error_increments_the_strike_counter():
    svc = SGPService(books=("novig",),
                     runners={"novig": _dead_book("novig", "events")})
    svc.refresh(TARGETS)
    assert svc._state["novig"].failures == 1


def test_three_consecutive_transport_errors_reinit_the_client():
    svc = SGPService(books=("caesars",),
                     runners={"caesars": _dead_book("caesars", "auth", None)})
    sentinel = object()
    svc._state["caesars"].client = sentinel
    svc.refresh(TARGETS)
    svc.refresh(TARGETS)
    assert svc._state["caesars"].client is sentinel      # not yet
    svc.refresh(TARGETS)
    assert svc._state["caesars"].client is None          # torn down for reinit


def test_transport_error_never_resets_the_strike_counter():
    """The whole bug: [] reset failures=0, so the reinit never fired."""
    state = {"n": 0}

    def flaky(_targets):
        state["n"] += 1
        raise BookTransportError("betmgm", "auth")

    svc = SGPService(books=("betmgm",), runners={"betmgm": flaky})
    svc.refresh(TARGETS)
    svc.refresh(TARGETS)
    assert svc._state["betmgm"].failures == 2
    assert state["n"] == 2


def test_valid_empty_slate_still_counts_as_success():
    svc = SGPService(books=("fanduel",), runners={"fanduel": lambda t: []})
    out = svc.refresh(TARGETS)
    assert out["fanduel"] == []
    assert svc._state["fanduel"].failures == 0
    assert svc._state["fanduel"].last_success is not None


def test_one_dead_book_does_not_affect_a_healthy_one():
    svc = SGPService(books=("draftkings", "prophetx"), runners={
        "draftkings": _dead_book(),
        "prophetx": lambda t: [_row("prophetx")],
    })
    out = svc.refresh(TARGETS)
    assert out["draftkings"] is None
    assert len(out["prophetx"]) == 1


# --------------------------------------------------------------------------- #
# sgp_cycle(): None preserves rows, [] clears them                            #
# --------------------------------------------------------------------------- #

class _FakeService:
    def __init__(self, results):
        self._results = results

    def refresh(self, _targets):
        return dict(self._results)


def _seed_rows(db_path, source):
    """Write one row through the real db layer so the schema matches what
    clear_source/ensure_table expect."""
    from mlb_sgp import db as sgp_db
    row = PricedRow(game_id="g1", combo="Home Spread + Over", period="FG",
                    spread_line=-1.5, total_line=8.5, bookmaker="DraftKings",
                    source=source, sgp_decimal=3.5, sgp_american=250,
                    fetch_time=CT)
    sgp_db.upsert_priced_rows([row], db_path=db_path)


def _count_rows(db_path, source):
    con = duckdb.connect(db_path, read_only=True)
    try:
        return con.execute(
            "SELECT count(*) FROM mlb_sgp_odds WHERE source = ?",
            [source]).fetchone()[0]
    finally:
        con.close()


def _run_cycle(monkeypatch, tmp_path, results):
    from kalshi_common import sgp_runner
    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets",
                        lambda both_teams=False: TARGETS)
    monkeypatch.setattr(sgp_runner, "write_target_lines",
                        lambda targets, db_path: None)
    db = str(tmp_path / "market.duckdb")
    _seed_rows(db, "draftkings_direct")
    counts = sgp_runner.sgp_cycle(bot_market_db=db,
                                  service=_FakeService(results))
    return db, counts


def test_sgp_cycle_preserves_rows_when_a_book_reports_transport_failure(
        monkeypatch, tmp_path):
    db, counts = _run_cycle(monkeypatch, tmp_path, {"draftkings": None})
    assert counts["draftkings"] == -1
    assert _count_rows(db, "draftkings_direct") == 1, \
        "a dead book must keep its previous rows, not have them deleted"


def test_sgp_cycle_clears_rows_on_a_genuinely_empty_slate(monkeypatch, tmp_path):
    db, counts = _run_cycle(monkeypatch, tmp_path, {"draftkings": []})
    assert counts["draftkings"] == 0
    assert _count_rows(db, "draftkings_direct") == 0, \
        "a valid empty slate must still clear stale rows (unchanged behavior)"


# --------------------------------------------------------------------------- #
# price_on_demand(): same strike path                                         #
# --------------------------------------------------------------------------- #

GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)
LEGS = [CanonicalLeg("KXMLBGAME-26JUL25NYYBOS", "spread", -1.5, "home"),
        CanonicalLeg("KXMLBGAME-26JUL25NYYBOS", "total", 8.5, "over")]


def _od_hooks(match_event=None, build_structure=None, price=None):
    return {
        "match_event": match_event or (lambda client, game: "EV"),
        "build_structure": build_structure or (lambda c, e, g: {"s": 1}),
        "resolve": lambda structure, legs, home, away: [
            ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(len(legs))],
        "price": price or (lambda client, refs, event: 3.0),
    }


def test_on_demand_transport_error_counts_a_strike():
    def boom(client, game):
        raise BookTransportError("draftkings", "events", status_code=403)

    svc = SGPService(books=("draftkings",),
                     on_demand_hooks={"draftkings": _od_hooks(match_event=boom)})
    assert svc.price_on_demand("draftkings", GAME, LEGS) is None
    assert svc._state["draftkings"].failures == 1


def test_on_demand_transport_error_reinits_client_after_three_strikes():
    def boom(client, refs, event):
        raise BookTransportError("caesars", "price", status_code=403)

    svc = SGPService(books=("caesars",),
                     on_demand_hooks={"caesars": _od_hooks(price=boom)})
    sentinel = object()
    svc._state["caesars"].client = sentinel
    svc.price_on_demand("caesars", GAME, LEGS)
    svc.price_on_demand("caesars", GAME, LEGS)
    assert svc._state["caesars"].client is sentinel
    svc.price_on_demand("caesars", GAME, LEGS)
    assert svc._state["caesars"].client is None


def test_on_demand_clean_miss_does_not_count_a_strike():
    """"This book doesn't list the game" is not a failure."""
    svc = SGPService(books=("novig",), on_demand_hooks={
        "novig": _od_hooks(match_event=lambda client, game: None)})
    assert svc.price_on_demand("novig", GAME, LEGS) is None
    assert svc._state["novig"].failures == 0


def test_on_demand_structure_miss_does_not_count_a_strike():
    svc = SGPService(books=("novig",), on_demand_hooks={
        "novig": _od_hooks(build_structure=lambda c, e, g: None)})
    assert svc.price_on_demand("novig", GAME, LEGS) is None
    assert svc._state["novig"].failures == 0
