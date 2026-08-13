"""On-demand per-fetch counters + parse tripwires (issue #35).

The on-demand path is the PRIORITY surface for this ticket: under epic #53
every maker quote is priced by a live at-quote fetch, so a book that quietly
stops resolving legs stops voting on quotes — strictly worse than a stale
sweep row, and until now equally invisible.

Hook-injection style mirrors test_sgp_service_on_demand.py.
"""
import logging

import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import FetchCountersSnapshot, GameRef, ResolvedLeg

EVT = "KXMLBGAME-26JUL09NYYBOS"
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)
LOGGER_NAME = "kalshi_common.sgp_service"


def _legs(n=2):
    return [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over"),
            CanonicalLeg(EVT, "ml", None, "home")][:n]


def _resolved(n=2):
    return [ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(n)]


# Four joint fairs summing to 1.0 with a uniform 1.2 overround.
CELL_FAIRS = [0.40, 0.20, 0.25, 0.15]
CELL_DECS = [1.0 / (1.2 * p) for p in CELL_FAIRS]


def _hooks(resolved, price_fn, event="EV", structure={"fake": "s"}):
    return {
        "match_event": lambda client, game: event,
        "build_structure": lambda client, ev, game: structure,
        "resolve": lambda structure, legs, home, away: resolved,
        "price": price_fn,
    }


def _svc(book, hooks):
    return SGPService(books=(book,), on_demand_hooks={book: hooks})


def _price_all(client, refs, event):
    return CELL_DECS[0]


# ------------------------------------------------------------------ #
# Counters ride on the result                                         #
# ------------------------------------------------------------------ #

def test_successful_on_demand_result_carries_counters():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        return CELL_DECS[len(calls) - 1]

    svc = _svc("novig", _hooks(_resolved(2), price))
    res = svc.price_on_demand("novig", GAME, _legs(2))

    assert res is not None
    assert isinstance(res.counters, FetchCountersSnapshot)
    assert res.counters.book == "novig"
    assert res.counters.path == "on_demand"
    assert res.counters.events_matched == 1
    assert res.counters.legs_attempted == 2
    assert res.counters.legs_resolved == 2
    # Route A prices 2^2 cells; every one came back.
    assert res.counters.targets_attempted == 4
    assert res.counters.prices_returned == 4


def test_counters_are_a_detached_snapshot():
    """The result is frozen and must not mutate after the fact."""
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: CELL_DECS[0]))
    res = svc.price_on_demand("novig", GAME, _legs(2))
    with pytest.raises(Exception):
        res.counters.prices_returned = 99


# ------------------------------------------------------------------ #
# Tripwires on the on-demand path                                     #
# ------------------------------------------------------------------ #

def test_resolve_miss_trips_the_wire(caplog):
    """The book listed the game and handed us structure, but not one leg
    resolved — the on-demand shape of the FanDuel regression."""
    svc = _svc("novig", _hooks(None, lambda c, r, e: 2.0))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None

    trips = [r for r in caplog.records if "PARSE TRIPWIRE" in r.getMessage()]
    assert len(trips) == 1
    assert trips[0].levelno == logging.WARNING
    msg = trips[0].getMessage()
    assert "book=novig" in msg and "path=on_demand" in msg
    assert "legs_unresolved" in msg
    assert "legs_attempted=2" in msg and "legs_resolved=0" in msg


def test_price_miss_trips_the_wire(caplog):
    """Legs resolved, every price call came back None, transport never
    complained — the book stopped pricing but stayed 'up'."""
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: None))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None

    trips = [r for r in caplog.records if "PARSE TRIPWIRE" in r.getMessage()]
    assert len(trips) == 1
    assert "prices_empty" in trips[0].getMessage()


def test_healthy_on_demand_fetch_has_no_tripwire(caplog):
    calls = []

    def price(client, refs, event):
        calls.append(refs)
        return CELL_DECS[len(calls) - 1]

    svc = _svc("novig", _hooks(_resolved(2), price))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is not None
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)


def test_event_miss_is_not_a_parse_tripwire(caplog):
    """A book that simply does not list the game is not broken — and with
    no events counted, nothing should fire."""
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: 2.0, event=None))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)


def test_structure_miss_is_not_a_parse_tripwire(caplog):
    """Matched the game but the structure fetch produced nothing: the legs
    were never attempted, so the legs tripwire must stay quiet (a structure
    fetch that FAILS raises BookTransportError — #33's story, not this one)."""
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: 2.0,
                               structure=None))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)


# ------------------------------------------------------------------ #
# Summary-line volume                                                 #
# ------------------------------------------------------------------ #

def test_on_demand_summary_is_debug_not_info(caplog):
    """price_on_demand runs per RFQ per book; an INFO summary each time
    would drown bot.log. The tripwire stays WARNING."""
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: CELL_DECS[0]))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        svc.price_on_demand("novig", GAME, _legs(2))
    summaries = [r for r in caplog.records
                 if r.getMessage().startswith("sgp fetch")]
    assert len(summaries) == 1
    assert summaries[0].levelno == logging.DEBUG


def test_summary_emitted_even_when_the_book_declines(caplog):
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: 2.0, event=None))
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        svc.price_on_demand("novig", GAME, _legs(2))
    assert any(r.getMessage().startswith("sgp fetch") for r in caplog.records)


def test_transport_failure_is_counted_and_summarized(caplog):
    """A BookTransportError out of a hook keeps #33's strike path AND still
    yields a summary — tallied as transport, never as a parse tripwire."""
    from mlb_sgp._shared import BookTransportError

    def boom(client, game):
        raise BookTransportError("novig", "events", status_code=503)

    hooks = _hooks(_resolved(2), lambda c, r, e: 2.0)
    hooks["match_event"] = boom
    svc = _svc("novig", hooks)
    with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None

    summaries = [r for r in caplog.records
                 if r.getMessage().startswith("sgp fetch")]
    assert len(summaries) == 1
    assert "transport_errors=1" in summaries[0].getMessage()
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)
    # #33's failure accounting must still have run.
    assert svc._state["novig"].failures == 1


# ------------------------------------------------------------------ #
# The real hooks bind counters into the book modules                  #
# ------------------------------------------------------------------ #

@pytest.mark.parametrize("book", ["draftkings", "fanduel", "prophetx",
                                  "novig", "betmgm", "caesars"])
def test_real_hooks_thread_counters_into_resolve_and_price(book, monkeypatch):
    """Not a shape test: the inner ``except`` blocks in resolve_legs /
    price_selection_set can only count a parse failure if the hook actually
    hands them the counter."""
    import inspect

    from kalshi_common.sgp_service import SGPService as S
    from mlb_sgp._shared import FetchCounters

    svc = S(books=(book,))
    counters = FetchCounters(book, "on_demand")
    monkeypatch.setattr(svc, "_ensure_client", lambda b: None)
    svc._state[book].client = object()
    hooks = svc._book_on_demand_hooks(book, counters=counters)

    seen = {}
    mod_name = {"draftkings": "draftkings", "fanduel": "fanduel",
                "prophetx": "prophetx", "novig": "novig",
                "betmgm": "betmgm", "caesars": "caesars"}[book]
    mod = __import__(f"mlb_sgp.{mod_name}", fromlist=["x"])

    def fake_resolve(*args, **kwargs):
        seen["resolve"] = kwargs.get("counters")
        return None

    def fake_price(*args, **kwargs):
        seen["price"] = kwargs.get("counters")
        return None

    monkeypatch.setattr(mod, "resolve_legs", fake_resolve)
    monkeypatch.setattr(mod, "price_selection_set", fake_price)
    hooks = svc._book_on_demand_hooks(book, counters=counters)

    hooks["resolve"]({"s": 1}, _legs(1), "home", "away")
    hooks["price"](object(), ["r0"], _FakeEvent())

    assert seen["resolve"] is counters, f"{book}: resolve lost its counters"
    assert seen["price"] is counters, f"{book}: price lost its counters"


def test_hooks_without_counters_still_work(monkeypatch):
    """``counters`` is optional on the factory — a caller that builds hooks
    the old way must not crash."""
    from kalshi_common.sgp_service import SGPService as S

    svc = S(books=("novig",))
    monkeypatch.setattr(svc, "_ensure_client", lambda b: None)
    svc._state["novig"].client = object()
    hooks = svc._book_on_demand_hooks("novig")
    assert set(hooks) == {"match_event", "build_structure", "resolve", "price"}


class _FakeEvent:
    event_id = "e1"
