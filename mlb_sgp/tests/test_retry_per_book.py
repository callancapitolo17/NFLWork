"""Issue #34 per-book regression: a transient blip must not cost a fetch.

Companion to test_retry_backoff.py, which pins the shared helper in
isolation. This file drives the SIX REAL CLIENTS through a fake session that
fails once and then succeeds, proving the helper is actually wired into each
book's wire calls — and, in the BetMGM case, that a blip on ONE event no
longer aborts the whole book's cycle.

Why that last case is the acceptance bar: issue #33 made transport failures
loud, and in doing so WIDENED their blast radius. Before #33 a 5xx on one
event's structure fetch returned [] and the cycle skipped that one game;
after #33 it raises, and the orchestrator loses every game. Retry is what
shrinks the window back down.

Every test here asserts the failure was RETRIED (call count), not merely
that the happy path works — a test that passes without the fix proves
nothing. Sleeps are recorded, not slept (conftest.py::recorded_sleeps).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

# The clients lazily import their legacy scraper module by TOP-LEVEL name,
# which only resolves with mlb_sgp/ itself on sys.path.
_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp._shared import BookTransportError, TargetLine  # noqa: E402


class FakeResponse:
    def __init__(self, status_code=200, body=None, text=None):
        self.status_code = status_code
        self._body = body if body is not None else {}
        self.text = text if text is not None else json.dumps(body or {})
        self.headers = {}

    def json(self):
        return self._body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise OSError(f"HTTP {self.status_code}")


SERVER_ERROR = FakeResponse(500, text="upstream hiccup")


class FlakySession:
    """Fails the first ``n_failures`` requests, then serves ``good``.

    ``failure`` defaults to a 500; pass an exception instance to model a
    connection reset / DNS blip instead.
    """

    def __init__(self, good, n_failures=1, failure=SERVER_ERROR,
                 only_url_containing=None):
        self.good = good
        self.n_failures = n_failures
        self.failure = failure
        self.only_url_containing = only_url_containing
        self.calls = 0
        self.urls: list[str] = []

    def _respond(self, url):
        self.calls += 1
        self.urls.append(url)
        targeted = (self.only_url_containing is None
                    or self.only_url_containing in str(url))
        if targeted and self.n_failures > 0:
            self.n_failures -= 1
            if isinstance(self.failure, BaseException):
                raise self.failure
            return self.failure
        return self.good(url) if callable(self.good) else self.good

    def get(self, url, **kwargs):
        return self._respond(url)

    def post(self, url, **kwargs):
        return self._respond(url)


# --------------------------------------------------------------------------- #
# ProphetX                                                                     #
# --------------------------------------------------------------------------- #

def _px_client(session):
    from mlb_sgp.prophetx_client import ProphetXClient
    c = ProphetXClient.__new__(ProphetXClient)
    c.session = session
    c.verbose = False
    return c


def test_px_list_events_recovers_from_a_transient_500(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"data": {"tournaments": []}}))
    assert _px_client(session).list_events() == []
    assert session.calls == 2                      # 1 failure + 1 recovery
    assert len(recorded_sleeps) == 1


def test_px_fetch_event_markets_recovers_from_a_connection_reset(
        recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"data": {"markets": []}}),
                           failure=ConnectionError("reset by peer"))
    assert _px_client(session).fetch_event_markets("e1") == []
    assert session.calls == 2


def test_px_persistent_500_still_raises_after_the_retries(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {}), n_failures=99)
    with pytest.raises(BookTransportError) as exc:
        _px_client(session).list_events()
    assert session.calls == 3                      # BACKGROUND max_attempts
    assert exc.value.status_code == 500


# --------------------------------------------------------------------------- #
# Novig                                                                        #
# --------------------------------------------------------------------------- #

def _nv_client(session):
    from mlb_sgp.novig_client import NovigClient
    c = NovigClient.__new__(NovigClient)
    c.session = session
    c.verbose = False
    return c


def test_nv_list_events_recovers_from_a_transient_dns_failure(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"data": {"event": []}}),
                           failure=ConnectionError("nodename nor servname"))
    assert _nv_client(session).list_events() == []
    assert session.calls == 2


def test_nv_fetch_event_legs_recovers_from_a_transient_500(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"data": {"event": [{}]}}))
    _nv_client(session).fetch_event_legs("e1")
    assert session.calls == 2


# --------------------------------------------------------------------------- #
# DraftKings                                                                   #
# --------------------------------------------------------------------------- #

def _dk_client(session):
    from mlb_sgp.dk_client import DraftKingsClient
    c = DraftKingsClient.__new__(DraftKingsClient)
    c.session = session
    c.verbose = False
    return c


def test_dk_fetch_event_selections_recovers_from_a_transient_500(
        recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"data": {"markets": []}}))
    assert _dk_client(session).fetch_event_selections("e1") == []
    assert session.calls == 2


def test_dk_403_is_not_retried(recorded_sleeps):
    """DK's production failure mode is a hard 403 — retrying it three times
    just triples the latency of a known-dead book."""
    session = FlakySession(FakeResponse(200, {}),
                           failure=FakeResponse(403, text="Access Denied"))
    with pytest.raises(BookTransportError) as exc:
        _dk_client(session).fetch_event_selections("e1")
    assert session.calls == 1
    assert exc.value.status_code == 403
    assert recorded_sleeps == []


# --------------------------------------------------------------------------- #
# FanDuel                                                                      #
# --------------------------------------------------------------------------- #

def _fd_client(session):
    from mlb_sgp.fd_client import FanDuelClient
    c = FanDuelClient.__new__(FanDuelClient)
    c.session = session
    c.verbose = False
    return c


def test_fd_fetch_event_page_recovers_from_a_transient_500(recorded_sleeps):
    markets, runners = _fd_client(
        FlakySession(FakeResponse(200, {"attachments": {}}))).fetch_event_page("e1")
    assert markets == [] and runners == []


def test_fd_404_is_not_retried_and_skips_just_that_event(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {}),
                           failure=FakeResponse(404, text="gone"))
    assert _fd_client(session).fetch_event_page("e1") == ([], [])
    assert session.calls == 1
    assert recorded_sleeps == []


# --------------------------------------------------------------------------- #
# BetMGM                                                                       #
# --------------------------------------------------------------------------- #

def _mgm_client(session, accessid="acc123"):
    from mlb_sgp.betmgm_client import BetMGMClient
    c = BetMGMClient.__new__(BetMGMClient)
    c.state = "pa"
    c.verbose = False
    c.session = session
    c._accessid = accessid
    return c


def test_mgm_accessid_recovers_from_a_transient_500(recorded_sleeps):
    """The auth gate is the highest-leverage retry at this book: every CDS
    read is gated on it, so one blip here loses the whole cycle."""
    good = FakeResponse(200, {}, text='{"publicAccessId":"harvested"}')
    session = FlakySession(good)
    assert _mgm_client(session, accessid="").accessid() == "harvested"
    assert session.calls == 2


def test_mgm_list_events_recovers_from_a_transient_500(recorded_sleeps):
    session = FlakySession(FakeResponse(200, {"fixtures": []}))
    assert _mgm_client(session).list_events() == []
    assert session.calls == 2


# --------------------------------------------------------------------------- #
# Caesars                                                                      #
# --------------------------------------------------------------------------- #

def _czr_client(session):
    from mlb_sgp.caesars_client import CaesarsClient
    c = CaesarsClient.__new__(CaesarsClient)
    c.state = "nj"
    c.verbose = False
    c.session = session
    c._token = "tok"
    c._device = "dev"
    c._cookies = {"aws-waf-token": "tok"}
    c._minted_at = 1e18          # far future: ensure_token short-circuits
    return c


def test_czr_recovers_from_a_transient_waf_challenge(recorded_sleeps):
    """The WAF challenge is HTTP 200 with an HTML body — only body_ok can
    see it, and asking again is how this book recovers."""
    good = FakeResponse(200, {"competitions": []},
                        text='{"competitions": []}')
    session = FlakySession(
        good, failure=FakeResponse(200, text="<html>waf challenge</html>"))
    assert _czr_client(session).list_events() == []
    assert session.calls == 2


def test_czr_403_is_no_longer_retried(recorded_sleeps):
    """Behavior change from #34: the old bespoke 4x1.5s loop burned ~6s on a
    hard 403 before admitting the book was dead. Now it fails immediately."""
    session = FlakySession(FakeResponse(200, {}),
                           failure=FakeResponse(403, text="Forbidden"))
    with pytest.raises(BookTransportError) as exc:
        _czr_client(session).list_events()
    assert session.calls == 1
    assert exc.value.status_code == 403
    assert recorded_sleeps == []


# --------------------------------------------------------------------------- #
# Price calls — RETRY_LIVE on BOTH paths                                       #
# --------------------------------------------------------------------------- #
# The on-demand quote path prices through price_selection_set at every book.
# DK and FD build that POST themselves rather than going through a client
# price method, so they need their own coverage — without it the highest-
# priority path (a live RFQ quote) would silently have no retry at 2 of 6
# books.

def test_dk_price_selection_set_retries_once_then_prices(recorded_sleeps):
    from mlb_sgp.draftkings import price_selection_set
    good = FakeResponse(200, {"bets": [
        {"trueOdds": 6.42, "selectionsMapped": ["a", "b"]}]})
    client = _dk_client(FlakySession(good))
    assert price_selection_set(client, ["a", "b"]) == 6.42
    assert len(recorded_sleeps) == 1
    assert sum(recorded_sleeps) <= 1.0          # RETRY_LIVE budget


def test_dk_price_selection_set_does_not_retry_a_422(recorded_sleeps):
    """422 = non-combinable. A real verdict, not a blip — retrying it would
    add latency to a quote that was never going to price."""
    from mlb_sgp.draftkings import price_selection_set
    session = FlakySession(FakeResponse(200, {}),
                           failure=FakeResponse(422, text="not combinable"))
    client = _dk_client(session)
    assert price_selection_set(client, ["a", "b"]) is None
    assert session.calls == 1
    assert recorded_sleeps == []


def test_fd_price_selection_set_retries_once_then_prices(recorded_sleeps):
    from mlb_sgp.fanduel import price_selection_set
    good = FakeResponse(200, {"betCombinations": [{
        "isSGM": True,
        # One entry per ref — FD's leg-coverage guard rejects an (N-1)-leg
        # combination, which would misprice the combo.
        "legCombinations": [{"n": 1}, {"n": 2}],
        "winAvgOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": 4.2}}},
    }]})
    client = _fd_client(FlakySession(good))
    assert price_selection_set(client, [("m1", "s1"), ("m2", "s2")]) == 4.2
    assert len(recorded_sleeps) == 1


def test_price_retries_stay_inside_the_live_budget_when_a_book_is_broken():
    """A fully-broken price endpoint must not blow the quote budget: the
    on-demand path bails after the failure, so the cost is one LIVE retry."""
    from mlb_sgp._shared import RETRY_LIVE
    from mlb_sgp.draftkings import price_selection_set
    slept: list[float] = []
    import mlb_sgp._shared as shared
    original = shared._sleep
    shared._sleep = slept.append
    try:
        client = _dk_client(FlakySession(FakeResponse(200, {}), n_failures=99))
        assert price_selection_set(client, ["a", "b"]) is None
    finally:
        shared._sleep = original
    assert sum(slept) <= RETRY_LIVE.max_total_delay_sec


# --------------------------------------------------------------------------- #
# THE ACCEPTANCE BAR — blast radius of a single transient failure              #
# --------------------------------------------------------------------------- #

def _mgm_market_payload():
    def _om(name, opts):
        return {"name": {"value": name}, "id": abs(hash(name)) % 100000,
                "isBetBuilder": True, "options": opts}

    def _opt(name, oid, odds, attr=None):
        return {"name": {"value": name}, "id": oid, "attr": attr,
                "price": {"odds": odds}}

    return [
        _om("Run Line Spread", [
            _opt("San Francisco Giants -1.5", 13, 2.4, "-1,5"),
            _opt("Athletics +1.5", 14, 1.6, "+1,5")]),
        _om("Totals", [_opt("Over 8.5", 15, 1.9),
                       _opt("Under 8.5", 16, 1.9)]),
    ]


def _mgm_two_game_session(failing_fixture_id: str | None, n_failures=1):
    """A session serving two games' fixture trees, optionally 500-ing on one.

    ``fixtureIds={id}`` is in the structure URL, so targeting one game's
    fetch is just a substring match. Only the structure and price calls go
    over this session — the test binds ``list_events`` on the instance so
    the fetch under test (``fetch_markets``) stays the real one.
    """
    def good(url):
        if "bettingoffer/picks" in str(url):        # the price call
            return FakeResponse(200, {
                "betBuilderPricingGroups": {
                    "g": {"odds": {"odds": 3.5, "americanOdds": 250},
                          "providerId": "angstrom"}}})
        return FakeResponse(200, {"fixtures": [
            {"optionMarkets": _mgm_market_payload()}]})

    return FlakySession(
        good, n_failures=n_failures,
        only_url_containing=(f"fixtureIds={failing_fixture_id}"
                             if failing_fixture_id else None))


def _mgm_two_targets():
    ct = datetime(2026, 6, 25, 23, 0, tzinfo=timezone.utc)
    return [
        TargetLine("g1", "San Francisco Giants", "Athletics", ct, "FG", -1.5, 8.5),
        TargetLine("g2", "San Francisco Giants", "Athletics", ct, "FG", -1.5, 8.5),
    ]


def _mgm_two_game_client(monkeypatch, session):
    """A real BetMGMClient whose only stubbed step is event LISTING.

    fetch_markets — the call whose retry is under test — stays the real
    implementation and goes over ``session``.
    """
    import betmgm
    from betmgm_client import Event
    evs = {
        "g1": Event(event_id="fix1", home_team="San Francisco Giants",
                    away_team="Athletics", start_time=""),
        "g2": Event(event_id="fix2", home_team="San Francisco Giants",
                    away_team="Athletics", start_time=""),
    }
    monkeypatch.setattr(betmgm, "_match_events", lambda events, targets: evs)
    client = _mgm_client(session)
    client.list_events = lambda *a, **k: list(evs.values())
    return client


def test_one_events_transient_500_no_longer_kills_the_whole_book(
        monkeypatch, recorded_sleeps):
    """THE regression this ticket exists to fix.

    Two games. Game 2's structure fetch 500s ONCE. Pre-#34 that
    BookTransportError propagated out of price_sgps and the book lost BOTH
    games (and, post-#54, its vote on every live quote). With retry, the blip
    is absorbed and both games price.
    """
    import betmgm
    session = _mgm_two_game_session(failing_fixture_id="fix2")
    client = _mgm_two_game_client(monkeypatch, session)

    rows = betmgm.price_sgps(_mgm_two_targets(), client=client, verbose=False)

    priced_games = {r.game_id for r in rows}
    assert priced_games == {"g1", "g2"}, (
        "a single transient 500 on one event still cost the book its cycle")
    assert len(recorded_sleeps) == 1               # exactly one retry needed


def test_a_persistently_broken_event_still_fails_the_book_loudly(
        monkeypatch, recorded_sleeps):
    """The complement, so the fix can't be mistaken for "swallow 5xx".

    #33 decided a 5xx is NOT a delisted event: if it survives every retry the
    book is genuinely broken and must fail loudly rather than silently
    reporting a short slate. Only the TRANSIENT window closed.
    """
    import betmgm
    session = _mgm_two_game_session(failing_fixture_id="fix2", n_failures=99)
    client = _mgm_two_game_client(monkeypatch, session)

    with pytest.raises(BookTransportError) as exc:
        betmgm.price_sgps(_mgm_two_targets(), client=client, verbose=False)
    assert exc.value.stage == "structure"
    assert exc.value.status_code == 500
