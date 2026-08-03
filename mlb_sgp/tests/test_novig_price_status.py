"""Regression tests for issue #40 — Novig's price declines carried no status.

THE TICKET'S PREMISE WAS STALE. #40 was filed as "api.novig.us unresolvable —
fix endpoint or retire the book". Verified 2026-07-30: the host resolves, the
GraphQL endpoint answers 200, the anonymous parlay endpoint answers coherently,
and the book prices for real. Production history agrees — over the last three
active days Novig returned rows on 69.8% of sweeps (median 136 rows), the third
best of the six books. The DNS line in bot.log was a transient blip. Novig is
NOT retired.

THE REAL DEFECT is diagnostic, and it is the same one issue #39 fixed for DK.
Novig's parlay endpoint has two very different non-200s:

    400 {"message": "Cannot price parlay"}   the book declines THIS combo —
                                             normal, most combos do this
    403 <html> ...                           we have been RATE-LIMITED

``submit_parlay`` collapsed both into ``{}`` + ``price_calls.record(False)``
and never called ``note_status``, so an all-failed cycle raised
``BookTransportError:price`` with ``status_code=None``. ``error_class`` in
``sgp_fetch_health`` therefore could not tell "Novig won't build these combos"
from "Novig has blocked us" — the single distinction that matters when the book
goes quiet.

The 403 is real and reproducible: two full sweeps back to back drove 100/100
price calls to failure, and the same combos recovered after a cooldown. The
production sweep's ~7-minute cadence is what has been hiding it. Epic #53's
live-at-RFQ pricing is exactly the burst pattern that will trigger it, so the
status must be visible before that lands.

No network: every test injects a fake session.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# The client lazily imports its legacy scraper module by TOP-LEVEL name.
_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp._shared import BookTransportError, PriceCallTally  # noqa: E402
from mlb_sgp.tests.test_book_transport_errors import (  # noqa: E402
    BoomSession, FakeResponse, FakeSession, _client)


RATE_LIMITED = FakeResponse(
    status_code=403,
    text="<HTML><HEAD><TITLE>403 Forbidden</TITLE></HEAD></HTML>")
CANNOT_PRICE = FakeResponse(
    status_code=400,
    body={"statusCode": 400, "message": "Cannot price parlay",
          "error": "Bad Request"})
PRICED = FakeResponse(
    status_code=201,
    body=[{"price": "0.35088", "status": "Unfilled",
           "legs": [{"vendor": "DRAFTKINGS"}]}])


def _novig(response):
    from mlb_sgp.novig_client import NovigClient
    return _client(NovigClient, FakeSession(response))


def _drive_to_verdict(client, response_desc=""):
    """Run enough price calls for the tally to reach an all-failed verdict."""
    tally = client.price_calls
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        client.submit_parlay(["outcome-a", "outcome-b"])
    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    return exc.value


# --------------------------------------------------------------------------- #
# 1. A rate-limit block must be reported as price:403                          #
# --------------------------------------------------------------------------- #

def test_rate_limited_price_calls_report_403():
    err = _drive_to_verdict(_novig(RATE_LIMITED))
    assert err.stage == "price"
    assert err.status_code == 403, (
        "a 403 from Novig's parlay endpoint means we are blocked, not that the "
        "book declined the combo — issue #40")


def test_rate_limit_surfaces_as_error_class_with_status():
    """The whole point: one string in ``sgp_fetch_health`` that says which."""
    from kalshi_common.sgp_health import transport_error_class
    err = _drive_to_verdict(_novig(RATE_LIMITED))
    assert transport_error_class(err) == "BookTransportError:price:403"


# --------------------------------------------------------------------------- #
# 2. A normal per-combo decline must be reported as price:400, NOT 403         #
# --------------------------------------------------------------------------- #

def test_cannot_price_parlay_reports_400_not_403():
    """Most Novig combos legitimately cannot be priced. A cycle of pure 400s
    means "the book won't build these combos", which must never be mistaken
    for the blockade that a 403 signals."""
    from kalshi_common.sgp_health import transport_error_class
    err = _drive_to_verdict(_novig(CANNOT_PRICE))
    assert err.status_code == 400
    assert transport_error_class(err) == "BookTransportError:price:400"


def test_a_decline_cycle_does_not_inherit_an_earlier_blockade_status():
    """The tally lives on a persistent client across sweeps. Yesterday's 403
    must not be reported against today's ordinary 400s."""
    from mlb_sgp.novig_client import NovigClient

    client = _client(NovigClient, FakeSession(RATE_LIMITED))
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        client.submit_parlay(["a", "b"])

    client.session = FakeSession(CANNOT_PRICE)          # a NEW cycle
    start = client.price_calls.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        client.submit_parlay(["a", "b"])
    with pytest.raises(BookTransportError) as exc:
        client.price_calls.verdict(start)
    assert exc.value.status_code == 400


# --------------------------------------------------------------------------- #
# 3. Transport-level failures (retries exhausted) keep their status            #
# --------------------------------------------------------------------------- #

def test_retry_exhaustion_still_carries_its_status():
    """``request_with_retry`` retries a 500 and finally raises with the status.
    ``submit_parlay`` swallows that into ``{}``; the status must survive."""
    err = _drive_to_verdict(_novig(FakeResponse(status_code=500, body={},
                                                text="oops")))
    assert err.status_code == 500


def test_connection_failure_reports_no_status_rather_than_inventing_one():
    """A DNS/connection error has no HTTP status. The verdict must say so
    instead of attributing the last status it happened to see."""
    from mlb_sgp.novig_client import NovigClient
    err = _drive_to_verdict(_client(NovigClient, BoomSession()))
    assert err.status_code is None


# --------------------------------------------------------------------------- #
# 4. Success must not poison the tally                                         #
# --------------------------------------------------------------------------- #

def test_a_successful_price_records_no_decline_status():
    client = _novig(PRICED)
    tally = client.price_calls
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        assert client.submit_parlay(["a", "b"])["decimal"] == pytest.approx(
            2.8500, rel=1e-3)
    tally.verdict(start)            # must NOT raise — every call succeeded


def test_a_mixed_cycle_does_not_reach_an_all_failed_verdict():
    """One good price is enough to prove the book is reachable."""
    from mlb_sgp.novig_client import NovigClient

    client = _client(NovigClient, FakeSession(RATE_LIMITED))
    tally = client.price_calls
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        client.submit_parlay(["a", "b"])
    client.session = FakeSession(PRICED)
    client.submit_parlay(["a", "b"])
    tally.verdict(start)            # must NOT raise


# --------------------------------------------------------------------------- #
# 5. The on-demand path shares the fix                                         #
# --------------------------------------------------------------------------- #

def test_on_demand_price_path_goes_through_submit_parlay():
    """``price_selection_set`` is the live-quote price call (epic #53). It must
    be fixed by the same change, or on-demand stays blind while the sweep sees
    the status."""
    import inspect

    from mlb_sgp import novig
    src = inspect.getsource(novig.price_selection_set)
    assert "submit_parlay" in src


def test_sweep_and_on_demand_use_DIFFERENT_submit_parlay_implementations():
    """The trap this ticket nearly shipped past. There are TWO submit_parlay
    functions: the client method (on-demand) and the legacy scraper function
    (sweep). Fixing only the client leaves every sweep verdict status-less."""
    import inspect

    from mlb_sgp.novig_client import NovigClient
    import scraper_novig_sgp

    assert inspect.isfunction(scraper_novig_sgp.submit_parlay)
    assert NovigClient.submit_parlay is not scraper_novig_sgp.submit_parlay


def test_legacy_submit_parlay_reports_a_403_block():
    import scraper_novig_sgp

    seen: list[int] = []
    priced, auth_failed = scraper_novig_sgp.submit_parlay(
        FakeSession(RATE_LIMITED), ["a", "b"], on_decline=seen.append)
    assert priced is None and auth_failed is True
    assert seen == [403]


def test_legacy_submit_parlay_reports_a_400_decline():
    import scraper_novig_sgp

    seen: list[int] = []
    priced, auth_failed = scraper_novig_sgp.submit_parlay(
        FakeSession(CANNOT_PRICE), ["a", "b"], on_decline=seen.append)
    assert priced is None and auth_failed is False
    assert seen == [400]


def test_legacy_submit_parlay_reports_nothing_on_success():
    import scraper_novig_sgp

    seen: list[int] = []
    priced, _ = scraper_novig_sgp.submit_parlay(
        FakeSession(PRICED), ["a", "b"], on_decline=seen.append)
    assert priced and seen == []


def test_legacy_submit_parlay_invents_no_status_for_a_connection_error():
    import scraper_novig_sgp

    seen: list[int] = []
    priced, _ = scraper_novig_sgp.submit_parlay(
        BoomSession(), ["a", "b"], on_decline=seen.append)
    assert priced is None and seen == []


def test_legacy_submit_parlay_on_decline_is_keyword_only():
    import inspect

    import scraper_novig_sgp
    sig = inspect.signature(scraper_novig_sgp.submit_parlay)
    with pytest.raises(TypeError):
        sig.bind(object(), ["a"], False, print)
    sig.bind(object(), ["a"], on_decline=print)


def test_legacy_submit_parlay_without_the_callback_is_unchanged():
    """Every pre-#40 call site passes no callback; that path must not raise."""
    import scraper_novig_sgp

    priced, auth_failed = scraper_novig_sgp.submit_parlay(
        FakeSession(RATE_LIMITED), ["a", "b"])
    assert priced is None and auth_failed is True


def test_orchestrator_binds_the_callback_to_the_real_price_fn():
    """``accepts_on_decline`` degrades silently by design, so pin the happy
    path: a rename that breaks the binding fails here rather than quietly
    dropping the status from every future sweep verdict."""
    from mlb_sgp._shared import accepts_on_decline
    import scraper_novig_sgp
    assert accepts_on_decline(scraper_novig_sgp.submit_parlay) is True


def test_orchestrator_end_to_end_reports_the_block_status():
    """The test that would have caught the two-implementations trap on its own:
    drive the REAL ``price_sgps`` against a 403 and assert the raised verdict
    names the status. Checking ``accepts_on_decline`` alone does not prove the
    binding survives the wrap and reaches the call site."""
    import datetime as dt
    from unittest.mock import patch

    from mlb_sgp._shared import TargetLine
    from mlb_sgp import novig
    from mlb_sgp.novig_client import Event

    ct = dt.datetime(2030, 7, 31, 2, 10, tzinfo=dt.timezone.utc)
    # Enough targets that the 4-combo fan-out clears MIN_ATTEMPTS_FOR_VERDICT;
    # the tally deliberately refuses to call a book dead on a tiny sample.
    totals = (7.5, 8.5, 9.5)
    targets = [TargetLine("g1", "Los Angeles Dodgers", "Seattle Mariners",
                          ct, "FG", -1.5, tot) for tot in totals]
    assert len(targets) * 4 >= PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT
    markets = [{"type": "SPREAD", "strike": -1.5, "is_consensus": True,
                "outcomes": []}] + [
               {"type": "TOTAL", "strike": tot, "is_consensus": True,
                "outcomes": []} for tot in totals]
    legs = {"fg": {"home_spread": {"id": "H"}, "away_spread": {"id": "A"},
                   "over": {"id": "O"}, "under": {"id": "U"}}, "f5": {}}

    class _Client:
        session = FakeSession(RATE_LIMITED)

        def list_events(self):
            return [Event("ev1", "Los Angeles Dodgers", "Seattle Mariners",
                          "LAD", "SEA", ct.isoformat())]

    def _match(nv_events, parlay_lines):
        return [{"game_id": "g1", "nv_event_id": "ev1", "nv_home_sym": "LAD",
                 "nv_away_sym": "SEA", "home_team": "Los Angeles Dodgers",
                 "away_team": "Seattle Mariners", "fg_spread_line": -1.5,
                 "fg_total_line": 8.5, "f5_spread_line": None,
                 "f5_total_line": None}]

    with patch("scraper_novig_sgp.fetch_event_legs",
               lambda s, g, v=False: (legs, markets)), \
         patch("scraper_novig_sgp.match_events", _match):
        with pytest.raises(BookTransportError) as exc:
            novig.price_sgps(targets, ("FG",), client=_Client())

    assert exc.value.stage == "price"
    assert exc.value.status_code == 403, (
        "the sweep orchestrator dropped the decline status — the legacy "
        "scraper_novig_sgp.submit_parlay is a SECOND implementation and needs "
        "the same on_decline wiring as the client method")


def test_orchestrator_tolerates_a_pre_issue_40_substitute():
    """Orchestrator tests monkeypatch this seam with stand-ins that predate the
    kwarg; binding it unconditionally would TypeError on every combo and turn
    the cycle into a silent zero-row result."""
    from mlb_sgp._shared import accepts_on_decline

    def legacy_stub(session, outcome_ids, verbose=False):
        return {"decimal": 2.0}, False

    assert accepts_on_decline(legacy_stub) is False


def test_on_demand_decline_is_not_filed_as_a_parse_failure():
    """Found by adversarial review. ``price_selection_set`` used to call
    ``float(priced.get("decimal"))`` on the ``{}`` that a declined/blocked
    combo returns. ``float(None)`` raises TypeError into the broad except and
    was counted as a PARSE failure — blaming our parser for a book that
    declined or rate-limited us, which is exactly the misdiagnosis issue #35
    removed elsewhere. caesars / betmgm / prophetx all guard this; Novig was
    the odd one out."""
    from mlb_sgp._shared import FetchCounters
    from mlb_sgp import novig

    counters = FetchCounters("novig", "on_demand")
    assert novig.price_selection_set(_novig(RATE_LIMITED), ["a", "b"],
                                     counters=counters) is None
    snap = counters.snapshot()
    assert snap.parse_failures == 0, (
        "a 403 block was filed as a parse failure — the parser is fine, the "
        "endpoint blocked us")


def test_on_demand_decline_on_a_400_is_not_a_parse_failure():
    from mlb_sgp._shared import FetchCounters
    from mlb_sgp import novig

    counters = FetchCounters("novig", "on_demand")
    assert novig.price_selection_set(_novig(CANNOT_PRICE), ["a", "b"],
                                     counters=counters) is None
    assert counters.snapshot().parse_failures == 0


def test_on_demand_still_prices_a_good_combo():
    """The guard must not swallow real prices."""
    from mlb_sgp import novig
    assert novig.price_selection_set(_novig(PRICED), ["a", "b"]) == pytest.approx(
        2.8500, rel=1e-3)


def test_on_demand_call_records_the_block_status():
    from mlb_sgp import novig

    client = _novig(RATE_LIMITED)
    tally = client.price_calls
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        assert novig.price_selection_set(client, ["a", "b"]) is None
    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    assert exc.value.status_code == 403
