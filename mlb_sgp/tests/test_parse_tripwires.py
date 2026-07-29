"""Sweep-path parse tripwires + per-fetch counters, all 6 books (issue #35).

The regression these lock down: a book's transport is healthy, its events
still match our slate, and our PARSER understands nothing it sends back. That
produced zero FanDuel alt-total rows for weeks with no error anywhere. After
this ticket the same shape emits a WARNING naming the book and the stage.
"""
from __future__ import annotations

import logging

import pytest

from mlb_sgp.tests.sweep_drives import DRIVES

BOOKS = tuple(DRIVES)


def _records(caplog, book):
    return [r for r in caplog.records
            if r.name == f"mlb_sgp.{book}"]


@pytest.mark.parametrize("book", BOOKS)
def test_healthy_fetch_emits_summary_and_no_tripwire(book, caplog, monkeypatch):
    with caplog.at_level(logging.DEBUG, logger=f"mlb_sgp.{book}"):
        rows = DRIVES[book](monkeypatch, healthy=True)

    assert rows, f"{book}: healthy drive produced no rows — fixture is wrong"
    summaries = [r for r in _records(caplog, book)
                 if r.getMessage().startswith("sgp fetch")]
    assert len(summaries) == 1, f"{book}: expected exactly one summary line"
    msg = summaries[0].getMessage()
    assert f"book={book}" in msg and "path=sweep" in msg
    assert "legs_resolved=1" in msg
    assert "prices_returned=" in msg
    assert not any("PARSE TRIPWIRE" in r.getMessage()
                   for r in _records(caplog, book)), f"{book}: false tripwire"


@pytest.mark.parametrize("book", BOOKS)
def test_zero_parse_on_nonzero_input_trips_the_wire(book, caplog, monkeypatch):
    """The canonical silent-death shape: events matched, nothing resolved."""
    with caplog.at_level(logging.DEBUG, logger=f"mlb_sgp.{book}"):
        rows = DRIVES[book](monkeypatch, healthy=False)

    assert rows == [], f"{book}: broken-parse drive should produce no rows"
    tripwires = [r for r in _records(caplog, book)
                 if "PARSE TRIPWIRE" in r.getMessage()]
    assert len(tripwires) == 1, f"{book}: expected exactly one tripwire line"
    assert tripwires[0].levelno == logging.WARNING
    msg = tripwires[0].getMessage()
    assert f"book={book}" in msg
    assert "legs_unresolved" in msg
    assert "legs_attempted=1" in msg and "legs_resolved=0" in msg


@pytest.mark.parametrize("book", BOOKS)
def test_summary_counts_events_matched(book, caplog, monkeypatch):
    with caplog.at_level(logging.DEBUG, logger=f"mlb_sgp.{book}"):
        DRIVES[book](monkeypatch, healthy=True)
    summary = next(r.getMessage() for r in _records(caplog, book)
                   if r.getMessage().startswith("sgp fetch"))
    assert "events_matched=1" in summary


@pytest.mark.parametrize("book", BOOKS)
def test_no_book_prints_to_stdout(book, capsys, monkeypatch):
    """#36: an in-process bot run must produce no stdout from the library."""
    DRIVES[book](monkeypatch, healthy=True)
    captured = capsys.readouterr()
    assert captured.out == "", f"{book} printed to stdout: {captured.out!r}"


def test_empty_slate_stays_quiet(caplog, monkeypatch):
    """An off-day (book lists no events) is NOT a parse failure."""
    from unittest.mock import MagicMock
    from mlb_sgp import caesars
    from mlb_sgp.tests.sweep_drives import target

    client = MagicMock()
    client.list_events.return_value = []
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        assert caesars.price_sgps([target()], periods=("FG",),
                                  client=client) == []
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)
    assert any("no events listed" in r.getMessage() for r in caplog.records)


def test_transport_error_is_counted_not_blamed_on_the_parser(caplog):
    """A dead book raises BookTransportError out of price_sgps; the summary
    still lands, and it must NOT be dressed up as a parse tripwire."""
    from unittest.mock import MagicMock
    from mlb_sgp import caesars
    from mlb_sgp._shared import BookTransportError
    from mlb_sgp.tests.sweep_drives import target

    client = MagicMock()
    client.list_events.side_effect = BookTransportError("caesars", "events",
                                                        status_code=403)
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        with pytest.raises(BookTransportError):
            caesars.price_sgps([target()], periods=("FG",), client=client)

    summaries = [r for r in caplog.records
                 if r.getMessage().startswith("sgp fetch")]
    assert len(summaries) == 1
    assert "transport_errors=1" in summaries[0].getMessage()
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)


def test_sanity_drops_are_counted_and_logged(caplog, monkeypatch):
    """A combo priced absurdly above its naive product is dropped (existing
    behavior) — now it is also COUNTED, so a mispricing storm is visible."""
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        # naive = 1.9 * 1.9 = 3.61; 99.0 is far past SANITY_MULT_RATIO * naive.
        rows = DRIVES["caesars"](monkeypatch, healthy=True,
                                 combo_decimal=99.0, leg_decimal=1.9)

    assert rows == []
    assert any("SANITY-DROP" in r.getMessage() and r.levelno == logging.DEBUG
               for r in caplog.records)
    summary = next(r.getMessage() for r in caplog.records
                   if r.getMessage().startswith("sgp fetch"))
    assert "sanity_drops=4" in summary
    # Every priced combo was dropped -> the price-side tripwire fires.
    assert any("prices_empty" in r.getMessage() for r in caplog.records)


def test_counters_do_not_double_count_retry_attempts(monkeypatch,
                                                     recorded_sleeps, caplog):
    """#34 retries a wire call up to N times; one LOGICAL fetch must still
    tick the counters once."""
    from unittest.mock import MagicMock
    from mlb_sgp import caesars
    from mlb_sgp.tests.sweep_drives import drive_caesars

    calls = {"n": 0}

    def flaky_price(legs):
        calls["n"] += 1
        return {"decimal": 2.5, "american": 150}

    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        drive_caesars(monkeypatch, healthy=True)
    summary = next(r.getMessage() for r in caplog.records
                   if r.getMessage().startswith("sgp fetch"))
    # 4 spread×total combos priced from ONE target -> 4 rows, 1 target.
    assert "targets_attempted=1" in summary
    assert "prices_returned=4" in summary


# ------------------------------------------------------------------ #
# A DEAD PRICE ENDPOINT MUST NOT READ AS A PARSE FAILURE              #
# ------------------------------------------------------------------ #
# Adversarial-review finding (2026-07-29): DK/FD/NV/PX price through
# module-level helpers whose request_with_retry(stage="price") RAISES
# BookTransportError, and that exception lands in the orchestrators' broad
# excepts. Filing it under parse_failures pointed a fixer at the parser while
# the real story was a 403 — and left the prices_empty tripwire armed on a
# thin slate, below PriceCallTally's 8-call verdict floor.

def test_dead_price_endpoint_counts_as_transport_not_parse(caplog, monkeypatch):
    from mlb_sgp._shared import BookTransportError
    from mlb_sgp.tests.sweep_drives import drive_fanduel

    dead = BookTransportError("fanduel", "price", status_code=403,
                              detail="failed after 2 attempt(s)")
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.fanduel"):
        rows = drive_fanduel(monkeypatch, healthy=True, price_exc=dead)

    assert rows == []
    summary = next(r.getMessage() for r in _records(caplog, "fanduel")
                   if r.getMessage().startswith("sgp fetch"))
    assert "parse_failures=0" in summary, "a 403 was blamed on the parser"
    assert "transport_errors=4" in summary
    # #33 owns this failure and is already loud about it; a PARSE tripwire on
    # top would send the fixer to the wrong file.
    assert not any("PARSE TRIPWIRE" in r.getMessage()
                   for r in _records(caplog, "fanduel"))


def test_client_swallowed_transport_still_leaves_the_wire_armed(caplog,
                                                                monkeypatch):
    """CZR/MGM clients swallow a price-stage BookTransportError and return
    None, so the orchestrator never sees an exception and has nothing to
    reclassify. prices_empty stays armed — correct, because from the
    orchestrator's view the book simply priced nothing, and no exception was
    misfiled as a parse bug on the way."""
    from mlb_sgp.tests.sweep_drives import drive_caesars

    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        rows = drive_caesars(monkeypatch, healthy=True, combo_decimal=None)

    assert rows == []
    summary = next(r.getMessage() for r in _records(caplog, "caesars")
                   if r.getMessage().startswith("sgp fetch"))
    # The client returned None cleanly: no exception, so nothing is counted
    # as either parse or transport.
    assert "parse_failures=0" in summary and "transport_errors=0" in summary
    assert "prices_returned=0" in summary
    assert any("prices_empty" in r.getMessage()
               for r in _records(caplog, "caesars"))


def test_malformed_price_payload_is_a_visible_parse_failure(caplog,
                                                            monkeypatch):
    """A price body with decimal=None used to blow up inside the sanity check
    and vanish into `except Exception: continue`. It is a genuine parse
    problem and must now be counted and named."""
    from mlb_sgp.tests.sweep_drives import drive_caesars

    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        rows = drive_caesars(monkeypatch, healthy=True,
                             price_exc=TypeError("decimal is None"))

    assert rows == []
    summary = next(r.getMessage() for r in _records(caplog, "caesars")
                   if r.getMessage().startswith("sgp fetch"))
    assert "parse_failures=4" in summary
    assert "stages=price_combo:4" in summary


def test_transport_exception_from_a_client_method_is_also_transport(
        caplog, monkeypatch):
    """If a client DOES let a price-stage BookTransportError escape (CZR's
    contract could change), the orchestrator must still bucket it as
    transport."""
    from mlb_sgp._shared import BookTransportError
    from mlb_sgp.tests.sweep_drives import drive_caesars

    dead = BookTransportError("caesars", "price", status_code=403)
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        rows = drive_caesars(monkeypatch, healthy=True, price_exc=dead)

    assert rows == []
    summary = next(r.getMessage() for r in _records(caplog, "caesars")
                   if r.getMessage().startswith("sgp fetch"))
    assert "parse_failures=0" in summary
    assert "transport_errors=" in summary and "transport_errors=0" not in summary
    assert not any("PARSE TRIPWIRE" in r.getMessage()
                   for r in _records(caplog, "caesars"))


def test_a_genuine_parser_bug_is_still_a_parse_failure(caplog, monkeypatch):
    """The reclassification must not swallow the case it exists to expose:
    a non-transport exception stays in parse_failures and still trips."""
    from mlb_sgp.tests.sweep_drives import drive_fanduel

    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.fanduel"):
        rows = drive_fanduel(monkeypatch, healthy=True,
                             price_exc=KeyError("selectionId"))

    assert rows == []
    summary = next(r.getMessage() for r in _records(caplog, "fanduel")
                   if r.getMessage().startswith("sgp fetch"))
    assert "parse_failures=4" in summary
    assert "transport_errors=0" in summary
    assert any("PARSE TRIPWIRE" in r.getMessage()
               for r in _records(caplog, "fanduel"))


def test_thick_slate_all_price_calls_failing_is_a_transport_verdict(caplog,
                                                                    monkeypatch):
    """Past PriceCallTally's 8-call floor, an all-failed price cycle raises
    BookTransportError out of price_sgps (#33). The summary must still land,
    tallied as transport, with no parse tripwire on top.

    Uses a stub client that owns a REAL PriceCallTally: ``price_tally_for``
    deliberately hands MagicMock clients a throwaway tally (#34), so a mock
    can never reach the verdict at CZR/MGM, where the client — not the
    orchestrator — does the recording.
    """
    from types import SimpleNamespace
    from mlb_sgp import caesars
    from mlb_sgp._shared import BookTransportError, PriceCallTally
    from mlb_sgp.tests.sweep_drives import (SPREAD, TOTAL, _czr_leg, target)
    from dataclasses import replace

    event = SimpleNamespace(event_id="e1", home_team="New York Yankees",
                            away_team="Boston Red Sox", start_time=None)
    monkeypatch.setattr(caesars, "_match_events",
                        lambda events, targets: {"g1": event})
    totals = {TOTAL + i: {"over": _czr_leg(1.9), "under": _czr_leg(1.9)}
              for i in range(3)}
    parsed = {"FG": {"spreads": {SPREAD: {"home": _czr_leg(1.9),
                                          "away": _czr_leg(1.9)}},
                     "totals": totals, "moneyline": None},
              "F5": {"spreads": {}, "totals": {}, "moneyline": None}}
    monkeypatch.setattr(caesars, "parse_markets", lambda ev: parsed)

    class DeadPriceClient:
        """Every price call fails the way a real client reports it: record
        the miss on the tally, return None."""
        def __init__(self):
            self.price_calls = PriceCallTally("caesars")

        def list_events(self):
            return [event]

        def fetch_event(self, event_id):
            return {"id": "e1"}

        def price_combo(self, legs):
            self.price_calls.record(False)
            return None

    targets = [replace(target(), total=TOTAL + i) for i in range(3)]
    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        with pytest.raises(BookTransportError) as excinfo:
            caesars.price_sgps(targets, periods=("FG",),
                               client=DeadPriceClient())

    assert excinfo.value.stage == "price"      # 12 calls, 0 successes
    summary = next(r.getMessage() for r in _records(caplog, "caesars")
                   if r.getMessage().startswith("sgp fetch"))
    assert "transport_errors=1" in summary
    assert not any("PARSE TRIPWIRE" in r.getMessage()
                   for r in _records(caplog, "caesars"))
