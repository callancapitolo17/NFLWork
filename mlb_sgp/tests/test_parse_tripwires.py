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
