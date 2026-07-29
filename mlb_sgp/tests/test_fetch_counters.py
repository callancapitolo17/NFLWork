"""Unit tests for the per-fetch drop counters + parse tripwires (issue #35).

The counters answer the question ``PriceCallTally`` (#33) cannot: a book's
transport can be perfectly healthy while OUR PARSER stopped understanding its
payload. The canonical casualty is the FanDuel alt-total paren regression —
weeks of zero rows with no error anywhere.
"""
from __future__ import annotations

import logging
import threading

import pytest

from mlb_sgp._shared import FetchCounters, FetchCountersSnapshot, fetch_counters


# ------------------------------------------------------------------ #
# Counter arithmetic                                                  #
# ------------------------------------------------------------------ #

def test_counters_start_at_zero():
    c = FetchCounters(book="fanduel", path="sweep")
    snap = c.snapshot()
    assert snap.book == "fanduel"
    assert snap.path == "sweep"
    assert (snap.events_seen, snap.events_matched, snap.targets_attempted,
            snap.legs_attempted, snap.legs_resolved, snap.prices_returned,
            snap.sanity_drops, snap.parse_failures,
            snap.transport_errors) == (0,) * 9
    assert snap.by_stage == ()


def test_bump_accumulates():
    c = FetchCounters(book="novig", path="on_demand")
    c.bump("legs_attempted", 3)
    c.bump("legs_resolved")
    c.bump("legs_resolved")
    snap = c.snapshot()
    assert snap.legs_attempted == 3
    assert snap.legs_resolved == 2


def test_bump_rejects_unknown_counter():
    """A typo'd counter name must fail loudly, not silently vanish."""
    c = FetchCounters(book="novig", path="sweep")
    with pytest.raises(AttributeError):
        c.bump("legs_resolvd")


def test_snapshot_is_immutable_and_detached():
    c = FetchCounters(book="caesars", path="sweep")
    c.bump("prices_returned", 4)
    snap = c.snapshot()
    c.bump("prices_returned", 4)          # keep counting after the snapshot
    assert snap.prices_returned == 4      # snapshot did not follow along
    with pytest.raises(Exception):        # frozen dataclass
        snap.prices_returned = 99


def test_bump_is_thread_safe():
    """Sweep pricing fans out over a thread pool, so += must be locked."""
    c = FetchCounters(book="draftkings", path="sweep")

    def worker():
        for _ in range(500):
            c.bump("prices_returned")

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert c.snapshot().prices_returned == 8 * 500


# ------------------------------------------------------------------ #
# Named parse failures                                                #
# ------------------------------------------------------------------ #

def test_parse_failure_records_named_stage(caplog):
    c = FetchCounters(book="prophetx", path="sweep")
    logger = logging.getLogger("test.parse_failure")
    with caplog.at_level(logging.DEBUG, logger="test.parse_failure"):
        c.record_parse_failure("price_combo", logger, exc=ValueError("boom"))
        c.record_parse_failure("price_combo", logger, exc=ValueError("boom2"))
        c.record_parse_failure("resolve_legs", logger, exc=KeyError("k"))

    snap = c.snapshot()
    assert snap.parse_failures == 3
    assert dict(snap.by_stage) == {"price_combo": 2, "resolve_legs": 1}


def test_first_parse_failure_per_stage_is_loud_rest_are_quiet(caplog):
    """Volume policy: an unexpected exception is a WARNING the FIRST time
    per (book, stage, fetch); repeats drop to DEBUG so one broken book
    cannot drown bot.log — the aggregate rides in the summary line."""
    c = FetchCounters(book="prophetx", path="sweep")
    logger = logging.getLogger("test.parse_volume")
    with caplog.at_level(logging.DEBUG, logger="test.parse_volume"):
        for i in range(5):
            c.record_parse_failure("price_combo", logger,
                                   exc=ValueError(f"boom{i}"))

    levels = [r.levelno for r in caplog.records]
    assert levels[0] == logging.WARNING
    assert levels[1:] == [logging.DEBUG] * 4


def test_record_parse_failure_returns_running_count():
    c = FetchCounters(book="betmgm", path="sweep")
    logger = logging.getLogger("test.parse_count")
    assert c.record_parse_failure("x", logger, exc=ValueError()) == 1
    assert c.record_parse_failure("x", logger, exc=ValueError()) == 2


# ------------------------------------------------------------------ #
# Tripwires                                                           #
# ------------------------------------------------------------------ #

def test_tripwire_fires_when_nothing_resolves():
    """THE FanDuel-alt-total case: the book offered lines, our parser
    understood none of them."""
    c = FetchCounters(book="fanduel", path="sweep")
    c.bump("events_seen", 12)
    c.bump("events_matched", 12)
    c.bump("legs_attempted", 40)
    assert "legs_unresolved" in c.tripwires()


def test_no_tripwire_when_some_legs_resolve():
    c = FetchCounters(book="fanduel", path="sweep")
    c.bump("legs_attempted", 40)
    c.bump("legs_resolved", 1)
    c.bump("targets_attempted", 1)
    c.bump("prices_returned", 1)
    assert c.tripwires() == []


def test_tripwire_fires_when_no_price_comes_back():
    c = FetchCounters(book="novig", path="sweep")
    c.bump("legs_attempted", 4)
    c.bump("legs_resolved", 4)
    c.bump("targets_attempted", 8)
    assert "prices_empty" in c.tripwires()


def test_no_price_tripwire_when_transport_already_failed():
    """A transport error is #33's story, already logged loudly there.
    Double-reporting it as a PARSE tripwire would send a fixer to the
    wrong file."""
    c = FetchCounters(book="novig", path="sweep")
    c.bump("legs_attempted", 4)
    c.bump("legs_resolved", 4)
    c.bump("targets_attempted", 8)
    c.bump("transport_errors", 1)
    assert "prices_empty" not in c.tripwires()


def test_tripwire_fires_when_events_listed_but_none_matched():
    """The book is up and listing games, but none map to our slate —
    a team-name/canonical-match regression, invisible until now."""
    c = FetchCounters(book="caesars", path="sweep")
    c.bump("events_seen", 14)
    assert "events_unmatched" in c.tripwires()


def test_no_event_tripwire_on_an_off_day():
    """Zero events listed is an off-day, not a regression."""
    c = FetchCounters(book="caesars", path="sweep")
    assert c.tripwires() == []


def test_idle_counters_trip_nothing():
    """A book skipped this cycle (nothing attempted) must stay silent."""
    c = FetchCounters(book="betmgm", path="on_demand")
    assert c.tripwires() == []


# ------------------------------------------------------------------ #
# Emission                                                            #
# ------------------------------------------------------------------ #

def test_emit_logs_one_summary_line(caplog):
    logger = logging.getLogger("test.emit_summary")
    c = FetchCounters(book="draftkings", path="sweep")
    c.bump("events_seen", 3)
    c.bump("events_matched", 3)
    c.bump("legs_attempted", 10)
    c.bump("legs_resolved", 9)
    c.bump("targets_attempted", 9)
    c.bump("prices_returned", 36)
    with caplog.at_level(logging.INFO, logger="test.emit_summary"):
        c.emit(logger)

    summaries = [r for r in caplog.records if "sgp fetch" in r.getMessage()]
    assert len(summaries) == 1
    msg = summaries[0].getMessage()
    assert summaries[0].levelno == logging.INFO
    for fragment in ("book=draftkings", "path=sweep", "events_matched=3",
                     "legs_attempted=10", "legs_resolved=9",
                     "prices_returned=36"):
        assert fragment in msg


def test_emit_logs_tripwire_warning(caplog):
    logger = logging.getLogger("test.emit_tripwire")
    c = FetchCounters(book="fanduel", path="sweep")
    c.bump("events_seen", 5)
    c.bump("events_matched", 5)
    c.bump("legs_attempted", 20)
    with caplog.at_level(logging.INFO, logger="test.emit_tripwire"):
        c.emit(logger)

    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "PARSE TRIPWIRE" in warnings[0].getMessage()
    assert "book=fanduel" in warnings[0].getMessage()
    assert "legs_unresolved" in warnings[0].getMessage()


def test_emit_summary_level_is_caller_chosen(caplog):
    """On-demand runs per RFQ per book, so its summary is DEBUG; the
    sweep's runs once per book per cycle and is INFO."""
    logger = logging.getLogger("test.emit_level")
    c = FetchCounters(book="novig", path="on_demand")
    c.bump("legs_attempted", 2)
    c.bump("legs_resolved", 2)
    with caplog.at_level(logging.DEBUG, logger="test.emit_level"):
        c.emit(logger, level=logging.DEBUG)
    assert [r.levelno for r in caplog.records] == [logging.DEBUG]


def test_emit_is_idempotent(caplog):
    """price_sgps emits in a finally; a nested emit must not double-log."""
    logger = logging.getLogger("test.emit_once")
    c = FetchCounters(book="novig", path="sweep")
    c.bump("legs_attempted", 2)
    with caplog.at_level(logging.INFO, logger="test.emit_once"):
        c.emit(logger)
        c.emit(logger)
    assert len([r for r in caplog.records if "sgp fetch" in r.getMessage()]) == 1


# ------------------------------------------------------------------ #
# Context manager                                                     #
# ------------------------------------------------------------------ #

def test_fetch_counters_emits_on_normal_exit(caplog):
    logger = logging.getLogger("test.ctx_ok")
    with caplog.at_level(logging.INFO, logger="test.ctx_ok"):
        with fetch_counters("caesars", "sweep", logger) as c:
            c.bump("legs_attempted", 4)
            c.bump("legs_resolved", 4)
            c.bump("targets_attempted", 4)
            c.bump("prices_returned", 16)
    assert any("sgp fetch" in r.getMessage() for r in caplog.records)


def test_fetch_counters_emits_on_early_return(caplog):
    """Every price_sgps has several ``return []`` short-circuits; the
    summary must survive all of them."""
    logger = logging.getLogger("test.ctx_return")

    def run():
        with fetch_counters("caesars", "sweep", logger) as c:
            c.bump("events_seen", 3)
            return []

    with caplog.at_level(logging.INFO, logger="test.ctx_return"):
        assert run() == []
    assert any("sgp fetch" in r.getMessage() for r in caplog.records)


def test_fetch_counters_counts_transport_error_and_reraises(caplog):
    """A BookTransportError still gets a summary on its way out, and is
    tallied so the price tripwire doesn't fire on top of it."""
    from mlb_sgp._shared import BookTransportError

    logger = logging.getLogger("test.ctx_transport")
    seen = {}
    with caplog.at_level(logging.INFO, logger="test.ctx_transport"):
        with pytest.raises(BookTransportError):
            with fetch_counters("draftkings", "sweep", logger) as c:
                seen["c"] = c
                c.bump("targets_attempted", 9)
                raise BookTransportError("draftkings", "price")

    assert seen["c"].snapshot().transport_errors == 1
    assert not any("PARSE TRIPWIRE" in r.getMessage() for r in caplog.records)


def test_fetch_counters_emits_on_unexpected_exception(caplog):
    logger = logging.getLogger("test.ctx_boom")
    with caplog.at_level(logging.INFO, logger="test.ctx_boom"):
        with pytest.raises(ValueError):
            with fetch_counters("novig", "sweep", logger) as c:
                c.bump("legs_attempted", 2)
                raise ValueError("boom")
    assert any("sgp fetch" in r.getMessage() for r in caplog.records)


def test_snapshot_type_is_the_shared_frozen_dataclass():
    """#38 will persist these; the type must be a stable frozen record."""
    c = FetchCounters(book="novig", path="sweep")
    assert isinstance(c.snapshot(), FetchCountersSnapshot)
