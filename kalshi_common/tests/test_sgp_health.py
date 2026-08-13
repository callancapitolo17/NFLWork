"""Per-book fetch-health history (issue #38).

The table this covers is the memory the epic was missing: the taker's feed
died 2026-06-28 and nobody noticed for four weeks because NOTHING recorded
that a book had stopped answering. #33 made a dead book raise, #35 made a
silent parser countable — this ticket persists both so a streak is visible
after the fact and #37 can alarm on it.

Two invariants dominate these tests:

  ONE LOGICAL FETCH = ONE ROW.  Not one per wire attempt (#34 retries inside
  a fetch), and not two when a deadline fires while the worker thread is
  still running.

  A HEALTH WRITE MUST NEVER REACH PRICING.  record() only buffers; flush()
  swallows. A locked DB costs us history, never a quote.
"""
from __future__ import annotations

import inspect
import json
import logging
import subprocess
import sys
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb
import pytest

from kalshi_common import health_queries, sgp_health
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_health import FetchHealthRecorder
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import (BookTransportError, FetchCounters, GameRef,
                             PricedRow, ResolvedLeg, TargetLine)

EVT = "KXMLBGAME-26JUL09NYYBOS"
CT = datetime(2026, 7, 9, 23, 0, tzinfo=timezone.utc)
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=CT)


# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #

def _target(game_id="g1", total=8.5):
    return TargetLine(game_id=game_id, home_team="Boston Red Sox",
                      away_team="New York Yankees", commence_time=CT,
                      period="FG", spread=-1.5, total=total)


def _row(combo="Home -1.5 + Over"):
    return PricedRow(game_id="g1", combo=combo, period="FG", spread_line=-1.5,
                     total_line=8.5, bookmaker="caesars", source="caesars_sgp",
                     sgp_decimal=2.5, sgp_american=150, fetch_time=CT)


def _legs(n=2):
    return [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over"),
            CanonicalLeg(EVT, "ml", None, "home")][:n]


def _resolved(n=2):
    return [ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(n)]


# Four joint fairs summing to 1.0 under a uniform 1.2 overround, so the
# partition devig actually succeeds and the fetch reads as 'ok'.
CELL_DECS = [1.0 / (1.2 * p) for p in (0.40, 0.20, 0.25, 0.15)]


def _hooks(resolved, price_fn, event="EV", structure={"fake": "s"}):
    return {"match_event": lambda client, game: event,
            "build_structure": lambda client, ev, game: structure,
            "resolve": lambda structure, legs, home, away: resolved,
            "price": price_fn}


def _read(db_path) -> list[dict]:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        cols = ["fetched_at", "book", "path", "outcome", "rows_or_prices",
                "duration_sec", "error_class", "counters_json"]
        rows = con.execute(
            f"SELECT {', '.join(cols)} FROM sgp_fetch_health "
            "ORDER BY fetched_at, book").fetchall()
        return [dict(zip(cols, r)) for r in rows]
    finally:
        con.close()


@pytest.fixture
def db(tmp_path):
    return tmp_path / "market.duckdb"


# ------------------------------------------------------------------ #
# Schema — TIMESTAMPTZ from day one (the #43 trap)                    #
# ------------------------------------------------------------------ #

def test_fetched_at_column_is_timestamptz(db):
    """mlb_sgp_odds.fetch_time being a naive TIMESTAMP is the whole reason
    #43 exists. Do not create a second instance of that bug."""
    rec = FetchHealthRecorder(str(db))
    rec.record(book="caesars", path="sweep", outcome="ok", rows_or_prices=4,
               duration_sec=1.0)
    rec.flush()

    con = duckdb.connect(str(db), read_only=True)
    try:
        types = dict(con.execute(
            "SELECT column_name, data_type FROM information_schema.columns "
            "WHERE table_name = 'sgp_fetch_health'").fetchall())
    finally:
        con.close()
    assert types["fetched_at"] == "TIMESTAMP WITH TIME ZONE"


def test_fetched_at_roundtrips_as_utc_aware(db):
    before = datetime.now(timezone.utc)
    rec = FetchHealthRecorder(str(db))
    rec.record(book="caesars", path="sweep", outcome="ok", rows_or_prices=4,
               duration_sec=1.0)
    rec.flush()
    after = datetime.now(timezone.utc)

    (row,) = _read(db)
    assert row["fetched_at"].tzinfo is not None
    assert before <= row["fetched_at"] <= after


# ------------------------------------------------------------------ #
# Recorder unit behavior                                              #
# ------------------------------------------------------------------ #

def test_record_buffers_until_flush(db):
    rec = FetchHealthRecorder(str(db))
    rec.record(book="novig", path="on_demand", outcome="ok",
               rows_or_prices=4, duration_sec=0.5)
    assert rec.pending() == 1
    assert not db.exists(), "record() must not touch the DB"
    assert rec.flush() == 1
    assert rec.pending() == 0
    assert len(_read(db)) == 1


def test_disabled_recorder_writes_nothing(db):
    """The dashboard CLI path and every test that doesn't ask for health
    must never create a DB file or grow a buffer."""
    rec = FetchHealthRecorder(None)
    assert not rec.enabled
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1)
    assert rec.pending() == 0
    assert rec.flush() == 0
    assert not db.exists()


def test_counters_are_serialized_as_json(db):
    counters = FetchCounters("fanduel", "sweep")
    counters.bump("events_seen", 15)
    counters.bump("events_matched", 14)
    counters.bump("legs_attempted", 120)
    counters.record_parse_failure("runners", logging.getLogger("t"),
                                  ValueError("bad"))
    rec = FetchHealthRecorder(str(db))
    rec.record(book="fanduel", path="sweep", outcome="empty",
               rows_or_prices=0, duration_sec=2.0,
               counters=counters.snapshot())
    rec.flush()

    payload = json.loads(_read(db)[0]["counters_json"])
    assert payload["events_seen"] == 15
    assert payload["events_matched"] == 14
    # Zeros are stored, not omitted: a query that has to COALESCE a missing
    # key is a query that will silently read 0 as "no data".
    assert payload["legs_resolved"] == 0
    assert payload["parse_failures"] == 1
    assert payload["by_stage"] == {"runners": 1}


def test_buffer_cap_drops_oldest(db):
    rec = FetchHealthRecorder(str(db), buffer_max=3)
    for i in range(5):
        rec.record(book=f"b{i}", path="sweep", outcome="ok",
                   rows_or_prices=i, duration_sec=0.1)
    assert rec.pending() == 3
    rec.flush()
    assert [r["book"] for r in _read(db)] == ["b2", "b3", "b4"]


def test_flush_failure_retains_buffer_and_never_raises(db, monkeypatch):
    rec = FetchHealthRecorder(str(db))
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1)

    def boom(*a, **k):
        raise duckdb.IOException("database is locked")

    monkeypatch.setattr(sgp_health.duckdb, "connect", boom)
    assert rec.flush() == 0          # no raise
    assert rec.pending() == 1        # retained for the next attempt

    monkeypatch.undo()
    assert rec.flush() == 1
    assert len(_read(db)) == 1


def test_record_survives_an_unserializable_payload(db):
    """A broken __repr__ inside a counters payload must not reach pricing."""
    class Hostile:
        def __repr__(self):
            raise RuntimeError("nope")

    rec = FetchHealthRecorder(str(db))
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1, counters=Hostile())
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=2,
               duration_sec=0.1)
    rec.flush()
    # The hostile row is dropped; the healthy one still lands.
    assert [r["rows_or_prices"] for r in _read(db)] == [2]


# ------------------------------------------------------------------ #
# Retention                                                           #
# ------------------------------------------------------------------ #

def test_prune_deletes_rows_past_retention(db):
    rec = FetchHealthRecorder(str(db), retention_days=30,
                              prune_interval_sec=0.0)
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1)
    rec.flush()
    con = duckdb.connect(str(db))
    try:
        con.execute(
            "INSERT INTO sgp_fetch_health " +
            "(fetched_at, book, path, outcome, rows_or_prices, duration_sec, error_class, counters_json) "
            "VALUES (?, 'novig', 'sweep', 'ok', 1, 0.1, NULL, NULL)",
            [datetime.now(timezone.utc) - timedelta(days=31)])
        con.execute(
            "INSERT INTO sgp_fetch_health " +
            "(fetched_at, book, path, outcome, rows_or_prices, duration_sec, error_class, counters_json) "
            "VALUES (?, 'novig', 'sweep', 'ok', 1, 0.1, NULL, NULL)",
            [datetime.now(timezone.utc) - timedelta(days=29)])
    finally:
        con.close()

    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1)
    rec.flush()
    assert len(_read(db)) == 3      # 2 fresh + the 29-day-old one


def test_prune_is_rate_limited(db):
    """flush() runs ~4x/sec in the maker loop; a DELETE scan on every one is
    pure waste."""
    clock = [1000.0]
    rec = FetchHealthRecorder(str(db), prune_interval_sec=3600.0,
                              now_fn=lambda: clock[0])
    pruned = []
    real_prune = rec._prune
    rec._prune = lambda con: pruned.append(clock[0]) or real_prune(con)

    for _ in range(3):
        rec.record(book="novig", path="sweep", outcome="ok",
                   rows_or_prices=1, duration_sec=0.1)
        rec.flush()
    assert len(pruned) == 1, "pruned more than once inside the interval"

    clock[0] += 3601.0
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=1,
               duration_sec=0.1)
    rec.flush()
    assert len(pruned) == 2


# ------------------------------------------------------------------ #
# Sweep path — outcome attribution                                    #
# ------------------------------------------------------------------ #

def _svc_sweep(db, runner, book="caesars", **kw):
    return SGPService(books=(book,), runners={book: runner},
                      health_db_path=str(db), **kw)


def test_sweep_rows_are_ok(db):
    svc = _svc_sweep(db, lambda t: [_row("a"), _row("b")])
    svc.refresh([_target()])
    (row,) = _read(db)
    assert row["book"] == "caesars"
    assert row["path"] == "sweep"
    assert row["outcome"] == "ok"
    assert row["rows_or_prices"] == 2
    assert row["duration_sec"] >= 0.0
    assert row["error_class"] is None


def test_sweep_empty_is_not_ok(db):
    """An empty slate and a healthy fetch are different facts. Conflating
    them is the meta-bug this epic exists to remove."""
    svc = _svc_sweep(db, lambda t: [])
    svc.refresh([_target()])
    (row,) = _read(db)
    assert row["outcome"] == "empty"
    assert row["rows_or_prices"] == 0


def test_sweep_transport_error_records_stage_and_status(db):
    def dead(targets):
        raise BookTransportError("caesars", "events", status_code=403)

    svc = _svc_sweep(db, dead)
    svc.refresh([_target()])
    (row,) = _read(db)
    assert row["outcome"] == "transport_error"
    assert row["error_class"] == "BookTransportError:events:403"


def test_sweep_transport_error_without_status(db):
    def dead(targets):
        raise BookTransportError("caesars", "price",
                                 detail="all 12 price calls failed")

    svc = _svc_sweep(db, dead)
    svc.refresh([_target()])
    assert _read(db)[0]["error_class"] == "BookTransportError:price"


def test_sweep_generic_exception_is_error_not_transport(db):
    """A TypeError in our own parse code is not a dead endpoint. Filing it
    as transport_error would point a fixer at the network — exactly the
    misdiagnosis #33 exists to remove."""
    def boom(targets):
        raise ValueError("our bug")

    svc = _svc_sweep(db, boom)
    svc.refresh([_target()])
    (row,) = _read(db)
    assert row["outcome"] == "error"
    assert row["error_class"] == "ValueError"


def test_book_skipped_by_min_refresh_writes_no_row(db):
    """A book that was never fetched has no health to report. Writing
    'empty' for it would fabricate a failure streak out of scheduling."""
    clock = [1000.0]
    svc = SGPService(books=("caesars",), runners={"caesars": lambda t: [_row()]},
                     min_refresh_sec={"caesars": 300}, health_db_path=str(db),
                     now_fn=lambda: clock[0])
    svc.refresh([_target()])          # runs
    clock[0] += 10.0
    svc.refresh([_target()])          # skipped — too soon
    assert len(_read(db)) == 1


def test_deadline_timeout_writes_exactly_one_row(db):
    """The worker thread is STILL RUNNING when the deadline fires. If the
    row were written by the worker, a slow book would land two rows per
    fetch (one 'timeout', one late 'ok') and every streak count would lie.
    """
    release = threading.Event()
    finished = threading.Event()

    def slow(targets):
        release.wait(10)
        finished.set()
        return [_row()]

    svc = SGPService(books=("caesars",), runners={"caesars": slow},
                     per_book_deadline_sec=0.2, health_db_path=str(db))
    results = svc.refresh([_target()])
    assert results["caesars"] is None

    release.set()
    assert finished.wait(10), "worker never finished"
    time.sleep(0.2)                  # give a stray writer a chance to appear
    svc.flush_health(force=True)

    rows = _read(db)
    assert len(rows) == 1, f"expected one row per fetch, got {rows}"
    assert rows[0]["outcome"] == "timeout"
    assert rows[0]["error_class"] == "FutureTimeout"


def test_duration_covers_the_whole_logical_fetch_including_retry_sleep(db):
    """#34 sleeps INSIDE a fetch. duration_sec must span the sleeps, and the
    retried fetch must still be one row, not one per attempt."""
    from mlb_sgp import _shared

    slept = []

    def flaky(targets):
        # Stand-in for request_with_retry's backoff between two attempts.
        _shared._sleep(0.05)
        slept.append(0.05)
        return [_row()]

    svc = _svc_sweep(db, flaky)
    svc.refresh([_target()])

    rows = _read(db)
    assert len(rows) == 1
    assert slept, "the fetch never retried — test is not measuring anything"
    assert rows[0]["duration_sec"] >= 0.05


def test_sweep_row_carries_the_orchestrator_counters(db, monkeypatch):
    """The sweep's FetchCounters lives inside price_sgps. If the service
    can't reach it, every sweep row loses its parse-side story — the
    FanDuel-regression signal #35 built."""
    import mlb_sgp.caesars as czr

    def fake_price_sgps(targets, periods=("FG",), client=None, verbose=False,
                        *, counters=None):
        counters.bump("events_seen", 7)
        counters.bump("legs_attempted", 3)
        return [_row()]

    monkeypatch.setattr(czr, "price_sgps", fake_price_sgps)

    svc = SGPService(books=("caesars",), health_db_path=str(db))
    svc._state["caesars"].client = object()
    monkeypatch.setattr(svc, "_ensure_client", lambda b: svc._state[b])
    svc.refresh([_target()])

    payload = json.loads(_read(db)[0]["counters_json"])
    assert payload["events_seen"] == 7
    assert payload["legs_attempted"] == 3


# ------------------------------------------------------------------ #
# On-demand path — the priority surface under epic #53                #
# ------------------------------------------------------------------ #

def _svc_od(db, hooks, book="novig"):
    return SGPService(books=(book,), on_demand_hooks={book: hooks},
                      health_db_path=str(db))


def test_on_demand_ok_row_carries_counters(db):
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        return CELL_DECS[len(calls) - 1]

    svc = _svc_od(db, _hooks(_resolved(2), price))
    assert svc.price_on_demand("novig", GAME, _legs(2)) is not None
    svc.flush_health(force=True)

    (row,) = _read(db)
    assert row["book"] == "novig"
    assert row["path"] == "on_demand"
    assert row["outcome"] == "ok"
    assert row["rows_or_prices"] == 4      # price calls that returned
    payload = json.loads(row["counters_json"])
    assert payload["targets_attempted"] == 4
    assert payload["prices_returned"] == 4


def test_on_demand_decline_is_empty(db):
    """The book simply doesn't list the game — not a failure, but not a
    successful price either."""
    svc = _svc_od(db, _hooks(_resolved(2), lambda c, r, e: 2.0, event=None))
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    svc.flush_health(force=True)
    (row,) = _read(db)
    assert row["outcome"] == "empty"
    assert row["error_class"] is None


def test_on_demand_transport_error_row(db):
    def boom(client, game):
        raise BookTransportError("novig", "events", status_code=503)

    hooks = _hooks(_resolved(2), lambda c, r, e: 2.0)
    hooks["match_event"] = boom
    svc = _svc_od(db, hooks)
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    svc.flush_health(force=True)
    (row,) = _read(db)
    assert row["outcome"] == "transport_error"
    assert row["error_class"] == "BookTransportError:events:503"


def test_on_demand_generic_exception_is_error(db):
    def boom(client, game):
        raise ValueError("our bug")

    hooks = _hooks(_resolved(2), lambda c, r, e: 2.0)
    hooks["match_event"] = boom
    svc = _svc_od(db, hooks)
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    svc.flush_health(force=True)
    (row,) = _read(db)
    assert row["outcome"] == "error"
    assert row["error_class"] == "ValueError"


def test_on_demand_no_fetch_writes_no_row(db):
    """Unknown book / empty leg set never touches a wire."""
    svc = _svc_od(db, _hooks(_resolved(2), lambda c, r, e: 2.0))
    assert svc.price_on_demand("novig", GAME, []) is None
    assert svc.price_on_demand("betmgm", GAME, _legs(2)) is None
    svc.flush_health(force=True)
    assert not db.exists() or _read(db) == []


def test_on_demand_duration_is_recorded(db):
    def slow_price(client, refs, event):
        time.sleep(0.02)
        return CELL_DECS[0]

    svc = _svc_od(db, _hooks(_resolved(1), slow_price))
    svc.price_on_demand("novig", GAME, _legs(1))
    svc.flush_health(force=True)
    assert _read(db)[0]["duration_sec"] >= 0.02


# ------------------------------------------------------------------ #
# Acceptance criteria from the ticket                                 #
# ------------------------------------------------------------------ #

def test_three_on_demand_fetches_plus_one_sweep_are_attributed(db):
    dead_calls = []

    def dead(client, game):
        dead_calls.append(1)
        raise BookTransportError("novig", "structure", status_code=403)

    dead_hooks = _hooks(_resolved(2), lambda c, r, e: 2.0)
    dead_hooks["match_event"] = dead

    svc = SGPService(books=("novig", "caesars"),
                     runners={"novig": lambda t: [],
                              "caesars": lambda t: [_row(), _row("b")]},
                     on_demand_hooks={"novig": dead_hooks},
                     health_db_path=str(db))

    for _ in range(3):
        svc.price_on_demand("novig", GAME, _legs(2))
    svc.refresh([_target()])
    svc.flush_health(force=True)

    rows = _read(db)
    assert len(rows) == 5, rows
    on_demand = [r for r in rows if r["path"] == "on_demand"]
    sweep = [r for r in rows if r["path"] == "sweep"]
    assert len(on_demand) == 3 and len(sweep) == 2

    # A dead book shows a transport_error streak, attributed to one book.
    assert {r["book"] for r in on_demand} == {"novig"}
    assert [r["outcome"] for r in on_demand] == ["transport_error"] * 3
    assert {r["error_class"] for r in on_demand} == {
        "BookTransportError:structure:403"}
    # ...while the healthy book on the same tick reads clean.
    by_book = {r["book"]: r for r in sweep}
    assert by_book["caesars"]["outcome"] == "ok"
    assert by_book["novig"]["outcome"] == "empty"


def test_dead_book_verdict_through_a_real_orchestrator(db, monkeypatch):
    """End-to-end shape of the DK-403 rot: transport is fine at the event
    level, every PRICE call fails, and #33's PriceCallTally verdict is what
    turns that into a transport verdict. A MagicMock client can never reach
    it (the CLIENT records at CZR), so the stub owns a real tally.
    """
    from dataclasses import replace
    from types import SimpleNamespace

    from mlb_sgp import caesars
    from mlb_sgp._shared import PriceCallTally

    event = SimpleNamespace(event_id="e1", home_team="New York Yankees",
                            away_team="Boston Red Sox", start_time=CT)
    monkeypatch.setattr(caesars, "_match_events",
                        lambda events, targets: {"g1": event})
    leg = {"price": {"decimal": 1.9}}
    totals = {8.5 + i: {"over": leg, "under": leg} for i in range(4)}
    monkeypatch.setattr(caesars, "parse_markets", lambda ev: {
        "FG": {"spreads": {-1.5: {"home": leg, "away": leg}},
               "totals": totals, "moneyline": None},
        "F5": {"spreads": {}, "totals": {}, "moneyline": None}})

    class DeadPriceClient:
        """Prices always fail — recorded on a REAL tally, as the client does."""
        def __init__(self):
            self.price_calls = PriceCallTally("caesars")

        def list_events(self, *a, **k):
            return [event]

        def fetch_event(self, *a, **k):
            return {"id": "e1"}

        def price_combo(self, *a, **k):
            self.price_calls.record(False)
            return None

    svc = SGPService(books=("caesars",), health_db_path=str(db))
    svc._state["caesars"].client = DeadPriceClient()
    monkeypatch.setattr(svc, "_ensure_client", lambda b: svc._state[b])

    targets = [replace(_target(), total=8.5 + i) for i in range(4)]
    svc.refresh(targets)

    (row,) = _read(db)
    assert row["outcome"] == "transport_error"
    assert row["error_class"].startswith("BookTransportError:price")


def test_held_write_lock_does_not_affect_pricing(db):
    """The health table shares mlb_sgp_odds' write lock. A locked DB must
    cost us history and nothing else."""
    holder = subprocess.Popen(
        [sys.executable, "-c",
         "import duckdb, sys, time\n"
         "con = duckdb.connect(sys.argv[1])\n"
         "con.execute('CREATE TABLE IF NOT EXISTS holder(i INTEGER)')\n"
         "print('locked', flush=True)\n"
         "time.sleep(60)\n",
         str(db)],
        stdout=subprocess.PIPE, text=True)
    try:
        assert holder.stdout.readline().strip() == "locked"

        seq = iter(CELL_DECS)        # deterministic 4-cell partition
        svc = SGPService(books=("novig",),
                         on_demand_hooks={"novig": _hooks(
                             _resolved(2), lambda c, r, e: next(seq))},
                         health_db_path=str(db))

        locked_result = svc.price_on_demand("novig", GAME, _legs(2))
        svc.flush_health(force=True)   # must not raise
        assert locked_result is not None
        assert svc.health_pending() == 1, "row must be retained for retry"
    finally:
        holder.kill()
        holder.wait(10)

    # Same pricing, no lock: identical fair, and now the history lands.
    seq2 = iter(CELL_DECS)
    svc2 = SGPService(books=("novig",),
                      on_demand_hooks={"novig": _hooks(
                          _resolved(2), lambda c, r, e: next(seq2))},
                      health_db_path=str(db))
    free_result = svc2.price_on_demand("novig", GAME, _legs(2))
    assert free_result.fair == pytest.approx(locked_result.fair)

    svc.flush_health(force=True)     # the retained row finally lands
    assert any(r["path"] == "on_demand" for r in _read(db))


# ------------------------------------------------------------------ #
# Signature guards — new params must be keyword-only                  #
# ------------------------------------------------------------------ #

def test_flush_health_is_coalesced(db):
    """Each flush holds a READ-WRITE handle on the market DB, and duckdb
    1.4.4 fails an in-process read-only connect for as long as one is open.
    The maker ticks 4x/sec; opening the pricing path's file that often for
    diagnostics is a window worth not having."""
    from kalshi_common import sgp_service

    svc = SGPService(books=("novig",),
                     on_demand_hooks={"novig": _hooks(
                         _resolved(1), lambda c, r, e: CELL_DECS[0])},
                     health_db_path=str(db))
    svc.price_on_demand("novig", GAME, _legs(1))
    assert svc.flush_health() == 1            # first tick drains

    svc.price_on_demand("novig", GAME, _legs(1))
    assert svc.flush_health() == 0, "flushed again inside the interval"
    assert svc.health_pending() == 1

    svc._last_health_flush -= sgp_service.HEALTH_FLUSH_MIN_INTERVAL_SEC + 1
    assert svc.flush_health() == 1


def test_refresh_flushes_sweep_rows_immediately(db):
    """The sweep runs once a minute on the tick thread that is about to take
    the same write lock anyway — no reason to make its rows wait."""
    svc = _svc_sweep(db, lambda t: [_row()])
    svc.refresh([_target()])
    assert _read(db)[0]["outcome"] == "ok"      # no flush_health() call
    assert svc.health_pending() == 0


def test_run_book_counters_is_keyword_only():
    """A third positional here would bind to `counters` — the same hazard
    #35 hit with BetMGM's fixture_id."""
    param = inspect.signature(SGPService._run_book).parameters["counters"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY


def test_health_db_path_is_keyword_only():
    param = inspect.signature(SGPService.__init__).parameters["health_db_path"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is None


def test_service_without_health_path_still_prices(db):
    """Every existing caller (dashboard, tests, the pre-#38 bots) builds the
    service with no health path and must be untouched."""
    svc = SGPService(books=("caesars",), runners={"caesars": lambda t: [_row()]})
    assert len(svc.refresh([_target()])["caesars"]) == 1
    svc.flush_health(force=True)
    assert not db.exists()


# ------------------------------------------------------------------ #
# The ready-made health query (fix-spec item 5)                       #
# ------------------------------------------------------------------ #

# The parser moved into kalshi_common.health_queries for #37 — the monitor
# needs the same named queries, and a second copy of the loader is a second
# place for the streak predicate to drift.
QUERY_FILE = health_queries.QUERY_FILE


def _named_query(name: str) -> str:
    return health_queries.named_query(name)


def _all_query_names() -> list[str]:
    return sorted(health_queries.all_queries())


def test_every_named_query_executes(db):
    """A shipped .sql file nobody runs is a broken .sql file. Seed one row
    of each shape and execute all of them."""
    rec = FetchHealthRecorder(str(db))
    rec.record(book="novig", path="sweep", outcome="ok", rows_or_prices=40,
               duration_sec=12.0, counters=FetchCounters("novig",
                                                         "sweep").snapshot())
    rec.record(book="novig", path="on_demand", outcome="transport_error",
               rows_or_prices=0, duration_sec=0.4,
               error_class="BookTransportError:price:403")
    rec.flush()

    names = _all_query_names()
    assert len(names) >= 4, names
    con = duckdb.connect(str(db), read_only=True)
    try:
        for name in names:
            con.execute(_named_query(name)).fetchall()   # must not raise
    finally:
        con.close()


def test_failure_streak_query_counts_only_the_current_run(db):
    """An old outage must not inflate today's streak, and a book that has
    since recovered must report nothing (#37 alerts off this)."""
    con = duckdb.connect(str(db))
    try:
        con.execute(sgp_health.SCHEMA_SQL)
        base = datetime.now(timezone.utc)
        seq = [
            # novig: two old failures, then a recovery, then 3 fresh failures
            ("novig", -50, "transport_error"), ("novig", -49, "transport_error"),
            ("novig", -40, "ok"),
            ("novig", -3, "transport_error"), ("novig", -2, "transport_error"),
            ("novig", -1, "transport_error"),
            # caesars: failed, then recovered — streak 0, must not appear
            ("caesars", -5, "timeout"), ("caesars", -1, "ok"),
        ]
        for book, mins, outcome in seq:
            con.execute(
                "INSERT INTO sgp_fetch_health " +
                "(fetched_at, book, path, outcome, rows_or_prices, duration_sec, error_class, counters_json) "
                "VALUES (?, ?, 'on_demand', ?, 0, 1.0, NULL, NULL)",
                [base + timedelta(minutes=mins), book, outcome])
    finally:
        con.close()

    con = duckdb.connect(str(db), read_only=True)
    try:
        rows = con.execute(_named_query("current_failure_streak")).fetchall()
        cols = [d[0] for d in con.description]
    finally:
        con.close()

    streaks = {r[cols.index("book")]: r[cols.index("failing_streak")]
               for r in rows}
    assert streaks == {"novig": 3}, streaks


def test_health_query_groups_by_book_and_path(db):
    """targets_attempted/prices_returned carry DIFFERENT UNITS per path
    (sweep ~1:4, on-demand 1:1 — see _shared.COUNTER_NAMES). A query that
    averages across paths reports a wrong ratio for every book."""
    rec = FetchHealthRecorder(str(db))
    for i in range(4):
        rec.record(book="novig", path="sweep", outcome="ok",
                   rows_or_prices=40, duration_sec=10.0 + i)
    for i in range(4):
        rec.record(book="novig", path="on_demand", outcome="ok",
                   rows_or_prices=4, duration_sec=1.0 + i)
    rec.record(book="novig", path="on_demand", outcome="transport_error",
               rows_or_prices=0, duration_sec=0.3,
               error_class="BookTransportError:price")
    rec.flush()

    con = duckdb.connect(str(db), read_only=True)
    try:
        rows = con.execute(_named_query("outcome_mix_24h")).fetchall()
        cols = [d[0] for d in con.description]
    finally:
        con.close()

    keyed = {(r[cols.index("book")], r[cols.index("path")]): r for r in rows}
    assert set(keyed) == {("novig", "sweep"), ("novig", "on_demand")}
    sweep = keyed[("novig", "sweep")]
    on_demand = keyed[("novig", "on_demand")]
    assert sweep[cols.index("n_fetches")] == 4
    assert on_demand[cols.index("n_fetches")] == 5
    assert on_demand[cols.index("n_transport_error")] == 1
    # p50/p95 must be computed within a path, not pooled.
    assert sweep[cols.index("p50_duration_sec")] > 5.0
    assert on_demand[cols.index("p50_duration_sec")] < 5.0
