"""Run-time book-health alerting (issue #37).

#38 gave us a memory of every fetch; this ticket is the alarm bell on top of
it. The bug it prevents is the one that already happened: the taker's feed
died 2026-06-28 and nobody noticed for four weeks.

Four invariants dominate these tests, and each of them is a way a naive
implementation goes wrong:

  'empty' IS NOT A FAILURE.  Five outcomes exist, not four. A book that
  answers and has no markets (off-day, thin slate, "won't price this combo"
  — which is the DEFAULT on-demand verdict) is healthy. Alerting on
  ``outcome != 'ok'`` pages you every quiet morning.

  ABSENCE IS NOT A FAILURE.  A book skipped by min_refresh_sec writes no
  health row and makes no observation. Streak logic that reads a gap as a
  failure alerts on scheduling.

  ONE NOTIFICATION PER INCIDENT.  Not per cycle. Re-armed only by a
  subsequent ok/empty.

  ZERO ALERTS WHEN THE BOTS ARE DOWN.  Which is their normal state for
  weeks at a time.
"""
from __future__ import annotations

import inspect
import logging
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import pytest

from kalshi_common.book_health import (HEALTHY_OUTCOMES, BookHealthAlerter,
                                       notify_macos)
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import (BookTransportError, GameRef, PricedRow,
                             ResolvedLeg, TargetLine)

EVT = "KXMLBGAME-26JUL09NYYBOS"
CT = datetime(2026, 7, 9, 23, 0, tzinfo=timezone.utc)
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=CT)
SIX_BOOKS = ("draftkings", "fanduel", "prophetx", "novig", "betmgm",
             "caesars")


# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #

class FakeNotifier:
    """Stand-in for the macOS notifier. Records (title, message)."""

    def __init__(self, raises: bool = False):
        self.calls: list[tuple[str, str]] = []
        self.raises = raises
        self._lock = threading.Lock()

    def __call__(self, title: str, message: str) -> None:
        with self._lock:
            self.calls.append((title, message))
        if self.raises:
            raise RuntimeError("osascript exploded")

    @property
    def n(self) -> int:
        with self._lock:
            return len(self.calls)

    @property
    def texts(self) -> str:
        with self._lock:
            return " || ".join(f"{t} {m}" for t, m in self.calls)


def _alerter(notifier=None, **kw):
    """An alerter wired for deterministic tests: notify inline, not on a
    daemon thread (thread behavior has its own test)."""
    kw.setdefault("label", "MM")
    kw.setdefault("notify_async", False)
    return BookHealthAlerter(notifier=notifier or FakeNotifier(), **kw)


def _fail(alerter, book="draftkings", path="sweep", n=1,
          outcome="transport_error", error_class="BookTransportError:events:403"):
    for _ in range(n):
        alerter.observe(book=book, path=path, outcome=outcome,
                        error_class=error_class)


def _target(game_id="g1", total=8.5):
    return TargetLine(game_id=game_id, home_team="Boston Red Sox",
                      away_team="New York Yankees", commence_time=CT,
                      period="FG", spread=-1.5, total=total)


def _row(combo="Home -1.5 + Over", book="caesars"):
    return PricedRow(game_id="g1", combo=combo, period="FG", spread_line=-1.5,
                     total_line=8.5, bookmaker=book, source=f"{book}_sgp",
                     sgp_decimal=2.5, sgp_american=150, fetch_time=CT)


def _legs(n=2):
    return [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over"),
            CanonicalLeg(EVT, "ml", None, "home")][:n]


def _resolved(n=2):
    return [ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(n)]


def _hooks(resolved, price_fn, event="EV", structure={"fake": "s"}):
    return {"match_event": lambda client, game: event,
            "build_structure": lambda client, ev, game: structure,
            "resolve": lambda structure, legs, home, away: resolved,
            "price": price_fn}


# ------------------------------------------------------------------ #
# THE TRAP: five outcomes, not four                                   #
# ------------------------------------------------------------------ #

def test_empty_is_never_a_failure():
    """A quiet morning must not page anyone. 'empty' means the book ANSWERED
    and had no markets — the correct state on an off-day, a thin slate, and
    (on-demand) 'this book won't price that combo', which is the default
    verdict for a plain None."""
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3)
    for _ in range(50):
        a.observe(book="draftkings", path="sweep", outcome="empty")
        a.dispatch()
    assert notifier.n == 0
    assert a.snapshot()["streaks"]["draftkings"]["sweep"] == 0


def test_healthy_outcomes_match_the_shipped_sql_predicate():
    """`current_failure_streak` in fetch_health_queries.sql already encodes
    the right rule. If someone edits one predicate and not the other, the
    monitor and the alerter will disagree about which books are dead."""
    sql = (Path(__file__).resolve().parents[1]
           / "fetch_health_queries.sql").read_text()
    block = sql.split("-- name: current_failure_streak", 1)[1]
    assert "outcome IN ('ok', 'empty')" in block
    assert HEALTHY_OUTCOMES == {"ok", "empty"}


def test_all_five_non_healthy_outcomes_count_as_failures():
    for outcome in ("transport_error", "timeout", "error"):
        notifier = FakeNotifier()
        a = _alerter(notifier, streak_threshold=2)
        _fail(a, n=2, outcome=outcome)
        a.dispatch()
        assert notifier.n == 1, outcome


# ------------------------------------------------------------------ #
# Rule A — book degraded: one notification per INCIDENT               #
# ------------------------------------------------------------------ #

def test_streak_alerts_once_and_only_once(caplog):
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3)

    with caplog.at_level(logging.ERROR, logger="kalshi_common.book_health"):
        for i in range(1, 4):
            _fail(a)
            alerts = a.dispatch()
            if i < 3:
                assert alerts == [], f"alerted after only {i} failures"
        assert notifier.n == 1
        # 20 more failures in the same incident: still one notification.
        _fail(a, n=20)
        a.dispatch()
        assert notifier.n == 1
    errors = [r for r in caplog.records if r.levelno == logging.ERROR]
    assert len(errors) == 1


def test_alert_text_carries_book_path_and_error_class():
    """`BookTransportError:events:403` is what tells DK's Akamai block apart
    from Novig's DNS death four weeks later."""
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3)
    _fail(a, book="draftkings", path="on_demand", n=3,
          error_class="BookTransportError:events:403")
    a.dispatch()
    text = notifier.texts
    assert "draftkings" in text
    assert "on_demand" in text
    assert "BookTransportError:events:403" in text
    assert "3" in text


def test_recovery_rearms_the_alert(caplog):
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3)

    _fail(a, n=3)
    a.dispatch()
    assert notifier.n == 1

    with caplog.at_level(logging.INFO, logger="kalshi_common.book_health"):
        a.observe(book="draftkings", path="sweep", outcome="ok")
        recovered = a.dispatch()
    assert [x.kind for x in recovered] == ["book_recovered"]
    assert notifier.n == 1, "recovery must be log-only, not a second toast"
    assert any("recover" in r.message.lower() for r in caplog.records)

    _fail(a, n=3)
    a.dispatch()
    assert notifier.n == 2, "a NEW incident must re-alert"


def test_empty_also_rearms():
    """'empty' is an answer. A book that comes back with nothing to price is
    reachable again, and the next outage is a new incident."""
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=2)
    _fail(a, n=2)
    a.dispatch()
    a.observe(book="draftkings", path="sweep", outcome="empty")
    a.dispatch()
    _fail(a, n=2)
    a.dispatch()
    assert notifier.n == 2


def test_streak_is_tracked_per_book_and_path():
    notifier = FakeNotifier()
    # min_healthy_books=1 isolates Rule A: with only two books observed here,
    # one degrading would legitimately also trip the capacity rule.
    a = _alerter(notifier, streak_threshold=3, min_healthy_books=1)
    _fail(a, book="draftkings", path="sweep", n=2)
    _fail(a, book="draftkings", path="on_demand", n=2)
    _fail(a, book="novig", path="sweep", n=2)
    a.dispatch()
    assert notifier.n == 0, "streaks leaked across (book, path)"

    _fail(a, book="draftkings", path="sweep", n=1)
    a.dispatch()
    assert notifier.n == 1
    assert "sweep" in notifier.texts and "on_demand" not in notifier.texts


# ------------------------------------------------------------------ #
# Absence is not failure                                              #
# ------------------------------------------------------------------ #

def test_skipped_book_neither_advances_nor_resets_the_streak():
    """min_refresh_sec skips write NO health row. Driven through the real
    service so this can't drift from how the sweep actually behaves."""
    clock = [1000.0]
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=3)

    calls = {"n": 0}

    def runner(targets):
        calls["n"] += 1
        raise BookTransportError("caesars", "events", status_code=403)

    svc = SGPService(books=("caesars",), runners={"caesars": runner},
                     min_refresh_sec={"caesars": 300},
                     now_fn=lambda: clock[0], book_health=alerter)
    # Two real attempts (min_refresh only skips after a SUCCESS, and there
    # are none here), then assert the streak reflects attempts, not ticks.
    svc.refresh([_target()])
    svc.refresh([_target()])
    svc.check_book_health()
    assert notifier.n == 0
    assert alerter.snapshot()["streaks"]["caesars"]["sweep"] == 2

    # 500 ticks in which the book is genuinely never fetched must not
    # manufacture failures out of scheduling.
    for _ in range(500):
        svc.check_book_health()
    assert notifier.n == 0
    assert alerter.snapshot()["streaks"]["caesars"]["sweep"] == 2
    assert calls["n"] == 2


def test_successful_book_skipped_by_min_refresh_stays_healthy():
    clock = [1000.0]
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=3, min_healthy_books=1)
    svc = SGPService(books=("caesars",), runners={"caesars": lambda t: [_row()]},
                     min_refresh_sec={"caesars": 300},
                     now_fn=lambda: clock[0], book_health=alerter)
    svc.refresh([_target()])
    for _ in range(10):
        svc.refresh([_target()])          # all skipped: no rows, no observes
        svc.check_book_health()
    assert notifier.n == 0
    snap = alerter.snapshot()
    assert snap["healthy_books"] == ["caesars"]
    assert snap["dark"] is False


# ------------------------------------------------------------------ #
# Idle bot                                                            #
# ------------------------------------------------------------------ #

def test_idle_bot_emits_nothing():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3, min_healthy_books=2)
    for _ in range(1000):
        assert a.dispatch() == []
    assert notifier.n == 0


def test_disabled_alerter_is_completely_inert():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=1, enabled=False)
    _fail(a, n=10)
    assert a.dispatch() == []
    assert notifier.n == 0


# ------------------------------------------------------------------ #
# Rule B — consensus capacity                                         #
# ------------------------------------------------------------------ #

def _observe_all(a, books, outcome="ok", path="sweep"):
    for b in books:
        a.observe(book=b, path=path, outcome=outcome)


def test_consensus_dark_when_healthy_drops_below_min():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3, min_healthy_books=2)
    _observe_all(a, SIX_BOOKS)
    a.dispatch()
    assert notifier.n == 0

    # Kill four of six: 2 healthy left == the minimum, still quotable.
    for b in SIX_BOOKS[:4]:
        _fail(a, book=b, n=3)
    a.dispatch()
    assert notifier.n == 4, "expected 4 book_degraded, no dark alert yet"
    assert "DARK" not in notifier.texts

    # The fifth death drops us to 1 healthy book — below MIN_AGREEING_BOOKS.
    _fail(a, book=SIX_BOOKS[4], n=3)
    alerts = a.dispatch()
    assert [x.kind for x in alerts] == ["book_degraded", "consensus_dark"]
    assert "DARK" in notifier.texts
    assert a.snapshot()["dark"] is True


def test_consensus_dark_alerts_once_then_restores():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=2, min_healthy_books=2)
    _observe_all(a, ("draftkings", "fanduel", "novig"))
    a.dispatch()
    for b in ("draftkings", "fanduel"):
        _fail(a, book=b, n=2)
    a.dispatch()
    n_after_dark = notifier.n
    assert a.snapshot()["dark"] is True

    for _ in range(50):                   # still dark, many ticks
        _fail(a, book="draftkings", n=1)
        a.dispatch()
    assert notifier.n == n_after_dark, "dark alert repeated per cycle"

    a.observe(book="draftkings", path="sweep", outcome="ok")
    alerts = a.dispatch()
    kinds = [x.kind for x in alerts]
    assert "consensus_restored" in kinds
    assert a.snapshot()["dark"] is False

    # And it can go dark again.
    _fail(a, book="draftkings", n=2)
    a.dispatch()
    assert a.snapshot()["dark"] is True


def test_first_observation_of_the_process_cannot_fire_dark():
    """Books are observed one at a time inside refresh(). If capacity were
    judged before the fleet is known, the very first ok would read as
    '1 healthy book < 2' and page on startup."""
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=3, min_healthy_books=2)
    a.observe(book="draftkings", path="sweep", outcome="ok")
    assert a.dispatch() == []
    assert notifier.n == 0
    a.observe(book="fanduel", path="sweep", outcome="ok")
    assert a.dispatch() == []
    assert notifier.n == 0


def test_book_degraded_on_any_considered_path_is_unhealthy():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=2, min_healthy_books=2)
    _observe_all(a, ("draftkings", "fanduel"))
    _observe_all(a, ("draftkings", "fanduel"), path="on_demand")
    a.dispatch()
    _fail(a, book="draftkings", path="on_demand", n=2)
    a.dispatch()
    snap = a.snapshot()
    assert snap["healthy_books"] == ["fanduel"]
    assert snap["dark"] is True


# ------------------------------------------------------------------ #
# Path scoping (fix-spec item 4 — #57 flips this in config, not code) #
# ------------------------------------------------------------------ #

def test_unconsidered_path_is_ignored_entirely():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=2, paths=("on_demand",))
    _fail(a, book="draftkings", path="sweep", n=50)
    a.dispatch()
    assert notifier.n == 0
    assert a.snapshot()["streaks"] == {}

    _fail(a, book="draftkings", path="on_demand", n=2)
    a.dispatch()
    assert notifier.n == 1


def test_unconsidered_path_does_not_count_toward_capacity():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=2, min_healthy_books=2,
                 paths=("on_demand",))
    _observe_all(a, SIX_BOOKS, path="sweep")     # sweep-only fleet
    assert a.dispatch() == []
    assert a.snapshot()["healthy_books"] == []


# ------------------------------------------------------------------ #
# The alerter can never break the bot                                 #
# ------------------------------------------------------------------ #

def test_notifier_exception_never_escapes():
    notifier = FakeNotifier(raises=True)
    a = _alerter(notifier, streak_threshold=1)
    _fail(a, n=1)
    a.dispatch()                       # must not raise
    assert notifier.n == 1
    _fail(a, book="novig", n=1)
    a.dispatch()


def test_observe_never_raises_on_junk():
    a = _alerter(FakeNotifier())
    a.observe(book="novig", path="sweep", outcome="not_a_real_outcome")
    a.observe(book=None, path=None, outcome=None)
    a.dispatch()


def test_notifier_runs_off_the_dispatch_thread():
    """The maker ticks 4x/sec. A 10s osascript timeout on the tick thread
    would stall quoting; the notifier does no DB work, so it is safe to
    hand off."""
    started = threading.Event()
    release = threading.Event()
    seen = {"thread": None}

    def slow(title, message):
        seen["thread"] = threading.current_thread()
        started.set()
        release.wait(5)

    a = BookHealthAlerter(label="MM", streak_threshold=1, notifier=slow)
    _fail(a, n=1)
    t0 = time.monotonic()
    a.dispatch()
    elapsed = time.monotonic() - t0
    assert started.wait(5), "notifier never ran"
    assert elapsed < 1.0, f"dispatch blocked {elapsed:.2f}s on the notifier"
    assert seen["thread"] is not threading.current_thread()
    release.set()


def test_concurrent_observes_are_consistent():
    notifier = FakeNotifier()
    a = _alerter(notifier, streak_threshold=100)

    def worker():
        for _ in range(100):
            a.observe(book="novig", path="on_demand",
                      outcome="transport_error", error_class="X")

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    a.dispatch()
    assert a.snapshot()["streaks"]["novig"]["on_demand"] == 800
    assert notifier.n == 1


def test_default_notifier_is_the_macos_one():
    a = BookHealthAlerter(label="MM")
    assert a.notifier is notify_macos


def test_notify_macos_swallows_a_broken_osascript(monkeypatch):
    def boom(*a, **kw):
        raise OSError("no osascript here")
    monkeypatch.setattr("kalshi_common.book_health.subprocess.run", boom)
    notify_macos("t", "m")            # must not raise


# ------------------------------------------------------------------ #
# Service integration — sweep and on-demand                           #
# ------------------------------------------------------------------ #

def test_sweep_transport_failures_reach_the_alerter():
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=3)

    def dead(targets):
        raise BookTransportError("caesars", "events", status_code=403)

    svc = SGPService(books=("caesars",), runners={"caesars": dead},
                     book_health=alerter)
    for _ in range(3):
        svc.refresh([_target()])
        svc.check_book_health()
    assert notifier.n == 1
    assert "BookTransportError:events:403" in notifier.texts
    assert "caesars" in notifier.texts


def test_sweep_success_keeps_the_book_healthy():
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=2, min_healthy_books=1)
    svc = SGPService(books=("caesars",),
                     runners={"caesars": lambda t: [_row()]},
                     book_health=alerter)
    for _ in range(10):
        svc.refresh([_target()])
        svc.check_book_health()
    assert notifier.n == 0
    assert alerter.snapshot()["healthy_books"] == ["caesars"]


def test_sweep_empty_result_keeps_the_book_healthy():
    """An off-day slate returns [] — the #33 fix means that is now genuinely
    'no markets', not 'book down'."""
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=2, min_healthy_books=1)
    svc = SGPService(books=("caesars",), runners={"caesars": lambda t: []},
                     book_health=alerter)
    for _ in range(10):
        svc.refresh([_target()])
        svc.check_book_health()
    assert notifier.n == 0


def test_on_demand_transport_failures_reach_the_alerter():
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=3)

    def boom(client, game):
        raise BookTransportError("novig", "events", status_code=503)

    hooks = _hooks(_resolved(2), lambda c, r, e: 2.0)
    hooks["match_event"] = boom
    svc = SGPService(books=("novig",), on_demand_hooks={"novig": hooks},
                     book_health=alerter)
    for _ in range(3):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None
        svc.check_book_health()
    assert notifier.n == 1
    assert "on_demand" in notifier.texts
    assert "BookTransportError:events:503" in notifier.texts


def test_on_demand_decline_is_not_a_failure():
    """'This book doesn't list that game' is the single most common
    on-demand outcome. If it counted, every RFQ on an unlisted game would
    march the streak toward an alert."""
    notifier = FakeNotifier()
    alerter = _alerter(notifier, streak_threshold=3)
    hooks = _hooks(_resolved(2), lambda c, r, e: 2.0, event=None)
    svc = SGPService(books=("novig",), on_demand_hooks={"novig": hooks},
                     book_health=alerter)
    for _ in range(30):
        assert svc.price_on_demand("novig", GAME, _legs(2)) is None
        svc.check_book_health()
    assert notifier.n == 0


def test_service_without_alerter_is_unchanged():
    """Every existing caller (dashboard CLI, all pre-#37 tests) builds the
    service with no alerter and must be untouched."""
    svc = SGPService(books=("caesars",),
                     runners={"caesars": lambda t: [_row()]})
    assert len(svc.refresh([_target()])["caesars"]) == 1
    assert svc.check_book_health() == []


def test_check_book_health_never_raises():
    class Exploding:
        def dispatch(self):
            raise RuntimeError("boom")

    svc = SGPService(books=("caesars",),
                     runners={"caesars": lambda t: [_row()]})
    svc.book_health = Exploding()
    svc.refresh([_target()])
    assert svc.check_book_health() == []


def test_book_health_is_keyword_only():
    """Repo rule: a new param must never bind positionally (here it would
    land on now_fn and silently break the sweep clock)."""
    sig = inspect.signature(SGPService.__init__)
    assert sig.parameters["book_health"].kind is inspect.Parameter.KEYWORD_ONLY
    with pytest.raises(TypeError):
        sig.bind(None, ("caesars",), 75.0, None, None, None, time.monotonic,
                 object())
