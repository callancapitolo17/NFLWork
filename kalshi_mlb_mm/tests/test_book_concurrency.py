"""Per-book on-demand concurrency (#101): each book's gate width.

Issue #40 pinned EVERY book to one in-flight pricing call because Novig
403s at ~26 rapid calls. That also capped total throughput near the
second-fastest book's serial rate. #101 makes the width per-book, so the
invariant these tests defend is no longer "one call per book" but:

    a book's simultaneous in-flight price_on_demand calls never exceed
    its configured lanes — and Novig's lanes are 1.
"""
import threading
import time

from kalshi_common import legset
from mlb_sgp._shared import GameRef, OnDemandBookResult
from kalshi_mlb_mm import config
from kalshi_mlb_mm.on_demand import OnDemandEngine

EVT = "KXMLBGAME-25JUN271905NYYBOS"
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)

# Enough overlapping jobs that a widened book WILL exceed one call unless
# its gate stops it; long enough per call that the overlap is not a race.
JOBS = 6
CALL_SEC = 0.15


def _legs(i: int):
    """A distinct leg set per job, so the engine never dedups them."""
    return [legset.CanonicalLeg(EVT, "spread", -1.5, "home"),
            legset.CanonicalLeg(EVT, "total", 7.5 + i, "over")]


class TrackingService:
    """Duck-typed SGPService recording peak concurrency per book."""

    def __init__(self, books):
        self.books = tuple(books)
        self.fair = 0.14
        self.on_demand_deadline_sec = 10.0
        self.lock = threading.Lock()
        self.active: dict[str, int] = {}
        self.max_active: dict[str, int] = {}

    def price_on_demand(self, book, game, legs):
        with self.lock:
            self.active[book] = self.active.get(book, 0) + 1
            self.max_active[book] = max(self.max_active.get(book, 0),
                                        self.active[book])
        time.sleep(CALL_SEC)
        with self.lock:
            self.active[book] -= 1
        return OnDemandBookResult(book=book, fair=self.fair,
                                  route="partition", n_cells_priced=8,
                                  latency_sec=CALL_SEC)


def _await_landing(eng, h, timeout=20.0):
    deadline = time.monotonic() + timeout
    while eng.landed_at(h) is None and time.monotonic() < deadline:
        time.sleep(0.02)
    return eng.landed_at(h) is not None


def _run_flood(books, book_concurrency):
    """Flood the engine with overlapping jobs; return peak concurrency."""
    svc = TrackingService(books)
    eng = OnDemandEngine(svc, book_concurrency=book_concurrency)
    hashes = []
    for i in range(JOBS):
        legs = _legs(i)
        h = legset.leg_set_hash(legs)
        hashes.append(h)
        eng.ensure_fetch(h, GAME, legs)
    assert len(set(hashes)) == JOBS
    for h in hashes:
        assert _await_landing(eng, h), "a job never landed"
    assert svc.max_active, "no calls recorded"
    return svc.max_active


def test_novig_stays_at_one_lane_under_a_job_flood():
    """The #40 guard itself: Novig must never see two calls at once, and
    the config must be what pins it — not an injected test value."""
    assert config.book_concurrency("novig") == 1
    peak = _run_flood(("novig",), book_concurrency=None)
    assert peak["novig"] == 1, peak


def test_widened_book_actually_runs_more_than_one_call_at_once():
    """The point of the ticket: a book given 3 lanes uses them. Without
    this the change is a no-op that still passes every #40 assertion."""
    peak = _run_flood(("fanduel",), book_concurrency={"fanduel": 3})
    assert peak["fanduel"] > 1, peak
    assert peak["fanduel"] <= 3, peak


def test_each_book_is_capped_at_its_own_width_not_a_shared_one():
    """Widths are independent: widening FanDuel must not leak lanes to
    Novig running in the same flight."""
    peak = _run_flood(("fanduel", "novig"),
                      book_concurrency={"fanduel": 3, "novig": 1})
    assert peak["novig"] == 1, peak
    assert peak["fanduel"] > 1, peak


def test_unknown_book_defaults_to_one_lane():
    """A book nobody has measured is never widened by accident."""
    assert config.book_concurrency("some_new_book") == 1
    peak = _run_flood(("some_new_book",), book_concurrency=None)
    assert peak["some_new_book"] == 1, peak


def test_zero_width_clamps_to_one_rather_than_deadlocking():
    """A 0 in config would be read as 'skip this book', but a
    Semaphore(0) never admits anyone: the flight would burn its whole
    deadline waiting instead. Clamp to 1 — skipping a book is the
    deadline-bounded acquire's job, not the width's."""
    peak = _run_flood(("novig",), book_concurrency={"novig": 0})
    assert peak["novig"] == 1, peak


def test_rollback_to_all_ones_restores_pre_101_behaviour(monkeypatch):
    """The named rollback in the ticket: every book back to one lane via
    config alone, no code change."""
    monkeypatch.setattr(config, "ON_DEMAND_BOOK_CONCURRENCY",
                        {b: 1 for b in config.ON_DEMAND_BOOK_CONCURRENCY})
    for book in ("fanduel", "draftkings", "betmgm", "prophetx", "novig",
                 "caesars"):
        assert config.book_concurrency(book) == 1
    peak = _run_flood(("fanduel", "draftkings"), book_concurrency=None)
    assert all(v == 1 for v in peak.values()), peak


# --------------------------------------------------------------------- #
# Env overrides must be loud (pre-merge review finding)                  #
# --------------------------------------------------------------------- #

def test_shipped_defaults_are_not_reported_as_widened():
    """A clean install announces nothing — the warning must mean something."""
    assert config.widened_books() == {}


def test_env_widening_a_book_is_reported(monkeypatch):
    """`ON_DEMAND_CONCURRENCY_NOVIG=3` silently re-opening #40 is exactly
    what this catches: config carries 'MUST remain 1' only as a comment,
    so the override has to announce itself at startup."""
    monkeypatch.setattr(config, "ON_DEMAND_BOOK_CONCURRENCY",
                        {**config.ON_DEMAND_BOOK_CONCURRENCY, "novig": 3})
    assert config.widened_books() == {"novig": (1, 3)}


def test_narrowing_a_book_is_not_reported(monkeypatch):
    """The rollback direction is safe and must stay quiet, or operators
    learn to ignore the warning."""
    monkeypatch.setattr(config, "ON_DEMAND_BOOK_CONCURRENCY",
                        {b: 1 for b in config.ON_DEMAND_BOOK_CONCURRENCY})
    assert config.widened_books() == {}
