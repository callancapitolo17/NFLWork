"""Tests for kalshi_common.sgp_service.SGPService.

All tests inject fake `runners` (the test seam) — no HTTP, no clients.
"""
import time
from datetime import datetime, timezone

from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import PricedRow, TargetLine

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)
TARGETS = [TargetLine(game_id="g1", home_team="H", away_team="A",
                      commence_time=CT, period="FG", spread=-1.5, total=8.5)]


def _row(book):
    return PricedRow(game_id="g1", combo="Home Spread + Over", period="FG",
                     spread_line=-1.5, total_line=8.5, bookmaker=book,
                     source=f"{book}_direct", sgp_decimal=3.5,
                     sgp_american=250, fetch_time=CT)


def test_refresh_returns_rows_per_successful_book():
    svc = SGPService(books=("draftkings", "novig"), runners={
        "draftkings": lambda t: [_row("draftkings")],
        "novig": lambda t: [],
    })
    out = svc.refresh(TARGETS)
    assert set(out) == {"draftkings", "novig"}
    assert len(out["draftkings"]) == 1
    assert out["novig"] == []          # success with zero rows != failure


def test_failed_book_is_none_and_does_not_affect_others():
    def boom(t):
        raise RuntimeError("dead session")
    svc = SGPService(books=("draftkings", "fanduel"), runners={
        "draftkings": boom,
        "fanduel": lambda t: [_row("fanduel")],
    })
    out = svc.refresh(TARGETS)
    assert out["draftkings"] is None
    assert len(out["fanduel"]) == 1


def test_client_reset_after_three_consecutive_failures():
    def boom(t):
        raise RuntimeError("dead session")
    svc = SGPService(books=("draftkings",), runners={"draftkings": boom})
    sentinel = object()
    svc._state["draftkings"].client = sentinel
    svc.refresh(TARGETS)
    svc.refresh(TARGETS)
    assert svc._state["draftkings"].client is sentinel   # not yet
    svc.refresh(TARGETS)                                 # 3rd failure
    assert svc._state["draftkings"].client is None       # torn down


def test_success_resets_failure_counter():
    state = {"fail": True}

    def flaky(t):
        if state["fail"]:
            raise RuntimeError("x")
        return [_row("draftkings")]
    svc = SGPService(books=("draftkings",), runners={"draftkings": flaky})
    svc._state["draftkings"].client = object()
    svc.refresh(TARGETS); svc.refresh(TARGETS)       # 2 failures
    state["fail"] = False
    svc.refresh(TARGETS)                             # success
    state["fail"] = True
    svc.refresh(TARGETS)                             # 1 failure (not 3rd)
    assert svc._state["draftkings"].client is not None


def test_per_book_deadline_times_out_slow_book_only():
    def slow(t):
        time.sleep(2.0)
        return [_row("prophetx")]
    svc = SGPService(books=("prophetx", "novig"),
                     per_book_deadline_sec=0.3,
                     runners={"prophetx": slow,
                              "novig": lambda t: [_row("novig")]})
    t0 = time.monotonic()
    out = svc.refresh(TARGETS)
    assert out["prophetx"] is None
    assert len(out["novig"]) == 1
    assert time.monotonic() - t0 < 1.5   # did not wait the full 2s sleep


def test_min_refresh_skips_book_until_due():
    clock = {"t": 1000.0}
    calls = {"n": 0}

    def px_runner(t):
        calls["n"] += 1
        return [_row("prophetx")]
    svc = SGPService(books=("prophetx",),
                     min_refresh_sec={"prophetx": 120},
                     now_fn=lambda: clock["t"],
                     runners={"prophetx": px_runner})
    out1 = svc.refresh(TARGETS)
    assert "prophetx" in out1 and calls["n"] == 1
    clock["t"] += 30                      # 30s later: not due
    out2 = svc.refresh(TARGETS)
    assert "prophetx" not in out2 and calls["n"] == 1
    clock["t"] += 200                     # well past 120s: due again
    out3 = svc.refresh(TARGETS)
    assert "prophetx" in out3 and calls["n"] == 2


def test_timeout_tears_down_client_immediately():
    """A single timeout must drop the book's client right away — the
    worker thread is still running and still holds the (non-thread-safe)
    curl_cffi session, so the next cycle must NOT reuse it."""
    import time as _t

    def slow(t):
        _t.sleep(2.0)
        return []
    svc = SGPService(books=("draftkings",), per_book_deadline_sec=0.2,
                     runners={"draftkings": slow})
    svc._state["draftkings"].client = object()
    out = svc.refresh(TARGETS)
    assert out["draftkings"] is None
    # Torn down after ONE timeout, not after 3 failures.
    assert svc._state["draftkings"].client is None
    assert svc._state["draftkings"].caches == {}
