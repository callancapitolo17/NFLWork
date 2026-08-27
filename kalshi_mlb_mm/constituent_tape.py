"""A short memory of Kalshi constituent single-leg prices (issue #99).

The leg surface (#96/#98) prices a cross-game combo from CACHED book fairs. A
row is only as good as the market it was built against, so before quoting we
ask the one real-time, book-independent source we have — Kalshi's own
constituent single-leg markets — whether the market moved AFTER the row was
built. Answering that needs a PAST Kalshi price, and the answer must cost
nothing: the maker already reads these markets three times over
(``main._leg_market_prices`` at quote time for #17/#23, the same call again in
the confirm last look, and ``singles.fetch_market_prices`` in the #23 risk
sweep). This module simply REMEMBERS those reads.

Inputs: raw ``{ticker: {"yes_bid": $, "yes_ask": $}}`` maps the bot already
holds. Outputs: the devigged P(YES) most recently observed at or before a
given instant. Side effects: NONE — no network, no DB, no API calls. Adding a
Kalshi request here would defeat the whole point.

Threading: every feed site runs on the maker's single tick thread, but the
lock is kept anyway so a future feeder (a WS price stream) cannot corrupt the
deques silently.
"""
from __future__ import annotations

import threading
from bisect import bisect_right
from collections import deque
from datetime import datetime, timedelta, timezone

from kalshi_mlb_mm.singles import devigged_yes


class ConstituentTape:
    """{ticker: [(observed_at, devigged P(YES))]}, newest last, bounded twice.

    Bounded by BOTH a retention window and a per-ticker point cap: the window
    is what the veto actually needs, and the cap is the backstop for a process
    that runs for days while one ticker is polled far more often than the
    window's worth of points would suggest.
    """

    def __init__(self, retention_sec: float, max_points: int):
        self._retention = timedelta(seconds=max(float(retention_sec), 0.0))
        self._max_points = max(int(max_points), 1)
        self._lock = threading.Lock()
        self._points: dict[str, deque] = {}

    # ---------------------------------------------------------------- #
    # write
    # ---------------------------------------------------------------- #
    def record(self, prices: dict | None, observed_at: datetime) -> int:
        """Remember one already-fetched batch of raw Kalshi bid/asks.

        Degenerate books (``yes_ask`` = 1.00, empty, crossed) devig to None
        and are simply not recorded — an unusable price is not an
        observation, and recording it would let the veto compare against a
        number ``devigged_yes`` itself refuses to produce.

        Returns how many tickers were recorded. Never raises: this is fed from
        inside the trading ticks, and a memory aid must never break one.
        """
        if not prices:
            return 0
        stamped = _as_aware(observed_at)
        n = 0
        with self._lock:
            for ticker, raw in prices.items():
                if not isinstance(raw, dict):
                    continue
                p = devigged_yes(raw.get("yes_bid"), raw.get("yes_ask"))
                if p is None:
                    continue
                points = self._points.get(str(ticker))
                if points is None:
                    points = self._points[str(ticker)] = deque(
                        maxlen=self._max_points)
                # Out-of-order arrivals would break the bisect in
                # price_at_or_before, and there is no useful ordering to
                # restore inside a deque — drop them instead. In practice the
                # feeders are one thread reading a monotonically later clock.
                if points and stamped < points[-1][0]:
                    continue
                points.append((stamped, p))
                n += 1
            self._prune(stamped)
        return n

    def _prune(self, now: datetime) -> None:
        """Drop points outside the retention window (caller holds the lock)."""
        cutoff = now - self._retention
        for ticker in list(self._points):
            points = self._points[ticker]
            while points and points[0][0] < cutoff:
                points.popleft()
            if not points:
                del self._points[ticker]

    # ---------------------------------------------------------------- #
    # read
    # ---------------------------------------------------------------- #
    def price_at_or_before(self, ticker: str,
                           when: datetime) -> tuple[float, datetime] | None:
        """(devigged P(YES), observed_at) for the NEWEST observation at or
        before `when`, or None if the tape never saw this ticker that far back.

        None means "no signal", not "no move" — the caller must fail OPEN and
        count it, exactly as ``singles.jumped_tickers`` treats an unreadable
        ticker. Anything else would turn a cold tape into a decline storm.
        """
        target = _as_aware(when)
        with self._lock:
            points = self._points.get(str(ticker))
            if not points:
                return None
            snapshot = list(points)
        idx = bisect_right([p[0] for p in snapshot], target)
        if idx == 0:
            return None
        observed_at, price = snapshot[idx - 1]
        return price, observed_at

    def ticker_count(self) -> int:
        with self._lock:
            return len(self._points)

    def point_count(self) -> int:
        with self._lock:
            return sum(len(p) for p in self._points.values())


def _as_aware(ts: datetime) -> datetime:
    """Naive timestamps are treated as session-local, mirroring
    ``risk._now_matching`` — the surface stamps ``built_at`` aware, but a
    caller passing a naive clock must not silently compare against UTC."""
    return ts if ts.tzinfo is not None else ts.astimezone(timezone.utc)
