"""Client for the runner's leg pricer (#282).

GET {RUNNER_URL}/rfq_price.json?rfq_id=...&legs=TICKER:side,TICKER:side
-> {rfq_id, priced_at_ms, legs: [{market_ticker, side, ok, reason, fair,
    fair_lo, fair_hi, bracketed, sigma, c, tau, line, kickoff_ms,
    fair_age_ms, kalshi_status, kalshi_yes_bid, kalshi_yes_ask,
    kalshi_no_ask, detail}]}   (legs in request order; prices 0-1 for the
    leg's side). Only `ok` says a leg may be quoted: a declined leg can still
    carry fair fields, which are logged but never multiplied.

Any transport error, non-200, or a response whose legs do not line up with
the request raises PricerError; the caller logs it as decline pricer_error.
"""
import http.client
import json
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass


class PricerError(Exception):
    pass


class PricerTimeout(PricerError):
    """The runner answered too slowly. One slow call (many uncached legs,
    each a Kalshi read on the runner) says little about the next, so only
    TIMEOUTS_BEFORE_BACKOFF in a row back off."""


@dataclass
class PriceResult:
    legs: list[dict]
    priced_at_ms: int | None
    elapsed_ms: float


def build_url(base_url: str, rfq_id: str, legs: list[tuple[str, str]]) -> str:
    q = urllib.parse.urlencode({
        "rfq_id": rfq_id,
        "legs": ",".join(f"{t}:{s}" for t, s in legs),
    }, safe=":,")
    return f"{base_url}/rfq_price.json?{q}"


def price(base_url: str, rfq_id: str, legs: list[tuple[str, str]], *,
          timeout: float, opener=urllib.request.urlopen) -> PriceResult:
    url = build_url(base_url, rfq_id, legs)
    t0 = time.monotonic()
    try:
        with opener(urllib.request.Request(url, method="GET"),
                    timeout=timeout) as r:
            body = json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        raise PricerError(f"http_{e.code}") from e
    except OSError as e:     # URLError and TimeoutError are OSErrors
        reason = getattr(e, "reason", e)
        if isinstance(e, TimeoutError) or isinstance(reason, TimeoutError):
            raise PricerTimeout("timeout") from e
        raise PricerError(f"transport: {e}") from e
    except http.client.HTTPException as e:   # e.g. IncompleteRead mid-body
        raise PricerError(f"transport: {type(e).__name__}") from e
    except (ValueError, UnicodeDecodeError) as e:
        raise PricerError(f"bad_json: {e}") from e
    elapsed_ms = (time.monotonic() - t0) * 1000.0
    out = body.get("legs") if isinstance(body, dict) else None
    if not isinstance(out, list) or len(out) != len(legs):
        raise PricerError("leg_count_mismatch")
    for (ticker, _side), got in zip(legs, out):
        if not isinstance(got, dict) or got.get("market_ticker") != ticker:
            raise PricerError(f"leg_order_mismatch at {ticker}")
    return PriceResult(out, body.get("priced_at_ms"), elapsed_ms)


class PricerBacklog(PricerError):
    """Not priced because the pricer queue was full; logged as its own
    reason, pricer_backlog. A runner stuck on timeouts backs off
    (CachingPricer), so its RFQs log pricer_error, not this."""


TIMEOUTS_BEFORE_BACKOFF = 3


class CachingPricer:
    """Leg-level cache + failure backoff in front of price().

    ~100 in-league RFQs/s re-ask the same legs (VM, 2026-11-07). Each
    (ticker, side) the runner priced ok is reused for ttl_sec, counted from
    when the call was SENT; a request asks the runner only for the rest.
    Declined legs are never cached (a one-off kalshi_unread must not decline
    every RFQ for ttl_sec). Every returned leg carries its own priced_at_ms
    and cache_age_ms, and fair_age_ms grows by the cache age, so the log
    never shows a reused price as fresh. A runner error, or
    TIMEOUTS_BEFORE_BACKOFF timeouts in a row, fails calls fast for
    backoff_sec.
    """

    def __init__(self, inner, *, ttl_sec: float, backoff_sec: float,
                 clock=time.monotonic):
        self._inner = inner            # inner(rfq_id, legs) -> PriceResult
        self._ttl = ttl_sec
        self._backoff = backoff_sec
        self._clock = clock
        # (ticker, side) -> (expiry, leg, runner priced_at_ms, stored_at)
        self._cache: dict[tuple[str, str], tuple[float, dict, int | None, float]] = {}
        self._backoff_until = 0.0
        self._last_error = ""
        self._timeouts_in_a_row = 0
        self._lock = threading.Lock()
        self.runner_calls = 0
        self.cache_hits = 0

    def __call__(self, rfq_id: str, legs: list[tuple[str, str]]) -> PriceResult:
        now = self._clock()
        with self._lock:
            if now < self._backoff_until:
                raise PricerError(f"backoff after {self._last_error}")
            hits = {k: v for k in legs
                    if (v := self._cache.get(k)) and v[0] > now}
            if len(self._cache) > 50_000:
                self._cache = {k: v for k, v in self._cache.items() if v[0] > now}
        missing = list(dict.fromkeys(k for k in legs if k not in hits))
        elapsed = 0.0
        fresh: dict = {}
        if missing:
            sent_at = self._clock()
            try:
                res = self._inner(rfq_id, missing)
            except PricerError as e:
                self._note_failure(e)
                raise
            elapsed = res.elapsed_ms
            expiry = sent_at + self._ttl
            fresh = {k: (expiry, leg, res.priced_at_ms, sent_at)
                     for k, leg in zip(missing, res.legs)}
            with self._lock:
                self.runner_calls += 1
                self._timeouts_in_a_row = 0
                self._cache.update({k: v for k, v in fresh.items() if v[1].get("ok")})
        with self._lock:
            self.cache_hits += len(legs) - len(missing)
        entries = {**hits, **fresh}
        now = self._clock()
        out = [self._aged_leg(*entries[k], now=now) for k in legs]
        stamps = [leg["priced_at_ms"] for leg in out if leg["priced_at_ms"] is not None]
        return PriceResult(out, min(stamps) if stamps else None, elapsed)

    @staticmethod
    def _aged_leg(_expiry, leg: dict, priced_at_ms, stored_at: float, *, now: float) -> dict:
        """A copy of a cached leg stamped with its true age."""
        age_ms = max(0.0, (now - stored_at) * 1000.0)
        aged = dict(leg)
        aged["priced_at_ms"] = priced_at_ms
        aged["cache_age_ms"] = age_ms
        if isinstance(leg.get("fair_age_ms"), (int, float)):
            aged["fair_age_ms"] = leg["fair_age_ms"] + age_ms
        return aged

    def _note_failure(self, e: PricerError) -> None:
        with self._lock:
            if isinstance(e, PricerTimeout):
                self._timeouts_in_a_row += 1
                if self._timeouts_in_a_row < TIMEOUTS_BEFORE_BACKOFF:
                    return
            self._backoff_until = self._clock() + self._backoff
            self._last_error = str(e)[:200]
