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
    """The runner answered too slowly. Not a reason to back off: a slow
    call with many uncached legs (each a Kalshi read on the runner) says
    nothing about the next call."""


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
    except TimeoutError as e:
        raise PricerTimeout("timeout") from e
    except urllib.error.URLError as e:
        if isinstance(e.reason, TimeoutError):
            raise PricerTimeout("timeout") from e
        raise PricerError(f"transport: {e}") from e
    except OSError as e:
        raise PricerError(f"transport: {e}") from e
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
    """Not priced because the pricer queue was full (logged as its own
    decline reason, pricer_backlog, so it is never mistaken for a runner
    failure)."""


class CachingPricer:
    """Leg-level cache + failure backoff in front of price().

    Measured 2026-11-07 on the VM: ~100 in-league RFQs/s, mostly the same
    legs re-asked. Each (ticker, side) is priced at most once per ttl_sec;
    a request sends the runner only the legs not in cache. After a runner
    failure, calls fail fast for backoff_sec instead of queueing more.
    """

    def __init__(self, inner, *, ttl_sec: float, backoff_sec: float,
                 clock=time.monotonic):
        self._inner = inner            # inner(rfq_id, legs) -> PriceResult
        self._ttl = ttl_sec
        self._backoff = backoff_sec
        self._clock = clock
        self._cache: dict[tuple[str, str], tuple[float, dict, int | None]] = {}
        self._backoff_until = 0.0
        self._last_error = ""
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
        elapsed, priced_at = 0.0, None
        if missing:
            try:
                res = self._inner(rfq_id, missing)
            except PricerTimeout:
                raise
            except PricerError as e:
                with self._lock:
                    self._backoff_until = self._clock() + self._backoff
                    self._last_error = str(e)[:200]
                raise
            elapsed, priced_at = res.elapsed_ms, res.priced_at_ms
            expiry = self._clock() + self._ttl
            fresh = {k: (expiry, leg, priced_at) for k, leg in zip(missing, res.legs)}
            with self._lock:
                self.runner_calls += 1
                self._cache.update(fresh)
            hits.update(fresh)
        with self._lock:
            self.cache_hits += len(legs) - len(missing)
        if priced_at is None:
            priced_at = min((v[2] for v in hits.values() if v[2] is not None), default=None)
        return PriceResult([dict(hits[k][1]) for k in legs], priced_at, elapsed)
