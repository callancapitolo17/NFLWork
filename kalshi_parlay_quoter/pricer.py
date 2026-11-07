"""Client for the runner's leg pricer (#282).

GET {RUNNER_URL}/rfq_price.json?rfq_id=...&legs=TICKER:side,TICKER:side
-> {rfq_id, priced_at_ms, legs: [{market_ticker, side, ok, reason, fair,
    fair_lo, fair_hi, bracketed, sigma, c, tau, line, kickoff_ms,
    fair_age_ms, kalshi_status, kalshi_yes_bid, kalshi_yes_ask,
    kalshi_no_ask}]}   (legs in request order; prices 0-1 for the leg's side)

Any transport error, non-200, or a response whose legs do not line up with
the request raises PricerError; the caller logs it as decline pricer_error.
"""
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass


class PricerError(Exception):
    pass


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
    except (urllib.error.URLError, TimeoutError, OSError) as e:
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
