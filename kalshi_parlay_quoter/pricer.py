"""Client for the runner's two leg pricers (#282). Reads only.

GET {RUNNER_URL}/rfq_price.json?rfq_id=...&legs=TICKER:side,TICKER:side
-> {rfq_id, priced_at_ms, legs: [{market_ticker, side, ok, reason, fair,
    fair_lo, fair_hi, bracketed, sigma, c, tau, line, kickoff_ms,
    fair_age_ms, kalshi_status, kalshi_yes_bid, kalshi_yes_ask,
    kalshi_no_ask, detail}]}   (legs in request order; prices 0-1 for the
    leg's side). Only `ok` says a leg may be quoted: a declined leg can still
    carry fair fields, which are logged but never multiplied. The runner
    reads Kalshi's market, book and start for every leg: price() below, used
    to confirm a would-quote parlay.

GET {RUNNER_URL}/rfq_fairs.json?legs=...  (at most FAIRS_CHUNK_LEGS legs)
-> the same leg shape priced off the board alone, no Kalshi read, the
    kalshi_* fields null: fairs() below, used to refresh the leg table.

Any transport error, non-200, or a response whose legs do not line up with
the request raises PricerError; the caller logs it as decline pricer_error.
"""
import http.client
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass


class PricerError(Exception):
    pass


class PricerTimeout(PricerError):
    """The runner answered too slowly (logged pricer_error, detail timeout)."""


class PricerBacklog(PricerError):
    """Not priced because too many RFQs were already waiting for a leg's
    first fair; logged as its own reason, pricer_backlog."""


class LegsNotReady(PricerError):
    """The leg table has no fresh fair for some leg yet."""


class PricerRejected(PricerError):
    """The runner refused the request itself (HTTP 4xx) or answered other
    legs than were asked: something in THIS request is wrong, so retrying it
    unchanged cannot work."""


# The runner's limit per /rfq_fairs.json call (rfqprice.js MAX_FAIR_LEGS).
FAIRS_CHUNK_LEGS = 200


@dataclass
class PriceResult:
    legs: list[dict]
    priced_at_ms: int | None
    elapsed_ms: float


@dataclass
class FairsResult(PriceResult):
    # League label -> when that league's snapshot landed on the runner (ms).
    league_loaded_at: dict[str, int]


def _legs_param(legs: list[tuple[str, str]]) -> str:
    return ",".join(f"{t}:{s}" for t, s in legs)


def build_url(base_url: str, rfq_id: str, legs: list[tuple[str, str]]) -> str:
    q = urllib.parse.urlencode({"rfq_id": rfq_id, "legs": _legs_param(legs)}, safe=":,")
    return f"{base_url}/rfq_price.json?{q}"


def build_fairs_url(base_url: str, legs: list[tuple[str, str]]) -> str:
    q = urllib.parse.urlencode({"legs": _legs_param(legs)}, safe=":,")
    return f"{base_url}/rfq_fairs.json?{q}"


def _get_json(url: str, *, timeout: float, opener) -> dict:
    """GET url as a JSON object. PricerRejected on HTTP 400, PricerTimeout
    on a timeout, PricerError on anything else that is not an answer."""
    try:
        with opener(urllib.request.Request(url, method="GET"),
                    timeout=timeout) as r:
            body = json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 400:    # the request's legs; 403/404/405 are about the runner
            raise PricerRejected(f"http_{e.code}") from e
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
    if not isinstance(body, dict):
        raise PricerError(f"bad_json: expected an object, got {type(body).__name__}")
    return body


def _get_legs(url: str, legs: list[tuple[str, str]], *, timeout: float, opener) -> tuple[PriceResult, dict]:
    """GET url and check the answer lines up with `legs`, leg for leg.
    Returns the legs and the whole answer."""
    t0 = time.monotonic()
    body = _get_json(url, timeout=timeout, opener=opener)
    elapsed_ms = (time.monotonic() - t0) * 1000.0
    out = body.get("legs")
    if not isinstance(out, list) or len(out) != len(legs):
        raise PricerRejected("leg_count_mismatch")
    for (ticker, side), got in zip(legs, out):
        if (not isinstance(got, dict) or got.get("market_ticker") != ticker
                or got.get("side") != side):
            raise PricerRejected(f"leg_order_mismatch at {ticker}")
    return PriceResult(out, body.get("priced_at_ms"), elapsed_ms), body


def _snapshot_times(value, where: str) -> dict[str, int]:
    if not isinstance(value, dict) or not all(isinstance(v, int) for v in value.values()):
        raise PricerError(f"expected {where} as {{league: ms}}, got {str(value)[:80]}")
    return value


def price(base_url: str, rfq_id: str, legs: list[tuple[str, str]], *,
          timeout: float, opener=urllib.request.urlopen) -> PriceResult:
    """Kalshi-checked prices for one RFQ's legs (the confirm call)."""
    return _get_legs(build_url(base_url, rfq_id, legs), legs, timeout=timeout, opener=opener)[0]


def fairs(base_url: str, legs: list[tuple[str, str]], *,
          timeout: float, opener=urllib.request.urlopen) -> FairsResult:
    """Board-only fairs for up to FAIRS_CHUNK_LEGS legs (the table refresh),
    with the snapshot times they were priced on."""
    if len(legs) > FAIRS_CHUNK_LEGS:
        raise ValueError(f"expected at most {FAIRS_CHUNK_LEGS} legs per fairs call, got {len(legs)}")
    result, body = _get_legs(build_fairs_url(base_url, legs), legs, timeout=timeout, opener=opener)
    return FairsResult(result.legs, result.priced_at_ms, result.elapsed_ms,
                       _snapshot_times(body.get("league_loaded_at"), "league_loaded_at"))


def league_snapshot_times(base_url: str, *, timeout: float,
                          opener=urllib.request.urlopen) -> dict[str, int]:
    """When each league's snapshot last landed on the runner, off its
    /health: {league label: ms}. The runner's main thread answers /health
    from what its board thread last posted, so this read costs the board
    thread (where the taker runs) nothing. PricerError when the runner is
    down, says it is not ok, or is a version without the field."""
    body = _get_json(f"{base_url}/health", timeout=timeout, opener=opener)
    if body.get("ok") is not True:
        raise PricerError("the runner's /health is not ok")
    rfq = body.get("rfq")
    return _snapshot_times(rfq.get("leagueLoadedAt") if isinstance(rfq, dict) else None,
                           "/health rfq.leagueLoadedAt")
