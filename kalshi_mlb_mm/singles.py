"""Kalshi constituent single-leg markets as an independent anchor (issue #23).

Every leg of an MLB combo IS its own 2-way Kalshi market (KXMLBSPREAD-*,
KXMLBTOTAL-*, KXMLBGAME-*), and those books move in REAL TIME while our
sportsbook SGP scrapes refresh every ~150-165s. That makes them the only
pricing input we have that is both independent of the books and fast.

Two consumers, both in main.py:
  * quote path -> `corr_sanity`: is the book-consensus combo fair even
    consistent with the live marginals? (Frechet bounds + correlation premium)
  * risk sweep -> `jumped_tickers`: did a constituent move after we quoted?

plus the confirm-window veto from issue #17 (`leg_market_prices` +
`singles_moved`), which lives here so all three share ONE implementation of
"read a Kalshi single and compare it".

Side effects: `leg_market_prices` and `fetch_market_prices` issue
`GET /markets/{ticker}`. Everything else is pure. No DB access.
"""
import logging
import math
import time
from dataclasses import dataclass

from kalshi_common import auth_client, fair_value

log = logging.getLogger(__name__)


def market_bid_ask(mkt: dict) -> tuple[float, float] | None:
    """(yes_bid, yes_ask) in DOLLARS from either Kalshi market response shape
    — int cents (`yes_bid`) or string-dollar (`yes_bid_dollars`); see the
    kalshi_price_gotchas note on *_dollars vs int shapes."""
    bid, ask = mkt.get("yes_bid"), mkt.get("yes_ask")
    if bid is not None and ask is not None:
        return float(bid) / 100.0, float(ask) / 100.0
    bid, ask = mkt.get("yes_bid_dollars"), mkt.get("yes_ask_dollars")
    if bid is not None and ask is not None:
        return float(bid), float(ask)
    return None


def leg_market_prices(legs: list[dict]) -> dict | None:
    """Raw Kalshi odds for every leg of a combo: {leg market_ticker:
    {"yes_bid": $, "yes_ask": $}}. Each leg IS its own Kalshi singles market,
    so this is one fast GET per leg (~<1s for a 2-3 leg combo). Returns None
    if ANY leg can't be read (fail-safe: no partial baselines)."""
    try:
        out = {}
        for leg in legs:
            mt = str(leg.get("market_ticker") or "")
            if not mt:
                return None
            if mt in out:
                continue
            _status, body, _hdrs = auth_client.api("GET", f"/markets/{mt}")
            mkt = body.get("market") if isinstance(body, dict) else None
            prices = market_bid_ask(mkt) if isinstance(mkt, dict) else None
            if prices is None:
                return None
            out[mt] = {"yes_bid": prices[0], "yes_ask": prices[1]}
        return out or None
    except Exception as e:
        log.warning("[leg_snapshot] fetch failed: %s", e)
        return None


def fetch_market_prices(tickers, budget_sec: float | None = None) -> dict[str, dict]:
    """Best-effort {ticker: {"yes_bid": $, "yes_ask": $}} for a SET of tickers.

    Unlike `leg_market_prices` this does NOT fail closed: unreadable tickers
    are simply absent from the result. That is deliberate and is the risk
    sweep's contract — a transient 500 must not look like a price move and
    flush every resting quote. The fill-blocking backstop is #17's confirm
    veto, which does fail closed.

    Cost: exactly one GET per DISTINCT ticker; callers dedup before calling.

    `budget_sec` stops the loop once that much wall clock is spent. The bot's
    ticks share ONE thread, so an unbounded poll here would delay the confirm
    tick — and Kalshi allows only 2s to confirm in High Volatility Markets.
    Polling fewer tickers is a far cheaper failure than missing a confirm
    window, and partial coverage is safe under the same contract that makes
    unreadable tickers safe. Callers rotate the iteration order so the tail
    is not starved across sweeps.
    """
    out = {}
    deadline = None if budget_sec is None else time.monotonic() + budget_sec
    for mt in tickers:
        if not mt:
            continue
        if deadline is not None and time.monotonic() >= deadline:
            log.warning("[constituent] poll budget %.2fs exhausted after %d ticker(s)",
                        budget_sec, len(out))
            break
        try:
            _status, body, _hdrs = auth_client.api("GET", f"/markets/{mt}")
            mkt = body.get("market") if isinstance(body, dict) else None
            prices = market_bid_ask(mkt) if isinstance(mkt, dict) else None
            if prices is not None:
                out[str(mt)] = {"yes_bid": prices[0], "yes_ask": prices[1]}
        except Exception as e:
            log.warning("[constituent] fetch failed for %s: %s", mt, e)
    return out


def singles_moved(snapshot: dict, fresh: dict) -> bool:
    """True if ANY leg's yes_bid or yes_ask differs between the quote-time
    snapshot and the fresh read. Zero tolerance by design: Kalshi prices move
    in 1c ticks, so 'moved at all' == 'moved >= one tick'. A leg-set mismatch
    counts as moved (fail-safe)."""
    if set(snapshot) != set(fresh):
        return True
    for mt, snap in snapshot.items():
        cur = fresh[mt]
        if (abs(float(snap["yes_bid"]) - float(cur["yes_bid"])) > 1e-9
                or abs(float(snap["yes_ask"]) - float(cur["yes_ask"])) > 1e-9):
            return True
    return False


def devigged_yes(yes_bid, yes_ask) -> float | None:
    """Probit-devigged fair P(YES) for one 2-way Kalshi market, from its
    top-of-book.

    Buying YES costs `yes_ask`; buying NO costs `1 - yes_bid` (the NO ask).
    Those two implied probabilities sum to more than 1 by exactly the quoted
    spread — the 2-way overround `devig_two_way` removes.

    Returns None, never raises, on every degenerate shape: `yes_ask = 1.00`
    (the known divide-by-zero gotcha), `yes_bid = 0`, an empty 0-100 book, or
    a crossed/locked book with no spread to devig. Callers treat None as
    "this leg has no usable anchor".
    """
    try:
        bid, ask = float(yes_bid), float(yes_ask)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(bid) and math.isfinite(ask)):
        return None
    no_ask = 1.0 - bid
    if not (0.0 < ask < 1.0) or not (0.0 < no_ask < 1.0):
        return None
    if bid >= ask:
        return None
    devigged = fair_value.devig_two_way(1.0 / ask, 1.0 / no_ask)
    if devigged is None:
        return None
    fair_yes = float(devigged[0])
    return fair_yes if 0.0 < fair_yes < 1.0 else None


def marginals_for_legs(snapshot: dict, legs: list[dict]) -> list[float] | None:
    """Devigged probability of the side EACH LEG ACTUALLY TAKES, in leg order.

    A leg resting on the NO side contributes 1 - P(YES). Returns None if any
    leg is missing from the snapshot or its book is degenerate: the premium
    gate needs every marginal or it has no baseline at all.
    """
    if not legs or not snapshot:
        return None
    out = []
    for leg in legs:
        prices = snapshot.get(str(leg.get("market_ticker") or ""))
        if not prices:
            return None
        p_yes = devigged_yes(prices.get("yes_bid"), prices.get("yes_ask"))
        if p_yes is None:
            return None
        p = p_yes if str(leg.get("side") or "yes").lower() == "yes" else 1.0 - p_yes
        if not (0.0 < p < 1.0):
            return None
        out.append(p)
    return out


@dataclass(frozen=True)
class CorrSanity:
    """Verdict of the correlation-premium / Frechet check on one combo fair."""
    ok: bool
    reason: str | None            # None | "frechet" | "premium"
    baseline_independent: float
    premium: float
    frechet_lo: float
    frechet_hi: float


def corr_sanity(combo_fair, marginals, premium_min: float,
                premium_max: float) -> CorrSanity | None:
    """Is `combo_fair` consistent with the live single-leg marginals?

    Two independent tests, most-certain first:

    1. FRECHET (parameter-free, always correct). For ANY correlation
       structure whatsoever, a joint probability obeys
       `max(0, sum(p) - (n-1)) <= p_joint <= min(p)`. A violation is not a
       judgement call — it is arithmetic saying our fair cannot be the
       probability of these legs jointly.
    2. PREMIUM BAND (heuristic, tunable). `premium = combo_fair / prod(p)` is
       the correlation multiplier the books are implying. Same-game legs
       legitimately run well above 1 (run line + moneyline ~2x), so the band
       is deliberately wide and ships log-only.

    Returns None when there is nothing to check (no marginals, unusable fair)
    so the caller can log the miss and carry on.
    """
    if not marginals:
        return None
    try:
        fair = float(combo_fair)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(fair) or not (0.0 < fair < 1.0):
        return None
    baseline = 1.0
    for p in marginals:
        baseline *= p
    if baseline <= 0.0:
        return None
    frechet_lo = max(0.0, sum(marginals) - (len(marginals) - 1))
    frechet_hi = min(marginals)
    premium = fair / baseline
    reason = None
    if not (frechet_lo <= fair <= frechet_hi):
        reason = "frechet"
    elif not (premium_min <= premium <= premium_max):
        reason = "premium"
    return CorrSanity(ok=reason is None, reason=reason,
                      baseline_independent=baseline, premium=premium,
                      frechet_lo=frechet_lo, frechet_hi=frechet_hi)


def jumped_tickers(baseline: dict, current: dict,
                   threshold: float) -> dict[str, float]:
    """{ticker: |delta devigged P(YES)|} for tickers that moved MORE than
    `threshold` between the two raw bid/ask snapshots.

    Tickers absent from either snapshot, or degenerate on either side, are
    SKIPPED — absence of signal is not a jump (see `fetch_market_prices`).
    Side is irrelevant here: |delta P(YES)| == |delta P(NO)|.
    """
    out = {}
    if not baseline or not current:
        return out
    for mt, base in baseline.items():
        cur = current.get(mt)
        if not cur or not isinstance(base, dict) or not isinstance(cur, dict):
            continue
        p_base = devigged_yes(base.get("yes_bid"), base.get("yes_ask"))
        p_cur = devigged_yes(cur.get("yes_bid"), cur.get("yes_ask"))
        if p_base is None or p_cur is None:
            continue
        delta = abs(p_cur - p_base)
        if delta > threshold:
            out[str(mt)] = delta
    return out
