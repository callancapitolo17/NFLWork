import sys
import time
import logging
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from kalshi_common import auth_client
from unabated_edge import config

log = logging.getLogger("unabated_edge")

def init():
    auth_client.configure(api_key_id=config.KALSHI_API_KEY_ID,
        private_key_path=config.KALSHI_PRIVATE_KEY_PATH, base_url=config.KALSHI_BASE_URL,
        project_root=str(Path(__file__).resolve().parent.parent.parent))

def list_events(series_ticker: str) -> list[dict]:
    status, body, _ = auth_client.api("GET", f"/events?series_ticker={series_ticker}&status=open&with_nested_markets=true")
    if status != 200 or not isinstance(body, dict):
        # A 429/500/auth failure silently returns []; that would zero out the
        # entire Kalshi side with no signal. Make it loud.
        log.warning("Kalshi list_events(%s) failed: status=%s", series_ticker, status)
        return []
    return body.get("events", [])

def _orderbook(market_ticker: str) -> dict | None:
    status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}/orderbook")
    if status != 200 or not isinstance(body, dict):
        return None
    return body


def _parse_levels(fp_levels, legacy_levels) -> list[tuple[float, float]]:
    """One side's bid levels normalized to [(price_dollars, qty), ...] sorted
    best (highest) first. Live API sends orderbook_fp *_dollars as STRING
    dollars, e.g. [["0.3500","73972.6"], ...] (verified 2026-06-28); legacy
    orderbook.* is integer cents."""
    if fp_levels:
        levels = [(float(p), float(q)) for p, q in fp_levels]
    elif legacy_levels:
        levels = [(l[0] / 100.0, float(l[1])) for l in legacy_levels]
    else:
        return []
    return sorted(levels, key=lambda pq: -pq[0])


def get_book(market_ticker: str) -> dict | None:
    """Full normalized orderbook {"yes_bids": [...], "no_bids": [...]} (each a
    best-first [(price, qty), ...] in dollars). None on fetch failure — the
    caller skips the market this tick rather than pricing off a missing book."""
    body = _orderbook(market_ticker)
    if body is None:
        return None
    fp = body.get("orderbook_fp") or {}
    legacy = body.get("orderbook") or {}
    return {
        "yes_bids": _parse_levels(fp.get("yes_dollars"), legacy.get("yes")),
        "no_bids": _parse_levels(fp.get("no_dollars"), legacy.get("no")),
    }


def yes_ask_from_book(book: dict | None) -> float | None:
    """Cost to BUY one YES contract = 1 - best (highest) no bid."""
    if book and book["no_bids"]:
        return round(1 - book["no_bids"][0][0], 2)
    return None


def no_ask_from_book(book: dict | None) -> float | None:
    """Cost to BUY one NO contract = 1 - best (highest) yes bid. Used to take the
    'under' side of a total (Kalshi lists Over-only rungs, so under = NO on Over-L)."""
    if book and book["yes_bids"]:
        return round(1 - book["yes_bids"][0][0], 2)
    return None


def best_yes_ask(market_ticker: str) -> float | None:
    return yes_ask_from_book(get_book(market_ticker))


def best_no_ask(market_ticker: str) -> float | None:
    return no_ask_from_book(get_book(market_ticker))


def recent_trades(market_ticker: str, lookback_sec: int) -> list[dict]:
    """Executed-trades tape for one market since now-lookback. Poll windows
    overlap on purpose — rows dedup on trade_id at insert. limit=100 caps a
    single poll; a >100-trade burst inside one window loses the excess (fine
    for dry-run capture, noted in README)."""
    min_ts = int(time.time()) - lookback_sec
    status, body, _ = auth_client.api(
        "GET", f"/markets/trades?ticker={market_ticker}&min_ts={min_ts}&limit=100")
    if status != 200 or not isinstance(body, dict):
        log.warning("Kalshi trades(%s) failed: status=%s", market_ticker, status)
        return []
    out = []
    for t in body.get("trades") or []:
        tid = t.get("trade_id")
        if not tid:
            continue
        # same *_dollars-vs-integer-cents duality as the orderbook payload
        yes_price = t.get("yes_price_dollars")
        if yes_price is not None:
            yes_price = float(yes_price)
        elif t.get("yes_price") is not None:
            yes_price = t["yes_price"] / 100.0
        count = t.get("count_fp") if t.get("count_fp") is not None else t.get("count")
        out.append({
            "trade_id": tid,
            "market_ticker": t.get("ticker") or market_ticker,
            "created_time": t.get("created_time"),
            "yes_price": yes_price,
            "count": float(count) if count is not None else None,
            "taker_side": t.get("taker_side"),
        })
    return out
