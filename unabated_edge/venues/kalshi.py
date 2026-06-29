import sys
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


def best_yes_ask(market_ticker: str) -> float | None:
    body = _orderbook(market_ticker)
    if body is None:
        return None
    # Live API returns fractional levels under orderbook_fp.no_dollars as STRING
    # dollars, e.g. [["0.3500","73972.6"], ...] (verified 2026-06-28). yes_ask =
    # 1 - best (highest) no bid. Fall back to legacy integer-cent orderbook.no.
    fp = (body.get("orderbook_fp") or {}).get("no_dollars")
    if fp:
        return round(1 - max(float(p) for p, _sz in fp), 2)
    no = (body.get("orderbook") or {}).get("no") or []
    if no:
        return round((100 - max(l[0] for l in no)) / 100.0, 2)
    return None


def best_no_ask(market_ticker: str) -> float | None:
    """Cost to BUY one NO contract = 1 - best (highest) yes bid. Used to take the
    'under' side of a total (Kalshi lists Over-only rungs, so under = NO on Over-L)."""
    body = _orderbook(market_ticker)
    if body is None:
        return None
    fp = (body.get("orderbook_fp") or {}).get("yes_dollars")
    if fp:
        return round(1 - max(float(p) for p, _sz in fp), 2)
    yes = (body.get("orderbook") or {}).get("yes") or []
    if yes:
        return round((100 - max(l[0] for l in yes)) / 100.0, 2)
    return None
