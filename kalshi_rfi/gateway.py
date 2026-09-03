"""Order-placement seam: shadow (default) vs live, YES bids only.

Minimal one-sided sibling of unabated_edge/maker/gateway.py — same v2
endpoint, same dead-man switch contract. Kept local rather than imported so
the two live bots stay decoupled (and this one can only ever place YES bids).
"""
import logging
import urllib.error
from abc import ABC, abstractmethod

from kalshi_common import auth_client

log = logging.getLogger(__name__)

# auth_client.api catches only HTTPError; connection resets, DNS failures,
# and timeouts raise through (adversarial review finding 2). The gateway is
# the containment boundary: a network exception must degrade to "this call
# failed", never crash the trading loop or the shutdown cancel sweep.
NETWORK_ERRORS = (urllib.error.URLError, OSError, TimeoutError)


class OrderGateway(ABC):
    is_live = False

    @abstractmethod
    def place_yes_bid(self, ticker: str, price_cents: int, count: int,
                      client_order_id: str,
                      exchange_index: int | None = None) -> str | None: ...

    @abstractmethod
    def cancel(self, order_id: str, *, ticker: str | None = None,
               exchange_index: int | None = None) -> bool: ...


class ShadowGateway(OrderGateway):
    """Places nothing; fabricates order ids. Every decision is logged to the
    DB regardless of gateway, so shadow mode is the free dataset."""

    def __init__(self):
        self._n = 0

    def place_yes_bid(self, ticker, price_cents, count, client_order_id,
                      exchange_index=None):
        self._n += 1
        return f"shadow-{self._n}"

    def cancel(self, order_id, *, ticker=None, exchange_index=None):
        return True


class LiveGateway(OrderGateway):
    is_live = True

    def place_yes_bid(self, ticker, price_cents, count, client_order_id,
                      exchange_index=None):
        # Kalshi v2 order API (v1 POST /portfolio/orders retired -> HTTP 410).
        # A YES buy is a bid at the YES price, in decimal-dollar strings.
        body = {"ticker": ticker, "side": "bid",
                "price": f"{price_cents / 100.0:.4f}",
                "count": f"{float(count):.2f}",
                "time_in_force": "good_till_canceled",
                "self_trade_prevention_type": "taker_at_cross",
                "client_order_id": client_order_id}
        if exchange_index is not None:
            # Exchange sharding (announced 2026-08-24): the ticker alone
            # auto-routes, but naming the shard skips that lookup's latency.
            # Never hardcode the index — baseball is 3 today, NFL is 0.
            body["exchange_index"] = int(exchange_index)
        try:
            status, resp, _ = auth_client.api(
                "POST", "/portfolio/events/orders", body)
        except NETWORK_ERRORS as e:
            # The POST may have landed despite the lost response — that
            # orphan is invisible to local state, which is exactly what the
            # periodic orphan sweep (ORPHAN_SWEEP_SEC) exists to cancel.
            log.error("place NETWORK FAILURE %s bid %dc x%d: %s",
                      ticker, price_cents, count, e)
            return None
        if status not in (200, 201) or not isinstance(resp, dict):
            log.warning("place failed %s bid %dc x%d: status=%s resp=%s",
                        ticker, price_cents, count, status, resp)
            return None
        return resp.get("order_id") or (resp.get("order") or {}).get("order_id")

    def cancel(self, order_id, *, ticker=None, exchange_index=None):
        """Cancel one resting order. Routing is not optional: DELETE carries
        no ticker in its path, so without `exchange_index` (or the
        `market_ticker` auto-route hint) Kalshi routes to shard 0 and
        returns 404 for a baseball order that lives on shard 3 — which the
        404-means-gone rule below then reported as a successful cancel while
        the order kept resting (bug 2, observed live 2026-08-27)."""
        path = f"/portfolio/events/orders/{order_id}"
        params = []
        if exchange_index is not None:
            params.append(f"exchange_index={int(exchange_index)}")
        if ticker:
            params.append(f"market_ticker={ticker}")
        if params:
            path = f"{path}?{'&'.join(params)}"
        try:
            status, _, _ = auth_client.api("DELETE", path)
        except NETWORK_ERRORS as e:
            log.error("cancel NETWORK FAILURE %s: %s", order_id, e)
            return False
        if status == 404:
            # A 404 is ambiguous: the order may be terminal (filled,
            # expired, auto-cancelled at close) or merely misrouted. Ask
            # Kalshi which it is — treating every 404 as "already gone" is
            # what silently stranded 13 live orders. Unknown fails CLOSED:
            # local state is held and the next cycle retries, so the phantom
            # entry the MM bug taught us about can only outlive a listing
            # endpoint that is itself down.
            still = self._order_still_resting(order_id)
            if still is None:
                log.error("cancel %s: 404 and could not verify — holding",
                          order_id)
                return False
            if still:
                log.error("cancel %s: 404 but STILL RESTING on Kalshi "
                          "(routing? exchange_index=%s) — holding",
                          order_id, exchange_index)
                return False
            log.info("cancel %s: verified gone (404) — treating as cancelled",
                     order_id)
            return True
        if status not in (200, 204):
            log.warning("cancel failed %s: status=%s", order_id, status)
            return False
        return True

    def _order_still_resting(self, order_id: str) -> bool | None:
        """True/False if Kalshi's resting listing settles it, None if the
        check itself failed. Omitting exchange_index lists ALL shards."""
        try:
            status, body, _ = auth_client.api(
                "GET", "/portfolio/orders?status=resting&limit=1000")
        except NETWORK_ERRORS as e:
            log.error("cancel verify network failure %s: %s", order_id, e)
            return None
        if status != 200 or not isinstance(body, dict):
            return None
        return any(o.get("order_id") == order_id
                   for o in body.get("orders") or [])


def make_gateway(mode: str | None, ack: str | None) -> OrderGateway | None:
    if not mode or mode == "off":
        return None
    if mode == "shadow":
        return ShadowGateway()
    if mode == "live":
        if ack != "1":
            raise SystemExit("RFI_MODE=live requires RFI_LIVE_ACK=1 (dead-man switch)")
        return LiveGateway()
    raise SystemExit(f"unknown RFI_MODE={mode!r}")
