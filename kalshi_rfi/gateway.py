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
                      client_order_id: str) -> str | None: ...

    @abstractmethod
    def cancel(self, order_id: str) -> bool: ...


class ShadowGateway(OrderGateway):
    """Places nothing; fabricates order ids. Every decision is logged to the
    DB regardless of gateway, so shadow mode is the free dataset."""

    def __init__(self):
        self._n = 0

    def place_yes_bid(self, ticker, price_cents, count, client_order_id):
        self._n += 1
        return f"shadow-{self._n}"

    def cancel(self, order_id):
        return True


class LiveGateway(OrderGateway):
    is_live = True

    def place_yes_bid(self, ticker, price_cents, count, client_order_id):
        # Kalshi v2 order API (v1 POST /portfolio/orders retired -> HTTP 410).
        # A YES buy is a bid at the YES price, in decimal-dollar strings.
        body = {"ticker": ticker, "side": "bid",
                "price": f"{price_cents / 100.0:.4f}",
                "count": f"{float(count):.2f}",
                "time_in_force": "good_till_canceled",
                "self_trade_prevention_type": "taker_at_cross",
                "client_order_id": client_order_id}
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

    def cancel(self, order_id):
        try:
            status, _, _ = auth_client.api(
                "DELETE", f"/portfolio/events/orders/{order_id}")
        except NETWORK_ERRORS as e:
            log.error("cancel NETWORK FAILURE %s: %s", order_id, e)
            return False
        if status == 404:
            # Explicit 404 = the order is already gone (filled out, expired,
            # or auto-cancelled on market close). Treating it as failure
            # wedges the ticker with a phantom resting entry forever — the
            # exact MM phantom-open-quotes bug (resolved 2026-08-12).
            log.info("cancel %s: already gone (404) — treating as cancelled",
                     order_id)
            return True
        if status not in (200, 204):
            log.warning("cancel failed %s: status=%s", order_id, status)
            return False
        return True


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
