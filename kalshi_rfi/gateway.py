"""Order-placement seam: shadow (default) vs live, YES bids only.

Minimal one-sided sibling of unabated_edge/maker/gateway.py — same v2
endpoint, same dead-man switch contract. Kept local rather than imported so
the two live bots stay decoupled (and this one can only ever place YES bids).
"""
import logging
from abc import ABC, abstractmethod

from kalshi_common import auth_client

log = logging.getLogger(__name__)


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
        status, resp, _ = auth_client.api("POST", "/portfolio/events/orders", body)
        if status not in (200, 201) or not isinstance(resp, dict):
            log.warning("place failed %s bid %dc x%d: status=%s resp=%s",
                        ticker, price_cents, count, status, resp)
            return None
        return resp.get("order_id") or (resp.get("order") or {}).get("order_id")

    def cancel(self, order_id):
        status, _, _ = auth_client.api(
            "DELETE", f"/portfolio/events/orders/{order_id}")
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
