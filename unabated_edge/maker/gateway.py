"""Order-placement seam. Engine talks only to QuoteGateway; swapping
shadow<->live changes no engine code (interface + implementations —
the kalshi_mlb_mm pattern)."""
import logging
from abc import ABC, abstractmethod

from kalshi_common import auth_client

log = logging.getLogger("unabated_edge")


class QuoteGateway(ABC):
    is_live = False

    @abstractmethod
    def place(self, ticker: str, side: str, price_cents: int, count: int,
              client_order_id: str) -> str | None: ...

    @abstractmethod
    def cancel(self, order_id: str) -> bool: ...


class ShadowGateway(QuoteGateway):
    """Places nothing; fabricates order ids. The engine logs every decision to
    the maker DB regardless of gateway, so shadow mode costs no code."""

    def __init__(self):
        self._n = 0

    def place(self, ticker, side, price_cents, count, client_order_id):
        self._n += 1
        return f"shadow-{self._n}"

    def cancel(self, order_id):
        return True


class LiveGateway(QuoteGateway):
    is_live = True

    def place(self, ticker, side, price_cents, count, client_order_id):
        # Kalshi v2 order API (v1 POST /portfolio/orders was retired -> HTTP 410).
        # side is bid/ask on the YES leg only: buy YES = bid @ price;
        # buy NO = sell YES (ask) @ (1 - price). price/count are decimal-DOLLAR
        # strings, not integer cents. Our logical price_cents is the price we
        # want to PAY on our logical side, so the ask price is the YES inverse.
        if side == "yes":
            v2_side, price_dollars = "bid", price_cents / 100.0
        else:
            v2_side, price_dollars = "ask", (100 - price_cents) / 100.0
        body = {"ticker": ticker, "side": v2_side,
                "price": f"{price_dollars:.4f}", "count": f"{float(count):.2f}",
                "time_in_force": "good_till_canceled",
                "self_trade_prevention_type": "taker_at_cross",
                "client_order_id": client_order_id}
        status, resp, _ = auth_client.api("POST", "/portfolio/events/orders", body)
        if status not in (200, 201) or not isinstance(resp, dict):
            log.warning("maker place failed %s %s(%s) $%.4f x%d: status=%s resp=%s",
                        ticker, side, v2_side, price_dollars, count, status, resp)
            return None
        # v2 returns order_id at top level; keep the nested fallback defensively.
        return resp.get("order_id") or (resp.get("order") or {}).get("order_id")

    def cancel(self, order_id):
        status, _, _ = auth_client.api("DELETE", f"/portfolio/events/orders/{order_id}")
        if status not in (200, 204):
            log.warning("maker cancel failed %s: status=%s", order_id, status)
            return False
        return True


def make_gateway(mode: str | None, ack: str | None) -> QuoteGateway | None:
    if not mode or mode == "off":
        return None
    if mode == "shadow":
        return ShadowGateway()
    if mode == "live":
        if ack != "1":
            raise SystemExit("MAKER_MODE=live requires MAKER_LIVE_ACK=1 (dead-man switch)")
        return LiveGateway()
    raise SystemExit(f"unknown MAKER_MODE={mode!r}")
