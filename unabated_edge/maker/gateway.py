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
        body = {"ticker": ticker, "client_order_id": client_order_id,
                "action": "buy", "side": side, "type": "limit", "count": int(count)}
        body["yes_price" if side == "yes" else "no_price"] = int(price_cents)
        status, resp, _ = auth_client.api("POST", "/portfolio/orders", body)
        if status not in (200, 201) or not isinstance(resp, dict):
            log.warning("maker place failed %s %s %dc x%d: status=%s resp=%s",
                        ticker, side, price_cents, count, status, resp)
            return None
        return (resp.get("order") or {}).get("order_id")

    def cancel(self, order_id):
        status, _, _ = auth_client.api("DELETE", f"/portfolio/orders/{order_id}")
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
