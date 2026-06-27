"""RFQ read-path: discover open RFQs + fetch a combo market's legs. v1 = REST poll;
a WsRFQSource can replace this behind the same interface later."""
from typing import Protocol

from kalshi_common import auth_client


class RFQSource(Protocol):
    def poll(self) -> list[dict]: ...
    def get_market(self, market_ticker: str) -> dict | None: ...


class RestRFQSource:
    def poll(self) -> list[dict]:
        # No creator filter → market-wide open list (others' RFQs we can quote on).
        status, body, _ = auth_client.api(
            "GET", "/communications/rfqs?status=open&limit=100")
        if status != 200 or not isinstance(body, dict):
            return []
        return body.get("rfqs") or []

    def get_market(self, market_ticker: str) -> dict | None:
        status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}")
        if status != 200 or not isinstance(body, dict):
            return None
        return body.get("market")
