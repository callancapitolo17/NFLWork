"""Quote write-path: submit / confirm / cancel. The pricing+risk code only ever
touches this interface, never raw HTTP — so a WS/FIX gateway can drop in later."""
from typing import Protocol

from kalshi_common import auth_client


class QuoteGateway(Protocol):
    def submit_quote(self, rfq_id: str, yes_bid: float, no_bid: float) -> str | None: ...
    def confirm(self, quote_id: str) -> bool: ...
    def cancel(self, quote_id: str) -> bool: ...


class RestQuoteGateway:
    def submit_quote(self, rfq_id: str, yes_bid: float, no_bid: float) -> str | None:
        body = {"rfq_id": rfq_id, "yes_bid": f"{yes_bid:.4f}",
                "no_bid": f"{no_bid:.4f}", "rest_remainder": False}
        status, resp, _ = auth_client.api("POST", "/communications/quotes", body=body)
        if status in (200, 201) and isinstance(resp, dict) and "id" in resp:
            return resp["id"]
        return None

    def confirm(self, quote_id: str) -> bool:
        # PUT /communications/quotes/{id}/confirm — no body per the OpenAPI spec.
        # Open item: verify the exact response/body on the first real accept.
        status, _, _ = auth_client.api("PUT", f"/communications/quotes/{quote_id}/confirm")
        return status in (200, 201, 204)

    def cancel(self, quote_id: str) -> bool:
        status, _, _ = auth_client.api("DELETE", f"/communications/quotes/{quote_id}")
        return status in (200, 204)
