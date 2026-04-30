"""Kalshi RFQ + Quote API client."""

from kalshi_mlb_rfq.auth_client import api


class KalshiAPIError(Exception):
    pass


def mint_combo_ticker(collection_ticker: str, selected_markets: list[dict]) -> tuple[str, str]:
    """Mint (or look up) a combo market by submitting selected_markets.

    Returns (combo_market_ticker, combo_event_ticker).
    Raises KalshiAPIError on non-200.
    """
    path = f"/multivariate_event_collections/{collection_ticker}"
    status, body, _ = api("POST", path, body={"selected_markets": selected_markets})
    if status != 200 or not isinstance(body, dict) or "market_ticker" not in body:
        raise KalshiAPIError(f"mint_combo_ticker failed: status={status} body={body}")
    return body["market_ticker"], body["event_ticker"]


def create_rfq(market_ticker: str, target_cost_dollars: float,
               replace_existing: bool = False) -> str:
    """Create an RFQ. Returns rfq_id.

    NOTE: replace_existing=True does NOT actually replace existing RFQs (recon-confirmed).
    The bot manages its own dedup via combo_cooldown and live_rfqs.
    """
    body = {
        "market_ticker": market_ticker,
        "rest_remainder": False,
        "target_cost_dollars": f"{target_cost_dollars:.2f}",
        "replace_existing": replace_existing,
    }
    status, resp, _ = api("POST", "/communications/rfqs", body=body)
    if status not in (200, 201) or not isinstance(resp, dict) or "id" not in resp:
        raise KalshiAPIError(f"create_rfq failed: status={status} body={resp}")
    return resp["id"]


def get_rfq(rfq_id: str) -> dict:
    """Get an RFQ object. Returns the inner 'rfq' dict, or raises KalshiAPIError."""
    status, body, _ = api("GET", f"/communications/rfqs/{rfq_id}")
    if status != 200 or not isinstance(body, dict) or "rfq" not in body:
        raise KalshiAPIError(f"get_rfq failed: status={status} body={body}")
    return body["rfq"]


def delete_rfq(rfq_id: str) -> bool:
    """DELETE an RFQ. Returns True if cancellation went through, False if already closed."""
    status, body, _ = api("DELETE", f"/communications/rfqs/{rfq_id}")
    if status == 204:
        return True
    # Common case: 400 'expired' = already cancelled or expired. Idempotent no-op.
    if status == 400 and isinstance(body, dict) and body.get("error", {}).get("code") == "expired":
        return False
    raise KalshiAPIError(f"delete_rfq failed: status={status} body={body}")


def poll_quotes(rfq_id: str, user_id: str) -> list[dict]:
    """Poll for quotes on an RFQ. Returns list of quote dicts."""
    path = f"/communications/quotes?rfq_id={rfq_id}&rfq_creator_user_id={user_id}"
    status, body, _ = api("GET", path)
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"poll_quotes failed: status={status} body={body}")
    return body.get("quotes") or []


def accept_quote(quote_id: str, contracts: int) -> dict | None:
    """Accept a quote for the given contract count.

    Returns the accept-response dict on success, or None if the quote walked /
    expired before our accept landed.
    """
    body = {"contracts": contracts}
    status, resp, _ = api("POST", f"/communications/quotes/{quote_id}/accept", body=body)
    # Kalshi proved on create_rfq they may return 201 even on action endpoints.
    # Accept both — the alternative is a silent failure with a real position
    # opened on Kalshi's side and no local tracking.
    if status in (200, 201):
        return resp if isinstance(resp, dict) else {}
    if status in (400, 409):
        # Common race: quote walked or expired. Not a hard error.
        return None
    raise KalshiAPIError(f"accept_quote failed: status={status} body={resp}")


def get_position_contracts(market_ticker: str) -> int:
    """Authoritative current position count for a ticker via /portfolio/positions.

    Returns 0 if no position. Raises KalshiAPIError on API failure.
    """
    status, body, _ = api("GET", f"/portfolio/positions?ticker={market_ticker}&limit=10")
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"get_position_contracts failed: status={status} body={body}")
    for p in body.get("market_positions") or []:
        if p.get("ticker") == market_ticker:
            return int(p.get("position", 0))
    return 0


def list_open_rfqs(user_id: str) -> list[dict]:
    """Return all open RFQs for the given user. Used at startup for phantom cleanup."""
    path = f"/communications/rfqs?status=open&creator_user_id={user_id}&limit=100"
    status, body, _ = api("GET", path)
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"list_open_rfqs failed: status={status} body={body}")
    return body.get("rfqs") or []
