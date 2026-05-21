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


def accept_quote(quote_id: str, accepted_side: str) -> tuple[dict | None, dict | str | None]:
    """Accept a quote on the given side ("yes" or "no").

    The Kalshi REST accept endpoint takes only `accepted_side` — sizing happens
    upstream when the RFQ is created via `target_cost_dollars`. Confirmed against
    docs.kalshi.com/openapi.yaml (PUT /communications/quotes/{id}/accept). The
    prior implementation sent `{"contracts": N}` which silently 400'd every
    accept with `invalid_parameters` for the lifetime of the bot.

    Returns (response, error_body):
      - success: ({}, None) — Kalshi returns 204 with no body
      - walked / expired race: (None, error_body) — error_body is whatever
        Kalshi returned (dict on JSON, str otherwise). Lets the caller
        distinguish 'quote_expired' / 'rfq_closed' / etc. in walk diagnostics.

    Raises KalshiAPIError on any other (unexpected) status.
    """
    body = {"accepted_side": accepted_side}
    # Kalshi accept endpoint requires PUT, not POST. POST returns a router-level
    # 404 ("404 page not found" plain text). Verified empirically 2026-05-02.
    status, resp, _ = api("PUT", f"/communications/quotes/{quote_id}/accept", body=body)
    # Success is 204 (empty body) per the OpenAPI spec. We also accept 200/201
    # defensively — Kalshi has historically returned 201 on action endpoints
    # (see create_rfq) and the cost of a silent missed-fill is higher than
    # over-accepting status codes.
    if status in (200, 201, 204):
        return (resp if isinstance(resp, dict) else {}), None
    if status in (400, 404, 409):
        # Common race: quote walked or expired. Also covers genuine bad
        # quote_ids. The caller uses the error body to disambiguate.
        return None, resp
    raise KalshiAPIError(f"accept_quote failed: status={status} body={resp}")


def get_rfq_safe(rfq_id: str) -> dict | None:
    """get_rfq variant that returns None on any failure. Used for post-walk
    diagnostics — must never crash the bot on a follow-up inspection."""
    try:
        return get_rfq(rfq_id)
    except Exception:
        return None


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
