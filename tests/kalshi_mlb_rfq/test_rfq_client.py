from unittest.mock import patch

import pytest

from kalshi_mlb_rfq import rfq_client


def test_mint_combo_returns_ticker():
    legs = [
        {"market_ticker": "KXMLBSPREAD-X-NYY2", "event_ticker": "KXMLBSPREAD-X", "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-X-8",     "event_ticker": "KXMLBTOTAL-X",  "side": "yes"},
    ]
    fake_response = (
        200,
        {"market_ticker": "KXMVECROSSCATEGORY-S-FOO", "event_ticker": "KXMVECROSSCATEGORY-S"},
        {},
    )
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake_response) as mock_api:
        ct, et = rfq_client.mint_combo_ticker("KXMVECROSSCATEGORY-R", legs)

    assert ct == "KXMVECROSSCATEGORY-S-FOO"
    assert et == "KXMVECROSSCATEGORY-S"
    call_args = mock_api.call_args
    assert call_args[0][0] == "POST"
    assert call_args[0][1] == "/multivariate_event_collections/KXMVECROSSCATEGORY-R"
    assert call_args[1]["body"] == {"selected_markets": legs}


def test_mint_combo_raises_on_error():
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(400, {"error": "bad"}, {})):
        with pytest.raises(rfq_client.KalshiAPIError):
            rfq_client.mint_combo_ticker("KXMVECROSSCATEGORY-R", [])


def test_create_rfq_returns_id_on_201():
    """Kalshi returns 201 Created on POST /communications/rfqs (REST convention)."""
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(201, {"id": "abc-123"}, {})):
        rid = rfq_client.create_rfq("KXMVECROSSCATEGORY-S-FOO", target_cost_dollars=0.50)
    assert rid == "abc-123"


def test_create_rfq_also_accepts_200():
    """Defensive: still accept 200 in case Kalshi normalizes to it."""
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(200, {"id": "abc-123"}, {})):
        rid = rfq_client.create_rfq("KXMVECROSSCATEGORY-S-FOO", target_cost_dollars=0.50)
    assert rid == "abc-123"


def test_get_rfq_returns_status():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(200, {"rfq": {"id": "abc", "status": "open"}}, {})):
        rfq = rfq_client.get_rfq("abc")
    assert rfq["status"] == "open"


def test_delete_rfq_handles_204():
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(204, {}, {})):
        ok = rfq_client.delete_rfq("abc")
    assert ok is True


def test_delete_rfq_handles_already_closed():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(400, {"error": {"code": "expired"}}, {})):
        ok = rfq_client.delete_rfq("abc")
    assert ok is False


def test_poll_quotes_returns_list():
    fake = (200, {"quotes": [
        {"id": "q1", "yes_bid_dollars": "0.13", "no_bid_dollars": "0.74", "status": "open"},
    ]}, {})
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake) as m:
        quotes = rfq_client.poll_quotes("rfq-1", user_id="user-uuid")

    assert len(quotes) == 1
    assert quotes[0]["id"] == "q1"
    args = m.call_args[0]
    assert "rfq_id=rfq-1" in args[1]
    assert "rfq_creator_user_id=user-uuid" in args[1]


def test_accept_quote_uses_put_verb():
    """Regression guard: accept_quote MUST call PUT, not POST. POST hits a
    router-level 404 ('404 page not found' plain text); PUT is the verb the
    Kalshi midland service routes correctly. Verified empirically 2026-05-02
    by probing both verbs against a fake quote_id."""
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(200, {"order": {"id": "ord-1"}}, {})) as mock_api:
        rfq_client.accept_quote("q1", contracts=10)
    method, path = mock_api.call_args.args[0], mock_api.call_args.args[1]
    assert method == "PUT", f"accept_quote called {method}, must be PUT"
    assert path == "/communications/quotes/q1/accept"


def test_accept_quote_returns_response():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(200, {"order": {"id": "ord-1"}}, {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp["order"]["id"] == "ord-1"


def test_accept_quote_also_accepts_201():
    """Defensive: Kalshi returned 201 on create_rfq despite REST convention
    saying 200 for action endpoints — accept here too to avoid stranded
    positions if they're inconsistent on the accept endpoint."""
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(201, {"order": {"id": "ord-2"}}, {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp["order"]["id"] == "ord-2"


def test_accept_quote_walked_returns_none():
    """400 (quote_walked) and 409 are race conditions — return None, no raise."""
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(400, {"error": {"code": "quote_walked"}}, {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp is None


def test_accept_quote_404_treated_as_walked():
    """If Kalshi returns 404 on a real quote_id (e.g. quote already accepted
    by someone else, or expired between our gate evaluation and the accept
    call), treat it as 'walked' — return None instead of raising. This avoids
    the bot's quote_poll loop logging a noisy KalshiAPIError when the outcome
    was simply that we lost the race."""
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(404, "404 page not found", {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp is None


def test_get_positions_for_combo():
    fake = (200, {"event_positions": [], "market_positions": [
        {"ticker": "KXMVECROSSCATEGORY-S-FOO", "position": 100},
        {"ticker": "OTHER", "position": 5},
    ]}, {})
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake):
        n = rfq_client.get_position_contracts("KXMVECROSSCATEGORY-S-FOO")
    assert n == 100


def test_list_open_rfqs_returns_list():
    fake = (200, {"rfqs": [{"id": "rfq-1"}, {"id": "rfq-2"}]}, {})
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake):
        result = rfq_client.list_open_rfqs("user-uuid")
    assert len(result) == 2
