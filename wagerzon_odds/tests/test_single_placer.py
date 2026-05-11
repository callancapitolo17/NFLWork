"""Tests for wagerzon_odds.single_placer

Mocks the WZ HTTP session and verifies:
- sel string for a 1-leg payload matches parlay placer's per-leg encoding
- preflight drift detection: returns price_moved when WZ shows different odds
- successful path: returns {status: 'placed', ticket_number: ...}
- error classification: auth_error / rejected / network_error / orphaned
- empty details from ConfirmWagerHelper → rejected (line pulled)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from unittest.mock import MagicMock, patch

import single_placer


def make_bet(**overrides):
    """Default 1-leg bet payload for tests."""
    bet = dict(
        bet_hash         = "abc123",
        bookmaker_key    = "wagerzon",
        idgm             = 100001,    # WZ game id
        play             = 1,         # 1 = home spread
        line             = -1.5,
        american_odds    = 110,
        kelly_bet        = 42,
        actual_size      = 42,
        wz_odds_at_place = 110,       # what client thought WZ was showing
    )
    bet.update(overrides)
    return bet


def _make_json_response(payload):
    """Return a MagicMock that acts like a JSON response with non-HTML content-type."""
    r = MagicMock()
    r.headers = {"content-type": "application/json"}
    r.json.return_value = payload
    return r


def test_sel_encoding_single_leg():
    bet = make_bet()
    sel = single_placer.build_sel_for_single(bet)
    # Same shape as the parlay placer's per-leg encoding: play_idgm_pts_odds
    assert sel == "1_100001_-1.5_110"


def test_preflight_drift_returns_price_moved():
    """If WZ's ConfirmWagerHelper returns Odds different from wz_odds_at_place,
    placer must NOT call MakeWagerHelper and must return price_moved."""
    bet = make_bet(wz_odds_at_place=110)
    fake_session = MagicMock()
    fake_session.post.return_value = _make_json_response({
        "result": {"details": [{"Win": 4400, "Risk": 4000, "Odds": 100}]}
        # Odds 100 differs from client's 110 -> drift
    })
    result = single_placer.place_single(account="primary", bet=bet,
                                        session=fake_session)
    assert result["status"] == "price_moved"
    assert result.get("ticket_number") is None
    # MakeWagerHelper should NOT have been called
    make_calls = [c for c in fake_session.post.call_args_list
                  if "MakeWagerHelper" in str(c)]
    assert len(make_calls) == 0


def test_successful_placement_returns_ticket():
    bet = make_bet(wz_odds_at_place=110)
    fake_session = MagicMock()
    # ConfirmWagerHelper: no drift
    confirm_response = _make_json_response({
        "result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}
    })
    # MakeWagerHelper: success with ticket
    make_response = _make_json_response({
        "result": {"WagerNumber": "T-9876", "AvailBalance": 153.50}
    })
    fake_session.post.side_effect = [confirm_response, make_response]

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-9876"
    assert result["balance_after"] == pytest.approx(153.50)


def test_wz_rejection_returns_rejected_with_reason():
    bet = make_bet()
    fake_session = MagicMock()
    confirm_response = _make_json_response({
        "result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}
    })
    make_response = _make_json_response({
        "result": {"ErrorCode": "INSUFFICIENT_BALANCE",
                   "ErrorMessage": "Not enough funds"}
    })
    fake_session.post.side_effect = [confirm_response, make_response]

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "rejected"
    assert "INSUFFICIENT_BALANCE" in (result.get("error_msg_key") or "") + \
           (result.get("error_msg") or "")


def test_network_error_returns_network_error():
    import requests
    bet = make_bet()
    fake_session = MagicMock()
    fake_session.post.side_effect = requests.ConnectionError("connection reset")

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "network_error"


def test_auth_error_returns_auth_error():
    """WZ returns an HTML page (session expired) — placer must return auth_error."""
    bet = make_bet()
    fake_session = MagicMock()
    # Simulate HTML response (non-JSON content-type) — same as parlay_placer detects
    auth_response = MagicMock()
    auth_response.headers = {"content-type": "text/html; charset=utf-8"}
    auth_response.json.side_effect = Exception("not JSON")
    fake_session.post.return_value = auth_response

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "auth_error"


def test_empty_details_treated_as_rejected():
    """ConfirmWagerHelper sometimes returns details=[] (line pulled) — mirror
    the parlay placer's behavior of converting this into 'rejected'."""
    bet = make_bet()
    fake_session = MagicMock()
    confirm_response = _make_json_response({"result": {"details": []}})
    fake_session.post.return_value = confirm_response

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "rejected"
    assert "details" in (result.get("error_msg") or "").lower() or \
           result.get("error_msg_key") == "empty_details"
