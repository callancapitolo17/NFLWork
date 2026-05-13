"""Tests for wagerzon_odds.single_placer

Mocks the WZ HTTP session and verifies:
- sel string for a 1-leg payload matches parlay placer's per-leg encoding
- preflight drift detection: returns price_moved when WZ shows different odds
- successful path: returns {status: 'placed', ticket_number: ...}
- error classification: auth_error / rejected / network_error / orphaned
- empty details from ConfirmWagerHelper → rejected (line pulled)
- confirmPassword included in PostWagerMultipleHelper payload (production path)
- both WZ error key conventions handled (ErrorMsgKey vs ErrorCode)
- orphaned when WagerNumber missing from make response
- missing_odds in preflight details → rejected
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
    placer must NOT call PostWagerMultipleHelper and must return price_moved."""
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
    # PostWagerMultipleHelper should NOT have been called
    make_calls = [c for c in fake_session.post.call_args_list
                  if "PostWagerMultipleHelper" in str(c)]
    assert len(make_calls) == 0


def test_successful_placement_returns_ticket():
    bet = make_bet(wz_odds_at_place=110)
    fake_session = MagicMock()
    # ConfirmWagerHelper: no drift
    confirm_response = _make_json_response({
        "result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}
    })
    # PostWagerMultipleHelper: success with ticket
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


def test_successful_placement_includes_confirm_password(monkeypatch):
    """Production path must send confirmPassword in PostWagerMultipleHelper.

    confirmPassword is required by WZ to authorize the placement. The test
    path (session= injected) skips it because mock sessions do not validate
    credentials.
    """
    from wagerzon_accounts import WagerzonAccount
    fake_account = WagerzonAccount(
        label="Wagerzon", suffix="", username="u", password="secret123"
    )
    monkeypatch.setattr("single_placer.get_account", lambda label: fake_account)

    fake_session = MagicMock()
    confirm = _make_json_response(
        {"result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}}
    )
    make = _make_json_response(
        {"result": {"WagerNumber": "T-1", "AvailBalance": 100.0}}
    )
    fake_session.post.side_effect = [confirm, make]

    monkeypatch.setattr("single_placer.wagerzon_auth.get_session",
                        lambda acct: fake_session)

    result = single_placer.place_single("Wagerzon", make_bet())

    # Inspect the second post call's data payload for confirmPassword
    make_call = fake_session.post.call_args_list[1]
    # post() is called as sess.post(URL, data=payload, ...) — keyword arg
    payload = make_call[1].get("data", {})
    assert payload.get("confirmPassword") == "secret123"
    assert result["status"] == "placed"


def test_wz_error_msg_key_convention_also_rejected():
    """WZ ConfirmWagerHelper may use ErrorMsgKey/ErrorMsg instead of
    ErrorCode/ErrorMessage. Both must be classified as 'rejected'."""
    bet = make_bet()
    fake_session = MagicMock()
    # ErrorMsgKey is the parlay_placer convention; single_placer must handle it too
    confirm_response = _make_json_response({
        "result": {
            "details": [],
            "ErrorMsgKey": "line_unavailable",
            "ErrorMsg": "Line is no longer available",
        }
    })
    fake_session.post.return_value = confirm_response

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "rejected"
    assert result.get("error_msg_key") == "line_unavailable"


def test_orphaned_when_make_succeeds_but_no_wager_number():
    """WZ confirmed submission but didn't return WagerNumber -> orphaned."""
    fake_session = MagicMock()
    confirm = _make_json_response(
        {"result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}}
    )
    make = _make_json_response(
        {"result": {"AvailBalance": 100.0}}  # no WagerNumber
    )
    fake_session.post.side_effect = [confirm, make]

    result = single_placer.place_single("primary", make_bet(), session=fake_session)
    assert result["status"] == "orphaned"
    assert result["error_msg_key"] == "no_ticket"


def test_missing_odds_in_preflight_details_returns_rejected():
    """Preflight details present but missing Odds field -> rejected."""
    fake_session = MagicMock()
    fake_response = _make_json_response({
        "result": {"details": [{"Win": 4620, "Risk": 4200}]}  # no Odds
    })
    fake_session.post.return_value = fake_response

    result = single_placer.place_single("primary", make_bet(), session=fake_session)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "missing_odds"
