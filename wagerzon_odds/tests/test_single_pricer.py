"""Unit tests for wagerzon_odds.single_pricer.get_single_price.

Mirrors test_parlay_placer.py pattern: pass a MagicMock session into
get_single_price so we never hit the live WZ API.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
from unittest.mock import MagicMock

import pytest

from single_pricer import get_single_price


def _mock_response(json_body=None, status=200, content_type="application/json", raise_exc=None):
    """Build a MagicMock requests.Response with the given JSON body."""
    if raise_exc is not None:
        sess = MagicMock()
        sess.post.side_effect = raise_exc
        return sess
    resp = MagicMock()
    resp.status_code = status
    resp.headers = {"content-type": content_type}
    resp.json.return_value = json_body or {}
    sess = MagicMock()
    sess.post.return_value = resp
    return sess


def _wz_confirm_body(win, odds, error_key=None, error_msg=None):
    """Build a fake ConfirmWagerHelper response shape for a single."""
    body = {
        "result": {
            "details": [
                {
                    "Amount": 25.0,
                    "Risk": 25.0,
                    "Win": win,
                    "details": [{"Odds": odds, "IsMLine": True}],
                }
            ],
            "ErrorMsgKey": error_key or "",
            "ErrorMsg":    error_msg or "",
        }
    }
    return body


@pytest.fixture
def bet():
    return {
        "idgm": 5632938,
        "play": 5,            # home ML per WZ play codes
        "line": 0.0,
        "american_odds": -140,
        "amount": 25.0,
        "pitcher": 0,
    }


def test_happy_path_returns_win_and_current_odds(bet):
    sess = _mock_response(_wz_confirm_body(win=17.86, odds=-140))
    out = get_single_price(sess, bet, amount=25)
    assert out["win"] == pytest.approx(17.86)
    assert out["current_wz_odds"] == -140
    assert out["error_msg"] == ""
    assert out["error_msg_key"] == ""


def test_error_msg_surfaced(bet):
    sess = _mock_response(_wz_confirm_body(
        win=0, odds=0,
        error_key="MAXMONEYLINEEXCEEDED",
        error_msg="Maximum money line risk exceeded",
    ))
    out = get_single_price(sess, bet, amount=10000)
    assert out["error_msg_key"] == "MAXMONEYLINEEXCEEDED"
    assert out["error_msg"] == "Maximum money line risk exceeded"
    assert out["win"] is None


def test_html_response_returns_auth_error(bet):
    sess = _mock_response(json_body={}, content_type="text/html")
    out = get_single_price(sess, bet, amount=25)
    assert out["error_msg_key"] == "session_expired"


def test_network_exception_returns_network_error(bet):
    import requests
    sess = _mock_response(raise_exc=requests.RequestException("connreset"))
    out = get_single_price(sess, bet, amount=25)
    assert out["error_msg_key"] == "network_error"
    assert "connreset" in out["error_msg"]
