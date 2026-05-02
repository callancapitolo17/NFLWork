"""Unit tests for wagerzon_balance."""
from __future__ import annotations

from datetime import datetime, timezone
import pytest
import requests_mock

from wagerzon_accounts import WagerzonAccount
import wagerzon_auth
import wagerzon_balance


@pytest.fixture(autouse=True)
def clean_caches():
    wagerzon_auth.clear_session_cache()
    yield
    wagerzon_auth.clear_session_cache()


@pytest.fixture
def acct():
    return WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


_LOGGED_IN_URL = wagerzon_auth.WAGERZON_BASE_URL.rstrip("/") + "/NewSchedule.aspx"


def _mock_login(m):
    """Helper: mock the login endpoints so get_session() succeeds."""
    LOGIN_HTML = (
        '<html><form>'
        '<input name="__VIEWSTATE" value="x"/>'
        '<input name="__VIEWSTATEGENERATOR" value="x"/>'
        '<input name="__EVENTVALIDATION" value="x"/>'
        '</form></html>'
    )
    m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_HTML)
    m.post(wagerzon_auth.WAGERZON_BASE_URL, text="",
           headers={"Location": _LOGGED_IN_URL},
           status_code=302)
    m.get(_LOGGED_IN_URL, text="logged in")


# Real WZ balance endpoint and a sample of its response shape.
BALANCE_URL = "https://backend.wagerzon.com/wager/PlayerInfoHelper.aspx"
MOCK_BALANCE_RESPONSE = {
    "result": {
        "AmountAtRisk": "715 ",
        "AvailBalance": "1,245.32 ",       # the gating number we display
        "BonusPoints": 0.0000,
        "CreditLimit": "2,000 ",
        "CurrentBalance": "1,300.00 ",     # the "cash" we expose for tooltip/debug
        "FreePlayAmount": "0 ",
        "RealAvailBalance": "3,245.32 ",
        # NOTE: real responses also include "Player" and "Password" fields
        # (yes, the password in plaintext). The parser MUST extract only the
        # numeric balance fields and never log the raw response body.
        "ErrorCode": {},
        "ErrorMsg": "",
    }
}


def test_fetch_returns_snapshot(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=MOCK_BALANCE_RESPONSE)
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.label == "Wagerzon"
        assert snap.available == 1245.32
        assert snap.cash == 1300.00
        assert snap.error is None
        assert snap.fetched_at.tzinfo == timezone.utc


def test_fetch_does_not_log_password_from_response(acct, caplog):
    """The WZ balance response includes the user's password in plaintext.
    The fetcher must never log the raw response body."""
    import logging
    response_with_password = {
        "result": {
            "AvailBalance": "100 ", "CurrentBalance": "100 ",
            "Player": "MYACCT", "Password": "supersecret-do-not-log",
            "ErrorCode": {}, "ErrorMsg": "",
        }
    }
    caplog.set_level(logging.DEBUG)
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=response_with_password)
        wagerzon_balance.fetch_available_balance(acct)
    assert "supersecret-do-not-log" not in caplog.text


def test_fetch_timeout_returns_error_snapshot(acct):
    import requests
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, exc=requests.exceptions.ConnectTimeout)
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error == "timeout"
        assert snap.available is None


def test_fetch_5xx_returns_wz_error_snapshot(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, status_code=503, text="oops")
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error == "wz_error"
        assert snap.available is None


def test_fetch_401_triggers_relogin_then_retries_once(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        # First balance call: 401. Second: success.
        m.get(BALANCE_URL, [
            {"status_code": 401, "text": "unauthorized"},
            {"json": MOCK_BALANCE_RESPONSE, "status_code": 200},
        ])
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error is None
        assert snap.available == 1245.32


def test_fetch_all_runs_in_parallel():
    a = WagerzonAccount(label="Wagerzon",  suffix="",  username="u1", password="p1")
    b = WagerzonAccount(label="WagerzonJ", suffix="J", username="u2", password="p2")
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=MOCK_BALANCE_RESPONSE)
        snaps = wagerzon_balance.fetch_all([a, b])
        labels = sorted(s.label for s in snaps)
        assert labels == ["Wagerzon", "WagerzonJ"]
