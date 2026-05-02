"""Unit tests for wagerzon_auth (login + session caching)."""
from __future__ import annotations

import re
import pytest
import requests_mock

from wagerzon_accounts import WagerzonAccount
import wagerzon_auth


@pytest.fixture
def acct():
    return WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


@pytest.fixture(autouse=True)
def reset_cache():
    wagerzon_auth.clear_session_cache()
    yield
    wagerzon_auth.clear_session_cache()


LOGIN_PAGE_HTML = (
    '<html><form>'
    '<input name="__VIEWSTATE" value="vs1"/>'
    '<input name="__VIEWSTATEGENERATOR" value="vsg1"/>'
    '<input name="__EVENTVALIDATION" value="ev1"/>'
    '</form></html>'
)
LOGGED_IN_URL = wagerzon_auth.WAGERZON_BASE_URL.rstrip("/") + "/NewSchedule.aspx"


def test_first_call_logs_in(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        m.get(LOGGED_IN_URL, text="logged in")
        sess = wagerzon_auth.get_session(acct)
        assert isinstance(sess, requests_mock.Adapter) is False  # got a Session
        assert m.call_count == 3  # GET login page + POST login form + GET redirect target


def test_second_call_returns_cached(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        m.get(LOGGED_IN_URL, text="logged in")
        s1 = wagerzon_auth.get_session(acct)
        s2 = wagerzon_auth.get_session(acct)
        assert s1 is s2
        assert m.call_count == 3  # GET login + POST + GET redirect; second call is cached


def test_already_logged_in_skips_form_post(acct):
    """If GET base URL bounces straight to NewSchedule, no POST required."""
    with requests_mock.Mocker() as m:
        # The cookie jar already has a valid session — GET redirects.
        m.get(wagerzon_auth.WAGERZON_BASE_URL,
              status_code=302, headers={"Location": LOGGED_IN_URL})
        m.get(LOGGED_IN_URL, text="hello")
        sess = wagerzon_auth.get_session(acct)
        # Only the redirect chain happened; no form POST.
        post_calls = [r for r in m.request_history if r.method == "POST"]
        assert len(post_calls) == 0


def test_two_accounts_get_separate_sessions():
    a = WagerzonAccount(label="Wagerzon",  suffix="",  username="u1", password="p1")
    b = WagerzonAccount(label="WagerzonJ", suffix="J", username="u2", password="p2")
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        m.get(LOGGED_IN_URL, text="logged in")
        s_a = wagerzon_auth.get_session(a)
        s_b = wagerzon_auth.get_session(b)
        assert s_a is not s_b


def test_clear_session_cache_forces_relogin(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        m.get(LOGGED_IN_URL, text="logged in")
        wagerzon_auth.get_session(acct)
        wagerzon_auth.clear_session_cache(label=acct.label)
        wagerzon_auth.get_session(acct)
        # Two full login sequences (GET + POST + GET redirect each) = 6 calls
        assert m.call_count == 6
