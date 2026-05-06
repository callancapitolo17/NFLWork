# wagerzon_odds/test_parlay_placer.py
"""Unit tests for parlay_placer. All Wagerzon HTTP calls are mocked."""
import pytest
import json
from unittest.mock import patch, MagicMock
from parlay_placer import Leg, ParlaySpec, encode_sel, encode_detail_data
from wagerzon_accounts import WagerzonAccount


@pytest.fixture
def primary_acct():
    """A throwaway WagerzonAccount used to satisfy the explicit-account
    contract on place_parlays / _confirm_preflight / _post_wagers. Tests
    never actually log in (they patch wagerzon_auth.get_session), so the
    credentials don't have to be real — but `password` IS used directly
    by _post_wagers as the `confirmPassword` field, so set it to
    something distinguishable in case a future test asserts on it."""
    return WagerzonAccount(
        label="Wagerzon", suffix="",
        username="test_user", password="test_pw",
    )


@pytest.fixture(autouse=True)
def reset_session():
    """Clear cached sessions in wagerzon_auth before/after each test to
    avoid leakage between cases. Phase 4 dropped the old module-level
    _CACHED_SESSION on parlay_placer; the cache now lives in wagerzon_auth."""
    import wagerzon_auth
    wagerzon_auth.clear_session_cache()
    yield
    wagerzon_auth.clear_session_cache()


def test_leg_dataclass_basic():
    leg = Leg(idgm=5632938, play=0, points=-1.5, odds=117, pitcher=0)
    assert leg.idgm == 5632938
    assert leg.play == 0


def test_parlay_spec_dataclass_basic():
    spec = ParlaySpec(
        parlay_hash="abc123",
        legs=[
            Leg(idgm=5632938, play=1, points=-1.5, odds=117),
            Leg(idgm=5632938, play=2, points=-7.5, odds=105),
        ],
        amount=15.0,
        expected_win=12.50,
        expected_risk=15.0,
    )
    assert len(spec.legs) == 2


def test_encode_sel_two_legs():
    legs = [
        Leg(idgm=5632938, play=1, points=-1.5, odds=117),
        Leg(idgm=5632938, play=2, points=-7.5, odds=105),
    ]
    assert encode_sel(legs) == "1_5632938_-1.5_117,2_5632938_-7.5_105"


def test_encode_sel_negative_odds():
    legs = [Leg(idgm=5632938, play=5, points=0, odds=-140)]
    assert encode_sel(legs) == "5_5632938_0_-140"


def test_encode_detail_data_shape():
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    out = encode_detail_data(legs, amount=15.0)
    # detail_data is a list of dicts ready for json.dumps
    assert out[0]["IdGame"] == 5632938
    assert out[0]["Play"] == 1
    assert out[0]["Amount"] == "15"
    assert out[0]["Points"]["selected"] is True


def test_encode_detail_data_non_integer_amount():
    """Kelly recommendations may be fractional dollars (e.g. $15.75)."""
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    out = encode_detail_data(legs, amount=15.75)
    assert out[0]["Amount"] == "15.75"


CONFIRM_OK_RESPONSE = {
    "result": {
        "details": [{"Risk": 15.0, "Win": 30.0, "WagerType": 1,
                     "WagerTypeDesc": "PARLAY (2 TEAMS)"}],
        "Confirm": True,
    }
}

POST_OK_RESPONSE = {
    "result": [{
        "WagerPostResult": {
            "details": [],
            "IDWT": 438173,
            "TicketNumber": "212147131",
            "Risk": 15.0,
            "Win": 30.0,
            "Confirm": True,
            "ErrorMsg": "",
            "ErrorMsgKey": "",
            "ErrorCode": {},
        }
    }]
}

POST_REJECTED_RESPONSE = {
    "result": [{
        "WagerPostResult": {
            "details": [],
            "Confirm": False,
            "ErrorMsg": "Insufficient funds",
            "ErrorMsgKey": "insufficient_funds",
            "ErrorCode": {"code": 42},
        }
    }]
}


def _make_session_with_post(json_response):
    """Build a mocked requests.Session whose .post returns the given JSON."""
    sess = MagicMock()
    r = MagicMock()
    r.json.return_value = json_response
    r.headers = {"content-type": "application/json"}
    r.status_code = 200
    r.text = json.dumps(json_response)
    sess.post.return_value = r
    return sess


def test_confirm_preflight_returns_win_risk(monkeypatch, primary_acct):
    sess = _make_session_with_post(CONFIRM_OK_RESPONSE)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    win, risk = parlay_placer._confirm_preflight(legs, amount=15.0,
                                                 account=primary_acct)
    assert win == 30.0
    assert risk == 15.0


def test_drift_check_within_penny_passes():
    import parlay_placer
    assert parlay_placer._drift_ok(expected=30.00, actual=30.005) is True
    assert parlay_placer._drift_ok(expected=30.00, actual=29.995) is True


def test_drift_check_beyond_penny_fails():
    import parlay_placer
    assert parlay_placer._drift_ok(expected=30.00, actual=29.50) is False
    assert parlay_placer._drift_ok(expected=30.00, actual=28.50) is False


def test_post_wagers_success(monkeypatch, primary_acct):
    sess = _make_session_with_post(POST_OK_RESPONSE)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    specs = [ParlaySpec(
        parlay_hash="h1",
        legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
        amount=15.0, expected_win=30.0, expected_risk=15.0,
    )]
    results = parlay_placer._post_wagers(specs, account=primary_acct)
    assert len(results) == 1
    r = results[0]
    assert r.status == "placed"
    assert r.ticket_number == "212147131"
    assert r.idwt == 438173
    assert r.actual_win == 30.0


def test_post_wagers_rejected(monkeypatch, primary_acct):
    sess = _make_session_with_post(POST_REJECTED_RESPONSE)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    specs = [ParlaySpec(
        parlay_hash="h1",
        legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
        amount=15.0, expected_win=30.0, expected_risk=15.0,
    )]
    results = parlay_placer._post_wagers(specs, account=primary_acct)
    assert results[0].status == "rejected"
    assert results[0].error_msg_key == "insufficient_funds"


# ---------------------------------------------------------------------------
# Integration tests for place_parlays (Task 8)
# ---------------------------------------------------------------------------

def _spec(hash="h", win=30.0, risk=15.0):
    return ParlaySpec(
        parlay_hash=hash,
        legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
        amount=risk, expected_win=win, expected_risk=risk,
    )


def _session_for_calls(*responses_per_call):
    """Return a session whose .post returns these responses in order."""
    sess = MagicMock()
    mocks = []
    for r in responses_per_call:
        mr = MagicMock()
        mr.json.return_value = r
        mr.headers = {"content-type": "application/json"}
        mr.text = json.dumps(r)
        mocks.append(mr)
    sess.post.side_effect = mocks
    return sess


def test_place_parlays_happy_path(monkeypatch, primary_acct):
    sess = _session_for_calls(CONFIRM_OK_RESPONSE, POST_OK_RESPONSE)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "placed"


def test_place_parlays_price_drift_aborts(monkeypatch, primary_acct):
    drifted_response = {
        "result": {"details": [{"Win": 28.50, "Risk": 15.0}], "Confirm": True}
    }
    sess = _session_for_calls(drifted_response)  # only the preflight call expected
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    results = parlay_placer.place_parlays([_spec(win=30.0)], primary_acct)
    assert results[0].status == "price_moved"
    # Drift message uses two-decimal format (per _drift_error_msg)
    assert "$28.50" in results[0].error_msg
    assert "$30.00" in results[0].error_msg
    # PostWager should NOT have been called
    assert sess.post.call_count == 1


def test_place_parlays_auth_retry_succeeds(monkeypatch, primary_acct):
    """First preflight returns HTML (session expired) -> re-login -> retry -> OK."""
    import parlay_placer
    html_resp = MagicMock()
    html_resp.headers = {"content-type": "text/html"}
    html_resp.text = "<html>login</html>"
    ok_resp = MagicMock()
    ok_resp.headers = {"content-type": "application/json"}
    ok_resp.json.return_value = CONFIRM_OK_RESPONSE
    ok_resp.text = json.dumps(CONFIRM_OK_RESPONSE)
    post_ok = MagicMock()
    post_ok.headers = {"content-type": "application/json"}
    post_ok.json.return_value = POST_OK_RESPONSE
    post_ok.text = json.dumps(POST_OK_RESPONSE)

    session1 = MagicMock(); session1.post.return_value = html_resp
    session2 = MagicMock(); session2.post.side_effect = [ok_resp, post_ok]
    call_count = [0]

    def _get_sess(account):
        call_count[0] += 1
        return session1 if call_count[0] == 1 else session2

    monkeypatch.setattr("wagerzon_auth.get_session", _get_sess)
    monkeypatch.setattr("wagerzon_auth.clear_session_cache",
                        lambda label=None: None)
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "placed"


def test_place_parlays_network_error_no_retry(monkeypatch, primary_acct):
    import parlay_placer
    sess = MagicMock()
    sess.post.side_effect = parlay_placer.requests.exceptions.Timeout()
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "network_error"
    assert sess.post.call_count == 1  # no retry


def test_place_parlays_post_auth_retry_succeeds(monkeypatch, primary_acct):
    """Preflight OK -> post raises AuthExpired -> re-login -> re-preflight OK
    -> post succeeds. No drift, no network error."""
    import parlay_placer

    def _json_resp(payload):
        r = MagicMock()
        r.headers = {"content-type": "application/json"}
        r.json.return_value = payload
        r.text = json.dumps(payload)
        return r

    html_resp = MagicMock()
    html_resp.headers = {"content-type": "text/html"}
    html_resp.text = "<html>login</html>"

    sess = MagicMock()
    # Calls in order: preflight (OK) -> post (HTML, AuthExpired) ->
    #                 re-preflight after re-login (OK) -> post retry (OK)
    sess.post.side_effect = [
        _json_resp(CONFIRM_OK_RESPONSE),
        html_resp,
        _json_resp(CONFIRM_OK_RESPONSE),
        _json_resp(POST_OK_RESPONSE),
    ]
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    monkeypatch.setattr("wagerzon_auth.clear_session_cache",
                        lambda label=None: None)

    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "placed"
    assert sess.post.call_count == 4


def test_place_parlays_empty_list(monkeypatch, primary_acct):
    """Empty input must return empty output without making any API calls."""
    import parlay_placer
    sess = MagicMock()
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    assert parlay_placer.place_parlays([], primary_acct) == []
    assert sess.post.call_count == 0


def test_place_parlays_duplicate_hash_raises(monkeypatch, primary_acct):
    """Duplicate parlay_hash in a batch is a programming error, raise loudly."""
    import parlay_placer
    a = _spec(hash="dup")
    b = _spec(hash="dup")
    with pytest.raises(ValueError, match="duplicate parlay_hash"):
        parlay_placer.place_parlays([a, b], primary_acct)


CONFIRM_REJECTED_INSUFFICIENT_FUNDS = {
    "result": {
        "details": [],
        "ErrorMsgKey": "insufficient_funds",
        "ErrorMsg": "Insufficient funds.",
        "Confirm": False,
    }
}

CONFIRM_REJECTED_LINE_UNAVAILABLE = {
    "result": {
        "details": [],
        "ErrorMsgKey": "line_unavailable",
        "ErrorMsg": "One or more lines are no longer available.",
        "Confirm": False,
    }
}

CONFIRM_REJECTED_NO_KEY = {
    "result": {
        "details": [],
        "Confirm": False,
    }
}


def test_place_parlays_preflight_insufficient_funds(monkeypatch, primary_acct):
    """Empty details + ErrorMsgKey=insufficient_funds → rejected, no IndexError."""
    sess = _session_for_calls(CONFIRM_REJECTED_INSUFFICIENT_FUNDS)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "rejected"
    assert results[0].error_msg_key == "insufficient_funds"
    assert "insufficient balance" in results[0].error_msg
    # PostWager must NOT be called when preflight rejects
    assert sess.post.call_count == 1


def test_place_parlays_preflight_line_unavailable(monkeypatch, primary_acct):
    """Empty details + ErrorMsgKey=line_unavailable → rejected as 'line pulled'."""
    sess = _session_for_calls(CONFIRM_REJECTED_LINE_UNAVAILABLE)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "rejected"
    assert results[0].error_msg_key == "line_unavailable"
    assert "line pulled" in results[0].error_msg


def test_place_parlays_preflight_rejected_no_key(monkeypatch, primary_acct):
    """Empty details with no error key → graceful 'rejected', no IndexError."""
    sess = _session_for_calls(CONFIRM_REJECTED_NO_KEY)
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    import parlay_placer
    results = parlay_placer.place_parlays([_spec()], primary_acct)
    assert results[0].status == "rejected"
    assert results[0].error_msg_key == ""
    assert "rejected" in results[0].error_msg.lower()


def test_place_parlays_post_retry_per_spec_drift(monkeypatch, primary_acct):
    """Multi-spec batch where post AuthExpires; on retry one spec drifts.
    Only the drifted spec gets price_moved; others should still be placed."""
    import parlay_placer

    def _json_resp(payload):
        r = MagicMock()
        r.headers = {"content-type": "application/json"}
        r.json.return_value = payload
        r.text = json.dumps(payload)
        return r

    html_resp = MagicMock()
    html_resp.headers = {"content-type": "text/html"}
    html_resp.text = "<html>login</html>"

    drifted_confirm = {
        "result": {"details": [{"Win": 28.50, "Risk": 15.0}], "Confirm": True}
    }
    # Post-OK response for the SINGLE surviving spec (B), so its result list has 1 entry
    post_ok_one = {
        "result": [{
            "WagerPostResult": {
                "details": [], "IDWT": 999001, "TicketNumber": "TKT-B",
                "Risk": 15.0, "Win": 30.0, "Confirm": True,
                "ErrorMsg": "", "ErrorMsgKey": "", "ErrorCode": {},
            }
        }]
    }

    sess = MagicMock()
    # Sequence:
    #   preflight A OK, preflight B OK -> post (HTML AuthExpired) ->
    #   re-preflight A drifted, re-preflight B OK -> post B only OK
    sess.post.side_effect = [
        _json_resp(CONFIRM_OK_RESPONSE),  # A preflight
        _json_resp(CONFIRM_OK_RESPONSE),  # B preflight
        html_resp,                        # post -> AuthExpired
        _json_resp(drifted_confirm),      # A re-preflight (drifted)
        _json_resp(CONFIRM_OK_RESPONSE),  # B re-preflight (OK)
        _json_resp(post_ok_one),          # post B only
    ]
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    monkeypatch.setattr("wagerzon_auth.clear_session_cache",
                        lambda label=None: None)

    a = _spec(hash="A", win=30.0)
    b = _spec(hash="B", win=30.0)
    results = parlay_placer.place_parlays([a, b], primary_acct)

    # Order preserved: A first, B second
    assert results[0].parlay_hash == "A"
    assert results[1].parlay_hash == "B"
    assert results[0].status == "price_moved"
    assert results[1].status == "placed"
    assert results[1].ticket_number == "TKT-B"


def test_place_parlays_post_retry_per_spec_preflight_rejected(monkeypatch, primary_acct):
    """Multi-spec batch where post AuthExpires; on retry one spec's line is
    pulled (re-preflight returns empty details). Only the rejected spec gets
    status='rejected'; the other should still be placed."""
    import parlay_placer

    def _json_resp(payload):
        r = MagicMock()
        r.headers = {"content-type": "application/json"}
        r.json.return_value = payload
        r.text = json.dumps(payload)
        return r

    html_resp = MagicMock()
    html_resp.headers = {"content-type": "text/html"}
    html_resp.text = "<html>login</html>"

    # Post-OK response for the SINGLE surviving spec (A) — A is placed
    # alone after B is rejected on re-preflight.
    post_ok_one = {
        "result": [{
            "WagerPostResult": {
                "details": [], "IDWT": 999002, "TicketNumber": "TKT-A",
                "Risk": 15.0, "Win": 30.0, "Confirm": True,
                "ErrorMsg": "", "ErrorMsgKey": "", "ErrorCode": {},
            }
        }]
    }

    sess = MagicMock()
    # Sequence:
    #   preflight A OK, preflight B OK -> post (HTML AuthExpired) ->
    #   re-preflight A OK, re-preflight B line_unavailable ->
    #   post A only OK
    sess.post.side_effect = [
        _json_resp(CONFIRM_OK_RESPONSE),                 # A preflight
        _json_resp(CONFIRM_OK_RESPONSE),                 # B preflight
        html_resp,                                       # post -> AuthExpired
        _json_resp(CONFIRM_OK_RESPONSE),                 # A re-preflight (OK)
        _json_resp(CONFIRM_REJECTED_LINE_UNAVAILABLE),   # B re-preflight (rejected)
        _json_resp(post_ok_one),                         # post A only
    ]
    monkeypatch.setattr("wagerzon_auth.get_session", lambda acct: sess)
    monkeypatch.setattr("wagerzon_auth.clear_session_cache",
                        lambda label=None: None)

    a = _spec(hash="A", win=30.0)
    b = _spec(hash="B", win=30.0)
    results = parlay_placer.place_parlays([a, b], primary_acct)

    assert results[0].parlay_hash == "A"
    assert results[1].parlay_hash == "B"
    assert results[0].status == "placed"
    assert results[0].ticket_number == "TKT-A"
    assert results[1].status == "rejected"
    assert results[1].error_msg_key == "line_unavailable"
    assert "line pulled" in results[1].error_msg
