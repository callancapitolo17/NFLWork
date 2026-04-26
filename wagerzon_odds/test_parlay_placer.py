# wagerzon_odds/test_parlay_placer.py
"""Unit tests for parlay_placer. All Wagerzon HTTP calls are mocked."""
import pytest
import json
from unittest.mock import patch, MagicMock
from parlay_placer import Leg, ParlaySpec, encode_sel, encode_detail_data


@pytest.fixture(autouse=True)
def reset_session():
    """Clear cached session before/after each test to avoid leakage."""
    import parlay_placer
    parlay_placer._clear_session()
    yield
    parlay_placer._clear_session()


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


def test_get_session_caches_within_call(monkeypatch):
    """Two calls to _get_session in the same place_parlays() invocation
    should reuse one session, not log in twice."""
    import parlay_placer
    monkeypatch.setenv("WAGERZON_USERNAME", "user")
    monkeypatch.setenv("WAGERZON_PASSWORD", "pass")

    with patch("parlay_placer.requests.Session") as MockSession:
        mock_session = MagicMock()
        # Simulate already-authenticated GET (URL contains NewSchedule)
        response = MagicMock(url="https://backend.wagerzon.com/wager/NewSchedule.aspx")
        mock_session.get.return_value = response
        MockSession.return_value = mock_session

        s1 = parlay_placer._get_session()
        s2 = parlay_placer._get_session()
        assert s1 is s2
        # Only one Session() instantiation
        assert MockSession.call_count == 1


def test_get_session_form_post_login_path(monkeypatch):
    """When GET lands on the login page, _get_session parses ASP.NET tokens
    out of the HTML and form-POSTs Account/Password."""
    import parlay_placer
    monkeypatch.setenv("WAGERZON_USERNAME", "user42")
    monkeypatch.setenv("WAGERZON_PASSWORD", "secret-pw")

    login_html = """
    <html><body><form>
      <input name="__VIEWSTATE" value="VS123" />
      <input name="__VIEWSTATEGENERATOR" value="VSG456" />
      <input name="__EVENTVALIDATION" value="EV789" />
      <input name="Account" />
      <input name="Password" />
    </form></body></html>
    """

    with patch("parlay_placer.requests.Session") as MockSession:
        mock_session = MagicMock()
        get_resp = MagicMock(url="https://backend.wagerzon.com/")
        get_resp.text = login_html
        mock_session.get.return_value = get_resp
        # Login POST returns a generic OK
        post_resp = MagicMock()
        mock_session.post.return_value = post_resp
        MockSession.return_value = mock_session

        sess = parlay_placer._get_session()
        assert sess is mock_session

        # Verify the POST happened with the parsed tokens + creds
        mock_session.post.assert_called_once()
        post_kwargs = mock_session.post.call_args.kwargs
        data = post_kwargs["data"]
        assert data["__VIEWSTATE"] == "VS123"
        assert data["__VIEWSTATEGENERATOR"] == "VSG456"
        assert data["__EVENTVALIDATION"] == "EV789"
        assert data["Account"] == "user42"
        assert data["Password"] == "secret-pw"


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


def test_confirm_preflight_returns_win_risk(monkeypatch):
    monkeypatch.setattr("parlay_placer._get_session",
                        lambda: _make_session_with_post(CONFIRM_OK_RESPONSE))
    import parlay_placer
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    win, risk = parlay_placer._confirm_preflight(legs, amount=15.0)
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


def test_post_wagers_success(monkeypatch):
    monkeypatch.setattr("parlay_placer._get_session",
                        lambda: _make_session_with_post(POST_OK_RESPONSE))
    monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
    import parlay_placer
    specs = [ParlaySpec(
        parlay_hash="h1",
        legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
        amount=15.0, expected_win=30.0, expected_risk=15.0,
    )]
    results = parlay_placer._post_wagers(specs)
    assert len(results) == 1
    r = results[0]
    assert r.status == "placed"
    assert r.ticket_number == "212147131"
    assert r.idwt == 438173
    assert r.actual_win == 30.0


def test_post_wagers_rejected(monkeypatch):
    monkeypatch.setattr("parlay_placer._get_session",
                        lambda: _make_session_with_post(POST_REJECTED_RESPONSE))
    monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
    import parlay_placer
    specs = [ParlaySpec(
        parlay_hash="h1",
        legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
        amount=15.0, expected_win=30.0, expected_risk=15.0,
    )]
    results = parlay_placer._post_wagers(specs)
    assert results[0].status == "rejected"
    assert results[0].error_msg_key == "insufficient_funds"
