from unittest.mock import patch, MagicMock

from kalshi_mlb_rfq import auth_client


def test_get_returns_parsed_json():
    fake_resp = MagicMock()
    fake_resp.status = 200
    fake_resp.read.return_value = b'{"foo": "bar"}'
    fake_resp.headers = {}
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False

    with patch("kalshi_mlb_rfq.auth_client.urllib.request.urlopen", return_value=fake_resp), \
         patch("kalshi_mlb_rfq.auth_client._sign", return_value=("sig", "1234")):
        status, body, _ = auth_client.api("GET", "/exchange/status")

    assert status == 200
    assert body == {"foo": "bar"}


def test_delete_handles_empty_body():
    fake_resp = MagicMock()
    fake_resp.status = 204
    fake_resp.read.return_value = b''
    fake_resp.headers = {}
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False

    with patch("kalshi_mlb_rfq.auth_client.urllib.request.urlopen", return_value=fake_resp), \
         patch("kalshi_mlb_rfq.auth_client._sign", return_value=("sig", "1234")):
        status, body, _ = auth_client.api("DELETE", "/communications/rfqs/abc")

    assert status == 204
    assert body == {}
