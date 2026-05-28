from kalshi_mlb_mm import quote_gateway
from kalshi_common import auth_client


def test_submit_returns_quote_id(monkeypatch):
    calls = {}
    def fake_api(method, path, body=None, timeout=30):
        calls["method"], calls["path"], calls["body"] = method, path, body
        return 201, {"id": "q123"}, {}
    monkeypatch.setattr(auth_client, "api", fake_api)
    gw = quote_gateway.RestQuoteGateway()
    qid = gw.submit_quote("rfq1", 0.52, 0.43)
    assert qid == "q123"
    assert calls["method"] == "POST" and calls["path"] == "/communications/quotes"
    assert calls["body"]["yes_bid"] == "0.5200" and calls["body"]["rest_remainder"] is False


def test_submit_none_on_error(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (400, {"error": "x"}, {}))
    assert quote_gateway.RestQuoteGateway().submit_quote("rfq1", 0.5, 0.4) is None


def test_confirm_true_on_204(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (204, {}, {}))
    assert quote_gateway.RestQuoteGateway().confirm("q1") is True


def test_cancel_true_on_204(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (204, {}, {}))
    assert quote_gateway.RestQuoteGateway().cancel("q1") is True
