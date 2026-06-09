from kalshi_mlb_mm import rfq_source
from kalshi_common import auth_client

def test_poll_returns_rfqs(monkeypatch):
    def fake_api(method, path, body=None, timeout=30):
        assert path.startswith("/communications/rfqs?status=open")
        return 200, {"rfqs": [{"id": "r1", "market_ticker": "MT", "contracts": 3}]}, {}
    monkeypatch.setattr(auth_client, "api", fake_api)
    out = rfq_source.RestRFQSource().poll()
    assert out[0]["id"] == "r1"

def test_poll_empty_on_error(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (500, "x", {}))
    assert rfq_source.RestRFQSource().poll() == []

def test_get_market_returns_inner(monkeypatch):
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (200, {"market": {"ticker": "MT", "mve_selected_legs": []}}, {}))
    m = rfq_source.RestRFQSource().get_market("MT")
    assert m["ticker"] == "MT"

def test_get_market_none_on_error(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (404, {"error": "x"}, {}))
    assert rfq_source.RestRFQSource().get_market("MT") is None
