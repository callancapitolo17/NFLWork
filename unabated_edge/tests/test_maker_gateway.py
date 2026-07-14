import pytest
from unabated_edge.maker import gateway


def test_make_gateway_modes():
    assert gateway.make_gateway("off", None) is None
    assert isinstance(gateway.make_gateway("shadow", None), gateway.ShadowGateway)
    with pytest.raises(SystemExit):
        gateway.make_gateway("live", None)         # dead-man switch
    assert isinstance(gateway.make_gateway("live", "1"), gateway.LiveGateway)
    with pytest.raises(SystemExit):
        gateway.make_gateway("bogus", None)


def test_shadow_place_and_cancel():
    gw = gateway.ShadowGateway()
    oid = gw.place("T-O25", "yes", 40, 30, "coid")
    assert oid == "shadow-1" and not gw.is_live
    assert gw.cancel(oid) is True


def test_live_place_payload_yes_and_no(monkeypatch):
    calls = []
    def fake_api(method, path, body=None, timeout=30):
        calls.append((method, path, body))
        return 201, {"order": {"order_id": "ord-1"}}, {}
    monkeypatch.setattr(gateway.auth_client, "api", fake_api)
    gw = gateway.LiveGateway()
    assert gw.place("T-O25", "yes", 40, 30, "c1") == "ord-1"
    assert gw.place("T-O25", "no", 57, 12, "c2") == "ord-1"
    (m1, p1, b1), (m2, p2, b2) = calls
    assert (m1, p1) == ("POST", "/portfolio/orders")
    assert b1 == {"ticker": "T-O25", "client_order_id": "c1", "action": "buy",
                  "side": "yes", "type": "limit", "count": 30, "yes_price": 40}
    assert b2["no_price"] == 57 and "yes_price" not in b2


def test_live_place_failure_returns_none(monkeypatch):
    monkeypatch.setattr(gateway.auth_client, "api", lambda m, p, body=None, timeout=30: (500, "err", {}))
    assert gateway.LiveGateway().place("T", "yes", 40, 30, "c") is None


def test_live_cancel(monkeypatch):
    seen = {}
    def fake_api(method, path, body=None, timeout=30):
        seen["call"] = (method, path)
        return 200, {}, {}
    monkeypatch.setattr(gateway.auth_client, "api", fake_api)
    assert gateway.LiveGateway().cancel("ord-9") is True
    assert seen["call"] == ("DELETE", "/portfolio/orders/ord-9")
