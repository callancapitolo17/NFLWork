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
    """v2: POST /portfolio/events/orders; buy YES = bid @ price; buy NO =
    sell YES (ask) @ (1 - price); price/count are decimal-dollar strings."""
    calls = []
    def fake_api(method, path, body=None, timeout=30):
        calls.append((method, path, body))
        return 201, {"order_id": "ord-1"}, {}          # v2: order_id at top level
    monkeypatch.setattr(gateway.auth_client, "api", fake_api)
    gw = gateway.LiveGateway()
    assert gw.place("T-O25", "yes", 40, 30, "c1") == "ord-1"
    assert gw.place("T-O25", "no", 57, 12, "c2") == "ord-1"
    (m1, p1, b1), (m2, p2, b2) = calls
    assert (m1, p1) == ("POST", "/portfolio/events/orders")
    # buy YES @ 40c -> bid @ "0.4000"
    assert b1["ticker"] == "T-O25" and b1["side"] == "bid" and b1["price"] == "0.4000"
    assert b1["count"] == "30.00" and b1["time_in_force"] == "good_till_canceled"
    assert b1["self_trade_prevention_type"] == "taker_at_cross" and b1["client_order_id"] == "c1"
    # buy NO @ 57c -> sell YES (ask) @ (1 - 0.57) = "0.4300"
    assert b2["side"] == "ask" and b2["price"] == "0.4300"


def test_live_place_reads_nested_order_id(monkeypatch):
    monkeypatch.setattr(gateway.auth_client, "api",
        lambda m, p, body=None, timeout=30: (201, {"order": {"order_id": "nested-1"}}, {}))
    assert gateway.LiveGateway().place("T", "yes", 40, 30, "c") == "nested-1"


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
    assert seen["call"] == ("DELETE", "/portfolio/events/orders/ord-9")
