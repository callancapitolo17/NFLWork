import datetime
import pytest
from unabated_edge import config
from unabated_edge.maker import state as mstate

_NOW = datetime.datetime(2026, 7, 11, 12, 0, tzinfo=datetime.timezone.utc)


def _state_with_order():
    s = mstate.MakerState()
    s.register_ticker("T-O25", "soccer", 1, 2.5)
    s.on_place("T-O25", "yes", "ord-1", 40, 30)
    return s


def test_exposure_includes_resting_as_fills():
    s = _state_with_order()
    s.fills[1] = [(2.5, "no", 10, 0.58)]
    exp = s.exposure_fills(1)
    assert (2.5, "yes", 30, 0.40) in exp and (2.5, "no", 10, 0.58) in exp
    assert s.exposure_fills(1, exclude=("T-O25", "yes")) == [(2.5, "no", 10, 0.58)]


def test_poll_fills_parses_updates_and_dedups(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    from unabated_edge.storage import connect
    store.init()
    s = _state_with_order()
    payload = {"fills": [
        {"trade_id": "tr1", "order_id": "ord-1", "ticker": "T-O25", "side": "yes",
         "count_fp": "12.5", "yes_price_dollars": "0.4000", "no_price_dollars": "0.6000", "fee": 4},
        {"trade_id": "tr2", "order_id": "someone-else", "ticker": "T-O25", "side": "yes",
         "count": 5, "yes_price": 40},
    ]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    new = mstate.poll_fills(s, _NOW)
    assert new == ["tr1"]                                  # not ours -> ignored
    assert s.fills[1] == [(2.5, "yes", 12.5, 0.40)]
    assert s.fills_by_ticker["T-O25"] == 12.5
    assert s.resting[("T-O25", "yes")]["count"] == 17.5    # 30 - 12.5
    assert mstate.poll_fills(s, _NOW) == []                # dedup on trade_id
    # Verify fee is stored as dollars (4 cents = 0.04 dollars)
    with connect(config.MAKER_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT fee FROM maker_fills WHERE trade_id='tr1'").fetchone()[0] == 0.04


def test_poll_fills_full_fill_clears_resting(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    payload = {"fills": [{"trade_id": "trX", "order_id": "ord-1", "ticker": "T-O25",
                          "side": "yes", "count": 30, "yes_price": 40}]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    mstate.poll_fills(s, _NOW)
    assert ("T-O25", "yes") not in s.resting


def test_poll_positions_detects_mismatch(monkeypatch):
    s = _state_with_order()
    s.fills_by_ticker["T-O25"] = 12.5
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position_fp": "12.50"}]}, {}))
    assert mstate.poll_positions(s) is True
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position": 3}]}, {}))
    assert mstate.poll_positions(s) is False


def test_poll_settlements_updates_daily_pnl(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    s.fills_by_ticker["T-O25"] = 12.5
    payload = {"settlements": [{"ticker": "T-O25", "revenue": 1250,        # cents
                                "yes_total_cost": 500, "no_total_cost": 0}]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    mstate.poll_settlements(s, _NOW)
    assert round(s.settled_pnl_today, 2) == 7.50
    assert "T-O25" in s.settled and "T-O25" not in s.fills_by_ticker
    mstate.poll_settlements(s, _NOW)                       # idempotent
    assert round(s.settled_pnl_today, 2) == 7.50


def test_roll_day_resets_pnl():
    s = mstate.MakerState()
    s.roll_day(_NOW)
    s.settled_pnl_today = -100.0
    s.roll_day(_NOW + datetime.timedelta(days=1))
    assert s.settled_pnl_today == 0.0


def test_money_real_payload_shapes():
    assert mstate._money({"revenue": 16807}, "revenue") == 168.07
    assert mstate._money({"fee_cost": "2.168600"}, "fee_cost") == 2.1686
    assert mstate._money({"yes_total_cost_dollars": "127.112500"}, "yes_total_cost") == 127.1125
    assert mstate._money({}, "revenue") == 0.0


def test_poll_settlements_real_shape_subtracts_fee(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    payload = {"settlements": [{
        "ticker": "T-O25", "market_result": "yes", "revenue": 16807,
        "yes_total_cost_dollars": "127.112500", "no_total_cost_dollars": "0.000000",
        "fee_cost": "2.168600", "yes_count_fp": "168.07", "no_count_fp": "0.00", "value": 100,
    }]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    mstate.poll_settlements(s, _NOW)
    assert s.settled_pnl_today == pytest.approx(168.07 - 127.1125 - 2.1686)


def test_poll_positions_orphan_in_series_position_trips(monkeypatch):
    s = mstate.MakerState()
    payload = {"market_positions": [{"ticker": "T-O25", "position_fp": "5.00"}]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    assert mstate.poll_positions(s, ("T-",)) is False
    # No series_prefixes passed -> orphan never checked (documents the need for prefixes)
    assert mstate.poll_positions(s) is True


def test_poll_positions_baseline_restart(monkeypatch):
    s = mstate.MakerState()
    s.position_baseline["T-O25"] = -223.0
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position_fp": "-223.00"}]}, {}))
    assert mstate.poll_positions(s) is True
    s.fills_by_ticker["T-O25"] = 10.0
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position_fp": "-213.00"}]}, {}))
    assert mstate.poll_positions(s) is True
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position_fp": "-223.00"}]}, {}))
    assert mstate.poll_positions(s) is False


def test_startup_sync_adopts_and_cancels(monkeypatch):
    s = mstate.MakerState()
    s.on_place("T-O26", "yes", "o-ours", 40, 10)   # register one of our own orders

    def fake_api(method, path, body=None, timeout=30):
        if "/portfolio/positions" in path:
            return 200, {"market_positions": [
                {"ticker": "T-O25", "position_fp": "7.00"},
                {"ticker": "X-1", "position_fp": "9.00"},
            ]}, {}
        if "/portfolio/orders" in path:
            return 200, {"orders": [
                {"order_id": "o-orphan", "ticker": "T-O26"},
                {"order_id": "o-ours", "ticker": "T-O26"},
                {"order_id": "o-x", "ticker": "X-1"},
            ]}, {}
        raise AssertionError(f"unexpected path {path}")

    monkeypatch.setattr(mstate.auth_client, "api", fake_api)

    class _FakeGateway:
        def __init__(self):
            self.cancels = []
        def cancel(self, order_id):
            self.cancels.append(order_id)
            return True

    gw = _FakeGateway()
    mstate.startup_sync(s, gw, ("T-",))
    assert s.position_baseline == {"T-O25": 7.0}
    assert gw.cancels == ["o-orphan"]


def test_sweep_orphan_orders_failure_returns_zero(monkeypatch):
    s = mstate.MakerState()
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (500, "err", {}))
    class _FakeGateway:
        def cancel(self, order_id):
            return True
    assert mstate.sweep_orphan_orders(s, _FakeGateway(), ("T-",)) == 0
