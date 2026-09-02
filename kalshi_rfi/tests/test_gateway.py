"""Unit tests for gateway failure semantics (no network — api monkeypatched)."""
import urllib.error

import pytest

from kalshi_common import auth_client
from kalshi_rfi.gateway import LiveGateway, make_gateway


class TestCancelRouting:
    """Bug 2 (observed live 2026-08-27): Kalshi sharded its exchange
    (~2026-08-24) and baseball moved to shard 3. DELETE carries no ticker in
    its path, so an unrouted cancel hits shard 0 and 404s — and the
    404-means-gone rule then reported success while 13 real orders kept
    resting."""

    def test_cancel_targets_the_market_shard(self, monkeypatch):
        seen = {}

        def api(method, path, *a, **k):
            seen["method"], seen["path"] = method, path
            return 200, {"reduced_by": "13.00"}, None

        monkeypatch.setattr(auth_client, "api", api)
        assert LiveGateway().cancel("oid-1", ticker="KXMLBRFI-X",
                                    exchange_index=3) is True
        assert seen["method"] == "DELETE"
        assert "exchange_index=3" in seen["path"]
        assert "market_ticker=KXMLBRFI-X" in seen["path"]

    def test_unknown_shard_falls_back_to_ticker_autorouting(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(auth_client, "api",
                            lambda m, path, *a, **k: (seen.update(path=path)
                                                      or (200, {}, None)))
        LiveGateway().cancel("oid-1", ticker="KXMLBRFI-X")
        assert seen["path"].endswith("?market_ticker=KXMLBRFI-X")
        assert "exchange_index" not in seen["path"]


class TestCancelSemantics:
    def test_404_verified_gone_is_success(self, monkeypatch):
        # MM phantom-open-quotes lesson (resolved 2026-08-12): an order that
        # really is terminal must clear local state, or the ticker wedges
        # with a phantom resting entry forever.
        def api(method, path, *a, **k):
            if method == "DELETE":
                return 404, {"error": "not_found"}, None
            return 200, {"orders": [{"order_id": "other"}]}, None

        monkeypatch.setattr(auth_client, "api", api)
        assert LiveGateway().cancel("oid-1") is True

    def test_404_while_still_resting_is_failure(self, monkeypatch):
        # The misrouted cancel: Kalshi says not_found, but the order is
        # right there in the resting listing. Reporting success here is
        # what stranded live orders.
        def api(method, path, *a, **k):
            if method == "DELETE":
                return 404, {"error": "not_found"}, None
            return 200, {"orders": [{"order_id": "oid-1",
                                     "ticker": "KXMLBRFI-X"}]}, None

        monkeypatch.setattr(auth_client, "api", api)
        assert LiveGateway().cancel("oid-1") is False

    def test_404_with_unverifiable_listing_fails_closed(self, monkeypatch):
        def api(method, path, *a, **k):
            return (404, {}, None) if method == "DELETE" else (503, {}, None)

        monkeypatch.setattr(auth_client, "api", api)
        assert LiveGateway().cancel("oid-1") is False

    def test_5xx_is_failure(self, monkeypatch):
        monkeypatch.setattr(auth_client, "api",
                            lambda *a, **k: (503, {}, None))
        assert LiveGateway().cancel("oid-1") is False

    def test_network_error_is_failure_not_crash(self, monkeypatch):
        def boom(*a, **k):
            raise urllib.error.URLError("connection refused")
        monkeypatch.setattr(auth_client, "api", boom)
        assert LiveGateway().cancel("oid-1") is False


class TestPlaceSemantics:
    def test_place_sends_exchange_index_in_the_body(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(auth_client, "api",
                            lambda m, p, body=None, **k: (seen.update(body=body)
                                                          or (201, {"order_id": "o1"}, None)))
        assert LiveGateway().place_yes_bid("T", 41, 10, "c1",
                                           exchange_index=3) == "o1"
        assert seen["body"]["exchange_index"] == 3

    def test_place_omits_exchange_index_when_unknown(self, monkeypatch):
        # Omitted, Kalshi auto-routes on the ticker — a hardcoded 3 would be
        # wrong for every non-baseball series.
        seen = {}
        monkeypatch.setattr(auth_client, "api",
                            lambda m, p, body=None, **k: (seen.update(body=body)
                                                          or (201, {"order_id": "o1"}, None)))
        LiveGateway().place_yes_bid("T", 41, 10, "c1")
        assert "exchange_index" not in seen["body"]

    def test_network_error_returns_none_not_crash(self, monkeypatch):
        def boom(*a, **k):
            raise TimeoutError("read timed out")
        monkeypatch.setattr(auth_client, "api", boom)
        assert LiveGateway().place_yes_bid("T", 41, 10, "c1") is None


class TestDeadManSwitch:
    def test_live_requires_ack(self):
        with pytest.raises(SystemExit):
            make_gateway("live", None)

    def test_off_returns_none_and_shadow_places_nothing(self):
        assert make_gateway("off", None) is None
        gw = make_gateway("shadow", None)
        assert gw.place_yes_bid("T", 41, 10, "c1").startswith("shadow-")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
