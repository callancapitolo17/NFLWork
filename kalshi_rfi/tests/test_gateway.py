"""Unit tests for gateway failure semantics (no network — api monkeypatched)."""
import urllib.error

import pytest

from kalshi_common import auth_client
from kalshi_rfi.gateway import LiveGateway, make_gateway


class TestCancelSemantics:
    def test_404_is_order_gone_not_failure(self, monkeypatch):
        # MM phantom-open-quotes lesson (resolved 2026-08-12): an explicit
        # 404 means the order is already terminal — treating it as failure
        # wedges the ticker with a phantom resting entry forever.
        monkeypatch.setattr(auth_client, "api",
                            lambda *a, **k: (404, {}, None))
        assert LiveGateway().cancel("oid-1") is True

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
