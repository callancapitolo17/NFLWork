"""Issue #102 — the sidecar transport for DK's price call.

``DK_PRICE_TRANSPORT=sidecar`` routes ``draftkings.price_selection_set``
through dk_price_sidecar instead of curl_cffi. The contract must be identical
to the HTTP path: float on a priced set, None on a decline, and a COUNTED
transport error when DK (or the sidecar itself) is refusing us. Default
transport stays "http" so nothing changes until the operator flips it.

No network: ``_sidecar_post`` is the seam.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp import draftkings                                   # noqa: E402
from mlb_sgp._shared import BookTransportError, FetchCounters    # noqa: E402

_REFS = ["0HC1N150_1", "0OU1O750_1"]


@pytest.fixture
def sidecar(monkeypatch):
    """Route through the sidecar and let each test script its reply."""
    monkeypatch.setattr(draftkings, "DK_PRICE_TRANSPORT", "sidecar")
    state = {"reply": None, "calls": 0}

    def fake_post(refs, *, timeout):
        state["calls"] += 1
        state["refs"] = list(refs)
        reply = state["reply"]
        if isinstance(reply, Exception):
            raise reply
        return reply

    monkeypatch.setattr(draftkings, "_sidecar_post", fake_post)
    return state


def _price(refs=_REFS):
    counters = FetchCounters("draftkings", "on_demand")
    counters.bump("targets_attempted")
    result = draftkings.price_selection_set(object(), refs, counters=counters)
    return result, counters.snapshot()


def test_default_transport_is_http_and_never_touches_the_sidecar(monkeypatch):
    called = []
    monkeypatch.setattr(draftkings, "DK_PRICE_TRANSPORT", "http")
    monkeypatch.setattr(draftkings, "_sidecar_post",
                        lambda *a, **k: called.append(1))

    class _Resp:
        status_code = 422
        text = ""
        def json(self): return {}

    class _Client:
        class session:
            @staticmethod
            def post(*a, **k): return _Resp()

    draftkings.price_selection_set(_Client(), _REFS)
    assert called == [], "http transport must not call the sidecar"


def test_priced_set_returns_the_correlated_decimal(sidecar):
    sidecar["reply"] = (200, {"status": 200, "true_odds": 7.5, "display": "+650"})
    result, snap = _price()
    assert result == 7.5
    assert snap.transport_errors == 0
    assert sidecar["refs"] == _REFS


def test_200_without_a_price_is_a_decline(sidecar):
    """Singles-only or restricted response: DK answered, combo not built."""
    sidecar["reply"] = (200, {"status": 200, "true_odds": None,
                              "restrictions": [{"restrictionType": "NonCombinableGroup"}]})
    result, snap = _price()
    assert result is None
    assert snap.transport_errors == 0


def test_422_is_a_decline_not_a_transport_error(sidecar):
    sidecar["reply"] = (422, {"status": 422, "error": "NonCombinable"})
    result, snap = _price()
    assert result is None
    assert snap.transport_errors == 0


def test_403_from_dk_via_sidecar_counts_a_transport_error(sidecar):
    """The whole point of #102's first fix, preserved on the new transport."""
    sidecar["reply"] = (403, {"status": 403, "error": "Access Denied"})
    result, snap = _price()
    assert result is None
    assert snap.transport_errors == 1


def test_sidecar_browser_failure_counts_a_transport_error(sidecar):
    """503 = the sidecar's Chrome is down. DK is unavailable, not declining."""
    sidecar["reply"] = (503, {"status": 503, "error": "browser: crashed"})
    result, snap = _price()
    assert result is None
    assert snap.transport_errors == 1


def test_sidecar_unreachable_counts_a_transport_error(sidecar):
    sidecar["reply"] = BookTransportError("draftkings", "price",
                                          detail="dk_price_sidecar unreachable")
    result, snap = _price()
    assert result is None
    assert snap.transport_errors == 1


def test_sidecar_path_never_raises_to_the_caller(sidecar):
    sidecar["reply"] = (403, {})
    assert draftkings.price_selection_set(object(), _REFS) is None
