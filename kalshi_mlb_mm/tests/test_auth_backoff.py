"""N6 — exponential backoff on transient 429/503 from Kalshi.

The shared kalshi_common.auth_client.api() retries 429/503 up to 3 times with
exponential backoff (200/400/800ms) + jitter. All other status codes return
immediately on the first attempt, so the taker's behavior is unchanged.
"""
import io
import json
from unittest.mock import patch

import urllib.error

from kalshi_common import auth_client


def _configure_stub():
    """Inject a stub signer so api() doesn't try to read a real key."""
    auth_client._API_KEY_ID = "stub-key"
    auth_client._PRIVATE_KEY_PATH = "/tmp/none"
    auth_client._BASE_URL = "https://example.invalid/trade-api/v2"
    auth_client._sign_request = lambda _pk, _ts, _m, _p: "fake-sig"


class _FakeResponse:
    def __init__(self, status: int, body: dict):
        self.status = status
        self._payload = json.dumps(body).encode()
        self.headers = {}

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _http_error(status: int, body: dict) -> urllib.error.HTTPError:
    payload = json.dumps(body).encode()
    return urllib.error.HTTPError(
        url="https://example.invalid",
        code=status,
        msg="x",
        hdrs=None,
        fp=io.BytesIO(payload),
    )


def test_api_retries_429_then_returns_200(monkeypatch):
    """Two 429 responses, then a 200. api() must retry and return the 200 body."""
    _configure_stub()
    sleeps = []
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: sleeps.append(s))

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise _http_error(429, {"error": "rate_limited"})
        return _FakeResponse(200, {"ok": True})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 200
    assert isinstance(body, dict) and body.get("ok") is True
    assert calls["n"] == 3, "expected exactly 3 attempts (2 retries + success)"
    # Two backoff sleeps for two retries. Each >= 200ms.
    assert len(sleeps) == 2
    assert sleeps[0] >= 0.2 and sleeps[1] >= 0.4


def test_api_retries_503_then_returns_200(monkeypatch):
    """503 is also retried."""
    _configure_stub()
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: None)

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        if calls["n"] == 1:
            raise _http_error(503, {"error": "unavailable"})
        return _FakeResponse(200, {"ok": True})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 200
    assert body == {"ok": True}
    assert calls["n"] == 2


def test_api_returns_immediately_on_200(monkeypatch):
    """A 200 first response means no retries and no sleeps."""
    _configure_stub()
    sleeps = []
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: sleeps.append(s))

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        return _FakeResponse(200, {"ok": True})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 200 and body == {"ok": True}
    assert calls["n"] == 1
    assert sleeps == [], "no sleep should happen on first-attempt success"


def test_api_does_not_retry_400(monkeypatch):
    """4xx other than 429 must NOT be retried — they are real client errors."""
    _configure_stub()
    sleeps = []
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: sleeps.append(s))

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        raise _http_error(400, {"error": "bad_request"})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 400
    assert calls["n"] == 1, "400 must not be retried"
    assert sleeps == []


def test_api_returns_last_attempt_when_all_retries_429(monkeypatch):
    """If every attempt returns 429, api() returns the final 429 (caller can handle)."""
    _configure_stub()
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: None)

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        raise _http_error(429, {"error": "rate_limited"})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 429
    assert calls["n"] == 4, "expected exactly 4 attempts (1 initial + 3 retries)"
