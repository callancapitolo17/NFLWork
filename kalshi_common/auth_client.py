"""Authenticated HTTP wrapper around Kalshi /trade-api/v2.

Config-agnostic: callers inject credentials via configure(). Defaults read from
env so ad-hoc scripts work without a configure() call.
"""
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

_API_KEY_ID = os.environ.get("KALSHI_API_KEY_ID")
_PRIVATE_KEY_PATH = os.environ.get("KALSHI_PRIVATE_KEY_PATH")
_BASE_URL = os.environ.get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
_sign_request = None


def configure(api_key_id, private_key_path, base_url, project_root):
    """Inject per-bot credentials + the repo root used to locate kalshi_draft/auth.py."""
    global _API_KEY_ID, _PRIVATE_KEY_PATH, _BASE_URL, _sign_request
    _API_KEY_ID = api_key_id
    _PRIVATE_KEY_PATH = private_key_path
    _BASE_URL = base_url
    candidates = [Path(project_root) / "kalshi_draft",
                  Path(project_root).parent.parent / "kalshi_draft"]
    for c in candidates:
        if (c / "auth.py").exists():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            break
    try:
        from auth import sign_request as sr
        _sign_request = sr
    except ImportError:
        # In tests we patch _sign_request anyway; provide a stub so import doesn't fail.
        def _stub(_pk, _ts, _method, _path):
            raise RuntimeError("sign_request unavailable; set KALSHI_PRIVATE_KEY_PATH and "
                               "ensure kalshi_draft/auth.py is reachable")
        _sign_request = _stub


def _sign(method: str, path_no_query: str) -> tuple[str, str]:
    """Returns (signature, timestamp_ms_str)."""
    if _sign_request is None:
        raise RuntimeError("auth_client.configure() not called")
    ts = str(int(datetime.now(timezone.utc).timestamp() * 1000))
    sig = _sign_request(_PRIVATE_KEY_PATH, ts, method, f"/trade-api/v2{path_no_query}")
    return sig, ts


_RETRY_STATUSES = (429, 503)
_MAX_RETRIES = 3  # 4 total attempts max
# N9: only retry on methods that are safe to repeat without side effects.
# POST is excluded — the first 429/503 may have already reached Kalshi and
# created a resource (e.g. a quote); retrying would produce a duplicate.
_IDEMPOTENT_METHODS = ("GET", "PUT", "DELETE")


def _attempt(method: str, path: str, body: dict | None,
             timeout: int) -> tuple[int, dict | str, dict]:
    """Single signed request; returns (status, body_or_text, headers)."""
    path_no_query = path.split("?", 1)[0]
    sig, ts = _sign(method, path_no_query)
    url = f"{_BASE_URL}{path}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("KALSHI-ACCESS-KEY", _API_KEY_ID or "")
    req.add_header("KALSHI-ACCESS-SIGNATURE", sig)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", ts)
    req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            text = r.read().decode()
            headers = dict(r.headers)
            if not text:
                return r.status, {}, headers
            try:
                return r.status, json.loads(text), headers
            except json.JSONDecodeError:
                return r.status, text, headers
    except urllib.error.HTTPError as e:
        text = e.read().decode()
        headers = dict(e.headers) if e.headers else {}
        try:
            return e.code, json.loads(text), headers
        except (json.JSONDecodeError, ValueError):
            return e.code, text, headers


def api(method: str, path: str, body: dict | None = None,
        timeout: int = 30) -> tuple[int, dict | str, dict]:
    """Make an authenticated request.

    Args:
        method: HTTP verb.
        path: Path under /trade-api/v2 (e.g., "/exchange/status"). Query string allowed;
              only the path portion is signed.
        body: Optional JSON body.
        timeout: Seconds.

    Returns:
        (status_code, parsed_body_or_text, response_headers_dict)

    N6: transient 429/503 errors are retried up to 3 times with exponential
    backoff + jitter (200/400/800ms + uniform(0, 100ms)). Each retry re-signs
    with a fresh timestamp. All other status codes return on the first attempt.

    N9: retry is ONLY performed for idempotent methods (GET, PUT, DELETE).
    POST returns on the first attempt regardless of status — a 429/503 on a
    POST may have already reached Kalshi and created the resource (e.g. a
    quote), so retrying would produce a duplicate.
    """
    last = _attempt(method, path, body, timeout)
    if method.upper() not in _IDEMPOTENT_METHODS:
        return last  # N9: non-idempotent — never retry
    for attempt in range(_MAX_RETRIES):
        status = last[0]
        if status not in _RETRY_STATUSES:
            return last
        # Backoff before next attempt: 200ms * 2^attempt + uniform(0, 100ms).
        time.sleep(0.2 * (2 ** attempt) + random.uniform(0, 0.1))
        last = _attempt(method, path, body, timeout)
    return last
