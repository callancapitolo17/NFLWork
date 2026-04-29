"""Authenticated HTTP wrapper around Kalshi /trade-api/v2.

Borrows the signing helper from kalshi_draft/auth.py (already in the main repo).
"""

import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

from kalshi_mlb_rfq.config import (
    KALSHI_API_KEY_ID,
    KALSHI_BASE_URL,
    KALSHI_PRIVATE_KEY_PATH,
    PROJECT_ROOT,
)

# kalshi_draft is gitignored; we have to find it on disk. It lives at the repo root,
# but PROJECT_ROOT differs between main checkout (= NFLWork) and a worktree
# (= NFLWork/.worktrees/<name>). Try both.
def _find_kalshi_draft() -> "str | None":
    candidates = [
        PROJECT_ROOT / "kalshi_draft",                # main repo
        PROJECT_ROOT.parent.parent / "kalshi_draft",  # inside a worktree
    ]
    for c in candidates:
        if (c / "auth.py").exists():
            return str(c)
    return None


_kalshi_draft_dir = _find_kalshi_draft()
if _kalshi_draft_dir and _kalshi_draft_dir not in sys.path:
    sys.path.insert(0, _kalshi_draft_dir)
try:
    from auth import sign_request as _sign_request
except ImportError:
    # In tests we patch _sign anyway; provide a stub so import doesn't fail.
    def _sign_request(_pk, _ts, _method, _path):
        raise RuntimeError("sign_request unavailable; set KALSHI_PRIVATE_KEY_PATH and "
                           "ensure kalshi_draft/auth.py is reachable")


def _sign(method: str, path_no_query: str) -> tuple[str, str]:
    """Returns (signature, timestamp_ms_str)."""
    ts = str(int(datetime.now(timezone.utc).timestamp() * 1000))
    sig = _sign_request(KALSHI_PRIVATE_KEY_PATH, ts, method,
                        f"/trade-api/v2{path_no_query}")
    return sig, ts


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
    """
    path_no_query = path.split("?", 1)[0]
    sig, ts = _sign(method, path_no_query)

    url = f"{KALSHI_BASE_URL}{path}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("KALSHI-ACCESS-KEY", KALSHI_API_KEY_ID or "")
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
