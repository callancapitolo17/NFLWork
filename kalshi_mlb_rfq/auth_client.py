"""Backward-compat shim. Logic lives in kalshi_common.auth_client.

Configures the shared client from the taker's config at import time so existing
`auth_client.api(...)` callers in kalshi_mlb_rfq work unchanged.
"""
from kalshi_common import auth_client as _impl
from kalshi_mlb_rfq.config import (
    KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_BASE_URL, PROJECT_ROOT,
)

_impl.configure(KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_BASE_URL, PROJECT_ROOT)

api = _impl.api
