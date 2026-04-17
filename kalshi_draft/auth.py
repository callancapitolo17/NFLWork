"""
Kalshi API Authentication
RSA-PSS SHA-256 signing for authenticated endpoints (portfolio, orders).
Extracted from kalshi_coaching.py for reuse across dashboards.
"""

import urllib.request
import json
import base64
import os
import time
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"

# Process-local throttle state for public_request.
# Kalshi's public rate limit is ~10 rps per IP; 600ms keeps us well under
# it while also protecting against bursts from the new portal pipeline.
_LAST_REQUEST_TS = 0.0
_MIN_INTERVAL = 0.6
_BACKOFFS = [1, 2, 4, 8, 16]


def load_credentials(env_dir=None):
    """Load Kalshi API credentials from .env file."""
    if env_dir is None:
        env_dir = Path(__file__).parent
    env_path = env_dir / ".env"
    creds = {"api_key_id": None, "private_key_path": None}

    if env_path.exists():
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("KALSHI_API_KEY_ID="):
                    creds["api_key_id"] = line.split("=", 1)[1].strip()
                elif line.startswith("KALSHI_PRIVATE_KEY_PATH="):
                    creds["private_key_path"] = line.split("=", 1)[1].strip()

    return creds


def sign_request(private_key_path: str, timestamp_str: str, method: str, path: str) -> str:
    """Sign a request using RSA-PSS with SHA-256 (Kalshi's required format)."""
    try:
        from cryptography.hazmat.primitives import hashes, serialization
        from cryptography.hazmat.primitives.asymmetric import padding
    except ImportError:
        print("  cryptography package required for authenticated requests")
        print("  Install with: pip install cryptography")
        return None

    with open(private_key_path, "rb") as f:
        private_key = serialization.load_pem_private_key(f.read(), password=None)

    path_without_query = path.split("?")[0]
    message = f"{timestamp_str}{method}{path_without_query}".encode()
    signature = private_key.sign(
        message,
        padding.PSS(
            mgf=padding.MGF1(hashes.SHA256()),
            salt_length=padding.PSS.DIGEST_LENGTH
        ),
        hashes.SHA256()
    )
    return base64.b64encode(signature).decode()


def authenticated_request(method: str, path: str, api_key_id: str, private_key_path: str):
    """Make an authenticated request to Kalshi API."""
    timestamp = datetime.now(timezone.utc)
    timestamp_str = str(int(timestamp.timestamp() * 1000))

    full_path = f"/trade-api/v2{path}"
    signature = sign_request(private_key_path, timestamp_str, method, full_path)
    if not signature:
        return None

    url = f"{BASE_URL}{path}"
    req = urllib.request.Request(url, method=method)
    req.add_header("KALSHI-ACCESS-KEY", api_key_id)
    req.add_header("KALSHI-ACCESS-SIGNATURE", signature)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", timestamp_str)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        print(f"  Auth request failed: {e.code} - {e.read().decode()}")
        return None
    except Exception as e:
        print(f"  Auth request error: {e}")
        return None


def public_request(path: str):
    """Make a public (unauthenticated) request to Kalshi API. Throttled.

    Enforces a minimum 600ms gap between successive calls (process-local)
    and retries on HTTP 429 with exponential backoff (1s, 2s, 4s, 8s, 16s).
    """
    global _LAST_REQUEST_TS
    elapsed = time.time() - _LAST_REQUEST_TS
    if elapsed < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - elapsed)
    url = f"{BASE_URL}{path}"
    for attempt, backoff in enumerate([0] + _BACKOFFS):
        if backoff:
            time.sleep(backoff)
        try:
            req = urllib.request.Request(url)
            req.add_header("Content-Type", "application/json")
            with urllib.request.urlopen(req) as response:
                _LAST_REQUEST_TS = time.time()
                return json.loads(response.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < len(_BACKOFFS):
                continue
            print(f"  Request failed: {e.code} - {e.read().decode()}")
            return None
        except Exception as e:
            print(f"  Request error: {e}")
            return None
