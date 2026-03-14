"""
Kalshi API order management.
Handles authenticated requests for placing, amending, and cancelling orders.
"""

import sys
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# Import sign_request from kalshi_draft/auth.py
# kalshi_draft is untracked, so look in the main repo via git
import subprocess as _sp
_repo_root = Path(_sp.check_output(
    ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
    cwd=Path(__file__).parent, text=True
).strip()).parent
_kalshi_draft = _repo_root / "kalshi_draft"
sys.path.insert(0, str(_kalshi_draft))
from auth import sign_request

from config import KALSHI_BASE_URL, KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH


def _authenticated_request(method, path, body=None):
    """Make an authenticated request to Kalshi API with optional JSON body."""
    if not KALSHI_API_KEY_ID or not KALSHI_PRIVATE_KEY_PATH:
        print("Error: Kalshi credentials not configured. Check .env file.")
        return None

    timestamp = datetime.now(timezone.utc)
    timestamp_str = str(int(timestamp.timestamp() * 1000))

    full_path = f"/trade-api/v2{path}"
    signature = sign_request(KALSHI_PRIVATE_KEY_PATH, timestamp_str, method, full_path)
    if not signature:
        return None

    url = f"{KALSHI_BASE_URL}{path}"

    if body:
        data = json.dumps(body).encode()
        req = urllib.request.Request(url, data=data, method=method)
    else:
        req = urllib.request.Request(url, method=method)

    req.add_header("KALSHI-ACCESS-KEY", KALSHI_API_KEY_ID)
    req.add_header("KALSHI-ACCESS-SIGNATURE", signature)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", timestamp_str)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        body_text = e.read().decode()
        print(f"  API error {e.code}: {body_text}")
        return None
    except Exception as e:
        print(f"  Request error: {e}")
        return None


def place_order(ticker, side, price, count, post_only=True):
    """Place a new limit order on Kalshi.

    Args:
        ticker: Market ticker (e.g., "KXNCAAMB1HSPREAD-26MAR10-DUKE-3")
        side: "yes" or "no"
        price: Price in cents (1-99)
        count: Number of contracts
        post_only: If True, reject if order would immediately fill

    Returns:
        Order response dict or None on failure.
    """
    body = {
        "ticker": ticker,
        "action": "buy",
        "side": side,
        "type": "limit",
        "count": count,
        "time_in_force": "good_till_canceled",
    }

    if side == "yes":
        body["yes_price"] = price
    else:
        body["no_price"] = price

    if post_only:
        body["post_only"] = True

    result = _authenticated_request("POST", "/portfolio/orders", body=body)
    if result and "order" in result:
        order = result["order"]
        print(f"  Placed {side} {count}x @ {price}c on {ticker} → order_id={order.get('order_id', '?')}")
        return order
    return None


def amend_order(order_id, side, price=None, count=None):
    """Amend a resting order's price and/or quantity.

    Args:
        order_id: The order to amend
        side: "yes" or "no" — determines which price field to send
        price: New price in cents (optional)
        count: New count (optional)

    Returns:
        Amended order dict or None.
    """
    body = {}
    if price is not None:
        if side == "yes":
            body["yes_price"] = price
        else:
            body["no_price"] = price
    if count is not None:
        body["count"] = count

    if not body:
        return None

    result = _authenticated_request("POST", f"/portfolio/orders/{order_id}/amend", body=body)
    if result and "order" in result:
        print(f"  Amended order {order_id}: {side}_price={price}, count={count}")
        return result["order"]
    return None


def cancel_order(order_id):
    """Cancel a single resting order."""
    result = _authenticated_request("DELETE", f"/portfolio/orders/{order_id}")
    if result is not None:
        print(f"  Cancelled order {order_id}")
        return True
    return False


def batch_cancel(order_ids):
    """Cancel multiple orders in one API call (0.2 rate limit units each).

    Args:
        order_ids: List of order IDs to cancel (max 20)
    """
    if not order_ids:
        return True

    # Batch cancel up to 20 at a time
    for i in range(0, len(order_ids), 20):
        batch = order_ids[i:i+20]
        body = {"ids": batch}
        result = _authenticated_request("DELETE", "/portfolio/orders/batched", body=body)
        if result is None:
            print(f"  Warning: batch cancel failed for {len(batch)} orders")
            return False
        print(f"  Batch cancelled {len(batch)} orders")
    return True


def get_resting_orders(ticker=None):
    """Get all resting orders, optionally filtered by ticker."""
    path = "/portfolio/orders?status=resting&limit=200"
    if ticker:
        path += f"&ticker={ticker}"

    result = _authenticated_request("GET", path)
    if result and "orders" in result:
        return result["orders"]
    return []


def get_positions():
    """Get all portfolio positions."""
    result = _authenticated_request("GET", "/portfolio/positions?limit=200")
    if result and "market_positions" in result:
        return result["market_positions"]
    return []


def cancel_all_orders():
    """Emergency: cancel ALL resting orders."""
    orders = get_resting_orders()
    if not orders:
        print("  No resting orders to cancel.")
        return True

    order_ids = [o["order_id"] for o in orders]
    print(f"  KILL SWITCH: Cancelling {len(order_ids)} orders...")
    return batch_cancel(order_ids)
