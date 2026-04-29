"""Fill / halt notifications: bot.log + optional webhook."""

import json
import urllib.error
import urllib.request
from datetime import datetime, timezone

from kalshi_mlb_rfq.config import LOG_PATH, NOTIFY_WEBHOOK_URL


def _append_log(line: str):
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a") as f:
        f.write(line + "\n")


def _post_webhook(payload: dict):
    if not NOTIFY_WEBHOOK_URL:
        return
    try:
        data = json.dumps(payload).encode()
        req = urllib.request.Request(NOTIFY_WEBHOOK_URL, data=data, method="POST")
        req.add_header("Content-Type", "application/json")
        with urllib.request.urlopen(req, timeout=5):
            pass
    except (urllib.error.URLError, TimeoutError, OSError):
        # Notification failure is non-fatal.
        pass


def fill(rfq_id: str, combo_market_ticker: str, contracts: float,
         price: float, ev_pct: float):
    ts = datetime.now(timezone.utc).isoformat()
    _append_log(f"{ts} [FILL] rfq={rfq_id} ticker={combo_market_ticker} "
                f"n={contracts} price={price} ev={ev_pct:.2%}")
    _post_webhook({"event": "fill", "rfq_id": rfq_id,
                    "ticker": combo_market_ticker, "contracts": contracts,
                    "price": price, "ev_pct": ev_pct, "ts": ts})


def halt(reason: str, detail: str | None = None):
    ts = datetime.now(timezone.utc).isoformat()
    _append_log(f"{ts} [HALT] reason={reason} detail={detail or ''}")
    _post_webhook({"event": "halt", "reason": reason,
                    "detail": detail, "ts": ts})
