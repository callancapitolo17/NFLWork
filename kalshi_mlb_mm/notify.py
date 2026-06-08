"""Minimal alerting: logger + optional webhook (mirrors taker notify)."""
import json
import logging
import urllib.request

from kalshi_mlb_mm.config import NOTIFY_WEBHOOK_URL

log = logging.getLogger("kalshi_mlb_mm.notify")


def _post(text: str):
    log.info(text)
    if not NOTIFY_WEBHOOK_URL:
        return
    try:
        req = urllib.request.Request(
            NOTIFY_WEBHOOK_URL, data=json.dumps({"text": text}).encode(),
            headers={"Content-Type": "application/json"}, method="POST")
        urllib.request.urlopen(req, timeout=5)
    except Exception:
        pass


def halt(reason: str, detail: str = ""):
    _post(f"[MM HALT] {reason} {detail}")


def fill(combo_market_ticker: str, side: str, contracts: int, price: float):
    _post(f"[MM FILL] {combo_market_ticker} {side} x{contracts} @ {price:.3f}")


def resume(reason: str, detail: str = ""):
    _post(f"[MM RESUME] {reason} {detail}")
