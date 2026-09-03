"""DraftKings price sidecar — a real Chrome that prices SGPs for the bots.

WHY THIS EXISTS (issue #102). Since ~2026-08-20 Akamai gates
``POST */api/wager/v1/calculateBets`` on two things: the CLIENT's TLS/HTTP2
fingerprint (curl_cffi fails on every profile; a real, non-headless Chrome
passes) and BURST SHAPE (concurrent calls are denied; sequential ones pass).
Nothing at the HTTP layer fixes that, so DK's price call has to be made by a
real browser, one call at a time. This process owns that browser so the bots
never have to: they POST ``{"selections": [...]}`` to ``/price`` on localhost
and get DK's answer back, exactly as if DK were an ordinary HTTP book.

Inputs:  HTTP on ``DK_SIDECAR_PORT`` (default 8095), loopback only.
Outputs: ``POST /price`` -> DK's own status (200 priced / 422 declined /
         403 blocked) with ``{"true_odds", "display", "restrictions"}``;
         ``GET /health`` -> liveness + counters.
Side effects: launches a persistent Chrome profile under
``DK_SIDECAR_PROFILE_DIR``; makes one DK request per /price call plus a
periodic page reload. No DuckDB, no writes anywhere else.

Serialization is BY CONSTRUCTION: ``http.server.HTTPServer`` is
single-threaded, so concurrent bot requests queue in the socket backlog and
reach DK one at a time, ``DK_SIDECAR_MIN_INTERVAL_SEC`` apart. That is the
no-burst rule enforced somewhere the bots cannot bypass it.

Run with a Python that has Playwright + Google Chrome installed:
    dk_price_sidecar/run.sh
"""
from __future__ import annotations

import json
import logging
import os
import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

logger = logging.getLogger("dk_price_sidecar")

# ---- configuration (env, read once at import; no hidden state below) ------
PORT = int(os.environ.get("DK_SIDECAR_PORT", "8095"))
PROFILE_DIR = Path(os.environ.get(
    "DK_SIDECAR_PROFILE_DIR", str(Path.home() / ".dk_price_sidecar" / "profile")))
# Measured 2026-09-02: 1.5s spacing -> 100% pass; 4 concurrent -> 5%.
MIN_INTERVAL_SEC = float(os.environ.get("DK_SIDECAR_MIN_INTERVAL_SEC", "1.0"))
# Reload the sportsbook page so Akamai's cookies stay fresh.
PAGE_RELOAD_SEC = float(os.environ.get("DK_SIDECAR_PAGE_RELOAD_SEC", "600"))
MINIMIZE_WINDOW = os.environ.get("DK_SIDECAR_MINIMIZE", "1") != "0"

SPORTSBOOK_URL = "https://sportsbook.draftkings.com/leagues/baseball/mlb"
# wagerBaseApiHost from DK's own client config, bare path (no /en/) — the
# form dkBetSlip.js builds. See mlb_sgp/README.md "DraftKings price host".
CALCULATE_BETS_URL = "https://gaming-us-wv.draftkings.com/api/wager/v1/calculateBets"
# DK's betslip header set, captured 2026-09-02. The request must originate
# cross-site from the sportsbook page with these; the browser adds
# origin/referer/sec-* itself.
DK_HEADERS = {
    "accept": "application/json",
    "content-type": "application/json",
    "clienttype": "Website",
    "x-api-features": json.dumps({"EnableFullSGPDrivenFlow": True}),
    "x-client-name": "web",
    "x-client-feature": "betslip",
    "x-client-page": "league",
    "x-client-version": "2636.2.1.11",
    "x-client-widget-name": "betslip",
    "x-client-widget-version": "2627.2.1",
}

# Runs INSIDE the page. Returns DK's status and body; a thrown fetch (the
# deterministic first-call throw on a fresh page, or a cross-origin 403 with
# no CORS headers) comes back as status 0 so the caller can tell it apart.
_FETCH_JS = """async ({url, headers, body}) => {
    try {
        const r = await fetch(url, {
            method: "POST",
            headers: {...headers, "x-request-client-timestamp": String(Date.now())},
            credentials: "include",
            body: JSON.stringify(body),
        });
        return {status: r.status, text: await r.text()};
    } catch (e) {
        return {status: 0, text: String(e)};
    }
}"""


def build_calculate_bets_body(selection_ids: list[str]) -> dict:
    """The production body mlb_sgp.draftkings sends — verified accepted
    unchanged on 2026-09-02 (YourBet trueOdds=7.5)."""
    return {
        "selections": [],
        "selectionsForYourBet": [{"id": sid, "yourBetGroup": 0}
                                 for sid in selection_ids],
        "selectionsForCombinator": [],
        "selectionsForProgressiveParlay": [],
        "oddsStyle": "american",
    }


def parse_calculate_bets(status: int, text: str, n_legs: int) -> dict:
    """Map DK's raw answer to the sidecar's response body.

    ``true_odds`` is the correlated decimal for the FULL leg set (the
    ``YourBet`` bet whose ``selectionsMapped`` covers every leg), or None.
    A 200 whose bets are only singles, or that carries
    ``combinabilityRestrictions``, is a decline, not a price — the same rule
    ``mlb_sgp.draftkings.price_selection_set`` applies.
    """
    out = {"status": status, "true_odds": None, "display": None,
           "restrictions": None}
    if status != 200:
        out["error"] = text[:200]
        return out
    try:
        data = json.loads(text)
    except ValueError:
        out["error"] = "non-JSON body"
        return out
    restrictions = data.get("combinabilityRestrictions") or None
    out["restrictions"] = restrictions
    if restrictions:
        return out
    for bet in data.get("bets", []):
        mapped = bet.get("selectionsMapped") or []
        if bet.get("trueOdds") and len(mapped) >= n_legs:
            out["true_odds"] = float(bet["trueOdds"])
            out["display"] = bet.get("displayOdds")
            break
    return out


class DkBrowser:
    """One persistent, non-headless Chrome that makes DK's own price call.

    ``price()`` is the only DK-facing method. It enforces the inter-call
    interval, reloads the sportsbook page on a timer, and relaunches the
    browser after any Playwright-level failure (a crashed Chrome must never
    take the sidecar down with it — the NEXT call relaunches).
    """

    def __init__(self):
        self._pw = None
        self._ctx = None
        self._page = None
        self._last_call_at = 0.0
        self._page_loaded_at = 0.0
        self.calls = 0
        self.last_status: int | None = None
        self.relaunches = 0

    # -- lifecycle ---------------------------------------------------------
    def launch(self) -> None:
        from playwright.sync_api import sync_playwright   # lazy: tests run without it
        self.close()
        PROFILE_DIR.mkdir(parents=True, exist_ok=True)
        self._pw = sync_playwright().start()
        # NOT headless: that is exactly what Akamai's rule detects. The
        # window is minimized instead (verified 3/3 on 2026-09-02).
        self._ctx = self._pw.chromium.launch_persistent_context(
            str(PROFILE_DIR), headless=False, channel="chrome",
            args=["--disable-blink-features=AutomationControlled"],
            ignore_default_args=["--enable-automation"],
            viewport={"width": 1280, "height": 800})
        self._page = self._ctx.pages[0] if self._ctx.pages else self._ctx.new_page()
        if MINIMIZE_WINDOW:
            try:
                cdp = self._ctx.new_cdp_session(self._page)
                wid = cdp.send("Browser.getWindowForTarget")["windowId"]
                cdp.send("Browser.setWindowBounds",
                         {"windowId": wid, "bounds": {"windowState": "minimized"}})
            except Exception as e:                       # cosmetic, never fatal
                logger.warning("could not minimize window: %s", e)
        self._load_page()
        # The first fetch on a fresh page throws deterministically; absorb it
        # here so no bot request pays for it.
        self._page.evaluate(_FETCH_JS, {"url": CALCULATE_BETS_URL,
                                        "headers": DK_HEADERS,
                                        "body": build_calculate_bets_body([])})
        logger.info("browser ready (profile=%s)", PROFILE_DIR)

    def _load_page(self) -> None:
        self._page.goto(SPORTSBOOK_URL, wait_until="domcontentloaded", timeout=60000)
        self._page.wait_for_timeout(5000)
        self._page_loaded_at = time.monotonic()

    def close(self) -> None:
        for closer in ((self._ctx.close if self._ctx else None),
                       (self._pw.stop if self._pw else None)):
            try:
                if closer:
                    closer()
            except Exception:
                pass
        self._pw = self._ctx = self._page = None

    def ready(self) -> bool:
        return self._page is not None

    # -- the one DK-facing call -------------------------------------------
    def price(self, selection_ids: list[str]) -> dict:
        if not self.ready():
            self.launch()
        if time.monotonic() - self._page_loaded_at > PAGE_RELOAD_SEC:
            self._load_page()
        wait = MIN_INTERVAL_SEC - (time.monotonic() - self._last_call_at)
        if wait > 0:
            time.sleep(wait)
        try:
            raw = self._page.evaluate(_FETCH_JS, {
                "url": CALCULATE_BETS_URL, "headers": DK_HEADERS,
                "body": build_calculate_bets_body(selection_ids)})
        except Exception as e:
            # Playwright-level failure = the browser is gone. Relaunch on the
            # next call; report THIS one as unavailable.
            logger.error("browser call failed, will relaunch: %s", e)
            self.relaunches += 1
            self.close()
            return {"status": 503, "true_odds": None, "display": None,
                    "restrictions": None, "error": f"browser: {e}"[:200]}
        finally:
            self._last_call_at = time.monotonic()
            self.calls += 1
        self.last_status = raw["status"]
        return parse_calculate_bets(raw["status"], raw["text"], len(selection_ids))


class _Handler(BaseHTTPRequestHandler):
    browser: DkBrowser = None    # set by serve()
    started_at = time.monotonic()

    def _send(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path != "/health":
            return self._send(404, {"error": "not found"})
        b = self.browser
        self._send(200, {"ok": b.ready(), "calls": b.calls,
                         "last_status": b.last_status, "relaunches": b.relaunches,
                         "uptime_sec": round(time.monotonic() - self.started_at)})

    def do_POST(self):
        if self.path != "/price":
            return self._send(404, {"error": "not found"})
        try:
            n = int(self.headers.get("Content-Length", "0"))
            req = json.loads(self.rfile.read(n) or b"{}")
            ids = req.get("selections")
            if not isinstance(ids, list) or not ids or not all(isinstance(s, str) for s in ids):
                return self._send(400, {"error": "selections must be a non-empty list of ids"})
        except (ValueError, TypeError) as e:
            return self._send(400, {"error": f"bad request: {e}"})
        result = self.browser.price(ids)
        # The sidecar's HTTP status IS DK's status (0 = fetch threw -> 502),
        # so the bot's existing 200/422/403 handling applies unchanged.
        http_status = result["status"] if result["status"] >= 100 else 502
        self._send(http_status, result)

    def log_message(self, fmt, *args):            # route to logging, not stderr
        logger.debug("%s " + fmt, self.address_string(), *args)


def serve() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    browser = DkBrowser()
    browser.launch()
    _Handler.browser = browser
    # Loopback only: this process speaks for a real browser session and
    # must never be reachable off the machine.
    httpd = HTTPServer(("127.0.0.1", PORT), _Handler)
    logger.info("dk_price_sidecar listening on http://127.0.0.1:%d "
                "(min_interval=%.1fs, reload=%.0fs)", PORT, MIN_INTERVAL_SEC, PAGE_RELOAD_SEC)
    try:
        httpd.serve_forever()
    finally:
        browser.close()


if __name__ == "__main__":
    serve()
