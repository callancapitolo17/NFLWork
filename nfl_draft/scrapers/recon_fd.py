#!/usr/bin/env python3
"""FanDuel NFL Draft markets recon.

Goal: capture the raw JSON that powers FD's NFL Draft futures page and save it
to nfl_draft/tests/fixtures/fanduel/draft_markets.json so the parser subagent
can work offline.

Key finding from live recon (2026-04-17):
    FD puts NFL Draft under the NFL content-managed-page as tab 391
    ('NFL Draft'). The single endpoint:

        GET https://api.sportsbook.fanduel.com/sbapi/content-managed-page
            ?page=CUSTOM&customPageId=nfl
            &_ak=<api_key>&timezone=America%2FNew_York

    returns ~700KB of JSON: 69 markets, 3 events, with `layout.tabs['391']`
    listing which markets (via cards->coupons) belong to the NFL Draft tab.

    The 3 FD required headers (x-application, x-sportsbook-region,
    x-px-context) from mlb_sgp/scraper_fanduel_sgp.py still apply - without
    them FD returns empty bodies.

Strategy:

Phase 1 - REST against content-managed-page (customPageId=nfl):
    Pull the whole NFL futures bundle. The parser will filter to tab 391.

Phase 2 - Headed Playwright fallback:
    Same-site, but useful when the PX token is stale (user logs in or waits
    for FD to issue fresh tokens).

Usage:
    python nfl_draft/scrapers/recon_fd.py              # REST
    python nfl_draft/scrapers/recon_fd.py --browser    # force browser
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from nfl_draft.scrapers._recon_util import (
    ensure_fixture_dirs,
    print_diagnostics,
    save_fixture,
)

from curl_cffi import requests as cffi_requests

# ---------------------------------------------------------------------------
# FD config - mirrors mlb_sgp/scraper_fanduel_sgp.py
# ---------------------------------------------------------------------------

FD_AK = "FhMFpcPWXMeyZxOx"
# PerimeterX visitor token from the MLB scraper. If FD starts returning 400s
# or empty bodies, refresh by opening sportsbook.fanduel.com in Chrome DevTools,
# find any event-page request, copy its x-px-context header value here.
FD_PX_CONTEXT = (
    "_pxvid=687ed2ad-33ac-11f1-ba08-95fbba040e2b;"
    "pxcts=6119da3c-33ae-11f1-b93d-4b2b160f34a7;"
)

FD_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/146.0.0.0 Safari/537.36"
    ),
    "Accept": "application/json",
    "Referer": "https://sportsbook.fanduel.com/",
    "Origin": "https://sportsbook.fanduel.com",
    "x-application": FD_AK,
    "x-sportsbook-region": "NJ",
    "x-px-context": FD_PX_CONTEXT,
}

# Confirmed endpoint: the content-managed-page for the whole NFL page includes
# NFL Draft markets as tab 391 ('NFL Draft'). FD doesn't expose a draft-only
# customPageId slug - the NFL bundle is the one source of truth.
FD_NFL_PAGE_URL = (
    "https://api.sportsbook.fanduel.com/sbapi/content-managed-page"
    f"?page=CUSTOM&customPageId=nfl&_ak={FD_AK}&timezone=America%2FNew_York"
)
# tab id 391 = "NFL Draft" tab under the NFL page (observed 2026-04-17).
# If FD renumbers, scan `layout.tabs` for title='NFL Draft'; the script does
# that discovery below so we don't hardcode 391 as a required value.
FD_NFL_DRAFT_TAB_TITLE = "NFL Draft"


def _find_draft_tab(data: dict) -> tuple[str | None, str | None]:
    """Scan layout.tabs for the NFL Draft tab. Returns (tab_id, title)."""
    if not isinstance(data, dict):
        return None, None
    tabs = (data.get("layout", {}) or {}).get("tabs", {}) or {}
    for tab_id, tab in tabs.items():
        title = (tab or {}).get("title", "") or ""
        if "draft" in title.lower():
            return tab_id, title
    return None, None


def run_rest_phase() -> tuple[dict | None, str | None, dict]:
    """Fetch the NFL content-managed-page; confirm it includes an NFL Draft tab."""
    print("Phase 1: REST fetch of FD NFL content-managed-page...")
    session = cffi_requests.Session(impersonate="chrome")

    try:
        resp = session.get(FD_NFL_PAGE_URL, headers=FD_HEADERS, timeout=20)
    except Exception as exc:
        print(f"  request failed: {exc!r}")
        return None, None, {}
    if resp.status_code != 200:
        print(f"  HTTP {resp.status_code}")
        if resp.status_code == 400:
            print("  If FD returns 400 with empty body, the x-px-context token")
            print("  in FD_PX_CONTEXT has likely expired. Refresh from a browser session.")
        return None, None, {}
    try:
        data = resp.json()
    except Exception:
        print(f"  non-JSON body ({len(resp.text)} chars)")
        return None, None, {}

    # Must have attachments with markets AND a tab titled 'NFL Draft'.
    attachments = data.get("attachments", {})
    markets = attachments.get("markets", {}) or {}
    if not markets:
        print(f"  bundle has no markets - FD may have rotated the draft off NFL page")
        return None, None, {}

    tab_id, tab_title = _find_draft_tab(data)
    if tab_id is None:
        print(f"  no 'NFL Draft' tab in layout - FD may have renamed it")
        print(f"  tabs present: {[t.get('title') for t in (data.get('layout', {}).get('tabs', {}) or {}).values()]}")
        return None, None, {}

    # Count how many markets belong to the draft tab for a sanity signal.
    draft_market_count = _count_draft_markets(data, tab_id)
    print(f"  found tab {tab_id} {tab_title!r} with {draft_market_count} "
          f"NFL-Draft markets (of {len(markets)} total on NFL page)")

    return data, FD_NFL_PAGE_URL, {
        "page_slug": "nfl",
        "draft_tab_id": tab_id,
        "draft_tab_title": tab_title,
        "draft_market_count": draft_market_count,
        "total_market_count": len(markets),
        "method": "rest_content_managed_page",
    }


def _count_draft_markets(data: dict, tab_id: str) -> int:
    """Count markets reachable from layout.tabs[tab_id].cards[*].coupons[*].marketIds.

    FD's page structure is:
        layout.tabs[id].cards[*] -> coupons[*] -> marketIds: [...]
    We walk that tree and count unique market IDs.
    """
    seen: set = set()
    tab = (data.get("layout", {}).get("tabs", {}).get(tab_id) or {})
    for card in tab.get("cards", []) or []:
        # Cards can reference coupons by id (lookup in layout.coupons) or have
        # them inline. Handle both.
        coupons = card.get("coupons", []) or []
        if coupons and isinstance(coupons[0], dict):
            coupon_iter = coupons
        else:
            coupon_map = (data.get("layout", {}).get("coupons", {}) or {})
            coupon_iter = [coupon_map.get(str(c)) for c in coupons if coupon_map.get(str(c))]
        for coup in coupon_iter:
            if not coup:
                continue
            for mid in coup.get("marketIds", []) or []:
                seen.add(str(mid))
    # Fallback: if the cards->coupons path yields nothing, count markets whose
    # name contains "draft".
    if not seen:
        markets = data.get("attachments", {}).get("markets", {}) or {}
        for mid, m in markets.items():
            if "draft" in (m.get("marketName") or "").lower():
                seen.add(mid)
    return len(seen)


# ---------------------------------------------------------------------------
# Phase 2: Playwright fallback
# ---------------------------------------------------------------------------

FD_DRAFT_PAGE_CANDIDATES = [
    "https://sportsbook.fanduel.com/navigation/nfl-draft",
    "https://sportsbook.fanduel.com/navigation/nfl?tab=nfl-draft",
    "https://sportsbook.fanduel.com/navigation/nfl",
]


def run_browser_phase() -> tuple[dict | None, str | None, dict]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("  Playwright not installed. Install with:")
        print("    pip install playwright && playwright install chromium")
        return None, None, {}

    print("Phase 2: opening headed Chrome to capture draft markets...")
    profile = Path(__file__).resolve().parent.parent / ".cookies" / "fd_recon_profile"
    profile.mkdir(parents=True, exist_ok=True)

    captured: list[tuple[str, bytes]] = []

    def on_response(response):
        url = response.url
        if "fanduel.com" not in url:
            return
        ct = (response.headers.get("content-type") or "").lower()
        if "json" not in ct:
            return
        lowered = url.lower()
        # event-page, content-managed-page, and search are the three big
        # payload endpoints. Anything else is usually pings/tracking.
        if not any(hint in lowered for hint in (
            "event-page", "content-managed-page", "search", "draft",
        )):
            return
        try:
            body = response.body()
        except Exception:
            return
        if len(body) < 500:
            return
        captured.append((url, body))

    with sync_playwright() as p:
        ctx = p.chromium.launch_persistent_context(
            user_data_dir=str(profile),
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = ctx.pages[0] if ctx.pages else ctx.new_page()
        page.on("response", on_response)

        for url in FD_DRAFT_PAGE_CANDIDATES:
            print(f"  navigating to {url}")
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as exc:
                print(f"    nav failed: {exc!r}")
                continue
            page.wait_for_timeout(8000)

        print("\n  If NFL Draft markets aren't on screen, navigate there now.")
        print("  Press ENTER when markets are visible to save the capture...")
        try:
            input()
        except (EOFError, KeyboardInterrupt):
            print("  (no stdin - continuing with whatever was captured)")

        ctx.close()

    if not captured:
        print("  No FD market responses captured.")
        return None, None, {}

    captured.sort(key=lambda t: len(t[1]), reverse=True)
    best_url, best_body = captured[0]
    try:
        data = json.loads(best_body.decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"  biggest response wasn't JSON: {exc!r}")
        return None, None, {}
    print(f"  captured {len(captured)} responses; saving largest ({len(best_body):,} bytes)")
    return data, best_url, {"method": "playwright_capture", "url": best_url}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="FD NFL Draft markets recon")
parser.add_argument(
    "--browser",
    action="store_true",
    help="Skip the REST probe and go straight to the headed browser capture.",
)
args = parser.parse_args()

print("=" * 60)
print("  FANDUEL NFL DRAFT RECON")
print("=" * 60)
ensure_fixture_dirs()

data = None
url = None
meta: dict = {}

if not args.browser:
    data, url, meta = run_rest_phase()

if data is None:
    data, url, meta = run_browser_phase()

if data is None:
    print("\nERROR: Could not capture FD NFL Draft markets.")
    print("  Next steps:")
    print("  - Re-run with --browser and navigate to the draft page manually")
    print("  - If REST returns 400: the x-px-context token in FD_PX_CONTEXT is stale.")
    print("    Open sportsbook.fanduel.com in Chrome DevTools, copy a fresh")
    print("    x-px-context value from any event-page request, paste it here.")
    sys.exit(1)

path = save_fixture("fanduel", data, meta=meta)
print_diagnostics("fanduel", path, url, data)
print(f"  method used: {meta.get('method', 'unknown')}")
if "page_slug" in meta:
    print(f"  page slug: {meta['page_slug']}")
if "draft_tab_id" in meta:
    print(f"  draft tab_id: {meta['draft_tab_id']}  <- filter markets to this tab")
if "draft_market_count" in meta:
    print(f"  {meta['draft_market_count']} draft markets / "
          f"{meta.get('total_market_count', '?')} total on NFL page")
