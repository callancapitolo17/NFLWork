#!/usr/bin/env python3
"""DraftKings NFL Draft markets recon.

Goal: capture the raw JSON that powers DK's NFL Draft futures page and save it
to nfl_draft/tests/fixtures/draftkings/draft_markets.json so the parser subagent
can work offline.

Key finding from live recon (2026-04-17):
    DK stores NFL Draft markets under NFL (leagueId=88808), categoryId=1803
    ('2026 Draft'). The category has ~20 subcategories (Pick Number, Player
    Drafted By, 1st Selected By Position, Draft Position, Top 5 Pick,
    Top 10 Pick, Round 1 Pick, Mr. Irrelevant, etc).

    The public endpoint is:
        https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent
        /dkusnj/v1/leagues/88808                              -> full index
        /leagues/88808/categories/1803                        -> default sub
        /leagues/88808/categories/1803/subcategories/<sub_id> -> that sub's markets

    No auth or cookies needed; curl_cffi impersonate='chrome' bypasses Akamai.

Strategy (two-phase):

Phase 1 - Pure REST (fast, no browser):
    1. Fetch /leagues/88808 to discover all subcategories under category 1803.
    2. Fetch each 1803 subcategory and merge into one fixture envelope.
    3. Save. Done.

Phase 2 - Headed Playwright fallback:
    If Akamai rate-limits the REST path, open Chrome to the NFL Draft URL and
    capture any sportscontent JSON response that fires during page load.

Usage:
    python nfl_draft/scrapers/recon_dk.py              # REST first
    python nfl_draft/scrapers/recon_dk.py --browser    # force browser fallback
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make the parent package importable when running as a plain script.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from nfl_draft.scrapers._recon_util import (
    ensure_fixture_dirs,
    print_diagnostics,
    save_fixture,
    fixture_path,
)

from curl_cffi import requests as cffi_requests

# ---------------------------------------------------------------------------
# DK public REST config
# ---------------------------------------------------------------------------

# Same REST pattern as mlb_sgp/scraper_draftkings_sgp.py. No auth, no cookies.
# Akamai bot detection is defeated by curl_cffi's Chrome TLS fingerprint.
DK_BASE = "https://sportsbook.draftkings.com"
DK_NASH = "https://sportsbook-nash.draftkings.com"

# leagueId=88808 is NFL. categoryId=1803 is "2026 Draft". Confirmed via live
# probe 2026-04-17. If DK bumps the season, run this script again to re-discover
# (the probe output will show a category named e.g. "2027 Draft").
DK_NFL_LEAGUE_ID = 88808
DK_NFL_DRAFT_CATEGORY_ID = 1803

# Headers that match a fresh Chrome session. curl_cffi sets TLS + JA3 to match;
# these application headers mirror what the site sends.
DK_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/146.0.0.0 Safari/537.36"
    ),
    "Accept": "application/json, text/plain, */*",
    "Referer": "https://sportsbook.draftkings.com/leagues/football/nfl",
    "Origin": "https://sportsbook.draftkings.com",
}


def _league_url(league_id: int) -> str:
    return f"{DK_NASH}/sites/US-SB/api/sportscontent/dkusnj/v1/leagues/{league_id}"


def _subcategory_url(league_id: int, category_id: int, subcategory_id: int) -> str:
    return (
        f"{DK_NASH}/sites/US-SB/api/sportscontent/dkusnj/v1/leagues/"
        f"{league_id}/categories/{category_id}/subcategories/{subcategory_id}"
    )


def _fetch_json(session: cffi_requests.Session, url: str) -> dict | None:
    """GET a DK API URL and return parsed JSON, or None on any failure.

    Wraps exceptions and status checks so the caller can just test truthiness.
    """
    try:
        resp = session.get(url, headers=DK_HEADERS, timeout=20)
    except Exception as exc:
        print(f"  GET {url[-60:]}: request error {exc!r}")
        return None
    if resp.status_code != 200:
        print(f"  GET {url[-60:]}: HTTP {resp.status_code}")
        return None
    try:
        return resp.json()
    except Exception:
        print(f"  GET {url[-60:]}: non-JSON body ({len(resp.text)} chars)")
        return None


def run_rest_phase() -> tuple[dict | None, str | None, dict]:
    """REST path: fetch the NFL league index, then every 1803 subcategory.

    Returns (envelope, source_url, meta) where envelope is a composite dict:
        {
            "league_index": <full NFL league response>,
            "draft_subcategories": {
                "17747": <Pick Number response>,
                "19620": <Draft Position response>,
                ...
            },
            "draft_category_name": "2026 Draft",
        }

    The parser subagent can then consume either the league_index (for the
    subcategory catalog) or any specific draft_subcategories entry for that
    sub's markets. All responses have the same shape: `markets`, `selections`,
    `events`, `subcategories`, etc.
    """
    print("Phase 1: REST fetch of DK NFL Draft markets...")
    session = cffi_requests.Session(impersonate="chrome")
    # Warm Akamai: hit the NFL page so the session has fresh bot-check cookies.
    # Without this, subsequent API calls sometimes 403.
    try:
        session.get(f"{DK_BASE}/leagues/football/nfl", headers=DK_HEADERS, timeout=20)
    except Exception:
        pass

    # Step 1: fetch the NFL league index to get the subcategory list.
    league_url = _league_url(DK_NFL_LEAGUE_ID)
    league_data = _fetch_json(session, league_url)
    if not league_data:
        print(f"  could not fetch league index for {DK_NFL_LEAGUE_ID}")
        return None, None, {}

    # Find category 1803 and its subcategories.
    categories = {c.get("id"): c for c in league_data.get("categories", [])}
    draft_cat = categories.get(DK_NFL_DRAFT_CATEGORY_ID)
    if not draft_cat:
        # Category ID changed - dump what we see so the user can eyeball.
        print(f"  category {DK_NFL_DRAFT_CATEGORY_ID} not in NFL index. Found:")
        for cid, c in categories.items():
            print(f"    {cid}: {c.get('name')!r}")
        return None, None, {}

    draft_cat_name = draft_cat.get("name", "?")
    print(f"  found category {DK_NFL_DRAFT_CATEGORY_ID}: {draft_cat_name!r}")

    # Subcategories are in the league index; filter by categoryId.
    draft_subs = [
        s for s in league_data.get("subcategories", [])
        if s.get("categoryId") == DK_NFL_DRAFT_CATEGORY_ID
    ]
    print(f"  {len(draft_subs)} subcategories under '{draft_cat_name}'")

    # Step 2: fetch each draft subcategory.
    draft_subcats: dict[str, dict] = {}
    total_markets = 0
    for sub in draft_subs:
        sid = sub.get("id")
        sname = sub.get("name")
        sub_url = _subcategory_url(DK_NFL_LEAGUE_ID, DK_NFL_DRAFT_CATEGORY_ID, sid)
        sub_data = _fetch_json(session, sub_url)
        if not sub_data:
            print(f"    sub {sid} ({sname!r}): skipped")
            continue
        m_count = len(sub_data.get("markets", []))
        s_count = len(sub_data.get("selections", []))
        draft_subcats[str(sid)] = sub_data
        total_markets += m_count
        print(f"    sub {sid} {sname!r}: {m_count} markets, {s_count} selections")

    if not draft_subcats:
        print("  no draft subcategories returned data - DK session may be blocked")
        return None, None, {}

    envelope = {
        "league_id": DK_NFL_LEAGUE_ID,
        "draft_category_id": DK_NFL_DRAFT_CATEGORY_ID,
        "draft_category_name": draft_cat_name,
        "league_index": league_data,
        "draft_subcategories": draft_subcats,
        "subcategory_catalog": [
            {"id": s.get("id"), "name": s.get("name"),
             "categoryId": s.get("categoryId")}
            for s in draft_subs
        ],
    }
    meta = {
        "league_id": DK_NFL_LEAGUE_ID,
        "draft_category_id": DK_NFL_DRAFT_CATEGORY_ID,
        "draft_category_name": draft_cat_name,
        "subcategory_count": len(draft_subcats),
        "total_markets": total_markets,
        "method": "rest_controldata",
    }
    return envelope, league_url, meta


# ---------------------------------------------------------------------------
# Phase 2: Playwright fallback
# ---------------------------------------------------------------------------

DK_DRAFT_PAGE_CANDIDATES = [
    "https://sportsbook.draftkings.com/leagues/football/nfl-draft",
    "https://sportsbook.draftkings.com/event/nfl-draft",
    "https://sportsbook.draftkings.com/leagues/football/nfl",
]


def run_browser_phase() -> tuple[dict | None, str | None, dict]:
    """Open headed Chrome and capture the biggest sportscontent JSON response.

    The user doesn't need to click anything — the page fires enough XHRs on
    load to expose the draft markets endpoint. We listen for responses whose
    URL looks like a DK markets API, collect them, and at the end pick the
    largest one. That heuristic worked for MLB SGP recon.
    """
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("  Playwright not installed in this venv.")
        print("  Install with: pip install playwright && playwright install chromium")
        return None, None, {}

    print("Phase 2: opening headed Chrome to capture draft markets...")
    # Use a profile dir under nfl_draft/.cookies/ so it's already gitignored.
    profile = Path(__file__).resolve().parent.parent / ".cookies" / "dk_recon_profile"
    profile.mkdir(parents=True, exist_ok=True)

    captured: list[tuple[str, bytes]] = []  # (url, body_bytes)

    def on_response(response):
        url = response.url
        # Only care about market/offer JSON responses on DK hosts.
        if "draftkings.com" not in url:
            return
        # Filter by content-type; DK returns HTML on the main page, JSON on APIs.
        ct = (response.headers.get("content-type") or "").lower()
        if "json" not in ct:
            return
        lowered = url.lower()
        if not any(hint in lowered for hint in (
            "sportscontent", "eventgroup", "offercategory", "markets",
            "categories", "subcategories",
        )):
            return
        try:
            # response.body() can fail on live streams; wrap defensively.
            body = response.body()
        except Exception:
            return
        if len(body) < 500:
            return  # skip pings, pixels, tiny responses
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

        for url in DK_DRAFT_PAGE_CANDIDATES:
            print(f"  navigating to {url}")
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as exc:
                print(f"    nav failed: {exc!r}")
                continue
            # Let XHRs fire. 8s is plenty for a futures page.
            page.wait_for_timeout(8000)

        # If the user needs to click something (e.g. NFL Draft sub-tab), give
        # them a chance before we close. This is the only manual step.
        print("\n  If the NFL Draft markets aren't on screen, navigate there now.")
        print("  Press ENTER when the markets are visible to save the capture...")
        try:
            input()
        except (EOFError, KeyboardInterrupt):
            print("  (no stdin — continuing with whatever was captured)")

        ctx.close()

    if not captured:
        print("  No DK markets responses captured. Session may have been blocked.")
        return None, None, {}

    # Pick the largest body — that's almost always the full markets payload.
    import json
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
# Main (no __main__ guard per spec — direct script invocation)
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DK NFL Draft markets recon")
    parser.add_argument(
        "--browser",
        action="store_true",
        help="Skip the REST probe and go straight to the headed browser capture.",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  DRAFTKINGS NFL DRAFT RECON")
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
        print("\nERROR: Could not capture DK NFL Draft markets.")
        print("  Next steps:")
        print("  - Re-run with --browser and navigate to the NFL Draft page manually")
        print("  - If the REST probe failed, Akamai may have rate-limited; wait 5 min")
        print("  - If browser capture is empty, DK may have removed the NFL Draft page")
        sys.exit(1)

    path = save_fixture("draftkings", data, meta=meta)
    print_diagnostics("draftkings", path, url, data)
    print(f"  method used: {meta.get('method', 'unknown')}")
    if "league_id" in meta:
        print(f"  NFL league_id: {meta['league_id']}")
    if "draft_category_id" in meta:
        print(f"  draft category_id: {meta['draft_category_id']} "
              f"({meta.get('draft_category_name')!r})")
    if "subcategory_count" in meta:
        print(f"  captured {meta['subcategory_count']} subcategories, "
              f"{meta.get('total_markets', '?')} markets total")
