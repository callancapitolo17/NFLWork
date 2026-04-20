#!/usr/bin/env python3
"""Bookmaker.eu NFL Draft markets recon.

Goal: capture the raw JSON that powers BM's NFL Draft futures (aka "NFL Props"
/ "NFL Odds To Win" / "NFL Specials") and save it to
nfl_draft/tests/fixtures/bookmaker/draft_markets.json.

This is HARDER than DK/FD because:
- BM's internal API is login-gated (Cloudflare + ASP.NET session cookies).
- Futures markets live under league IDs we don't know yet. The existing
  bookmaker_odds/scraper.py hardcodes league IDs for CBB/NBA/MLB game lines
  but has no entry for NFL Draft.

Strategy (two-phase):

Phase 1 - Probe known league ID ranges using saved cookies:
    BM's /gateway/BetslipProxy.aspx/GetSchedule accepts any league ID and
    returns a Schedule skeleton. If the ID is used, Schedule.Data.Leagues is
    non-empty and the League name tells us what the ID represents.

    The existing scraper proves IDs 3-6 are main-sport game lines. Futures /
    props typically live in the 4000-7000 range. We probe a curated list of
    candidates and any ID whose league name contains "draft", "specials",
    "props", or "to win" is a hit.

Phase 2 - Headed Playwright capture:
    If probes fail (session expired, or BM moved the IDs), open Chrome to BM's
    NFL Specials page, let the user log in if needed, and capture the
    GetSchedule response body that loads when the page shows draft markets.

Usage:
    python nfl_draft/scrapers/recon_bm.py                 # probe first, browser if needed
    python nfl_draft/scrapers/recon_bm.py --browser       # go straight to browser
    python nfl_draft/scrapers/recon_bm.py --probe 4029    # probe a specific league ID
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from nfl_draft.scrapers._recon_util import (
    _main_repo_root,
    ensure_fixture_dirs,
    fixture_dir,
    load_env,
    print_diagnostics,
    save_fixture,
)

from curl_cffi import requests as cffi_requests

# ---------------------------------------------------------------------------
# BM config - mirrors bookmaker_odds/scraper.py
# ---------------------------------------------------------------------------

BM_SITE_URL = "https://be.bookmaker.eu/en/sports/"
BM_API_BASE = "https://be.bookmaker.eu/gateway/BetslipProxy.aspx"
# Reuse the cookie file the production scraper already maintains so we don't
# have to log in twice. Resolve against the main repo root so it works from a
# worktree too (bookmaker_odds/ lives in the main working tree).
BM_COOKIE_PATH = _main_repo_root() / "bookmaker_odds" / ".bookmaker_cookies.json"


def _load_cookies(session: cffi_requests.Session) -> bool:
    """Load saved BM cookies from the bookmaker_odds dir. Returns True if found."""
    if not BM_COOKIE_PATH.exists():
        return False
    try:
        cookies = json.loads(BM_COOKIE_PATH.read_text())
        for name, value in cookies.items():
            session.cookies.set(name, value, domain=".bookmaker.eu")
        return bool(cookies)
    except Exception:
        return False


def _login(session: cffi_requests.Session, username: str, password: str) -> bool:
    """Same login flow as bookmaker_odds/scraper.py."""
    body = {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "Player": username,
                    "Password": password,
                    "loginKey": "",
                }
            }
        }
    }
    try:
        resp = session.post(f"{BM_API_BASE}/Login", json=body, timeout=15)
        return resp.status_code == 200
    except Exception:
        return False


def _fetch_schedule(session: cffi_requests.Session, league_id: str | int) -> dict | None:
    """Single GetSchedule call for one league ID. Returns parsed JSON or None."""
    body = {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "LeaguesIdList": str(league_id),
                    "LanguageId": "0",
                    "LineStyle": "E",
                    "ScheduleType": "american",
                    "LinkDeriv": "true",
                }
            }
        }
    }
    try:
        resp = session.post(f"{BM_API_BASE}/GetSchedule", json=body, timeout=15)
    except Exception:
        return None
    if resp.status_code != 200:
        return None
    try:
        return resp.json()
    except Exception:
        return None


def _looks_like_draft_league(league_entry: dict) -> bool:
    """Heuristic: does this league name/description suggest active NFL Draft futures?

    We look at BM's `Description` (or `desc`) field. To avoid matching old years,
    we reject any name with 2019-2025. Active draft = 2026 or unversioned.
    """
    name = (league_entry.get("Description") or league_entry.get("desc")
            or league_entry.get("name") or "").upper()
    if not name:
        return False
    # Reject stale years.
    if any(str(y) in name for y in range(2019, 2026)):
        return False
    # Must mention both NFL and DRAFT to count as an NFL Draft market.
    if "NFL" in name and "DRAFT" in name:
        return True
    return False


def _count_games_and_extract_leagues(data: dict) -> list[dict]:
    """Pull the League[] array out of a GetSchedule response."""
    if not isinstance(data, dict):
        return []
    return (
        data.get("Schedule", {})
        .get("Data", {})
        .get("Leagues", {})
        .get("League", [])
        or []
    )


# League ID candidates for BM. Live probing (2026-04-17) showed BM uses small
# dense IDs: 1-18 are main sports (NFL=1, CBB=4, MLB=5, NHL=7, etc). The 100-
# 1500 range was fully parallel-scanned for this account and found only:
#   104/105 = NFL 1ST/2ND HALVES, 117/208/211 = TNT "ODDS TO WIN" (all empty)
# ZERO leagues matched "NFL Draft" — BM likely gates draft futures behind a
# funded account or has them offline until closer to draft day.
#
# The REST probe here exists mostly for quick sanity. The authoritative path
# is the browser fallback, which captures the GetSchedule POST body's
# LeaguesIdList value from the live page.
def _candidate_league_ids() -> list[int]:
    return [
        # Main + known extras (shape probes)
        1, 11, 18, 104, 105, 202, 205, 503, 13554,
        # Observed TNT 'ODDS TO WIN' futures buckets — NFL Draft might land here
        # when it goes live.
        117, 208, 211,
    ]


def run_rest_phase(probe_only: int | None = None) -> tuple[dict | None, str | None, dict]:
    """Probe BM league IDs via the authenticated GetSchedule endpoint.

    If probe_only is set, hit just that one ID and return whatever it returns
    (useful for debugging when the user already knows the ID).
    """
    print("Phase 1: REST probe of BM NFL Draft league IDs...")
    load_env()
    username = os.getenv("BOOKMAKER_USERNAME")
    password = os.getenv("BOOKMAKER_PASSWORD")

    session = cffi_requests.Session(impersonate="chrome")
    loaded = _load_cookies(session)
    if loaded:
        print(f"  loaded cookies from {BM_COOKIE_PATH.name}")
    else:
        print("  no saved BM cookies; will need to refresh via recon_bookmaker.py")

    # Hit the site to refresh Cloudflare tokens. A 403 means cookies are stale.
    try:
        site_resp = session.get(BM_SITE_URL, timeout=15)
        if site_resp.status_code == 403:
            print("  BM returned 403 — Cloudflare cookies expired.")
            print("  Run: python bookmaker_odds/recon_bookmaker.py to refresh, then retry.")
            return None, None, {}
    except Exception as exc:
        print(f"  cloudflare warmup failed: {exc!r}")
        return None, None, {}

    # Re-authenticate if we have creds (cheap; avoids stale-session surprises).
    if username and password:
        _login(session, username, password)

    # Probe candidates.
    candidates = [probe_only] if probe_only else _candidate_league_ids()

    draft_hits: list[tuple[int, dict, dict]] = []  # (league_id, response, league_entry)

    for lid in candidates:
        data = _fetch_schedule(session, lid)
        leagues = _count_games_and_extract_leagues(data)
        if not leagues:
            continue
        # Any league entry in this response that looks like NFL Draft = hit.
        for le in leagues:
            if _looks_like_draft_league(le):
                name = le.get("desc") or le.get("Description") or le.get("name")
                print(f"  HIT lg={lid}: {name!r}")
                draft_hits.append((lid, data, le))
                break
        else:
            # Not a draft match, but worth logging the first few we find so the
            # user can eyeball misses. Only log when probing a specific ID.
            if probe_only and leagues:
                name = leagues[0].get("desc") or leagues[0].get("Description")
                print(f"  lg={lid} returned league: {name!r} (not draft)")

    if not draft_hits:
        return None, None, {}

    # Pick the hit with the biggest raw JSON — more games = more likely to be
    # the main Draft markets vs a one-off special.
    draft_hits.sort(key=lambda t: len(json.dumps(t[1])), reverse=True)
    winner_id, winner_data, winner_entry = draft_hits[0]
    # Build meta with ALL discovered IDs so the real scraper can use them all.
    all_ids = [h[0] for h in draft_hits]
    url = f"{BM_API_BASE}/GetSchedule  (lg={winner_id})"
    return winner_data, url, {
        "league_id": winner_id,
        "league_name": winner_entry.get("desc") or winner_entry.get("Description"),
        "all_draft_league_ids": all_ids,
        "method": "rest_probe",
    }


# ---------------------------------------------------------------------------
# Phase 2: Playwright capture
# ---------------------------------------------------------------------------

BM_DRAFT_PAGE_CANDIDATES = [
    # BM's sport tree is filter-driven by JS; we land on the sports root and
    # ask the user to navigate to NFL Draft / NFL Futures. The user click
    # triggers the GetSchedule POST that we capture.
    "https://be.bookmaker.eu/en/sports/",
    "https://be.bookmaker.eu/en/sports/football/",
]


def run_browser_phase() -> tuple[dict | None, str | None, dict]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("  Playwright not installed.")
        return None, None, {}

    print("Phase 2: opening headed Chrome to capture BM draft markets...")
    # Reuse the profile from the main BM recon so cookies/Cloudflare session carry over.
    # Resolve against main repo root (not a worktree-relative path).
    profile = _main_repo_root() / "bookmaker_odds" / ".bookmaker_profile"
    profile.mkdir(parents=True, exist_ok=True)

    captured: list[tuple[str, bytes, str | None]] = []  # (url, body, post_data)

    def on_response(response):
        url = response.url
        if "bookmaker.eu" not in url:
            return
        # We specifically want GetSchedule responses — skip GetConfig, Login,
        # GetBetslip, GetPlayerInfo, and all the other RPCs on BetslipProxy.
        if "GetSchedule" not in url:
            return
        try:
            body = response.body()
        except Exception:
            return
        # GetSchedule with no data returns ~80 bytes (empty League list).
        # A real draft response is 10KB+. Filter aggressively so the biggest
        # capture is always a real payload.
        if len(body) < 2000:
            return
        try:
            post = response.request.post_data
        except Exception:
            post = None
        # Verify the body actually has a non-empty League array — otherwise
        # it's a GetSchedule for an ID that returned nothing.
        try:
            j = json.loads(body.decode("utf-8", errors="replace"))
        except Exception:
            return
        leagues = (j.get("Schedule", {}).get("Data", {})
                    .get("Leagues", {}).get("League", []) or [])
        if isinstance(leagues, dict):
            leagues = [leagues]
        if not leagues:
            return
        captured.append((url, body, post))

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

        for url in BM_DRAFT_PAGE_CANDIDATES:
            print(f"  navigating to {url}")
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as exc:
                print(f"    nav failed: {exc!r}")
                continue
            page.wait_for_timeout(5000)

        print("\n  Log in if needed, then navigate to NFL Draft / NFL Props / NFL Futures.")
        print("  We listen for GetSchedule responses with real league data >2KB.")
        print("  Press ENTER when the draft markets are visible to save the capture...")
        try:
            input()
        except (EOFError, KeyboardInterrupt):
            pass

        ctx.close()

    if not captured:
        print("  No GetSchedule responses with league data captured.")
        print("  - Did you reach the NFL Draft markets page?")
        print("  - Did the page actually show odds (vs an empty card)?")
        print("  - GetConfig / Login / empty GetSchedule responses are filtered out;")
        print("    if you ONLY saw those, BM didn't serve the draft markets.")
        return None, None, {}

    # Pick the biggest GetSchedule body. Try to extract the league ID from the
    # POST body (LeaguesIdList).
    captured.sort(key=lambda t: len(t[1]), reverse=True)
    best_url, best_body, best_post = captured[0]
    try:
        data = json.loads(best_body.decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"  biggest response wasn't JSON: {exc!r}")
        return None, None, {}

    meta: dict = {"method": "playwright_capture", "url": best_url}
    if best_post:
        # Parse the LeaguesIdList value out of the POST body.
        try:
            post_json = json.loads(best_post)
            lid = (post_json.get("o", {}).get("BORequestData", {})
                   .get("BOParameters", {}).get("LeaguesIdList"))
            if lid:
                meta["league_id"] = lid
                print(f"  discovered league_id from POST body: {lid}")
        except Exception:
            pass
    leagues = _count_games_and_extract_leagues(data)
    if leagues:
        meta["league_name"] = (leagues[0].get("desc")
                               or leagues[0].get("Description"))
    print(f"  captured {len(captured)} responses; saving largest ({len(best_body):,} bytes)")
    return data, best_url, meta


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BM NFL Draft markets recon")
    parser.add_argument("--browser", action="store_true",
                        help="Skip probe; go straight to browser capture.")
    parser.add_argument("--probe", type=int, default=None,
                        help="Probe a single specific league ID (for debugging).")
    args = parser.parse_args()

    print("=" * 60)
    print("  BOOKMAKER NFL DRAFT RECON")
    print("=" * 60)
    ensure_fixture_dirs()

    data = None
    url = None
    meta: dict = {}

    if not args.browser:
        data, url, meta = run_rest_phase(probe_only=args.probe)

    if data is None and not args.probe:
        data, url, meta = run_browser_phase()

    if data is None:
        print("\nERROR: Could not capture BM NFL Draft markets.")
        print("  Next steps:")
        print("  - Run: python bookmaker_odds/recon_bookmaker.py  (to refresh Cloudflare cookies)")
        print("  - Then re-run this script")
        print("  - If still failing, run with --browser and manually navigate to NFL Draft")
        sys.exit(1)

    path = save_fixture("bookmaker", data, meta=meta)
    print_diagnostics("bookmaker", path, url, data)
    print(f"  method used: {meta.get('method', 'unknown')}")
    if "league_id" in meta:
        print(f"  league_id: {meta['league_id']}  <- hardcode in the real scraper")
    if "league_name" in meta:
        print(f"  league_name: {meta['league_name']}")
    if "all_draft_league_ids" in meta and len(meta["all_draft_league_ids"]) > 1:
        print(f"  ALL matching IDs: {meta['all_draft_league_ids']}")
