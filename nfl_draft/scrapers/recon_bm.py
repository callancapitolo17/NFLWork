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


def _save_cookies(session: cffi_requests.Session) -> bool:
    """Persist session cookies to disk for future runs.

    Mirrors bookmaker_odds/scraper.py._save_cookies so both paths share the
    same cookie jar file. Returns True on success, False on any IO error.
    """
    try:
        cookies = dict(session.cookies)
        BM_COOKIE_PATH.parent.mkdir(parents=True, exist_ok=True)
        BM_COOKIE_PATH.write_text(json.dumps(cookies))
        return True
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


def _auto_relogin(session: cffi_requests.Session) -> bool:
    """Attempt a pure-HTTP re-login using saved credentials.

    Returns True iff creds were present AND _login returned 200. The caller
    decides whether to retry the probe. Never prints credential VALUES — only
    the outcome. Mirrors the refresh-on-auth-fail pattern from
    bet_logger/scraper_betonline.py.refresh_access_token: one retry, no loop.
    """
    load_env()
    username = os.getenv("BOOKMAKER_USERNAME")
    password = os.getenv("BOOKMAKER_PASSWORD")
    if not username or not password:
        print("  auto-re-login: BOOKMAKER_USERNAME / BOOKMAKER_PASSWORD missing from bet_logger/.env")
        return False
    ok = _login(session, username, password)
    if ok:
        print("  auto-re-login: login succeeded")
    else:
        print("  auto-re-login: login failed (non-200 response)")
    return ok


def _looks_like_expired_cookies(probe_results: list[dict | None]) -> bool:
    """Heuristic: does this set of probe responses indicate expired cookies?

    BM returns a few telltale shapes when the ASP.NET session cookie is stale:
      - HTTP 401/403 at the transport layer -> _fetch_schedule returns None
      - A guest / 'IsDummyPlayer: true' marker embedded anywhere in the response
      - Every GetSchedule returns an empty Schedule.Data.Leagues.League[]

    We accept a list of probe responses (one per league ID we tried). If EVERY
    probe came back None or empty, we flag it as expired. A single non-empty
    response means the session is fine and we shouldn't risk a login attempt.
    """
    if not probe_results:
        return False

    saw_any_leagues = False
    for data in probe_results:
        if data is None:
            # Transport-level failure (403/non-200/exception). Could be expiry;
            # keep looking for evidence of a real response.
            continue
        if not isinstance(data, dict):
            continue
        # Explicit guest/dummy marker: if BM hints we're not logged in, that's
        # the strongest signal. Walk the common locations.
        if _has_dummy_player_marker(data):
            return True
        leagues = _count_games_and_extract_leagues(data)
        if leagues:
            saw_any_leagues = True
            break

    # If NOTHING came back with leagues, treat as expired so caller can try a
    # fresh login. Harmless if the account really does have no draft markets:
    # we attempt ONE login, retry ONCE, and bail.
    return not saw_any_leagues


def _has_dummy_player_marker(data: dict) -> bool:
    """Scan a BM response for the `IsDummyPlayer: true` guest marker.

    BM's GetConfig / GetSchedule sometimes nest player status at variable depths.
    A shallow recursive scan is cheap and keeps this tolerant to shape drift.
    """
    stack: list = [data]
    depth = 0
    while stack and depth < 6:  # cap recursion depth
        current = stack.pop()
        if isinstance(current, dict):
            # Key can appear as 'IsDummyPlayer' or 'isDummyPlayer' depending on endpoint.
            for k, v in current.items():
                if isinstance(k, str) and k.lower() == "isdummyplayer":
                    if v is True or (isinstance(v, str) and v.lower() == "true"):
                        return True
                if isinstance(v, (dict, list)):
                    stack.append(v)
        elif isinstance(current, list):
            stack.extend(current)
        depth += 1
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


def _probe_candidates(session: cffi_requests.Session,
                      candidates: list[int]) -> tuple[list[tuple[int, dict, dict]],
                                                      list[dict | None]]:
    """Probe a list of league IDs; return (draft_hits, raw_responses).

    `draft_hits` is the subset that matched the NFL-Draft league heuristic.
    `raw_responses` is every probe's raw JSON (or None if the call failed) so
    a caller can feed it into _looks_like_expired_cookies without re-probing.
    """
    draft_hits: list[tuple[int, dict, dict]] = []
    raw_responses: list[dict | None] = []
    for lid in candidates:
        data = _fetch_schedule(session, lid)
        raw_responses.append(data)
        leagues = _count_games_and_extract_leagues(data)
        if not leagues:
            continue
        for le in leagues:
            if _looks_like_draft_league(le):
                name = le.get("desc") or le.get("Description") or le.get("name")
                print(f"  HIT lg={lid}: {name!r}")
                draft_hits.append((lid, data, le))
                break
        else:
            # When probing a single ID, echo the non-match so the user can see
            # what BM returned. Skip this noise for the bulk scan.
            if len(candidates) == 1 and leagues:
                name = leagues[0].get("desc") or leagues[0].get("Description")
                print(f"  lg={lid} returned league: {name!r} (not draft)")
    return draft_hits, raw_responses


def _pick_winner(draft_hits: list[tuple[int, dict, dict]]
                 ) -> tuple[dict, str, dict] | tuple[None, None, dict]:
    """Choose the biggest hit and build the return tuple for run_rest_phase."""
    if not draft_hits:
        return None, None, {}
    # Biggest raw JSON wins — more markets = more likely the main draft board
    # rather than an incidental one-off special.
    draft_hits.sort(key=lambda t: len(json.dumps(t[1])), reverse=True)
    winner_id, winner_data, winner_entry = draft_hits[0]
    all_ids = [h[0] for h in draft_hits]
    url = f"{BM_API_BASE}/GetSchedule  (lg={winner_id})"
    return winner_data, url, {
        "league_id": winner_id,
        "league_name": winner_entry.get("desc") or winner_entry.get("Description"),
        "all_draft_league_ids": all_ids,
        "method": "rest_probe",
    }


def run_rest_phase(probe_only: int | None = None) -> tuple[dict | None, str | None, dict]:
    """Probe BM league IDs via the authenticated GetSchedule endpoint.

    Uses the BetOnline-style refresh pattern: if the first round of probes
    looks like expired cookies (empty everywhere / guest response), try ONE
    auto-re-login via saved credentials, save fresh cookies, and retry the
    probes exactly once. Never loops.

    If probe_only is set, hit just that one ID (useful for debugging).
    """
    print("Phase 1: REST probe of BM NFL Draft league IDs...")

    session = cffi_requests.Session(impersonate="chrome")
    loaded = _load_cookies(session)
    if loaded:
        print(f"  loaded cookies from {BM_COOKIE_PATH.name}")
    else:
        print("  no saved BM cookies; will attempt auto-re-login")

    # Hit the site to refresh Cloudflare tokens. A 403 means cf_clearance is
    # gone — pure HTTP can't mint it, so we bail to the Playwright fallback.
    try:
        site_resp = session.get(BM_SITE_URL, timeout=15)
        if site_resp.status_code == 403:
            print("  BM returned 403 — Cloudflare cookies expired.")
            print("  Run: python bookmaker_odds/recon_bookmaker.py to refresh, then retry.")
            return None, None, {}
    except Exception as exc:
        print(f"  cloudflare warmup failed: {exc!r}")
        return None, None, {}

    candidates = [probe_only] if probe_only else _candidate_league_ids()

    # First probe pass — uses whatever cookies we have on disk.
    draft_hits, raw_responses = _probe_candidates(session, candidates)

    # If we got real hits, ship them. Don't touch login; don't rotate cookies.
    if draft_hits:
        return _pick_winner(draft_hits)

    # No hits. Heuristic: does it look like the session is expired, or did BM
    # simply not post draft markets today? If expired-looking, try ONE login
    # and retry ONCE. This mirrors bet_logger/scraper_betonline.py's
    # refresh_access_token -> authenticated_call retry shape.
    if not _looks_like_expired_cookies(raw_responses):
        print("  probes returned no draft hits, but responses don't look expired.")
        print("  BM may not be posting NFL Draft markets for this account right now.")
        return None, None, {}

    print("  cookies appear stale; attempting auto-re-login...")
    if not _auto_relogin(session):
        print("  falling back: run `python bookmaker_odds/recon_bookmaker.py` "
              "to refresh cookies manually.")
        return None, None, {}

    # Persist fresh cookies so the next recon / scraper run picks them up.
    if _save_cookies(session):
        print(f"  saved refreshed cookies to {BM_COOKIE_PATH.name}")

    # Retry probes exactly once. If still empty, it's genuinely no-draft.
    draft_hits, raw_responses = _probe_candidates(session, candidates)
    if draft_hits:
        return _pick_winner(draft_hits)

    if _looks_like_expired_cookies(raw_responses):
        print("  still no data after re-login — Cloudflare may also need refreshing.")
        print("  Run: python bookmaker_odds/recon_bookmaker.py to refresh manually.")
    else:
        print("  auto-re-login succeeded but BM isn't posting NFL Draft markets "
              "for this account right now. Will auto-refresh when they appear.")
    return None, None, {}


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
