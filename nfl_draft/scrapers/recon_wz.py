#!/usr/bin/env python3
"""Wagerzon NFL Draft markets recon.

Goal: capture the raw JSON that powers WZ's NFL Draft futures page and save it
to nfl_draft/tests/fixtures/wagerzon/draft_markets.json.

WZ is an ASP.NET site. Authentication: GET the login page to capture
__VIEWSTATE + friends, POST credentials, then the ASP.NET_SessionId cookie
carries auth. Data comes from NewScheduleHelper.aspx?WT=0&lg=<comma-separated-ids>
as JSON.

The NFL config in wagerzon_odds/config.py uses lg=4029,430,432,... for game
lines. We don't know which ID (if any of those) carries the NFL Draft markets.

Strategy (two-phase):

Phase 1 - Probe lg values:
    Authenticate via the standard scraper_v2.py pattern, then hit
    NewScheduleHelper.aspx with each candidate lg ID. Any response with a
    Competition name like "NFL Draft" / "NFL Futures" / "NFL Odds to Win"
    is a hit.

Phase 2 - Headed Playwright fallback:
    Open WZ in Chrome, let user log in + navigate to NFL Draft page, and
    capture the NewScheduleHelper.aspx request that fires. The lg query-string
    parameter tells us the answer.

Usage:
    python nfl_draft/scrapers/recon_wz.py                 # probe first
    python nfl_draft/scrapers/recon_wz.py --browser       # go to browser
    python nfl_draft/scrapers/recon_wz.py --probe 4029    # single-ID probe
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from nfl_draft.scrapers._recon_util import (
    _main_repo_root,
    ensure_fixture_dirs,
    load_env,
    print_diagnostics,
    save_fixture,
)

import requests

# ---------------------------------------------------------------------------
# WZ config - mirrors wagerzon_odds/scraper_v2.py + config.py
# ---------------------------------------------------------------------------

WZ_BASE_URL = "https://backend.wagerzon.com"
WZ_HELPER_URL = f"{WZ_BASE_URL}/wager/NewScheduleHelper.aspx"


def _wz_login(session: requests.Session, username: str, password: str) -> bool:
    """Login via the ASP.NET form POST pattern from wagerzon_odds/scraper_v2.py."""
    resp = session.get(WZ_BASE_URL, timeout=15)
    if resp.status_code != 200:
        return False
    # If already redirected to schedule, we're already authenticated.
    if "NewSchedule" in resp.url:
        return True

    html = resp.text
    fields: dict[str, str] = {}
    for name in ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"):
        m = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if m:
            fields[name] = m.group(1)
    if "__VIEWSTATE" not in fields:
        return False

    fields["Account"] = username
    fields["Password"] = password
    fields["BtnSubmit"] = ""
    try:
        resp = session.post(WZ_BASE_URL, data=fields, timeout=15)
    except Exception:
        return False
    return resp.status_code == 200


def _wz_fetch(session: requests.Session, lg: str | int) -> dict | None:
    """Single NewScheduleHelper.aspx fetch for one (or many, comma-sep) lg IDs."""
    url = f"{WZ_HELPER_URL}?WT=0&lg={lg}"
    try:
        resp = session.get(url, timeout=20, headers={
            "Accept": "application/json, text/plain, */*",
            "X-Requested-With": "XMLHttpRequest",
        })
    except Exception:
        return None
    if resp.status_code != 200:
        return None
    try:
        return resp.json()
    except Exception:
        return None


def _find_competitions(data: dict) -> list[dict]:
    """Extract league entries from a NewScheduleHelper response.

    WZ's JSON shape (confirmed 2026-04-17):
        { "result": { "listLeagues": [ [<league>, ...], [<league>, ...], ... ] } }

    Each entry has: Description, IdLeague, IdSport, Games, ... — what WZ calls
    a "league" we treat as a competition bucket for filtering. We flatten the
    nested list so callers get one list of dicts.
    """
    if not isinstance(data, dict):
        return []
    result = data.get("result") or data.get("NewSchedule", {}).get("Data", {})
    if not isinstance(result, dict):
        return []
    buckets = result.get("listLeagues") or []
    flat: list[dict] = []
    for bucket in buckets:
        if isinstance(bucket, list):
            for lg in bucket:
                if isinstance(lg, dict):
                    flat.append(lg)
        elif isinstance(bucket, dict):
            flat.append(bucket)
    return flat


def _looks_like_draft_competition(comp: dict) -> bool:
    """Heuristic match for ACTIVE (2026) NFL Draft in a league entry.

    We reject old years (2020-2025) to avoid capturing stale leagues that WZ
    leaves on the books. The active season is detected by 'NFL DRAFT 2026' or
    a generic 'NFL DRAFT' with no year (WZ also uses unversioned labels like
    'NFL DRAFT TOTALS' / 'NFL DRAFT DRAFT SPECIALS').
    """
    name = (comp.get("Description") or comp.get("desc")
            or comp.get("Name") or comp.get("name") or "").upper()
    if not name:
        return False
    # Any old-year tag? Reject.
    if any(str(y) in name for y in range(2019, 2026)):
        return False
    # Must mention NFL + DRAFT somewhere.
    if "NFL" not in name or "DRAFT" not in name:
        return False
    return True


def _candidate_lg_ids() -> list[int]:
    """Confirmed (live 2026-04-17) lg IDs that carry active NFL Draft markets.

    Discovered via wide sweep against backend.wagerzon.com. These all returned
    Descriptions matching 'NFL DRAFT 2026 - ...' or unversioned 'NFL DRAFT ...'.
    If WZ rotates IDs between draft years, rerun this script's --browser mode
    to capture the actual lg= from the URL, then re-seed this list.

    Active overall-pick leagues:
        1270 (#1), 2867 (#2), 1370 (#3), 1867 (#4), 2611 (#5),
        3699-3703 (#6-10), 2614 (#11-15), 2618 (Mr. Irrelevant)
    Position/player props:
        2579 (Player Selected By Position), 478 (Players Draft Position),
        590 (1st Round Props), 977 (Other Props), 1231 (Totals), 3237
    Team-to-draft-player props:
        2516, 2560, 2561, 2562, 2563, 2564, 2565, 2580, 2581, 4537
    Prop bundles:
        4242 (Top 5), 4243 (Top 10), 4245 (1st Round), 4246 (Draft Specials)
    Live:
        1269 (Live NFL Draft)
    """
    return [
        # Overall pick props (2026)
        1270, 2867, 1370, 1867, 2611, 3699, 3700, 3701, 3702, 3703,
        2614, 2618,
        # Position / player selection props
        2579, 478, 590, 977, 1231, 3237,
        # Team to draft specific players
        2516, 2560, 2561, 2562, 2563, 2564, 2565, 2580, 2581, 4537,
        # Prop bundles
        4242, 4243, 4245, 4246,
        # Live
        1269,
    ]


def run_rest_phase(probe_only: int | None = None) -> tuple[dict | None, str | None, dict]:
    """Fetch all active NFL-Draft leagues in a single NewScheduleHelper call.

    WZ's endpoint accepts a comma-separated lg list, so one request returns all
    draft leagues at once. That's faster than probing each ID individually and
    matches how the real scraper will be wired (single multi-ID fetch).
    """
    print("Phase 1: REST fetch of WZ NFL Draft leagues...")
    load_env()
    username = os.getenv("WAGERZON_USERNAME")
    password = os.getenv("WAGERZON_PASSWORD")
    if not username or not password:
        print("  WAGERZON_USERNAME / WAGERZON_PASSWORD missing from .env")
        return None, None, {}

    session = requests.Session()
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
            "(KHTML, like Gecko) Chrome/146.0.0.0 Safari/537.36"
        ),
    })

    if not _wz_login(session, username, password):
        print("  WZ login failed. Credentials may be wrong or login page changed.")
        return None, None, {}
    print("  logged in")

    if probe_only:
        # Single-ID debug mode: just fetch that one, dump shape, return it.
        data = _wz_fetch(session, probe_only)
        comps = _find_competitions(data)
        if comps:
            for comp in comps[:5]:
                name = (comp.get("Description") or comp.get("desc"))
                print(f"  lg={probe_only} -> {name!r}")
        if not data:
            return None, None, {}
        return data, f"{WZ_HELPER_URL}?WT=0&lg={probe_only}", {
            "lg": probe_only,
            "method": "rest_probe_single",
        }

    # Build one comma-separated request for all candidate IDs.
    candidates = _candidate_lg_ids()
    lg_param = ",".join(str(x) for x in candidates)
    data = _wz_fetch(session, lg_param)
    if not data:
        print(f"  multi-ID fetch failed for {len(candidates)} lg values")
        return None, None, {}

    comps = _find_competitions(data)
    print(f"  response contains {len(comps)} league buckets")

    # Count matching vs total — logs which IDs actually came back with draft data.
    matches = []
    for comp in comps:
        if _looks_like_draft_competition(comp):
            matches.append(comp)
    print(f"  {len(matches)} leagues matched as active NFL-Draft markets:")
    for comp in matches:
        name = comp.get("Description") or comp.get("desc")
        idlg = comp.get("IdLeague")
        n_games = len(comp.get("Games", []) or [])
        print(f"    IdLeague={idlg} games={n_games} {name!r}")

    if not matches:
        # Still return the raw bundle — caller can inspect why no matches.
        print("  no leagues matched the draft heuristic. "
              "Raw response saved for manual inspection.")

    all_ids = [m.get("IdLeague") for m in matches if m.get("IdLeague") is not None]
    url = f"{WZ_HELPER_URL}?WT=0&lg={lg_param}"
    return data, url, {
        "lg_list": candidates,
        "lg_param": lg_param,
        "draft_league_ids": all_ids,
        "matched_league_count": len(matches),
        "total_league_count": len(comps),
        "method": "rest_multi_fetch",
    }


# ---------------------------------------------------------------------------
# Phase 2: Playwright capture
# ---------------------------------------------------------------------------

WZ_DRAFT_PAGE_CANDIDATES = [
    # We don't know the exact lg for draft, so we land on NFL and the user
    # clicks into the futures tab.
    "https://backend.wagerzon.com/wager/NewSchedule.aspx?WT=0&lg=4029,430,432,433",
    "https://backend.wagerzon.com/",
]


def run_browser_phase() -> tuple[dict | None, str | None, dict]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("  Playwright not installed.")
        return None, None, {}

    print("Phase 2: opening headed Chrome to capture WZ draft markets...")
    # Resolve against main repo root so the profile survives worktree swaps.
    profile = _main_repo_root() / "wagerzon_odds" / ".wagerzon_profile"
    profile.mkdir(parents=True, exist_ok=True)

    captured: list[tuple[str, bytes]] = []

    def on_response(response):
        url = response.url
        if "wagerzon.com" not in url:
            return
        if "NewScheduleHelper" not in url:
            return
        try:
            body = response.body()
        except Exception:
            return
        if len(body) < 300:
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

        for url in WZ_DRAFT_PAGE_CANDIDATES:
            print(f"  navigating to {url}")
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as exc:
                print(f"    nav failed: {exc!r}")
                continue
            page.wait_for_timeout(4000)

        print("\n  Log in if needed, then navigate to NFL Draft / NFL Props / NFL Futures.")
        print("  Press ENTER when draft markets are visible to save the capture...")
        try:
            input()
        except (EOFError, KeyboardInterrupt):
            pass

        ctx.close()

    if not captured:
        print("  No NewScheduleHelper responses captured.")
        return None, None, {}

    captured.sort(key=lambda t: len(t[1]), reverse=True)
    best_url, best_body = captured[0]
    try:
        data = json.loads(best_body.decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"  biggest response wasn't JSON: {exc!r}")
        return None, None, {}

    meta: dict = {"method": "playwright_capture", "url": best_url}
    # Extract lg from the captured URL's query string.
    m = re.search(r"[?&]lg=([^&]+)", best_url)
    if m:
        meta["lg"] = m.group(1)
        print(f"  discovered lg={m.group(1)} from captured URL")
    comps = _find_competitions(data)
    if comps:
        meta["competition_name"] = (comps[0].get("desc")
                                    or comps[0].get("Description")
                                    or comps[0].get("Name"))
    print(f"  captured {len(captured)} responses; saving largest ({len(best_body):,} bytes)")
    return data, best_url, meta


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WZ NFL Draft markets recon")
    parser.add_argument("--browser", action="store_true",
                        help="Skip probe; go straight to browser capture.")
    parser.add_argument("--probe", type=int, default=None,
                        help="Probe a single specific lg ID (for debugging).")
    args = parser.parse_args()

    print("=" * 60)
    print("  WAGERZON NFL DRAFT RECON")
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
        print("\nERROR: Could not capture WZ NFL Draft markets.")
        print("  Next steps:")
        print("  - Confirm WAGERZON_USERNAME / WAGERZON_PASSWORD are in bet_logger/.env")
        print("  - Re-run with --browser to navigate manually and capture the lg value")
        sys.exit(1)

    path = save_fixture("wagerzon", data, meta=meta)
    print_diagnostics("wagerzon", path, url, data)
    print(f"  method used: {meta.get('method', 'unknown')}")
    if "lg" in meta:
        print(f"  lg: {meta['lg']}  <- single-ID probe result")
    if "lg_param" in meta:
        print(f"  multi-ID request used: lg={meta['lg_param']}")
    if "draft_league_ids" in meta:
        print(f"  ACTIVE draft IdLeagues: {meta['draft_league_ids']}")
        print(f"  -> use these in wagerzon_odds/config.py for the draft sport")
    if "matched_league_count" in meta:
        print(f"  {meta['matched_league_count']} / {meta.get('total_league_count', '?')} "
              "buckets matched the active-draft heuristic")
