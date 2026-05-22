#!/usr/bin/env python3
"""
Wagerzon Place-Single (F-Period) Recon

Goal: capture the EXACT HTTP request Wagerzon's browser sends when "Place Bet"
is clicked on a confirmed SINGLE bet at an F-period market (F3 / F7 total or
spread). The parlay recon (recon_place_parlay.py) gave us the parlay shape,
but F-period singles fail with GAMELINECHANGE because F3/F7 share an idgm
with the FG market — we need to see what extra field/play-code/header the
browser sends to disambiguate "I want the F3 total, not the FG total."

What we're trying to learn:
    - Does play=3 mean "FG under" only, with F3/F7 under needing a different
      play code (e.g. play=33, 43, 53)?
    - Does detailData need an extra field (Period, IdGameType, etc.) we're not
      sending today?
    - Is the disambiguation in sel's encoding (a fifth segment?) or
      elsewhere entirely?
    - The wire form for CreateWagerHelper / ConfirmWagerHelper / PostWager-
      MultipleHelper for an F-period leg.

Usage:
    cd /Users/callancapitolo/NFLWork
    source wagerzon_odds/venv/bin/activate
    python3 wagerzon_odds/recon_place_single_fperiod.py

Walkthrough:
    1. Browser opens to Wagerzon (persistent profile in .wagerzon_profile)
    2. Log in if needed, navigate to MLB schedule
    3. Build a SINGLE bet on a FIRST-3-INNINGS or FIRST-7-INNINGS total
       (or spread) — minimum stake (typically $1).
    4. Get to the page where "Place Bet" is the next click.
    5. Press ENTER in terminal → phase boundary marker.
    6. Click "Place Bet" → wait for confirmation → press ENTER.
    7. Script saves recon_place_single_fperiod.json with full request /
       response data plus a summary of the three key POSTs we care about:
       CreateWagerHelper, ConfirmWagerHelper, PostWagerMultipleHelper.

NOTE: This places a REAL bet. Use Wagerzon's minimum stake. Choose an
F-period market deliberately (the whole point is to capture that wire shape).
"""

import json
import os
from playwright.sync_api import sync_playwright

WAGERZON_URL = "https://backend.wagerzon.com"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")
OUTPUT_PATH = os.path.join(os.path.dirname(__file__),
                            "recon_place_single_fperiod.json")

# Endpoints whose behavior we already understand. Any POST hitting an endpoint
# NOT in this list during the click_place_bet phase is a placement candidate.
KNOWN_ENDPOINTS = (
    "NewScheduleHelper",
    "ConfirmWagerHelper",
    "CreateWagerHelper",
    "HistoryHelper",
    "Login",
    "Welcome",
    "Logout",
    "Root.aspx",
    ".css",
    ".js",
    ".woff",
    ".png",
    ".jpg",
    ".gif",
    ".svg",
    "google",
    "doubleclick",
    "analytics",
)

# The three POSTs we specifically care about; we'll pretty-print these in
# the final summary even though they're "known" endpoints.
KEY_PLACEMENT_POSTS = (
    "CreateWagerHelper",
    "ConfirmWagerHelper",
    "PostWagerMultipleHelper",
)


def is_placement_candidate(url: str) -> bool:
    """Unknown POST → potential placement endpoint we hadn't documented."""
    return not any(k in url for k in KNOWN_ENDPOINTS)


def is_key_post(url: str) -> bool:
    return any(k in url for k in KEY_PLACEMENT_POSTS)


def run_recon():
    current_phase = {"name": "startup", "requests": []}
    all_phases = []

    def handle_request(request):
        if "wagerzon.com" not in request.url:
            return
        if request.resource_type in ("image", "stylesheet", "font", "media"):
            return

        current_phase["requests"].append({
            "url": request.url,
            "method": request.method,
            "resource_type": request.resource_type,
            "headers": dict(request.headers),
            "post_data": request.post_data,
        })

    def handle_response(response):
        if "wagerzon.com" not in response.url:
            return
        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                entry["response_headers"] = dict(response.headers)
                try:
                    ct = response.headers.get("content-type", "")
                    if any(t in ct for t in ("json", "javascript", "text", "html")):
                        body = response.text()
                        entry["response_body"] = body[:15000]
                        entry["response_size"] = len(body)
                except Exception:
                    pass
                break

    def start_phase(name: str):
        nonlocal current_phase
        if current_phase["requests"]:
            all_phases.append(current_phase)
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing phase: {name} ---")

    def show_phase_summary(highlight: bool = False):
        reqs = [r for r in current_phase["requests"] if r["method"] != "GET"]
        if not reqs:
            print("  (no POST requests captured this phase)")
            return
        for req in reqs:
            url = req["url"]
            status = req.get("status", "?")
            mark = ""
            if highlight:
                if is_key_post(url):
                    mark = "  *** KEY POST ***"
                elif is_placement_candidate(url):
                    mark = "  *** UNKNOWN PLACEMENT CANDIDATE ***"
            print(f"  [{status}] {req['method']} {url[:140]}{mark}")
            if req.get("post_data"):
                print(f"       BODY: {req['post_data'][:600]}")
            if req.get("response_body"):
                print(f"       RESP: {req['response_body'][:400]}")

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()
        page.on("request", handle_request)
        page.on("response", handle_response)

        # ---- Phase 1: setup ----
        start_phase("page_load")
        print(f"Opening {WAGERZON_URL} ...")
        page.goto(WAGERZON_URL, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(2000)

        print("\n" + "=" * 78)
        print("STEP 1: Log in if needed. Navigate to the MLB schedule.")
        print("        Build a SINGLE bet at an F-PERIOD market:")
        print("           - 1st 3 innings total (e.g. Under 2.5)  ← preferred")
        print("           - 1st 7 innings total (e.g. Over 6.5)")
        print("           - 1st 3 / 1st 7 spread (also acceptable)")
        print("        Use Wagerzon's MINIMUM stake ($1).")
        print("        Stop when you are STARING AT 'Place Bet' as the next click.")
        print()
        print("Press ENTER when ready to capture the placement click.")
        print("=" * 78)
        input()
        show_phase_summary()

        # ---- Phase 2: the click ----
        start_phase("click_place_bet")
        print("\n" + "=" * 78)
        print("STEP 2: Click 'Place Bet' NOW. Wait for the confirmation page.")
        print("        DO NOT navigate away. Just wait until you see a result.")
        print("        If WZ refuses ('line moved' etc.), that's still useful data —")
        print("        the request/response is what we want to see.")
        print()
        print("Press ENTER after the bet is placed (or rejected).")
        print("=" * 78)
        input()
        page.wait_for_timeout(2000)
        show_phase_summary(highlight=True)

        # ---- Phase 3: post-confirmation noise ----
        start_phase("post_confirmation")
        print("\n" + "=" * 78)
        print("STEP 3: Capturing any follow-up requests (history refresh etc.)")
        print("Press ENTER when done.")
        print("=" * 78)
        input()
        show_phase_summary()

        # Save full capture
        all_phases.append(current_phase)
        with open(OUTPUT_PATH, "w") as f:
            json.dump(all_phases, f, indent=2, default=str)
        print(f"\nFull capture saved to: {OUTPUT_PATH}")

        # Final highlight of the three KEY POSTs we actually care about
        print("\n" + "=" * 78)
        print("KEY F-PERIOD PLACEMENT POSTs (CreateWagerHelper / ConfirmWagerHelper")
        print("/ PostWagerMultipleHelper, with full request body + response):")
        print("=" * 78)
        for phase in all_phases:
            if phase["name"] != "click_place_bet":
                continue
            for req in phase["requests"]:
                if req["method"] == "POST" and is_key_post(req["url"]):
                    endpoint = req["url"].rsplit("/", 1)[-1]
                    print(f"\n--- {endpoint} ---")
                    print(f"    Status: {req.get('status', '?')}")
                    print(f"    Body:   {(req.get('post_data') or '')[:1000]}")
                    if req.get("response_body"):
                        print(f"    Resp:   {req['response_body'][:1000]}")

        # Also surface anything we DIDN'T expect (e.g. a new endpoint specific
        # to F-period singles)
        unknown = [
            req for phase in all_phases
            if phase["name"] == "click_place_bet"
            for req in phase["requests"]
            if req["method"] == "POST" and is_placement_candidate(req["url"])
        ]
        if unknown:
            print("\n" + "=" * 78)
            print("UNKNOWN POSTs (endpoints we hadn't documented before):")
            print("=" * 78)
            for req in unknown:
                print(f"\n  URL:    {req['url']}")
                print(f"  Body:   {(req.get('post_data') or '')[:500]}")
                if req.get("response_body"):
                    print(f"  Resp:   {req['response_body'][:500]}")

        print("\n" + "=" * 78)
        print(f"Done. Send {OUTPUT_PATH} to Claude for analysis.")
        print("Press ENTER to close browser.")
        print("=" * 78)
        input()
        context.close()


if __name__ == "__main__":
    run_recon()
