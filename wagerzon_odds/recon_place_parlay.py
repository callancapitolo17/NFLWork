#!/usr/bin/env python3
"""
Wagerzon Place-Bet Recon

Goal: capture the HTTP request Wagerzon sends when "Place Bet" / "Submit Wager"
is clicked on a confirmed parlay slip. This is the missing piece for fully-API
parlay placement (the placement endpoint, headers, body, and response shape).

Usage:
    cd /Users/callancapitolo/NFLWork
    source wagerzon_odds/venv/bin/activate
    python3 wagerzon_odds/recon_place_parlay.py

Walkthrough:
    1. Browser opens to Wagerzon (persistent profile in .wagerzon_profile)
    2. Log in if needed, navigate to MLB schedule
    3. Build a small parlay (use Wagerzon's minimum, typically $1)
    4. Get to the page where "Place Bet" is the next click
    5. Press ENTER in terminal -> phase boundary marker
    6. Click "Place Bet" -> wait for confirmation -> press ENTER
    7. Script saves recon_place_parlay.json with full request/response data
       and prints the placement candidates

NOTE: This places a REAL bet. Use Wagerzon's minimum stake.
"""

import json
import os
from playwright.sync_api import sync_playwright

WAGERZON_URL = "https://backend.wagerzon.com"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")
OUTPUT_PATH = os.path.join(os.path.dirname(__file__), "recon_place_parlay.json")

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


def is_placement_candidate(url: str) -> bool:
    """Unknown POST -> potential placement endpoint."""
    return not any(k in url for k in KNOWN_ENDPOINTS)


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
                        entry["response_body"] = body[:10000]
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
            if highlight and is_placement_candidate(url):
                mark = "  *** PLACEMENT CANDIDATE ***"
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

        print("\n" + "=" * 70)
        print("STEP 1: Log in if needed. Navigate to MLB schedule.")
        print("        Build a small parlay (use Wagerzon's minimum stake).")
        print("        Stop when you are STARING AT 'Place Bet' / 'Submit Wager'")
        print("        as the next click.")
        print("Press ENTER when ready to capture the placement click.")
        print("=" * 70)
        input()
        show_phase_summary()

        # ---- Phase 2: the click ----
        start_phase("click_place_bet")
        print("\n" + "=" * 70)
        print("STEP 2: Click 'Place Bet' NOW. Wait for the confirmation page.")
        print("        DO NOT navigate away. Just wait until you see a result.")
        print("Press ENTER after the bet is placed.")
        print("=" * 70)
        input()
        page.wait_for_timeout(2000)
        show_phase_summary(highlight=True)

        # ---- Phase 3: post-confirmation noise ----
        start_phase("post_confirmation")
        print("\n" + "=" * 70)
        print("STEP 3: Capturing any follow-up requests (history refresh etc.)")
        print("Press ENTER when done.")
        print("=" * 70)
        input()
        show_phase_summary()

        # Save full capture
        all_phases.append(current_phase)
        with open(OUTPUT_PATH, "w") as f:
            json.dump(all_phases, f, indent=2, default=str)
        print(f"\nFull capture saved to: {OUTPUT_PATH}")

        # Final highlight
        print("\n" + "=" * 70)
        print("PLACEMENT CANDIDATES (unknown POSTs from click_place_bet):")
        print("=" * 70)
        candidates = []
        for phase in all_phases:
            if phase["name"] != "click_place_bet":
                continue
            for req in phase["requests"]:
                if req["method"] == "POST" and is_placement_candidate(req["url"]):
                    candidates.append(req)

        if not candidates:
            print("\n  (none found - check the saved JSON manually)")
        else:
            for i, req in enumerate(candidates, 1):
                print(f"\n[{i}] URL:    {req['url']}")
                print(f"    Status: {req.get('status', '?')}")
                print(f"    Body:   {(req.get('post_data') or '')[:500]}")
                if req.get("response_body"):
                    print(f"    Resp:   {req['response_body'][:500]}")
                # Show non-trivial headers (omit user-agent etc.)
                interesting = {
                    k: v for k, v in (req.get("headers") or {}).items()
                    if k.lower() in (
                        "x-requested-with", "x-csrf-token", "x-xsrf-token",
                        "__requestverificationtoken", "authorization",
                        "content-type", "origin", "referer",
                    )
                }
                if interesting:
                    print(f"    Headers (interesting): {interesting}")

        print("\n" + "=" * 70)
        print("Done. Press ENTER to close browser.")
        print("=" * 70)
        input()
        context.close()


if __name__ == "__main__":
    run_recon()
