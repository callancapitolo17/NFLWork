#!/usr/bin/env python3
"""
Wagerzon Parlay Slip Recon
Captures network traffic when adding legs to bet slip and entering a bet amount.
Goal: find the API call that calculates parlay payout so we can replicate it.

Usage:
    python3 wagerzon_odds/recon_parlay_slip.py

Steps:
    1. Browser opens to MLB schedule (auto-login if profile exists)
    2. Manually add a spread leg to the bet slip
    3. Manually add a total leg (same game) to the bet slip
    4. Select "Parlay" in the slip
    5. Type a bet amount (e.g. $100)
    6. Press ENTER in terminal at each prompt
    7. Script shows all API calls captured during each step
"""

from playwright.sync_api import sync_playwright
import os
import json

WAGERZON_URL = "https://backend.wagerzon.com"
MLB_SCHEDULE_URL = f"{WAGERZON_URL}/wager/NewSchedule.aspx?WT=0&lg=417,1280"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")


def run_recon():
    # Collect requests per phase
    current_phase = {"name": "startup", "requests": []}
    all_phases = []

    def handle_request(request):
        if "wagerzon.com" not in request.url:
            return
        # Skip static assets
        if request.resource_type in ("image", "stylesheet", "font"):
            return

        entry = {
            "url": request.url,
            "method": request.method,
            "resource_type": request.resource_type,
            "post_data": request.post_data,
        }
        current_phase["requests"].append(entry)

    def handle_response(response):
        if "wagerzon.com" not in response.url:
            return
        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                try:
                    content_type = response.headers.get("content-type", "")
                    if "json" in content_type or "javascript" in content_type:
                        body = response.text()
                        entry["response_preview"] = body[:5000]
                        entry["response_size"] = len(body)
                    elif "html" in content_type:
                        body = response.text()
                        entry["response_size"] = len(body)
                        # Look for bet slip / payout related content
                        for keyword in ["payout", "parlay", "towin", "to_win",
                                        "wager", "betslip", "slip", "PayOut",
                                        "calcul", "TotalRisk", "TotalWin"]:
                            if keyword.lower() in body.lower():
                                # Extract surrounding context
                                idx = body.lower().find(keyword.lower())
                                snippet = body[max(0, idx-100):idx+200]
                                if "snippets" not in entry:
                                    entry["snippets"] = []
                                entry["snippets"].append({
                                    "keyword": keyword,
                                    "context": snippet
                                })
                except Exception:
                    pass
                break

    def start_phase(name):
        nonlocal current_phase
        if current_phase["requests"]:
            all_phases.append(current_phase)
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing: {name} ---")

    def show_phase_results():
        reqs = current_phase["requests"]
        if not reqs:
            print("  (no requests captured)")
            return

        for req in reqs:
            rtype = req.get("resource_type", "?")
            status = req.get("status", "?")
            method = req["method"]
            url = req["url"]

            # Highlight interesting requests
            is_api = rtype in ("xhr", "fetch")
            is_post = method == "POST"
            has_json = "response_preview" in req
            has_snippets = "snippets" in req

            prefix = "  "
            if is_api or is_post:
                prefix = "  *** "

            print(f"{prefix}[{status}] {method} {rtype} {url[:120]}")

            if req.get("post_data"):
                post = req["post_data"]
                if len(post) > 500:
                    print(f"       POST data: {post[:500]}...")
                else:
                    print(f"       POST data: {post}")

            if has_json:
                print(f"       JSON response ({req['response_size']} chars):")
                print(f"       {req['response_preview'][:500]}")

            if has_snippets:
                print(f"       SLIP-RELATED KEYWORDS FOUND:")
                for s in req["snippets"][:3]:
                    print(f"         '{s['keyword']}': ...{s['context'][:150]}...")

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

        # Phase 1: Navigate to MLB schedule
        start_phase("page_load")
        print(f"Navigating to {MLB_SCHEDULE_URL}...")
        page.goto(MLB_SCHEDULE_URL, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(3000)

        print("\n" + "=" * 60)
        print("If you need to log in, do so now.")
        print("Once you see the MLB schedule, press ENTER.")
        print("=" * 60)
        input()
        show_phase_results()

        # Phase 2: Add first leg (spread)
        start_phase("add_spread_leg")
        print("=" * 60)
        print("STEP 1: Click on a SPREAD bet to add it to the slip.")
        print("        (e.g., click a team's run line -1.5)")
        print("Press ENTER after clicking.")
        print("=" * 60)
        input()
        page.wait_for_timeout(1000)
        show_phase_results()

        # Phase 3: Add second leg (total)
        start_phase("add_total_leg")
        print("\n" + "=" * 60)
        print("STEP 2: Click on a TOTAL bet (same game) to add it.")
        print("        (e.g., click Over 8.5)")
        print("Press ENTER after clicking.")
        print("=" * 60)
        input()
        page.wait_for_timeout(1000)
        show_phase_results()

        # Phase 4: Select parlay
        start_phase("select_parlay")
        print("\n" + "=" * 60)
        print("STEP 3: Select 'Parlay' in the bet slip.")
        print("        (look for a Parlay/Teaser toggle)")
        print("Press ENTER after selecting.")
        print("=" * 60)
        input()
        page.wait_for_timeout(1000)
        show_phase_results()

        # Phase 5: Enter bet amount
        start_phase("enter_amount")
        print("\n" + "=" * 60)
        print("STEP 4: Type a bet amount (e.g., 100).")
        print("        Watch for the 'To Win' / payout to appear.")
        print("Press ENTER after entering the amount.")
        print("=" * 60)
        input()
        page.wait_for_timeout(2000)
        show_phase_results()

        # Phase 6: Check the slip HTML for payout info
        start_phase("inspect_slip")
        print("\n" + "=" * 60)
        print("Inspecting bet slip DOM for payout values...")
        print("=" * 60)

        # Try to find bet slip elements
        try:
            slip_html = page.evaluate("""() => {
                // Look for common bet slip containers
                const selectors = [
                    '#betslip', '.betslip', '.bet-slip',
                    '#slip', '.slip', '.wager-slip',
                    '[class*="slip"]', '[id*="slip"]',
                    '[class*="Slip"]', '[id*="Slip"]',
                    '[class*="wager"]', '[id*="wager"]',
                    '[class*="Wager"]', '[id*="Wager"]',
                    '[class*="parlay"]', '[id*="parlay"]',
                    '.betCard', '#betCard',
                ];
                for (const sel of selectors) {
                    const el = document.querySelector(sel);
                    if (el) {
                        return {selector: sel, html: el.innerHTML.substring(0, 5000)};
                    }
                }
                // Fallback: look for "to win" or "payout" text anywhere
                const body = document.body.innerHTML;
                const payoutIdx = body.toLowerCase().indexOf('to win');
                if (payoutIdx >= 0) {
                    return {selector: 'body (to win)', html: body.substring(payoutIdx - 200, payoutIdx + 500)};
                }
                const parlayIdx = body.toLowerCase().indexOf('parlay');
                if (parlayIdx >= 0) {
                    return {selector: 'body (parlay)', html: body.substring(parlayIdx - 200, parlayIdx + 500)};
                }
                return {selector: 'none found', html: ''};
            }""")
            print(f"  Found slip element: {slip_html['selector']}")
            if slip_html['html']:
                print(f"  HTML preview:\n{slip_html['html'][:2000]}")
        except Exception as e:
            print(f"  Error inspecting DOM: {e}")

        # Also try to find any JavaScript that calculates payouts
        try:
            js_calcs = page.evaluate("""() => {
                const scripts = document.querySelectorAll('script');
                const relevant = [];
                for (const s of scripts) {
                    const text = s.textContent || '';
                    if (text.match(/payout|parlay|towin|multiply|decimal|odds/i)) {
                        relevant.push(text.substring(0, 2000));
                    }
                }
                return relevant;
            }""")
            if js_calcs:
                print(f"\n  Found {len(js_calcs)} scripts with payout-related code:")
                for i, js in enumerate(js_calcs):
                    print(f"\n  --- Script {i+1} ---")
                    print(f"  {js[:1000]}")
        except Exception as e:
            print(f"  Error scanning scripts: {e}")

        # Save all phases
        all_phases.append(current_phase)
        out_path = os.path.join(os.path.dirname(__file__), "recon_parlay_slip.json")
        with open(out_path, "w") as f:
            json.dump(all_phases, f, indent=2, default=str)
        print(f"\nSaved all captured data to {out_path}")

        print("\n" + "=" * 60)
        print("Done! Press ENTER to close browser.")
        print("=" * 60)
        input()
        context.close()


if __name__ == "__main__":
    run_recon()
