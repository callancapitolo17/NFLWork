"""
BFA Recon Script
- Uses real Chrome to capture API requests and auth tokens
- Persistent profile so cookies survive between runs
- Captures all XHR/fetch to api.bfagaming.com and auth.bfagaming.com
- Saves cookies (including Keycloak refresh token) for browser-free scraper

Usage:
    python3 bfa_odds/recon_bfa.py

After running:
    1. Log in manually when prompted
    2. Navigate to each sport (NFL, NBA, CBB) and click into a game detail
    3. Press ENTER at each prompt to advance
    4. Cookies and API captures are saved to bfa_odds/
"""

from playwright.sync_api import sync_playwright
import os
import json

BFA_URL = "https://bfagaming.com"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".bfa_profile")
OUT_DIR = os.path.dirname(__file__)


def run_recon():
    captured_requests = []

    def handle_request(request):
        """Capture API/XHR requests to bfagaming.com domains."""
        if request.resource_type in ("xhr", "fetch"):
            url = request.url
            # Only capture relevant API calls
            if "bfagaming.com" not in url:
                return
            entry = {
                "url": url,
                "method": request.method,
                "headers": dict(request.headers),
                "post_data": request.post_data,
                "resource_type": request.resource_type,
            }
            captured_requests.append(entry)
            # Highlight auth vs odds requests
            if "token" in url or "auth" in url:
                print(f"  [AUTH] {request.method} {url[:120]}")
            elif "oddsservice" in url:
                print(f"  [ODDS] {request.method} {url[:120]}")
            else:
                print(f"  [API]  {request.method} {url[:120]}")

    def handle_response(response):
        """Capture response details for API calls."""
        if response.request.resource_type not in ("xhr", "fetch"):
            return
        if "bfagaming.com" not in response.url:
            return

        for entry in captured_requests:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                entry["response_headers"] = dict(response.headers)
                try:
                    content_type = response.headers.get("content-type", "")
                    if "json" in content_type:
                        body = response.text()
                        entry["response_body_preview"] = body[:10000]
                        entry["response_body_size"] = len(body)
                except Exception:
                    pass
                break

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()

        # Step 1: Navigate to BFA
        print("Navigating to BFA...")
        page.goto(BFA_URL, wait_until="domcontentloaded", timeout=60000)

        print("\n" + "=" * 60)
        print("1. Log in if needed (Keycloak login page)")
        print("2. Once you see the main BFA site, press ENTER here")
        print("=" * 60)
        input()

        # Step 2: Attach network listeners
        print("\nAttaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 3: Navigate to sportsbook
        print("Navigating to sportsbook...")
        page.goto(f"{BFA_URL}/sports", wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(5000)

        print(f"\nCaptured {len(captured_requests)} requests so far")
        print("\n" + "=" * 60)
        print("Now click on a sport (e.g., NCAA(B)) and then click")
        print("'More wagers' on a game to load the detail page.")
        print("This captures the events/popular and event/{id} API calls.")
        print("\nClick through multiple sports if you want to discover all sport keys.")
        print("Press ENTER when done exploring.")
        print("=" * 60)
        input()

        page.wait_for_timeout(3000)

        # Step 4: Save captured requests
        api_path = os.path.join(OUT_DIR, "recon_bfa_api.json")
        with open(api_path, "w", encoding="utf-8") as f:
            json.dump(captured_requests, f, indent=2, default=str)
        print(f"\nSaved {len(captured_requests)} API requests to {api_path}")

        # Print summary
        print(f"\n{'=' * 60}")
        print("API REQUEST SUMMARY")
        print(f"{'=' * 60}")

        auth_requests = []
        odds_requests = []
        other_requests = []

        for req in captured_requests:
            url = req["url"]
            if "token" in url or "auth" in url:
                auth_requests.append(req)
            elif "oddsservice" in url:
                odds_requests.append(req)
            else:
                other_requests.append(req)

        if auth_requests:
            print("\n--- AUTH REQUESTS ---")
            for req in auth_requests:
                status = req.get("status", "?")
                print(f"  [{status}] {req['method']} {req['url'][:120]}")
                if req.get("post_data"):
                    print(f"     POST: {req['post_data'][:200]}")

        if odds_requests:
            print("\n--- ODDS API REQUESTS ---")
            for req in odds_requests:
                status = req.get("status", "?")
                size = req.get("response_body_size", "?")
                print(f"  [{status}] {req['method']} {req['url']}")
                if size != "?":
                    print(f"     Response: {size} chars")

        if other_requests:
            print(f"\n--- OTHER ({len(other_requests)} requests) ---")
            for req in other_requests[:10]:
                status = req.get("status", "?")
                print(f"  [{status}] {req['method']} {req['url'][:120]}")

        # Save cookies
        cookies = context.cookies()
        cookie_path = os.path.join(OUT_DIR, "recon_bfa_cookies.json")
        with open(cookie_path, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)
        print(f"\nSaved {len(cookies)} cookies to {cookie_path}")

        # Extract key info for the scraper
        print(f"\n{'=' * 60}")
        print("KEY INFO FOR SCRAPER")
        print(f"{'=' * 60}")

        # Find token endpoint and client_id
        for req in auth_requests:
            if req.get("post_data") and "refresh_token" in req["post_data"]:
                print(f"\nToken URL: {req['url']}")
                print(f"POST data: {req['post_data'][:300]}")

        # Find sport keys
        sport_keys = set()
        for req in odds_requests:
            if "events/popular/" in req["url"]:
                key = req["url"].split("events/popular/")[1].split("?")[0]
                sport_keys.add(key)
                print(f"\nSport key: {key}")
                print(f"  URL: {req['url']}")

        # Find player/agent IDs
        for req in odds_requests[:1]:
            url = req["url"]
            if "playerId=" in url:
                pid = url.split("playerId=")[1].split("&")[0]
                print(f"\nplayerId: {pid}")
            if "agentId=" in url:
                aid = url.split("agentId=")[1].split("&")[0]
                print(f"agentId: {aid}")

        print("\nPress ENTER to close browser.")
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
