"""
Hoop88 Recon Script
- Uses real Chrome to capture all network requests during login + navigation
- Discovers: auth flow, JSON API endpoints, data format, cookie/token mechanics
- Persistent profile so cookies survive between runs

Usage:
    hoop88_odds/venv/bin/python hoop88_odds/recon_hoop88.py

After running:
    1. Log in manually when prompted
    2. Click BASKETBALL → NCAA Basketball in the sidebar
    3. Click "1st Half" period tab
    4. Press ENTER at each prompt to advance
    5. Captured requests + cookies saved to hoop88_odds/
"""

from playwright.sync_api import sync_playwright
import os
import json

HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".hoop88_profile")
OUT_DIR = os.path.dirname(__file__)


def run_recon():
    captured_requests = []

    def handle_request(request):
        """Capture all requests to discover API endpoints."""
        url = request.url
        # Filter to the site's domain (skip CDN, analytics, etc.)
        if not any(domain in url for domain in ["hoop88", "h88"]):
            # Also capture if it looks like an API call from the site
            if request.resource_type not in ("xhr", "fetch"):
                return

        entry = {
            "url": url,
            "method": request.method,
            "headers": dict(request.headers),
            "post_data": request.post_data,
            "resource_type": request.resource_type,
        }
        captured_requests.append(entry)

        # Categorize for live output
        rtype = request.resource_type
        if rtype == "document":
            print(f"  [PAGE] {request.method} {url[:120]}")
            if request.method == "POST" and request.post_data:
                post = request.post_data
                if "customerID" in post or "Password" in post:
                    print(f"     LOGIN POST detected")
                print(f"     POST data length: {len(post)} chars")
        elif rtype in ("xhr", "fetch"):
            content_type = request.headers.get("content-type", "")
            accept = request.headers.get("accept", "")
            if "json" in accept or "json" in content_type:
                print(f"  [API]  {request.method} {url[:120]}")
            else:
                print(f"  [XHR]  {request.method} {url[:120]}")
            # Show POST body for API calls
            if request.method == "POST" and request.post_data:
                preview = request.post_data[:300]
                print(f"     Body: {preview}")
        elif rtype == "script":
            print(f"  [JS]   {url[:80]}")
        elif rtype == "websocket":
            print(f"  [WS]   {url[:120]}")

    def handle_response(response):
        """Capture response details, especially JSON bodies."""
        url = response.url
        # Match same filter as requests
        if not any(domain in url for domain in ["hoop88", "h88"]):
            # Check if this was captured as xhr/fetch
            matching = [e for e in captured_requests if e["url"] == url and "status" not in e]
            if not matching:
                return

        for entry in captured_requests:
            if entry["url"] == url and "status" not in entry:
                entry["status"] = response.status
                entry["response_headers"] = dict(response.headers)
                try:
                    content_type = response.headers.get("content-type", "")
                    if "json" in content_type:
                        body = response.text()
                        entry["response_body_preview"] = body[:10000]
                        entry["response_body_size"] = len(body)
                    elif "html" in content_type and entry["resource_type"] == "document":
                        body = response.text()
                        entry["response_body_size"] = len(body)
                        # Check for login form
                        if "customerID" in body and "Password" in body:
                            entry["page_type"] = "login"
                        elif "data-panel" in body:
                            entry["page_type"] = "app"
                            entry["has_structured_data"] = True
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

        # Attach network listeners BEFORE navigating
        print("Attaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 1: Navigate to Hoop88
        print(f"\nNavigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until="domcontentloaded", timeout=60000)

        print("\n" + "=" * 60)
        print("1. If you see a login form, log in manually")
        print("2. Once logged in, press ENTER")
        print("=" * 60)
        input()

        current_url = page.url
        print(f"\nCurrent URL: {current_url}")
        print(f"Captured {len(captured_requests)} requests so far")

        # Step 2: Navigate to Basketball
        print("\n" + "=" * 60)
        print("Now click BASKETBALL in the sidebar.")
        print("Wait for games to load, then press ENTER.")
        print("=" * 60)
        input()

        print(f"\nCaptured {len(captured_requests)} requests so far")

        # Step 3: Navigate to NCAA Basketball
        print("\n" + "=" * 60)
        print("Now click 'NCAA Basketball' in the sidebar.")
        print("Wait for games to load, then press ENTER.")
        print("=" * 60)
        input()

        print(f"\nCaptured {len(captured_requests)} requests so far")

        # Step 4: Click 1st Half period tab
        print("\n" + "=" * 60)
        print("Now click '1st Half' in the sidebar (under NCAA Basketball).")
        print("Wait for games to load, then press ENTER.")
        print("=" * 60)
        input()

        print(f"\nCaptured {len(captured_requests)} requests so far")

        # Step 5: Try other navigation (NBA, different periods)
        print("\n" + "=" * 60)
        print("Explore freely: try NBA, different periods, team totals, etc.")
        print("Each click may trigger new API calls.")
        print("Press ENTER when done exploring.")
        print("=" * 60)
        input()

        page.wait_for_timeout(2000)

        # Save captured requests
        api_path = os.path.join(OUT_DIR, "recon_hoop88_api.json")
        with open(api_path, "w", encoding="utf-8") as f:
            json.dump(captured_requests, f, indent=2, default=str)
        print(f"\nSaved {len(captured_requests)} requests to {api_path}")

        # Print summary
        print(f"\n{'=' * 60}")
        print("REQUEST SUMMARY")
        print(f"{'=' * 60}")

        page_requests = []
        xhr_requests = []
        other_requests = []

        for req in captured_requests:
            rtype = req.get("resource_type", "")
            if rtype == "document":
                page_requests.append(req)
            elif rtype in ("xhr", "fetch"):
                xhr_requests.append(req)
            else:
                other_requests.append(req)

        if page_requests:
            print("\n--- PAGE REQUESTS ---")
            for req in page_requests:
                status = req.get("status", "?")
                page_type = req.get("page_type", "")
                size = req.get("response_body_size", "?")
                print(f"  [{status}] {req['method']} {req['url'][:120]}")
                if page_type:
                    print(f"     Page type: {page_type}")
                if size != "?":
                    print(f"     Size: {size} chars")

        if xhr_requests:
            print(f"\n--- XHR/FETCH REQUESTS ({len(xhr_requests)} total) ---")
            print("*** These are the API endpoints! ***")
            for req in xhr_requests:
                status = req.get("status", "?")
                size = req.get("response_body_size", "?")
                print(f"  [{status}] {req['method']} {req['url']}")
                if req.get("post_data"):
                    preview = req["post_data"][:200]
                    print(f"     Request body: {preview}")
                if size != "?":
                    print(f"     Response: {size} chars")
                if req.get("response_body_preview"):
                    preview = req["response_body_preview"][:300]
                    print(f"     Response body: {preview}")
                # Show auth headers
                auth = req.get("headers", {}).get("authorization", "")
                if auth:
                    print(f"     Auth header: {auth[:80]}...")
                token = req.get("headers", {}).get("x-token", "")
                if token:
                    print(f"     X-Token: {token[:80]}...")
        else:
            print("\n--- No XHR/fetch requests detected ---")
            print("Site may use WebSockets or server-rendered HTML")

        if other_requests:
            print(f"\n--- OTHER ({len(other_requests)} requests) ---")
            ws = [r for r in other_requests if "websocket" in r.get("resource_type", "")]
            if ws:
                print(f"  WebSocket connections: {len(ws)}")
                for r in ws:
                    print(f"    {r['url'][:120]}")
            script_count = sum(1 for r in other_requests if r.get("resource_type") == "script")
            print(f"  {script_count} scripts, {len(other_requests) - script_count - len(ws)} other")

        # Save cookies
        cookies = context.cookies()
        cookie_path = os.path.join(OUT_DIR, "recon_hoop88_cookies.json")
        with open(cookie_path, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)
        print(f"\nSaved {len(cookies)} cookies to {cookie_path}")

        # Key findings
        print(f"\n{'=' * 60}")
        print("KEY FINDINGS FOR SCRAPER")
        print(f"{'=' * 60}")

        # Auth analysis
        auth_cookies = [c for c in cookies if any(k in c["name"].lower()
                        for k in ["session", "auth", "token", "jwt", "sid"])]
        if auth_cookies:
            print("\nAuth-related cookies:")
            for c in auth_cookies:
                val = c["value"][:60] + "..." if len(c["value"]) > 60 else c["value"]
                print(f"  {c['name']} = {val} (domain: {c['domain']})")
        else:
            print("\nAll cookies:")
            for c in cookies:
                val = c["value"][:60] + "..." if len(c["value"]) > 60 else c["value"]
                print(f"  {c['name']} = {val} (domain: {c['domain']})")

        # API endpoint summary
        if xhr_requests:
            unique_urls = set()
            for req in xhr_requests:
                # Strip query params for grouping
                base_url = req["url"].split("?")[0]
                unique_urls.add(base_url)
            print(f"\nUnique API endpoints: {len(unique_urls)}")
            for url in sorted(unique_urls):
                count = sum(1 for r in xhr_requests if r["url"].split("?")[0] == url)
                print(f"  [{count}x] {url}")

            # Check for auth headers in API calls
            auth_headers = set()
            for req in xhr_requests:
                headers = req.get("headers", {})
                for key in ["authorization", "x-token", "x-auth", "x-session"]:
                    if key in headers:
                        auth_headers.add(key)
            if auth_headers:
                print(f"\nAuth headers used in API calls: {auth_headers}")
            else:
                print("\nNo auth headers found — API may use cookie-based auth")

        print("\nPress ENTER to close browser.")
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
