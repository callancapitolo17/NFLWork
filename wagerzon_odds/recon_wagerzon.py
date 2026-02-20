"""
Wagerzon Recon Script
- Uses real Chrome to capture all network requests during login + page load
- Discovers: ASP.NET hidden form fields, cookies, any hidden XHR/AJAX APIs
- Persistent profile so cookies survive between runs

Usage:
    python3 wagerzon_odds/recon_wagerzon.py

After running:
    1. Log in manually when prompted
    2. Schedule page loads automatically (CBB)
    3. Press ENTER at each prompt to advance
    4. Captured requests + cookies saved to wagerzon_odds/
"""

from playwright.sync_api import sync_playwright
import os
import json

WAGERZON_URL = "https://backend.wagerzon.com"
SCHEDULE_URL = f"{WAGERZON_URL}/wager/NewSchedule.aspx?WT=0&lg=43,403,45"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")
OUT_DIR = os.path.dirname(__file__)


def run_recon():
    captured_requests = []

    def handle_request(request):
        """Capture ALL requests (form POSTs, XHR, fetch, document loads)."""
        url = request.url
        if "wagerzon.com" not in url:
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
            # For POST requests, show form data (login, ASP.NET postbacks)
            if request.method == "POST" and request.post_data:
                # Truncate but show key fields
                post = request.post_data
                if "__VIEWSTATE" in post:
                    print(f"     ASP.NET form POST (has __VIEWSTATE)")
                if "Account" in post or "Password" in post or "BtnSubmit" in post:
                    print(f"     LOGIN form POST detected")
                print(f"     POST data length: {len(post)} chars")
        elif rtype in ("xhr", "fetch"):
            if "json" in request.headers.get("accept", ""):
                print(f"  [API]  {request.method} {url[:120]}")
            else:
                print(f"  [XHR]  {request.method} {url[:120]}")
        elif rtype == "script":
            print(f"  [JS]   {request.method} {url[:80]}")

    def handle_response(response):
        """Capture response details."""
        if "wagerzon.com" not in response.url:
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
                    elif "html" in content_type and entry["resource_type"] == "document":
                        body = response.text()
                        entry["response_body_size"] = len(body)
                        # Check for key ASP.NET form fields
                        aspnet_fields = []
                        for field in ["__VIEWSTATE", "__EVENTVALIDATION",
                                      "__VIEWSTATEGENERATOR", "__EVENTTARGET"]:
                            if field in body:
                                aspnet_fields.append(field)
                        if aspnet_fields:
                            entry["aspnet_fields_found"] = aspnet_fields
                        # Check if it's a login page or schedule page
                        if "Account" in body and "Password" in body:
                            entry["page_type"] = "login"
                        elif "Competition" in body and "GameRow" in body:
                            entry["page_type"] = "schedule"
                            entry["has_odds_data"] = True
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

        # Attach network listeners BEFORE navigating (capture login flow)
        print("Attaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 1: Navigate to Wagerzon login/schedule
        print(f"\nNavigating to {SCHEDULE_URL}...")
        page.goto(SCHEDULE_URL, wait_until="domcontentloaded", timeout=60000)

        print("\n" + "=" * 60)
        print("1. If redirected to login, log in manually")
        print("2. Once you see the schedule page with games, press ENTER")
        print("=" * 60)
        input()

        # Check if we're on the schedule page now
        current_url = page.url
        print(f"\nCurrent URL: {current_url}")
        print(f"Captured {len(captured_requests)} requests so far")

        # Step 2: If already on schedule, reload to capture a clean page load
        print("\n" + "=" * 60)
        print("Now reloading the schedule page to capture a clean request flow...")
        print("=" * 60)
        page.goto(SCHEDULE_URL, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(5000)

        print(f"\nTotal captured: {len(captured_requests)} requests")

        # Step 3: Try navigating to a different sport to see URL patterns
        print("\n" + "=" * 60)
        print("Try clicking different leagues/sports in the sidebar if available.")
        print("Or navigate to NBA: NewSchedule.aspx?WT=0&lg=3,301,303")
        print("Press ENTER when done exploring.")
        print("=" * 60)
        input()

        page.wait_for_timeout(2000)

        # Step 4: Save captured requests
        api_path = os.path.join(OUT_DIR, "recon_wagerzon_api.json")
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

        # Document/page requests (most important for ASP.NET)
        if page_requests:
            print("\n--- PAGE REQUESTS (document loads + form POSTs) ---")
            for req in page_requests:
                status = req.get("status", "?")
                method = req["method"]
                page_type = req.get("page_type", "")
                aspnet = req.get("aspnet_fields_found", [])
                has_odds = req.get("has_odds_data", False)
                size = req.get("response_body_size", "?")

                print(f"  [{status}] {method} {req['url'][:120]}")
                if page_type:
                    print(f"     Page type: {page_type}")
                if aspnet:
                    print(f"     ASP.NET fields: {aspnet}")
                if has_odds:
                    print(f"     Contains odds data: YES")
                if size != "?":
                    print(f"     Response size: {size} chars")
                if method == "POST" and req.get("post_data"):
                    # Show form field names (not values - may contain credentials)
                    post = req["post_data"]
                    fields = []
                    for part in post.split("&"):
                        if "=" in part:
                            key = part.split("=")[0]
                            fields.append(key)
                    print(f"     Form fields: {fields[:20]}")

        # XHR/fetch requests (hidden APIs!)
        if xhr_requests:
            print(f"\n--- XHR/FETCH REQUESTS ({len(xhr_requests)} total) ---")
            print("*** These indicate hidden APIs! ***")
            for req in xhr_requests:
                status = req.get("status", "?")
                size = req.get("response_body_size", "?")
                print(f"  [{status}] {req['method']} {req['url']}")
                if size != "?":
                    print(f"     Response: {size} chars")
                if req.get("response_body_preview"):
                    preview = req["response_body_preview"][:200]
                    print(f"     Body: {preview}")
        else:
            print("\n--- No XHR/fetch requests detected ---")
            print("This confirms Wagerzon is pure server-rendered HTML (ASP.NET)")

        if other_requests:
            print(f"\n--- OTHER ({len(other_requests)} requests: scripts, css, etc.) ---")
            script_count = sum(1 for r in other_requests if r.get("resource_type") == "script")
            print(f"  {script_count} script loads, {len(other_requests) - script_count} other")

        # Save cookies
        cookies = context.cookies()
        cookie_path = os.path.join(OUT_DIR, "recon_wagerzon_cookies.json")
        with open(cookie_path, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)
        print(f"\nSaved {len(cookies)} cookies to {cookie_path}")

        # Key findings for scraper implementation
        print(f"\n{'=' * 60}")
        print("KEY FINDINGS FOR SCRAPER")
        print(f"{'=' * 60}")

        # Auth cookies
        auth_cookies = [c for c in cookies if "session" in c["name"].lower()
                        or "auth" in c["name"].lower()
                        or "asp" in c["name"].lower()]
        if auth_cookies:
            print("\nAuth-related cookies:")
            for c in auth_cookies:
                print(f"  {c['name']} = {c['value'][:50]}... (domain: {c['domain']})")
        else:
            print("\nAll cookies:")
            for c in cookies:
                print(f"  {c['name']} = {c['value'][:50]}... (domain: {c['domain']})")

        # Login form analysis
        login_posts = [r for r in page_requests
                       if r["method"] == "POST" and r.get("post_data")
                       and ("Account" in r["post_data"] or "Password" in r["post_data"])]
        if login_posts:
            print("\nLogin POST detected:")
            for req in login_posts:
                print(f"  URL: {req['url']}")
                fields = []
                for part in req["post_data"].split("&"):
                    if "=" in part:
                        key = part.split("=")[0]
                        fields.append(key)
                print(f"  Fields: {fields}")

        # Schedule page analysis
        schedule_pages = [r for r in page_requests
                          if r.get("has_odds_data")]
        if schedule_pages:
            print(f"\nSchedule pages with odds data: {len(schedule_pages)}")
            for req in schedule_pages:
                print(f"  URL: {req['url'][:120]}")
                print(f"  Size: {req.get('response_body_size', '?')} chars")

        # Hidden API summary
        if xhr_requests:
            print(f"\nHIDDEN APIs FOUND: {len(xhr_requests)}")
            print("These can potentially replace HTML scraping entirely!")
        else:
            print("\nNo hidden APIs - scraper should use requests + BeautifulSoup")
            print("Login via form POST, then GET schedule pages")

        print("\nPress ENTER to close browser.")
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
