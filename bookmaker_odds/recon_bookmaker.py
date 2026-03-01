"""
Bookmaker.eu Recon Script
- Uses real Chrome (not Playwright Chromium) to bypass Cloudflare
- Persistent profile so cookies survive between runs
- Captures all network API requests on odds pages
- Saves requests to recon_bookmaker_api.json for building a browser-free scraper

Usage:
    python recon_bookmaker.py

Steps:
    1. Browser opens to be.bookmaker.eu
    2. Cloudflare challenge resolves automatically
    3. You can log in if needed (odds may be visible without login)
    4. Press ENTER — script navigates to sports pages and captures API calls
    5. Results saved to recon_bookmaker_api.json
"""

from playwright.sync_api import sync_playwright
import os
import json

# be.bookmaker.eu is the Angular SPA where odds are loaded via API
SITE_URL = "https://be.bookmaker.eu/en/sports/"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".bookmaker_profile")
OUT_DIR = os.path.dirname(__file__)

# Pages to visit — each should trigger API calls for odds data
ODDS_PAGES = [
    ("CBB", "https://be.bookmaker.eu/en/sports/basketball/college-basketball/"),
    ("NBA", "https://be.bookmaker.eu/en/sports/basketball/nba/"),
    ("NHL", "https://be.bookmaker.eu/en/sports/ice-hockey/nhl/"),
]


def run_recon():
    captured_requests = []
    seen_urls = set()  # Deduplicate requests across pages

    def handle_request(request):
        """Capture all API/XHR/fetch requests."""
        if request.resource_type in ("xhr", "fetch"):
            url = request.url
            # Skip analytics, tracking, chat
            skip = ("google", "hotjar", "liveperson", "cloudflare", "sentry",
                    "segment", "amplitude", "datadog", "newrelic")
            if any(s in url.lower() for s in skip):
                return

            dedup_key = f"{request.method}:{url}"
            if dedup_key in seen_urls:
                return
            seen_urls.add(dedup_key)

            entry = {
                "url": url,
                "method": request.method,
                "headers": dict(request.headers),
                "post_data": request.post_data,
                "resource_type": request.resource_type,
            }
            captured_requests.append(entry)
            print(f"  >> {request.method} {url[:140]}")

    def handle_response(response):
        """Capture response details for API calls."""
        if response.request.resource_type in ("xhr", "fetch"):
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

        # Step 1: Navigate to Bookmaker
        print("Navigating to Bookmaker.eu...")
        page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)

        # Step 2: Wait for Cloudflare + optional login
        print("\n" + "=" * 60)
        print("1. Wait for Cloudflare check to pass")
        print("2. Log in if needed (odds may work without login)")
        print("3. Once the sports page loads, press ENTER here")
        print("=" * 60)
        input()

        # Step 3: Attach network listeners
        print("\nAttaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 4: Visit each odds page to capture API calls
        for label, url in ODDS_PAGES:
            print(f"\n--- Navigating to {label}: {url} ---")
            page.goto(url, wait_until="domcontentloaded", timeout=60000)
            # Wait for Angular to bootstrap and API calls to fire
            page.wait_for_timeout(8000)
            print(f"    Captured {len(captured_requests)} total requests so far")

        # Step 5: Let user explore more pages if desired
        print(f"\n{'=' * 60}")
        print(f"Captured {len(captured_requests)} API requests so far.")
        print("You can now navigate to other pages (1H lines, props, etc.)")
        print("to capture more endpoints. Press ENTER when done.")
        print(f"{'=' * 60}")
        input()

        page.wait_for_timeout(3000)

        # Step 6: Save captured requests
        out_path = os.path.join(OUT_DIR, "recon_bookmaker_api.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(captured_requests, f, indent=2, default=str)
        print(f"\nSaved {len(captured_requests)} API requests to {out_path}")

        # Print summary
        print(f"\n{'=' * 60}")
        print("API REQUEST SUMMARY")
        print(f"{'=' * 60}")
        for i, req in enumerate(captured_requests, 1):
            status = req.get("status", "?")
            size = req.get("response_body_size", "?")
            print(f"  {i}. [{status}] {req['method']} {req['url'][:120]}")
            if req.get("post_data"):
                print(f"     POST: {req['post_data'][:200]}")
            if size != "?":
                print(f"     Response: {size:,} chars")

        # Highlight likely odds endpoints (large JSON responses)
        print(f"\n{'=' * 60}")
        print("LIKELY ODDS ENDPOINTS (large JSON responses)")
        print(f"{'=' * 60}")
        odds_candidates = [
            r for r in captured_requests
            if r.get("response_body_size", 0) > 500
            and "json" in r.get("response_headers", {}).get("content-type", "")
        ]
        for i, req in enumerate(odds_candidates, 1):
            print(f"\n  {i}. [{req.get('status')}] {req['method']} {req['url']}")
            print(f"     Size: {req.get('response_body_size', 0):,} chars")
            preview = req.get("response_body_preview", "")[:300]
            print(f"     Preview: {preview}")

        # Save cookies
        cookies = context.cookies()
        cookie_path = os.path.join(OUT_DIR, "recon_bookmaker_cookies.json")
        with open(cookie_path, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)
        print(f"\nSaved {len(cookies)} cookies to {cookie_path}")

        print("\nPress ENTER to close browser.")
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
