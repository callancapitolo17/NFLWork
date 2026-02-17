"""
BetOnline Recon Script
- Uses real Chrome (not Playwright Chromium) to bypass Cloudflare
- Persistent profile so cookies survive between runs
- Captures all network API requests on the bet history page
- Saves requests to recon_betonline_api.json for building a browser-free scraper
"""

from playwright.sync_api import sync_playwright
import os
import json

BET_HISTORY_URL = "https://www.betonline.ag/my-account/bet-history"
LOGIN_URL = "https://www.betonline.ag"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".betonline_profile")
OUT_DIR = os.path.dirname(__file__)


def run_recon():
    captured_requests = []

    def handle_request(request):
        """Capture all API/XHR requests."""
        if request.resource_type in ("xhr", "fetch"):
            entry = {
                "url": request.url,
                "method": request.method,
                "headers": dict(request.headers),
                "post_data": request.post_data,
                "resource_type": request.resource_type,
            }
            captured_requests.append(entry)
            print(f"  📡 {request.method} {request.url[:120]}")

    def handle_response(response):
        """Capture response details for API calls."""
        if response.request.resource_type in ("xhr", "fetch"):
            # Find matching request and add response info
            for entry in captured_requests:
                if entry["url"] == response.url and "status" not in entry:
                    entry["status"] = response.status
                    entry["response_headers"] = dict(response.headers)
                    # Try to capture response body for small JSON responses
                    try:
                        content_type = response.headers.get("content-type", "")
                        if "json" in content_type:
                            body = response.text()
                            # Only save first 5000 chars to keep file manageable
                            entry["response_body_preview"] = body[:5000]
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

        # Step 1: Go to BetOnline
        print("Navigating to BetOnline...")
        page.goto(LOGIN_URL, wait_until="domcontentloaded", timeout=60000)

        # Step 2: Wait for Cloudflare + manual login
        print("\n" + "=" * 60)
        print("1. Wait for Cloudflare check to pass (should be automatic)")
        print("2. Log in manually if needed")
        print("3. Once you see the main site, press ENTER here")
        print("=" * 60)
        input()

        # Step 3: Attach network listeners BEFORE navigating to bet history
        print("\nAttaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 4: Navigate to bet history
        print("Navigating to bet history...\n")
        page.goto(BET_HISTORY_URL, wait_until="domcontentloaded", timeout=60000)

        # Wait for page to fully render and API calls to complete
        print("\nWaiting for API calls to complete...")
        page.wait_for_timeout(10000)

        print(f"\n{'=' * 60}")
        print(f"Captured {len(captured_requests)} API requests")
        print(f"{'=' * 60}")

        # Step 5: Try changing date filter to capture that API call too
        print("\nNow try changing the date filter (e.g., click '30 days').")
        print("This will capture the filter API call.")
        print("Press ENTER when done.")
        input()

        # Wait for any new API calls
        page.wait_for_timeout(3000)

        # Step 6: Save captured requests
        out_path = os.path.join(OUT_DIR, "recon_betonline_api.json")
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
            print(f"  {i}. [{status}] {req['method']} {req['url'][:100]}")
            if req.get("post_data"):
                print(f"     POST: {req['post_data'][:100]}")
            if size != "?":
                print(f"     Response: {size} chars")

        # Also save cookies for potential use with requests library
        cookies = context.cookies()
        cookie_path = os.path.join(OUT_DIR, "recon_betonline_cookies.json")
        with open(cookie_path, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)
        print(f"\nSaved {len(cookies)} cookies to {cookie_path}")

        print("\nPress ENTER to close browser.")
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
