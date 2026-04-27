"""
Wagerzon MLB Specials Recon Script

Sibling of recon_wagerzon.py. Launches Chromium with the existing
.wagerzon_profile, navigates to the MLB specials page (league 4899),
captures every network request + response, and dumps to recon_specials.json.

Usage:
    python3 wagerzon_odds/recon_specials.py

Prerequisites:
    - Playwright installed: pip install playwright && playwright install chromium
    - .wagerzon_profile directory populated (log in via recon_wagerzon.py once
      if starting from scratch).

Workflow:
    1. Script launches persistent-context Chromium.
    2. Script navigates to the candidate specials URL (see CANDIDATE_URLS).
    3. You manually verify the specials page rendered correctly — if the URL
       is wrong, close the browser, fix CANDIDATE_URLS, re-run.
    4. Press ENTER in the terminal when the page is fully loaded with all
       triple-play / grand-slam rows visible.
    5. Captured traffic is saved to wagerzon_odds/recon_specials.json.
    6. Inspect the JSON to determine: JSON API endpoint? HTML-only? URL
       structure? Section-header format?

The follow-up plan uses this recon output to design the production scraper.
"""

from playwright.sync_api import sync_playwright
import json
import os

WAGERZON_URL = "https://backend.wagerzon.com"

# Candidate URLs for MLB specials. If the first fails, adjust based on what
# you see in the browser's devtools Network tab when you manually navigate
# to the specials page.
CANDIDATE_URLS = [
    f"{WAGERZON_URL}/wager/NewSchedule.aspx?WT=0&lg=4899",
]

PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")
OUT_FILE = os.path.join(os.path.dirname(__file__), "recon_specials.json")


def run_recon():
    captured = []

    def handle_request(request):
        if "wagerzon.com" not in request.url:
            return
        captured.append({
            "phase": "request",
            "url": request.url,
            "method": request.method,
            "headers": dict(request.headers),
            "post_data": request.post_data,
            "resource_type": request.resource_type,
        })
        rtype = request.resource_type
        tag = {"document": "[PAGE]", "xhr": "[API ]", "fetch": "[API ]"}.get(rtype, "[____]")
        print(f"  {tag} {request.method} {request.url[:130]}")

    def handle_response(response):
        if "wagerzon.com" not in response.url:
            return
        try:
            content_type = response.headers.get("content-type", "")
            body = None
            if "json" in content_type:
                body = response.text()[:20000]
            elif "html" in content_type and response.request.resource_type == "document":
                body = response.text()[:40000]
            captured.append({
                "phase": "response",
                "url": response.url,
                "status": response.status,
                "content_type": content_type,
                "body_preview": body,
                "body_length": len(response.text()) if body is not None else None,
            })
        except Exception as e:
            captured.append({
                "phase": "response",
                "url": response.url,
                "status": response.status,
                "error": str(e),
            })

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",  # use installed Chrome — bundled Chromium crashes on this machine
            headless=False,
            viewport={"width": 1440, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        # Reuse the first existing tab so login traffic is captured by our handlers
        page = context.pages[0] if context.pages else context.new_page()
        page.on("request", handle_request)
        page.on("response", handle_response)

        # Step 1: navigate to candidate specials URL. If redirected to login,
        # the user logs in in the SAME tab; handlers capture the auth flow.
        for url in CANDIDATE_URLS:
            print(f"\n=== Navigating to candidate URL: {url} ===")
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as e:
                print(f"Navigation failed: {e}")

        print("\n" + "=" * 60)
        print("1. If redirected to login, log in manually IN THIS BROWSER TAB")
        print("2. After login, navigate to the MLB specials page")
        print("3. Once specials are fully visible, press ENTER here")
        print("=" * 60)
        input()

        # Step 2: reload the specials page to get a clean post-auth capture
        print(f"\nReloading {CANDIDATE_URLS[0]} to capture a clean specials request...")
        try:
            page.goto(CANDIDATE_URLS[0], wait_until="domcontentloaded", timeout=60000)
            page.wait_for_timeout(5000)
        except Exception as e:
            print(f"Reload failed: {e}")
        print(f"Total events captured so far: {len(captured)}")

        input("\nPress ENTER once the page is fully loaded a second time (or to skip)...")

        out = {
            "candidate_urls": CANDIDATE_URLS,
            "captured_events": captured,
        }
        with open(OUT_FILE, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nSaved {len(captured)} captured events to {OUT_FILE}")

        input("Press ENTER to close the browser...")
        context.close()


if __name__ == "__main__":
    run_recon()
