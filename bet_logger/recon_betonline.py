"""
BetOnline Recon Script — Automated Token Capture

Uses real Chrome with persistent profile to:
1. Pass Cloudflare challenge (automatic with saved cookies)
2. Log in via Keycloak (automatic with saved credentials)
3. Navigate to bet history to trigger API token generation
4. Capture krefresh token and cookies for the headless scraper

By default, runs fully automated — same Cloudflare/login patterns as
bet_placer/navigator_betonline.py. Use --interactive for manual control.
"""

import os
import sys
import json
import time
import argparse

from playwright.sync_api import sync_playwright
from dotenv import load_dotenv

load_dotenv()

BET_HISTORY_URL = "https://www.betonline.ag/my-account/bet-history"
LOGIN_URL = "https://www.betonline.ag"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".betonline_profile")
OUT_DIR = os.path.dirname(__file__)
COOKIES_FILE = os.path.join(OUT_DIR, "recon_betonline_cookies.json")

BETONLINE_USERNAME = os.getenv("BETONLINE_USERNAME")
BETONLINE_PASSWORD = os.getenv("BETONLINE_PASSWORD")


# ── Cloudflare + Auth (same patterns as navigator_betonline.py) ──


def wait_for_cloudflare(page, max_wait: int = 60) -> bool:
    """Wait for Cloudflare challenge to pass.

    With a persistent Chrome profile, Cloudflare usually passes instantly
    because the cf_clearance cookie is already saved.
    """
    for i in range(max_wait):
        title = page.title().lower()
        url = page.url.lower()

        if "just a moment" in title or "challenge" in title:
            if i == 0:
                print("  Waiting for Cloudflare...")
            time.sleep(1)
            continue

        if "betonline" in title or "betonline" in url:
            if i > 0:
                print(f"  Cloudflare passed ({i}s)")
            return True

        time.sleep(1)

    print("  Cloudflare may not have passed — continuing anyway")
    return False


def is_logged_in(page) -> bool:
    """Check if the user session is authenticated."""
    try:
        balance = page.locator('[data-testid="balance"], .balance, [class*="balance"]')
        if balance.count() > 0 and balance.first.is_visible(timeout=5000):
            return True
    except Exception:
        pass

    try:
        login_btn = page.locator('button:has-text("Log In"), a:has-text("Log In"), [data-testid="login"]')
        if login_btn.count() > 0 and login_btn.first.is_visible(timeout=3000):
            return False
    except Exception:
        pass

    try:
        deposit = page.locator('button:has-text("Deposit"), a:has-text("Deposit")')
        if deposit.count() > 0 and deposit.first.is_visible(timeout=3000):
            return True
    except Exception:
        pass

    # Default to logged in (persistent profile usually works)
    return True


def do_login(page):
    """Log in via Keycloak. Same approach as navigator_betonline._login."""
    if not BETONLINE_USERNAME or not BETONLINE_PASSWORD:
        raise RuntimeError(
            "BETONLINE_USERNAME/PASSWORD not set in .env. "
            "Cannot auto-login. Run with --interactive for manual login."
        )

    print("  Logging in...")
    try:
        # Check if we're already on the Keycloak login form.
        # The persistent profile sometimes redirects straight to Keycloak
        # instead of the main BetOnline site. In that case we skip the
        # "Log In" button click and fill credentials directly.
        on_keycloak = (
            "auth" in page.url.lower() or
            "login" in page.url.lower() or
            page.locator('#kc-login, input[name="username"]').count() > 0
        )

        if not on_keycloak:
            # On the main site — click the Log In button to navigate to Keycloak
            login_btn = page.locator('[data-testid="login"], a:has-text("Log In")')
            if login_btn.count() > 0:
                login_btn.first.click()
                time.sleep(3)

        # Now we should be on the Keycloak form — fill credentials
        user_field = page.locator('input[name="username"]')
        if user_field.count() == 0:
            raise RuntimeError(f"No username field found. URL: {page.url}")

        user_field.first.click()
        page.keyboard.type(BETONLINE_USERNAME, delay=50)

        pass_field = page.locator('input[name="password"]').first
        pass_field.click()
        page.keyboard.type(BETONLINE_PASSWORD, delay=50)

        # Click submit — the #kc-login button enables after fields are filled
        time.sleep(0.5)
        submit = page.locator('#kc-login, input[type="submit"], button[type="submit"]')
        if submit.count() > 0:
            submit.first.click()
        else:
            page.keyboard.press("Enter")

        time.sleep(5)
        print(f"  Login submitted. URL: {page.url}")
        return

    except Exception as e:
        raise RuntimeError(f"Auto-login failed: {e}")


# ── Main recon flow ──────────────────────────────────────────────


def run_recon(interactive: bool = False):
    """Run recon to capture BetOnline auth tokens.

    Automated mode (default): Cloudflare + Keycloak login handled
    automatically using persistent Chrome profile + .env credentials.

    Interactive mode (--interactive): pauses for manual steps at each stage.
    """
    captured_requests = []

    def handle_request(request):
        """Capture API/XHR requests."""
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
            for entry in captured_requests:
                if entry["url"] == response.url and "status" not in entry:
                    entry["status"] = response.status
                    entry["response_headers"] = dict(response.headers)
                    try:
                        content_type = response.headers.get("content-type", "")
                        if "json" in content_type:
                            body = response.text()
                            entry["response_body_preview"] = body[:5000]
                            entry["response_body_size"] = len(body)
                    except Exception:
                        pass
                    break

    print("=" * 60)
    print("BETONLINE RECON — Token Capture")
    print("=" * 60)

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()

        # Step 1: Navigate to BetOnline
        print("\nNavigating to BetOnline...")
        page.goto(LOGIN_URL, wait_until="domcontentloaded", timeout=60000)

        # Step 2: Handle Cloudflare
        wait_for_cloudflare(page)

        if interactive:
            print("\nCloudflare/login: handle manually if needed, then press ENTER.")
            input()

        # Step 3: Login if needed
        if not is_logged_in(page):
            if interactive:
                print("Not logged in. Please log in manually, then press ENTER.")
                input()
            else:
                do_login(page)
                wait_for_cloudflare(page)
        else:
            print("  Already authenticated")

        # Step 4: Attach network listeners, then navigate to bet history
        print("\nAttaching network listeners...")
        page.on("request", handle_request)
        page.on("response", handle_response)

        print("Navigating to bet history...")
        page.goto(BET_HISTORY_URL, wait_until="domcontentloaded", timeout=60000)

        # Wait for API calls to complete (token exchange + data load)
        print("Waiting for API calls to complete...")
        page.wait_for_timeout(10000)

        print(f"\nCaptured {len(captured_requests)} API requests")

        if interactive:
            print("\nOptional: change date filters to capture more API calls.")
            print("Press ENTER when done.")
            input()
            page.wait_for_timeout(3000)

        # Step 5: Save captured requests
        out_path = os.path.join(OUT_DIR, "recon_betonline_api.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(captured_requests, f, indent=2, default=str)
        print(f"Saved {len(captured_requests)} API requests to {out_path}")

        # Step 6: Save cookies (includes krefresh token)
        cookies = context.cookies()
        with open(COOKIES_FILE, "w", encoding="utf-8") as f:
            json.dump(cookies, f, indent=2)

        # Verify we got the krefresh token
        krefresh = None
        for c in cookies:
            if c['name'] == 'krefresh':
                krefresh = c['value']
                break

        print(f"\n{'=' * 60}")
        if krefresh:
            print(f"Saved {len(cookies)} cookies to {COOKIES_FILE}")
            print(f"krefresh token captured ({len(krefresh)} chars)")
            print("The headless scraper should now work.")
        else:
            print(f"WARNING: Saved {len(cookies)} cookies but krefresh NOT found.")
            print("Try running with --interactive to debug.")
        print(f"{'=' * 60}")

        if interactive:
            print("\nPress ENTER to close browser.")
            input()

        context.close()
        print("Done.")

    return bool(krefresh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Capture BetOnline auth tokens')
    parser.add_argument('--interactive', action='store_true',
                        help='Pause for manual steps (for debugging)')
    args = parser.parse_args()

    success = run_recon(interactive=args.interactive)
    sys.exit(0 if success else 1)
