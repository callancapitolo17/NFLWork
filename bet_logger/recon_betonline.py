"""
BetOnline Recon Script
- Uses real Chrome (not Playwright Chromium) to bypass Cloudflare
- Persistent profile so cookies survive between runs
- Waits for manual login, then saves bet history HTML
"""

from playwright.sync_api import sync_playwright
import os

BET_HISTORY_URL = "https://www.betonline.ag/my-account/bet-history"
LOGIN_URL = "https://www.betonline.ag"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".betonline_profile")


def run_recon():
    with sync_playwright() as p:
        # Use real Chrome with a persistent profile — looks like a real user
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",          # Use installed Chrome, not Chromium
            headless=False,
            slow_mo=300,
            viewport={"width": 1400, "height": 900},
            args=[
                "--disable-blink-features=AutomationControlled",
            ],
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

        # Step 3: Navigate to bet history
        print("Navigating to bet history...")
        page.goto(BET_HISTORY_URL, wait_until="domcontentloaded", timeout=60000)

        # Give JS time to render
        page.wait_for_timeout(5000)

        # Step 4: Save the full page HTML
        html = page.content()
        out_path = os.path.join(os.path.dirname(__file__), "recon_betonline.html")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"\nSaved page HTML to {out_path} ({len(html):,} chars)")

        # Step 5: Keep browser open for DevTools inspection
        print("\n" + "=" * 60)
        print("Browser is still open. Use DevTools (Cmd+Option+I) to inspect.")
        print("Look for:")
        print("  - Login form field IDs/names")
        print("  - Bet history table/card selectors")
        print("  - Date filter controls")
        print("  - How parlays display legs")
        print("When done, press ENTER to close.")
        print("=" * 60)
        input()

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
