#!/usr/bin/env python3
"""
Test Pikkit odds scraping - Broncos +4.5 + Under 43.5 parlay.
"""

import json
from pathlib import Path
from playwright.sync_api import sync_playwright

PIKKIT_URL = "https://app.pikkit.com"
SESSION_FILE = Path(__file__).parent / ".pikkit_session.json"


def decimal_to_american(decimal_odds: float) -> int:
    if decimal_odds >= 2.0:
        return int(round((decimal_odds - 1) * 100))
    else:
        return int(round(-100 / (decimal_odds - 1)))


def load_session():
    if SESSION_FILE.exists():
        try:
            with open(SESSION_FILE, 'r') as f:
                return json.load(f)
        except:
            return None
    return None


def test_broncos_parlay(visible: bool = True):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=not visible)

        saved_session = load_session()
        if not saved_session:
            print("No saved session")
            browser.close()
            return

        context = browser.new_context(
            viewport={'width': 1920, 'height': 1080},
            storage_state=saved_session
        )

        page = context.new_page()

        # Capture betslip API
        betslip_responses = []

        def capture_betslip(response):
            if 'prod-website.pikkit.app/betslip' in response.url:
                try:
                    data = response.json()
                    if data and isinstance(data, dict):
                        betslip_responses.append(data)
                        parlay = data.get('parlay', {})
                        if parlay and parlay.get('options'):
                            options = parlay['options']
                            print(f"\n=== CAPTURED PARLAY: {len(options)} books ===")
                            for opt in options:
                                book = opt.get('institution_name', 'Unknown')
                                decimal_odds = opt.get('odds', 0)
                                try:
                                    decimal_odds = float(decimal_odds)
                                    american = decimal_to_american(decimal_odds)
                                    print(f"  {book}: {decimal_odds:.3f} = {american:+d}")
                                except:
                                    pass
                except:
                    pass

        page.on('response', capture_betslip)

        # Navigate to the Broncos game directly
        print("Navigating to Broncos game...")
        page.goto(f"{PIKKIT_URL}/event/696d767b8ab2c5f0c9d744ae", wait_until='networkidle')
        page.wait_for_timeout(3000)

        current_url = page.url
        print(f"Current URL: {current_url}")

        page.screenshot(path="debug_game_page.png")
        print("Saved debug_game_page.png")

        # Click on spread cell for Broncos
        print("\n=== Building parlay ===")
        print("Clicking spread (+4)...")

        # Use Playwright locator to click on the element containing "+4" followed by odds
        # The structure shows: +4 (top line) -110 (bottom line)
        try:
            # Try to find and click the spread cell - look for the "-110" odds near "+4"
            spread_cell = page.locator('text=-110').first
            if spread_cell.count() > 0:
                spread_cell.click()
                page.wait_for_timeout(1500)
                print("Clicked spread via -110 odds")
            else:
                print("Could not find -110 odds")
        except Exception as e:
            print(f"Spread click error: {e}")

        page.screenshot(path="debug_after_spread.png")

        # Click on under total
        print("Clicking under (u43.5)...")
        try:
            # Look for u43.5 text
            under_cell = page.locator('text=u43.5').first
            if under_cell.count() > 0:
                under_cell.click()
                page.wait_for_timeout(1500)
                print("Clicked under via u43.5 text")
            else:
                # Try clicking on "-115" for the under
                under_odds = page.locator('text=-115').first
                if under_odds.count() > 0:
                    under_odds.click()
                    page.wait_for_timeout(1500)
                    print("Clicked under via -115 odds")
                else:
                    print("Could not find u43.5 or -115")
        except Exception as e:
            print(f"Under click error: {e}")

        page.screenshot(path="debug_after_under.png")

        # Check betslip state
        betslip_count = page.evaluate('''() => {
            const text = document.body.innerText;
            const match = text.match(/(\d+)\s*Place Bets/);
            return match ? parseInt(match[1]) : 0;
        }''')
        print(f"Betslip count: {betslip_count}")

        # Click Parlay tab if we have 2 selections
        if betslip_count >= 2:
            print("\nClicking Parlay tab...")
            try:
                page.get_by_text("Parlay", exact=True).first.click()
                page.wait_for_timeout(2000)
                print("Clicked Parlay tab")
            except:
                print("Could not find Parlay tab")

        # Wait for parlay info
        parlay_info = None
        for i in range(10):
            page.wait_for_timeout(500)
            parlay_info = page.evaluate(r'''() => {
                const text = document.body.innerText;
                const match = text.match(/(\d+)-Leg Same Game Parlay\s+([+-]\d+)/);
                if (match) {
                    return {legs: parseInt(match[1]), best_odds: match[2]};
                }
                return null;
            }''')
            if parlay_info:
                print(f"Parlay detected: {parlay_info}")
                break

        page.screenshot(path="debug_final_parlay.png")
        print("\nSaved debug_final_parlay.png")

        # Results
        print(f"\n{'='*60}")
        print("FINAL RESULTS")
        print(f"{'='*60}")

        if parlay_info:
            print(f"Parlay: {parlay_info['legs']}-Leg at {parlay_info['best_odds']}")

        if betslip_responses:
            print(f"\nCaptured {len(betslip_responses)} betslip responses")

            # Find the parlay response
            for resp in reversed(betslip_responses):
                parlay = resp.get('parlay', {})
                if parlay and parlay.get('options'):
                    options = parlay['options']
                    print(f"\nParlay from {len(options)} books:")

                    best_american = 0
                    best_book = ""
                    all_books = []

                    for opt in options:
                        book = opt.get('institution_name', 'Unknown')
                        decimal_odds = opt.get('odds', 0)
                        try:
                            decimal_odds = float(decimal_odds)
                            american = decimal_to_american(decimal_odds)
                            all_books.append(f"{book}: {american:+d}")
                            if american > best_american:
                                best_american = american
                                best_book = book
                        except:
                            pass

                    print(f"\nAll books: {all_books}")
                    print(f"\nBest odds: {best_american:+d} at {best_book}")
                    break
        else:
            print("\nNo parlay API responses captured!")

        browser.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--headless', action='store_true')
    args = parser.parse_args()

    test_broncos_parlay(visible=not args.headless)
