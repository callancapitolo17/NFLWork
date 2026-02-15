#!/usr/bin/env python3
"""Test building Rams +2.5 + Under 45.5 parlay on Pikkit."""

from scraper_pikkit import load_session, is_logged_in, decimal_to_american
from playwright.sync_api import sync_playwright
import re

RAMS_GAME_URL = "https://app.pikkit.com/event/696da1508ab2c5f0c9e021ca"

def test_rams_parlay():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        session = load_session()
        context = browser.new_context(
            viewport={"width": 1920, "height": 1080},
            storage_state=session
        )
        page = context.new_page()

        # Capture API
        api_data = [None]
        def capture(response):
            if "prod-website.pikkit.app/betslip" in response.url:
                try:
                    api_data[0] = response.json()
                except:
                    pass
        page.on("response", capture)

        print("Navigating to Rams vs Seahawks...")
        page.goto(RAMS_GAME_URL, wait_until="networkidle")
        page.wait_for_timeout(5000)

        print("Clicking Odds tab...")
        # Click Odds tab via JS
        page.evaluate('''() => {
            const els = document.querySelectorAll("button, div, span, a");
            for (const el of els) {
                const text = (el.innerText || "").trim();
                if (text === "Odds") {
                    el.click();
                    return true;
                }
            }
            return false;
        }''')

        # Wait for content to load (check for "Spread...Picks" text)
        print("Waiting for Odds content...")
        for i in range(30):
            page.wait_for_timeout(1000)
            loaded = page.evaluate('''() => {
                const text = document.body.innerText;
                return text.includes("Spread") && text.includes("Picks") && !text.includes("Loading");
            }''')
            if loaded:
                print(f"  Loaded after {i+1}s")
                break

        page.screenshot(path="debug_odds_content.png")

        # Expand Spread section
        print("Expanding Spread...")
        page.evaluate('''() => {
            const els = document.querySelectorAll("*");
            for (const el of els) {
                const text = (el.innerText || "").trim();
                if (text.startsWith("Spread") && text.includes("Picks") && text.length < 30) {
                    el.click();
                    return true;
                }
            }
            return false;
        }''')
        page.wait_for_timeout(3000)

        # Check available spreads
        spreads = page.evaluate(r'''() => {
            const found = [];
            const els = document.querySelectorAll("*");
            for (const el of els) {
                const text = (el.innerText || "").trim();
                if (/^[+-]\d+\.?\d*$/.test(text)) {
                    const rect = el.getBoundingClientRect();
                    if (rect.y > 300 && rect.height > 10 && rect.width > 10) {
                        found.push(text);
                    }
                }
            }
            return [...new Set(found)];
        }''')
        print(f"Available spreads: {spreads}")

        # Click +2.5 if available
        has_25 = "+2.5" in spreads
        if has_25:
            page.evaluate('''() => {
                const els = document.querySelectorAll("*");
                for (const el of els) {
                    const text = (el.innerText || "").trim();
                    if (text === "+2.5") {
                        const rect = el.getBoundingClientRect();
                        if (rect.y > 300 && rect.height > 10) {
                            const parent = el.parentElement;
                            if (parent) parent.click();
                            else el.click();
                            return true;
                        }
                    }
                }
                return false;
            }''')
            print("Clicked +2.5")
            page.wait_for_timeout(2000)
        else:
            print("+2.5 NOT available - this explains odds mismatch!")

        # Expand Total section
        print("Expanding Total...")
        page.evaluate('''() => {
            const els = document.querySelectorAll("*");
            for (const el of els) {
                const text = (el.innerText || "").trim();
                if (text.startsWith("Total") && !text.startsWith("Team") && text.includes("Picks") && text.length < 30) {
                    el.click();
                    return true;
                }
            }
            return false;
        }''')
        page.wait_for_timeout(3000)

        # Check available totals
        totals = page.evaluate(r'''() => {
            const found = [];
            const els = document.querySelectorAll("*");
            for (const el of els) {
                const text = (el.innerText || "").trim();
                if (/^\d{2}\.5$/.test(text)) {
                    const rect = el.getBoundingClientRect();
                    if (rect.y > 300 && rect.height > 10) {
                        found.push(text);
                    }
                }
            }
            return [...new Set(found)].sort();
        }''')
        print(f"Available totals: {totals}")

        # Find and click Under 45.5
        has_455 = "45.5" in totals
        if has_455:
            # First scroll to 45.5
            page.evaluate('''() => {
                const els = document.querySelectorAll("*");
                for (const el of els) {
                    if ((el.innerText || "").trim() === "45.5") {
                        el.scrollIntoView({block: "center"});
                        return true;
                    }
                }
                return false;
            }''')
            page.wait_for_timeout(1000)

            # Now find and click the Under cell (right side, +103 in screenshot)
            clicked = page.evaluate(r'''() => {
                const els = document.querySelectorAll("*");

                // Find 45.5 row Y position (need to re-query after scroll)
                let rowY = null;
                let totalEl = null;
                for (const el of els) {
                    if ((el.innerText || "").trim() === "45.5") {
                        const rect = el.getBoundingClientRect();
                        if (rect.height > 10 && rect.width > 20 && rect.y > 100 && rect.y < 800) {
                            rowY = rect.y;
                            totalEl = el;
                            break;
                        }
                    }
                }
                if (!rowY) return {error: "45.5 row not found after scroll"};

                // Find clickable elements ONLY on the 45.5 row (tight Y tolerance)
                const cells = [];
                for (const el of els) {
                    const text = (el.innerText || "").trim();
                    // Match odds like -104, +103, -110, etc.
                    if (/^[+-]\d{2,4}$/.test(text)) {
                        const rect = el.getBoundingClientRect();
                        // Very tight Y tolerance (within 10px) to get ONLY 45.5 row
                        if (Math.abs(rect.y - rowY) < 10 && rect.x > 300) {
                            cells.push({el, text, x: rect.x, y: rect.y});
                        }
                    }
                }

                if (cells.length === 0) {
                    return {error: "No odds cells found on 45.5 row", rowY: rowY};
                }

                // Sort by x position (left to right): Over is left, Under is right
                cells.sort((a, b) => a.x - b.x);

                // Under should be the RIGHTMOST cell (highest x value)
                // Over is left (-104), Under is right (+103)
                const underCell = cells[cells.length - 1];  // Last cell = rightmost = Under

                // Click the cell or its parent
                const clickTarget = underCell.el.closest("button") ||
                                   underCell.el.closest("div[role='button']") ||
                                   underCell.el.parentElement ||
                                   underCell.el;
                clickTarget.click();

                return {
                    clicked: true,
                    odds: underCell.text,
                    cellCount: cells.length,
                    allCells: cells.map(c => c.text)
                };
            }''')
            print(f"Under 45.5: {clicked}")
            page.wait_for_timeout(2000)
        else:
            print("45.5 NOT available!")

        # Click Parlay tab
        page.evaluate('''() => {
            const els = document.querySelectorAll("*");
            for (const el of els) {
                if ((el.innerText || "").trim() === "Parlay") {
                    el.click();
                    return true;
                }
            }
            return false;
        }''')
        page.wait_for_timeout(3000)

        page.screenshot(path="debug_rams_final.png")

        # Check results
        betslip_count = page.evaluate(r'''() => {
            const match = document.body.innerText.match(/(\d+)\s*Place Bets/);
            return match ? parseInt(match[1]) : 0;
        }''')
        print(f"\nBetslip count: {betslip_count}")

        best_odds = page.evaluate(r'''() => {
            const match = document.body.innerText.match(/(\d+)-Leg.*?([+-]\d{3})/);
            return match ? match[2] : null;
        }''')
        print(f"Best odds displayed: {best_odds}")

        # Print API data
        if api_data[0] and isinstance(api_data[0], dict):
            parlay = api_data[0].get("parlay") or {}
            options = parlay.get("options") or []
            if options:
                print(f"\nAPI ({len(options)} books):")
                for opt in options:
                    book = opt.get("institution_name", "?")
                    dec = float(opt.get("odds", 0))
                    amer = decimal_to_american(dec)
                    print(f"  {book}: {amer:+d}")
        else:
            print("\nNo parlay API data captured")

        browser.close()


if __name__ == "__main__":
    test_rams_parlay()
