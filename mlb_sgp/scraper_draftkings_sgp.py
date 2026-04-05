#!/usr/bin/env python3
"""
DraftKings MLB SGP Scraper

Fetches Same Game Parlay (SGP) odds from DraftKings for MLB spread+total combos.

DraftKings' SGP pricing is behind Akamai bot protection. The only way to get
correlation-adjusted parlay prices is to connect to a real Chrome browser via
Chrome DevTools Protocol (CDP) and programmatically click spread + total legs
in the SGP builder, intercepting the calculateBets API response.

Setup:
    1. Start Chrome with remote debugging:
       /Applications/Google\\ Chrome.app/Contents/MacOS/Google\\ Chrome --remote-debugging-port=9222

    2. Navigate to DraftKings MLB in that Chrome window:
       https://sportsbook.draftkings.com/leagues/baseball/mlb

    3. Run this scraper:
       python scraper_draftkings_sgp.py --verbose

How it works:
    - Public REST API fetches game list and market IDs (no browser needed)
    - CDP connects to your Chrome, clicks into each game's SGP builder
    - Clicks the Run Line + Total odds cells
    - Intercepts the calculateBets API response with trueOdds (correlation-adjusted)
    - Stores results in mlb_sgp_odds table in mlb.duckdb
"""

import argparse
import re
import time
import requests
from playwright.sync_api import sync_playwright

from db import ensure_table, upsert_sgp_odds

# ---------------------------------------------------------------------------
# DK public REST API (no browser needed)
# ---------------------------------------------------------------------------

DK_LEAGUE_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/controldata/"
    "league/leagueSubcategory/v1/markets"
)
DK_MLB_LEAGUE_ID = "84240"
DK_BASE_URL = "https://sportsbook.draftkings.com"

REST_HEADERS = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "application/json",
    "Referer": "https://sportsbook.draftkings.com/",
}


def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    else:
        return int(round(-100 / (dec - 1)))


def fetch_dk_events() -> list[dict]:
    """Fetch today's MLB events from DK public API."""
    params = {
        "isBatchable": "false",
        "templateVars": DK_MLB_LEAGUE_ID,
        "eventsQuery": (
            f"$filter=leagueId eq '{DK_MLB_LEAGUE_ID}' "
            f"AND clientMetadata/Subcategories/any(s: s/Id eq '4519')"
        ),
        "marketsQuery": (
            "$filter=clientMetadata/subCategoryId eq '4519' "
            "AND tags/all(t: t ne 'SportcastBetBuilder')"
        ),
        "include": "Events",
        "entity": "events",
    }
    resp = requests.get(DK_LEAGUE_URL, params=params, headers=REST_HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    events = []
    for evt in data.get("events", []):
        participants = evt.get("participants", [])
        home = next((p for p in participants if p.get("venueRole") == "Home"), {})
        away = next((p for p in participants if p.get("venueRole") == "Away"), {})
        events.append({
            "dk_event_id": evt["id"],
            "name": evt.get("name", ""),
            "seo_id": evt.get("seoIdentifier", ""),
            "home_team": home.get("name", ""),
            "away_team": away.get("name", ""),
            "start_time": evt.get("startEventDate", ""),
        })
    return events


# ---------------------------------------------------------------------------
# SGP scraping via CDP click-and-capture
# ---------------------------------------------------------------------------

def clear_betslip(page):
    """Clear all selections from the DK betslip."""
    try:
        clear = page.locator('text="Clear All"')
        if clear.count() > 0 and clear.first.is_visible():
            clear.first.click(force=True)
            page.wait_for_timeout(1000)
            return
    except Exception:
        pass
    # Fallback: click X buttons
    page.evaluate(r'''() => {
        document.querySelectorAll('[aria-label*="emove"], [data-testid*="remove"]')
            .forEach(b => b.click());
    }''')
    page.wait_for_timeout(500)


def find_odds_cell(page, label: str, line_value: str):
    """
    Find an odds cell in the DK SGP builder by matching its text content.

    DK renders cells with label + line + odds in separate spans, but the
    parent element's textContent concatenates them (e.g., "O7.5-110").

    Returns {text, x, y} or None.
    """
    return page.evaluate(f'''() => {{
        const label = '{label}';
        const lineVal = '{line_value}';
        const target = label + lineVal;
        const all = document.querySelectorAll('td, div, button, span');
        for (const el of all) {{
            const raw = el.textContent || '';
            // Remove all whitespace for matching
            const t = raw.replace(/\\s+/g, '');
            if (t.includes(target) && /\\d{{3}}/.test(t) && t.length < 25) {{
                const r = el.getBoundingClientRect();
                if (r.width > 20 && r.width < 250 && r.height > 10 && r.height < 100 && r.y > 0) {{
                    return {{text: raw.replace(/\\s+/g, ' ').trim(), x: r.x + r.width/2, y: r.y + r.height/2}};
                }}
            }}
        }}
        return null;
    }}''')


def navigate_to_game_sgp(page, dk_event_id: str, verbose: bool = False) -> bool:
    """
    Navigate to a game's SGP builder.

    Tries SPA link click first (from MLB schedule), falls back to
    direct URL + reload.
    """
    clear_betslip(page)

    # Try clicking a game link on the current page (SPA navigation)
    clicked = page.evaluate(f'''() => {{
        const links = document.querySelectorAll('a');
        for (const a of links) {{
            if (a.href && a.href.includes('{dk_event_id}')) {{
                a.click();
                return true;
            }}
        }}
        return false;
    }}''')

    if not clicked:
        if verbose:
            print("    No game link found, using direct URL")
        # Direct navigation — may not render in SPA
        try:
            page.goto(
                f"{DK_BASE_URL}/leagues/baseball/mlb",
                wait_until="commit", timeout=15000,
            )
        except Exception:
            pass
        page.wait_for_timeout(3000)
        # Retry link click after schedule loads
        page.evaluate(f'''() => {{
            const links = document.querySelectorAll('a');
            for (const a of links) {{
                if (a.href && a.href.includes('{dk_event_id}')) {{
                    a.click();
                    return true;
                }}
            }}
            return false;
        }}''')

    # Wait for game content
    for i in range(20):
        page.wait_for_timeout(1000)
        has_content = page.evaluate(
            "() => document.body.innerText.includes('Run Line') || "
            "document.body.innerText.includes('SGP')"
        )
        if has_content:
            if verbose:
                print(f"    Game loaded in {i+1}s")
            break
    else:
        return False

    # Click SGP tab
    for sel in ['text="SGP"', 'a:has-text("SGP")']:
        try:
            loc = page.locator(sel).first
            if loc.count() > 0:
                loc.click(force=True, timeout=3000)
                page.wait_for_timeout(2000)
                break
        except Exception:
            continue

    # Verify Run Line + Total visible
    for _ in range(10):
        if page.evaluate(
            "() => document.body.innerText.includes('Run Line') && "
            "document.body.innerText.includes('Total')"
        ):
            return True
        page.wait_for_timeout(500)

    return page.evaluate("() => document.body.innerText.includes('Run Line')")


def read_game_lines(page) -> dict:
    """Read the available spread and total lines from the SGP builder."""
    text = page.evaluate("() => document.body.innerText")

    # Find spreads (e.g., -1.5, +1.5)
    spreads = list(set(re.findall(r'[+-]1\.5', text)))

    # Find totals — DK renders as "O 7.5" or "U 8" (with space between O/U and number)
    total_matches = re.findall(r'[OU]\s+(\d+\.?\d*)', text)
    total_value = total_matches[0] if total_matches else None

    return {"spreads": spreads, "total": total_value}


def scrape_combo(page, spread: str, ou: str, total: str, combo_name: str,
                 verbose: bool = False) -> dict | None:
    """
    Click spread + total in the SGP builder and capture the calculateBets response.

    Returns {combo_name, trueOdds, displayOdds, selections} or None.
    """
    # Set up response capture
    calc_responses = []

    def on_resp(response):
        if "calculateBets" in response.url or "sportsdata/v2/sgp" in response.url:
            try:
                calc_responses.append(response.json())
            except Exception:
                pass

    page.on("response", on_resp)

    # Clear previous
    clear_betslip(page)
    page.wait_for_timeout(800)

    # Click spread cell
    spread_cell = find_odds_cell(page, spread, "")
    if not spread_cell:
        if verbose:
            print(f"      Spread '{spread}' not found")
        page.remove_listener("response", on_resp)
        return None

    if verbose:
        print(f"      Spread: '{spread_cell['text']}' at ({spread_cell['x']:.0f}, {spread_cell['y']:.0f})")
    page.mouse.click(spread_cell['x'], spread_cell['y'])
    page.wait_for_timeout(2000)

    # Click total cell
    total_cell = find_odds_cell(page, ou, total)
    if not total_cell:
        if verbose:
            print(f"      Total '{ou}{total}' not found")
        page.remove_listener("response", on_resp)
        clear_betslip(page)
        return None

    if verbose:
        print(f"      Total: '{total_cell['text']}' at ({total_cell['x']:.0f}, {total_cell['y']:.0f})")
    page.mouse.click(total_cell['x'], total_cell['y'])

    # Wait for 2-leg SGP response
    for _ in range(10):
        page.wait_for_timeout(500)
        for resp_body in calc_responses:
            for bet in resp_body.get("bets", []):
                mapped = bet.get("selectionsMapped", [])
                if bet.get("trueOdds") and len(mapped) >= 2:
                    page.remove_listener("response", on_resp)
                    return {
                        "combo_name": combo_name,
                        "trueOdds": bet["trueOdds"],
                        "displayOdds": bet.get("displayOdds"),
                        "selections": [m.get("id") for m in mapped],
                    }

    page.remove_listener("response", on_resp)
    return None


def scrape_game(page, event: dict, verbose: bool = False) -> list[dict]:
    """Scrape all 4 SGP combos for a game."""
    results = []

    if not navigate_to_game_sgp(page, event["dk_event_id"], verbose):
        print("    Could not load SGP builder")
        return results

    lines = read_game_lines(page)
    if verbose:
        print(f"    Lines: spreads={lines['spreads']}, total={lines['total']}")

    if not lines["total"]:
        print("    No total line found")
        return results

    total = lines["total"]
    combos = [
        ("-1.5", "O", total, "Home Spread + Over"),
        ("-1.5", "U", total, "Home Spread + Under"),
        ("+1.5", "O", total, "Away Spread + Over"),
        ("+1.5", "U", total, "Away Spread + Under"),
    ]

    for spread, ou, tot, name in combos:
        if verbose:
            print(f"\n    {name}: {spread} + {ou}{tot}")

        result = scrape_combo(page, spread, ou, tot, name, verbose)
        if result:
            odds = result["trueOdds"]
            am = decimal_to_american(odds)
            print(f"    {name}: {odds:.4f} ({am:+d})")
            results.append(result)
        elif verbose:
            print(f"    {name}: not captured")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scrape_dk_sgp(verbose: bool = False):
    """Main: connect to Chrome via CDP, scrape SGP odds for all MLB games."""
    print("Fetching DraftKings MLB events...")
    dk_events = fetch_dk_events()
    print(f"  Found {len(dk_events)} events")

    if not dk_events:
        return []

    ensure_table()
    all_rows = []

    with sync_playwright() as p:
        try:
            browser = p.chromium.connect_over_cdp("http://localhost:9222")
        except Exception as e:
            print(f"\nCould not connect to Chrome.")
            print("Start Chrome with remote debugging first:")
            print('  /Applications/Google\\ Chrome.app/Contents/MacOS/Google\\ Chrome --remote-debugging-port=9222')
            print("Then navigate to: https://sportsbook.draftkings.com/leagues/baseball/mlb")
            return []

        context = browser.contexts[0]
        page = context.pages[0] if context.pages else context.new_page()

        # Ensure we start on the MLB schedule (needed for SPA game link clicks)
        if "leagues/baseball/mlb" not in page.url:
            print("  Navigate to DK MLB page in your Chrome first!")
            print("  URL: https://sportsbook.draftkings.com/leagues/baseball/mlb")
            browser.close()
            return []

        for event in dk_events:
            print(f"\n{'='*60}")
            print(f"  {event['name']}")
            print(f"  Event: {event['dk_event_id']}")
            print(f"{'='*60}")

            game_results = scrape_game(page, event, verbose)

            for gr in game_results:
                all_rows.append({
                    "game_id": event["dk_event_id"],
                    "combo": gr["combo_name"],
                    "period": "FG",
                    "bookmaker": "draftkings",
                    "sgp_decimal": round(gr["trueOdds"], 4),
                    "sgp_american": decimal_to_american(gr["trueOdds"]),
                    "source": "draftkings_direct",
                })

            # Navigate back to MLB schedule for next game
            try:
                page.go_back()
                page.wait_for_timeout(2000)
            except Exception:
                try:
                    page.goto(f"{DK_BASE_URL}/leagues/baseball/mlb",
                              wait_until="commit", timeout=15000)
                except Exception:
                    pass
                page.wait_for_timeout(3000)

        browser.close()

    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} DK SGP odds to mlb_sgp_odds")
        print(f"{'='*60}")
    else:
        print("\nNo SGP odds collected.")

    return all_rows


def main():
    parser = argparse.ArgumentParser(description="DraftKings MLB SGP Scraper (CDP)")
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    args = parser.parse_args()

    print("=" * 60)
    print("  DRAFTKINGS MLB SGP SCRAPER")
    print("=" * 60)
    print("  Requires Chrome with --remote-debugging-port=9222")
    print("  Navigate to DK MLB page before running")
    print()

    scrape_dk_sgp(verbose=args.verbose)


if __name__ == "__main__":
    main()
