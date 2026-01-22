#!/usr/bin/env python3
"""
DraftKings Same-Game Parlay (SGP) Odds Scraper
Scrapes SGP odds for correlated spread+total parlays to compare with hoop88
"""

import os
import re
import duckdb
from datetime import datetime
from dotenv import load_dotenv

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

load_dotenv()

# DraftKings doesn't require login to view SGP odds
DK_BASE_URL = "https://sportsbook.draftkings.com"
DK_NFL_URL = f"{DK_BASE_URL}/leagues/football/nfl"


def scrape_dk_sgp_for_game(page, away_team: str, home_team: str, spread: float, total: float) -> dict:
    """
    Build SGP on DraftKings for a specific game and get the correlated odds.

    Args:
        page: Playwright page object
        away_team: Away team name
        home_team: Home team name
        spread: Underdog spread value (positive number)
        total: Game total

    Returns:
        Dict with SGP odds for both correlated combinations
    """
    result = {
        'underdog_under_odds': None,
        'favorite_over_odds': None,
        'error': None
    }

    try:
        # Find the game card - search by team names
        # DraftKings uses team names in their game cards
        game_found = False

        # Look for game card containing both team names
        game_cards = page.locator('[data-testid="event-card"], .sportsbook-event-accordion__wrapper')
        card_count = game_cards.count()
        print(f"  Found {card_count} game cards")

        for i in range(card_count):
            card = game_cards.nth(i)
            card_text = card.text_content()

            # Check if both teams are in this card (fuzzy match)
            away_match = away_team.lower() in card_text.lower()
            home_match = home_team.lower() in card_text.lower()

            if away_match or home_match:
                print(f"  Found potential match: {card_text[:100]}...")
                game_found = True

                # Click to expand the game and access SGP
                try:
                    card.click()
                    page.wait_for_timeout(2000)
                except:
                    pass

                break

        if not game_found:
            result['error'] = f"Game not found: {away_team} @ {home_team}"
            return result

        # Look for Same Game Parlay tab/section
        sgp_selectors = [
            'text=Same Game Parlay',
            'text=SGP',
            '[data-testid="sgp-tab"]',
            'button:has-text("Same Game")',
        ]

        for selector in sgp_selectors:
            try:
                elem = page.locator(selector)
                if elem.count() > 0 and elem.first.is_visible():
                    elem.first.click()
                    page.wait_for_timeout(2000)
                    print(f"  Clicked SGP tab using: {selector}")
                    break
            except:
                continue

        # Now we need to:
        # 1. Select the underdog spread
        # 2. Select the Under total
        # 3. Read the SGP odds
        # 4. Clear and repeat for favorite spread + Over

        # This is complex because DK's SGP builder is highly dynamic
        # For now, we'll try to capture the structure

        page.screenshot(path="/Users/callancapitolo/NFLWork/hoop88_correlation/debug_dk_sgp.png")

        # Get the SGP section HTML for analysis
        html = page.content()
        with open("/Users/callancapitolo/NFLWork/hoop88_correlation/debug_dk_sgp.html", "w") as f:
            f.write(html)

        print("  Saved DK SGP page for analysis")

    except Exception as e:
        result['error'] = str(e)

    return result


def scrape_draftkings_sgp(games: list, headless: bool = True) -> list:
    """
    Scrape DraftKings SGP odds for a list of games.

    Args:
        games: List of game dicts from hoop88 scraper
        headless: Run browser in headless mode

    Returns:
        List of games with DK SGP odds added
    """
    if not games:
        print("No games provided")
        return []

    results = []

    with sync_playwright() as p:
        # Use chromium with stealth-like settings
        browser = p.chromium.launch(
            headless=headless,
            args=['--disable-blink-features=AutomationControlled']
        )
        context = browser.new_context(
            viewport={'width': 1920, 'height': 1080},
            user_agent='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        )
        page = context.new_page()

        # Remove webdriver property that bots often have
        page.add_init_script("""
            Object.defineProperty(navigator, 'webdriver', {
                get: () => undefined
            });
        """)

        print(f"Navigating to DraftKings NFL...")
        try:
            page.goto(DK_NFL_URL, timeout=90000, wait_until='domcontentloaded')
            page.wait_for_timeout(5000)
        except Exception as e:
            print(f"Error navigating to DK: {e}")
            print("DraftKings may be blocking automated access or not available in your region.")
            browser.close()
            return games  # Return original games without DK data

        # Take initial screenshot
        page.screenshot(path="/Users/callancapitolo/NFLWork/hoop88_correlation/debug_dk_initial.png")
        print("Saved initial DK page screenshot")

        # Check if we need to handle any popups/modals
        try:
            close_modal = page.locator('[data-testid="modal-close"], .modal-close, button:has-text("Close")')
            if close_modal.count() > 0 and close_modal.first.is_visible():
                close_modal.first.click()
                page.wait_for_timeout(1000)
        except:
            pass

        # For each game, try to find it and get SGP odds
        for game in games:
            print(f"\nLooking for: {game['away_team']} @ {game['home_team']}")

            sgp_result = scrape_dk_sgp_for_game(
                page,
                game['away_team'],
                game['home_team'],
                game.get('underdog_spread', 0),
                game.get('total', 0)
            )

            game_with_dk = game.copy()
            game_with_dk['dk_underdog_under'] = sgp_result.get('underdog_under_odds')
            game_with_dk['dk_favorite_over'] = sgp_result.get('favorite_over_odds')
            game_with_dk['dk_error'] = sgp_result.get('error')

            results.append(game_with_dk)

        browser.close()

    return results


def save_dk_odds_to_db(games: list):
    """Save DraftKings SGP odds to DuckDB for tracking."""
    conn = duckdb.connect('/Users/callancapitolo/NFLWork/hoop88_correlation/dk_sgp.duckdb')

    conn.execute("""
        CREATE TABLE IF NOT EXISTS dk_sgp_odds (
            fetch_time TIMESTAMP,
            league VARCHAR,
            away_team VARCHAR,
            home_team VARCHAR,
            underdog VARCHAR,
            underdog_spread FLOAT,
            total FLOAT,
            hoop88_underdog_under INTEGER,
            hoop88_favorite_over INTEGER,
            dk_underdog_under INTEGER,
            dk_favorite_over INTEGER,
            edge_underdog_under FLOAT,
            edge_favorite_over FLOAT
        )
    """)

    fetch_time = datetime.now()

    for game in games:
        # Calculate edge if we have both prices
        edge_uu = None
        edge_fo = None

        if game.get('dk_underdog_under') and game.get('parlay_underdog_under'):
            # Edge = (hoop88_decimal - dk_decimal) / dk_decimal
            hoop88_dec = american_to_decimal(game['parlay_underdog_under'])
            dk_dec = american_to_decimal(game['dk_underdog_under'])
            edge_uu = (hoop88_dec - dk_dec) / dk_dec * 100

        if game.get('dk_favorite_over') and game.get('parlay_favorite_over'):
            hoop88_dec = american_to_decimal(game['parlay_favorite_over'])
            dk_dec = american_to_decimal(game['dk_favorite_over'])
            edge_fo = (hoop88_dec - dk_dec) / dk_dec * 100

        conn.execute("""
            INSERT INTO dk_sgp_odds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            fetch_time,
            game.get('league', ''),
            game['away_team'],
            game['home_team'],
            game.get('underdog', ''),
            game.get('underdog_spread', 0),
            game.get('total', 0),
            game.get('parlay_underdog_under'),
            game.get('parlay_favorite_over'),
            game.get('dk_underdog_under'),
            game.get('dk_favorite_over'),
            edge_uu,
            edge_fo
        ])

    conn.close()
    print(f"Saved {len(games)} games to DuckDB")


def american_to_decimal(american: int) -> float:
    """Convert American odds to decimal."""
    if american == 0:
        return 1.91
    if american > 0:
        return (american / 100) + 1
    else:
        return (100 / abs(american)) + 1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape DraftKings SGP odds')
    parser.add_argument('--visible', action='store_true', help='Show browser window')
    args = parser.parse_args()

    print("=" * 60)
    print("DRAFTKINGS SGP ODDS SCRAPER")
    print("=" * 60)

    # Test with sample games
    test_games = [
        {
            'league': 'NFL',
            'away_team': 'Patriots',
            'home_team': 'Broncos',
            'underdog': 'Broncos',
            'underdog_spread': 5.5,
            'total': 40.5,
            'parlay_underdog_under': 264,
            'parlay_favorite_over': 264,
        },
        {
            'league': 'NFL',
            'away_team': 'Rams',
            'home_team': 'Seahawks',
            'underdog': 'Rams',
            'underdog_spread': 3.0,
            'total': 47.0,
            'parlay_underdog_under': 250,
            'parlay_favorite_over': 282,
        }
    ]

    results = scrape_draftkings_sgp(test_games, headless=not args.visible)

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)

    for game in results:
        print(f"\n{game['away_team']} @ {game['home_team']}")
        print(f"  Hoop88 {game['underdog']} + Under: {game.get('parlay_underdog_under', 'N/A'):+d}")
        print(f"  DK SGP {game['underdog']} + Under: {game.get('dk_underdog_under', 'N/A')}")
        if game.get('dk_error'):
            print(f"  Error: {game['dk_error']}")
