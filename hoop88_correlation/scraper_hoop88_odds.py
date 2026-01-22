#!/usr/bin/env python3
"""
Hoop88 Football Odds Scraper
Scrapes spreads and totals for NFL and NCAAF across all periods,
builds correlated parlays to get parlay odds from bet slip.
"""

import os
import re
from datetime import datetime
from bs4 import BeautifulSoup
from dotenv import load_dotenv

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

load_dotenv()

HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")

# Available leagues and periods on hoop88
LEAGUES = {
    'NFL': {
        'selector': 'NFL',
        'periods': ['NFL', '1st Half', '1st Quarter', '2nd Quarter', '3rd Quarter', '4th Quarter']
    },
    'NCAAF': {
        'selector': 'NCAA Football',
        'periods': ['NCAA Football', '1st Half', '1st Quarter']
    }
}

# Short names for display
PERIOD_SHORT_NAMES = {
    'NFL': 'FG',
    'NCAA Football': 'FG',
    '1st Half': '1H',
    '1st Quarter': '1Q',
    '2nd Quarter': '2Q',
    '3rd Quarter': '3Q',
    '4th Quarter': '4Q',
}


def parse_odds(odds_text: str) -> int:
    """Extract American odds from text like '-110' or '+105'."""
    match = re.search(r'([+-]\d+)', odds_text)
    return int(match.group(1)) if match else -110


def parse_spread(spread_text: str) -> float:
    """Extract spread value from text like '-4½' or '+3'."""
    text = spread_text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
    match = re.search(r'([+-]?\d+\.?\d*)', text)
    return float(match.group(1)) if match else 0.0


def parse_total(total_text: str) -> float:
    """Extract total value from text like 'O 45½' or 'U 32'."""
    text = total_text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
    match = re.search(r'(\d+\.?\d*)', text)
    return float(match.group(1)) if match else 0.0


def parse_spread_odds(raw_text: str) -> tuple:
    """
    Parse spread and odds from raw text like '-5½ -110' or '+3 -120' or '-½ +105'.
    Returns (spread_float, odds_int)
    """
    if not raw_text:
        return (0.0, -110)

    # Replace fractions with decimal equivalents
    # Use regex to handle cases like -½ -> -0.5 and -5½ -> -5.5
    text = raw_text
    text = re.sub(r'(\d)½', r'\g<1>.5', text)  # 5½ -> 5.5
    text = re.sub(r'(\d)¼', r'\g<1>.25', text)  # 5¼ -> 5.25
    text = re.sub(r'(\d)¾', r'\g<1>.75', text)  # 5¾ -> 5.75
    text = re.sub(r'([+-])½', r'\g<1>0.5', text)  # -½ -> -0.5, +½ -> +0.5
    text = re.sub(r'([+-])¼', r'\g<1>0.25', text)
    text = re.sub(r'([+-])¾', r'\g<1>0.75', text)
    text = text.replace('½', '0.5').replace('¼', '0.25').replace('¾', '0.75')  # standalone

    # Match patterns like "-5.5 -110" or "+3 -120" or "-0.5 +105"
    # Pattern: optional sign, digits with optional decimal, space, sign + digits for odds
    match = re.search(r'([+-]?\d+\.?\d*)\s*([+-]\d+)', text)
    if match:
        spread = float(match.group(1))
        odds = int(match.group(2))
        return (spread, odds)

    return (0.0, -110)


def parse_total_odds(raw_text: str) -> tuple:
    """
    Parse total and odds from raw text like 'O 40½ -110' or 'U 47 -110'.
    Returns (total_float, odds_int)
    """
    if not raw_text:
        return (0.0, -110)

    # Replace fractions
    text = raw_text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')

    # Match patterns like "O 40.5 -110" or "U 47 -110"
    match = re.search(r'[OU]\s*(\d+\.?\d*)\s*([+-]\d+)', text)
    if match:
        total = float(match.group(1))
        odds = int(match.group(2))
        return (total, odds)

    return (0.0, -110)


def american_to_decimal(american: int) -> float:
    """Convert American odds to decimal."""
    if american == 0:
        return 1.91  # Default to -110 equivalent
    if american > 0:
        return (american / 100) + 1
    else:
        return (100 / abs(american)) + 1


def calculate_parlay_odds(leg1_american: int, leg2_american: int) -> int:
    """Calculate parlay odds from two legs (uncorrelated multiplication)."""
    dec1 = american_to_decimal(leg1_american)
    dec2 = american_to_decimal(leg2_american)
    combined_decimal = dec1 * dec2

    # Convert back to American
    if combined_decimal >= 2:
        return int(round((combined_decimal - 1) * 100))
    else:
        return int(round(-100 / (combined_decimal - 1)))


def scrape_hoop88_1h_odds(headless: bool = True):
    """
    Scrape 1H football odds from hoop88 and build correlated parlays.

    Returns list of games with:
    - team names
    - 1H spreads and odds
    - 1H totals and odds
    - Calculated parlay odds for correlated combinations
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    games = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Login if needed
        username_field = page.locator('input[name="customerID"]')
        if username_field.count() > 0 and username_field.is_visible():
            print("Logging in...")
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')
            page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1500)
            print("✅ Login successful")

        # Navigate to Football section
        print("Navigating to Football...")
        page.click('text=FOOTBALL')
        page.wait_for_timeout(2000)

        # Wait for game data to load - look for line elements
        print("Waiting for game data to load...")
        try:
            page.wait_for_selector('[data-panel="line"]', timeout=10000)
        except:
            print("No [data-panel='line'] found, trying alternatives...")

        # Wait a bit more for dynamic content
        page.wait_for_timeout(2000)

        # Take screenshot for debugging
        page.screenshot(path="/Users/callancapitolo/NFLWork/hoop88_correlation/debug_schedule.png")

        # Now navigate to 1st Half via the sidebar menu
        print("Looking for 1st Half option...")
        try:
            # The 1st Half option is in the sidebar under FOOTBALL
            # Look for the checkbox or link with "1st Half" text
            # Try multiple selectors
            selectors_to_try = [
                'label:has-text("1st Half")',
                'span:has-text("1st Half")',
                'a:has-text("1st Half")',
                '[data-sub-type*="1H"]',
                '[data-sport-sub-type*="1H"]',
            ]

            clicked = False
            for selector in selectors_to_try:
                try:
                    elem = page.locator(selector)
                    if elem.count() > 0 and elem.first.is_visible():
                        elem.first.click()
                        page.wait_for_timeout(2000)
                        print(f"✅ Clicked 1st Half using: {selector}")
                        clicked = True
                        break
                except:
                    continue

            if not clicked:
                print("Could not find 1st Half option - continuing with current view")

        except Exception as e:
            print(f"Error navigating to 1st Half: {e}")

        # Take another screenshot after navigation
        page.screenshot(path="/Users/callancapitolo/NFLWork/hoop88_correlation/debug_1h.png")

        # Get the page content
        html = page.content()

        # Save HTML for analysis
        with open("/Users/callancapitolo/NFLWork/hoop88_correlation/debug_schedule.html", "w") as f:
            f.write(html)

        print("Parsing game data using JavaScript evaluation...")

        # Use a single JavaScript evaluation to extract all game data
        game_data = page.evaluate('''() => {
            const games = [];
            const linePanels = document.querySelectorAll('[data-panel="line"].GAME');

            linePanels.forEach(panel => {
                try {
                    const firstTeam = panel.querySelector('[data-field="first-team"]');
                    const secondTeam = panel.querySelector('[data-field="second-team"]');
                    if (!firstTeam || !secondTeam) return;

                    const leagueName = panel.querySelector('[data-field="league-name"]');
                    const lineRows = panel.querySelectorAll('.lines .line');
                    if (lineRows.length < 2) return;

                    const awayRow = lineRows[0];
                    const homeRow = lineRows[1];

                    // Get text from the first span inside each group (contains line + odds)
                    function getGroupText(row, groupClass) {
                        const group = row.querySelector('.' + groupClass);
                        if (!group) return '';
                        const span = group.querySelector('span');
                        if (span) return span.textContent.trim();
                        // Fallback to first text node
                        const text = group.textContent || '';
                        const lines = text.split('\\n').filter(l => l.trim());
                        return lines[0]?.trim() || '';
                    }

                    games.push({
                        league: leagueName ? leagueName.textContent.trim() : '',
                        away_team: firstTeam.textContent.trim(),
                        home_team: secondTeam.textContent.trim(),
                        away_spread_raw: getGroupText(awayRow, 'group-2'),
                        home_spread_raw: getGroupText(homeRow, 'group-2'),
                        over_raw: getGroupText(awayRow, 'group-4'),
                        under_raw: getGroupText(homeRow, 'group-4'),
                        away_ml_raw: getGroupText(awayRow, 'group-3'),
                        home_ml_raw: getGroupText(homeRow, 'group-3'),
                    });
                } catch (e) {
                    console.error('Error parsing panel:', e);
                }
            });

            return games;
        }''')

        print(f"\nFound {len(game_data)} games with raw data:")
        for game in game_data:
            print(f"  {game['away_team']} @ {game['home_team']}")
            print(f"    Spread raw: '{game['away_spread_raw']}' / '{game['home_spread_raw']}'")
            print(f"    Total raw: '{game['over_raw']}' / '{game['under_raw']}'")
            print(f"    ML raw: '{game['away_ml_raw']}' / '{game['home_ml_raw']}'")

            # Parse the raw values
            away_spread, away_spread_odds = parse_spread_odds(game['away_spread_raw'])
            home_spread, home_spread_odds = parse_spread_odds(game['home_spread_raw'])
            total, over_odds = parse_total_odds(game['over_raw'])
            _, under_odds = parse_total_odds(game['under_raw'])
            away_ml = parse_odds(game['away_ml_raw'])
            home_ml = parse_odds(game['home_ml_raw'])

            # Skip games with missing data
            if not away_spread or not home_spread:
                print(f"    ⚠️ Skipping - missing spread data")
                continue
            if over_odds == 0 or under_odds == 0:
                print(f"    ⚠️ Skipping - missing total odds")
                continue

            # Determine favorite/underdog
            if away_spread > 0:
                underdog = game['away_team']
                underdog_spread = away_spread
                underdog_spread_odds = away_spread_odds
                favorite = game['home_team']
                favorite_spread = home_spread
                favorite_spread_odds = home_spread_odds
            else:
                underdog = game['home_team']
                underdog_spread = home_spread
                underdog_spread_odds = home_spread_odds
                favorite = game['away_team']
                favorite_spread = away_spread
                favorite_spread_odds = away_spread_odds

            # Calculate correlated parlay odds (uncorrelated multiplication)
            parlay1_odds = calculate_parlay_odds(underdog_spread_odds, under_odds)
            parlay2_odds = calculate_parlay_odds(favorite_spread_odds, over_odds)

            game_result = {
                'league': game['league'],
                'away_team': game['away_team'],
                'home_team': game['home_team'],
                'away_spread': away_spread,
                'away_spread_odds': away_spread_odds,
                'home_spread': home_spread,
                'home_spread_odds': home_spread_odds,
                'total': total,
                'over_odds': over_odds,
                'under_odds': under_odds,
                'away_ml': away_ml,
                'home_ml': home_ml,
                'underdog': underdog,
                'underdog_spread': underdog_spread,
                'underdog_spread_odds': underdog_spread_odds,
                'favorite': favorite,
                'favorite_spread': favorite_spread,
                'favorite_spread_odds': favorite_spread_odds,
                'parlay_underdog_under': parlay1_odds,
                'parlay_favorite_over': parlay2_odds,
            }

            print(f"    Parsed: Spread {away_spread}({away_spread_odds:+d})/{home_spread}({home_spread_odds:+d})")
            print(f"    Parsed: Total {total} O({over_odds:+d})/U({under_odds:+d})")
            print(f"    Correlated Parlays:")
            print(f"      {underdog} {underdog_spread:+.1f} + Under {total}: {parlay1_odds:+d}")
            print(f"      {favorite} {favorite_spread:+.1f} + Over {total}: {parlay2_odds:+d}")

            games.append(game_result)

        browser.close()

    return games


def navigate_to_period(page, league: str, period: str) -> bool:
    """
    Navigate to a specific league and period on hoop88.

    The sidebar has multiple sections with same period names (e.g., multiple "1st Half").
    We need to find the period that belongs to the correct league section.

    Args:
        league: 'NFL' or 'NCAAF'
        period: Full period name like 'NFL', 'NCAA Football', '1st Half', '1st Quarter'

    Returns True if successful.
    """
    try:
        league_selector = LEAGUES[league]['selector']  # 'NFL' or 'NCAA Football'

        # For full game, just click the league
        if period == league_selector:
            clicked = page.evaluate(f'''() => {{
                const elements = document.querySelectorAll('[data-sport-sub-type="{period}"]');
                for (const elem of elements) {{
                    elem.click();
                    return true;
                }}
                return false;
            }}''')
            if not clicked:
                print(f"  League '{period}' not found")
                return False
            page.wait_for_timeout(2500)
            return True

        # For periods (1st Half, 1st Quarter, etc.), we need to find the one
        # that belongs to our league section. The sidebar structure is ordered:
        # NFL -> 1st Half (NFL) -> 1st Quarter (NFL) -> NCAA Football -> 1st Half (NCAAF) -> ...
        #
        # So we find all sport-sub-type elements, locate our league, then find
        # the period that comes after it (before the next league).

        clicked = page.evaluate(f'''() => {{
            const allElements = Array.from(document.querySelectorAll('[data-sport-sub-type]'));
            const leagueSelector = "{league_selector}";
            const periodSelector = "{period}";

            // Find the index of our league
            let leagueIdx = -1;
            for (let i = 0; i < allElements.length; i++) {{
                if (allElements[i].getAttribute('data-sport-sub-type') === leagueSelector) {{
                    leagueIdx = i;
                    break;
                }}
            }}

            if (leagueIdx === -1) {{
                console.log('League not found: ' + leagueSelector);
                return false;
            }}

            // Now find the period after the league but before the next major section
            // Major sections are: NFL, NCAA Football, NBA, NHL, etc.
            const majorSections = ['NFL', 'NCAA Football', 'NBA', 'NHL', 'NCAA Basketball'];

            for (let i = leagueIdx + 1; i < allElements.length; i++) {{
                const subType = allElements[i].getAttribute('data-sport-sub-type');

                // Stop if we hit another major section
                if (majorSections.includes(subType) && subType !== leagueSelector) {{
                    break;
                }}

                // Found our period!
                if (subType === periodSelector) {{
                    allElements[i].click();
                    return true;
                }}
            }}

            console.log('Period not found in league section: ' + periodSelector);
            return false;
        }}''')

        if not clicked:
            print(f"  Period '{period}' not found in {league} section")
            return False

        page.wait_for_timeout(2500)
        return True

    except Exception as e:
        print(f"  Error navigating to {league}/{period}: {e}")
        return False


def get_parlay_odds_from_betslip(page, risk_amount: int = 100) -> int:
    """
    Read the actual parlay odds from the bet slip.
    Assumes legs have already been added to bet slip.
    Returns American odds.
    """
    # Click on PARLAY tab via JavaScript
    page.evaluate('''() => {
        const tabs = document.querySelectorAll('[data-wager-type="Parlay"]');
        for (const tab of tabs) {
            if (tab.offsetParent !== null) {
                tab.click();
                return true;
            }
        }
        if (tabs.length > 0) tabs[0].click();
        return true;
    }''')
    page.wait_for_timeout(1000)

    # Enter risk amount via JS
    page.evaluate(f'''() => {{
        const inputs = document.querySelectorAll('[data-wager-panel="Parlay"] input[data-field="risk"]');
        for (const input of inputs) {{
            if (input.offsetParent !== null) {{
                input.value = '{risk_amount}';
                input.dispatchEvent(new Event('input', {{ bubbles: true }}));
                input.dispatchEvent(new Event('change', {{ bubbles: true }}));
                input.dispatchEvent(new Event('keyup', {{ bubbles: true }}));
                return true;
            }}
        }}
        return false;
    }}''')
    page.wait_for_timeout(1500)

    # Read the To Win value
    to_win_value = page.evaluate('''() => {
        const inputs = document.querySelectorAll('[data-wager-panel="Parlay"] input[data-field="win"]');
        for (const input of inputs) {
            if (input.offsetParent !== null && input.value) {
                return input.value;
            }
        }
        for (const input of inputs) {
            if (input.value) return input.value;
        }
        return '';
    }''')

    if to_win_value and to_win_value.strip():
        try:
            win = float(to_win_value.replace('$', '').replace(',', ''))
            if win >= risk_amount:
                return int(round(win))  # Plus odds
            else:
                return int(round(-100 * risk_amount / win))  # Minus odds
        except:
            return 0
    return 0


def clear_betslip(page):
    """Clear all bets from the bet slip."""
    page.evaluate('''() => {
        const clearBtn = document.querySelector('[data-action="clear-all"]');
        if (clearBtn) clearBtn.click();
    }''')
    page.wait_for_timeout(800)


def scrape_hoop88_with_actual_odds(headless: bool = True):
    """
    Scrape hoop88 football odds AND get ACTUAL parlay odds from bet slip.

    Returns list of games with:
    - team names, spreads and odds, totals and odds
    - ACTUAL parlay odds from hoop88 (not calculated)
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    games = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Login
        username_field = page.locator('input[name="customerID"]')
        if username_field.count() > 0 and username_field.is_visible():
            print("Logging in...")
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')
            page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1500)
            print("✅ Login successful")

        # Navigate to Football
        print("Navigating to Football...")
        page.click('text=FOOTBALL')
        page.wait_for_timeout(3000)

        try:
            page.wait_for_selector('[data-panel="line"]', timeout=10000)
        except:
            pass

        # Extract game data using JS
        game_data = page.evaluate('''() => {
            const games = [];
            const linePanels = document.querySelectorAll('[data-panel="line"].GAME');

            linePanels.forEach((panel, idx) => {
                try {
                    const firstTeam = panel.querySelector('[data-field="first-team"]');
                    const secondTeam = panel.querySelector('[data-field="second-team"]');
                    if (!firstTeam || !secondTeam) return;

                    const leagueName = panel.querySelector('[data-field="league-name"]');
                    const lineRows = panel.querySelectorAll('.lines .line');
                    if (lineRows.length < 2) return;

                    const awayRow = lineRows[0];
                    const homeRow = lineRows[1];

                    function getGroupText(row, groupClass) {
                        const group = row.querySelector('.' + groupClass);
                        if (!group) return '';
                        const span = group.querySelector('span');
                        if (span) return span.textContent.trim();
                        return '';
                    }

                    games.push({
                        panelIndex: idx,
                        league: leagueName ? leagueName.textContent.trim() : '',
                        away_team: firstTeam.textContent.trim(),
                        home_team: secondTeam.textContent.trim(),
                        away_spread_raw: getGroupText(awayRow, 'group-2'),
                        home_spread_raw: getGroupText(homeRow, 'group-2'),
                        over_raw: getGroupText(awayRow, 'group-4'),
                        under_raw: getGroupText(homeRow, 'group-4'),
                    });
                } catch (e) {}
            });

            return games;
        }''')

        print(f"Found {len(game_data)} games")

        for game_info in game_data:
            print(f"\nProcessing: {game_info['away_team']} @ {game_info['home_team']}")

            # Parse spreads and totals
            away_spread, away_spread_odds = parse_spread_odds(game_info['away_spread_raw'])
            home_spread, home_spread_odds = parse_spread_odds(game_info['home_spread_raw'])
            total, over_odds = parse_total_odds(game_info['over_raw'])
            _, under_odds = parse_total_odds(game_info['under_raw'])

            if not away_spread or not home_spread:
                print(f"  Skipping - missing spread data")
                continue

            # Determine favorite/underdog
            if away_spread > 0:
                underdog, underdog_spread, underdog_spread_odds = game_info['away_team'], away_spread, away_spread_odds
                favorite, favorite_spread, favorite_spread_odds = game_info['home_team'], home_spread, home_spread_odds
                underdog_is_away = True
            else:
                underdog, underdog_spread, underdog_spread_odds = game_info['home_team'], home_spread, home_spread_odds
                favorite, favorite_spread, favorite_spread_odds = game_info['away_team'], away_spread, away_spread_odds
                underdog_is_away = False

            panel = page.locator('[data-panel="line"].GAME').nth(game_info['panelIndex'])

            # === Get Underdog + Under parlay odds ===
            clear_betslip(page)
            page.wait_for_timeout(500)

            # Click underdog spread
            if underdog_is_away:
                spread_btn = panel.locator('.lines .line').first.locator('.group-2 .line-play').first
            else:
                spread_btn = panel.locator('.lines .line').nth(1).locator('.group-2 .line-play').first
            spread_btn.click()
            page.wait_for_timeout(800)

            # Click under total (always second row)
            under_btn = panel.locator('.lines .line').nth(1).locator('.group-4 .line-play').first
            under_btn.click()
            page.wait_for_timeout(800)

            parlay_underdog_under = get_parlay_odds_from_betslip(page)
            print(f"  {underdog} {underdog_spread:+.1f} + Under {total}: +{parlay_underdog_under}")

            # === Get Favorite + Over parlay odds ===
            clear_betslip(page)
            page.wait_for_timeout(500)

            # Click favorite spread
            if underdog_is_away:
                spread_btn = panel.locator('.lines .line').nth(1).locator('.group-2 .line-play').first
            else:
                spread_btn = panel.locator('.lines .line').first.locator('.group-2 .line-play').first
            spread_btn.click()
            page.wait_for_timeout(800)

            # Click over total (always first row)
            over_btn = panel.locator('.lines .line').first.locator('.group-4 .line-play').first
            over_btn.click()
            page.wait_for_timeout(800)

            parlay_favorite_over = get_parlay_odds_from_betslip(page)
            print(f"  {favorite} {favorite_spread:+.1f} + Over {total}: +{parlay_favorite_over}")

            # Calculate what we would have calculated (for comparison)
            calc_underdog_under = calculate_parlay_odds(underdog_spread_odds, under_odds)
            calc_favorite_over = calculate_parlay_odds(favorite_spread_odds, over_odds)

            game_result = {
                'league': game_info['league'],
                'away_team': game_info['away_team'],
                'home_team': game_info['home_team'],
                'away_spread': away_spread,
                'away_spread_odds': away_spread_odds,
                'home_spread': home_spread,
                'home_spread_odds': home_spread_odds,
                'total': total,
                'over_odds': over_odds,
                'under_odds': under_odds,
                'underdog': underdog,
                'underdog_spread': underdog_spread,
                'underdog_spread_odds': underdog_spread_odds,
                'favorite': favorite,
                'favorite_spread': favorite_spread,
                'favorite_spread_odds': favorite_spread_odds,
                'parlay_underdog_under': parlay_underdog_under,
                'parlay_favorite_over': parlay_favorite_over,
                'calc_underdog_under': calc_underdog_under,
                'calc_favorite_over': calc_favorite_over,
            }

            games.append(game_result)

        clear_betslip(page)
        browser.close()

    return games


def scrape_hoop88_all_markets(leagues: list = None, periods: list = None, headless: bool = True):
    """
    Scrape hoop88 football odds for multiple leagues and periods.

    Args:
        leagues: List of leagues to scrape ['NFL', 'NCAAF']. Default: both.
        periods: List of periods to scrape ['FG', '1H', '1Q', '2Q', '3Q', '4Q']. Default: all.
        headless: Run browser in headless mode.

    Returns:
        List of games with actual parlay odds from bet slip.
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    if leagues is None:
        leagues = ['NFL', 'NCAAF']
    if periods is None:
        periods = ['FG', '1H', '1Q', '2Q', '3Q', '4Q']

    # Map short names to full period names
    period_map = {
        'FG': None,  # Full game - use league selector
        '1H': '1st Half',
        '1Q': '1st Quarter',
        '2Q': '2nd Quarter',
        '3Q': '3rd Quarter',
        '4Q': '4th Quarter',
    }

    all_games = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Login
        username_field = page.locator('input[name="customerID"]')
        if username_field.count() > 0 and username_field.is_visible():
            print("Logging in...")
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')
            page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1500)
            print("✅ Login successful")

        # Navigate to Football first
        print("Navigating to Football...")
        page.click('text=FOOTBALL')
        page.wait_for_timeout(3000)

        for league in leagues:
            if league not in LEAGUES:
                print(f"Unknown league: {league}")
                continue

            league_config = LEAGUES[league]
            available_periods = league_config['periods']

            for period_short in periods:
                period_full = period_map.get(period_short)

                # For full game, use league selector
                if period_full is None:
                    period_full = league_config['selector']

                # Check if this period is available for this league
                if period_full not in available_periods:
                    continue

                print(f"\n{'='*60}")
                print(f"Scraping {league} {period_short}")
                print(f"{'='*60}")

                # Navigate to the period
                if not navigate_to_period(page, league, period_full):
                    print(f"  Skipping - could not navigate to {period_full}")
                    continue

                page.wait_for_timeout(2000)

                # Extract game data
                game_data = page.evaluate('''() => {
                    const games = [];
                    const linePanels = document.querySelectorAll('[data-panel="line"].GAME');

                    linePanels.forEach((panel, idx) => {
                        try {
                            const firstTeam = panel.querySelector('[data-field="first-team"]');
                            const secondTeam = panel.querySelector('[data-field="second-team"]');
                            if (!firstTeam || !secondTeam) return;

                            const leagueName = panel.querySelector('[data-field="league-name"]');
                            const lineRows = panel.querySelectorAll('.lines .line');
                            if (lineRows.length < 2) return;

                            const awayRow = lineRows[0];
                            const homeRow = lineRows[1];

                            function getGroupText(row, groupClass) {
                                const group = row.querySelector('.' + groupClass);
                                if (!group) return '';
                                const span = group.querySelector('span');
                                if (span) return span.textContent.trim();
                                return '';
                            }

                            games.push({
                                panelIndex: idx,
                                league: leagueName ? leagueName.textContent.trim() : '',
                                away_team: firstTeam.textContent.trim(),
                                home_team: secondTeam.textContent.trim(),
                                away_spread_raw: getGroupText(awayRow, 'group-2'),
                                home_spread_raw: getGroupText(homeRow, 'group-2'),
                                over_raw: getGroupText(awayRow, 'group-4'),
                                under_raw: getGroupText(homeRow, 'group-4'),
                            });
                        } catch (e) {}
                    });

                    return games;
                }''')

                print(f"Found {len(game_data)} games")

                for game_info in game_data:
                    print(f"\nProcessing: {game_info['away_team']} @ {game_info['home_team']}")

                    # Parse spreads and totals
                    away_spread, away_spread_odds = parse_spread_odds(game_info['away_spread_raw'])
                    home_spread, home_spread_odds = parse_spread_odds(game_info['home_spread_raw'])
                    total, over_odds = parse_total_odds(game_info['over_raw'])
                    _, under_odds = parse_total_odds(game_info['under_raw'])

                    if not away_spread or not home_spread:
                        print(f"  Skipping - missing spread data")
                        continue

                    # Determine favorite/underdog
                    if away_spread > 0:
                        underdog, underdog_spread, underdog_spread_odds = game_info['away_team'], away_spread, away_spread_odds
                        favorite, favorite_spread, favorite_spread_odds = game_info['home_team'], home_spread, home_spread_odds
                        underdog_is_away = True
                    else:
                        underdog, underdog_spread, underdog_spread_odds = game_info['home_team'], home_spread, home_spread_odds
                        favorite, favorite_spread, favorite_spread_odds = game_info['away_team'], away_spread, away_spread_odds
                        underdog_is_away = False

                    panel = page.locator('[data-panel="line"].GAME').nth(game_info['panelIndex'])

                    # === Get Underdog + Under parlay odds ===
                    clear_betslip(page)
                    page.wait_for_timeout(500)

                    try:
                        # Click underdog spread
                        if underdog_is_away:
                            spread_btn = panel.locator('.lines .line').first.locator('.group-2 .line-play').first
                        else:
                            spread_btn = panel.locator('.lines .line').nth(1).locator('.group-2 .line-play').first
                        spread_btn.click()
                        page.wait_for_timeout(800)

                        # Click under total (always second row)
                        under_btn = panel.locator('.lines .line').nth(1).locator('.group-4 .line-play').first
                        under_btn.click()
                        page.wait_for_timeout(800)

                        parlay_underdog_under = get_parlay_odds_from_betslip(page)
                        print(f"  {underdog} {underdog_spread:+.1f} + Under {total}: +{parlay_underdog_under}")

                        # === Get Favorite + Over parlay odds ===
                        clear_betslip(page)
                        page.wait_for_timeout(500)

                        # Click favorite spread
                        if underdog_is_away:
                            spread_btn = panel.locator('.lines .line').nth(1).locator('.group-2 .line-play').first
                        else:
                            spread_btn = panel.locator('.lines .line').first.locator('.group-2 .line-play').first
                        spread_btn.click()
                        page.wait_for_timeout(800)

                        # Click over total (always first row)
                        over_btn = panel.locator('.lines .line').first.locator('.group-4 .line-play').first
                        over_btn.click()
                        page.wait_for_timeout(800)

                        parlay_favorite_over = get_parlay_odds_from_betslip(page)
                        print(f"  {favorite} {favorite_spread:+.1f} + Over {total}: +{parlay_favorite_over}")

                    except Exception as e:
                        print(f"  Error getting parlay odds: {e}")
                        parlay_underdog_under = 0
                        parlay_favorite_over = 0

                    # Calculate what we would have calculated (for comparison)
                    calc_underdog_under = calculate_parlay_odds(underdog_spread_odds, under_odds)
                    calc_favorite_over = calculate_parlay_odds(favorite_spread_odds, over_odds)

                    game_result = {
                        'league': league,
                        'period': period_short,
                        'away_team': game_info['away_team'],
                        'home_team': game_info['home_team'],
                        'away_spread': away_spread,
                        'away_spread_odds': away_spread_odds,
                        'home_spread': home_spread,
                        'home_spread_odds': home_spread_odds,
                        'total': total,
                        'over_odds': over_odds,
                        'under_odds': under_odds,
                        'underdog': underdog,
                        'underdog_spread': underdog_spread,
                        'underdog_spread_odds': underdog_spread_odds,
                        'favorite': favorite,
                        'favorite_spread': favorite_spread,
                        'favorite_spread_odds': favorite_spread_odds,
                        'parlay_underdog_under': parlay_underdog_under,
                        'parlay_favorite_over': parlay_favorite_over,
                        'calc_underdog_under': calc_underdog_under,
                        'calc_favorite_over': calc_favorite_over,
                    }

                    all_games.append(game_result)

        clear_betslip(page)
        browser.close()

    return all_games


def scrape_with_betslip_interaction(headless: bool = False):
    """
    Full scrape with bet slip interaction to get ACTUAL parlay odds from hoop88.

    This version:
    1. Scrapes the schedule for games
    2. For each game, clicks spread + total to add to bet slip
    3. Reads the actual parlay odds from bet slip
    4. Compares to calculated (uncorrelated) parlay odds
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    games = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Login if needed
        username_field = page.locator('input[name="customerID"]')
        if username_field.count() > 0 and username_field.is_visible():
            print("Logging in...")
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')
            page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1500)
            print("✅ Login successful")

        # Navigate to Football
        print("Navigating to Football...")
        page.click('text=FOOTBALL')
        page.wait_for_timeout(2000)

        # Wait for games to load
        try:
            page.wait_for_selector('[data-panel="line"]', timeout=10000)
        except:
            pass
        page.wait_for_timeout(2000)

        # First, let's extract game data
        game_data = page.evaluate('''() => {
            const games = [];
            const linePanels = document.querySelectorAll('[data-panel="line"].GAME');

            linePanels.forEach((panel, idx) => {
                const firstTeam = panel.querySelector('[data-field="first-team"]');
                const secondTeam = panel.querySelector('[data-field="second-team"]');
                if (!firstTeam || !secondTeam) return;

                const lineRows = panel.querySelectorAll('.lines .line');
                if (lineRows.length < 2) return;

                games.push({
                    panelIndex: idx,
                    away_team: firstTeam.textContent.trim(),
                    home_team: secondTeam.textContent.trim(),
                });
            });
            return games;
        }''')

        print(f"Found {len(game_data)} games")

        # For each game, test the parlay by clicking on bet slip
        for game_info in game_data:
            print(f"\n Testing {game_info['away_team']} @ {game_info['home_team']}...")

            # Click on away spread to add to bet slip
            try:
                panel = page.locator(f'[data-panel="line"].GAME').nth(game_info['panelIndex'])
                away_spread_btn = panel.locator('.lines .line').first.locator('.group-2 .line-play').first

                if away_spread_btn.count() > 0:
                    away_spread_btn.click()
                    page.wait_for_timeout(1000)
                    print("  Clicked away spread")

                # Click on over total
                over_btn = panel.locator('.lines .line').first.locator('.group-4 .line-play').first
                if over_btn.count() > 0:
                    over_btn.click()
                    page.wait_for_timeout(1000)
                    print("  Clicked over total")

                # Now read the parlay odds from bet slip
                page.wait_for_timeout(1500)

                # Click on PARLAY tab to see combined parlay odds
                try:
                    # Try multiple selectors for the PARLAY tab
                    parlay_clicked = False
                    parlay_selectors = [
                        'button:has-text("PARLAY")',
                        'text=PARLAY',
                        '.betting-stage button >> text=PARLAY',
                        '[data-type="parlay"]',
                    ]
                    for sel in parlay_selectors:
                        try:
                            tab = page.locator(sel)
                            if tab.count() > 0:
                                tab.first.click()
                                page.wait_for_timeout(500)
                                parlay_clicked = True
                                print(f"  Clicked PARLAY tab using: {sel}")
                                break
                        except:
                            continue

                    if not parlay_clicked:
                        # Try JavaScript click
                        page.evaluate('''() => {
                            const buttons = document.querySelectorAll('button');
                            for (const btn of buttons) {
                                if (btn.textContent.includes('PARLAY')) {
                                    btn.click();
                                    return true;
                                }
                            }
                            return false;
                        }''')
                        print("  Tried JS click for PARLAY")

                    page.wait_for_timeout(1000)

                    # Enter a bet amount to trigger payout calculation
                    bet_input = page.locator('input[type="text"], input[type="number"]').filter(has_text='').first
                    if bet_input.count() > 0:
                        bet_input.fill('100')
                        page.wait_for_timeout(500)
                        print("  Entered $100 bet amount")

                except Exception as e:
                    print(f"  Could not interact with parlay: {e}")

                page.screenshot(path=f"/Users/callancapitolo/NFLWork/hoop88_correlation/debug_betslip_{game_info['panelIndex']}.png")

                # Try to find parlay odds in bet slip
                parlay_info = page.evaluate('''() => {
                    // Look for the parlay odds display
                    // After clicking PARLAY tab, look for the combined odds
                    const result = {
                        raw_text: '',
                        odds: null
                    };

                    // The bet slip panel
                    const betSlip = document.querySelector('.slide-up[data-panel="bets"], .betting-stage');
                    if (betSlip) {
                        result.raw_text = betSlip.textContent.substring(0, 1000);

                        // Look for odds patterns like "+264" or "-110"
                        const oddsMatch = result.raw_text.match(/To Win[\\s\\S]*?([+-]\\d+)/);
                        if (oddsMatch) {
                            result.odds = oddsMatch[1];
                        }

                        // Also try to find "Parlay" section odds
                        const parlaySection = betSlip.querySelector('.parlay, [data-type="parlay"]');
                        if (parlaySection) {
                            result.parlay_section = parlaySection.textContent;
                        }
                    }

                    return result;
                }''')

                print(f"  Parlay info: {parlay_info}")

                # Clear the bet slip
                try:
                    clear_btn = page.locator('text=Clear All, [data-action="clear-all"]')
                    if clear_btn.count() > 0:
                        clear_btn.first.click()
                        page.wait_for_timeout(1000)
                except:
                    pass

            except Exception as e:
                print(f"  Error: {e}")

        # Keep browser open for inspection if not headless
        if not headless:
            print("\nBrowser left open for inspection. Press Enter to close...")
            input()

        browser.close()

    return games


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape hoop88 football odds')
    parser.add_argument('--visible', action='store_true', help='Show browser window')
    parser.add_argument('--all-markets', action='store_true', help='Scrape all leagues and periods')
    parser.add_argument('--leagues', nargs='+', default=['NFL', 'NCAAF'], help='Leagues to scrape (default: NFL NCAAF)')
    parser.add_argument('--periods', nargs='+', default=['FG', '1H', '1Q'], help='Periods to scrape (default: FG 1H 1Q)')
    args = parser.parse_args()

    print("=" * 60)
    print("HOOP88 FOOTBALL ODDS SCRAPER")
    print("=" * 60)

    if args.all_markets:
        print(f"Leagues: {args.leagues}")
        print(f"Periods: {args.periods}")
        games = scrape_hoop88_all_markets(
            leagues=args.leagues,
            periods=args.periods,
            headless=not args.visible
        )
    else:
        games = scrape_hoop88_with_actual_odds(headless=not args.visible)

    print(f"\n{'='*60}")
    print(f"SUMMARY: Scraped {len(games)} games")
    print(f"{'='*60}")

    for game in games:
        period = game.get('period', 'FG')
        print(f"\n{game['league']} {period} | {game['away_team']} @ {game['home_team']}")
        print(f"  {game['underdog']} {game['underdog_spread']:+.1f} + Under {game['total']}: +{game['parlay_underdog_under']}")
        print(f"  {game['favorite']} {game['favorite_spread']:+.1f} + Over {game['total']}: +{game['parlay_favorite_over']}")
