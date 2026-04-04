#!/usr/bin/env python3
"""
Pikkit Pro common functions — sport-agnostic session management,
API response parsing, and SGP odds building.

Extracted from hoop88_correlation/scraper_pikkit.py for reuse across sports.
"""

import os
import json
import re
from pathlib import Path
from dotenv import load_dotenv

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

load_dotenv()

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PIKKIT_URL = os.getenv("PIKKIT_URL", "https://app.pikkit.com")

# Session storage for persistent login — defaults to mlb_sgp/ directory
SESSION_FILE = Path(os.getenv("PIKKIT_SESSION_FILE", str(Path(__file__).parent / ".pikkit_session.json")))

# Period configuration per sport.
# Each sport maps human-readable period names to Pikkit <select> values,
# and provides a function that returns the filter category pill to click
# (e.g. "Halves", "Quarters") before selecting from the dropdown.
PERIOD_CONFIGS = {
    "NFL": {
        "period_to_select": {
            '1st Half': 'half_1', '2nd Half': 'half_2',
            '1st Quarter': 'quarter_1',
            '2nd Quarter': 'other', '3rd Quarter': 'other', '4th Quarter': 'other',
        },
        "filter_category": lambda period: "Halves" if "Half" in period else "Quarters" if "Quarter" in period else None,
    },
    "NCAAFB": {
        "period_to_select": {
            '1st Half': 'half_1', '2nd Half': 'half_2',
            '1st Quarter': 'quarter_1',
            '2nd Quarter': 'other', '3rd Quarter': 'other', '4th Quarter': 'other',
        },
        "filter_category": lambda period: "Halves" if "Half" in period else "Quarters" if "Quarter" in period else None,
    },
    "NBA": {
        "period_to_select": {
            '1st Half': 'half_1', '2nd Half': 'half_2',
            '1st Quarter': 'quarter_1',
            '2nd Quarter': 'other', '3rd Quarter': 'other', '4th Quarter': 'other',
        },
        "filter_category": lambda period: "Halves" if "Half" in period else "Quarters" if "Quarter" in period else None,
    },
    "MLB": {
        # MLB may not have sub-periods on Pikkit; Full Game is default.
        # F5 (first 5 innings) mapping TBD after manual recon.
        "period_to_select": {},
        "filter_category": lambda period: None,
    },
}


# ---------------------------------------------------------------------------
# Session management
# ---------------------------------------------------------------------------

def save_session(context):
    """Save browser session cookies and storage for reuse."""
    storage = context.storage_state()
    with open(SESSION_FILE, 'w') as f:
        json.dump(storage, f)
    print(f"Session saved to {SESSION_FILE}")


def load_session():
    """Load saved session if available."""
    if SESSION_FILE.exists():
        try:
            with open(SESSION_FILE, 'r') as f:
                storage = json.load(f)
            return storage
        except Exception as e:
            print(f"Warning: could not load session: {e}")
            return None
    return None


def is_logged_in(page) -> bool:
    """Check if user is logged in by looking for profile/balance indicators."""
    return page.evaluate(r'''() => {
        // Check for balance indicator (logged-in users see their balance)
        const balanceText = document.body.innerText;
        const hasBalance = /\$[\d,]+\.\d{2}/.test(balanceText);

        // Check for login button or login link (means not logged in)
        const allLinks = Array.from(document.querySelectorAll('a, button'));
        const loginBtn = allLinks.find(el => {
            const text = el.innerText.toLowerCase();
            const href = el.href || '';
            return text.includes('log in') || text.includes('sign in') || href.includes('login');
        });

        return hasBalance && !loginBtn;
    }''')


# ---------------------------------------------------------------------------
# Betslip utilities
# ---------------------------------------------------------------------------

def clear_betslip(page):
    """Clear all selections from the Pikkit betslip by clicking selected items to deselect."""
    # Simply reload or click Straight tab to clear - most reliable approach
    try:
        straight_tab = page.get_by_text("Straight", exact=True)
        if straight_tab.count() > 0:
            straight_tab.first.click()
            page.wait_for_timeout(500)
    except Exception:
        pass

    # Try clicking any blue/selected cells to deselect
    for _ in range(3):
        cleared = page.evaluate(r'''() => {
            let cleared = 0;
            // Find divs with blue background (selected cells)
            document.querySelectorAll('div').forEach(div => {
                const style = window.getComputedStyle(div);
                const bg = style.backgroundColor;
                // Blue background indicates selection (rgb values for blue)
                if (bg && (bg.includes('59, 130, 246') || bg.includes('37, 99, 235') ||
                           bg.includes('96, 165, 250'))) {
                    div.click();
                    cleared++;
                }
            });
            return cleared;
        }''')
        if cleared == 0:
            break
        page.wait_for_timeout(300)

    page.wait_for_timeout(300)


def click_parlay_tab(page) -> bool:
    """Click the Parlay tab in the Place Bets section."""
    clicked = page.evaluate('''() => {
        const buttons = document.querySelectorAll('button, div[role="button"]');
        for (const btn of buttons) {
            if (btn.innerText.trim() === 'Parlay') {
                btn.click();
                return true;
            }
        }
        return false;
    }''')

    if clicked:
        page.wait_for_timeout(1000)
    return clicked


# ---------------------------------------------------------------------------
# Odds math
# ---------------------------------------------------------------------------

def decimal_to_american(decimal_odds: float) -> int:
    """Convert decimal odds to American odds."""
    if decimal_odds >= 2.0:
        return int(round((decimal_odds - 1) * 100))
    else:
        return int(round(-100 / (decimal_odds - 1)))


# ---------------------------------------------------------------------------
# API data extraction
# ---------------------------------------------------------------------------

def extract_parlay_odds_from_api(betslip_data: dict) -> dict:
    """
    Extract SGP parlay odds from the Pikkit /betslip API response.

    Args:
        betslip_data: The JSON response from /betslip API

    Returns:
        Dict with:
        - best_odds: Best available odds (int, American)
        - all_odds: List of odds (int, American) from all books
        - book_odds: List of {book: str, odds: int} for filtering by book
        - leg_count: Number of legs in parlay
    """
    result = {
        'best_odds': None,
        'all_odds': [],
        'book_odds': [],
        'leg_count': 0
    }

    if not betslip_data:
        return result

    parlay = betslip_data.get('parlay')
    if not parlay:
        return result

    options = parlay.get('options', [])
    if not options:
        return result

    # Count legs from first option's outcomes
    if options[0].get('outcomes'):
        result['leg_count'] = len(options[0]['outcomes'])

    # Extract book names and odds from each option
    # Filter out books with adjusted/altered legs
    best_decimal = 0
    for opt in options:
        institution_name = opt.get('institution_name', '')
        decimal_odds = opt.get('odds')

        if not institution_name or not decimal_odds:
            continue

        # Check if any leg has been adjusted (different line than requested)
        outcomes = opt.get('outcomes', [])
        has_adjusted_leg = any(outcome.get('adjusted') for outcome in outcomes)
        if has_adjusted_leg:
            # Skip this book - their line doesn't match what we requested
            continue

        # Ensure decimal_odds is a float
        try:
            decimal_odds = float(decimal_odds)
        except (TypeError, ValueError):
            continue

        # Convert decimal to American
        american_odds = decimal_to_american(decimal_odds)

        result['all_odds'].append(american_odds)
        result['book_odds'].append({
            'book': institution_name,
            'odds': american_odds
        })

        if decimal_odds > best_decimal:
            best_decimal = decimal_odds
            result['best_odds'] = american_odds

    return result


# ---------------------------------------------------------------------------
# DOM fallback
# ---------------------------------------------------------------------------

def extract_parlay_odds(page) -> dict:
    """
    Extract SGP parlay odds from all sportsbooks shown in Pikkit.
    Falls back to DOM scraping if API data not available.

    Returns:
        Dict with:
        - best_odds: Best available odds (int)
        - all_odds: List of odds (int) from all books
        - book_odds: List of {book: str, odds: int} for filtering by book
        - leg_count: Number of legs in parlay
    """
    # Fallback to DOM scraping for leg count and best odds display
    result = page.evaluate(r'''() => {
        const result = {
            best_odds: null,
            all_odds: [],
            book_odds: [],
            leg_count: 0
        };

        const bodyText = document.body.innerText;

        // Find "X-Leg Same Game Parlay +XXX" to get leg count and best odds
        const legMatch = bodyText.match(/(\d+)-Leg Same Game Parlay\s+([+-]\d+)/);
        if (legMatch) {
            result.leg_count = parseInt(legMatch[1]);
            result.best_odds = parseInt(legMatch[2]);
        }

        // Find all odds in the parlay section
        const placeIdx = bodyText.indexOf('Place Bets');
        const bookIdx = bodyText.indexOf('Book', placeIdx);

        if (bookIdx > 0) {
            const section = bodyText.substring(bookIdx, bookIdx + 2000);
            const lines = section.split('\n').map(l => l.trim()).filter(l => l);

            for (const line of lines) {
                if (/^[+-]\d{3}$/.test(line)) {
                    const odds = parseInt(line);
                    result.all_odds.push(odds);
                }
            }
        }

        // Deduplicate all_odds
        const seen = new Set();
        result.all_odds = result.all_odds.filter(o => {
            if (seen.has(o)) return false;
            seen.add(o);
            return true;
        });

        return result;
    }''')

    return result


# ---------------------------------------------------------------------------
# Navigation helpers
# ---------------------------------------------------------------------------

def navigate_to_events(page, sport: str = "NFL") -> bool:
    """Navigate to the Events page and filter by sport."""
    print(f"Navigating to Events ({sport})...")

    # Go directly to events page
    page.goto(f"{PIKKIT_URL}/events", wait_until='networkidle')
    page.wait_for_timeout(2000)

    # Click sport filter - look for text containing the sport name
    clicked = page.evaluate(f'''() => {{
        // Get all clickable elements
        const elements = document.querySelectorAll('button, div, span, p');
        for (const elem of elements) {{
            const text = elem.innerText.trim();
            // Match exact sport name or with icon
            if (text === "{sport}" || text === "🏈 {sport}") {{
                elem.click();
                return "exact: " + text;
            }}
        }}
        // Try partial match
        for (const elem of elements) {{
            const text = elem.innerText.trim();
            if (text.length < 20 && text.toUpperCase().includes("{sport}")) {{
                elem.click();
                return "partial: " + text;
            }}
        }}
        return false;
    }}''')

    if clicked:
        print(f"  Clicked filter: {clicked}")
        page.wait_for_timeout(2000)
    else:
        print(f"  Could not find {sport} filter - checking available filters")
        # List what filters are available
        filters = page.evaluate('''() => {
            const result = [];
            const elements = document.querySelectorAll('button, div[role="button"]');
            elements.forEach(el => {
                const text = el.innerText.trim();
                if (text.length > 0 && text.length < 30) {
                    result.push(text);
                }
            });
            return result.slice(0, 20);
        }''')
        print(f"  Available: {filters}")

    return True


def navigate_to_date(page, target_date: str) -> bool:
    """Navigate to a specific date on the Events page.

    Args:
        target_date: Date string like "Jan 25" or "20" (day number)
    """
    print(f"Navigating to date: {target_date}")

    clicked = page.evaluate(f'''() => {{
        const elements = document.querySelectorAll('button, div, span');
        for (const elem of elements) {{
            const text = elem.innerText.trim();
            if (text.includes("{target_date}")) {{
                elem.click();
                return text;
            }}
        }}
        return false;
    }}''')

    if clicked:
        print(f"  Clicked date: {clicked}")
        page.wait_for_timeout(2000)
        return True
    return False


# ---------------------------------------------------------------------------
# Game finding
# ---------------------------------------------------------------------------

def find_game_on_pikkit(page, away_team: str, home_team: str, sport: str = "NFL") -> str:
    """
    Find a game on Pikkit and return the game URL.

    Args:
        away_team: Away team name (or partial, e.g., "Miami")
        home_team: Home team name (or partial, e.g., "Indiana")
        sport: Sport filter ("NFL", "NCAAF", "NCAAFB", "NBA", "MLB")

    Returns:
        Game URL if found, None otherwise
    """
    # Map common sport names to Pikkit's naming
    sport_map = {
        "NCAAF": "NCAAFB",  # Pikkit uses NCAAFB
        "CFB": "NCAAFB",
        "MLB": "MLB",
    }
    pikkit_sport = sport_map.get(sport.upper(), sport.upper())

    # Go to events page
    page.goto(f"{PIKKIT_URL}/events", wait_until='networkidle')
    page.wait_for_timeout(2000)

    # Click sport filter pill at the top using Playwright locator (more reliable)
    try:
        sport_btn = page.locator(f'text={pikkit_sport}').first
        if sport_btn.count() > 0:
            sport_btn.click()
            page.wait_for_timeout(2000)
    except Exception as e:
        print(f"  Could not click sport filter {pikkit_sport}: {e}")

    # Extract key words from team names for flexible matching
    # "Miami FL" -> "miami", "New England Patriots" -> "patriots", etc.
    def get_search_terms(team_name):
        name = team_name.lower().strip()

        # Remove state abbreviations and common suffixes
        for suffix in [' fl', ' oh', ' tx', ' ca', ' ny', ' state', ' st']:
            name = name.replace(suffix, '')

        words = name.split()
        terms = [name]  # Full cleaned name

        # Add individual words (for matching "Patriots" in "New England Patriots")
        for word in words:
            if len(word) > 3:  # Skip short words like "the", "at", etc.
                terms.append(word)

        # Add last word (often the team nickname)
        if words:
            terms.append(words[-1])

        return list(set(terms))  # Deduplicate

    away_terms = get_search_terms(away_team)
    home_terms = get_search_terms(home_team)

    print(f"  Searching for: {away_terms} or {home_terms}")

    # Check if our teams are on the page using flexible matching
    teams_found = page.evaluate(f'''() => {{
        const text = document.body.innerText.toLowerCase();
        const awayTerms = {away_terms};
        const homeTerms = {home_terms};

        const awayFound = awayTerms.some(term => text.includes(term));
        const homeFound = homeTerms.some(term => text.includes(term));

        return awayFound || homeFound;
    }}''')

    if not teams_found:
        print(f"  Teams not found on events page: {away_team}, {home_team}")
        return None

    # Find and click "More wagers" near our team names
    # Pikkit uses JavaScript routing, so we need to click and get URL after navigation
    try:
        # Find "More wagers" link in a game card with BOTH our teams
        # Key: check Level 1 parent (the game card) has our teams, not higher levels
        clicked = page.evaluate(f'''() => {{
            const awayTerms = {away_terms};
            const homeTerms = {home_terms};

            // Find all "More wagers" elements
            const moreWagersLinks = Array.from(document.querySelectorAll('*'))
                .filter(el => el.innerText && el.innerText.trim() === 'More wagers →');

            for (const link of moreWagersLinks) {{
                // Game card is at grandparent level (parent of parent)
                // Link -> parent (small) -> grandparent (game card ~120-130 chars)
                let container = link.parentElement?.parentElement;
                if (!container) continue;

                const text = (container.innerText || '').toLowerCase();
                const textLen = text.length;

                // Game cards are typically 100-200 chars
                if (textLen < 80 || textLen > 200) continue;

                // Check if this container has BOTH our teams
                const hasAway = awayTerms.some(term => text.includes(term));
                const hasHome = homeTerms.some(term => text.includes(term));

                if (hasAway && hasHome) {{
                    // Found the correct game card
                    link.click();
                    return true;
                }}
            }}

            return false;
        }}''')

        if clicked:
            page.wait_for_timeout(3000)
            current_url = page.url
            if '/event/' in current_url:
                return current_url

        # Alternative: Try clicking on team name directly to navigate
        # This is safer than clicking random "More wagers" links
        for term in home_terms + away_terms:
            if len(term) > 4:  # Skip short terms
                try:
                    team_link = page.locator(f'text={term}').first
                    if team_link.count() > 0:
                        team_link.click()
                        page.wait_for_timeout(3000)
                        if '/event/' in page.url:
                            return page.url
                except Exception as e:
                    print(f"Warning: could not click team link '{term}': {e}")

        # DO NOT fall back to clicking first "More wagers" - this causes wrong game selection
        print(f"  Could not find specific game for {away_team} vs {home_team}")

    except Exception as e:
        print(f"  Error finding game: {e}")

    return None


# ---------------------------------------------------------------------------
# SGP odds builder
# ---------------------------------------------------------------------------

def get_sgp_odds_for_parlay(page, period: str, spread_team: str, spread_value: float,
                            total_value: float, over_under: str, sport: str = "NFL") -> dict:
    """
    Build a parlay on Pikkit using exact line matching via Odds tab.

    Navigates to the Odds tab, expands Spread and Total sections to find
    the exact lines matching the target, then builds the parlay.

    Uses the /betslip API response to get accurate book names and odds.

    Args:
        page: Playwright page on a Pikkit game
        period: "Full", "1st Half", "1st Quarter", etc.
        spread_team: Team name for spread leg (used to identify underdog vs favorite)
        spread_value: Target spread value (e.g., 5.0 for +5.0, -3.0 for favorite)
        total_value: Target total value (e.g., 42.5)
        over_under: "Over" or "Under"
        sport: Sport key for period config lookup (default "NFL")

    Returns:
        Dict with best_odds, all_odds, book_odds, pikkit_spread, pikkit_total
    """
    result = {'best_odds': None, 'all_odds': [], 'book_odds': [], 'pikkit_spread': None, 'pikkit_total': None}

    # Set up API response capture for /betslip endpoint
    # Only keep responses with 2+ leg parlays to avoid using 1-leg data
    betslip_data = [None]  # Use list to allow mutation in closure

    def capture_betslip(response):
        if 'prod-website.pikkit.app/betslip' in response.url:
            try:
                data = response.json()
                parlay = data.get('parlay', {})
                options = parlay.get('options', [])
                # Only keep if it has 2+ legs
                if options and len(options[0].get('outcomes', [])) >= 2:
                    betslip_data[0] = data
            except Exception as e:
                print(f"Warning: could not parse betslip response: {e}")

    page.on('response', capture_betslip)

    # Format the target spread string (e.g., "+5" or "-3.5")
    if spread_value > 0:
        target_spread = f"+{spread_value:g}"
    else:
        target_spread = f"{spread_value:g}"

    # Format target total
    target_total = f"{total_value:g}"
    ou_prefix = 'o' if over_under.lower() == 'over' else 'u'
    is_over = over_under.lower() == 'over'

    # Clear any existing betslip selections before building new parlay
    clear_betslip(page)
    page.wait_for_timeout(500)

    # Click Odds tab and wait for content to actually load
    try:
        odds_tab = page.locator('text=Odds').first
        if odds_tab.count() > 0:
            odds_tab.click()
            # Wait for filter pills to appear (content varies by sport)
            for _ in range(20):  # Up to 10 seconds
                page.wait_for_timeout(500)
                loaded = page.evaluate(r'''() => {
                    const text = document.body.innerText;
                    // "Popular" and "Spread" appear across all sports once odds load
                    return text.includes('Popular') && text.includes('Spread');
                }''')
                if loaded:
                    break
    except Exception as e:
        print(f"    Could not click Odds tab: {e}")
        return result

    # Select period if not Full Game
    # Look up period config for this sport
    sport_config = PERIOD_CONFIGS.get(sport.upper(), PERIOD_CONFIGS.get("NFL"))
    period_to_select_value = sport_config["period_to_select"]
    filter_category_fn = sport_config["filter_category"]

    if period != "Full":
        try:
            # Use the sport-specific filter category (e.g., "Halves", "Quarters")
            filter_cat = filter_category_fn(period)
            if filter_cat:
                page.locator(f'text={filter_cat}').first.click()
                page.wait_for_timeout(1500)

            # Use the <select> dropdown to pick the specific period
            select_val = period_to_select_value.get(period)
            if select_val:
                page.select_option('select', value=select_val)
                page.wait_for_timeout(1500)
        except Exception as e:
            print(f"Warning: could not select period '{period}': {e}")

    # Build a period prefix for matching section headers (e.g., "2nd Quarter Spread")
    # For "Full" period, sections are just "Spread" and "Total"
    # For quarters, sections are "1st Quarter Spread", "2nd Quarter Total", etc.
    section_prefix = period if period != "Full" else ""

    # === STEP 1: Expand Spread section ===
    spread_clicked = False
    try:
        spread_box = page.evaluate(r'''(prefix) => {
            const elements = document.querySelectorAll('*');
            for (const el of elements) {
                const text = el.innerText;
                if (!text || !text.includes('Spread') || !text.includes('Picks')) continue;
                if (text.includes('Moneyline') || text.includes('Total')) continue;
                // Match period-specific section (e.g., "2nd Quarter Spread") or generic "Spread"
                if (prefix && !text.includes(prefix)) continue;
                const rect = el.getBoundingClientRect();
                if (rect.width > 100 && rect.height > 20 && rect.height < 100) {
                    return {x: rect.x, y: rect.y, width: rect.width, height: rect.height};
                }
            }
            return null;
        }''', section_prefix)

        if spread_box:
            page.mouse.click(spread_box['x'] + spread_box['width']/2, spread_box['y'] + spread_box['height']/2)
            page.wait_for_timeout(2000)

            spread_val = target_spread

            # Scroll down incrementally to find the target spread (virtual scrolling)
            # The spread list is long and only visible items are rendered
            found_spread = False
            for scroll_attempt in range(15):  # Try scrolling up to 15 times
                found = page.evaluate(r'''(spreadVal) => {
                    const elements = document.querySelectorAll('*');
                    for (const el of elements) {
                        const text = (el.innerText || '').trim();
                        if (text.length > 20) continue;
                        const firstToken = text.split(/[\s\n]+/)[0];
                        if (firstToken === spreadVal) {
                            el.scrollIntoView({behavior: 'instant', block: 'center'});
                            return true;
                        }
                    }
                    return false;
                }''', spread_val)

                if found:
                    found_spread = True
                    page.wait_for_timeout(500)
                    break

                # Scroll down using mouse wheel in the spread area
                page.mouse.move(spread_box['x'] + spread_box['width']/2, spread_box['y'] + 200)
                page.mouse.wheel(0, 400)  # Scroll down more aggressively
                page.wait_for_timeout(400)

            if not found_spread:
                # Try using keyboard to scroll
                page.keyboard.press('End')
                page.wait_for_timeout(500)
                page.keyboard.press('Home')
                page.wait_for_timeout(500)
                # One more attempt after keyboard scroll
                page.evaluate(r'''(spreadVal) => {
                    const elements = document.querySelectorAll('*');
                    for (const el of elements) {
                        const text = (el.innerText || '').trim();
                        if (text.length > 20) continue;
                        const firstToken = text.split(/[\s\n]+/)[0];
                        if (firstToken === spreadVal) {
                            el.scrollIntoView({behavior: 'instant', block: 'center'});
                            return true;
                        }
                    }
                    return false;
                }''', spread_val)
                page.wait_for_timeout(500)

            # Find and click spread cell with target value and reasonable odds
            cell_pos = page.evaluate(r'''(spreadVal) => {
                const elements = document.querySelectorAll('*');
                const candidates = [];

                for (const el of elements) {
                    const text = (el.innerText || '').trim();
                    if (text.length > 20) continue;

                    // Extract parts: split by whitespace/newlines
                    const parts = text.split(/[\s\n]+/);
                    // Exact match on first token (spread value)
                    if (parts[0] !== spreadVal) continue;
                    if (parts.length < 2) continue;

                    const oddsStr = parts[parts.length - 1];
                    const odds = parseInt(oddsStr);
                    if (isNaN(odds)) continue;

                    const rect = el.getBoundingClientRect();
                    if (rect.width < 20 || rect.width > 500) continue;
                    if (rect.height < 20 || rect.height > 100) continue;

                    candidates.push({
                        el: el,
                        x: rect.x + rect.width/2,
                        y: rect.y + rect.height/2,
                        odds: odds,
                        absOdds: Math.abs(odds)
                    });
                }

                if (candidates.length === 0) return null;

                // Pick odds closest to -110 (typical vig)
                candidates.sort((a, b) => {
                    return Math.abs(a.absOdds - 110) - Math.abs(b.absOdds - 110);
                });

                const best = candidates[0];
                best.el.scrollIntoView({behavior: 'instant', block: 'center'});
                return {x: best.x, y: best.y, odds: best.odds};
            }''', spread_val)
            page.wait_for_timeout(500)

            if cell_pos:
                # Re-get position after scroll and click
                cell_pos = page.evaluate(r'''(spreadVal) => {
                    const elements = document.querySelectorAll('*');
                    for (const el of elements) {
                        const text = (el.innerText || '').trim();
                        if (text.length > 20) continue;

                        const parts = text.split(/[\s\n]+/);
                        if (parts[0] !== spreadVal) continue;
                        if (parts.length < 2) continue;
                        const odds = parseInt(parts[parts.length - 1]);
                        if (isNaN(odds)) continue;
                        if (Math.abs(Math.abs(odds) - 110) > 100) continue;

                        const rect = el.getBoundingClientRect();
                        if (rect.width < 20 || rect.width > 500) continue;
                        if (rect.height < 20 || rect.height > 100) continue;
                        if (rect.y < 0 || rect.y > 800) continue;

                        return {x: rect.x + rect.width/2, y: rect.y + rect.height/2, odds: odds};
                    }
                    return null;
                }''', spread_val)
                page.wait_for_timeout(300)

            if cell_pos:
                page.mouse.click(cell_pos['x'], cell_pos['y'])
                spread_clicked = True
                result['pikkit_spread'] = target_spread
                page.wait_for_timeout(1000)
    except Exception as e:
        print(f"    Error clicking spread: {e}")

    if not spread_clicked:
        print(f"    Could not find exact spread {target_spread} on Pikkit")
        return result

    # === STEP 2: Expand Total section and click exact total ===
    total_clicked = False
    try:
        # First scroll the Total row into view (it may be far down after Spread expansion)
        # Match period-specific "Total" section but NOT "Team Total"
        page.evaluate(r'''(prefix) => {
            const elements = document.querySelectorAll('*');
            for (const el of elements) {
                const text = el.innerText;
                if (!text || !text.includes('Total') || !text.includes('Picks')) continue;
                if (text.includes('Team') || text.includes('Moneyline') || text.includes('Spread')) continue;
                if (prefix && !text.includes(prefix)) continue;
                const rect = el.getBoundingClientRect();
                if (rect.width > 100 && rect.height > 20 && rect.height < 100) {
                    el.scrollIntoView({behavior: 'instant', block: 'center'});
                    return true;
                }
            }
            return false;
        }''', section_prefix)
        page.wait_for_timeout(500)

        # Get Total row position and click to expand
        total_box = page.evaluate(r'''(prefix) => {
            const elements = document.querySelectorAll('*');
            for (const el of elements) {
                const text = el.innerText;
                if (!text || !text.includes('Total') || !text.includes('Picks')) continue;
                if (text.includes('Team') || text.includes('Moneyline') || text.includes('Spread')) continue;
                if (prefix && !text.includes(prefix)) continue;
                const rect = el.getBoundingClientRect();
                if (rect.width > 100 && rect.height > 20 && rect.height < 100) {
                    return {x: rect.x + rect.width/2, y: rect.y + rect.height/2};
                }
            }
            return null;
        }''', section_prefix)

        if total_box:
            page.mouse.click(total_box['x'], total_box['y'])
            page.wait_for_timeout(2000)

            # Scroll to the exact total value row
            page.evaluate(f'''() => {{
                const elements = document.querySelectorAll('*');
                for (const el of elements) {{
                    const text = (el.innerText || '').trim();
                    if (text === '{target_total}') {{
                        el.scrollIntoView({{behavior: 'instant', block: 'center'}});
                        return true;
                    }}
                }}
                return false;
            }}''')
            page.wait_for_timeout(500)

            # Find the Over or Under odds cell for this total
            # Over is first odds column (+101), Under is second (-107)
            # Need to find the correct odds value on the same row as the total
            ou_cell = page.evaluate(f'''() => {{
                const targetTotal = '{target_total}';
                const isOver = {str(is_over).lower()};
                const elements = document.querySelectorAll('*');

                // Find the total value element first
                let totalEl = null;
                for (const el of elements) {{
                    const text = (el.innerText || '').trim();
                    if (text === targetTotal) {{
                        const rect = el.getBoundingClientRect();
                        if (rect.width > 0 && rect.height > 0) {{
                            totalEl = {{x: rect.x, y: rect.y, centerY: rect.y + rect.height/2}};
                            break;
                        }}
                    }}
                }}

                if (!totalEl) return {{error: 'total not found'}};

                // Now find the Over and Under cells on the same row
                // They should be to the right of the total value
                const oddsCells = [];
                for (const el of elements) {{
                    const text = (el.innerText || '').trim();
                    // Match odds pattern like +101, -107, +225, etc.
                    if (/^[+-]\\d{{2,3}}$/.test(text)) {{
                        const rect = el.getBoundingClientRect();
                        // Check if on same row (within 40px vertically) and to the right
                        if (Math.abs(rect.y - totalEl.y) < 40 && rect.x > totalEl.x) {{
                            oddsCells.push({{
                                x: rect.x + rect.width/2,
                                y: rect.y + rect.height/2,
                                xPos: rect.x,
                                text: text
                            }});
                        }}
                    }}
                }}

                // Sort by x position (left to right)
                oddsCells.sort((a, b) => a.xPos - b.xPos);

                // Over is first (index 0), Under is second (index 1)
                const targetIdx = isOver ? 0 : 1;
                if (oddsCells.length > targetIdx) {{
                    return oddsCells[targetIdx];
                }}
                return {{error: 'odds cells not found', count: oddsCells.length}};
            }}''')

            if ou_cell and ou_cell.get('x'):
                page.mouse.click(ou_cell['x'], ou_cell['y'])
                total_clicked = True
                result['pikkit_total'] = f"{ou_prefix}{target_total}"
                page.wait_for_timeout(1000)
    except Exception as e:
        print(f"    Error clicking total: {e}")

    if not total_clicked:
        print(f"    Could not find exact total {over_under} {total_value} on Pikkit")
        return result

    # === STEP 3: Click Parlay tab and extract odds ===
    try:
        parlay_tab = page.get_by_text("Parlay", exact=True)
        if parlay_tab.count() > 0:
            parlay_tab.first.click()
    except Exception:
        click_parlay_tab(page)

    # Wait for parlay odds to load
    for _ in range(10):
        page.wait_for_timeout(500)
        parlay_info = page.evaluate(r'''() => {
            const text = document.body.innerText;
            const match = text.match(/(\d+)-Leg Same Game Parlay\s+([+-]\d+)/);
            if (match) {
                return {legs: parseInt(match[1]), odds: parseInt(match[2])};
            }
            return null;
        }''')
        if parlay_info and parlay_info.get('legs') == 2:
            break

    page.wait_for_timeout(500)

    # Verify we have exactly 2 legs
    leg_check = page.evaluate(r'''() => {
        const text = document.body.innerText;
        const match = text.match(/(\d+)-Leg Same Game Parlay/);
        return match ? parseInt(match[1]) : 0;
    }''')

    if leg_check != 2:
        print(f"    Warning: Expected 2-leg parlay, got {leg_check} legs")
        # Remove the response handler before returning
        page.remove_listener('response', capture_betslip)
        return result

    # Extract parlay odds - prefer API data, fallback to DOM scraping
    if betslip_data[0]:
        odds = extract_parlay_odds_from_api(betslip_data[0])
    else:
        print(f"    No API data, using DOM scraping")
        odds = extract_parlay_odds(page)

    result['best_odds'] = odds.get('best_odds')
    result['all_odds'] = odds.get('all_odds', [])
    result['book_odds'] = odds.get('book_odds', [])

    # Remove the response handler
    page.remove_listener('response', capture_betslip)

    return result
