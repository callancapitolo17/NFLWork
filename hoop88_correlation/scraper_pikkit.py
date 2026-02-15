#!/usr/bin/env python3
"""
Pikkit Pro SGP Odds Scraper
Scrapes correlation-adjusted parlay odds from Pikkit Pro.

Pikkit is an odds comparison tool that aggregates lines across sportsbooks.
It shows SGP (Same Game Parlay) odds from multiple books when you build a parlay.

Workflow:
1. Navigate to game page
2. Select period (Full, 1st Half, 1st Quarter)
3. Click spread selection (e.g., +4.5)
4. Click total selection (e.g., o23 or u24)
5. Click "Parlay" tab
6. Extract SGP odds from all listed sportsbooks

Pikkit uses phone number + SMS verification for login.
After first login, session is saved to avoid repeated SMS verification.
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

PIKKIT_URL = os.getenv("PIKKIT_URL", "https://app.pikkit.com")

# Session storage for persistent login
SESSION_FILE = Path(__file__).parent / ".pikkit_session.json"
DEBUG_DIR = Path(__file__).parent


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


def select_period(page, period: str) -> bool:
    """
    Select a period filter on the Pikkit game page.

    Args:
        period: "Full", "1st Half", or "1st Quarter"

    Returns:
        True if period was found and clicked
    """
    clicked = page.evaluate(f'''() => {{
        const buttons = document.querySelectorAll('button, div[role="button"]');
        for (const btn of buttons) {{
            const text = btn.innerText.trim();
            if (text === "{period}") {{
                btn.click();
                return true;
            }}
        }}
        return false;
    }}''')

    if clicked:
        page.wait_for_timeout(1500)
    return clicked


def click_spread(page, team_name: str, spread_value: str) -> bool:
    """
    Click on a spread selection for a team.

    Args:
        team_name: Team name (e.g., "Miami Hurricanes")
        spread_value: Spread value (e.g., "+4.5" or "-3.5")

    Returns:
        True if spread was clicked
    """
    # Find and click the spread button
    clicked = page.evaluate(f'''() => {{
        // Find the team row
        const rows = document.querySelectorAll('div, tr');
        for (const row of rows) {{
            const text = row.innerText;
            if (text.includes("{team_name}")) {{
                // Find spread button in this row
                const buttons = row.querySelectorAll('button, div[class*="odds"]');
                for (const btn of buttons) {{
                    const btnText = btn.innerText.trim();
                    // Match spread like +4.5 or -3.5
                    if (btnText.includes("{spread_value}") ||
                        (btnText.startsWith("{spread_value.split()[0]}") && /[+-]\\d/.test(btnText))) {{
                        btn.click();
                        return true;
                    }}
                }}
            }}
        }}
        return false;
    }}''')

    if clicked:
        page.wait_for_timeout(800)
    return clicked


def click_total(page, over_under: str, total_value: str) -> bool:
    """
    Click on a total (over/under) selection.

    Args:
        over_under: "o" for over, "u" for under
        total_value: Total value (e.g., "23" or "40.5")

    Returns:
        True if total was clicked
    """
    prefix = over_under.lower()[0]  # 'o' or 'u'

    clicked = page.evaluate(f'''() => {{
        const buttons = document.querySelectorAll('button, div[class*="odds"]');
        for (const btn of buttons) {{
            const text = btn.innerText.trim().toLowerCase();
            // Match patterns like o23, u24, o40.5
            if (text.startsWith("{prefix}{total_value}") ||
                text.includes("{prefix}{total_value}")) {{
                btn.click();
                return true;
            }}
        }}
        return false;
    }}''')

    if clicked:
        page.wait_for_timeout(800)
    return clicked


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


def decimal_to_american(decimal_odds: float) -> int:
    """Convert decimal odds to American odds."""
    if decimal_odds >= 2.0:
        return int(round((decimal_odds - 1) * 100))
    else:
        return int(round(-100 / (decimal_odds - 1)))


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
        except Exception:
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


def find_games_on_events_page(page) -> list:
    """Find games on the Events page and extract odds data."""
    print("Looking for games...")

    games = page.evaluate(r'''() => {
        const games = [];
        const pageText = document.body.innerText;

        // Parse games from page text
        // Pattern: Team1, odds, spread, total, Team2, odds, spread, total, time, More wagers
        const lines = pageText.split('\n').map(l => l.trim()).filter(l => l);

        // Find "More wagers" markers which end each game section
        for (let i = 0; i < lines.length; i++) {
            if (lines[i].includes('More wagers')) {
                // Look backwards to find game data
                // Typical pattern:
                // Team1 Name
                // +/-odds (moneyline)
                // +/-spread
                // odds
                // o/u total
                // odds
                // Team2 Name
                // ...similar...
                // time
                // More wagers

                // Find the time (e.g., "1/25 12:00 PM")
                let timeIdx = i - 1;
                while (timeIdx > 0 && !lines[timeIdx].match(/\d+\/\d+\s+\d+:\d+/)) {
                    timeIdx--;
                }

                if (timeIdx > 10) {
                    // Extract chunk of text for this game
                    const chunk = lines.slice(timeIdx - 15, i).join('\n');

                    // Find team names (lines that look like team names)
                    const teamPattern = /^[A-Z][a-z]+\s+[A-Z][a-z]+/;
                    const potentialTeams = lines.slice(timeIdx - 15, timeIdx)
                        .filter(l => teamPattern.test(l) && l.length > 5 && l.length < 40);

                    if (potentialTeams.length >= 2) {
                        const game = {
                            away: potentialTeams[0],
                            home: potentialTeams[1],
                            time: lines[timeIdx],
                            raw: chunk
                        };

                        // Extract odds using regex
                        const spreadMatch = chunk.match(/([+-]\d+\.?\d*)\s+(-?\d+)/g);
                        const totalMatch = chunk.match(/([ou])(\d+\.?\d*)\s+(-?\d+)/gi);

                        if (spreadMatch) game.spreads = spreadMatch;
                        if (totalMatch) game.totals = totalMatch;

                        games.push(game);
                    }
                }
            }
        }

        // Also find "More wagers" links for navigation
        const moreWagersLinks = Array.from(document.querySelectorAll('a')).filter(a =>
            a.innerText.includes('More wagers')
        );

        moreWagersLinks.forEach((link, idx) => {
            if (games[idx]) {
                games[idx].href = link.href;
            }
        });

        return games;
    }''')

    print(f"Found {len(games)} potential games")
    return games

    print(f"Found {len(games)} potential games")
    return games


def search_for_game(page, team1: str, team2: str) -> str:
    """Search for a specific game and return the URL."""
    print(f"Searching for {team1} vs {team2}...")

    # Navigate to search
    page.goto(f"{PIKKIT_URL}/search", wait_until='networkidle')
    page.wait_for_timeout(2000)

    # Type search query
    search_input = page.query_selector('input[type="search"], input[placeholder*="Search"]')
    if search_input:
        search_input.fill(team1)
        page.wait_for_timeout(1500)

        # Look for game result
        game_url = page.evaluate(f'''() => {{
            const results = document.querySelectorAll('a[href*="/events/"]');
            for (const link of results) {{
                const text = link.innerText.toLowerCase();
                if (text.includes("{team1.lower()}") || text.includes("{team2.lower()}")) {{
                    return link.href;
                }}
            }}
            return null;
        }}''')

        return game_url

    return None


def extract_game_odds(page, period: str = "1st Half") -> dict:
    """
    Extract spread and total odds from a game page.

    Args:
        page: Playwright page on a game detail view
        period: "Full Game", "1st Half", "1st Quarter", etc.

    Returns:
        Dict with spread and total odds from various books
    """
    print(f"Extracting {period} odds...")

    # Click appropriate filter (Halves, Quarters, etc.)
    if "Half" in period:
        page.evaluate('''() => {
            const btns = document.querySelectorAll('button, div[role="button"]');
            for (const btn of btns) {
                if (btn.innerText.includes('Halves')) {
                    btn.click();
                    return true;
                }
            }
            return false;
        }''')
        page.wait_for_timeout(1500)
    elif "Quarter" in period:
        page.evaluate('''() => {
            const btns = document.querySelectorAll('button, div[role="button"]');
            for (const btn of btns) {
                if (btn.innerText.includes('Quarters')) {
                    btn.click();
                    return true;
                }
            }
            return false;
        }''')
        page.wait_for_timeout(1500)

    # Extract the odds data
    odds_data = page.evaluate(r'''() => {
        const result = {
            spreads: [],
            totals: [],
            raw_text: document.body.innerText.substring(0, 5000)
        };

        // Find all odds buttons/elements
        const allBtns = Array.from(document.querySelectorAll('[class*="odds"], button'));
        const oddsElements = allBtns.filter(el => {
            const text = el.innerText;
            return text.includes('+') || text.includes('-');
        });

        oddsElements.forEach(elem => {
            const text = elem.innerText.trim();
            // Match odds patterns like -110, +150, etc.
            const oddsMatch = text.match(/([+-]\d+)/);
            if (oddsMatch) {
                // Try to determine if this is spread or total
                const parentText = elem.parentElement?.innerText || '';
                result.spreads.push({
                    odds: oddsMatch[1],
                    context: parentText.substring(0, 50)
                });
            }
        });

        return result;
    }''')

    return odds_data


def get_sgp_odds(page, spread_team: str, spread_line: float, total_line: float, over_under: str) -> dict:
    """
    Try to get SGP odds for a spread + total combination.

    Args:
        page: Playwright page on game detail
        spread_team: Team name for spread bet
        spread_line: Spread value (e.g., +4.5)
        total_line: Total value (e.g., 41.5)
        over_under: "Over" or "Under"

    Returns:
        Dict with SGP odds from different books
    """
    print(f"Getting SGP odds for {spread_team} {spread_line:+.1f} + {over_under} {total_line}...")

    # This would involve:
    # 1. Finding and clicking the spread selection
    # 2. Finding and clicking the total selection
    # 3. Reading the combined parlay odds from the betslip

    # For now, we'll extract available odds and let find_edge.py handle comparison
    sgp_odds = page.evaluate(f'''() => {{
        const result = {{
            found: false,
            odds: {{}},
            debug: []
        }};

        // Look for spread and total sections
        const allText = document.body.innerText;

        // Find spread odds near the team name
        const spreadRegex = new RegExp("{spread_team}.*?([+-]\\d+)", "i");
        const spreadMatch = allText.match(spreadRegex);
        if (spreadMatch) {{
            result.debug.push("Found spread: " + spreadMatch[0]);
        }}

        // Find total odds
        const totalRegex = new RegExp("{over_under}\\s*{total_line}.*?([+-]\\d+)", "i");
        const totalMatch = allText.match(totalRegex);
        if (totalMatch) {{
            result.debug.push("Found total: " + totalMatch[0]);
        }}

        return result;
    }}''')

    return sgp_odds


def explore_game_page(page, game_url: str):
    """Explore a specific game page to understand its structure."""
    print(f"\nNavigating to game: {game_url}")
    if page.url != game_url:
        page.goto(game_url, wait_until='networkidle')
        page.wait_for_timeout(3000)

    # Take screenshot
    page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_game.png"))

    # Get page structure
    structure = page.evaluate('''() => {
        const result = {
            title: document.title,
            url: window.location.href,
            tabs: [],
            filters: [],
            markets: [],
            betslip: null
        };

        // Find tab elements (For You, Odds, Injuries, etc.)
        document.querySelectorAll('[role="tab"], button').forEach(tab => {
            const text = tab.innerText.trim();
            if (text && text.length < 30) {
                result.tabs.push(text);
            }
        });

        // Find filter pills (Popular, Player Props, Halves, etc.)
        document.querySelectorAll('[class*="pill"], [class*="filter"], [class*="chip"]').forEach(pill => {
            const text = pill.innerText.trim();
            if (text && text.length < 30) {
                result.filters.push(text);
            }
        });

        // Find market sections
        document.querySelectorAll('h2, h3, [class*="market"], [class*="heading"]').forEach(heading => {
            const text = heading.innerText.trim();
            if (text && text.length < 50 && !result.markets.includes(text)) {
                result.markets.push(text);
            }
        });

        // Check betslip
        const betslip = document.querySelector('[class*="betslip"], [class*="slip"]');
        if (betslip) {
            result.betslip = betslip.innerText.substring(0, 200);
        }

        return result;
    }''')

    print(f"\nGame page structure:")
    print(f"  Title: {structure['title']}")
    print(f"  Tabs: {structure['tabs'][:10]}")
    print(f"  Filters: {structure['filters'][:10]}")
    print(f"  Markets: {structure['markets'][:10]}")
    print(f"  Betslip: {structure['betslip']}")

    # Save HTML
    html = page.content()
    with open(DEBUG_DIR / "debug_pikkit_game.html", "w") as f:
        f.write(html)
    print(f"\nSaved game page HTML to debug_pikkit_game.html")

    return structure


def extract_game_odds_from_page(page, period: str = "Full") -> dict:
    """
    Extract spread and total odds from the current game page.

    Args:
        page: Playwright page on a game detail
        period: "Full", "1st Half", "1st Quarter"

    Returns:
        Dict with teams, spreads, totals, and odds
    """
    print(f"\nExtracting odds for period: {period}")

    # Click on period filter if not Full
    if period != "Full":
        clicked = page.evaluate(f'''() => {{
            const buttons = document.querySelectorAll('button, div[role="button"]');
            for (const btn of buttons) {{
                if (btn.innerText.trim() === "{period}") {{
                    btn.click();
                    return true;
                }}
            }}
            return false;
        }}''')
        if clicked:
            print(f"  Clicked '{period}' filter")
            page.wait_for_timeout(2000)
        else:
            print(f"  Could not find '{period}' filter")

    page.screenshot(path=str(DEBUG_DIR / f"debug_pikkit_{period.replace(' ', '_')}.png"))

    # Extract odds data from the page
    odds_data = page.evaluate(r'''() => {
        const result = {
            teams: [],
            spreads: [],
            totals: [],
            period: '',
            raw_text: ''
        };

        // Get the main content area
        const content = document.body.innerText;
        result.raw_text = content.substring(0, 3000);

        // Find team names
        const teamElements = document.querySelectorAll('p, span');
        for (const el of teamElements) {
            const text = el.innerText.trim();
            if (text.includes('Patriots') || text.includes('Broncos') ||
                text.includes('Seahawks') || text.includes('Rams')) {
                if (!result.teams.includes(text) && text.length < 50) {
                    result.teams.push(text);
                }
            }
        }

        // Find clickable odds buttons
        const oddsButtons = [];
        document.querySelectorAll('button, div[class*="css-"]').forEach(el => {
            const text = el.innerText.trim();
            // Match odds patterns like -4.5, +6.5, -106, +222, o40.5, u42
            if (/^[+-]?\d+\.?\d*$/.test(text) || /^[ou]\d+/.test(text.toLowerCase())) {
                oddsButtons.push({
                    text: text,
                    rect: el.getBoundingClientRect()
                });
            }
        });
        result.odds_buttons = oddsButtons.slice(0, 20);

        return result;
    }''')

    print(f"  Teams found: {odds_data.get('teams', [])[:4]}")
    print(f"  Odds buttons: {[b['text'] for b in odds_data.get('odds_buttons', [])[:10]]}")

    return odds_data


def test_betslip_parlay(page) -> dict:
    """
    Test if clicking multiple odds creates a parlay in the betslip.

    Returns:
        Dict with betslip contents after adding selections
    """
    print("\nTesting betslip parlay functionality...")

    # Try to click on a spread odds button
    clicked_spread = page.evaluate(r'''() => {
        // Find the spread row (-4.5 or +6.5 type)
        const allElements = document.querySelectorAll('div, button');
        for (const el of allElements) {
            const text = el.innerText.trim();
            // Look for spread values like -4.5, +6.5
            if (/^[+-]\d+\.5$/.test(text)) {
                el.click();
                return text;
            }
        }
        return null;
    }''')

    if clicked_spread:
        print(f"  Clicked spread: {clicked_spread}")
        page.wait_for_timeout(1500)

    # Take screenshot after first click
    page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_after_spread_click.png"))

    # Try to click on a total odds button
    clicked_total = page.evaluate(r'''() => {
        // Find total values like o40.5, u42
        const allElements = document.querySelectorAll('div, button, p');
        for (const el of allElements) {
            const text = el.innerText.trim().toLowerCase();
            if (/^[ou]\d+/.test(text)) {
                el.click();
                return text;
            }
        }
        return null;
    }''')

    if clicked_total:
        print(f"  Clicked total: {clicked_total}")
        page.wait_for_timeout(1500)

    # Take screenshot after second click
    page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_after_total_click.png"))

    # Check the betslip
    betslip_content = page.evaluate(r'''() => {
        const result = {
            text: '',
            selections: 0,
            parlay_odds: null
        };

        // Find betslip area
        const betslip = document.body.innerText;

        // Look for "Place Bets" section
        const placeIdx = betslip.indexOf('Place Bets');
        if (placeIdx > 0) {
            result.text = betslip.substring(placeIdx, placeIdx + 500);
        }

        // Count selections
        const selMatch = result.text.match(/(\d+)\s*(?:selection|bet|pick)/i);
        if (selMatch) {
            result.selections = parseInt(selMatch[1]);
        }

        // Look for parlay odds
        const oddsMatch = result.text.match(/parlay[^+\-]*([+-]\d+)/i);
        if (oddsMatch) {
            result.parlay_odds = oddsMatch[1];
        }

        return result;
    }''')

    print(f"  Betslip selections: {betslip_content.get('selections', 0)}")
    print(f"  Betslip text: {betslip_content.get('text', '')[:200]}")

    return betslip_content


def explore_pikkit_structure(headless: bool = False):
    """
    Explore Pikkit's website structure to understand how to scrape SGP odds.
    Run with --visible to see the browser and manually login if needed.
    """
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)

        # Try to load saved session
        saved_session = load_session()
        if saved_session:
            context = browser.new_context(
                viewport={'width': 1920, 'height': 1080},
                storage_state=saved_session
            )
            print("Loaded saved session")
        else:
            context = browser.new_context(viewport={'width': 1920, 'height': 1080})

        page = context.new_page()

        print(f"Navigating to {PIKKIT_URL}...")
        page.goto(PIKKIT_URL, wait_until='networkidle')
        page.wait_for_timeout(3000)

        # Check if logged in
        logged_in = is_logged_in(page)
        print(f"Logged in: {logged_in}")

        if not logged_in and not headless:
            print("\n" + "="*60)
            print("NOT LOGGED IN - Please login manually")
            print("Browser is open for manual interaction.")
            print("Session will be saved after 120 seconds (2 minutes).")
            print("="*60)
            page.wait_for_timeout(120000)
            save_session(context)
            logged_in = is_logged_in(page)

        if logged_in:
            # Take screenshot of logged-in state
            page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_main.png"))

            # Try to search for a specific game instead of browsing events
            print("\nSearching for Patriots game...")
            page.goto(f"{PIKKIT_URL}/search", wait_until='networkidle')
            page.wait_for_timeout(2000)

            # Find and fill search input
            search_filled = page.evaluate('''() => {
                const inputs = document.querySelectorAll('input');
                for (const input of inputs) {
                    if (input.placeholder?.toLowerCase().includes('search') ||
                        input.type === 'search' ||
                        input.type === 'text') {
                        input.value = 'Patriots';
                        input.dispatchEvent(new Event('input', { bubbles: true }));
                        return true;
                    }
                }
                return false;
            }''')

            if search_filled:
                print("  Typed 'Patriots' in search")
                page.wait_for_timeout(2000)
                page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_search.png"))

                # Look for search results
                results = page.evaluate('''() => {
                    const links = Array.from(document.querySelectorAll('a'));
                    return links
                        .filter(a => a.href.includes('/events/'))
                        .map(a => ({
                            text: a.innerText.substring(0, 100),
                            href: a.href
                        }))
                        .slice(0, 10);
                }''')

                print(f"  Search results: {len(results)}")
                for r in results[:5]:
                    print(f"    {r['text'][:50]} -> {r['href']}")

                # If we found a game, explore it
                if results:
                    print(f"\nExploring first result...")
                    explore_game_page(page, results[0]['href'])

            # Also try the Events page with NFL filter
            navigate_to_events(page, sport="NFL")
            page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_events.png"))

            # Try navigating to Jan 25 for NFL games
            navigate_to_date(page, "25")
            page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_jan25.png"))

            # Find games
            games = find_games_on_events_page(page)
            print(f"\nGames found on events page:")
            for g in games[:10]:
                print(f"  {g.get('away', '?')} @ {g.get('home', '?')}")
                print(f"    Time: {g.get('time', '?')}")
                print(f"    Spreads: {g.get('spreads', [])}")
                print(f"    Totals: {g.get('totals', [])}")
                if g.get('href'):
                    print(f"    Link: {g.get('href')}")

            # Click "More wagers" to get to game detail page
            # Try using Playwright's native locator and click
            print(f"\n{'='*60}")
            print("Clicking 'More wagers' to explore game detail...")
            print(f"{'='*60}")

            try:
                # Use Playwright's text locator
                more_wagers = page.get_by_text("More wagers", exact=False).first
                if more_wagers:
                    print("Found 'More wagers' element, clicking...")
                    more_wagers.click()
                    page.wait_for_timeout(3000)

                    new_url = page.url
                    print(f"Current URL after click: {new_url}")
                    page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_game.png"))

                    # Explore the game page if we navigated there
                    if '/event/' in new_url:
                        explore_game_page(page, new_url)

                        # Extract odds for different periods
                        print("\n" + "="*60)
                        print("EXTRACTING ODDS FROM PIKKIT GAME PAGE")
                        print("="*60)
                        for period in ["Full", "1st Half", "1st Quarter"]:
                            odds = extract_game_odds_from_page(page, period)

                        # Test betslip parlay functionality
                        print("\n" + "="*60)
                        print("TESTING BETSLIP PARLAY")
                        print("="*60)
                        betslip = test_betslip_parlay(page)
                    else:
                        print("Did not navigate to game page, trying team name click...")
                        # Try clicking on the team name instead
                        patriots = page.get_by_text("New England Patriots").first
                        if patriots:
                            patriots.click()
                            page.wait_for_timeout(3000)
                            new_url = page.url
                            print(f"URL after Patriots click: {new_url}")
                            page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_game.png"))
                            if '/event/' in new_url:
                                explore_game_page(page, new_url)
                                # Extract odds and test betslip
                                for period in ["Full", "1st Half"]:
                                    extract_game_odds_from_page(page, period)
                                test_betslip_parlay(page)
            except Exception as e:
                print(f"Error clicking: {e}")

                # Save HTML for debugging
                html = page.content()
                with open(DEBUG_DIR / "debug_pikkit_click_error.html", "w") as f:
                    f.write(html)

            # If we found a game, explore it
            if games and games[0].get('href'):
                structure = explore_game_page(page, games[0]['href'])

                # Try to extract odds
                odds = extract_game_odds(page, "1st Half")
                print(f"\nExtracted odds data:")
                print(f"  Spreads found: {len(odds.get('spreads', []))}")
        else:
            print("\nNot logged in. Run with --visible to login manually.")
            page.screenshot(path=str(DEBUG_DIR / "debug_pikkit_not_logged_in.png"))

        # Save HTML for detailed analysis
        html = page.content()
        with open(DEBUG_DIR / "debug_pikkit.html", "w") as f:
            f.write(html)
        print("\nSaved current page HTML to debug_pikkit.html")

        browser.close()


def find_game_on_pikkit(page, away_team: str, home_team: str, sport: str = "NFL") -> str:
    """
    Find a game on Pikkit and return the game URL.

    Args:
        away_team: Away team name (or partial, e.g., "Miami")
        home_team: Home team name (or partial, e.g., "Indiana")
        sport: Sport filter ("NFL", "NCAAF", "NCAAFB", "NBA")

    Returns:
        Game URL if found, None otherwise
    """
    # Map common sport names to Pikkit's naming
    sport_map = {
        "NCAAF": "NCAAFB",  # Pikkit uses NCAAFB
        "CFB": "NCAAFB",
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
                except:
                    pass

        # DO NOT fall back to clicking first "More wagers" - this causes wrong game selection
        print(f"  Could not find specific game for {away_team} vs {home_team}")

    except Exception as e:
        print(f"  Error finding game: {e}")

    return None

    # Look for the game by team names
    # Try to find "More wagers" link near the team names
    game_url = page.evaluate(f'''() => {{
        const links = document.querySelectorAll('a');
        const awayLower = "{away_team}".toLowerCase();
        const homeLower = "{home_team}".toLowerCase();

        for (const link of links) {{
            const href = link.href || '';
            const text = link.innerText.toLowerCase();
            const parentText = link.parentElement?.innerText?.toLowerCase() || '';
            const grandparentText = link.parentElement?.parentElement?.innerText?.toLowerCase() || '';

            // Check if this link or its context mentions our teams
            const contextText = text + ' ' + parentText + ' ' + grandparentText;

            // Match partial team names (e.g., "Patriots" matches "New England Patriots")
            const awayParts = awayLower.split(' ');
            const homeParts = homeLower.split(' ');

            const hasAway = awayParts.some(part => part.length > 3 && contextText.includes(part));
            const hasHome = homeParts.some(part => part.length > 3 && contextText.includes(part));

            if ((hasAway || hasHome) && href.includes('/event/')) {{
                return href;
            }}
        }}
        return null;
    }}''')

    return game_url


def get_sgp_odds_for_parlay(page, period: str, spread_team: str, spread_value: float,
                            total_value: float, over_under: str) -> dict:
    """
    Build a parlay on Pikkit using exact line matching via Odds tab.

    Navigates to the Odds tab, expands Spread and Total sections to find
    the exact lines matching hoop88, then builds the parlay.

    Uses the /betslip API response to get accurate book names and odds.

    Args:
        page: Playwright page on a Pikkit game
        period: "Full", "1st Half", "1st Quarter"
        spread_team: Team name for spread leg (used to identify underdog vs favorite)
        spread_value: Target spread value (e.g., 5.0 for +5.0, -3.0 for favorite)
        total_value: Target total value (e.g., 42.5)
        over_under: "Over" or "Under"

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
            except:
                pass

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
            # Wait for filter pills (Popular, Quarters, etc.) to appear
            for _ in range(20):  # Up to 10 seconds
                page.wait_for_timeout(500)
                loaded = page.evaluate(r'''() => {
                    const text = document.body.innerText;
                    return text.includes('Popular') && text.includes('Quarters');
                }''')
                if loaded:
                    break
    except Exception as e:
        print(f"    Could not click Odds tab: {e}")
        return result

    # Select period if not Full Game
    # Pikkit UI: filter pills (Popular, Halves, Quarters) + <select> dropdown for specific period
    # Quarters dropdown: "1st Quarter" = quarter_1, everything else = "Other"
    period_to_select_value = {
        '1st Half': 'half_1', '2nd Half': 'half_2',
        '1st Quarter': 'quarter_1',
        '2nd Quarter': 'other', '3rd Quarter': 'other', '4th Quarter': 'other',
    }
    if period != "Full":
        try:
            if "Half" in period:
                page.locator('text=Halves').first.click()
                page.wait_for_timeout(1500)
            elif "Quarter" in period:
                page.locator('text=Quarters').first.click()
                page.wait_for_timeout(1500)

            # Use the <select> dropdown to pick the specific period
            select_val = period_to_select_value.get(period)
            if select_val:
                page.select_option('select', value=select_val)
                page.wait_for_timeout(1500)
        except Exception:
            pass

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
            over_odds_pattern = r"[+-]\\d{{2,3}}"  # Over odds (first)
            under_odds_pattern = r"-\\d{{2,3}}"  # Under odds typically negative

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


def get_sgp_odds_for_parlay_fallback(page, period: str, spread_team: str, spread_value: float,
                            total_value: float, over_under: str) -> dict:
    """
    Fallback: Build a parlay using the For You tab with closest available lines.
    Used when exact line matching fails.
    """
    result = {'best_odds': None, 'all_odds': [], 'book_odds': [], 'pikkit_spread': None, 'pikkit_total': None}

    # Set up API response capture for /betslip endpoint
    betslip_data = [None]

    def capture_betslip(response):
        if 'prod-website.pikkit.app/betslip' in response.url:
            try:
                betslip_data[0] = response.json()
            except:
                pass

    page.on('response', capture_betslip)

    # Clear any existing selections
    clear_betslip(page)
    page.wait_for_timeout(500)

    # Select period using Playwright locator
    if period != "Full":
        try:
            period_btn = page.get_by_text(period, exact=True)
            if period_btn.count() > 0:
                period_btn.first.click()
                page.wait_for_timeout(1500)
        except Exception:
            pass

    # Get all available lines on the page
    available_lines = page.evaluate(r'''() => {
        const spreads = [];
        const totals = [];

        document.querySelectorAll('p').forEach(p => {
            const text = p.innerText.trim();
            if (/^[+-]\d+\.?\d*$/.test(text)) {
                const val = parseFloat(text);
                if (Math.abs(val) <= 50) {
                    spreads.push({text: text, value: val});
                }
            }
            if (/^[ou]\d+\.?\d*$/i.test(text)) {
                const val = parseFloat(text.substring(1));
                totals.push({text: text, value: val, type: text[0].toLowerCase()});
            }
        });

        return {spreads, totals};
    }''')

    # Find closest spread
    matching_spreads = [s for s in available_lines['spreads']
                       if (s['value'] > 0) == (spread_value > 0)]
    if not matching_spreads:
        return result

    closest_spread = min(matching_spreads, key=lambda s: abs(s['value'] - spread_value))
    spread_str = closest_spread['text']

    # Find closest total
    ou_type = 'o' if over_under.lower() == 'over' else 'u'
    matching_totals = [t for t in available_lines['totals'] if t['type'] == ou_type]
    if not matching_totals:
        return result

    closest_total = min(matching_totals, key=lambda t: abs(t['value'] - total_value))
    total_str = closest_total['text']

    # Click spread and total
    spread_clicked = page.evaluate(f'''() => {{
        const pElements = document.querySelectorAll('p');
        for (const p of pElements) {{
            if (p.innerText.trim() === "{spread_str}") {{
                const cell = p.parentElement;
                if (cell && cell.tagName === 'DIV') {{
                    cell.click();
                    return true;
                }}
            }}
        }}
        return false;
    }}''')

    if not spread_clicked:
        return result
    page.wait_for_timeout(800)

    total_clicked = page.evaluate(f'''() => {{
        const pElements = document.querySelectorAll('p');
        for (const p of pElements) {{
            if (p.innerText.trim().toLowerCase() === "{total_str.lower()}") {{
                const cell = p.parentElement;
                if (cell && cell.tagName === 'DIV') {{
                    cell.click();
                    return true;
                }}
            }}
        }}
        return false;
    }}''')

    if not total_clicked:
        print(f"    Could not click total: {total_str}")
        return {'best_odds': None, 'all_odds': [], 'book_odds': []}

    page.wait_for_timeout(800)

    # Click Parlay tab
    try:
        parlay_tab = page.get_by_text("Parlay", exact=True)
        if parlay_tab.count() > 0:
            parlay_tab.first.click()
    except Exception:
        click_parlay_tab(page)

    # Wait for parlay odds to load (watch for "2-Leg Same Game Parlay" specifically)
    for _ in range(10):  # Try up to 10 times (5 seconds)
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

    page.wait_for_timeout(500)  # Extra buffer

    # Verify we have exactly 2 legs
    leg_check = page.evaluate(r'''() => {
        const text = document.body.innerText;
        const match = text.match(/(\d+)-Leg Same Game Parlay/);
        return match ? parseInt(match[1]) : 0;
    }''')

    if leg_check != 2:
        print(f"    Warning: Expected 2-leg parlay, got {leg_check} legs")
        page.remove_listener('response', capture_betslip)
        return {'best_odds': None, 'all_odds': [], 'book_odds': [], 'pikkit_spread': spread_str, 'pikkit_total': total_str}

    # Extract parlay odds - prefer API data, fallback to DOM scraping
    if betslip_data[0]:
        odds = extract_parlay_odds_from_api(betslip_data[0])
    else:
        odds = extract_parlay_odds(page)

    # Add the actual Pikkit lines used
    odds['pikkit_spread'] = spread_str
    odds['pikkit_total'] = total_str

    # Remove the response handler
    page.remove_listener('response', capture_betslip)

    return odds


def scrape_pikkit_sgp(games: list, headless: bool = True) -> list:
    """
    Get SGP odds from Pikkit for a list of games.

    Args:
        games: List of game dicts from hoop88 scraper
        headless: Run browser in headless mode

    Returns:
        Games list with pikkit_sgp_odds added
    """
    if not games:
        return []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)

        # Try to load saved session
        saved_session = load_session()
        if saved_session:
            context = browser.new_context(
                viewport={'width': 1920, 'height': 1080},
                storage_state=saved_session
            )
            print("Loaded saved Pikkit session")
        else:
            context = browser.new_context(viewport={'width': 1920, 'height': 1080})
            print("No saved session - run with --visible to login")
            browser.close()
            return games

        page = context.new_page()

        print(f"Navigating to {PIKKIT_URL}...")
        page.goto(PIKKIT_URL, wait_until='networkidle')
        page.wait_for_timeout(3000)

        if not is_logged_in(page):
            print("Not logged in - run 'python scraper_pikkit.py --visible' first")
            browser.close()
            return games

        # Process each game
        for game in games:
            away_team = game.get('away_team', '')
            home_team = game.get('home_team', '')
            league = game.get('league', 'NFL')
            period = game.get('period', 'FG')

            # Map period short names
            period_map = {
                'FG': 'Full', '1H': '1st Half', '2H': '2nd Half',
                '1Q': '1st Quarter', '2Q': '2nd Quarter',
                '3Q': '3rd Quarter', '4Q': '4th Quarter',
            }
            period_full = period_map.get(period, period)

            print(f"\nProcessing: {away_team} @ {home_team} ({league} {period})")

            # Find the game on Pikkit
            sport = "NCAAF" if "NCAA" in league or league == "NCAAF" else "NFL"
            game_url = find_game_on_pikkit(page, away_team, home_team, sport)

            if not game_url:
                print(f"  Could not find game on Pikkit")
                continue

            page.goto(game_url, wait_until='networkidle')
            page.wait_for_timeout(2000)

            # Get SGP odds for underdog + under
            underdog = game.get('underdog', '')
            underdog_spread = game.get('underdog_spread', 0)
            total = game.get('total', 0)

            print(f"  Getting SGP for {underdog} {underdog_spread:+.1f} + Under {total}...")
            sgp_uu = get_sgp_odds_for_parlay(
                page, period_full, underdog, underdog_spread, total, "Under"
            )
            game['pikkit_underdog_under'] = sgp_uu
            if sgp_uu.get('best_odds'):
                print(f"    Best SGP odds: {sgp_uu['best_odds']:+d}")

            # Get SGP odds for favorite + over
            favorite = game.get('favorite', '')
            favorite_spread = game.get('favorite_spread', 0)

            print(f"  Getting SGP for {favorite} {favorite_spread:+.1f} + Over {total}...")
            sgp_fo = get_sgp_odds_for_parlay(
                page, period_full, favorite, favorite_spread, total, "Over"
            )
            game['pikkit_favorite_over'] = sgp_fo
            if sgp_fo.get('best_odds'):
                print(f"    Best SGP odds: {sgp_fo['best_odds']:+d}")

        browser.close()

    return games


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pikkit Pro SGP odds scraper')
    parser.add_argument('--visible', action='store_true', help='Show browser window')
    parser.add_argument('--explore', action='store_true', help='Explore Pikkit structure')
    parser.add_argument('--game', type=str, help='Specific game URL to explore')
    parser.add_argument('--test-sgp', action='store_true', help='Test SGP scraping with a sample game')
    args = parser.parse_args()

    print("=" * 60)
    print("PIKKIT PRO SGP SCRAPER")
    print("=" * 60)

    if args.test_sgp:
        # Test SGP scraping with a sample parlay
        print("\nTesting SGP scraping...")

        with sync_playwright() as p:
            browser = p.chromium.launch(headless=not args.visible)

            saved_session = load_session()
            if saved_session:
                context = browser.new_context(
                    viewport={'width': 1920, 'height': 1080},
                    storage_state=saved_session
                )
                print("Loaded saved session")
            else:
                context = browser.new_context(viewport={'width': 1920, 'height': 1080})
                print("No saved session - will need to login")

            page = context.new_page()
            page.goto(PIKKIT_URL, wait_until='networkidle')
            page.wait_for_timeout(2000)

            if not is_logged_in(page):
                print("Not logged in. Running with --visible to allow manual login...")
                if not args.visible:
                    print("Please run with --visible flag to login")
                    browser.close()
                else:
                    print("Please login manually. Session will be saved after 120 seconds.")
                    page.wait_for_timeout(120000)
                    save_session(context)
            else:
                print("Logged in!")

                # Find an NBA game (try common teams)
                game_url = None
                nba_teams = [
                    ("Lakers", "Celtics"),
                    ("Warriors", "Nuggets"),
                    ("Knicks", "Heat"),
                    ("Suns", "Mavericks"),
                    ("Bucks", "76ers"),
                ]

                for away, home in nba_teams:
                    game_url = find_game_on_pikkit(page, away, home, "NBA")
                    if game_url:
                        print(f"Found NBA game: {away} vs {home}")
                        break

                # Fallback: just navigate to NBA events and click first game
                if not game_url:
                    print("Trying to find any NBA game...")
                    page.goto(f"{PIKKIT_URL}/events", wait_until='networkidle')
                    page.wait_for_timeout(2000)
                    try:
                        nba_btn = page.locator('text=NBA').first
                        if nba_btn.count() > 0:
                            nba_btn.click()
                            page.wait_for_timeout(2000)
                            # Click first "More wagers" link
                            more_wagers = page.get_by_text("More wagers", exact=False).first
                            if more_wagers.count() > 0:
                                more_wagers.click()
                                page.wait_for_timeout(3000)
                                if '/event/' in page.url:
                                    game_url = page.url
                                    print(f"Found NBA game via events page: {game_url}")
                    except Exception as e:
                        print(f"Error finding NBA game: {e}")

                if game_url:
                    print(f"Found game: {game_url}")
                    if page.url != game_url:
                        page.goto(game_url, wait_until='networkidle')
                        page.wait_for_timeout(2000)

                    # Extract team names and available lines from the page
                    page_info = page.evaluate(r'''() => {
                        const text = document.body.innerText;
                        const spreads = [];
                        const totals = [];

                        // Find spread values
                        const spreadMatches = text.match(/[+-]\d+\.?\d*/g) || [];
                        spreadMatches.forEach(s => {
                            const val = parseFloat(s);
                            if (Math.abs(val) <= 20 && Math.abs(val) >= 1) {
                                spreads.push(val);
                            }
                        });

                        // Find total values (NBA totals are typically 200-250)
                        const totalMatches = text.match(/[ou](\d{3}\.?\d*)/gi) || [];
                        totalMatches.forEach(t => {
                            totals.push(t);
                        });

                        return {spreads: [...new Set(spreads)].slice(0, 5), totals: [...new Set(totals)].slice(0, 5)};
                    }''')
                    print(f"  Available spreads: {page_info.get('spreads', [])}")
                    print(f"  Available totals: {page_info.get('totals', [])}")

                    # Test getting SGP odds for Full Game spread + under
                    # Use first available spread and a reasonable NBA total
                    test_spread = page_info.get('spreads', [5.5])[0] if page_info.get('spreads') else 5.5
                    test_total = 220.5  # Typical NBA total

                    print(f"\nTesting Full Game parlay: spread {test_spread:+.1f} + Under {test_total}...")
                    sgp = get_sgp_odds_for_parlay(page, "Full", "test", test_spread, test_total, "Under")
                    print(f"  Result: {sgp}")

                    page.screenshot(path=str(DEBUG_DIR / "debug_sgp_test.png"))
                    print(f"\nScreenshot saved to debug_sgp_test.png")
                else:
                    print("Could not find any NBA games")

            if args.visible:
                print("\nBrowser open - press Ctrl+C to exit")
                try:
                    page.wait_for_timeout(300000)
                except KeyboardInterrupt:
                    pass

            browser.close()

    elif args.game:
        # Explore a specific game
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=not args.visible)
            saved_session = load_session()
            if saved_session:
                context = browser.new_context(
                    viewport={'width': 1920, 'height': 1080},
                    storage_state=saved_session
                )
            else:
                context = browser.new_context(viewport={'width': 1920, 'height': 1080})
            page = context.new_page()
            explore_game_page(page, args.game)
            if not args.visible:
                browser.close()
            else:
                print("\nBrowser open - press Ctrl+C to exit")
                page.wait_for_timeout(300000)
    else:
        explore_pikkit_structure(headless=not args.visible)
