#!/usr/bin/env python3
"""
BFA Gaming Bet History Scraper
Scrapes bet history from https://bfagaming.com and uploads to Google Sheets
"""

import os
import sys
import re
from datetime import datetime, timedelta
from bs4 import BeautifulSoup
from dotenv import load_dotenv

# Playwright import
try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

from sheets import append_bets_to_sheet

load_dotenv()

# Configuration
BFA_URL = os.getenv("BFA_URL", "https://bfagaming.com")
BFA_USERNAME = os.getenv("BFA_USERNAME")
BFA_PASSWORD = os.getenv("BFA_PASSWORD")
BFA_HISTORY_URL = "https://bfagaming.com/my-bets"


def calculate_american_odds_from_line(line_text: str) -> int:
    """Extract American odds from line text like '-110' or '+120'."""
    match = re.search(r'([+-]\d+)', line_text)
    if match:
        return int(match.group(1))
    return -110  # Default juice


def calculate_decimal_odds_from_american(american_odds: int) -> float:
    """Convert American odds to decimal odds format."""
    if american_odds > 0:
        return (american_odds / 100) + 1
    elif american_odds < 0:
        return (100 / abs(american_odds)) + 1
    return 1.0


def parse_status(status_text: str) -> str:
    """Extract result from status text."""
    status_upper = status_text.upper().strip()
    if 'WIN' in status_upper or 'WON' in status_upper:
        return 'W'
    elif 'LOST' in status_upper or 'LOSE' in status_upper:
        return 'L'
    elif 'PUSH' in status_upper or 'TIE' in status_upper or 'CANCELLED' in status_upper:
        return 'X'
    return ''


def parse_sport(description: str) -> str:
    """Determine sport from bet description text."""
    text_upper = description.upper()

    # NFL teams
    nfl_teams = [
        'BEARS', 'LIONS', 'PACKERS', 'VIKINGS',  # NFC North
        'COWBOYS', 'EAGLES', 'GIANTS', 'COMMANDERS',  # NFC East
        'FALCONS', 'PANTHERS', 'SAINTS', 'BUCCANEERS',  # NFC South
        'CARDINALS', 'RAMS', 'SEAHAWKS', '49ERS',  # NFC West
        'BENGALS', 'BROWNS', 'RAVENS', 'STEELERS',  # AFC North
        'BILLS', 'DOLPHINS', 'PATRIOTS', 'JETS',  # AFC East
        'COLTS', 'JAGUARS', 'TEXANS', 'TITANS',  # AFC South
        'BRONCOS', 'CHARGERS', 'CHIEFS', 'RAIDERS'  # AFC West
    ]

    # NBA teams
    nba_teams = [
        'CELTICS', 'NETS', 'KNICKS', '76ERS', 'RAPTORS',  # Atlantic
        'BULLS', 'CAVALIERS', 'PISTONS', 'PACERS', 'BUCKS',  # Central
        'HAWKS', 'HORNETS', 'HEAT', 'MAGIC', 'WIZARDS',  # Southeast
        'NUGGETS', 'TIMBERWOLVES', 'THUNDER', 'BLAZERS', 'JAZZ',  # Northwest
        'WARRIORS', 'CLIPPERS', 'LAKERS', 'SUNS', 'KINGS',  # Pacific
        'MAVERICKS', 'ROCKETS', 'GRIZZLIES', 'PELICANS', 'SPURS'  # Southwest
    ]

    # NHL teams
    nhl_teams = [
        'BRUINS', 'SABRES', 'RED WINGS', 'PANTHERS', 'CANADIENS',
        'SENATORS', 'LIGHTNING', 'MAPLE LEAFS', 'HURRICANES', 'BLUE JACKETS',
        'DEVILS', 'ISLANDERS', 'RANGERS', 'FLYERS', 'PENGUINS',
        'CAPITALS', 'BLACKHAWKS', 'AVALANCHE', 'STARS', 'WILD',
        'PREDATORS', 'BLUES', 'JETS', 'DUCKS', 'COYOTES',
        'FLAMES', 'OILERS', 'KINGS', 'SHARKS', 'KRAKEN', 'GOLDEN KNIGHTS'
    ]

    # Check NFL first (most common based on your bets)
    for team in nfl_teams:
        if team in text_upper:
            return 'NFL'

    # Check NBA
    for team in nba_teams:
        if team in text_upper:
            return 'NBA'

    # Check NHL
    for team in nhl_teams:
        if team in text_upper:
            return 'NHL'

    # Check for explicit sport mentions
    if 'NFL' in text_upper or 'FOOTBALL' in text_upper:
        return 'NFL'
    elif 'NBA' in text_upper or 'BASKETBALL' in text_upper:
        return 'NBA'
    elif 'NHL' in text_upper or 'HOCKEY' in text_upper:
        return 'NHL'
    elif 'MLB' in text_upper or 'BASEBALL' in text_upper:
        return 'MLB'
    elif 'NCAAF' in text_upper or 'COLLEGE FOOTBALL' in text_upper:
        return 'NCAAF'
    elif 'NCAAM' in text_upper or 'COLLEGE BASKETBALL' in text_upper:
        return 'NCAAM'

    return ''


def parse_bet_type(description: str) -> str:
    """Extract bet type from description header."""
    desc_upper = description.upper()

    if 'TEASER' in desc_upper:
        return 'Teaser'
    elif 'PARLAY' in desc_upper:
        return 'Parlay'
    elif 'STRAIGHT' in desc_upper:
        return 'Straight'
    elif 'PROP' in desc_upper:
        return 'Prop'

    # Check if it's a single bet (no multi-leg indicators)
    if 'TEAM' not in desc_upper:
        return 'Straight'

    return 'Parlay'


def parse_date(date_str: str) -> str:
    """Convert BFA date string to MM/DD/YY format."""
    # Input format: "Jan 10, 2:23PM PST" or "Jan 12, 9:51AM PST"
    try:
        # Remove timezone
        date_part = re.sub(r'\s*(PST|EST|CST|MST|PDT|EDT|CDT|MDT)\s*$', '', date_str.strip())

        # Try parsing with year
        try:
            dt = datetime.strptime(date_part, "%b %d, %Y %I:%M%p")
        except ValueError:
            # No year in date string - prepend current year to avoid deprecation warning
            now = datetime.now()
            date_with_year = f"{date_part}, {now.year}"
            dt = datetime.strptime(date_with_year, "%b %d, %I:%M%p, %Y")
            if dt > now + timedelta(days=1):  # If more than 1 day in future, use last year
                dt = dt.replace(year=now.year - 1)

        return dt.strftime("%-m/%-d/%y")
    except Exception as e:
        print(f"  Warning: Could not parse date '{date_str}': {e}")
        return date_str


def parse_risk_win(risk_win_text: str) -> tuple:
    """Parse risk/win text into separate values."""
    # Format: "200.00/360.00"
    try:
        parts = risk_win_text.replace('$', '').replace(',', '').split('/')
        if len(parts) == 2:
            risk = float(parts[0].strip())
            win = float(parts[1].strip())
            return risk, win
    except Exception as e:
        print(f"  Warning: Could not parse risk/win '{risk_win_text}': {e}")
    return 0.0, 0.0


def parse_leg(leg_text: str) -> dict:
    """
    Parse a single bet leg from BFA format.
    Example: "[378] CHICAGO BEARS +8-110 (B+6)"
    Returns: {'team': 'CHICAGO BEARS', 'line': '+8', 'odds': -110, 'game_id': '378'}
    """
    leg_text = leg_text.strip()

    # Extract game ID if present: [378]
    game_id = ''
    game_match = re.match(r'\[(\d+)\]\s*', leg_text)
    if game_match:
        game_id = game_match.group(1)
        leg_text = leg_text[game_match.end():]

    # Remove teaser points at end: (B+6)
    leg_text = re.sub(r'\s*\([^)]+\)\s*$', '', leg_text)

    # Extract odds from end: -110 or +120
    odds = -110
    odds_match = re.search(r'([+-]\d+)\s*$', leg_text)
    if odds_match:
        odds = int(odds_match.group(1))
        leg_text = leg_text[:odds_match.start()].strip()

    # Extract line/spread: +8, -3.5, O 45.5, U 42
    line = ''
    # Look for spread pattern
    spread_match = re.search(r'\s+([+-]?\d+\.?\d*)\s*$', leg_text)
    if spread_match:
        line = spread_match.group(1)
        if not line.startswith('+') and not line.startswith('-'):
            line = '+' + line if float(line) > 0 else line
        leg_text = leg_text[:spread_match.start()].strip()

    # Over/Under pattern
    ou_match = re.search(r'\s+([OU])\s*(\d+\.?\d*)\s*$', leg_text, re.IGNORECASE)
    if ou_match:
        direction = 'Over' if ou_match.group(1).upper() == 'O' else 'Under'
        line = f"{direction} {ou_match.group(2)}"
        leg_text = leg_text[:ou_match.start()].strip()

    team = leg_text.strip()

    return {
        'team': team,
        'line': line,
        'odds': odds,
        'game_id': game_id
    }


def clean_description(description_lines: list) -> str:
    """Create clean bet description from parsed legs."""
    # First line is the bet type header, skip it
    if not description_lines:
        return ''

    legs = []
    for line in description_lines[1:]:  # Skip header line
        line = line.strip()
        if not line:
            continue

        leg = parse_leg(line)
        if leg['team']:
            leg_desc = leg['team']
            if leg['line']:
                leg_desc += f" {leg['line']}"
            legs.append(leg_desc)

    return ' | '.join(legs) if legs else description_lines[0]


def parse_bets_from_html(html_content: str) -> list:
    """Parse bet data from BFA Gaming HTML table."""
    soup = BeautifulSoup(html_content, 'html.parser')

    bets = []

    # Find the table - BFA uses a standard table structure
    table = soup.find('table')
    if not table:
        print("Could not find table element")
        return []

    # Get all rows
    rows = table.find_all('tr')
    print(f"Found {len(rows)} table rows")

    # Skip header row
    for row in rows[1:]:
        try:
            cells = row.find_all('td')
            if len(cells) < 6:
                continue

            # Column mapping from screenshot:
            # 0: Ticket (ID)
            # 1: Description
            # 2: Risk/Win
            # 3: W/L (profit/loss amount)
            # 4: Result
            # 5: Placed Time
            # 6: Settled Time (optional)

            ticket_text = cells[0].get_text(strip=True)
            description_text = cells[1].get_text(separator='\n', strip=True)
            risk_win_text = cells[2].get_text(strip=True)
            wl_amount = cells[3].get_text(strip=True)
            result_text = cells[4].get_text(strip=True)
            placed_time = cells[5].get_text(strip=True)

            # Skip non-bet rows (transfers, deposits, etc.)
            if 'Transfer' in description_text or 'Deposit' in description_text:
                print(f"  Skipping non-bet row: {description_text[:50]}...")
                continue

            # Skip rows without proper bet data
            if not risk_win_text or '/' not in risk_win_text:
                continue

            # Parse description lines
            desc_lines = [line.strip() for line in description_text.split('\n') if line.strip()]
            if not desc_lines:
                continue

            # Get bet type from first line
            bet_type = parse_bet_type(desc_lines[0])

            # Build clean description
            description = clean_description(desc_lines)
            if not description:
                description = description_text.replace('\n', ' | ')

            # Parse risk/win
            risk, win = parse_risk_win(risk_win_text)

            # Parse date
            date = parse_date(placed_time)

            # Parse result
            result = parse_status(result_text)

            # Detect sport
            sport = parse_sport(description_text)

            # Calculate odds from risk/win ratio
            if risk > 0 and win > 0:
                american_odds = int(round((win / risk) * 100))
            else:
                american_odds = 0

            decimal_odds = calculate_decimal_odds_from_american(american_odds)

            # Extract first leg's line for the 'line' column
            first_line = ''
            if len(desc_lines) > 1:
                first_leg = parse_leg(desc_lines[1])
                first_line = first_leg.get('line', '')

            bet = {
                'date': date,
                'platform': 'BFA',
                'sport': sport,
                'description': description,
                'bet_type': bet_type,
                'line': first_line,
                'odds': american_odds,
                'bet_amount': risk,
                'dec': decimal_odds,
                'result': 'win' if result == 'W' else 'loss' if result == 'L' else 'push' if result == 'X' else ''
            }

            bets.append(bet)
            print(f"  ✓ {date} - {bet_type} - {description[:60]}... - ${risk:.2f} - {result_text}")

        except Exception as e:
            print(f"  ✗ Error parsing bet row: {e}")
            import traceback
            traceback.print_exc()
            continue

    return bets


def scrape_bfa(weeks_back: int = 1, headless: bool = True) -> list:
    """
    Log into BFA Gaming and scrape bet history.

    Args:
        weeks_back: Number of weeks back to fetch (0=current week, 1=last week, etc.)
        headless: Run browser in headless mode (no visible window)

    Returns:
        List of parsed bet dictionaries
    """
    if not BFA_USERNAME or not BFA_PASSWORD:
        raise ValueError("BFA_USERNAME and BFA_PASSWORD must be set in .env file")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {BFA_HISTORY_URL}...")
        page.goto(BFA_HISTORY_URL, wait_until='networkidle', timeout=60000)
        page.wait_for_timeout(3000)  # Wait for Blazor app to initialize

        # Check if we need to login by looking for "You must log in" message
        print("Checking login status...")
        must_login = page.locator('text=You must log in').count() > 0
        login_button = page.locator('button:has-text("Log In"), button:has-text("Log in")')

        if must_login or login_button.count() > 0:
            print("Login required. Clicking login button...")
            try:
                # Click the Log In button in the header - this redirects to Keycloak auth
                header_login = page.locator('.header-auth-container button:has-text("Log In")')
                if header_login.count() > 0:
                    header_login.click()
                else:
                    login_button.first.click()

                # Wait for redirect to Keycloak auth page
                print("Waiting for Keycloak login page...")
                page.wait_for_url('**/realms/**', timeout=15000)
                page.wait_for_timeout(2000)

                # Keycloak uses specific IDs for the login form
                print("Filling Keycloak login form...")

                # Fill username using Keycloak's #username field
                username_field = page.locator('#username')
                if username_field.count() > 0:
                    username_field.fill(BFA_USERNAME)
                    print("✅ Filled username")
                else:
                    # Fallback to name attribute
                    page.fill('input[name="username"]', BFA_USERNAME)
                    print("✅ Filled username (fallback)")

                page.wait_for_timeout(500)

                # Fill password using Keycloak's #password field
                password_field = page.locator('#password')
                if password_field.count() > 0:
                    password_field.fill(BFA_PASSWORD)
                    print("✅ Filled password")
                else:
                    page.fill('input[name="password"]', BFA_PASSWORD)
                    print("✅ Filled password (fallback)")

                page.wait_for_timeout(500)

                # Click the Keycloak login button
                login_submit = page.locator('#kc-login')
                if login_submit.count() > 0:
                    login_submit.click()
                    print("✅ Clicked Keycloak login button")
                else:
                    page.click('input[type="submit"]')
                    print("✅ Clicked submit button (fallback)")

                # Wait for redirect back to BFA after successful login
                print("Waiting for redirect back to BFA...")
                page.wait_for_url('**/bfagaming.com/**', timeout=30000)
                page.wait_for_load_state('networkidle', timeout=30000)
                page.wait_for_timeout(3000)

                print("✅ Login successful! Redirected back to BFA")

            except Exception as e:
                print(f"Auto-login failed: {e}")
                page.screenshot(path="debug_bfa_login_error.png")
                if not headless:
                    print("Please log in manually (60 seconds)...")
                    page.wait_for_timeout(60000)
                else:
                    raise RuntimeError(f"Login failed in headless mode: {e}")
        else:
            print("Already logged in")

        # Click on "Settled" tab to see bet history
        print("Clicking 'Settled' tab...")
        try:
            settled_tab = page.locator('.scroller-item:has-text("Settled"), div:has-text("Settled"):not(:has(*))')
            if settled_tab.count() > 0:
                settled_tab.first.click()
                page.wait_for_load_state('networkidle', timeout=30000)
                page.wait_for_timeout(2000)
                print("✅ Clicked Settled tab")
            else:
                print("Settled tab not found")
        except Exception as e:
            print(f"Could not click Settled tab: {e}")

        # Set date filter - BFA uses HTML5 date inputs requiring ISO format (YYYY-MM-DD)
        # weeks_back: 0 = current week (Monday to today), 1 = last week (Mon-Sun), 2 = 2 weeks ago, etc.
        week_label = "current week" if weeks_back == 0 else f"{weeks_back} week(s) back"
        print(f"Setting date filter for {week_label}...")
        try:
            today = datetime.now()

            if weeks_back == 0:
                # Current week: this Monday to today
                days_since_monday = today.weekday()  # Monday=0
                start_date = today - timedelta(days=days_since_monday)
                end_date = today
            else:
                # Previous week(s): Monday to Sunday
                # Find last Sunday (end of most recent complete week)
                days_since_sunday = (today.weekday() + 1) % 7  # Monday=0, Sunday=6
                if days_since_sunday == 0:
                    days_since_sunday = 7  # If today is Sunday, go back to previous Sunday
                last_sunday = today - timedelta(days=days_since_sunday)

                # Go back additional weeks if needed
                end_date = last_sunday - timedelta(weeks=weeks_back - 1)
                start_date = end_date - timedelta(days=6)

            # HTML5 date inputs require ISO format: YYYY-MM-DD
            start_str = start_date.strftime("%Y-%m-%d")
            end_str = end_date.strftime("%Y-%m-%d")

            # Fill the from-date input
            from_date = page.locator('#from-date')
            if from_date.count() > 0:
                from_date.fill(start_str)
                print(f"✅ Set from-date: {start_str} ({start_date.strftime('%A')})")
            else:
                print("from-date input not found")

            page.wait_for_timeout(500)

            # Fill the to-date input
            to_date = page.locator('#to-date')
            if to_date.count() > 0:
                to_date.fill(end_str)
                print(f"✅ Set to-date: {end_str} ({end_date.strftime('%A')})")
            else:
                print("to-date input not found")

            page.wait_for_timeout(500)

            # Click Apply Filter button
            apply_btn = page.locator('.apply-filter-btn, button:has-text("Apply Filter")')
            if apply_btn.count() > 0:
                apply_btn.first.click()
                print("✅ Clicked Apply Filter")
                page.wait_for_load_state('networkidle', timeout=30000)
                page.wait_for_timeout(3000)  # Wait for table to update
            else:
                print("Apply Filter button not found")

        except Exception as e:
            print(f"Could not set date filter: {e}")

        # Wait for bet table to load
        print("Waiting for bet table to load...")
        try:
            page.wait_for_selector('table', timeout=15000)
            page.wait_for_timeout(2000)
            print("✅ Table loaded")
        except:
            print("Warning: Table not found, checking for data...")

        # Get the page HTML
        page_html = page.content()

        # Save for debugging
        with open('debug_bfa_page.html', 'w', encoding='utf-8') as f:
            f.write(page_html)
        print("Saved page HTML to debug_bfa_page.html")

        browser.close()
        return parse_bets_from_html(page_html)


def scrape_from_file(filepath: str) -> list:
    """Parse bets from a saved HTML file (for testing)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    return parse_bets_from_html(content)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from BFA Gaming')
    parser.add_argument('--test', action='store_true', help='Parse from saved HTML file instead of live scrape')
    parser.add_argument('--weeks', type=int, default=1, help='Weeks back to fetch (0=current week, 1=last week, default: 1)')
    parser.add_argument('--visible', action='store_true', help='Show browser window (default is headless)')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    args = parser.parse_args()

    if args.test:
        print("TEST MODE: Parsing from saved HTML file")
        test_file = "debug_bfa_page.html"

        if not os.path.exists(test_file):
            print(f"Error: {test_file} not found")
            print("Run the scraper first to generate this file, or provide your own HTML")
            sys.exit(1)

        bets = scrape_from_file(test_file)
        print(f"\nParsed {len(bets)} bets:")
        for i, bet in enumerate(bets, 1):
            print(f"\n{i}. {bet['date']} - {bet['bet_type']}")
            print(f"   Sport: {bet['sport']}")
            print(f"   Description: {bet['description']}")
            print(f"   Line: {bet['line']}")
            print(f"   Odds: {bet['odds']:+d} (Decimal: {bet['dec']:.2f})")
            print(f"   Amount: ${bet['bet_amount']:.2f}")
            print(f"   Result: {bet['result']}")
    else:
        print("=" * 60)
        print("BFA GAMING BET HISTORY SCRAPER")
        print("=" * 60)

        try:
            bets = scrape_bfa(weeks_back=args.weeks, headless=not args.visible)

            print(f"\n{'=' * 60}")
            print(f"Successfully scraped {len(bets)} bets from BFA Gaming")
            print(f"{'=' * 60}\n")

            if bets and not args.dry_run:
                print("Uploading to Google Sheets...")
                result = append_bets_to_sheet(bets)

                if result['status'] == 'success':
                    print(f"\n✅ SUCCESS! Added {result['rows_added']} new bets to sheet")
                    print(f"   Rows {result['start_row']} to {result['end_row']}")
                elif result['status'] == 'skipped':
                    print(f"\n⚠️  {result['message']}")
                else:
                    print(f"\n❌ Error uploading to sheets: {result.get('message', 'Unknown error')}")
            elif args.dry_run:
                print("Dry run - skipping upload to Google Sheets")
            else:
                print("No bets found to upload")

        except Exception as e:
            print(f"\n❌ Error: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    sys.exit(0)
