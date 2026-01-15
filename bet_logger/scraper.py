#!/usr/bin/env python3
"""
Wagerzon Bet History Scraper
Scrapes bet history from Wagerzon and uploads to Google Sheets
"""

import re
import os
import sys
from datetime import datetime
from bs4 import BeautifulSoup
from dotenv import load_dotenv

# Only import playwright when needed (not for --test mode)
if "--test" not in sys.argv:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sync_playwright = None

load_dotenv()

# Configuration
WAGERZON_URL = os.getenv("WAGERZON_URL", "https://www.wagerzon.com")  # Update if different
WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")


def calculate_american_odds(risk: float, win: float) -> int:
    """Convert risk/win to American odds format."""
    if win >= risk:
        # Underdog: positive odds
        return int(round((win / risk) * 100))
    else:
        # Favorite: negative odds
        return int(round(-(risk / win) * 100))


def calculate_decimal_odds_from_american(american_odds: int) -> float:
    """Convert American odds to decimal odds format."""
    if american_odds > 0:
        return (american_odds / 100) + 1
    else:
        return (100 / abs(american_odds)) + 1


def parse_bet_type(description_strong: str, description_spans: list) -> str:
    """Extract bet type from description."""
    strong_lower = description_strong.lower()

    if "prop builder" in strong_lower:
        # Check if it's a parlay or straight from the span content
        if description_spans:
            first_span = description_spans[0].lower()
            if "[parlay:" in first_span:
                return "Parlay"
            elif "[straight:" in first_span:
                return "Straight"
        return "Prop"
    elif "parlay" in strong_lower:
        return "Parlay"
    elif "straight" in strong_lower:
        return "Straight"
    elif "prop:" in strong_lower:
        return "Prop"
    elif "trifecta" in strong_lower:
        return "Trifecta"
    elif "superfecta" in strong_lower:
        return "Superfecta"
    else:
        return ""


def parse_sport(description: str) -> str:
    """Try to determine sport from description."""
    description_upper = description.upper()

    # NFL indicators
    if any(x in description_upper for x in ["NFL", "PANTHERS", "RAMS", "BEARS", "PACKERS",
                                             "CHIEFS", "EAGLES", "COWBOYS", "49ERS", "BILLS",
                                             "RAVENS", "LIONS", "BRONCOS", "JETS", "GIANTS"]):
        return "NFL"

    # NHL indicators
    if any(x in description_upper for x in ["SHARKS", "BLUE JACKETS", "BRUINS", "RANGERS",
                                             "MAPLE LEAFS", "CANADIENS", "PENGUINS", "FLYERS",
                                             "BLACKHAWKS", "RED WINGS", "AVALANCHE", "KNIGHTS"]):
        return "NHL"

    # NBA indicators
    if any(x in description_upper for x in ["LAKERS", "CELTICS", "WARRIORS", "NETS", "KNICKS",
                                             "HEAT", "BUCKS", "SUNS", "MAVERICKS", "CLIPPERS"]):
        return "NBA"

    # NCAAF indicators
    if any(x in description_upper for x in ["NCAAF", "COLLEGE FOOTBALL"]):
        return "NCAAF"

    # NCAAM/NCAAB indicators
    if any(x in description_upper for x in ["NCAAM", "NCAAB", "COLLEGE BASKETBALL"]):
        return "NCAAM"

    return ""


def parse_line(description_spans: list) -> str:
    """Extract line/spread from description if present."""
    for span in description_spans:
        # Look for patterns like "Over 0.5", "+3.5", "-7.5", "Under 45.5"
        over_under = re.search(r'(Over|Under)\s+([\d.]+)', span, re.IGNORECASE)
        if over_under:
            return f"{over_under.group(1)} {over_under.group(2)}"

        # Look for spread patterns with Points - avoid matching standalone odds like +350
        spread_with_points = re.search(r'([+-][\d.]+)\s*Points?', span, re.IGNORECASE)
        if spread_with_points:
            return spread_with_points.group(1)

    return ""


def clean_rtf_encoding(text: str) -> str:
    """Clean RTF encoding artifacts from text."""
    # Replace common RTF escape sequences
    replacements = {
        r"\'bd": "½",
        r"\'bc": "¼",
        r"\'be": "¾",
        "&amp;": "&",
    }
    for pattern, replacement in replacements.items():
        text = text.replace(pattern, replacement)
    return text


def clean_description(description_strong: str, description_spans: list) -> str:
    """Create clean bet description."""
    parts = []

    # Add the main title if it's informative
    if description_strong and "prop builder" not in description_strong.lower():
        parts.append(description_strong)

    # Add span content, cleaning up the bracket codes
    for span in description_spans:
        # Remove the [Parlay: xxx] or [Straight: xxx] prefixes
        cleaned = re.sub(r'\[(Parlay|Straight):\s*\d+\]\s*', '', span)
        # Remove [xxx] number codes
        cleaned = re.sub(r'\[\d+\]\s*', '', cleaned)
        cleaned = cleaned.strip()
        if cleaned:
            parts.append(cleaned)

    result = " | ".join(parts) if parts else description_strong
    return clean_rtf_encoding(result)


def parse_bets_from_html(html_content: str) -> list:
    """Parse bet data from HTML table."""
    soup = BeautifulSoup(html_content, 'html.parser')
    table = soup.find('table', id='wager-results-table')

    if not table:
        print("Could not find wager-results-table")
        return []

    bets = []
    rows = table.find_all('tr', class_='history-GameRow')

    for row in rows:
        cells = row.find_all('td')
        if len(cells) < 6:
            continue

        # Skip transfer rows
        description_div = cells[2].find('div')
        if description_div:
            strong = description_div.find('strong')
            if strong and 'transfer' in strong.get_text().lower():
                continue

        try:
            # Parse Placed (date)
            placed_text = cells[0].get_text(separator=' ').strip()
            date_match = re.search(r'(\d{2}/\d{2}/\d{4})', placed_text)
            if date_match:
                date_str = date_match.group(1)
                date_obj = datetime.strptime(date_str, '%m/%d/%Y')
                formatted_date = date_obj.strftime('%-m/%-d/%y')  # Format as "1/6/26" to match existing data
            else:
                formatted_date = ""

            # Parse Description
            description_strong = ""
            description_spans = []

            if description_div:
                strong = description_div.find('strong')
                if strong:
                    description_strong = strong.get_text().strip()

                spans = description_div.find_all('span')
                for span in spans:
                    # Skip result/italic spans
                    if 'italic' in span.get('style', ''):
                        continue
                    span_text = span.get_text(separator=' ').strip()
                    if span_text:
                        description_spans.append(span_text)

            # Parse Risk/Win
            risk_win_text = cells[3].get_text().strip()
            risk_win_match = re.search(r'([\d.]+)\s*/\s*([\d.]+)', risk_win_text)
            if risk_win_match:
                risk = float(risk_win_match.group(1))
                win = float(risk_win_match.group(2))
                american_odds = calculate_american_odds(risk, win)
                decimal_odds = calculate_decimal_odds_from_american(american_odds)
                bet_amount = risk
            else:
                american_odds = None
                decimal_odds = None
                bet_amount = None

            # Parse Result
            result_text = cells[4].get_text().strip().upper()
            if 'WIN' in result_text and 'LOSE' not in result_text:
                result = "win"
            elif 'LOSE' in result_text:
                result = "loss"
            else:
                result = ""

            # Build bet record
            bet = {
                'date': formatted_date,
                'platform': 'Wagerzon',
                'sport': parse_sport(' '.join(description_spans)),
                'description': clean_description(description_strong, description_spans),
                'bet_type': parse_bet_type(description_strong, description_spans),
                'line': parse_line(description_spans),
                'odds': american_odds,
                'bet_amount': bet_amount,
                'dec': decimal_odds,
                'result': result
            }

            bets.append(bet)

        except Exception as e:
            print(f"Error parsing row: {e}")
            continue

    return bets


def scrape_wagerzon(days_back: int = 7) -> list:
    """
    Log into Wagerzon and scrape bet history.

    Args:
        days_back: Number of days of history to fetch (not implemented yet)

    Returns:
        List of parsed bet dictionaries
    """
    if not WAGERZON_USERNAME or not WAGERZON_PASSWORD:
        raise ValueError("WAGERZON_USERNAME and WAGERZON_PASSWORD must be set in .env file")

    with sync_playwright() as p:
        # Launch browser in headless mode by default
        browser = p.chromium.launch(headless=True)
        context = browser.new_context()
        page = context.new_page()

        print(f"Navigating to {WAGERZON_URL}...")
        page.goto(WAGERZON_URL)

        # Wait for page to load
        page.wait_for_load_state('networkidle')

        # TODO: Update these selectors based on actual Wagerzon login page
        # These are placeholders - you'll need to inspect the actual login form
        print("Attempting login...")

        # Check if already on history page (auto-logged in via cookies)
        current_url = page.url
        print(f"Initial URL: {current_url}")

        if "History" in current_url:
            print("Already on history page (session restored)!")
        else:
            # Try to navigate directly to history page first
            print("Navigating to history page...")
            page.goto("https://backend.wagerzon.com/wager/History.aspx")
            page.wait_for_load_state('networkidle')

            current_url = page.url
            print(f"After navigation: {current_url}")

            # If redirected to login page, try to log in
            if "Login" in current_url or "History" not in current_url:
                print("Login required. Looking for login form...")
                try:
                    # Try common login form selectors
                    username_field = page.locator('input[type="text"], input[name="username"], input[name="txtUsername"]').first
                    username_field.fill(WAGERZON_USERNAME, timeout=10000)

                    password_field = page.locator('input[type="password"], input[name="password"], input[name="txtPassword"]').first
                    password_field.fill(WAGERZON_PASSWORD)

                    login_button = page.locator('input[type="submit"], button[type="submit"], input[value="Login"], button:has-text("Login")').first
                    login_button.click()

                    page.wait_for_load_state('networkidle')
                    print("Login submitted!")

                    # Navigate to history after login
                    import time
                    time.sleep(2)
                    if "History" not in page.url:
                        page.goto("https://backend.wagerzon.com/wager/History.aspx")
                        page.wait_for_load_state('networkidle')

                except Exception as e:
                    print(f"Auto-login failed: {e}")
                    print("Please log in manually in the browser window...")
                    print("You have 60 seconds to log in and navigate to bet history.")
                    try:
                        page.wait_for_url("**/History*", timeout=60000)
                    except:
                        print("Trying to continue anyway...")

        # Debug: print current URL
        print(f"Current URL: {page.url}")

        # Try to change dropdown to show more history (e.g., "Last Week" or "This Month")
        try:
            dropdown = page.locator('select').first
            # Get available options
            options = dropdown.locator('option').all_text_contents()
            print(f"Available time periods: {options}")

            # Try to select a broader time range
            for preferred in ['All', 'This Month', 'Last Month', 'Last Week']:
                if preferred in options:
                    dropdown.select_option(label=preferred)
                    print(f"Selected time period: {preferred}")
                    break

            # Wait for table to reload
            page.wait_for_load_state('networkidle')
            import time
            time.sleep(2)  # Extra wait for data to load
        except Exception as e:
            print(f"Could not change time period: {e}")

        # Wait for actual bet rows to load (not just the table header)
        print("Waiting for bet data to load...")
        try:
            # Wait for either bet rows or "no wagers" message
            page.wait_for_selector('.history-GameRow, :text("No wagers")', timeout=30000)
        except:
            pass  # Continue anyway and see what we get


        # Wait for the table to load
        try:
            page.wait_for_selector('#wager-results-table', timeout=30000)
        except Exception as e:
            print(f"Could not find table: {e}")
            # Save HTML for debugging
            with open("/Users/callancapitolo/bet_logger/debug_page.html", "w") as f:
                f.write(page.content())
            print("Saved page HTML to debug_page.html")
            browser.close()
            return []

        # Get the table HTML
        table_html = page.locator('#wager-results-table').evaluate('el => el.outerHTML')

        # Also get the full container in case we need it
        try:
            container_html = page.locator('.table-responsive').evaluate('el => el.outerHTML')
        except:
            container_html = table_html

        browser.close()

        # Parse the bets
        bets = parse_bets_from_html(container_html)

        return bets


def scrape_from_file(filepath: str) -> list:
    """Parse bets from a saved HTML file (for testing)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Handle RTF wrapper if present
    if content.startswith('{\\rtf'):
        # Extract HTML from RTF
        html_match = re.search(r'<div class="table-responsive">.*</div>', content, re.DOTALL)
        if html_match:
            content = html_match.group(0)

    return parse_bets_from_html(content)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        # Test mode: parse from the saved RTF file
        test_file = "/Users/callancapitolo/Desktop/wagerzon.rtf"
        print(f"Testing parser with {test_file}...")
        bets = scrape_from_file(test_file)

        print(f"\nFound {len(bets)} bets:\n")
        for bet in bets:
            print(f"Date: {bet['date']}")
            print(f"Platform: {bet['platform']}")
            print(f"Sport: {bet['sport']}")
            print(f"Description: {bet['description']}")
            print(f"Bet Type: {bet['bet_type']}")
            print(f"Line: {bet['line']}")
            print(f"Odds: {bet['odds']}")
            print(f"Bet Amount: ${bet['bet_amount']}")
            print(f"Result: {bet['result']}")
            print("-" * 50)
    else:
        # Production mode: scrape from website
        print("Starting Wagerzon scraper...")
        bets = scrape_wagerzon()
        print(f"\nScraped {len(bets)} bets")

        # Import and use sheets module
        from sheets import append_bets_to_sheet
        append_bets_to_sheet(bets)

    sys.exit(0)
