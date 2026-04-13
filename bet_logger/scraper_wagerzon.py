#!/usr/bin/env python3
"""
Wagerzon Bet History Scraper
Scrapes bet history from Wagerzon and uploads to Google Sheets.

Auth: ASP.NET form POST (same approach as wagerzon_odds/scraper_v2.py).
No browser/Playwright required — pure HTTP via requests.Session.
"""

import re
import os
import sys
from datetime import datetime
from bs4 import BeautifulSoup
from dotenv import load_dotenv
import requests
from utils import (
    calculate_american_odds,
    calculate_decimal_odds_from_american,
    parse_sport,
)

load_dotenv()

# Configuration
WAGERZON_BASE_URL = "https://backend.wagerzon.com"
WAGERZON_HISTORY_URL = f"{WAGERZON_BASE_URL}/wager/History.aspx"
WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")


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


def login(session: requests.Session):
    """Login to Wagerzon via ASP.NET form POST.

    Same auth approach as wagerzon_odds/scraper_v2.py:
    1. GET the login page to capture __VIEWSTATE and other hidden fields
    2. POST credentials with those hidden fields
    3. Session cookie (ASP.NET_SessionId) maintains auth for subsequent requests
    """
    resp = session.get(WAGERZON_BASE_URL, timeout=15)
    resp.raise_for_status()

    # If already redirected to a logged-in page, session is still valid
    if "History" in resp.url or "NewSchedule" in resp.url:
        print("Already authenticated")
        return

    html = resp.text

    # Extract ASP.NET hidden fields from the login form.
    # ASP.NET uses these to validate form submissions — they're generated
    # server-side and must be posted back exactly as received.
    fields = {}
    for name in ["__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"]:
        match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if match:
            fields[name] = match.group(1)

    if "__VIEWSTATE" not in fields:
        raise RuntimeError("Could not find __VIEWSTATE on login page — page structure may have changed")

    fields["Account"] = WAGERZON_USERNAME
    fields["Password"] = WAGERZON_PASSWORD
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
    resp.raise_for_status()
    print("Logged in successfully")


def scrape_wagerzon(days_back: int = 7) -> list:
    """
    Log into Wagerzon via HTTP and scrape bet history.

    Uses requests.Session instead of a browser — faster, more reliable,
    and doesn't break on headless click timeouts.

    Returns:
        List of parsed bet dictionaries
    """
    if not WAGERZON_USERNAME or not WAGERZON_PASSWORD:
        raise ValueError("WAGERZON_USERNAME and WAGERZON_PASSWORD must be set in .env file")

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    })

    # Step 1: Authenticate
    print("Logging in to Wagerzon...")
    login(session)

    # Step 2: Fetch the bet history page HTML
    print("Fetching bet history page...")
    resp = session.get(WAGERZON_HISTORY_URL, timeout=15)
    resp.raise_for_status()

    html = resp.text

    # Check if we got redirected back to login (auth failed silently)
    if "BtnSubmit" in html and "Account" in html and "wager-results-table" not in html:
        raise RuntimeError("Auth session invalid — redirected back to login page")

    # Step 3: Parse the HTML using the existing parser
    bets = parse_bets_from_html(html)

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
        try:
            bets = scrape_wagerzon()
        except Exception as e:
            print(f"\n❌ Error: {e}")
            sys.exit(1)

        print(f"\nScraped {len(bets)} bets")

        if not bets:
            print("No bets found — exiting with error")
            sys.exit(1)

        # Import and use sheets module
        from sheets import append_bets_to_sheet
        append_bets_to_sheet(bets)
