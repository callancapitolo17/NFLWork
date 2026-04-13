#!/usr/bin/env python3
"""
Hoop88 Bet History Scraper
Scrapes bet history from Hoop88 and uploads to Google Sheets.

Auth: JWT via /cloud/api/System/authenticateCustomer (same approach as hoop88_odds/scraper.py).
No browser/Playwright required — pure HTTP via requests.Session.
"""

import os
import sys
import re
from datetime import datetime
from bs4 import BeautifulSoup
from dotenv import load_dotenv
import requests
from utils import (
    calculate_american_odds,
    calculate_american_odds_from_line,
    calculate_decimal_odds_from_american,
    parse_status,
    parse_sport,
    parse_risk_win,
)

from sheets import append_bets_to_sheet

load_dotenv()

# Configuration
HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")
API_BASE = f"{HOOP88_URL}/cloud/api"


def parse_bet_type(type_label: str) -> str:
    """Extract bet type from type label."""
    type_lower = type_label.lower()
    
    if 'parlay' in type_lower:
        return 'Parlay'
    elif 'straight' in type_lower:
        return 'Straight'
    elif 'teaser' in type_lower:
        return 'Teaser'
    elif 'prop' in type_lower:
        return 'Prop'
    
    return ''


def parse_line(line_selected: str, period: str) -> str:
    """Extract line/spread from line text."""
    # Examples: "-½  +130", "O 7½  +105", "U 19  -110", "+1  -110"
    
    # Over/Under
    over_under = re.search(r'([OU])\s+([\d½¼¾.]+)', line_selected)
    if over_under:
        direction = 'Over' if over_under.group(1) == 'O' else 'Under'
        value = over_under.group(2).replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
        return f"{direction} {value}"
    
    # Spread
    spread = re.search(r'([+-][\d½¼¾.]+)', line_selected)
    if spread:
        value = spread.group(1).replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
        return value
    
    return ''


def clean_description(team_name: str, line_selected: str, period: str) -> str:
    """Create clean bet description."""
    # Clean up the line text to remove odds
    line_clean = re.sub(r'\s*[+-]\d+\s*$', '', line_selected).strip()
    
    if period:
        period = period.strip()
        return f"{team_name} {line_clean} - {period}"
    return f"{team_name} {line_clean}"


def parse_date(date_str: str) -> str:
    """Convert date string to MM/DD/YY format."""
    # Input format: "Jan 12, 2026 4:30 PM - (PST)"
    try:
        # Remove timezone part
        date_part = date_str.split(' - ')[0]
        # Parse date
        dt = datetime.strptime(date_part, "%b %d, %Y %I:%M %p")
        # Format as MM/DD/YY
        return dt.strftime("%-m/%-d/%y")
    except:
        return date_str


def parse_contest_bet(row) -> dict:
    """Parse a type-C contest/props bet which has a different HTML structure."""
    try:
        # Contest bets have the description in span.font-normal
        font_normal = row.find('span', class_='font-normal')
        if not font_normal:
            return None

        description = font_normal.get_text(strip=True)

        # The header-bet is inside content-descripcion for type-C bets
        content_desc = row.find('div', class_='content-descripcion')
        if not content_desc:
            return None

        header_bet = content_desc.find('div', class_='header-bet')
        if not header_bet:
            return None

        # Extract key information
        risk_win_text = ''
        accepted_date = ''
        status_text = ''

        for ln_div in header_bet.find_all('div', class_='ln'):
            d_l = ln_div.find('div', class_='d-l')
            d_r = ln_div.find('div', class_='d-r')

            if d_l and d_r:
                label = d_l.get_text(strip=True)
                value = d_r.get_text(strip=True)

                if label == 'Risk/Win':
                    risk_win_text = value
                elif label == 'Accepted':
                    accepted_date = value
                elif label == 'Status':
                    status_text = value

        risk, win = parse_risk_win(risk_win_text)
        date = parse_date(accepted_date)

        # Determine result from status text or wager-status class
        final_result = None
        if 'lost' in status_text.lower():
            final_result = 'loss'
        elif 'win' in status_text.lower() or 'won' in status_text.lower():
            final_result = 'win'
        elif 'push' in status_text.lower() or 'cancel' in status_text.lower():
            final_result = 'push'

        # Fallback to wager-status class in the amount cell
        if not final_result:
            result_cell = row.find('td', class_='to-risk')
            if result_cell:
                status_span = result_cell.find('span', class_=re.compile(r'wager-status-'))
                if status_span:
                    for cls in status_span.get('class', []):
                        if 'wager-status-W' in cls:
                            final_result = 'win'
                            break
                        elif 'wager-status-L' in cls:
                            final_result = 'loss'
                            break
                        elif 'wager-status-X' in cls:
                            final_result = 'push'
                            break

        # Detect sport from description
        sport = parse_sport(description, '')

        # Calculate odds from risk/win
        american_odds = calculate_american_odds(risk, win)
        decimal_odds = calculate_decimal_odds_from_american(american_odds)

        return {
            'date': date,
            'platform': 'Hoop88',
            'sport': sport,
            'description': description,
            'bet_type': 'Contest',
            'line': '',
            'odds': american_odds,
            'bet_amount': risk,
            'dec': decimal_odds,
            'result': final_result or ''
        }

    except Exception as e:
        print(f"  ✗ Error parsing contest bet: {e}")
        return None


def parse_bets_from_html(html_content: str) -> list:
    """Parse bet data from Hoop88 HTML."""
    soup = BeautifulSoup(html_content, 'html.parser')
    tbody = soup.find('tbody')
    
    if not tbody:
        print("Could not find tbody element")
        return []
    
    bets = []
    rows = tbody.find_all('tr', {'data-ticket': True, 'data-type': True})
    
    print(f"Found {len(rows)} bet rows")
    
    for row in rows:
        try:
            ticket_num = row.get('data-ticket', '')
            bet_data_type = row.get('data-type', '')

            # Handle type-C bets (contest/props bets) differently
            if bet_data_type == 'C':
                bet = parse_contest_bet(row)
                if bet:
                    bets.append(bet)
                    print(f"  ✓ {bet['date']} - {bet['bet_type']} - {bet['description'][:60]}... - ${bet['bet_amount']:.2f}")
                continue

            # Find the type label
            type_label_elem = row.find('label', class_='type-magic')
            type_label = type_label_elem.get_text(strip=True) if type_label_elem else ''
            bet_type = parse_bet_type(type_label)

            # Find header bet info
            header_bet = row.find('div', class_='header-bet')
            if not header_bet:
                continue
            
            # Extract key information from header
            risk_win_text = ''
            accepted_date = ''
            note = ''
            game = ''
            selection = ''
            
            for ln_div in header_bet.find_all('div', class_='ln'):
                d_l = ln_div.find('div', class_='d-l')
                d_r = ln_div.find('div', class_='d-r')
                
                if d_l and d_r:
                    label = d_l.get_text(strip=True)
                    value = d_r.get_text(strip=True)
                    
                    if label == 'Risk/Win':
                        risk_win_text = value
                    elif label == 'Accepted':
                        accepted_date = value
                    elif label == 'Note':
                        note = value
                    elif label == 'Game':
                        game = value
                    elif label == 'Selection':
                        selection = value
            
            risk, win = parse_risk_win(risk_win_text)
            date = parse_date(accepted_date)
            sport = parse_sport(selection, note)
            
            # Parse individual legs of parlay or straight bet
            wager_simples = row.find_all('div', class_='wager-simple')
            
            if len(wager_simples) == 0:
                continue
            
            # For parlays with multiple legs, we'll create one combined description
            leg_descriptions = []
            leg_data = []  # Track odds and result for each leg

            for wager_simple in wager_simples:
                # Extract team name
                choosen = wager_simple.find('span', class_='choosen')
                team_name = choosen.get_text(strip=True) if choosen else ''

                # Extract line and odds
                line_selected = wager_simple.find('span', class_='line-selected')
                line_text = line_selected.get_text(strip=True) if line_selected else ''

                # Extract period
                period_desc = wager_simple.find('span', class_='period-description-s')
                period = period_desc.get_text(strip=True) if period_desc else ''

                # Extract status
                status_span = wager_simple.find('span', class_=re.compile(r'status-'))
                leg_result = status_span.get_text(strip=True) if status_span else ''

                # Parse odds from line
                odds = calculate_american_odds_from_line(line_text)
                decimal_odds = calculate_decimal_odds_from_american(odds)

                leg_data.append({
                    'result': leg_result,
                    'american_odds': odds,
                    'decimal_odds': decimal_odds
                })

                # Build leg description
                line_value = parse_line(line_text, period)
                leg_desc = clean_description(team_name, line_text, period)
                leg_descriptions.append(leg_desc)

            # Extract results for determining outcome
            leg_results = [leg['result'] for leg in leg_data]

            # Determine parlay result from leg results:
            # - If ANY leg loses -> parlay loses
            # - If ALL legs push -> parlay pushes
            # - If some legs win and others push (no losses) -> parlay wins (at reduced odds)
            # - If ALL legs win -> parlay wins
            if 'L' in leg_results:
                final_result = 'loss'
            elif all(r == 'X' for r in leg_results):
                final_result = 'push'
            elif all(r in ('W', 'X') for r in leg_results) and 'W' in leg_results:
                final_result = 'win'
            else:
                final_result = None
            
            # Combine all legs into description
            description = ' | '.join(leg_descriptions)

            # Calculate odds - for parlays with pushed legs, recalculate from winning legs only
            has_push = 'X' in leg_results
            if bet_type == 'Parlay' and has_push and final_result == 'win':
                # Multiply decimal odds of only the winning legs (exclude pushed legs)
                combined_decimal = 1.0
                for leg in leg_data:
                    if leg['result'] == 'W':
                        combined_decimal *= leg['decimal_odds']
                decimal_odds = combined_decimal
                # Convert back to American odds
                if decimal_odds >= 2.0:
                    american_odds = int(round((decimal_odds - 1) * 100))
                else:
                    american_odds = int(round(-100 / (decimal_odds - 1)))
            elif bet_type == 'Parlay' and risk > 0:
                # Calculate from potential win (original parlay odds)
                american_odds = int(round((win / risk) * 100))
            else:
                # For straight bets, use the single leg odds
                american_odds = leg_data[0]['american_odds'] if leg_data else 0
            
            # Determine final result for the entire bet from the wager-status class
            result_cell = row.find('td', class_='to-risk')
            if result_cell:
                # Look for the wager-status span with class like wager-status-W, wager-status-L, wager-status-X
                status_span = result_cell.find('span', class_=re.compile(r'wager-status-'))
                if status_span:
                    status_class = status_span.get('class', [])
                    # Find the class that contains the status
                    for cls in status_class:
                        if 'wager-status-W' in cls:
                            final_result = 'win'
                            break
                        elif 'wager-status-L' in cls:
                            final_result = 'loss'
                            break
                        elif 'wager-status-X' in cls:
                            final_result = 'push'
                            break
                
                # Fallback: check text content if class detection didn't work
                if not final_result:
                    result_text = result_cell.get_text(strip=True)
                    if result_text.startswith('+'):
                        final_result = 'win'
                    elif result_text.startswith('-') and result_text != '-':
                        final_result = 'loss'
                    elif '$0' in result_text:
                        final_result = 'push'
            
            # Only recalculate decimal_odds if we didn't already set it for pushed parlays
            if not (bet_type == 'Parlay' and has_push and final_result == 'win'):
                decimal_odds = calculate_decimal_odds_from_american(american_odds)

            # Extract line for the sheet (first leg's line)
            first_line = ''
            if wager_simples:
                first_line_elem = wager_simples[0].find('span', class_='line-selected')
                if first_line_elem:
                    first_line = parse_line(first_line_elem.get_text(strip=True), '')
            
            bet = {
                'date': date,
                'platform': 'Hoop88',
                'sport': sport,
                'description': description,
                'bet_type': bet_type,
                'line': first_line,
                'odds': american_odds,
                'bet_amount': risk,
                'dec': decimal_odds,
                'result': final_result or ''
            }
            
            bets.append(bet)
            print(f"  ✓ {date} - {bet_type} - {description[:60]}... - ${risk:.2f}")
            
        except Exception as e:
            print(f"  ✗ Error parsing bet row: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    return bets


def login(session: requests.Session) -> str:
    """Login to Hoop88 via JWT and return the auth token.

    Same auth approach as hoop88_odds/scraper.py:
    POST /cloud/api/System/authenticateCustomer with credentials.
    Returns JWT on success (HTTP 200). HTTP 204 = bad credentials.
    """
    domain = HOOP88_URL.replace("https://", "").replace("http://", "")
    resp = session.post(f"{API_BASE}/System/authenticateCustomer", data={
        "customerID": HOOP88_USERNAME,
        "password": HOOP88_PASSWORD,
        "state": "true",
        "multiaccount": "1",
        "response_type": "code",
        "client_id": HOOP88_USERNAME,
        "domain": domain,
        "redirect_uri": domain,
        "operation": "authenticateCustomer",
        "RRO": "1",
    }, timeout=15)

    if resp.status_code == 204:
        raise RuntimeError("Login failed — bad credentials (HTTP 204)")
    resp.raise_for_status()

    data = resp.json()
    token = data.get("code")
    if not token:
        raise RuntimeError(f"Login succeeded but no token in response: {list(data.keys())}")

    print("✅ Login successful")
    return token


def fetch_figures(session: requests.Session, token: str, weeks_back: int = 1) -> str:
    """Fetch the weekly figures HTML from Hoop88's figures API.

    The figures page contains a table with bet details for each day.
    The API endpoint mirrors what the JS frontend calls when you click
    on the Balance box and select a week.
    """
    headers = {"Authorization": f"Bearer {token}"}

    # Fetch the weekly figures — this returns HTML that includes bet rows
    resp = session.post(f"{API_BASE}/Figures/get-figure", data={
        "week": str(weeks_back),
    }, headers=headers, timeout=15)
    resp.raise_for_status()

    return resp.text


def fetch_figure_detail(session: requests.Session, token: str, weeks_back: int = 1) -> str:
    """Fetch the bet detail HTML for a given week.

    This is the equivalent of clicking the 'Week' total in the figures
    table — it returns the expanded bet-level detail rows.
    """
    headers = {"Authorization": f"Bearer {token}"}

    # Try the detail endpoint — this returns the expanded bet rows
    resp = session.post(f"{API_BASE}/Figures/get-figure-detail", data={
        "week": str(weeks_back),
        "index": "10",  # "10" = Week total (same as the old data-index="10")
    }, headers=headers, timeout=15)
    resp.raise_for_status()

    return resp.text


def scrape_hoop88(weeks_back: int = 1) -> list:
    """
    Log into Hoop88 via JWT and scrape bet history.

    Uses requests.Session instead of a browser — faster, more reliable,
    and doesn't break on Playwright selector changes.

    Args:
        weeks_back: Number of weeks back to fetch (0=This Week, 1=Last Week, etc.)

    Returns:
        List of parsed bet dictionaries
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env file")

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    })

    # Step 1: Authenticate via JWT
    print("Logging in to Hoop88...")
    token = login(session)

    # Step 2: Fetch bet detail HTML
    week_label = 'This Week' if weeks_back == 0 else 'Last Week' if weeks_back == 1 else f'{weeks_back} Weeks ago'
    print(f"Fetching bet details for {week_label}...")

    try:
        html = fetch_figure_detail(session, token, weeks_back)
    except requests.HTTPError as e:
        # If the detail endpoint doesn't work, try the full figures page
        print(f"Detail endpoint failed ({e}), trying figures page...")
        html = fetch_figures(session, token, weeks_back)

    if not html or len(html.strip()) < 50:
        print("No bet data returned from API")
        return []

    # Step 3: Parse with existing parser
    bets = parse_bets_from_html(html)

    return bets


def scrape_from_file(filepath: str) -> list:
    """Parse bets from a saved HTML file (for testing)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    return parse_bets_from_html(content)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from Hoop88')
    parser.add_argument('--test', action='store_true', help='Parse from saved HTML file instead of live scrape')
    parser.add_argument('--weeks', type=int, default=1, help='Weeks back to fetch (0=This Week, 1=Last Week, default: 1)')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    args = parser.parse_args()

    if args.test:
        print("TEST MODE: Parsing from saved HTML file")
        test_file = "hoop88_test.html"

        if not os.path.exists(test_file):
            print(f"Error: {test_file} not found")
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
        print("HOOP88 BET HISTORY SCRAPER")
        print("=" * 60)

        try:
            bets = scrape_hoop88(weeks_back=args.weeks)

            print(f"\n{'=' * 60}")
            print(f"Successfully scraped {len(bets)} bets from Hoop88")
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
