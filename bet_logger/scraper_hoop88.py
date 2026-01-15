#!/usr/bin/env python3
"""
Hoop88 Bet History Scraper
Scrapes bet history from Hoop88 and uploads to Google Sheets
"""

import os
import sys
import re
from datetime import datetime
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
HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")
HOOP88_HISTORY_URL = "https://hoop88.com/sports.html?v=1768432347363"


def calculate_american_odds_from_line(line_text: str) -> int:
    """Extract American odds from line text like '+120' or '-110'."""
    match = re.search(r'([+-]\d+)', line_text)
    if match:
        return int(match.group(1))
    return 0


def calculate_decimal_odds_from_american(american_odds: int) -> float:
    """Convert American odds to decimal odds format."""
    if american_odds > 0:
        return (american_odds / 100) + 1
    elif american_odds < 0:
        return (100 / abs(american_odds)) + 1
    return 1.0


def parse_status(status_text: str) -> str:
    """Extract result from status text."""
    status_upper = status_text.upper()
    if 'WIN' in status_upper or 'WON' in status_upper:
        return 'W'
    elif 'LOST' in status_upper or 'LOSS' in status_upper:
        return 'L'
    elif 'CANCELLED' in status_upper or 'PUSHED' in status_upper or 'PUSH' in status_upper:
        return 'X'
    return ''


def parse_sport(selection_text: str, note_text: str = '') -> str:
    """Determine sport from selection or note text."""
    combined_text = (selection_text + ' ' + note_text).upper()
    
    if 'NFL' in combined_text:
        return 'NFL'
    elif 'NBA' in combined_text:
        return 'NBA'
    elif 'NHL' in combined_text:
        return 'NHL'
    elif 'MLB' in combined_text:
        return 'MLB'
    elif 'NCAAF' in combined_text or 'COLLEGE FOOTBALL' in combined_text:
        return 'NCAAF'
    elif 'COLLEGE' in combined_text or 'BOWL' in combined_text or 'PEACH' in combined_text:
        # Generic college detection - includes bowl games
        return 'NCAAF'
    elif 'NCAAB' in combined_text or 'NCAAM' in combined_text or 'COLLEGE BASKETBALL' in combined_text:
        return 'NCAAM'
    
    return ''


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


def parse_risk_win(risk_win_text: str) -> tuple:
    """Parse risk/win text into separate values."""
    # Format: "$ 125.00 / $ 325.00" or "$ 0.00 / $ 603.48"
    try:
        parts = risk_win_text.split('/')
        if len(parts) == 2:
            risk = float(parts[0].replace('$', '').replace(',', '').strip())
            win = float(parts[1].replace('$', '').replace(',', '').strip())
            return risk, win
    except:
        pass
    return 0.0, 0.0


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
            all_odds = []
            final_result = None
            
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
                all_odds.append(odds)
                
                # Build leg description
                line_value = parse_line(line_text, period)
                leg_desc = clean_description(team_name, line_text, period)
                leg_descriptions.append(leg_desc)
                
                # Track result - if any leg loses, parlay loses
                if leg_result == 'L':
                    final_result = 'loss'
                elif leg_result == 'X' and final_result != 'loss':
                    final_result = 'push'
                elif leg_result == 'W' and final_result is None:
                    final_result = 'win'
            
            # Combine all legs into description
            description = ' | '.join(leg_descriptions)
            
            # For parlay, calculate combined odds from risk/win
            if bet_type == 'Parlay' and risk > 0:
                # Calculate from potential win
                american_odds = int(round((win / risk) * 100))
            else:
                # For straight bets, use the single leg odds
                american_odds = all_odds[0] if all_odds else 0
            
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
            
            decimal_odds = calculate_decimal_odds_from_american(american_odds)
            
            # Extract line for the sheet (first leg's line)
            first_line = parse_line(wager_simples[0].find('span', class_='line-selected').get_text(strip=True), '') if wager_simples else ''
            
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


def scrape_hoop88(weeks_back: int = 1, headless: bool = True) -> list:
    """
    Log into Hoop88 and scrape bet history.

    Args:
        weeks_back: Number of weeks back to fetch (0=This Week, 1=Last Week, etc.)
        headless: Run browser in headless mode (no visible window)

    Returns:
        List of parsed bet dictionaries
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env file")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Check if login form is present
        print("Checking login status...")
        username_field = page.locator('input[name="customerID"]')
        login_needed = username_field.count() > 0 and username_field.is_visible()

        if login_needed:
            print("Logging in...")
            try:
                page.fill('input[name="customerID"]', HOOP88_USERNAME)
                page.fill('input[name="Password"]', HOOP88_PASSWORD)
                page.click('button[data-action="login"]')
                page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
                page.wait_for_load_state('networkidle')
                page.wait_for_timeout(1500)  # Brief pause for page to stabilize after login
                print("✅ Login successful")
            except Exception as e:
                print(f"Auto-login failed: {e}")
                if not headless:
                    print("Please log in manually (30 seconds)...")
                    page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=30000)
                else:
                    raise RuntimeError("Login failed in headless mode")
        else:
            print("Already logged in")

        # Click the Balance box to open bet history
        print("Opening bet history...")
        try:
            page.evaluate('document.querySelector(\'div[data-action="get-figure"]\').click()')
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1000)
            print("✅ Opened bet history")
        except Exception as e:
            if not headless:
                print(f"Could not click Balance box: {e}")
                print("Please navigate to bet history manually (15 seconds)...")
                page.wait_for_timeout(15000)
            else:
                raise RuntimeError(f"Could not open bet history in headless mode: {e}")

        # Select the desired week
        print(f"Selecting week {weeks_back}...")
        try:
            dropdown = page.locator('select[data-list="week"]')
            dropdown.wait_for(state='attached', timeout=10000)
            dropdown.select_option(value=str(weeks_back))
            page.wait_for_load_state('networkidle')
            week_label = 'This Week' if weeks_back == 0 else 'Last Week' if weeks_back == 1 else f'{weeks_back} Weeks ago'
            print(f"✅ Selected: {week_label}")
        except Exception as e:
            print(f"Could not change week filter: {e}")

        # Wait for the weekly figures table to load (contains Carry, Mon, Tue, etc.)
        print("Waiting for weekly figures to load...")
        try:
            # Wait for the "Week" row to appear in the table
            page.wait_for_selector('text=Week', timeout=15000)
            page.wait_for_timeout(1000)  # Brief pause after table loads
            print("✅ Weekly figures loaded")
        except Exception as e:
            print(f"Warning: Could not confirm table loaded: {e}")

        # Click on Week value (data-index="10") to expand bet details
        print("Expanding bet details...")
        try:
            page.evaluate('document.querySelector(\'span[data-trigger="true"][data-index="10"]\').click()')
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(2000)
            print("✅ Expanded bet details")
        except Exception as e:
            print(f"Could not expand bet details: {e}")

        # Wait for bet rows to load
        print("Loading bet data...")
        try:
            page.wait_for_selector('tr[data-ticket]', timeout=15000)
            print("✅ Bet data loaded")
        except:
            print("Warning: No bet rows found")

        # Get the table HTML
        table_html = None
        for selector in ['#DataTables_Table_0', 'table.dataTable', 'table:has(tr[data-ticket])']:
            try:
                locator = page.locator(selector)
                if locator.count() > 0:
                    table_html = locator.first.evaluate('el => el.outerHTML')
                    break
            except:
                continue

        if not table_html:
            print("❌ Could not get table HTML")
            page.screenshot(path="debug_hoop88_error.png")
            browser.close()
            return []

        browser.close()
        return parse_bets_from_html(table_html)


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
    parser.add_argument('--visible', action='store_true', help='Show browser window (default is headless)')
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
            bets = scrape_hoop88(weeks_back=args.weeks, headless=not args.visible)

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

    sys.exit(0)
