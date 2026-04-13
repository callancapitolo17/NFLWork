#!/usr/bin/env python3
"""
BetOnline Bet History Scraper
Fetches bet history via BetOnline's REST API (no browser needed).

Auth flow:
  1. Read saved cookies/tokens from recon_betonline_cookies.json
  2. Use refresh token to get a fresh access token
  3. Call the bet history API endpoint
  4. Save rotated refresh token for next run

If the refresh token has expired (>3 days without use), re-run recon_betonline.py
to capture fresh tokens via the browser.
"""

import os
import sys
import re
import json
import requests
from datetime import datetime, timezone, timedelta
from dotenv import load_dotenv
from utils import (
    calculate_american_odds,
    calculate_american_odds_from_line,
    calculate_decimal_odds_from_american,
    parse_sport,
)

from sheets import append_bets_to_sheet, get_sheets_service, SPREADSHEET_ID, SHEET_NAME

load_dotenv()

# Paths
SCRIPT_DIR = os.path.dirname(__file__)
COOKIES_FILE = os.path.join(SCRIPT_DIR, "recon_betonline_cookies.json")

# API endpoints
TOKEN_URL = "https://api.betonline.ag/api/auth/realms/betonline/protocol/openid-connect/token"
BET_HISTORY_URL = "https://api.betonline.ag/report/api/report/get-bet-history"

# Shared request headers
BASE_HEADERS = {
    "Accept": "application/json, text/plain, */*",
    "Origin": "https://www.betonline.ag",
    "Referer": "https://www.betonline.ag/",
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
    ),
    "gsetting": "bolnasite",
    "contests": "na",
    "gmt-offset": "-8",
    "utc-offset": "480",
}


def load_cookies() -> list:
    """Load saved cookies from recon output."""
    if not os.path.exists(COOKIES_FILE):
        raise FileNotFoundError(
            f"Cookie file not found: {COOKIES_FILE}\n"
            "Run recon_betonline.py first to capture auth tokens."
        )
    with open(COOKIES_FILE, 'r') as f:
        return json.load(f)


def build_session(cookies: list) -> requests.Session:
    """Build a requests session with BetOnline Cloudflare cookies."""
    session = requests.Session()
    for c in cookies:
        if 'betonline.ag' in c.get('domain', ''):
            session.cookies.set(
                c['name'], c['value'],
                domain=c['domain'],
                path=c.get('path', '/'),
            )
    return session


def refresh_access_token(session: requests.Session, cookies: list) -> str:
    """
    Use the saved refresh token to get a fresh access token.
    Saves the rotated refresh token back to the cookies file.

    Returns the new access token.
    Raises RuntimeError if the refresh token has expired.
    """
    krefresh = None
    for c in cookies:
        if c['name'] == 'krefresh':
            krefresh = c['value']
            break

    if not krefresh:
        raise RuntimeError("No refresh token (krefresh) found in cookies. Run recon_betonline.py.")

    resp = session.post(
        TOKEN_URL,
        data={
            "grant_type": "refresh_token",
            "refresh_token": krefresh,
            "client_id": "betonline-web",
        },
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            **{k: v for k, v in BASE_HEADERS.items() if k in ("User-Agent", "Origin", "Referer")},
        },
        timeout=15,
    )

    if resp.status_code != 200:
        raise RuntimeError(
            f"Token refresh failed (HTTP {resp.status_code}). "
            "Refresh token likely expired (>3 days). Run recon_betonline.py to re-authenticate."
        )

    token_data = resp.json()
    new_access = token_data["access_token"]
    new_refresh = token_data.get("refresh_token")

    # Save rotated refresh token back to cookies file
    if new_refresh and new_refresh != krefresh:
        for c in cookies:
            if c['name'] == 'krefresh':
                c['value'] = new_refresh
                break
        with open(COOKIES_FILE, 'w') as f:
            json.dump(cookies, f, indent=2)

    return new_access


def build_api_headers(access_token: str) -> dict:
    """Build the full header set for BetOnline API calls."""
    now = datetime.now(timezone.utc)
    ms = int(now.timestamp() * 1000)

    return {
        "Authorization": f"Bearer {access_token}",
        "Content-Type": "application/json",
        **BASE_HEADERS,
        "actual-time": str(ms),
        "iso-time": now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{ms % 1000:03d}Z",
        "utc-time": now.strftime("%a, %d %b %Y %H:%M:%S GMT"),
    }


def fetch_bet_history(session: requests.Session, headers: dict,
                      start_date: str, end_date: str) -> list:
    """
    Fetch bet history from the API.

    Args:
        start_date: ISO date string (YYYY-MM-DD), converted to midnight UTC
        end_date: ISO date string (YYYY-MM-DD), converted to midnight UTC

    Returns:
        List of bet dicts from the API's Data array
    """
    all_bets = []
    start_pos = 0
    page_size = 100

    while True:
        resp = session.post(
            BET_HISTORY_URL,
            headers=headers,
            json={
                "Id": None,
                "StartDate": f"{start_date}T00:00:00.000Z",
                "EndDate": f"{end_date}T00:00:00.000Z",
                "Status": None,
                "Product": None,
                "WagerType": None,
                "FreePlayFlag": None,
                "StartPosition": start_pos,
                "TotalPerPage": page_size,
                "IsDailyFigureReport": False,
            },
            timeout=15,
        )

        if resp.status_code != 200:
            print(f"API error: HTTP {resp.status_code}")
            break

        data = resp.json()
        bets = data.get("Data", [])
        total_rows = data.get("TotalRows", 0)

        all_bets.extend(bets)

        if len(all_bets) >= total_rows or not bets:
            break

        start_pos += page_size

    return all_bets


# ── Parsing helpers ──────────────────────────────────────────────


def parse_date(date_str: str) -> str:
    """Convert API date (ISO 8601) to sheet format (M/D/YY)."""
    try:
        dt = datetime.fromisoformat(date_str.replace("Z", "+00:00"))
        return dt.strftime("%-m/%-d/%y")
    except (ValueError, TypeError):
        return date_str


def parse_description_odds(desc: str) -> int:
    """Extract American odds from description (e.g., '... under 63 -110')."""
    match = re.search(r'([+-]\d+)\s*(?:for\s|$)', desc.strip())
    if match:
        return int(match.group(1))
    # Fallback: last +/- number
    match = re.search(r'([+-]\d+)\s*$', desc.strip())
    if match:
        return int(match.group(1))
    return 0


def parse_description_line(desc: str, wager_type: str) -> str:
    """Extract line/spread from description."""
    if wager_type == 'Money Line':
        return ''

    ou_match = re.search(r'(under|over)\s+([\d½¼¾.]+)', desc, re.IGNORECASE)
    if ou_match:
        direction = ou_match.group(1).capitalize()
        value = ou_match.group(2).replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
        return f"{direction} {value}"

    spread_match = re.search(r'([+-][\d½¼¾.]+)\s+[+-]\d+', desc)
    if spread_match:
        value = spread_match.group(1).replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
        return value

    return ''


def clean_description(desc: str, wager_type: str) -> str:
    """Clean the API description for the Google Sheet."""
    # Remove 'Desktop - ' or 'Mobile - ' prefix
    desc = re.sub(r'^(Desktop|Mobile)\s*-\s*', '', desc)

    if wager_type == 'Same Game Parlay':
        match = re.match(r'[A-Z]+\s*-\s*[A-Z]+\s*-\s*', desc)
        if match:
            return desc[match.end():].strip()
        match = re.match(r'[A-Z]+\s*-\s*', desc)
        if match:
            return desc[match.end():].strip()
        return desc.strip()

    # Straight bet: SPORT - NUM Description
    match = re.match(r'[A-Z]+\s*-\s*\d+\s+', desc)
    if match:
        return desc[match.end():].strip()

    return desc.strip()



def map_bet_type(wager_type: str) -> str:
    """Map API WagerType to standard format."""
    mapping = {
        'Spread': 'Straight',
        'Total': 'Straight',
        'Money Line': 'Straight',
        'Same Game Parlay': 'Parlay',
    }
    return mapping.get(wager_type, wager_type)


def map_status(wager_status: str) -> str:
    """Map API WagerStatus to standard format."""
    s = wager_status.strip().lower()
    if s == 'won':
        return 'win'
    elif s == 'lost':
        return 'loss'
    elif s in ('cancelled', 'push', 'void'):
        return 'push'
    return ''


def parse_api_bets(api_bets: list) -> list:
    """Convert API bet objects to the standard bet dict format."""
    bets = []

    for b in api_bets:
        try:
            raw_desc = b.get('Description', '')
            wager_type = b.get('WagerType', '')
            wager_status = b.get('WagerStatus', '')
            risk = b.get('Risk', 0) or 0
            to_win = b.get('ToWin', 0) or 0

            # Skip pending bets
            if wager_status.lower() == 'pending':
                continue

            # Calculate odds
            if wager_type == 'Same Game Parlay':
                american_odds = calculate_american_odds(risk, to_win) if to_win > 0 and risk > 0 else 0
            else:
                american_odds = parse_description_odds(raw_desc)
                if american_odds == 0 and to_win > 0 and risk > 0:
                    american_odds = calculate_american_odds(risk, to_win)

            decimal_odds = calculate_decimal_odds_from_american(american_odds) if american_odds != 0 else 0

            bet = {
                'date': parse_date(b.get('Date', '')),
                'platform': 'BetOnline',
                'sport': parse_sport(raw_desc),
                'description': clean_description(raw_desc, wager_type),
                'bet_type': map_bet_type(wager_type),
                'line': parse_description_line(raw_desc, wager_type),
                'odds': american_odds,
                'bet_amount': risk,
                'dec': decimal_odds,
                'result': map_status(wager_status),
            }

            bets.append(bet)
            desc_preview = bet['description'][:60]
            print(f"  ✓ {bet['date']} - {wager_type} - {desc_preview} - ${risk:.2f} - {bet['result']}")

        except Exception as e:
            print(f"  ✗ Error parsing bet: {e}")
            import traceback
            traceback.print_exc()
            continue

    return bets


# ── Sheet lookup ─────────────────────────────────────────────────


def get_last_betonline_date() -> str:
    """
    Query Google Sheets for the most recent BetOnline bet date.
    Returns date as YYYY-MM-DD for the API, or '' if none found.
    """
    try:
        service = get_sheets_service()
        result = service.spreadsheets().values().get(
            spreadsheetId=SPREADSHEET_ID,
            range=f"{SHEET_NAME}!A:B"
        ).execute()

        values = result.get('values', [])
        latest = None

        for row in values[1:]:
            if len(row) >= 2 and row[1].strip() == 'BetOnline':
                date_str = row[0].strip()
                dt = None
                for fmt in ("%m/%d/%Y", "%m/%d/%y"):
                    try:
                        dt = datetime.strptime(date_str, fmt)
                        break
                    except ValueError:
                        continue
                if dt and (latest is None or dt > latest):
                    latest = dt

        if latest:
            return latest.strftime("%Y-%m-%d")
        return ''

    except Exception as e:
        print(f"Could not fetch last BetOnline date from sheet: {e}")
        return ''


# ── Main scraper ─────────────────────────────────────────────────


def scrape_betonline(days_back: int = 7, from_date: str = None) -> list:
    """
    Fetch bet history from BetOnline's API.

    Args:
        days_back: Number of days of history to fetch (default: 7).
                   Ignored if from_date is provided.
        from_date: Explicit start date as YYYY-MM-DD.

    Returns:
        List of parsed bet dictionaries.
    """
    # Load auth
    cookies = load_cookies()
    session = build_session(cookies)

    print("Refreshing access token...")
    access_token = refresh_access_token(session, cookies)
    print("✅ Token refreshed")

    # Build date range (use local time — BetOnline rejects future EndDates in UTC)
    today = datetime.now().strftime("%Y-%m-%d")
    if from_date:
        start = from_date
    else:
        start = (datetime.now() - timedelta(days=days_back)).strftime("%Y-%m-%d")

    print(f"Fetching bets from {start} to {today}...")
    headers = build_api_headers(access_token)
    api_bets = fetch_bet_history(session, headers, start, today)

    print(f"\nAPI returned {len(api_bets)} bets")
    return parse_api_bets(api_bets)


def scrape_from_file(filepath: str) -> list:
    """Parse bets from a saved API JSON response (for testing)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    api_bets = data if isinstance(data, list) else data.get('Data', [])
    return parse_api_bets(api_bets)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from BetOnline')
    parser.add_argument('--test', action='store_true', help='Parse from saved JSON file instead of live API')
    parser.add_argument('--days', type=int, default=7, help='Days of history to fetch (default: 7)')
    parser.add_argument('--since-last', action='store_true', help='Scrape from the date of the last BetOnline bet in Google Sheets')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    parser.add_argument('--refresh-only', action='store_true', help='Just refresh the auth token and exit (keeps token alive)')
    args = parser.parse_args()

    if args.refresh_only:
        try:
            cookies = load_cookies()
            session = build_session(cookies)
            refresh_access_token(session, cookies)
            print("Token refreshed successfully")
        except Exception as e:
            print(f"Token refresh failed: {e}")
            sys.exit(1)
        sys.exit(0)

    if args.test:
        print("TEST MODE: Parsing from saved JSON file")
        test_file = os.path.join(SCRIPT_DIR, "recon_betonline_api_response.json")

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
            odds_str = f"{bet['odds']:+d}" if bet['odds'] else '?'
            dec_str = f"{bet['dec']:.2f}" if bet['dec'] else '?'
            print(f"   Odds: {odds_str} (Decimal: {dec_str})")
            print(f"   Amount: ${bet['bet_amount']:.2f}")
            print(f"   Result: {bet['result']}")
    else:
        print("=" * 60)
        print("BETONLINE BET HISTORY SCRAPER")
        print("=" * 60)

        try:
            from_date = None
            if args.since_last:
                print("Looking up last BetOnline bet date from Google Sheets...")
                from_date = get_last_betonline_date()
                if from_date:
                    print(f"Last BetOnline bet: {from_date}")
                    # Start from day after last bet to avoid re-fetching existing bets
                    from_date = (datetime.strptime(from_date, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
                    print(f"Scraping from {from_date} to today")
                else:
                    print("No existing BetOnline bets found in sheet, using default range")

            bets = scrape_betonline(days_back=args.days, from_date=from_date)

            print(f"\n{'=' * 60}")
            print(f"Successfully scraped {len(bets)} bets from BetOnline")
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
