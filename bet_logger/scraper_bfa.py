#!/usr/bin/env python3
"""
BFA Gaming Bet History Scraper
Fetches bet history via BFA's REST API (no browser needed).

Auth flow:
  1. Read saved auth from recon_bfa_auth.json
  2. Use refresh token to get a fresh access token from Keycloak
  3. Call the GetPlayerHistory API endpoint with pagination
  4. Save rotated refresh token for next run

If the refresh token has expired (~1 hour without use), re-run recon_bfa.py
to capture fresh tokens via the browser.
"""

import os
import sys
import re
import json
import requests
from datetime import datetime, timedelta
from dotenv import load_dotenv
from utils import (
    calculate_american_odds,
    calculate_decimal_odds_from_american,
    parse_sport,
)

from sheets import append_bets_to_sheet, get_sheets_service, SPREADSHEET_ID, SHEET_NAME

load_dotenv()

# Paths
SCRIPT_DIR = os.path.dirname(__file__)
AUTH_FILE = os.path.join(SCRIPT_DIR, "recon_bfa_auth.json")

# API endpoints
TOKEN_URL = "https://auth.bfagaming.com/realms/players_realm/protocol/openid-connect/token"
BET_HISTORY_URL = "https://api.bfagaming.com/history/api/GetPlayerHistory"

# Shared request headers
BASE_HEADERS = {
    "Accept": "application/json",
    "Origin": "https://bfagaming.com",
    "Referer": "https://bfagaming.com/",
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36"
    ),
}

RECORDS_PER_PAGE = 100


def load_auth() -> dict:
    """Load saved auth data (refresh_token, player_id)."""
    if not os.path.exists(AUTH_FILE):
        raise FileNotFoundError(
            f"Auth file not found: {AUTH_FILE}\n"
            "Run recon_bfa.py first to capture auth tokens."
        )
    with open(AUTH_FILE, 'r') as f:
        return json.load(f)


def refresh_access_token(auth: dict) -> str:
    """
    Use the saved refresh token to get a fresh access token.
    Saves the rotated refresh token back to the auth file.

    Returns the new access token.
    Raises RuntimeError if the refresh token has expired.
    """
    refresh_token = auth.get('refresh_token')
    if not refresh_token:
        raise RuntimeError("No refresh_token in auth file. Run recon_bfa.py.")

    resp = requests.post(
        TOKEN_URL,
        data={
            "grant_type": "refresh_token",
            "refresh_token": refresh_token,
            "client_id": "bfagaming",
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
            "Refresh token likely expired (~1 hour). Run recon_bfa.py to re-authenticate."
        )

    token_data = resp.json()
    new_access = token_data["access_token"]
    new_refresh = token_data.get("refresh_token")

    # Save rotated refresh token
    if new_refresh and new_refresh != refresh_token:
        auth['refresh_token'] = new_refresh
        with open(AUTH_FILE, 'w') as f:
            json.dump(auth, f, indent=2)

    return new_access


def fetch_bet_history(access_token: str, player_id: str,
                      start_date: str, end_date: str) -> list:
    """
    Fetch all bet history pages from the API.

    Args:
        start_date: YYYY-MM-DD
        end_date: YYYY-MM-DD

    Returns:
        List of wager dicts from the API
    """
    all_wagers = []
    page = 0

    headers = {
        "Authorization": f"Bearer {access_token}",
        **BASE_HEADERS,
    }

    while True:
        resp = requests.get(
            BET_HISTORY_URL,
            params={
                "playerId": player_id,
                "startDate": start_date,
                "endDate": end_date,
                "page": page,
                "recordsByPage": RECORDS_PER_PAGE,
            },
            headers=headers,
            timeout=15,
        )

        if resp.status_code != 200:
            print(f"API error: HTTP {resp.status_code}")
            break

        data = resp.json()
        wagers = data.get("wagers", [])
        total_records = data.get("totalRecords", 0)

        all_wagers.extend(wagers)

        if not wagers or len(all_wagers) >= total_records:
            break

        page += 1

    return all_wagers


# ── Parsing helpers ──────────────────────────────────────────────


def parse_date(date_str: str) -> str:
    """Convert API date (ISO 8601) to sheet format (M/D/YY)."""
    try:
        dt = datetime.fromisoformat(date_str)
        return dt.strftime("%-m/%-d/%y")
    except (ValueError, TypeError):
        return date_str


def parse_description(raw_desc: str) -> dict:
    """
    Parse BFA bet description into structured components.

    Formats seen:
      [1601] TOTAL o67½+116 \\r(WASHINGTON 1H vrs RUTGERS 1H)
      [96072] TOTAL u38½-115 \\r(SAINT LOUIS 1H TEAM PTS vrs ...)
      [1306551] GRAMBLING 1H +168
      [1306564] LAMAR 1H -2-110

    Returns:
        {'clean_desc': str, 'line': str, 'odds': int}
    """
    # Normalize unicode fractions and whitespace
    desc = raw_desc.replace('\u00bd', '.5').replace('\u00bc', '.25').replace('\u00be', '.75')
    desc = desc.replace('\r', ' ').replace('\n', ' ')
    desc = re.sub(r'^\[\d+\]\s*', '', desc).strip()
    desc = re.sub(r'\s+', ' ', desc)

    # Pattern 1: TOTAL o/u LINE ODDS (CONTEXT)
    # Odds can be +/-NNN or "EV" (even money = +100)
    total_match = re.match(
        r'TOTAL\s+([ou])([\d.]+)((?:[+-]\d+)|EV)\s*\((.+?)\)\s*$',
        desc, re.IGNORECASE
    )
    if total_match:
        direction = 'Over' if total_match.group(1).lower() == 'o' else 'Under'
        line_val = total_match.group(2)
        odds_str = total_match.group(3)
        odds = 100 if odds_str.upper() == 'EV' else int(odds_str)
        context = total_match.group(4).strip()

        is_team_total = 'TEAM PTS' in context.upper()

        if is_team_total:
            # "SAINT LOUIS 1H TEAM PTS vrs SAINT LOUIS 1H TEAM PTS"
            team_match = re.match(r'(.+?)\s+\d+H\s+TEAM PTS', context, re.IGNORECASE)
            team = team_match.group(1).strip() if team_match else context.split(' vrs ')[0].strip()
            clean = f"{team} 1H Team Total {direction} {line_val}"
        else:
            # "WASHINGTON 1H vrs RUTGERS 1H"
            clean = f"{context} {direction} {line_val}"

        return {'clean_desc': clean, 'line': f"{direction} {line_val}", 'odds': odds}

    # Pattern 2: TEAM PERIOD SPREAD ODDS (e.g., "LAMAR 1H -2-110")
    spread_match = re.match(r'(.+?)\s+(\d+H)\s+([+-][\d.]+)([+-]\d+)\s*$', desc)
    if spread_match:
        team = spread_match.group(1).strip()
        period = spread_match.group(2)
        spread = spread_match.group(3)
        odds = int(spread_match.group(4))
        return {'clean_desc': f"{team} {period} {spread}", 'line': spread, 'odds': odds}

    # Pattern 3: TEAM PERIOD ODDS (moneyline, e.g., "GRAMBLING 1H +168")
    ml_match = re.match(r'(.+?)\s+(\d+H)\s+([+-]\d+)\s*$', desc)
    if ml_match:
        team = ml_match.group(1).strip()
        period = ml_match.group(2)
        odds = int(ml_match.group(3))
        return {'clean_desc': f"{team} {period} {'+' if odds > 0 else ''}{odds}", 'line': '', 'odds': odds}

    # Pattern 4: Full game spread (TEAM SPREAD ODDS, no period)
    spread_fg = re.match(r'(.+?)\s+([+-][\d.]+)([+-]\d+)\s*$', desc)
    if spread_fg:
        team = spread_fg.group(1).strip()
        spread = spread_fg.group(2)
        odds = int(spread_fg.group(3))
        return {'clean_desc': f"{team} {spread}", 'line': spread, 'odds': odds}

    # Pattern 5: Full game ML (TEAM ODDS, no period)
    ml_fg = re.match(r'(.+?)\s+([+-]\d+)\s*$', desc)
    if ml_fg:
        team = ml_fg.group(1).strip()
        odds = int(ml_fg.group(2))
        return {'clean_desc': team, 'line': '', 'odds': odds}

    # Fallback
    odds = 0
    odds_match = re.search(r'([+-]\d+)\s*$', desc)
    if odds_match:
        odds = int(odds_match.group(1))

    return {'clean_desc': desc, 'line': '', 'odds': odds}


def parse_api_bets(wagers: list) -> list:
    """Convert API wager objects to the standard bet dict format."""
    bets = []

    for w in wagers:
        try:
            raw_desc = w.get('description', '')
            wager_type = w.get('type', '')
            result = w.get('result', '')
            risk = w.get('risk', 0) or 0
            win_amount = w.get('win', 0) or 0

            # Skip pending/open bets
            if not result or result.upper() in ('PENDING', 'OPEN'):
                continue

            # Parse description
            parsed = parse_description(raw_desc)

            # If parsing didn't extract odds, calculate from risk/win
            american_odds = parsed['odds']
            if american_odds == 0 and risk > 0 and win_amount > 0:
                american_odds = calculate_american_odds(risk, win_amount)

            decimal_odds = calculate_decimal_odds_from_american(american_odds) if american_odds != 0 else 0

            # Map bet type
            bet_type = 'Straight'
            wt = wager_type.upper()
            if 'PARLAY' in wt:
                bet_type = 'Parlay'
            elif 'TEASER' in wt:
                bet_type = 'Teaser'

            # Map result
            result_mapped = ''
            r = result.upper()
            if r == 'WIN':
                result_mapped = 'win'
            elif r == 'LOSE':
                result_mapped = 'loss'
            elif r in ('PUSH', 'CANCELLED', 'VOID'):
                result_mapped = 'push'

            # Sport detection — BFA is primarily CBB, default to NCAAM
            sport = parse_sport(raw_desc)
            if not sport:
                sport = 'NCAAM'

            bet = {
                'date': parse_date(w.get('placedDate', '')),
                'platform': 'BFA',
                'sport': sport,
                'description': parsed['clean_desc'].strip(),
                'bet_type': bet_type,
                'line': parsed['line'],
                'odds': american_odds,
                'bet_amount': risk,
                'dec': decimal_odds,
                'result': result_mapped,
            }

            bets.append(bet)
            desc_preview = bet['description'][:60]
            print(f"  {bet['date']} - {bet_type} - {desc_preview} - ${risk:.2f} - {bet['result']}")

        except Exception as e:
            print(f"  Error parsing wager: {e}")
            import traceback
            traceback.print_exc()
            continue

    return bets


# ── Sheet lookup ─────────────────────────────────────────────────


def get_last_bfa_date() -> str:
    """
    Query Google Sheets for the most recent BFA bet date.
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
            if len(row) >= 2 and row[1].strip() == 'BFA':
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
        print(f"Could not fetch last BFA date from sheet: {e}")
        return ''


# ── Main scraper ─────────────────────────────────────────────────


def scrape_bfa(days_back: int = 7, from_date: str = None) -> list:
    """
    Fetch bet history from BFA's API.

    Args:
        days_back: Number of days of history to fetch (default: 7).
                   Ignored if from_date is provided.
        from_date: Explicit start date as YYYY-MM-DD.

    Returns:
        List of parsed bet dictionaries.
    """
    # Load auth
    auth = load_auth()
    player_id = auth.get('player_id')
    if not player_id:
        raise RuntimeError("No player_id in auth file. Run recon_bfa.py.")

    print("Refreshing access token...")
    access_token = refresh_access_token(auth)
    print("Token refreshed")

    # Build date range
    today = datetime.now().strftime("%Y-%m-%d")
    if from_date:
        start = from_date
    else:
        start = (datetime.now() - timedelta(days=days_back)).strftime("%Y-%m-%d")

    print(f"Fetching bets from {start} to {today}...")
    wagers = fetch_bet_history(access_token, player_id, start, today)

    print(f"\nAPI returned {len(wagers)} wagers")
    return parse_api_bets(wagers)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from BFA Gaming')
    parser.add_argument('--days', type=int, default=7, help='Days of history to fetch (default: 7)')
    parser.add_argument('--since-last', action='store_true', help='Scrape from the date of the last BFA bet in Google Sheets')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    parser.add_argument('--refresh-only', action='store_true', help='Just refresh the auth token and exit')
    args = parser.parse_args()

    if args.refresh_only:
        try:
            auth = load_auth()
            refresh_access_token(auth)
            print("Token refreshed successfully")
        except Exception as e:
            print(f"Token refresh failed: {e}")
            sys.exit(1)
        sys.exit(0)

    print("=" * 60)
    print("BFA GAMING BET HISTORY SCRAPER")
    print("=" * 60)

    try:
        from_date = None
        if args.since_last:
            print("Looking up last BFA bet date from Google Sheets...")
            from_date = get_last_bfa_date()
            if from_date:
                print(f"Last BFA bet: {from_date}")
                # Start from day after last bet to avoid re-fetching
                from_date = (datetime.strptime(from_date, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
                print(f"Scraping from {from_date} to today")
            else:
                print("No existing BFA bets found in sheet, using default range")

        bets = scrape_bfa(days_back=args.days, from_date=from_date)

        print(f"\n{'=' * 60}")
        print(f"Successfully scraped {len(bets)} bets from BFA Gaming")
        print(f"{'=' * 60}\n")

        if bets and not args.dry_run:
            print("Uploading to Google Sheets...")
            result = append_bets_to_sheet(bets)

            if result['status'] == 'success':
                print(f"\nSUCCESS! Added {result['rows_added']} new bets to sheet")
                print(f"   Rows {result['start_row']} to {result['end_row']}")
            elif result['status'] == 'skipped':
                print(f"\n{result['message']}")
            else:
                print(f"\nError uploading to sheets: {result.get('message', 'Unknown error')}")
        elif args.dry_run:
            print("Dry run - skipping upload to Google Sheets")
        else:
            print("No bets found to upload")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    sys.exit(0)
