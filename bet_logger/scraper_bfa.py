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

# Account configurations
ACCOUNTS = {
    'default': {
        'auth_file': os.path.join(SCRIPT_DIR, 'recon_bfa_auth.json'),
        'platform': 'Betfastaction',
        'bet_adjustment': 0,
        'shared_sheet': None,
    },
    'j': {
        'auth_file': os.path.join(SCRIPT_DIR, 'recon_bfaj_auth.json'),
        'platform': 'BFAJ',
        'bet_adjustment': -15,
        'shared_sheet': 'Shared',
    },
}

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


def load_auth(auth_file: str = None) -> dict:
    """Load saved auth data (refresh_token, player_id)."""
    path = auth_file or AUTH_FILE
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Auth file not found: {path}\n"
            "Run recon_bfa.py first to capture auth tokens."
        )
    with open(path, 'r') as f:
        data = json.load(f)
    data['_auth_file'] = path  # store path for token rotation save-back
    return data


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
        save_path = auth.get('_auth_file', AUTH_FILE)
        save_data = {k: v for k, v in auth.items() if not k.startswith('_')}
        with open(save_path, 'w') as f:
            json.dump(save_data, f, indent=2)

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


def parse_api_bets(wagers: list, platform: str = 'Betfastaction',
                   bet_adjustment: float = 0) -> list:
    """Convert API wager objects to the standard bet dict format.

    Args:
        platform: Platform label for the sheet (e.g. 'Betfastaction', 'BFAJ')
        bet_adjustment: Amount to add to each bet_amount (negative to subtract)
    """
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

            adjusted_risk = risk + bet_adjustment if bet_adjustment else risk

            bet = {
                'date': parse_date(w.get('placedDate', '')),
                'platform': platform,
                'sport': sport,
                'description': parsed['clean_desc'].strip(),
                'bet_type': bet_type,
                'line': parsed['line'],
                'odds': american_odds,
                'bet_amount': adjusted_risk,
                'dec': decimal_odds,
                'result': result_mapped,
                '_raw_risk': risk,  # original amount for verification
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


def get_last_bfa_date(platform: str = 'Betfastaction') -> str:
    """
    Query Google Sheets for the most recent bet date for a given platform.
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
            if len(row) >= 2 and row[1].strip() == platform:
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
        print(f"Could not fetch last {platform} date from sheet: {e}")
        return ''


# ── Main scraper ─────────────────────────────────────────────────


def scrape_bfa(days_back: int = 7, from_date: str = None,
               account_name: str = 'default') -> list:
    """
    Fetch bet history from BFA's API.

    Args:
        days_back: Number of days of history to fetch (default: 7).
                   Ignored if from_date is provided.
        from_date: Explicit start date as YYYY-MM-DD.
        account_name: Account key ('default' or 'j').

    Returns:
        List of parsed bet dictionaries.
    """
    acct = ACCOUNTS[account_name]

    # Load auth
    auth = load_auth(acct['auth_file'])
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
    return parse_api_bets(wagers, platform=acct['platform'],
                          bet_adjustment=acct['bet_adjustment'])


if __name__ == "__main__":
    import argparse
    import copy

    parser = argparse.ArgumentParser(description='Scrape bet history from BFA Gaming')
    parser.add_argument('--days', type=int, default=7, help='Days of history to fetch (default: 7)')
    parser.add_argument('--since-last', action='store_true', help='Scrape from the date of the last bet in Google Sheets')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    parser.add_argument('--refresh-only', action='store_true', help='Just refresh the auth token and exit')
    parser.add_argument('--account', choices=list(ACCOUNTS.keys()), default='default',
                        help='Which BFA account to scrape (default or j)')
    args = parser.parse_args()

    acct = ACCOUNTS[args.account]

    if args.refresh_only:
        try:
            auth = load_auth(acct['auth_file'])
            refresh_access_token(auth)
            print("Token refreshed successfully")
        except Exception as e:
            print(f"Token refresh failed: {e}")
            sys.exit(1)
        sys.exit(0)

    print("=" * 60)
    print(f"BFA GAMING BET HISTORY SCRAPER — {acct['platform']}")
    print("=" * 60)

    try:
        from_date = None
        if args.since_last:
            print(f"Looking up last {acct['platform']} bet date from Google Sheets...")
            from_date = get_last_bfa_date(platform=acct['platform'])
            if from_date:
                print(f"Last {acct['platform']} bet: {from_date}")
                from_date = (datetime.strptime(from_date, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
                print(f"Scraping from {from_date} to today")
            else:
                print(f"No existing {acct['platform']} bets found in sheet, using default range")

        bets = scrape_bfa(days_back=args.days, from_date=from_date,
                          account_name=args.account)

        print(f"\n{'=' * 60}")
        print(f"Successfully scraped {len(bets)} bets from {acct['platform']}")
        print(f"{'=' * 60}\n")

        if bets and not args.dry_run:
            # Upload adjusted bets to main sheet (Sheet1)
            print("Uploading to Google Sheets (main)...")
            result = append_bets_to_sheet(bets)

            if result['status'] == 'success':
                print(f"\nSUCCESS! Added {result['rows_added']} new bets to Sheet1")
                print(f"   Rows {result['start_row']} to {result['end_row']}")
            elif result['status'] == 'skipped':
                print(f"\n{result['message']}")
            else:
                print(f"\nError uploading to sheets: {result.get('message', 'Unknown error')}")

            # For accounts with shared_sheet, also upload raw (unadjusted) bets
            if acct['shared_sheet']:
                print(f"\nUploading raw bets to '{acct['shared_sheet']}' tab...")
                raw_bets = []
                for bet in bets:
                    raw = copy.copy(bet)
                    raw['bet_amount'] = raw.pop('_raw_risk', raw['bet_amount'])
                    raw_bets.append(raw)
                raw_result = append_bets_to_sheet(raw_bets, sheet_name=acct['shared_sheet'])
                if raw_result['status'] == 'success':
                    print(f"Added {raw_result['rows_added']} raw bets to {acct['shared_sheet']}")
                elif raw_result['status'] == 'skipped':
                    print(f"{raw_result['message']}")

        elif args.dry_run:
            print("Dry run - skipping upload to Google Sheets")
            if acct['bet_adjustment']:
                print(f"\nBet adjustment: ${acct['bet_adjustment']:+.0f} per bet")
                raw_total = sum(b.get('_raw_risk', 0) for b in bets)
                adj_total = sum(b['bet_amount'] for b in bets)
                print(f"Raw total wagered:      ${raw_total:,.2f}")
                print(f"Adjusted total wagered: ${adj_total:,.2f}")
        else:
            print("No bets found to upload")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    sys.exit(0)
