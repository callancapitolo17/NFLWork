#!/usr/bin/env python3
"""
Wagerzon Bet History Scraper
Scrapes bet history from Wagerzon and uploads to Google Sheets.

Auth: ASP.NET form POST (same approach as wagerzon_odds/scraper_v2.py).
Data: HistoryHelper.aspx JSON API — returns structured bet data with
sport IDs, risk/win amounts, and individual leg details.

No browser/Playwright required — pure HTTP via requests.Session.
"""

import re
import os
import sys
from datetime import datetime
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
WAGERZON_HISTORY_URL = f"{WAGERZON_BASE_URL}/wager/HistoryHelper.aspx"
WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")


# ── Auth ────────────────────────────────────────────────────────


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
    if "History" in resp.url or "NewSchedule" in resp.url or "Welcome" in resp.url:
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


# ── JSON API ────────────────────────────────────────────────────


def fetch_history_json(session: requests.Session, week: int = 0) -> dict:
    """Fetch bet history JSON from HistoryHelper.aspx.

    The React widget on History.aspx calls this endpoint to load data.
    Returns structured JSON with daily breakdowns and individual wagers.

    Args:
        week: Week offset (0 = current week, 1 = last week, etc.)
    """
    resp = session.get(
        WAGERZON_HISTORY_URL,
        params={"week": week},
        headers={
            "X-Requested-With": "XMLHttpRequest",
            "Accept": "application/json",
        },
        timeout=15,
    )
    resp.raise_for_status()

    data = resp.json()
    result = data.get("result", {})

    if "ErrorMessage" in result:
        raise RuntimeError(f"API error: {result['ErrorMessage']}")

    return result


# ── Parsing ─────────────────────────────────────────────────────


def clean_desc(raw: str) -> str:
    """Clean a DetailDesc string from the API.

    The API returns descriptions like:
      "[967] TOTAL o7½-120 (CHI CUBS vrs TB RAYS)<BR>( JAMESON TAILLON - R / SHANE MCCLANAHAN - L )"

    We strip the rotation number, convert unicode fractions, and clean HTML.
    """
    desc = raw.replace('\u00bd', '.5').replace('\u00bc', '.25').replace('\u00be', '.75')
    desc = desc.replace('<BR>', ' ').replace('<br>', ' ')
    desc = re.sub(r'^\[\d+\]\s*', '', desc)
    desc = re.sub(r'\s+', ' ', desc).strip()
    return desc


def parse_line_from_desc(desc: str) -> str:
    """Extract line/spread from a cleaned description."""
    # Over/Under: "TOTAL o7.5-120 ..."
    ou = re.search(r'TOTAL\s+([ou])([\d.]+)', desc, re.IGNORECASE)
    if ou:
        direction = 'Over' if ou.group(1).lower() == 'o' else 'Under'
        return f"{direction} {ou.group(2)}"

    # Spread: "TEAM -3.5-110"
    spread = re.search(r'([+-][\d.]+)(?:[+-]\d+|EV)\s*(?:\(|$)', desc)
    if spread:
        return spread.group(1)

    return ''


def parse_odds_from_desc(desc: str) -> int:
    """Extract American odds from a cleaned description.

    Patterns: "-110", "+150", "EV" (even = +100)
    """
    # TOTAL o7.5-120 or TEAM -3.5+150
    match = re.search(r'[\d.]+([+-]\d+)\s*(?:\(|$)', desc)
    if match:
        return int(match.group(1))

    # EV (even money)
    if re.search(r'[\d.]+EV\s*(?:\(|$)', desc, re.IGNORECASE):
        return 100

    # Moneyline: "TEAM +168"
    match = re.search(r'\s([+-]\d+)\s*(?:\(|$)', desc)
    if match:
        return int(match.group(1))

    return 0


def map_sport(id_sport: str) -> str:
    """Map Wagerzon's IdSport field to our standard sport labels.

    The API returns sport IDs like "MLB", "NFL", "NBA", etc.
    This is the primary sport detection — no guessing from team names.
    """
    mapping = {
        'MLB': 'MLB',
        'NFL': 'NFL',
        'NBA': 'NBA',
        'NHL': 'NHL',
        'CFB': 'NCAAF',
        'CBK': 'NCAAM',
        'SOC': 'Soccer',
        'TNS': 'Tennis',
        'MMA': 'MMA',
        'BOX': 'Boxing',
        'GLF': 'Golf',
        'NASCAR': 'NASCAR',
        'WNBA': 'WNBA',
    }
    return mapping.get(id_sport, id_sport)


def map_result(result_str: str) -> str:
    """Map API result to standard format."""
    r = result_str.upper()
    if r == 'WIN':
        return 'win'
    elif r == 'LOSE':
        return 'loss'
    elif r in ('PUSH', 'CANCELLED', 'VOID', 'NO ACTION'):
        return 'push'
    return ''


def parse_bet_type(header_desc: str) -> str:
    """Extract bet type from the wager's HeaderDesc field.

    Examples: "STRAIGHT", "PARLAY (2 TEAMS)", "TEASER (3 TEAMS)"
    """
    h = header_desc.upper()
    if 'PARLAY' in h:
        return 'Parlay'
    elif 'TEASER' in h:
        return 'Teaser'
    elif 'STRAIGHT' in h:
        return 'Straight'
    elif 'PROP' in h:
        return 'Prop'
    elif 'IF BET' in h:
        return 'If Bet'
    return 'Straight'


def parse_date(date_str: str) -> str:
    """Convert API date (MM/DD/YYYY) to sheet format (M/D/YY)."""
    try:
        dt = datetime.strptime(date_str, '%m/%d/%Y')
        return dt.strftime('%-m/%-d/%y')
    except (ValueError, TypeError):
        return date_str


def parse_api_bets(history: dict) -> list:
    """Convert HistoryHelper JSON to the standard bet dict format.

    The JSON has a 'details' list where each item is a day with a
    list of wagers. Each wager has legs in its own 'details' list.
    """
    bets = []
    days = history.get('details', [])

    for day in days:
        for wager in day.get('wager', []):
            try:
                # Skip non-wager transactions (cash in/out, transfers)
                if wager.get('WagerOrTrans') != 'WAGER':
                    continue

                risk = float(wager.get('RiskAmount', 0))
                win = float(wager.get('WinAmount', 0))
                result = map_result(wager.get('Result', ''))

                # Skip pending/open bets
                if not result:
                    continue

                # Calculate odds from risk/win
                american_odds = calculate_american_odds(risk, win) if risk > 0 and win > 0 else 0
                decimal_odds = calculate_decimal_odds_from_american(american_odds) if american_odds else 0

                # Build description from legs
                legs = wager.get('details', [])
                leg_descs = []
                line = ''
                sport = ''

                for leg in legs:
                    raw = leg.get('DetailDesc', '')
                    cleaned = clean_desc(raw)
                    leg_descs.append(cleaned)

                    # Use first leg for line and sport
                    if not line:
                        line = parse_line_from_desc(cleaned)
                    if not sport:
                        # IdSport from the API is the primary source
                        id_sport = leg.get('IdSport', '')
                        sport = map_sport(id_sport) if id_sport else ''

                # Fallback: parse sport from description text if API didn't provide it
                if not sport:
                    sport = parse_sport(' '.join(leg_descs))

                description = ' | '.join(leg_descs)
                bet_type = parse_bet_type(wager.get('HeaderDesc', ''))

                # For straight bets, try to get more precise odds from the description
                if bet_type == 'Straight' and len(legs) == 1:
                    desc_odds = parse_odds_from_desc(clean_desc(legs[0].get('DetailDesc', '')))
                    if desc_odds:
                        american_odds = desc_odds
                        decimal_odds = calculate_decimal_odds_from_american(american_odds)

                bet = {
                    'date': parse_date(wager.get('PlacedDate', '')),
                    'platform': 'Wagerzon',
                    'sport': sport,
                    'description': description,
                    'bet_type': bet_type,
                    'line': line,
                    'odds': american_odds,
                    'bet_amount': risk,
                    'dec': decimal_odds,
                    'result': result,
                }

                bets.append(bet)
                print(f"  {bet['date']} - {bet_type} - {description[:60]} - ${risk:.2f} - {result}")

            except Exception as e:
                print(f"  Error parsing wager: {e}")
                continue

    return bets


# ── Main ────────────────────────────────────────────────────────


def scrape_wagerzon(weeks_back: int = 1) -> list:
    """
    Log into Wagerzon via HTTP and fetch bet history from the JSON API.

    Uses HistoryHelper.aspx — the same endpoint the React frontend calls.
    Returns structured data with sport IDs, so no HTML parsing or team
    name guessing needed.

    Args:
        weeks_back: Which week to fetch (0 = current, 1 = last week, etc.)

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

    # Step 2: Fetch bet history JSON
    print(f"Fetching bet history (week {weeks_back})...")
    history = fetch_history_json(session, week=weeks_back)

    days = history.get('details', [])
    total_wagers = sum(len(d.get('wager', [])) for d in days)
    print(f"API returned {total_wagers} wagers across {len(days)} days\n")

    # Step 3: Parse into standard bet format
    bets = parse_api_bets(history)

    return bets


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from Wagerzon')
    parser.add_argument('--weeks', type=int, default=1,
                        help='Weeks back to fetch (0=This Week, 1=Last Week, default: 1)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Scrape but do not upload to Google Sheets')
    args = parser.parse_args()

    print("=" * 60)
    print("WAGERZON BET HISTORY SCRAPER")
    print("=" * 60)

    try:
        bets = scrape_wagerzon(weeks_back=args.weeks)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"Successfully scraped {len(bets)} bets from Wagerzon")
    print(f"{'=' * 60}\n")

    if bets and not args.dry_run:
        from sheets import append_bets_to_sheet
        result = append_bets_to_sheet(bets)

        if result['status'] == 'success':
            print(f"\n✅ SUCCESS! Added {result['rows_added']} new bets to sheet")
            print(f"   Rows {result['start_row']} to {result['end_row']}")
        elif result['status'] == 'skipped':
            print(f"\n⚠️  {result['message']}")
        else:
            print(f"\n❌ Error uploading to sheets: {result.get('message', 'Unknown error')}")
    elif args.dry_run:
        print("Dry run — skipping upload to Google Sheets")
    elif not bets:
        print("No bets found to upload")
        sys.exit(1)
