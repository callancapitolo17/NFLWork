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
# Account configurations.
# Wagerzon supports two accounts using the same credentials flow but
# different sheet labels and stake adjustments. The primary uses a
# multiplier of 1.0 (full risk attributed to user). WagerzonJ is a
# partner account where the user holds 87.5% of the risk; raw bets
# also land on the Shared tab for week-over-week reconciliation.
ACCOUNTS = {
    'default': {
        'username_env': 'WAGERZON_USERNAME',
        'password_env': 'WAGERZON_PASSWORD',
        'platform': 'Wagerzon',
        'bet_multiplier': 1.0,
        'shared_sheet': None,
    },
    'j': {
        'username_env': 'WAGERZONJ_USERNAME',
        'password_env': 'WAGERZONJ_PASSWORD',
        'platform': 'WagerzonJ',
        'bet_multiplier': 0.875,
        'shared_sheet': 'Shared',
    },
}
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

    We strip the rotation number and clean HTML, but keep unicode fractions
    (½ ¼ ¾) to match the existing sheet format.
    """
    desc = re.sub(r'<[^>]+>', ' ', raw)  # Strip all HTML tags (<BR>, <em>, etc.)
    desc = re.sub(r'^\[\d+\]\s*', '', desc)
    desc = re.sub(r'\s+', ' ', desc).strip()
    return desc


def _to_decimal(text: str) -> str:
    """Convert unicode fractions to decimal for numeric parsing only."""
    return text.replace('\u00bd', '.5').replace('\u00bc', '.25').replace('\u00be', '.75')


def parse_line_from_desc(desc: str) -> str:
    """Extract line/spread from a cleaned description."""
    desc = _to_decimal(desc)  # Convert ½ → .5 for numeric matching
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
    desc = _to_decimal(desc)  # Convert ½ → .5 for numeric matching
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
    return mapping.get(id_sport, '')


def map_result(result_str: str) -> str:
    """Map API result to standard format."""
    r = result_str.upper()
    if r == 'WIN':
        return 'win'
    elif r == 'LOSE':
        return 'loss'
    elif r in ('PUSH', 'CANCELLED', 'VOID', 'NO ACTION', 'NO BET'):
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


def parse_api_bets(history: dict, platform: str = 'Wagerzon',
                   bet_multiplier: float = 1.0) -> list:
    """Convert HistoryHelper JSON to the standard bet dict format.

    The JSON has a 'details' list where each item is a day with a
    list of wagers. Each wager has legs in its own 'details' list.

    Args:
        history: The 'result' dict from HistoryHelper.aspx.
        platform: Sheet platform label (e.g. 'Wagerzon', 'WagerzonJ').
        bet_multiplier: Fraction of the wagered amount attributable to
                        the user. 1.0 for the primary account; 0.875
                        for the partner-shared WagerzonJ account.
                        The original risk is preserved in '_raw_risk'
                        so it can be uploaded to a verification tab.
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
                win_potential = float(wager.get('WinAmount', 0))
                win_loss = float(wager.get('WinLoss', 0))
                result = map_result(wager.get('Result', ''))

                # Skip pending/open bets
                if not result:
                    continue

                # Use the actual settlement (WinLoss) to calculate odds.
                # WinLoss reflects the real payout — e.g. parlays with pushed
                # legs pay reduced odds, and WinLoss captures that.
                if result == 'win' and risk > 0 and win_loss > 0:
                    american_odds = calculate_american_odds(risk, win_loss)
                elif risk > 0 and win_potential > 0:
                    # For losses and pushes, use potential win to show the
                    # odds when the bet was placed.
                    american_odds = calculate_american_odds(risk, win_potential)
                else:
                    american_odds = 0
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

                # Prepend HeaderDesc to match old format:
                # "PARLAY (2 TEAMS) | leg1 | leg2" or "STRAIGHT BET | desc"
                header = wager.get('HeaderDesc', '').strip()
                leg_text = ' | '.join(leg_descs)
                description = f"{header} | {leg_text}" if header else leg_text
                bet_type = parse_bet_type(header)

                # Apply per-account stake adjustment.
                # bet_amount = user's share for Sheet1; _raw_risk = original
                # for the Shared verification tab.
                adjusted_risk = round(risk * bet_multiplier, 2)

                bet = {
                    'date': parse_date(wager.get('PlacedDate', '')),
                    'platform': platform,
                    'sport': sport,
                    'description': description,
                    'bet_type': bet_type,
                    'line': line,
                    'odds': american_odds,
                    'bet_amount': adjusted_risk,
                    'dec': decimal_odds,
                    'result': result,
                    '_raw_risk': risk,
                }

                bets.append(bet)
                print(f"  {bet['date']} - {bet_type} - {description[:60]} - ${adjusted_risk:.2f} (raw ${risk:.2f}) - {result}")

            except Exception as e:
                print(f"  Error parsing wager: {e}")
                continue

    return bets


# ── Main ────────────────────────────────────────────────────────


def scrape_wagerzon(weeks_back: int = 1, all_weeks: bool = False) -> list:
    """
    Log into Wagerzon via HTTP and fetch bet history from the JSON API.

    Uses HistoryHelper.aspx — the same endpoint the React frontend calls.
    Returns structured data with sport IDs, so no HTML parsing or team
    name guessing needed.

    Args:
        weeks_back: Which week to fetch (0 = current, 1 = last week, etc.)
        all_weeks: If True, fetch weeks 0 through N until an empty week is
                   found. Useful for catching up after missed weeks.

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

    # Step 2: Determine which weeks to fetch
    if all_weeks:
        weeks_to_fetch = range(0, 10)  # cap at 10 weeks back
    else:
        weeks_to_fetch = [weeks_back]

    # Step 3: Fetch and parse each week
    all_bets = []
    for week in weeks_to_fetch:
        week_label = 'This Week' if week == 0 else f'Week {week}'
        print(f"Fetching bet history ({week_label})...")
        history = fetch_history_json(session, week=week)

        days = history.get('details', [])
        total_wagers = sum(len(d.get('wager', [])) for d in days)
        print(f"  {total_wagers} wagers across {len(days)} days")

        if total_wagers == 0 and all_weeks and week > 0:
            print(f"  Empty week — stopping.\n")
            break

        bets = parse_api_bets(history)
        all_bets.extend(bets)
        print()

    return all_bets


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from Wagerzon')
    parser.add_argument('--weeks', type=int, default=1,
                        help='Weeks back to fetch (0=This Week, 1=Last Week, default: 1)')
    parser.add_argument('--all-weeks', action='store_true',
                        help='Fetch all weeks until an empty one (catches up after missed weeks)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Scrape but do not upload to Google Sheets')
    args = parser.parse_args()

    print("=" * 60)
    print("WAGERZON BET HISTORY SCRAPER")
    print("=" * 60)

    try:
        bets = scrape_wagerzon(weeks_back=args.weeks, all_weeks=args.all_weeks)
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
