#!/usr/bin/env python3
"""
BetOnline Bet History -> Google Sheet

Reads BetOnline bets from the Unabated Ticket bets service on the Oracle VM
(GET <BETS_SERVICE_URL>/bets.json, over Tailscale) instead of logging in to
BetOnline. Since the 2026-10-14 phone cutover the BetOnline refresh token lives
only on the VM; a second login from this Mac would rotate it under the VM and
Keycloak's reuse detection would kill both, so this script holds no BetOnline
credentials and never calls BetOnline.

Inputs:  BETS_SERVICE_URL (bet_logger/.env, default DEFAULT_BETS_SERVICE_URL);
         the VM stores each report row verbatim under record["raw"]
         (unabated_ticket/bets_service/sources/betonline.py), which is rebuilt
         into the report's own shape and parsed exactly as before.
Outputs: settled bets appended to the Google Sheet (sheets.append_bets_to_sheet,
         which skips rows already there). Pending bets are left out, as before.
Fails loudly when the VM is unreachable or its last BetOnline poll failed, so
the Monday runner reports BetOnline as FAILED rather than logging stale data.
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

SCRIPT_DIR = os.path.dirname(__file__)

# The VM's tailnet name (tailscale serve fronts the bets service with HTTPS).
DEFAULT_BETS_SERVICE_URL = "https://mlb-stack.tail87ab2d.ts.net"
BETS_SERVICE_SOURCE = "betonline"
BETS_SERVICE_TIMEOUT_SEC = 30
# /bets.json windows on the settle time; a bet settles on or after the day it was
# placed, so a window this much wider than the placed-date range holds every bet.
SETTLE_WINDOW_MARGIN_DAYS = 2
# The report fields the sheet parser reads, as the VM keeps them under record["raw"].
RAW_TO_REPORT_FIELDS = {
    "description": "Description",
    "wagerType": "WagerType",
    "wagerStatus": "WagerStatus",
    "risk": "Risk",
    "toWin": "ToWin",
    "date": "Date",
}


def fetch_bets_payload(service_url: str, days: int) -> dict:
    """GET /bets.json?days=N from the bets service. Raises with the reason when the
    VM cannot be reached, answers non-200, or its last BetOnline poll failed."""
    url = f"{service_url.rstrip('/')}/bets.json"
    try:
        resp = requests.get(url, params={"days": days}, timeout=BETS_SERVICE_TIMEOUT_SEC)
    except requests.RequestException as e:
        raise RuntimeError(
            f"could not reach the bets service at {url} ({e}). "
            "Is Tailscale up on this Mac?") from e
    if resp.status_code != 200:
        raise RuntimeError(f"bets service {url} answered HTTP {resp.status_code}: {resp.text[:200]}")
    payload = resp.json()
    status = payload.get("sources", {}).get(BETS_SERVICE_SOURCE)
    if status is None:
        raise RuntimeError(f"bets service has no '{BETS_SERVICE_SOURCE}' source; sources seen: "
                           f"{sorted(payload.get('sources', {}))}")
    if not status.get("ok"):
        raise RuntimeError(
            f"the VM's last BetOnline poll failed ({status.get('error')}); last good poll "
            f"{status.get('fetchedAt')}. Its BetOnline token may need a new login on the VM.")
    return payload


def report_rows_from_payload(payload: dict, from_date: str) -> list:
    """The BetOnline report rows placed on or after from_date (YYYY-MM-DD, the
    report's own clock), one per wager. A parlay is stored as one record per leg,
    each carrying the same raw row, so rows are de-duplicated on the wager id."""
    rows_by_wager_id = {}
    for record in payload.get("bets", []):
        if record.get("venue") != BETS_SERVICE_SOURCE:
            continue
        raw = record.get("raw") or {}
        wager_id = raw.get("nativeId")
        if not wager_id:
            raise RuntimeError(f"BetOnline record {record.get('id')} carries no raw.nativeId")
        placed_date = str(raw.get("date") or "")[:10]
        if not placed_date:
            raise RuntimeError(f"BetOnline record {record.get('id')} carries no raw.date")
        if placed_date < from_date:
            continue
        rows_by_wager_id[wager_id] = {
            report_field: raw.get(raw_field) for raw_field, report_field in RAW_TO_REPORT_FIELDS.items()
        }
    return sorted(rows_by_wager_id.values(), key=lambda row: str(row["Date"] or ""))


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


def scrape_betonline(days_back: int = 7, from_date: str = None,
                     service_url: str = DEFAULT_BETS_SERVICE_URL) -> list:
    """
    Settled BetOnline bets placed on or after from_date (else the last days_back
    days), read from the bets service on the VM. Read-only: no BetOnline call,
    no token, no DB write.
    """
    today = datetime.now()
    if from_date:
        start = from_date
    else:
        start = (today - timedelta(days=days_back)).strftime("%Y-%m-%d")
    days_in_range = (today - datetime.strptime(start, "%Y-%m-%d")).days
    window_days = max(days_in_range, 0) + SETTLE_WINDOW_MARGIN_DAYS

    print(f"Reading BetOnline bets placed since {start} from {service_url}...")
    payload = fetch_bets_payload(service_url, window_days)
    report_rows = report_rows_from_payload(payload, start)

    print(f"\nBets service returned {len(report_rows)} BetOnline bets")
    return parse_api_bets(report_rows)


def scrape_from_file(filepath: str) -> list:
    """Parse bets from a saved API JSON response (for testing)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    api_bets = data if isinstance(data, list) else data.get('Data', [])
    return parse_api_bets(api_bets)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Scrape bet history from BetOnline')
    parser.add_argument('--test', action='store_true', help='Parse from a saved report JSON file instead of the bets service')
    parser.add_argument('--days', type=int, default=7, help='Days of history to fetch (default: 7)')
    parser.add_argument('--since-last', action='store_true', help='Scrape from the date of the last BetOnline bet in Google Sheets')
    parser.add_argument('--dry-run', action='store_true', help='Scrape but do not upload to Google Sheets')
    args = parser.parse_args()
    service_url = os.getenv('BETS_SERVICE_URL') or DEFAULT_BETS_SERVICE_URL

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

            bets = scrape_betonline(days_back=args.days, from_date=from_date, service_url=service_url)

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
