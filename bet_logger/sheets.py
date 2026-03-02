#!/usr/bin/env python3
"""
Google Sheets integration for bet logging.
Uploads scraped bet data to the specified Google Sheet.
"""

import os
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from dotenv import load_dotenv

load_dotenv()

# Google Sheets configuration
SPREADSHEET_ID = os.getenv("GOOGLE_SHEET_ID", "1t9_7HmsrQAu34HI_gMIDPErC2kThz2MHBwwzZIYU6H4")
CREDENTIALS_FILE = os.getenv("GOOGLE_CREDENTIALS_FILE", "credentials.json")
SCOPES = ['https://www.googleapis.com/auth/spreadsheets']

# Sheet configuration - adjust these based on your sheet structure
SHEET_NAME = os.getenv("SHEET_NAME", "Sheet1")  # Name of the sheet tab
START_ROW = int(os.getenv("START_ROW", "383"))  # First row to write to

# Column mapping (A=1, B=2, etc.)
# Your columns: Date, Platform, Sports, Bet Description, Bet Type, Line, Odds, Bet Amount, Dec (auto), Result
COLUMN_ORDER = ['date', 'platform', 'sport', 'description', 'bet_type', 'line', 'odds', 'bet_amount', 'dec', 'result']


def get_sheets_service():
    """Initialize and return Google Sheets API service."""
    if not os.path.exists(CREDENTIALS_FILE):
        raise FileNotFoundError(
            f"Credentials file '{CREDENTIALS_FILE}' not found. "
            "Please download your service account credentials from Google Cloud Console."
        )

    creds = Credentials.from_service_account_file(CREDENTIALS_FILE, scopes=SCOPES)
    service = build('sheets', 'v4', credentials=creds)
    return service


def get_next_empty_row(service, sheet_name: str = None) -> int:
    """Find the next empty row in the sheet."""
    tab = sheet_name or SHEET_NAME
    try:
        # Get all values in column A to find last row with data
        result = service.spreadsheets().values().get(
            spreadsheetId=SPREADSHEET_ID,
            range=f"{tab}!A:A"
        ).execute()

        values = result.get('values', [])
        return len(values) + 1

    except HttpError as e:
        print(f"Error finding next row: {e}")
        return START_ROW


def _normalize_date(date_str: str) -> str:
    """Normalize date string to M/D/YYYY for consistent comparison."""
    from datetime import datetime as dt
    date_str = date_str.strip()
    for fmt in ("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d"):
        try:
            parsed = dt.strptime(date_str, fmt)
            return parsed.strftime("%-m/%-d/%Y")
        except ValueError:
            continue
    return date_str


def _normalize_amount(amount_str: str) -> str:
    """Normalize bet amount to plain float string for consistent comparison."""
    amount_str = amount_str.strip().replace("$", "").replace(",", "")
    try:
        return f"{float(amount_str):.2f}"
    except (ValueError, TypeError):
        return amount_str


def get_existing_bets(service, sheet_name: str = None) -> dict:
    """
    Get existing bets from sheet for duplicate detection.
    Returns a dict mapping (date, platform, description, amount) -> count,
    so legitimate duplicate bets (same game placed twice) are tracked.
    """
    tab = sheet_name or SHEET_NAME
    try:
        result = service.spreadsheets().values().get(
            spreadsheetId=SPREADSHEET_ID,
            range=f"{tab}!A:H"
        ).execute()

        values = result.get('values', [])
        existing = {}

        for row in values[1:]:  # Skip header row
            if len(row) >= 8:
                date = _normalize_date(row[0]) if row[0] else ''
                platform = row[1].strip() if len(row) > 1 and row[1] else ''
                description = row[3].strip() if len(row) > 3 and row[3] else ''
                bet_amount = _normalize_amount(row[7]) if len(row) > 7 and row[7] else ''

                if date and description and bet_amount:
                    key = (date, platform, description, bet_amount)
                    existing[key] = existing.get(key, 0) + 1

        return existing

    except HttpError as e:
        print(f"Error fetching existing bets: {e}")
        return {}


def filter_duplicates(bets: list, existing: dict) -> list:
    """
    Filter out bets that already exist in the sheet.
    Uses count-based tracking so legitimate duplicate bets
    (same game placed twice) aren't wrongly filtered.
    """
    new_bets = []
    # Track how many times each key appears in the incoming batch
    seen_counts = {}

    for bet in bets:
        date = _normalize_date(bet.get('date', ''))
        platform = bet.get('platform', '')
        description = bet.get('description', '').strip()
        bet_amount = _normalize_amount(f"{bet.get('bet_amount', 0):.2f}") if bet.get('bet_amount') else ''

        key = (date, platform, description, bet_amount)
        seen_counts[key] = seen_counts.get(key, 0) + 1

        existing_count = existing.get(key, 0)

        if seen_counts[key] <= existing_count:
            print(f"  Skipping duplicate: {date} - {description[:50]}...")
        else:
            new_bets.append(bet)

    return new_bets


def format_bet_row(bet: dict) -> list:
    """Convert bet dictionary to row values matching spreadsheet columns."""
    row = []
    for col in COLUMN_ORDER:
        value = bet.get(col, '')

        # Format specific columns
        if col == 'odds' and value is not None:
            # Format American odds with +/- prefix
            if isinstance(value, (int, float)):
                value = f"+{int(value)}" if value > 0 else str(int(value))
        elif col == 'bet_amount' and value is not None:
            # Format as currency
            value = f"${value:.2f}"
        elif col == 'dec' and value is not None:
            # Pass raw number - let Google Sheets handle display formatting
            pass
        elif value is None:
            value = ''

        row.append(value)

    return row


def ensure_sheet_exists(service, sheet_name: str):
    """Create a sheet tab if it doesn't already exist."""
    try:
        meta = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
        existing = [s['properties']['title'] for s in meta.get('sheets', [])]
        if sheet_name in existing:
            return
        service.spreadsheets().batchUpdate(
            spreadsheetId=SPREADSHEET_ID,
            body={'requests': [{'addSheet': {'properties': {'title': sheet_name}}}]}
        ).execute()
        print(f"Created new sheet tab: {sheet_name}")
    except HttpError as e:
        print(f"Warning: could not ensure sheet '{sheet_name}' exists: {e}")


def append_bets_to_sheet(bets: list, start_row: int = None, sheet_name: str = None) -> dict:
    """
    Append bet data to Google Sheet.

    Args:
        bets: List of bet dictionaries from scraper
        start_row: Optional specific row to start at (otherwise finds next empty)

    Returns:
        dict with status and details
    """
    if not bets:
        return {"status": "skipped", "message": "No bets to upload"}

    tab = sheet_name or SHEET_NAME

    try:
        service = get_sheets_service()

        # Create tab if it doesn't exist (only for non-default tabs)
        if sheet_name:
            ensure_sheet_exists(service, sheet_name)

        # Check for duplicates
        print("Checking for duplicates...")
        existing = get_existing_bets(service, sheet_name=tab)
        print(f"Found {len(existing)} existing bets in {tab}")

        original_count = len(bets)
        bets = filter_duplicates(bets, existing)

        if not bets:
            print(f"All {original_count} bets already exist in {tab}. Nothing to upload.")
            return {"status": "skipped", "message": "All bets already exist", "duplicates": original_count}

        print(f"Uploading {len(bets)} new bets ({original_count - len(bets)} duplicates skipped)")

        # Find next empty row if not specified
        if start_row is None:
            start_row = get_next_empty_row(service, sheet_name=tab)

        # Format all bets as rows
        rows = [format_bet_row(bet) for bet in bets]

        # Calculate range (A through I for 9 columns)
        end_col = chr(ord('A') + len(COLUMN_ORDER) - 1)  # 'I' for 9 columns
        range_str = f"{tab}!A{start_row}:{end_col}{start_row + len(rows) - 1}"

        print(f"Uploading {len(rows)} bets to {range_str}...")

        # Upload to sheet
        body = {'values': rows}
        result = service.spreadsheets().values().update(
            spreadsheetId=SPREADSHEET_ID,
            range=range_str,
            valueInputOption='USER_ENTERED',  # This allows formulas and formatting
            body=body
        ).execute()

        updated_cells = result.get('updatedCells', 0)
        print(f"Successfully updated {updated_cells} cells!")

        return {
            "status": "success",
            "rows_added": len(rows),
            "start_row": start_row,
            "end_row": start_row + len(rows) - 1,
            "cells_updated": updated_cells
        }

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return {"status": "error", "message": str(e)}

    except HttpError as e:
        print(f"Google Sheets API error: {e}")
        return {"status": "error", "message": str(e)}


def test_connection():
    """Test the Google Sheets connection."""
    try:
        service = get_sheets_service()

        # Try to read sheet metadata
        result = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
        title = result.get('properties', {}).get('title', 'Unknown')
        print(f"Successfully connected to sheet: {title}")

        # Get sheet names
        sheets = result.get('sheets', [])
        print(f"Available sheets: {[s['properties']['title'] for s in sheets]}")

        return True

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False

    except HttpError as e:
        print(f"Google Sheets API error: {e}")
        if "403" in str(e):
            print("\nPermission denied. Make sure you've shared the spreadsheet with your service account email.")
        return False


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        print("Testing Google Sheets connection...")
        test_connection()
    else:
        print("Use --test to test the connection")
        print("Or import this module and use append_bets_to_sheet(bets)")
