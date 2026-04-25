#!/usr/bin/env python3
"""
Create the MLB Summary tab in the bet-logging Google Sheet.

One-shot setup script. Writes pure spreadsheet formulas that auto-update
from Sheet1 as new Wagerzon bets are scraped in. Also creates two embedded
charts (equity curve + weekly P&L bars). No helper sheet. No cron.

Filter: sport=MLB + bet_type=Parlay + description starts with
"PARLAY (2 TEAMS)" + description has a spread token + has a total token.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sheets import (  # noqa: E402
    get_sheets_service,
    SPREADSHEET_ID,
    SHEET_NAME,
    ensure_sheet_exists,
)
from googleapiclient.errors import HttpError  # noqa: E402


# ── Constants ─────────────────────────────────────────────────────────────
TAB = "MLB Summary"
SRC = SHEET_NAME  # "Sheet1"
N = 10000         # Max row the formulas scan in Sheet1


# ── Filter regex ──────────────────────────────────────────────────────────
# A spread leg contains a signed half-integer like "+1.5" or "-2.5". The
# `.5` requirement excludes American odds (which are integers like -110,
# +140 — never half-numbers). Lookarounds prevent matching inside numbers
# that happen to contain ".5" digits — e.g., we don't want to match the
# "+1.5" inside hypothetical "+1.50" odds, though Wagerzon doesn't use that.
SPREAD_REGEX = re.compile(r"(?<![\d.])[+\-]\d+\.5(?!\d)")

# A total leg contains "Over N.N" or "Under N.N". We accept the word forms
# only (the bare "O 8.5" / "U 8.5" shorthand does not appear in Wagerzon's
# MLB output — confirmed against the scraper's clean_desc output).
TOTAL_REGEX = re.compile(r"\b(Over|Under)\s+\d+(?:\.\d+)?\b", re.IGNORECASE)

# F5 markers in leg descriptions.
F5_REGEX = re.compile(r"\b(?:1st\s*5|F5|First\s*5)\b", re.IGNORECASE)


def is_mlb_correlated_parlay(description: str, sport: str, bet_type: str) -> bool:
    """True iff this Sheet1 row is a tracked MLB correlated parlay."""
    if sport != "MLB":
        return False
    if bet_type != "Parlay":
        return False
    if not description.startswith("PARLAY (2 TEAMS)"):
        return False
    if not SPREAD_REGEX.search(description):
        return False
    if not TOTAL_REGEX.search(description):
        return False
    return True


def is_f5_parlay(description: str) -> bool:
    """True iff this parlay is a First-5-Innings parlay."""
    return bool(F5_REGEX.search(description))


def main():
    # Filled in by later tasks
    raise NotImplementedError("main() not implemented yet — see later tasks")


if __name__ == "__main__":
    main()
