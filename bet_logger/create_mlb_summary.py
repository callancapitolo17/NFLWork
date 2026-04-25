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


# ── Formula builders ──────────────────────────────────────────────────────
# Column map in Sheet1:
#   A=Date  B=Platform  C=Sport  D=Description  E=Bet Type
#   F=Line  G=Odds      H=Stake  I=Dec Odds     J=Result

S = f"'{SRC}'"
COL_DATE = f"{S}!A$2:A${N}"
COL_SPORT = f"{S}!C$2:C${N}"
COL_DESC = f"{S}!D$2:D${N}"
COL_TYPE = f"{S}!E$2:E${N}"
COL_STAKE = f"{S}!H$2:H${N}"
COL_DEC = f"{S}!I$2:I${N}"
COL_RESULT = f"{S}!J$2:J${N}"

# Filter mask as a Sheets expression (no leading `=`). Five conditions,
# multiplied together (boolean AND). Leaves a 1/0 array suitable for
# SUMPRODUCT. Mirrors is_mlb_correlated_parlay() exactly.
FILTER_MASK = (
    f'({COL_SPORT}="MLB")'
    f'*({COL_TYPE}="Parlay")'
    f'*(LEFT({COL_DESC},16)="PARLAY (2 TEAMS)")'
    f'*IFERROR(REGEXMATCH({COL_DESC},'
    r'"(?<![\d.])[+\-]\d+\.5(?!\d)"'
    f'),FALSE)'
    f'*IFERROR(REGEXMATCH({COL_DESC},'
    r'"\b(?i)(Over|Under)\s+\d+(?:\.\d+)?\b"'
    f'),FALSE)'
)

# The FG/F5 modifiers — multiplied INTO the filter mask to split the set.
F5_COND = (
    f'IFERROR(REGEXMATCH({COL_DESC},'
    r'"\b(?i)(?:1st\s*5|F5|First\s*5)\b"'
    f'),FALSE)'
)
FG_ADD = f"*(1-{F5_COND})"   # NOT F5
F5_ADD = f"*{F5_COND}"       # IS F5

SETTLED = f'({COL_RESULT}<>"")'
WIN = f'({COL_RESULT}="win")'
LOSS = f'({COL_RESULT}="loss")'
PUSH = f'({COL_RESULT}="push")'


def placed_f(mask=FILTER_MASK, extra=""):
    """Count of qualifying bets (settled + pending)."""
    return f"=SUMPRODUCT({mask}{extra})"


def settled_f(mask=FILTER_MASK, extra=""):
    """Count of qualifying settled bets."""
    return f"=SUMPRODUCT({mask}{extra}*{SETTLED})"


def wagered_f(mask=FILTER_MASK, extra=""):
    """Sum of stake across settled qualifying bets."""
    return f"=SUMPRODUCT({mask}{extra}*{SETTLED}*{COL_STAKE})"


def _profit_expr(mask=FILTER_MASK, extra=""):
    """Profit expression WITHOUT leading `=` (for embedding)."""
    wins_payout = (
        f"SUMPRODUCT({mask}{extra}*{WIN}*({COL_DEC}<>\"\")"
        f"*{COL_STAKE}*({COL_DEC}-1))"
    )
    losses = f"SUMPRODUCT({mask}{extra}*{LOSS}*{COL_STAKE})"
    return f"{wins_payout}-{losses}"


def profit_f(mask=FILTER_MASK, extra=""):
    return f"={_profit_expr(mask, extra)}"


def wins_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{WIN})"


def losses_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{LOSS})"


def pushes_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{PUSH})"


def record_f(mask=FILTER_MASK, extra=""):
    w = f"SUMPRODUCT({mask}{extra}*{WIN})"
    l = f"SUMPRODUCT({mask}{extra}*{LOSS})"
    p = f"SUMPRODUCT({mask}{extra}*{PUSH})"
    return f'={w}&"-"&{l}&"-"&{p}'


def avg_odds_f(mask=FILTER_MASK, extra=""):
    num = f"SUMPRODUCT({mask}{extra}*{SETTLED}*{COL_DEC})"
    den = f"SUMPRODUCT({mask}{extra}*{SETTLED})"
    return f"=IF({den}=0,0,{num}/{den})"


def main():
    # Filled in by later tasks
    raise NotImplementedError("main() not implemented yet — see later tasks")


if __name__ == "__main__":
    main()
