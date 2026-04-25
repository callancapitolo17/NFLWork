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

# Cap on the dynamic trend tables — bounds the MMULT triangle so cumulative
# formulas don't materialize a ~1000x1000 matrix on every recalc.
DAILY_ROWS_MAX = 2000
WEEKLY_ROWS_MAX = 500


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
FG_ADD = f"*NOT({F5_COND})"  # NOT F5
F5_ADD = f"*{F5_COND}"       # IS F5

SETTLED = f'({COL_RESULT}<>"")'
WIN = f'({COL_RESULT}="win")'
LOSS = f'({COL_RESULT}="loss")'
PUSH = f'({COL_RESULT}="push")'


# All builders below take `(mask=FILTER_MASK, extra="")`. By convention
# `extra` is a string that begins with `*` (e.g., `FG_ADD`, `F5_ADD`,
# or per-row date conditions) — it gets concatenated into the SUMPRODUCT
# expression as an additional multiplicative term. Passing an `extra`
# without a leading `*` will produce a malformed Sheets expression.
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


# Internal "expression" variants — return SUMPRODUCT(...) WITHOUT leading `=`
# so they can be embedded inside compound formulas (win rate, record, etc.).

def _wins_expr(mask=FILTER_MASK, extra=""):
    return f"SUMPRODUCT({mask}{extra}*{WIN})"


def _losses_expr(mask=FILTER_MASK, extra=""):
    return f"SUMPRODUCT({mask}{extra}*{LOSS})"


def _pushes_expr(mask=FILTER_MASK, extra=""):
    return f"SUMPRODUCT({mask}{extra}*{PUSH})"


def winrate_f(mask=FILTER_MASK, extra=""):
    """Win rate = wins / (wins + losses). Pushes excluded from denominator."""
    w = _wins_expr(mask, extra)
    lo = _losses_expr(mask, extra)
    return f"=IF(({w}+{lo})=0,0,{w}/({w}+{lo}))"


def record_f(mask=FILTER_MASK, extra=""):
    w = _wins_expr(mask, extra)
    lo = _losses_expr(mask, extra)
    p = _pushes_expr(mask, extra)
    return f'={w}&"-"&{lo}&"-"&{p}'


def avg_odds_f(mask=FILTER_MASK, extra=""):
    num = f"SUMPRODUCT({mask}{extra}*{SETTLED}*{COL_DEC})"
    den = f"SUMPRODUCT({mask}{extra}*{SETTLED})"
    return f"=IF({den}=0,0,{num}/{den})"


# ── Row builder ───────────────────────────────────────────────────────────

# Plain Python lists do not allow attribute assignment (e.g.
# `lst.foo = 1` raises AttributeError). We need to attach anchor row
# numbers (`__anchor_daily_start__`, `__anchor_weekly_start__`) to the
# return value of build_rows() so downstream tasks (charts, wiring) can
# locate the dynamic blocks without re-deriving their positions. A trivial
# subclass of `list` accepts attributes on its instances while behaving
# identically to a list everywhere else.
class _RowList(list):
    """list subclass that permits arbitrary attribute assignment."""
    pass


def build_rows():
    """Return the 2D list of cell values for the MLB Summary tab.

    Layout (row indices 1-based):
      1  Title
      2  (blank)
      3  "OVERALL STATS" header
      4  column labels
      5-13 overall metrics (anchor row numbers used by ROI/win-rate formulas)
      14 (blank)
      15 "BY GAME WINDOW" header
      16 column labels
      17 FG row
      18 F5 row
      19 (blank)
      20 "EQUITY CURVE (Daily P&L)" header
      21 column labels ("Date","Bets","Wagered","P&L","Cumulative P&L")
      22+ daily data rows (dynamic, auto-extend)
      ...
      (a later section for weekly — filled by Task 5)
    """
    rows = _RowList()

    # Row 1: title
    rows.append(["MLB CORRELATED PARLAYS — WAGERZON"])
    # Row 2: blank
    rows.append([])

    # Row 3: section header
    rows.append(["OVERALL STATS"])
    # Row 4: column labels
    rows.append(["", "Value"])
    # Row 5: Total bets placed
    rows.append(["Total bets placed", placed_f()])
    # Row 6: Settled bets
    rows.append(["Settled bets (W+L+P)", settled_f()])
    # Row 7: Total wagered
    rows.append(["Total wagered", wagered_f()])
    # Row 8: Net P&L
    rows.append(["Net P&L", profit_f()])
    # Row 9: ROI = Net P&L / Total wagered (anchors derived from current
    # row positions so inserting/removing metrics above won't silently break)
    r_wagered = len(rows) - 1   # the just-appended "Total wagered" row
    r_pnl = len(rows)           # the just-appended "Net P&L" row
    rows.append(["ROI", f"=IF(B{r_wagered}=0,0,B{r_pnl}/B{r_wagered})"])
    # Row 10: Win rate = wins / (wins + losses), NOT wins / settled
    rows.append(["Win rate", winrate_f()])
    # Row 11: Record
    rows.append(["Record (W-L-P)", record_f()])
    # Row 12: Avg decimal odds
    rows.append(["Avg decimal odds", avg_odds_f()])
    # Row 13: blank spacer
    rows.append([])

    # Row 14: "BY GAME WINDOW" header
    rows.append(["BY GAME WINDOW"])
    # Row 15: column labels
    rows.append([
        "Window", "Bets", "Settled", "Wagered",
        "P&L", "ROI", "Win Rate", "Record", "Avg Odds",
    ])
    # Rows 16-17: FG and F5 rows
    for label, extra in (("Full Game (FG)", FG_ADD), ("First 5 (F5)", F5_ADD)):
        r = len(rows) + 1
        rows.append([
            label,
            placed_f(extra=extra),
            settled_f(extra=extra),
            wagered_f(extra=extra),
            profit_f(extra=extra),
            f"=IF(D{r}=0,0,E{r}/D{r})",     # ROI
            winrate_f(extra=extra),
            record_f(extra=extra),
            avg_odds_f(extra=extra),
        ])

    # Row 18: blank spacer (Task 5 continues from here)
    rows.append([])

    # ── Block 3a: Daily P&L (dynamic, feeds the equity curve chart) ──
    rows.append(["DAILY P&L"])
    rows.append(["Date", "Bets", "Wagered", "P&L", "Cumulative P&L"])
    daily_header_row = len(rows)          # 1-based row of column labels
    daily_start_row = daily_header_row + 1
    daily_end = daily_start_row + DAILY_ROWS_MAX - 1

    # Column A: a single dynamic formula. The spilled dates fill downward.
    # We put the formula in the top cell; rows below it will be filled by
    # the ARRAYFORMULA spill when the sheet is opened.
    rows.append([
        (
            f"=IFERROR(SORT(UNIQUE(FILTER({COL_DATE},{FILTER_MASK}=1))),\"\")"
        ),
        # The other four columns use ARRAYFORMULA over the spilled date column
        # in A — bounded to daily_end to avoid materializing a ~1000x1000 MMULT.
        f"=ARRAYFORMULA(IF(A{daily_start_row}:A{daily_end}=\"\",\"\","
        f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A{daily_end})*{FILTER_MASK}*{SETTLED})))",
        f"=ARRAYFORMULA(IF(A{daily_start_row}:A{daily_end}=\"\",\"\","
        f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A{daily_end})*{FILTER_MASK}*{SETTLED}*{COL_STAKE})))",
        # P&L per day — wins payout minus losses stake
        (
            f"=ARRAYFORMULA(IF(A{daily_start_row}:A{daily_end}=\"\",\"\","
            f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A{daily_end})*{FILTER_MASK}*{WIN}"
            f"*({COL_DEC}<>\"\")*{COL_STAKE}*({COL_DEC}-1))"
            f"-SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A{daily_end})*{FILTER_MASK}*{LOSS}*{COL_STAKE})))"
        ),
        # Cumulative P&L — running sum of column D starting from daily_start_row.
        # ARRAYFORMULA of MMULT gives a running sum without per-row formulas.
        # Bounded to daily_end rows to avoid materializing a huge MMULT triangle.
        (
            f"=ARRAYFORMULA(IF(A{daily_start_row}:A{daily_end}=\"\",\"\","
            f"MMULT(--(ROW(D{daily_start_row}:D{daily_end})>=TRANSPOSE(ROW(D{daily_start_row}:D{daily_end}))),"
            f"IFERROR(D{daily_start_row}:D{daily_end}*1,0))))"
        ),
    ])
    rows.append([])  # spacer after daily block

    # ── Block 3b: Weekly P&L (dynamic, feeds the bar chart) ──
    rows.append(["WEEKLY P&L"])
    rows.append(["Week of (Mon)", "Bets", "Wagered", "P&L", "Cumulative P&L"])
    weekly_header_row = len(rows)
    weekly_start_row = weekly_header_row + 1
    weekly_end = weekly_start_row + WEEKLY_ROWS_MAX - 1

    # Week-start bucket = A - WEEKDAY(A, 2) + 1  (Monday of that date's week).
    # Wrap the FILTER in ARRAYFORMULA so the transform is applied element-wise.
    week_of_date = (
        f"({COL_DATE}-WEEKDAY({COL_DATE},2)+1)"
    )
    rows.append([
        (
            f"=IFERROR(SORT(UNIQUE("
            f"FILTER(ARRAYFORMULA({week_of_date}),{FILTER_MASK}=1))),\"\")"
        ),
        f"=ARRAYFORMULA(IF(A{weekly_start_row}:A{weekly_end}=\"\",\"\","
        f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A{weekly_end})*{FILTER_MASK}*{SETTLED})))",
        f"=ARRAYFORMULA(IF(A{weekly_start_row}:A{weekly_end}=\"\",\"\","
        f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A{weekly_end})*{FILTER_MASK}*{SETTLED}*{COL_STAKE})))",
        (
            f"=ARRAYFORMULA(IF(A{weekly_start_row}:A{weekly_end}=\"\",\"\","
            f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A{weekly_end})*{FILTER_MASK}*{WIN}"
            f"*({COL_DEC}<>\"\")*{COL_STAKE}*({COL_DEC}-1))"
            f"-SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A{weekly_end})*{FILTER_MASK}*{LOSS}*{COL_STAKE})))"
        ),
        (
            f"=ARRAYFORMULA(IF(A{weekly_start_row}:A{weekly_end}=\"\",\"\","
            f"MMULT(--(ROW(D{weekly_start_row}:D{weekly_end})>=TRANSPOSE(ROW(D{weekly_start_row}:D{weekly_end}))),"
            f"IFERROR(D{weekly_start_row}:D{weekly_end}*1,0))))"
        ),
    ])

    # Also stash the anchor row numbers so Task 6 (charts) can reference
    # them without re-computing. Attach as an attribute on the returned list.
    rows._anchor_daily_start = daily_start_row
    rows._anchor_weekly_start = weekly_start_row

    return rows


# ── Chart builders ────────────────────────────────────────────────────────

def get_sheet_id(service, tab_name):
    meta = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
    for s in meta.get("sheets", []):
        if s["properties"]["title"] == tab_name:
            return s["properties"]["sheetId"]
    return None


def build_chart_requests(sheet_id, daily_start_row, weekly_start_row, last_row):
    """Return Sheets API batchUpdate requests for the two charts.

    Both charts read a generous range (through `last_row`) so they pick up
    new dates automatically as the dynamic tables extend downward.
    """
    # Convert 1-based row numbers to API's 0-based indices.
    # The daily data block occupies columns A-E (indices 0-5 exclusive).
    daily_header_idx = daily_start_row - 2       # 0-based index of the header row
    daily_data_start = daily_start_row - 1       # 0-based index of first data row

    weekly_header_idx = weekly_start_row - 2
    weekly_data_start = weekly_start_row - 1

    # EQUITY CURVE (line chart) — X: Date (col A), Y: Cumulative P&L (col E)
    equity_chart = {
        "addChart": {
            "chart": {
                "spec": {
                    "title": "Equity Curve — Cumulative P&L",
                    "basicChart": {
                        "chartType": "LINE",
                        "legendPosition": "NO_LEGEND",
                        "axis": [
                            {"position": "BOTTOM_AXIS", "title": "Date"},
                            {"position": "LEFT_AXIS", "title": "Cumulative P&L"},
                        ],
                        "domains": [{
                            "domain": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": daily_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 0, "endColumnIndex": 1,
                            }]}}
                        }],
                        "series": [{
                            "series": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": daily_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 4, "endColumnIndex": 5,
                            }]}},
                            "targetAxis": "LEFT_AXIS",
                        }],
                        "headerCount": 0,
                    },
                },
                "position": {
                    "overlayPosition": {
                        "anchorCell": {
                            "sheetId": sheet_id,
                            "rowIndex": daily_header_idx - 12,  # above DAILY P&L block
                            "columnIndex": 6,                   # column G
                        },
                        "widthPixels": 600,
                        "heightPixels": 300,
                    }
                },
            }
        }
    }

    # WEEKLY BARS — X: Week of (col A), Y: P&L (col D)
    weekly_chart = {
        "addChart": {
            "chart": {
                "spec": {
                    "title": "Weekly P&L",
                    "basicChart": {
                        "chartType": "COLUMN",
                        "legendPosition": "NO_LEGEND",
                        "axis": [
                            {"position": "BOTTOM_AXIS", "title": "Week of"},
                            {"position": "LEFT_AXIS", "title": "P&L"},
                        ],
                        "domains": [{
                            "domain": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": weekly_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 0, "endColumnIndex": 1,
                            }]}}
                        }],
                        "series": [{
                            "series": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": weekly_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 3, "endColumnIndex": 4,
                            }]}},
                            "targetAxis": "LEFT_AXIS",
                        }],
                        "headerCount": 0,
                    },
                },
                "position": {
                    "overlayPosition": {
                        "anchorCell": {
                            "sheetId": sheet_id,
                            "rowIndex": weekly_header_idx - 6,
                            "columnIndex": 6,                   # column G, below equity curve
                        },
                        "widthPixels": 600,
                        "heightPixels": 300,
                    }
                },
            }
        }
    }

    return [equity_chart, weekly_chart]


def main():
    # Filled in by later tasks
    raise NotImplementedError("main() not implemented yet — see later tasks")


if __name__ == "__main__":
    main()
