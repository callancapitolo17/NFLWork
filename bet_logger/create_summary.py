#!/usr/bin/env python3
"""
Create CBB Answer Key Summary tab with auto-updating formulas.
Re-run to add rows for new dates/weeks.
"""

import os, sys, re
from datetime import datetime, timedelta
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheets import get_sheets_service, SPREADSHEET_ID, SHEET_NAME, ensure_sheet_exists
from googleapiclient.errors import HttpError

TAB = "CBB Summary"
HLP = "_cbb_ak"  # Hidden helper sheet (avoids ARRAYFORMULA spill issues)
SRC = SHEET_NAME  # "Sheet1"
N = 5000
# Helper column references (on a separate hidden sheet)
AK = f"'{HLP}'!A$2:A${N}"
MKT = f"'{HLP}'!B$2:B${N}"
S = f"'{SRC}'"


# ── Formula builders ──────────────────────────────────────────────────────

def profit_f(extra=""):
    """SUMPRODUCT profit = wins_payout - losses_wagered."""
    return (
        f'=SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="win")*({S}!I$2:I${N}<>"")*{S}!H$2:H${N}*({S}!I$2:I${N}-1))'
        f'-SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="loss")*{S}!H$2:H${N})'
    )


def wagered_f(extra=""):
    return f'=SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}<>"")*{S}!H$2:H${N})'


def bets_f(extra=""):
    return f'=SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}<>""))'


def wins_f(extra=""):
    return f'=SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="win"))'


def _wins_expr(extra=""):
    """Wins expression WITHOUT leading = (for embedding in other formulas)."""
    return f'SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="win"))'


def _profit_expr(extra=""):
    """Profit expression WITHOUT leading = (for embedding in other formulas)."""
    return (
        f'SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="win")*({S}!I$2:I${N}<>"")*{S}!H$2:H${N}*({S}!I$2:I${N}-1))'
        f'-SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="loss")*{S}!H$2:H${N})'
    )


def record_f(extra=""):
    w = f'SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="win"))'
    l = f'SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="loss"))'
    p = f'SUMPRODUCT({extra}({AK}="AK")*({S}!J$2:J${N}="push"))'
    return f'={w}&"-"&{l}&"-"&{p}'


def row_formulas(label, extra, r):
    """Full breakdown row: label, bets, wagered, profit, roi, winrate, record."""
    return [
        label,
        bets_f(extra),
        wagered_f(extra),
        profit_f(extra),
        f'=IF(C{r}=0,0,D{r}/C{r})',
        f'=IF(B{r}=0,0,{_wins_expr(extra)}/B{r})',
        record_f(extra),
    ]


def plat_extra(platforms):
    """Filter for one or more platforms (OR logic via addition)."""
    if isinstance(platforms, str):
        platforms = [platforms]
    if len(platforms) == 1:
        return f'({S}!B$2:B${N}="{platforms[0]}")*'
    # OR: (B="X")+(B="Y") wrapped in sign() to get 0/1
    parts = "+".join(f'({S}!B$2:B${N}="{p}")' for p in platforms)
    return f'SIGN({parts})*'


def mkt_extra(market_type):
    return f'({MKT}="{market_type}")*'


def date_extra(d):
    return f'({S}!A$2:A${N}=DATE({d.year},{d.month},{d.day}))*'


def week_extra(monday, sunday):
    return (
        f'({S}!A$2:A${N}>=DATE({monday.year},{monday.month},{monday.day}))*'
        f'({S}!A$2:A${N}<=DATE({sunday.year},{sunday.month},{sunday.day}))*'
    )


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    service = get_sheets_service()
    ensure_sheet_exists(service, TAB)
    ensure_sheet_exists(service, HLP)

    # Clear existing content
    for tab in (TAB, HLP):
        try:
            service.spreadsheets().values().clear(
                spreadsheetId=SPREADSHEET_ID, range=f"'{tab}'!A:Z"
            ).execute()
        except HttpError:
            pass

    # Write helper ARRAYFORMULAs FIRST (separate sheet avoids spill interference)
    print("Writing helper sheet...")
    ak_formula = (
        f'=ARRAYFORMULA(IF({S}!A2:A{N}="","",'
        f'IF(({S}!C2:C{N}="NCAAM")*'
        f'IFERROR(REGEXMATCH({S}!D2:D{N},"(?i)(1H|2H|1st Half|2nd Half)"),FALSE),'
        f'"AK","")))'
    )
    mkt_formula = (
        f'=ARRAYFORMULA(IF(A2:A{N}<>"AK","",'
        f'IF(IFERROR(REGEXMATCH({S}!D2:D{N},"(?i)team total"),FALSE),"Team Total",'
        f'IF(IFERROR(REGEXMATCH({S}!D2:D{N},'
        f'"(?i)(/.*( over| under| O | U )|vrs.*(over|under)|TOTAL)"),FALSE),"Total",'
        f'IF(IFERROR(REGEXMATCH({S}!D2:D{N},'
        r'"(?i)(½|[+-]\d*\.5|[+-]\d+\s+-\d+.*for)"),'
        f'FALSE),"Spread",'
        f'"Moneyline")))))'
    )
    service.spreadsheets().values().update(
        spreadsheetId=SPREADSHEET_ID,
        range=f"'{HLP}'!A1:B2",
        valueInputOption="USER_ENTERED",
        body={"values": [["is_ak", "market_type"], [ak_formula, mkt_formula]]},
    ).execute()

    # Read data
    print("Reading data...")
    result = service.spreadsheets().values().get(
        spreadsheetId=SPREADSHEET_ID,
        range=f"'{SRC}'!A2:J",
        valueRenderOption="FORMATTED_VALUE",
    ).execute()
    data = result.get("values", [])

    # Filter to AK bets for discovering unique values
    ak_re = re.compile(r"1H|2H|1st Half|2nd Half", re.IGNORECASE)
    ak_bets = [
        r for r in data
        if len(r) > 3 and r[2].strip() == "NCAAM" and ak_re.search(r[3])
    ]
    print(f"  {len(ak_bets)} answer key bets found")

    # Unique platforms (sorted by frequency), merging BFA + BFAJ
    BFA_NAMES = {"Betfastaction", "BFAJ"}
    raw_counts = Counter(r[1].strip() for r in ak_bets if len(r) > 1)
    # Merge BFA variants into one entry
    merged = Counter()
    for p, c in raw_counts.items():
        merged["BFA" if p in BFA_NAMES else p] += c
    # Build platform list: each entry is (label, [sheet_values])
    platform_map = {}
    for p in raw_counts:
        key = "BFA" if p in BFA_NAMES else p
        platform_map.setdefault(key, []).append(p)
    platforms = [(label, platform_map[label]) for label, _ in merged.most_common()]
    print(f"  Platforms: {[l for l,_ in platforms]}")

    # Find earliest date, then generate every day through end of CBB season
    first_dates = set()
    for r in ak_bets:
        if r and r[0].strip():
            for fmt in ("%m/%d/%Y", "%m/%d/%y"):
                try:
                    first_dates.add(datetime.strptime(r[0].strip(), fmt))
                    break
                except ValueError:
                    continue
    start_date = min(first_dates)
    end_date = datetime(start_date.year, 4, 15)  # Through April 15 (covers March Madness)
    dates = []
    d = start_date
    while d <= end_date:
        dates.append(d)
        d += timedelta(days=1)
    print(f"  Dates: {dates[0].strftime('%m/%d')} - {dates[-1].strftime('%m/%d')} ({len(dates)} days pre-populated)")

    # Weeks (Monday-based) through end date
    weeks = []
    monday = start_date - timedelta(days=start_date.weekday())
    while monday <= end_date:
        weeks.append(monday)
        monday += timedelta(days=7)
    print(f"  Weeks: {len(weeks)} pre-populated")

    # ── Build rows ──
    rows = []

    # Title
    rows.append(["CBB ANSWER KEY PERFORMANCE"])  # 1
    rows.append([])  # 2

    # === OVERALL STATS ===
    rows.append(["OVERALL STATS"])  # 3
    rows.append(["", "Value"])  # 4
    rows.append(["Total Bets", f'=COUNTIF({AK},"AK")'])  # 5  (uses helper sheet)
    rows.append(["Settled Bets", bets_f()])  # 6
    rows.append(["Pending Bets", "=B5-B6"])  # 7
    rows.append(["Total Wagered", wagered_f()])  # 8
    rows.append(["Profit/Loss", profit_f()])  # 9
    rows.append(["ROI", "=IF(B8=0,0,B9/B8)"])  # 10
    rows.append(["Win Rate", f'=IF(B6=0,0,{_wins_expr()}/B6)'])  # 11
    rows.append(["Record (W-L-P)", record_f()])  # 12
    rows.append(["Units +/-", "=B9/100"])  # 13
    rows.append(["Avg Decimal Odds",
                  f'=IF(B6=0,0,SUMPRODUCT(({AK}="AK")*({S}!J$2:J${N}<>"")*{S}!I$2:I${N})/B6)'])  # 14
    rows.append([])  # 15

    # === BY PLATFORM ===
    rows.append(["BY PLATFORM"])  # 16
    rows.append(["Platform", "Bets", "Wagered", "Profit", "ROI", "Win Rate", "Record"])
    for label, sheet_vals in platforms:
        r = len(rows) + 1
        rows.append(row_formulas(label, plat_extra(sheet_vals), r))
    rows.append([])

    # === BY MARKET TYPE ===
    rows.append(["BY MARKET TYPE"])
    rows.append(["Market", "Bets", "Wagered", "Profit", "ROI", "Win Rate", "Record"])
    for m in ["Team Total", "Total", "Spread", "Moneyline"]:
        r = len(rows) + 1
        rows.append(row_formulas(m, mkt_extra(m), r))
    rows.append([])

    # === BY PLATFORM & MARKET TYPE ===
    rows.append(["BY PLATFORM & MARKET TYPE"])
    rows.append(["Platform / Market", "Bets", "Wagered", "Profit", "ROI", "Win Rate", "Record"])
    for label, sheet_vals in platforms:
        for m in ["Team Total", "Total", "Spread", "Moneyline"]:
            r = len(rows) + 1
            combined = plat_extra(sheet_vals) + mkt_extra(m)
            rows.append(row_formulas(f"{label} — {m}", combined, r))
    rows.append([])

    # === DAILY P&L ===
    rows.append(["DAILY P&L"])
    rows.append(["Date", "Bets", "Wagered", "Profit", "ROI", "Win Rate", "Cumulative P&L"])
    daily_start = len(rows) + 1
    for d in dates:
        r = len(rows) + 1
        de = date_extra(d)
        bets_expr = f'SUMPRODUCT({de}({AK}="AK")*({S}!J$2:J${N}<>""))'
        rows.append([
            d.strftime("%-m/%-d/%y"),
            f'=IF({bets_expr}=0,"",{bets_expr})',
            f'=IF(B{r}="","",SUMPRODUCT({de}({AK}="AK")*({S}!J$2:J${N}<>"")*{S}!H$2:H${N}))',
            f'=IF(B{r}="","",{_profit_expr(de)})',
            f'=IF(B{r}="","",D{r}/C{r})',
            f'=IF(B{r}="","",{_wins_expr(de)}/B{r})',
            f'=IF(B{r}="","",SUM(D${daily_start}:D{r}))',
        ])
    rows.append([])

    # === WEEKLY P&L ===
    rows.append(["WEEKLY P&L"])
    rows.append(["Week", "Bets", "Wagered", "Profit", "ROI", "Win Rate", "Cumulative P&L"])
    weekly_start = len(rows) + 1
    for monday in weeks:
        r = len(rows) + 1
        sunday = monday + timedelta(days=6)
        we = week_extra(monday, sunday)
        rows.append([
            f'Week of {monday.strftime("%-m/%-d")}',
            bets_f(we), wagered_f(we), profit_f(we),
            f'=IF(C{r}=0,0,D{r}/C{r})',
            f'=IF(B{r}=0,0,{_wins_expr(we)}/B{r})',
            f'=SUM(D${weekly_start}:D{r})',
        ])

    # ── Write to sheet ──
    print(f"\nWriting {len(rows)} rows...")
    service.spreadsheets().values().update(
        spreadsheetId=SPREADSHEET_ID,
        range=f"'{TAB}'!A1:G{len(rows)}",
        valueInputOption="USER_ENTERED",
        body={"values": rows},
    ).execute()

    print("Formulas written!")

    # ── Formatting ──
    apply_formatting(service, rows)
    print("Done! Open the sheet to see 'CBB Summary' tab.")


# ── Formatting helpers ────────────────────────────────────────────────────

def get_sheet_id(service, name):
    meta = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
    for s in meta.get("sheets", []):
        if s["properties"]["title"] == name:
            return s["properties"]["sheetId"]
    return None


def _text(sid, r1, r2, c1, c2, bold=False, size=None):
    tf = {"bold": bold}
    if size:
        tf["fontSize"] = size
    return {
        "repeatCell": {
            "range": {"sheetId": sid, "startRowIndex": r1, "endRowIndex": r2,
                      "startColumnIndex": c1, "endColumnIndex": c2},
            "cell": {"userEnteredFormat": {"textFormat": tf}},
            "fields": "userEnteredFormat.textFormat",
        }
    }


def _bg(sid, r1, r2, c1, c2, red, green, blue):
    return {
        "repeatCell": {
            "range": {"sheetId": sid, "startRowIndex": r1, "endRowIndex": r2,
                      "startColumnIndex": c1, "endColumnIndex": c2},
            "cell": {"userEnteredFormat": {
                "backgroundColor": {"red": red, "green": green, "blue": blue}
            }},
            "fields": "userEnteredFormat.backgroundColor",
        }
    }


def _num(sid, r1, r2, c1, c2, ntype, pattern):
    return {
        "repeatCell": {
            "range": {"sheetId": sid, "startRowIndex": r1, "endRowIndex": r2,
                      "startColumnIndex": c1, "endColumnIndex": c2},
            "cell": {"userEnteredFormat": {
                "numberFormat": {"type": ntype, "pattern": pattern}
            }},
            "fields": "userEnteredFormat.numberFormat",
        }
    }


def _cond(sid, r1, r2, c1, c2, cond_type, value, color):
    return {
        "addConditionalFormatRule": {
            "rule": {
                "ranges": [{"sheetId": sid, "startRowIndex": r1, "endRowIndex": r2,
                            "startColumnIndex": c1, "endColumnIndex": c2}],
                "booleanRule": {
                    "condition": {"type": cond_type,
                                  "values": [{"userEnteredValue": value}]},
                    "format": {"textFormat": {
                        "foregroundColorStyle": {"rgbColor": color}
                    }},
                },
            },
            "index": 0,
        }
    }


def apply_formatting(service, rows):
    sid = get_sheet_id(service, TAB)
    if sid is None:
        return

    reqs = []

    # Identify special rows
    sections = {
        "CBB ANSWER KEY PERFORMANCE", "OVERALL STATS", "BY PLATFORM",
        "BY MARKET TYPE", "DAILY P&L", "WEEKLY P&L",
    }
    section_idxs = [i for i, r in enumerate(rows) if r and r[0] in sections]
    header_idxs = [i for i, r in enumerate(rows) if r and len(r) > 1 and r[1] in ("Value", "Bets")]

    # Data section ranges (rows between header and next blank)
    data_ranges = []
    for i, r in enumerate(rows):
        if r and len(r) > 1 and r[1] == "Bets":
            s = i + 1
            e = s
            while e < len(rows) and rows[e] and rows[e][0]:
                e += 1
            if e > s:
                data_ranges.append((s, e))

    # Title
    reqs.append(_text(sid, 0, 1, 0, 7, bold=True, size=16))

    # Section headers
    for r in section_idxs:
        reqs.append(_bg(sid, r, r + 1, 0, 7, 0.85, 0.92, 1.0))
        reqs.append(_text(sid, r, r + 1, 0, 7, bold=True, size=12))

    # Column headers
    for r in header_idxs:
        reqs.append(_bg(sid, r, r + 1, 0, 7, 0.93, 0.93, 0.93))
        reqs.append(_text(sid, r, r + 1, 0, 7, bold=True))

    # Overall section formats (rows 8-14 = indices 7-13, column B = index 1)
    reqs.append(_num(sid, 7, 8, 1, 2, "CURRENCY", "$#,##0.00"))     # Wagered
    reqs.append(_num(sid, 8, 9, 1, 2, "CURRENCY", "$#,##0.00"))     # Profit
    reqs.append(_num(sid, 9, 10, 1, 2, "PERCENT", "0.00%"))         # ROI
    reqs.append(_num(sid, 10, 11, 1, 2, "PERCENT", "0.0%"))         # Win Rate
    reqs.append(_num(sid, 12, 13, 1, 2, "NUMBER", "+#,##0.00;-#,##0.00"))  # Units
    reqs.append(_num(sid, 13, 14, 1, 2, "NUMBER", "0.00"))          # Avg Odds

    # Data sections: C=Wagered($), D=Profit($), E=ROI(%), F=WinRate(%)
    for s, e in data_ranges:
        reqs.append(_num(sid, s, e, 2, 3, "CURRENCY", "$#,##0.00"))
        reqs.append(_num(sid, s, e, 3, 4, "CURRENCY", "$#,##0.00"))
        reqs.append(_num(sid, s, e, 4, 5, "PERCENT", "0.00%"))
        reqs.append(_num(sid, s, e, 5, 6, "PERCENT", "0.0%"))
        # Conditional formatting on Profit column (D)
        green = {"red": 0.13, "green": 0.55, "blue": 0.13}
        red = {"red": 0.80, "green": 0.13, "blue": 0.13}
        reqs.append(_cond(sid, s, e, 3, 4, "NUMBER_GREATER", "0", green))
        reqs.append(_cond(sid, s, e, 3, 4, "NUMBER_LESS", "0", red))

    # Cumulative P&L column (G) for daily/weekly sections
    for i, r in enumerate(rows):
        if r and len(r) > 6 and r[6] == "Cumulative P&L":
            s = i + 1
            e = s
            while e < len(rows) and rows[e] and rows[e][0]:
                e += 1
            if e > s:
                reqs.append(_num(sid, s, e, 6, 7, "CURRENCY", "$#,##0.00"))
                reqs.append(_cond(sid, s, e, 6, 7, "NUMBER_GREATER", "0", green))
                reqs.append(_cond(sid, s, e, 6, 7, "NUMBER_LESS", "0", red))

    # Column widths
    for i, w in enumerate([155, 80, 110, 115, 70, 80, 125]):
        reqs.append({
            "updateDimensionProperties": {
                "range": {"sheetId": sid, "dimension": "COLUMNS",
                          "startIndex": i, "endIndex": i + 1},
                "properties": {"pixelSize": w}, "fields": "pixelSize",
            }
        })

    # Hide helper sheet
    hlp_id = get_sheet_id(service, HLP)
    if hlp_id is not None:
        reqs.append({
            "updateSheetProperties": {
                "properties": {"sheetId": hlp_id, "hidden": True},
                "fields": "hidden",
            }
        })

    service.spreadsheets().batchUpdate(
        spreadsheetId=SPREADSHEET_ID,
        body={"requests": reqs},
    ).execute()
    print("Formatting applied!")


if __name__ == "__main__":
    main()
