"""Kalshi ticker grammar -> league / bet type / period / teams / side.

Port of the parsing half of extension/bets.js (the node tests on
tests/fixtures/bets/kalshi_fixture.json pin the semantics). Pure functions.

Facts (recon 2026-09-11, plan § Recon):
  event ticker   <SERIES>-<YYMMMDD>[HHMM]<AWAY><HOME>[G1|G2]; football suffixes
                 carry the date only (Eastern), MLB suffixes carry HHMM too.
  event.title    "Missouri St. vs Texas A&M: Spread" / "PIT Steelers vs NE
                 Patriots" — first team away, second home; the ": Market" tail
                 is absent on moneyline events.
  event.sub_title "MOSU vs TXAM (Sep 5)" — the codes in away/home order; a
                 spread or moneyline strike names its YES team by code.
  market         spread "-TXAM39" -> floor_strike 38.5 -> YES = TXAM -38.5,
                 NO = other team +38.5; total "-52" -> floor 51.5 -> YES = Over,
                 NO = Under; moneyline YES = the team, NO = the other team OR A
                 TIE in leagues that can end level (approx flag).
"""
import re

from unabated_ticket.bets_service.normalize import eastern_wall_clock_to_iso

TIE_CAVEAT = "kalshi_no_side_includes_tie"
# Leagues whose games can end level, so a NO on a team market also wins on a tie.
LEAGUES_WITH_TIES = frozenset({"nfl", "cfb", "soccer"})

# Series whose tickers this module reads (same table as extension/bets.js).
# `fixed_points` overrides floor_strike (RFI is "1st-inning runs >= 1" = Over 0.5).
GAME_SERIES: dict[str, dict] = {
    "KXNFLGAME": {"league": "nfl", "betType": "moneyline", "period": "FG"},
    "KXNFLSPREAD": {"league": "nfl", "betType": "spread", "period": "FG"},
    "KXNFLTOTAL": {"league": "nfl", "betType": "total", "period": "FG"},
    "KXNFL1HSPREAD": {"league": "nfl", "betType": "spread", "period": "1H"},
    "KXNFL1HTOTAL": {"league": "nfl", "betType": "total", "period": "1H"},
    "KXNCAAFGAME": {"league": "cfb", "betType": "moneyline", "period": "FG"},
    "KXNCAAFSPREAD": {"league": "cfb", "betType": "spread", "period": "FG"},
    "KXNCAAFTOTAL": {"league": "cfb", "betType": "total", "period": "FG"},
    "KXNCAAF1HSPREAD": {"league": "cfb", "betType": "spread", "period": "1H"},
    "KXNCAAF1HTOTAL": {"league": "cfb", "betType": "total", "period": "1H"},
    "KXMLBGAME": {"league": "mlb", "betType": "moneyline", "period": "FG"},
    "KXMLBSPREAD": {"league": "mlb", "betType": "spread", "period": "FG"},
    "KXMLBTOTAL": {"league": "mlb", "betType": "total", "period": "FG"},
    "KXMLBF5": {"league": "mlb", "betType": "moneyline", "period": "F5"},
    "KXMLBF5SPREAD": {"league": "mlb", "betType": "spread", "period": "F5"},
    "KXMLBF5TOTAL": {"league": "mlb", "betType": "total", "period": "F5"},
    "KXMLBRFI": {"league": "mlb", "betType": "total", "period": "I1", "fixed_points": 0.5},
}
# Futures, props and the bots' combos seen in the account: shown, never matched.
NON_GAME_SERIES = frozenset({
    "KXJOINCLUB", "KXWCAWARD", "KXNEXTTEAMNFL", "KXNFLOROTY", "KXSTARTINGQBWEEK1",
    "KXWCGOALLEADER", "KXMVECROSSCATEGORY",
})
REASON_NOT_GAME = "not a game market"
REASON_UNKNOWN_SERIES = "unknown Kalshi series"

_MONTHS = {"JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
           "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12}
_EVENT_SUFFIX_RE = re.compile(r"^(\d{2})([A-Z]{3})(\d{2})(\d{4})?([A-Z0-9]+?)(G[12])?$")
_EVENT_TITLE_RE = re.compile(r"^(.+?) vs (.+?)(?:: .+)?$")
_SUB_TITLE_CODES_RE = re.compile(r"^([A-Z0-9]+) vs ([A-Z0-9]+)")
_SPREAD_STRIKE_RE = re.compile(r"^([A-Z]+[A-Z0-9]*?)(\d+)$")


def series_of(ticker: str) -> str:
    return ticker.split("-")[0]


def parse_event_suffix(event_ticker: str) -> dict | None:
    """{eventDate, eventStart, gameNumber} from the event-ticker suffix, or None.
    eventDate is the Eastern date "YYYY-MM-DD"; eventStart is set only when the
    suffix carries HHMM (MLB)."""
    suffix = event_ticker[event_ticker.find("-") + 1:]
    match = _EVENT_SUFFIX_RE.match(suffix)
    if not match or match.group(2) not in _MONTHS:
        return None
    year = 2000 + int(match.group(1))
    month = _MONTHS[match.group(2)]
    day = int(match.group(3))
    event_start = None
    if match.group(4):
        hour = int(match.group(4)[:2])
        minute = int(match.group(4)[2:])
        event_start = eastern_wall_clock_to_iso(year, month, day, hour, minute)
    game_number = int(match.group(6)[1]) if match.group(6) else None
    return {"eventDate": f"{year:04d}-{month:02d}-{day:02d}", "eventStart": event_start,
            "gameNumber": game_number}


def parse_event_teams(event: dict) -> dict | None:
    """{awayTeam, homeTeam, codes: [awayCode, homeCode]} from the public event payload."""
    title = _EVENT_TITLE_RE.match(event.get("title") or "")
    codes = _SUB_TITLE_CODES_RE.match(event.get("sub_title") or "")
    if not title or not codes:
        return None
    return {"awayTeam": title.group(1), "homeTeam": title.group(2),
            "codes": [codes.group(1), codes.group(2)]}


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def parse_strike(market_ticker: str, market: dict, bet_type: str, fixed_points: float | None,
                 codes: list[str]) -> dict | str:
    """Which side (0 away / 1 home) a spread or moneyline strike names as YES and
    the strike's own number: {yesIndex, points}, or an error string."""
    strike = market_ticker[len(market["event_ticker"]) + 1:]
    if bet_type == "total":
        points = fixed_points if fixed_points is not None else market.get("floor_strike")
        if not _is_number(points):
            return f"no floor_strike on {market_ticker}"
        return {"yesIndex": None, "points": points}
    if bet_type == "moneyline":
        if strike not in codes:
            return f"strike {strike} not in {'/'.join(codes)}"
        return {"yesIndex": codes.index(strike), "points": None}
    match = _SPREAD_STRIKE_RE.match(strike)
    if not match or match.group(1) not in codes:
        return f"strike {strike} not in {'/'.join(codes)}"
    if not _is_number(market.get("floor_strike")):
        return f"no floor_strike on {market_ticker}"
    return {"yesIndex": codes.index(match.group(1)), "points": market["floor_strike"]}


def side_of_contract(bet_type: str, strike: dict, contract_side: str, league: str) -> dict:
    """Side + points for one (market, yes|no) in the record's convention: points
    is the side's own number. {side, points, approx}."""
    approx: list[str] = []
    if bet_type == "total":
        return {"side": "over" if contract_side == "yes" else "under",
                "points": strike["points"], "approx": approx}
    team_index = strike["yesIndex"] if contract_side == "yes" else 1 - strike["yesIndex"]
    side = "away" if team_index == 0 else "home"
    if bet_type == "moneyline":
        if contract_side == "no" and league in LEAGUES_WITH_TIES:
            approx.append(TIE_CAVEAT)
        return {"side": side, "points": None, "approx": approx}
    points = -strike["points"] if contract_side == "yes" else strike["points"]
    return {"side": side, "points": points, "approx": approx}
