"""MLB run-totals adapter: Unabated bt3 anchor -> Kalshi MLB total ladder.

Constants below (league prefix, series) were live-verified in Task 1 recon
(plan Appendix A) - update there if they drift.

Task 1 found Kalshi's MLB event titles are CITY-based with short
disambiguator suffixes ("New York Y vs Philadelphia: Total Runs"), not
nickname-based like soccer's World Cup titles ("Colombia vs Ghana: ...").
Free-text title parsing is fragile against that format (e.g. the lone
suffix letter distinguishing "Chicago C" (Cubs) from "Chicago WS" (White
Sox), or "A's" replacing a city entirely for the Athletics). event_teams
instead reads the two club codes out of the event's `event_ticker`, e.g.
"KXMLBTOTAL-26JUL251805NYYPHI" (NYY=Yankees, PHI=Phillies), through an
explicit code table, and fails closed (empty frozenset) on anything that
doesn't parse into two known codes."""
import re

from unabated_edge.sports.totals import TotalsLadderAdapter

MLB_TOTAL_SERIES = "KXMLBTOTAL"   # Task 1 Appendix A

# Unabated's resolved team-name spelling is "{City} {Nickname}" (Task 1
# Appendix A sample events, e.g. "Pittsburgh Pirates", "Chicago Cubs") ->
# canonical nickname, so the Unabated side collapses onto the same key the
# ticker-derived Kalshi side produces.
_MLB_TEAMS = {
    "arizona diamondbacks": "Diamondbacks", "atlanta braves": "Braves",
    "baltimore orioles": "Orioles", "boston red sox": "Red Sox",
    "chicago cubs": "Cubs", "chicago white sox": "White Sox",
    "cincinnati reds": "Reds", "cleveland guardians": "Guardians",
    "colorado rockies": "Rockies", "detroit tigers": "Tigers",
    "houston astros": "Astros", "kansas city royals": "Royals",
    "los angeles angels": "Angels", "los angeles dodgers": "Dodgers",
    "miami marlins": "Marlins", "milwaukee brewers": "Brewers",
    "minnesota twins": "Twins", "new york mets": "Mets",
    "new york yankees": "Yankees", "athletics": "Athletics",
    "oakland athletics": "Athletics", "philadelphia phillies": "Phillies",
    "pittsburgh pirates": "Pirates", "san diego padres": "Padres",
    "san francisco giants": "Giants", "seattle mariners": "Mariners",
    "st. louis cardinals": "Cardinals", "st louis cardinals": "Cardinals",
    "tampa bay rays": "Rays", "texas rangers": "Rangers",
    "toronto blue jays": "Blue Jays", "washington nationals": "Nationals",
}
# Nicknames map to themselves so a bare Kalshi-derived nickname also passes
# through canon unchanged.
_MLB_TEAMS.update({v.lower(): v for v in list(_MLB_TEAMS.values())})

# Kalshi event_ticker club codes -> canonical nickname. Codes actually
# observed live in Task 1 recon (NYY, PHI, CHC, PIT, ATH, MIN, LAD, NYM,
# LAA, SF) plus AZ (observed in Task 4's own live pairing smoke against
# today's slate — Arizona uses "AZ", not the "ARI" this file first guessed;
# kept both since they collapse to the same nickname). The remaining codes
# follow the same standard 2/3-letter convention (including the city-less
# "ATH" for the Athletics) but are not yet individually verified against a
# live ticker. A few other common alternate spellings (CWS/CHW, WSH/WAS,
# KC/KAN) are included defensively — harmless duplicates since they map to
# the same nickname, not a correctness risk. If a future live pairing smoke
# (Task 4 brief Step 5) shows a low match rate, diff against a fuller slate
# and correct entries here.
_MLB_CLUB_CODES = {
    "ARI": "Diamondbacks", "AZ": "Diamondbacks", "ATL": "Braves",
    "BAL": "Orioles", "BOS": "Red Sox",
    "CHC": "Cubs", "CWS": "White Sox", "CHW": "White Sox", "CIN": "Reds",
    "CLE": "Guardians", "COL": "Rockies", "DET": "Tigers", "HOU": "Astros",
    "KC": "Royals", "KAN": "Royals", "LAA": "Angels", "LAD": "Dodgers",
    "MIA": "Marlins", "MIL": "Brewers", "MIN": "Twins", "NYM": "Mets",
    "NYY": "Yankees", "ATH": "Athletics", "PHI": "Phillies", "PIT": "Pirates",
    "SD": "Padres", "SF": "Giants", "SEA": "Mariners", "STL": "Cardinals",
    "TB": "Rays", "TEX": "Rangers", "TOR": "Blue Jays",
    "WSH": "Nationals", "WAS": "Nationals",
}

# event_ticker shape verified live (Task 1): "KXMLBTOTAL-26JUL251805NYYPHI"
# = series + "-" + an 11-char date/time block (2-digit year + 3-letter
# month + 2-digit day + 4-digit HHMM) + the two club codes concatenated
# with no separator between them (each code is 2 or 3 letters).
_TICKER_RE = re.compile(r"^KXMLBTOTAL-\d{2}[A-Z]{3}\d{6}([A-Z]+)$")


class Mlb(TotalsLadderAdapter):
    sport = "mlb"
    league_prefix = "lg5"   # Task 1 Appendix A

    def canon_team(self, name: str) -> str:
        key = (name or "").strip().lower()
        return _MLB_TEAMS.get(key, (name or "").strip())

    def kalshi_series(self) -> str:
        return MLB_TOTAL_SERIES

    def event_teams(self, kalshi_event: dict) -> frozenset:
        """Parse the two club codes out of event_ticker (see module
        docstring for why title parsing is unsafe for MLB). Fails closed
        (empty frozenset) when the ticker doesn't match the expected shape
        or its trailing letters don't split into two known club codes."""
        ticker = (kalshi_event.get("event_ticker") or "")
        m = _TICKER_RE.match(ticker)
        if not m:
            return frozenset()
        codes = m.group(1)
        for away_len in (2, 3):
            home_len = len(codes) - away_len
            if home_len not in (2, 3):
                continue
            away, home = codes[:away_len], codes[away_len:]
            if away in _MLB_CLUB_CODES and home in _MLB_CLUB_CODES:
                return frozenset({_MLB_CLUB_CODES[away], _MLB_CLUB_CODES[home]})
        return frozenset()
