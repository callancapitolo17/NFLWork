"""Kalshi MLB leg-typing helpers shared between taker and maker bots.

Converts raw Kalshi market-ticker dicts to typed fair_value.SpreadLeg /
fair_value.TotalLeg instances, and extracts canonical spread / total line
values from a legs list.
"""
from kalshi_common import fair_value

# 3-letter Kalshi team code → mlb_parlay_lines.home_team / away_team canonical name.
# Kalshi uses 3-letter codes; mlb_parlay_lines stores Odds-API canonical names.
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "ATH": "Athletics",
    "AZ": "Arizona Diamondbacks", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


def _parse_event_suffix(suffix: str) -> tuple[str | None, str | None]:
    """Split a KXMLB* event suffix into (away_code, home_code).

    Format: YYMMMDDHHMM{AwayCode}{HomeCode}. Date prefix is fixed at 11 chars.
    Each team code is 2 or 3 letters (KC/SF/SD/TB/AZ are 2-letter; the rest
    are 3-letter). Probes 3- then 2-letter home splits and returns the first
    where both codes are valid in _MLB_CODE_TO_TEAM. Returns (None, None) if
    no split matches — caller drops the event.
    """
    if len(suffix) < 11 + 4:  # date prefix + at least 2+2 team chars
        return None, None
    team_block = suffix[11:]
    for home_len in (3, 2):
        if len(team_block) <= home_len:
            continue
        home = team_block[-home_len:]
        away = team_block[:-home_len]
        if home in _MLB_CODE_TO_TEAM and away in _MLB_CODE_TO_TEAM:
            return away, home
    return None, None


def _home_code_from_event_ticker(event_ticker: str) -> str | None:
    """Parse the home-team code from a Kalshi event ticker (2- or 3-letter)."""
    if "-" not in event_ticker:
        return None
    suffix = event_ticker.rsplit("-", 1)[-1]
    _, home = _parse_event_suffix(suffix)
    return home


def _leg_dict_to_typed(leg: dict, game_id: str):
    """Convert {market_ticker, event_ticker, side} to fair_value typed leg.

    Determines team_is_home by parsing the home code from the event_ticker
    (no DB lookup needed — the ticker self-encodes the home/away convention).
    """
    mt = leg["market_ticker"]
    et = leg.get("event_ticker", "")
    side = leg["side"]
    if mt.startswith("KXMLBSPREAD-"):
        suffix = mt.rsplit("-", 1)[-1]
        n_chars = "".join(c for c in suffix if c.isdigit())
        team_chars = "".join(c for c in suffix if not c.isdigit())
        if not n_chars or not team_chars:
            return None
        n = int(n_chars)
        home_code = _home_code_from_event_ticker(et)
        team_is_home = (home_code is not None and team_chars == home_code)
        return fair_value.SpreadLeg(team_is_home=team_is_home, line_n=n, side=side)
    if mt.startswith("KXMLBTOTAL-"):
        try:
            n = int(mt.rsplit("-", 1)[-1])
        except ValueError:
            return None
        return fair_value.TotalLeg(line_n=n, side=side)
    return None


def _spread_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBSPREAD-"):
            suffix = l["market_ticker"].rsplit("-", 1)[-1]
            digits = "".join(c for c in suffix if c.isdigit())
            if digits:
                n = int(digits)
                return -(n - 0.5)
    return 0.0


def _total_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBTOTAL-"):
            try:
                n = int(l["market_ticker"].rsplit("-", 1)[-1])
                return n - 0.5
            except ValueError:
                continue
    return 0.0
