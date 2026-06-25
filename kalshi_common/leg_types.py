"""Kalshi MLB leg-typing helpers shared between taker and maker bots.

Converts raw Kalshi market-ticker dicts to typed fair_value.SpreadLeg /
fair_value.TotalLeg instances, and extracts canonical spread / total line
values from a legs list.
"""
from dataclasses import dataclass

from kalshi_common import fair_value

# The 4-cell devig families stored in mlb_sgp_odds. Each family partitions the
# outcome space, so devig_book n-way-devigs the four cells together.
SPREAD_TOTAL_FAMILY = (
    "Home Spread + Over", "Home Spread + Under",
    "Away Spread + Over", "Away Spread + Under",
)
ML_TOTAL_FAMILY = (
    "Home ML + Over", "Home ML + Under",
    "Away ML + Over", "Away ML + Under",
)

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


def _moneyline_side(leg: dict) -> tuple[bool, str] | None:
    """(team_is_home, side) for a KXMLBGAME (moneyline) leg, or None.

    KXMLBGAME-{event_suffix}-{TEAMCODE}: YES = that team wins. team_is_home iff
    the ticker's team code equals the home code parsed from the event ticker.
    """
    mt = str(leg.get("market_ticker", ""))
    team = mt.rsplit("-", 1)[-1]
    home = _home_code_from_event_ticker(str(leg.get("event_ticker", "")))
    if not team or home is None:
        return None
    return (team == home, leg.get("side", "yes"))


def _event_codes_from_legs(legs: list[dict]) -> tuple[str | None, str | None]:
    """(away_code, home_code) parsed from any leg's event ticker. All legs of a
    game share the same event suffix, so the first resolvable one wins."""
    for l in legs:
        et = str(l.get("event_ticker", ""))
        if "-" in et:
            away, home = _parse_event_suffix(et.rsplit("-", 1)[-1])
            if away and home:
                return away, home
    return None, None


@dataclass(frozen=True)
class ComboDescriptor:
    """How to look a 2-leg Kalshi combo up in the mlb_sgp_odds grid.

    kind         : "spread_total" | "ml_total"
    spread_line  : home-perspective spread (None for ml_total — there's no spread)
    total_line   : the Over/Under line
    target_combo : the single grid cell matching the legs' STATED sides
                   (e.g. "Away ML + Under")
    combo_family : the 4 cells to n-way devig together
    away_code    : Kalshi away team code (for game resolution)
    home_code    : Kalshi home team code
    """
    kind: str
    spread_line: float | None
    total_line: float | None
    target_combo: str
    combo_family: tuple
    away_code: str
    home_code: str


def combo_descriptor(legs: list[dict]) -> "ComboDescriptor | None":
    """Classify a 2-leg combo into a grid-lookup descriptor, or None if it isn't
    a shape the SGP grid carries ({spread,total} or {moneyline,total}).

    The target cell is derived from the ACTUAL leg sides (not hardcoded), which
    is also the fix for the old maker's "Home Spread + Over"-always bug.
    """
    if not legs or len(legs) != 2:
        return None

    def _pfx(l):
        return str(l.get("market_ticker", "")).split("-")[0]

    by: dict[str, list[dict]] = {}
    for l in legs:
        by.setdefault(_pfx(l), []).append(l)

    away, home = _event_codes_from_legs(legs)
    if not away or not home:
        return None

    total_legs = by.get("KXMLBTOTAL")
    if not total_legs or len(total_legs) != 1:
        return None
    total = _leg_dict_to_typed(total_legs[0], "")        # TotalLeg
    if total is None:
        return None
    total_line = total.line_n - 0.5
    over_part = "Over" if total.side == "yes" else "Under"

    if len(by.get("KXMLBSPREAD", [])) == 1:
        spread = _leg_dict_to_typed(by["KXMLBSPREAD"][0], "")   # SpreadLeg
        if spread is None:
            return None
        # At a half-point line exactly one team covers, so spread 'no' == other.
        home_covers = ((spread.team_is_home and spread.side == "yes")
                       or (not spread.team_is_home and spread.side == "no"))
        part = "Home" if home_covers else "Away"
        return ComboDescriptor(
            kind="spread_total",
            spread_line=-(spread.line_n - 0.5),
            total_line=total_line,
            target_combo=f"{part} Spread + {over_part}",
            combo_family=SPREAD_TOTAL_FAMILY,
            away_code=away, home_code=home)

    if len(by.get("KXMLBGAME", [])) == 1:
        ml = _moneyline_side(by["KXMLBGAME"][0])
        if ml is None:
            return None
        team_is_home, side = ml
        home_ml = ((team_is_home and side == "yes")
                   or (not team_is_home and side == "no"))
        part = "Home" if home_ml else "Away"
        return ComboDescriptor(
            kind="ml_total",
            spread_line=None,
            total_line=total_line,
            target_combo=f"{part} ML + {over_part}",
            combo_family=ML_TOTAL_FAMILY,
            away_code=away, home_code=home)

    return None
