"""Team-name resolution: Wagerzon nicknames/abbreviations and Odds API full
names -> nflverse team abbreviations (the codes in nfl_game_outcomes)."""
from __future__ import annotations

# nflverse abbreviation -> (Odds API full name, Wagerzon nickname)
TEAMS: dict[str, tuple[str, str]] = {
    "ARI": ("Arizona Cardinals", "CARDINALS"),
    "ATL": ("Atlanta Falcons", "FALCONS"),
    "BAL": ("Baltimore Ravens", "RAVENS"),
    "BUF": ("Buffalo Bills", "BILLS"),
    "CAR": ("Carolina Panthers", "PANTHERS"),
    "CHI": ("Chicago Bears", "BEARS"),
    "CIN": ("Cincinnati Bengals", "BENGALS"),
    "CLE": ("Cleveland Browns", "BROWNS"),
    "DAL": ("Dallas Cowboys", "COWBOYS"),
    "DEN": ("Denver Broncos", "BRONCOS"),
    "DET": ("Detroit Lions", "LIONS"),
    "GB": ("Green Bay Packers", "PACKERS"),
    "HOU": ("Houston Texans", "TEXANS"),
    "IND": ("Indianapolis Colts", "COLTS"),
    "JAX": ("Jacksonville Jaguars", "JAGUARS"),
    "KC": ("Kansas City Chiefs", "CHIEFS"),
    "LV": ("Las Vegas Raiders", "RAIDERS"),
    "LAC": ("Los Angeles Chargers", "CHARGERS"),
    "LA": ("Los Angeles Rams", "RAMS"),
    "MIA": ("Miami Dolphins", "DOLPHINS"),
    "MIN": ("Minnesota Vikings", "VIKINGS"),
    "NE": ("New England Patriots", "PATRIOTS"),
    "NO": ("New Orleans Saints", "SAINTS"),
    "NYG": ("New York Giants", "GIANTS"),
    "NYJ": ("New York Jets", "JETS"),
    "PHI": ("Philadelphia Eagles", "EAGLES"),
    "PIT": ("Pittsburgh Steelers", "STEELERS"),
    "SF": ("San Francisco 49ers", "49ERS"),
    "SEA": ("Seattle Seahawks", "SEAHAWKS"),
    "TB": ("Tampa Bay Buccaneers", "BUCCANEERS"),
    "TEN": ("Tennessee Titans", "TITANS"),
    "WAS": ("Washington Commanders", "COMMANDERS"),
}

# Wagerzon's own abbreviations where they differ from nflverse's.
WZ_ABBREVIATION_ALIASES = {
    "LVR": "LV", "NOS": "NO", "LAR": "LA", "WSH": "WAS", "JAC": "JAX",
    "GNB": "GB", "KAN": "KC", "NWE": "NE", "SFO": "SF", "TAM": "TB",
    "BUCS": "TB", "NINERS": "SF",
}

ODDS_API_NAME_TO_ABBR = {full: abbr for abbr, (full, _) in TEAMS.items()}
NICKNAME_TO_ABBR = {nick: abbr for abbr, (_, nick) in TEAMS.items()}


def resolve_wz_team(token: str) -> str | None:
    """'RAIDERS' / 'LVR' / 'LV' -> 'LV'; None when the token is not a team."""
    token = token.strip().upper()
    if token in NICKNAME_TO_ABBR:
        return NICKNAME_TO_ABBR[token]
    if token in WZ_ABBREVIATION_ALIASES:
        return WZ_ABBREVIATION_ALIASES[token]
    if token in TEAMS:
        return token
    return None
