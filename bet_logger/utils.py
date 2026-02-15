"""
Shared utility functions for bet scrapers.
Common odds calculations, sport detection, and parsing helpers.
"""

import re


def calculate_american_odds(risk: float, win: float) -> int:
    """Convert risk/win amounts to American odds format."""
    if risk <= 0 or win <= 0:
        return 0
    if win >= risk:
        return int(round((win / risk) * 100))
    else:
        return int(round(-(risk / win) * 100))


def calculate_american_odds_from_line(line_text: str, default: int = 0) -> int:
    """Extract American odds from line text like '+120' or '-110'."""
    match = re.search(r'([+-]\d+)', line_text)
    if match:
        return int(match.group(1))
    return default


def calculate_decimal_odds_from_american(american_odds: int) -> float:
    """Convert American odds to decimal odds format."""
    if american_odds > 0:
        return (american_odds / 100) + 1
    elif american_odds < 0:
        return (100 / abs(american_odds)) + 1
    return 1.0


def parse_status(status_text: str) -> str:
    """
    Extract result from status text.
    Returns: 'W', 'L', 'X', or ''
    """
    status_upper = status_text.upper().strip()
    if 'WIN' in status_upper or 'WON' in status_upper:
        return 'W'
    elif 'LOST' in status_upper or 'LOSE' in status_upper or 'LOSS' in status_upper:
        return 'L'
    elif 'PUSH' in status_upper or 'PUSHED' in status_upper or 'TIE' in status_upper or 'CANCELLED' in status_upper:
        return 'X'
    return ''


def parse_risk_win(risk_win_text: str) -> tuple:
    """
    Parse risk/win text into (risk, win) floats.
    Handles formats: "$ 125.00 / $ 325.00", "200.00/360.00", "125/325"
    """
    try:
        parts = risk_win_text.replace('$', '').replace(',', '').split('/')
        if len(parts) == 2:
            risk = float(parts[0].strip())
            win = float(parts[1].strip())
            return risk, win
    except (ValueError, IndexError):
        pass
    return 0.0, 0.0


# Complete team lists for sport detection
NFL_TEAMS = [
    'BEARS', 'LIONS', 'PACKERS', 'VIKINGS',
    'COWBOYS', 'EAGLES', 'GIANTS', 'COMMANDERS',
    'FALCONS', 'PANTHERS', 'SAINTS', 'BUCCANEERS',
    'CARDINALS', 'RAMS', 'SEAHAWKS', '49ERS',
    'BENGALS', 'BROWNS', 'RAVENS', 'STEELERS',
    'BILLS', 'DOLPHINS', 'PATRIOTS', 'JETS',
    'COLTS', 'JAGUARS', 'TEXANS', 'TITANS',
    'BRONCOS', 'CHARGERS', 'CHIEFS', 'RAIDERS',
]

NBA_TEAMS = [
    'CELTICS', 'NETS', 'KNICKS', '76ERS', 'RAPTORS',
    'BULLS', 'CAVALIERS', 'PISTONS', 'PACERS', 'BUCKS',
    'HAWKS', 'HORNETS', 'HEAT', 'MAGIC', 'WIZARDS',
    'NUGGETS', 'TIMBERWOLVES', 'THUNDER', 'BLAZERS', 'JAZZ',
    'WARRIORS', 'CLIPPERS', 'LAKERS', 'SUNS', 'KINGS',
    'MAVERICKS', 'ROCKETS', 'GRIZZLIES', 'PELICANS', 'SPURS',
]

NHL_TEAMS = [
    'BRUINS', 'SABRES', 'RED WINGS', 'CANADIENS',
    'SENATORS', 'LIGHTNING', 'MAPLE LEAFS', 'HURRICANES', 'BLUE JACKETS',
    'DEVILS', 'ISLANDERS', 'RANGERS', 'FLYERS', 'PENGUINS',
    'CAPITALS', 'BLACKHAWKS', 'AVALANCHE', 'STARS', 'WILD',
    'PREDATORS', 'BLUES', 'DUCKS', 'COYOTES',
    'FLAMES', 'OILERS', 'SHARKS', 'KRAKEN', 'GOLDEN KNIGHTS',
]


def parse_sport(*text_args: str) -> str:
    """
    Determine sport from one or more text fields.
    Accepts variable number of strings (description, note, selection, etc.)
    and combines them for matching.
    """
    combined = ' '.join(text_args).upper()

    # Check team names first (most reliable)
    for team in NFL_TEAMS:
        if team in combined:
            return 'NFL'

    for team in NBA_TEAMS:
        if team in combined:
            return 'NBA'

    for team in NHL_TEAMS:
        if team in combined:
            return 'NHL'

    # Keyword fallbacks (check college before generic 'FOOTBALL'/'BASKETBALL')
    if 'NCAAF' in combined or 'COLLEGE FOOTBALL' in combined:
        return 'NCAAF'
    if 'NCAAB' in combined or 'NCAAM' in combined or 'COLLEGE BASKETBALL' in combined:
        return 'NCAAM'
    if 'COLLEGE' in combined or 'BOWL' in combined:
        return 'NCAAF'
    if 'NFL' in combined or 'FOOTBALL' in combined:
        return 'NFL'
    if 'NBA' in combined or 'BASKETBALL' in combined:
        return 'NBA'
    if 'NHL' in combined or 'HOCKEY' in combined:
        return 'NHL'
    if 'MLB' in combined or 'BASEBALL' in combined:
        return 'MLB'

    return ''
