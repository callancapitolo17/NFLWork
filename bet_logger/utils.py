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


# ── Sport detection ──────────────────────────────────────────────
#
# Strategy (in order of priority):
#   1. Sport prefix — some APIs prepend the sport category to descriptions
#      (e.g. BetOnline: "Football - 123 Team...", "Soccer - MLS - ...")
#   2. League/keyword patterns — matches league names, acronyms, and
#      sport-specific terms (e.g. "Premier League", "UFC", "ATP")
#   3. Team name fallback — hardcoded team lists for the 4 major US leagues
#      (NFL, NBA, NHL, MLB). Used when the API doesn't provide a sport field.

# Maps a description prefix to the sport label returned.
# BetOnline descriptions start with these; other books might too.
# The prefix is stripped after detection, so "Football" here means the
# API literally starts the description with "Football - ...".
SPORT_PREFIXES = {
    'FOOTBALL': '_football',   # Could be NFL or NCAAF — needs further check
    'BASKETBALL': '_basketball',  # Could be NBA or NCAAM — needs further check
    'HOCKEY': 'NHL',
    'BASEBALL': '_baseball',   # Could be MLB or college — needs further check
    'SOCCER': 'Soccer',
    'TENNIS': 'Tennis',
    'FIGHTING': 'MMA',
    'MARTIAL ARTS': 'MMA',
    'GOLF': 'Golf',
    'MOTOR SPORTS': 'NASCAR',
    'MOTORSPORTS': 'NASCAR',
    'BOXING': 'Boxing',
    'RUGBY': 'Rugby',
    'CRICKET': 'Cricket',
    'TABLE TENNIS': 'Table Tennis',
    'ESPORTS': 'Esports',
    'ENTERTAINMENT': 'Entertainment',
    'POLITICS': 'Politics',
}

# League name patterns → sport label.
# Checked via substring match against the combined text.
LEAGUE_KEYWORDS = [
    # Soccer leagues
    ('PREMIER LEAGUE', 'Soccer'), ('LA LIGA', 'Soccer'), ('SERIE A', 'Soccer'),
    ('BUNDESLIGA', 'Soccer'), ('LIGUE 1', 'Soccer'), ('CHAMPIONS LEAGUE', 'Soccer'),
    ('EUROPA LEAGUE', 'Soccer'), ('MLS', 'Soccer'), ('LIGA MX', 'Soccer'),
    ('WORLD CUP', 'Soccer'), ('CONCACAF', 'Soccer'), ('EPL', 'Soccer'),
    # Tennis
    ('ATP', 'Tennis'), ('WTA', 'Tennis'), ('ROLAND GARROS', 'Tennis'),
    ('WIMBLEDON', 'Tennis'), ('US OPEN TENNIS', 'Tennis'),
    # MMA/Boxing
    ('UFC', 'MMA'), ('BELLATOR', 'MMA'), ('PFL', 'MMA'),
    ('BOXING', 'Boxing'),
    # Golf
    ('PGA', 'Golf'), ('LPGA', 'Golf'), ('LIV GOLF', 'Golf'),
    ('MASTERS', 'Golf'), ('RYDER CUP', 'Golf'),
    # Racing
    ('NASCAR', 'NASCAR'), ('FORMULA 1', 'F1'), ('F1', 'F1'),
    ('INDYCAR', 'IndyCar'),
    # Other US leagues
    ('WNBA', 'WNBA'), ('CFL', 'CFL'), ('XFL', 'XFL'), ('UFL', 'UFL'),
    ('USFL', 'USFL'), ('NLL', 'NLL'), ('MLL', 'MLL'),
]

# Team name lists — fallback when no prefix or keyword matches.
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
    'PREDATORS', 'BLUES', 'DUCKS', 'COYOTES', 'UTAH HOCKEY CLUB',
    'FLAMES', 'OILERS', 'SHARKS', 'KRAKEN', 'GOLDEN KNIGHTS',
]

MLB_TEAMS = [
    'YANKEES', 'RED SOX', 'BLUE JAYS', 'RAYS', 'ORIOLES',
    'WHITE SOX', 'GUARDIANS', 'TWINS', 'ROYALS', 'TIGERS',
    'ASTROS', 'MARINERS', 'ANGELS', 'ATHLETICS', 'RANGERS',
    'METS', 'BRAVES', 'PHILLIES', 'MARLINS', 'NATIONALS',
    'CUBS', 'BREWERS', 'REDS', 'PIRATES', 'CARDINALS',
    'DODGERS', 'PADRES', 'DIAMONDBACKS', 'ROCKIES', 'GIANTS',
]

# Teams that appear in multiple leagues need context to disambiguate.
# "RANGERS" is in both NHL and MLB; "GIANTS" in NFL and MLB; "CARDINALS" in NFL and MLB.
# We check in league order: NFL > NBA > NHL > MLB, so NFL wins for shared names.
# This is imperfect but matches the existing behavior. The prefix-based detection
# (layer 1) resolves this correctly when available.


def _detect_from_prefix(text_upper: str) -> str:
    """Layer 1: Check if text starts with a known sport prefix (e.g. 'Football - ...')."""
    for prefix, sport in SPORT_PREFIXES.items():
        if text_upper.startswith(prefix):
            return sport
    return ''


def _resolve_ambiguous(sport_hint: str, text_upper: str) -> str:
    """Resolve ambiguous sport hints like '_football' into NFL vs NCAAF.

    When the API tells us the sport category (e.g. "Football") but not the
    specific league, we check for pro team names to distinguish NFL vs NCAAF.
    """
    if sport_hint == '_football':
        if 'NFL' in text_upper:
            return 'NFL'
        for team in NFL_TEAMS:
            if team in text_upper:
                return 'NFL'
        return 'NCAAF'

    if sport_hint == '_basketball':
        if 'NBA' in text_upper:
            return 'NBA'
        for team in NBA_TEAMS:
            if team in text_upper:
                return 'NBA'
        return 'NCAAM'

    if sport_hint == '_baseball':
        if 'MLB' in text_upper:
            return 'MLB'
        for team in MLB_TEAMS:
            if team in text_upper:
                return 'MLB'
        return 'Baseball'

    return sport_hint


def parse_sport(*text_args: str) -> str:
    """
    Determine sport from one or more text fields.

    Uses a 3-layer approach:
      1. Sport prefix (e.g. BetOnline's "Football - ...")
      2. League/keyword patterns (e.g. "Premier League", "UFC", "ATP")
      3. Team name fallback (NFL, NBA, NHL, MLB team lists)

    Accepts variable number of strings (description, note, selection, etc.)
    and combines them for matching.
    """
    combined = ' '.join(text_args).upper()
    if not combined.strip():
        return ''

    # Strip device prefix (BetOnline prepends "Desktop - " or "Mobile - ")
    combined = re.sub(r'^(DESKTOP|MOBILE)\s*-\s*', '', combined)

    # Layer 1: Sport prefix — most reliable when available.
    # This catches BetOnline-style descriptions that start with the sport.
    sport = _detect_from_prefix(combined)
    if sport:
        return _resolve_ambiguous(sport, combined)

    # Layer 2: League/keyword patterns — catches named leagues and acronyms.
    # Check college-specific keywords before generic sport keywords to avoid
    # mis-labeling "College Basketball" as "NBA" just because "BASKETBALL" appears.
    if 'NCAAF' in combined or 'COLLEGE FOOTBALL' in combined:
        return 'NCAAF'
    if 'NCAAB' in combined or 'NCAAM' in combined or 'COLLEGE BASKETBALL' in combined:
        return 'NCAAM'

    for keyword, sport in LEAGUE_KEYWORDS:
        if keyword in combined:
            return sport

    # Generic sport keywords (after college check to avoid misclassification)
    if 'NFL' in combined:
        return 'NFL'
    if 'NBA' in combined:
        return 'NBA'
    if 'NHL' in combined or 'HOCKEY' in combined:
        return 'NHL'
    if 'MLB' in combined or 'BASEBALL' in combined:
        return 'MLB'
    if 'FOOTBALL' in combined:
        # Bare "FOOTBALL" without NFL/NCAAF — check for pro teams
        for team in NFL_TEAMS:
            if team in combined:
                return 'NFL'
        return 'NCAAF'
    if 'BASKETBALL' in combined:
        for team in NBA_TEAMS:
            if team in combined:
                return 'NBA'
        return 'NCAAM'
    if 'SOCCER' in combined:
        return 'Soccer'
    if 'TENNIS' in combined:
        return 'Tennis'

    # Layer 3: Team name fallback — scan for known team names.
    for team in NFL_TEAMS:
        if team in combined:
            return 'NFL'
    for team in NBA_TEAMS:
        if team in combined:
            return 'NBA'
    for team in NHL_TEAMS:
        if team in combined:
            return 'NHL'
    for team in MLB_TEAMS:
        if team in combined:
            return 'MLB'

    return ''
