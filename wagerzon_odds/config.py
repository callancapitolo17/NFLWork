"""
Wagerzon Scraper Configuration
Sport-agnostic configuration for scraping different leagues
"""

# Sport configurations
SPORTS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "url_params": "lg=4029,430,432,433,2239,2140,98,2767,124",
        "table_name": "nfl_odds",
        # Section header patterns -> market names
        "section_markets": {
            "SUPER BOWL": "spreads",           # Full game
            "1ST HALVES": "spreads_h1",
            "QUARTERS": "spreads_quarters",    # Will be split by rotation prefix
            "TEAM TOTALS": "team_totals",
            "1ST HALF (3WAY)": "h2h_h1_3way",
            "1ST QUARTER (3-WAY)": "h2h_q1_3way",
        },
        # Rotation number prefixes -> periods
        "rotation_periods": {
            "1": "h1",   # 1xxx = 1st half
            "2": "h2",   # 2xxx = 2nd half (if available)
            "3": "q1",   # 3xxx = 1st quarter
            "4": "q2",   # 4xxx = 2nd quarter
            "5": "q3",   # 5xxx = 3rd quarter
            "6": "q4",   # 6xxx = 4th quarter
        },
        # Team name cleanup patterns
        "team_prefix_patterns": ["1H ", "2H ", "1Q ", "2Q ", "3Q ", "4Q "],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "url_params": "lg=43,403,45",
        "table_name": "cbb_odds",
        "section_markets": {
            "COLLEGE BASKETBALL": "spreads",
            "1H COLLEGE BASKETBALL": "spreads_h1",
        },
        "rotation_periods": {
            "1": "h1",
        },
        "team_prefix_patterns": ["1H "],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "url_params": "lg=3,301,303",  # Placeholder - needs discovery
        "table_name": "nba_odds",
        "section_markets": {
            "NBA": "spreads",
            "1ST HALVES": "spreads_h1",
            "QUARTERS": "spreads_quarters",
        },
        "rotation_periods": {
            "1": "h1",
            "3": "q1",
            "4": "q2",
            "5": "q3",
            "6": "q4",
        },
        "team_prefix_patterns": ["1H ", "1Q ", "2Q ", "3Q ", "4Q "],
    },
}

# Markets to skip (props, futures, etc. that don't fit standard format)
SKIP_SECTIONS = [
    "COIN TOSS",
    "TD SCORER",
    "MVP",
    "WINNER",
    "FUTURES",
    "PROPS",
    "TEAM TOTAL",  # Different structure - handle separately later
    "3WAY",        # 3-way markets have different structure
    "3-WAY",
]

# Wagerzon backend URL
WAGERZON_BASE_URL = "https://backend.wagerzon.com"
WAGERZON_SCHEDULE_URL = f"{WAGERZON_BASE_URL}/wager/NewSchedule.aspx?WT=0&"


def get_sport_url(sport: str) -> str:
    """Get the full URL for a sport."""
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    return f"{WAGERZON_SCHEDULE_URL}{SPORTS[sport]['url_params']}"


def get_sport_config(sport: str) -> dict:
    """Get configuration for a sport."""
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    return SPORTS[sport]
