"""
Wagerzon Scraper Configuration
Sport-agnostic configuration for scraping different leagues via JSON API.
"""

# Sport configurations
SPORTS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "url_params": "lg=4029,430,432,433,2239,2140,98,2767,124",
        "table_name": "nfl_odds",
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "url_params": "lg=43,403,45",
        "table_name": "cbb_odds",
    },
    "nba": {
        "sport_key": "basketball_nba",
        "url_params": "lg=535,443,626,2917,761",
        "table_name": "nba_odds",
    },
    "college_baseball": {
        "sport_key": "baseball_ncaa",
        "url_params": "lg=762,1554,4321",  # 762=College Baseball, 1554=NCAA Game, 4321=1st 5 Innings
        "table_name": "college_baseball_odds",
    },
    "mlb": {
        "sport_key": "baseball_mlb",
        "url_params": "lg=121,2904,4320",  # 121=MLB, 2904=MLB Game, 4320=1st 5 Innings
        "table_name": "mlb_odds",
    },
}

# Wagerzon backend URLs
WAGERZON_BASE_URL = "https://backend.wagerzon.com"
WAGERZON_HELPER_URL = f"{WAGERZON_BASE_URL}/wager/NewScheduleHelper.aspx"
WAGERZON_SCHEDULE_URL = f"{WAGERZON_BASE_URL}/wager/NewSchedule.aspx?WT=0&"


def get_sport_url(sport: str) -> str:
    """Get the full schedule page URL for a sport."""
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    return f"{WAGERZON_SCHEDULE_URL}{SPORTS[sport]['url_params']}"


def get_sport_config(sport: str) -> dict:
    """Get configuration for a sport."""
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    return SPORTS[sport]
