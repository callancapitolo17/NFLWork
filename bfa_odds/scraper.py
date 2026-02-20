#!/usr/bin/env python3
"""
BFA Gaming Odds Scraper for Answer Key Pipeline
Fetches derivative markets (1H, 1Q, 2Q, 3Q, 4Q) via REST API and saves to DuckDB.

No auth required — the BFA odds API is publicly accessible.
"""

import sys
import duckdb
import requests
from datetime import datetime, timezone
from pathlib import Path

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

DB_PATH = Path(__file__).parent / "bfa.duckdb"

# API configuration
API_BASE = "https://api.bfagaming.com"
API_PARAMS = {"playerId": 0, "agentId": 0, "fixtureType": 0}
API_HEADERS = {
    "Authorization": "Bearer",
    "User-Agent": "bfa-client/1.0.0",
    "Referer": "https://bfagaming.com/",
}

# Sport configurations
SPORT_CONFIGS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "api_slug": "nfl",
        "league_id": 889,
        "periods": ["1H", "1Q", "2Q", "3Q", "4Q"],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "api_slug": "nba",
        "league_id": 487,
        "periods": ["1H", "1Q", "2Q", "3Q", "4Q"],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "api_slug": "ncaa_b",
        "league_id": 493,
        "periods": ["1H"],
    },
}

# Map period short names to display names (must match API period "name" field)
PERIOD_DISPLAY = {
    "1H": "1st Half",
    "1Q": "1st Quarter",
    "2Q": "2nd Quarter",
    "3Q": "3rd Quarter",
    "4Q": "4th Quarter",
}

# Map to standard period names for R
PERIOD_TO_STANDARD = {
    "1H": "Half1",
    "1Q": "1",
    "2Q": "2",
    "3Q": "3",
    "4Q": "4",
}

# Map period short names to market suffix (matching Wagerzon/API format)
PERIOD_TO_MARKET_SUFFIX = {
    "1H": "h1",
    "1Q": "q1",
    "2Q": "q2",
    "3Q": "q3",
    "4Q": "q4",
}

# NFL team name mapping (short names -> Odds API full names)
NFL_TEAM_NAMES = {
    "49ers": "San Francisco 49ers",
    "Bears": "Chicago Bears",
    "Bengals": "Cincinnati Bengals",
    "Bills": "Buffalo Bills",
    "Broncos": "Denver Broncos",
    "Browns": "Cleveland Browns",
    "Buccaneers": "Tampa Bay Buccaneers",
    "Cardinals": "Arizona Cardinals",
    "Chargers": "Los Angeles Chargers",
    "Chiefs": "Kansas City Chiefs",
    "Colts": "Indianapolis Colts",
    "Commanders": "Washington Commanders",
    "Cowboys": "Dallas Cowboys",
    "Dolphins": "Miami Dolphins",
    "Eagles": "Philadelphia Eagles",
    "Falcons": "Atlanta Falcons",
    "Giants": "New York Giants",
    "Jaguars": "Jacksonville Jaguars",
    "Jets": "New York Jets",
    "Lions": "Detroit Lions",
    "Packers": "Green Bay Packers",
    "Panthers": "Carolina Panthers",
    "Patriots": "New England Patriots",
    "Raiders": "Las Vegas Raiders",
    "Rams": "Los Angeles Rams",
    "Ravens": "Baltimore Ravens",
    "Saints": "New Orleans Saints",
    "Seahawks": "Seattle Seahawks",
    "Steelers": "Pittsburgh Steelers",
    "Texans": "Houston Texans",
    "Titans": "Tennessee Titans",
    "Vikings": "Minnesota Vikings",
}


def normalize_team_name(name: str, sport: str) -> str:
    """Normalize team name to match Odds API format."""
    name = name.strip()
    if sport == "nfl":
        if name in NFL_TEAM_NAMES:
            return NFL_TEAM_NAMES[name]
        last_word = name.split()[-1] if name else ""
        if last_word in NFL_TEAM_NAMES:
            return NFL_TEAM_NAMES[last_word]
        for short, full in NFL_TEAM_NAMES.items():
            if short.lower() in name.lower():
                return full
    return name


# =============================================================================
# API CLIENT
# =============================================================================


def fetch_games(session: requests.Session, sport: str) -> list:
    """Fetch game list from BFA API.

    Tries the popular endpoint first, falls back to league-based endpoint.
    Returns list of game objects (type=1 regular games only).
    """
    config = SPORT_CONFIGS[sport]

    # Try popular endpoint first
    url = f"{API_BASE}/oddsservice/events/popular/{config['api_slug']}"
    resp = session.get(url, params={**API_PARAMS, "set": "Auto"}, timeout=15)
    resp.raise_for_status()
    data = resp.json()
    games = [g for g in data.get("games", []) if g.get("type") == 1]

    # Fall back to league endpoint if popular returns empty
    if not games and config.get("league_id"):
        url = f"{API_BASE}/oddsservice/events/leagues"
        resp = session.get(
            url,
            params={**API_PARAMS, "id": config["league_id"], "set": "Auto"},
            timeout=15,
        )
        resp.raise_for_status()
        data = resp.json()
        games = [g for g in data.get("games", []) if g.get("type") == 1]

    return games


def fetch_game_detail(session: requests.Session, game_id: int) -> dict:
    """Fetch full game detail (all markets, periods, alt lines)."""
    url = f"{API_BASE}/oddsservice/event/{game_id}"
    resp = session.get(url, params=API_PARAMS, timeout=15)
    resp.raise_for_status()
    return resp.json()


# =============================================================================
# JSON PARSING
# =============================================================================


def parse_game_odds(game: dict, config: dict, game_index: int,
                    away_team: str, home_team: str, fetch_time) -> list:
    """Parse API game detail JSON into DuckDB record format.

    Returns list of dicts matching the existing 18-column DuckDB schema.
    """
    records = []

    fixture = game.get("fixtures", [{}])[0]
    contestants = fixture.get("contestants", [])
    if len(contestants) < 2:
        return records

    # Determine away/home by contestant order (first listed = away, matching Playwright scraper)
    away_id = contestants[0]["id"]
    home_id = contestants[1]["id"]

    # Parse game date/time from ISO format
    game_dt = fixture.get("date", "")
    game_date = game_dt.split("T")[0] if game_dt else ""
    game_time = game_dt.split("T")[1].replace("Z", "") if "T" in game_dt else ""

    # Build period name → number mapping from the game data
    period_name_to_num = {}
    for p in game.get("periods", []):
        period_name_to_num[p["name"]] = p["number"]

    markets = game.get("markets", [])

    for period_short in config["periods"]:
        period_display = PERIOD_DISPLAY[period_short]
        period_number = period_name_to_num.get(period_display)
        if period_number is None:
            continue

        period_standard = PERIOD_TO_STANDARD[period_short]
        market_suffix = PERIOD_TO_MARKET_SUFFIX[period_short]

        # Filter markets for this period
        period_markets = [m for m in markets if m["periodNumber"] == period_number]

        # Find market objects by type
        spread_market = next((m for m in period_markets if m["type"] == 2), None)
        total_market = next((m for m in period_markets if m["type"] == 3), None)
        ml_market = next((m for m in period_markets if m["type"] == 1), None)

        # Parse main spread (index=0)
        away_spread = home_spread = away_spread_price = home_spread_price = None
        if spread_market:
            for odd in spread_market.get("odds", []):
                if odd["index"] == 0 and odd.get("status", 0) == 0:
                    if odd.get("contestantId") == away_id:
                        away_spread = odd["line"]
                        away_spread_price = odd["price"]
                    elif odd.get("contestantId") == home_id:
                        home_spread = odd["line"]
                        home_spread_price = odd["price"]

        # Parse main total (index=0)
        total = over_price = under_price = None
        if total_market:
            for odd in total_market.get("odds", []):
                if odd["index"] == 0 and odd.get("status", 0) == 0:
                    if odd.get("side") == 4:  # over
                        total = odd["line"]
                        over_price = odd["price"]
                    elif odd.get("side") == 5:  # under
                        under_price = odd["price"]

        # Parse main moneyline (index=0)
        away_ml = home_ml = None
        if ml_market:
            for odd in ml_market.get("odds", []):
                if odd["index"] == 0 and odd.get("status", 0) == 0:
                    if odd.get("contestantId") == away_id:
                        away_ml = odd["price"]
                    elif odd.get("contestantId") == home_id:
                        home_ml = odd["price"]

        game_base = {
            "fetch_time": fetch_time,
            "sport_key": config["sport_key"],
            "game_date": game_date,
            "game_time": game_time,
            "away_team": away_team,
            "home_team": home_team,
            "period": period_standard,
        }

        # Add main record
        if away_spread is not None or total is not None or away_ml is not None:
            records.append({
                **game_base,
                "game_id": f"bfa-{game_index}-{period_short}",
                "market": f"spreads_{market_suffix}",
                "away_spread": away_spread,
                "away_spread_price": away_spread_price,
                "home_spread": home_spread,
                "home_spread_price": home_spread_price,
                "total": total,
                "over_price": over_price,
                "under_price": under_price,
                "away_ml": away_ml,
                "home_ml": home_ml,
            })

        # Alt spreads (index > 0)
        if spread_market:
            alt_indices = sorted(set(
                o["index"] for o in spread_market.get("odds", []) if o["index"] > 0
            ))
            alt_count = 0
            for idx in alt_indices:
                alt_away_spread = alt_away_price = alt_home_spread = alt_home_price = None
                for odd in spread_market["odds"]:
                    if odd["index"] == idx and odd.get("status", 0) == 0:
                        if odd.get("contestantId") == away_id:
                            alt_away_spread = odd["line"]
                            alt_away_price = odd["price"]
                        elif odd.get("contestantId") == home_id:
                            alt_home_spread = odd["line"]
                            alt_home_price = odd["price"]

                # Sanity check: alt spread within 15 pts of main
                if alt_away_spread is not None and away_spread is not None:
                    if abs(alt_away_spread - away_spread) > 15:
                        continue

                if alt_away_spread is not None:
                    records.append({
                        **game_base,
                        "game_id": f"bfa-{game_index}-{period_short}-alt-s{alt_count}",
                        "market": f"alternate_spreads_{market_suffix}",
                        "away_spread": alt_away_spread,
                        "away_spread_price": alt_away_price,
                        "home_spread": alt_home_spread,
                        "home_spread_price": alt_home_price,
                        "total": None,
                        "over_price": None,
                        "under_price": None,
                        "away_ml": None,
                        "home_ml": None,
                    })
                    alt_count += 1

        # Alt totals (index > 0)
        if total_market:
            alt_indices = sorted(set(
                o["index"] for o in total_market.get("odds", []) if o["index"] > 0
            ))
            alt_count = 0
            for idx in alt_indices:
                alt_total = alt_over_price = alt_under_price = None
                for odd in total_market["odds"]:
                    if odd["index"] == idx and odd.get("status", 0) == 0:
                        if odd.get("side") == 4:
                            alt_total = odd["line"]
                            alt_over_price = odd["price"]
                        elif odd.get("side") == 5:
                            alt_under_price = odd["price"]

                # Sanity check: alt total within 25 pts of main
                if alt_total is not None and total is not None:
                    if abs(alt_total - total) > 25:
                        continue

                if alt_total is not None:
                    records.append({
                        **game_base,
                        "game_id": f"bfa-{game_index}-{period_short}-alt-t{alt_count}",
                        "market": f"alternate_totals_{market_suffix}",
                        "away_spread": None,
                        "away_spread_price": None,
                        "home_spread": None,
                        "home_spread_price": None,
                        "total": alt_total,
                        "over_price": alt_over_price,
                        "under_price": alt_under_price,
                        "away_ml": None,
                        "home_ml": None,
                    })
                    alt_count += 1

        # Team totals: type 4 = home, type 5 = away
        home_tt_market = next((m for m in period_markets if m["type"] == 4), None)
        away_tt_market = next((m for m in period_markets if m["type"] == 5), None)

        for tt_market, tt_side in [(home_tt_market, "home"), (away_tt_market, "away")]:
            if tt_market is None:
                continue
            tt_line = tt_over_price = tt_under_price = None
            for odd in tt_market.get("odds", []):
                if odd["index"] == 0 and odd.get("status", 0) == 0:
                    if odd.get("side") == 4:  # over
                        tt_line = odd["line"]
                        tt_over_price = odd["price"]
                    elif odd.get("side") == 5:  # under
                        tt_under_price = odd["price"]

            if tt_line is not None:
                records.append({
                    **game_base,
                    "game_id": f"bfa-{game_index}-{period_short}-tt-{tt_side}",
                    "market": f"team_totals_{tt_side}_{market_suffix}",
                    "away_spread": None,
                    "away_spread_price": None,
                    "home_spread": None,
                    "home_spread_price": None,
                    "total": tt_line,
                    "over_price": tt_over_price,
                    "under_price": tt_under_price,
                    "away_ml": None,
                    "home_ml": None,
                })

    return records


# =============================================================================
# DATABASE
# =============================================================================


def init_database(sport: str):
    """Initialize DuckDB with the odds table for a sport."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            fetch_time TIMESTAMP,
            sport_key VARCHAR,
            game_id VARCHAR,
            game_date VARCHAR,
            game_time VARCHAR,
            away_team VARCHAR,
            home_team VARCHAR,
            market VARCHAR,
            period VARCHAR,
            away_spread FLOAT,
            away_spread_price INTEGER,
            home_spread FLOAT,
            home_spread_price INTEGER,
            total FLOAT,
            over_price INTEGER,
            under_price INTEGER,
            away_ml INTEGER,
            home_ml INTEGER
        )
    """)
    conn.close()


def save_to_database(sport: str, odds_data: list):
    """Save scraped odds to DuckDB."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    columns = [
        "fetch_time", "sport_key", "game_id", "game_date", "game_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml"
    ]

    placeholders = ", ".join(["?" for _ in columns])

    # Clear old data before inserting fresh scrape
    conn.execute(f"DELETE FROM {table_name}")

    conn.executemany(f"""
        INSERT INTO {table_name} ({", ".join(columns)})
        VALUES ({placeholders})
    """, [
        tuple(d[col] for col in columns)
        for d in odds_data
    ])

    result = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()
    print(f"Database now has {result[0]} total records in {table_name}")

    conn.close()


# =============================================================================
# MAIN SCRAPER
# =============================================================================


def scrape_bfa(sport: str, headless: bool = True):
    """Scrape BFA derivative odds via REST API and save to DuckDB.

    The headless parameter is kept for backward compatibility but is ignored —
    no browser is needed.
    """
    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    fetch_time = datetime.now(timezone.utc)
    session = requests.Session()
    session.headers.update(API_HEADERS)

    # Step 1: Fetch game list
    print(f"Fetching {sport.upper()} games from BFA API...")
    games = fetch_games(session, sport)
    print(f"Found {len(games)} games")

    if not games:
        print("No games found.")
        return []

    # Load team name resolution
    team_dict = load_team_dict(sport) if sport != "nfl" else {}
    canonical_games = load_canonical_games(sport) if sport != "nfl" else []

    all_odds = []

    for i, game in enumerate(games):
        fixture = game.get("fixtures", [{}])[0]
        contestants = fixture.get("contestants", [])
        if len(contestants) < 2:
            continue

        # Get raw team names (first listed = away, second = home)
        away_raw = contestants[0]["name"]
        home_raw = contestants[1]["name"]

        # Resolve team names
        if team_dict or canonical_games:
            away_team, home_team = resolve_team_names(
                away_raw, home_raw, team_dict, canonical_games
            )
        else:
            away_team = normalize_team_name(away_raw, sport)
            home_team = normalize_team_name(home_raw, sport)

        print(f"  [{i+1}/{len(games)}] {away_team} @ {home_team}", end="")

        try:
            # Step 2: Fetch game detail for all markets
            detail = fetch_game_detail(session, game["id"])
            records = parse_game_odds(detail, config, i, away_team, home_team, fetch_time)
            all_odds.extend(records)

            # Count main vs alt records
            main_count = sum(1 for r in records if "alt" not in r["game_id"])
            alt_count = len(records) - main_count
            print(f" — {main_count} main, {alt_count} alts")

        except Exception as e:
            print(f" — Error: {e}")

    # Save to DuckDB
    if all_odds:
        save_to_database(sport, all_odds)

    print(f"\nScraped {len(all_odds)} total records for {sport.upper()}")
    return all_odds


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    scrape_bfa(sport)
