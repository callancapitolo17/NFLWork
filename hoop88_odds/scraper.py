#!/usr/bin/env python3
"""
Hoop88 Odds Scraper - REST API
Fetches odds via /cloud/api/Lines/Get_LeagueLines2 JSON endpoint.
No browser required — auth via JWT, data via POST requests.
"""

import os
import sys
import duckdb
import requests
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from zoneinfo import ZoneInfo

from dotenv import load_dotenv

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Hoop88's API returns naive Pacific wall-clock in GameDateTime (confirmed
# by tools/TZ_AUDIT_FINDINGS.md — modal Δ vs sharp UTC was +7.02h = PDT).
# Prior code assumed ET, which silently shifted every game by 3 hours.
# We now parse as America/Los_Angeles and store the UTC instant as
# game_start_time TIMESTAMPTZ so downstream code never has to guess the TZ.
_HOOP88_TZ = ZoneInfo("America/Los_Angeles")


def _hoop88_game_start_time(game_dt: str) -> Optional[datetime]:
    """Parse H88's naive PT GameDateTime string into a UTC-aware datetime.

    game_dt: "YYYY-MM-DD HH:MM" or "YYYY-MM-DD HH:MM:SS"
    Returns: tz-aware UTC datetime, or None if malformed/empty.
    """
    if not game_dt or len(game_dt) < 10:
        return None
    try:
        date_part, _, time_part = game_dt.partition(" ")
        if not date_part or not time_part:
            return None
        y, m, d = map(int, date_part.split("-"))
        time_part = time_part[:5]  # truncate seconds if present
        hh, mm = map(int, time_part.split(":"))
    except (ValueError, AttributeError):
        return None
    naive_pt = datetime(y, m, d, hh, mm, tzinfo=_HOOP88_TZ)
    return naive_pt.astimezone(timezone.utc)

# Load .env from bet_logger directory
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")

DB_PATH = Path(__file__).parent / "hoop88.duckdb"

API_BASE = f"{HOOP88_URL}/cloud/api"

# Sport configurations
SPORT_CONFIGS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "sport_type": "FOOTBALL",
        "sport_sub_type": "NFL",
        "periods": [
            {"short": "1H", "api_name": "1st Half", "api_num": 1, "suffix": "h1", "standard": "Half1"},
            {"short": "1Q", "api_name": "1st Quarter", "api_num": 3, "suffix": "q1", "standard": "1"},
            {"short": "2Q", "api_name": "2nd Quarter", "api_num": 4, "suffix": "q2", "standard": "2"},
            {"short": "3Q", "api_name": "3rd Quarter", "api_num": 5, "suffix": "q3", "standard": "3"},
            {"short": "4Q", "api_name": "4th Quarter", "api_num": 6, "suffix": "q4", "standard": "4"},
        ],
    },
    "ncaaf": {
        "sport_key": "americanfootball_ncaaf",
        "table_name": "ncaaf_odds",
        "sport_type": "FOOTBALL",
        "sport_sub_type": "NCAA Football",
        "periods": [
            {"short": "1H", "api_name": "1st Half", "api_num": 1, "suffix": "h1", "standard": "Half1"},
            {"short": "1Q", "api_name": "1st Quarter", "api_num": 3, "suffix": "q1", "standard": "1"},
        ],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "sport_type": "BASKETBALL",
        "sport_sub_type": "NCAA",
        "periods": [
            {"short": "1H", "api_name": "1st Half", "api_num": 1, "suffix": "h1", "standard": "Half1"},
        ],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "sport_type": "BASKETBALL",
        "sport_sub_type": "NBA",
        "periods": [
            {"short": "1H", "api_name": "1st Half", "api_num": 1, "suffix": "h1", "standard": "Half1"},
            {"short": "1Q", "api_name": "1st Quarter", "api_num": 3, "suffix": "q1", "standard": "1"},
        ],
    },
    "college_baseball": {
        "sport_key": "baseball_ncaa",
        "table_name": "college_baseball_odds",
        "sport_type": "BASEBALL",
        "sport_sub_type": "",  # empty returns all baseball; filter by rot num to get college only
        "periods": [
            {"short": "FG", "api_name": "Game", "api_num": 0, "suffix": "fg", "standard": "Full"},
        ],
        "rot_range": (3001, 3999),  # college baseball rotation numbers are in 3xxx range
    },
    "mlb": {
        "sport_key": "baseball_mlb",
        "table_name": "mlb_odds",
        "sport_type": "BASEBALL",
        "sport_sub_type": "MLB",
        "periods": [
            {"short": "FG", "api_name": "Game", "api_num": 0, "suffix": "fg", "standard": "Full"},
            {"short": "F5", "api_name": "1st 5 Innings", "api_num": 1, "suffix": "f5", "standard": "F5"},
        ],
        "rot_range": (901, 999),  # MLB rotation numbers
    },
}

# NFL team name mapping (Hoop88 short names -> Odds API full names)
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
    if sport == "nfl":
        return NFL_TEAM_NAMES.get(name, name)
    return name


# =============================================================================
# AUTH
# =============================================================================


def login(session: requests.Session) -> str:
    """Login to Hoop88 and return JWT token.

    POST /cloud/api/System/authenticateCustomer with credentials.
    Returns JWT on success (status 200). Status 204 = bad credentials.
    """
    resp = session.post(f"{API_BASE}/System/authenticateCustomer", data={
        "customerID": HOOP88_USERNAME,
        "password": HOOP88_PASSWORD,
        "state": "true",
        "multiaccount": "1",
        "response_type": "code",
        "client_id": HOOP88_USERNAME,
        "domain": HOOP88_URL.replace("https://", "").replace("http://", ""),
        "redirect_uri": HOOP88_URL.replace("https://", "").replace("http://", ""),
        "operation": "authenticateCustomer",
        "RRO": "1",
    }, timeout=15)

    if resp.status_code == 204:
        raise RuntimeError("Login failed — bad credentials (HTTP 204)")
    resp.raise_for_status()

    data = resp.json()
    token = data.get("code")
    if not token:
        raise RuntimeError(f"Login succeeded but no token found in response: {list(data.keys())}")

    return token


# =============================================================================
# API CLIENT
# =============================================================================


def fetch_lines(session: requests.Session, config: dict, period: dict) -> list:
    """Fetch lines for a sport/period from Get_LeagueLines2."""
    resp = session.post(f"{API_BASE}/Lines/Get_LeagueLines2", data={
        "sportType": config["sport_type"],
        "sportSubType": config["sport_sub_type"],
        "period": period["api_name"],
        "periodNumber": period["api_num"],
        "propDescription": "Game",
        "wagerType": "Straight",
        "office": "COINDEVIL",
        "customerID": HOOP88_USERNAME,
        "hourFilter": "0",
        "keyword": "",
        "correlationID": "",
        "grouping": "",
        "periods": "0",
        "rotOrder": "0",
        "placeLateFlag": "false",
        "RRO": "1",
        "agentSite": "0",
        "operation": "Get_LeagueLines2",
    }, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    return data.get("Lines", [])


# =============================================================================
# JSON PARSING
# =============================================================================


def parse_lines(lines: list, config: dict, period: dict,
                team_dict: dict, canonical_games: list,
                fetch_time: str, sport: str) -> list[dict]:
    """Parse Get_LeagueLines2 JSON response into DuckDB records."""
    records = []
    sport_key = config["sport_key"]
    suffix = period["suffix"]
    period_standard = period["standard"]

    rot_range = config.get("rot_range")

    for game in lines:
        if game.get("Status") != "O":
            continue

        # Filter by rotation number range (e.g., college baseball = 3xxx)
        if rot_range:
            rot1 = game.get("Team1RotNum", 0) or 0
            rot2 = game.get("Team2RotNum", 0) or 0
            if not (rot_range[0] <= rot1 <= rot_range[1] and
                    rot_range[0] <= rot2 <= rot_range[1]):
                continue

        # Team names: Team1 = away, Team2 = home
        away_raw = (game.get("Team1ID") or "").strip()
        home_raw = (game.get("Team2ID") or "").strip()
        if not away_raw or not home_raw:
            continue

        # Resolve team names
        if team_dict or canonical_games:
            away_team, home_team = resolve_team_names(
                away_raw, home_raw, team_dict, canonical_games
            )
        else:
            away_team = normalize_team_name(away_raw, sport)
            home_team = normalize_team_name(home_raw, sport)

        # Game start time — parse H88's naive PT wall-clock into UTC instant.
        game_start_time = _hoop88_game_start_time(game.get("GameDateTime", ""))

        away_rot = str(game.get("Team1RotNum", ""))
        home_rot = str(game.get("Team2RotNum", ""))
        game_id = f"{away_rot}-{home_rot}"

        # Spread: Spread field is always positive, FavoredTeamID tells direction
        spread_val = game.get("Spread")
        favored = (game.get("FavoredTeamID") or "").strip()
        spread_adj1 = game.get("SpreadAdj1")  # Team1 spread juice
        spread_adj2 = game.get("SpreadAdj2")  # Team2 spread juice

        away_spread = home_spread = None
        away_spread_price = home_spread_price = None
        if spread_val is not None and spread_val != 0:
            if favored == away_raw:
                away_spread = -abs(spread_val)
                home_spread = abs(spread_val)
            else:
                away_spread = abs(spread_val)
                home_spread = -abs(spread_val)
            away_spread_price = spread_adj1
            home_spread_price = spread_adj2

        # Moneyline
        away_ml = game.get("MoneyLine1")
        home_ml = game.get("MoneyLine2")

        # Total
        total = game.get("TotalPoints")
        over_price = game.get("TtlPtsAdj1")
        under_price = game.get("TtlPtsAdj2")

        base = {
            "fetch_time": fetch_time,
            "sport_key": sport_key,
            "game_start_time": game_start_time,
            "away_team": away_team,
            "home_team": home_team,
        }

        # Main line record (spread + total + ML)
        if away_spread is not None or total is not None or away_ml is not None:
            records.append({
                **base,
                "game_id": game_id,
                "market": f"spreads_{suffix}",
                "period": period_standard,
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

        # Team totals (away = Team1, home = Team2)
        t1_total = game.get("Team1TotalPoints")
        t2_total = game.get("Team2TotalPoints")

        if t1_total is not None:
            records.append({
                **base,
                "game_id": f"{game_id}-tt-away",
                "market": f"team_totals_away_{suffix}",
                "period": period_standard,
                "away_spread": None,
                "away_spread_price": None,
                "home_spread": None,
                "home_spread_price": None,
                "total": t1_total,
                "over_price": game.get("Team1TtlPtsAdj1"),
                "under_price": game.get("Team1TtlPtsAdj2"),
                "away_ml": None,
                "home_ml": None,
            })

        if t2_total is not None:
            records.append({
                **base,
                "game_id": f"{game_id}-tt-home",
                "market": f"team_totals_home_{suffix}",
                "period": period_standard,
                "away_spread": None,
                "away_spread_price": None,
                "home_spread": None,
                "home_spread_price": None,
                "total": t2_total,
                "over_price": game.get("Team2TtlPtsAdj1"),
                "under_price": game.get("Team2TtlPtsAdj2"),
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
            fetch_time TIMESTAMPTZ,
            sport_key VARCHAR,
            game_id VARCHAR,
            game_start_time TIMESTAMPTZ,
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
    """Save scraped odds to DuckDB atomically (stage-and-swap)."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    columns = [
        "fetch_time", "sport_key", "game_id", "game_start_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml"
    ]

    # Stage into a TEMP table cloned from the live schema, then atomically
    # swap. CREATE OR REPLACE TABLE ... AS SELECT is DuckDB's closest
    # equivalent of an atomic rename — readers see either the entire old
    # snapshot or the entire new one, never a half-written state. On an
    # empty scrape, leave the prior snapshot in place so consumers don't
    # see a momentarily empty table (matches BFA/DK/FD/WZ pattern).
    if odds_data:
        conn.execute(
            f"CREATE OR REPLACE TEMP TABLE {table_name}_new AS "
            f"SELECT * FROM {table_name} LIMIT 0"
        )
        placeholders = ", ".join(["?" for _ in columns])
        tuples = [tuple(d.get(c) for c in columns) for d in odds_data]
        conn.executemany(
            f"INSERT INTO {table_name}_new VALUES ({placeholders})", tuples
        )
        conn.execute(
            f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM {table_name}_new"
        )
        conn.execute(f"DROP TABLE IF EXISTS {table_name}_new")
    else:
        print(f"[hoop88] empty scrape — leaving prior snapshot in place", flush=True)

    result = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()
    print(f"Database now has {result[0]} total records in {table_name}")

    conn.close()


# =============================================================================
# MAIN
# =============================================================================


def scrape_hoop88(sport: str, headless: bool = True):
    """Scrape Hoop88 odds via REST API and save to DuckDB.

    The headless parameter is kept for backward compatibility but is ignored —
    no browser is needed.
    """
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    # Aware UTC datetime — DuckDB stores it directly as TIMESTAMPTZ. A
    # naive strftime string would be interpreted as local-TZ on insert
    # and silently shift (same bug DK/FD hit pre-migration).
    fetch_time = datetime.now(timezone.utc)

    # Load team name resolution
    team_dict = load_team_dict(sport) if sport != "nfl" else {}
    canonical_games = load_canonical_games(sport) if sport != "nfl" else []

    session = requests.Session()
    session.headers.update({
        "X-Requested-With": "XMLHttpRequest",
        "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
    })

    # Step 1: Get Cloudflare cookie by loading the page first
    print("Getting session cookies...")
    session.get(HOOP88_URL, timeout=15)

    # Step 2: Authenticate and get JWT
    print("Logging in to Hoop88...")
    token = login(session)
    session.headers["Authorization"] = f"Bearer {token}"
    print("Logged in successfully")

    # Step 3: Fetch lines for each period
    all_odds = []

    for period in config["periods"]:
        print(f"Fetching {sport.upper()} {period['short']}...")
        lines = fetch_lines(session, config, period)
        print(f"  Found {len(lines)} games")

        records = parse_lines(
            lines, config, period,
            team_dict, canonical_games,
            fetch_time, sport
        )
        all_odds.extend(records)

    # Summary by market
    market_counts = {}
    for rec in all_odds:
        m = rec["market"]
        market_counts[m] = market_counts.get(m, 0) + 1
    for m, c in sorted(market_counts.items()):
        print(f"  {m}: {c} records")

    # Step 4: Save to DuckDB (always save, even if empty, to clear stale data)
    save_to_database(sport, all_odds)

    print(f"\nScraped {len(all_odds)} total records for {sport.upper()}")
    return all_odds


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    print(f"Starting Hoop88 {sport.upper()} odds scraper (REST API)...")
    scrape_hoop88(sport)
