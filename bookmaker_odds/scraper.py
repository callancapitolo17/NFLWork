#!/usr/bin/env python3
"""
Bookmaker.eu Odds Scraper
Fetches odds from Bookmaker.eu's internal API via curl_cffi (Chrome TLS fingerprint).
No browser needed — uses saved cookies from recon_bookmaker.py.

Usage:
    python scraper.py cbb
    python scraper.py nba

If cookies expire, run recon_bookmaker.py to refresh them.
"""

import json
import os
import sys
import duckdb
from datetime import datetime, timezone
from pathlib import Path

from curl_cffi import requests as cffi_requests
from dotenv import load_dotenv

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Load .env from bet_logger directory
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

BOOKMAKER_USERNAME = os.getenv("BOOKMAKER_USERNAME")
BOOKMAKER_PASSWORD = os.getenv("BOOKMAKER_PASSWORD")

SITE_URL = "https://be.bookmaker.eu/en/sports/"
API_BASE = "https://be.bookmaker.eu/gateway/BetslipProxy.aspx"
COOKIE_PATH = Path(__file__).parent / ".bookmaker_cookies.json"
DB_PATH = Path(__file__).parent / "bookmaker.duckdb"

# Sport configurations — league IDs discovered via recon
SPORT_CONFIGS = {
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "leagues": [
            {"id": "4", "name": "Game Lines", "market": "spreads", "period": "fg"},
            {"id": "13554", "name": "Extra Games", "market": "spreads", "period": "fg"},
            {"id": "205", "name": "1st Halves", "market": "spreads_h1", "period": "Half1"},
        ],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "leagues": [
            {"id": "3", "name": "Game Lines", "market": "spreads", "period": "fg"},
            {"id": "202", "name": "1st Halves", "market": "spreads_h1", "period": "Half1"},
        ],
    },
    "mlb": {
        "sport_key": "baseball_mlb",
        "table_name": "mlb_odds",
        "leagues": [
            {"id": "5", "name": "Game Lines", "market": "spreads", "period": "fg"},
            {"id": "6", "name": "1st 5 Innings", "market": "spreads_f5", "period": "F5"},
            {"id": "503", "name": "1st 3 Innings", "market": "spreads_f3", "period": "F3"},
        ],
    },
}


# =============================================================================
# API CLIENT (curl_cffi)
# =============================================================================


def _make_schedule_body(league_id: str) -> dict:
    return {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "LeaguesIdList": league_id,
                    "LanguageId": "0",
                    "LineStyle": "E",
                    "ScheduleType": "american",
                    "LinkDeriv": "true",
                }
            }
        }
    }


def _save_cookies(session):
    """Save session cookies to disk for reuse."""
    cookies = dict(session.cookies)
    with open(COOKIE_PATH, "w") as f:
        json.dump(cookies, f)


def _load_cookies(session):
    """Load saved cookies into session."""
    if not COOKIE_PATH.exists():
        return False
    try:
        with open(COOKIE_PATH) as f:
            cookies = json.load(f)
        for name, value in cookies.items():
            session.cookies.set(name, value, domain=".bookmaker.eu")
        return bool(cookies)
    except (json.JSONDecodeError, OSError):
        return False


def _create_session() -> cffi_requests.Session:
    """Create a curl_cffi session that impersonates Chrome."""
    session = cffi_requests.Session(impersonate="chrome")
    _load_cookies(session)
    return session


def fetch_schedule(session, league_id: str) -> dict | None:
    """Fetch GetSchedule via curl_cffi. Returns None if blocked."""
    try:
        resp = session.post(
            f"{API_BASE}/GetSchedule",
            json=_make_schedule_body(league_id),
            timeout=15,
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        if "Schedule" not in data:
            return None
        return data
    except Exception:
        return None


def login(session, username: str, password: str) -> bool:
    """Login via curl_cffi. Returns True on success."""
    body = {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "Player": username,
                    "Password": password,
                    "loginKey": "",
                }
            }
        }
    }
    try:
        resp = session.post(f"{API_BASE}/Login", json=body, timeout=15)
        return resp.status_code == 200
    except Exception:
        return False


def _has_games(data: dict | None) -> bool:
    """Check if GetSchedule response contains games."""
    if not data:
        return False
    return bool(
        data.get("Schedule", {})
        .get("Data", {})
        .get("Leagues", {})
        .get("League", [])
    )


def refresh_cookies():
    """Run recon_bookmaker.py to get fresh Cloudflare cookies."""
    import subprocess
    recon_script = Path(__file__).parent / "recon_bookmaker.py"
    python = sys.executable
    print(f"\nRunning recon script to refresh cookies...")
    print(f"  {python} {recon_script}\n")
    subprocess.run([python, str(recon_script)], check=True)

    # Extract cookies from recon output into our format
    recon_cookies = Path(__file__).parent / "recon_bookmaker_cookies.json"
    if recon_cookies.exists():
        with open(recon_cookies) as f:
            pw_cookies = json.load(f)
        cookie_dict = {}
        for c in pw_cookies:
            if "bookmaker.eu" in c.get("domain", ""):
                cookie_dict[c["name"]] = c["value"]
        with open(COOKIE_PATH, "w") as f:
            json.dump(cookie_dict, f)
        print(f"Saved {len(cookie_dict)} cookies from recon.")


# =============================================================================
# JSON PARSING
# =============================================================================


def parse_schedule(schedule_data: dict, league_config: dict, sport_key: str,
                   team_dict: dict, canonical_games: list,
                   fetch_time: str) -> list[dict]:
    """Parse GetSchedule response into 18-column DuckDB records.

    Response structure: Schedule.Data.Leagues.League[].dateGroup[].game[]
    Each game has Derivatives.line[] where index "0" is the main line.
    """
    records = []
    market = league_config["market"]
    period = league_config["period"]

    leagues = (schedule_data
               .get("Schedule", {})
               .get("Data", {})
               .get("Leagues", {})
               .get("League", []))

    for league in leagues:
        for date_group in league.get("dateGroup", []):
            for game in date_group.get("game", []):
                record = _parse_game(
                    game, market, period, sport_key,
                    team_dict, canonical_games, fetch_time
                )
                if record:
                    records.append(record)

    return records


def _parse_game(game: dict, market: str, period: str, sport_key: str,
                team_dict: dict, canonical_games: list,
                fetch_time: str) -> dict | None:
    """Parse a single game object into a DuckDB record."""
    if game.get("Stat") != "O":
        return None

    away_raw = (game.get("vtm") or "").strip()
    home_raw = (game.get("htm") or "").strip()
    if not away_raw or not home_raw:
        return None

    # Resolve team names via shared canonical_match
    if team_dict or canonical_games:
        away_team, home_team = resolve_team_names(
            away_raw, home_raw, team_dict, canonical_games
        )
    else:
        away_team, home_team = away_raw, home_raw

    # Validate resolved teams against canonical games to filter out
    # cross-sport contamination (e.g., UFC fighters in MLB league 206)
    if canonical_games:
        matched = any(
            (away_team == cg["away_team"] and home_team == cg["home_team"]) or
            (away_team == cg["home_team"] and home_team == cg["away_team"])
            for cg in canonical_games
        )
        if not matched:
            return None

    # Game date (YYYYMMDD → MM/DD) and time (HH:MM:SS → HH:MM)
    gmdt = game.get("gmdt", "")
    gmtm = game.get("gmtm", "")
    game_date = f"{gmdt[4:6]}/{gmdt[6:8]}" if len(gmdt) >= 8 else ""
    game_time = gmtm[:5] if gmtm else ""

    game_id = str(game.get("idgm", ""))

    # Get main line (index == "0")
    lines = game.get("Derivatives", {}).get("line", [])
    main_line = None
    for line in lines:
        if str(line.get("index")) == "0":
            main_line = line
            break

    if not main_line:
        return None

    # Parse spread (check s_sp status flag)
    away_spread = away_spread_price = home_spread = home_spread_price = None
    if main_line.get("s_sp") == 1:
        try:
            away_spread = float(main_line["vsprdt"])
            home_spread = float(main_line["hsprdt"])
            away_spread_price = int(main_line["vsprdoddst"])
            home_spread_price = int(main_line["hsprdoddst"])
        except (KeyError, ValueError, TypeError):
            pass

    # Parse total (check s_tot status flag)
    total = over_price = under_price = None
    if main_line.get("s_tot") == 1:
        try:
            total = float(main_line["ovt"])
            over_price = int(main_line["ovoddst"])
            under_price = int(main_line["unoddst"])
        except (KeyError, ValueError, TypeError):
            pass

    # Parse moneyline (check s_ml status flag)
    away_ml = home_ml = None
    if main_line.get("s_ml") == 1:
        try:
            away_ml = int(main_line["voddst"])
            home_ml = int(main_line["hoddst"])
        except (KeyError, ValueError, TypeError):
            pass

    # Skip if no odds available at all
    if away_spread is None and total is None and away_ml is None:
        return None

    return {
        "fetch_time": fetch_time,
        "sport_key": sport_key,
        "game_id": game_id,
        "game_date": game_date,
        "game_time": game_time,
        "away_team": away_team,
        "home_team": home_team,
        "market": market,
        "period": period,
        "away_spread": away_spread,
        "away_spread_price": away_spread_price,
        "home_spread": home_spread,
        "home_spread_price": home_spread_price,
        "total": total,
        "over_price": over_price,
        "under_price": under_price,
        "away_ml": away_ml,
        "home_ml": home_ml,
    }


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
    """Save scraped odds to DuckDB (replaces previous scrape)."""
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
# MAIN
# =============================================================================


def scrape_bookmaker(sport: str):
    """Scrape Bookmaker.eu odds via curl_cffi and save to DuckDB.

    Uses saved cookies from recon_bookmaker.py for Cloudflare bypass.
    If cookies are expired, launches recon script to refresh them.
    """
    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    # Load team name resolution
    team_dict = load_team_dict(sport)
    canonical_games = load_canonical_games(sport)

    session = _create_session()

    # Hit the site to refresh Cloudflare cookies
    try:
        resp = session.get(SITE_URL, timeout=15)
        blocked = resp.status_code == 403
    except Exception:
        blocked = True

    # Test with first league
    first_league = config["leagues"][0]
    test_data = None if blocked else fetch_schedule(session, first_league["id"])

    # If blocked or no data, try refreshing cookies via recon
    if not _has_games(test_data):
        if blocked:
            print("Cloudflare blocked request — cookies expired.")
        else:
            # Not blocked but no games — try login first
            if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
                print("No games returned, logging in...")
                if login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD):
                    test_data = fetch_schedule(session, first_league["id"])

        # Still no games? Refresh cookies via recon
        if not _has_games(test_data):
            refresh_cookies()
            session = _create_session()
            resp = session.get(SITE_URL, timeout=15)
            if resp.status_code == 403:
                print("ERROR: Still blocked after recon. Check Cloudflare manually.")
                return []
            if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
                login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
            test_data = fetch_schedule(session, first_league["id"])
            if not _has_games(test_data):
                print("ERROR: No games returned after recon + login.")
                return []

    # Fetch all leagues
    all_odds = []
    for i, league in enumerate(config["leagues"]):
        print(f"Fetching {league['name']} (league {league['id']})...")
        data = test_data if i == 0 else fetch_schedule(session, league["id"])
        if data is None:
            print(f"  Failed to fetch league {league['id']}")
            continue

        records = parse_schedule(
            data, league, config["sport_key"],
            team_dict, canonical_games, fetch_time
        )
        all_odds.extend(records)
        print(f"  Parsed {len(records)} records")

    # Save cookies for next run
    _save_cookies(session)

    # Summary by market
    market_counts = {}
    for rec in all_odds:
        m = rec["market"]
        market_counts[m] = market_counts.get(m, 0) + 1
    for m, c in sorted(market_counts.items()):
        print(f"  {m}: {c} records")

    # Save to DuckDB
    if all_odds:
        save_to_database(sport, all_odds)

    print(f"\nScraped {len(all_odds)} total records for {sport.upper()}")
    return all_odds


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "cbb"
    print(f"Starting Bookmaker.eu {sport.upper()} odds scraper...")
    scrape_bookmaker(sport)
