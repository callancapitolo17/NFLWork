#!/usr/bin/env python3
"""
Hoop88 Odds Scraper for Answer Key Pipeline
Scrapes derivative markets (1H, 1Q, 2Q, 3Q, 4Q) and saves to DuckDB
"""

import os
import re
import sys
import duckdb
import sys
from datetime import datetime, timezone
from pathlib import Path
from dotenv import load_dotenv

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

# Load .env from bet_logger directory
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")

DB_PATH = Path(__file__).parent / "hoop88.duckdb"

# Sport configurations
SPORT_CONFIGS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "sport_category": "FOOTBALL",
        "league_selector": "NFL",
        "periods": ["1H", "1Q", "2Q", "3Q", "4Q"],  # Derivatives only
    },
    "ncaaf": {
        "sport_key": "americanfootball_ncaaf",
        "table_name": "ncaaf_odds",
        "sport_category": "FOOTBALL",
        "league_selector": "NCAA Football",
        "periods": ["1H", "1Q"],  # NCAAF has fewer periods on Hoop88
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "sport_category": "BASKETBALL",
        "league_selector": "NCAA Basketball",
        "periods": ["1H"],
    },
}

# Map period short names to Hoop88 UI names
PERIOD_MAP = {
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


def parse_spread_odds(raw_text: str) -> tuple:
    """Parse spread and odds from text like '-5½ -110'."""
    if not raw_text:
        return (None, None)

    text = raw_text
    text = re.sub(r'(\d)½', r'\g<1>.5', text)
    text = re.sub(r'(\d)¼', r'\g<1>.25', text)
    text = re.sub(r'(\d)¾', r'\g<1>.75', text)
    text = re.sub(r'([+-])½', r'\g<1>0.5', text)
    text = re.sub(r'([+-])¼', r'\g<1>0.25', text)
    text = re.sub(r'([+-])¾', r'\g<1>0.75', text)
    text = text.replace('½', '0.5').replace('¼', '0.25').replace('¾', '0.75')

    match = re.search(r'([+-]?\d+\.?\d*)\s*([+-]\d+)', text)
    if match:
        return (float(match.group(1)), int(match.group(2)))
    return (None, None)


def parse_total_odds(raw_text: str) -> tuple:
    """Parse total and odds from text like 'O 45½ -110'."""
    if not raw_text:
        return (None, None)

    text = raw_text
    text = text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')

    total_match = re.search(r'(\d+\.?\d*)', text)
    odds_match = re.search(r'([+-]\d+)', text)

    if total_match and odds_match:
        return (float(total_match.group(1)), int(odds_match.group(1)))
    return (None, None)


def parse_ml_odds(raw_text: str) -> int:
    """Parse moneyline odds from text like '+150' or '-200'."""
    if not raw_text:
        return None
    match = re.search(r'([+-]\d+)', raw_text)
    return int(match.group(1)) if match else None


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


def navigate_to_period(page, league_selector: str, period_name: str) -> bool:
    """Navigate to a specific period in Hoop88 using data-sport-sub-type selectors."""
    try:
        # For periods (1st Half, 1st Quarter, etc.), find the one belonging to our league
        # The sidebar is ordered: NFL -> 1st Half (NFL) -> 1st Quarter (NFL) -> NCAA Football -> ...
        clicked = page.evaluate(f'''() => {{
            const allElements = Array.from(document.querySelectorAll('[data-sport-sub-type]'));
            const leagueSelector = "{league_selector}";
            const periodSelector = "{period_name}";

            // Find the index of our league
            let leagueIdx = -1;
            for (let i = 0; i < allElements.length; i++) {{
                if (allElements[i].getAttribute('data-sport-sub-type') === leagueSelector) {{
                    leagueIdx = i;
                    break;
                }}
            }}

            if (leagueIdx === -1) {{
                console.log('League not found: ' + leagueSelector);
                return false;
            }}

            // Find the period after the league but before the next major section
            const majorSections = ['NFL', 'NCAA Football', 'NBA', 'NHL', 'NCAA Basketball'];

            for (let i = leagueIdx + 1; i < allElements.length; i++) {{
                const subType = allElements[i].getAttribute('data-sport-sub-type');

                // Stop if we hit another major section
                if (majorSections.includes(subType) && subType !== leagueSelector) {{
                    break;
                }}

                // Found our period!
                if (subType === periodSelector) {{
                    allElements[i].click();
                    return true;
                }}
            }}

            console.log('Period not found in league section: ' + periodSelector);
            return false;
        }}''')

        if not clicked:
            print(f"  Period '{period_name}' not found in {league_selector} section")
            return False

        page.wait_for_timeout(2500)
        return True

    except Exception as e:
        print(f"  Error navigating to {period_name}: {e}")
        return False


def scrape_hoop88(sport: str, headless: bool = True):
    """Scrape Hoop88 derivative odds and save to DuckDB."""
    if not HOOP88_USERNAME or not HOOP88_PASSWORD:
        raise ValueError("HOOP88_USERNAME and HOOP88_PASSWORD must be set in .env")

    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    fetch_time = datetime.now(tz=timezone.utc)
    all_odds = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {HOOP88_URL}...")
        page.goto(HOOP88_URL, wait_until='networkidle')

        # Login
        username_field = page.locator('input[name="customerID"]')
        if username_field.count() > 0 and username_field.is_visible():
            print("Logging in...")
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')
            page.wait_for_selector('input[name="customerID"]', state='hidden', timeout=15000)
            page.wait_for_load_state('networkidle')
            page.wait_for_timeout(1500)
            print("Login successful")

        # Navigate to sport category
        sport_cat = config["sport_category"]
        print(f"Navigating to {sport_cat}...")
        page.click(f'text={sport_cat}')
        page.wait_for_timeout(3000)

        # Load team name resolution data for non-NFL sports
        team_dict = load_team_dict(sport) if sport != "nfl" else {}
        canonical_games = load_canonical_games(sport) if sport != "nfl" else []

        for period_short in config["periods"]:
            period_full = PERIOD_MAP[period_short]
            period_standard = PERIOD_TO_STANDARD[period_short]

            print(f"Scraping {sport.upper()} {period_short}...")

            if not navigate_to_period(page, config["league_selector"], period_full):
                print(f"  Skipping - could not navigate to {period_full}")
                continue

            # Extract game data
            game_data = page.evaluate('''() => {
                const games = [];
                const linePanels = document.querySelectorAll('[data-panel="line"].GAME');

                linePanels.forEach((panel, idx) => {
                    try {
                        const firstTeam = panel.querySelector('[data-field="first-team"]');
                        const secondTeam = panel.querySelector('[data-field="second-team"]');
                        if (!firstTeam || !secondTeam) return;

                        const lineRows = panel.querySelectorAll('.lines .line');
                        if (lineRows.length < 2) return;

                        const awayRow = lineRows[0];
                        const homeRow = lineRows[1];

                        function getGroupText(row, groupClass) {
                            const group = row.querySelector('.' + groupClass);
                            if (!group) return '';
                            const span = group.querySelector('span');
                            if (span) return span.textContent.trim();
                            return '';
                        }

                        // Get game date/time from panel
                        const dateField = panel.querySelector('[data-field="date"]');
                        const timeField = panel.querySelector('[data-field="time"]');

                        games.push({
                            away_team: firstTeam.textContent.trim(),
                            home_team: secondTeam.textContent.trim(),
                            game_date: dateField ? dateField.textContent.trim() : '',
                            game_time: timeField ? timeField.textContent.trim() : '',
                            away_spread_raw: getGroupText(awayRow, 'group-2'),
                            home_spread_raw: getGroupText(homeRow, 'group-2'),
                            over_raw: getGroupText(awayRow, 'group-4'),
                            under_raw: getGroupText(homeRow, 'group-4'),
                            away_ml_raw: getGroupText(awayRow, 'group-1'),
                            home_ml_raw: getGroupText(homeRow, 'group-1'),
                        });
                    } catch (e) {}
                });

                return games;
            }''')

            print(f"  Found {len(game_data)} games")

            for i, game in enumerate(game_data):
                away_spread, away_spread_price = parse_spread_odds(game['away_spread_raw'])
                home_spread, home_spread_price = parse_spread_odds(game['home_spread_raw'])
                total, over_price = parse_total_odds(game['over_raw'])
                _, under_price = parse_total_odds(game['under_raw'])
                away_ml = parse_ml_odds(game['away_ml_raw'])
                home_ml = parse_ml_odds(game['home_ml_raw'])

                # Resolve team names
                if team_dict or canonical_games:
                    away, home = resolve_team_names(
                        game['away_team'], game['home_team'],
                        team_dict, canonical_games
                    )
                else:
                    away = normalize_team_name(game['away_team'], sport)
                    home = normalize_team_name(game['home_team'], sport)

                # Create game_id from teams
                game_id = f"{period_short}-{i}"

                # Market name matching Wagerzon/API pattern
                market_suffix = PERIOD_TO_MARKET_SUFFIX[period_short]
                market = f"spreads_{market_suffix}"

                all_odds.append({
                    "fetch_time": fetch_time,
                    "sport_key": config["sport_key"],
                    "game_id": game_id,
                    "game_date": game['game_date'],
                    "game_time": game['game_time'],
                    "away_team": away,
                    "home_team": home,
                    "market": market,
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

        browser.close()

    # Save to DuckDB
    if all_odds:
        save_to_database(sport, all_odds)

    print(f"\nScraped {len(all_odds)} total records for {sport.upper()}")
    return all_odds


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


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    scrape_hoop88(sport)
