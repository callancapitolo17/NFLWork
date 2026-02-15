#!/usr/bin/env python3
"""
BFA Gaming Odds Scraper for Answer Key Pipeline
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

BFA_URL = os.getenv("BFA_URL", "https://bfagaming.com")
BFA_USERNAME = os.getenv("BFA_USERNAME")
BFA_PASSWORD = os.getenv("BFA_PASSWORD")

DB_PATH = Path(__file__).parent / "bfa.duckdb"

# Sport configurations
SPORT_CONFIGS = {
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "nav_text": "NFL",
        "periods": ["1H", "1Q", "2Q", "3Q", "4Q"],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "nav_text": "NBA",
        "periods": ["1H", "1Q", "2Q", "3Q", "4Q"],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "nav_text": "NCAA(B)",
        "periods": ["1H"],
    },
}

# Map period short names to display names (for finding in UI)
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
    # Full names map to themselves
    "San Francisco 49ers": "San Francisco 49ers",
    "Chicago Bears": "Chicago Bears",
    "Cincinnati Bengals": "Cincinnati Bengals",
    "Buffalo Bills": "Buffalo Bills",
    "Denver Broncos": "Denver Broncos",
    "Cleveland Browns": "Cleveland Browns",
    "Tampa Bay Buccaneers": "Tampa Bay Buccaneers",
    "Arizona Cardinals": "Arizona Cardinals",
    "Los Angeles Chargers": "Los Angeles Chargers",
    "Kansas City Chiefs": "Kansas City Chiefs",
    "Indianapolis Colts": "Indianapolis Colts",
    "Washington Commanders": "Washington Commanders",
    "Dallas Cowboys": "Dallas Cowboys",
    "Miami Dolphins": "Miami Dolphins",
    "Philadelphia Eagles": "Philadelphia Eagles",
    "Atlanta Falcons": "Atlanta Falcons",
    "New York Giants": "New York Giants",
    "Jacksonville Jaguars": "Jacksonville Jaguars",
    "New York Jets": "New York Jets",
    "Detroit Lions": "Detroit Lions",
    "Green Bay Packers": "Green Bay Packers",
    "Carolina Panthers": "Carolina Panthers",
    "New England Patriots": "New England Patriots",
    "Las Vegas Raiders": "Las Vegas Raiders",
    "Los Angeles Rams": "Los Angeles Rams",
    "Baltimore Ravens": "Baltimore Ravens",
    "New Orleans Saints": "New Orleans Saints",
    "Seattle Seahawks": "Seattle Seahawks",
    "Pittsburgh Steelers": "Pittsburgh Steelers",
    "Houston Texans": "Houston Texans",
    "Tennessee Titans": "Tennessee Titans",
    "Minnesota Vikings": "Minnesota Vikings",
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


def parse_spread(text: str) -> float:
    """Extract spread value from text like '-3', '+4.5'."""
    if not text:
        return None
    text = text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
    text = text.replace('−', '-')
    match = re.search(r'([+-]?\d+\.?\d*)', text)
    return float(match.group(1)) if match else None


def parse_american_odds(text: str) -> int:
    """Extract American odds from text like '-110' or '+150' or 'EV'."""
    if not text:
        return None
    if text.upper() == 'EV':
        return 100  # Even money
    match = re.search(r'([+-]\d+)', text.replace('−', '-'))
    return int(match.group(1)) if match else None


def parse_total(text: str) -> float:
    """Extract total value from text like 'o22.5' or 'u45'."""
    if not text:
        return None
    text = text.replace('½', '.5').replace('¼', '.25').replace('¾', '.75')
    match = re.search(r'(\d+\.?\d*)', text)
    return float(match.group(1)) if match else None


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


def login_bfa(page) -> bool:
    """Login to BFA using Keycloak auth flow."""
    print("Checking login status...")

    must_login = page.locator('text=You must log in').count() > 0
    login_button = page.locator('button:has-text("Log In"), button:has-text("Log in")')

    if not must_login and login_button.count() == 0:
        print("Already logged in")
        return True

    print("Login required. Clicking login button...")
    try:
        header_login = page.locator('.header-auth-container button:has-text("Log In")')
        if header_login.count() > 0:
            header_login.click()
        else:
            login_button.first.click()

        print("Waiting for Keycloak login page...")
        page.wait_for_url('**/realms/**', timeout=15000)
        page.wait_for_timeout(2000)

        username_field = page.locator('#username')
        if username_field.count() > 0:
            username_field.fill(BFA_USERNAME)
        else:
            page.fill('input[name="username"]', BFA_USERNAME)

        page.wait_for_timeout(500)

        password_field = page.locator('#password')
        if password_field.count() > 0:
            password_field.fill(BFA_PASSWORD)
        else:
            page.fill('input[name="password"]', BFA_PASSWORD)

        page.wait_for_timeout(500)

        login_submit = page.locator('#kc-login')
        if login_submit.count() > 0:
            login_submit.click()
        else:
            page.click('input[type="submit"]')

        print("Waiting for redirect back to BFA...")
        page.wait_for_url('**bfagaming.com**', timeout=30000)
        page.wait_for_load_state('networkidle', timeout=30000)
        page.wait_for_timeout(3000)

        print("Login successful")
        return True

    except Exception as e:
        print(f"Login failed: {e}")
        return False


def extract_games_from_main_page(page) -> list:
    """Extract game info from the main sport page."""
    return page.evaluate('''() => {
        const games = [];
        const fixtures = document.querySelectorAll('.threeopt-fixture-component');

        fixtures.forEach((fixture, idx) => {
            try {
                const teamNames = fixture.querySelectorAll('.team-name-value');
                if (teamNames.length < 2) return;

                const awayTeam = teamNames[0].textContent.trim();
                const homeTeam = teamNames[1].textContent.trim();

                // Get "More wagers" link (single ID: /details/fixtureId) which has period tabs
                // NOT the Props link (two IDs: /details/eventId/fixtureId)
                const allLinks = fixture.closest('.fixture-view-component')?.querySelectorAll('a[href*="/details/"]');
                let detailsUrl = null;
                allLinks?.forEach(link => {
                    const href = link.getAttribute('href');
                    const text = link.textContent.toLowerCase();
                    // Main markets URL has single ID (e.g., /details/2648238)
                    // Props URL has two IDs (e.g., /details/2648233/2648238)
                    if (href && text.includes('more wagers')) {
                        detailsUrl = href;
                    }
                });
                // Fallback: get the shortest details URL (single ID)
                if (!detailsUrl && allLinks?.length > 0) {
                    let shortest = null;
                    allLinks.forEach(link => {
                        const href = link.getAttribute('href');
                        if (href && (!shortest || href.length < shortest.length)) {
                            shortest = href;
                        }
                    });
                    detailsUrl = shortest;
                }

                // Get game date/time
                const footer = fixture.closest('.fixture-view-component')?.querySelector('.fixture-live-info span');
                const gameDateTime = footer ? footer.textContent.trim() : '';

                games.push({
                    away_team: awayTeam,
                    home_team: homeTeam,
                    details_url: detailsUrl,
                    game_date_time: gameDateTime
                });
            } catch (e) {}
        });

        return games;
    }''')


def extract_period_odds(page) -> dict:
    """Extract odds from the current period tab.

    BFA structure: first 6 odd-buttons are the main lines:
    - Away spread + price
    - Away ML
    - Over total + price
    - Home spread + price
    - Home ML
    - Under total + price
    """
    return page.evaluate('''() => {
        const odds = {
            away_spread: null, away_spread_price: null,
            home_spread: null, home_spread_price: null,
            total: null, over_price: null, under_price: null,
            away_ml: null, home_ml: null
        };

        // Get all odd buttons
        const buttons = document.querySelectorAll('.odd-button');
        if (buttons.length < 6) return odds;

        function parseButton(button) {
            const spans = button.querySelectorAll('span');
            const values = [];
            spans.forEach(s => {
                const text = s.textContent.trim();
                if (text && !text.includes('<') && text.length < 20) {
                    values.push(text);
                }
            });
            return values;
        }

        // Button 0: Away spread (e.g., ["-3-110", "-3", "-110"])
        const awaySpreadVals = parseButton(buttons[0]);
        if (awaySpreadVals.length >= 2) {
            odds.away_spread = awaySpreadVals[1];  // "-3"
            odds.away_spread_price = awaySpreadVals[2] || awaySpreadVals[1].match(/[+-]\\d+$/)?.[0];  // "-110"
        }

        // Button 1: Away ML (e.g., ["-184", "-184"])
        const awayMLVals = parseButton(buttons[1]);
        if (awayMLVals.length >= 1) {
            odds.away_ml = awayMLVals[1] || awayMLVals[0];
        }

        // Button 2: Over total (e.g., ["o22.5-110", "o22.5", "-110"])
        const overVals = parseButton(buttons[2]);
        if (overVals.length >= 2) {
            odds.total = overVals[1];  // "o22.5"
            odds.over_price = overVals[2] || overVals[1].match(/[+-]\\d+$/)?.[0];  // "-110"
        }

        // Button 3: Home spread (e.g., ["+3-110", "+3", "-110"])
        const homeSpreadVals = parseButton(buttons[3]);
        if (homeSpreadVals.length >= 2) {
            odds.home_spread = homeSpreadVals[1];  // "+3"
            odds.home_spread_price = homeSpreadVals[2] || homeSpreadVals[1].match(/[+-]\\d+$/)?.[0];  // "-110"
        }

        // Button 4: Home ML (e.g., ["+155", "+155"])
        const homeMLVals = parseButton(buttons[4]);
        if (homeMLVals.length >= 1) {
            odds.home_ml = homeMLVals[1] || homeMLVals[0];
        }

        // Button 5: Under total (e.g., ["u22.5-110", "u22.5", "-110"])
        const underVals = parseButton(buttons[5]);
        if (underVals.length >= 2) {
            odds.under_price = underVals[2] || underVals[1].match(/[+-]\\d+$/)?.[0];  // "-110"
        }

        return odds;
    }''')


def scrape_bfa(sport: str, headless: bool = True):
    """Scrape BFA derivative odds and save to DuckDB."""
    if not BFA_USERNAME or not BFA_PASSWORD:
        raise ValueError("BFA_USERNAME and BFA_PASSWORD must be set in .env")

    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    fetch_time = datetime.now(timezone.utc)
    all_odds = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        print(f"Navigating to {BFA_URL}...")
        page.goto(BFA_URL, wait_until='networkidle', timeout=60000)
        page.wait_for_timeout(3000)

        if not login_bfa(page):
            print("Could not login to BFA")
            browser.close()
            return []

        # Navigate to sportsbook
        print("Navigating to sportsbook...")
        page.goto(f"{BFA_URL}/sports", wait_until='networkidle', timeout=60000)

        # Wait for content
        print("Waiting for sports content...")
        for i in range(10):
            page.wait_for_timeout(2000)
            content = page.content()
            if config['nav_text'] in content and 'threeopt-fixture-component' in content:
                print(f"  Sports content loaded after {(i+1)*2} seconds")
                break
            print(f"  Waiting... ({(i+1)*2}s)")

        # Navigate to sport
        print(f"Navigating to {config['nav_text']}...")
        sport_nav = page.locator(f'.main-nav-item:has-text("{config["nav_text"]}")')
        if sport_nav.count() > 0:
            sport_nav.first.click()
            page.wait_for_timeout(3000)
            print(f"  Clicked {config['nav_text']} in navigation")

        page.wait_for_load_state('networkidle', timeout=30000)
        page.wait_for_timeout(2000)

        # Extract games
        print("Extracting games from main page...")
        games = extract_games_from_main_page(page)
        print(f"Found {len(games)} games")

        for game in games[:3]:
            print(f"  {game['away_team']} @ {game['home_team']}")
            if game.get('details_url'):
                print(f"    Details: {game['details_url']}")

        # Load team name resolution data for non-NFL sports
        team_dict = load_team_dict(sport) if sport != "nfl" else {}
        canonical_games = load_canonical_games(sport) if sport != "nfl" else []

        # Process each game
        for i, game in enumerate(games):
            if team_dict or canonical_games:
                away_team, home_team = resolve_team_names(
                    game['away_team'], game['home_team'],
                    team_dict, canonical_games
                )
            else:
                away_team = normalize_team_name(game['away_team'], sport)
                home_team = normalize_team_name(game['home_team'], sport)

            print(f"\nProcessing game {i+1}/{len(games)}: {away_team} @ {home_team}")

            if not game.get('details_url'):
                print("  No details URL, skipping")
                continue

            # Navigate to game details
            details_url = f"{BFA_URL}{game['details_url']}"
            print(f"  Navigating to: {details_url}")

            try:
                page.goto(details_url, wait_until='networkidle', timeout=30000)
                page.wait_for_timeout(3000)  # Wait for Blazor content to load

                # Extract odds for each period
                for period_short in config['periods']:
                    period_display = PERIOD_DISPLAY[period_short]
                    period_standard = PERIOD_TO_STANDARD[period_short]
                    market_suffix = PERIOD_TO_MARKET_SUFFIX[period_short]

                    # Click the period tab
                    tab = page.locator(f'.scroller-item:has-text("{period_display}")')
                    if tab.count() == 0:
                        print(f"    {period_short}: Tab not found")
                        continue

                    tab.first.click()
                    page.wait_for_timeout(1500)

                    # Extract odds
                    raw_odds = extract_period_odds(page)

                    # Parse values
                    away_spread = parse_spread(raw_odds.get('away_spread'))
                    away_spread_price = parse_american_odds(raw_odds.get('away_spread_price'))
                    home_spread = parse_spread(raw_odds.get('home_spread'))
                    home_spread_price = parse_american_odds(raw_odds.get('home_spread_price'))
                    total = parse_total(raw_odds.get('total'))
                    over_price = parse_american_odds(raw_odds.get('over_price'))
                    under_price = parse_american_odds(raw_odds.get('under_price'))
                    away_ml = parse_american_odds(raw_odds.get('away_ml'))
                    home_ml = parse_american_odds(raw_odds.get('home_ml'))

                    if away_spread is not None or total is not None:
                        print(f"    {period_short}: Spread {away_spread}/{home_spread}, Total {total}, ML {away_ml}/{home_ml}")

                        all_odds.append({
                            "fetch_time": fetch_time,
                            "sport_key": config["sport_key"],
                            "game_id": f"bfa-{i}-{period_short}",
                            "game_date": game.get('game_date_time', '').split(',')[0].strip() if game.get('game_date_time') else '',
                            "game_time": game.get('game_date_time', '').split(',')[1].strip() if ',' in game.get('game_date_time', '') else '',
                            "away_team": away_team,
                            "home_team": home_team,
                            "market": f"spreads_{market_suffix}",
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
                    else:
                        print(f"    {period_short}: No odds found")

            except Exception as e:
                print(f"  Error: {e}")

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
    headless = "--headless" in sys.argv or len(sys.argv) <= 2
    scrape_bfa(sport, headless=headless)
