#!/usr/bin/env python3
"""
Wagerzon Odds Scraper v2 - Sport-Agnostic
Scrapes odds from Wagerzon and stores in format compatible with The Odds API
"""

import os
import re
import time
import duckdb
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from bs4 import BeautifulSoup
from dotenv import load_dotenv
from playwright.sync_api import sync_playwright

import sys

from config import (
    get_sport_url,
    get_sport_config,
    SKIP_SECTIONS,
    WAGERZON_BASE_URL,
)
from team_mapping import normalize_team_name

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Load environment
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")

DB_PATH = Path(__file__).parent / "wagerzon.duckdb"


def init_database(sport: str):
    """Initialize DuckDB with the odds table for a sport."""
    config = get_sport_config(sport)
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    # Create table matching The Odds API format
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
            -- Spread fields
            away_spread FLOAT,
            away_spread_price INTEGER,
            home_spread FLOAT,
            home_spread_price INTEGER,
            -- Total fields
            total FLOAT,
            over_price INTEGER,
            under_price INTEGER,
            -- Moneyline fields
            away_ml INTEGER,
            home_ml INTEGER
        )
    """)
    conn.close()


def login_to_wagerzon(page, sport_url: str):
    """Log into Wagerzon and navigate to sport page."""
    print(f"Navigating to {sport_url}...")
    page.goto(sport_url)
    page.wait_for_load_state('networkidle')

    # Check if we need to login
    if "NewSchedule" not in page.url:
        print("Login required...")

        # Fill credentials
        account_field = page.locator('input#Account')
        account_field.fill(WAGERZON_USERNAME)

        password_field = page.locator('input#Password')
        password_field.fill(WAGERZON_PASSWORD)

        # Submit
        submit_btn = page.locator('button[name="BtnSubmit"]')
        submit_btn.click()

        page.wait_for_load_state('networkidle')
        time.sleep(2)

        # Navigate to sport page after login
        if "NewSchedule" not in page.url:
            page.goto(sport_url)
            page.wait_for_load_state('networkidle')

    print(f"Current URL: {page.url}")
    return page


def parse_spread(text: str) -> tuple[Optional[float], Optional[int]]:
    """Parse spread text like '+1½-110' into (spread, odds)."""
    if not text:
        return None, None
    text = text.strip().replace('½', '.5')

    # Match patterns like "+1.5-110", "-3-110", "PK-110", "PK+100", "-.5-110", "+.5-110"
    # The spread can be: PK, or [+-]?[digits]?[.digits]? (allowing -.5 without leading 0)
    match = re.match(r'(PK|[+-]?\d*\.?\d+)\s*([+-]\d+)', text)
    if match:
        spread_str = match.group(1)
        odds_str = match.group(2)
        spread = 0.0 if spread_str == 'PK' else float(spread_str)
        odds = int(odds_str)
        return spread, odds

    return None, None


def parse_total(text: str) -> tuple[Optional[float], Optional[int]]:
    """Parse total text like 'o137-110' into (total, odds)."""
    if not text:
        return None, None
    text = text.strip().replace('½', '.5')

    # Match patterns like "o137-110", "u137-110", "o7½+101"
    match = re.match(r'[ou](\d+\.?\d*)\s*([+-]?\d+|EV)', text)
    if match:
        total = float(match.group(1))
        odds_str = match.group(2)
        odds = 100 if odds_str == 'EV' else int(odds_str)
        return total, odds

    return None, None


def parse_moneyline(text: str) -> Optional[int]:
    """Parse moneyline text like '+105' or '-125'."""
    if not text:
        return None
    text = text.strip()

    if text == 'EV':
        return 100

    try:
        # Handle both "+150" and "150" formats
        return int(text.replace('+', ''))
    except ValueError:
        return None


def determine_period_from_rotation(rotation: str, config: dict) -> str:
    """Determine period from rotation number prefix."""
    if not rotation:
        return "fg"

    # 3-digit rotations are full game
    if len(rotation) <= 3:
        return "fg"

    # 4+ digit rotations: first digit indicates period
    prefix = rotation[0]
    return config.get("rotation_periods", {}).get(prefix, "fg")


def clean_team_name(team: str, config: dict) -> str:
    """Remove period prefixes from team names."""
    for prefix in config.get("team_prefix_patterns", []):
        if team.startswith(prefix):
            team = team[len(prefix):]
    return team.strip()


def should_skip_section(section_header: str) -> bool:
    """Check if a section should be skipped (props, futures, etc.)."""
    upper = section_header.upper()
    return any(skip in upper for skip in SKIP_SECTIONS)


def determine_market_from_section(section_header: str, config: dict) -> Optional[str]:
    """Determine market type from section header."""
    upper = section_header.upper()

    for pattern, market in config.get("section_markets", {}).items():
        if pattern in upper:
            return market

    return None


def scrape_odds(page, sport: str) -> list[dict]:
    """Scrape odds from the current page."""
    config = get_sport_config(sport)

    time.sleep(3)  # Wait for dynamic content

    html = page.content()
    soup = BeautifulSoup(html, 'html.parser')

    # Save for debugging
    debug_path = Path(__file__).parent / f"debug_{sport}.html"
    with open(debug_path, 'w') as f:
        f.write(html)

    # Load team name resolution data for non-NFL sports
    team_dict = load_team_dict(sport) if sport != "nfl" else {}
    canonical_games = load_canonical_games(sport) if sport != "nfl" else []

    odds_data = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    sport_key = config["sport_key"]

    competitions = soup.find_all('div', class_='Competition')
    print(f"Found {len(competitions)} competitions")

    for comp in competitions:
        # Find section header for this competition
        section_h1 = comp.find_previous('h1')
        section_header = section_h1.get_text().strip() if section_h1 else ""

        # Skip props and futures
        if should_skip_section(section_header):
            continue

        # Determine base market from section
        base_market = determine_market_from_section(section_header, config)
        if not base_market:
            continue

        game_rows = comp.find_all('div', class_='GameRow')
        if len(game_rows) < 2:
            continue

        away_row = game_rows[0]
        home_row = game_rows[1]

        # Get team names
        away_team_span = away_row.find('span', class_='Team')
        home_team_span = home_row.find('span', class_='Team')

        if not away_team_span or not home_team_span:
            continue

        away_team_raw = away_team_span.get_text().strip()
        home_team_raw = home_team_span.get_text().strip()

        # Clean team names (remove period prefixes) and normalize to standard names
        away_team_cleaned = clean_team_name(away_team_raw, config)
        home_team_cleaned = clean_team_name(home_team_raw, config)
        if team_dict or canonical_games:
            away_team, home_team = resolve_team_names(
                away_team_cleaned, home_team_cleaned,
                team_dict, canonical_games
            )
        else:
            away_team = normalize_team_name(away_team_cleaned, sport)
            home_team = normalize_team_name(home_team_cleaned, sport)

        # Get rotation numbers
        away_cols = away_row.find_all('div', class_='col-xs-1')
        home_cols = home_row.find_all('div', class_='col-xs-1')

        away_rot = ""
        game_date = ""
        game_time = ""

        for col in away_cols:
            text = col.get_text().strip()
            if text.isdigit() and len(text) >= 3:
                away_rot = text
            elif text and not game_date:
                game_date = text

        for col in home_cols:
            text = col.get_text().strip()
            if not text.isdigit() and text and not game_time:
                game_time = text

        # Determine period from rotation number
        period = determine_period_from_rotation(away_rot, config)

        # Build market name with period
        if base_market == "spreads_quarters":
            market = f"spreads_{period}"
        elif base_market.endswith("_h1") or base_market.endswith("_q1"):
            market = base_market
        else:
            market = base_market if period == "fg" else f"{base_market}_{period}"

        # Parse odds from buttons
        away_btns = away_row.find_all('div', class_='btn-odds')
        home_btns = home_row.find_all('div', class_='btn-odds')

        # Away team odds
        away_spread, away_spread_price = None, None
        total, over_price = None, None
        away_ml = None

        if len(away_btns) >= 1:
            away_spread, away_spread_price = parse_spread(away_btns[0].get_text())
        if len(away_btns) >= 2:
            total, over_price = parse_total(away_btns[1].get_text())
        if len(away_btns) >= 3:
            away_ml = parse_moneyline(away_btns[2].get_text())

        # Home team odds
        home_spread, home_spread_price = None, None
        _, under_price = None, None
        home_ml = None

        if len(home_btns) >= 1:
            home_spread, home_spread_price = parse_spread(home_btns[0].get_text())
        if len(home_btns) >= 2:
            _, under_price = parse_total(home_btns[1].get_text())
        if len(home_btns) >= 3:
            home_ml = parse_moneyline(home_btns[2].get_text())

        # Create game ID from rotation
        home_rot = ""
        for col in home_cols:
            text = col.get_text().strip()
            if text.isdigit() and len(text) >= 3:
                home_rot = text
        game_id = f"{away_rot}-{home_rot}"

        record = {
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

        odds_data.append(record)
        print(f"  {market}: {away_team} @ {home_team} | {away_spread}/{home_spread} | {total}")

    return odds_data


def save_odds(odds_data: list[dict], sport: str):
    """Save odds to DuckDB."""
    if not odds_data:
        print("No odds data to save")
        return

    config = get_sport_config(sport)
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    # Insert records
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
    print(f"\nDatabase now has {result[0]} total records in {table_name}")

    conn.close()


def scrape_wagerzon(sport: str, headless: bool = True):
    """Main function to scrape odds for a sport."""
    if not WAGERZON_USERNAME or not WAGERZON_PASSWORD:
        raise ValueError("WAGERZON_USERNAME and WAGERZON_PASSWORD must be set in bet_logger/.env")

    sport_url = get_sport_url(sport)
    init_database(sport)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context()
        page = context.new_page()

        try:
            login_to_wagerzon(page, sport_url)
            time.sleep(3)

            odds_data = scrape_odds(page, sport)

            if odds_data:
                save_odds(odds_data, sport)
                print(f"\nSaved {len(odds_data)} records for {sport}")
            else:
                print("No odds data scraped")

            return odds_data

        finally:
            browser.close()


if __name__ == "__main__":
    import sys

    # Default to NFL, can pass sport as argument
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    headless = "--no-headless" not in sys.argv

    print(f"Starting Wagerzon {sport.upper()} odds scraper...")
    scrape_wagerzon(sport, headless=headless)
