#!/usr/bin/env python3
"""
Wagerzon College Basketball Odds Scraper
Scrapes CBB spreads and totals (full game + 1H) and stores in DuckDB
"""

import os
import re
import time
import duckdb
from datetime import datetime, timezone
from pathlib import Path
from bs4 import BeautifulSoup
from dotenv import load_dotenv
from playwright.sync_api import sync_playwright

# Load environment from bet_logger's .env
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

# Configuration
WAGERZON_URL = "https://backend.wagerzon.com"
CBB_URL = "https://backend.wagerzon.com/wager/NewSchedule.aspx?WT=0&lg=43,403,45"
WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")

DB_PATH = Path(__file__).parent / "wagerzon_cbb.duckdb"


def init_database():
    """Initialize DuckDB with the odds table."""
    conn = duckdb.connect(str(DB_PATH))
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cbb_odds (
            fetch_time TIMESTAMP,
            game_date VARCHAR,
            game_time VARCHAR,
            game_id VARCHAR,
            away_team VARCHAR,
            home_team VARCHAR,
            period VARCHAR,
            away_spread FLOAT,
            away_spread_odds INTEGER,
            home_spread FLOAT,
            home_spread_odds INTEGER,
            total FLOAT,
            over_odds INTEGER,
            under_odds INTEGER,
            away_ml INTEGER,
            home_ml INTEGER
        )
    """)
    conn.close()


def login_to_wagerzon(page):
    """Log into Wagerzon."""
    print(f"Navigating to {WAGERZON_URL}...")
    page.goto(WAGERZON_URL)
    page.wait_for_load_state('networkidle')

    # Try to navigate directly to CBB page
    print("Navigating to CBB odds page...")
    page.goto(CBB_URL)
    page.wait_for_load_state('networkidle')

    current_url = page.url

    # If redirected to login page, log in
    if "Login" in current_url or "NewSchedule" not in current_url:
        print("Login required...")
        try:
            username_field = page.locator('input[type="text"], input[name="username"], input[name="txtUsername"]').first
            username_field.fill(WAGERZON_USERNAME, timeout=10000)

            password_field = page.locator('input[type="password"]').first
            password_field.fill(WAGERZON_PASSWORD)

            login_button = page.locator('input[type="submit"], button[type="submit"]').first
            login_button.click()

            page.wait_for_load_state('networkidle')
            print("Login submitted!")

            time.sleep(2)

            if "NewSchedule" not in page.url:
                page.goto(CBB_URL)
                page.wait_for_load_state('networkidle')

        except Exception as e:
            print(f"Auto-login failed: {e}")
            raise

    print(f"Current URL: {page.url}")
    return page


def parse_spread_odds(text: str) -> tuple:
    """Parse spread text like '+1½-110' into (spread, odds)."""
    if not text:
        return None, None
    text = text.strip()

    # Replace ½ with .5
    text = text.replace('½', '.5')

    # Match patterns like "+1.5-110", "-3-110", "PK-110"
    match = re.match(r'(PK|[+-]?\d+\.?\d*)\s*([+-]\d+)', text)
    if match:
        spread_str = match.group(1)
        odds_str = match.group(2)

        spread = 0.0 if spread_str == 'PK' else float(spread_str)
        odds = int(odds_str)
        return spread, odds

    return None, None


def parse_total_odds(text: str) -> tuple:
    """Parse total text like 'o137-110' into (total, odds)."""
    if not text:
        return None, None
    text = text.strip()

    # Replace ½ with .5
    text = text.replace('½', '.5')

    # Match patterns like "o137-110", "u137-110"
    match = re.match(r'[ou](\d+\.?\d*)\s*([+-]\d+)', text)
    if match:
        total = float(match.group(1))
        odds = int(match.group(2))
        return total, odds

    return None, None


def parse_moneyline(text: str) -> int:
    """Parse moneyline text like '+105' or '-125'."""
    if not text:
        return None
    text = text.strip()

    try:
        return int(text.replace('+', ''))
    except ValueError:
        return None


def scrape_odds_from_page(page) -> list:
    """Scrape odds from the current page."""
    time.sleep(3)  # Wait for dynamic content

    html_content = page.content()
    soup = BeautifulSoup(html_content, 'html.parser')

    # Save HTML for debugging
    debug_path = Path(__file__).parent / "debug_cbb.html"
    with open(debug_path, 'w') as f:
        f.write(html_content)

    odds_data = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    # Find all Competition containers (each contains a matchup)
    competitions = soup.find_all('div', class_='Competition')
    print(f"Found {len(competitions)} competitions")

    current_period = 'FG'  # Default to full game

    # Track section headers to determine period
    for h1 in soup.find_all('h1'):
        h1_text = h1.get_text().strip().upper()
        if '1H COLLEGE BASKETBALL' in h1_text:
            # Mark where 1H section starts
            pass

    # Process each competition
    for comp in competitions:
        # Check if this is in 1H section by looking at team names
        game_rows = comp.find_all('div', class_='GameRow')
        if len(game_rows) < 2:
            continue  # Need at least 2 rows (away + home)

        # Get the first team name to check if it's 1H
        team_span = game_rows[0].find('span', class_='Team')
        if not team_span:
            continue

        first_team = team_span.get_text().strip()
        period = '1H' if first_team.startswith('1H ') else 'FG'

        # Skip team totals
        if 'TEAM TOTAL' in first_team.upper():
            continue

        # Parse away team row
        away_row = game_rows[0]
        home_row = game_rows[1]

        # Get date/time
        date_col = away_row.find('div', class_='col-xs-1')
        game_date = date_col.get_text().strip() if date_col else ''

        time_col = home_row.find('div', class_='col-xs-1')
        game_time = time_col.get_text().strip() if time_col else ''

        # Get rotation numbers
        rot_cols = away_row.find_all('div', class_='col-xs-1')
        away_rot = rot_cols[1].get_text().strip() if len(rot_cols) > 1 else ''

        rot_cols_home = home_row.find_all('div', class_='col-xs-1')
        home_rot = rot_cols_home[1].get_text().strip() if len(rot_cols_home) > 1 else ''

        # Get team names
        away_team_span = away_row.find('span', class_='Team')
        home_team_span = home_row.find('span', class_='Team')

        away_team = away_team_span.get_text().strip() if away_team_span else ''
        home_team = home_team_span.get_text().strip() if home_team_span else ''

        # Remove "1H " prefix if present
        away_team = away_team.replace('1H ', '').strip()
        home_team = home_team.replace('1H ', '').strip()

        # Get odds buttons
        away_odds_btns = away_row.find_all('div', class_='btn-odds')
        home_odds_btns = home_row.find_all('div', class_='btn-odds')

        # Parse spread, total, moneyline for away team
        away_spread, away_spread_odds = None, None
        total, over_odds = None, None
        away_ml = None

        if len(away_odds_btns) >= 1:
            away_spread, away_spread_odds = parse_spread_odds(away_odds_btns[0].get_text())
        if len(away_odds_btns) >= 2:
            total, over_odds = parse_total_odds(away_odds_btns[1].get_text())
        if len(away_odds_btns) >= 3:
            away_ml = parse_moneyline(away_odds_btns[2].get_text())

        # Parse for home team
        home_spread, home_spread_odds = None, None
        _, under_odds = None, None
        home_ml = None

        if len(home_odds_btns) >= 1:
            home_spread, home_spread_odds = parse_spread_odds(home_odds_btns[0].get_text())
        if len(home_odds_btns) >= 2:
            _, under_odds = parse_total_odds(home_odds_btns[1].get_text())
        if len(home_odds_btns) >= 3:
            home_ml = parse_moneyline(home_odds_btns[2].get_text())

        # Create game ID from rotation numbers
        game_id = f"{away_rot}-{home_rot}"

        odds_record = {
            'fetch_time': fetch_time,
            'game_date': game_date,
            'game_time': game_time,
            'game_id': game_id,
            'away_team': away_team,
            'home_team': home_team,
            'period': period,
            'away_spread': away_spread,
            'away_spread_odds': away_spread_odds,
            'home_spread': home_spread,
            'home_spread_odds': home_spread_odds,
            'total': total,
            'over_odds': over_odds,
            'under_odds': under_odds,
            'away_ml': away_ml,
            'home_ml': home_ml
        }

        odds_data.append(odds_record)
        print(f"  {period}: {away_team} @ {home_team} | Spread: {away_spread}/{home_spread} | Total: {total}")

    return odds_data


def save_odds_to_db(odds_data: list):
    """Save odds to DuckDB."""
    if not odds_data:
        print("No odds data to save")
        return

    conn = duckdb.connect(str(DB_PATH))

    conn.executemany("""
        INSERT INTO cbb_odds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, [
        (
            d["fetch_time"], d["game_date"], d["game_time"], d["game_id"],
            d["away_team"], d["home_team"], d["period"],
            d["away_spread"], d["away_spread_odds"],
            d["home_spread"], d["home_spread_odds"],
            d["total"], d["over_odds"], d["under_odds"],
            d["away_ml"], d["home_ml"]
        )
        for d in odds_data
    ])

    result = conn.execute("SELECT COUNT(*) FROM cbb_odds").fetchone()
    print(f"\nDatabase now has {result[0]} total records")

    conn.close()


def scrape_wagerzon_cbb(headless: bool = True):
    """Main function to scrape CBB odds."""
    if not WAGERZON_USERNAME or not WAGERZON_PASSWORD:
        raise ValueError("WAGERZON_USERNAME and WAGERZON_PASSWORD must be set in bet_logger/.env")

    init_database()

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context()
        page = context.new_page()

        try:
            login_to_wagerzon(page)
            time.sleep(3)

            odds_data = scrape_odds_from_page(page)

            if odds_data:
                save_odds_to_db(odds_data)
                print(f"\nSaved {len(odds_data)} games")
            else:
                print("No odds data scraped")

        finally:
            browser.close()

    return odds_data


if __name__ == "__main__":
    import sys

    headless = "--headless" in sys.argv or "-h" not in sys.argv

    print("Starting Wagerzon CBB odds scraper...")
    scrape_wagerzon_cbb(headless=True)
