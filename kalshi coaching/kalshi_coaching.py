#!/usr/bin/env python3
"""
Kalshi NFL Coaching Odds Fetcher
Fetches current odds for NFL head coaching hire markets and stores in DuckDB.
Run with: ./venv/bin/python kalshi_coaching.py
"""

import urllib.request
import json
import duckdb
import time
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"

# NFL coaching series tickers
NFL_COACHING_SERIES = [
    "KXNEWCOACHDAL",      # Dallas Cowboys
    "KXNEWCOACHJAX",      # Jacksonville Jaguars
    "KXNEWCOACHNO",       # New Orleans Saints
    "KXNEWCOACHLV",       # Las Vegas Raiders
    "KXNEWCOACHCHIBEARS", # Chicago Bears
    "KXNEWCOACHNYJ",      # New York Jets
    "KXNEWCOACHNE",       # New England Patriots
    "KXTENNCOACH",        # Tennessee Titans
    "KXNEXTNFLCOACH",     # General next NFL coach markets
    "KXNFLHIRECOACH",     # NFL hire coach
    "KXNYGCOACH",         # NY Giants
    "KXATLCOACH",         # Atlanta Falcons
]

# Map series tickers to team names
SERIES_TO_TEAM = {
    "KXNEWCOACHDAL": "Dallas Cowboys",
    "KXNEWCOACHJAX": "Jacksonville Jaguars",
    "KXNEWCOACHNO": "New Orleans Saints",
    "KXNEWCOACHLV": "Las Vegas Raiders",
    "KXNEWCOACHCHIBEARS": "Chicago Bears",
    "KXNEWCOACHNYJ": "New York Jets",
    "KXNEWCOACHNE": "New England Patriots",
    "KXTENNCOACH": "Tennessee Titans",
    "KXNYGCOACH": "New York Giants",
    "KXATLCOACH": "Atlanta Falcons",
}


def fetch_markets_for_series(series_ticker: str) -> list:
    """Fetch all open markets for a given series ticker."""
    markets = []
    cursor = None

    while True:
        params = f"series_ticker={series_ticker}&status=open&limit=100"
        if cursor:
            params += f"&cursor={cursor}"

        url = f"{BASE_URL}/markets?{params}"

        try:
            with urllib.request.urlopen(url) as response:
                if response.status != 200:
                    print(f"Error fetching {series_ticker}: {response.status}")
                    break
                data = json.loads(response.read().decode())
        except urllib.error.HTTPError as e:
            print(f"Error fetching {series_ticker}: {e.code}")
            break
        except Exception as e:
            print(f"Error fetching {series_ticker}: {e}")
            break

        markets.extend(data.get("markets", []))

        cursor = data.get("cursor")
        if not cursor or len(data.get("markets", [])) == 0:
            break

    return markets


def extract_team_from_title(title: str) -> str:
    """Extract team name from market title."""
    title_lower = title.lower()

    team_keywords = {
        "dallas": "Dallas Cowboys",
        "cowboys": "Dallas Cowboys",
        "jacksonville": "Jacksonville Jaguars",
        "jaguars": "Jacksonville Jaguars",
        "new orleans": "New Orleans Saints",
        "saints": "New Orleans Saints",
        "las vegas": "Las Vegas Raiders",
        "raiders": "Las Vegas Raiders",
        "chicago": "Chicago Bears",
        "bears": "Chicago Bears",
        "jets": "New York Jets",
        "new england": "New England Patriots",
        "patriots": "New England Patriots",
        "tennessee": "Tennessee Titans",
        "titans": "Tennessee Titans",
        "giants": "New York Giants",
        "atlanta": "Atlanta Falcons",
        "falcons": "Atlanta Falcons",
        "miami": "Miami Dolphins",
        "dolphins": "Miami Dolphins",
        "arizona": "Arizona Cardinals",
        "cardinals": "Arizona Cardinals",
        "pittsburgh": "Pittsburgh Steelers",
        "steelers": "Pittsburgh Steelers",
        "baltimore": "Baltimore Ravens",
        "ravens": "Baltimore Ravens",
    }

    for keyword, team in team_keywords.items():
        if keyword in title_lower:
            return team
    return "Unknown"


def fetch_all_coaching_odds() -> list:
    """Fetch all NFL coaching odds from Kalshi."""
    all_odds = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    for i, series in enumerate(NFL_COACHING_SERIES):
        if i > 0:
            time.sleep(0.1)  # Small delay to avoid rate limiting
        print(f"Fetching {series}...")
        markets = fetch_markets_for_series(series)

        for market in markets:
            # Extract candidate name
            candidate = (
                market.get("yes_sub_title") or
                market.get("custom_strike", {}).get("Person") or
                market.get("no_sub_title") or
                "Unknown"
            )

            # Skip if no real candidate
            if candidate in ["Unknown", "", None]:
                continue

            # Extract team from series or title
            team = SERIES_TO_TEAM.get(series)
            if not team:
                team = extract_team_from_title(market.get("title", ""))

            # Skip non-NFL coaching markets that might have slipped in
            if team == "Unknown":
                continue

            odds_data = {
                "fetch_time": fetch_time,
                "team": team,
                "candidate": candidate,
                "ticker": market.get("ticker", ""),
                "yes_bid": market.get("yes_bid", 0),
                "yes_ask": market.get("yes_ask", 0),
                "last_price": market.get("last_price", 0),
                "volume": market.get("volume", 0),
                "volume_24h": market.get("volume_24h", 0),
                "liquidity": market.get("liquidity", 0),
                "open_interest": market.get("open_interest", 0),
                "series_ticker": series,
                "event_ticker": market.get("event_ticker", ""),
            }
            all_odds.append(odds_data)

    return all_odds


def save_to_duckdb(odds_data: list, db_path: str = None):
    """Save odds data to DuckDB, appending to existing table."""
    if db_path is None:
        db_path = Path(__file__).parent / "kalshi_coaching.duckdb"

    conn = duckdb.connect(str(db_path))

    # Create table if not exists
    conn.execute("""
        CREATE TABLE IF NOT EXISTS coaching_odds (
            fetch_time TIMESTAMP,
            team VARCHAR,
            candidate VARCHAR,
            ticker VARCHAR,
            yes_bid INTEGER,
            yes_ask INTEGER,
            last_price INTEGER,
            volume BIGINT,
            volume_24h BIGINT,
            liquidity BIGINT,
            open_interest INTEGER,
            series_ticker VARCHAR,
            event_ticker VARCHAR
        )
    """)

    # Insert data
    if odds_data:
        conn.executemany("""
            INSERT INTO coaching_odds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            (
                d["fetch_time"], d["team"], d["candidate"], d["ticker"],
                d["yes_bid"], d["yes_ask"], d["last_price"], d["volume"],
                d["volume_24h"], d["liquidity"], d["open_interest"],
                d["series_ticker"], d["event_ticker"]
            )
            for d in odds_data
        ])

    # Show table stats
    result = conn.execute("SELECT COUNT(*) as total, COUNT(DISTINCT fetch_time) as snapshots FROM coaching_odds").fetchone()
    print(f"\nDatabase: {result[0]:,} total records across {result[1]} snapshots")

    conn.close()
    print(f"Saved {len(odds_data)} records to {db_path}")
    return db_path


def main():
    """Main function to fetch and save coaching odds."""
    print("Fetching NFL coaching odds from Kalshi...")
    odds = fetch_all_coaching_odds()

    print(f"\nFound {len(odds)} candidate markets")

    # Print summary
    teams = {}
    for o in odds:
        team = o["team"]
        if team not in teams:
            teams[team] = []
        teams[team].append(o)

    print("\n" + "="*60)
    for team, candidates in sorted(teams.items()):
        print(f"\n{team}:")
        # Sort by last_price descending
        sorted_candidates = sorted(candidates, key=lambda x: x["last_price"], reverse=True)
        for c in sorted_candidates[:5]:  # Top 5
            prob = c["last_price"]
            print(f"  {c['candidate']}: {prob}% (vol: {c['volume']:,})")

    # Save to DuckDB
    db_path = save_to_duckdb(odds)

    return odds, db_path


if __name__ == "__main__":
    main()
