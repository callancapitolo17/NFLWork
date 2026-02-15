#!/usr/bin/env python3
"""
Kalshi NFL Coaching Odds Fetcher
Fetches current odds for NFL head coaching hire markets and stores in DuckDB.
Also fetches positions and resting orders when API credentials are configured.
Run with: ./venv/bin/python kalshi_coaching.py
"""

import urllib.request
import json
import duckdb
import time
import os
import base64
import hashlib
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"

# Load credentials from .env if available
def load_credentials():
    """Load Kalshi API credentials from .env file."""
    env_path = Path(__file__).parent / ".env"
    creds = {"api_key_id": None, "private_key_path": None}

    if env_path.exists():
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("KALSHI_API_KEY_ID="):
                    creds["api_key_id"] = line.split("=", 1)[1].strip()
                elif line.startswith("KALSHI_PRIVATE_KEY_PATH="):
                    creds["private_key_path"] = line.split("=", 1)[1].strip()

    return creds


def sign_request(private_key_path: str, timestamp_str: str, method: str, path: str) -> str:
    """Sign a request using RSA-PSS with SHA-256 (Kalshi's required format)."""
    try:
        from cryptography.hazmat.primitives import hashes, serialization
        from cryptography.hazmat.primitives.asymmetric import padding
    except ImportError:
        print("  cryptography package required for authenticated requests")
        print("  Install with: pip install cryptography")
        return None

    with open(private_key_path, "rb") as f:
        private_key = serialization.load_pem_private_key(f.read(), password=None)

    # Kalshi signature format: timestamp + method + path (without query params)
    # Strip query params if present
    path_without_query = path.split("?")[0]
    message = f"{timestamp_str}{method}{path_without_query}".encode()
    signature = private_key.sign(
        message,
        padding.PSS(
            mgf=padding.MGF1(hashes.SHA256()),
            salt_length=padding.PSS.DIGEST_LENGTH  # Kalshi uses DIGEST_LENGTH, not MAX_LENGTH
        ),
        hashes.SHA256()
    )
    return base64.b64encode(signature).decode()


def authenticated_request(method: str, path: str, api_key_id: str, private_key_path: str):
    """Make an authenticated request to Kalshi API."""
    timestamp = datetime.now(timezone.utc)
    timestamp_str = str(int(timestamp.timestamp() * 1000))

    # Full path for signing must include /trade-api/v2 prefix
    full_path = f"/trade-api/v2{path}"

    signature = sign_request(private_key_path, timestamp_str, method, full_path)
    if not signature:
        return None

    url = f"{BASE_URL}{path}"
    req = urllib.request.Request(url, method=method)
    req.add_header("KALSHI-ACCESS-KEY", api_key_id)
    req.add_header("KALSHI-ACCESS-SIGNATURE", signature)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", timestamp_str)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        print(f"  Auth request failed: {e.code} - {e.read().decode()}")
        return None
    except Exception as e:
        print(f"  Auth request error: {e}")
        return None


def fetch_positions(api_key_id: str, private_key_path: str) -> list:
    """Fetch all open positions from Kalshi portfolio."""
    positions = []
    cursor = None

    while True:
        path = "/portfolio/positions?limit=100"
        if cursor:
            path += f"&cursor={cursor}"

        data = authenticated_request("GET", path, api_key_id, private_key_path)
        if not data:
            break

        positions.extend(data.get("market_positions", []))
        cursor = data.get("cursor")
        if not cursor:
            break

    return positions


def fetch_orders(api_key_id: str, private_key_path: str) -> list:
    """Fetch all resting orders from Kalshi portfolio."""
    orders = []
    cursor = None

    while True:
        path = "/portfolio/orders?status=resting&limit=100"
        if cursor:
            path += f"&cursor={cursor}"

        data = authenticated_request("GET", path, api_key_id, private_key_path)
        if not data:
            break

        orders.extend(data.get("orders", []))
        cursor = data.get("cursor")
        if not cursor:
            break

    return orders


def fetch_market_details(ticker: str) -> dict:
    """Fetch market details including expiration time."""
    url = f"{BASE_URL}/markets/{ticker}"
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            return data.get("market", {})
    except Exception:
        return {}

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
        "buffalo": "Buffalo Bills",
        "bills": "Buffalo Bills",
        "cleveland": "Cleveland Browns",
        "browns": "Cleveland Browns",
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
                "no_bid": market.get("no_bid", 0),
                "no_ask": market.get("no_ask", 0),
                "last_price": market.get("last_price", 0),
                "volume": market.get("volume", 0),
                "volume_24h": market.get("volume_24h", 0),
                "liquidity": int(float(market.get("liquidity_dollars") or 0)),
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

    # Create table if not exists (v2 with no_bid/no_ask)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS coaching_odds_v2 (
            fetch_time TIMESTAMP,
            team VARCHAR,
            candidate VARCHAR,
            ticker VARCHAR,
            yes_bid INTEGER,
            yes_ask INTEGER,
            no_bid INTEGER,
            no_ask INTEGER,
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
            INSERT INTO coaching_odds_v2 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            (
                d["fetch_time"], d["team"], d["candidate"], d["ticker"],
                d["yes_bid"], d["yes_ask"], d["no_bid"], d["no_ask"],
                d["last_price"], d["volume"], d["volume_24h"], d["liquidity"],
                d["open_interest"], d["series_ticker"], d["event_ticker"]
            )
            for d in odds_data
        ])

    # Show table stats
    result = conn.execute("SELECT COUNT(*) as total, COUNT(DISTINCT fetch_time) as snapshots FROM coaching_odds_v2").fetchone()
    print(f"\nDatabase: {result[0]:,} total records across {result[1]} snapshots")

    conn.close()
    print(f"Saved {len(odds_data)} records to {db_path}")
    return db_path


def save_positions_to_duckdb(positions: list, db_path: str = None):
    """Save positions to DuckDB and return changes from last snapshot."""
    if db_path is None:
        db_path = Path(__file__).parent / "kalshi_coaching.duckdb"

    conn = duckdb.connect(str(db_path))
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    # Create positions table if not exists
    conn.execute("""
        CREATE TABLE IF NOT EXISTS positions (
            fetch_time TIMESTAMP,
            ticker VARCHAR,
            position INTEGER,
            market_exposure DOUBLE,
            realized_pnl DOUBLE,
            resting_orders_count INTEGER,
            total_traded DOUBLE,
            fees_paid DOUBLE
        )
    """)

    # Get previous positions for comparison
    try:
        prev_positions = conn.execute("""
            SELECT ticker, position, market_exposure
            FROM positions
            WHERE fetch_time = (SELECT MAX(fetch_time) FROM positions)
        """).fetchall()

        # Check if previous values are in cents (sum of exposure > 50000 suggests cents)
        total_exp = sum(abs(float(row[2])) for row in prev_positions)
        is_cents = total_exp > 50000  # $500+ exposure = 50000 cents minimum

        prev_map = {}
        for row in prev_positions:
            exp = float(row[2])
            if is_cents:
                exp = exp / 100.0  # Convert cents to dollars
            prev_map[row[0]] = {"position": row[1], "exposure": exp}
    except Exception:
        prev_map = {}

    # Calculate changes (including closed positions)
    changes = []
    current_tickers = set()

    for p in positions:
        ticker = p.get("ticker", "")
        current_tickers.add(ticker)
        current_pos = p.get("position", 0)
        # Use _dollars fields from API (convert to float, fall back to cents/100)
        raw_exp = p.get("market_exposure_dollars")
        if raw_exp is not None:
            current_exp = float(raw_exp)
        else:
            current_exp = p.get("market_exposure", 0) / 100.0

        prev = prev_map.get(ticker, {"position": 0, "exposure": 0.0})
        pos_change = current_pos - prev["position"]
        exp_change = current_exp - float(prev["exposure"])

        if pos_change != 0 or abs(exp_change) > 0.001:
            changes.append({
                "ticker": ticker,
                "pos_change": pos_change,
                "exp_change": exp_change,
                "current_pos": current_pos,
                "current_exp": current_exp
            })

    # Check for closed positions (existed before, not in current)
    for ticker, prev in prev_map.items():
        if ticker not in current_tickers and prev["position"] != 0:
            changes.append({
                "ticker": ticker,
                "pos_change": -prev["position"],
                "exp_change": -float(prev["exposure"]),
                "current_pos": 0,
                "current_exp": 0.0,
                "closed": True
            })

    # Helper to convert API dollar values (may be string or None)
    def to_dollars(p, dollar_key, cents_key):
        val = p.get(dollar_key)
        if val is not None:
            return float(val)
        cents = p.get(cents_key, 0)
        return float(cents) / 100.0 if cents else 0.0

    # Insert current positions (using _dollars fields from API)
    if positions:
        conn.executemany("""
            INSERT INTO positions VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            (
                fetch_time,
                p.get("ticker", ""),
                p.get("position", 0),
                to_dollars(p, "market_exposure_dollars", "market_exposure"),
                to_dollars(p, "realized_pnl_dollars", "realized_pnl"),
                p.get("resting_orders_count", 0),
                to_dollars(p, "total_traded_dollars", "total_traded"),
                to_dollars(p, "fees_paid_dollars", "fees_paid")
            )
            for p in positions
        ])

    conn.close()
    return changes


def save_orders_to_duckdb(orders: list, db_path: str = None):
    """Save resting orders to DuckDB."""
    if db_path is None:
        db_path = Path(__file__).parent / "kalshi_coaching.duckdb"

    conn = duckdb.connect(str(db_path))
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    # Create orders table if not exists
    conn.execute("""
        CREATE TABLE IF NOT EXISTS resting_orders (
            fetch_time TIMESTAMP,
            order_id VARCHAR,
            ticker VARCHAR,
            side VARCHAR,
            type VARCHAR,
            yes_price INTEGER,
            no_price INTEGER,
            remaining_count INTEGER,
            created_time TIMESTAMP,
            expiration_time TIMESTAMP
        )
    """)

    # Insert current orders
    if orders:
        conn.executemany("""
            INSERT INTO resting_orders VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            (
                fetch_time,
                o.get("order_id", ""),
                o.get("ticker", ""),
                o.get("side", ""),
                o.get("type", ""),
                o.get("yes_price", 0),
                o.get("no_price", 0),
                o.get("remaining_count", 0),
                o.get("created_time"),
                o.get("expiration_time")
            )
            for o in orders
        ])

    conn.close()


def get_expiration_for_tickers(tickers: list, db_path: str = None) -> dict:
    """Fetch expiration times for a list of tickers and cache in DB."""
    if db_path is None:
        db_path = Path(__file__).parent / "kalshi_coaching.duckdb"

    conn = duckdb.connect(str(db_path))

    # Create market_info table if not exists
    conn.execute("""
        CREATE TABLE IF NOT EXISTS market_info (
            ticker VARCHAR PRIMARY KEY,
            title VARCHAR,
            subtitle VARCHAR,
            expiration_time TIMESTAMP,
            close_time TIMESTAMP,
            updated_at TIMESTAMP
        )
    """)

    # Check which tickers we already have cached
    existing = set()
    try:
        rows = conn.execute("SELECT ticker FROM market_info").fetchall()
        existing = {row[0] for row in rows}
    except Exception:
        pass

    expirations = {}
    updated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    for ticker in tickers:
        # Fetch from API if not cached
        if ticker not in existing:
            market = fetch_market_details(ticker)
            if market:
                info = {
                    "expiration_time": market.get("expiration_time"),
                    "close_time": market.get("close_time"),
                    "title": market.get("title", ""),
                    "subtitle": market.get("yes_sub_title", "")
                }
                expirations[ticker] = info

                # Save to cache
                conn.execute("""
                    INSERT OR REPLACE INTO market_info VALUES (?, ?, ?, ?, ?, ?)
                """, [
                    ticker,
                    info["title"],
                    info["subtitle"],
                    info["expiration_time"],
                    info["close_time"],
                    updated_at
                ])
            time.sleep(0.05)  # Small delay
        else:
            # Load from cache
            row = conn.execute(
                "SELECT title, subtitle, expiration_time, close_time FROM market_info WHERE ticker = ?",
                [ticker]
            ).fetchone()
            if row:
                expirations[ticker] = {
                    "title": row[0],
                    "subtitle": row[1],
                    "expiration_time": row[2],
                    "close_time": row[3]
                }

    conn.close()
    return expirations


def display_portfolio(positions: list, orders: list, changes: list, expirations: dict):
    """Display portfolio positions, orders, and changes."""
    print("\n" + "="*70)
    print("PORTFOLIO SUMMARY")
    print("="*70)

    # Helper to get dollar value from API response
    def get_dollars(p, dollar_key, cents_key):
        val = p.get(dollar_key)
        if val is not None:
            return float(val)
        cents = p.get(cents_key, 0)
        return float(cents) / 100.0 if cents else 0.0

    # Display positions
    if positions:
        print(f"\nOPEN POSITIONS ({len(positions)}):")
        print("-"*70)
        total_exposure = 0
        for p in sorted(positions, key=lambda x: abs(get_dollars(x, "market_exposure_dollars", "market_exposure")), reverse=True):
            ticker = p.get("ticker", "")
            pos = p.get("position", 0)
            # Use _dollars field from API (cents values also available)
            exp = get_dollars(p, "market_exposure_dollars", "market_exposure")
            total_exposure += exp

            # Get expiration
            exp_info = expirations.get(ticker, {})
            exp_time = exp_info.get("expiration_time", "")
            subtitle = exp_info.get("subtitle", ticker)

            # Format expiration
            exp_str = ""
            if exp_time:
                try:
                    exp_dt = datetime.fromisoformat(exp_time.replace("Z", "+00:00"))
                    days_left = (exp_dt - datetime.now(timezone.utc)).days
                    exp_str = f" | Expires: {exp_dt.strftime('%b %d')} ({days_left}d)"
                except Exception:
                    exp_str = ""

            # Check for changes
            change_str = ""
            for c in changes:
                if c["ticker"] == ticker:
                    if c["pos_change"] != 0:
                        sign = "+" if c["pos_change"] > 0 else ""
                        change_str = f" [CHG: {sign}{c['pos_change']} contracts]"
                    break

            side = "YES" if pos > 0 else "NO"
            print(f"  {subtitle[:35]:<35} | {abs(pos):>4} {side} | ${exp:>8.2f}{exp_str}{change_str}")

        print(f"\n  Total Exposure: ${total_exposure:,.2f}")
    else:
        print("\nNo open positions")

    # Display changes since last run
    if changes:
        print(f"\nCHANGES SINCE LAST RUN:")
        print("-"*70)
        for c in changes:
            exp_info = expirations.get(c["ticker"], {})
            subtitle = exp_info.get("subtitle", c["ticker"])
            sign = "+" if c["pos_change"] > 0 else ""
            exp_sign = "+" if c["exp_change"] > 0 else ""
            print(f"  {subtitle[:40]:<40} | Pos: {sign}{c['pos_change']:>3} | Exp: {exp_sign}${c['exp_change']:.2f}")

    # Display resting orders
    if orders:
        print(f"\nRESTING ORDERS ({len(orders)}):")
        print("-"*70)
        for o in orders:
            ticker = o.get("ticker", "")
            side = o.get("side", "").upper()
            remaining = o.get("remaining_count", 0)
            price = o.get("yes_price", 0) if side == "YES" else o.get("no_price", 0)

            exp_info = expirations.get(ticker, {})
            subtitle = exp_info.get("subtitle", ticker)

            # Order expiration
            order_exp = o.get("expiration_time")
            order_exp_str = ""
            if order_exp:
                try:
                    exp_dt = datetime.fromisoformat(order_exp.replace("Z", "+00:00"))
                    order_exp_str = f" | Order exp: {exp_dt.strftime('%b %d %H:%M')}"
                except Exception:
                    pass

            print(f"  {subtitle[:35]:<35} | {remaining:>3}x {side} @ {price}c{order_exp_str}")
    else:
        print("\nNo resting orders")


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

    # Check for API credentials and fetch portfolio
    creds = load_credentials()
    if creds["api_key_id"] and creds["private_key_path"]:
        print("\n" + "="*60)
        print("Fetching portfolio (positions & orders)...")

        positions = fetch_positions(creds["api_key_id"], creds["private_key_path"])
        orders = fetch_orders(creds["api_key_id"], creds["private_key_path"])

        if positions or orders:
            # Get all tickers we need expiration info for
            tickers = set()
            for p in positions:
                tickers.add(p.get("ticker", ""))
            for o in orders:
                tickers.add(o.get("ticker", ""))

            print(f"  Found {len(positions)} positions, {len(orders)} resting orders")
            print(f"  Fetching market info for {len(tickers)} tickers...")
            expirations = get_expiration_for_tickers(list(tickers), db_path)

            # Save and get changes
            changes = save_positions_to_duckdb(positions, db_path)
            save_orders_to_duckdb(orders, db_path)

            # Display portfolio
            display_portfolio(positions, orders, changes, expirations)
        else:
            print("  No positions or orders found (or auth failed)")
    else:
        print("\n[Portfolio tracking disabled - no API credentials]")
        print("  To enable: copy .env.example to .env and add your Kalshi credentials")

    return odds, db_path


if __name__ == "__main__":
    main()
