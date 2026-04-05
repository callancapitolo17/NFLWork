"""
Kalshi NFL Draft Market Fetcher
Auto-discovers all draft-related series and fetches odds + portfolio data.
"""

import time
import duckdb
from datetime import datetime, timezone
from pathlib import Path

from auth import (
    load_credentials, authenticated_request, public_request, BASE_URL
)
from db import DB_PATH, get_connection, init_schema


def discover_draft_series():
    """Auto-discover all NFL Draft related series from Kalshi API.

    Queries the /series endpoint and filters for draft-related series.
    Also tries known prefixes as fallback.
    """
    discovered = []

    # Strategy 1: Search series endpoint and filter for NFL Draft
    nfl_keywords = {"nfl", "pro football", "football draft"}
    cursor = None
    while True:
        path = "/series?limit=100"
        if cursor:
            path += f"&cursor={cursor}"

        data = public_request(path)
        if not data:
            break

        for s in data.get("series", []):
            ticker = s.get("ticker", "")
            title = (s.get("title") or "").lower()

            # Must be NFL-specific draft market
            is_nfl_draft = (
                (ticker.startswith("KXNFL") and "draft" in ticker.lower())
                or ("draft" in title and any(k in title for k in nfl_keywords))
            )
            if is_nfl_draft:
                discovered.append({
                    "series_ticker": ticker,
                    "title": s.get("title", ""),
                })

        cursor = data.get("cursor")
        if not cursor:
            break

        time.sleep(0.1)

    # Strategy 2: Try known prefixes as fallback
    known_prefixes = [
        "KXNFLDRAFT1",      # #1 overall pick (by player)
        "KXNFLDRAFT1ST",    # #1 pick (by team)
        "KXNFLDRAFTTOP5",   # Top 5 pick
        "KXNFLDRAFTTOP10",  # Top 10 pick
        "KXNFLDRAFTQB",     # First QB drafted
        "KXNFLDRAFTWR",     # First WR drafted
        "KXNFLDRAFTRB",     # First RB drafted
        "KXNFLDRAFTPOS",    # Draft position
    ]

    discovered_tickers = {d["series_ticker"] for d in discovered}

    for prefix in known_prefixes:
        if prefix not in discovered_tickers:
            # Check if this series has any open markets
            data = public_request(f"/markets?series_ticker={prefix}&status=open&limit=1")
            if data and data.get("markets"):
                discovered.append({
                    "series_ticker": prefix,
                    "title": data["markets"][0].get("title", prefix),
                })
                discovered_tickers.add(prefix)
                time.sleep(0.1)

    return discovered


def classify_series(series_ticker, title):
    """Classify a series into a category based on ticker/title patterns."""
    ticker_lower = series_ticker.lower()
    title_lower = (title or "").lower()

    if "1st" in ticker_lower or "team" in title_lower or "make the" in title_lower:
        return "team_pick"
    elif "top5" in ticker_lower or "top 5" in title_lower:
        return "range_top5"
    elif "top10" in ticker_lower or "top 10" in title_lower:
        return "range_top10"
    elif "qb" in ticker_lower or "quarterback" in title_lower:
        return "position_qb"
    elif "wr" in ticker_lower or "wide receiver" in title_lower:
        return "position_wr"
    elif "rb" in ticker_lower or "running back" in title_lower:
        return "position_rb"
    elif "draft1" in ticker_lower or "first pick" in title_lower or "#1" in title_lower:
        return "player_pick_1"
    elif "draft" in ticker_lower:
        return "player_general"
    return "unknown"


def fetch_markets_for_series(series_ticker):
    """Fetch all open markets for a given series ticker with pagination."""
    markets = []
    cursor = None

    while True:
        path = f"/markets?series_ticker={series_ticker}&status=open&limit=100"
        if cursor:
            path += f"&cursor={cursor}"

        data = public_request(path)
        if not data:
            break

        markets.extend(data.get("markets", []))

        cursor = data.get("cursor")
        if not cursor or len(data.get("markets", [])) == 0:
            break

        time.sleep(0.1)

    return markets


def extract_candidate(market):
    """Extract candidate name (player or team) from market data."""
    return (
        market.get("yes_sub_title")
        or market.get("custom_strike", {}).get("Person")
        or market.get("custom_strike", {}).get("Team")
        or market.get("no_sub_title")
        or "Unknown"
    )


def fetch_all_draft_odds():
    """Discover series, fetch all markets, return structured odds list."""
    print("Discovering NFL Draft series...")
    series_list = discover_draft_series()
    print(f"  Found {len(series_list)} draft series")

    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    all_odds = []

    # Save discovered series to DB
    con = get_connection()
    for s in series_list:
        category = classify_series(s["series_ticker"], s["title"])
        con.execute("""
            INSERT OR REPLACE INTO draft_series VALUES (?, ?, ?, ?)
        """, [s["series_ticker"], s["title"], category, fetch_time])
    con.close()

    # Fetch markets for each series
    for i, s in enumerate(series_list):
        if i > 0:
            time.sleep(0.1)

        series_ticker = s["series_ticker"]
        print(f"  Fetching {series_ticker} ({s['title']})...")
        markets = fetch_markets_for_series(series_ticker)

        for market in markets:
            candidate = extract_candidate(market)
            if candidate in ("Unknown", "", None):
                continue

            all_odds.append({
                "fetch_time": fetch_time,
                "series_ticker": series_ticker,
                "event_ticker": market.get("event_ticker", ""),
                "ticker": market.get("ticker", ""),
                "market_title": market.get("title", ""),
                "candidate": candidate,
                "yes_bid": market.get("yes_bid", 0),
                "yes_ask": market.get("yes_ask", 0),
                "no_bid": market.get("no_bid", 0),
                "no_ask": market.get("no_ask", 0),
                "last_price": market.get("last_price", 0),
                "volume": market.get("volume", 0),
                "volume_24h": market.get("volume_24h", 0),
                "liquidity": int(float(market.get("liquidity_dollars") or 0)),
                "open_interest": market.get("open_interest", 0),
            })

    return all_odds


def save_odds_snapshot(odds_data):
    """Append odds snapshot to draft_odds table."""
    if not odds_data:
        print("No odds data to save")
        return

    con = get_connection()
    con.executemany("""
        INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, [
        (
            d["fetch_time"], d["series_ticker"], d["event_ticker"],
            d["ticker"], d["market_title"], d["candidate"],
            d["yes_bid"], d["yes_ask"], d["no_bid"], d["no_ask"],
            d["last_price"], d["volume"], d["volume_24h"],
            d["liquidity"], d["open_interest"]
        )
        for d in odds_data
    ])

    result = con.execute("""
        SELECT COUNT(*) as total, COUNT(DISTINCT fetch_time) as snapshots
        FROM draft_odds
    """).fetchone()
    print(f"\nDatabase: {result[0]:,} total records across {result[1]} snapshots")
    con.close()


def cache_market_info(odds_data):
    """Cache market metadata for portfolio display."""
    con = get_connection()
    updated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    for d in odds_data:
        con.execute("""
            INSERT OR REPLACE INTO market_info (ticker, title, subtitle, series_ticker, updated_at)
            VALUES (?, ?, ?, ?, ?)
        """, [d["ticker"], d["market_title"], d["candidate"], d["series_ticker"], updated_at])

    con.close()


def fetch_portfolio():
    """Fetch positions and resting orders (requires auth)."""
    creds = load_credentials()
    if not creds["api_key_id"] or not creds["private_key_path"]:
        print("[Portfolio tracking disabled - no API credentials]")
        return [], []

    print("Fetching portfolio...")
    api_key = creds["api_key_id"]
    key_path = creds["private_key_path"]

    # Fetch positions
    positions = []
    cursor = None
    while True:
        path = "/portfolio/positions?limit=100"
        if cursor:
            path += f"&cursor={cursor}"
        data = authenticated_request("GET", path, api_key, key_path)
        if not data:
            break
        positions.extend(data.get("market_positions", []))
        cursor = data.get("cursor")
        if not cursor:
            break

    # Fetch resting orders
    orders = []
    cursor = None
    while True:
        path = "/portfolio/orders?status=resting&limit=100"
        if cursor:
            path += f"&cursor={cursor}"
        data = authenticated_request("GET", path, api_key, key_path)
        if not data:
            break
        orders.extend(data.get("orders", []))
        cursor = data.get("cursor")
        if not cursor:
            break

    print(f"  {len(positions)} positions, {len(orders)} resting orders")
    return positions, orders


def save_portfolio(positions, orders):
    """Save portfolio data to DuckDB."""
    con = get_connection()
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    def to_dollars(p, dollar_key, cents_key):
        val = p.get(dollar_key)
        if val is not None:
            return float(val)
        cents = p.get(cents_key, 0)
        return float(cents) / 100.0 if cents else 0.0

    if positions:
        con.executemany("""
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

    if orders:
        con.executemany("""
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

    con.close()


def run():
    """Main fetch pipeline: discover series, fetch odds, fetch portfolio."""
    init_schema()

    # Fetch odds
    odds = fetch_all_draft_odds()
    print(f"\nFound {len(odds)} markets across {len(set(d['series_ticker'] for d in odds))} series")

    # Print summary
    by_series = {}
    for o in odds:
        by_series.setdefault(o["series_ticker"], []).append(o)

    for series, markets in sorted(by_series.items()):
        top = sorted(markets, key=lambda x: x["last_price"], reverse=True)[:3]
        print(f"\n{series} ({len(markets)} markets):")
        for m in top:
            print(f"  {m['candidate']}: {m['last_price']}% (vol: {m['volume']:,})")

    # Save to DB
    save_odds_snapshot(odds)
    cache_market_info(odds)

    # Portfolio
    positions, orders = fetch_portfolio()
    if positions or orders:
        save_portfolio(positions, orders)

    return odds


if __name__ == "__main__":
    run()
