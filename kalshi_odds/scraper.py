#!/usr/bin/env python3
"""
Kalshi CBB 1H Odds Scraper for Answer Key Pipeline
Fetches 1st half spreads, totals, and moneylines from Kalshi's public API.

Kalshi is a prediction market exchange — prices are probabilities (0-100 cents).
This scraper converts to American odds with the taker fee baked in so EV
calculations reflect actual executable cost.

No auth required for market data.
"""

import sys
import re
import json
import time
import urllib.request
import duckdb
from datetime import datetime, timezone
from pathlib import Path
from collections import defaultdict

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

KALSHI_BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"


def public_request(path):
    """Make a public (unauthenticated) request to Kalshi API."""
    url = f"{KALSHI_BASE_URL}{path}"
    try:
        req = urllib.request.Request(url)
        req.add_header("Content-Type", "application/json")
        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        print(f"  Request failed: {e.code} - {e.read().decode()}")
        return None
    except Exception as e:
        print(f"  Request error: {e}")
        return None

DB_PATH = Path(__file__).parent / "kalshi.duckdb"

# Series tickers for CBB 1H markets
SERIES_TICKERS = {
    "spreads": "KXNCAAMB1HSPREAD",
    "totals": "KXNCAAMB1HTOTAL",
    "moneyline": "KXNCAAMB1HWINNER",
}

TAKER_FEE_RATE = 0.07


# =============================================================================
# PRICE CONVERSION
# =============================================================================


def cents_to_american(ask_cents):
    """Convert Kalshi ask price (cents 0-100) to American odds with taker fee.

    The taker fee is 7% * P * (1-P) per contract. We add this to the ask price
    to get the effective cost, then convert to American odds.

    Returns (american_odds, effective_cents) tuple.
    """
    if ask_cents is None or ask_cents <= 0 or ask_cents >= 100:
        return None, None
    p = ask_cents / 100.0
    fee = TAKER_FEE_RATE * p * (1 - p)
    effective_p = p + fee
    if effective_p >= 1.0:
        return None, None
    effective_cents = round(effective_p * 100, 2)
    if effective_p > 0.5:
        return round(-100 * effective_p / (1 - effective_p)), effective_cents
    else:
        return round(100 * (1 - effective_p) / effective_p), effective_cents


# =============================================================================
# MARKET FETCHING
# =============================================================================


def fetch_markets(series_ticker):
    """Fetch all open markets for a series with cursor pagination."""
    markets = []
    cursor = None

    while True:
        path = f"/markets?series_ticker={series_ticker}&status=open&limit=200"
        if cursor:
            path += f"&cursor={cursor}"

        data = public_request(path)
        if not data:
            break

        batch = data.get("markets", [])
        markets.extend(batch)

        cursor = data.get("cursor")
        if not cursor or len(batch) == 0:
            break

        time.sleep(0.1)

    return markets


# =============================================================================
# TEAM NAME PARSING
# =============================================================================


def parse_matchup_title(title):
    """Parse game title to extract team names.

    Handles formats:
    - "Away vs Home: First Half Winner?"
    - "Away at Home: ..."
    - "Will TeamName win the 1H by over X.5 points?"
    - "Away vs Home: First Half Total?"
    """
    # Remove trailing question mark and common suffixes
    clean = title.strip().rstrip("?")
    clean = re.sub(r":\s*(First Half (Winner|Total))$", "", clean)

    # "Away vs Home" or "Away at Home"
    for sep in [" vs ", " at "]:
        if sep in clean:
            parts = clean.split(sep, 1)
            return parts[0].strip(), parts[1].strip()

    return None, None


def parse_spread_team(title):
    """Extract team name from spread market title.

    Format: "Will TeamName win the 1H by over X.5 points?"
    or: "TeamName wins by over X.5 Points?"
    """
    m = re.match(r"Will (.+?) win the 1H by over", title)
    if m:
        return m.group(1).strip()
    m = re.match(r"(.+?) wins (?:the 1H )?by over", title)
    if m:
        return m.group(1).strip()
    return None


def resolve_home_away(away_raw, home_raw, team_dict, canonical_games):
    """Resolve team names and determine home/away using canonical games."""
    return resolve_team_names(away_raw, home_raw, team_dict, canonical_games)


# =============================================================================
# MARKET PARSERS
# =============================================================================


def parse_spread_records(markets, team_dict, canonical_games, fetch_time):
    """Parse 1H spread markets into 18-column records.

    Each Kalshi spread contract like "Team wins by >X" independently provides
    BOTH sides of a traditional spread:
      - YES = Team -X (team covers)   → use yes_ask
      - NO  = Opponent +X (team doesn't cover) → use no_ask (100 - yes_bid)

    We process each contract individually — no pairing needed.
    Multiple contracts at the same strike for the same event will produce
    multiple records; the R pipeline picks the best line.
    """
    records = []
    by_event = defaultdict(list)
    for m in markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Collect all team names mentioned across contracts
        team_names = set()
        for m in event_markets:
            team = parse_spread_team(m.get("title", ""))
            if team:
                team_names.add(team)

        if len(team_names) != 2:
            continue

        team_list = sorted(team_names)

        # Resolve teams
        away_resolved, home_resolved = resolve_home_away(
            team_list[0], team_list[1], team_dict, canonical_games
        )

        # Map raw Kalshi names to home/away
        raw_to_side = {}
        for raw_name in team_list:
            if _fuzzy_team_match(raw_name.lower(), away_resolved.lower()):
                raw_to_side[raw_name] = "away"
            elif _fuzzy_team_match(raw_name.lower(), home_resolved.lower()):
                raw_to_side[raw_name] = "home"

        if len(raw_to_side) != 2:
            raw_to_side = {team_list[0]: "away", team_list[1]: "home"}

        close_time = event_markets[0].get("close_time", "")
        game_date, game_time_str = _parse_datetime(close_time)

        for m in event_markets:
            strike = m.get("floor_strike")
            if strike is None:
                continue

            if not _is_liquid(m):
                continue

            contract_team = parse_spread_team(m.get("title", ""))
            if not contract_team or contract_team not in raw_to_side:
                continue

            side = raw_to_side[contract_team]

            # "Team wins by >X" → YES = team -X, NO = opponent +X
            yes_ask = m.get("yes_ask", 0)
            # NO ask = 100 - yes_bid (what you pay to bet NO)
            yes_bid = m.get("yes_bid", 0)
            no_ask = 100 - yes_bid if yes_bid > 0 else 0

            cover_odds, cover_eff = cents_to_american(yes_ask)     # team -X
            opponent_odds, opponent_eff = cents_to_american(no_ask)    # opponent +X

            if cover_odds is None or opponent_odds is None:
                continue

            if side == "home":
                # Contract is for home team: home -X, away +X
                record = {
                    "home_spread": -strike,
                    "home_spread_price": cover_odds,
                    "home_spread_cents": cover_eff,
                    "away_spread": strike,
                    "away_spread_price": opponent_odds,
                    "away_spread_cents": opponent_eff,
                }
            else:
                # Contract is for away team: away -X, home +X
                record = {
                    "away_spread": -strike,
                    "away_spread_price": cover_odds,
                    "away_spread_cents": cover_eff,
                    "home_spread": strike,
                    "home_spread_price": opponent_odds,
                    "home_spread_cents": opponent_eff,
                }

            record.update({
                "fetch_time": fetch_time,
                "sport_key": "basketball_ncaab",
                "game_id": f"kalshi-{event_ticker}-{m['ticker']}",
                "game_date": game_date,
                "game_time": game_time_str,
                "away_team": away_resolved,
                "home_team": home_resolved,
                "market": "spreads_h1",
                "period": "Half1",
                "total": None,
                "over_price": None,
                "over_cents": None,
                "under_price": None,
                "under_cents": None,
                "away_ml": None,
                "away_ml_cents": None,
                "home_ml": None,
                "home_ml_cents": None,
            })

            records.append(record)

    return records


def parse_total_records(markets, team_dict, canonical_games, fetch_time):
    """Parse 1H total markets into 18-column records.

    Each event has multiple total lines. Each line has 1 contract:
    YES = over, NO = under.
    """
    records = []
    by_event = defaultdict(list)
    for m in markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Parse teams from event title
        title = event_markets[0].get("title", "")
        away_raw, home_raw = parse_matchup_title(title)
        if not away_raw or not home_raw:
            continue

        away_resolved, home_resolved = resolve_home_away(
            away_raw, home_raw, team_dict, canonical_games
        )

        close_time = event_markets[0].get("close_time", "")
        game_date, game_time_str = _parse_datetime(close_time)

        for m in event_markets:
            strike = m.get("floor_strike")
            if strike is None:
                continue

            if not _is_liquid(m):
                continue

            over_ask = m.get("yes_ask", 0)
            under_ask = m.get("no_ask", 0)

            over_odds, over_eff = cents_to_american(over_ask)
            under_odds, under_eff = cents_to_american(under_ask)

            if over_odds is None or under_odds is None:
                continue

            records.append({
                "fetch_time": fetch_time,
                "sport_key": "basketball_ncaab",
                "game_id": f"kalshi-{event_ticker}-{int(strike)}",
                "game_date": game_date,
                "game_time": game_time_str,
                "away_team": away_resolved,
                "home_team": home_resolved,
                "market": "totals_h1",
                "period": "Half1",
                "away_spread": None,
                "away_spread_price": None,
                "away_spread_cents": None,
                "home_spread": None,
                "home_spread_price": None,
                "home_spread_cents": None,
                "total": strike,
                "over_price": over_odds,
                "over_cents": over_eff,
                "under_price": under_odds,
                "under_cents": under_eff,
                "away_ml": None,
                "away_ml_cents": None,
                "home_ml": None,
                "home_ml_cents": None,
            })

    return records


def parse_moneyline_records(markets, team_dict, canonical_games, fetch_time):
    """Parse 1H winner (3-way) markets into 18-column records.

    Each event has 3 contracts: team A, team B, and tie.
    We only output the two team contracts as a moneyline pair.
    """
    records = []
    by_event = defaultdict(list)
    for m in markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Separate tie from team contracts
        team_contracts = []
        for m in event_markets:
            sub = (m.get("yes_sub_title") or "").lower()
            if "tie" in sub:
                continue
            team_contracts.append(m)

        if len(team_contracts) != 2:
            continue

        # Parse teams from event title
        title = event_markets[0].get("title", "")
        away_raw, home_raw = parse_matchup_title(title)
        if not away_raw or not home_raw:
            continue

        away_resolved, home_resolved = resolve_home_away(
            away_raw, home_raw, team_dict, canonical_games
        )

        close_time = event_markets[0].get("close_time", "")
        game_date, game_time_str = _parse_datetime(close_time)

        # Match contracts to home/away
        # The yes_sub_title contains the team name
        home_contract = None
        away_contract = None
        for m in team_contracts:
            sub = (m.get("yes_sub_title") or "").strip()
            # Try matching against raw names
            sub_lower = sub.lower()
            if _fuzzy_team_match(sub_lower, home_raw):
                home_contract = m
            elif _fuzzy_team_match(sub_lower, away_raw):
                away_contract = m

        # Fallback: try order from title (first = away, second = home)
        if not home_contract or not away_contract:
            # Just assign by title order
            c0_sub = (team_contracts[0].get("yes_sub_title") or "").strip().lower()
            if _fuzzy_team_match(c0_sub, away_raw):
                away_contract = team_contracts[0]
                home_contract = team_contracts[1]
            else:
                home_contract = team_contracts[0]
                away_contract = team_contracts[1]

        if not home_contract or not away_contract:
            continue

        # Skip illiquid
        if not _is_liquid(home_contract) or not _is_liquid(away_contract):
            continue

        home_ask = home_contract.get("yes_ask", 0)
        away_ask = away_contract.get("yes_ask", 0)

        home_odds, home_eff = cents_to_american(home_ask)
        away_odds, away_eff = cents_to_american(away_ask)

        if home_odds is None or away_odds is None:
            continue

        records.append({
            "fetch_time": fetch_time,
            "sport_key": "basketball_ncaab",
            "game_id": f"kalshi-{event_ticker}",
            "game_date": game_date,
            "game_time": game_time_str,
            "away_team": away_resolved,
            "home_team": home_resolved,
            "market": "h2h_h1",
            "period": "Half1",
            "away_spread": None,
            "away_spread_price": None,
            "away_spread_cents": None,
            "home_spread": None,
            "home_spread_price": None,
            "home_spread_cents": None,
            "total": None,
            "over_price": None,
            "over_cents": None,
            "under_price": None,
            "under_cents": None,
            "away_ml": away_odds,
            "away_ml_cents": away_eff,
            "home_ml": home_odds,
            "home_ml_cents": home_eff,
        })

    return records


# =============================================================================
# HELPERS
# =============================================================================


def _is_liquid(market, max_spread=20):
    """Check if a market has sufficient liquidity to be actionable."""
    ask = market.get("yes_ask", 0)
    bid = market.get("yes_bid", 0)
    if ask <= 0:
        return False
    if bid <= 0:
        return False
    if (ask - bid) > max_spread:
        return False
    return True


def _parse_datetime(iso_str):
    """Parse ISO datetime string to game_date (MM/DD) and game_time (HH:MM UTC)."""
    if not iso_str:
        return "", ""
    try:
        dt = datetime.fromisoformat(iso_str.replace("Z", "+00:00"))
        return dt.strftime("%m/%d"), dt.strftime("%H:%M")
    except (ValueError, TypeError):
        return "", ""


def _fuzzy_team_match(name1, name2):
    """Simple fuzzy check: one name is a substring of the other (case-insensitive)."""
    if not name1 or not name2:
        return False
    n1 = re.sub(r"['\.\-]", "", name1).lower().strip()
    n2 = re.sub(r"['\.\-]", "", name2).lower().strip()
    return n1 in n2 or n2 in n1


# =============================================================================
# DATABASE
# =============================================================================


def init_database():
    """Initialize DuckDB with the cbb_odds table."""
    conn = duckdb.connect(str(DB_PATH))
    try:
        # Recreate table to pick up schema changes (data is ephemeral — cleared each run)
        conn.execute("DROP TABLE IF EXISTS cbb_odds")
        conn.execute("""
            CREATE TABLE cbb_odds (
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
                away_spread_cents FLOAT,
                home_spread FLOAT,
                home_spread_price INTEGER,
                home_spread_cents FLOAT,
                total FLOAT,
                over_price INTEGER,
                over_cents FLOAT,
                under_price INTEGER,
                under_cents FLOAT,
                away_ml INTEGER,
                away_ml_cents FLOAT,
                home_ml INTEGER,
                home_ml_cents FLOAT
            )
        """)
    finally:
        conn.close()


def save_to_database(odds_data):
    """Save scraped odds to DuckDB (clear + insert)."""
    conn = duckdb.connect(str(DB_PATH))
    try:
        columns = [
            "fetch_time", "sport_key", "game_id", "game_date", "game_time",
            "away_team", "home_team", "market", "period",
            "away_spread", "away_spread_price", "away_spread_cents",
            "home_spread", "home_spread_price", "home_spread_cents",
            "total", "over_price", "over_cents", "under_price", "under_cents",
            "away_ml", "away_ml_cents", "home_ml", "home_ml_cents"
        ]

        placeholders = ", ".join(["?" for _ in columns])

        conn.execute("DELETE FROM cbb_odds")

        conn.executemany(f"""
            INSERT INTO cbb_odds ({", ".join(columns)})
            VALUES ({placeholders})
        """, [
            tuple(d[col] for col in columns)
            for d in odds_data
        ])

        result = conn.execute("SELECT COUNT(*) FROM cbb_odds").fetchone()
        print(f"Database now has {result[0]} total records in cbb_odds")
    finally:
        conn.close()


# =============================================================================
# MAIN SCRAPER
# =============================================================================


def scrape_kalshi(sport="cbb"):
    """Scrape Kalshi CBB 1H odds and save to DuckDB."""
    if sport != "cbb":
        print(f"Kalshi scraper only supports CBB currently, got: {sport}")
        return []

    init_database()
    fetch_time = datetime.now(timezone.utc)

    # Load team name resolution
    team_dict = load_team_dict("cbb")
    canonical_games = load_canonical_games("cbb")

    all_records = []

    # Fetch and parse each market type
    for market_type, series_ticker in SERIES_TICKERS.items():
        print(f"Fetching Kalshi {market_type} ({series_ticker})...")
        markets = fetch_markets(series_ticker)
        print(f"  Found {len(markets)} open markets")

        if not markets:
            continue

        if market_type == "spreads":
            records = parse_spread_records(markets, team_dict, canonical_games, fetch_time)
        elif market_type == "totals":
            records = parse_total_records(markets, team_dict, canonical_games, fetch_time)
        elif market_type == "moneyline":
            records = parse_moneyline_records(markets, team_dict, canonical_games, fetch_time)
        else:
            continue

        all_records.extend(records)
        print(f"  Parsed {len(records)} records")

    # Save to DuckDB
    if all_records:
        save_to_database(all_records)

    print(f"\nScraped {len(all_records)} total Kalshi 1H records")
    return all_records


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "cbb"
    scrape_kalshi(sport)
