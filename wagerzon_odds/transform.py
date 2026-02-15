"""
Transform Wagerzon odds data to The Odds API format
"""

import duckdb
from pathlib import Path
from typing import Optional
from datetime import datetime

from team_mapping import normalize_team_name

DB_PATH = Path(__file__).parent / "wagerzon.duckdb"


def get_wagerzon_odds(sport: str, table_name: Optional[str] = None) -> list[dict]:
    """
    Retrieve Wagerzon odds from database.

    Args:
        sport: Sport key (nfl, nba, cbb)
        table_name: Override table name (default: {sport}_odds)

    Returns:
        List of odds records
    """
    if table_name is None:
        table_name = f"{sport}_odds"

    conn = duckdb.connect(str(DB_PATH), read_only=True)

    try:
        rows = conn.execute(f"SELECT * FROM {table_name}").fetchall()
        columns = [desc[0] for desc in conn.description]

        return [dict(zip(columns, row)) for row in rows]
    finally:
        conn.close()


def transform_to_odds_api_format(
    wagerzon_odds: list[dict],
    sport: str = "nfl"
) -> list[dict]:
    """
    Transform Wagerzon odds to The Odds API format.

    The Odds API format has these key columns:
    - bookmaker_key: "wagerzon"
    - market: "spreads", "spreads_h1", "totals_q1", etc.
    - line: the spread or total value
    - odds_away, odds_home: American odds for spreads
    - odds_over, odds_under: American odds for totals
    - home_team, away_team: Full team names
    - market_type: "spreads" or "totals"

    Args:
        wagerzon_odds: Raw odds from Wagerzon database
        sport: Sport key for team name normalization

    Returns:
        List of records in Odds API format
    """
    transformed = []

    for record in wagerzon_odds:
        # Normalize team names
        home_team = normalize_team_name(record["home_team"], sport)
        away_team = normalize_team_name(record["away_team"], sport)

        # Extract market info
        market = record["market"]
        period = record["period"]

        # Create spread record if we have spread data
        if record["away_spread"] is not None:
            spread_record = {
                "bookmaker_key": "wagerzon",
                "sport_key": record["sport_key"],
                "market": market,
                "market_type": "spreads",
                "line": record["home_spread"],  # Line is typically from home perspective
                "odds_away": record["away_spread_price"],
                "odds_home": record["home_spread_price"],
                "odds_over": None,
                "odds_under": None,
                "home_team": home_team,
                "away_team": away_team,
                "game_date": record["game_date"],
                "game_time": record["game_time"],
                "game_id": record["game_id"],
                "fetch_time": record["fetch_time"],
                "period": period,
                # Store both spread values for reference
                "away_spread": record["away_spread"],
                "home_spread": record["home_spread"],
            }
            transformed.append(spread_record)

        # Create totals record if we have totals data
        if record["total"] is not None:
            # Determine totals market name
            totals_market = market.replace("spreads", "totals")

            totals_record = {
                "bookmaker_key": "wagerzon",
                "sport_key": record["sport_key"],
                "market": totals_market,
                "market_type": "totals",
                "line": record["total"],
                "odds_away": None,
                "odds_home": None,
                "odds_over": record["over_price"],
                "odds_under": record["under_price"],
                "home_team": home_team,
                "away_team": away_team,
                "game_date": record["game_date"],
                "game_time": record["game_time"],
                "game_id": record["game_id"],
                "fetch_time": record["fetch_time"],
                "period": period,
            }
            transformed.append(totals_record)

        # Create moneyline record if we have ML data
        if record["away_ml"] is not None:
            # Determine ML market name
            ml_market = market.replace("spreads", "h2h")

            ml_record = {
                "bookmaker_key": "wagerzon",
                "sport_key": record["sport_key"],
                "market": ml_market,
                "market_type": "h2h",
                "line": None,
                "odds_away": record["away_ml"],
                "odds_home": record["home_ml"],
                "odds_over": None,
                "odds_under": None,
                "home_team": home_team,
                "away_team": away_team,
                "game_date": record["game_date"],
                "game_time": record["game_time"],
                "game_id": record["game_id"],
                "fetch_time": record["fetch_time"],
                "period": period,
            }
            transformed.append(ml_record)

    return transformed


def get_transformed_odds(sport: str = "nfl") -> list[dict]:
    """
    Get Wagerzon odds transformed to Odds API format.

    Args:
        sport: Sport key

    Returns:
        List of transformed odds records
    """
    raw_odds = get_wagerzon_odds(sport)
    return transform_to_odds_api_format(raw_odds, sport)


def print_odds_comparison(sport: str = "nfl"):
    """Print Wagerzon odds in a readable format for comparison."""
    odds = get_transformed_odds(sport)

    print(f"\n=== Wagerzon {sport.upper()} Odds (Odds API Format) ===\n")

    # Group by game
    games = {}
    for record in odds:
        key = f"{record['away_team']} @ {record['home_team']}"
        if key not in games:
            games[key] = []
        games[key].append(record)

    for game, records in games.items():
        print(f"\n{game}")
        print("-" * 60)

        for r in sorted(records, key=lambda x: (x["market_type"], x["market"])):
            if r["market_type"] == "spreads":
                print(f"  {r['market']:15} | Line: {r['line']:+.1f} | Away: {r['odds_away']:+d} Home: {r['odds_home']:+d}")
            elif r["market_type"] == "totals":
                print(f"  {r['market']:15} | Total: {r['line']:.1f} | Over: {r['odds_over']:+d} Under: {r['odds_under']:+d}")
            elif r["market_type"] == "h2h":
                print(f"  {r['market']:15} | ML | Away: {r['odds_away']:+d} Home: {r['odds_home']:+d}")


if __name__ == "__main__":
    import sys
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    print_odds_comparison(sport)
