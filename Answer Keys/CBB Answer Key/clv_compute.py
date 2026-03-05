#!/usr/bin/env python3
"""
CLV (Closing Line Value) Computation

End-of-day batch job that computes CLV for all placed bets whose games
have completed. Fetches Pinnacle closing odds via Odds API historical
endpoint (market CLV) and uses offshore book snapshots (book CLV).

Usage: python clv_compute.py
"""
import json
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import duckdb
import requests
from scipy.stats import norm

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).parent.resolve()
DASHBOARD_DB = BASE_DIR.parent / "CBB Dashboard" / "cbb_dashboard.duckdb"
SPORT_KEY = "basketball_ncaab"
API_BASE = "https://api.the-odds-api.com/v4"

# Sharp books to use for market CLV (in priority order)
SHARP_BOOKS = ["pinnacle", "draftkings", "fanduel", "betmgm"]

# Markets we bet on and need closing odds for
MARKETS = "spreads_h1,totals_h1,team_totals_h1,spreads,totals,team_totals"

API_DELAY = 0.6  # seconds between API calls to respect rate limits

# Standard deviations by market type for normal CDF re-pricing
SIGMA = {
    "totals_h1": 6.5,
    "totals": 11.0,
    "team_totals_h1": 4.5,
    "team_totals_home_h1": 4.5,
    "team_totals_away_h1": 4.5,
    "team_totals": 7.0,
    "spreads_h1": 6.0,
    "spreads": 10.0,
}


# =============================================================================
# PROBABILITY MATH
# =============================================================================

def odds_to_prob(american: int) -> float:
    """Convert American odds to implied probability."""
    if american < 0:
        return (-american) / (-american + 100)
    return 100 / (american + 100)


def devig(odds_side: int, odds_other: int) -> float:
    """Remove vig and return true probability for the side."""
    p1 = odds_to_prob(odds_side)
    p2 = odds_to_prob(odds_other)
    return p1 / (p1 + p2)


def reprice_at_line(closing_line, closing_novig_prob, placement_line,
                    market_type, bet_on):
    """Re-price using normal CDF when line has moved.

    Given closing odds at one line, compute the implied probability
    at a different line using normal distribution assumptions.

    For totals: X = total points. Over: P(X > line), Under: P(X < line).
    For spreads: X = home_margin. Home: P(X > -line) since line is stored
    as negative for favorites. Away: P(X < line).
    """
    sigma = SIGMA.get(market_type, 8.0)

    # Home spread lines are stored negative; flip to get the threshold
    if bet_on == "home":
        eff_close = -closing_line
        eff_place = -placement_line
    else:
        eff_close = closing_line
        eff_place = placement_line

    delta = (eff_place - eff_close) / sigma

    # "Over-type" bets: P(X > threshold), prob decreases as threshold rises
    if bet_on in ("over", "home"):
        z = norm.ppf(1 - closing_novig_prob)
        return 1 - norm.cdf(delta + z)

    # "Under-type" bets: P(X < threshold), prob increases as threshold rises
    # Covers: under, away
    z = norm.ppf(closing_novig_prob)
    return norm.cdf(delta + z)


# =============================================================================
# ODDS API
# =============================================================================

def fetch_closing_odds(event_id: str, commence_time: str, api_key: str) -> dict | None:
    """Fetch historical odds from Odds API at 15 min before game time.

    Returns dict of {(market, bet_on, line): (odds, counter_odds, bookmaker)} or None.
    """
    # Parse commence_time and snapshot at T-15min
    try:
        if isinstance(commence_time, str):
            ct = datetime.fromisoformat(commence_time.replace("Z", "+00:00"))
        else:
            ct = commence_time
        snapshot = (ct - timedelta(minutes=15)).strftime("%Y-%m-%dT%H:%M:%SZ")
    except Exception:
        return None

    resp = requests.get(
        f"{API_BASE}/historical/sports/{SPORT_KEY}/events/{event_id}/odds",
        params={
            "apiKey": api_key,
            "date": snapshot,
            "regions": "us,us2,eu",
            "markets": MARKETS,
            "oddsFormat": "american",
            "dateFormat": "iso",
        },
        timeout=30,
    )

    if resp.status_code != 200:
        return None

    data = resp.json().get("data")
    if not data or "bookmakers" not in data:
        return None

    # Parse bookmakers, preferring sharp books
    closing = {}
    for bm in data["bookmakers"]:
        bm_key = bm["key"]
        for mkt in bm.get("markets", []):
            market_key = mkt["key"]
            outcomes = mkt.get("outcomes", [])

            # Group outcomes into pairs for devigging
            pairs = {}
            for o in outcomes:
                name = o["name"].lower()
                point = o.get("point")
                price = o["price"]

                # Determine bet_on
                if name == "over":
                    bet_on = "over"
                elif name == "under":
                    bet_on = "under"
                elif name == data.get("home_team", ""):
                    bet_on = "home"
                elif name == data.get("away_team", ""):
                    bet_on = "away"
                else:
                    bet_on = name

                key = (market_key, point)
                if key not in pairs:
                    pairs[key] = {}
                pairs[key][bet_on] = price

            # Store each side with its counter
            for (mk, pt), sides in pairs.items():
                for side, price in sides.items():
                    # Find counter side
                    counter_map = {"over": "under", "under": "over",
                                   "home": "away", "away": "home"}
                    counter = sides.get(counter_map.get(side))

                    lookup = (mk, side, pt)
                    # Prefer sharper books
                    existing_book = closing.get(lookup, (None, None, None))[2]
                    if existing_book is None or (
                        bm_key in SHARP_BOOKS and
                        (existing_book not in SHARP_BOOKS or
                         SHARP_BOOKS.index(bm_key) < SHARP_BOOKS.index(existing_book))
                    ):
                        closing[lookup] = (price, counter, bm_key)

    return closing


# =============================================================================
# CLV COMPUTATION
# =============================================================================

def compute_clv_for_bet(bet, market_closing, book_closing):
    """Compute market and book CLV for a single bet.

    Returns dict with CLV fields to insert into bet_clv table.
    """
    result = {
        "bet_hash": bet["bet_hash"],
        "game_id": bet["game_id"],
        "game_time": bet["game_time"],
        "bookmaker": bet["bookmaker"],
        "market": bet["market"],
        "bet_on": bet["bet_on"],
        "placement_line": bet["line"],
        "placement_odds": bet["odds"],
        "placement_novig_prob": None,
        "market_closing_line": None,
        "market_closing_odds": None,
        "market_closing_counter_odds": None,
        "market_closing_novig_prob": None,
        "market_clv": None,
        "book_closing_line": None,
        "book_closing_odds": None,
        "book_closing_counter_odds": None,
        "book_closing_novig_prob": None,
        "book_clv": None,
        "line_moved": None,
        "clv_method": None,
        "sigma_used": None,
    }

    bet_on = bet["bet_on"]
    placement_line = bet["line"]
    placement_odds = bet["odds"]
    market_type = bet["market"]

    # --- Market CLV ---
    if market_closing:
        # Try exact line match first, then any line for the same market/side
        mc = market_closing.get((market_type, bet_on, placement_line))
        method = "direct"
        closing_line = placement_line

        if mc is None:
            # Look for same market/side at any line
            candidates = {k: v for k, v in market_closing.items()
                          if k[0] == market_type and k[1] == bet_on}
            if candidates:
                # Pick the one closest to our line
                closest_key = min(candidates.keys(),
                                  key=lambda k: abs((k[2] or 0) - (placement_line or 0)))
                mc = candidates[closest_key]
                closing_line = closest_key[2]
                method = "normal_cdf"

        if mc:
            close_odds, close_counter, _ = mc
            result["market_closing_line"] = closing_line
            result["market_closing_odds"] = close_odds
            result["market_closing_counter_odds"] = close_counter

            if close_counter is not None:
                closing_novig = devig(close_odds, close_counter)
                result["market_closing_novig_prob"] = closing_novig

                if method == "direct":
                    # Same line: direct comparison
                    # We need counter odds for placement too — estimate from model_prob
                    placement_novig = bet.get("model_prob", odds_to_prob(placement_odds))
                    result["placement_novig_prob"] = placement_novig
                    result["market_clv"] = placement_novig - closing_novig
                    result["line_moved"] = False
                    result["clv_method"] = "direct"
                else:
                    # Line moved: re-price closing odds at our placement line
                    sigma = SIGMA.get(market_type, 8.0)
                    adj_prob = reprice_at_line(
                        closing_line, closing_novig, placement_line,
                        market_type, bet_on
                    )
                    placement_novig = bet.get("model_prob", odds_to_prob(placement_odds))
                    result["placement_novig_prob"] = placement_novig
                    result["market_clv"] = placement_novig - adj_prob
                    result["line_moved"] = True
                    result["clv_method"] = "normal_cdf"
                    result["sigma_used"] = sigma

    # --- Book CLV ---
    if book_closing:
        bc = book_closing.get((market_type, bet_on, placement_line))
        method = "direct"
        closing_line = placement_line

        if bc is None:
            candidates = {k: v for k, v in book_closing.items()
                          if k[0] == market_type and k[1] == bet_on}
            if candidates:
                closest_key = min(candidates.keys(),
                                  key=lambda k: abs((k[2] or 0) - (placement_line or 0)))
                bc = candidates[closest_key]
                closing_line = closest_key[2]
                method = "normal_cdf"

        if bc:
            close_odds, close_counter = bc
            result["book_closing_line"] = closing_line
            result["book_closing_odds"] = close_odds
            result["book_closing_counter_odds"] = close_counter

            if close_counter is not None:
                closing_novig = devig(close_odds, close_counter)
                result["book_closing_novig_prob"] = closing_novig

                if result["placement_novig_prob"] is None:
                    result["placement_novig_prob"] = bet.get("model_prob", odds_to_prob(placement_odds))

                if method == "direct":
                    result["book_clv"] = result["placement_novig_prob"] - closing_novig
                    if result["line_moved"] is None:
                        result["line_moved"] = False
                    if result["clv_method"] is None:
                        result["clv_method"] = "direct"
                else:
                    sigma = SIGMA.get(market_type, 8.0)
                    adj_prob = reprice_at_line(
                        closing_line, closing_novig, placement_line,
                        market_type, bet_on
                    )
                    result["book_clv"] = result["placement_novig_prob"] - adj_prob
                    if result["line_moved"] is None:
                        result["line_moved"] = True
                    if result["clv_method"] is None:
                        result["clv_method"] = "normal_cdf"
                        result["sigma_used"] = sigma

    return result


def load_book_closing(game_id: str, bookmaker: str, game_time) -> dict | None:
    """Load the latest closing snapshot for a specific book/game before tipoff.

    Returns dict of {(market, bet_on, line): (odds, counter_odds)} or None.
    """
    try:
        con = duckdb.connect(str(DASHBOARD_DB), read_only=True)
        rows = con.execute("""
            SELECT market, bet_on, line, odds, counter_odds
            FROM closing_snapshots
            WHERE bookmaker = ?
              AND snapshot_time = (
                  SELECT MAX(snapshot_time) FROM closing_snapshots
                  WHERE bookmaker = ? AND snapshot_time < ?
              )
        """, [bookmaker, bookmaker, game_time]).fetchall()
        con.close()
    except Exception:
        return None

    if not rows:
        return None

    closing = {}
    for market, bet_on, line, odds, counter_odds in rows:
        closing[(market, bet_on, line)] = (odds, counter_odds)

    return closing


# =============================================================================
# MAIN
# =============================================================================

def main():
    api_key = os.environ.get("ODDS_API_KEY")
    if not api_key:
        print("Error: ODDS_API_KEY environment variable not set")
        return 1

    if not DASHBOARD_DB.exists():
        print(f"Error: Dashboard DB not found at {DASHBOARD_DB}")
        return 1

    # Find placed bets for completed games that don't have CLV yet
    con = duckdb.connect(str(DASHBOARD_DB), read_only=True)
    bets = con.execute("""
        SELECT p.*
        FROM placed_bets p
        LEFT JOIN bet_clv c ON p.bet_hash = c.bet_hash
        WHERE p.game_time < CURRENT_TIMESTAMP
          AND c.bet_hash IS NULL
        ORDER BY p.game_time
    """).fetchdf()
    con.close()

    if bets.empty:
        print("No new bets to compute CLV for.")
        return 0

    print(f"Computing CLV for {len(bets)} bets across "
          f"{bets['game_id'].nunique()} games...")

    # Group bets by game for efficient API calls
    games = bets.groupby("game_id").first()[["game_time"]].to_dict("index")

    # Fetch market closing odds for each game
    market_closings = {}
    for i, (game_id, info) in enumerate(games.items()):
        print(f"  Fetching market closing odds for game {i+1}/{len(games)}: {game_id}")
        closing = fetch_closing_odds(game_id, str(info["game_time"]), api_key)
        if closing:
            market_closings[game_id] = closing
            print(f"    Got {len(closing)} closing lines")
        else:
            print(f"    No closing data available")

        if i < len(games) - 1:
            time.sleep(API_DELAY)

    # Compute CLV for each bet
    results = []
    for _, bet in bets.iterrows():
        bet_dict = bet.to_dict()
        market_closing = market_closings.get(bet_dict["game_id"])
        book_closing = load_book_closing(
            bet_dict["game_id"], bet_dict["bookmaker"], bet_dict["game_time"]
        )

        clv = compute_clv_for_bet(bet_dict, market_closing, book_closing)
        results.append(clv)

    # Save results
    if results:
        con = duckdb.connect(str(DASHBOARD_DB))
        for r in results:
            con.execute("""
                INSERT INTO bet_clv (
                    bet_hash, game_id, game_time, bookmaker, market, bet_on,
                    placement_line, placement_odds, placement_novig_prob,
                    market_closing_line, market_closing_odds, market_closing_counter_odds,
                    market_closing_novig_prob, market_clv,
                    book_closing_line, book_closing_odds, book_closing_counter_odds,
                    book_closing_novig_prob, book_clv,
                    line_moved, clv_method, sigma_used, computed_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
                ON CONFLICT (bet_hash) DO NOTHING
            """, [
                r["bet_hash"], r["game_id"], r["game_time"], r["bookmaker"],
                r["market"], r["bet_on"],
                r["placement_line"], r["placement_odds"], r["placement_novig_prob"],
                r["market_closing_line"], r["market_closing_odds"],
                r["market_closing_counter_odds"], r["market_closing_novig_prob"],
                r["market_clv"],
                r["book_closing_line"], r["book_closing_odds"],
                r["book_closing_counter_odds"], r["book_closing_novig_prob"],
                r["book_clv"],
                r["line_moved"], r["clv_method"], r["sigma_used"],
            ])
        con.close()

    # Summary
    has_market = sum(1 for r in results if r["market_clv"] is not None)
    has_book = sum(1 for r in results if r["book_clv"] is not None)
    avg_mclv = (sum(r["market_clv"] for r in results if r["market_clv"] is not None) / has_market
                if has_market else 0)
    avg_bclv = (sum(r["book_clv"] for r in results if r["book_clv"] is not None) / has_book
                if has_book else 0)

    print(f"\n=== CLV RESULTS ===")
    print(f"  Bets processed:  {len(results)}")
    print(f"  Market CLV data: {has_market}/{len(results)} bets")
    print(f"  Book CLV data:   {has_book}/{len(results)} bets")
    if has_market:
        print(f"  Avg market CLV:  {avg_mclv:+.3f} ({avg_mclv*100:+.1f}%)")
    if has_book:
        print(f"  Avg book CLV:    {avg_bclv:+.3f} ({avg_bclv*100:+.1f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
