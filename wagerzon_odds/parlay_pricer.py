#!/usr/bin/env python3
"""
Wagerzon Parlay Pricer
Fetches exact parlay payouts from Wagerzon's ConfirmWagerHelper API.

For each MLB game with spread + total odds (FG and F5), builds 4 parlay combos
(home spread+over, home spread+under, away spread+over, away spread+under)
and queries the real payout. Stores results in DuckDB.

Usage:
    python3 wagerzon_odds/parlay_pricer.py [sport]
    python3 wagerzon_odds/parlay_pricer.py mlb  (default)
"""

import json
import sys
import duckdb
from datetime import datetime, timezone
from pathlib import Path

from scraper_v2 import login, DB_PATH
from config import WAGERZON_BASE_URL
import requests

CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"

# Play codes (from Wagerzon React bundle main.db15c074.js):
#   0 = away spread, 1 = home spread
#   2 = over total,  3 = under total
#   4 = away ML,     5 = home ML
PLAY_AWAY_SPREAD = 0
PLAY_HOME_SPREAD = 1
PLAY_OVER = 2
PLAY_UNDER = 3


def get_parlay_price(session: requests.Session, idgm: int, legs: list[dict],
                     amount: int = 100) -> dict | None:
    """Call ConfirmWagerHelper to get exact parlay payout.

    Args:
        session: Authenticated requests session
        idgm: Wagerzon internal game ID
        legs: List of dicts with {play, points, odds}
        amount: Bet amount for price query (default $100)

    Returns:
        Dict with {win, decimal, american} or None on error
    """
    sel = ",".join(
        f"{l['play']}_{idgm}_{l['points']}_{l['odds']}" for l in legs
    )
    detail_data = [
        {
            "Amount": str(amount),
            "RiskWin": 0,
            "TeaserPointsPurchased": 0,
            "IdGame": idgm,
            "Play": l["play"],
            "Pitcher": 3,
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
        for l in legs
    ]

    try:
        resp = session.post(
            CONFIRM_URL,
            data={
                "IDWT": "0",
                "WT": "1",
                "amountType": "0",
                "open": "0",
                "sameAmount": "false",
                "sameAmountNumber": "0",
                "useFreePlayAmount": "false",
                "sel": sel,
                "detailData": json.dumps(detail_data),
            },
            timeout=15,
        )
        resp.raise_for_status()
        result = resp.json().get("result", {})

        # Check for errors
        if result.get("ErrorMsgKey"):
            print(f"  API error: {result.get('ErrorMsg', result['ErrorMsgKey'])}")
            return None

        details = result.get("details", [])
        if not details:
            return None

        win = details[0].get("Win", 0)
        if win <= 0:
            return None

        decimal_odds = 1 + win / amount
        american = round(win / amount * 100) if win > amount else round(-amount / win * 100)

        return {"win": win, "decimal": round(decimal_odds, 4), "american": american}

    except Exception as e:
        print(f"  Request error: {e}")
        return None


def price_mlb_parlays(session: requests.Session):
    """Price all MLB FG + F5 spread+total parlay combos."""
    conn = duckdb.connect(str(DB_PATH))
    try:
        fg_results = _price_mlb_parlays_inner(session, conn, period="fg")
        f5_results = _price_mlb_parlays_inner(session, conn, period="f5")
        all_results = fg_results + f5_results
        _save_parlay_prices(conn, all_results)
    finally:
        conn.close()


def _price_mlb_parlays_inner(session: requests.Session, conn, period: str = "fg"):
    """Price parlay combos for a given period. Returns list of result tuples.

    Args:
        session: Authenticated requests session
        conn: DuckDB connection
        period: "fg" for full game, "f5" for first 5 innings

    Returns:
        List of tuples: (fetch_time, home, away, idgm, combo, period, decimal, american, win)
    """
    # Period-specific query filters and combo name prefix
    if period == "f5":
        where_clause = "period = 'h1' AND market = 'spreads_h1'"
        combo_prefix = "F5 "
        label = "F5"
    else:
        where_clause = "period = 'fg' AND market = 'spreads'"
        combo_prefix = ""
        label = "FG"

    games = conn.execute(f"""
        SELECT DISTINCT home_team, away_team, idgm,
               away_spread, away_spread_price, home_spread, home_spread_price,
               total, over_price, under_price
        FROM mlb_odds
        WHERE {where_clause}
          AND idgm IS NOT NULL
          AND away_spread IS NOT NULL AND total IS NOT NULL
    """).fetchall()

    cols = ["home_team", "away_team", "idgm",
            "away_spread", "away_spread_price", "home_spread", "home_spread_price",
            "total", "over_price", "under_price"]

    if not games:
        print(f"No MLB {label} games with idgm found. Run scraper first.")
        return []

    print(f"\nPricing {label} parlays for {len(games)} MLB games...")

    results = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    for row in games:
        g = dict(zip(cols, row))
        idgm = g["idgm"]
        game_label = f"{g['away_team']} @ {g['home_team']}"

        combos = [
            {
                "combo": f"{combo_prefix}Home Spread + Over",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": g["home_spread_price"]},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": g["over_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Home Spread + Under",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": g["home_spread_price"]},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": g["under_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Away Spread + Over",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": g["away_spread_price"]},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": g["over_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Away Spread + Under",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": g["away_spread_price"]},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": g["under_price"]},
                ],
            },
        ]

        for c in combos:
            price = get_parlay_price(session, idgm, c["legs"])
            if price:
                results.append((
                    fetch_time, g["home_team"], g["away_team"], idgm,
                    c["combo"], period, price["decimal"], price["american"], price["win"]
                ))
                print(f"  {game_label} | {c['combo']}: +{price['american']} (${price['win']} on $100)")
            else:
                print(f"  {game_label} | {c['combo']}: FAILED")

    print(f"{label}: {len(results)} prices fetched")
    return results


def _save_parlay_prices(conn, results: list):
    """Save combined FG + F5 parlay prices to DuckDB."""
    if not results:
        print("No prices fetched")
        return

    conn.execute("DROP TABLE IF EXISTS mlb_parlay_prices")
    conn.execute("""
        CREATE TABLE mlb_parlay_prices (
            fetch_time TIMESTAMP,
            home_team VARCHAR,
            away_team VARCHAR,
            idgm INTEGER,
            combo VARCHAR,
            period VARCHAR,
            wz_decimal DOUBLE,
            wz_american INTEGER,
            wz_win DOUBLE
        )
    """)
    conn.executemany(
        "INSERT INTO mlb_parlay_prices VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        results,
    )
    print(f"\nSaved {len(results)} parlay prices to mlb_parlay_prices")


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "mlb"

    if sport != "mlb":
        print(f"Parlay pricer currently only supports MLB (got: {sport})")
        sys.exit(1)

    session = requests.Session()
    print("Logging in to Wagerzon...")
    login(session)

    price_mlb_parlays(session)
