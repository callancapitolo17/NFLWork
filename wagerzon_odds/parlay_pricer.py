#!/usr/bin/env python3
"""
Wagerzon Parlay Pricer
Fetches exact parlay payouts from Wagerzon's ConfirmWagerHelper API.

For each MLB game with FG spread + total odds, builds 4 parlay combos
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
    """Price all MLB FG spread+total parlay combos."""
    conn = duckdb.connect(str(DB_PATH))

    # Read FG spreads with idgm
    games = conn.execute("""
        SELECT DISTINCT home_team, away_team, idgm,
               away_spread, away_spread_price, home_spread, home_spread_price,
               total, over_price, under_price
        FROM mlb_odds
        WHERE period = 'fg' AND market = 'spreads'
          AND idgm IS NOT NULL
          AND away_spread IS NOT NULL AND total IS NOT NULL
    """).fetchall()

    cols = ["home_team", "away_team", "idgm",
            "away_spread", "away_spread_price", "home_spread", "home_spread_price",
            "total", "over_price", "under_price"]

    if not games:
        print("No MLB FG games with idgm found. Run scraper first.")
        conn.close()
        return

    print(f"Pricing parlays for {len(games)} MLB games...")

    # Build 4 combos per game
    results = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    for row in games:
        g = dict(zip(cols, row))
        idgm = g["idgm"]
        game_label = f"{g['away_team']} @ {g['home_team']}"

        combos = [
            {
                "combo": "Home Spread + Over",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": abs(g["home_spread_price"])},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": abs(g["over_price"])},
                ],
            },
            {
                "combo": "Home Spread + Under",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": abs(g["home_spread_price"])},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": abs(g["under_price"])},
                ],
            },
            {
                "combo": "Away Spread + Over",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": abs(g["away_spread_price"])},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": abs(g["over_price"])},
                ],
            },
            {
                "combo": "Away Spread + Under",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": abs(g["away_spread_price"])},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": abs(g["under_price"])},
                ],
            },
        ]

        for c in combos:
            price = get_parlay_price(session, idgm, c["legs"])
            if price:
                results.append((
                    fetch_time, g["home_team"], g["away_team"], idgm,
                    c["combo"], price["decimal"], price["american"], price["win"]
                ))
                print(f"  {game_label} | {c['combo']}: +{price['american']} (${price['win']} on $100)")
            else:
                print(f"  {game_label} | {c['combo']}: FAILED")

    # Save to DuckDB
    if results:
        conn.execute("DROP TABLE IF EXISTS mlb_parlay_prices")
        conn.execute("""
            CREATE TABLE mlb_parlay_prices (
                fetch_time TIMESTAMP,
                home_team VARCHAR,
                away_team VARCHAR,
                idgm INTEGER,
                combo VARCHAR,
                wz_decimal DOUBLE,
                wz_american INTEGER,
                wz_win DOUBLE
            )
        """)
        conn.executemany(
            "INSERT INTO mlb_parlay_prices VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            results,
        )
        print(f"\nSaved {len(results)} parlay prices to mlb_parlay_prices")
    else:
        print("No prices fetched")

    conn.close()


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "mlb"

    if sport != "mlb":
        print(f"Parlay pricer currently only supports MLB (got: {sport})")
        sys.exit(1)

    session = requests.Session()
    print("Logging in to Wagerzon...")
    login(session)

    price_mlb_parlays(session)
