#!/usr/bin/env python3
"""
Kalshi Orderbook Scanner — Find mispriced contracts vs fair value.

Loads predictions from the answer key, fetches live Kalshi markets,
and prints a sorted table of +EV opportunities above a threshold.

Usage:
    python scan_book.py              # Default 5% EV threshold
    python scan_book.py --min-ev 3   # Custom threshold (percent)
"""

import sys
from pathlib import Path

# Add kalshi_mm to path for imports
sys.path.insert(0, str(Path(__file__).parent))

import config
import db

# Add kalshi_odds to path for market fetching
sys.path.insert(0, str(config.PROJECT_ROOT / "kalshi_odds"))
from scraper import fetch_markets, parse_spread_team

# Add Answer Keys for team name resolution
sys.path.insert(0, str(config.ANSWER_KEYS_DIR))
from canonical_match import load_team_dict, load_canonical_games

from main import match_kalshi_to_predictions, fetch_all_markets, match_all_markets


def scan(min_ev_pct=5.0):
    """Scan orderbook for mispriced contracts and print opportunities."""

    # 1. Load predictions
    predictions, updated_at = db.load_predictions()
    if not predictions:
        print("No predictions available. Run the CBB pipeline first:")
        print('  cd "Answer Keys" && python run.py --sport cbb')
        return

    age_str = "unknown"
    if updated_at:
        from datetime import datetime, timezone
        ts = updated_at
        if ts.tzinfo is None:
            ts = ts.replace(tzinfo=timezone.utc)
        age_sec = (datetime.now(timezone.utc) - ts).total_seconds()
        age_str = f"{age_sec / 60:.0f}m ago"

    print(f"Predictions: {len(predictions)} loaded ({age_str})")

    # 2. Fetch Kalshi markets (all enabled types)
    enabled = ", ".join(sorted(config.ENABLED_MARKET_TYPES))
    print(f"Fetching Kalshi markets (enabled: {enabled})...")
    all_kalshi = fetch_all_markets()
    if not all_kalshi:
        print("No open Kalshi markets found.")
        return
    total_contracts = sum(len(v) for v in all_kalshi.values())
    print(f"  {total_contracts} open contracts")

    # 3. Match to predictions
    team_dict = load_team_dict("cbb")
    canonical_games = load_canonical_games("cbb")
    quotable = match_all_markets(all_kalshi, predictions, team_dict, canonical_games)
    print(f"  {len(quotable)} matched to predictions\n")

    if not quotable:
        print("No markets matched predictions. Check team name resolution.")
        return

    # 4. Calculate EV for each opportunity
    opportunities = []
    for m in quotable:
        fair_cents = m["fair_prob"] * 100
        yes_bid = m["book_bid"]
        yes_ask = m["book_ask"]
        mtype = m.get("market_type", "spreads")

        # Format contract description based on market type
        if mtype == "totals":
            contract_desc = f"O/U {m['strike']}"
        elif mtype == "moneyline":
            contract_desc = f"ML {m.get('contract_team', '?')}"
        else:
            contract_desc = f"{m.get('contract_team', '?')} by >{m.get('strike', '?')}"

        # Buy YES at yes_ask
        if yes_ask > 0 and fair_cents > yes_ask:
            ev_pct = (fair_cents - yes_ask) / yes_ask * 100
            if ev_pct >= min_ev_pct:
                opportunities.append({
                    "ticker": m["ticker"],
                    "game": f"{m['away_team']} @ {m['home_team']}",
                    "contract": contract_desc,
                    "market_type": mtype,
                    "fair_prob": m["fair_prob"],
                    "yes_bid": yes_bid,
                    "yes_ask": yes_ask,
                    "side": "BUY YES",
                    "price": yes_ask,
                    "ev_pct": ev_pct,
                })

        # Buy NO at (100 - yes_bid)
        if yes_bid > 0:
            no_price = 100 - yes_bid
            fair_no = 100 - fair_cents
            if fair_no > no_price and no_price > 0:
                ev_pct = (fair_no - no_price) / no_price * 100
                if ev_pct >= min_ev_pct:
                    opportunities.append({
                        "ticker": m["ticker"],
                        "game": f"{m['away_team']} @ {m['home_team']}",
                        "contract": contract_desc,
                        "market_type": mtype,
                        "fair_prob": m["fair_prob"],
                        "yes_bid": yes_bid,
                        "yes_ask": yes_ask,
                        "side": "BUY NO",
                        "price": no_price,
                        "ev_pct": ev_pct,
                    })

    # 5. Sort and print
    opportunities.sort(key=lambda x: x["ev_pct"], reverse=True)

    if not opportunities:
        print(f"No opportunities above {min_ev_pct:.0f}% EV found.")
        print(f"  ({len(quotable)} markets scanned)")
        return

    # Header
    print(f"{'TICKER':<35} {'TYPE':<5} {'FAIR':>5} {'BID':>4} {'ASK':>4} {'SIDE':<8} {'PRICE':>5} {'EV%':>7}")
    print("-" * 85)

    for o in opportunities:
        tabbr = o['market_type'][:3].upper()
        print(f"{o['ticker']:<35} {tabbr:<5} {o['fair_prob']:>5.1%} {o['yes_bid']:>4} {o['yes_ask']:>4} "
              f"{o['side']:<8} {o['price']:>5} {o['ev_pct']:>6.1f}%")

    print(f"\n{len(opportunities)} opportunities found (min EV: {min_ev_pct:.0f}%)")
    print(f"Scanned {len(quotable)} matched markets")


if __name__ == "__main__":
    threshold = 5.0
    for i, arg in enumerate(sys.argv):
        if arg == "--min-ev" and i + 1 < len(sys.argv):
            threshold = float(sys.argv[i + 1])
    scan(min_ev_pct=threshold)
