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
                     amount: int = 10000) -> dict | None:
    """Call ConfirmWagerHelper to get exact parlay payout.

    Args:
        session: Authenticated requests session
        idgm: Wagerzon internal game ID
        legs: List of dicts with {play, points, odds}
        amount: Bet amount for price query. Default $10000 maximises decimal
            precision (±$0.005 band vs ±$0.50 at $100). MAXPARLAYRISKEXCEED at
            $10000 triggers a fallback to $100. Pass a specific amount when
            you need the exact integer payout at that stake (e.g. Stage 2
            empirical nudge).

    Returns:
        Dict with {win, decimal, american, amount} or None on error
    """
    sel = ",".join(
        f"{l['play']}_{idgm}_{l['points']}_{l['odds']}" for l in legs
    )
    # RiskWin="2" marks this as a preview quote. With RiskWin=0 (real-bet mode),
    # WZ balance-checks and returns BALANCEEXCEED whenever account balance < amount.
    # With RiskWin="2" the balance check is skipped and we get a price at any amount.
    detail_data = [
        {
            "Amount": str(amount),
            "RiskWin": "2",
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

        details = result.get("details", [])
        if not details:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        # In preview mode (RiskWin="2"), ErrorMsgKey can fire with useful data still
        # populated — e.g. MINWAGERONLINE returns the correct Risk/Win for the queried
        # amount but flags that it wouldn't be placeable as a real bet. Only treat the
        # response as failed when Risk and Win are both zero (e.g. MAXPARLAYRISKEXCEED).
        win = details[0].get("Win", 0)
        risk = details[0].get("Risk", 0)
        if win <= 0 or risk <= 0:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        decimal_odds = 1 + win / amount
        american = round(win / amount * 100) if win > amount else round(-amount / win * 100)

        return {
            "win": win,
            "decimal": round(decimal_odds, 4),
            "american": american,
            "amount": amount,
        }

    except Exception as e:
        print(f"  Request error: {e}")
        return None


def get_combined_parlay_price(session: requests.Session, legs: list[dict],
                               amount: int = 10000) -> dict | None:
    """Call ConfirmWagerHelper for a cross-game parlay (legs from multiple games).

    Differs from get_parlay_price() in that each leg carries its own idgm.

    Args:
        session: Authenticated requests session
        legs: List of dicts with {idgm, play, points, odds} per leg
        amount: Bet amount for price query (defaults to 10000 for precision;
                falls back to 100 if MAXPARLAYRISKEXCEED — caller's responsibility)

    Returns:
        Dict with {win, decimal, american, amount} or None on error
    """
    sel = ",".join(
        f"{l['play']}_{l['idgm']}_{l['points']}_{l['odds']}" for l in legs
    )
    detail_data = [
        {
            "Amount": str(amount),
            "RiskWin": "2",
            "TeaserPointsPurchased": 0,
            "IdGame": l["idgm"],
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

        details = result.get("details", [])
        if not details:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        win = details[0].get("Win", 0)
        risk = details[0].get("Risk", 0)
        if win <= 0 or risk <= 0:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        decimal_odds = 1 + win / amount
        american = round(win / amount * 100) if win > amount else round(-amount / win * 100)

        return {
            "win": win,
            "decimal": round(decimal_odds, 4),
            "american": american,
            "amount": amount,
        }
    except Exception as e:
        print(f"  Request error: {e}")
        return None


def get_parlay_price_with_fallback(session: requests.Session, idgm: int,
                                    legs: list[dict]) -> dict | None:
    """Stage 1 price lookup with descending amount ladder.

    Tries $10000 first for maximum decimal precision. If WZ returns
    MAXPARLAYRISKEXCEED (max parlay risk for this combo is below $10000), falls
    back to $100 which is always inside WZ's per-parlay max.

    Use this for initial pricing (populating mlb_parlay_prices). For exact
    payout at a specific stake, call get_parlay_price(..., amount=stake).
    """
    for amount in (10000, 100):
        result = get_parlay_price(session, idgm, legs, amount=amount)
        if result is not None:
            return result
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
            price = get_parlay_price_with_fallback(session, idgm, c["legs"])
            if price:
                results.append((
                    fetch_time, g["home_team"], g["away_team"], idgm,
                    c["combo"], period, price["decimal"], price["american"], price["win"]
                ))
                print(f"  {game_label} | {c['combo']}: +{price['american']} "
                      f"(${price['win']} on ${price['amount']})")
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


# ─────────────────────────────────────────────────────────────────────────────
# STAGE 2: Empirical nudge + exact payout
# ─────────────────────────────────────────────────────────────────────────────


NUDGE_RANGE = 2  # query W-2..W+2 around Kelly-ideal wager

# Path to the answer-keys DuckDB that owns mlb_parlay_opportunities.
# parlay_pricer.py lives in wagerzon_odds/; the answer-keys DB lives in
# Answer Keys/mlb.duckdb relative to the repo root.
_MLB_DB_PATH = Path(__file__).resolve().parent.parent / "Answer Keys" / "mlb.duckdb"


def _ensure_exact_columns(conn):
    """Ensure optional exact-payout columns exist on mlb_parlay_opportunities.

    DuckDB lacks `ADD COLUMN IF NOT EXISTS`, so we probe information_schema.
    Columns added here are nullable so prior rows (written by an older
    mlb_correlated_parlay.R) continue to load cleanly.
    """
    existing = {
        r[0] for r in conn.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name = 'mlb_parlay_opportunities'"
        ).fetchall()
    }
    if "exact_to_win" not in existing:
        conn.execute("ALTER TABLE mlb_parlay_opportunities ADD COLUMN exact_to_win INTEGER")
    if "exact_wager" not in existing:
        conn.execute("ALTER TABLE mlb_parlay_opportunities ADD COLUMN exact_wager INTEGER")


def _combo_to_legs(combo: str, idgm: int, spread_line: float, total_line: float,
                    spread_price: int, total_price: int) -> list[dict]:
    """Reconstruct ConfirmWagerHelper legs from the combo row."""
    is_home = "Home" in combo
    is_over = "Over" in combo
    # Spread points: away uses raw spread_line (positive dog), home uses signed spread_line.
    # mlb_correlated_parlay.R stores spread_line as the signed number for the picked side.
    spread_leg = {
        "play": PLAY_HOME_SPREAD if is_home else PLAY_AWAY_SPREAD,
        "points": spread_line,
        "odds": spread_price,
    }
    # Total points: WZ expects -total for over, +total for under (recon-confirmed).
    total_leg = {
        "play": PLAY_OVER if is_over else PLAY_UNDER,
        "points": -abs(total_line) if is_over else total_line,
        "odds": total_price,
    }
    return [spread_leg, total_leg]


def compute_exact_payouts(session: requests.Session):
    """Stage 2: empirical nudge + exact payout per sized parlay.

    For each row in mlb_parlay_opportunities with kelly_bet > 0:
        1. Query ConfirmWagerHelper at stake ∈ [kelly_bet - 2, kelly_bet + 2]
        2. Pick the stake that maximises WZ's returned win / stake ratio
        3. Write (exact_wager, exact_to_win) back to the row

    Replaces the math-based nudge in mlb_correlated_parlay.R. Uses WZ's real
    rounding behaviour instead of predicting it from our stored decimal.
    """
    conn = duckdb.connect(str(_MLB_DB_PATH))
    try:
        _ensure_exact_columns(conn)

        sized = conn.execute("""
            SELECT parlay_hash, game, combo, idgm,
                   spread_line, total_line, spread_price, total_price,
                   CAST(kelly_bet AS INTEGER) AS kelly_bet
            FROM mlb_parlay_opportunities
            WHERE kelly_bet > 0 AND idgm IS NOT NULL
        """).fetchall()

        if not sized:
            print("No sized parlays with idgm — nothing to price.")
            return

        print(f"Pricing {len(sized)} sized parlays at exact stakes...")
        updates = []
        for (parlay_hash, game, combo, idgm,
             spread_line, total_line, spread_price, total_price,
             kelly_bet) in sized:
            legs = _combo_to_legs(combo, idgm, spread_line, total_line,
                                  spread_price, total_price)

            candidates = []
            for stake in range(max(1, kelly_bet - NUDGE_RANGE),
                               kelly_bet + NUDGE_RANGE + 1):
                price = get_parlay_price(session, idgm, legs, amount=stake)
                if price is None:
                    continue  # MAXPARLAYRISKEXCEED etc.
                candidates.append((stake, int(price["win"])))

            if not candidates:
                print(f"  {game} | {combo}: no valid candidate stakes — skipping")
                continue

            # Pick the wager with the highest win/stake ratio; break ties by
            # taking the stake closest to Kelly-ideal to keep growth on track.
            best = max(candidates,
                       key=lambda sw: (sw[1] / sw[0], -abs(sw[0] - kelly_bet)))
            best_stake, best_win = best
            updates.append((best_stake, best_win, parlay_hash))
            print(f"  {game} | {combo}: Kelly={kelly_bet} → "
                  f"nudged={best_stake} (to_win=${best_win})")

        if not updates:
            print("No updates to apply.")
            return

        conn.executemany(
            "UPDATE mlb_parlay_opportunities "
            "SET exact_wager = ?, exact_to_win = ? "
            "WHERE parlay_hash = ?",
            updates,
        )
        print(f"Saved {len(updates)} exact payouts to mlb_parlay_opportunities.")
    finally:
        conn.close()


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = "price"
    sport = "mlb"
    for a in args:
        if a == "--exact-payouts":
            mode = "exact"
        elif not a.startswith("--"):
            sport = a

    if sport != "mlb":
        print(f"Parlay pricer currently only supports MLB (got: {sport})")
        sys.exit(1)

    session = requests.Session()
    print("Logging in to Wagerzon...")
    login(session)

    if mode == "exact":
        compute_exact_payouts(session)
    else:
        price_mlb_parlays(session)
