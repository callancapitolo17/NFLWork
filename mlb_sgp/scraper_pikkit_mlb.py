#!/usr/bin/env python3
"""
MLB Pikkit Pro SGP Scraper

Fetches SGP (Same Game Parlay) odds from Pikkit Pro for MLB games.
Pikkit aggregates SGP prices from multiple books (FanDuel, DraftKings,
Novig, ProphetX) — so one scrape gives us multi-book pricing.

The scraper reads today's games from mlb_consensus_temp, finds each game
on Pikkit, builds spread+total parlays, and intercepts the /betslip API
response to extract correlation-adjusted odds from each book.

Results are stored in the mlb_sgp_odds table in mlb.duckdb.

Usage:
    cd mlb_sgp
    python scraper_pikkit_mlb.py               # headless, all games
    python scraper_pikkit_mlb.py --visible      # show browser (for debugging)
    python scraper_pikkit_mlb.py --game-id XYZ  # single game only
"""

import argparse
import sys
import duckdb
from pathlib import Path
from datetime import datetime

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    print("Playwright not installed. Run: pip install playwright && playwright install chromium")
    sys.exit(1)

from pikkit_common import (
    PIKKIT_URL,
    save_session,
    load_session,
    is_logged_in,
    clear_betslip,
    decimal_to_american,
    extract_parlay_odds_from_api,
    find_game_on_pikkit,
    get_sgp_odds_for_parlay,
)
from db import ensure_table, upsert_sgp_odds, MLB_DB

# ---------------------------------------------------------------------------
# Combo definitions — must match mlb_correlated_parlay.R exactly
# ---------------------------------------------------------------------------
# FG combos: "Home Spread + Over", "Home Spread + Under",
#            "Away Spread + Over", "Away Spread + Under"
# F5 combos: prefix with "F5 " (e.g., "F5 Home Spread + Over")

# Pikkit period names — Full Game is default for MLB.
# F5 (first 5 innings) support is TBD pending manual recon of Pikkit's MLB UI.
MLB_PERIODS = {
    "FG": "Full",
    # "F5": "1st 5 Innings",  # uncomment once confirmed on Pikkit
}


def load_todays_games(db_path: str = None) -> list[dict]:
    """
    Load today's MLB games from the consensus table.

    Returns list of dicts with:
        id, home_team, away_team, commence_time,
        home_spread, away_spread, total_line
    """
    db_path = db_path or str(MLB_DB)
    con = duckdb.connect(db_path, read_only=True)
    try:
        # Check that the required tables exist
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_consensus_temp" not in tables:
            print("No mlb_consensus_temp table — run the MLB pipeline first.")
            return []

        # Pull games with spread and total lines from consensus + Wagerzon
        # We need the spread/total to know what legs to build on Pikkit
        games = con.execute("""
            SELECT
                c.id,
                c.home_team,
                c.away_team,
                c.commence_time
            FROM mlb_consensus_temp c
            ORDER BY c.commence_time
        """).fetchdf()

        if games.empty:
            print("No games found in mlb_consensus_temp.")
            return []

        return games.to_dict("records")
    finally:
        con.close()


def load_parlay_lines(db_path: str = None) -> dict:
    """
    Load spread and total lines from mlb_parlay_opportunities.

    Returns dict keyed by game_id -> {spread_line, total_line, home_team, away_team}
    """
    db_path = db_path or str(MLB_DB)
    con = duckdb.connect(db_path, read_only=True)
    try:
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_parlay_opportunities" not in tables:
            print("No mlb_parlay_opportunities table — run mlb_correlated_parlay.R first.")
            return {}

        # Get distinct game lines (one row per game per period)
        lines = con.execute("""
            SELECT DISTINCT
                game_id,
                combo,
                spread_line,
                total_line
            FROM mlb_parlay_opportunities
            WHERE combo LIKE '%Home Spread + Over%'
        """).fetchdf()

        result = {}
        for _, row in lines.iterrows():
            gid = row["game_id"]
            period = "F5" if row["combo"].startswith("F5 ") else "FG"
            key = f"{gid}|{period}"
            result[key] = {
                "spread_line": row["spread_line"],
                "total_line": row["total_line"],
            }
        return result
    finally:
        con.close()


def build_combos(spread_line: float, total_line: float, period_code: str) -> list[dict]:
    """
    Build the 4 parlay combos for a game + period.

    Args:
        spread_line: Home team spread (e.g., -1.5 means home is favorite)
        total_line: Game total (e.g., 8.5)
        period_code: "FG" or "F5"

    Returns:
        List of combo dicts with: combo_name, side ("home"/"away"), spread_value,
        total_value, over_under ("Over"/"Under")
    """
    prefix = f"{period_code} " if period_code != "FG" else ""

    # spread_line is from home's perspective
    home_spread = spread_line
    away_spread = -spread_line

    return [
        {
            "combo_name": f"{prefix}Home Spread + Over",
            "side": "home",
            "spread_value": home_spread,
            "total_value": total_line,
            "over_under": "Over",
        },
        {
            "combo_name": f"{prefix}Home Spread + Under",
            "side": "home",
            "spread_value": home_spread,
            "total_value": total_line,
            "over_under": "Under",
        },
        {
            "combo_name": f"{prefix}Away Spread + Over",
            "side": "away",
            "spread_value": away_spread,
            "total_value": total_line,
            "over_under": "Over",
        },
        {
            "combo_name": f"{prefix}Away Spread + Under",
            "side": "away",
            "spread_value": away_spread,
            "total_value": total_line,
            "over_under": "Under",
        },
    ]


def scrape_mlb_sgp(headless: bool = True, game_id_filter: str = None,
                    trusted_books: list[str] = None):
    """
    Main scraping function. Opens Pikkit, finds MLB games, builds parlays,
    and stores SGP odds in DuckDB.

    Args:
        headless: Run browser headless (False = visible for debugging)
        game_id_filter: Only scrape this game_id (None = all games)
        trusted_books: Only keep odds from these books (None = keep all)
    """
    if trusted_books is None:
        trusted_books = ["fanduel", "draftkings", "novig", "prophetx"]

    # Load game data from DuckDB
    games = load_todays_games()
    if not games:
        return

    parlay_lines = load_parlay_lines()
    if not parlay_lines:
        return

    if game_id_filter:
        games = [g for g in games if g["id"] == game_id_filter]
        if not games:
            print(f"Game {game_id_filter} not found in consensus.")
            return

    print(f"Found {len(games)} MLB games to process")
    print(f"Trusted books: {', '.join(trusted_books)}")

    # Ensure the output table exists
    ensure_table()

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)

        # Load saved Pikkit session
        saved_session = load_session()
        if saved_session:
            context = browser.new_context(
                viewport={"width": 1920, "height": 1080},
                storage_state=saved_session,
            )
            print("Loaded saved Pikkit session")
        else:
            context = browser.new_context(viewport={"width": 1920, "height": 1080})
            print("No saved session — run with --visible to login manually")
            browser.close()
            return

        page = context.new_page()

        # Navigate to Pikkit and verify login
        page.goto(PIKKIT_URL, wait_until="networkidle")
        page.wait_for_timeout(3000)

        if not is_logged_in(page):
            if not headless:
                print("\n" + "=" * 60)
                print("NOT LOGGED IN — please login manually.")
                print("Session will be saved after 120 seconds.")
                print("=" * 60)
                page.wait_for_timeout(120000)
                save_session(context)
            else:
                print("Not logged in. Run with --visible to login first.")
                browser.close()
                return

        print("\nPikkit logged in!")
        all_rows = []

        for game in games:
            game_id = game["id"]
            home = game["home_team"]
            away = game["away_team"]

            print(f"\n{'='*60}")
            print(f"  {away} @ {home}")
            print(f"  game_id: {game_id}")
            print(f"{'='*60}")

            # Find the game on Pikkit
            game_url = find_game_on_pikkit(page, away, home, sport="MLB")
            if not game_url:
                print("  Could not find game on Pikkit — skipping")
                continue

            # Navigate to the game page
            page.goto(game_url, wait_until="networkidle")
            page.wait_for_timeout(2000)

            # Process each period (FG, maybe F5)
            for period_code, pikkit_period in MLB_PERIODS.items():
                key = f"{game_id}|{period_code}"
                lines = parlay_lines.get(key)
                if not lines:
                    print(f"  No {period_code} lines found for this game — skipping")
                    continue

                spread_line = lines["spread_line"]
                total_line = lines["total_line"]

                print(f"\n  {period_code}: spread={spread_line:+.1f}, total={total_line}")

                combos = build_combos(spread_line, total_line, period_code)

                for combo in combos:
                    # Determine which team name to use for the spread leg
                    spread_team = home if combo["side"] == "home" else away

                    print(f"\n    {combo['combo_name']}")
                    print(f"    {spread_team} {combo['spread_value']:+.1f} + {combo['over_under']} {combo['total_value']}")

                    # Build the parlay on Pikkit and get SGP odds
                    sgp = get_sgp_odds_for_parlay(
                        page,
                        period=pikkit_period,
                        spread_team=spread_team,
                        spread_value=combo["spread_value"],
                        total_value=combo["total_value"],
                        over_under=combo["over_under"],
                        sport="MLB",
                    )

                    if not sgp.get("book_odds"):
                        print("    No SGP odds returned")
                        continue

                    # Convert book_odds to DB rows
                    for bo in sgp["book_odds"]:
                        book_name = bo["book"].lower().strip()

                        # Filter to trusted books
                        if not any(tb in book_name for tb in trusted_books):
                            continue

                        american = bo["odds"]
                        # Convert American back to decimal for storage
                        if american > 0:
                            dec_odds = 1 + american / 100
                        else:
                            dec_odds = 1 + 100 / abs(american)

                        all_rows.append({
                            "game_id": game_id,
                            "combo": combo["combo_name"],
                            "period": period_code,
                            "bookmaker": book_name,
                            "sgp_decimal": round(dec_odds, 4),
                            "sgp_american": american,
                            "source": "pikkit",
                        })

                    # Print summary for this combo
                    books_str = ", ".join(
                        f"{bo['book']}:{bo['odds']:+d}" for bo in sgp["book_odds"]
                        if any(tb in bo["book"].lower() for tb in trusted_books)
                    )
                    if sgp.get("best_odds"):
                        print(f"    Best: {sgp['best_odds']:+d} | Books: {books_str}")

                    # Re-navigate to game page before next combo (resets betslip)
                    page.goto(game_url, wait_until="networkidle")
                    page.wait_for_timeout(1500)

        browser.close()

    # Write all collected rows to DuckDB
    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} SGP odds rows to mlb_sgp_odds")
        print(f"{'='*60}")
    else:
        print("\nNo SGP odds collected.")


def login_only():
    """Just open Pikkit for login and save session. No DB reads needed."""
    from playwright.sync_api import sync_playwright as sp

    print("Opening Pikkit for login...")
    print("Login with your phone number + SMS code.")
    print("Session will save automatically after 120 seconds.\n")

    with sp() as p:
        browser = p.chromium.launch(headless=False)
        saved = load_session()
        if saved:
            context = browser.new_context(
                viewport={"width": 1920, "height": 1080},
                storage_state=saved,
            )
            print("Loaded existing session — checking if still valid...")
        else:
            context = browser.new_context(viewport={"width": 1920, "height": 1080})
            print("No existing session found.")

        page = context.new_page()
        page.goto(PIKKIT_URL, wait_until="networkidle")
        page.wait_for_timeout(3000)

        if is_logged_in(page):
            print("Already logged in! Session is valid.")
            # Re-save to refresh the session file
            save_session(context)
        else:
            print("\nNot logged in — please login manually in the browser.")
            print("Waiting 120 seconds for you to complete SMS login...")
            page.wait_for_timeout(120000)
            save_session(context)
            if is_logged_in(page):
                print("Login successful! Session saved.")
            else:
                print("Still not logged in. Try again.")

        browser.close()


def main():
    parser = argparse.ArgumentParser(description="MLB Pikkit Pro SGP Scraper")
    parser.add_argument("--visible", action="store_true", help="Show browser window")
    parser.add_argument("--login", action="store_true", help="Just login and save session (no scraping)")
    parser.add_argument("--game-id", type=str, help="Only scrape this game_id")
    parser.add_argument(
        "--books",
        nargs="+",
        default=["fanduel", "draftkings", "novig", "prophetx"],
        help="Trusted books to keep (default: fanduel draftkings novig prophetx)",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  MLB PIKKIT SGP SCRAPER")
    print("=" * 60)

    if args.login:
        login_only()
        return

    scrape_mlb_sgp(
        headless=not args.visible,
        game_id_filter=args.game_id,
        trusted_books=args.books,
    )


if __name__ == "__main__":
    main()
