#!/usr/bin/env python3
"""
DraftKings MLB SGP Scraper (Pure REST API)

Fetches Same Game Parlay (SGP) odds from DraftKings for MLB spread+total combos.
No browser needed — uses curl_cffi with Chrome TLS impersonation to bypass Akamai.

Usage:
    cd mlb_sgp
    source venv/bin/activate
    python scraper_draftkings_sgp.py           # all games
    python scraper_draftkings_sgp.py --verbose  # show details
"""

import argparse
import sys
import time
import duckdb
from pathlib import Path
from curl_cffi import requests as cffi_requests

# Resolve repo root dynamically (works from worktrees too)
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
# If we're in a worktree, the repo root is the worktree root's parent's parent
# e.g., .worktrees/mlb-sgp/mlb_sgp -> .worktrees/mlb-sgp -> NFLWork
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))
_ANSWER_KEYS = _REPO_ROOT / "Answer Keys"

sys.path.insert(0, str(_ANSWER_KEYS))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

from db import ensure_table, upsert_sgp_odds, MLB_DB

# ---------------------------------------------------------------------------
# DK API config
# ---------------------------------------------------------------------------

DK_BASE_URL = "https://sportsbook.draftkings.com"

# Public REST (no Akamai)
DK_LEAGUE_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "controldata/league/leagueSubcategory/v1/markets"
)
DK_EVENT_MARKETS_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "controldata/event/eventSubcategory/v1/markets"
)

# SGP pricing (Akamai — bypassed by curl_cffi)
DK_CALCULATE_BETS_URL = (
    "https://gaming-us-nj.draftkings.com/api/wager/v1/calculateBets"
)

DK_MLB_LEAGUE_ID = "84240"

# Over suffix is always _1, Under is always _3.
# Spread suffix (_1 vs _3 for home/away) varies per game — discovered at runtime.
OVER_SUFFIX = "_1"
UNDER_SUFFIX = "_3"

# Track consecutive failures to detect session expiry
_consecutive_failures = 0
_MAX_CONSECUTIVE_FAILURES = 3


def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    else:
        return int(round(-100 / (dec - 1)))


# ---------------------------------------------------------------------------
# Step 1: Init session
# ---------------------------------------------------------------------------

def init_session() -> cffi_requests.Session:
    """Create a curl_cffi session with Chrome TLS impersonation and DK cookies."""
    session = cffi_requests.Session(impersonate="chrome")
    session.get(DK_BASE_URL + "/leagues/baseball/mlb", timeout=30)
    return session


# ---------------------------------------------------------------------------
# Step 2: Fetch DK events + market IDs (public REST, no Akamai)
# ---------------------------------------------------------------------------

def fetch_dk_events(session: cffi_requests.Session) -> list[dict]:
    """Fetch today's MLB events with home/away teams."""
    resp = session.get(DK_LEAGUE_URL, params={
        "isBatchable": "false",
        "templateVars": DK_MLB_LEAGUE_ID,
        "eventsQuery": (
            f"$filter=leagueId eq '{DK_MLB_LEAGUE_ID}' "
            f"AND clientMetadata/Subcategories/any(s: s/Id eq '4519')"
        ),
        "marketsQuery": (
            "$filter=clientMetadata/subCategoryId eq '4519' "
            "AND tags/all(t: t ne 'SportcastBetBuilder')"
        ),
        "include": "Events",
        "entity": "events",
    }, timeout=30)
    resp.raise_for_status()

    events = []
    for evt in resp.json().get("events", []):
        participants = evt.get("participants", [])
        home = next((p for p in participants if p.get("venueRole") == "Home"), {})
        away = next((p for p in participants if p.get("venueRole") == "Away"), {})
        events.append({
            "dk_event_id": evt["id"],
            "name": evt.get("name", ""),
            "dk_home": home.get("name", ""),
            "dk_away": away.get("name", ""),
            "start_time": evt.get("startEventDate", ""),
        })
    return events


def fetch_market_num(session: cffi_requests.Session, dk_event_id: str) -> str | None:
    """Fetch the Run Line market number for a DK event."""
    resp = session.get(DK_EVENT_MARKETS_URL, params={
        "isBatchable": "false",
        "templateVars": dk_event_id,
        "marketsQuery": (
            f"$filter=eventId eq '{dk_event_id}' "
            f"AND clientMetadata/subCategoryId eq '4519' "
            f"AND tags/all(t: t ne 'SportcastBetBuilder')"
        ),
        "include": "MarketSplits",
        "entity": "markets",
    }, timeout=15)

    if resp.status_code != 200:
        return None

    for m in resp.json().get("markets", []):
        if m.get("name") == "Run Line":
            return m["id"].split("_")[-1] if "_" in m["id"] else m["id"]
    return None


# ---------------------------------------------------------------------------
# Step 3: Load game lines from DuckDB
# ---------------------------------------------------------------------------

def load_parlay_lines() -> dict:
    """
    Load spread + total lines from mlb_parlay_opportunities.
    Returns {canonical_game_id: {spread_line, total_line, home_team, away_team}}
    """
    con = duckdb.connect(str(MLB_DB), read_only=True)
    try:
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_parlay_opportunities" not in tables:
            print("  No mlb_parlay_opportunities table — run the MLB pipeline first.")
            return {}

        rows = con.execute("""
            SELECT DISTINCT game_id, spread_line, total_line, home_team, away_team
            FROM mlb_parlay_opportunities
            WHERE combo = 'Home Spread + Over'
        """).fetchall()

        return {
            row[0]: {
                "spread_line": row[1],
                "total_line": row[2],
                "home_team": row[3],
                "away_team": row[4],
            }
            for row in rows
        }
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Step 4: Match DK events to our game_ids via canonical_match
# ---------------------------------------------------------------------------

def match_events(dk_events: list[dict], parlay_lines: dict) -> list[dict]:
    """Match DK events to our game_ids using canonical_match.py."""
    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")

    matched = []
    for dk_evt in dk_events:
        resolved = resolve_team_names(
            dk_evt["dk_away"], dk_evt["dk_home"],
            team_dict, canonical_games,
        )
        if not resolved or not resolved[0] or not resolved[1]:
            continue

        canon_away, canon_home = resolved

        for game_id, lines in parlay_lines.items():
            if lines["home_team"] == canon_home and lines["away_team"] == canon_away:
                matched.append({
                    "game_id": game_id,
                    "dk_event_id": dk_evt["dk_event_id"],
                    "home_team": canon_home,
                    "away_team": canon_away,
                    "dk_name": dk_evt["name"],
                    "spread_line": lines["spread_line"],
                    "total_line": lines["total_line"],
                })
                break

    return matched


# ---------------------------------------------------------------------------
# Step 5: Build selection IDs and call calculateBets
# ---------------------------------------------------------------------------

def build_selection_id(market_num: str, market_type: str, sign: str,
                       line_value: float, suffix: str) -> str:
    """Build a DK selection ID (e.g., 0HC84191347N150_1)."""
    encoded_line = str(int(abs(line_value) * 100))
    return f"0{market_type}{market_num}{sign}{encoded_line}{suffix}"


def calculate_sgp(session: cffi_requests.Session,
                  spread_sel: str, total_sel: str) -> dict | None:
    """
    Call DK's calculateBets endpoint with two selections.
    Returns {trueOdds, displayOdds, legs} or None.
    """
    resp = session.post(DK_CALCULATE_BETS_URL, json={
        "selections": [],
        "selectionsForYourBet": [
            {"id": spread_sel, "yourBetGroup": 0},
            {"id": total_sel, "yourBetGroup": 0},
        ],
        "selectionsForCombinator": [],
        "selectionsForProgressiveParlay": [],
        "oddsStyle": "american",
    }, headers={"Content-Type": "application/json"}, timeout=10)

    if resp.status_code != 200:
        return None

    data = resp.json()
    for bet in data.get("bets", []):
        mapped = bet.get("selectionsMapped", [])
        if bet.get("trueOdds") and len(mapped) >= 2:
            legs = []
            for s in data.get("selectionsForYourBet", []) + data.get("selections", []):
                if s.get("trueOdds"):
                    legs.append({
                        "id": s["id"],
                        "points": s.get("points"),
                        "odds": s["trueOdds"],
                        "display": s.get("displayOdds", ""),
                    })
            return {
                "trueOdds": bet["trueOdds"],
                "displayOdds": bet.get("displayOdds", ""),
                "legs": legs,
            }
    return None


def discover_suffixes(session: cffi_requests.Session, market_num: str,
                      home_sign: str, spread: float, total: float,
                      verbose: bool = False) -> tuple | None:
    """
    Discover which spread suffix and total line work for this game.

    Tries the most common pattern (_1/_3) first, then reverses.
    The Wagerzon total should always match DK since DK offers many totals,
    but we try ±0.5 as a safety net.

    Returns (home_suffix, away_suffix, actual_total) or None.
    """
    for home_suf, away_suf in [("_1", "_3"), ("_3", "_1")]:
        for total_offset in [0, 0.5, -0.5]:
            test_total = total + total_offset
            spread_sel = build_selection_id(market_num, "HC", home_sign, spread, home_suf)
            total_sel = build_selection_id(market_num, "OU", "O", test_total, OVER_SUFFIX)
            sgp = calculate_sgp(session, spread_sel, total_sel)
            if sgp:
                if verbose and total_offset != 0:
                    print(f"    (DK total is {test_total}, not {total})")
                if verbose:
                    print(f"    (Suffixes: home={home_suf}, away={away_suf})")
                return home_suf, away_suf, test_total

    return None


def scrape_game(session: cffi_requests.Session, game: dict,
                market_num: str, verbose: bool = False) -> list[dict]:
    """Scrape all 4 SGP combos for one game."""
    global _consecutive_failures
    results = []
    spread_line = game["spread_line"]
    total = game["total_line"]

    if spread_line < 0:
        home_sign, away_sign = "N", "P"
    else:
        home_sign, away_sign = "P", "N"

    spread = abs(spread_line)

    discovery = discover_suffixes(session, market_num, home_sign, spread, total, verbose)
    if not discovery:
        _consecutive_failures += 1
        if verbose:
            print("    Could not find valid suffix/total combination")
        return results

    _consecutive_failures = 0
    home_suf, away_suf, actual_total = discovery

    combos = [
        ("Home Spread + Over",  home_sign, home_suf, "O", OVER_SUFFIX),
        ("Home Spread + Under", home_sign, home_suf, "U", UNDER_SUFFIX),
        ("Away Spread + Over",  away_sign, away_suf, "O", OVER_SUFFIX),
        ("Away Spread + Under", away_sign, away_suf, "U", UNDER_SUFFIX),
    ]

    for combo_name, spread_sign, spread_suf, ou_sign, total_suf in combos:
        spread_sel = build_selection_id(market_num, "HC", spread_sign, spread, spread_suf)
        total_sel = build_selection_id(market_num, "OU", ou_sign, actual_total, total_suf)

        if verbose:
            print(f"    {combo_name}: {spread_sel} + {total_sel}")

        sgp = calculate_sgp(session, spread_sel, total_sel)

        if sgp:
            odds = sgp["trueOdds"]
            am = decimal_to_american(odds)
            total_note = f" [total={actual_total}]" if actual_total != total else ""
            print(f"    {combo_name}: {odds:.4f} ({am:+d}){total_note}")
            results.append({
                "combo_name": combo_name,
                "trueOdds": odds,
                "displayOdds": sgp["displayOdds"],
                "actual_total": actual_total,
            })
        elif verbose:
            print(f"    {combo_name}: no price")

        time.sleep(0.1)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scrape_dk_sgp(verbose: bool = False):
    """Main: fetch DK SGP odds for all MLB games via pure REST API."""
    global _consecutive_failures

    print("Loading parlay lines from DuckDB...")
    parlay_lines = load_parlay_lines()
    if not parlay_lines:
        return

    print(f"  {len(parlay_lines)} games with lines")

    print("Initializing DK session...")
    session = init_session()

    print("Fetching DraftKings MLB events...")
    dk_events = fetch_dk_events(session)
    print(f"  {len(dk_events)} DK events")

    print("Matching teams via canonical_match...")
    matched = match_events(dk_events, parlay_lines)
    print(f"  {len(matched)} matched games")

    if not matched:
        print("  No matches found.")
        return

    ensure_table()
    all_rows = []
    _consecutive_failures = 0

    for game in matched:
        # Re-init session if multiple consecutive games fail (session expiry)
        if _consecutive_failures >= _MAX_CONSECUTIVE_FAILURES:
            print("  Re-initializing session (possible expiry)...")
            session = init_session()
            _consecutive_failures = 0

        print(f"\n{'='*60}")
        print(f"  {game['away_team']} @ {game['home_team']}")
        print(f"  Spread: {game['spread_line']:+.1f} | Total: {game['total_line']}")
        print(f"{'='*60}")

        market_num = fetch_market_num(session, game["dk_event_id"])
        if not market_num:
            print("  No Run Line market — skipping")
            continue

        if verbose:
            print(f"  Market num: {market_num}")

        game_results = scrape_game(session, game, market_num, verbose)

        for gr in game_results:
            all_rows.append({
                "game_id": game["game_id"],
                "combo": gr["combo_name"],
                "period": "FG",
                "bookmaker": "draftkings",
                "sgp_decimal": round(gr["trueOdds"], 4),
                "sgp_american": decimal_to_american(gr["trueOdds"]),
                "source": "draftkings_direct",
            })

    # Batch write to DuckDB
    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} DK SGP odds to mlb_sgp_odds")
        print(f"{'='*60}")
    else:
        print("\nNo SGP odds collected.")

    return all_rows


def main():
    parser = argparse.ArgumentParser(description="DraftKings MLB SGP Scraper")
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    args = parser.parse_args()

    print("=" * 60)
    print("  DRAFTKINGS MLB SGP SCRAPER (Pure REST)")
    print("=" * 60)

    scrape_dk_sgp(verbose=args.verbose)


if __name__ == "__main__":
    main()
