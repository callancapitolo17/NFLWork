#!/usr/bin/env python3
"""
DraftKings MLB SGP Scraper (Pure REST API)

Fetches Same Game Parlay (SGP) odds from DraftKings for MLB spread+total combos.
No browser needed — uses curl_cffi with Chrome TLS impersonation to bypass Akamai.

How it works:
1. Fetch DK events via public REST API
2. For each game, fetch the SGP parlays data (curl_cffi) — returns ALL selection IDs
3. Look up the exact selection IDs for our spread + total lines
4. Call calculateBets with those IDs to get correlation-adjusted SGP price

Usage:
    cd mlb_sgp
    source venv/bin/activate
    python scraper_draftkings_sgp.py           # all games
    python scraper_draftkings_sgp.py --verbose  # show details
"""

import argparse
import re
import sys
import time
import duckdb
from pathlib import Path
from curl_cffi import requests as cffi_requests

# Resolve repo root dynamically (works from worktrees too)
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
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

# Public REST — event listing (no Akamai)
DK_LEAGUE_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "controldata/league/leagueSubcategory/v1/markets"
)

# SGP parlays — full selection ID listing (curl_cffi bypasses Akamai)
DK_SGP_PARLAYS_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "parlays/v1/sgp/events"
)

# SGP pricing — correlation-adjusted odds (curl_cffi bypasses Akamai)
DK_CALCULATE_BETS_URL = (
    "https://gaming-us-nj.draftkings.com/api/wager/v1/calculateBets"
)

DK_MLB_LEAGUE_ID = "84240"

MAX_CONSECUTIVE_FAILURES = 3


def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    else:
        return int(round(-100 / (dec - 1)))


# ---------------------------------------------------------------------------
# Step 1: Init session
# ---------------------------------------------------------------------------

def init_session() -> cffi_requests.Session:
    """Create a curl_cffi session with Chrome TLS impersonation."""
    session = cffi_requests.Session(impersonate="chrome")
    session.get(DK_BASE_URL + "/leagues/baseball/mlb", timeout=30)
    return session


# ---------------------------------------------------------------------------
# Step 2: Fetch DK events (public REST)
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


# ---------------------------------------------------------------------------
# Step 3: Load game lines from DuckDB
# ---------------------------------------------------------------------------

def load_parlay_lines() -> dict:
    """Load spread + total lines from mlb_parlay_opportunities."""
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
# Step 4: Match DK events to our game_ids
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
# Step 5: Fetch exact selection IDs from SGP parlays endpoint
# ---------------------------------------------------------------------------

def fetch_main_market_num(session: cffi_requests.Session, dk_event_id: str) -> str | None:
    """Fetch the main Run Line market number from the public API (subcategory 4519)."""
    resp = session.get(
        "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
        "controldata/event/eventSubcategory/v1/markets",
        params={
            "isBatchable": "false",
            "templateVars": dk_event_id,
            "marketsQuery": (
                f"$filter=eventId eq '{dk_event_id}' "
                f"AND clientMetadata/subCategoryId eq '4519' "
                f"AND tags/all(t: t ne 'SportcastBetBuilder')"
            ),
            "include": "MarketSplits",
            "entity": "markets",
        },
        timeout=15,
    )
    if resp.status_code != 200:
        return None
    for m in resp.json().get("markets", []):
        if m.get("name") == "Run Line":
            return m["id"].split("_")[-1] if "_" in m["id"] else m["id"]
    return None


def fetch_selection_ids(session: cffi_requests.Session, dk_event_id: str,
                        main_market_num: str | None = None,
                        verbose: bool = False) -> dict:
    """
    Fetch all SGP selection IDs for a game from the parlays endpoint.

    The 2MB+ response contains selection IDs for every market (main, alt,
    innings, props). We parse it to find the IDs for spreads and totals.

    For spreads: we prefer IDs from the main market (full game Run Line).
    For totals: the main market only has one total, so we also collect
    alt totals from other markets to match any Wagerzon total.

    Returns {
        'spreads': {('N', 1.5): 'sel_id', ('P', 1.5): 'sel_id', ...},
        'totals': {('O', 7.5): 'sel_id', ('U', 7.5): 'sel_id', ...},
    }
    """
    resp = session.get(
        f"{DK_SGP_PARLAYS_URL}/{dk_event_id}",
        timeout=60,
    )

    if resp.status_code != 200:
        return {"spreads": {}, "totals": {}}

    text = resp.text

    # Parse all selection IDs from the response
    spread_matches = re.findall(r'0HC(\d+)([NP])(\d+)(_\d+)', text)
    total_matches = re.findall(r'0OU(\d+)([OU])(\d+)(_\d+)', text)

    # Find market numbers that have BOTH spreads and totals — these are
    # game line markets (main + alt). Inning-only markets have one or the other.
    # Real market numbers are 8+ digits; filter out short spurious matches.
    spread_mnums = set(mnum for mnum, _, _, _ in spread_matches if len(mnum) >= 8)
    total_mnums = set(mnum for mnum, _, _, _ in total_matches if len(mnum) >= 8)
    game_line_mnums = spread_mnums & total_mnums
    if main_market_num:
        game_line_mnums.add(main_market_num)

    # Build spread lookup — only from game line markets, prefer main market
    spreads = {}
    for mnum, sign, line, suf in spread_matches:
        if len(mnum) < 8 or mnum not in game_line_mnums:
            continue
        line_val = int(line) / 100
        key = (sign, line_val)
        sel_id = f"0HC{mnum}{sign}{line}{suf}"
        # Main market overwrites alt; alt only fills gaps
        if key not in spreads or mnum == main_market_num:
            spreads[key] = sel_id

    # Build total lookup — only from game line markets, prefer main market
    totals = {}
    for mnum, ou, line, suf in total_matches:
        if len(mnum) < 8 or mnum not in game_line_mnums:
            continue
        line_val = int(line) / 100
        key = (ou, line_val)
        sel_id = f"0OU{mnum}{ou}{line}{suf}"
        if key not in totals or mnum == main_market_num:
            totals[key] = sel_id

    if verbose:
        spread_lines = sorted(set(k[1] for k in spreads))
        total_lines = sorted(set(k[1] for k in totals if k[0] == 'O'))
        print(f"    Spreads available: {spread_lines}")
        print(f"    Totals available (O): {total_lines}")

    return {"spreads": spreads, "totals": totals}


# ---------------------------------------------------------------------------
# Step 6: Call calculateBets for SGP pricing
# ---------------------------------------------------------------------------

def calculate_sgp(session: cffi_requests.Session,
                  spread_sel: str, total_sel: str,
                  verbose: bool = False) -> dict | None:
    """Call DK's calculateBets with two selections. Returns {trueOdds, displayOdds} or None."""
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

    if resp.status_code == 422:
        if verbose:
            try:
                err = resp.json()
                code = err.get("statusCode", "")
                desc = err.get("description", "")
                print(f"      422: {code} — {desc[:100]}")
            except Exception:
                pass
        return None

    if resp.status_code != 200:
        if verbose:
            print(f"      HTTP {resp.status_code}")
        return None

    data = resp.json()

    # Check for combinability restrictions (cross-market rejection)
    restrictions = data.get("combinabilityRestrictions", [])
    if restrictions:
        if verbose:
            print(f"      NonCombinable: selections can't be combined in SGP")
        return None

    for bet in data.get("bets", []):
        mapped = bet.get("selectionsMapped", [])
        if bet.get("trueOdds") and len(mapped) >= 2:
            return {
                "trueOdds": bet["trueOdds"],
                "displayOdds": bet.get("displayOdds", ""),
            }
    return None


# ---------------------------------------------------------------------------
# Step 7: Scrape all 4 combos for a game
# ---------------------------------------------------------------------------

def scrape_game(session: cffi_requests.Session, game: dict,
                sel_ids: dict, verbose: bool = False) -> list[dict] | None:
    """
    Scrape all 4 SGP combos using exact selection IDs.
    Returns list of results, or None if selection IDs not found (signals failure).
    """
    results = []
    spread_line = game["spread_line"]
    total = game["total_line"]

    # Determine home/away spread signs
    if spread_line < 0:
        home_sign, away_sign = "N", "P"
    else:
        home_sign, away_sign = "P", "N"

    spread = abs(spread_line)

    # Look up exact selection IDs
    home_spread_sel = sel_ids["spreads"].get((home_sign, spread))
    away_spread_sel = sel_ids["spreads"].get((away_sign, spread))
    over_sel = sel_ids["totals"].get(("O", total))
    under_sel = sel_ids["totals"].get(("U", total))

    if not home_spread_sel or not away_spread_sel:
        print(f"    Spread ±{spread} not found in DK selection IDs")
        return None

    if not over_sel or not under_sel:
        print(f"    Total {total} not found in DK selection IDs")
        return None

    combos = [
        ("Home Spread + Over",  home_spread_sel, over_sel),
        ("Home Spread + Under", home_spread_sel, under_sel),
        ("Away Spread + Over",  away_spread_sel, over_sel),
        ("Away Spread + Under", away_spread_sel, under_sel),
    ]

    for combo_name, spread_sel, total_sel in combos:
        if verbose:
            print(f"    {combo_name}: {spread_sel} + {total_sel}")

        sgp = calculate_sgp(session, spread_sel, total_sel, verbose)

        if sgp:
            odds = sgp["trueOdds"]
            am = decimal_to_american(odds)
            print(f"    {combo_name}: {odds:.4f} ({am:+d})")
            results.append({
                "combo_name": combo_name,
                "trueOdds": odds,
                "displayOdds": sgp["displayOdds"],
            })
        else:
            print(f"    {combo_name}: no price")

        time.sleep(0.1)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scrape_dk_sgp(verbose: bool = False):
    """Main: fetch DK SGP odds for all MLB games via pure REST API."""

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

    # Deduplicate: keep only the earliest event per team matchup
    # DK returns both today's and tomorrow's games
    seen_matchups = set()
    deduped_events = []
    for evt in sorted(dk_events, key=lambda e: e["start_time"]):
        matchup = (evt["dk_home"], evt["dk_away"])
        if matchup not in seen_matchups:
            seen_matchups.add(matchup)
            deduped_events.append(evt)

    print("Matching teams via canonical_match...")
    matched = match_events(deduped_events, parlay_lines)
    print(f"  {len(matched)} matched games")

    if not matched:
        print("  No matches found.")
        return

    ensure_table()
    all_rows = []
    consecutive_failures = 0

    for game in matched:
        if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
            print("  Re-initializing session (possible expiry)...")
            session = init_session()
            consecutive_failures = 0

        print(f"\n{'='*60}")
        print(f"  {game['away_team']} @ {game['home_team']}")
        print(f"  Spread: {game['spread_line']:+.1f} | Total: {game['total_line']}")
        print(f"{'='*60}")

        # Fetch main market number (for prioritizing full-game IDs)
        main_market_num = fetch_main_market_num(session, game["dk_event_id"])

        # Fetch all selection IDs for this game
        sel_ids = fetch_selection_ids(
            session, game["dk_event_id"], main_market_num, verbose,
        )
        if not sel_ids["spreads"]:
            print("  No selection IDs found — skipping")
            consecutive_failures += 1
            continue

        game_results = scrape_game(session, game, sel_ids, verbose)

        # Retry once if all combos failed (DK may temporarily reject)
        if game_results is not None and len(game_results) == 0:
            time.sleep(2)
            game_results = scrape_game(session, game, sel_ids, verbose)

        if game_results is None:
            consecutive_failures += 1
            continue

        consecutive_failures = 0

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
