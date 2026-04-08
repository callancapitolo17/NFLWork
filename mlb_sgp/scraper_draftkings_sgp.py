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
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    """Load FG and F5 spread + total lines from mlb_parlay_opportunities.

    Returns dict keyed by game_id with both fg_* and f5_* lines.
    """
    con = duckdb.connect(str(MLB_DB), read_only=True)
    try:
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_parlay_opportunities" not in tables:
            print("  No mlb_parlay_opportunities table — run the MLB pipeline first.")
            return {}

        # FG combos — the R pipeline writes only the highest-edge combo per
        # game, so any one of the 4 FG combo names may be present. The
        # spread/total lines are identical across all 4 combos for a given
        # game, so we collapse with ANY_VALUE.
        fg_rows = con.execute("""
            SELECT game_id,
                   ANY_VALUE(spread_line) AS spread_line,
                   ANY_VALUE(total_line)  AS total_line,
                   ANY_VALUE(home_team)   AS home_team,
                   ANY_VALUE(away_team)   AS away_team
            FROM mlb_parlay_opportunities
            WHERE combo IN ('Home Spread + Over', 'Home Spread + Under',
                            'Away Spread + Over', 'Away Spread + Under')
            GROUP BY game_id
        """).fetchall()

        f5_rows = con.execute("""
            SELECT game_id,
                   ANY_VALUE(spread_line) AS spread_line,
                   ANY_VALUE(total_line)  AS total_line
            FROM mlb_parlay_opportunities
            WHERE combo IN ('F5 Home Spread + Over', 'F5 Home Spread + Under',
                            'F5 Away Spread + Over', 'F5 Away Spread + Under')
            GROUP BY game_id
        """).fetchall()

        result = {}
        for row in fg_rows:
            result[row[0]] = {
                "fg_spread_line": row[1],
                "fg_total_line": row[2],
                "home_team": row[3],
                "away_team": row[4],
                "f5_spread_line": None,
                "f5_total_line": None,
            }
        for row in f5_rows:
            if row[0] in result:
                result[row[0]]["f5_spread_line"] = row[1]
                result[row[0]]["f5_total_line"] = row[2]

        return result
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
                    "fg_spread_line": lines["fg_spread_line"],
                    "fg_total_line": lines["fg_total_line"],
                    "f5_spread_line": lines["f5_spread_line"],
                    "f5_total_line": lines["f5_total_line"],
                })
                break

    return matched


# ---------------------------------------------------------------------------
# Step 5: Fetch exact selection IDs from SGP parlays endpoint
# ---------------------------------------------------------------------------

DK_EVENT_MARKETS_URL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "controldata/event/eventSubcategory/v1/markets"
)


def _fetch_subcat_markets(session, dk_event_id, subcat_id):
    """Fetch all markets in a given subcategory. Returns list of (id, name)."""
    resp = session.get(DK_EVENT_MARKETS_URL, params={
        "isBatchable": "false",
        "templateVars": dk_event_id,
        "marketsQuery": (
            f"$filter=eventId eq '{dk_event_id}' "
            f"AND clientMetadata/subCategoryId eq '{subcat_id}' "
            f"AND tags/all(t: t ne 'SportcastBetBuilder')"
        ),
        "include": "MarketSplits",
        "entity": "markets",
    }, timeout=15)
    if resp.status_code != 200:
        return []
    return [(m["id"], m.get("name", "")) for m in resp.json().get("markets", [])]


def _strip_prefix(mid: str) -> str:
    return mid.split("_")[-1] if "_" in mid else mid


def fetch_main_market_nums(session: cffi_requests.Session, dk_event_id: str) -> dict:
    """Fetch main Run Line + Total market numbers for both FG and F5.

    Returns {
        "fg": {"run_line": "...", "total": "..."},
        "f5": {"run_line": "...", "total": "..."},
    }
    Any field may be None if the market is unavailable.
    """
    out = {
        "fg": {"run_line": None, "total": None},
        "f5": {"run_line": None, "total": None},
    }
    for m_id, name in _fetch_subcat_markets(session, dk_event_id, "4519"):
        if name == "Run Line":
            out["fg"]["run_line"] = _strip_prefix(m_id)
        elif name == "Total":
            out["fg"]["total"] = _strip_prefix(m_id)
    for m_id, name in _fetch_subcat_markets(session, dk_event_id, "15628"):
        if name == "Run Line - 1st 5 Innings":
            out["f5"]["run_line"] = _strip_prefix(m_id)
        elif name == "Total Runs - 1st 5 Innings":
            out["f5"]["total"] = _strip_prefix(m_id)
    return out


def fetch_selection_ids(session: cffi_requests.Session, dk_event_id: str,
                        main_market_nums: dict | None = None,
                        verbose: bool = False) -> dict:
    """
    Fetch all SGP selection IDs for a game from the parlays endpoint, split
    by period (FG vs F5).

    Returns {
        'fg': {'spreads': {...}, 'totals': {...}},
        'f5': {'spreads': {...}, 'totals': {...}},
    }

    Period assignment for game line markets:
    - Main market_num for each period is provided as a seed
    - For other game line markets (alts), we partition by total line range:
      F5 totals are <= 5.5, FG totals are >= 6.0
    """
    main_market_nums = main_market_nums or {
        "fg": {"run_line": None, "total": None},
        "f5": {"run_line": None, "total": None},
    }
    fg_rl = main_market_nums["fg"].get("run_line")
    fg_tot = main_market_nums["fg"].get("total")
    f5_rl = main_market_nums["f5"].get("run_line")
    f5_tot = main_market_nums["f5"].get("total")

    resp = session.get(
        f"{DK_SGP_PARLAYS_URL}/{dk_event_id}",
        timeout=60,
    )

    empty = {"fg": {"spreads": {}, "totals": {}},
             "f5": {"spreads": {}, "totals": {}}}
    if resp.status_code != 200:
        return empty

    text = resp.text

    spread_matches = re.findall(r'0HC(\d+)([NP])(\d+)(_\d+)', text)
    total_matches = re.findall(r'0OU(\d+)([OU])(\d+)(_\d+)', text)

    # Compute spread/total range per market_num
    spread_lines_per_mnum = {}
    for mnum, _, line, _ in spread_matches:
        if len(mnum) < 8:
            continue
        spread_lines_per_mnum.setdefault(mnum, []).append(int(line) / 100)
    total_lines_per_mnum = {}
    for mnum, _, line, _ in total_matches:
        if len(mnum) < 8:
            continue
        total_lines_per_mnum.setdefault(mnum, []).append(int(line) / 100)

    # Period classification:
    # - main markets are explicit seeds
    # - alt spread markets: classify by max spread line (F5 ≤ 1.5, FG up to 5+)
    # - alt total markets: classify by max total line (F5 ≤ 5.5, FG ≥ 6.0)
    # - markets that don't look like game lines (innings/props) are excluded
    fg_spread_mnums = {fg_rl} if fg_rl else set()
    f5_spread_mnums = {f5_rl} if f5_rl else set()
    fg_total_mnums = {fg_tot} if fg_tot else set()
    f5_total_mnums = {f5_tot} if f5_tot else set()

    for mnum, lines in spread_lines_per_mnum.items():
        if mnum in fg_spread_mnums or mnum in f5_spread_mnums:
            continue
        # Only classify alt markets that ALSO appear as totals (game-line alt
        # markets pair both, e.g. 84220463). Pure-spread inning props skipped.
        if mnum not in total_lines_per_mnum:
            continue
        max_sp = max(lines)
        max_tot = max(total_lines_per_mnum[mnum])
        # F5 alts: small spreads AND small totals
        if max_sp <= 2.0 and max_tot <= 5.5:
            f5_spread_mnums.add(mnum)
            f5_total_mnums.add(mnum)
        else:
            fg_spread_mnums.add(mnum)
            fg_total_mnums.add(mnum)

    # Standalone total markets (no spreads in same mnum) — classify by total range
    for mnum, lines in total_lines_per_mnum.items():
        if mnum in fg_total_mnums or mnum in f5_total_mnums:
            continue
        if mnum in spread_lines_per_mnum:
            continue  # already handled above
        max_tot = max(lines)
        min_tot = min(lines)
        # Skip inning prop markets (very small totals like 0.5/1.5 only)
        if max_tot <= 2.0:
            continue
        if max_tot <= 5.5 and min_tot >= 2.5:
            f5_total_mnums.add(mnum)
        elif min_tot >= 4.0:
            fg_total_mnums.add(mnum)

    out = {"fg": {"spreads": {}, "totals": {}},
           "f5": {"spreads": {}, "totals": {}}}

    # Each (sign, line) key maps to a LIST of candidate sel_ids. DK returns
    # multiple suffix variants (_1, _3) per selection — one may be "closed"
    # while another is live. We'll try them in order, preferring the main
    # market and deduping. Pricing code retries through the list.
    def assign_spread(per_mnums, per_main, per_key):
        bucket = out[per_key]["spreads"]
        for mnum, sign, line, suf in spread_matches:
            if mnum not in per_mnums:
                continue
            line_val = int(line) / 100
            key = (sign, line_val)
            sel_id = f"0HC{mnum}{sign}{line}{suf}"
            lst = bucket.setdefault(key, [])
            if sel_id in lst:
                continue
            # Main market goes to front; alts to back
            if mnum == per_main:
                lst.insert(0, sel_id)
            else:
                lst.append(sel_id)

    def assign_total(per_mnums, per_main, per_key):
        bucket = out[per_key]["totals"]
        for mnum, ou, line, suf in total_matches:
            if mnum not in per_mnums:
                continue
            line_val = int(line) / 100
            key = (ou, line_val)
            sel_id = f"0OU{mnum}{ou}{line}{suf}"
            lst = bucket.setdefault(key, [])
            if sel_id in lst:
                continue
            if mnum == per_main:
                lst.insert(0, sel_id)
            else:
                lst.append(sel_id)

    assign_spread(fg_spread_mnums, fg_rl, "fg")
    assign_spread(f5_spread_mnums, f5_rl, "f5")
    assign_total(fg_total_mnums, fg_tot, "fg")
    assign_total(f5_total_mnums, f5_tot, "f5")

    if verbose:
        for per in ("fg", "f5"):
            sp = sorted(set(k[1] for k in out[per]["spreads"]))
            to = sorted(set(k[1] for k in out[per]["totals"] if k[0] == 'O'))
            print(f"    [{per.upper()}] spreads: {sp}  totals(O): {to}")

    return out


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
# Main
# ---------------------------------------------------------------------------

def scrape_dk_sgp(verbose: bool = False):
    """Main: fetch DK SGP odds for all MLB games via pure REST API.

    Uses parallel fetching for speed:
    - Phase 1: Fetch market nums + selection IDs for all games in parallel
    - Phase 2: Price all spread+total combos in parallel
    - Phase 3: Batch write to DuckDB
    """

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

    # ── Phase 1: Parallel fetch of market nums + selection IDs ──
    print("Fetching selection IDs (parallel)...")
    t0 = time.time()

    game_data = {}  # game_id -> sel_ids (per-period dict)

    def fetch_one_game(game):
        dk_eid = game["dk_event_id"]
        main_nums = fetch_main_market_nums(session, dk_eid)
        sel_ids = fetch_selection_ids(session, dk_eid, main_nums, verbose)
        return game["game_id"], sel_ids

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(fetch_one_game, g): g for g in matched}
        for future in as_completed(futures):
            try:
                game_id, sel_ids = future.result()
                game_data[game_id] = sel_ids
            except Exception as e:
                game = futures[future]
                print(f"  Error fetching {game['dk_name']}: {e}")

    print(f"  Fetched {len(game_data)} games in {time.time() - t0:.1f}s")

    # ── Phase 2: Build combo pairs and price in parallel ──
    print("Pricing SGP combos (parallel)...")
    t1 = time.time()

    # Build all (game_id, period, combo_name, spread_sel, total_sel) tuples
    combo_items = []
    for game in matched:
        gid = game["game_id"]
        if gid not in game_data:
            continue

        sel_ids_per_period = game_data[gid]

        for period in ("fg", "f5"):
            spread_line = game[f"{period}_spread_line"]
            total = game[f"{period}_total_line"]
            if spread_line is None or total is None:
                continue

            sel_ids = sel_ids_per_period.get(period, {"spreads": {}, "totals": {}})
            if not sel_ids["spreads"]:
                continue

            if spread_line < 0:
                home_sign, away_sign = "N", "P"
            else:
                home_sign, away_sign = "P", "N"

            spread = abs(spread_line)

            home_spread_sels = sel_ids["spreads"].get((home_sign, spread)) or []
            away_spread_sels = sel_ids["spreads"].get((away_sign, spread)) or []
            over_sels = sel_ids["totals"].get(("O", total)) or []
            under_sels = sel_ids["totals"].get(("U", total)) or []

            if not home_spread_sels or not away_spread_sels or not over_sels or not under_sels:
                continue

            prefix = "" if period == "fg" else "F5 "
            for combo_name, sp_sels, tot_sels in [
                ("Home Spread + Over",  home_spread_sels, over_sels),
                ("Home Spread + Under", home_spread_sels, under_sels),
                ("Away Spread + Over",  away_spread_sels, over_sels),
                ("Away Spread + Under", away_spread_sels, under_sels),
            ]:
                combo_items.append((gid, period, prefix + combo_name, sp_sels, tot_sels))

    # Price all combos in parallel
    pricing_results = []  # (game_id, period, combo_name, sgp_result)

    def price_one(item):
        gid, period, combo_name, sp_sels, tot_sels = item
        # Try every (spread_variant, total_variant) combo until one returns a
        # price. DK has stale "_1"/"_3" suffix variants — one may be closed
        # while another is live.
        sgp = None
        for sp in sp_sels:
            for to in tot_sels:
                sgp = calculate_sgp(session, sp, to, verbose=verbose)
                if sgp:
                    break
            if sgp:
                break
        return gid, period, combo_name, sgp

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(price_one, item) for item in combo_items]
        for future in as_completed(futures):
            try:
                pricing_results.append(future.result())
            except Exception:
                pass

    print(f"  Priced {len(pricing_results)} combos in {time.time() - t1:.1f}s")

    # ── Phase 2b: Retry failed games ──
    # Find games where all 4 combos failed — retry once
    priced_by_game = {}
    for gid, period, combo_name, sgp in pricing_results:
        if sgp:
            priced_by_game.setdefault(gid, []).append((period, combo_name, sgp))

    failed_games = set()
    for game in matched:
        gid = game["game_id"]
        if gid in game_data and gid not in priced_by_game:
            failed_games.add(gid)

    if failed_games:
        retry_items = [item for item in combo_items if item[0] in failed_games]
        if retry_items:
            print(f"  Retrying {len(failed_games)} failed games...")
            time.sleep(1)
            with ThreadPoolExecutor(max_workers=6) as pool:
                futures = [pool.submit(price_one, item) for item in retry_items]
                for future in as_completed(futures):
                    try:
                        gid, period, combo_name, sgp = future.result()
                        if sgp:
                            priced_by_game.setdefault(gid, []).append((period, combo_name, sgp))
                    except Exception:
                        pass

    # ── Phase 3: Collect results and write to DuckDB ──
    ensure_table()
    all_rows = []

    # Build game lookup for printing
    game_lookup = {g["game_id"]: g for g in matched}

    for gid, combos in sorted(priced_by_game.items()):
        game = game_lookup[gid]
        print(f"\n  {game['away_team']} @ {game['home_team']}")
        for period, combo_name, sgp in combos:
            odds = sgp["trueOdds"]
            am = decimal_to_american(odds)
            print(f"    {combo_name}: {odds:.4f} ({am:+d})")
            all_rows.append({
                "game_id": gid,
                "combo": combo_name,
                "period": "FG" if period == "fg" else "F5",
                "bookmaker": "draftkings",
                "sgp_decimal": round(odds, 4),
                "sgp_american": am,
                "source": "draftkings_direct",
            })

    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} DK SGP odds in {time.time() - t0:.1f}s total")
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
