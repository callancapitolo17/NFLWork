#!/usr/bin/env python3
"""
FanDuel MLB SGP Scraper (Pure REST API)

Fetches Same Game Parlay (SGP) odds from FanDuel for MLB FG spread+total combos.
Uses curl_cffi with Chrome TLS impersonation. Three required headers:
  - x-application       (FD's API key, also passed as ?_ak= query param)
  - x-sportsbook-region (geo: NJ — same value works from any state, FD just routes)
  - x-px-context        (PerimeterX visitor token — captured from a real Chrome session)

How it works:
1. List today's MLB events from the scan endpoint (competitionId 11196870)
2. Match each event to our internal game_ids via canonical_match
3. For each game, GET event-page to extract Run Line + Total Runs runner IDs
4. For each game, build 4 combos (Home+Over, Home+Under, Away+Over, Away+Under)
   matching the same logic as wagerzon_odds/parlay_pricer.py
5. POST each combo to implyBets, parse the betCombinations entry where isSGM=true,
   that's the correlated SGP price
6. Upsert into mlb_sgp_odds with bookmaker='fanduel'

Notes:
- F5 (1st 5 Innings) SGPs not yet supported — FD does not appear to expose F5
  markets via event-page during the early-season window. Add when available.
- Only main lines (no alts), matching the Wagerzon parlay pricer.
- The PerimeterX token is hardcoded for v1. If it stops working we'll add a
  Playwright bootstrap step that fetches a fresh token from a real browser.

Usage:
    cd mlb_sgp
    source venv/bin/activate
    python3 scraper_fanduel_sgp.py
    python3 scraper_fanduel_sgp.py --verbose
"""

import argparse
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import duckdb
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
# FD API config
# ---------------------------------------------------------------------------

FD_AK = "FhMFpcPWXMeyZxOx"  # FD's "application" token, sent as both query param + header

# Captured from a real Chrome session against fanduel.com. PerimeterX visitor
# token. Long-lived in practice (days/weeks) but not forever. If FD starts
# returning 400s with empty bodies or implyBets stops returning the isSGM entry,
# refresh this from a fresh browser session via the recon script.
FD_PX_CONTEXT = "_pxvid=687ed2ad-33ac-11f1-ba08-95fbba040e2b;pxcts=6119da3c-33ae-11f1-b93d-4b2b160f34a7;"

FD_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/146.0.0.0 Safari/537.36",
    "Accept": "application/json",
    "Referer": "https://sportsbook.fanduel.com/",
    "Origin": "https://sportsbook.fanduel.com",
    "x-application": FD_AK,
    "x-sportsbook-region": "NJ",
    "x-px-context": FD_PX_CONTEXT,
}

FD_SCAN_URL = "https://scan.nj.sportsbook.fanduel.com/api/sports/navigation/facet/v1.0/search"
FD_EVENT_PAGE_URL = "https://api.sportsbook.fanduel.com/sbapi/event-page"
FD_IMPLY_BETS_URL = (
    "https://sib.nj.sportsbook.fanduel.com"
    "/api/sports/fixedodds/transactional/v1/implyBets?pricePolicy=SUGGESTED"
)

FD_MLB_COMPETITION_ID = 11196870  # MLB regular season


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1)))


def init_session() -> cffi_requests.Session:
    """curl_cffi session with Chrome TLS impersonation."""
    return cffi_requests.Session(impersonate="chrome")


# ---------------------------------------------------------------------------
# Step 1: Fetch FD MLB events
# ---------------------------------------------------------------------------

def fetch_fd_events(session: cffi_requests.Session) -> list[dict]:
    """List today's MLB events. Name format: 'Away (P Pitcher) @ Home (P Pitcher)'."""
    body = {
        "filter": {
            "competitionIds": [FD_MLB_COMPETITION_ID],
            "contentGroup": {"language": "en", "regionCode": "NAMERICA"},
            "marketLevels": ["AVB_EVENT"],
            "maxResults": 100,
            "productTypes": ["SPORTSBOOK"],
            "selectBy": "FIRST_TO_START",
        },
        "facets": [{"type": "EVENT"}],
        "currencyCode": "USD",
    }
    resp = session.post(
        FD_SCAN_URL,
        headers={**FD_HEADERS, "Content-Type": "application/json"},
        data=json.dumps(body),
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    # Events live in attachments.events; the scan endpoint returns them via tree walk
    events = {}

    def walk(o):
        if isinstance(o, dict):
            if "eventId" in o and "openDate" in o and "name" in o:
                events[o["eventId"]] = o
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for it in o:
                walk(it)

    walk(data)

    out = []
    for e in events.values():
        away, home = _parse_event_name(e.get("name", ""))
        if not away or not home:
            continue
        out.append({
            "fd_event_id": str(e["eventId"]),
            "name": e["name"],
            "fd_away": away,
            "fd_home": home,
            "open_date": e.get("openDate", ""),
        })
    return out


def _parse_event_name(name: str) -> tuple[str, str]:
    """'Away Team (P Pitcher) @ Home Team (P Pitcher)' -> ('Away Team', 'Home Team')."""
    if " @ " not in name:
        return "", ""
    away_raw, home_raw = name.split(" @ ", 1)
    # Strip pitcher in parens
    away = re.sub(r"\s*\([^)]*\)\s*$", "", away_raw).strip()
    home = re.sub(r"\s*\([^)]*\)\s*$", "", home_raw).strip()
    return away, home


# ---------------------------------------------------------------------------
# Step 2: Load FG parlay lines from DuckDB
# ---------------------------------------------------------------------------

def load_parlay_lines() -> dict:
    """Load FG spread + total lines from mlb_parlay_opportunities. FG only.

    Returns dict keyed by game_id with fg_spread_line, fg_total_line, home_team, away_team.
    Mirrors load_parlay_lines() in scraper_draftkings_sgp.py but skips F5 (FD doesn't
    expose F5 SGP markets at this time).
    """
    con = duckdb.connect(str(MLB_DB), read_only=True)
    try:
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_parlay_opportunities" not in tables:
            print("  No mlb_parlay_opportunities table — run the MLB pipeline first.")
            return {}

        rows = con.execute("""
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

        return {
            r[0]: {
                "fg_spread_line": r[1],
                "fg_total_line": r[2],
                "home_team": r[3],
                "away_team": r[4],
            }
            for r in rows
        }
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Step 3: Match FD events to our game_ids
# ---------------------------------------------------------------------------

def match_events(fd_events: list[dict], parlay_lines: dict) -> list[dict]:
    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")

    matched = []
    for ev in fd_events:
        resolved = resolve_team_names(
            ev["fd_away"], ev["fd_home"], team_dict, canonical_games,
        )
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        canon_away, canon_home = resolved

        for game_id, lines in parlay_lines.items():
            if lines["home_team"] == canon_home and lines["away_team"] == canon_away:
                matched.append({
                    "game_id": game_id,
                    "fd_event_id": ev["fd_event_id"],
                    "home_team": canon_home,
                    "away_team": canon_away,
                    "fd_name": ev["name"],
                    "fg_spread_line": lines["fg_spread_line"],
                    "fg_total_line": lines["fg_total_line"],
                })
                break
    return matched


# ---------------------------------------------------------------------------
# Step 4: Fetch event-page and extract runner IDs
# ---------------------------------------------------------------------------

def fetch_event_runners(session: cffi_requests.Session, fd_event_id: str) -> dict:
    """For one event, return {'spreads': {(sign, line): (marketId, selectionId)},
                              'totals':  {(O|U, line): (marketId, selectionId)}}.

    sign: 'N' for negative spread, 'P' for positive. line is abs value.
    """
    url = (f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={fd_event_id}"
           f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true")
    resp = session.get(url, headers=FD_HEADERS, timeout=20)
    if resp.status_code != 200:
        return {"spreads": {}, "totals": {}}

    data = resp.json()
    markets = []

    def walk(o):
        if isinstance(o, dict):
            if "marketId" in o and "runners" in o and "marketName" in o:
                markets.append(o)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for it in o:
                walk(it)

    walk(data)
    # Dedupe by marketId
    seen = {m["marketId"]: m for m in markets}

    spreads = {}
    totals = {}
    for mid, m in seen.items():
        name = m.get("marketName", "")
        if name == "Run Line":
            for run in m.get("runners", []):
                hc = run.get("handicap")
                sid = run.get("selectionId")
                if hc is None or sid is None:
                    continue
                sign = "N" if hc < 0 else "P"
                spreads[(sign, abs(hc))] = (mid, sid)
        elif name == "Total Runs":
            for run in m.get("runners", []):
                hc = run.get("handicap")
                sid = run.get("selectionId")
                rn = (run.get("runnerName") or "").lower()
                if hc is None or sid is None:
                    continue
                ou = "O" if rn.startswith("over") else "U" if rn.startswith("under") else None
                if ou is None:
                    continue
                totals[(ou, hc)] = (mid, sid)

    return {"spreads": spreads, "totals": totals}


# ---------------------------------------------------------------------------
# Step 5: Price one combo via implyBets
# ---------------------------------------------------------------------------

def price_combo(session: cffi_requests.Session,
                spread_market: str, spread_sel: int,
                total_market: str, total_sel: int,
                verbose: bool = False) -> dict | None:
    """POST implyBets, return the SGP combined entry (isSGM=true) or None."""
    body = {"betLegs": [
        {"legType": "SIMPLE_SELECTION",
         "betRunners": [{"runner": {"marketId": spread_market, "selectionId": spread_sel}}]},
        {"legType": "SIMPLE_SELECTION",
         "betRunners": [{"runner": {"marketId": total_market, "selectionId": total_sel}}]},
    ]}
    resp = session.post(
        FD_IMPLY_BETS_URL,
        headers={**FD_HEADERS, "Content-Type": "application/json"},
        data=json.dumps(body),
        timeout=15,
    )
    if resp.status_code != 200:
        if verbose:
            print(f"      HTTP {resp.status_code}")
        return None

    try:
        data = resp.json()
    except Exception:
        return None

    for bc in data.get("betCombinations", []):
        if not bc.get("isSGM"):
            continue
        odds = bc.get("winAvgOdds", {}).get("trueOdds", {}).get("decimalOdds", {}).get("decimalOdds")
        am = bc.get("winAvgOdds", {}).get("americanDisplayOdds", {}).get("americanOddsInt")
        if odds is None or am is None:
            continue
        return {"decimal": float(odds), "american": int(am)}

    if verbose:
        print(f"      No isSGM entry in response (singles only)")
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scrape_fd_sgp(verbose: bool = False):
    print("Loading parlay lines from DuckDB...")
    parlay_lines = load_parlay_lines()
    if not parlay_lines:
        return
    print(f"  {len(parlay_lines)} games with lines")

    print("Initializing FD session...")
    session = init_session()

    print("Fetching FanDuel MLB events...")
    fd_events = fetch_fd_events(session)
    print(f"  {len(fd_events)} FD events")

    # Filter out events that have already started — those return live (in-game)
    # markets with running-score adjusted handicaps that don't match our
    # pre-game lines. Compare openDate (ISO UTC) against now.
    from datetime import datetime, timezone
    now_utc = datetime.now(timezone.utc).isoformat()
    fd_events = [e for e in fd_events if e["open_date"] > now_utc]

    # Dedupe by matchup, prefer earliest upcoming game
    seen = set()
    deduped = []
    for ev in sorted(fd_events, key=lambda e: e["open_date"]):
        key = (ev["fd_home"], ev["fd_away"])
        if key not in seen:
            seen.add(key)
            deduped.append(ev)

    print("Matching teams via canonical_match...")
    matched = match_events(deduped, parlay_lines)

    # Dedupe by canonical game_id — multiple FD events can resolve to the same
    # game (e.g. a pre-game listing and a separate live in-progress listing of
    # the same matchup). Events are already sorted by open_date, so keeping
    # the first occurrence keeps the pre-game one (with main ±1.5 lines)
    # over the live one (with running-score adjusted lines).
    seen_gids = set()
    deduped_matches = []
    for m in matched:
        if m["game_id"] in seen_gids:
            continue
        seen_gids.add(m["game_id"])
        deduped_matches.append(m)
    matched = deduped_matches

    print(f"  {len(matched)} matched games")
    if not matched:
        return

    # ── Phase 1: parallel fetch of runner IDs ──
    print("Fetching event runners (parallel)...")
    t0 = time.time()
    game_runners = {}

    def fetch_one(game):
        return game["game_id"], fetch_event_runners(session, game["fd_event_id"])

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(fetch_one, g): g for g in matched}
        for f in as_completed(futures):
            try:
                gid, runners = f.result()
                game_runners[gid] = runners
            except Exception as e:
                print(f"  Error: {e}")

    print(f"  Fetched {len(game_runners)} games in {time.time() - t0:.1f}s")

    # ── Phase 2: build combo items ──
    combo_items = []
    for game in matched:
        gid = game["game_id"]
        if gid not in game_runners:
            continue
        runners = game_runners[gid]
        if not runners["spreads"] or not runners["totals"]:
            continue

        spread_line = game["fg_spread_line"]
        total_line = game["fg_total_line"]
        if spread_line is None or total_line is None:
            continue

        # Determine signs (mirrors DK scraper logic)
        if spread_line < 0:
            home_sign, away_sign = "N", "P"
        else:
            home_sign, away_sign = "P", "N"
        sp = abs(spread_line)

        home_spread = runners["spreads"].get((home_sign, sp))
        away_spread = runners["spreads"].get((away_sign, sp))
        over = runners["totals"].get(("O", total_line))
        under = runners["totals"].get(("U", total_line))

        if not (home_spread and away_spread and over and under):
            if verbose:
                missing = []
                if not home_spread: missing.append(f"home {home_sign}{sp}")
                if not away_spread: missing.append(f"away {away_sign}{sp}")
                if not over: missing.append(f"over {total_line}")
                if not under: missing.append(f"under {total_line}")
                print(f"  {game['away_team']} @ {game['home_team']}: missing {missing}")
            continue

        for combo_name, sp_pair, tot_pair in [
            ("Home Spread + Over",  home_spread, over),
            ("Home Spread + Under", home_spread, under),
            ("Away Spread + Over",  away_spread, over),
            ("Away Spread + Under", away_spread, under),
        ]:
            combo_items.append((gid, combo_name, sp_pair, tot_pair))

    print(f"Pricing {len(combo_items)} SGP combos (parallel)...")
    t1 = time.time()

    def price_one(item):
        gid, name, sp, tot = item
        result = price_combo(session, sp[0], sp[1], tot[0], tot[1], verbose=verbose)
        return gid, name, result

    pricing_results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(price_one, it) for it in combo_items]
        for f in as_completed(futures):
            try:
                pricing_results.append(f.result())
            except Exception as e:
                if verbose:
                    print(f"  pricing error: {e}")

    print(f"  Priced {sum(1 for _,_,r in pricing_results if r)}/{len(pricing_results)} "
          f"combos in {time.time() - t1:.1f}s")

    # ── Phase 3: write to DuckDB ──
    ensure_table()
    rows = []
    game_lookup = {g["game_id"]: g for g in matched}
    by_game = {}
    for gid, name, result in pricing_results:
        if not result:
            continue
        by_game.setdefault(gid, []).append((name, result))

    for gid, items in sorted(by_game.items()):
        g = game_lookup[gid]
        print(f"\n  {g['away_team']} @ {g['home_team']}")
        for name, r in items:
            print(f"    {name}: {r['decimal']:.4f} ({r['american']:+d})")
            rows.append({
                "game_id": gid,
                "combo": name,
                "period": "FG",
                "bookmaker": "fanduel",
                "sgp_decimal": round(r["decimal"], 4),
                "sgp_american": r["american"],
                "source": "fanduel_direct",
            })

    if rows:
        upsert_sgp_odds(rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(rows)} FD SGP odds in {time.time() - t0:.1f}s total")
        print(f"{'='*60}")
    else:
        print("\nNo SGP odds collected.")

    return rows


def main():
    parser = argparse.ArgumentParser(description="FanDuel MLB SGP Scraper")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    print("=" * 60)
    print("  FANDUEL MLB SGP SCRAPER (Pure REST)")
    print("=" * 60)
    scrape_fd_sgp(verbose=args.verbose)


if __name__ == "__main__":
    main()
