#!/usr/bin/env python3
"""
FanDuel MLB SGP Scraper (Pure REST API)

Fetches Same Game Parlay (SGP) odds from FanDuel for MLB spread+total combos.
Supports FG (full game) and F5 (first 5 innings), main + alt spreads and totals.
Uses curl_cffi with Chrome TLS impersonation. Three required headers:
  - x-application       (FD's API key, also passed as ?_ak= query param)
  - x-sportsbook-region (geo: NJ — same value works from any state, FD just routes)
  - x-px-context        (PerimeterX visitor token — captured from a real Chrome session)

How it works:
1. List today's MLB events from the scan endpoint (competitionId 11196870)
2. Match each event to our internal game_ids via canonical_match
3. For each game, GET event-page with SGP tab to extract ALL runner IDs
   (main + alt spreads and totals, FG + F5)
4. Match exact Wagerzon lines — if FD doesn't have the exact line, skip the game
5. For each matched combo, POST to implyBets and parse the isSGM=true entry
6. Upsert into mlb_sgp_odds with bookmaker='fanduel'

The PerimeterX token is hardcoded for v1. If FD starts returning 400s with empty
bodies or implyBets stops returning the isSGM entry, refresh from a real browser
session (see README.md for instructions).

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

from db import ensure_table, upsert_sgp_odds, clear_source, MLB_DB

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

# Market names on FD's SGP tab, mapped to (period, type)
_MARKET_MAP = {
    # FG main
    "Run Line":                              ("fg", "spreads", "main"),
    "Total Runs":                            ("fg", "totals",  "main"),
    # FG alt
    "Alternate Run Lines":                   ("fg", "spreads", "alt"),
    "Alternate Total Runs":                  ("fg", "totals",  "alt"),
    # F5 main
    "First 5 Innings Run Line":              ("f5", "spreads", "main"),
    "First 5 Innings Total Runs":            ("f5", "totals",  "main"),
    # F5 alt
    "First 5 Innings Alternate Run Lines":   ("f5", "spreads", "alt"),
    "First 5 Innings Alternate Total Runs":  ("f5", "totals",  "alt"),
}


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
    away = re.sub(r"\s*\([^)]*\)\s*$", "", away_raw).strip()
    home = re.sub(r"\s*\([^)]*\)\s*$", "", home_raw).strip()
    return away, home


# ---------------------------------------------------------------------------
# Step 2: Load FG + F5 parlay lines from DuckDB
# ---------------------------------------------------------------------------

def load_parlay_lines() -> dict:
    """Load FG and F5 spread + total lines from mlb_parlay_lines staging table.

    The staging table is written by mlb_correlated_parlay.R before calling
    this scraper, breaking the old circular dependency where scrapers read
    from mlb_parlay_opportunities (the output of the script that calls us).

    Returns dict keyed by game_id with fg_*, f5_* lines plus home/away teams.
    Mirrors load_parlay_lines() in scraper_draftkings_sgp.py.
    """
    con = duckdb.connect(str(MLB_DB), read_only=True)
    try:
        tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
        if "mlb_parlay_lines" not in tables:
            print("  No mlb_parlay_lines table — run the MLB pipeline first.")
            return {}

        rows = con.execute("""
            SELECT game_id, home_team, away_team,
                   fg_spread, fg_total, f5_spread, f5_total,
                   commence_time
            FROM mlb_parlay_lines
        """).fetchall()

        result = {}
        for row in rows:
            result[row[0]] = {
                "fg_spread_line": row[3],
                "fg_total_line": row[4],
                "home_team": row[1],
                "away_team": row[2],
                "f5_spread_line": row[5],
                "f5_total_line": row[6],
                "commence_time": row[7],  # UTC ISO string for doubleheader matching
            }

        return result
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Step 3: Match FD events to our game_ids
# ---------------------------------------------------------------------------

def _utc_hour(ts) -> str:
    """Extract the UTC hour ("HH") from a timestamp.

    Accepts either an ISO string (from FD API responses) or a
    datetime object (from DuckDB TIMESTAMP columns). Returns "" if the
    input is empty/None so callers can fall back to team-only matching.
    """
    if not ts:
        return ""
    # DuckDB returns datetime objects for TIMESTAMP columns
    if hasattr(ts, "hour"):
        return f"{ts.hour:02d}"
    # API responses are ISO strings like "2026-06-15T17:05:00Z"
    return ts[11:13] if len(ts) >= 13 else ""


def match_events(fd_events: list[dict], parlay_lines: dict) -> list[dict]:
    """Match FD events to our game_ids using canonical_match.py.

    Matching uses team names AND start-time hour (UTC) so that doubleheaders
    (two games between the same teams on the same day) map to the correct
    game_id instead of both collapsing onto game 1.
    """
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
        fd_hour = _utc_hour(ev.get("open_date", ""))

        for game_id, lines in parlay_lines.items():
            if lines["home_team"] == canon_home and lines["away_team"] == canon_away:
                pl_hour = _utc_hour(lines.get("commence_time", ""))
                # If commence_time is available, require hour match to
                # distinguish doubleheaders. If missing (backward compat),
                # fall back to team-only matching.
                if pl_hour and fd_hour and pl_hour != fd_hour:
                    continue
                matched.append({
                    "game_id": game_id,
                    "fd_event_id": ev["fd_event_id"],
                    "home_team": canon_home,
                    "away_team": canon_away,
                    "fd_home": ev["fd_home"],
                    "fd_away": ev["fd_away"],
                    "fd_name": ev["name"],
                    "fg_spread_line": lines["fg_spread_line"],
                    "fg_total_line": lines["fg_total_line"],
                    "f5_spread_line": lines["f5_spread_line"],
                    "f5_total_line": lines["f5_total_line"],
                })
                break
    return matched


# ---------------------------------------------------------------------------
# Step 4: Fetch event-page (SGP tab) and extract ALL runner IDs
# ---------------------------------------------------------------------------

def fetch_event_runners(session: cffi_requests.Session, fd_event_id: str,
                        fd_home: str, fd_away: str) -> dict:
    """Fetch all SGP-eligible runners for one event (main + alt, FG + F5).

    Returns:
        {'fg': {'spreads': {('home'|'away', line): (marketId, selectionId), ...},
                'totals':  {('O'|'U', line): (marketId, selectionId), ...}},
         'f5': {'spreads': {...}, 'totals': {...}}}

    For spreads, key is ('home', abs_line) or ('away', abs_line) — determined
    by matching the runner's team name against fd_home/fd_away. This handles
    alt markets where both teams appear at each line.

    For totals, key is ('O', line) or ('U', line) — unambiguous from runner name.
    """
    # Use SGP tab to get all 156+ markets (main + alts + F5)
    url = (f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={fd_event_id}"
           f"&tab=same-game-parlay-"
           f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true")
    resp = session.get(url, headers=FD_HEADERS, timeout=20)

    empty = {"fg": {"spreads": {}, "totals": {}},
             "f5": {"spreads": {}, "totals": {}}}
    if resp.status_code != 200:
        return empty

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
    seen = {m["marketId"]: m for m in markets}

    out = {"fg": {"spreads": {}, "totals": {}},
           "f5": {"spreads": {}, "totals": {}}}

    for mid, m in seen.items():
        name = m.get("marketName", "")
        if name not in _MARKET_MAP:
            continue
        period, mtype, main_or_alt = _MARKET_MAP[name]
        bucket = out[period][mtype]

        for run in m.get("runners", []):
            sid = run.get("selectionId")
            rn = run.get("runnerName") or ""
            hc = run.get("handicap")
            if sid is None:
                continue

            if mtype == "spreads":
                parsed = _parse_spread_runner(rn, hc, main_or_alt, fd_home, fd_away)
                if parsed is None:
                    continue
                side, line = parsed
                key = (side, line)
                # Main takes priority over alt
                if key not in bucket or main_or_alt == "main":
                    bucket[key] = (mid, sid)

            elif mtype == "totals":
                parsed = _parse_total_runner(rn, hc, main_or_alt)
                if parsed is None:
                    continue
                ou, line = parsed
                key = (ou, line)
                if key not in bucket or main_or_alt == "main":
                    bucket[key] = (mid, sid)

    return out


def _parse_spread_runner(rn: str, hc, main_or_alt: str,
                         fd_home: str, fd_away: str) -> tuple[str, float] | None:
    """Parse a spread runner into ('home'|'away', abs_line).

    Main market: runnerName='Cincinnati Reds', handicap=-1.5
    Alt market:  runnerName='Cincinnati Reds +3.5', handicap=0
    """
    if main_or_alt == "main":
        if hc is None:
            return None
        side = "home" if rn == fd_home else "away" if rn == fd_away else None
        if side is None:
            return None
        return (side, abs(hc))
    else:
        # Alt: parse 'Cincinnati Reds +3.5' -> team, signed line
        m = re.match(r'(.+?)\s*([+-]\d+\.?\d*)$', rn)
        if not m:
            return None
        team = m.group(1).strip()
        line = abs(float(m.group(2)))
        side = "home" if team == fd_home else "away" if team == fd_away else None
        if side is None:
            return None
        return (side, line)


def _parse_total_runner(rn: str, hc, main_or_alt: str) -> tuple[str, float] | None:
    """Parse a total runner into ('O'|'U', line).

    Main market: runnerName='Over', handicap=8
    Alt market:  runnerName='Over (8.5)', handicap=0
    """
    rn_lower = rn.lower()
    if main_or_alt == "main":
        if hc is None:
            return None
        ou = "O" if rn_lower.startswith("over") else "U" if rn_lower.startswith("under") else None
        if ou is None:
            return None
        return (ou, float(hc))
    else:
        # Alt: parse 'Over (8.5)' or 'Under (7.5)'
        m = re.match(r'(Over|Under)\s*\((\d+\.?\d*)\)', rn, re.IGNORECASE)
        if not m:
            return None
        ou = "O" if m.group(1).lower() == "over" else "U"
        line = float(m.group(2))
        return (ou, line)


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
    # Wipe all previous FD SGP prices so only this run's results exist.
    clear_source("fanduel_direct")

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

    # Dedupe by matchup + start hour, prefer earliest per slot.
    # This filters out live-game duplicates (same matchup, same hour) while
    # preserving both games of a doubleheader (same teams, different hours).
    seen = set()
    deduped = []
    for ev in sorted(fd_events, key=lambda e: e["open_date"]):
        key = (ev["fd_home"], ev["fd_away"], ev["open_date"][:13])
        if key not in seen:
            seen.add(key)
            deduped.append(ev)

    print("Matching teams via canonical_match...")
    matched = match_events(deduped, parlay_lines)

    # Dedupe by canonical game_id
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
        return game["game_id"], fetch_event_runners(
            session, game["fd_event_id"], game["fd_home"], game["fd_away"])

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

        sel_ids_per_period = game_runners[gid]

        for period in ("fg", "f5"):
            spread_line = game[f"{period}_spread_line"]
            total_line = game[f"{period}_total_line"]
            if spread_line is None or total_line is None:
                continue

            sel = sel_ids_per_period.get(period, {"spreads": {}, "totals": {}})
            if not sel["spreads"]:
                continue

            sp = abs(spread_line)

            home_spread = sel["spreads"].get(("home", sp))
            away_spread = sel["spreads"].get(("away", sp))
            over = sel["totals"].get(("O", total_line))
            under = sel["totals"].get(("U", total_line))

            if not (home_spread and away_spread and over and under):
                if verbose:
                    missing = []
                    if not home_spread: missing.append(f"home {sp}")
                    if not away_spread: missing.append(f"away {sp}")
                    if not over: missing.append(f"over {total_line}")
                    if not under: missing.append(f"under {total_line}")
                    print(f"  {game['away_team']} @ {game['home_team']} [{period.upper()}]: missing {missing}")
                continue

            prefix = "" if period == "fg" else "F5 "
            for combo_name, sp_pair, tot_pair in [
                ("Home Spread + Over",  home_spread, over),
                ("Home Spread + Under", home_spread, under),
                ("Away Spread + Over",  away_spread, over),
                ("Away Spread + Under", away_spread, under),
            ]:
                combo_items.append((gid, period, prefix + combo_name, sp_pair, tot_pair))

    print(f"Pricing {len(combo_items)} SGP combos (parallel)...")
    t1 = time.time()

    def price_one(item):
        gid, period, name, sp, tot = item
        result = price_combo(session, sp[0], sp[1], tot[0], tot[1], verbose=verbose)
        return gid, period, name, result

    pricing_results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(price_one, it) for it in combo_items]
        for f in as_completed(futures):
            try:
                pricing_results.append(f.result())
            except Exception as e:
                if verbose:
                    print(f"  pricing error: {e}")

    print(f"  Priced {sum(1 for _,_,_,r in pricing_results if r)}/{len(pricing_results)} "
          f"combos in {time.time() - t1:.1f}s")

    # ── Phase 3: write to DuckDB ──
    ensure_table()
    rows = []
    game_lookup = {g["game_id"]: g for g in matched}
    by_game = {}
    for gid, period, name, result in pricing_results:
        if not result:
            continue
        by_game.setdefault(gid, []).append((period, name, result))

    for gid, items in sorted(by_game.items()):
        g = game_lookup[gid]
        print(f"\n  {g['away_team']} @ {g['home_team']}")
        for period, name, r in items:
            print(f"    {name}: {r['decimal']:.4f} ({r['american']:+d})")
            rows.append({
                "game_id": gid,
                "combo": name,
                "period": "FG" if period == "fg" else "F5",
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
