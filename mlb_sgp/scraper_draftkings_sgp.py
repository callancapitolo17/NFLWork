#!/usr/bin/env python3
"""DraftKings MLB SGP Scraper — thin shim invoking the draftkings library.

Reads target_lines from MLB_DB (default: mlb_mm.duckdb; bot overrides via
MLB_SGP_DB_PATH env var). Calls mlb_sgp.draftkings.price_sgps() with the
periods configured via MLB_SGP_PERIODS env var (default: FG,F5).
Writes PricedRow results back to MLB_DB via mlb_sgp.db.upsert_priced_rows.

Legacy helpers (fetch_dk_events, match_events, calculate_sgp, etc.) are
preserved in this file because the new draftkings orchestrator + dk_client
import them lazily. They stay here during the transition; a follow-up
refactor can lift them into the library module.
"""

import os
import re
import sys
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

from db import MLB_DB, _connect_with_retry
from integer_line_derivation import is_integer_line, derive_fair_probs

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
    """Load FG and F5 spread + total lines from mlb_parlay_lines staging table.

    The staging table is written by mlb_correlated_parlay.R before calling
    this scraper, breaking the old circular dependency where scrapers read
    from mlb_parlay_opportunities (the output of the script that calls us).

    Returns dict keyed by game_id with both fg_* and f5_* lines.
    """
    # Use _connect_with_retry so parallel scraper runs don't lose the read
    # to another scraper's brief write lock on mlb_mm.duckdb.
    con = _connect_with_retry(str(MLB_DB), read_only=True)
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
# Integer-line fallback helper
# ---------------------------------------------------------------------------

def try_integer_fallback_dk(
    session,
    sel_ids: dict,
    spread_line: float,
    total_line: float,
    canonical: set,
    verbose: bool = False,
):
    """Try to derive 4 fair_probs at an integer total line from adjacent
    half-point alts. Returns the derive_fair_probs result dict or None.

    Mirrors the spread/total-sel lookup the main loop does for the matched
    line, but for the two adjacent half-point alts.
    """
    if not is_integer_line(total_line):
        return None

    lo_total = total_line - 0.5
    hi_total = total_line + 0.5

    # Sign convention same as main loop
    if spread_line < 0:
        home_sign, away_sign = "N", "P"
    else:
        home_sign, away_sign = "P", "N"
    spread = abs(spread_line)

    home_spread_sels = sel_ids["spreads"].get((home_sign, spread, "1")) or []
    away_spread_sels = sel_ids["spreads"].get((away_sign, spread, "3")) or []

    over_lo_sels  = sel_ids["totals"].get(("O", lo_total)) or []
    under_lo_sels = sel_ids["totals"].get(("U", lo_total)) or []
    over_hi_sels  = sel_ids["totals"].get(("O", hi_total)) or []
    under_hi_sels = sel_ids["totals"].get(("U", hi_total)) or []

    if not (home_spread_sels and away_spread_sels and
            over_lo_sels and under_lo_sels and over_hi_sels and under_hi_sels):
        if verbose:
            print(f"      integer fallback: missing alts for {total_line}")
        return None

    # Helper: price one canonical-canonical combo, return decimal or None
    def _price(sp_sels, tot_sels):
        for sp in sp_sels:
            sp_mnum = _market_num(sp)
            if sp_mnum not in canonical:
                continue
            for to in tot_sels:
                to_mnum = _market_num(to)
                if to_mnum not in canonical:
                    continue
                sgp = calculate_sgp(session, sp, to, verbose=verbose)
                if sgp:
                    return sgp["trueOdds"]
        return None

    # 8 calls: 4 combos × 2 alts
    decimals_lo = {
        "home_over":  _price(home_spread_sels, over_lo_sels),
        "home_under": _price(home_spread_sels, under_lo_sels),
        "away_over":  _price(away_spread_sels, over_lo_sels),
        "away_under": _price(away_spread_sels, under_lo_sels),
    }
    decimals_hi = {
        "home_over":  _price(home_spread_sels, over_hi_sels),
        "home_under": _price(home_spread_sels, under_hi_sels),
        "away_over":  _price(away_spread_sels, over_hi_sels),
        "away_under": _price(away_spread_sels, under_hi_sels),
    }

    if any(d is None for d in decimals_lo.values()) or any(d is None for d in decimals_hi.values()):
        if verbose:
            print(f"      integer fallback: pricing call failed for {total_line}")
        return None

    return derive_fair_probs(decimals_lo, decimals_hi)


# ---------------------------------------------------------------------------
# Step 4: Match DK events to our game_ids
# ---------------------------------------------------------------------------

def _utc_bucket(ts) -> str:
    """Extract a UTC "YYYY-MM-DDTHH" bucket string from a timestamp.

    Used as a match key between WZ parlay_lines commence_time and DK event
    start_time. Date + hour granularity is required because hour alone would
    collide two games on different days that happen to start at the same
    UTC hour (e.g. a 2-game series between the same teams on consecutive
    days, both scheduled at 01:40 UTC).

    Accepts either an ISO string (from DK/FD API responses) or a datetime
    object (from DuckDB TIMESTAMP columns). Returns "" if the input is
    empty/None so callers can fall back to team-only matching.
    """
    if not ts:
        return ""
    # DuckDB returns datetime objects for TIMESTAMP columns
    if hasattr(ts, "year"):
        return f"{ts.year:04d}-{ts.month:02d}-{ts.day:02d}T{ts.hour:02d}"
    # API responses are ISO strings like "2026-06-15T17:05:00Z"
    return ts[:13] if len(ts) >= 13 else ""


def match_events(dk_events: list[dict], parlay_lines: dict) -> list[dict]:
    """Match DK events to our game_ids using canonical_match.py.

    Matching uses team names AND a UTC date+hour bucket so that (a)
    doubleheaders (two games same teams, same day, different hours) map to
    the correct game_id, and (b) consecutive-day series games between the
    same teams at the same hour don't collapse onto each other.
    """
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
        dk_bucket = _utc_bucket(dk_evt["start_time"])

        for game_id, lines in parlay_lines.items():
            if lines["home_team"] == canon_home and lines["away_team"] == canon_away:
                pl_bucket = _utc_bucket(lines.get("commence_time", ""))
                # If commence_time is available, require date+hour match so
                # same-teams games on different days (or in different hours
                # of the same day) don't collide. If missing (backward
                # compat), fall back to team-only matching.
                if pl_bucket and dk_bucket and pl_bucket != dk_bucket:
                    continue
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


# Selection IDs look like 0HC84252605P150_3 (spread) or 0OU84265271U750_3 (total).
# The 8+ digits immediately after the 0HC/0OU prefix are the DK market_num. A
# "same-market" selection pair means both legs were priced from the same DK
# market; cross-market pairs can look acceptable to calculateBets but return
# nonsense prices that don't match what DK actually offers in the UI.
_MARKET_NUM_RE = re.compile(r"^0(?:HC|OU)(\d+)")


def _market_num(sel_id: str) -> str:
    m = _MARKET_NUM_RE.match(sel_id)
    return m.group(1) if m else ""


# Moneyline selection ids look like 0ML85265198_1 (home, participant 1) or
# 0ML85265198_3 (away, participant 3) — same shape as spreads/totals but no
# line. The 8+ digits after 0ML are the DK market_num.
_ML_RE = re.compile(r"0ML(\d+)(_\d+)")


def _extract_moneyline_selections(text: str, ml_mnum: str | None) -> dict:
    """Pull home/away moneyline selection ids for one market_num out of the raw
    parlays response text. Returns ``{"1": [home_sels], "3": [away_sels]}``
    (participant "1" = home, "3" = away). Empty lists when ml_mnum is falsy or
    no selections match — keeps the per-book sign/participant convention
    independently testable, mirroring _extract_offered_lines_dk.
    """
    out = {"1": [], "3": []}
    if not ml_mnum:
        return out
    for mnum, suf in _ML_RE.findall(text):
        if mnum != ml_mnum:
            continue
        participant = suf[1:]          # "_1" -> "1"
        sel_id = f"0ML{mnum}{suf}"
        if participant in out and sel_id not in out[participant]:
            out[participant].append(sel_id)
    return out


def fetch_main_market_nums(session: cffi_requests.Session, dk_event_id: str) -> dict:
    """Fetch main Run Line + Total + Moneyline market numbers for both FG and F5.

    Returns {
        "fg": {"run_line": "...", "total": "...", "moneyline": "..."},
        "f5": {"run_line": "...", "total": "...", "moneyline": "..."},
    }
    Any field may be None if the market is unavailable.
    """
    out = {
        "fg": {"run_line": None, "total": None, "moneyline": None},
        "f5": {"run_line": None, "total": None, "moneyline": None},
    }
    for m_id, name in _fetch_subcat_markets(session, dk_event_id, "4519"):
        if name == "Run Line":
            out["fg"]["run_line"] = _strip_prefix(m_id)
        elif name == "Total":
            out["fg"]["total"] = _strip_prefix(m_id)
        elif name == "Moneyline":
            out["fg"]["moneyline"] = _strip_prefix(m_id)
    for m_id, name in _fetch_subcat_markets(session, dk_event_id, "15628"):
        if name == "Run Line - 1st 5 Innings":
            out["f5"]["run_line"] = _strip_prefix(m_id)
        elif name == "Total Runs - 1st 5 Innings":
            out["f5"]["total"] = _strip_prefix(m_id)
        elif name == "Moneyline - 1st 5 Innings":
            out["f5"]["moneyline"] = _strip_prefix(m_id)
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

    empty = {"fg": {"spreads": {}, "totals": {}, "moneyline": {"1": [], "3": []}},
             "f5": {"spreads": {}, "totals": {}, "moneyline": {"1": [], "3": []}}}
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

    # Each (sign, line, participant) key maps to a LIST of candidate sel_ids.
    # The suffix _1 = home team (participant 1), _3 = away team (participant 3).
    # These represent DIFFERENT TEAMS, not open/closed variants.  Keying by
    # participant ensures the pricing code selects the correct team's spread.
    def assign_spread(per_mnums, per_main, per_key):
        bucket = out[per_key]["spreads"]
        for mnum, sign, line, suf in spread_matches:
            if mnum not in per_mnums:
                continue
            line_val = int(line) / 100
            participant = suf[1:]  # "_1" -> "1", "_3" -> "3"
            key = (sign, line_val, participant)
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

    # Canonical markets per period = main RL + main total + any market that has
    # BOTH spreads AND totals (primary alt markets). Cross-pairs between two
    # canonical markets are allowed during pricing — e.g. main spread from M1
    # plus alt total from primary alt M2 is a legitimate SGP combo that DK's
    # UI supports. Secondary markets (totals-only standalones like DK's
    # 84267310 with non-canonical prices) are excluded, so we never pair a
    # main leg against those weird prices.
    out["fg"]["canonical"] = {m for m in (fg_rl, fg_tot) if m} | (
        fg_spread_mnums & fg_total_mnums
    )
    out["f5"]["canonical"] = {m for m in (f5_rl, f5_tot) if m} | (
        f5_spread_mnums & f5_total_mnums
    )

    # Moneyline selections (home/away) per period, keyed off the main-market
    # moneyline market_num. ML has no line, so it lives in its own bucket rather
    # than the (sign, line, participant) spread map.
    fg_ml = main_market_nums["fg"].get("moneyline")
    f5_ml = main_market_nums["f5"].get("moneyline")
    out["fg"]["moneyline"] = _extract_moneyline_selections(text, fg_ml)
    out["f5"]["moneyline"] = _extract_moneyline_selections(text, f5_ml)

    if verbose:
        for per in ("fg", "f5"):
            sp = sorted(set(k[1] for k in out[per]["spreads"]))
            to = sorted(set(k[1] for k in out[per]["totals"] if k[0] == 'O'))
            print(f"    [{per.upper()}] spreads: {sp}  totals(O): {to}  "
                  f"canonical: {sorted(out[per]['canonical'])}")

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
# Shim entry point
# ---------------------------------------------------------------------------
#
# Everything below replaces the legacy scrape_dk_sgp() + main() pair. The shim
# delegates orchestration to mlb_sgp.draftkings.price_sgps(), which composes
# the helpers above. Helpers stay defined at module level because
# mlb_sgp/draftkings.py and mlb_sgp/dk_client.py import them lazily.


def main():
    """Entry point: load targets, price SGPs, write rows.

    Path / import gymnastics: this script can be launched in three modes —
    (1) `python mlb_sgp/scraper_draftkings_sgp.py` from the repo root
        (test_sgp_regression's shim test path),
    (2) `python scraper_draftkings_sgp.py` from inside mlb_sgp/
        (legacy direct-invocation path, also how kalshi_mlb_rfq spawns it),
    (3) `python -m mlb_sgp.scraper_draftkings_sgp` (package mode).

    For (1) and (2) we need `mlb_sgp` importable as a top-level package so
    `mlb_sgp.draftkings` resolves. Make that work by ensuring the repo root
    is on sys.path before importing the orchestrator + db modules.
    """
    if str(_REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(_REPO_ROOT))

    from mlb_sgp import db
    from mlb_sgp._shared import load_target_lines

    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)

    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0

    periods_raw = os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(",")
    periods = tuple(p.strip() for p in periods_raw if p.strip())

    from mlb_sgp import draftkings
    print(f"  DK shim: {len(targets)} target lines, periods={periods}")
    rows = draftkings.price_sgps(targets, periods=periods, verbose=False)
    print(f"  DK shim: priced {len(rows)} rows")

    # Wipe both source labels so stale rows from a previous run never linger.
    # The orchestrator tags _direct vs _interpolated rows; clearing both
    # preserves the "fresh prices only" invariant per source.
    db.clear_source("draftkings_direct", db_path=db_path)
    db.clear_source("draftkings_interpolated", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
