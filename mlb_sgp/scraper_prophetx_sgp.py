#!/usr/bin/env python3
"""
ProphetX MLB SGP Scraper (Pure REST API)

Fetches Same Game Parlay (SGP) odds from ProphetX for MLB spread+total combos,
both full-game (FG) and first-5-innings (F5). Uses curl_cffi to post to the
RFQ endpoint and writes results to mlb_sgp_odds with source='prophetx_direct'.

Architectural notes (see recon_prophetx_sgp.py + memory notes for details):
    - ProphetX is a P2P exchange with an RFQ-style parlay pricer. The endpoint
      /parlay/public/api/v1/user/request returns an array of offers at
      different stake tiers. We use offers[0].odds (tightest available) as
      "the price", consistent with DK/FD one-price semantics.
    - All API calls go to www.prophetx.co — no separate backend domain.
    - Market tree has explicit names for FG ('Run Line', 'Total Runs') vs F5
      ('1st-5th Inning Spread', '1st-5th Inning Total Runs') — no heuristic
      period classification needed.
    - Selection IDs are structured (marketId, outcomeId, lineId, line), not
      regex-parsed strings.
    - LineIds are hex fingerprints that change when lines move — we fetch the
      market tree and POST the RFQ in quick succession.

Usage:
    cd mlb_sgp
    source venv/bin/activate
    python scraper_prophetx_sgp.py           # all games
    python scraper_prophetx_sgp.py --verbose  # show details
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

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

from db import ensure_table, upsert_sgp_odds, clear_source, MLB_DB, _connect_with_retry
from integer_line_derivation import is_integer_line, derive_fair_probs

# ---------------------------------------------------------------------------
# ProphetX API config
# ---------------------------------------------------------------------------
PX_BASE = "https://www.prophetx.co"
PX_TOURNAMENTS = f"{PX_BASE}/trade/public/api/v1/tournaments"
PX_EVENT_MARKETS = f"{PX_BASE}/trade/public/api/v2/events"
PX_PARLAY_REQUEST = f"{PX_BASE}/parlay/public/api/v1/user/request"

# Persistent profile from the recon run — cookies live here if auth is needed
PX_PROFILE = _THIS_DIR / ".prophetx_profile"

# Market-name mapping (captured verbatim from live market tree on 2026-04-24)
MARKET_NAMES = {
    "fg": {
        "spread": "Run Line",
        "total":  "Total Runs",
    },
    "f5": {
        # F5 market names discovered via recon (may vary; see NAME_ALIASES)
        "spread": "1st-5th Inning Spread",
        "total":  "1st-5th Inning Total Runs",
    },
}
# Accept a few likely variants in case ProphetX renames a market; first hit wins
NAME_ALIASES = {
    "1st-5th Inning Spread": ["1st-5th Inning Run Line", "1st 5 Innings Spread", "F5 Run Line"],
    "1st-5th Inning Total Runs": ["1st 5 Innings Total Runs", "F5 Total Runs"],
}

# RFQ behaviour
RFQ_STAKE = 1            # small probe stake — response returns full offer ladder regardless
RFQ_TIMEOUT = 10
MARKETS_TIMEOUT = 15

# Minimum stake threshold for offer selection. ProphetX's offer ladder
# interleaves low-stake "teaser" tiers (e.g. $50) with real-liquidity tiers
# ($2000+). offers[0] is often the teaser; the UI filters these out. We pick
# the best-odds offer whose available stake meets this threshold, matching
# the price a user would actually see on the site. Falls back to offers[0]
# if nothing clears the threshold (thin markets).
MIN_OFFER_STAKE = 150

# Sanity filter: drop combos where the parlay decimal exceeds this multiplier
# times the naive independent multiply of the two legs' single-line decimals.
# Legitimate anti-correlated combos legitimately reach ~1.15×; 1.5× gives
# headroom for stronger correlation while still catching the systematic
# ProphetX F5-Over bug (decimals 5-7× naive multiply). Tune if needed.
SANITY_MULT_RATIO = 1.5

PARALLEL_MARKETS = 4     # conservative vs DK's 6 — RFQ limits unknown
PARALLEL_PRICING = 4


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1)))


def american_to_decimal(am: int) -> float:
    if am >= 0:
        return 1 + am / 100.0
    return 1 + 100.0 / (-am)


def _utc_bucket(ts) -> str:
    """UTC YYYY-MM-DDTHH bucket for doubleheader matching (same idiom as DK)."""
    if not ts:
        return ""
    if hasattr(ts, "year"):
        return f"{ts.year:04d}-{ts.month:02d}-{ts.day:02d}T{ts.hour:02d}"
    return ts[:13] if len(ts) >= 13 else ""


# ---------------------------------------------------------------------------
# Session + auth
# ---------------------------------------------------------------------------
def init_session() -> cffi_requests.Session:
    """Create a curl_cffi session with Chrome TLS impersonation and attempt to
    load cookies from the Playwright profile if present. ProphetX's /public/
    paths *may* work anonymously; we load cookies preemptively so a single
    unauthenticated 401 doesn't derail the whole run.
    """
    session = cffi_requests.Session(impersonate="chrome")
    # Warm up + seed any anti-bot cookies
    try:
        session.get(PX_BASE, timeout=15)
    except Exception:
        pass

    cookies_loaded = _load_profile_cookies(session)
    if cookies_loaded:
        print(f"  Loaded {cookies_loaded} cookies from .prophetx_profile")
    else:
        print("  No cookies loaded (trying anonymous first; will fall back if 401)")
    return session


def _load_profile_cookies(session) -> int:
    """Extract cookies from the Playwright persistent context SQLite.

    Playwright stores cookies in .prophetx_profile/Default/Cookies (sqlite3).
    We read *without* locking the file. Returns number of cookies loaded.
    """
    db_path = PX_PROFILE / "Default" / "Cookies"
    if not db_path.exists():
        return 0
    try:
        import sqlite3
        # Copy-read is safer than reading a possibly-locked DB; use URI read-only
        con = sqlite3.connect(f"file:{db_path}?mode=ro&immutable=1", uri=True)
        cur = con.cursor()
        cur.execute(
            "SELECT name, value, host_key FROM cookies "
            "WHERE host_key LIKE '%prophetx%' OR host_key LIKE '%betprophet%'"
        )
        count = 0
        for name, value, host in cur.fetchall():
            session.cookies.set(name, value, domain=host.lstrip("."))
            count += 1
        con.close()
        return count
    except Exception as e:
        print(f"  (cookie load failed: {e} — continuing anonymously)")
        return 0


# ---------------------------------------------------------------------------
# Step 1: Event discovery
# ---------------------------------------------------------------------------
def fetch_prophetx_mlb_events(session) -> list[dict]:
    """List upcoming MLB events from the tournaments-with-events endpoint.

    Filters: sport.name == 'Baseball' AND tournament.name contains 'MLB' or
    'Major League'. If ProphetX renames tournaments, widen this filter.
    """
    resp = session.get(PX_TOURNAMENTS, params={
        "expand": "events", "type": "highlight", "limit": 150,
    }, timeout=15)
    resp.raise_for_status()

    events = []
    for t in resp.json().get("data", {}).get("tournaments", []):
        tname = t.get("name", "")
        is_mlb = ("MLB" in tname) or ("Major League" in tname)
        for evt in t.get("sportEvents", []) or []:
            sport_name = (evt.get("sport") or {}).get("name", "")
            if sport_name != "Baseball":
                continue
            if not is_mlb:
                continue
            competitors = evt.get("competitors") or []
            if len(competitors) < 2:
                continue
            # ProphetX uses seq 0=home, 1=away in recent data; fall back safely
            by_seq = {c.get("seq"): c for c in competitors}
            home = by_seq.get(0) or competitors[0]
            away = by_seq.get(1) or competitors[1]
            events.append({
                "px_event_id": evt.get("id"),
                "px_home": home.get("name", ""),
                "px_away": away.get("name", ""),
                "px_home_competitor_id": home.get("id"),
                "px_away_competitor_id": away.get("id"),
                "scheduled": evt.get("scheduled", ""),
                "name": evt.get("name", ""),
                "tournament": tname,
            })
    return events


# ---------------------------------------------------------------------------
# Step 2: Match events to our canonical game_ids
# ---------------------------------------------------------------------------
def load_parlay_lines() -> dict:
    """Read FG+F5 spread/total lines + canonical team names from mlb_parlay_lines."""
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
                   fg_spread, fg_total, f5_spread, f5_total, commence_time
            FROM mlb_parlay_lines
        """).fetchall()
        out = {}
        for row in rows:
            out[row[0]] = {
                "home_team": row[1], "away_team": row[2],
                "fg_spread_line": row[3], "fg_total_line": row[4],
                "f5_spread_line": row[5], "f5_total_line": row[6],
                "commence_time": row[7],
            }
        return out
    finally:
        con.close()


def match_events(px_events: list[dict], parlay_lines: dict) -> list[dict]:
    """Map ProphetX events to our game_ids via team-name canonicalization +
    UTC-hour bucket (handles doubleheaders). Same pattern as DK scraper."""
    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")
    matched = []
    for evt in px_events:
        resolved = resolve_team_names(
            evt["px_away"], evt["px_home"], team_dict, canonical_games,
        )
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        canon_away, canon_home = resolved
        px_bucket = _utc_bucket(evt["scheduled"])
        for gid, lines in parlay_lines.items():
            if lines["home_team"] != canon_home or lines["away_team"] != canon_away:
                continue
            pl_bucket = _utc_bucket(lines.get("commence_time", ""))
            if pl_bucket and px_bucket and pl_bucket != px_bucket:
                continue
            matched.append({
                "game_id": gid,
                "px_event_id": evt["px_event_id"],
                "home_team": canon_home, "away_team": canon_away,
                "px_home_competitor_id": evt["px_home_competitor_id"],
                "px_away_competitor_id": evt["px_away_competitor_id"],
                "fg_spread_line": lines["fg_spread_line"],
                "fg_total_line":  lines["fg_total_line"],
                "f5_spread_line": lines["f5_spread_line"],
                "f5_total_line":  lines["f5_total_line"],
            })
            break
    return matched


# ---------------------------------------------------------------------------
# Step 3: Fetch market tree + locate the specific legs we need
# ---------------------------------------------------------------------------
def _find_market(markets: list[dict], target_name: str) -> dict | None:
    """Locate a market by exact name, with alias fallback."""
    names_to_try = [target_name] + NAME_ALIASES.get(target_name, [])
    for name in names_to_try:
        for m in markets:
            if m.get("name") == name:
                return m
    return None


def _pick_selection(market: dict, predicate) -> dict | None:
    """Scan market.marketLines[].selections[][] for the first selection
    matching `predicate(selection_dict)`. Returns the raw selection dict
    so caller can extract outcomeId (selection['id']), lineID, line."""
    for line_grp in market.get("marketLines") or []:
        for side in line_grp.get("selections") or []:
            for sel in side:
                if predicate(sel):
                    return sel
    # Some markets (moneyline) have flat selections at the top level
    for side in market.get("selections") or []:
        for sel in side:
            if predicate(sel):
                return sel
    return None


def _verify_competitor_ids(markets, home_id, away_id) -> bool:
    """Sanity check: tournaments endpoint gives us (home_id, away_id); this
    verifies they match the moneyline outcomes in the markets endpoint. If
    the two endpoints use different ID namespaces (the reviewer's C1 concern)
    we'd silently get zero legs — this detects and loudly skips instead.
    """
    for m in markets:
        if m.get("name") == "Moneyline":
            outs = m.get("outcomes") or []
            ids_seen = {o.get("competitorId") for o in outs if o.get("competitorId") is not None}
            if ids_seen and ids_seen == {home_id, away_id}:
                return True
            # Namespace mismatch: log what we saw so debugging is painless
            return False
    return True  # no moneyline market — can't verify; don't block


def fetch_event_legs(session, game: dict, verbose: bool = False) -> dict:
    """For one game, fetch its market tree and extract the 4 canonical legs
    per period: home spread, away spread, over, under.

    Returns:
        {
            "fg": {
                "home_spread": {"marketId", "outcomeId", "lineId", "line"} | None,
                "away_spread": {...} | None,
                "over":        {...} | None,
                "under":       {...} | None,
            },
            "f5": {...}
        }
    Any leg may be None if its market isn't available for this event.
    """
    url = f"{PX_EVENT_MARKETS}/{game['px_event_id']}/markets"
    resp = session.get(url, timeout=MARKETS_TIMEOUT)
    out = {"fg": {"home_spread": None, "away_spread": None, "over": None, "under": None},
           "f5": {"home_spread": None, "away_spread": None, "over": None, "under": None}}
    if resp.status_code != 200:
        if verbose:
            print(f"      markets HTTP {resp.status_code} for {game['px_event_id']}")
        return out, []

    markets = resp.json().get("data", {}).get("markets", [])
    home_id = game["px_home_competitor_id"]
    away_id = game["px_away_competitor_id"]

    # C1 guard: if tournaments/markets use different competitorId namespaces,
    # every competitorId-based lookup below silently fails. Detect + skip loudly.
    if not _verify_competitor_ids(markets, home_id, away_id):
        print(f"  WARN: {game.get('home_team')} / {game.get('away_team')} "
              f"(event {game['px_event_id']}): competitorId mismatch between "
              f"tournaments (home={home_id}, away={away_id}) and markets — skipping. "
              "If this persists, switch to name-based outcome matching.")
        return out, markets

    def _leg_from_sel(sel, market_id):
        if not sel:
            return None
        return {
            "marketId":  market_id,
            "outcomeId": sel.get("id"),
            "lineId":    sel.get("lineID"),
            "line":      sel.get("line"),
            "leg_am":    sel.get("odds"),  # single-leg American odds for sanity check
        }

    for period in ("fg", "f5"):
        target_spread_line = game[f"{period}_spread_line"]
        target_total_line  = game[f"{period}_total_line"]
        if target_spread_line is None or target_total_line is None:
            continue

        # --- Spread market ---
        spread_mkt = _find_market(markets, MARKET_NAMES[period]["spread"])
        if spread_mkt:
            mid = spread_mkt.get("id")
            # Home side: competitor == home AND line == +target (if home is underdog)
            # ProphetX stores line from each outcome's perspective.
            # If parlay_lines reports home_spread = +1.5, home outcome has line == +1.5.
            # If home_spread = -1.5, home outcome has line == -1.5.
            home_sel = _pick_selection(
                spread_mkt,
                lambda s, lv=target_spread_line, hid=home_id: (
                    s.get("competitorId") == hid
                    and _line_eq(s.get("line"), lv)
                ),
            )
            away_sel = _pick_selection(
                spread_mkt,
                lambda s, lv=-target_spread_line, aid=away_id: (
                    s.get("competitorId") == aid
                    and _line_eq(s.get("line"), lv)
                ),
            )
            out[period]["home_spread"] = _leg_from_sel(home_sel, mid)
            out[period]["away_spread"] = _leg_from_sel(away_sel, mid)

        # --- Total market ---
        total_mkt = _find_market(markets, MARKET_NAMES[period]["total"])
        if total_mkt:
            mid = total_mkt.get("id")
            over_sel = _pick_selection(
                total_mkt,
                lambda s, lv=target_total_line: (
                    (s.get("name") or "").lower().startswith("over")
                    and _line_eq(s.get("line"), lv)
                ),
            )
            under_sel = _pick_selection(
                total_mkt,
                lambda s, lv=target_total_line: (
                    (s.get("name") or "").lower().startswith("under")
                    and _line_eq(s.get("line"), lv)
                ),
            )
            out[period]["over"]  = _leg_from_sel(over_sel,  mid)
            out[period]["under"] = _leg_from_sel(under_sel, mid)

        if verbose:
            missing = [k for k, v in out[period].items() if v is None]
            if missing:
                print(f"      [{period.upper()}] {game['game_id'][:8]}: "
                      f"missing {missing}  (spread={target_spread_line}, total={target_total_line})")

    return out, markets


def _line_eq(a, b, eps=1e-6) -> bool:
    """Float-safe equality for line values (they're halves or quarters)."""
    if a is None or b is None:
        return False
    try:
        return abs(float(a) - float(b)) < eps
    except (TypeError, ValueError):
        return False


# ---------------------------------------------------------------------------
# Integer-line fallback helper
# ---------------------------------------------------------------------------

def try_integer_fallback_px(
    session,
    markets: list,
    home_id,
    away_id,
    px_event_id,
    spread_line: float,
    total_line: float,
    period: str,
    verbose: bool = False,
):
    """PX: when WZ quotes an integer total that PX's market tree doesn't have
    at the exact strike, look up the two adjacent half-point alt total markets
    (total_line - 0.5 and total_line + 0.5), price all 8 combos via the shared
    submit_parlay_rfq (which already applies MIN_OFFER_STAKE filtering), and
    call derive_fair_probs to recover the integer-line fair probabilities.

    Args:
        markets: raw market list from fetch_event_legs.
        home_id / away_id: PX competitor IDs for spread leg lookup.
        px_event_id: ProphetX event ID — required for RFQ leg dicts.
        spread_line: target spread from mlb_parlay_lines (home perspective).
        total_line: the integer total line from mlb_parlay_lines.
        period: 'fg' or 'f5' — selects the correct market name strings.

    Returns a derive_fair_probs result dict or None.
    """
    if not is_integer_line(total_line):
        return None

    lo, hi = total_line - 0.5, total_line + 0.5

    def _leg_from_sel(sel, market_id):
        if not sel:
            return None
        return {
            "marketId":  market_id,
            "outcomeId": sel.get("id"),
            "lineId":    sel.get("lineID"),
            "line":      sel.get("line"),
            "leg_am":    sel.get("odds"),
        }

    # --- Spread legs at the exact spread_line ---
    spread_mkt = _find_market(markets, MARKET_NAMES[period]["spread"])
    if spread_mkt is None:
        if verbose:
            print(f"      integer fallback: no spread market for {period} spread={spread_line}")
        return None
    mid_spread = spread_mkt.get("id")
    home_sel = _pick_selection(
        spread_mkt,
        lambda s, lv=spread_line, hid=home_id: (
            s.get("competitorId") == hid and _line_eq(s.get("line"), lv)
        ),
    )
    away_sel = _pick_selection(
        spread_mkt,
        lambda s, lv=-spread_line, aid=away_id: (
            s.get("competitorId") == aid and _line_eq(s.get("line"), lv)
        ),
    )
    home_leg = _leg_from_sel(home_sel, mid_spread)
    away_leg = _leg_from_sel(away_sel, mid_spread)
    if not (home_leg and away_leg):
        if verbose:
            print(f"      integer fallback: spread legs missing at spread={spread_line}")
        return None

    # --- Adjacent total markets at lo and hi ---
    total_mkt = _find_market(markets, MARKET_NAMES[period]["total"])
    if total_mkt is None:
        if verbose:
            print(f"      integer fallback: no total market for {period}")
        return None
    mid_total = total_mkt.get("id")

    def _get_total_leg(name_prefix, strike):
        sel = _pick_selection(
            total_mkt,
            lambda s, p=name_prefix, lv=strike: (
                (s.get("name") or "").lower().startswith(p)
                and _line_eq(s.get("line"), lv)
            ),
        )
        return _leg_from_sel(sel, mid_total)

    over_lo  = _get_total_leg("over",  lo)
    under_lo = _get_total_leg("under", lo)
    over_hi  = _get_total_leg("over",  hi)
    under_hi = _get_total_leg("under", hi)

    if not all([over_lo, under_lo, over_hi, under_hi]):
        if verbose:
            missing = []
            if not over_lo:  missing.append(f"over {lo}")
            if not under_lo: missing.append(f"under {lo}")
            if not over_hi:  missing.append(f"over {hi}")
            if not under_hi: missing.append(f"under {hi}")
            print(f"      integer fallback: missing adjacent alt legs {missing} for total={total_line}")
        return None

    def _build_rfq_legs(sp_leg, tot_leg):
        return [
            {"sportEventId": px_event_id,
             "marketId": sp_leg["marketId"], "outcomeId": sp_leg["outcomeId"],
             "lineId":   sp_leg["lineId"],   "line": sp_leg["line"]},
            {"sportEventId": px_event_id,
             "marketId": tot_leg["marketId"], "outcomeId": tot_leg["outcomeId"],
             "lineId":   tot_leg["lineId"],   "line": tot_leg["line"]},
        ]

    def _price(sp_leg, tot_leg):
        priced, _ = submit_parlay_rfq(session, _build_rfq_legs(sp_leg, tot_leg), verbose=verbose)
        return priced["decimal"] if priced else None

    decimals_lo = {
        "home_over":  _price(home_leg, over_lo),
        "home_under": _price(home_leg, under_lo),
        "away_over":  _price(away_leg, over_lo),
        "away_under": _price(away_leg, under_lo),
    }
    decimals_hi = {
        "home_over":  _price(home_leg, over_hi),
        "home_under": _price(home_leg, under_hi),
        "away_over":  _price(away_leg, over_hi),
        "away_under": _price(away_leg, under_hi),
    }

    if any(d is None for d in decimals_lo.values()) or any(d is None for d in decimals_hi.values()):
        if verbose:
            print(f"      integer fallback: pricing failed for total={total_line}")
        return None

    return derive_fair_probs(decimals_lo, decimals_hi)


# ---------------------------------------------------------------------------
# Step 4: RFQ pricing
# ---------------------------------------------------------------------------
def submit_parlay_rfq(session, legs: list[dict], verbose: bool = False) -> tuple[dict | None, bool]:
    """POST two legs to the RFQ endpoint.

    Returns (priced, auth_failed) where:
      - priced is {'decimal','american','offers'} on success, else None
      - auth_failed is True iff we got a 401/403 (caller tracks for M2 alert)
    """
    payload = {"marketLines": legs, "stake": RFQ_STAKE}
    try:
        resp = session.post(
            PX_PARLAY_REQUEST, json=payload,
            headers={"Content-Type": "application/json"},
            timeout=RFQ_TIMEOUT,
        )
    except Exception as e:
        if verbose:
            print(f"      RFQ error: {e}")
        return None, False

    if resp.status_code == 401 or resp.status_code == 403:
        if verbose:
            print(f"      RFQ {resp.status_code} — auth required; cookies may be stale")
        return None, True
    if resp.status_code != 200:
        if verbose:
            body_preview = resp.text[:200] if resp.text else ""
            print(f"      RFQ HTTP {resp.status_code}: {body_preview}")
        return None, False

    try:
        data = resp.json()
    except Exception:
        return None, False

    if not data.get("success"):
        if verbose:
            print(f"      RFQ success=false: {data.get('error') or data}")
        return None, False

    offers = (data.get("data") or {}).get("offers") or []
    if not offers:
        return None, False

    # Pick the best-odds offer with a liquidity stake >= MIN_OFFER_STAKE.
    # ProphetX's array is loosely best-first but includes low-stake teaser
    # tiers that the UI hides. Filtering by stake matches the price a retail
    # user sees. Fall back to offers[0] only if nothing clears the threshold.
    qualifying = [o for o in offers
                  if (o.get("stake") or 0) >= MIN_OFFER_STAKE and o.get("odds") is not None]
    chosen = qualifying[0] if qualifying else offers[0]
    am = chosen.get("odds")
    if am is None:
        return None, False
    return {
        "decimal": round(american_to_decimal(int(am)), 4),
        "american": int(am),
        "offers": offers,  # kept for logging / future tier-aware pricing
    }, False


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def scrape_prophetx_sgp(verbose: bool = False):
    # Clear previous prices — missing combo means NA, not stale.
    clear_source("prophetx_direct")
    clear_source("prophetx_interpolated")

    print("Loading parlay lines from DuckDB...")
    parlay_lines = load_parlay_lines()
    if not parlay_lines:
        return
    print(f"  {len(parlay_lines)} games with lines")

    print("Initializing ProphetX session...")
    session = init_session()

    print("Fetching ProphetX MLB events...")
    try:
        px_events = fetch_prophetx_mlb_events(session)
    except Exception as e:
        print(f"  Failed to fetch events: {e}")
        return
    print(f"  {len(px_events)} ProphetX MLB events")

    # Drop live / already-started games
    from datetime import datetime, timezone
    now_utc = datetime.now(timezone.utc).isoformat()
    px_events = [e for e in px_events if (e.get("scheduled") or "") > now_utc]

    # Dedup: one event per (home, away, start-hour) — preserves doubleheaders
    seen = set()
    deduped = []
    for e in sorted(px_events, key=lambda x: x.get("scheduled") or ""):
        key = (e["px_home"], e["px_away"], (e.get("scheduled") or "")[:13])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(e)

    print("Matching teams via canonical_match...")
    matched = match_events(deduped, parlay_lines)
    print(f"  {len(matched)} matched games")
    if not matched:
        print("  No matches found.")
        return

    # ── Phase 1: fetch market trees + extract legs (parallel) ──
    print("Fetching market trees (parallel)...")
    t0 = time.time()
    legs_by_game = {}
    markets_by_game = {}  # raw market lists for integer-line fallback

    def _fetch(g):
        return g["game_id"], fetch_event_legs(session, g, verbose)

    with ThreadPoolExecutor(max_workers=PARALLEL_MARKETS) as pool:
        futures = {pool.submit(_fetch, g): g for g in matched}
        for fut in as_completed(futures):
            try:
                gid, (legs, markets) = fut.result()
                legs_by_game[gid] = legs
                markets_by_game[gid] = markets
            except Exception as e:
                g = futures[fut]
                print(f"  Error fetching markets for {g.get('home_team')}/{g.get('away_team')}: {e}")
    print(f"  Fetched {len(legs_by_game)} games in {time.time() - t0:.1f}s")

    # ── Phase 2: build combo items and RFQ (parallel) ──
    print("Pricing SGP combos (parallel RFQ)...")
    t1 = time.time()
    combo_items = []
    interpolated_results = []   # (game_id, period, combo_name, fair_prob)
    for g in matched:
        gid = g["game_id"]
        legs = legs_by_game.get(gid)
        if not legs:
            continue
        for period in ("fg", "f5"):
            p = legs.get(period, {})
            home = p.get("home_spread")
            away = p.get("away_spread")
            over = p.get("over")
            under = p.get("under")
            if not (home and away and over and under):
                # Exact-line lookup miss — try integer-line fallback.
                total_line = g[f"{period}_total_line"]
                spread_line = g[f"{period}_spread_line"]
                if total_line is None or spread_line is None:
                    continue
                markets = markets_by_game.get(gid, [])
                fallback = try_integer_fallback_px(
                    session, markets,
                    g["px_home_competitor_id"], g["px_away_competitor_id"],
                    g["px_event_id"],
                    spread_line, total_line, period,
                    verbose=verbose,
                )
                if fallback is None:
                    continue
                prefix = "" if period == "fg" else "F5 "
                COMBO_DISPLAY = {
                    "home_over":  "Home Spread + Over",
                    "home_under": "Home Spread + Under",
                    "away_over":  "Away Spread + Over",
                    "away_under": "Away Spread + Under",
                }
                for k, name in COMBO_DISPLAY.items():
                    interpolated_results.append((
                        gid, period, prefix + name, fallback["fair_probs"][k]
                    ))
                continue
            prefix = "" if period == "fg" else "F5 "
            for combo_name, sp, to in [
                ("Home Spread + Over",  home, over),
                ("Home Spread + Under", home, under),
                ("Away Spread + Over",  away, over),
                ("Away Spread + Under", away, under),
            ]:
                combo_items.append((gid, period, prefix + combo_name, sp, to))

    pricing_results = []
    event_id_by_game = {g["game_id"]: g["px_event_id"] for g in matched}
    auth_failures = {"count": 0}

    def _price(item):
        gid, period, combo_name, sp, to = item
        eid = event_id_by_game.get(gid)
        legs = [
            {"sportEventId": eid,
             "marketId": sp["marketId"], "outcomeId": sp["outcomeId"],
             "lineId":   sp["lineId"],   "line": sp["line"]},
            {"sportEventId": eid,
             "marketId": to["marketId"], "outcomeId": to["outcomeId"],
             "lineId":   to["lineId"],   "line": to["line"]},
        ]
        priced, auth_failed = submit_parlay_rfq(session, legs, verbose=verbose)
        if auth_failed:
            auth_failures["count"] += 1

        # Sanity filter: reject combos that exceed naive leg-product by > SANITY_MULT_RATIO.
        # Catches ProphetX's systematic F5-Over bug where the parlay pricer quotes
        # decimals 5-7x larger than independent multiply suggests.
        if priced is not None and sp.get("leg_am") is not None and to.get("leg_am") is not None:
            sp_dec = american_to_decimal(int(sp["leg_am"]))
            to_dec = american_to_decimal(int(to["leg_am"]))
            naive = sp_dec * to_dec
            if priced["decimal"] > naive * SANITY_MULT_RATIO:
                if verbose:
                    print(f"      SANITY-DROP {combo_name[:30]:<30} "
                          f"parlay={priced['decimal']:.2f} > naive×{SANITY_MULT_RATIO}="
                          f"{naive*SANITY_MULT_RATIO:.2f}")
                priced = None
        return gid, period, combo_name, priced

    with ThreadPoolExecutor(max_workers=PARALLEL_PRICING) as pool:
        futures = {pool.submit(_price, it): it for it in combo_items}
        for fut in as_completed(futures):
            try:
                pricing_results.append(fut.result())
            except Exception as e:
                it = futures[fut]
                print(f"  Error pricing {it[2]} ({it[0][:8]}): {e}")
                pricing_results.append((it[0], it[1], it[2], None))

    print(f"  Priced {len(pricing_results)} combos in {time.time() - t1:.1f}s")

    # ── Phase 2b: retry failed combos once ──
    failed = [r for r in pricing_results if r[3] is None]
    if failed:
        by_key = {(it[0], it[2]): it for it in combo_items}
        retry_items = [by_key[(r[0], r[2])] for r in failed if (r[0], r[2]) in by_key]
        if retry_items:
            print(f"  Retrying {len(retry_items)} failed combos...")
            time.sleep(1)
            with ThreadPoolExecutor(max_workers=PARALLEL_PRICING) as pool:
                fmap = {pool.submit(_price, it): it for it in retry_items}
                # Dict-based merge: latest result per (game, combo) wins.
                results_by_key = {(r[0], r[2]): r for r in pricing_results}
                for fut in as_completed(fmap):
                    try:
                        res = fut.result()
                        if res[3] is not None:
                            results_by_key[(res[0], res[2])] = res
                    except Exception as e:
                        it = fmap[fut]
                        print(f"  Retry error {it[2]} ({it[0][:8]}): {e}")
                pricing_results = list(results_by_key.values())

    # ── Phase 3: write ──
    ensure_table()
    all_rows = []
    game_lookup = {g["game_id"]: g for g in matched}
    priced_by_game = {}
    for gid, period, combo_name, priced in pricing_results:
        if priced:
            priced_by_game.setdefault(gid, []).append((period, combo_name, priced))

    for gid, combos in sorted(priced_by_game.items()):
        g = game_lookup[gid]
        print(f"\n  {g['away_team']} @ {g['home_team']}")
        # Per-game vig instrumentation — split by period
        vig_by_period = {"fg": [], "f5": []}
        for period, combo_name, priced in combos:
            dec = priced["decimal"]
            am = priced["american"]
            print(f"    {combo_name}: {dec:.4f} ({am:+d})  "
                  f"[{len(priced['offers'])} tiers]")
            vig_by_period[period].append(dec)
            all_rows.append({
                "game_id": gid, "combo": combo_name,
                "period": "FG" if period == "fg" else "F5",
                "bookmaker": "prophetx",
                "sgp_decimal": dec, "sgp_american": am,
                "source": "prophetx_direct",
            })
        for period, decs in vig_by_period.items():
            if len(decs) == 4:
                vig = sum(1 / d for d in decs)
                label = "FG" if period == "fg" else "F5"
                # Flag when measured vig materially exceeds our default (1.10).
                # Consistent readings >1.12 OR <1.08 mean PROPHETX_SGP_VIG_DEFAULT
                # needs re-tuning in either direction.
                flag = ""
                if vig > 1.12:
                    flag = "  <-- high vs default 1.10"
                elif vig < 1.08:
                    flag = "  <-- low vs default 1.10"
                print(f"    [{label} vig: {vig:.4f}]{flag}")

    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} ProphetX SGP odds in {time.time() - t0:.1f}s total")
        print(f"{'='*60}")
    else:
        # M2: don't fail silently. Loud diagnostic when zero rows written.
        reason = ""
        if auth_failures["count"] > 0:
            reason = (f" — {auth_failures['count']} RFQ calls returned 401/403 "
                      "(auth). Cookies from .prophetx_profile may be stale; log in "
                      "again via recon_prophetx_sgp.py to refresh.")
        elif not matched:
            reason = " — zero events matched; check MLB tournament name filter."
        else:
            reason = " — events matched but no RFQ returned priced offers. Try --verbose."
        print(f"\n!! ERROR: No ProphetX SGP odds collected{reason}")

    # Append rows for integer-line interpolation fallback
    interp_rows = []
    for gid, period, combo_name, fair_prob in interpolated_results:
        decimal = round(1.0 / fair_prob, 4)
        american = decimal_to_american(decimal)
        game = game_lookup.get(gid, {})
        print(f"\n  {game.get('away_team', '?')} @ {game.get('home_team', '?')} [interpolated]")
        print(f"    {combo_name}: {decimal:.4f} ({american:+d})")
        interp_rows.append({
            "game_id":      gid,
            "combo":        combo_name,
            "period":       "FG" if period == "fg" else "F5",
            "bookmaker":    "prophetx",
            "sgp_decimal":  decimal,
            "sgp_american": american,
            "source":       "prophetx_interpolated",
        })

    if interp_rows:
        upsert_sgp_odds(interp_rows)
        print(f"  Wrote {len(interp_rows)} interpolated ProphetX SGP odds")

    return all_rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="ProphetX MLB SGP Scraper")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    print("=" * 60)
    print("  PROPHETX MLB SGP SCRAPER (RFQ)")
    print("=" * 60)
    scrape_prophetx_sgp(verbose=args.verbose)


if __name__ == "__main__":
    main()
