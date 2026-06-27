#!/usr/bin/env python3
"""ProphetX MLB SGP Scraper — thin shim invoking the prophetx library.

Reads target_lines from MLB_DB (default: mlb_mm.duckdb; bot overrides via
MLB_SGP_DB_PATH env var). Calls mlb_sgp.prophetx.price_sgps() with the
periods configured via MLB_SGP_PERIODS env var (default: FG,F5).
Writes PricedRow results back to MLB_DB via mlb_sgp.db.upsert_priced_rows.

Legacy helpers (init_session, fetch_prophetx_mlb_events, match_events,
_find_market, _pick_selection, _verify_competitor_ids, fetch_event_legs,
submit_parlay_rfq, try_integer_fallback_px, MARKET_NAMES / NAME_ALIASES /
MIN_OFFER_STAKE / SANITY_MULT_RATIO constants) are preserved in this file
because the new prophetx orchestrator + prophetx_client import them lazily.
They stay here during the transition; a follow-up refactor can lift them
into the library module.

PX-specific SANITY_MULT_RATIO defense (F5-Over bug) lives in the prophetx.py
orchestrator now. The integer-fallback path is retained in this module for
backward compat with legacy callers but is no longer wired to any emitter —
the orchestrator does not invoke it because PX has produced zero
`prophetx_interpolated` rows in practice (Task 11 design decision).

Note on the legacy live-event filter
------------------------------------
The previous orchestrator filtered ``scheduled > now_utc`` and deduped
events by ``(home, away, hour)`` to suppress live (in-progress) games.
The shim intentionally drops both: ``mlb_target_lines`` / ``mlb_parlay_lines``
are written upstream and already contain only the games we want to price,
so a second client-side filter would be redundant. Matches DK/FD's shim
shape and behavior.
"""

import os
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
        "moneyline": "Moneyline",
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
    so caller can extract outcomeId (selection['id']), lineID, line.

    Some PX market lines have `selections: [null, null]` and put the actual
    outcomes under a separate `outcomes` field — skip None sides defensively.
    """
    for line_grp in market.get("marketLines") or []:
        for side in line_grp.get("selections") or []:
            if side is None:
                continue
            for sel in side:
                if sel is None:
                    continue
                if predicate(sel):
                    return sel
        # Fallback: this line_grp uses `outcomes` (flat list of selections)
        # instead of `selections` (nested per-side). Try that shape too.
        for sel in line_grp.get("outcomes") or []:
            if sel is None:
                continue
            if predicate(sel):
                return sel
    # Some markets (moneyline) have flat selections at the top level
    for side in market.get("selections") or []:
        if side is None:
            continue
        for sel in side:
            if sel is None:
                continue
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
# Shim entry point
# ---------------------------------------------------------------------------
#
# Everything below replaces the legacy scrape_prophetx_sgp() + main() pair.
# The shim delegates orchestration to mlb_sgp.prophetx.price_sgps(), which
# composes the helpers above. Helpers stay defined at module level because
# mlb_sgp/prophetx.py and mlb_sgp/prophetx_client.py import them lazily.


def main():
    """Entry point: load targets, price SGPs, write rows.

    Path / import gymnastics: this script can be launched in three modes —
    (1) `python mlb_sgp/scraper_prophetx_sgp.py` from the repo root
        (test_sgp_regression's shim test path),
    (2) `python scraper_prophetx_sgp.py` from inside mlb_sgp/
        (legacy direct-invocation path, also how kalshi_mlb_rfq spawns it),
    (3) `python -m mlb_sgp.scraper_prophetx_sgp` (package mode).

    For (1) and (2) we need `mlb_sgp` importable as a top-level package so
    `mlb_sgp.prophetx` resolves. Make that work by ensuring the repo root
    is on sys.path before importing the orchestrator + db modules.

    Live-event filtering: the legacy ``scrape_prophetx_sgp`` post-filtered
    PX events by ``scheduled > now_utc`` and deduped by ``(home, away, hour)``.
    The shim deliberately drops both — ``mlb_target_lines`` /
    ``mlb_parlay_lines`` is written upstream and only contains the games
    we want to price, so the post-filter would be redundant.
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

    from mlb_sgp import prophetx
    print(f"  PX shim: {len(targets)} target lines, periods={periods}", flush=True)

    # Incremental batched writes (2026-06-10). PX prices line-by-line over the
    # RFQ endpoint (~250-350s for a full slate) and used to write ONCE at the
    # end — so when an orchestrator's scrape deadline killed it, ZERO rows
    # landed (FD/Novig write incrementally and survive the same kill). Batch
    # the targets and upsert after each batch so a deadline kill keeps
    # everything priced so far. One client is built up front and reused so
    # batching adds no session overhead.
    #
    # Wipe both source labels ONCE up front so stale rows from a previous run
    # never linger ("fresh prices only" invariant — partial fresh > stale).
    # The orchestrator currently emits only _direct rows; clearing
    # _interpolated too preserves the invariant in case the legacy
    # integer-fallback path is ever re-wired.
    db.clear_source("prophetx_direct", db_path=db_path)
    db.clear_source("prophetx_interpolated", db_path=db_path)

    batch_size = int(os.environ.get("PX_WRITE_BATCH", "50"))
    client = prophetx.ProphetXClient(verbose=False)
    total = 0
    for i in range(0, len(targets), batch_size):
        batch = targets[i:i + batch_size]
        rows = prophetx.price_sgps(batch, periods=periods, client=client,
                                   verbose=False)
        if rows:
            db.upsert_priced_rows(rows, db_path=db_path)
            total += len(rows)
        print(f"  PX shim: batch {i // batch_size + 1} — {len(rows)} rows "
              f"(total {total})", flush=True)
    print(f"  PX shim: priced {total} rows", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
