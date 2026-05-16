#!/usr/bin/env python3
"""Novig MLB SGP Scraper — thin shim invoking the novig library.

Reads target_lines from MLB_DB (default: mlb_mm.duckdb; bot overrides via
MLB_SGP_DB_PATH env var). Calls mlb_sgp.novig.price_sgps() with the
periods configured via MLB_SGP_PERIODS env var (default: FG,F5).
Writes PricedRow results back to MLB_DB via mlb_sgp.db.upsert_priced_rows.

Legacy helpers (init_session, fetch_novig_mlb_events, match_events,
fetch_event_legs, submit_parlay, try_integer_fallback_nv,
_load_event_markets_query, load_parlay_lines, _find_outcome_in_spread /
_total, _empty_legs, _utc_bucket, _float_eq, decimal_to_american,
prob_to_decimal, _gql, and the MLB_EVENTS_QUERY / EVENT_MARKETS_QUERY /
EVENT_MARKETS_PATH / SPREAD_TYPE / TOTAL_TYPE / NOVIG_GRAPHQL /
NOVIG_PARLAY / RFQ_TIMEOUT / GQL_TIMEOUT / EVENT_WINDOW_HOURS /
SANITY_MULT_RATIO constants) are preserved in this file because
mlb_sgp/novig.py and mlb_sgp/novig_client.py import them lazily.
They stay here during the transition; a follow-up refactor can lift
them into the library module.

Final scraper in the DK/FD/PX/NV shim sequence — all four now route
through their library orchestrators with the same shim contract:
load_target_lines, price_sgps, upsert_priced_rows.

Note on dropped legacy orchestration
------------------------------------
The previous orchestrator (1) post-filtered by ``scheduled > now_utc``
and deduped events by ``(home, away, hour)``, (2) ran phase-1 market
fetches in a ``ThreadPoolExecutor(max_workers=PARALLEL_MARKETS)`` and
phase-2 RFQ pricing in another pool of size ``PARALLEL_PRICING``, and
(3) retried failed RFQ combos once after an 0.8-1.6s jitter. The shim
intentionally drops all three:

* Live-event filter: ``mlb_target_lines`` / ``mlb_parlay_lines`` is
  written upstream and already contains only the games we want to
  price, so a second client-side filter would be redundant.
  ``match_events`` (preserved) still maps Novig events to canonical
  game_ids, so any natural duplicates resolve to the same canonical
  game.
* Parallel pricing: the orchestrator is serial. For a 60s bot cadence
  this is acceptable; if smoke-testing reveals throughput issues we
  can add parallelization to the orchestrator later. Same call made
  for the DK/FD/PX shims.
* Retry-on-failure pass: the bot's cadence loop already invokes the
  scraper once per cycle, so the next cycle is the retry. Same call
  made for the DK/FD/PX shims.
"""

import os
import json
import sys
from datetime import datetime, timezone, timedelta
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

# Preserved helpers below (load_parlay_lines) still read mlb_parlay_lines via
# the legacy db.py helpers. The shim's main() uses mlb_sgp.db separately.
from db import MLB_DB, _connect_with_retry
from integer_line_derivation import is_integer_line, derive_fair_probs

# ---------------------------------------------------------------------------
# Novig API config
# ---------------------------------------------------------------------------
NOVIG_GRAPHQL = "https://api.novig.us/v1/graphql"
NOVIG_PARLAY  = "https://api.novig.us/nbx/v1/parlay/request/unauthenticated"

# Hasura query to list upcoming MLB events — written specifically for this scraper
# to get proper league-filtered results (LiveEventTicker_Query's upcoming branch
# has no league filter).
MLB_EVENTS_QUERY = """
query MLBEvents($start_gte: timestamptz!, $start_lte: timestamptz!) {
  event(where: {
    league: {_eq: "MLB"}, type: {_eq: "Game"},
    status: {_eq: "OPEN_PREGAME"},
    scheduled_start: {_gte: $start_gte, _lte: $start_lte}
  }, order_by: {scheduled_start: asc}) {
    id scheduled_start
    game {
      homeTeam { name symbol short_name }
      awayTeam { name symbol short_name }
    }
  }
}
"""

# Reuse Novig's own EventMarkets_Query — captured from recon and committed
# to disk so the scraper bootstraps cleanly on a fresh checkout. Returns the
# full market tree (with outcomes + fragments). We only use the markets[] array.
EVENT_MARKETS_QUERY = None  # populated at runtime from EVENT_MARKETS_PATH
EVENT_MARKETS_PATH = _THIS_DIR / "novig_event_markets_query.json"
# Legacy cache path (written by older scraper versions); kept for compat.
_LEGACY_CACHE_PATH = _THIS_DIR / ".novig_event_markets_query.json"

# Market type names (Novig uses SPREAD_1H / TOTAL_1H for F5)
SPREAD_TYPE = {"fg": "SPREAD",    "f5": "SPREAD_1H"}
TOTAL_TYPE  = {"fg": "TOTAL",     "f5": "TOTAL_1H"}

PARALLEL_MARKETS = 4
PARALLEL_PRICING = 4
RFQ_TIMEOUT = 15
GQL_TIMEOUT = 20

EVENT_WINDOW_HOURS = 48   # how far ahead to look for upcoming games

# Sanity filter: drop combos where the parlay decimal exceeds this multiplier
# times the naive independent multiply of the two legs. Legitimate anti-correlation
# tops out around ~1.15×; 1.5× gives headroom while still catching systematic
# mispricings (5×+ ratios observed on ProphetX — same pattern could occur here).
SANITY_MULT_RATIO = 1.5


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1)))


def prob_to_decimal(p: float) -> float:
    return 1.0 / p if p > 0 else float("inf")


def _utc_bucket(ts) -> str:
    """UTC YYYY-MM-DDTHH bucket for doubleheader matching."""
    if not ts:
        return ""
    if hasattr(ts, "year"):
        return f"{ts.year:04d}-{ts.month:02d}-{ts.day:02d}T{ts.hour:02d}"
    return ts[:13] if len(ts) >= 13 else ""


# ---------------------------------------------------------------------------
# GraphQL loading: EventMarkets_Query
# ---------------------------------------------------------------------------
def _load_event_markets_query() -> str:
    """Load the captured EventMarkets_Query JSON.

    Tries (in order):
      1. The committed canonical file `novig_event_markets_query.json`
         — what ships on main and works on a fresh checkout.
      2. The legacy hidden cache `.novig_event_markets_query.json`
         — written by older scraper versions; kept for compatibility.
      3. The recon JSON `recon_novig_sgp.json` (if available locally),
         extracting the captured GraphQL post_data.
    """
    global EVENT_MARKETS_QUERY
    if EVENT_MARKETS_QUERY is not None:
        return EVENT_MARKETS_QUERY

    for path in (EVENT_MARKETS_PATH, _LEGACY_CACHE_PATH):
        if path.exists():
            try:
                EVENT_MARKETS_QUERY = path.read_text()
                return EVENT_MARKETS_QUERY
            except Exception:
                pass

    # Last-resort fallback: extract from recon JSON if it's around
    recon_path = _THIS_DIR / "recon_novig_sgp.json"
    if recon_path.exists():
        try:
            phases = json.loads(recon_path.read_text())
            for ph in phases:
                for r in ph.get("requests", []):
                    if "/v1/graphql" not in r.get("url", ""):
                        continue
                    pd = r.get("post_data") or ""
                    if "EventMarkets_Query" in pd[:120]:
                        EVENT_MARKETS_QUERY = pd
                        return EVENT_MARKETS_QUERY
        except Exception as e:
            print(f"  (could not parse recon JSON: {e})")

    raise RuntimeError(
        "EventMarkets_Query text unavailable. Expected at "
        f"{EVENT_MARKETS_PATH}. The committed canonical query file is "
        "missing — recover from /tmp/novig_EventMarkets_Query.json or "
        "re-run recon_novig_sgp.py."
    )


# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
def init_session() -> cffi_requests.Session:
    """No auth required — plain session with Chrome impersonation for safety."""
    session = cffi_requests.Session(impersonate="chrome")
    # Warm up: fetch marketing page to seed any anti-bot cookies (Cloudflare
    # __cf_bm, etc.) before we start hammering the API.
    try:
        session.get("https://www.novig.us", timeout=10)
    except Exception:
        pass  # non-fatal — API is the point, not the landing page
    return session


def _gql(session, body: str) -> dict:
    """POST a raw GraphQL payload string. Returns parsed JSON (or raises)."""
    resp = session.post(
        NOVIG_GRAPHQL, data=body,
        headers={"Content-Type": "application/json"},
        timeout=GQL_TIMEOUT,
    )
    resp.raise_for_status()
    return resp.json()


# ---------------------------------------------------------------------------
# Event discovery
# ---------------------------------------------------------------------------
def fetch_novig_mlb_events(session) -> list[dict]:
    """List upcoming MLB events within the event window. Returns list of
    {nv_event_id, nv_home, nv_away, nv_home_sym, nv_away_sym, scheduled}."""
    now = datetime.now(timezone.utc)
    cutoff = (now + timedelta(hours=EVENT_WINDOW_HOURS)).isoformat()
    payload = json.dumps({
        "query": MLB_EVENTS_QUERY,
        "variables": {"start_gte": now.isoformat(), "start_lte": cutoff},
    })
    data = _gql(session, payload)
    events_raw = (data.get("data") or {}).get("event") or []

    out = []
    for e in events_raw:
        g = e.get("game") or {}
        ht = g.get("homeTeam") or {}
        at = g.get("awayTeam") or {}
        if not (ht.get("name") and at.get("name")):
            continue
        out.append({
            "nv_event_id": e.get("id"),
            "nv_home":     ht.get("name"),
            "nv_away":     at.get("name"),
            "nv_home_sym": ht.get("symbol") or ht.get("short_name"),
            "nv_away_sym": at.get("symbol") or at.get("short_name"),
            "scheduled":   e.get("scheduled_start", ""),
        })
    return out


# ---------------------------------------------------------------------------
# Canonical game matching
# ---------------------------------------------------------------------------
def load_parlay_lines() -> dict:
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
        return {
            row[0]: {
                "home_team": row[1], "away_team": row[2],
                "fg_spread_line": row[3], "fg_total_line": row[4],
                "f5_spread_line": row[5], "f5_total_line": row[6],
                "commence_time": row[7],
            }
            for row in rows
        }
    finally:
        con.close()


def match_events(nv_events: list[dict], parlay_lines: dict) -> list[dict]:
    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")
    matched = []
    for evt in nv_events:
        resolved = resolve_team_names(
            evt["nv_away"], evt["nv_home"], team_dict, canonical_games,
        )
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        canon_away, canon_home = resolved
        nv_bucket = _utc_bucket(evt["scheduled"])
        for gid, lines in parlay_lines.items():
            if lines["home_team"] != canon_home or lines["away_team"] != canon_away:
                continue
            pl_bucket = _utc_bucket(lines.get("commence_time", ""))
            if pl_bucket and nv_bucket and pl_bucket != nv_bucket:
                continue
            matched.append({
                "game_id": gid,
                "nv_event_id": evt["nv_event_id"],
                "home_team": canon_home, "away_team": canon_away,
                "nv_home_sym": evt["nv_home_sym"], "nv_away_sym": evt["nv_away_sym"],
                "fg_spread_line": lines["fg_spread_line"],
                "fg_total_line":  lines["fg_total_line"],
                "f5_spread_line": lines["f5_spread_line"],
                "f5_total_line":  lines["f5_total_line"],
            })
            break
    return matched


# ---------------------------------------------------------------------------
# Market tree -> per-period legs
# ---------------------------------------------------------------------------
def _find_outcome_in_spread(market, home_sym, away_sym) -> tuple[dict | None, dict | None]:
    """Return (home_leg, away_leg) from a SPREAD or SPREAD_1H market.
    Each leg is {"id": uuid, "available": implied_prob}. `available` used for
    the naive-multiply sanity check in _price. Identifies via competitor.symbol.
    """
    home_leg = None
    away_leg = None
    for o in market.get("outcomes") or []:
        comp = o.get("competitor") or {}
        sym = comp.get("symbol")
        leg = {"id": o.get("id"), "available": o.get("available")}
        if sym == home_sym:
            home_leg = leg
        elif sym == away_sym:
            away_leg = leg
    return home_leg, away_leg


def _find_outcome_in_total(market) -> tuple[dict | None, dict | None]:
    """Return (over_leg, under_leg) from a TOTAL market.
    Each leg is {"id": uuid, "available": implied_prob}."""
    over_leg = None
    under_leg = None
    for o in market.get("outcomes") or []:
        desc = (o.get("description") or "").lower()
        leg = {"id": o.get("id"), "available": o.get("available")}
        if desc.startswith("over "):
            over_leg = leg
        elif desc.startswith("under "):
            under_leg = leg
    return over_leg, under_leg


def _float_eq(a, b, eps=1e-6) -> bool:
    if a is None or b is None:
        return False
    try:
        return abs(float(a) - float(b)) < eps
    except (TypeError, ValueError):
        return False


def fetch_event_legs(session, game: dict, verbose: bool = False) -> tuple[dict, list]:
    """Fetch market tree for one event, extract the 4 outcome UUIDs per period.

    Returns (legs, markets) where legs is the per-period dict and markets is the
    raw market list from the API (needed by the integer-line fallback helper).
    """
    query_text = _load_event_markets_query()
    # The captured post_data is already a full JSON blob including operationName,
    # variables, query — we need to rewrite the eventId variable.
    q_obj = json.loads(query_text)
    q_obj["variables"]["eventId"] = game["nv_event_id"]
    try:
        data = _gql(session, json.dumps(q_obj))
    except Exception as e:
        if verbose:
            print(f"      EventMarkets error for {game['nv_event_id']}: {e}")
        return _empty_legs(), []

    ev = ((data.get("data") or {}).get("event") or [{}])[0]
    markets = ev.get("markets") or []

    home_sym = game["nv_home_sym"]
    away_sym = game["nv_away_sym"]
    out = _empty_legs()

    for period in ("fg", "f5"):
        target_spread = game[f"{period}_spread_line"]
        target_total  = game[f"{period}_total_line"]
        if target_spread is None or target_total is None:
            continue

        # --- Spread market (type = SPREAD or SPREAD_1H) ---
        # Prefer is_consensus=true when multiple markets share the same strike
        # (review C2). For a given (type, strike), a non-consensus duplicate
        # could be a stale/promo market we don't want.
        spread_type = SPREAD_TYPE[period]
        spread_matches = [m for m in markets
                          if m.get("type") == spread_type
                          and _float_eq(m.get("strike"), target_spread)]
        spread_mkt = next((m for m in spread_matches if m.get("is_consensus") is True),
                          spread_matches[0] if spread_matches else None)
        if spread_mkt is None and period == "f5" and verbose:
            # Diagnostic only: Novig exposes only the main F5 line, so an unmatched
            # target means our fg_pipeline's F5 line moved relative to Novig's.
            cand = [m for m in markets if m.get("type") == spread_type]
            if cand:
                print(f"      [F5] {game['game_id'][:8]}: target spread "
                      f"{target_spread} != Novig's F5 line {cand[0].get('strike')} "
                      "— skipping for strict match")
        if spread_mkt:
            home_leg, away_leg = _find_outcome_in_spread(spread_mkt, home_sym, away_sym)
            if home_leg and away_leg:
                out[period]["home_spread"] = home_leg
                out[period]["away_spread"] = away_leg

        # --- Total market (type = TOTAL or TOTAL_1H) ---
        total_type = TOTAL_TYPE[period]
        total_matches = [m for m in markets
                         if m.get("type") == total_type
                         and _float_eq(m.get("strike"), target_total)]
        total_mkt = next((m for m in total_matches if m.get("is_consensus") is True),
                         total_matches[0] if total_matches else None)
        if total_mkt:
            over_leg, under_leg = _find_outcome_in_total(total_mkt)
            if over_leg and under_leg:
                out[period]["over"]  = over_leg
                out[period]["under"] = under_leg

        if verbose:
            missing = [k for k, v in out[period].items() if v is None]
            if missing:
                print(f"      [{period.upper()}] {game['game_id'][:8]}: "
                      f"missing {missing}  (target spread={target_spread}, total={target_total})")

    return out, markets


def _empty_legs():
    return {
        "fg": {"home_spread": None, "away_spread": None, "over": None, "under": None},
        "f5": {"home_spread": None, "away_spread": None, "over": None, "under": None},
    }


# ---------------------------------------------------------------------------
# Integer-line fallback helper
# ---------------------------------------------------------------------------

def try_integer_fallback_nv(
    session,
    markets: list,
    home_sym: str,
    away_sym: str,
    spread_line: float,
    total_line: float,
    period: str,
    verbose: bool = False,
):
    """NV: when WZ quotes an integer total that Novig's market tree doesn't have
    at the exact strike, look up the two adjacent half-point alt total markets
    (total_line - 0.5 and total_line + 0.5), price all 8 combos, and call
    derive_fair_probs to recover the integer-line fair probabilities.

    Args:
        markets: raw market list from fetch_event_legs (the full API response).
        home_sym / away_sym: Novig competitor symbols for spread leg lookup.
        spread_line: target spread from mlb_parlay_lines (home perspective).
        total_line: the integer total line from mlb_parlay_lines.
        period: 'fg' or 'f5' — selects the correct market type strings.

    Returns a derive_fair_probs result dict or None.
    """
    if not is_integer_line(total_line):
        return None

    lo, hi = total_line - 0.5, total_line + 0.5

    # --- Spread legs: reuse the exact-match spread that was already found ---
    spread_type = SPREAD_TYPE[period]
    spread_matches = [m for m in markets
                      if m.get("type") == spread_type
                      and _float_eq(m.get("strike"), spread_line)]
    spread_mkt = next((m for m in spread_matches if m.get("is_consensus") is True),
                      spread_matches[0] if spread_matches else None)
    if spread_mkt is None:
        if verbose:
            print(f"      integer fallback: no spread market at {spread_line}")
        return None
    home_leg, away_leg = _find_outcome_in_spread(spread_mkt, home_sym, away_sym)
    if not (home_leg and away_leg):
        if verbose:
            print(f"      integer fallback: spread legs missing for {spread_line}")
        return None

    # --- Adjacent total markets ---
    total_type = TOTAL_TYPE[period]

    def _get_total_legs(strike):
        mkt_list = [m for m in markets
                    if m.get("type") == total_type and _float_eq(m.get("strike"), strike)]
        mkt = next((m for m in mkt_list if m.get("is_consensus") is True),
                   mkt_list[0] if mkt_list else None)
        if mkt is None:
            return None, None
        return _find_outcome_in_total(mkt)

    over_lo, under_lo = _get_total_legs(lo)
    over_hi, under_hi = _get_total_legs(hi)

    if not all([over_lo, under_lo, over_hi, under_hi]):
        if verbose:
            missing = []
            if not over_lo:  missing.append(f"over {lo}")
            if not under_lo: missing.append(f"under {lo}")
            if not over_hi:  missing.append(f"over {hi}")
            if not under_hi: missing.append(f"under {hi}")
            print(f"      integer fallback: missing adjacent alt legs {missing} for total={total_line}")
        return None

    def _price(sp_leg, tot_leg):
        priced, _ = submit_parlay(session, [sp_leg["id"], tot_leg["id"]], verbose=verbose)
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
# Parlay RFQ
# ---------------------------------------------------------------------------
def submit_parlay(session, outcome_ids: list[str],
                  verbose: bool = False) -> tuple[dict | None, bool]:
    """POST 2+ outcome UUIDs to the Novig parlay endpoint.

    Returns (priced, auth_failed) where priced is
    {'decimal','american','price_str','status','raw_offers'} or None.
    """
    payload = {"outcomes": [{"id": oid} for oid in outcome_ids], "boostId": None}
    try:
        resp = session.post(
            NOVIG_PARLAY, json=payload,
            headers={"Content-Type": "application/json"},
            timeout=RFQ_TIMEOUT,
        )
    except Exception as e:
        if verbose:
            print(f"      parlay error: {e}")
        return None, False

    if resp.status_code in (401, 403):
        if verbose:
            print(f"      parlay {resp.status_code} — unexpected for /unauthenticated endpoint")
        return None, True
    if resp.status_code not in (200, 201):
        if verbose:
            preview = resp.text[:200] if resp.text else ""
            print(f"      parlay HTTP {resp.status_code}: {preview}")
        return None, False

    try:
        offers = resp.json()
    except Exception:
        return None, False

    if not offers or not isinstance(offers, list):
        return None, False

    top = offers[0]
    price_str = top.get("price")
    if price_str is None:
        return None, False
    try:
        p = float(price_str)
    except (TypeError, ValueError):
        return None, False
    if not (0 < p < 1):
        return None, False

    dec = round(prob_to_decimal(p), 4)
    return {
        "decimal":    dec,
        "american":   decimal_to_american(dec),
        "price_str":  price_str,
        "status":     top.get("status"),
        "raw_offers": offers,
    }, False


# ---------------------------------------------------------------------------
# Shim entry point
# ---------------------------------------------------------------------------
#
# Everything below replaces the legacy scrape_novig_sgp() + main() pair.
# The shim delegates orchestration to mlb_sgp.novig.price_sgps(), which
# composes the helpers above. Helpers stay defined at module level because
# mlb_sgp/novig.py and mlb_sgp/novig_client.py import them lazily.


def main():
    """Entry point: load targets, price SGPs, write rows.

    Path / import gymnastics: this script can be launched in three modes —
    (1) `python mlb_sgp/scraper_novig_sgp.py` from the repo root
        (test_sgp_regression's shim test path),
    (2) `python scraper_novig_sgp.py` from inside mlb_sgp/
        (legacy direct-invocation path, also how kalshi_mlb_rfq spawns it),
    (3) `python -m mlb_sgp.scraper_novig_sgp` (package mode).

    For (1) and (2) we need `mlb_sgp` importable as a top-level package so
    `mlb_sgp.novig` resolves. Make that work by ensuring the repo root
    is on sys.path before importing the orchestrator + db modules.

    Live-event filter / parallel pricing / retry pass: the legacy
    ``scrape_novig_sgp`` post-filtered Novig events by ``scheduled > now_utc``,
    deduped events by ``(home, away, hour)``, fanned out market fetches and
    RFQ pricing over a ThreadPoolExecutor, and retried failed combos once
    after an 0.8-1.6s jitter. The shim deliberately drops all three:
    ``mlb_target_lines`` / ``mlb_parlay_lines`` is written upstream and
    only contains the games we want to price (and ``match_events`` still
    handles canonical resolution so doubleheaders / dupes collapse to the
    same game_id); the orchestrator is serial because the bot's 60s cadence
    loop tolerates it and the next cycle is the natural retry. Matches
    DK/FD/PX shim shape.
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

    from mlb_sgp import novig
    print(f"  NV shim: {len(targets)} target lines, periods={periods}")
    rows = novig.price_sgps(targets, periods=periods, verbose=False)
    print(f"  NV shim: priced {len(rows)} rows")

    # Wipe both source labels so stale rows from a previous run never linger.
    # novig.price_sgps emits both ``novig_direct`` (main RFQ path) and
    # ``novig_interpolated`` (integer-line fallback path) rows; clearing
    # both preserves the "fresh prices only" invariant.
    db.clear_source("novig_direct", db_path=db_path)
    db.clear_source("novig_interpolated", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
