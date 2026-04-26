#!/usr/bin/env python3
"""
Novig MLB SGP Scraper

Fetches Same Game Parlay (SGP) prices from Novig's GraphQL + REST API, both
full-game and first-5-innings (Novig types F5 as *_1H). Writes to mlb_sgp_odds
with source='novig_direct'.

Architectural notes (see plan doc + recon memory for details):
    - All API calls go to api.novig.us. No auth required — the parlay endpoint
      is literally named /unauthenticated.
    - Novig pulls leg prices from DraftKings (every observed leg returned
      vendor="DRAFTKINGS"). The combined parlay price IS correlation-adjusted
      (empirically 1.09x-1.12x more generous than naive independent multiply
      for fav+over combos). What Novig adds to the blend = the delta between
      Novig's correlation model and DK's calculateBets, applied to the same
      underlying leg prices.
    - Market types: SPREAD (FG), TOTAL (FG), SPREAD_1H (F5 spread),
      TOTAL_1H (F5 total). Novig offers only the main F5 line, no alts.
    - strike field is signed from home's perspective — matches our
      mlb_parlay_lines.fg_spread convention exactly.
    - Response prices are implied-probability strings like "0.26600"; we
      convert to decimal via 1/p and to American odds for storage.

Usage:
    cd mlb_sgp
    venv/bin/python3 scraper_novig_sgp.py           # all games
    venv/bin/python3 scraper_novig_sgp.py --verbose  # show details
"""

import argparse
import json
import sys
import time
from datetime import datetime, timezone, timedelta
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

from db import ensure_table, upsert_sgp_odds, clear_source, MLB_DB

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
    con = duckdb.connect(str(MLB_DB), read_only=True)
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


def fetch_event_legs(session, game: dict, verbose: bool = False) -> dict:
    """Fetch market tree for one event, extract the 4 outcome UUIDs per period."""
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
        return _empty_legs()

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

    return out


def _empty_legs():
    return {
        "fg": {"home_spread": None, "away_spread": None, "over": None, "under": None},
        "f5": {"home_spread": None, "away_spread": None, "over": None, "under": None},
    }


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
# Main
# ---------------------------------------------------------------------------
def scrape_novig_sgp(verbose: bool = False):
    clear_source("novig_direct")

    print("Loading parlay lines from DuckDB...")
    parlay_lines = load_parlay_lines()
    if not parlay_lines:
        return
    print(f"  {len(parlay_lines)} games with lines")

    print("Initializing Novig session (anonymous)...")
    session = init_session()

    print("Fetching Novig MLB events...")
    try:
        nv_events = fetch_novig_mlb_events(session)
    except Exception as e:
        print(f"  Failed to fetch events: {e}")
        return
    print(f"  {len(nv_events)} upcoming MLB events from Novig")

    # De-dupe by (home, away, start-hour) — preserves doubleheaders
    seen = set()
    deduped = []
    for e in sorted(nv_events, key=lambda x: x.get("scheduled") or ""):
        key = (e["nv_home"], e["nv_away"], (e.get("scheduled") or "")[:13])
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

    # ── Phase 1: fetch market trees + extract outcome UUIDs ──
    print("Fetching market trees (parallel)...")
    t0 = time.time()
    legs_by_game = {}

    def _fetch(g):
        return g["game_id"], fetch_event_legs(session, g, verbose)

    with ThreadPoolExecutor(max_workers=PARALLEL_MARKETS) as pool:
        futures = {pool.submit(_fetch, g): g for g in matched}
        for fut in as_completed(futures):
            try:
                gid, legs = fut.result()
                legs_by_game[gid] = legs
            except Exception as e:
                g = futures[fut]
                print(f"  Error fetching markets for "
                      f"{g.get('home_team')}/{g.get('away_team')}: {e}")
    print(f"  Fetched {len(legs_by_game)} games in {time.time() - t0:.1f}s")

    # ── Phase 2: build combo items + parlay RFQ ──
    print("Pricing SGP combos (parallel)...")
    t1 = time.time()
    combo_items = []
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
    auth_failures = {"count": 0}

    def _price(item):
        gid, period, combo_name, sp, to = item
        priced, auth_failed = submit_parlay(session, [sp["id"], to["id"]], verbose=verbose)
        if auth_failed:
            auth_failures["count"] += 1

        # Sanity filter: reject if parlay decimal > naive leg-product × SANITY_MULT_RATIO.
        # Same logic as ProphetX — catches systematic mispricings where a book's
        # correlation model returns decimals multiples larger than independent-multiply.
        if priced is not None:
            sp_av = sp.get("available")
            to_av = to.get("available")
            if sp_av and to_av and sp_av > 0 and to_av > 0:
                naive = (1.0 / sp_av) * (1.0 / to_av)
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
            # Jitter between 0.8s-1.6s to avoid a thundering 3-worker burst exactly
            # 1s after the original batch, which is worst-case for IP-based rate limits.
            import random
            time.sleep(random.uniform(0.8, 1.6))
            with ThreadPoolExecutor(max_workers=PARALLEL_PRICING) as pool:
                fmap = {pool.submit(_price, it): it for it in retry_items}
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
        vig_by_period = {"fg": [], "f5": []}
        for period, combo_name, priced in combos:
            dec = priced["decimal"]
            am = priced["american"]
            print(f"    {combo_name}: {dec:.4f} ({am:+d})  "
                  f"p={priced['price_str']}  {priced.get('status','?')}")
            vig_by_period[period].append(dec)
            all_rows.append({
                "game_id":      gid,
                "combo":        combo_name,
                "period":       "FG" if period == "fg" else "F5",
                "bookmaker":    "novig",
                "sgp_decimal":  dec,
                "sgp_american": am,
                "source":       "novig_direct",
            })
        for period, decs in vig_by_period.items():
            if len(decs) == 4:
                vig = sum(1 / d for d in decs)
                label = "FG" if period == "fg" else "F5"
                flag = ""
                if vig > 1.12:
                    flag = "  <-- high vs default 1.10"
                elif vig < 1.08:
                    flag = "  <-- low vs default 1.10"
                print(f"    [{label} vig: {vig:.4f}]{flag}")

    if all_rows:
        upsert_sgp_odds(all_rows)
        print(f"\n{'='*60}")
        print(f"  Wrote {len(all_rows)} Novig SGP odds in {time.time() - t0:.1f}s total")
        print(f"{'='*60}")
    else:
        if auth_failures["count"] > 0:
            reason = (f" — {auth_failures['count']} RFQ calls returned 401/403. "
                      "This is unexpected for /unauthenticated — Novig may be "
                      "rate-limiting or requiring auth now.")
        elif not matched:
            reason = " — zero events matched; check MLB league filter."
        else:
            reason = " — events matched but no parlay prices returned. Try --verbose."
        print(f"\n!! ERROR: No Novig SGP odds collected{reason}")

    return all_rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Novig MLB SGP Scraper")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    print("=" * 60)
    print("  NOVIG MLB SGP SCRAPER (GraphQL + parlay RFQ)")
    print("=" * 60)
    scrape_novig_sgp(verbose=args.verbose)


if __name__ == "__main__":
    main()
