"""Kalshi MLB RFQ Bot — autonomous taker daemon."""

import argparse
import json
import os
import signal
import subprocess
import threading
import time
import uuid
from datetime import datetime, timedelta, timezone

import duckdb
import numpy as np
import pandas as pd

from kalshi_mlb_rfq import (
    auth_client, combo_enumerator, config, db, ev_calc,
    fair_value, kelly, notify, rfq_client, risk,
)
from kalshi_mlb_rfq.config import (
    ANSWER_KEY_DB, KILL_FILE, MAX_BOOK_STALENESS_SEC,
)

VERSION = "0.1.0"

_running = threading.Event()
_running.set()
ACCEPT_LOCK = threading.Lock()

# Test-only override; otherwise use config.MAX_LIVE_RFQS.
MAX_LIVE_RFQS_OVERRIDE: int | None = None

# Track positions-API health across calls.
_POSITIONS_API_FAIL_COUNT = 0


def _signal_handler(_sig, _frame):
    _running.clear()


def _max_live_rfqs() -> int:
    return MAX_LIVE_RFQS_OVERRIDE if MAX_LIVE_RFQS_OVERRIDE is not None else config.MAX_LIVE_RFQS


def _record_positions_api_result(success: bool):
    global _POSITIONS_API_FAIL_COUNT
    if success:
        _POSITIONS_API_FAIL_COUNT = 0
    else:
        _POSITIONS_API_FAIL_COUNT += 1


# ------------------------------------------------------------------------ #
# Pipeline-refresh subprocess                                              #
# ------------------------------------------------------------------------ #

def _run_pipeline():
    """Trigger the MLB R answer-key pipeline. Notify on failure; do not halt RFQs."""
    cmd = ["Rscript", str(ANSWER_KEY_DB.parent / "MLB Answer Key" / "MLB.R")]
    try:
        proc = subprocess.run(cmd, cwd=str(ANSWER_KEY_DB.parent),
                              timeout=240, capture_output=True, text=True)
        if proc.returncode != 0:
            notify.halt("pipeline_refresh_failed",
                        detail=f"exit={proc.returncode}; stderr={proc.stderr[:500]}")
    except subprocess.TimeoutExpired:
        notify.halt("pipeline_refresh_failed", detail="timeout")
    except Exception as e:
        notify.halt("pipeline_refresh_failed", detail=str(e))


# ------------------------------------------------------------------------ #
# Phantom RFQ cleanup at startup                                           #
# ------------------------------------------------------------------------ #

def _phantom_rfq_cleanup():
    try:
        kalshi_open = rfq_client.list_open_rfqs(config.KALSHI_USER_ID or "")
    except Exception as e:
        print(f"  startup: list_open_rfqs failed: {e}", flush=True)
        return

    if not kalshi_open:
        print("  startup: no open RFQs on Kalshi", flush=True)
        return

    with db.connect(read_only=True) as con:
        ours = {r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()}

    for rfq in kalshi_open:
        rid = rfq.get("id")
        if rid and rid not in ours:
            try:
                rfq_client.delete_rfq(rid)
                print(f"  startup: cancelled phantom rfq {rid}", flush=True)
            except Exception as e:
                print(f"  startup: failed to cancel phantom {rid}: {e}", flush=True)


# ------------------------------------------------------------------------ #
# Answer-key DB helpers (read-only)                                        #
# ------------------------------------------------------------------------ #

def _samples_generated_at() -> datetime | None:
    if not ANSWER_KEY_DB.exists():
        return None
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT generated_at FROM mlb_samples_meta ORDER BY generated_at DESC LIMIT 1"
        ).fetchone()
    except duckdb.CatalogException:
        return None
    finally:
        con.close()
    return row[0] if row else None


def _commence_time_for_game(game_id: str) -> datetime | None:
    """Pull commence_time from mlb_parlay_lines (the canonical source for it)."""
    if not ANSWER_KEY_DB.exists():
        return None
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT commence_time FROM mlb_parlay_lines WHERE game_id=? LIMIT 1",
            [game_id]
        ).fetchone()
    except duckdb.CatalogException:
        return None
    finally:
        con.close()
    return row[0] if row else None


# 3-letter Kalshi team code → mlb_parlay_lines.home_team / away_team canonical name.
# Kalshi uses 3-letter codes; mlb_parlay_lines stores Odds-API canonical names.
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


def _load_samples_for_game(game_id: str) -> pd.DataFrame | None:
    if not ANSWER_KEY_DB.exists():
        return None
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        df = con.execute(
            "SELECT home_margin, total_final_score, home_margin_f5, total_f5, "
            "home_scored_first FROM mlb_game_samples WHERE game_id=?", [game_id]
        ).fetchdf()
    except duckdb.CatalogException:
        return None
    finally:
        con.close()
    return df if not df.empty else None


def _load_book_fairs(game_id: str, spread_line: float, total_line: float) -> dict[str, float]:
    """Pull mlb_sgp_odds rows and devig per book at the given (spread, total)."""
    if not ANSWER_KEY_DB.exists():
        return {}
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        rows = con.execute(
            "SELECT bookmaker, combo, sgp_decimal "
            "FROM mlb_sgp_odds WHERE game_id=? AND period='FG' "
            "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
            [game_id, MAX_BOOK_STALENESS_SEC]
        ).fetchdf()
    except duckdb.CatalogException:
        return {}
    finally:
        con.close()

    out: dict[str, float] = {}
    for book in rows["bookmaker"].unique() if not rows.empty else []:
        sub = rows[rows["bookmaker"] == book].copy()
        fair_per_book = fair_value.devig_book(
            sub, combo="Home Spread + Over",
            vig_fallback=_vig_fallback(book),
        )
        if fair_per_book is not None:
            out[book] = fair_per_book
    return out


def _vig_fallback(book: str) -> float:
    return {
        "draftkings": config.DK_VIG_FALLBACK,
        "fanduel": config.FD_VIG_FALLBACK,
        "prophetx": config.PX_VIG_FALLBACK,
        "novig": config.NOVIG_VIG_FALLBACK,
    }.get(book, 0.10)


def _home_code_from_event_ticker(event_ticker: str) -> str | None:
    """Parse the 3-letter home-team code from a Kalshi event ticker.

    Event suffix format: YYMMMDDHHMM{AwayCode}{HomeCode}, with the AwayCode/HomeCode
    being 2-3 letters each. Convention: away first, home last. We take the trailing
    3 letters as the home code (matches every code in _MLB_CODE_TO_TEAM).
    """
    if "-" not in event_ticker:
        return None
    suffix = event_ticker.rsplit("-", 1)[-1]
    if len(suffix) < 3:
        return None
    return suffix[-3:]


def _current_book_lines_for_combo(game_id: str) -> dict | None:
    """Latest spread/total snapshot from mlb_sgp_odds for line-move detection."""
    if not ANSWER_KEY_DB.exists():
        return None
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT spread_line, total_line FROM mlb_sgp_odds "
            "WHERE game_id=? AND period='FG' "
            "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND "
            "ORDER BY fetch_time DESC LIMIT 1",
            [game_id, MAX_BOOK_STALENESS_SEC]
        ).fetchone()
    except duckdb.CatalogException:
        return None
    finally:
        con.close()
    if not row:
        return None
    return {"spread": float(row[0]), "total": float(row[1])}


# ------------------------------------------------------------------------ #
# Leg-typing helpers                                                       #
# ------------------------------------------------------------------------ #

def _leg_dict_to_typed(leg: dict, game_id: str):
    """Convert {market_ticker, event_ticker, side} to fair_value typed leg.

    Determines team_is_home by parsing the home code from the event_ticker
    (no DB lookup needed — the ticker self-encodes the home/away convention).
    """
    mt = leg["market_ticker"]
    et = leg.get("event_ticker", "")
    side = leg["side"]
    if mt.startswith("KXMLBSPREAD-"):
        suffix = mt.rsplit("-", 1)[-1]
        n_chars = "".join(c for c in suffix if c.isdigit())
        team_chars = "".join(c for c in suffix if not c.isdigit())
        if not n_chars or not team_chars:
            return None
        n = int(n_chars)
        home_code = _home_code_from_event_ticker(et)
        team_is_home = (home_code is not None and team_chars == home_code)
        return fair_value.SpreadLeg(team_is_home=team_is_home, line_n=n, side=side)
    if mt.startswith("KXMLBTOTAL-"):
        try:
            n = int(mt.rsplit("-", 1)[-1])
        except ValueError:
            return None
        return fair_value.TotalLeg(line_n=n, side=side)
    return None


def _spread_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBSPREAD-"):
            suffix = l["market_ticker"].rsplit("-", 1)[-1]
            digits = "".join(c for c in suffix if c.isdigit())
            if digits:
                n = int(digits)
                return -(n - 0.5)
    return 0.0


def _total_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBTOTAL-"):
            try:
                n = int(l["market_ticker"].rsplit("-", 1)[-1])
                return n - 0.5
            except ValueError:
                continue
    return 0.0


# ------------------------------------------------------------------------ #
# Fair-value provider + gate aggregator + Kelly sizing                     #
# ------------------------------------------------------------------------ #

def _fresh_blended_fair(combo_market_ticker: str) -> float | None:
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json, game_id FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    if not row:
        return None
    legs_json, game_id = row
    legs = json.loads(legs_json)

    samples = _load_samples_for_game(game_id)
    if samples is None:
        return None

    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed):
        return None

    model = fair_value.model_fair(samples, typed)
    spread_line = _spread_line_from_legs(legs)
    total_line = _total_line_from_legs(legs)
    book_fairs = _load_book_fairs(game_id, spread_line, total_line)
    return fair_value.blend(model, book_fairs)


def _all_per_accept_gates_pass(quote: dict, fair: float,
                                combo_meta: dict) -> tuple[bool, str]:
    """Run every per-accept gate. Returns (pass, decision_label)."""
    # Prediction staleness
    gen_at = _samples_generated_at()
    if gen_at is None or not risk.staleness_ok(gen_at, config.MAX_PREDICTION_STALENESS_SEC):
        return False, "declined_stale_predictions"

    # Tipoff window
    ct = _commence_time_for_game(combo_meta["game_id"])
    if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
        return False, "declined_tipoff"

    # Line-move check (uses reference_lines snapshot from RFQ submission)
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        ref_row = con.execute(
            "SELECT lines_json FROM reference_lines WHERE rfq_id=?", [rfq_id]
        ).fetchone()
    if ref_row:
        ref_lines = json.loads(ref_row[0])
        current_lines = _current_book_lines_for_combo(combo_meta["game_id"])
        if current_lines and not risk.line_move_ok(
                ref_lines, current_lines, config.LINE_MOVE_THRESHOLD):
            return False, "declined_line_move"

    # Positions API health
    if _POSITIONS_API_FAIL_COUNT >= config.POSITIONS_HEALTH_RETRIES:
        return False, "declined_positions_unhealthy"

    # Fair bounds
    if not risk.fair_in_bounds(fair, config.MIN_FAIR_PROB, config.MAX_FAIR_PROB):
        return False, "declined_kelly_zero"

    # Sanity bound
    no_bid = float(quote.get("no_bid_dollars") or 0)
    quote_implied = 1 - no_bid
    if not risk.sanity_bound_ok(quote_implied, fair, config.MAX_QUOTE_DEVIATION):
        return False, "declined_sanity"

    # Kill switch
    if not risk.kill_switch_ok():
        return False, "declined_killswitch"

    # Cooldown
    leg_set_hash = combo_meta["leg_set_hash"]
    with db.connect(read_only=True) as con:
        rows = con.execute("SELECT leg_set_hash, cooled_until FROM combo_cooldown").fetchall()
    cd_map = {h: u for h, u in rows}
    if not risk.cooldown_ok(leg_set_hash, cd_map):
        return False, "declined_cooldown"

    # Inverse-combo
    legs = json.loads(combo_meta.get("legs_json") or "[]")
    if legs:
        with db.connect(read_only=True) as con:
            held = {h for (h,) in con.execute(
                "SELECT cc.leg_set_hash FROM combo_cache cc "
                "JOIN positions p ON p.combo_market_ticker = cc.combo_market_ticker "
                "WHERE p.net_contracts > 0").fetchall()}
        if not risk.inverse_combo_ok(legs, held):
            return False, "declined_inverse_lock"

    # Per-game cap + daily cap
    today_start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        today_fills = [
            {"game_id": r[0], "contracts": r[1], "price_dollars": r[2]}
            for r in con.execute(
                "SELECT game_id, contracts, price_dollars FROM fills WHERE filled_at >= ?",
                [today_start]
            ).fetchall()
        ]
    if not risk.per_game_cap_ok(combo_meta["game_id"], today_fills,
                                 config.BANKROLL, config.MAX_GAME_EXPOSURE_PCT):
        return False, "declined_per_game_cap"
    if not risk.daily_cap_ok(today_fills, config.DAILY_EXPOSURE_CAP_USD):
        return False, "declined_daily_cap"

    # Fill-ratio halt
    with db.connect(read_only=True) as con:
        window = con.execute(
            "SELECT decision FROM quote_log WHERE decision IN "
            "('accepted', 'failed_quote_walked') ORDER BY observed_at DESC LIMIT ?",
            [config.FILL_RATIO_WINDOW]
        ).fetchall()
    if not risk.fill_ratio_ok([{"decision": d} for (d,) in window],
                                config.MIN_FILL_RATIO):
        return False, "halted_low_fill_ratio"

    return True, "passed"


def _kelly_size_for_quote(quote: dict, fair: float) -> int:
    """Conditional Kelly sizing using mlb_game_samples + existing positions."""
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        meta = con.execute(
            "SELECT game_id, combo_market_ticker FROM live_rfqs WHERE rfq_id=?",
            [rfq_id]
        ).fetchone()
    if not meta:
        return 0
    game_id, combo_market_ticker = meta

    samples = _load_samples_for_game(game_id)
    if samples is None or samples.empty:
        return 0

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    if not row:
        return 0
    legs = json.loads(row[0])
    typed_legs = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed_legs):
        return 0
    mask = pd.Series([True] * len(samples), index=samples.index)
    for leg in typed_legs:
        mask &= fair_value._hit_mask(samples, leg)
    outcome_vec = mask.astype(int).values

    # Existing positions on the same game.
    with db.connect(read_only=True) as con:
        pos_rows = con.execute(
            "SELECT cc.legs_json, p.net_contracts, p.weighted_price "
            "FROM positions p JOIN combo_cache cc "
            "ON cc.combo_market_ticker = p.combo_market_ticker "
            "WHERE p.game_id = ? AND p.net_contracts > 0",
            [game_id]
        ).fetchall()
    existing = []
    for legs_json, n, price in pos_rows:
        leg_objs = [_leg_dict_to_typed(l, game_id) for l in json.loads(legs_json)]
        if any(l is None for l in leg_objs):
            continue
        sub_mask = pd.Series([True] * len(samples), index=samples.index)
        for leg in leg_objs:
            sub_mask &= fair_value._hit_mask(samples, leg)
        existing.append({
            "outcome_vec": sub_mask.astype(int).values,
            "contracts": float(n),
            "effective_price": float(price),
        })

    no_bid = float(quote.get("no_bid_dollars") or 0)
    yes_ask = 1 - no_bid
    fee = ev_calc.fee_per_contract(yes_ask)
    effective_price = yes_ask + fee

    return kelly.kelly_size_combo(
        outcome_vec=np.asarray(outcome_vec),
        existing_positions=existing,
        effective_price=effective_price,
        bankroll=config.BANKROLL,
        kelly_fraction=config.KELLY_FRACTION,
    )


# ------------------------------------------------------------------------ #
# Quote logging + accept                                                   #
# ------------------------------------------------------------------------ #

def _log_quote_decision(quote: dict, fair: float | None,
                         decision: str, reason: str | None = None,
                         post_fee_ev: float | None = None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_log (quote_id, rfq_id, combo_market_ticker, "
            "creator_id, yes_bid_dollars, no_bid_dollars, blended_fair_at_eval, "
            "post_fee_ev_pct, decision, reason_detail, observed_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
            [quote["id"], quote["rfq_id"], quote.get("market_ticker"),
             quote.get("creator_id"),
             float(quote["yes_bid_dollars"]) if quote.get("yes_bid_dollars") else None,
             float(quote["no_bid_dollars"]) if quote.get("no_bid_dollars") else None,
             fair, post_fee_ev, decision, reason,
             datetime.now(timezone.utc)],
        )


def _evaluate_quote(quote: dict, dry_run: bool):
    """Per-quote: evaluate gates, accept if all pass, log decision."""
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT combo_market_ticker, leg_set_hash, game_id "
            "FROM live_rfqs WHERE rfq_id=?", [rfq_id]
        ).fetchone()
    if not row:
        return
    combo_market_ticker, leg_set_hash, game_id = row
    quote = {**quote, "market_ticker": combo_market_ticker}

    with db.connect(read_only=True) as con:
        cc_row = con.execute(
            "SELECT legs_json FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    legs_json = cc_row[0] if cc_row else "[]"

    with ACCEPT_LOCK:
        fair = _fresh_blended_fair(combo_market_ticker)
        if fair is None:
            _log_quote_decision(quote, None, "declined_ev", reason="no_fresh_fair")
            return

        passed, decision = _all_per_accept_gates_pass(
            quote, fair, {"leg_set_hash": leg_set_hash, "game_id": game_id,
                          "legs_json": legs_json})
        if not passed:
            _log_quote_decision(quote, fair, decision)
            return

        no_bid = float(quote.get("no_bid_dollars") or 0)
        ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)
        if ev_pct < config.MIN_EV_PCT:
            _log_quote_decision(quote, fair, "declined_ev", post_fee_ev=ev_pct)
            return

        if dry_run:
            _log_quote_decision(quote, fair, "declined_dry_run", post_fee_ev=ev_pct)
            return

        contracts = _kelly_size_for_quote(quote, fair)
        if contracts <= 0:
            _log_quote_decision(quote, fair, "declined_kelly_zero",
                                 post_fee_ev=ev_pct)
            return

        resp = rfq_client.accept_quote(quote["id"], contracts=contracts)
        if resp is None:
            _log_quote_decision(quote, fair, "failed_quote_walked",
                                 post_fee_ev=ev_pct)
            return

        # Post-accept fill reconciliation via /portfolio/positions
        try:
            actual = rfq_client.get_position_contracts(combo_market_ticker)
            _record_positions_api_result(True)
        except Exception:
            _record_positions_api_result(False)
            actual = contracts

        yes_ask = 1.0 - no_bid
        fee = ev_calc.fee_per_contract(yes_ask)
        with db.connect() as con:
            con.execute(
                "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                "game_id, side, contracts, price_dollars, fee_dollars, "
                "blended_fair_at_fill, expected_ev_dollars, filled_at, raw_response) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [str(uuid.uuid4()), quote["id"], rfq_id, combo_market_ticker, game_id,
                 "yes", actual, yes_ask, fee, fair, ev_dollars,
                 datetime.now(timezone.utc), str(resp)],
            )
            con.execute(
                "INSERT INTO positions (combo_market_ticker, game_id, "
                "net_contracts, weighted_price, legs_json, updated_at) VALUES "
                "(?, ?, ?, ?, ?, ?) ON CONFLICT (combo_market_ticker) DO UPDATE SET "
                "net_contracts = positions.net_contracts + EXCLUDED.net_contracts, "
                "weighted_price = (positions.weighted_price * positions.net_contracts + "
                "                  EXCLUDED.weighted_price * EXCLUDED.net_contracts) / "
                "                 (positions.net_contracts + EXCLUDED.net_contracts), "
                "updated_at = EXCLUDED.updated_at",
                [combo_market_ticker, game_id, actual, yes_ask, legs_json,
                 datetime.now(timezone.utc)],
            )
            cooled_until = datetime.now(timezone.utc) + timedelta(seconds=config.COMBO_COOLDOWN_SEC)
            con.execute(
                "INSERT INTO combo_cooldown (leg_set_hash, game_id, cooled_until, reason) "
                "VALUES (?, ?, ?, ?) ON CONFLICT (leg_set_hash) DO UPDATE SET "
                "cooled_until = EXCLUDED.cooled_until",
                [leg_set_hash, game_id, cooled_until, "post_accept"],
            )
        _log_quote_decision(quote, fair, "accepted", post_fee_ev=ev_pct)
        notify.fill(rfq_id=rfq_id, combo_market_ticker=combo_market_ticker,
                    contracts=actual, price=yes_ask, ev_pct=ev_pct)


# ------------------------------------------------------------------------ #
# RFQ refresh — continuous priority-queue pipeline                          #
# ------------------------------------------------------------------------ #

def mint_and_create_rfq(candidate: combo_enumerator.ComboCandidate,
                         target_cost_dollars: float = 1.0) -> tuple[str, str]:
    """Mint combo ticker (or cache hit) + create RFQ. Returns (rfq_id, combo_ticker)."""
    legs = list(candidate.legs)

    with db.connect(read_only=True) as con:
        cached = con.execute(
            "SELECT combo_market_ticker FROM combo_cache WHERE leg_set_hash=?",
            [candidate.leg_set_hash],
        ).fetchone()

    if cached:
        combo_ticker = cached[0]
    else:
        combo_ticker, combo_event = rfq_client.mint_combo_ticker(
            config.MVE_COLLECTION_TICKER, legs)
        with db.connect() as con:
            con.execute(
                "INSERT INTO combo_cache (leg_set_hash, collection_ticker, "
                "combo_market_ticker, combo_event_ticker, legs_json, game_id) "
                "VALUES (?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                [candidate.leg_set_hash, config.MVE_COLLECTION_TICKER,
                 combo_ticker, combo_event, json.dumps(legs), candidate.game_id],
            )

    rfq_id = rfq_client.create_rfq(combo_ticker, target_cost_dollars=target_cost_dollars)
    return rfq_id, combo_ticker


def _refresh_rfqs(candidates: list[combo_enumerator.ComboCandidate],
                   fair_scores: dict[str, tuple[float, float]],
                   dry_run: bool):
    """Continuous priority-queue pipeline.

    fair_scores: leg_set_hash → (blended_fair, kalshi_ref).
    """
    scored = [(c, *fair_scores[c.leg_set_hash])
              for c in candidates if c.leg_set_hash in fair_scores]
    ranked = combo_enumerator.rank_by_edge(scored)
    target = ranked[: _max_live_rfqs()]
    target_hashes = {c.leg_set_hash for c in target}

    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT rfq_id, leg_set_hash FROM live_rfqs WHERE status='open'"
        ).fetchall()
    live_hashes = {h: rid for rid, h in live}

    # Drop: in DB but not in target.
    for h, rid in list(live_hashes.items()):
        if h not in target_hashes:
            try:
                rfq_client.delete_rfq(rid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_rfqs SET status='cancelled', closed_at=?, "
                        "cancellation_reason='out_of_top_n' WHERE rfq_id=?",
                        [datetime.now(timezone.utc), rid],
                    )
            except Exception as e:
                print(f"  drop {rid} failed: {e}", flush=True)

    # Add: in target but not in DB.
    for c in target:
        if c.leg_set_hash in live_hashes:
            continue
        try:
            rid, combo_ticker = mint_and_create_rfq(c)
            blended_fair, kalshi_ref = fair_scores[c.leg_set_hash]
            edge = combo_enumerator.edge_score(blended_fair, kalshi_ref)
            with db.connect() as con:
                con.execute(
                    "INSERT INTO live_rfqs (rfq_id, combo_market_ticker, leg_set_hash, "
                    "game_id, blended_fair_at_submit, kalshi_ref_at_submit, "
                    "edge_at_submit, status, submitted_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [rid, combo_ticker, c.leg_set_hash, c.game_id,
                     blended_fair, kalshi_ref, edge, "open",
                     datetime.now(timezone.utc)],
                )
                # Snapshot reference lines for line-move detection.
                lines_now = _current_book_lines_for_combo(c.game_id)
                if lines_now:
                    con.execute(
                        "INSERT INTO reference_lines (rfq_id, lines_json, snapped_at) "
                        "VALUES (?, ?, ?) ON CONFLICT (rfq_id) DO NOTHING",
                        [rid, json.dumps(lines_now), datetime.now(timezone.utc)],
                    )
        except Exception as e:
            print(f"  add {c.leg_set_hash[:8]} failed: {e}", flush=True)


# ------------------------------------------------------------------------ #
# Quote poll loop                                                          #
# ------------------------------------------------------------------------ #

def _poll_all_live_rfqs(dry_run: bool):
    with db.connect(read_only=True) as con:
        live = [r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()]
    for rid in live:
        try:
            quotes = rfq_client.poll_quotes(rid, user_id=config.KALSHI_USER_ID or "")
        except Exception:
            continue
        for q in quotes:
            if q.get("status") != "open":
                continue
            with db.connect(read_only=True) as con:
                already = con.execute(
                    "SELECT 1 FROM quote_log WHERE quote_id=?", [q["id"]]
                ).fetchone()
            if already:
                continue
            _evaluate_quote(q, dry_run=dry_run)


# ------------------------------------------------------------------------ #
# Risk sweep (cleanup, not safety — accept gates are the primary defense)  #
# ------------------------------------------------------------------------ #

def _risk_sweep():
    """Tipoff cancel for live RFQs whose game starts within TIPOFF_CANCEL_MIN."""
    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT rfq_id, game_id FROM live_rfqs WHERE status='open'"
        ).fetchall()

    now = datetime.now(timezone.utc)
    for rid, game_id in live:
        ct = _commence_time_for_game(game_id)
        if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN, now=now):
            try:
                rfq_client.delete_rfq(rid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_rfqs SET status='cancelled', closed_at=?, "
                        "cancellation_reason='tipoff' WHERE rfq_id=?",
                        [now, rid],
                    )
            except Exception:
                pass


# ------------------------------------------------------------------------ #
# Per-cycle enumerate + score across all open MLB games                    #
# ------------------------------------------------------------------------ #

def _kalshi_available_spreads(suffix: str, home_code: str, away_code: str
                                 ) -> list[tuple[float, str]]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBSPREAD-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        prefix = f"KXMLBSPREAD-{suffix}-"
        if not ticker.startswith(prefix):
            continue
        spread_part = ticker[len(prefix):]
        digits = "".join(c for c in spread_part if c.isdigit())
        team_chars = "".join(c for c in spread_part if not c.isdigit())
        if not digits or not team_chars:
            continue
        n = int(digits)
        line = -(n - 0.5)
        who = "home" if team_chars == home_code else "away"
        out.append((line, who))
    return out


def _kalshi_available_totals(suffix: str) -> list[float]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBTOTAL-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        try:
            n = int(ticker.rsplit("-", 1)[-1])
            out.append(n - 0.5)
        except ValueError:
            continue
    return out


def _kalshi_last_price(market_ticker: str) -> float:
    status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}")
    if status != 200 or not isinstance(body, dict):
        return 0.0
    m = body.get("market") if isinstance(body, dict) else None
    if not m:
        return 0.0
    try:
        return float(m.get("last_price_dollars", "0") or 0.0)
    except ValueError:
        return 0.0


def _resolve_game_id(home_code: str, away_code: str) -> str | None:
    """Map Kalshi 3-letter codes to mlb.duckdb game_id via mlb_parlay_lines."""
    home = _MLB_CODE_TO_TEAM.get(home_code)
    away = _MLB_CODE_TO_TEAM.get(away_code)
    if not home or not away:
        return None
    if not ANSWER_KEY_DB.exists():
        return None
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT game_id FROM mlb_parlay_lines "
            "WHERE home_team=? AND away_team=? "
            "AND commence_time > NOW() AND commence_time < NOW() + INTERVAL '24 HOUR' "
            "ORDER BY commence_time LIMIT 1",
            [home, away]
        ).fetchone()
    except duckdb.CatalogException:
        return None
    finally:
        con.close()
    return row[0] if row else None


def _enumerate_and_score_all_games() -> tuple[list[combo_enumerator.ComboCandidate],
                                                dict[str, tuple[float, float]]]:
    status, body, _ = auth_client.api(
        "GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    events = body.get("events", []) if status == 200 and isinstance(body, dict) else []

    candidates_all: list[combo_enumerator.ComboCandidate] = []
    fair_scores: dict[str, tuple[float, float]] = {}

    for ev in events:
        event_ticker = ev.get("event_ticker", "")
        if not event_ticker.startswith("KXMLBGAME-"):
            continue
        suffix = event_ticker.replace("KXMLBGAME-", "")
        if len(suffix) < 6:
            continue
        away_code = suffix[-6:-3]
        home_code = suffix[-3:]

        avail_spreads = _kalshi_available_spreads(suffix, home_code, away_code)
        avail_totals = _kalshi_available_totals(suffix)
        if not avail_spreads or not avail_totals:
            continue

        game_id = _resolve_game_id(home_code, away_code)
        if game_id is None:
            continue

        samples = _load_samples_for_game(game_id)
        if samples is None:
            continue

        for cand in combo_enumerator.enumerate_2leg(
                game_id=game_id, event_suffix=suffix,
                home_code=home_code, away_code=away_code,
                available_spreads=avail_spreads,
                available_totals=avail_totals):
            typed = [_leg_dict_to_typed(dict(l), game_id) for l in cand.legs]
            if any(l is None for l in typed):
                continue
            model = fair_value.model_fair(samples, typed)
            spread_line = _spread_line_from_legs([dict(l) for l in cand.legs])
            total_line = _total_line_from_legs([dict(l) for l in cand.legs])
            books = _load_book_fairs(game_id, spread_line, total_line)
            blended = fair_value.blend(model, books)
            if blended is None:
                continue
            if not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
                continue
            kalshi_ref = _kalshi_last_price(cand.legs[0]["market_ticker"])
            candidates_all.append(cand)
            fair_scores[cand.leg_set_hash] = (blended, kalshi_ref)

    return candidates_all, fair_scores


# ------------------------------------------------------------------------ #
# Main loop                                                                #
# ------------------------------------------------------------------------ #

def main_loop(dry_run: bool):
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run, version=VERSION)
    print(f"=== Kalshi MLB RFQ Bot — session {sid} (dry_run={dry_run}) ===", flush=True)
    _phantom_rfq_cleanup()

    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0

    try:
        while _running.is_set():
            now = time.time()

            if KILL_FILE.exists():
                notify.halt("kill_switch")
                time.sleep(config.RISK_SWEEP_SEC)
                continue

            if now - last_rfq_refresh >= config.RFQ_REFRESH_SEC:
                try:
                    candidates, fair_scores = _enumerate_and_score_all_games()
                    _refresh_rfqs(candidates, fair_scores, dry_run=dry_run)
                except Exception as e:
                    print(f"  rfq_refresh error: {e}", flush=True)
                last_rfq_refresh = now

            if now - last_quote_poll >= config.QUOTE_POLL_SEC:
                try:
                    _poll_all_live_rfqs(dry_run=dry_run)
                except Exception as e:
                    print(f"  quote_poll error: {e}", flush=True)
                last_quote_poll = now

            if now - last_risk_sweep >= config.RISK_SWEEP_SEC:
                try:
                    _risk_sweep()
                except Exception as e:
                    print(f"  risk_sweep error: {e}", flush=True)
                last_risk_sweep = now

            if now - last_pipeline >= config.PIPELINE_REFRESH_SEC:
                _run_pipeline()
                last_pipeline = now

            if now - last_heartbeat >= 60:
                print(f"  [HB] {datetime.now(timezone.utc).isoformat()} alive", flush=True)
                last_heartbeat = now

            time.sleep(0.5)
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT rfq_id FROM live_rfqs WHERE status='open'"
            ).fetchall()]
        for rid in live:
            try:
                rfq_client.delete_rfq(rid)
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Run full loop without calling accept_quote.")
    args = parser.parse_args()

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
