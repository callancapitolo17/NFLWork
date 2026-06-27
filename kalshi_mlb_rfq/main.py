"""Kalshi MLB RFQ Bot — autonomous taker daemon."""

import argparse
import json
import logging
import os
import random as _random
import signal
import statistics
import subprocess
import threading
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from kalshi_mlb_rfq import (
    auth_client, combo_enumerator, config, correlation, db, ev_calc,
    fair_value, kelly, notify, research, rfq_client, risk, sgp_runner,
)
from kalshi_mlb_rfq.config import (
    ANSWER_KEY_DB, KILL_FILE, MAX_BOOK_STALENESS_SEC,
)
from kalshi_mlb_rfq.log_setup import setup_logging
from kalshi_common.leg_types import (
    _MLB_CODE_TO_TEAM, _parse_event_suffix, _home_code_from_event_ticker,
    _leg_dict_to_typed, _spread_line_from_legs, _total_line_from_legs,
)

log = logging.getLogger("kalshi_mlb_rfq")

VERSION = "0.1.0"

_running = threading.Event()
_running.set()
ACCEPT_LOCK = threading.Lock()
_CACHE_LOCK = threading.Lock()

# Test-only override; otherwise use config.MAX_LIVE_RFQS.
MAX_LIVE_RFQS_OVERRIDE: int | None = None

# Track positions-API health across calls.
_POSITIONS_API_FAIL_COUNT = 0

# ------------------------------------------------------------------------ #
# In-memory caches of answer-key data.                                     #
# Refreshed on startup and on each PIPELINE_REFRESH_SEC tick. Decouples    #
# the bot's hot path from DuckDB lock contention with the R pipeline +     #
# mlb_dashboard_server. ~10 MB resident — negligible.                      #
# ------------------------------------------------------------------------ #
_SAMPLES_CACHE: dict[str, pd.DataFrame] = {}              # game_id → samples df
_SGP_ODDS_CACHE: pd.DataFrame | None = None                # full mlb_sgp_odds (last hour, FG period)
_PARLAY_LINES_CACHE: dict[str, dict] = {}                  # game_id → {home, away, commence_time}
_SAMPLES_META_GENERATED_AT: datetime | None = None
_CACHE_LOADED_AT: datetime | None = None


def _signal_handler(_sig, _frame):
    _running.clear()


def _max_live_rfqs() -> int:
    return MAX_LIVE_RFQS_OVERRIDE if MAX_LIVE_RFQS_OVERRIDE is not None else config.MAX_LIVE_RFQS


def _research_sample() -> bool:
    """Return True if the next candidate event should be logged.
    Driven by RESEARCH_CANDIDATE_SAMPLING (default 1.0 = log all). Always
    safe to call; short-circuits with zero RNG cost when the knob is 1.0.
    """
    s = config.RESEARCH_CANDIDATE_SAMPLING
    return s >= 1.0 or _random.random() < s


def _emit_candidate_event(outcome: str, *, game_id: str, cand,
                          spread_line: float, total_line: float,
                          model: float | None = None,
                          books: dict | None = None,
                          blended: float | None = None,
                          kalshi_ref: float | None = None,
                          kelly: tuple | None = None) -> None:
    """Emit one candidate_evaluated event (gated by _research_sample).

    Defined at module level (NOT as a per-iter closure in the enumeration
    loop) so we don't allocate a new closure per candidate per tick — and
    so a future deferred-emit refactor can't trip the closure-capture
    footgun. The combo's own market_ticker is intentionally omitted
    (it doesn't exist until RFQ minting); leg_set_hash is the join key
    for analysis against live_rfqs / quote_log.
    """
    if not _research_sample():
        return
    research.emit(
        "candidate_evaluated",
        game_id=game_id,
        leg_set_hash=cand.leg_set_hash,
        spread_line=spread_line, total_line=total_line,
        outcome=outcome, model_fair=model, book_fairs=books,
        n_books=(len(books) if books else 0),
        blended_fair=blended, kalshi_ref=kalshi_ref,
        kelly_yes_n=(kelly[0] if kelly else None),
        kelly_no_n=(kelly[1] if kelly else None),
        worst_yes_ask=(kelly[2] if kelly else None),
        worst_no_ask=(kelly[3] if kelly else None),
    )


def _record_positions_api_result(success: bool):
    global _POSITIONS_API_FAIL_COUNT
    if success:
        _POSITIONS_API_FAIL_COUNT = 0
    else:
        _POSITIONS_API_FAIL_COUNT += 1


# ------------------------------------------------------------------------ #
# Cache refresh — populate _SAMPLES_CACHE / _SGP_ODDS_CACHE etc. in one    #
# DuckDB session with retry-with-backoff for the brief lock window when   #
# the R pipeline is writing.                                              #
# ------------------------------------------------------------------------ #

def _refresh_caches(retries: int = 5) -> bool:
    """Reload all answer-key data into memory. Returns True if successful, False otherwise.

    Tries up to `retries` times with exponential backoff to handle the brief lock
    window when the R pipeline (or dashboard_server) is writing.
    """
    # _SGP_ODDS_CACHE is intentionally NOT in this global list — see
    # _refresh_sgp_cache, which owns it. _refresh_caches must not touch it.
    global _SAMPLES_CACHE, _PARLAY_LINES_CACHE
    global _SAMPLES_META_GENERATED_AT, _CACHE_LOADED_AT

    if not ANSWER_KEY_DB.exists():
        log.warning("cache_refresh: %s does not exist", ANSWER_KEY_DB)
        return False

    last_err = None
    for attempt in range(retries):
        try:
            con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
        except duckdb.IOException as e:
            last_err = e
            wait = 1.0 * (2 ** attempt)
            log.warning("cache_refresh: lock conflict (attempt %d/%d); retrying in %.1fs",
                        attempt + 1, retries, wait)
            time.sleep(wait)
            continue

        try:
            t0 = time.time()

            # mlb_game_samples and mlb_samples_meta are only loaded when the model
            # is active. In book-only mode (USE_MODEL=False) these tables may not
            # exist or be up to date; skipping them avoids DuckDB catalog errors
            # and removes the dependency on the R pipeline entirely.
            if config.USE_MODEL:
                samples_df = con.execute(
                    "SELECT game_id, sim_idx, home_margin, total_final_score, "
                    "home_margin_f5, total_f5 FROM mlb_game_samples"
                ).fetchdf()
                samples_by_game = {gid: g.reset_index(drop=True)
                                    for gid, g in samples_df.groupby("game_id")}
                meta_row = con.execute(
                    "SELECT generated_at FROM mlb_samples_meta "
                    "ORDER BY generated_at DESC LIMIT 1"
                ).fetchone()
                generated_at = meta_row[0] if meta_row else None
            else:
                samples_by_game = {}
                generated_at = None

            # mlb_sgp_odds is owned exclusively by _refresh_sgp_cache (reads
            # the bot market DB on its own 60s cadence). Do NOT touch
            # _SGP_ODDS_CACHE in this function — _refresh_caches runs on the
            # 10-min PIPELINE_REFRESH_SEC tick and would wipe the fresh
            # bot-DB rows _refresh_sgp_cache loaded seconds ago.

            # mlb_parlay_lines replaced by bot DB::mlb_target_lines cache.
            # Schedule (team names + commence_time) now lives directly in
            # mlb_target_lines rows, written by sgp_runner from the Odds API.
            parlay_lines = _build_parlay_lines_cache(
                bot_db=str(config.BOT_MARKET_DB),
            )
        except duckdb.CatalogException as e:
            # Schema mismatch / table missing — typically the R pipeline is
            # mid-atomic-replace (DROP + CREATE). Retry with backoff parallel
            # to the IOException handling above instead of waiting the full
            # PIPELINE_REFRESH_SEC tick.
            last_err = e
            try:
                con.close()
            except Exception:
                pass
            wait = 1.0 * (2 ** attempt)
            log.warning("cache_refresh: schema mismatch (attempt %d/%d); retrying in %.1fs — %s",
                        attempt + 1, retries, wait, e)
            time.sleep(wait)
            continue
        finally:
            try:
                con.close()
            except Exception:
                pass

        # Atomic swap into the caches under the lock.
        # _SGP_ODDS_CACHE intentionally NOT touched here — owned by
        # _refresh_sgp_cache.
        with _CACHE_LOCK:
            _SAMPLES_CACHE = samples_by_game
            _PARLAY_LINES_CACHE = parlay_lines
            _SAMPLES_META_GENERATED_AT = generated_at
            _CACHE_LOADED_AT = datetime.now(timezone.utc)

        elapsed = time.time() - t0
        log.info("cache_refresh: %d games, %d parlay_lines, samples gen_at=%s (%.1fs)",
                 len(samples_by_game), len(parlay_lines), generated_at, elapsed)
        return True

    log.warning("cache_refresh: gave up after %d attempts; last error: %s", retries, last_err)
    return False


def _refresh_sgp_cache(retries: int = 3) -> bool:
    """Narrow variant of _refresh_caches — reloads only mlb_sgp_odds from
    the bot market DB. Atomic swap under _CACHE_LOCK.

    Returns False on missing DB / table; True on success.
    """
    global _SGP_ODDS_CACHE
    bot_db = str(config.BOT_MARKET_DB)
    if not Path(bot_db).exists():
        return False
    last_err = None
    for attempt in range(retries):
        try:
            con = duckdb.connect(bot_db, read_only=True)
        except duckdb.IOException as e:
            last_err = e
            time.sleep(0.5 * (2 ** attempt))
            continue
        try:
            tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
            if "mlb_sgp_odds" not in tables:
                return False
            sgp_df = con.execute(
                "SELECT game_id, combo, period, bookmaker, sgp_decimal, fetch_time, "
                "spread_line, total_line "
                "FROM mlb_sgp_odds WHERE period='FG' "
                "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
                [config.MAX_BOOK_STALENESS_SEC],
            ).fetchdf()
        finally:
            con.close()

        with _CACHE_LOCK:
            _SGP_ODDS_CACHE = sgp_df
        log.info("sgp_cache_refresh: %d rows from bot_market_db", len(sgp_df))
        return True

    log.warning("sgp_cache_refresh: gave up after %d attempts; %s", retries, last_err)
    return False


def _build_parlay_lines_cache(bot_db: str) -> dict[str, dict]:
    """Build the bot's parlay_lines cache from bot DB::mlb_target_lines only.

    Returns {game_id: {home_team, away_team, commence_time, fg_lines: list[(spread, total)]}}.
    Drops games with no target lines.

    Schedule (game_id, team names, commence_time) is sourced from
    mlb_target_lines rows directly — written by sgp_runner.write_target_lines
    from Odds API data on each cycle. The old mlb.duckdb::mlb_odds_temp
    read path was removed in Task 29 to decouple the bot from the R-locked
    dashboard DB.
    """
    from pathlib import Path as _Path

    out: dict[str, dict] = {}
    if not _Path(bot_db).exists():
        return out
    con = duckdb.connect(bot_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_target_lines" not in tables:
            return out
        rows = con.execute(
            "SELECT game_id, home_team, away_team, commence_time, spread, total "
            "FROM mlb_target_lines WHERE period = 'FG' "
            "ORDER BY game_id, spread, total"
        ).fetchall()
    finally:
        con.close()
    for game_id, home, away, ct, spread, total in rows:
        if game_id not in out:
            out[game_id] = {
                "home_team": home, "away_team": away,
                "commence_time": ct, "fg_lines": [],
            }
        out[game_id]["fg_lines"].append((spread, total))
    return out


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
    """Cancel any RFQs on Kalshi that are OURS (cross-category combo) but missing
    from live_rfqs. CRITICAL: never touches RFQs that belong to other bots / user
    actions — those are identified by their market_ticker NOT starting with our
    combo prefix (e.g., NFL Draft RFQs from kalshi_draft are off-limits).
    Also requires KALSHI_USER_ID to be set so we only see our own RFQs (without
    it, Kalshi returns the global open-RFQ list, which we have no business
    touching).
    """
    if not config.KALSHI_USER_ID:
        log.warning("startup: KALSHI_USER_ID not set — skipping phantom cleanup "
                    "(safety: would otherwise fetch RFQs from all users)")
        return
    try:
        kalshi_open = rfq_client.list_open_rfqs(config.KALSHI_USER_ID)
    except Exception as e:
        log.warning("startup: list_open_rfqs failed: %s", e)
        return

    if not kalshi_open:
        log.info("startup: no open RFQs on Kalshi")
        return

    with db.connect(read_only=True) as con:
        ours = {r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()}

    # Filter to our combo namespace. Our combos are minted via the cross-category
    # MVE collection so the resulting market_ticker starts with the collection's
    # series ticker (KXMVECROSSCATEGORY-S-... per recon). Anything else — NFL Draft
    # RFQs, Mention RFQs, single-market RFQs from other bots — belongs to someone
    # else. Skip them entirely.
    our_combo_prefix = config.MVE_COLLECTION_TICKER.replace("-R", "-S")  # "KXMVECROSSCATEGORY-S"

    skipped_other_bot = 0
    cancelled = 0
    failed = 0
    for rfq in kalshi_open:
        rid = rfq.get("id")
        market_ticker = rfq.get("market_ticker", "") or ""
        if not rid or rid in ours:
            continue
        if not market_ticker.startswith(our_combo_prefix):
            skipped_other_bot += 1
            continue
        try:
            rfq_client.delete_rfq(rid)
            cancelled += 1
            log.info("startup: cancelled phantom rfq %s", rid)
        except Exception as e:
            failed += 1
            log.warning("startup: failed to cancel phantom %s: %s", rid, e)

    log.info("startup: phantom cleanup — cancelled=%d failed=%d skipped_other_bot=%d",
             cancelled, failed, skipped_other_bot)


# ------------------------------------------------------------------------ #
# Answer-key DB helpers (read-only)                                        #
# ------------------------------------------------------------------------ #

def _samples_generated_at() -> datetime | None:
    """Read from in-memory cache, populated by _refresh_caches()."""
    return _SAMPLES_META_GENERATED_AT


def _commence_time_for_game(game_id: str) -> datetime | None:
    """Read from in-memory cache, populated by _refresh_caches()."""
    row = _PARLAY_LINES_CACHE.get(game_id)
    return row["commence_time"] if row else None



def _load_samples_for_game(game_id: str) -> pd.DataFrame | None:
    """Read from in-memory cache, populated by _refresh_caches()."""
    df = _SAMPLES_CACHE.get(game_id)
    return df if df is not None and not df.empty else None


def _load_book_fairs(game_id: str, spread_line: float, total_line: float,
                     spread_side: str = "home", total_side: str = "over") -> dict[str, float]:
    """Per-line book lookup with N>=2 gate.

    Returns {book -> devigged_fair} only if at least MIN_BOOK_COUNT_FOR_BLEND
    books priced the matching (game, spread, total) tuple. Empty dict
    otherwise — caller treats as 'no book signal, drop candidate'.

    spread_side/total_side select WHICH of the four SGP quadrants to devig
    ("Home Spread + Over" etc.). The cache stores one home-perspective
    spread_line per (game, total) with all four quadrant rows, so the row
    filter is the same for every quadrant — only the devig target differs.
    Defaults preserve the legacy "Home Spread + Over" behaviour for callers
    that don't specify a quadrant.
    """
    if _SGP_ODDS_CACHE is None or _SGP_ODDS_CACHE.empty:
        return {}
    if "spread_line" not in _SGP_ODDS_CACHE.columns or "total_line" not in _SGP_ODDS_CACHE.columns:
        # Transition state — cache still has legacy schema. Drop candidates
        # until Task 26 wires _refresh_sgp_cache (new-schema reader).
        return {}
    rows = _SGP_ODDS_CACHE[
        (_SGP_ODDS_CACHE["game_id"] == game_id)
        & (_SGP_ODDS_CACHE["spread_line"].astype(float).round(2) == round(float(spread_line), 2))
        & (_SGP_ODDS_CACHE["total_line"].astype(float).round(2) == round(float(total_line), 2))
    ]
    if rows.empty:
        return {}
    label = f"{spread_side.title()} Spread + {total_side.title()}"
    out: dict[str, float] = {}
    for book in rows["bookmaker"].unique():
        sub = rows[rows["bookmaker"] == book].copy()
        fair_per_book = fair_value.devig_book(
            sub, combo=label,
            vig_fallback=_vig_fallback(book),
        )
        if fair_per_book is not None:
            out[book] = fair_per_book
    if len(out) < config.MIN_BOOK_COUNT_FOR_BLEND:
        return {}
    return out


def _vig_fallback(book: str) -> float:
    return {
        "draftkings": config.DK_VIG_FALLBACK,
        "fanduel": config.FD_VIG_FALLBACK,
        "prophetx": config.PX_VIG_FALLBACK,
        "novig": config.NOVIG_VIG_FALLBACK,
        "betmgm": config.BETMGM_VIG_FALLBACK,
        "caesars": config.CAESARS_VIG_FALLBACK,
    }.get(book, 0.10)


# ------------------------------------------------------------------------ #
# Fair-value provider + gate aggregator + Kelly sizing                     #
# ------------------------------------------------------------------------ #

def _book_only_fair(book_fairs: dict) -> float | None:
    """Median of book fairs (>=2-book floor enforced upstream). None if empty."""
    vals = [v for v in book_fairs.values() if v is not None]
    return statistics.median(vals) if vals else None


def _staleness_gate_ok() -> tuple[bool, str]:
    """Model-prediction staleness — only meaningful when the model prices.

    When USE_MODEL is False the bot never loads mlb_game_samples / mlb_samples_meta,
    so _SAMPLES_META_GENERATED_AT will always be None. Gating on it would decline
    every quote; skip the check entirely in book-only mode.
    """
    if not config.USE_MODEL:
        return True, "passed"
    gen_at = _samples_generated_at()
    if gen_at is None or not risk.staleness_ok(gen_at, config.MAX_PREDICTION_STALENESS_SEC):
        return False, "declined_stale_predictions"
    return True, "passed"


def _grid_lookup(game_id: str):
    """Return a closure: (spread_line, total_line, spread_side, total_side) -> median
    book probability for that SGP quadrant, or None if <MIN_BOOK_COUNT_FOR_BLEND
    books priced it.

    The closure devigged each book's row the same way _load_book_fairs does,
    but accepts the quadrant label dynamically so all four combos are priceable.
    """
    def lookup(spread_line, total_line, spread_side, total_side):
        if _SGP_ODDS_CACHE is None or _SGP_ODDS_CACHE.empty:
            return None
        label = f"{spread_side.title()} Spread + {total_side.title()}"
        rows = _SGP_ODDS_CACHE[
            (_SGP_ODDS_CACHE["game_id"] == game_id)
            & (_SGP_ODDS_CACHE["spread_line"].astype(float).round(2) == round(float(spread_line), 2))
            & (_SGP_ODDS_CACHE["total_line"].astype(float).round(2) == round(float(total_line), 2))
        ]
        if rows.empty:
            return None
        out = []
        for book in rows["bookmaker"].unique():
            sub = rows[rows["bookmaker"] == book]
            f = fair_value.devig_book(sub, combo=label, vig_fallback=_vig_fallback(book))
            if f is not None:
                out.append(f)
        if len(out) < config.MIN_BOOK_COUNT_FOR_BLEND:
            return None
        return statistics.median(out)
    return lookup


def _combo_region_from_legs(typed_legs: list) -> "correlation.ComboRegion | None":
    """Map a list of typed legs (SpreadLeg / TotalLeg) to a ComboRegion.

    Returns None if the combo lacks exactly one spread and one total leg
    (e.g. player-prop cross-category legs that can't be priced on the grid).

    Spread-side convention: ComboRegion.spread_side is the COVERING team
    ("home" when the home team covers, "away" when the away team covers),
    independent of which Kalshi side (yes/no) the bet is on.

    SpreadLeg mapping:
      team_is_home=True,  side="yes"  → home covers    → "home"
      team_is_home=True,  side="no"   → away covers    → "away"
      team_is_home=False, side="yes"  → away covers    → "away"
      team_is_home=False, side="no"   → home covers    → "home"
    ⟹ spread_side = "home" if (team_is_home == (side == "yes")) else "away"

    spread_line sign depends on team_is_home: home-margin leg → -(line_n-0.5)
    (negative grid); away-margin leg → +(line_n-0.5) (positive grid).
    spread_side selects the cell within the grid (which team covers).

    TotalLeg mapping:
      side="yes" → over;  side="no" → under
      total_line = line_n - 0.5  (same formula as _total_line_from_legs)
    """
    spread_leg = None
    total_leg = None
    for leg in typed_legs:
        if isinstance(leg, fair_value.SpreadLeg):
            if spread_leg is not None:
                return None  # multiple spread legs → not a grid combo
            spread_leg = leg
        elif isinstance(leg, fair_value.TotalLeg):
            if total_leg is not None:
                return None  # multiple total legs → not a grid combo
            total_leg = leg
    if spread_leg is None or total_leg is None:
        return None  # missing one of the two legs

    spread_side = "home" if (spread_leg.team_is_home == (spread_leg.side == "yes")) else "away"
    # Home-perspective signed line: home margin → -(n-0.5); away margin → +(n-0.5).
    # The sign selects the grid (home-favorite vs away-favorite); spread_side
    # selects the cell within it.
    spread_line = (-(spread_leg.line_n - 0.5) if spread_leg.team_is_home
                   else (spread_leg.line_n - 0.5))
    total_side = "over" if total_leg.side == "yes" else "under"
    total_line = total_leg.line_n - 0.5         # e.g. n=8 → 7.5

    return correlation.ComboRegion(
        spread_side=spread_side,
        spread_line=spread_line,
        total_side=total_side,
        total_line=total_line,
    )


def _combo_fair_for_region(game_id: str, r: "correlation.ComboRegion") -> "float | None":
    """Book-implied probability for a ComboRegion on a specific game."""
    return _grid_lookup(game_id)(r.spread_line, r.total_line, r.spread_side, r.total_side)


def _book_implied_cov(game_id: str,
                       new_region: "correlation.ComboRegion",
                       new_price: float,
                       pos_region: "correlation.ComboRegion",
                       pos_price: float) -> float:
    """Return-space covariance between a new bet and an existing position,
    using book-implied joint probability from the SGP grid.

    Falls back to ρ=1 (joint = min(P_a, P_b)) when the grid cannot price the
    joint — conservative because it maximally penalises correlated positions.
    Returns 0.0 if either single-leg probability is missing (can't price one
    side at all → skip the correlation term rather than guess).
    """
    p_a = _combo_fair_for_region(game_id, new_region)
    p_b = _combo_fair_for_region(game_id, pos_region)
    if p_a is None or p_b is None:
        return 0.0  # cannot price one leg book-only → no correlation term
    if not (0 < new_price < 1) or not (0 < pos_price < 1):
        return 0.0  # degenerate price would cause ZeroDivisionError in cov_returns
    j = correlation.joint_prob(new_region, pos_region, p_a, p_b, _grid_lookup(game_id))
    if j is None:
        j = min(p_a, p_b)  # ρ=1 conservative fallback
    return correlation.cov_returns(p_a, p_b, j, new_price, pos_price)


def _fresh_blended_fair(combo_market_ticker: str) -> tuple[float | None, dict]:
    """Return (blended_fair, book_fairs).

    blended_fair is None when no fresh fair can be formed; book_fairs is the
    per-book devigged dict (empty when none). Returning book_fairs lets the
    caller log the pricing components without a second _load_book_fairs pass.
    """
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json, game_id FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    if not row:
        return None, {}
    legs_json, game_id = row
    legs = json.loads(legs_json)

    total_line = _total_line_from_legs(legs)
    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed):
        return None, {}
    region = _combo_region_from_legs(typed)
    if region is None:
        return None, {}
    book_fairs = _load_book_fairs(game_id, region.spread_line, region.total_line,
                                  region.spread_side, region.total_side)

    if not config.USE_MODEL:
        return _book_only_fair(book_fairs), book_fairs

    samples = _load_samples_for_game(game_id)
    if samples is None:
        return None, {}
    model = fair_value.model_fair(samples, typed)
    return fair_value.blend(model, book_fairs), book_fairs


def _all_per_accept_gates_pass(quote: dict, fair: float,
                                combo_meta: dict) -> tuple[bool, str]:
    """Run every per-accept gate. Returns (pass, decision_label)."""
    # Prediction staleness (skipped when USE_MODEL is False)
    ok, decision = _staleness_gate_ok()
    if not ok:
        return False, decision

    # Tipoff window
    ct = _commence_time_for_game(combo_meta["game_id"])
    if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
        return False, "declined_tipoff"

    # Line-move gate removed (D7 of line-source pivot): each candidate's
    # (spread, total) is baked into its Kalshi ticker — line can't move
    # per-candidate. Drift handled by RFQ refresh re-scoring every 30s.

    # Positions API health
    if _POSITIONS_API_FAIL_COUNT >= config.POSITIONS_HEALTH_RETRIES:
        return False, "declined_positions_unhealthy"

    # Fair bounds
    if not risk.fair_in_bounds(fair, config.MIN_FAIR_PROB, config.MAX_FAIR_PROB):
        return False, "declined_kelly_zero"

    # Kill switch
    if not risk.kill_switch_ok():
        return False, "declined_killswitch"

    # Cooldown is now side-aware and is checked in _evaluate_quote AFTER side
    # selection (we need the chosen side to look up the right cooldown row).

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


def _outcome_vec_for_legs(samples: pd.DataFrame, typed_legs: list,
                           side: str = "yes") -> np.ndarray | None:
    """AND the per-leg hit masks into one 0/1 outcome vector. For NO the
    vector inverts (we win when the combo MISSES). Returns None on bad legs."""
    if any(l is None for l in typed_legs):
        return None
    mask = pd.Series([True] * len(samples), index=samples.index)
    for leg in typed_legs:
        mask &= fair_value._hit_mask(samples, leg)
    hit_vec = mask.astype(int).values
    return hit_vec if side == "yes" else (1 - hit_vec)


def _load_existing_positions_for_game(game_id: str,
                                       samples: pd.DataFrame) -> list[dict]:
    """Existing combo positions on this game, with outcome vectors for
    conditional Kelly. Each row's outcome_vec is the combo-hit mask if
    the held side is YES, or its complement if NO — so a long-NO position
    correctly pays off when the combo misses. Filters out positions whose
    legs no longer type-check against current samples."""
    with db.connect(read_only=True) as con:
        pos_rows = con.execute(
            "SELECT cc.legs_json, p.net_contracts, p.weighted_price, p.side "
            "FROM positions p JOIN combo_cache cc "
            "ON cc.combo_market_ticker = p.combo_market_ticker "
            "WHERE p.game_id = ? AND p.net_contracts > 0",
            [game_id]
        ).fetchall()
    existing = []
    for legs_json, n, price, pos_side in pos_rows:
        leg_objs = [_leg_dict_to_typed(l, game_id) for l in json.loads(legs_json)]
        vec = _outcome_vec_for_legs(samples, leg_objs, side=pos_side)
        if vec is None:
            continue
        existing.append({
            "outcome_vec": vec,
            "contracts": float(n),
            "effective_price": float(price),
        })
    return existing


def _load_existing_positions_book(game_id: str,
                                   new_region: "correlation.ComboRegion | None",
                                   new_price: float,
                                   side: str) -> list[dict]:
    """Book-only equivalent of _load_existing_positions_for_game.

    Reads held same-game positions from the positions table (leg-sets via
    combo_cache), types each position's legs, derives its ComboRegion, and
    builds the cov_return via the book-implied correlation engine.

    Positions whose legs can't be mapped to a ComboRegion (non-grid combos,
    e.g. player-prop cross-category legs) are skipped — they contribute 0
    covariance by omission, which is equivalent to assuming independence.

    If new_region is None (the new bet itself is non-grid), we can't compute
    any covariance, so return an empty list (no correlation adjustment).
    """
    if new_region is None:
        return []
    with db.connect(read_only=True) as con:
        pos_rows = con.execute(
            "SELECT cc.legs_json, p.net_contracts, p.weighted_price, p.side "
            "FROM positions p JOIN combo_cache cc "
            "ON cc.combo_market_ticker = p.combo_market_ticker "
            "WHERE p.game_id = ? AND p.net_contracts > 0",
            [game_id]
        ).fetchall()
    existing = []
    for legs_json, n, price, pos_side in pos_rows:
        typed_pos = [_leg_dict_to_typed(l, game_id) for l in json.loads(legs_json)]
        if any(t is None for t in typed_pos):
            continue
        pos_region = _combo_region_from_legs(typed_pos)
        if pos_region is None:
            continue  # non-grid position; skip correlation term
        cov = _book_implied_cov(game_id, new_region, new_price, pos_region, float(price))
        existing.append({
            "cov_return": cov,
            "contracts": float(n),
            "effective_price": float(price),
        })
    return existing


def _kelly_size_for_quote(quote: dict, fair: float, side: str = "yes") -> int:
    """Conditional Kelly sizing at quote-evaluation time. The maker's
    no_bid/yes_bid is known, so effective_price is exact (no estimation).

    `side` is the side we're considering BUYING ("yes" or "no"). For NO,
    the outcome_vec inverts (we win when the combo misses), effective_price
    uses no_ask, and blended_fair flips to (1 - fair). Existing positions
    inherit their stored side via _load_existing_positions_for_game.

    Used as a go/no-go gate before accept — REST accept is all-or-nothing,
    so this can't downsize the fill. A Kelly==0 result means even accepting
    the maker's full quote would be -EV after correlation with held positions.
    """
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        meta = con.execute(
            "SELECT game_id, combo_market_ticker FROM live_rfqs WHERE rfq_id=?",
            [rfq_id]
        ).fetchone()
    if not meta:
        return 0
    game_id, combo_market_ticker = meta

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    if not row:
        return 0
    legs = json.loads(row[0])
    typed_legs = [_leg_dict_to_typed(l, game_id) for l in legs]

    no_bid  = float(quote.get("no_bid_dollars")  or 0)
    yes_bid = float(quote.get("yes_bid_dollars") or 0)
    if side == "yes":
        ask = 1 - no_bid
    else:
        ask = 1 - yes_bid
    fee = ev_calc.fee_per_contract(ask)
    effective_price = ask + fee

    if config.USE_MODEL:
        samples = _load_samples_for_game(game_id)
        if samples is None or samples.empty:
            return 0
        outcome_vec = _outcome_vec_for_legs(samples, typed_legs, side=side)
        if outcome_vec is None:
            return 0
        existing = _load_existing_positions_for_game(game_id, samples)
    else:
        outcome_vec = None
        new_region = _combo_region_from_legs(typed_legs)
        existing = _load_existing_positions_book(game_id, new_region, effective_price, side)

    return kelly.kelly_size_combo(
        outcome_vec=outcome_vec,
        blended_fair=fair if side == "yes" else (1 - fair),
        existing_positions=existing,
        effective_price=effective_price,
        bankroll=config.BANKROLL,
        kelly_fraction=config.KELLY_FRACTION,
    )


def _worst_acceptable_ask(blended_fair: float, side: str,
                            ev_floor_pct: float) -> float:
    """Highest ask on `side` whose post-fee EV % still meets ev_floor_pct.

    Symmetric over both sides via ev_calc:
      - side="yes": ev = post_fee_ev_buy_yes(blended_fair, no_bid),
                    yes_ask = 1 - no_bid; max yes_ask is bounded by fair.
      - side="no":  ev = post_fee_ev_buy_no(blended_fair, yes_bid),
                    no_ask  = 1 - yes_bid; max no_ask is bounded by (1-fair).

    Binary search because the fee is itself a function of price, so closed-
    form inversion is messy. 30 iters → sub-cent precision.

    Returns 0.0 if no acceptable price exists in the side's range (typically
    when the fee floor exceeds the available edge for that side).
    """
    if side == "yes":
        hi = max(0.0001, blended_fair - 1e-6)
        def ev_at(ask: float) -> float:
            _, ev_pct = ev_calc.post_fee_ev_buy_yes(blended_fair, 1.0 - ask)
            return ev_pct
    else:
        hi = max(0.0001, (1.0 - blended_fair) - 1e-6)
        def ev_at(ask: float) -> float:
            _, ev_pct = ev_calc.post_fee_ev_buy_no(blended_fair, 1.0 - ask)
            return ev_pct

    lo = 0.0001
    if hi <= lo:
        return 0.0
    for _ in range(30):
        mid = (lo + hi) / 2
        if ev_at(mid) >= ev_floor_pct:
            lo = mid   # this ask is acceptable; try a higher one
        else:
            hi = mid
    return lo if ev_at(lo) >= ev_floor_pct else 0.0


def _kelly_size_for_candidate(game_id: str,
                                typed_legs: list,
                                samples: "pd.DataFrame | None",
                                blended_fair: float,
                                leg_set_hash: str | None = None,
                                ) -> tuple[int, int, float, float]:
    """Kelly counts and worst-acceptable prices for BOTH sides at create time.

    Returns (kelly_yes_n, kelly_no_n, worst_yes_ask, worst_no_ask). Any side
    with no acceptable price returns (n=0, ask=0.0) for that side.

    When USE_MODEL is True: loads game samples + builds outcome_vec per side.
    When USE_MODEL is False: uses book-implied correlation for existing positions;
      samples is expected to be None and is not consulted.

    Existing positions on the same game are loaded once and reused:
      - model path: via _load_existing_positions_for_game (outcome_vec based)
      - book path:  via _load_existing_positions_book (cov_return based)
    """
    if config.USE_MODEL:
        if samples is None or samples.empty:
            return 0, 0, 0.0, 0.0
        existing = _load_existing_positions_for_game(game_id, samples)
    else:
        # Pre-compute the new combo's region once; reused for both sides.
        new_region = _combo_region_from_legs(typed_legs)

    def kelly_for_side(side: str) -> tuple[int, float]:
        ask = _worst_acceptable_ask(
            blended_fair, side, config.KELLY_CREATE_EV_FLOOR_PCT)
        if ask <= 0:
            return 0, 0.0
        fee = ev_calc.fee_per_contract(ask)
        effective_price = ask + fee
        p = blended_fair if side == "yes" else (1.0 - blended_fair)

        if config.USE_MODEL:
            outcome_vec = _outcome_vec_for_legs(samples, typed_legs, side=side)
            if outcome_vec is None:
                return 0, 0.0
            existing_for_side = existing
        else:
            outcome_vec = None
            existing_for_side = _load_existing_positions_book(
                game_id, new_region, effective_price, side)

        n = kelly.kelly_size_combo(
            outcome_vec=outcome_vec,
            blended_fair=p,
            existing_positions=existing_for_side,
            effective_price=effective_price,
            bankroll=config.BANKROLL,
            kelly_fraction=config.KELLY_FRACTION,
        )
        return n, ask

    yes_n, yes_ask = kelly_for_side("yes")
    no_n,  no_ask  = kelly_for_side("no")
    existing_for_log = existing if config.USE_MODEL else []
    research.emit("kelly_sized", game_id=game_id,
                  leg_set_hash=leg_set_hash,
                  blended_fair=blended_fair,
                  kelly_yes_n=yes_n, kelly_no_n=no_n,
                  worst_yes_ask=yes_ask, worst_no_ask=no_ask,
                  n_existing_positions=len(existing_for_log),
                  bankroll=config.BANKROLL,
                  kelly_fraction=config.KELLY_FRACTION)
    return yes_n, no_n, yes_ask, no_ask


# ------------------------------------------------------------------------ #
# Quote logging + accept                                                   #
# ------------------------------------------------------------------------ #

def _log_quote_decision(quote: dict, fair: float | None,
                         decision: str, reason: str | None = None,
                         post_fee_ev: float | None = None,
                         diag: dict | None = None):
    """Write one quote evaluation outcome to quote_log.

    `diag` carries walk-diagnostic context (competitor state, latency timestamps,
    accept-error body, rfq terminal status) plus side-selection context
    (chosen_side, ev_yes_pct, ev_no_pct) and the hedge-formation tracker
    (hedge_added, hedge_original_side/price, hedge_new_price, hedge_current_fair,
    hedge_projected_net). All keys are optional; missing keys are stored as NULL.
    See db.py SCHEMA_SQL + MIGRATE_SQL for the column set.
    """
    d = diag or {}
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_log (quote_id, rfq_id, combo_market_ticker, "
            "creator_id, yes_bid_dollars, no_bid_dollars, blended_fair_at_eval, "
            "post_fee_ev_pct, decision, reason_detail, observed_at, "
            "competitor_count, best_competitor_no_bid_dollars, "
            "accept_response_body, rfq_terminal_status, "
            "quote_first_seen_at, accept_attempted_at, accept_response_at, "
            "chosen_side, ev_yes_pct, ev_no_pct, "
            "hedge_added, hedge_original_side, hedge_original_price, "
            "hedge_new_price, hedge_current_fair, hedge_projected_net) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
            "?, ?, ?, ?, ?, ?, ?, ?, ?) "
            "ON CONFLICT DO NOTHING",
            [quote["id"], quote["rfq_id"], quote.get("market_ticker"),
             quote.get("creator_id"),
             float(quote["yes_bid_dollars"]) if quote.get("yes_bid_dollars") else None,
             float(quote["no_bid_dollars"]) if quote.get("no_bid_dollars") else None,
             fair, post_fee_ev, decision, reason,
             datetime.now(timezone.utc),
             d.get("competitor_count"),
             d.get("best_competitor_no_bid_dollars"),
             d.get("accept_response_body"),
             d.get("rfq_terminal_status"),
             d.get("quote_first_seen_at"),
             d.get("accept_attempted_at"),
             d.get("accept_response_at"),
             d.get("chosen_side"),
             d.get("ev_yes_pct"),
             d.get("ev_no_pct"),
             d.get("hedge_added"),
             d.get("hedge_original_side"),
             d.get("hedge_original_price"),
             d.get("hedge_new_price"),
             d.get("hedge_current_fair"),
             d.get("hedge_projected_net")],
        )


def _evaluate_quote(quote: dict, dry_run: bool,
                    *,
                    competitor_count: int | None = None,
                    best_competitor_no_bid_dollars: float | None = None,
                    quote_first_seen_at: datetime | None = None):
    """Per-quote: evaluate gates, accept if all pass, log decision.

    The kwargs starting with `competitor_count` carry walk-diagnostic context
    computed by the caller from the full quotes-on-RFQ list. They are stored on
    every decision row (declined and walked alike) so we can later distinguish
    'we were too slow' from 'we were too cheap' on walks. See _log_quote_decision.
    """
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT combo_market_ticker, leg_set_hash, game_id, intended_side "
            "FROM live_rfqs WHERE rfq_id=?", [rfq_id]
        ).fetchone()
    if not row:
        return
    combo_market_ticker, leg_set_hash, game_id, intended_side = row
    quote = {**quote, "market_ticker": combo_market_ticker}

    with db.connect(read_only=True) as con:
        cc_row = con.execute(
            "SELECT legs_json FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    legs_json = cc_row[0] if cc_row else "[]"

    diag: dict = {
        "competitor_count": competitor_count,
        "best_competitor_no_bid_dollars": best_competitor_no_bid_dollars,
        "quote_first_seen_at": quote_first_seen_at,
    }

    with ACCEPT_LOCK:
        fair, book_fairs = _fresh_blended_fair(combo_market_ticker)
        if fair is None:
            _log_quote_decision(quote, None, "declined_ev", reason="no_fresh_fair", diag=diag)
            return
        # quote_priced: reuse the per-book fairs already computed by
        # _fresh_blended_fair — no redundant pandas recompute per quote eval.
        research.emit("quote_priced", game_id=game_id,
                      combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                      quote_id=quote.get("id"), blended_fair=fair,
                      book_fairs=book_fairs)

        passed, decision = _all_per_accept_gates_pass(
            quote, fair, {"leg_set_hash": leg_set_hash, "game_id": game_id,
                          "legs_json": legs_json})
        research.emit("gate_evaluated", game_id=game_id,
                      combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                      quote_id=quote.get("id"),
                      decision=decision, passed=passed,
                      blended_fair=fair)
        if not passed:
            _log_quote_decision(quote, fair, decision, diag=diag)
            return

        no_bid  = float(quote.get("no_bid_dollars")  or 0)
        yes_bid = float(quote.get("yes_bid_dollars") or 0)

        # Symmetric EV: evaluate both sides against current fair after fees.
        ev_yes_dollars, ev_yes_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)
        ev_no_dollars,  ev_no_pct  = ev_calc.post_fee_ev_buy_no (fair, yes_bid)
        yes_ok = ev_yes_pct >= config.MIN_EV_PCT and no_bid  > 0
        no_ok  = ev_no_pct  >= config.MIN_EV_PCT and yes_bid > 0

        diag["ev_yes_pct"] = ev_yes_pct
        diag["ev_no_pct"]  = ev_no_pct

        # Defensive math-invariant check. For any LP making money on the
        # spread, yes_ask + no_ask + fees > 1, which means both gates can't
        # pass simultaneously against the same fair. If this ever fires the
        # fee model or fair-value pipeline has drifted — alert + decline.
        if yes_ok and no_ok:
            log.warning(
                "[MATH_INVARIANT_BROKEN] both sides +EV from one LP: "
                "yes_ask=%.4f no_ask=%.4f fair=%.4f ev_yes_pct=%.4f ev_no_pct=%.4f",
                1 - no_bid, 1 - yes_bid, fair, ev_yes_pct, ev_no_pct,
            )
            _log_quote_decision(quote, fair, "declined_math_invariant",
                                 diag=diag)
            return

        if not (yes_ok or no_ok):
            # Log against the (better of the two) EV so historic queries that
            # filter on post_fee_ev_pct still see meaningful numbers.
            _log_quote_decision(quote, fair, "declined_ev",
                                 post_fee_ev=max(ev_yes_pct, ev_no_pct),
                                 diag=diag)
            return

        chosen = "yes" if yes_ok else "no"
        chosen_ev_pct     = ev_yes_pct     if chosen == "yes" else ev_no_pct
        chosen_ev_dollars = ev_yes_dollars if chosen == "yes" else ev_no_dollars
        diag["chosen_side"] = chosen

        # Per-RFQ intended_side gate. Each RFQ is created sized for one
        # specific side (via two-RFQ Kelly create-time sizing). If the
        # maker's +EV side doesn't match this RFQ's intent, decline — the
        # sibling RFQ on the other side will handle that scenario at its
        # correct size. NULL intended_side = legacy pre-deploy RFQ; treat
        # as "no constraint" for backward compatibility.
        if intended_side is not None and intended_side != chosen:
            _log_quote_decision(quote, fair, "declined_side_mismatch",
                                 post_fee_ev=chosen_ev_pct, diag=diag)
            return

        # Side-aware cooldown — moved here from the pre-eval gate sweep so we
        # can look up the right (leg_set_hash, side) row. The other side
        # remains eligible for its own RFQ.
        with db.connect(read_only=True) as con:
            cd_rows = con.execute(
                "SELECT leg_set_hash, side, cooled_until FROM combo_cooldown"
            ).fetchall()
        cd_map = {(h, s): u for h, s, u in cd_rows}
        if not risk.cooldown_ok(leg_set_hash, chosen, cd_map):
            _log_quote_decision(quote, fair, "declined_cooldown",
                                 post_fee_ev=chosen_ev_pct, diag=diag)
            return

        # Hedge-formation diagnostic (non-blocking). If we already hold the
        # opposite side on this combo, tag this fill so we can monitor whether
        # across-time hedges net to a real loss pattern. Forward-looking +EV
        # math already says taking this side is correct given the held side;
        # we just want to *measure* the cumulative effect.
        opposite = "no" if chosen == "yes" else "yes"
        with db.connect(read_only=True) as con:
            held = con.execute(
                "SELECT net_contracts, weighted_price FROM positions "
                "WHERE combo_market_ticker=? AND side=? "
                "AND net_contracts > 0 LIMIT 1",
                [combo_market_ticker, opposite]
            ).fetchone()
        if held:
            held_n, held_price = held
            new_ask = (1 - no_bid) if chosen == "yes" else (1 - yes_bid)
            diag["hedge_added"]          = True
            diag["hedge_original_side"]  = opposite
            diag["hedge_original_price"] = float(held_price)
            diag["hedge_new_price"]      = new_ask
            diag["hedge_current_fair"]   = fair
            # Combined net P&L at settlement on the matched contracts, before
            # fees: $1 payout minus the sum of asks paid on the two sides.
            diag["hedge_projected_net"]  = 1.0 - held_price - new_ask

        if dry_run:
            _log_quote_decision(quote, fair, "declined_dry_run",
                                 post_fee_ev=chosen_ev_pct, diag=diag)
            return

        # Kelly is a go/no-go gate at $1 RFQs (sizing happens upstream via
        # target_cost_dollars). Pass `chosen` so we evaluate the right side.
        contracts = _kelly_size_for_quote(quote, fair, side=chosen)
        if contracts <= 0:
            _log_quote_decision(quote, fair, "declined_kelly_zero",
                                 post_fee_ev=chosen_ev_pct, diag=diag)
            return

        # Side semantics (INVERTED FROM INTUITION — verified empirically
        # 2026-05-21 from a real fill):
        #   accepted_side="yes" → Kalshi makes us BUY NO at no_ask
        #   accepted_side="no"  → Kalshi makes us BUY YES at yes_ask
        # The field names the side of the LP's two-sided quote we're
        # accepting (the side the LP is bidding on); by accepting it we SELL
        # them that side, leaving us long the opposite. So to open a LONG
        # YES position we send accepted_side="no"; for LONG NO we send
        # accepted_side="yes". See README "Accept semantics".
        accepted_side = "no" if chosen == "yes" else "yes"
        diag["accept_attempted_at"] = datetime.now(timezone.utc)
        resp, err_body = rfq_client.accept_quote(quote["id"], accepted_side=accepted_side)
        diag["accept_response_at"] = datetime.now(timezone.utc)
        if resp is None:
            # Walked. Capture Kalshi's error body + a fresh look at the RFQ's
            # terminal state. get_rfq_safe swallows any errors so diagnostics
            # never crash the bot.
            try:
                diag["accept_response_body"] = (
                    json.dumps(err_body) if err_body is not None else None
                )
            except (TypeError, ValueError):
                diag["accept_response_body"] = str(err_body)
            rfq_obj = rfq_client.get_rfq_safe(rfq_id)
            diag["rfq_terminal_status"] = (
                rfq_obj.get("status") if rfq_obj else None
            )
            research.emit("walk_diagnosed", game_id=game_id,
                          combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                          quote_id=quote.get("id"),
                          accept_response_body=diag.get("accept_response_body"),
                          rfq_terminal_status=diag.get("rfq_terminal_status"))
            _log_quote_decision(quote, fair, "failed_quote_walked",
                                 post_fee_ev=chosen_ev_pct, diag=diag)
            return

        # Post-accept fill reconciliation via /portfolio/positions.
        # `get_position_contracts` returns a SIGNED int: positive = long YES,
        # negative = long NO (short YES). The absolute value is our actual
        # position size; the sign is a sanity check against `chosen`.
        # On API failure we fall back to the offered-size field on the side
        # we filled against — for chosen="yes" we accepted="no" so the LP's
        # NO-side offer is what filled (no_contracts_fp), and vice versa.
        try:
            signed = rfq_client.get_position_contracts(combo_market_ticker)
            _record_positions_api_result(True)
            if (chosen == "yes" and signed < 0) or (chosen == "no" and signed > 0):
                log.warning(
                    "[position_direction_mismatch] chosen=%s but signed_position=%s on %s",
                    chosen, signed, combo_market_ticker,
                )
            actual = abs(signed)
        except Exception:
            _record_positions_api_result(False)
            if chosen == "yes":
                fallback_fp = quote.get("no_contracts_fp") or quote.get("contracts_fp")
            else:
                fallback_fp = quote.get("yes_contracts_fp") or quote.get("contracts_fp")
            actual = int(float(fallback_fp)) if fallback_fp else 0

        ask = (1.0 - no_bid) if chosen == "yes" else (1.0 - yes_bid)
        fee = ev_calc.fee_per_contract(ask)
        with db.connect() as con:
            con.execute(
                "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                "game_id, side, contracts, price_dollars, fee_dollars, "
                "blended_fair_at_fill, expected_ev_dollars, filled_at, raw_response) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [str(uuid.uuid4()), quote["id"], rfq_id, combo_market_ticker, game_id,
                 chosen, actual, ask, fee, fair, chosen_ev_dollars,
                 datetime.now(timezone.utc), str(resp)],
            )
            con.execute(
                "INSERT INTO positions (combo_market_ticker, side, game_id, "
                "net_contracts, weighted_price, legs_json, updated_at) VALUES "
                "(?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT (combo_market_ticker, side) DO UPDATE SET "
                "net_contracts = positions.net_contracts + EXCLUDED.net_contracts, "
                "weighted_price = (positions.weighted_price * positions.net_contracts + "
                "                  EXCLUDED.weighted_price * EXCLUDED.net_contracts) / "
                "                 (positions.net_contracts + EXCLUDED.net_contracts), "
                "updated_at = EXCLUDED.updated_at",
                [combo_market_ticker, chosen, game_id, actual, ask, legs_json,
                 datetime.now(timezone.utc)],
            )
            cooled_until = datetime.now(timezone.utc) + timedelta(seconds=config.COMBO_COOLDOWN_SEC)
            con.execute(
                "INSERT INTO combo_cooldown (leg_set_hash, side, game_id, cooled_until, reason) "
                "VALUES (?, ?, ?, ?, ?) "
                "ON CONFLICT (leg_set_hash, side) DO UPDATE SET "
                "cooled_until = EXCLUDED.cooled_until",
                [leg_set_hash, chosen, game_id, cooled_until, "post_accept"],
            )
        # Persist post-fill position snapshot to the research firehose
        # (positions row is mutated in place, so this captures the history).
        # Wrapped: the fill is already committed — a read-back error must not
        # bubble past the accept and skip _log_quote_decision / notify.fill.
        try:
            with db.connect(read_only=True) as con:
                after = con.execute(
                    "SELECT net_contracts, weighted_price FROM positions "
                    "WHERE combo_market_ticker=? AND side=?",
                    [combo_market_ticker, diag["chosen_side"]]).fetchone()
        except Exception:
            after = None
        research.emit("position_snapshot", game_id=game_id,
                      combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                      quote_id=quote.get("id"),
                      side=diag["chosen_side"], contracts_added=actual,
                      net_contracts_after=(after[0] if after else None),
                      weighted_price_after=(after[1] if after else None))
        _log_quote_decision(quote, fair, "accepted",
                             post_fee_ev=chosen_ev_pct, diag=diag)
        notify.fill(rfq_id=rfq_id, combo_market_ticker=combo_market_ticker,
                    contracts=actual, price=ask, ev_pct=chosen_ev_pct)


# ------------------------------------------------------------------------ #
# RFQ refresh — continuous priority-queue pipeline                          #
# ------------------------------------------------------------------------ #

def mint_and_create_rfq(candidate: combo_enumerator.ComboCandidate,
                         target_cost_dollars: float,
                         replace_existing: bool = False) -> tuple[str, str]:
    """Mint combo ticker (or cache hit) + create RFQ with a dollar budget.
    Returns (rfq_id, combo_ticker).

    Sizes via `target_cost_dollars` (= Kelly's contract count × worst-
    acceptable price for the intended side). When the maker quotes better
    than worst-case, the dollar budget naturally translates into more
    contracts at the lower price, closer to ideal Kelly. The MIN_EV_PCT
    gate still rejects quotes worse than the EV cliff, and the per-RFQ
    intended_side gate (in _evaluate_quote) declines quotes that favor
    the opposite side.

    `replace_existing` bypasses Kalshi's per-(user, market_ticker) RFQ
    dedup — required for the second-side RFQ on the same combo ticker.
    See rfq_client.create_rfq for the (misnamed-flag) details.
    """
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

    rfq_id = rfq_client.create_rfq(
        combo_ticker,
        target_cost_dollars=target_cost_dollars,
        replace_existing=replace_existing,
    )
    return rfq_id, combo_ticker


def _refresh_rfqs(candidates: list[combo_enumerator.ComboCandidate],
                   fair_scores: dict[str, tuple[float, float]],
                   kelly_sizes: dict[str, tuple[int, int, float, float]],
                   dry_run: bool):
    """Continuous priority-queue pipeline. Per candidate, up to TWO RFQs
    are sent — one per side with positive Kelly — each tagged with its
    intended_side.

    Args:
      fair_scores: leg_set_hash → (blended_fair, kalshi_ref).
      kelly_sizes: leg_set_hash → (kelly_yes_n, kelly_no_n, worst_yes_ask,
        worst_no_ask). Sides with 0 Kelly count or 0 worst-ask are skipped.
    """
    scored = [(c, *fair_scores[c.leg_set_hash])
              for c in candidates if c.leg_set_hash in fair_scores]
    ranked = combo_enumerator.rank_by_edge(scored)
    target = ranked[: _max_live_rfqs()]
    target_hashes = {c.leg_set_hash for c in target}

    # Dedup is now keyed by (leg_set_hash, intended_side). NULL legacy rows
    # (pre-deploy RFQs) are treated as having no constraint and are kept
    # alive — we don't reissue against them.
    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT rfq_id, leg_set_hash, intended_side FROM live_rfqs "
            "WHERE status='open'"
        ).fetchall()
    live_pairs: dict[tuple[str, str | None], str] = {
        (h, s): rid for rid, h, s in live
    }

    # Drop: any live RFQ whose leg_set_hash isn't in the new target list.
    # Iterate over a copy because we mutate the dict for retries.
    for (h, _side), rid in list(live_pairs.items()):
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
                log.warning("drop %s failed: %s", rid, e)

    # Add: for each target candidate, emit up to 2 RFQs (one per side with
    # positive Kelly). Skip any (leg_set_hash, side) pair already live.
    skipped_kelly_zero = 0
    skipped_too_small = 0
    added = 0
    for c in target:
        ks = kelly_sizes.get(c.leg_set_hash)
        if ks is None:
            continue
        yes_n, no_n, yes_ask, no_ask = ks
        sides: list[tuple[str, int, float]] = []
        if yes_n > 0 and yes_ask > 0:
            sides.append(("yes", yes_n, yes_ask))
        if no_n > 0 and no_ask > 0:
            sides.append(("no", no_n, no_ask))
        if not sides:
            skipped_kelly_zero += 1
            continue
        blended_fair, kalshi_ref = fair_scores[c.leg_set_hash]
        edge = combo_enumerator.edge_score(blended_fair, kalshi_ref)
        # Track whether ANY side already has a live RFQ on this combo
        # ticker — the second side we send must use replace_existing=True
        # to bypass Kalshi's per-(user, market_ticker) RFQ dedup. Probed
        # 2026-05-24: with replace_existing=True both RFQs coexist as
        # status=open; with False the second 409s.
        combo_has_live_rfq = any(
            (c.leg_set_hash, _s) in live_pairs for _s in ("yes", "no")
        )
        sent_this_cycle_for_combo = combo_has_live_rfq
        for side, n, ask in sides:
            if (c.leg_set_hash, side) in live_pairs:
                continue
            target_cost = round(n * ask, 2)
            if target_cost < 0.01:
                skipped_too_small += 1
                continue
            try:
                rid, combo_ticker = mint_and_create_rfq(
                    c,
                    target_cost_dollars=target_cost,
                    replace_existing=sent_this_cycle_for_combo,
                )
                sent_this_cycle_for_combo = True
                with db.connect() as con:
                    con.execute(
                        "INSERT INTO live_rfqs (rfq_id, combo_market_ticker, "
                        "leg_set_hash, game_id, blended_fair_at_submit, "
                        "kalshi_ref_at_submit, edge_at_submit, intended_side, "
                        "kelly_yes_n_at_submit, kelly_no_n_at_submit, "
                        "worst_yes_ask_at_submit, worst_no_ask_at_submit, "
                        "target_cost_dollars_at_submit, "
                        "status, submitted_at) "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        [rid, combo_ticker, c.leg_set_hash, c.game_id,
                         blended_fair, kalshi_ref, edge, side,
                         int(yes_n), int(no_n),
                         float(yes_ask), float(no_ask),
                         float(target_cost),
                         "open", datetime.now(timezone.utc)],
                    )
                added += 1
            except Exception as e:
                log.warning("add %s/%s failed: %s", c.leg_set_hash[:8], side, e)
                research.emit("rfq_submit_failed",
                              game_id=getattr(c, "game_id", None),
                              leg_set_hash=c.leg_set_hash,
                              side=side, error=str(e))
    if skipped_kelly_zero or skipped_too_small or added:
        log.info("rfq_refresh: add=%d skipped_kelly_zero=%d skipped_too_small=%d",
                 added, skipped_kelly_zero, skipped_too_small)


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
        # Stamp "we first saw this batch" once per poll_quotes call. All quotes
        # in `quotes` arrived in the same response, so they share the same
        # observed-arrival time from the bot's perspective.
        first_seen_at = datetime.now(timezone.utc)
        open_quotes = [q for q in quotes if q.get("status") == "open"]
        for q in open_quotes:
            # Competitor stats relative to *this* quote: every OTHER open quote
            # on the same RFQ. Used in walk diagnostics to tell whether a
            # cheaper/better quote was sitting right there when we walked on
            # this one.
            others = [oq for oq in open_quotes if oq.get("id") != q.get("id")]
            competitor_count = len(others)
            other_no_bids = [
                float(oq["no_bid_dollars"]) for oq in others
                if oq.get("no_bid_dollars") is not None
            ]
            # Higher no_bid_dollars ⇒ lower yes_ask ⇒ better for the buyer
            # (this bot is the RFQ creator buying YES). "Best" = max.
            best_competitor_no_bid = max(other_no_bids) if other_no_bids else None

            with db.connect(read_only=True) as con:
                already = con.execute(
                    "SELECT 1 FROM quote_log WHERE quote_id=?", [q["id"]]
                ).fetchone()
            if already:
                continue
            _evaluate_quote(
                q,
                dry_run=dry_run,
                competitor_count=competitor_count,
                best_competitor_no_bid_dollars=best_competitor_no_bid,
                quote_first_seen_at=first_seen_at,
            )


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
    """Map Kalshi 3-letter codes to game_id via in-memory mlb_parlay_lines cache."""
    home = _MLB_CODE_TO_TEAM.get(home_code)
    away = _MLB_CODE_TO_TEAM.get(away_code)
    if not home or not away:
        return None
    now = datetime.now(timezone.utc)
    horizon = now + timedelta(hours=24)
    for game_id, row in _PARLAY_LINES_CACHE.items():
        if row["home_team"] != home or row["away_team"] != away:
            continue
        ct = row["commence_time"]
        if ct is None:
            continue
        if ct.tzinfo is None:
            ct = ct.replace(tzinfo=timezone.utc)
        if now < ct < horizon:
            return game_id
    return None


def _enumerate_and_score_all_games() -> tuple[list[combo_enumerator.ComboCandidate],
                                                dict[str, tuple[float, float]],
                                                dict[str, tuple[int, int, float, float]]]:
    """Returns (candidates, fair_scores, kelly_sizes) where:
      - fair_scores[leg_set_hash] = (blended_fair, kalshi_ref) for ranking
      - kelly_sizes[leg_set_hash] = (kelly_yes_n, kelly_no_n,
        worst_yes_ask, worst_no_ask) — per-side Kelly state for two-RFQ
        creation in _refresh_rfqs. Candidates with both sides zero are
        still ranked (so the cycle counts what was considered), but
        _refresh_rfqs will skip them at create time.
    """
    status, body, _ = auth_client.api(
        "GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    events = body.get("events", []) if status == 200 and isinstance(body, dict) else []

    candidates_all: list[combo_enumerator.ComboCandidate] = []
    fair_scores: dict[str, tuple[float, float]] = {}
    kelly_sizes: dict[str, tuple[int, int, float, float]] = {}

    for ev in events:
        event_ticker = ev.get("event_ticker", "")
        if not event_ticker.startswith("KXMLBGAME-"):
            continue
        suffix = event_ticker.replace("KXMLBGAME-", "")
        away_code, home_code = _parse_event_suffix(suffix)
        if away_code is None or home_code is None:
            continue

        avail_spreads = _kalshi_available_spreads(suffix, home_code, away_code)
        avail_totals = _kalshi_available_totals(suffix)
        if not avail_spreads or not avail_totals:
            continue

        game_id = _resolve_game_id(home_code, away_code)
        if game_id is None:
            continue

        # (book-only) do not skip games without samples
        if config.USE_MODEL:
            samples = _load_samples_for_game(game_id)
            if samples is None:
                continue
        else:
            samples = None

        for cand in combo_enumerator.enumerate_2leg(
                game_id=game_id, event_suffix=suffix,
                home_code=home_code, away_code=away_code,
                available_spreads=avail_spreads,
                available_totals=avail_totals):
            typed = [_leg_dict_to_typed(dict(l), game_id) for l in cand.legs]
            spread_line = _spread_line_from_legs([dict(l) for l in cand.legs])
            total_line = _total_line_from_legs([dict(l) for l in cand.legs])

            if any(l is None for l in typed):
                _emit_candidate_event("rejected_no_mapping", game_id=game_id,
                                      cand=cand, spread_line=spread_line,
                                      total_line=total_line)
                continue
            region = _combo_region_from_legs(typed)
            if region is None:
                _emit_candidate_event("rejected_no_mapping", game_id=game_id,
                                      cand=cand, spread_line=spread_line,
                                      total_line=total_line)
                continue
            books = _load_book_fairs(game_id, region.spread_line, region.total_line,
                                     region.spread_side, region.total_side)
            if config.USE_MODEL:
                model = fair_value.model_fair(samples, typed)
                blended = fair_value.blend(model, books)
            else:
                model = None
                blended = _book_only_fair(books)
            if blended is None:
                _emit_candidate_event("rejected_no_book_data", game_id=game_id,
                                      cand=cand, spread_line=spread_line,
                                      total_line=total_line,
                                      model=model, books=books)
                continue
            if not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
                _emit_candidate_event("rejected_fair_oob", game_id=game_id,
                                      cand=cand, spread_line=spread_line,
                                      total_line=total_line,
                                      model=model, books=books, blended=blended)
                continue
            kalshi_ref = _kalshi_last_price(cand.legs[0]["market_ticker"])
            yes_n, no_n, yes_ask, no_ask = _kelly_size_for_candidate(
                game_id, typed, samples, blended,
                leg_set_hash=cand.leg_set_hash)
            _emit_candidate_event("submitted", game_id=game_id, cand=cand,
                                  spread_line=spread_line, total_line=total_line,
                                  model=model, books=books, blended=blended,
                                  kalshi_ref=kalshi_ref,
                                  kelly=(yes_n, no_n, yes_ask, no_ask))
            candidates_all.append(cand)
            fair_scores[cand.leg_set_hash] = (blended, kalshi_ref)
            kelly_sizes[cand.leg_set_hash] = (yes_n, no_n, yes_ask, no_ask)

    return candidates_all, fair_scores, kelly_sizes


# ------------------------------------------------------------------------ #
# Main loop                                                                #
# ------------------------------------------------------------------------ #

def main_loop(dry_run: bool):
    setup_logging()
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run, version=VERSION)
    research.init_research_db()
    research.set_session(sid)
    log.info("=== Kalshi MLB RFQ Bot — session %s (dry_run=%s) ===", sid, dry_run)
    # Initial cache load — bot is useless until this succeeds.
    if not _refresh_caches():
        log.warning("startup: cache_refresh failed; bot will retry on pipeline-refresh tick")
    _phantom_rfq_cleanup()

    # Synchronous SGP warm-up: run one full SGP cycle before entering the
    # loop so _SGP_ODDS_CACHE is non-empty when the first RFQ refresh fires.
    # Blocks ~60-90s; without this the first RFQ refresh would have zero
    # book fairs and burn Kalshi quota on candidates that can't pass the
    # N>=2 books gate.
    from kalshi_common.sgp_service import SGPService
    sgp_service = SGPService(per_book_deadline_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
    log.info("startup: warming SGP cache (one synchronous scrape tick)...")
    try:
        rcs = sgp_runner.sgp_cycle(
            bot_market_db=str(config.BOT_MARKET_DB),
            service=sgp_service,
            both_teams=True,
        )
        log.info("startup: SGP warm-up done — return codes %s", rcs)
    except Exception as e:
        log.warning("startup: SGP warm-up failed (%s); bot will retry on first cadence tick", e)
    _refresh_sgp_cache()
    # Rebuild parlay_lines cache now that target_lines has rows
    _refresh_caches()

    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0
    last_sgp_cycle = time.time()  # warm-up just finished

    try:
        while _running.is_set():
            now = time.time()

            if KILL_FILE.exists():
                notify.halt("kill_switch")
                time.sleep(config.RISK_SWEEP_SEC)
                continue

            if now - last_rfq_refresh >= config.RFQ_REFRESH_SEC:
                t_ref = time.time()
                try:
                    candidates, fair_scores, kelly_sizes = _enumerate_and_score_all_games()
                    _refresh_rfqs(candidates, fair_scores, kelly_sizes, dry_run=dry_run)
                    log.info("rfq_refresh: %d candidates (%.1fs)",
                             len(candidates), time.time() - t_ref)
                except Exception as e:
                    log.warning("rfq_refresh error: %s", e)
                last_rfq_refresh = now

            if now - last_quote_poll >= config.QUOTE_POLL_SEC:
                try:
                    _poll_all_live_rfqs(dry_run=dry_run)
                except Exception as e:
                    log.warning("quote_poll error: %s", e)
                last_quote_poll = now

            if now - last_risk_sweep >= config.RISK_SWEEP_SEC:
                try:
                    _risk_sweep()
                except Exception as e:
                    log.warning("risk_sweep error: %s", e)
                last_risk_sweep = now

            if now - last_sgp_cycle >= config.SGP_REFRESH_SEC:
                t_sgp = time.time()
                try:
                    rcs = sgp_runner.sgp_cycle(
                        bot_market_db=str(config.BOT_MARKET_DB),
                        service=sgp_service,
                        both_teams=True,
                    )
                    _refresh_sgp_cache()
                    _refresh_caches()  # parlay_lines_cache may have new (spread, total) tuples
                    log.info("sgp_cycle: rcs=%s (%.1fs)", rcs, time.time() - t_sgp)
                except Exception as e:
                    log.warning("sgp_cycle error: %s", e)
                last_sgp_cycle = now

            if now - last_pipeline >= config.PIPELINE_REFRESH_SEC:
                _run_pipeline()
                # Reload caches after the pipeline has had a chance to refresh data.
                # Retries inside _refresh_caches handle any brief lock window.
                _refresh_caches()
                last_pipeline = now

            if now - last_heartbeat >= 60:
                log.info("[HB] %s alive", datetime.now(timezone.utc).isoformat())
                last_heartbeat = now

            research.flush()   # persist this tick's buffered research events
            time.sleep(0.5)
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT rfq_id FROM live_rfqs WHERE status='open'"
            ).fetchall()]
        cancelled = 0
        failed = 0
        for rid in live:
            try:
                rfq_client.delete_rfq(rid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_rfqs SET status='cancelled', closed_at=?, "
                        "cancellation_reason='shutdown' WHERE rfq_id=?",
                        [datetime.now(timezone.utc), rid],
                    )
                cancelled += 1
            except Exception as e:
                failed += 1
                log.warning("shutdown: cancel %s failed: %s", rid, e)
        if live:
            log.info("shutdown: drained live RFQs — cancelled=%d failed=%d",
                     cancelled, failed)
        research.flush()   # persist the final tick's buffered research events
        db.end_session(sid)
        log.info("=== shutdown complete ===")


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
