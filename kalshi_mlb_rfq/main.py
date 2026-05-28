"""Kalshi MLB RFQ Bot — autonomous taker daemon."""

import argparse
import json
import logging
import os
import random as _random
import signal
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
    auth_client, combo_enumerator, config, db, ev_calc,
    fair_value, kelly, notify, research, rfq_client, risk, sgp_runner,
)
from kalshi_mlb_rfq.config import (
    ANSWER_KEY_DB, KILL_FILE, MAX_BOOK_STALENESS_SEC,
)
from kalshi_mlb_rfq.log_setup import setup_logging

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
            # mlb_game_samples
            t0 = time.time()
            samples_df = con.execute(
                "SELECT game_id, sim_idx, home_margin, total_final_score, "
                "home_margin_f5, total_f5 FROM mlb_game_samples"
            ).fetchdf()
            samples_by_game = {gid: g.reset_index(drop=True)
                                for gid, g in samples_df.groupby("game_id")}

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

            # samples meta — generated_at
            meta_row = con.execute(
                "SELECT generated_at FROM mlb_samples_meta "
                "ORDER BY generated_at DESC LIMIT 1"
            ).fetchone()
            generated_at = meta_row[0] if meta_row else None
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


# 3-letter Kalshi team code → mlb_parlay_lines.home_team / away_team canonical name.
# Kalshi uses 3-letter codes; mlb_parlay_lines stores Odds-API canonical names.
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "ATH": "Athletics",
    "AZ": "Arizona Diamondbacks", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


def _parse_event_suffix(suffix: str) -> tuple[str | None, str | None]:
    """Split a KXMLB* event suffix into (away_code, home_code).

    Format: YYMMMDDHHMM{AwayCode}{HomeCode}. Date prefix is fixed at 11 chars.
    Each team code is 2 or 3 letters (KC/SF/SD/TB/AZ are 2-letter; the rest
    are 3-letter). Probes 3- then 2-letter home splits and returns the first
    where both codes are valid in _MLB_CODE_TO_TEAM. Returns (None, None) if
    no split matches — caller drops the event.
    """
    if len(suffix) < 11 + 4:  # date prefix + at least 2+2 team chars
        return None, None
    team_block = suffix[11:]
    for home_len in (3, 2):
        if len(team_block) <= home_len:
            continue
        home = team_block[-home_len:]
        away = team_block[:-home_len]
        if home in _MLB_CODE_TO_TEAM and away in _MLB_CODE_TO_TEAM:
            return away, home
    return None, None


def _load_samples_for_game(game_id: str) -> pd.DataFrame | None:
    """Read from in-memory cache, populated by _refresh_caches()."""
    df = _SAMPLES_CACHE.get(game_id)
    return df if df is not None and not df.empty else None


def _load_book_fairs(game_id: str, spread_line: float, total_line: float) -> dict[str, float]:
    """Per-line book lookup with N>=2 gate.

    Returns {book -> devigged_fair} only if at least MIN_BOOK_COUNT_FOR_BLEND
    books priced the matching (game, spread, total) tuple. Empty dict
    otherwise — caller treats as 'no book signal, drop candidate'.
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
    out: dict[str, float] = {}
    for book in rows["bookmaker"].unique():
        sub = rows[rows["bookmaker"] == book].copy()
        fair_per_book = fair_value.devig_book(
            sub, combo="Home Spread + Over",
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
    }.get(book, 0.10)


def _home_code_from_event_ticker(event_ticker: str) -> str | None:
    """Parse the home-team code from a Kalshi event ticker (2- or 3-letter)."""
    if "-" not in event_ticker:
        return None
    suffix = event_ticker.rsplit("-", 1)[-1]
    _, home = _parse_event_suffix(suffix)
    return home


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
    outcome_vec = _outcome_vec_for_legs(samples, typed_legs, side=side)
    if outcome_vec is None:
        return 0

    existing = _load_existing_positions_for_game(game_id, samples)

    no_bid  = float(quote.get("no_bid_dollars")  or 0)
    yes_bid = float(quote.get("yes_bid_dollars") or 0)
    if side == "yes":
        ask = 1 - no_bid
    else:
        ask = 1 - yes_bid
    fee = ev_calc.fee_per_contract(ask)
    effective_price = ask + fee

    return kelly.kelly_size_combo(
        outcome_vec=np.asarray(outcome_vec),
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
                                samples: pd.DataFrame,
                                blended_fair: float,
                                ) -> tuple[int, int, float, float]:
    """Kelly counts and worst-acceptable prices for BOTH sides at create time.

    Returns (kelly_yes_n, kelly_no_n, worst_yes_ask, worst_no_ask). Any side
    with no acceptable price returns (n=0, ask=0.0) for that side.

    Sizing rationale per side:
      - Compute worst-acceptable ask via binary search at MIN_EV_PCT floor.
      - Build the side's outcome_vec (combo hits for YES, misses for NO).
      - Kelly with the side's effective_price + (blended_fair for YES,
        1-blended_fair for NO).
    Existing positions on the same game are loaded once and reused; each
    row's outcome_vec is inverted if its stored side is "no" (handled by
    _load_existing_positions_for_game).
    """
    if samples is None or samples.empty:
        return 0, 0, 0.0, 0.0

    existing = _load_existing_positions_for_game(game_id, samples)

    def kelly_for_side(side: str) -> tuple[int, float]:
        outcome_vec = _outcome_vec_for_legs(samples, typed_legs, side=side)
        if outcome_vec is None:
            return 0, 0.0
        ask = _worst_acceptable_ask(
            blended_fair, side, config.KELLY_CREATE_EV_FLOOR_PCT)
        if ask <= 0:
            return 0, 0.0
        fee = ev_calc.fee_per_contract(ask)
        effective_price = ask + fee
        p = blended_fair if side == "yes" else (1.0 - blended_fair)
        n = kelly.kelly_size_combo(
            outcome_vec=np.asarray(outcome_vec),
            blended_fair=p,
            existing_positions=existing,
            effective_price=effective_price,
            bankroll=config.BANKROLL,
            kelly_fraction=config.KELLY_FRACTION,
        )
        return n, ask

    yes_n, yes_ask = kelly_for_side("yes")
    no_n,  no_ask  = kelly_for_side("no")
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
        fair = _fresh_blended_fair(combo_market_ticker)
        if fair is None:
            _log_quote_decision(quote, None, "declined_ev", reason="no_fresh_fair", diag=diag)
            return

        passed, decision = _all_per_accept_gates_pass(
            quote, fair, {"leg_set_hash": leg_set_hash, "game_id": game_id,
                          "legs_json": legs_json})
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

        samples = _load_samples_for_game(game_id)
        if samples is None:
            continue

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
            model = fair_value.model_fair(samples, typed)
            books = _load_book_fairs(game_id, spread_line, total_line)
            blended = fair_value.blend(model, books)
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
                game_id, typed, samples, blended)
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
    log.info("startup: warming SGP cache (one synchronous scrape tick)...")
    try:
        rcs = sgp_runner.sgp_cycle(
            bot_market_db=str(config.BOT_MARKET_DB),
            scraper_dir=str(config.MLB_SGP_DIR),
            venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
            timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC,
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
                        scraper_dir=str(config.MLB_SGP_DIR),
                        venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                        timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC,
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
