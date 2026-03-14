"""
Risk controls for the market maker.
Enforces position limits, staleness checks, and kill switch logic.
"""

import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import db
from config import (
    MAX_POSITION_PER_MARKET, MAX_TOTAL_EXPOSURE_DOLLARS, MAX_MARKETS,
    MAX_TOTAL_DIRECTIONAL, MAX_STALENESS_SEC, LINE_MOVE_THRESHOLD,
    TIPOFF_PULLBACK_MIN, BOOKMAKER_SCRAPER, BET105_SCRAPER
)


def check_position_limit(ticker, proposed_delta):
    """Check if a proposed fill would exceed position limits.

    Returns True if the trade is allowed, False if it would breach limits.
    """
    pos = db.get_position(ticker)
    new_net = pos["net_yes"] + proposed_delta
    return abs(new_net) <= MAX_POSITION_PER_MARKET


def check_exposure_limit():
    """Check if total exposure is within limits.

    Returns (is_ok, current_exposure).
    """
    exposure = db.compute_total_exposure()
    return exposure < MAX_TOTAL_EXPOSURE_DOLLARS, exposure


def check_directional_limit():
    """Check if global net directional exposure is within limits.

    Prevents correlated risk across events: a sharp filling YES bids
    on 10 different games simultaneously creates massive directional
    exposure that per-event limits don't catch.

    Returns:
        (can_go_long, can_go_short) — True if allowed in that direction.
    """
    net = db.get_global_net_position()
    return net <= MAX_TOTAL_DIRECTIONAL, net >= -MAX_TOTAL_DIRECTIONAL


def check_market_count(current_count):
    """Check if we're quoting too many markets."""
    return current_count < MAX_MARKETS


def check_staleness(prediction_updated_at):
    """Check if predictions are too stale to quote.

    Args:
        prediction_updated_at: datetime when predictions were last generated

    Returns:
        (is_fresh, age_seconds)
    """
    if prediction_updated_at is None:
        return False, float("inf")

    # Normalize both timestamps to the same basis before comparing.
    # CBB.R may write naive local time or UTC — handle both.
    if prediction_updated_at.tzinfo is not None:
        # Timezone-aware: compare in UTC
        now = datetime.now(timezone.utc)
        pred_ts = prediction_updated_at
    else:
        # Naive (assumed local time from R): compare in local time
        now = datetime.now()
        pred_ts = prediction_updated_at

    age = (now - pred_ts).total_seconds()

    # Sanity: negative age means clock skew or timezone mismatch.
    # Treat as stale rather than trading on potentially wrong data.
    if age < -60:
        return False, age

    return age < MAX_STALENESS_SEC, age


def check_tipoff_proximity(commence_time):
    """Check if a game is too close to tipoff to quote.

    Args:
        commence_time: Game start time (datetime or ISO string)

    Returns:
        True if safe to quote, False if too close to tipoff.
    """
    if commence_time is None:
        return True  # Unknown tipoff, allow quoting

    now = datetime.now(timezone.utc)

    if isinstance(commence_time, str):
        from datetime import datetime as dt
        try:
            commence_time = dt.fromisoformat(commence_time.replace("Z", "+00:00"))
        except (ValueError, TypeError):
            return True

    if commence_time.tzinfo is None:
        commence_time = commence_time.replace(tzinfo=timezone.utc)

    minutes_to_tip = (commence_time - now).total_seconds() / 60
    return minutes_to_tip > TIPOFF_PULLBACK_MIN


def run_line_monitor():
    """Run Bookmaker + Bet105 scrapers and return current 1H lines.

    Returns list of dicts with home_team, away_team, market, line_value per game.
    """
    current_lines = []

    for scraper_path, book_name in [
        (BOOKMAKER_SCRAPER, "bookmaker"),
        (BET105_SCRAPER, "bet105"),
    ]:
        if not scraper_path.exists():
            continue

        try:
            # Run scraper silently
            result = subprocess.run(
                [sys.executable, str(scraper_path), "cbb"],
                capture_output=True, text=True, timeout=60,
                cwd=str(scraper_path.parent)
            )
            if result.returncode != 0:
                print(f"  Warning: {book_name} scraper failed: {result.stderr[:200]}")
                continue

            # Read 1H spreads from the scraper's DuckDB
            import duckdb
            db_path = str(scraper_path.parent / f"{book_name}.duckdb")
            if not Path(db_path).exists():
                continue

            conn = duckdb.connect(db_path, read_only=True)
            try:
                rows = conn.execute("""
                    SELECT home_team, away_team, market, home_spread as line_value
                    FROM cbb_odds
                    WHERE market = 'spreads_h1' AND home_spread IS NOT NULL
                """).fetchall()
                for row in rows:
                    current_lines.append({
                        "home_team": row[0],
                        "away_team": row[1],
                        "market": row[2],
                        "line_value": row[3],
                        "source": book_name,
                    })
            finally:
                conn.close()

        except subprocess.TimeoutExpired:
            print(f"  Warning: {book_name} scraper timed out")
        except Exception as e:
            print(f"  Warning: {book_name} monitor error: {e}")

    return current_lines


def detect_line_moves(current_lines, reference_lines):
    """Compare current lines to reference lines and detect moves.

    Args:
        current_lines: List of current line dicts from run_line_monitor()
        reference_lines: List of reference line dicts from db.get_reference_lines()

    Returns:
        List of (home_team, away_team) tuples where lines moved > threshold.
    """
    moved_games = []

    # Build reference lookup: (home_team, away_team) -> line_value
    ref_lookup = {}
    for ref in reference_lines:
        key = (ref["home_team"], ref["away_team"])
        ref_lookup[key] = ref["line_value"]

    for line in current_lines:
        key = (line["home_team"], line["away_team"])
        if key in ref_lookup:
            delta = abs(line["line_value"] - ref_lookup[key])
            if delta > LINE_MOVE_THRESHOLD:
                moved_games.append(key)
                print(f"  LINE MOVE: {key[1]} vs {key[0]}: "
                      f"{ref_lookup[key]} → {line['line_value']} (delta={delta:.1f})")

    return list(set(moved_games))
