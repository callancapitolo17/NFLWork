"""
Kalshi CBB 1H Taker — Crosses the spread when EV exceeds threshold.

Can run standalone (python taker.py --dry-run) when the MM is not running,
or be imported by main.py to run in the same process (required because
DuckDB only allows one process to access a file at a time).

Usage (standalone — only when MM is NOT running):
    python taker.py              # Live mode
    python taker.py --dry-run    # Scan only, don't place orders
"""

import sys
import time
import signal
import uuid
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import config
import db
import orders
import risk

# --- Global state ---
# Cooldown: {ticker: prediction_updated_at} — skip ticker until predictions refresh
TAKEN_THIS_CYCLE = {}

# Prediction cache to avoid re-reading DB every 2s
PRED_CACHE = {"predictions": [], "updated_at": None, "checked_at": 0}


def compute_taker_fee(price_cents):
    """Compute Kalshi taker fee in cents for a given price.

    Fee = 7% * P * (1-P) * 100 (in cents per contract).
    """
    p = price_cents / 100.0
    return config.TAKER_FEE_RATE * p * (1 - p) * 100


def compute_take_ev(fair_cents, book_price, side):
    """Compute EV% for taking a side, after taker fees.

    Returns (ev_pct, fee_cents) or (None, None) if not +EV enough.
    """
    if book_price <= 0 or book_price >= 100:
        return None, None

    fee = compute_taker_fee(book_price)

    if side == "yes":
        edge = fair_cents - book_price - fee
        ev_pct = edge / book_price
    else:
        fair_no = 100 - fair_cents
        edge = fair_no - book_price - fee
        ev_pct = edge / book_price

    if ev_pct < config.MIN_TAKE_EV_PCT:
        return None, None

    return ev_pct, fee


def is_on_cooldown(ticker, current_pred_ts):
    """Check if ticker was already taken on this prediction cycle."""
    if ticker not in TAKEN_THIS_CYCLE:
        return False
    return TAKEN_THIS_CYCLE[ticker] == current_pred_ts


def set_cooldown(ticker, pred_ts):
    """Mark ticker as taken for this prediction cycle."""
    TAKEN_THIS_CYCLE[ticker] = pred_ts


def clear_cooldowns():
    """Clear all cooldowns (called when predictions refresh)."""
    TAKEN_THIS_CYCLE.clear()


def refresh_predictions():
    """Reload predictions from cbb.duckdb if changed. Cache for efficiency."""
    now = time.time()
    # Only check DB every 10 seconds (predictions refresh every ~10min)
    if now - PRED_CACHE["checked_at"] < 10:
        return PRED_CACHE["predictions"], PRED_CACHE["updated_at"]

    PRED_CACHE["checked_at"] = now
    predictions, updated_at = db.load_predictions()

    if updated_at != PRED_CACHE["updated_at"] and predictions:
        PRED_CACHE["predictions"] = predictions
        PRED_CACHE["updated_at"] = updated_at
        # Clear all cooldowns — fresh predictions mean fresh signals
        clear_cooldowns()
        print(f"  [TAKER] Predictions refreshed ({len(predictions)} rows, cooldowns cleared)")

    return PRED_CACHE["predictions"], PRED_CACHE["updated_at"]


def execute_take(market, side, price, ev_pct, fee_cents, pred_ts, dry_run=False):
    """Place a taker order (post_only=False), then cancel unfilled remainder.

    Returns number of contracts actually filled (0 if failed).
    """
    ticker = market["ticker"]
    take_id = f"take-{uuid.uuid4().hex[:8]}"

    print(f"  [TAKER] TAKE: {side.upper()} {config.TAKE_CONTRACT_SIZE}x @ {price}c on {ticker} "
          f"(fair={market['fair_prob']*100:.1f}c, EV={ev_pct:.1%}, fee={fee_cents:.1f}c)")

    if dry_run:
        set_cooldown(ticker, pred_ts)
        return 0

    # Place order with post_only=False to cross the spread
    result = orders.place_order(
        ticker, side, price, config.TAKE_CONTRACT_SIZE, post_only=False
    )

    if not result:
        return 0

    order_id = result.get("order_id", "")
    remaining = result.get("remaining_count", 0)
    filled = config.TAKE_CONTRACT_SIZE - remaining

    # Cancel unfilled remainder immediately — we're a taker, not a maker
    if remaining > 0:
        orders.cancel_order(order_id)
        print(f"    Cancelled {remaining} unfilled (filled {filled}/{config.TAKE_CONTRACT_SIZE})")
    elif filled > 0:
        print(f"    Fully filled {filled} contracts")

    # Record fill and update position
    if filled > 0:
        db.update_position(
            ticker, side, price, filled,
            event_ticker=market.get("event_ticker"),
            home_team=market.get("home_team"),
            away_team=market.get("away_team"),
        )
        db.record_fill(
            fill_id=f"{take_id}-fill",
            ticker=ticker, side=side, action="buy",
            price=price, count=filled,
            fee_cents=fee_cents * filled,
            order_id=order_id,
        )
        db.log_take(
            take_id=take_id, ticker=ticker,
            event_ticker=market.get("event_ticker", ""),
            side=side, price=price,
            fair_cents=market["fair_prob"] * 100,
            ev_pct=ev_pct, fee_cents=fee_cents,
            count_requested=config.TAKE_CONTRACT_SIZE,
            count_filled=filled, order_id=order_id,
        )

    # Cooldown regardless of fill count — don't re-attempt until predictions refresh
    set_cooldown(ticker, pred_ts)
    return filled


def run_take_cycle(quotable_markets, prediction_updated_at, dry_run=False):
    """Scan all markets for +EV takes and execute."""
    # Global risk checks
    is_fresh, pred_age = risk.check_staleness(prediction_updated_at)
    if not is_fresh:
        return  # Silent — don't spam logs every 2s

    exposure_ok, exposure = risk.check_exposure_limit()
    if not exposure_ok:
        return

    takes_this_cycle = 0

    for market in quotable_markets:
        ticker = market["ticker"]
        event_ticker = market.get("event_ticker", "")
        fair_cents = market["fair_prob"] * 100

        # Skip if already taken on this prediction cycle
        if is_on_cooldown(ticker, prediction_updated_at):
            continue

        # Skip tipoff proximity
        if not risk.check_tipoff_proximity(market.get("commence_time")):
            continue

        # Get current position
        pos = db.get_position(ticker)
        net = pos["net_yes"]
        event_net = db.get_event_net_position(event_ticker)

        yes_ask = market.get("book_ask", 0)
        yes_bid = market.get("book_bid", 0)

        # --- Evaluate YES take (buy YES at ask) ---
        if yes_ask > 0:
            proposed_net = net + config.TAKE_CONTRACT_SIZE
            if (abs(proposed_net) <= config.MAX_POSITION_PER_MARKET
                    and abs(event_net + config.TAKE_CONTRACT_SIZE) <= config.MAX_POSITION_PER_EVENT):
                ev_pct, fee = compute_take_ev(fair_cents, yes_ask, "yes")
                if ev_pct is not None:
                    filled = execute_take(market, "yes", yes_ask, ev_pct, fee,
                                          prediction_updated_at, dry_run=dry_run)
                    if filled > 0:
                        takes_this_cycle += 1
                        continue  # Don't also check NO for same ticker

        # --- Evaluate NO take (buy NO at 100-bid) ---
        if yes_bid > 0:
            no_price = 100 - yes_bid
            proposed_net = net - config.TAKE_CONTRACT_SIZE
            if (abs(proposed_net) <= config.MAX_POSITION_PER_MARKET
                    and abs(event_net - config.TAKE_CONTRACT_SIZE) <= config.MAX_POSITION_PER_EVENT):
                ev_pct, fee = compute_take_ev(fair_cents, no_price, "no")
                if ev_pct is not None:
                    filled = execute_take(market, "no", no_price, ev_pct, fee,
                                          prediction_updated_at, dry_run=dry_run)
                    if filled > 0:
                        takes_this_cycle += 1

    if takes_this_cycle > 0:
        print(f"  [TAKER] Executed {takes_this_cycle} takes this cycle")


# --- Standalone mode (only when MM is not running) ---

def main():
    """Run taker as standalone process. Only use when MM is NOT running."""
    from scraper import fetch_markets
    from canonical_match import load_team_dict, load_canonical_games
    from main import match_kalshi_to_predictions, enforce_monotonicity, refresh_book_data

    dry_run = "--dry-run" in sys.argv
    running = True

    def _signal_handler(sig, frame):
        nonlocal running
        running = False

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)

    print("=" * 60)
    print(f"  Kalshi CBB 1H Taker Bot (standalone)")
    print(f"  {'DRY RUN' if dry_run else 'LIVE'}")
    print(f"  Min Take EV: {config.MIN_TAKE_EV_PCT:.0%} (after {config.TAKER_FEE_RATE:.0%} fee)")
    print(f"  Take Size: {config.TAKE_CONTRACT_SIZE} contracts")
    print(f"  Poll: {config.TAKE_POLL_SEC}s | Cooldown: per prediction refresh")
    print(f"  Max position: {config.MAX_POSITION_PER_MARKET} | Max exposure: ${config.MAX_TOTAL_EXPOSURE_DOLLARS}")
    print("=" * 60)

    db.init_database()
    db.init_taker_tables()

    print("\nLoading predictions...")
    predictions, prediction_updated_at = db.load_predictions()
    if not predictions:
        print("No predictions. Run pipeline first:")
        print('  cd "Answer Keys" && python run.py --sport cbb')
        return

    pred_age = "?"
    if prediction_updated_at:
        ts = prediction_updated_at
        if ts.tzinfo is None:
            ts = ts.replace(tzinfo=timezone.utc)
        age_sec = (datetime.now(timezone.utc) - ts).total_seconds()
        pred_age = f"{age_sec:.0f}s"
    print(f"  Loaded {len(predictions)} predictions (age: {pred_age})")

    PRED_CACHE.update({
        "predictions": predictions,
        "updated_at": prediction_updated_at,
        "checked_at": time.time(),
    })

    print("Loading team name resolution...")
    team_dict = load_team_dict("cbb")
    canonical_games = load_canonical_games("cbb")

    print(f"Fetching Kalshi 1H spread markets ({config.SPREAD_SERIES})...")
    kalshi_markets = fetch_markets(config.SPREAD_SERIES)
    print(f"  Found {len(kalshi_markets)} open markets")

    quotable = match_kalshi_to_predictions(
        kalshi_markets, predictions, team_dict, canonical_games
    )
    enforce_monotonicity(quotable)
    print(f"  {len(quotable)} quotable markets")

    if not quotable:
        print("No markets matched predictions.")
        return

    last_full_refresh = time.time()

    print(f"\nStarting taker loop (poll every {config.TAKE_POLL_SEC}s)...\n")

    try:
        while running:
            now = time.time()

            try:
                predictions, prediction_updated_at = refresh_predictions()

                if now - last_full_refresh >= 60:
                    fresh_markets = fetch_markets(config.SPREAD_SERIES)
                    if fresh_markets:
                        kalshi_markets = fresh_markets
                    quotable = match_kalshi_to_predictions(
                        kalshi_markets, predictions, team_dict, canonical_games
                    )
                    enforce_monotonicity(quotable)
                    last_full_refresh = now
                else:
                    refresh_book_data(quotable)

                run_take_cycle(quotable, prediction_updated_at, dry_run=dry_run)

            except Exception as e:
                print(f"  Cycle error (will retry): {e}")

            time.sleep(config.TAKE_POLL_SEC)

    except KeyboardInterrupt:
        pass
    finally:
        print("\nTaker bot shut down.")


if __name__ == "__main__":
    main()
