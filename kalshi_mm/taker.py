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
import kelly
import orders
import risk

# --- Global state ---
# Per-ticker tactical cooldown: {ticker: last_take_timestamp}
# Prevents hammering the same contract. Expires after 10s.
# Cross-market correlation is handled by Kelly conditional sizing.
TICKER_COOLDOWN = {}
_TICKER_COOLDOWN_SEC = 10

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


def is_on_cooldown(ticker):
    """Check if ticker was recently taken (10s tactical cooldown)."""
    last_take = TICKER_COOLDOWN.get(ticker)
    if last_take and (time.time() - last_take) < _TICKER_COOLDOWN_SEC:
        return True
    return False


def set_cooldown(ticker):
    """Mark ticker as recently taken."""
    TICKER_COOLDOWN[ticker] = time.time()


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
        print(f"  [TAKER] Predictions refreshed ({len(predictions)} rows)")

    return PRED_CACHE["predictions"], PRED_CACHE["updated_at"]


# Cached book snapshot to avoid multiple API calls per cycle
_BOOK_CACHE = {"data": {}, "fetched_at": 0}
_BOOK_CACHE_TTL = 1.5  # seconds — fresh enough for taker, saves rate limit


def _confirm_book_price(ticker, side):
    """Re-fetch the book for a single ticker to confirm price before taking.

    Returns (yes_bid, yes_ask) in cents, or (0, 0) if fetch fails.
    Uses a short-lived cache so multiple takes in one cycle share one API call.
    """
    now = time.time()
    if now - _BOOK_CACHE["fetched_at"] > _BOOK_CACHE_TTL:
        try:
            from scraper import fetch_markets
            new_data = {}
            for mtype, series in config.MARKET_SERIES.items():
                if mtype not in config.ENABLED_MARKET_TYPES:
                    continue
                fresh = fetch_markets(series)
                if fresh:
                    for m in fresh:
                        new_data[m["ticker"]] = (
                            int(round(float(m.get("yes_bid_dollars", 0)) * 100)),
                            int(round(float(m.get("yes_ask_dollars", 0)) * 100))
                        )
            if new_data:
                _BOOK_CACHE["data"] = new_data
                _BOOK_CACHE["fetched_at"] = now
        except Exception:
            pass

    result = _BOOK_CACHE["data"].get(ticker)
    return result if result else (0, 0)


def _reconcile_fill_count(order_id, total_requested):
    """Query the API for the actual fill count on an order.

    Handles the gap between place_order response and cancel: additional
    fills can happen in that window. This is the source of truth.
    """
    order = orders.get_order(order_id)
    if not order:
        return 0
    remaining = order.get("remaining_count", 0)
    # For cancelled orders, remaining_count reflects what was left at cancel time.
    # filled = original - remaining
    return total_requested - remaining


def execute_take(market, side, price, ev_pct, fee_cents, dry_run=False,
                 take_size=None):
    """Place a taker order (post_only=False), then cancel unfilled remainder.

    Returns number of contracts actually filled (0 if failed).
    """
    if take_size is None:
        take_size = config.TAKE_CONTRACT_SIZE
    ticker = market["ticker"]
    event_ticker = market.get("event_ticker", "")
    take_id = f"take-{uuid.uuid4().hex[:8]}"

    home_team = market.get("home_team")
    away_team = market.get("away_team")
    market_type = market.get("market_type", "spreads")

    if dry_run:
        print(f"  [TAKER] TAKE (dry): {side.upper()} {take_size}x @ {price}c on {ticker} "
              f"(fair={market['fair_prob']*100:.1f}c, EV={ev_pct:.1%}, fee={fee_cents:.1f}c) [{market_type}]")
        set_cooldown(ticker)
        return 0

    # CONFIRM: re-fetch book to make sure the price is still there
    live_bid, live_ask = _confirm_book_price(ticker, side)
    if side == "yes":
        if live_ask <= 0 or live_ask > price:
            # Ask disappeared or moved up — opportunity gone
            return 0
        # Use the live price (might be even better)
        price = live_ask
    else:
        live_no_price = 100 - live_bid if live_bid > 0 else 0
        if live_no_price <= 0 or live_no_price > price:
            return 0
        price = live_no_price

    # Re-check EV with confirmed price
    fair_cents = market["fair_prob"] * 100
    ev_pct, fee_cents = compute_take_ev(fair_cents, price, side)
    if ev_pct is None:
        return 0  # No longer +EV at confirmed price

    print(f"  [TAKER] TAKE: {side.upper()} {take_size}x @ {price}c on {ticker} "
          f"(fair={fair_cents:.1f}c, EV={ev_pct:.1%}, fee={fee_cents:.1f}c)")

    # Place order with post_only=False to cross the spread
    result = orders.place_order(
        ticker, side, price, take_size, post_only=False
    )

    if not result:
        return 0

    order_id = result.get("order_id", "")
    remaining = result.get("remaining_count", 0)

    # Cancel unfilled remainder immediately — we're a taker, not a maker.
    # Retry up to 3 times: a failed cancel leaves a toxic resting order
    # at a price chosen to cross the spread (easy pickings for sharps).
    if remaining > 0:
        cancelled = False
        for attempt in range(3):
            if orders.cancel_order(order_id):
                cancelled = True
                break
            time.sleep(0.1 * (attempt + 1))
        if not cancelled:
            print(f"    WARNING: Failed to cancel remainder for {order_id} after 3 attempts!")

    # RECONCILE: re-check the order via API to get actual final fill count.
    # Fills can happen between place_order response and cancel_order.
    # The place_order response's remaining_count may be stale.
    filled = _reconcile_fill_count(order_id, take_size)

    if filled > 0:
        print(f"    Filled {filled}/{take_size} contracts")
    else:
        print(f"    No fills (order cancelled or rejected)")

    # Record fill and update position
    if filled > 0:
        # Compute line_value in home-spread convention for position tracking
        strike = market.get("strike")
        if market_type == "spreads" and strike is not None:
            fill_line = -strike if market.get("is_home_contract") else strike
        elif market_type == "totals" and strike is not None:
            fill_line = strike
        else:
            fill_line = None

        db.update_position(
            ticker, side, price, filled,
            event_ticker=event_ticker,
            home_team=market.get("home_team"),
            away_team=market.get("away_team"),
            market_type=market_type,
            fair_prob=market.get("fair_prob"),
            line_value=fill_line,
            contract_team=market.get("contract_team"),
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
            event_ticker=event_ticker,
            side=side, price=price,
            fair_cents=fair_cents,
            ev_pct=ev_pct, fee_cents=fee_cents,
            count_requested=take_size,
            count_filled=filled, order_id=order_id,
        )

    # Tactical cooldown — prevent hammering the same contract
    set_cooldown(ticker)
    return filled


def run_take_cycle(quotable_markets, prediction_updated_at, dry_run=False,
                   resting_by_ticker=None):
    """Scan all markets for +EV takes and execute.

    Args:
        quotable_markets: Markets with fair values
        prediction_updated_at: Timestamp of current predictions
        dry_run: If True, scan only
        resting_by_ticker: MM's resting order state — used to avoid doubling up
    """
    if resting_by_ticker is None:
        resting_by_ticker = {}

    # Clear per-cycle caches — pass quotable_markets for Kalshi position mapping
    kelly.clear_positions_cache(current_markets=quotable_markets)
    risk.clear_exposure_cache()

    # SAFETY: don't take if Kalshi positions API failed
    if not kelly.positions_api_healthy():
        return

    # Staleness check (safety — don't trade on stale predictions)
    is_fresh, pred_age = risk.check_staleness(prediction_updated_at)
    if not is_fresh:
        return  # Silent — don't spam logs every 2s

    takes_this_cycle = 0

    for market in quotable_markets:
        ticker = market["ticker"]
        event_ticker = market.get("event_ticker", "")
        home_team = market.get("home_team")
        away_team = market.get("away_team")
        market_type = market.get("market_type", "spreads")
        fair_cents = market["fair_prob"] * 100

        # Skip if recently taken (10s tactical cooldown)
        if is_on_cooldown(ticker):
            continue

        # Skip tipoff proximity
        if not risk.check_tipoff_proximity(market.get("commence_time")):
            continue

        # Hard exposure cap — compute remaining room
        allowed, cur_exp, max_exp, _ = risk.check_game_type_exposure(
            home_team or "", away_team or "", market_type)
        exposure_room = max_exp - cur_exp if allowed else 0
        if exposure_room <= 0:
            continue

        game_key = market.get("game_key", (home_team or "", away_team or ""))

        # Check what resting orders the MM has on this ticker
        mm_resting = resting_by_ticker.get(ticker, {})

        yes_ask = market.get("book_ask", 0)
        yes_bid = market.get("book_bid", 0)

        # Compute Kelly sizes for this take (separate for YES and NO sides)
        # Uses actual execution price + taker fee (not bid/ask for resting orders)
        if config.USE_KELLY_SIZING:
            game_id = market.get("game_id")
            placed = []
            if game_id and game_key[0]:
                placed = kelly.get_placed_positions_for_game(game_key, game_id, current_markets=quotable_markets)

            if yes_ask > 0:
                fee = compute_taker_fee(yes_ask)
                yes_take_size = kelly.kelly_size_for_take(
                    market, "yes", yes_ask, fee,
                    placed_positions=placed,
                    bankroll=config.BANKROLL, kelly_mult=config.KELLY_FRACTION
                )
            else:
                yes_take_size = 0

            if yes_bid > 0:
                no_price = 100 - yes_bid
                fee = compute_taker_fee(no_price)
                no_take_size = kelly.kelly_size_for_take(
                    market, "no", no_price, fee,
                    placed_positions=placed,
                    bankroll=config.BANKROLL, kelly_mult=config.KELLY_FRACTION
                )
            else:
                no_take_size = 0
        else:
            yes_take_size = config.TAKE_CONTRACT_SIZE
            no_take_size = config.TAKE_CONTRACT_SIZE

        # Clamp take sizes to remaining exposure room
        if yes_ask > 0 and yes_take_size > 0:
            max_contracts = int(exposure_room / (yes_ask / 100))
            yes_take_size = min(yes_take_size, max_contracts)
        if yes_bid > 0 and no_take_size > 0:
            no_price_for_cap = 100 - yes_bid
            max_contracts = int(exposure_room / (no_price_for_cap / 100)) if no_price_for_cap > 0 else 0
            no_take_size = min(no_take_size, max_contracts)

        # --- Evaluate YES take (buy YES at ask) ---
        if yes_ask > 0 and yes_take_size > 0:
            # Don't take YES if MM already has a resting YES bid (would double up)
            if not mm_resting.get("bid_order_id"):
                ev_pct, fee = compute_take_ev(fair_cents, yes_ask, "yes")
                if ev_pct is not None:
                    filled = execute_take(market, "yes", yes_ask, ev_pct, fee,
                                          dry_run=dry_run,
                                          take_size=yes_take_size)
                    if filled > 0:
                        takes_this_cycle += 1
                        kelly.clear_positions_cache()  # Invalidate after fill
                        continue  # Don't also check NO for same ticker

        # --- Evaluate NO take (buy NO at 100-bid) ---
        if yes_bid > 0 and no_take_size > 0:
            # Don't take NO if MM already has a resting NO ask (would double up)
            if not mm_resting.get("ask_order_id"):
                no_price = 100 - yes_bid
                ev_pct, fee = compute_take_ev(fair_cents, no_price, "no")
                if ev_pct is not None:
                    filled = execute_take(market, "no", no_price, ev_pct, fee,
                                          dry_run=dry_run,
                                          take_size=no_take_size)
                    if filled > 0:
                        takes_this_cycle += 1
                        kelly.clear_positions_cache()  # Invalidate after fill

    if takes_this_cycle > 0:
        print(f"  [TAKER] Executed {takes_this_cycle} takes this cycle")


# --- Standalone mode (only when MM is not running) ---

def main():
    """Run taker as standalone process. Only use when MM is NOT running."""
    from scraper import fetch_markets
    from canonical_match import load_team_dict, load_canonical_games
    from main import (match_kalshi_to_predictions, enforce_monotonicity,
                      refresh_book_data, fetch_all_markets, match_all_markets)

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
    size_str = "Kelly-sized" if config.USE_KELLY_SIZING else f"{config.TAKE_CONTRACT_SIZE} contracts"
    print(f"  Take Size: {size_str}")
    print(f"  Poll: {config.TAKE_POLL_SEC}s | Cooldown: 10s per ticker")
    print(f"  Risk: Kelly-only (no hard position limits)")
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

    enabled = ", ".join(sorted(config.ENABLED_MARKET_TYPES))
    print(f"Fetching Kalshi 1H markets (enabled: {enabled})...")
    all_kalshi = fetch_all_markets()
    total_markets = sum(len(v) for v in all_kalshi.values())
    print(f"  Found {total_markets} open markets")

    quotable = match_all_markets(all_kalshi, predictions, team_dict, canonical_games)
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
                    fresh_all = fetch_all_markets()
                    if fresh_all:
                        all_kalshi = fresh_all
                    quotable = match_all_markets(
                        all_kalshi, predictions, team_dict, canonical_games
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
