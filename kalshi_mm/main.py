#!/usr/bin/env python3
"""
Kalshi CBB 1H Market Maker — Low-Stakes Prototype

Posts resting orders on Kalshi 1H spread markets using the answer key's
sample-based predictions as fair value. Monitors Bookmaker + Bet105 for
line moves and pulls stale quotes.

Usage:
    python main.py              # Run market maker
    python main.py --dry-run    # Compute quotes but don't place orders
"""

import sys
import time
import signal
import uuid
from datetime import datetime, timezone
from pathlib import Path

# Add kalshi_mm to path for imports
sys.path.insert(0, str(Path(__file__).parent))

import config
import db
import orders
import quoter
import risk

# Add kalshi_odds to path for market fetching
sys.path.insert(0, str(config.PROJECT_ROOT / "kalshi_odds"))
from scraper import fetch_markets, parse_spread_team, parse_matchup_title

# Add Answer Keys for team name resolution
sys.path.insert(0, str(config.ANSWER_KEYS_DIR))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names


# Global state
SESSION_ID = str(uuid.uuid4())[:8]
RUNNING = True
DRY_RUN = "--dry-run" in sys.argv
TOTAL_FILLS = 0
TOTAL_FEES = 0.0


def signal_handler(sig, frame):
    """Graceful shutdown on SIGINT/SIGTERM."""
    global RUNNING
    print("\n\nShutdown signal received. Cancelling all orders...")
    RUNNING = False


signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)


def match_kalshi_to_predictions(kalshi_markets, predictions, team_dict, canonical_games):
    """Match Kalshi spread markets to answer key predictions.

    Returns list of quotable markets with their fair values.
    """
    from collections import defaultdict
    from scraper import _fuzzy_team_match

    quotable = []

    # Group Kalshi markets by event
    by_event = defaultdict(list)
    for m in kalshi_markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Extract team names from contracts
        team_names = set()
        for m in event_markets:
            team = parse_spread_team(m.get("title", ""))
            if team:
                team_names.add(team)

        if len(team_names) != 2:
            continue

        team_list = sorted(team_names)

        # Resolve to canonical names
        away_resolved, home_resolved = resolve_team_names(
            team_list[0], team_list[1], team_dict, canonical_games
        )

        # Find matching predictions
        matching_preds = [
            p for p in predictions
            if p["market"] == "spreads_h1"
            and _fuzzy_team_match(p["home_team"].lower(), home_resolved.lower())
            and _fuzzy_team_match(p["away_team"].lower(), away_resolved.lower())
        ]

        if not matching_preds:
            continue

        # For each contract in this event, check if we have a prediction at that strike
        for m in event_markets:
            strike = m.get("floor_strike")
            if strike is None:
                continue

            team = parse_spread_team(m.get("title", ""))
            if not team:
                continue

            # Determine if this contract's team is home or away
            is_home = _fuzzy_team_match(team.lower(), home_resolved.lower())
            is_away = _fuzzy_team_match(team.lower(), away_resolved.lower())
            if not is_home and not is_away:
                continue

            # "Team wins by >X" → team is -X. Match to prediction at that spread.
            # If contract team is home: home_spread = -strike
            # The prediction has book_home_spread, so we look for -strike
            target_spread = -strike if is_home else strike

            pred = None
            for p in matching_preds:
                if abs(p["line_value"] - target_spread) < 0.1:
                    pred = p
                    break

            if pred is None:
                continue

            # Fair probability for YES on this contract ("team wins by >strike")
            # pred has prob_side1 = home_cover_prob at line_value = home_spread
            if is_home:
                # Home team -strike: home_spread = -strike, prob_side1 = P(home covers)
                fair_prob = pred["prob_side1"]
            else:
                # Away team -strike: home_spread = +strike (away is favored)
                # P(away wins by >strike) = P(home margin < -strike) = 1 - P(home covers +strike)
                fair_prob = 1 - pred["prob_side1"]

            quotable.append({
                "ticker": m["ticker"],
                "event_ticker": event_ticker,
                "home_team": home_resolved,
                "away_team": away_resolved,
                "contract_team": team,
                "strike": strike,
                "is_home_contract": is_home,
                "fair_prob": fair_prob,
                "commence_time": pred.get("commence_time"),
            })

    return quotable


def run_quote_cycle(quotable_markets, resting_by_ticker, prediction_updated_at):
    """Compute quotes and place/amend/cancel orders for all quotable markets.

    Args:
        quotable_markets: List of markets with fair values
        resting_by_ticker: Dict of ticker -> resting order info
        prediction_updated_at: Timestamp of when predictions were last updated

    Returns:
        Updated resting_by_ticker dict.
    """
    global TOTAL_FILLS

    # Check overall risk
    is_fresh, pred_age = risk.check_staleness(prediction_updated_at)
    if not is_fresh:
        print(f"  STALE PREDICTIONS ({pred_age:.0f}s old). Pulling all quotes.")
        if not DRY_RUN:
            orders.cancel_all_orders()
            resting_by_ticker.clear()
        return resting_by_ticker

    exposure_ok, exposure = risk.check_exposure_limit()
    if not exposure_ok:
        print(f"  EXPOSURE LIMIT (${exposure:.2f}). Pulling all quotes.")
        if not DRY_RUN:
            orders.cancel_all_orders()
            resting_by_ticker.clear()
        return resting_by_ticker

    quoted_count = 0
    for market in quotable_markets:
        ticker = market["ticker"]

        # Check tipoff proximity
        if not risk.check_tipoff_proximity(market.get("commence_time")):
            # Cancel if we have resting orders for this ticker
            if ticker in resting_by_ticker:
                if not DRY_RUN:
                    for side_key in ["bid_order_id", "ask_order_id"]:
                        oid = resting_by_ticker[ticker].get(side_key)
                        if oid:
                            orders.cancel_order(oid)
                    db.remove_resting_order(resting_by_ticker[ticker].get("bid_order_id", ""))
                    db.remove_resting_order(resting_by_ticker[ticker].get("ask_order_id", ""))
                del resting_by_ticker[ticker]
            continue

        # Check market count limit
        if quoted_count >= config.MAX_MARKETS and ticker not in resting_by_ticker:
            continue

        # Get current position
        pos = db.get_position(ticker)
        net = pos["net_yes"]

        # Check position limit — skip bid side if at max long, skip ask side if at max short
        at_max_long = net >= config.MAX_POSITION_PER_MARKET
        at_max_short = net <= -config.MAX_POSITION_PER_MARKET

        # Compute quote
        quote = quoter.compute_quotes(market["fair_prob"], net)
        if quote is None:
            # Cancel existing orders if any
            if ticker in resting_by_ticker:
                if not DRY_RUN:
                    for side_key in ["bid_order_id", "ask_order_id"]:
                        oid = resting_by_ticker[ticker].get(side_key)
                        if oid:
                            orders.cancel_order(oid)
                del resting_by_ticker[ticker]
            continue

        # Log the quote
        db.log_quote(ticker, market["fair_prob"], quote["bid_yes"], quote["ask_yes"],
                     net, quote["skew"], pred_age)

        if DRY_RUN:
            print(quoter.format_quote_summary(ticker, quote))
            quoted_count += 1
            continue

        # Manage bid (buy YES)
        existing = resting_by_ticker.get(ticker, {})

        if not at_max_long:
            bid_oid = existing.get("bid_order_id")
            if bid_oid:
                # Amend if price changed
                if quoter.should_amend(existing.get("bid_price", 0), quote["bid_yes"]):
                    result = orders.amend_order(bid_oid, price=quote["bid_yes"], count=quote["size"])
                    if result:
                        existing["bid_price"] = quote["bid_yes"]
                        db.save_resting_order(bid_oid, ticker, "yes", "buy",
                                              quote["bid_yes"], quote["size"])
            else:
                # Place new bid
                result = orders.place_order(ticker, "yes", quote["bid_yes"], quote["size"])
                if result:
                    bid_oid = result.get("order_id")
                    existing["bid_order_id"] = bid_oid
                    existing["bid_price"] = quote["bid_yes"]
                    db.save_resting_order(bid_oid, ticker, "yes", "buy",
                                          quote["bid_yes"], quote["size"])
        elif existing.get("bid_order_id"):
            orders.cancel_order(existing["bid_order_id"])
            db.remove_resting_order(existing["bid_order_id"])
            existing.pop("bid_order_id", None)

        # Manage ask (buy NO, which is selling YES)
        if not at_max_short:
            ask_oid = existing.get("ask_order_id")
            no_price = 100 - quote["ask_yes"]  # Buying NO at this price = selling YES at ask_yes
            if ask_oid:
                if quoter.should_amend(existing.get("ask_price", 0), quote["ask_yes"]):
                    result = orders.amend_order(ask_oid, price=no_price, count=quote["size"])
                    if result:
                        existing["ask_price"] = quote["ask_yes"]
                        db.save_resting_order(ask_oid, ticker, "no", "buy",
                                              no_price, quote["size"])
            else:
                result = orders.place_order(ticker, "no", no_price, quote["size"])
                if result:
                    ask_oid = result.get("order_id")
                    existing["ask_order_id"] = ask_oid
                    existing["ask_price"] = quote["ask_yes"]
                    db.save_resting_order(ask_oid, ticker, "no", "buy",
                                          no_price, quote["size"])
        elif existing.get("ask_order_id"):
            orders.cancel_order(existing["ask_order_id"])
            db.remove_resting_order(existing["ask_order_id"])
            existing.pop("ask_order_id", None)

        resting_by_ticker[ticker] = existing
        quoted_count += 1

    print(f"  Quoting {quoted_count} markets (exposure: ${exposure:.2f})")
    return resting_by_ticker


def poll_for_fills(resting_by_ticker):
    """Check resting orders for fills and update positions."""
    global TOTAL_FILLS, TOTAL_FEES

    if DRY_RUN:
        return

    api_orders = orders.get_resting_orders()
    api_order_map = {o["order_id"]: o for o in api_orders}

    for ticker, info in list(resting_by_ticker.items()):
        for side_key, side in [("bid_order_id", "yes"), ("ask_order_id", "no")]:
            oid = info.get(side_key)
            if not oid:
                continue

            api_order = api_order_map.get(oid)

            if api_order is None:
                # Order no longer resting — filled, cancelled, or expired.
                # Clean up our local state so we don't leak stale entries.
                print(f"  Order {oid} no longer resting (filled/cancelled externally)")
                db.remove_resting_order(oid)
                info.pop(side_key, None)
                price_key = "bid_price" if side == "yes" else "ask_price"
                info.pop(price_key, None)
                continue

            remaining = api_order.get("remaining_count", api_order.get("count", 0))
            original = api_order.get("count", 0)
            filled_count = original - remaining

            if filled_count > 0:
                # Partial or full fill
                price = api_order.get("yes_price", 0) if side == "yes" else api_order.get("no_price", 0)

                print(f"  FILL: {side} {filled_count}x @ {price}c on {ticker}"
                      f"{' (partial)' if remaining > 0 else ''}")
                TOTAL_FILLS += filled_count

                # Update position
                db.update_position(ticker, side, price, filled_count)
                db.record_fill(
                    fill_id=f"{oid}-fill-{filled_count}",
                    ticker=ticker, side=side, action="buy",
                    price=price, count=filled_count, fee_cents=0,
                    order_id=oid
                )

                if remaining == 0:
                    # Fully filled — remove
                    db.remove_resting_order(oid)
                    info.pop(side_key, None)
                    price_key = "bid_price" if side == "yes" else "ask_price"
                    info.pop(price_key, None)

    # Clean up empty ticker entries to prevent state drift
    for ticker in [t for t, info in resting_by_ticker.items()
                   if not info.get("bid_order_id") and not info.get("ask_order_id")]:
        del resting_by_ticker[ticker]


def main():
    global RUNNING

    print("=" * 60)
    print(f"  Kalshi CBB 1H Market Maker — Session {SESSION_ID}")
    print(f"  {'DRY RUN' if DRY_RUN else 'LIVE'}")
    print(f"  Half-spread: {config.HALF_SPREAD_CENTS}c | Size: {config.CONTRACT_SIZE}")
    print(f"  Max position: {config.MAX_POSITION_PER_MARKET} | Max exposure: ${config.MAX_TOTAL_EXPOSURE_DOLLARS}")
    print(f"  API: {config.KALSHI_BASE_URL}")
    print("=" * 60)

    # Init
    db.init_database()
    db.start_session(SESSION_ID)

    # Load predictions
    print("\nLoading predictions from answer key...")
    predictions, prediction_updated_at = db.load_predictions()
    if not predictions:
        print("No predictions available. Run the CBB pipeline first:")
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

    # Load team resolution
    print("Loading team name resolution...")
    team_dict = load_team_dict("cbb")
    canonical_games = load_canonical_games("cbb")

    # Fetch Kalshi 1H spread markets
    print(f"Fetching Kalshi 1H spread markets ({config.SPREAD_SERIES})...")
    kalshi_markets = fetch_markets(config.SPREAD_SERIES)
    print(f"  Found {len(kalshi_markets)} open markets")

    if not kalshi_markets:
        print("No open Kalshi 1H spread markets found.")
        return

    # Match to predictions
    print("Matching Kalshi markets to predictions...")
    quotable = match_kalshi_to_predictions(kalshi_markets, predictions, team_dict, canonical_games)
    print(f"  {len(quotable)} quotable markets found")

    if not quotable:
        print("No markets matched predictions. Check team name resolution.")
        return

    # Snapshot reference lines for monitoring
    print("Snapshotting reference lines from Bookmaker/Bet105...")
    ref_lines = risk.run_line_monitor()
    if ref_lines:
        db.save_reference_lines(ref_lines)
        print(f"  Saved {len(ref_lines)} reference lines")
    else:
        print("  No reference lines available (scrapers may have no 1H data)")

    # Main loop
    print(f"\nStarting main loop (quote every {config.QUOTE_CYCLE_SEC}s, "
          f"monitor every {config.MONITOR_CYCLE_SEC}s)...\n")

    resting_by_ticker = {}
    last_quote_time = 0
    last_fill_poll = 0
    last_monitor_time = 0
    last_refresh_time = time.time()

    try:
        while RUNNING:
            now = time.time()

            # Quote cycle
            if now - last_quote_time >= config.QUOTE_CYCLE_SEC:
                print(f"\n--- Quote cycle @ {datetime.now().strftime('%H:%M:%S')} ---")
                resting_by_ticker = run_quote_cycle(quotable, resting_by_ticker, prediction_updated_at)
                last_quote_time = now

            # Fill polling
            if now - last_fill_poll >= config.FILL_POLL_SEC:
                poll_for_fills(resting_by_ticker)
                last_fill_poll = now

            # Line-move monitoring
            if now - last_monitor_time >= config.MONITOR_CYCLE_SEC:
                print(f"\n--- Line monitor @ {datetime.now().strftime('%H:%M:%S')} ---")
                current_lines = risk.run_line_monitor()
                ref_lines = db.get_reference_lines()

                if current_lines and ref_lines:
                    moved = risk.detect_line_moves(current_lines, ref_lines)
                    if moved:
                        print(f"  {len(moved)} games with line moves. Pulling quotes.")
                        # Pull quotes for moved games
                        for ticker, info in list(resting_by_ticker.items()):
                            market = next((m for m in quotable if m["ticker"] == ticker), None)
                            if market and (market["home_team"], market["away_team"]) in moved:
                                if not DRY_RUN:
                                    for side_key in ["bid_order_id", "ask_order_id"]:
                                        oid = info.get(side_key)
                                        if oid:
                                            orders.cancel_order(oid)
                                del resting_by_ticker[ticker]
                last_monitor_time = now

            # Periodic prediction refresh (every 15 min)
            if now - last_refresh_time >= 900:
                print("\n--- Refreshing predictions ---")
                new_preds, new_ts = db.load_predictions()
                if new_preds and new_ts != prediction_updated_at:
                    predictions = new_preds
                    prediction_updated_at = new_ts
                    # Re-match
                    quotable = match_kalshi_to_predictions(
                        kalshi_markets, predictions, team_dict, canonical_games
                    )
                    print(f"  Refreshed: {len(quotable)} quotable markets")
                    # Update reference lines
                    ref_lines = risk.run_line_monitor()
                    if ref_lines:
                        db.save_reference_lines(ref_lines)
                last_refresh_time = now

            time.sleep(1)

    except Exception as e:
        print(f"\n  UNHANDLED ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Kill switch: cancel everything — retry up to 3 times
        print("\n\nShutting down...")
        if not DRY_RUN and resting_by_ticker:
            print("Cancelling all resting orders...")
            for attempt in range(3):
                if orders.cancel_all_orders():
                    break
                print(f"  Kill switch retry {attempt + 1}/3...")
                time.sleep(1)

        # Log session summary
        positions = db.get_all_positions()
        print(f"\n{'=' * 60}")
        print(f"  Session {SESSION_ID} complete")
        print(f"  Total fills: {TOTAL_FILLS}")
        print(f"  Open positions: {len(positions)}")
        for p in positions:
            print(f"    {p['ticker']}: net_yes={p['net_yes']}, avg={p['avg_entry_price']:.1f}c")
        print(f"{'=' * 60}")

        db.end_session(SESSION_ID, TOTAL_FILLS, 0.0, TOTAL_FEES,
                       len(resting_by_ticker))


if __name__ == "__main__":
    main()
