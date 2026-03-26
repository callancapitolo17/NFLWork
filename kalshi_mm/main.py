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

import duckdb
import random
import subprocess
import sys
import time
import signal
import uuid
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Add kalshi_mm to path for imports
sys.path.insert(0, str(Path(__file__).parent))

import config
import db
import kelly
import orders
import quoter
import risk
import taker

# Add kalshi_odds to path for market fetching
sys.path.insert(0, str(config.PROJECT_ROOT / "kalshi_odds"))
from scraper import (fetch_markets, parse_spread_team, parse_matchup_title,
                     parse_race_to_10_team, _fuzzy_team_match)

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

def _nuke_and_adopt(resting_by_ticker):
    """Cancel all bot orders and adopt any survivors to prevent duplicates.

    Used when predictions are stale, pipeline fails, etc. After cancelling,
    re-queries Kalshi to catch bot orders that survived (rate limit, partial failure).
    Manual (user-placed) orders are preserved and not adopted into bot tracking.
    """
    orders.cancel_all_bot_orders()
    db.clear_all_resting_orders()
    resting_by_ticker.clear()
    survivors = orders.get_resting_orders()
    # Only adopt bot orders — manual orders should not be managed by the bot
    bot_survivors = [o for o in survivors if orders.is_bot_order(o)]
    manual_count = len(survivors) - len(bot_survivors)
    if manual_count:
        print(f"  Skipping {manual_count} manual order(s) during adopt")
    if bot_survivors:
        print(f"  WARNING: {len(bot_survivors)} bot orders survived cancel — adopting")
        for o in bot_survivors:
            ticker = o.get("ticker", "")
            oid = o.get("order_id", "")
            side = o.get("side", "")
            price = int(round(float(o.get("yes_price_dollars", 0)) * 100))
            existing = resting_by_ticker.get(ticker, {})
            if side == "yes":
                existing["bid_order_id"] = oid
                existing["bid_price"] = price
            elif side == "no":
                existing["ask_order_id"] = oid
                existing["ask_price"] = 100 - price
            resting_by_ticker[ticker] = existing


# --- Background pipeline runner ---
_pipeline_proc = None
_pipeline_start_time = None
_pipeline_backoff = config.PIPELINE_REFRESH_SEC


def start_pipeline():
    """Launch prediction pipeline as background subprocess. Non-blocking."""
    global _pipeline_proc, _pipeline_start_time
    if _pipeline_proc and _pipeline_proc.poll() is None:
        return  # Already running
    print(f"\n--- Starting pipeline @ {datetime.now().strftime('%H:%M:%S')} ---")
    _pipeline_proc = subprocess.Popen(
        [sys.executable, "run.py", "cbb"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        cwd=str(config.ANSWER_KEYS_DIR)
    )
    _pipeline_start_time = time.time()


def check_pipeline_completion(resting_by_ticker):
    """Check if background pipeline finished. Returns True if new predictions ready."""
    global _pipeline_proc, _pipeline_start_time, _pipeline_backoff
    if _pipeline_proc is None:
        return False

    rc = _pipeline_proc.poll()
    if rc is None:
        # Still running — check timeout (120s)
        if time.time() - _pipeline_start_time > 120:
            _pipeline_proc.kill()
            _pipeline_proc = None
            print("  Pipeline timed out (120s) — killed. Pulling all quotes.")
            _pipeline_backoff = min(_pipeline_backoff * 2, 2400)
            if not DRY_RUN:
                _nuke_and_adopt(resting_by_ticker)
            return False
        return False

    elapsed = time.time() - _pipeline_start_time
    _pipeline_proc = None

    if rc == 0:
        print(f"  Pipeline complete ({elapsed:.0f}s)")
        _pipeline_backoff = config.PIPELINE_REFRESH_SEC
        return True
    else:
        print(f"  Pipeline FAILED (exit {rc}, {elapsed:.0f}s). Pulling all quotes. Retry in {_pipeline_backoff}s")
        _pipeline_backoff = min(_pipeline_backoff * 2, 2400)
        if not DRY_RUN:
            _nuke_and_adopt(resting_by_ticker)
        return False


def enforce_monotonicity(quotable_markets):
    """Clamp fair probabilities so strike ordering is consistent.

    Spreads: P(team wins by >N) >= P(team wins by >N+k) — non-increasing.
    Totals: P(over N) >= P(over N+k) — higher total = lower over prob.
    Moneylines: No strikes, skip.

    Without this, a sharp can arb mispriced strikes against each other.
    """
    from collections import defaultdict

    groups = defaultdict(list)
    for m in quotable_markets:
        mtype = m.get("market_type", "spreads")
        if mtype in ("moneyline", "race_to_10"):
            continue  # No strikes to enforce
        if mtype == "totals":
            # Totals: group by event_ticker alone (no team direction)
            key = ("totals", m["event_ticker"])
        else:
            # Spreads: group by (event_ticker, is_home_contract)
            key = ("spreads", m["event_ticker"], m["is_home_contract"])
        groups[key].append(m)

    for key, markets in groups.items():
        # Sort by strike ascending — probability should decrease
        markets.sort(key=lambda m: m["strike"])
        for i in range(1, len(markets)):
            if markets[i]["fair_prob"] > markets[i - 1]["fair_prob"]:
                markets[i]["fair_prob"] = markets[i - 1]["fair_prob"]


def refresh_book_data(quotable_markets):
    """Fetch fresh yes_bid/yes_ask from all enabled series and update in-place."""
    book = {}
    for mtype, series in config.MARKET_SERIES.items():
        if mtype not in config.ENABLED_MARKET_TYPES:
            continue
        fresh = fetch_markets(series)
        if fresh:
            for m in fresh:
                book[m["ticker"]] = (
                    int(round(float(m.get("yes_bid_dollars", 0)) * 100)),
                    int(round(float(m.get("yes_ask_dollars", 0)) * 100))
                )
    # Timestamp AFTER all API calls complete so staleness check is accurate
    now = time.time()
    for market in quotable_markets:
        bid_ask = book.get(market["ticker"])
        if bid_ask:
            market["book_bid"], market["book_ask"] = bid_ask
            market["book_fetched_at"] = now


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

        # Resolve to canonical names (resolve_team_names preserves input order,
        # so we don't know home/away yet — try both orderings against predictions)
        resolved_a, resolved_b = resolve_team_names(
            team_list[0], team_list[1], team_dict, canonical_games
        )

        # Try both orderings to find matching predictions
        matching_preds = [
            p for p in predictions
            if p["market"] == "spreads_h1"
            and _fuzzy_team_match(p["home_team"].lower(), resolved_b.lower())
            and _fuzzy_team_match(p["away_team"].lower(), resolved_a.lower())
        ]
        home_resolved, away_resolved = resolved_b, resolved_a

        if not matching_preds:
            # Try swapped ordering
            matching_preds = [
                p for p in predictions
                if p["market"] == "spreads_h1"
                and _fuzzy_team_match(p["home_team"].lower(), resolved_a.lower())
                and _fuzzy_team_match(p["away_team"].lower(), resolved_b.lower())
            ]
            home_resolved, away_resolved = resolved_a, resolved_b

        if not matching_preds:
            # Fallback: try raw Kalshi names directly against predictions
            # (handles cases where resolve_team_names fails but "Indiana" still
            # fuzzy-matches "Indiana Hoosiers" in predictions)
            for raw_a, raw_b in [(team_list[0], team_list[1]), (team_list[1], team_list[0])]:
                matching_preds = [
                    p for p in predictions
                    if p["market"] == "spreads_h1"
                    and _fuzzy_team_match(p["home_team"].lower(), raw_b.lower())
                    and _fuzzy_team_match(p["away_team"].lower(), raw_a.lower())
                ]
                if matching_preds:
                    home_resolved, away_resolved = raw_b, raw_a
                    break

        if not matching_preds:
            continue

        # Now home_resolved/away_resolved match the prediction's home/away

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
            if is_home and is_away:
                continue  # Ambiguous match (e.g. both contain "State") — skip
            if not is_home and not is_away:
                continue

            # "Team wins by >X" → team is -X. Match to prediction at that spread.
            # If contract team is home: home_spread = -strike
            # The prediction has book_home_spread, so we look for -strike
            target_spread = -strike if is_home else strike

            pred = None
            for p in matching_preds:
                if p["line_value"] is not None and abs(p["line_value"] - target_spread) < 0.1:
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
                "home_team": pred["home_team"],
                "away_team": pred["away_team"],
                "game_key": (pred["home_team"], pred["away_team"]),
                "game_id": pred.get("id"),
                "market_type": "spreads",
                "contract_team": team,
                "strike": strike,
                "is_home_contract": is_home,
                "fair_prob": fair_prob,
                "commence_time": pred.get("commence_time"),
                "book_bid": int(round(float(m.get("yes_bid_dollars", 0)) * 100)),
                "book_ask": int(round(float(m.get("yes_ask_dollars", 0)) * 100)),
            })

    return quotable


def fetch_all_markets():
    """Fetch markets from all enabled series. Returns {market_type: [markets]}."""
    result = {}
    for mtype, series in config.MARKET_SERIES.items():
        if mtype not in config.ENABLED_MARKET_TYPES:
            continue
        markets = fetch_markets(series)
        if markets:
            result[mtype] = markets
    return result


def match_total_markets(kalshi_markets, predictions, team_dict, canonical_games):
    """Match Kalshi 1H total markets to answer key predictions.

    Each total contract: YES = over, NO = under. Strike is the total line.
    """
    from collections import defaultdict

    quotable = []
    by_event = defaultdict(list)
    for m in kalshi_markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Parse teams from event title
        title = event_markets[0].get("title", "")
        away_raw, home_raw = parse_matchup_title(title)
        if not away_raw or not home_raw:
            continue

        away_resolved, home_resolved = resolve_team_names(
            away_raw, home_raw, team_dict, canonical_games
        )

        # Find matching totals predictions — try both orderings
        matching_preds = [
            p for p in predictions
            if p["market"] == "totals_h1"
            and _fuzzy_team_match(p["home_team"].lower(), home_resolved.lower())
            and _fuzzy_team_match(p["away_team"].lower(), away_resolved.lower())
        ]
        if not matching_preds:
            matching_preds = [
                p for p in predictions
                if p["market"] == "totals_h1"
                and _fuzzy_team_match(p["home_team"].lower(), away_resolved.lower())
                and _fuzzy_team_match(p["away_team"].lower(), home_resolved.lower())
            ]
            if matching_preds:
                home_resolved, away_resolved = away_resolved, home_resolved

        if not matching_preds:
            # Fallback: try raw Kalshi names
            for raw_h, raw_a in [(home_raw, away_raw), (away_raw, home_raw)]:
                matching_preds = [
                    p for p in predictions
                    if p["market"] == "totals_h1"
                    and _fuzzy_team_match(p["home_team"].lower(), raw_h.lower())
                    and _fuzzy_team_match(p["away_team"].lower(), raw_a.lower())
                ]
                if matching_preds:
                    home_resolved, away_resolved = raw_h, raw_a
                    break

        if not matching_preds:
            continue

        for m in event_markets:
            strike = m.get("floor_strike")
            if strike is None:
                continue

            # Find prediction at this total line
            pred = None
            for p in matching_preds:
                if p["line_value"] is not None and abs(p["line_value"] - strike) < 0.1:
                    pred = p
                    break
            if pred is None:
                continue

            # Use prediction's team names for game_key — guarantees consistency
            # across matchers even when raw/resolved names differ
            pred_home = pred["home_team"]
            pred_away = pred["away_team"]

            # fair_prob for YES = over probability = prob_side1
            quotable.append({
                "ticker": m["ticker"],
                "event_ticker": event_ticker,
                "home_team": pred_home,
                "away_team": pred_away,
                "game_key": (pred_home, pred_away),
                "game_id": pred.get("id"),
                "market_type": "totals",
                "contract_team": None,
                "strike": strike,
                "is_home_contract": None,
                "fair_prob": pred["prob_side1"],
                "commence_time": pred.get("commence_time"),
                "book_bid": int(round(float(m.get("yes_bid_dollars", 0)) * 100)),
                "book_ask": int(round(float(m.get("yes_ask_dollars", 0)) * 100)),
                "book_fetched_at": time.time(),
            })

    return quotable


def match_moneyline_markets(kalshi_markets, predictions, team_dict, canonical_games):
    """Match Kalshi 1H winner (3-way) markets to answer key predictions.

    Each event has 3 contracts: home, away, tie. Each is an independent binary market.
    Uses 3-way probabilities from the answer key (home/away/tie sum to 1.0).
    """
    from collections import defaultdict

    quotable = []
    by_event = defaultdict(list)
    for m in kalshi_markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        # Separate tie from team contracts
        team_contracts = []
        tie_contract = None
        for m in event_markets:
            sub = (m.get("yes_sub_title") or "").lower()
            if "tie" in sub:
                tie_contract = m
            else:
                team_contracts.append(m)

        if len(team_contracts) != 2:
            continue

        # Parse teams from event title
        title = event_markets[0].get("title", "")
        away_raw, home_raw = parse_matchup_title(title)
        if not away_raw or not home_raw:
            continue

        away_resolved, home_resolved = resolve_team_names(
            away_raw, home_raw, team_dict, canonical_games
        )

        # Find matching h2h_h1 prediction (one per game, no line matching)
        matching_preds = [
            p for p in predictions
            if p["market"] == "h2h_h1"
            and _fuzzy_team_match(p["home_team"].lower(), home_resolved.lower())
            and _fuzzy_team_match(p["away_team"].lower(), away_resolved.lower())
        ]
        if not matching_preds:
            matching_preds = [
                p for p in predictions
                if p["market"] == "h2h_h1"
                and _fuzzy_team_match(p["home_team"].lower(), away_resolved.lower())
                and _fuzzy_team_match(p["away_team"].lower(), home_resolved.lower())
            ]
            if matching_preds:
                home_resolved, away_resolved = away_resolved, home_resolved

        if not matching_preds:
            # Fallback: raw names
            for raw_h, raw_a in [(home_raw, away_raw), (away_raw, home_raw)]:
                matching_preds = [
                    p for p in predictions
                    if p["market"] == "h2h_h1"
                    and _fuzzy_team_match(p["home_team"].lower(), raw_h.lower())
                    and _fuzzy_team_match(p["away_team"].lower(), raw_a.lower())
                ]
                if matching_preds:
                    home_resolved, away_resolved = raw_h, raw_a
                    break

        if not matching_preds:
            continue

        pred = matching_preds[0]
        home_prob = pred["prob_side1"]

        # Use prediction's team names for game_key — guarantees consistency
        # across matchers even when raw/resolved names differ
        pred_home = pred["home_team"]
        pred_away = pred["away_team"]

        # Require 3-way probabilities from the answer key.
        # 2-way probs (ties excluded) systematically overvalue both team contracts
        # by the tie probability (~5-10% for CBB halves), giving sharps free edge.
        if pred.get("prob_side2") is None or pred.get("prob_tie") is None:
            continue  # Refuse to quote — pipeline must export 3-way probs

        away_prob = pred["prob_side2"]
        tie_prob = pred["prob_tie"]

        # Match team contracts to home/away
        home_contract = None
        away_contract = None
        for m in team_contracts:
            sub = (m.get("yes_sub_title") or "").strip()
            sub_lower = sub.lower()
            if _fuzzy_team_match(sub_lower, home_raw):
                home_contract = m
            elif _fuzzy_team_match(sub_lower, away_raw):
                away_contract = m

        if not home_contract or not away_contract:
            # Fallback: assign by title order
            c0_sub = (team_contracts[0].get("yes_sub_title") or "").strip().lower()
            if _fuzzy_team_match(c0_sub, away_raw):
                away_contract = team_contracts[0]
                home_contract = team_contracts[1]
            else:
                home_contract = team_contracts[0]
                away_contract = team_contracts[1]

        if not home_contract or not away_contract:
            continue

        now_ts = time.time()

        # Home YES contract
        quotable.append({
            "ticker": home_contract["ticker"],
            "event_ticker": event_ticker,
            "home_team": pred_home,
            "away_team": pred_away,
            "game_key": (pred_home, pred_away),
            "game_id": pred.get("id"),
            "market_type": "moneyline",
            "contract_team": pred_home,
            "strike": None,
            "is_home_contract": True,
            "fair_prob": home_prob,
            "commence_time": pred.get("commence_time"),
            "book_bid": int(round(float(home_contract.get("yes_bid_dollars", 0)) * 100)),
            "book_ask": int(round(float(home_contract.get("yes_ask_dollars", 0)) * 100)),
            "book_fetched_at": now_ts,
        })

        # Away YES contract
        quotable.append({
            "ticker": away_contract["ticker"],
            "event_ticker": event_ticker,
            "home_team": pred_home,
            "away_team": pred_away,
            "game_key": (pred_home, pred_away),
            "game_id": pred.get("id"),
            "market_type": "moneyline",
            "contract_team": pred_away,
            "strike": None,
            "is_home_contract": False,
            "fair_prob": away_prob,
            "commence_time": pred.get("commence_time"),
            "book_bid": int(round(float(away_contract.get("yes_bid_dollars", 0)) * 100)),
            "book_ask": int(round(float(away_contract.get("yes_ask_dollars", 0)) * 100)),
            "book_fetched_at": now_ts,
        })

        # Tie contract (if exists and we have tie probability)
        if tie_contract and tie_prob is not None and tie_prob > 0:
            quotable.append({
                "ticker": tie_contract["ticker"],
                "event_ticker": event_ticker,
                "home_team": pred_home,
                "away_team": pred_away,
                "game_key": (pred_home, pred_away),
                "game_id": pred.get("id"),
                "market_type": "moneyline",
                "contract_team": "Tie",
                "strike": None,
                "is_home_contract": None,
                "fair_prob": tie_prob,
                "commence_time": pred.get("commence_time"),
                "book_bid": int(round(float(tie_contract.get("yes_bid_dollars", 0)) * 100)),
                "book_ask": int(round(float(tie_contract.get("yes_ask_dollars", 0)) * 100)),
                "book_fetched_at": now_ts,
            })

    return quotable


def match_race_to_10_markets(kalshi_markets, predictions, team_dict, canonical_games):
    """Match Kalshi race-to-10 markets to answer key predictions.

    Each event has 2 contracts: one per team. YES = that team reaches 10 first.
    Fair value comes from race_to_10 prediction rows.
    """
    from collections import defaultdict

    quotable = []
    by_event = defaultdict(list)
    for m in kalshi_markets:
        by_event[m["event_ticker"]].append(m)

    for event_ticker, event_markets in by_event.items():
        if len(event_markets) != 2:
            continue

        # Extract team names from contract titles
        team_names = []
        for m in event_markets:
            team = parse_race_to_10_team(m.get("title", ""))
            if team:
                team_names.append(team)
        if len(team_names) != 2:
            continue

        # Resolve to canonical names — try both orderings
        resolved_a, resolved_b = resolve_team_names(
            team_names[0], team_names[1], team_dict, canonical_games
        )

        # Find matching race_to_10 prediction
        matching_preds = None
        pred_home = pred_away = None
        for try_home, try_away in [(resolved_b, resolved_a), (resolved_a, resolved_b)]:
            matches = [
                p for p in predictions
                if p["market"] == "race_to_10"
                and _fuzzy_team_match(p["home_team"].lower(), try_home.lower())
                and _fuzzy_team_match(p["away_team"].lower(), try_away.lower())
            ]
            if matches:
                matching_preds = matches
                pred_home, pred_away = try_home, try_away
                break

        # Fallback: raw names
        if not matching_preds:
            for raw_h, raw_a in [(team_names[0], team_names[1]),
                                 (team_names[1], team_names[0])]:
                matches = [
                    p for p in predictions
                    if p["market"] == "race_to_10"
                    and _fuzzy_team_match(p["home_team"].lower(), raw_h.lower())
                    and _fuzzy_team_match(p["away_team"].lower(), raw_a.lower())
                ]
                if matches:
                    matching_preds = matches
                    pred_home, pred_away = raw_h, raw_a
                    break

        if not matching_preds:
            continue

        pred = matching_preds[0]
        home_first10_prob = pred["prob_side1"]  # prob home reaches 10 first

        now_ts = time.time()
        for m in event_markets:
            contract_team = parse_race_to_10_team(m.get("title", ""))
            if not contract_team:
                continue

            # Determine if this contract is for home or away team
            is_home = _fuzzy_team_match(contract_team.lower(), pred_home.lower())
            fair_prob = home_first10_prob if is_home else (1 - home_first10_prob)

            quotable.append({
                "ticker": m["ticker"],
                "event_ticker": event_ticker,
                "home_team": pred["home_team"],
                "away_team": pred["away_team"],
                "game_key": (pred["home_team"], pred["away_team"]),
                "game_id": pred.get("id"),
                "market_type": "race_to_10",
                "contract_team": pred["home_team"] if is_home else pred["away_team"],
                "strike": None,
                "is_home_contract": is_home,
                "fair_prob": fair_prob,
                "commence_time": pred.get("commence_time"),
                "book_bid": int(round(float(m.get("yes_bid_dollars", 0)) * 100)),
                "book_ask": int(round(float(m.get("yes_ask_dollars", 0)) * 100)),
                "book_fetched_at": now_ts,
            })

    return quotable


def match_all_markets(all_kalshi, predictions, team_dict, canonical_games):
    """Match all enabled market types to predictions. Returns unified quotable list."""
    quotable = []

    if "spreads" in all_kalshi:
        spread_q = match_kalshi_to_predictions(
            all_kalshi["spreads"], predictions, team_dict, canonical_games
        )
        for m in spread_q:
            m.setdefault("book_fetched_at", time.time())
        quotable.extend(spread_q)

    if "totals" in all_kalshi:
        quotable.extend(match_total_markets(
            all_kalshi["totals"], predictions, team_dict, canonical_games
        ))

    if "moneyline" in all_kalshi:
        quotable.extend(match_moneyline_markets(
            all_kalshi["moneyline"], predictions, team_dict, canonical_games
        ))

    if "race_to_10" in all_kalshi:
        quotable.extend(match_race_to_10_markets(
            all_kalshi["race_to_10"], predictions, team_dict, canonical_games
        ))

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
    # Clear per-cycle caches — pass quotable_markets so Kalshi positions
    # can be mapped to game metadata (home_team, market_type, etc.)
    kelly.clear_positions_cache(current_markets=quotable_markets)
    risk.clear_exposure_cache()

    # SAFETY: verify Kalshi positions loaded successfully before sizing.
    # If API fails, refuse to quote — sizing with no position data causes
    # unbounded accumulation (the root cause of the Tennessee disaster).
    if not kelly.positions_api_healthy():
        print("  WARNING: Kalshi positions API failed — skipping quote cycle")
        return resting_by_ticker

    # Check overall risk
    is_fresh, pred_age = risk.check_staleness(prediction_updated_at)
    if not is_fresh:
        print(f"  STALE PREDICTIONS ({pred_age:.0f}s old). Pulling all quotes.")
        if not DRY_RUN:
            _nuke_and_adopt(resting_by_ticker)
        return resting_by_ticker

    exposure = db.compute_total_exposure()  # For logging only

    quoted_count = 0
    quoted_events = set()
    pending_placements = []   # New orders to batch-place after the loop
    pending_cancel_ids = []   # Order IDs to batch-cancel after the loop
    cycle_ts = time.time()

    # Pre-compute maker sizes using budget-based approach.
    # Pass 1: Run quoter on all markets to get posting prices.
    # Pass 2: Kelly sizes at the quoter's actual price (not book price),
    # so that contract_count × posting_price never exceeds the dollar budget.
    kelly_sizes = {}  # ticker -> {"bid_size": int, "ask_size": int}
    # Pre-pass: run quoter on all markets to get posting prices, then
    # pass those prices to Kelly so dollar→contract conversion uses the
    # actual posting price (not the stale book price).
    # Cache full quote dicts to avoid running the quoter twice per ticker.
    # NOTE: Do NOT call _detect_penny_loop here — it mutates state.
    # The main loop calls it; the rare case where penny-loop fires will
    # re-run the quoter with zeroed book prices.
    _quote_cache = {}  # ticker -> quote dict (from pre-pass, without penny-loop)
    if config.USE_KELLY_SIZING:
        quoter_prices = {}  # ticker -> {"bid_yes": int, "ask_yes": int}
        for mkt in quotable_markets:
            tk = mkt["ticker"]
            if not risk.check_tipoff_proximity(mkt.get("commence_time")):
                continue
            pre_quote = quoter.compute_quotes(
                mkt["fair_prob"], kelly.get_net_position(tk),
                book_bid=mkt.get("book_bid", 0),
                book_ask=mkt.get("book_ask", 0))
            if pre_quote:
                quoter_prices[tk] = {
                    "bid_yes": pre_quote["bid_yes"],
                    "ask_yes": pre_quote["ask_yes"],
                }
                _quote_cache[tk] = pre_quote

        positions_map = {p.get("ticker", ""): p.get("net_yes", 0)
                         for p in kelly.get_cached_positions() or []}
        kelly_sizes = kelly.compute_maker_sizes(
            quotable_markets, positions_map,
            config.BANKROLL, config.KELLY_FRACTION,
            config.MAX_GAME_TYPE_EXPOSURE_PCT,
            quoter_prices)

    for market in quotable_markets:
        ticker = market["ticker"]
        event_ticker = market.get("event_ticker", "")

        # Check tipoff proximity
        if not risk.check_tipoff_proximity(market.get("commence_time")):
            # Cancel if we have resting orders for this ticker
            if ticker in resting_by_ticker:
                if not DRY_RUN:
                    for side_key in ["bid_order_id", "ask_order_id"]:
                        oid = resting_by_ticker[ticker].get(side_key)
                        if oid:
                            pending_cancel_ids.append(oid)
                del resting_by_ticker[ticker]
            continue

        # Skip markets with stale book data (use cycle start time, not wall clock,
        # since order placement for earlier markets can take 10+ seconds)
        book_age = cycle_ts - market.get("book_fetched_at", 0)
        if book_age > config.MAX_BOOK_STALENESS_SEC:
            continue

        # Get current position for inventory skew
        net = kelly.get_net_position(ticker)
        home_team = market.get("home_team", "")
        away_team = market.get("away_team", "")
        game_key = market.get("game_key", (home_team, away_team))

        # Budget-based sizing: look up pre-computed sizes (budget per game+type)
        if config.USE_KELLY_SIZING:
            ks = kelly_sizes.get(ticker, {"bid_size": 0, "ask_size": 0})
            bid_size = ks["bid_size"]
            ask_size = ks["ask_size"]
        else:
            bid_size = config.CONTRACT_SIZE
            ask_size = config.CONTRACT_SIZE

        # Kelly says don't quote either side → skip entirely
        if bid_size <= 0 and ask_size <= 0:
            if ticker in resting_by_ticker:
                if not DRY_RUN:
                    for side_key in ["bid_order_id", "ask_order_id"]:
                        oid = resting_by_ticker[ticker].get(side_key)
                        if oid:
                            pending_cancel_ids.append(oid)
                del resting_by_ticker[ticker]
            continue

        # Compute quote (orderbook-aware, with anti-penny-loop).
        # Reuse pre-pass cache when possible; re-run only if penny-loop fires
        # (which zeros book prices, producing different quotes).
        book_bid = market.get("book_bid", 0)
        book_ask = market.get("book_ask", 0)
        penny_loop = quoter._detect_penny_loop(ticker, book_bid, book_ask)
        if penny_loop:
            book_bid, book_ask = 0, 0

        cached = _quote_cache.get(ticker)
        if cached and not penny_loop:
            # Reuse pre-pass quote, just update sizes
            quote = dict(cached)
            quote["bid_size"] = bid_size
            quote["ask_size"] = ask_size
        else:
            quote = quoter.compute_quotes(
                market["fair_prob"], net,
                book_bid=book_bid, book_ask=book_ask,
                bid_size=bid_size, ask_size=ask_size
            )
        if quote is None:
            # Cancel existing orders if any
            if ticker in resting_by_ticker:
                if not DRY_RUN:
                    for side_key in ["bid_order_id", "ask_order_id"]:
                        oid = resting_by_ticker[ticker].get(side_key)
                        if oid:
                            pending_cancel_ids.append(oid)
                del resting_by_ticker[ticker]
            continue

        # Log the quote
        db.log_quote(ticker, market["fair_prob"], quote["bid_yes"], quote["ask_yes"],
                     net, quote["skew"], pred_age,
                     book_bid=market.get("book_bid", 0),
                     book_ask=market.get("book_ask", 0))

        if DRY_RUN:
            print(quoter.format_quote_summary(ticker, quote))
            quoted_count += 1
            quoted_events.add(event_ticker)
            continue

        # Manage bid and ask: amend existing, collect new placements, collect cancels
        existing = resting_by_ticker.get(ticker, {})

        if bid_size > 0:
            bid_oid = existing.get("bid_order_id")
            if bid_oid:
                # Amend if price changed
                if quoter.should_amend(existing.get("bid_price", 0), quote["bid_yes"]):
                    result = orders.amend_order(bid_oid, ticker, "yes", price=quote["bid_yes"], count=quote["bid_size"])
                    if result:
                        existing["bid_price"] = quote["bid_yes"]
                        db.save_resting_order(bid_oid, ticker, "yes", "buy",
                                              quote["bid_yes"], quote["bid_size"])
            else:
                # Defer new bid to batch placement
                pending_placements.append({
                    "ticker": ticker, "side": "yes",
                    "price": quote["bid_yes"], "count": quote["bid_size"],
                    "post_only": True,
                    "_side_key": "bid",  # internal: which side to wire up
                })
        elif existing.get("bid_order_id"):
            pending_cancel_ids.append(existing["bid_order_id"])
            existing.pop("bid_order_id", None)

        # Manage ask (buy NO, which is selling YES)
        if ask_size > 0:
            ask_oid = existing.get("ask_order_id")
            no_price = 100 - quote["ask_yes"]  # Buying NO at this price = selling YES at ask_yes
            if ask_oid:
                if quoter.should_amend(existing.get("ask_price", 0), quote["ask_yes"]):
                    result = orders.amend_order(ask_oid, ticker, "no", price=no_price, count=quote["ask_size"])
                    if result:
                        existing["ask_price"] = quote["ask_yes"]
                        db.save_resting_order(ask_oid, ticker, "no", "buy",
                                              no_price, quote["ask_size"])
            else:
                pending_placements.append({
                    "ticker": ticker, "side": "no",
                    "price": no_price, "count": quote["ask_size"],
                    "post_only": True,
                    "_side_key": "ask",
                    "_ask_yes": quote["ask_yes"],  # store for metadata
                })
        elif existing.get("ask_order_id"):
            pending_cancel_ids.append(existing["ask_order_id"])
            existing.pop("ask_order_id", None)

        # Store market metadata so poll_for_fills can find it even after
        # quotable list refreshes (prevents ghost positions with NULL home/away)
        existing["_home_team"] = market.get("home_team")
        existing["_away_team"] = market.get("away_team")
        existing["_market_type"] = market.get("market_type")
        existing["_event_ticker"] = event_ticker
        existing["_fair_prob"] = market.get("fair_prob")
        existing["_bid_size"] = quote["bid_size"]
        existing["_ask_size"] = quote["ask_size"]
        existing["_contract_team"] = market.get("contract_team")
        # Compute line_value in home-spread convention for position tracking
        mt = market.get("market_type", "spreads")
        strike = market.get("strike")
        if mt == "spreads" and strike is not None:
            existing["_line_value"] = -strike if market.get("is_home_contract") else strike
        elif mt == "totals" and strike is not None:
            existing["_line_value"] = strike
        else:
            existing["_line_value"] = None
        existing["_commence_time"] = market.get("commence_time")
        resting_by_ticker[ticker] = existing
        quoted_count += 1
        quoted_events.add(event_ticker)

    # --- Batch cancel collected order IDs ---
    if pending_cancel_ids:
        cancelled = orders.batch_cancel(pending_cancel_ids)
        db.batch_remove_resting_orders(cancelled)

    # --- Batch place collected new orders ---
    if pending_placements:
        results = orders.batch_place(pending_placements)
        db_specs = []
        for spec, result in zip(pending_placements, results):
            if not result:
                continue
            oid = result.get("order_id")
            if not oid:
                continue
            ticker = spec["ticker"]
            existing = resting_by_ticker.get(ticker, {})
            if spec["_side_key"] == "bid":
                existing["bid_order_id"] = oid
                existing["bid_price"] = spec["price"]
                db_specs.append({"order_id": oid, "ticker": ticker, "side": "yes",
                                 "action": "buy", "price": spec["price"], "count": spec["count"]})
            else:  # ask
                existing["ask_order_id"] = oid
                existing["ask_price"] = spec.get("_ask_yes", 100 - spec["price"])
                db_specs.append({"order_id": oid, "ticker": ticker, "side": "no",
                                 "action": "buy", "price": spec["price"], "count": spec["count"]})
            resting_by_ticker[ticker] = existing
        db.batch_save_resting_orders(db_specs)

    # Flush buffered quote logs in one write
    db.flush_quote_log()

    print(f"  Quoting {quoted_count} tickers across {len(quoted_events)} games (exposure: ${exposure:.2f})")
    if pending_placements:
        print(f"  Batch placed {len([r for r in results if r])} / {len(pending_placements)} new orders")
    return resting_by_ticker


def sweep_tipoff_cancel(resting_by_ticker):
    """Cancel all resting orders for games within tipoff proximity.

    Runs every quote cycle, independent of quotable_markets.
    Uses _commence_time stored on each resting order entry.
    Only removes tracking after confirmed cancel to prevent orphaned orders.
    """
    cancel_ids = []
    oid_to_ticker = {}  # map order ID → ticker for partial-success cleanup
    for ticker, info in list(resting_by_ticker.items()):
        ct = info.get("_commence_time")
        if not risk.check_tipoff_proximity(ct):
            print(f"  TIPOFF CANCEL: {ticker} (commence={ct})")
            for side_key in ["bid_order_id", "ask_order_id"]:
                oid = info.get(side_key)
                if oid:
                    cancel_ids.append(oid)
                    oid_to_ticker[oid] = ticker

    if not cancel_ids:
        return

    if DRY_RUN:
        for ticker in set(oid_to_ticker.values()):
            del resting_by_ticker[ticker]
        print(f"  Tipoff-cancelled {len(oid_to_ticker)} orders (dry run)")
        return

    cancelled = orders.batch_cancel(cancel_ids)
    db.batch_remove_resting_orders(cancelled)
    # Remove ticker only if ALL its orders were cancelled
    cancelled_tickers = set()
    failed_tickers = set()
    for oid, ticker in oid_to_ticker.items():
        if oid in cancelled:
            cancelled_tickers.add(ticker)
        else:
            failed_tickers.add(ticker)
    for ticker in cancelled_tickers - failed_tickers:
        del resting_by_ticker[ticker]
    if cancelled:
        print(f"  Tipoff-cancelled {len(cancelled_tickers - failed_tickers)} tickers")
    if failed_tickers:
        print(f"  WARNING: Tipoff cancel failed for {len(failed_tickers)} tickers — will retry next sweep")


def poll_for_fills(resting_by_ticker, quotable_markets_ref=None):
    """Check resting orders for fills and update positions."""
    global TOTAL_FILLS, TOTAL_FEES

    if DRY_RUN:
        return

    if quotable_markets_ref is None:
        quotable_markets_ref = []

    api_orders = orders.get_resting_orders()
    api_order_map = {o["order_id"]: o for o in api_orders}

    # SAFETY: Cancel any orders we don't recognize (phantom orders from
    # place_order responses that were lost — order went through on Kalshi
    # but we never stored the order_id, so we'd re-place every cycle)
    tracked_oids = set()
    for info in resting_by_ticker.values():
        for key in ["bid_order_id", "ask_order_id"]:
            oid = info.get(key)
            if oid:
                tracked_oids.add(oid)

    phantom_ids = []
    for api_order in api_orders:
        oid = api_order["order_id"]
        if oid not in tracked_oids:
            # Skip manual (user-placed) orders — only cancel bot phantoms
            if not orders.is_bot_order(api_order):
                continue
            print(f"  PHANTOM ORDER detected: {oid}")
            phantom_ids.append(oid)
    if phantom_ids:
        print(f"  Batch cancelling {len(phantom_ids)} phantom orders...")
        orders.batch_cancel(phantom_ids)

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
            cumulative_filled = original - remaining

            # Track previously seen fills to compute incremental delta
            prev_filled_key = f"{side_key}_filled"
            prev_filled = info.get(prev_filled_key, 0)
            new_fills = cumulative_filled - prev_filled

            if new_fills > 0:
                price = api_order.get("yes_price", 0) if side == "yes" else api_order.get("no_price", 0)

                print(f"  FILL: {side} {new_fills}x @ {price}c on {ticker}"
                      f"{' (partial)' if remaining > 0 else ''}")
                TOTAL_FILLS += new_fills

                # Find market metadata — prefer resting_by_ticker (survives quotable
                # refreshes) over quotable_markets_ref (may not contain this ticker
                # after a refresh, causing ghost positions with NULL home/away)
                resting_info = resting_by_ticker.get(ticker, {})
                event_ticker = resting_info.get("_event_ticker")
                fill_home = resting_info.get("_home_team")
                fill_away = resting_info.get("_away_team")
                fill_market_type = resting_info.get("_market_type")
                fill_fair_prob = resting_info.get("_fair_prob")
                fill_line_value = resting_info.get("_line_value")
                fill_contract_team = resting_info.get("_contract_team")

                # Fallback to quotable_markets_ref if resting metadata missing
                if not event_ticker:
                    for m in quotable_markets_ref:
                        if m["ticker"] == ticker:
                            event_ticker = m.get("event_ticker")
                            fill_home = fill_home or m.get("home_team")
                            fill_away = fill_away or m.get("away_team")
                            fill_market_type = fill_market_type or m.get("market_type")
                            fill_fair_prob = fill_fair_prob or m.get("fair_prob")
                            fill_contract_team = fill_contract_team or m.get("contract_team")
                            # Compute line_value from market if not in resting metadata
                            if fill_line_value is None:
                                mt = m.get("market_type", "spreads")
                                strike = m.get("strike")
                                if mt == "spreads" and strike is not None:
                                    fill_line_value = -strike if m.get("is_home_contract") else strike
                                elif mt == "totals" and strike is not None:
                                    fill_line_value = strike
                            break

                # Compute maker fee: fee_rate * P * (1-P) * 100 per contract
                p = price / 100.0
                maker_fee_per = config.MAKER_FEE_RATE * p * (1 - p) * 100
                total_fee = maker_fee_per * new_fills

                # Update position with incremental fills only
                db.update_position(ticker, side, price, new_fills,
                                   event_ticker=event_ticker,
                                   home_team=fill_home,
                                   away_team=fill_away,
                                   market_type=fill_market_type,
                                   fair_prob=fill_fair_prob,
                                   line_value=fill_line_value,
                                   contract_team=fill_contract_team)
                db.record_fill(
                    fill_id=f"{oid}-fill-{cumulative_filled}",
                    ticker=ticker, side=side, action="buy",
                    price=price, count=new_fills, fee_cents=total_fee,
                    order_id=oid
                )
                info[prev_filled_key] = cumulative_filled

                if remaining == 0:
                    # Fully filled — remove
                    db.remove_resting_order(oid)
                    info.pop(side_key, None)
                    info.pop(prev_filled_key, None)
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
    if config.USE_KELLY_SIZING:
        size_str = f"Kelly (f={config.KELLY_FRACTION}, B=${config.BANKROLL:.0f})"
    else:
        size_str = str(config.CONTRACT_SIZE)
    print(f"  Maker: EV>{config.MIN_EV_PCT:.0%} | Size: {size_str}")
    take_size_str = "Kelly-sized" if config.USE_KELLY_SIZING else str(config.TAKE_CONTRACT_SIZE)
    print(f"  Taker: EV>{config.MIN_TAKE_EV_PCT:.0%} (after {config.TAKER_FEE_RATE:.0%} fee) | Size: {take_size_str}")
    print(f"  Markets: {', '.join(sorted(config.ENABLED_MARKET_TYPES))}")
    print(f"  Risk: Kelly-only (no hard position limits)")
    print(f"  API: {config.KALSHI_BASE_URL}")
    print("=" * 60)

    # Init
    db.init_database()
    db.init_taker_tables()

    # SAFETY: Cancel any resting orders from previous sessions (crash recovery)
    # Then adopt any survivors into resting_by_ticker to prevent duplicates.
    resting_by_ticker = {}
    if not DRY_RUN:
        print("Checking for stale orders from previous sessions...")
        stale = orders.get_resting_orders()
        if stale:
            print(f"  Found {len(stale)} stale resting orders — cancelling all.")
            _nuke_and_adopt(resting_by_ticker)
        else:
            print("  No stale orders found.")

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
        _, age_sec = risk.check_staleness(prediction_updated_at)
        pred_age = f"{age_sec:.0f}s"
    print(f"  Loaded {len(predictions)} predictions (age: {pred_age})")

    # Load team resolution
    print("Loading team name resolution...")
    team_dict = load_team_dict("cbb")
    canonical_games = load_canonical_games("cbb")

    # Fetch Kalshi 1H markets (all enabled types)
    enabled = ", ".join(sorted(config.ENABLED_MARKET_TYPES))
    print(f"Fetching Kalshi 1H markets (enabled: {enabled})...")
    all_kalshi = fetch_all_markets()
    total_markets = sum(len(v) for v in all_kalshi.values())
    for mtype, markets in all_kalshi.items():
        print(f"  {mtype}: {len(markets)} open contracts")
    print(f"  Total: {total_markets} open markets")

    if not all_kalshi:
        print("No open Kalshi 1H markets found.")
        return

    # Match to predictions
    print("Matching Kalshi markets to predictions...")
    quotable = match_all_markets(all_kalshi, predictions, team_dict, canonical_games)
    enforce_monotonicity(quotable)
    # Count by type
    by_type = {}
    for m in quotable:
        t = m.get("market_type", "spreads")
        by_type[t] = by_type.get(t, 0) + 1
    type_summary = ", ".join(f"{t}={c}" for t, c in sorted(by_type.items()))
    print(f"  {len(quotable)} quotable markets ({type_summary})")

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

    # Pre-warm Kelly sample cache (avoids ~15 min cold start on first cycle)
    print("Pre-warming Kelly sample cache...")
    n_games = kelly.prewarm_sample_cache()
    print(f"  Cached samples for {n_games} games")

    # Main loop
    print(f"\nStarting main loop (quote every {config.QUOTE_CYCLE_SEC}s, "
          f"monitor every {config.MONITOR_CYCLE_SEC}s)...\n")

    # resting_by_ticker initialized above (startup cancel + adopt survivors)
    last_quote_time = 0
    last_fill_poll = 0
    last_monitor_time = 0
    last_pred_check = 0
    # If predictions are already stale at startup, trigger pipeline immediately
    _, pred_age = risk.check_staleness(prediction_updated_at)
    if pred_age > config.MAX_STALENESS_SEC:
        last_pipeline_time = 0
        print(f"  Predictions stale at startup ({pred_age:.0f}s) — triggering immediate pipeline refresh")
    else:
        last_pipeline_time = time.time()

    try:
        while RUNNING:
            now = time.time()

            # Check if background pipeline finished
            try:
                if check_pipeline_completion(resting_by_ticker):
                    new_preds, new_ts = db.load_predictions()
                    if new_preds and new_ts != prediction_updated_at:
                        predictions = new_preds
                        prediction_updated_at = new_ts
                        kelly.clear_sample_cache()  # Reload samples with new predictions

                    # Re-fetch Kalshi markets + re-match (all enabled types)
                    fresh_all = fetch_all_markets()
                    if fresh_all:
                        all_kalshi = fresh_all
                    new_quotable = match_all_markets(
                        all_kalshi, predictions, team_dict, canonical_games
                    )
                    enforce_monotonicity(new_quotable)
                    print(f"  Refreshed: {len(new_quotable)} quotable markets")

                    # Cancel orders for orphaned tickers
                    new_tickers = {m["ticker"] for m in new_quotable}
                    for ticker, info in list(resting_by_ticker.items()):
                        if ticker not in new_tickers:
                            print(f"  Orphaned ticker {ticker} — cancelling orders")
                            if not DRY_RUN:
                                for side_key in ["bid_order_id", "ask_order_id"]:
                                    oid = info.get(side_key)
                                    if oid:
                                        orders.cancel_order(oid)
                                        db.remove_resting_order(oid)
                            del resting_by_ticker[ticker]

                    quotable = new_quotable

                    # Update reference lines
                    ref_lines = risk.run_line_monitor()
                    if ref_lines:
                        db.save_reference_lines(ref_lines)
                    last_pipeline_time = now
            except Exception as e:
                print(f"  Pipeline refresh error (will retry): {e}")

            # Periodic prediction check — catches external pipeline runs
            if now - last_pred_check >= 30:
                last_pred_check = now
                try:
                    meta_ts = db.get_prediction_timestamp()
                    if meta_ts and meta_ts != prediction_updated_at:
                        new_preds, new_ts = db.load_predictions()
                        if new_preds and new_ts:
                            predictions = new_preds
                            prediction_updated_at = new_ts
                            kelly.clear_sample_cache()
                            fresh_all = fetch_all_markets()
                            if fresh_all:
                                all_kalshi = fresh_all
                            new_quotable = match_all_markets(
                                all_kalshi, predictions, team_dict, canonical_games
                            )
                            enforce_monotonicity(new_quotable)

                            # Cancel orders for orphaned tickers
                            new_tickers = {m["ticker"] for m in new_quotable}
                            for ticker, info in list(resting_by_ticker.items()):
                                if ticker not in new_tickers:
                                    print(f"  Orphaned ticker {ticker} — cancelling orders")
                                    if not DRY_RUN:
                                        for side_key in ["bid_order_id", "ask_order_id"]:
                                            oid = info.get(side_key)
                                            if oid:
                                                orders.cancel_order(oid)
                                                db.remove_resting_order(oid)
                                    del resting_by_ticker[ticker]

                            quotable = new_quotable
                            ref_lines = risk.run_line_monitor()
                            if ref_lines:
                                db.save_reference_lines(ref_lines)
                            print(f"  External predictions detected — refreshed {len(quotable)} markets, cache rewarmed")
                except Exception as e:
                    print(f"  Prediction check error: {e}")

            # Scheduled pipeline refresh (with backoff on failure)
            try:
                if now - last_pipeline_time >= _pipeline_backoff:
                    start_pipeline()
                    last_pipeline_time = now
            except Exception as e:
                print(f"  Pipeline start error: {e}")

            # Quote cycle — always poll fills first so position/skew is current
            if now - last_quote_time >= config.QUOTE_CYCLE_SEC:
                try:
                    poll_for_fills(resting_by_ticker, quotable)
                    last_fill_poll = now
                    refresh_book_data(quotable)  # Fresh orderbook each cycle
                    sweep_tipoff_cancel(resting_by_ticker)  # Cancel orders near tipoff
                    print(f"\n--- Quote cycle @ {datetime.now().strftime('%H:%M:%S')} ---")
                    resting_by_ticker = run_quote_cycle(quotable, resting_by_ticker, prediction_updated_at)
                except duckdb.Error as e:
                    print(f"  Quote cycle DB error — resetting connection: {e}")
                    db._reset_write_conn()
                except Exception as e:
                    print(f"  Quote cycle error (will retry next cycle): {e}")
                last_quote_time = now

            # Line-move monitoring
            if now - last_monitor_time >= config.MONITOR_CYCLE_SEC:
                try:
                    print(f"\n--- Line monitor @ {datetime.now().strftime('%H:%M:%S')} ---")
                    current_lines = risk.run_line_monitor()
                    ref_lines = db.get_reference_lines()

                    if current_lines and ref_lines:
                        moved = risk.detect_line_moves(current_lines, ref_lines)
                        if moved:
                            moved_set = set(moved)
                            print(f"  {len(moved_set)} games with line moves — pulling quotes, triggering refresh")
                            # Batch cancel orders for moved games
                            cancel_ids = []
                            oid_to_ticker = {}
                            for ticker, info in list(resting_by_ticker.items()):
                                game = (info.get("_home_team"), info.get("_away_team"))
                                if game in moved_set:
                                    for side_key in ["bid_order_id", "ask_order_id"]:
                                        oid = info.get(side_key)
                                        if oid:
                                            cancel_ids.append(oid)
                                            oid_to_ticker[oid] = ticker
                            if cancel_ids and not DRY_RUN:
                                cancelled = orders.batch_cancel(cancel_ids)
                                db.batch_remove_resting_orders(cancelled)
                                cancelled_tickers = set()
                                failed_tickers = set()
                                for oid, ticker in oid_to_ticker.items():
                                    if oid in cancelled:
                                        cancelled_tickers.add(ticker)
                                    else:
                                        failed_tickers.add(ticker)
                                for ticker in cancelled_tickers - failed_tickers:
                                    del resting_by_ticker[ticker]
                                if failed_tickers:
                                    print(f"  WARNING: Line-move cancel failed for {len(failed_tickers)} tickers — will retry")
                            elif DRY_RUN:
                                for ticker in set(oid_to_ticker.values()):
                                    del resting_by_ticker[ticker]
                            # Trigger immediate pipeline refresh
                            start_pipeline()
                except Exception as e:
                    print(f"  Line monitor error (will retry next cycle): {e}")
                last_monitor_time = now

            # Taker scan — runs every iteration (~1s) for fast execution
            try:
                taker.run_take_cycle(quotable, prediction_updated_at,
                                     dry_run=DRY_RUN,
                                     resting_by_ticker=resting_by_ticker)
            except Exception as e:
                print(f"  [TAKER] Cycle error (will retry): {e}")

            time.sleep(random.uniform(0.5, 1.5))  # Jitter to avoid predictable timing

    except KeyboardInterrupt:
        print("\n  Received interrupt signal, shutting down...")
    except Exception as e:
        print(f"\n  FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Kill background pipeline if running
        if _pipeline_proc and _pipeline_proc.poll() is None:
            _pipeline_proc.kill()
            print("  Killed background pipeline subprocess.")

        # Kill switch: always cancel ALL orders on Kalshi (catches phantom orders
        # from failed API responses and orders we lost track of)
        print("\n\nShutting down...")
        if not DRY_RUN:
            print("Cancelling all bot orders (manual orders preserved)...")
            for attempt in range(3):
                if orders.cancel_all_bot_orders():
                    break
                print(f"  Kill switch retry {attempt + 1}/3...")
                time.sleep(1)

        # Log session summary
        try:
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
        except Exception as e:
            print(f"  Shutdown DB error (non-fatal): {e}")
            print(f"  Session {SESSION_ID}: {TOTAL_FILLS} fills, ${TOTAL_FEES:.2f} fees")
        finally:
            db.close_write_conn()


if __name__ == "__main__":
    main()
