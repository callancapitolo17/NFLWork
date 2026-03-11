#!/usr/bin/env python3
"""
Bet auto-placer: Opens visible browser, navigates to bookmaker, pre-fills bet(s).
Called as subprocess by the dashboard server.

Usage:
    # Single bet (backward compat):
    python placer.py '{"bookmaker": "bfa", "home_team": "Duke", ...}'

    # Batch (all bets for one book):
    python placer.py '[{"bookmaker": "bfa", ...}, {"bookmaker": "bfa", ...}]'
"""

import sys
import json

from base_navigator import update_bet_status


def get_navigator(bookmaker: str):
    """Import and return the navigator class for a bookmaker."""
    if bookmaker == "wagerzon":
        from navigator_wagerzon import WagerzonNavigator
        return WagerzonNavigator()
    elif bookmaker == "hoop88":
        from navigator_hoop88 import Hoop88Navigator
        return Hoop88Navigator()
    elif bookmaker == "bfa":
        from navigator_bfa import BFANavigator
        return BFANavigator()
    else:
        return None


def main():
    if len(sys.argv) < 2:
        print("Usage: python placer.py '<bet_json_or_array>'")
        sys.exit(1)

    raw = json.loads(sys.argv[1])

    # Normalize to list
    if isinstance(raw, dict):
        bets = [raw]
    elif isinstance(raw, list):
        bets = raw
    else:
        print("Invalid input: expected JSON object or array")
        sys.exit(1)

    if not bets:
        print("No bets to place")
        sys.exit(0)

    bookmaker = bets[0].get("bookmaker", "")
    print(f"Auto-placing {len(bets)} bet(s) on {bookmaker}...", flush=True)

    navigator = get_navigator(bookmaker)
    if not navigator:
        print(f"No navigator for bookmaker: {bookmaker}", flush=True)
        for bet in bets:
            if bet.get("bet_hash"):
                update_bet_status(bet["bet_hash"], "nav_error")
        sys.exit(1)

    try:
        for bet in bets:
            if bet.get("bet_hash"):
                update_bet_status(bet["bet_hash"], "navigating")
        navigator.place_bets(bets)
    except Exception as e:
        print(f"Navigation failed: {e}", flush=True)
        for bet in bets:
            if bet.get("bet_hash"):
                update_bet_status(bet["bet_hash"], "nav_error")
        sys.exit(1)


if __name__ == "__main__":
    main()
