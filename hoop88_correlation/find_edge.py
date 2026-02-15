#!/usr/bin/env python3
"""
Correlated Parlay Edge Finder

This tool:
1. Scrapes hoop88 for football 1H/FG odds
2. Calculates uncorrelated parlay odds for spread+total combinations
3. Compares against Pikkit SGP odds (correlation-adjusted) from multiple books
4. Displays the edge (positive edge = hoop88 is offering +EV)
5. Calculates Kelly criterion bet sizing

Usage:
    python find_edge.py                        # Scrape hoop88, manual SGP entry
    python find_edge.py --pikkit               # Auto-fetch SGP odds from Pikkit
    python find_edge.py --pikkit --visible     # Show browser windows
    python find_edge.py --no-compare           # Just show hoop88 odds
"""

import argparse
from scraper_hoop88_odds import scrape_hoop88_1h_odds, scrape_hoop88_with_actual_odds, scrape_hoop88_all_markets

# Try to import Pikkit scraper
try:
    from scraper_pikkit import (
        load_session, is_logged_in, find_game_on_pikkit,
        get_sgp_odds_for_parlay, extract_parlay_odds, PIKKIT_URL
    )
    from playwright.sync_api import sync_playwright
    PIKKIT_AVAILABLE = True
except ImportError:
    PIKKIT_AVAILABLE = False

# Defaults for Kelly criterion
DEFAULT_BANKROLL = 26000
DEFAULT_KELLY_FRACTION = 0.25

# Defaults for SGP odds filtering and adjustment
DEFAULT_TRUSTED_BOOKS = ['prophetx', 'novig', 'fanduel', 'draftkings']
DEFAULT_MULTIPLIER = 1.1  # Conservative adjustment on decimal odds

# Line matching tolerance - prompt for manual entry if lines differ by more than this
LINE_MATCH_TOLERANCE = 0

# Hoop88 rounds parlay odds down to the nearest 5
# e.g., calculated +264 becomes +260, +315 becomes +310
# This is how they take their margin on parlays


def lines_match(hoop88_spread: float, pikkit_spread_str: str, hoop88_total: float, pikkit_total_str: str, tolerance: float = LINE_MATCH_TOLERANCE) -> bool:
    """
    Check if Pikkit lines are close enough to hoop88 for valid comparison.

    Args:
        hoop88_spread: Spread value from hoop88 (e.g., 3.0 or -5.0)
        pikkit_spread_str: Spread string from Pikkit (e.g., "+3.5" or "-4.5")
        hoop88_total: Total value from hoop88 (e.g., 20.5)
        pikkit_total_str: Total string from Pikkit (e.g., "u21.5" or "o42")
        tolerance: Max allowed difference (default 0.5)

    Returns:
        True if both spread and total are within tolerance
    """
    try:
        # Parse Pikkit spread (e.g., "+3.5" or "-2.5")
        pikkit_spread = float(pikkit_spread_str)
        spread_diff = abs(hoop88_spread - pikkit_spread)

        # Parse Pikkit total (e.g., "u21.5" or "o42")
        pikkit_total_val = float(pikkit_total_str[1:])  # Remove o/u prefix
        total_diff = abs(hoop88_total - pikkit_total_val)

        return spread_diff <= tolerance and total_diff <= tolerance
    except (ValueError, IndexError):
        return False


def american_to_decimal(american: int) -> float:
    """Convert American odds to decimal."""
    if american == 0:
        return 1.91
    if american > 0:
        return (american / 100) + 1
    else:
        return (100 / abs(american)) + 1


def decimal_to_american(decimal: float) -> int:
    """Convert decimal odds to American."""
    if decimal >= 2:
        return int(round((decimal - 1) * 100))
    else:
        return int(round(-100 / (decimal - 1)))


def get_adjusted_sgp_odds(book_odds: list, all_odds: list = None, trusted_books: list = None, multiplier: float = None) -> dict:
    """
    Filter to trusted books, average their odds, apply conservative multiplier.

    If book_odds is empty or has no trusted books, falls back to averaging all_odds.
    Note: Pikkit only shows the name of the best book - other books appear as logos.

    Args:
        book_odds: List of {'book': str, 'odds': int} from Pikkit (often only best book)
        all_odds: List of all odds (int) from all books (fallback)
        trusted_books: List of book names to filter to (case-insensitive)
        multiplier: Conservative adjustment to apply to decimal odds (e.g., 1.1)

    Returns:
        Dict with:
        - adjusted_odds: Final American odds after filtering and adjustment
        - avg_odds: Average American odds before adjustment
        - filtered_books: List of books that were included (or 'all' if using fallback)
        - filtered_odds: List of odds from filtered books
        - used_fallback: True if we used all_odds instead of book_odds
    """
    if trusted_books is None:
        trusted_books = DEFAULT_TRUSTED_BOOKS
    if multiplier is None:
        multiplier = DEFAULT_MULTIPLIER

    result = {'adjusted_odds': None, 'avg_odds': None, 'filtered_books': [], 'filtered_odds': [], 'used_fallback': False}

    # First try: filter book_odds to trusted books
    if book_odds:
        trusted_lower = [b.lower() for b in trusted_books]
        filtered = [o for o in book_odds if o.get('book', '').lower() in trusted_lower]

        if filtered:
            decimals = [american_to_decimal(o['odds']) for o in filtered]
            avg_decimal = sum(decimals) / len(decimals)
            avg_american = decimal_to_american(avg_decimal)
            adjusted_decimal = avg_decimal * multiplier
            adjusted_american = decimal_to_american(adjusted_decimal)

            return {
                'adjusted_odds': adjusted_american,
                'avg_odds': avg_american,
                'filtered_books': [o['book'] for o in filtered],
                'filtered_odds': [o['odds'] for o in filtered],
                'used_fallback': False
            }

    # Fallback: use average of all_odds
    # Pikkit only shows book names for the best odds, so we fall back to averaging
    # all available odds as a proxy for "market consensus"
    if all_odds:
        decimals = [american_to_decimal(o) for o in all_odds]
        avg_decimal = sum(decimals) / len(decimals)
        avg_american = decimal_to_american(avg_decimal)
        adjusted_decimal = avg_decimal * multiplier
        adjusted_american = decimal_to_american(adjusted_decimal)

        return {
            'adjusted_odds': adjusted_american,
            'avg_odds': avg_american,
            'filtered_books': ['all'],  # Indicates we used all books
            'filtered_odds': all_odds,
            'used_fallback': True
        }

    return result


def apply_hoop88_shave(calculated_odds: int) -> int:
    """
    Apply hoop88's odds rounding to get estimated actual parlay odds.
    Hoop88 rounds down to nearest 5 for plus odds.
    e.g., +264 -> +260, +315 -> +310
    """
    if calculated_odds > 0:
        # Round down to nearest 5
        return (calculated_odds // 5) * 5
    else:
        # For minus odds, round to make more negative (worse for bettor)
        # e.g., -245 becomes -250
        return ((calculated_odds - 4) // 5) * 5


def calculate_edge(hoop88_american: int, dk_american: int) -> float:
    """
    Calculate edge percentage.
    Positive edge means hoop88 is offering better odds than DK's correlated price.
    """
    hoop88_dec = american_to_decimal(hoop88_american)
    dk_dec = american_to_decimal(dk_american)

    # Edge = how much better hoop88 pays vs DK's "true" price
    edge = (hoop88_dec - dk_dec) / dk_dec * 100
    return edge


def implied_probability(american: int) -> float:
    """Convert American odds to implied probability."""
    decimal = american_to_decimal(american)
    return 1 / decimal


def calculate_kelly(hoop88_american: int, dk_american: int, bankroll: float, kelly_fraction: float) -> dict:
    """
    Calculate Kelly criterion bet size.

    Uses DK SGP odds as the "true" probability (since they account for correlation).
    Hoop88 odds represent what we're actually betting at.

    Kelly % = (bp - q) / b
    where:
        b = decimal odds - 1 (profit per unit wagered)
        p = true win probability (from DK)
        q = 1 - p

    Returns dict with kelly details.
    """
    hoop88_dec = american_to_decimal(hoop88_american)
    dk_dec = american_to_decimal(dk_american)

    # True probability from DK (they price correlation correctly)
    true_prob = 1 / dk_dec

    # b = profit per unit if we win at hoop88 odds
    b = hoop88_dec - 1
    p = true_prob
    q = 1 - p

    # Full Kelly percentage
    full_kelly = (b * p - q) / b if b > 0 else 0

    # Apply fractional Kelly
    fractional_kelly = full_kelly * kelly_fraction

    # Don't bet if Kelly is negative (no edge)
    if fractional_kelly <= 0:
        return {
            'full_kelly_pct': 0,
            'fractional_kelly_pct': 0,
            'recommended_bet': 0,
            'bankroll_pct': 0,
            'has_edge': False
        }

    recommended_bet = bankroll * fractional_kelly

    return {
        'full_kelly_pct': full_kelly * 100,
        'fractional_kelly_pct': fractional_kelly * 100,
        'recommended_bet': recommended_bet,
        'bankroll_pct': fractional_kelly * 100,
        'has_edge': True
    }


def display_games(games: list, actual_odds: bool = False):
    """Display scraped games with parlay odds."""
    print("\n" + "=" * 70)
    print("HOOP88 CORRELATED PARLAY OPPORTUNITIES")
    print("=" * 70)
    if actual_odds:
        print("(Using ACTUAL parlay odds from hoop88 bet slip)")
    else:
        print("(Note: Hoop88 rounds down to nearest 5, e.g. +264 -> +260)")

    for game in games:
        period = game.get('period', 'FG')
        print(f"\n{game['league']} {period} | {game['away_team']} @ {game['home_team']}")
        print("-" * 50)

        # Display raw odds
        print(f"  Spread: {game['away_team']} {game['away_spread']:+.1f} ({game['away_spread_odds']:+d}) | "
              f"{game['home_team']} {game['home_spread']:+.1f} ({game['home_spread_odds']:+d})")
        print(f"  Total: {game['total']} | Over ({game['over_odds']:+d}) | Under ({game['under_odds']:+d})")

        print(f"\n  CORRELATED PARLAYS:")

        # Underdog + Under (positive correlation)
        actual_uu = game['parlay_underdog_under']
        calc_uu = game.get('calc_underdog_under', actual_uu)
        print(f"    {game['underdog']} {game['underdog_spread']:+.1f} + Under {game['total']}")
        if actual_odds:
            print(f"       Actual: {actual_uu:+d} (calc: {calc_uu:+d}) "
                  f"({implied_probability(actual_uu)*100:.1f}% implied)")
        else:
            est_uu = apply_hoop88_shave(actual_uu)
            print(f"       Calculated: {actual_uu:+d} -> Est. Actual: {est_uu:+d} "
                  f"({implied_probability(est_uu)*100:.1f}% implied)")

        # Favorite + Over (positive correlation)
        actual_fo = game['parlay_favorite_over']
        calc_fo = game.get('calc_favorite_over', actual_fo)
        print(f"    {game['favorite']} {game['favorite_spread']:+.1f} + Over {game['total']}")
        if actual_odds:
            print(f"       Actual: {actual_fo:+d} (calc: {calc_fo:+d}) "
                  f"({implied_probability(actual_fo)*100:.1f}% implied)")
        else:
            est_fo = apply_hoop88_shave(actual_fo)
            print(f"       Calculated: {actual_fo:+d} -> Est. Actual: {est_fo:+d} "
                  f"({implied_probability(est_fo)*100:.1f}% implied)")

    print("\n" + "=" * 70)


def pikkit_compare(games: list, bankroll: float = DEFAULT_BANKROLL, kelly_fraction: float = DEFAULT_KELLY_FRACTION,
                   headless: bool = True, trusted_books: list = None, multiplier: float = None,
                   tolerance: float = LINE_MATCH_TOLERANCE):
    """
    Automatically fetch SGP odds from Pikkit and compare with hoop88.

    Args:
        games: List of games from hoop88 scraper
        bankroll: Bankroll for Kelly sizing
        kelly_fraction: Fraction of Kelly to use
        headless: Run browser in headless mode
        trusted_books: List of book names to filter SGP odds to
        multiplier: Conservative multiplier for decimal odds
    """
    if trusted_books is None:
        trusted_books = DEFAULT_TRUSTED_BOOKS
    if multiplier is None:
        multiplier = DEFAULT_MULTIPLIER
    if not PIKKIT_AVAILABLE:
        print("Pikkit scraper not available. Install playwright and check scraper_pikkit.py")
        return

    print("\n" + "=" * 70)
    print("FETCHING SGP ODDS FROM PIKKIT")
    print("=" * 70)

    edges_found = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)

        saved_session = load_session()
        if not saved_session:
            print("No Pikkit session found. Run: python scraper_pikkit.py --visible")
            browser.close()
            return

        context = browser.new_context(viewport={'width': 1920, 'height': 1080}, storage_state=saved_session)
        page = context.new_page()

        page.goto(PIKKIT_URL, wait_until='networkidle')
        page.wait_for_timeout(2000)

        if not is_logged_in(page):
            print("Pikkit session expired. Run: python scraper_pikkit.py --visible")
            browser.close()
            return

        print("Pikkit logged in!")

        # Group games by matchup to avoid navigating to same game multiple times
        game_groups = {}
        for game in games:
            key = f"{game['away_team']}@{game['home_team']}"
            if key not in game_groups:
                game_groups[key] = []
            game_groups[key].append(game)

        for matchup_key, matchup_games in game_groups.items():
            first_game = matchup_games[0]
            away_team = first_game['away_team']
            home_team = first_game['home_team']
            league = first_game.get('league', 'NFL')

            print(f"\n{'='*60}")
            print(f"{away_team} @ {home_team}")
            print(f"{'='*60}")

            # Find game on Pikkit
            game_url = find_game_on_pikkit(page, away_team, home_team, league)

            if not game_url:
                print(f"  Could not find game on Pikkit")
                continue

            # Navigate to game page
            if page.url != game_url:
                page.goto(game_url, wait_until='networkidle')
                page.wait_for_timeout(2000)

            # Process each period for this matchup
            for game in matchup_games:
                period = game.get('period', 'FG')
                period_map = {'FG': 'Full', '1H': '1st Half', '1Q': '1st Quarter',
                              '2Q': '2nd Quarter', '3Q': '3rd Quarter', '4Q': '4th Quarter'}
                pikkit_period = period_map.get(period, 'Full')

                print(f"\n  {league} {period}")
                print(f"  {'-'*40}")

                # Reload page to ensure clean state for each period
                page.goto(game_url, wait_until='networkidle')
                page.wait_for_timeout(2000)

                # Get hoop88 odds
                hoop88_uu = game['parlay_underdog_under']
                hoop88_fo = game['parlay_favorite_over']

                underdog = game['underdog']
                underdog_spread = game['underdog_spread']
                favorite = game['favorite']
                favorite_spread = game['favorite_spread']
                total = game['total']

                # === Underdog + Under ===
                parlay1_desc = f"{underdog} {underdog_spread:+.1f} + Under {total}"
                print(f"\n    {parlay1_desc}")
                print(f"    Hoop88: {hoop88_uu:+d}")

                # Get Pikkit SGP odds
                sgp_uu = get_sgp_odds_for_parlay(page, pikkit_period, underdog, underdog_spread, total, "Under")
                book_odds_uu = sgp_uu.get('book_odds', [])
                all_odds = sgp_uu.get('all_odds', [])

                pikkit_spread = sgp_uu.get('pikkit_spread', '')
                pikkit_total = sgp_uu.get('pikkit_total', '')

                # Get adjusted odds (filter to trusted books, average, apply multiplier)
                # Falls back to averaging all_odds if no trusted books found
                adjusted_result = get_adjusted_sgp_odds(book_odds_uu, all_odds, trusted_books, multiplier)
                pikkit_uu = adjusted_result['adjusted_odds']

                if pikkit_uu and pikkit_spread and pikkit_total:
                    # Check if lines match within tolerance
                    if lines_match(underdog_spread, pikkit_spread, total, pikkit_total, tolerance):
                        # Lines match - show adjusted odds with book breakdown
                        if adjusted_result['used_fallback']:
                            # Used all_odds average
                            print(f"    Pikkit SGP: {pikkit_uu:+d} (adj) [avg of all: {adjusted_result['avg_odds']:+d}]")
                        else:
                            # Used trusted books
                            filtered_str = ', '.join([f"{b}:{o:+d}" for b, o in zip(adjusted_result['filtered_books'], adjusted_result['filtered_odds'])])
                            print(f"    Pikkit SGP: {pikkit_uu:+d} (adj) [avg:{adjusted_result['avg_odds']:+d} | {filtered_str}]")
                    else:
                        # Lines don't match - prompt for manual entry
                        print(f"    ⚠️  LINE MISMATCH:")
                        print(f"       Hoop88: {underdog_spread:+.1f} + Under {total}")
                        print(f"       Pikkit: {pikkit_spread} + {pikkit_total}")
                        manual_odds = input(f"    Enter SGP odds (or Enter to skip): ").strip()
                        if manual_odds:
                            try:
                                pikkit_uu = int(manual_odds.replace('+', ''))
                                print(f"    Manual SGP: {pikkit_uu:+d}")
                            except ValueError:
                                print("    Invalid input, skipping...")
                                pikkit_uu = None
                        else:
                            pikkit_uu = None

                if pikkit_uu:
                    edge = calculate_edge(hoop88_uu, pikkit_uu)
                    kelly = calculate_kelly(hoop88_uu, pikkit_uu, bankroll, kelly_fraction)

                    if edge > 0:
                        print(f"    --> EDGE: {edge:+.1f}%  |  Kelly: ${kelly['recommended_bet']:,.0f}")
                        edges_found.append({
                            'game': f"{away_team} @ {home_team}",
                            'period': period,
                            'parlay': parlay1_desc,
                            'hoop88': hoop88_uu,
                            'sgp': pikkit_uu,
                            'sgp_avg': adjusted_result.get('avg_odds'),
                            'filtered_books': adjusted_result.get('filtered_books', []),
                            'pikkit_lines': f"{pikkit_spread} + {pikkit_total}",
                            'all_sgp': all_odds,
                            'edge': edge,
                            'kelly': kelly
                        })
                    else:
                        print(f"    --> No edge: {edge:+.1f}%")
                else:
                    print(f"    Pikkit SGP: Not available")

                # === Favorite + Over ===
                # Reload page to clear betslip (clicking cells to deselect doesn't work reliably)
                page.goto(game_url, wait_until='networkidle')
                page.wait_for_timeout(1500)

                parlay2_desc = f"{favorite} {favorite_spread:+.1f} + Over {total}"
                print(f"\n    {parlay2_desc}")
                print(f"    Hoop88: {hoop88_fo:+d}")

                sgp_fo = get_sgp_odds_for_parlay(page, pikkit_period, favorite, favorite_spread, total, "Over")
                book_odds_fo = sgp_fo.get('book_odds', [])
                all_odds_fo = sgp_fo.get('all_odds', [])

                pikkit_spread_fo = sgp_fo.get('pikkit_spread', '')
                pikkit_total_fo = sgp_fo.get('pikkit_total', '')

                # Get adjusted odds (filter to trusted books, average, apply multiplier)
                # Falls back to averaging all_odds if no trusted books found
                adjusted_result_fo = get_adjusted_sgp_odds(book_odds_fo, all_odds_fo, trusted_books, multiplier)
                pikkit_fo = adjusted_result_fo['adjusted_odds']

                if pikkit_fo and pikkit_spread_fo and pikkit_total_fo:
                    # Check if lines match within tolerance
                    if lines_match(favorite_spread, pikkit_spread_fo, total, pikkit_total_fo, tolerance):
                        # Lines match - show adjusted odds with book breakdown
                        if adjusted_result_fo['used_fallback']:
                            # Used all_odds average
                            print(f"    Pikkit SGP: {pikkit_fo:+d} (adj) [avg of all: {adjusted_result_fo['avg_odds']:+d}]")
                        else:
                            # Used trusted books
                            filtered_str = ', '.join([f"{b}:{o:+d}" for b, o in zip(adjusted_result_fo['filtered_books'], adjusted_result_fo['filtered_odds'])])
                            print(f"    Pikkit SGP: {pikkit_fo:+d} (adj) [avg:{adjusted_result_fo['avg_odds']:+d} | {filtered_str}]")
                    else:
                        # Lines don't match - prompt for manual entry
                        print(f"    ⚠️  LINE MISMATCH:")
                        print(f"       Hoop88: {favorite_spread:+.1f} + Over {total}")
                        print(f"       Pikkit: {pikkit_spread_fo} + {pikkit_total_fo}")
                        manual_odds = input(f"    Enter SGP odds (or Enter to skip): ").strip()
                        if manual_odds:
                            try:
                                pikkit_fo = int(manual_odds.replace('+', ''))
                                print(f"    Manual SGP: {pikkit_fo:+d}")
                            except ValueError:
                                print("    Invalid input, skipping...")
                                pikkit_fo = None
                        else:
                            pikkit_fo = None

                if pikkit_fo:
                    edge = calculate_edge(hoop88_fo, pikkit_fo)
                    kelly = calculate_kelly(hoop88_fo, pikkit_fo, bankroll, kelly_fraction)

                    if edge > 0:
                        print(f"    --> EDGE: {edge:+.1f}%  |  Kelly: ${kelly['recommended_bet']:,.0f}")
                        edges_found.append({
                            'game': f"{away_team} @ {home_team}",
                            'period': period,
                            'parlay': parlay2_desc,
                            'hoop88': hoop88_fo,
                            'sgp': pikkit_fo,
                            'sgp_avg': adjusted_result_fo.get('avg_odds'),
                            'filtered_books': adjusted_result_fo.get('filtered_books', []),
                            'pikkit_lines': f"{pikkit_spread_fo} + {pikkit_total_fo}",
                            'all_sgp': all_odds_fo,
                            'edge': edge,
                            'kelly': kelly
                        })
                    else:
                        print(f"    --> No edge: {edge:+.1f}%")
                else:
                    print(f"    Pikkit SGP: Not available")

        browser.close()

    # Summary of edges found
    if edges_found:
        print("\n" + "=" * 70)
        print("RECOMMENDED BETS - SORTED BY EDGE")
        print("=" * 70)
        print(f"Bankroll: ${bankroll:,.0f} | Kelly: {kelly_fraction:.0%}")
        edges_found.sort(key=lambda x: x['edge'], reverse=True)

        total_recommended = sum(e['kelly']['recommended_bet'] for e in edges_found)

        for e in edges_found:
            print(f"\n  {e['edge']:+.1f}% EDGE -> Bet ${e['kelly']['recommended_bet']:,.0f}")
            print(f"    {e['game']} ({e['period']})")
            print(f"    {e['parlay']}")
            print(f"    Hoop88: {e['hoop88']:+d} vs Pikkit: {e['sgp']:+d}")

        print(f"\n" + "-" * 70)
        print(f"  TOTAL RECOMMENDED ACTION: ${total_recommended:,.0f} ({total_recommended/bankroll*100:.1f}% of bankroll)")
    else:
        print("\n" + "=" * 70)
        print("No edges found.")
        print("=" * 70)


def interactive_compare(games: list, bankroll: float = DEFAULT_BANKROLL, kelly_fraction: float = DEFAULT_KELLY_FRACTION, actual_odds: bool = False):
    """Interactive mode to compare with DraftKings SGP odds (manual entry)."""
    print("\n" + "=" * 70)
    print("COMPARE WITH DRAFTKINGS SGP")
    print("=" * 70)
    print(f"\nBankroll: ${bankroll:,.0f} | Kelly Fraction: {kelly_fraction:.0%}")
    if actual_odds:
        print("(Using ACTUAL Hoop88 parlay odds from bet slip)")
    else:
        print("(Using estimated actual Hoop88 odds - rounded down to nearest 5)")
    print("Enter the DraftKings SGP odds for each parlay (or 'skip' to skip):\n")

    edges_found = []

    for game in games:
        period = game.get('period', 'FG')
        print(f"\n{game['league']} {period} | {game['away_team']} @ {game['home_team']}")
        print("-" * 50)

        # Parlay 1: Underdog + Under
        if actual_odds:
            hoop88_uu = game['parlay_underdog_under']
            calc_uu = game.get('calc_underdog_under', hoop88_uu)
        else:
            calc_uu = game['parlay_underdog_under']
            hoop88_uu = apply_hoop88_shave(calc_uu)
        parlay1_desc = f"{game['underdog']} {game['underdog_spread']:+.1f} + Under {game['total']}"
        print(f"\n  {parlay1_desc}")
        print(f"  Hoop88: {hoop88_uu:+d} (calc: {calc_uu:+d})")

        dk_input = input(f"  Enter DK SGP odds (e.g., +220): ").strip()
        if dk_input.lower() != 'skip' and dk_input:
            try:
                dk_odds = int(dk_input.replace('+', ''))
                edge = calculate_edge(hoop88_uu, dk_odds)
                kelly = calculate_kelly(hoop88_uu, dk_odds, bankroll, kelly_fraction)

                if edge > 0:
                    print(f"  ✅ EDGE: {edge:+.1f}%")
                    print(f"  💰 Kelly: ${kelly['recommended_bet']:,.0f} ({kelly['bankroll_pct']:.2f}% of bankroll)")
                    edges_found.append({
                        'game': f"{game['away_team']} @ {game['home_team']}",
                        'parlay': parlay1_desc,
                        'hoop88': hoop88_uu,
                        'dk': dk_odds,
                        'edge': edge,
                        'kelly': kelly
                    })
                else:
                    print(f"  ❌ No edge: {edge:+.1f}%")
            except ValueError:
                print(f"  Could not parse odds: {dk_input}")

        # Parlay 2: Favorite + Over
        if actual_odds:
            hoop88_fo = game['parlay_favorite_over']
            calc_fo = game.get('calc_favorite_over', hoop88_fo)
        else:
            calc_fo = game['parlay_favorite_over']
            hoop88_fo = apply_hoop88_shave(calc_fo)
        parlay2_desc = f"{game['favorite']} {game['favorite_spread']:+.1f} + Over {game['total']}"
        print(f"\n  {parlay2_desc}")
        print(f"  Hoop88: {hoop88_fo:+d} (calc: {calc_fo:+d})")

        dk_input = input(f"  Enter DK SGP odds (e.g., +220): ").strip()
        if dk_input.lower() != 'skip' and dk_input:
            try:
                dk_odds = int(dk_input.replace('+', ''))
                edge = calculate_edge(hoop88_fo, dk_odds)
                kelly = calculate_kelly(hoop88_fo, dk_odds, bankroll, kelly_fraction)

                if edge > 0:
                    print(f"  ✅ EDGE: {edge:+.1f}%")
                    print(f"  💰 Kelly: ${kelly['recommended_bet']:,.0f} ({kelly['bankroll_pct']:.2f}% of bankroll)")
                    edges_found.append({
                        'game': f"{game['away_team']} @ {game['home_team']}",
                        'parlay': parlay2_desc,
                        'hoop88': hoop88_fo,
                        'dk': dk_odds,
                        'edge': edge,
                        'kelly': kelly
                    })
                else:
                    print(f"  ❌ No edge: {edge:+.1f}%")
            except ValueError:
                print(f"  Could not parse odds: {dk_input}")

    # Summary of edges found
    if edges_found:
        print("\n" + "=" * 70)
        print("RECOMMENDED BETS - SORTED BY EDGE")
        print("=" * 70)
        print(f"Bankroll: ${bankroll:,.0f} | Kelly: {kelly_fraction:.0%}")
        edges_found.sort(key=lambda x: x['edge'], reverse=True)

        total_recommended = sum(e['kelly']['recommended_bet'] for e in edges_found)

        for e in edges_found:
            print(f"\n  {e['edge']:+.1f}% EDGE -> Bet ${e['kelly']['recommended_bet']:,.0f}")
            print(f"    {e['game']}")
            print(f"    {e['parlay']}")
            print(f"    Hoop88: {e['hoop88']:+d} vs DK: {e['dk']:+d}")

        print(f"\n" + "-" * 70)
        print(f"  TOTAL RECOMMENDED ACTION: ${total_recommended:,.0f} ({total_recommended/bankroll*100:.1f}% of bankroll)")
    else:
        print("\n  No edges found.")


def main():
    parser = argparse.ArgumentParser(description='Find correlated parlay edge vs SGP odds')
    parser.add_argument('--pikkit', action='store_true', help='Auto-fetch SGP odds from Pikkit')
    parser.add_argument('--no-compare', action='store_true', help='Skip SGP comparison (just show hoop88 odds)')
    parser.add_argument('--fast', action='store_true', help='Use estimated odds instead of actual (faster but less precise)')
    parser.add_argument('--visible', action='store_true', help='Show browser window')
    parser.add_argument('--leagues', nargs='+', default=['NFL', 'NCAAF'],
                        help='Leagues to scrape (default: NFL NCAAF)')
    parser.add_argument('--periods', nargs='+', default=['FG', '1H', '1Q'],
                        help='Periods to scrape: FG, 1H, 1Q, 2Q, 3Q, 4Q (default: FG 1H 1Q)')
    parser.add_argument('--bankroll', type=float, default=DEFAULT_BANKROLL,
                        help=f'Your bankroll for Kelly sizing (default: ${DEFAULT_BANKROLL:,})')
    parser.add_argument('--kelly-fraction', type=float, default=DEFAULT_KELLY_FRACTION,
                        help=f'Fraction of Kelly to use (default: {DEFAULT_KELLY_FRACTION})')
    parser.add_argument('--books', nargs='+', default=DEFAULT_TRUSTED_BOOKS,
                        help=f'Trusted books to filter SGP odds to (default: {" ".join(DEFAULT_TRUSTED_BOOKS)})')
    parser.add_argument('--multiplier', type=float, default=DEFAULT_MULTIPLIER,
                        help=f'Conservative multiplier for decimal odds (default: {DEFAULT_MULTIPLIER})')
    parser.add_argument('--tolerance', type=float, default=LINE_MATCH_TOLERANCE,
                        help=f'Line match tolerance for spread/total comparison (default: {LINE_MATCH_TOLERANCE})')
    args = parser.parse_args()

    # Defaults: actual=True, compare=True
    use_actual = not args.fast
    do_compare = not args.no_compare

    print("=" * 70)
    print("CORRELATED PARLAY EDGE FINDER")
    print("=" * 70)
    print(f"Leagues: {', '.join(args.leagues)}")
    print(f"Periods: {', '.join(args.periods)}")
    if args.pikkit:
        print("SGP Source: Pikkit (auto)")
        print(f"Trusted Books: {', '.join(args.books)}")
        print(f"Multiplier: {args.multiplier}x")
        print(f"Line Tolerance: {args.tolerance}")
    else:
        print("SGP Source: Manual entry")

    if use_actual:
        print("\nScraping hoop88 for ACTUAL parlay odds (interacting with bet slip)...")
        games = scrape_hoop88_all_markets(
            leagues=args.leagues,
            periods=args.periods,
            headless=not args.visible
        )
    else:
        print("\nScraping hoop88 for current football odds (fast mode - estimated odds)...")
        games = scrape_hoop88_1h_odds(headless=not args.visible)

    if not games:
        print("No games found. Check if there are any football games available.")
        return

    display_games(games, actual_odds=use_actual)

    if args.pikkit:
        # Auto-fetch SGP odds from Pikkit
        pikkit_compare(games, bankroll=args.bankroll, kelly_fraction=args.kelly_fraction,
                       headless=not args.visible, trusted_books=args.books, multiplier=args.multiplier,
                       tolerance=args.tolerance)
    elif do_compare:
        # Manual SGP entry
        interactive_compare(games, bankroll=args.bankroll, kelly_fraction=args.kelly_fraction, actual_odds=use_actual)
    else:
        print(f"\nKelly defaults: Bankroll=${args.bankroll:,.0f}, Fraction={args.kelly_fraction:.0%}")
        print("Run with --pikkit to auto-fetch SGP odds, or without --no-compare for manual entry.")


if __name__ == "__main__":
    main()
