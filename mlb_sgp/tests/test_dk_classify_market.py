"""Tests classify_market — DK market-name → (period, market_type) classifier."""
import pytest

from mlb_sgp.scraper_draftkings_singles import classify_market


# ---------------------------------------------------------------------------
# Real game-level markets — should classify
# ---------------------------------------------------------------------------

def test_fg_main_markets_classify():
    assert classify_market("Total Runs") == ("FG", "main")
    assert classify_market("Run Line") == ("FG", "main")
    assert classify_market("Moneyline") == ("FG", "main")


def test_fg_alternate_markets_classify():
    assert classify_market("Total Alternate") == ("FG", "alternate_totals")
    assert classify_market("Run Line Alternate") == ("FG", "alternate_spreads")


def test_f5_markets_classify():
    assert classify_market("Total Runs - 1st 5 Innings") == ("F5", "main")
    assert classify_market("Run Line - 1st 5 Innings") == ("F5", "main")


def test_f7_markets_classify():
    assert classify_market("Total Runs - 1st 7 Innings") == ("F7", "main")
    assert classify_market("Run Line - 1st 7 Innings") == ("F7", "main")


# ---------------------------------------------------------------------------
# F3 (newly enabled) — DK posts these and the scraper now collects them
# ---------------------------------------------------------------------------

def test_f3_totals_classify():
    """DK exposes 'Total Runs - 1st 3 Innings' bundling 1.5 / 2.5 / 3.5."""
    assert classify_market("Total Runs - 1st 3 Innings") == ("F3", "main")


def test_f3_spreads_classify():
    assert classify_market("Run Line - 1st 3 Innings") == ("F3", "main")


# ---------------------------------------------------------------------------
# Team-specific markets — MUST all return None
#
# Regression: pre-fix, 'Alternate <TEAM> Total Runs' slipped through the
# 'team total' substring filter and collided with the real game-totals
# market in parse_selections_to_wide_rows's (period, market_type, line)
# bucket, overwriting game-total odds with one team's team-total odds.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    # Caught by the original 'team total' substring filter
    "HOU Astros: Team Total Runs",
    "CHI Cubs: Team Total Runs - 1st 5 Innings",
    "SD Padres: Team Total Runs - 1st 3 Innings",
    # NEW: caught by the team-prefix filter
    "Alternate HOU Astros Total Runs",
    "Alternate CHI Cubs Total Runs",
    "Alternate SD Padres Total Runs",
    "HOU Astros Extra Base Hits",
    "CHI Cubs Home Runs",
    "CHI Cubs Walks",
    "HOU Astros Walks",
    "ATL Braves Hits",
])
def test_team_specific_markets_filtered(name):
    assert classify_market(name) is None


# ---------------------------------------------------------------------------
# Single-inning markets (3rd Inning, 5th Inning, etc.) — still filtered.
# Note "1st N Innings" is plural and stays.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "Total Runs - 5th Inning",
    "Run Line - 1st Inning",
    "7th Inning (3 Way)",
    "3rd Inning Over/Under 1.5 Runs",
    "Runs - 2nd Inning",
])
def test_single_inning_markets_filtered(name):
    assert classify_market(name) is None


# ---------------------------------------------------------------------------
# Player props and other out-of-scope keywords
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "Nico Hoerner Hits O/U",
    "Alex Bregman RBIs O/U",
    "Player To Record 2+ Hits + Runs + RBIs",
    "Correct Score",
    "Race To 3 Runs",
    "Total Bases",
    "Odd/Even Total Runs",
])
def test_props_and_oos_filtered(name):
    assert classify_market(name) is None
