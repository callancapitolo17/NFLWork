"""Regression tests for scraper_fanduel_singles.parse_runners_to_wide_rows.

Companion to test_fd_singles_parser.py. This file focuses on the
opposite-direction-alt-spread regression where abs(line) bucketing collapsed
KC -2.5/CHW +2.5 with KC +2.5/CHW -2.5 into the same row.
"""
from datetime import datetime

import pytest

from mlb_sgp.fd_client import Event, Runner
from mlb_sgp.scraper_fanduel_singles import classify_market, parse_runners_to_wide_rows


def test_parse_runners_alt_spread_emits_all_lines_signed_bucket():
    """FD alt-spread runners spanning both directions of the same magnitude
    (KC -2.5/+2.5 vs CHW -2.5/+2.5) must produce SEPARATE rows.

    Regression for the bug where abs(line) collapsed both directions into
    one bucket and the second pair silently overwrote the first."""
    event = Event(
        event_id="fd-test-altspread",
        away_team="Kansas City Royals",
        home_team="Chicago White Sox",
        start_time=datetime(2026, 5, 13, 19, 0),
    )
    # Two bet-pairs at magnitude 2.5: opposite favored directions.
    # FD encodes alt-spread lines in the runner NAME (line=None on the Runner)
    # so we mirror that here — the parser regexes the name to recover the line.
    runners = [
        # Pair A: KC favored at -2.5
        Runner(runner_id="rA1", market_id="mkt-alt",
               name="Kansas City Royals -2.5", line=None, american_odds=+190),
        Runner(runner_id="rA2", market_id="mkt-alt",
               name="Chicago White Sox +2.5", line=None, american_odds=-250),
        # Pair B: CHW favored at -2.5
        Runner(runner_id="rB1", market_id="mkt-alt",
               name="Kansas City Royals +2.5", line=None, american_odds=-620),
        Runner(runner_id="rB2", market_id="mkt-alt",
               name="Chicago White Sox -2.5", line=None, american_odds=+400),
    ]
    market_meta = {"mkt-alt": ("FG", "alternate_spreads")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                      fetch_time=datetime(2026, 5, 13))

    # Two distinct rows expected — one per signed home_team line.
    assert len(rows) == 2, f"expected 2 rows (signed dirs), got {len(rows)}"

    by_home_line = {r["home_spread"]: r for r in rows}
    assert sorted(by_home_line) == [-2.5, 2.5], \
        f"home_spread values should be -2.5 and 2.5, got {sorted(by_home_line)}"

    # Pair A: KC -2.5 / CHW +2.5  => home_spread = +2.5 (CHW is home, on +2.5)
    pair_a = by_home_line[2.5]
    assert pair_a["home_spread"] == 2.5
    assert pair_a["away_spread"] == -2.5
    assert pair_a["home_spread_price"] == -250    # CHW +2.5
    assert pair_a["away_spread_price"] == +190    # KC -2.5

    # Pair B: KC +2.5 / CHW -2.5  => home_spread = -2.5 (CHW is home, on -2.5)
    pair_b = by_home_line[-2.5]
    assert pair_b["home_spread"] == -2.5
    assert pair_b["away_spread"] == 2.5
    assert pair_b["home_spread_price"] == +400    # CHW -2.5
    assert pair_b["away_spread_price"] == -620    # KC +2.5


def test_classify_market_existing_whitelist_unchanged():
    """Sanity check that the signed-bucket fix didn't disturb classification."""
    assert classify_market("Run Line") == ("FG", "main")
    assert classify_market("Alternate Run Lines") == ("FG", "alternate_spreads")
    assert classify_market("First 5 Innings Alternate Run Lines") == ("F5", "alternate_spreads")
    assert classify_market("Random Garbage Market") is None


# (home, away) used for team-total exclusion cases
_H, _A = "San Francisco Giants", "Chicago White Sox"

@pytest.mark.parametrize("name,home,away,expected", [
    # --- ACCEPTED: real game lines, every period x type FD posts ---
    ("Run Line", None, None, ("FG", "main")),
    ("Total Runs", None, None, ("FG", "main")),
    ("Moneyline", None, None, ("FG", "main")),
    ("Alternate Run Lines", None, None, ("FG", "alternate_spreads")),
    ("Alternate Total Runs", None, None, ("FG", "alternate_totals")),
    ("First 5 Innings Run Line", None, None, ("F5", "main")),
    ("First 5 Innings Total Runs", None, None, ("F5", "main")),
    ("First 5 Innings Money Line", None, None, ("F5", "main")),  # NOTE the space
    ("First 5 Innings Alternate Run Lines", None, None, ("F5", "alternate_spreads")),
    ("First 5 Innings Alternate Total Runs", None, None, ("F5", "alternate_totals")),
    ("First 7 Innings Total Runs", None, None, ("F7", "main")),
    ("First 7 Innings Run Line", None, None, ("F7", "main")),
    ("First 3 Innings Total Runs", None, None, ("F3", "main")),
    ("First 3 Innings Run Line", None, None, ("F3", "main")),
    # --- REJECTED: one example per junk family ---
    ("First 5 Innings Run Line / Total Runs Parlay", None, None, None),  # parlay
    ("Line / Total Parlay 7", None, None, None),                         # parlay
    ("Total Runs (Bands)", None, None, None),                            # bands
    ("Moneyline Away Listed", None, None, None),                         # listed
    ("First 7 Innings Result", None, None, None),                        # 3-way result
    ("First 6 Innings Result", None, None, None),                        # 3-way result
    ("7th Inning Total Runs", None, None, None),                         # single inning
    ("7th Inning Run Line", None, None, None),                           # single inning
    ("First 5 Innings Winning Margin (5-Way)", None, None, None),        # winning margin
    ("Race To 7 Runs", None, None, None),                                # race to
    ("Tri-Bet", None, None, None),                                       # no line keyword
    ("Chicago White Sox Total Runs", _H, _A, None),                      # team total (away)
    ("San Francisco Giants Alt. Total Runs", _H, _A, None),              # team total (home)
    ("Random Garbage Market", None, None, None),
])
def test_classify_market_keyword(name, home, away, expected):
    from mlb_sgp.scraper_fanduel_singles import classify_market
    assert classify_market(name, home, away) == expected
