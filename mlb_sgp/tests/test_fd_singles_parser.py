"""Tests parse_runners_to_wide_rows — converts FD runners to mlb_odds rows.

Mirror of test_dk_singles_parser.py — same logic, FD's Runner dataclass instead
of DK's Selection. Verifies bucketing (one row for main, one per line for
alts, abs(line) coalesces both sides of an alt spread) and period isolation.
"""
from datetime import datetime
from mlb_sgp.fd_client import Event, Runner
from mlb_sgp.scraper_fanduel_singles import parse_runners_to_wide_rows


def test_main_game_lines_produce_one_row():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_main_spread", "Yankees -1.5", -1.5, -120),
        Runner("r2", "m_main_spread", "Red Sox +1.5", 1.5, 110),
        Runner("r3", "m_main_total", "Over 9.5", 9.5, -105),
        Runner("r4", "m_main_total", "Under 9.5", 9.5, -115),
        Runner("r5", "m_ml", "Yankees", None, -150),
        Runner("r6", "m_ml", "Red Sox", None, 130),
    ]
    market_meta = {
        "m_main_spread": ("FG", "main"),
        "m_main_total":  ("FG", "main"),
        "m_ml":          ("FG", "main"),
    }
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["home_spread"] == -1.5
    assert r["total"] == 9.5
    assert r["home_ml"] == -150


def test_alt_totals_emit_separate_rows():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_alt_total", "Over 8.5",  8.5, 120),
        Runner("r2", "m_alt_total", "Under 8.5", 8.5, -150),
        Runner("r3", "m_alt_total", "Over 10.5", 10.5, 180),
        Runner("r4", "m_alt_total", "Under 10.5", 10.5, -220),
    ]
    market_meta = {"m_alt_total": ("FG", "alternate_totals")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    totals = sorted([r["total"] for r in rows])
    assert totals == [8.5, 10.5]
    assert all(r["market"] == "alternate_totals" for r in rows)


def test_alt_spreads_coalesce_home_and_away():
    """Both sides of an alt-spread line share one row, keyed by abs(line)."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_alt_spread", "Yankees -2.5", -2.5, 120),
        Runner("r2", "m_alt_spread", "Red Sox +2.5", 2.5, -140),
    ]
    market_meta = {"m_alt_spread": ("FG", "alternate_spreads")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["home_spread"] == -2.5
    assert r["home_spread_price"] == 120
    assert r["away_spread"] == 2.5
    assert r["away_spread_price"] == -140


def test_fd_alt_spread_line_parsed_from_name():
    """FD alt-spread runners have line=None and the line embedded in the name
    (e.g. 'Athletics +3.5'). Parser must regex the name to recover the line —
    otherwise alt rows collapse to one bucket and overwrite each other.
    """
    event = Event("e1", "Athletics", "St. Louis Cardinals", "2026-05-12T22:00:00Z")
    runners = [
        # FD's actual runner format for "Alternate Run Lines":
        Runner("r1", "m_alt_spread", "St. Louis Cardinals +1.5", None, -162),
        Runner("r2", "m_alt_spread", "Athletics -1.5",          None,  130),
        Runner("r3", "m_alt_spread", "St. Louis Cardinals +2.5", None, -265),
        Runner("r4", "m_alt_spread", "Athletics -2.5",          None,  200),
    ]
    market_meta = {"m_alt_spread": ("FG", "alternate_spreads")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    # Sort by abs(home_spread) so we get deterministic order
    rows_sorted = sorted(rows, key=lambda r: abs(r["home_spread"] or 0))
    r1, r2 = rows_sorted
    assert r1["home_spread"] == -1.5 and r1["home_spread_price"] == 130
    assert r1["away_spread"] == 1.5 and r1["away_spread_price"] == -162
    assert r2["home_spread"] == -2.5 and r2["home_spread_price"] == 200
    assert r2["away_spread"] == 2.5 and r2["away_spread_price"] == -265
    # Critically: ML fields stay empty — alt-spread runners must NOT spill
    # into the moneyline branch even though their line was None on the
    # Runner dataclass.
    assert all(r["home_ml"] is None and r["away_ml"] is None for r in rows)


def test_fd_alt_total_line_parsed_from_name():
    """FD alt-total runners have line=None and the line in parens in the name
    (e.g. 'Over (8.5)'). Parser must regex the name and split by line.
    """
    event = Event("e1", "Athletics", "St. Louis Cardinals", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_alt_total", "Over (8.5)",  None, -194),
        Runner("r2", "m_alt_total", "Under (8.5)", None,  150),
        Runner("r3", "m_alt_total", "Over (10.5)", None,  114),
        Runner("r4", "m_alt_total", "Under (10.5)", None, -144),
    ]
    market_meta = {"m_alt_total": ("FG", "alternate_totals")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    totals = sorted([r["total"] for r in rows])
    assert totals == [8.5, 10.5]
    by_total = {r["total"]: r for r in rows}
    assert by_total[8.5]["over_price"] == -194
    assert by_total[8.5]["under_price"] == 150
    assert by_total[10.5]["over_price"] == 114
    assert by_total[10.5]["under_price"] == -144


def test_classify_market_whitelist_excludes_parlay_markets():
    """FD posts dozens of '<X> / Total Runs Parlay' / 'Line / Total Parlay N'
    markets that match naive 'total'/'run line' keywords. Whitelist ensures
    they're all skipped.
    """
    from mlb_sgp.scraper_fanduel_singles import classify_market
    # These should ALL return None — they look like main markets but aren't
    assert classify_market("Run Line / Total Runs Parlay") is None
    assert classify_market("Line / Total Parlay 9") is None
    assert classify_market("Money Line / Total Runs Parlay") is None
    assert classify_market("Home Run / Moneyline Parlay") is None
    # Team totals and bands
    assert classify_market("Athletics Total Runs") is None
    assert classify_market("St. Louis Cardinals Alt. Total Runs") is None
    assert classify_market("Total Runs (Bands)") is None
    # Moneyline-Listed variants (pitcher-conditional ML)
    assert classify_market("Moneyline Away Listed") is None
    assert classify_market("Moneyline Home Listed") is None
    # Single-inning markets and Race-To
    assert classify_market("5th Inning Run Line") is None
    assert classify_market("Race To 5 Runs") is None
    # First X Innings Result (3-way ML)
    assert classify_market("First 7 Innings Result") is None
    # Core in-scope names DO match
    assert classify_market("Run Line") == ("FG", "main")
    assert classify_market("Total Runs") == ("FG", "main")
    assert classify_market("Moneyline") == ("FG", "main")
    assert classify_market("Alternate Run Lines") == ("FG", "alternate_spreads")
    assert classify_market("Alternate Total Runs") == ("FG", "alternate_totals")
    assert classify_market("First 5 Innings Run Line") == ("F5", "main")
    assert classify_market("First 5 Innings Total Runs") == ("F5", "main")
    assert classify_market("First 5 Innings Money Line") == ("F5", "main")
    assert classify_market("First 5 Innings Alternate Run Lines") == ("F5", "alternate_spreads")
    assert classify_market("First 5 Innings Alternate Total Runs") == ("F5", "alternate_totals")


def test_f5_period_isolated_from_fg():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_fg_total",  "Over 9.5",  9.5,  -110),
        Runner("r2", "m_fg_total",  "Under 9.5", 9.5,  -110),
        Runner("r3", "m_f5_total",  "Over 4.5",  4.5,  -115),
        Runner("r4", "m_f5_total",  "Under 4.5", 4.5,  -105),
    ]
    market_meta = {
        "m_fg_total": ("FG", "main"),
        "m_f5_total": ("F5", "main"),
    }
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    periods = sorted([r["period"] for r in rows])
    assert periods == ["F5", "FG"]
