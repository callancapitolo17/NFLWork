"""Tests parse_selections_to_wide_rows — converts DK selections to mlb_odds rows."""
from datetime import datetime
from mlb_sgp.dk_client import Event, Selection
from mlb_sgp.scraper_draftkings_singles import parse_selections_to_wide_rows


def test_main_game_lines_produce_one_row():
    """FG main spread + total + ML for one game → one row marked 'main'."""
    event = Event(event_id="e1", home_team="Yankees", away_team="Red Sox",
                  start_time="2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_main_spread", "Yankees -1.5", -1.5, -120),
        Selection("s2", "m_main_spread", "Red Sox +1.5", 1.5, 110),
        Selection("s3", "m_main_total", "Over 9.5", 9.5, -105),
        Selection("s4", "m_main_total", "Under 9.5", 9.5, -115),
        Selection("s5", "m_ml", "Yankees", None, -150),
        Selection("s6", "m_ml", "Red Sox", None, 130),
    ]
    market_meta = {
        "m_main_spread": ("FG", "main"),
        "m_main_total":  ("FG", "main"),
        "m_ml":          ("FG", "main"),
    }
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["period"] == "FG"
    assert r["market"] == "main"
    assert r["home_spread"] == -1.5
    assert r["home_spread_price"] == -120
    assert r["away_spread"] == 1.5
    assert r["away_spread_price"] == 110
    assert r["total"] == 9.5
    assert r["over_price"] == -105
    assert r["under_price"] == -115
    assert r["home_ml"] == -150
    assert r["away_ml"] == 130


def test_alternate_spreads_emit_separate_rows():
    """Each alt-spread line gets its own row marked 'alternate_spreads'."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_alt_spread", "Yankees -2.5", -2.5, 120),
        Selection("s2", "m_alt_spread", "Red Sox +2.5", 2.5, -140),
        Selection("s3", "m_alt_spread", "Yankees -0.5", -0.5, -250),
        Selection("s4", "m_alt_spread", "Red Sox +0.5", 0.5, 200),
    ]
    market_meta = {"m_alt_spread": ("FG", "alternate_spreads")}
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    lines = sorted([r["home_spread"] for r in rows])
    assert lines == [-2.5, -0.5]
    assert all(r["market"] == "alternate_spreads" for r in rows)
    assert all(r["period"] == "FG" for r in rows)


def test_alternate_totals_emit_separate_rows():
    """Each alt-total line gets its own row marked 'alternate_totals'."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_alt_total", "Over 8.5",  8.5,  120),
        Selection("s2", "m_alt_total", "Under 8.5", 8.5,  -150),
        Selection("s3", "m_alt_total", "Over 10.5", 10.5, 180),
        Selection("s4", "m_alt_total", "Under 10.5", 10.5, -220),
    ]
    market_meta = {"m_alt_total": ("FG", "alternate_totals")}
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    totals = sorted([r["total"] for r in rows])
    assert totals == [8.5, 10.5]
    assert all(r["market"] == "alternate_totals" for r in rows)


def test_f5_period_isolated_from_fg():
    """F5 markets live on different rows than FG."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_fg_total",  "Over 9.5",  9.5,  -110),
        Selection("s2", "m_fg_total",  "Under 9.5", 9.5,  -110),
        Selection("s3", "m_f5_total",  "Over 4.5",  4.5,  -115),
        Selection("s4", "m_f5_total",  "Under 4.5", 4.5,  -105),
    ]
    market_meta = {
        "m_fg_total": ("FG", "main"),
        "m_f5_total": ("F5", "main"),
    }
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    periods = sorted([r["period"] for r in rows])
    assert periods == ["F5", "FG"]


def test_selection_with_unknown_market_skipped():
    """Selections whose market_id isn't in market_meta are silently skipped."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_main", "Over 9.5", 9.5, -110),
        Selection("s2", "m_main", "Under 9.5", 9.5, -110),
        Selection("s3", "m_unknown_prop", "Player X 2+ Hits", 2.5, 150),
    ]
    market_meta = {"m_main": ("FG", "main")}  # m_unknown_prop deliberately absent
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1  # only the main row, no row from the unknown market
