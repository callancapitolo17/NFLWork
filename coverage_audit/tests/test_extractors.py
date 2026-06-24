from coverage_audit.detectors import latest_run, trailing_baseline
from coverage_audit.tests.conftest import days_ago


def test_latest_run_picks_max_fetch_time(book_db):
    con = book_db([
        (days_ago(2), [("spreads", "FG"), ("totals", "FG")]),
        (days_ago(0), [("spreads", "FG")]),
    ])
    lr = latest_run(con, "mlb_odds")
    assert lr["markets"] == {("spreads", "FG")}
    assert lr["row_count"] == 1


def test_latest_run_empty_table_returns_none(book_db):
    con = book_db([])
    assert latest_run(con, "mlb_odds") is None


def test_trailing_baseline_fractions(book_db):
    con = book_db([
        (days_ago(3), [("spreads", "FG"), ("totals", "FG")]),
        (days_ago(2), [("spreads", "FG"), ("totals", "FG")]),
        (days_ago(1), [("spreads", "FG")]),  # totals missing once
    ])
    b = trailing_baseline(con, "mlb_odds", days=7)
    assert b["run_count"] == 3
    assert b["market_run_fraction"][("spreads", "FG")] == 1.0
    assert abs(b["market_run_fraction"][("totals", "FG")] - 2/3) < 1e-9
