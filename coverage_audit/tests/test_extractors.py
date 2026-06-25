import duckdb
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


def _specials_db(tmp_path):
    """A wagerzon_specials-like table: ts col is `scraped_at`, no market/period."""
    p = tmp_path / "specials.duckdb"
    c = duckdb.connect(str(p))
    c.execute("CREATE TABLE wagerzon_specials(scraped_at TIMESTAMPTZ, description VARCHAR)")
    run_ts = days_ago(0)  # one run = one timestamp shared by both rows
    c.execute("INSERT INTO wagerzon_specials VALUES (?, 'TRIPLE-PLAY')", [run_ts])
    c.execute("INSERT INTO wagerzon_specials VALUES (?, 'GRAND-SLAM')", [run_ts])
    c.close()
    return duckdb.connect(str(p), read_only=True)


def test_latest_run_custom_ts_col_no_market_period(tmp_path):
    con = _specials_db(tmp_path)
    lr = latest_run(con, "wagerzon_specials", ts_col="scraped_at")
    assert lr["markets"] == set()      # no market/period columns → empty, no crash
    assert lr["row_count"] == 2


def test_trailing_baseline_custom_ts_col_no_market_period(tmp_path):
    con = _specials_db(tmp_path)
    b = trailing_baseline(con, "wagerzon_specials", days=7, ts_col="scraped_at")
    assert b["run_count"] == 1
    assert b["market_run_fraction"] == {}
    assert b["median_row_count"] == 2
