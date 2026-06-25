import duckdb
from datetime import datetime, timezone, timedelta
from coverage_audit.detectors import screen_bookmakers, detect_screen_absence
from coverage_audit.registry import get_books

def test_screen_bookmakers_window_includes_recent_excludes_stale(tmp_path):
    # Each book has its own fetch_time; presence is a recent-window check, not a
    # single MAX instant. Recent book counts; a book stale beyond the window does not.
    p = tmp_path / "mm.duckdb"; c = duckdb.connect(str(p))
    c.execute("CREATE TABLE mlb_bets_book_prices(bookmaker VARCHAR, fetch_time TIMESTAMPTZ)")
    now = datetime.now(timezone.utc)
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('fanduel', ?)", [now - timedelta(hours=1)])
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('wagerzon', ?)", [now - timedelta(hours=40)])
    c.close()
    con = duckdb.connect(str(p), read_only=True)
    assert screen_bookmakers(con, within_hours=26.0) == {"fanduel"}
    con.close()

def test_screen_bookmakers_returns_all_books_in_window(tmp_path):
    # Books written seconds apart in the same run must ALL count (the old MAX-instant
    # bug returned only the single newest book).
    p = tmp_path / "mm.duckdb"; c = duckdb.connect(str(p))
    c.execute("CREATE TABLE mlb_bets_book_prices(bookmaker VARCHAR, fetch_time TIMESTAMPTZ)")
    now = datetime.now(timezone.utc)
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('wagerzon', ?)", [now])
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('fanduel', ?)", [now - timedelta(seconds=10)])
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('bet105', ?)", [now - timedelta(seconds=11)])
    c.close()
    con = duckdb.connect(str(p), read_only=True)
    assert screen_bookmakers(con) == {"wagerzon", "fanduel", "bet105"}
    con.close()

def test_detect_screen_absence_flags_expected_missing():
    gaps = detect_screen_absence(get_books(), screen_labels={"wagerzon", "fanduel"})
    absent = {g.book for g in gaps}
    # bet105, bookmaker, draftkings expected but missing; hoop88/bfa/kalshi never expected
    assert "draftkings_singles" in absent and "bet105" in absent
    assert "hoop88" not in absent and "kalshi" not in absent
