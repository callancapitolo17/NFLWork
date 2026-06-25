import duckdb
from datetime import datetime, timezone, timedelta
from coverage_audit.detectors import screen_bookmakers, detect_screen_absence
from coverage_audit.registry import get_books

def test_screen_bookmakers_at_latest_only(tmp_path):
    p = tmp_path / "mm.duckdb"; c = duckdb.connect(str(p))
    c.execute("CREATE TABLE mlb_bets_book_prices(bookmaker VARCHAR, fetch_time TIMESTAMPTZ)")
    now = datetime.now(timezone.utc)
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('wagerzon', ?)", [now - timedelta(days=1)])
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('fanduel', ?)", [now])
    c.close()
    con = duckdb.connect(str(p), read_only=True)
    assert screen_bookmakers(con) == {"fanduel"}
    con.close()

def test_detect_screen_absence_flags_expected_missing():
    gaps = detect_screen_absence(get_books(), screen_labels={"wagerzon", "fanduel"})
    absent = {g.book for g in gaps}
    # bet105, bookmaker, draftkings expected but missing; hoop88/bfa/kalshi never expected
    assert "draftkings_singles" in absent and "bet105" in absent
    assert "hoop88" not in absent and "kalshi" not in absent
