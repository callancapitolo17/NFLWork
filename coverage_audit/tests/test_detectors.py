from datetime import datetime, timezone, timedelta
from coverage_audit.detectors import Gap, detect_regressions, detect_freshness, detect_rowcount

NOW = datetime.now(timezone.utc)

def test_regression_flags_disappeared_market():
    latest = {"fetch_time": NOW, "markets": {("spreads","FG")}, "row_count": 10}
    baseline = {"run_count": 5, "market_run_fraction": {
        ("spreads","FG"): 1.0, ("alternate_totals","FG"): 1.0}, "median_row_count": 10}
    gaps = detect_regressions("fanduel_singles", latest, baseline)
    assert len(gaps) == 1
    g = gaps[0]
    assert (g.market, g.period) == ("alternate_totals","FG")
    assert g.gap_type == "regression" and g.severity == "alert"

def test_regression_ignores_rarely_seen_market():
    latest = {"fetch_time": NOW, "markets": set(), "row_count": 0}
    baseline = {"run_count": 10, "market_run_fraction": {("odd_even","FG"): 0.2}, "median_row_count": 1}
    assert detect_regressions("wagerzon", latest, baseline) == []

def test_freshness_flags_stale_and_missing():
    stale = {"fetch_time": NOW - timedelta(hours=30), "markets": set(), "row_count": 0}
    assert detect_freshness("bfa", stale)[0].gap_type == "freshness"
    assert detect_freshness("bfa", None)[0].severity == "alert"

def test_freshness_handles_naive_timestamp():
    # wagerzon_specials.scraped_at is a naive TIMESTAMP — must not raise on subtraction.
    naive_fresh = {"fetch_time": datetime.now().replace(tzinfo=None), "markets": set(), "row_count": 1}
    assert detect_freshness("wagerzon_specials", naive_fresh) == []   # fresh, no gap, no crash
    naive_stale = {"fetch_time": (datetime.now() - timedelta(hours=30)).replace(tzinfo=None),
                   "markets": set(), "row_count": 1}
    assert detect_freshness("wagerzon_specials", naive_stale)[0].gap_type == "freshness"

def test_rowcount_flags_drop():
    latest = {"fetch_time": NOW, "markets": set(), "row_count": 40}
    baseline = {"run_count": 7, "market_run_fraction": {}, "median_row_count": 600}
    g = detect_rowcount("fanduel_singles", latest, baseline)
    assert g[0].gap_type == "rowcount" and g[0].metric_value == 40
