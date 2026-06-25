import json, subprocess
from datetime import datetime, timezone
import coverage_audit.audit as audit

def test_notify_disabled_is_noop(monkeypatch):
    called = {"n": 0}
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.__setitem__("n", called["n"]+1))
    audit.notify("t", "m", enabled=False)
    assert called["n"] == 0

def test_notify_swallows_errors(monkeypatch):
    def boom(*a, **k): raise OSError("no osascript")
    monkeypatch.setattr(subprocess, "run", boom)
    audit.notify("t", "m", enabled=True)  # must not raise

def test_run_audit_returns_summary_and_writes(monkeypatch, tmp_path):
    # Point coverage DB at tmp and stub the per-book + screen scans to fixed gaps.
    import coverage_audit.db as dbmod
    monkeypatch.setattr(dbmod, "COVERAGE_DB_PATH", tmp_path / "cov.duckdb")
    from coverage_audit.detectors import Gap
    monkeypatch.setattr(audit, "scan_books", lambda: [Gap("bfa","freshness","alert","stale")])
    monkeypatch.setattr(audit, "scan_screen", lambda: [])
    monkeypatch.setattr(audit, "notify", lambda *a, **k: None)
    ts = datetime.now(timezone.utc)
    summary = audit.run_audit(audit_ts=ts, notify_enabled=False)
    assert summary["total_gaps"] == 1
    assert summary["new_gaps"] == 1
    assert json.dumps(summary)  # serializable


def test_scan_books_does_not_crash_on_query_error(monkeypatch, tmp_path):
    # A readable DB whose table lacks the expected columns must yield a gap, not raise.
    import duckdb
    from coverage_audit.registry import Book
    bad = tmp_path / "bad.duckdb"
    c = duckdb.connect(str(bad)); c.execute("CREATE TABLE mlb_odds(wrong_col INT)"); c.close()
    book = Book("bookbad", bad, "mlb_odds", None, False, None)  # queries fetch_time → BinderException
    monkeypatch.setattr(audit, "get_books", lambda: [book])
    gaps = audit.scan_books()  # must not raise
    assert len(gaps) == 1
    assert gaps[0].book == "bookbad" and gaps[0].gap_type == "freshness"
    assert "detection failed" in gaps[0].detail


def test_scan_screen_does_not_crash_on_bad_schema(monkeypatch, tmp_path):
    # Screen DB opens but the expected table is missing → gap, not a raised exception.
    import duckdb
    bad = tmp_path / "mm_bad.duckdb"
    c = duckdb.connect(str(bad)); c.execute("CREATE TABLE unrelated(x INT)"); c.close()
    monkeypatch.setattr(audit, "SCREEN_DB", bad)
    gaps = audit.scan_screen()  # must not raise
    assert len(gaps) == 1
    assert gaps[0].book == "(screen)" and gaps[0].gap_type == "book_absent"
    assert "unreadable" in gaps[0].detail
