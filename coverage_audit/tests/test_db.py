from pathlib import Path
import duckdb
from coverage_audit.db import connect_readonly, connect_coverage

def test_connect_readonly_missing_file_returns_none(tmp_path):
    assert connect_readonly(tmp_path / "nope.duckdb") is None

def test_connect_readonly_opens_existing(tmp_path):
    p = tmp_path / "x.duckdb"
    c = duckdb.connect(str(p)); c.execute("CREATE TABLE t(a INT)"); c.close()
    con = connect_readonly(p)
    assert con is not None
    assert con.execute("SELECT count(*) FROM t").fetchone()[0] == 0
    con.close()

def test_connect_coverage_creates_db(tmp_path, monkeypatch):
    import coverage_audit.db as dbmod
    monkeypatch.setattr(dbmod, "COVERAGE_DB_PATH", tmp_path / "coverage.duckdb")
    con = connect_coverage()
    assert (tmp_path / "coverage.duckdb").exists()
    con.close()
