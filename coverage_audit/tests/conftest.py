import duckdb
import pytest
from datetime import datetime, timezone, timedelta


@pytest.fixture
def book_db(tmp_path):
    """Build a tiny mlb_odds DB. runs = list of (fetch_time, [(market,period),...])."""
    def _build(runs, table="mlb_odds"):
        p = tmp_path / "book.duckdb"
        c = duckdb.connect(str(p))
        c.execute(f"CREATE TABLE {table}(fetch_time TIMESTAMPTZ, market VARCHAR, period VARCHAR)")
        for ft, pairs in runs:
            for m, per in pairs:
                c.execute(f"INSERT INTO {table} VALUES (?, ?, ?)", [ft, m, per])
        c.close()
        return duckdb.connect(str(p), read_only=True)
    return _build


def days_ago(n):
    return datetime.now(timezone.utc) - timedelta(days=n)
