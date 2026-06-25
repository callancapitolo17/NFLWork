import duckdb
from datetime import datetime, timezone, timedelta
from coverage_audit.detectors import Gap
from coverage_audit.store import ensure_schema, write_gaps, previous_gap_keys, new_gap_keys, gap_key

def _con(tmp_path):
    return duckdb.connect(str(tmp_path / "cov.duckdb"))

def test_write_and_previous_keys(tmp_path):
    con = _con(tmp_path); ensure_schema(con)
    t1 = datetime.now(timezone.utc) - timedelta(days=1)
    write_gaps(con, [Gap("bfa","freshness","alert","stale")], t1)
    assert previous_gap_keys(con) == {("bfa","freshness",None,None)}

def test_previous_keys_uses_only_latest_audit(tmp_path):
    con = _con(tmp_path); ensure_schema(con)
    t1 = datetime.now(timezone.utc) - timedelta(days=2)
    t2 = datetime.now(timezone.utc) - timedelta(days=1)
    write_gaps(con, [Gap("bfa","freshness","alert","old")], t1)
    write_gaps(con, [Gap("hoop88","rowcount","warn","drop")], t2)
    assert previous_gap_keys(con) == {("hoop88","rowcount",None,None)}

def test_new_gap_keys_detects_only_fresh():
    prev = {("bfa","freshness",None,None)}
    gaps = [Gap("bfa","freshness","alert","still"), Gap("fanduel_singles","regression","alert","x",market="totals",period="FG")]
    assert new_gap_keys(prev, gaps) == {("fanduel_singles","regression","totals","FG")}
