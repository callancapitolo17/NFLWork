def ensure_schema(con):
    con.execute("""
        CREATE TABLE IF NOT EXISTS coverage_gaps(
            audit_ts TIMESTAMPTZ, book VARCHAR, gap_type VARCHAR, severity VARCHAR,
            market VARCHAR, period VARCHAR, detail VARCHAR,
            metric_value DOUBLE, baseline_value DOUBLE)
    """)

def gap_key(g):
    return (g.book, g.gap_type, g.market, g.period)

def previous_gap_keys(con):
    row = con.execute("SELECT MAX(audit_ts) FROM coverage_gaps").fetchone()
    if row is None or row[0] is None:
        return set()
    keys = con.execute(
        "SELECT book, gap_type, market, period FROM coverage_gaps WHERE audit_ts = ?", [row[0]]
    ).fetchall()
    return {(b, t, m, p) for b, t, m, p in keys}

def write_gaps(con, gaps, audit_ts):
    for g in gaps:
        con.execute(
            "INSERT INTO coverage_gaps VALUES (?,?,?,?,?,?,?,?,?)",
            [audit_ts, g.book, g.gap_type, g.severity, g.market, g.period,
             g.detail, g.metric_value, g.baseline_value])

def new_gap_keys(prev_keys, gaps):
    return {gap_key(g) for g in gaps} - set(prev_keys)
