def latest_run(con, table):
    """Extract the latest run: fetch_time, markets, row_count."""
    row = con.execute(f"SELECT MAX(fetch_time) FROM {table}").fetchone()
    if row is None or row[0] is None:
        return None
    max_ft = row[0]
    pairs = con.execute(
        f"SELECT DISTINCT market, period FROM {table} WHERE fetch_time = ?", [max_ft]
    ).fetchall()
    count = con.execute(
        f"SELECT COUNT(*) FROM {table} WHERE fetch_time = ?", [max_ft]
    ).fetchone()[0]
    return {"fetch_time": max_ft, "markets": {(m, p) for m, p in pairs}, "row_count": count}


def trailing_baseline(con, table, days):
    """Compute baseline metrics over the last N days: run_count, market_run_fraction, median_row_count."""
    runs = con.execute(
        f"SELECT DISTINCT fetch_time FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days'"
    ).fetchall()
    run_times = [r[0] for r in runs]
    run_count = len(run_times)
    if run_count == 0:
        return {"run_count": 0, "market_run_fraction": {}, "median_row_count": None}
    # fraction of runs each (market,period) appears in
    rows = con.execute(
        f"SELECT market, period, COUNT(DISTINCT fetch_time) AS runs_present FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days' GROUP BY 1,2"
    ).fetchall()
    fraction = {(m, p): present / run_count for m, p, present in rows}
    median_rc = con.execute(
        f"SELECT median(rc) FROM (SELECT fetch_time, COUNT(*) AS rc FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days' GROUP BY fetch_time)"
    ).fetchone()[0]
    return {"run_count": run_count, "market_run_fraction": fraction, "median_row_count": median_rc}
