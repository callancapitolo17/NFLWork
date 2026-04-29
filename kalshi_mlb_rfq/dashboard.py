"""CLI status tool for the Kalshi MLB RFQ bot."""

from kalshi_mlb_rfq import db


def main():
    print("=" * 60)
    print("  Kalshi MLB RFQ Bot — status")
    print("=" * 60)

    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT COUNT(*) FROM live_rfqs WHERE status='open'"
        ).fetchone()[0]
        today_fills = con.execute(
            "SELECT COUNT(*), COALESCE(SUM(contracts*price_dollars), 0) "
            "FROM fills WHERE filled_at >= CURRENT_DATE"
        ).fetchone()
        decisions = con.execute(
            "SELECT decision, COUNT(*) FROM quote_log "
            "WHERE observed_at >= NOW() - INTERVAL '1 HOUR' GROUP BY decision"
        ).fetchall()
        recent_session = con.execute(
            "SELECT session_id, started_at, ended_at, dry_run "
            "FROM sessions ORDER BY started_at DESC LIMIT 1"
        ).fetchone()

    if recent_session:
        sid, started, ended, dry = recent_session
        status = "running" if ended is None else "ended"
        print(f"  Last session   : {sid[:8]}  ({status}, dry_run={dry}, started {started})")
    print(f"  Live RFQs      : {live}")
    print(f"  Today fills    : {today_fills[0]} totaling ${float(today_fills[1]):.2f}")
    print()
    print("  Quote-log decisions in last hour:")
    if not decisions:
        print("    (none)")
    for d, n in decisions:
        print(f"    {d:<35} {n}")


if __name__ == "__main__":
    main()
