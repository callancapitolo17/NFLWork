import argparse, json, subprocess, sys
from dataclasses import asdict
from datetime import datetime, timezone
from coverage_audit.registry import get_books, REPO_ROOT
from coverage_audit.db import connect_readonly, connect_coverage
from coverage_audit import detectors as D
from coverage_audit.store import ensure_schema, previous_gap_keys, write_gaps, new_gap_keys

SCREEN_DB = (REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb").resolve()
BASELINE_DAYS = 7

def notify(title, message, enabled=True):
    if not enabled:
        return
    try:
        subprocess.run(
            ["osascript", "-e", f'display notification "{message}" with title "{title}"'],
            check=False, capture_output=True, timeout=10)
    except Exception:
        pass

def scan_books():
    gaps = []
    for b in get_books():
        con = connect_readonly(b.db_path)
        if con is None:
            gaps.append(D.Gap(book=b.name, gap_type="freshness", severity="alert",
                              detail=f"DB unreadable or missing at {b.db_path}"))
            continue
        try:
            latest = D.latest_run(con, b.table, b.ts_col)
            baseline = D.trailing_baseline(con, b.table, BASELINE_DAYS, b.ts_col)
            gaps += D.detect_freshness(b.name, latest)
            gaps += D.detect_regressions(b.name, latest, baseline)
            gaps += D.detect_rowcount(b.name, latest, baseline)
        except Exception as e:
            # Never let one book's query error abort the whole audit (global constraint).
            gaps.append(D.Gap(book=b.name, gap_type="freshness", severity="alert",
                              detail=f"detection failed for table '{b.table}': {e}"))
        finally:
            con.close()
    return gaps

def scan_screen():
    con = connect_readonly(SCREEN_DB)
    if con is None:
        return [D.Gap(book="(screen)", gap_type="book_absent", severity="alert",
                      detail=f"screen DB unreadable at {SCREEN_DB}")]
    try:
        labels = D.screen_bookmakers(con)
    except Exception as e:
        # Readable DB but missing/renamed table → degrade to a gap, never abort the audit.
        return [D.Gap(book="(screen)", gap_type="book_absent", severity="alert",
                      detail=f"screen table unreadable in {SCREEN_DB}: {e}")]
    finally:
        con.close()
    return D.detect_screen_absence(get_books(), labels)

def run_audit(audit_ts=None, notify_enabled=True):
    audit_ts = audit_ts or datetime.now(timezone.utc)
    gaps = scan_books() + scan_screen()
    con = connect_coverage()
    try:
        ensure_schema(con)
        prev = previous_gap_keys(con)
        write_gaps(con, gaps, audit_ts)
    finally:
        con.close()
    fresh = new_gap_keys(prev, gaps)
    by_sev = {}
    for g in gaps:
        by_sev[g.severity] = by_sev.get(g.severity, 0) + 1
    summary = {
        "audit_ts": audit_ts.isoformat(),
        "total_gaps": len(gaps),
        "new_gaps": len(fresh),
        "by_severity": by_sev,
        "gaps": [asdict(g) for g in gaps],
    }
    if fresh:
        notify("MLB Coverage Audit",
               f"{len(fresh)} new gap(s), {len(gaps)} total", enabled=notify_enabled)
    return summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-notify", action="store_true")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args()
    summary = run_audit(notify_enabled=not args.no_notify)
    print(json.dumps(summary, indent=None if args.json_only else 2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
