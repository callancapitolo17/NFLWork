from dataclasses import dataclass
from datetime import datetime, timezone


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


@dataclass
class Gap:
    book: str
    gap_type: str          # regression | freshness | rowcount | book_absent
    severity: str          # alert | warn
    detail: str
    market: str | None = None
    period: str | None = None
    metric_value: float | None = None
    baseline_value: float | None = None


def detect_regressions(book_name, latest, baseline, presence_threshold=0.8):
    if latest is None or baseline["run_count"] == 0:
        return []
    gaps = []
    for (m, p), frac in baseline["market_run_fraction"].items():
        if frac >= presence_threshold and (m, p) not in latest["markets"]:
            gaps.append(Gap(
                book=book_name, gap_type="regression", severity="alert",
                market=m, period=p, metric_value=0.0, baseline_value=frac,
                detail=f"{m}/{p} present in {frac:.0%} of last runs but absent today "
                       f"(possible parse break or auth expiry — check recon)"))
    return gaps


def detect_freshness(book_name, latest, max_age_hours=26.0):
    if latest is None:
        return [Gap(book=book_name, gap_type="freshness", severity="alert",
                    detail="no rows at all in book DB (scraper never wrote / DB unreadable)")]
    age_h = (datetime.now(timezone.utc) - latest["fetch_time"]).total_seconds() / 3600.0
    if age_h > max_age_hours:
        return [Gap(book=book_name, gap_type="freshness", severity="alert",
                    metric_value=round(age_h, 1), baseline_value=max_age_hours,
                    detail=f"latest run is {age_h:.1f}h old (> {max_age_hours}h threshold)")]
    return []


def detect_rowcount(book_name, latest, baseline, min_ratio=0.5):
    med = baseline.get("median_row_count")
    if latest is None or not med:
        return []
    if latest["row_count"] < min_ratio * med:
        return [Gap(book=book_name, gap_type="rowcount", severity="warn",
                    metric_value=float(latest["row_count"]), baseline_value=float(med),
                    detail=f"row count {latest['row_count']} < {min_ratio:.0%} of trailing "
                           f"median {med:.0f} (possible partial parse break)")]
    return []


def screen_bookmakers(con, table="mlb_bets_book_prices"):
    row = con.execute(f"SELECT MAX(fetch_time) FROM {table}").fetchone()
    if row is None or row[0] is None:
        return set()
    labels = con.execute(
        f"SELECT DISTINCT bookmaker FROM {table} WHERE fetch_time = ?", [row[0]]
    ).fetchall()
    return {l[0] for l in labels}


def detect_screen_absence(books, screen_labels):
    gaps = []
    for b in books:
        if b.expected_on_screen and b.screen_name not in screen_labels:
            gaps.append(Gap(
                book=b.name, gap_type="book_absent", severity="alert",
                detail=f"book '{b.screen_name}' expected on odds screen but absent "
                       f"from latest mlb_bets_book_prices run"))
    return gaps
