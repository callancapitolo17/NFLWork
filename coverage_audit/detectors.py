from dataclasses import dataclass
from datetime import datetime, timezone


def _table_columns(con, table):
    """Set of column names in a table (lowercased), via PRAGMA table_info."""
    return {r[1].lower() for r in con.execute(f"PRAGMA table_info('{table}')").fetchall()}


def latest_run(con, table, ts_col="fetch_time"):
    """Extract the latest run: fetch_time, markets, row_count.

    `ts_col` is the per-book timestamp column. If the table lacks market/period
    columns (e.g. wagerzon_specials), `markets` is an empty set and only
    freshness/row-count are meaningful — market-regression simply finds nothing.
    """
    row = con.execute(f"SELECT MAX({ts_col}) FROM {table}").fetchone()
    if row is None or row[0] is None:
        return None
    max_ft = row[0]
    cols = _table_columns(con, table)
    if "market" in cols and "period" in cols:
        pairs = con.execute(
            f"SELECT DISTINCT market, period FROM {table} WHERE {ts_col} = ?", [max_ft]
        ).fetchall()
        markets = {(m, p) for m, p in pairs}
    else:
        markets = set()
    count = con.execute(
        f"SELECT COUNT(*) FROM {table} WHERE {ts_col} = ?", [max_ft]
    ).fetchone()[0]
    return {"fetch_time": max_ft, "markets": markets, "row_count": count}


def trailing_baseline(con, table, days, ts_col="fetch_time"):
    """Compute baseline metrics over the last N days: run_count, market_run_fraction, median_row_count.

    `market_run_fraction` is empty when the table has no market/period columns.
    """
    runs = con.execute(
        f"SELECT DISTINCT {ts_col} FROM {table} "
        f"WHERE {ts_col} > now() - INTERVAL '{int(days)} days'"
    ).fetchall()
    run_count = len(runs)
    if run_count == 0:
        return {"run_count": 0, "market_run_fraction": {}, "median_row_count": None}
    cols = _table_columns(con, table)
    if "market" in cols and "period" in cols:
        # fraction of runs each (market,period) appears in
        rows = con.execute(
            f"SELECT market, period, COUNT(DISTINCT {ts_col}) AS runs_present FROM {table} "
            f"WHERE {ts_col} > now() - INTERVAL '{int(days)} days' GROUP BY 1,2"
        ).fetchall()
        fraction = {(m, p): present / run_count for m, p, present in rows}
    else:
        fraction = {}
    median_rc = con.execute(
        f"SELECT median(rc) FROM (SELECT {ts_col}, COUNT(*) AS rc FROM {table} "
        f"WHERE {ts_col} > now() - INTERVAL '{int(days)} days' GROUP BY {ts_col})"
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
    ft = latest["fetch_time"]
    if ft.tzinfo is None:
        # Naive timestamps (e.g. wagerzon_specials.scraped_at is plain TIMESTAMP) → assume UTC.
        ft = ft.replace(tzinfo=timezone.utc)
    age_h = (datetime.now(timezone.utc) - ft).total_seconds() / 3600.0
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


def screen_bookmakers(con, table="mlb_bets_book_prices", within_hours=26.0):
    """Bookmakers present on the odds screen within the recent window.

    Each book's rows carry their OWN scraper fetch_time (they are not written at
    one shared instant), so a single MAX(fetch_time) filter would return only the
    one book with the newest row. A window-based check is the correct "is this
    book currently on the screen" test.
    """
    labels = con.execute(
        f"SELECT DISTINCT bookmaker FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(within_hours)} hours'"
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
