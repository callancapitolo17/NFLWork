"""Daily state-of-the-maker report (issue #14).

Aggregates the maker bot's three DuckDBs into a printed markdown report:
RFQ funnel + on-demand coverage, quotable universe, staleness, demand curve,
settlement P&L, health.

Inputs:  kalshi_mlb_mm.duckdb (state), kalshi_mlb_mm_research.duckdb
         (firehose), kalshi_mlb_mm_market.duckdb (SGP odds) — all opened
         READ-ONLY with the monitor's lock-retry pattern.
Outputs: markdown on stdout.
Side effects: NONE. This module never writes any database.

Usage:
    python -m kalshi_mlb_mm.report              # live DBs (run from repo root)
    python -m kalshi_mlb_mm.report --state-db X --research-db Y --market-db Z
"""

import argparse
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb

from kalshi_mlb_mm import config


class _Sentinel:
    """Falsy marker so `if rows` skips both LOCKED and MISSING naturally."""

    def __init__(self, name: str):
        self._name = name

    def __bool__(self) -> bool:
        return False

    def __repr__(self) -> str:
        return self._name


LOCKED = _Sentinel("LOCKED")
MISSING = _Sentinel("MISSING")

_RETRY_SLEEP_SEC = 0.06  # checkpoint locks clear fast (monitor pattern)


def _read(db_path: str, sql: str, params=None, retries: int = 3):
    """Run one read-only query → list of row tuples, or LOCKED / MISSING.

    Missing file → MISSING (fresh install / bot never ran here).
    Missing table (fresh DB) → [] — reads as "no data yet", not an error.
    Lock/IO errors → retry briefly (the live bot holds sub-second checkpoint
    locks), then LOCKED — callers must render a note, never fake emptiness.
    """
    if not Path(db_path).exists():
        return MISSING
    for attempt in range(retries):
        con = None
        try:
            con = duckdb.connect(db_path, read_only=True)
            return con.execute(sql, params or []).fetchall()
        except duckdb.CatalogException:
            return []
        except Exception:
            if attempt < retries - 1:
                time.sleep(_RETRY_SLEEP_SEC)
                continue
            return LOCKED
        finally:
            if con is not None:
                try:
                    con.close()
                except Exception:
                    pass


def _unavailable_note(rows, db_label: str) -> str | None:
    """Markdown note for a failed read, or None if rows are usable."""
    if rows is MISSING:
        return f"_{db_label} DB not found — has the bot run on this machine?_"
    if rows is LOCKED:
        return f"_{db_label} DB locked (bot busy) — re-run in a few seconds._"
    return None


QUOTED_DECISIONS = ("quoted", "dry_run_quote")
# Post-accept outcomes: an accept happened iff one of these was logged.
ACCEPT_DECISIONS = ("confirmed", "voided_no_legs", "voided_singles_moved",
                    "voided_no_fresh_books", "voided_last_look")


def _sql_in(values: tuple) -> str:
    """Render a constant IN-list — duckdb params don't expand tuples."""
    return "(" + ", ".join(f"'{v}'" for v in values) + ")"


def _usable(rows) -> bool:
    return rows is not LOCKED and rows is not MISSING


# ---------------------------------------------------------------------------
# Sections (each returns a markdown fragment; none may raise)
# ---------------------------------------------------------------------------

def funnel_counts(state_db: str, since: datetime) -> dict:
    """RFQ funnel + path split for one window. Zeros on empty/unreadable DB.

    Path classification is rfq-level: out-of-scope (never priceable),
    on-demand (hit the Phase-2 live-fetch path at least once), else grid.
    Approximation: an on-demand RFQ whose fetch lands within its first 2s
    discovery tick never logs `on_demand_pending` and counts as grid.
    """
    out = {
        "seen": 0, "in_scope": 0, "quoted": 0, "accepted": 0,
        "confirmed": 0, "filled": 0, "decisions": [],
        "paths": {p: {"rfqs": 0, "quoted": 0}
                  for p in ("out_of_scope", "on_demand", "grid")},
        "unavailable": None,
    }

    seen = _read(state_db,
                 "SELECT count(*), count(*) FILTER (in_scope) "
                 "FROM seen_rfqs WHERE first_seen_at >= ?", [since])
    if not _usable(seen):
        out["unavailable"] = _unavailable_note(seen, "state")
        return out
    if seen:
        out["seen"], out["in_scope"] = seen[0]

    stages = _read(
        state_db,
        f"SELECT count(DISTINCT rfq_id) FILTER (decision IN {_sql_in(QUOTED_DECISIONS)}), "
        f"       count(DISTINCT rfq_id) FILTER (decision IN {_sql_in(ACCEPT_DECISIONS)}), "
        f"       count(DISTINCT rfq_id) FILTER (decision = 'confirmed') "
        "FROM quote_decisions WHERE observed_at >= ?", [since])
    if _usable(stages) and stages:
        out["quoted"], out["accepted"], out["confirmed"] = stages[0]

    filled = _read(state_db,
                   "SELECT count(*) FROM fills "
                   "WHERE contracts > 0 AND filled_at >= ?", [since])
    if _usable(filled) and filled:
        out["filled"] = filled[0][0]

    decisions = _read(state_db,
                      "SELECT decision, coalesce(reason, ''), count(*) "
                      "FROM quote_decisions WHERE observed_at >= ? "
                      "GROUP BY 1, 2 ORDER BY 3 DESC, 1, 2", [since])
    if _usable(decisions):
        out["decisions"] = [tuple(r) for r in decisions]

    paths = _read(
        state_db,
        "WITH win AS (SELECT rfq_id, in_scope FROM seen_rfqs "
        "             WHERE first_seen_at >= ?), "
        "od AS (SELECT DISTINCT rfq_id FROM quote_decisions "
        "       WHERE reason = 'on_demand_pending' AND observed_at >= ?), "
        f"q AS (SELECT DISTINCT rfq_id FROM quote_decisions "
        f"      WHERE decision IN {_sql_in(QUOTED_DECISIONS)} "
        "       AND observed_at >= ?) "
        "SELECT CASE WHEN NOT w.in_scope THEN 'out_of_scope' "
        "            WHEN od.rfq_id IS NOT NULL THEN 'on_demand' "
        "            ELSE 'grid' END AS path, "
        "       count(*), count(q.rfq_id) "
        "FROM win w "
        "LEFT JOIN od ON w.rfq_id = od.rfq_id "
        "LEFT JOIN q ON w.rfq_id = q.rfq_id "
        "GROUP BY 1", [since, since, since])
    if _usable(paths):
        for path, rfqs, quoted in paths:
            out["paths"][path] = {"rfqs": rfqs, "quoted": quoted}
    return out


def _percentile(values: list[float], q: float) -> float | None:
    """Linear-interpolated percentile (q in [0,1]); None on empty input."""
    if not values:
        return None
    s = sorted(values)
    pos = (len(s) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def on_demand_stats(research_db: str, since: datetime) -> dict:
    """Phase-2 on-demand pricing coverage: fetch success rate + latency.

    A fetch flight is `on_demand_requested`; a landing is `on_demand_result`
    with the same leg_set_hash. A hash with a request but no result = the
    whole fetch failed (all 6 books errored/timed out — main.py re-enqueues).
    Per-book latency_sec comes from the result payload's books dict; per-book
    failures show up as absence (book missing from the dict).
    """
    import json

    out = {"flights_requested": 0, "hashes_requested": 0, "hashes_landed": 0,
           "success_rate": None, "wall_p50": None, "wall_p95": None,
           "per_book": [], "unavailable": None}

    requests = _read(research_db,
                     "SELECT json_extract_string(payload, '$.leg_set_hash'), ts "
                     "FROM events WHERE event_type = 'on_demand_requested' "
                     "AND ts >= ?", [since])
    if not _usable(requests):
        out["unavailable"] = _unavailable_note(requests, "research")
        return out
    results = _read(research_db,
                    "SELECT payload, ts FROM events "
                    "WHERE event_type = 'on_demand_result' AND ts >= ?",
                    [since])
    if not _usable(results):
        out["unavailable"] = _unavailable_note(results, "research")
        return out

    out["flights_requested"] = len(requests)
    req_ts_by_hash: dict[str, list] = {}
    for leg_hash, ts in requests:
        req_ts_by_hash.setdefault(leg_hash, []).append(ts)
    out["hashes_requested"] = len(req_ts_by_hash)

    landed_hashes = set()
    wall_latencies = []
    book_latencies: dict[str, list[float]] = {}
    for payload_str, result_ts in results:
        payload = json.loads(payload_str)
        leg_hash = payload.get("leg_set_hash")
        landed_hashes.add(leg_hash)
        prior = [t for t in req_ts_by_hash.get(leg_hash, ())
                 if t <= result_ts]
        if prior:
            wall_latencies.append((result_ts - max(prior)).total_seconds())
        for book, info in (payload.get("books") or {}).items():
            latency = info.get("latency_sec")
            if latency is not None:
                book_latencies.setdefault(book, []).append(float(latency))

    out["hashes_landed"] = len(landed_hashes & set(req_ts_by_hash))
    if out["hashes_requested"]:
        out["success_rate"] = out["hashes_landed"] / out["hashes_requested"]
    out["wall_p50"] = _percentile(wall_latencies, 0.50)
    out["wall_p95"] = _percentile(wall_latencies, 0.95)
    out["per_book"] = sorted(
        (book, len(lats), _percentile(lats, 0.50), _percentile(lats, 0.95))
        for book, lats in book_latencies.items())
    return out


def _fmt(value, spec: str = ".1f") -> str:
    return "—" if value is None else format(value, spec)


def _render_on_demand(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    if stats["flights_requested"] == 0:
        return "_no on-demand fetches in window._"
    lines = [
        f"- fetch flights: {stats['flights_requested']} "
        f"({stats['hashes_requested']} distinct leg sets)",
        f"- landed: {stats['hashes_landed']} / {stats['hashes_requested']} "
        f"(success rate {_fmt(stats['success_rate'], '.0%')})",
        f"- wall latency request→result: p50 {_fmt(stats['wall_p50'])}s, "
        f"p95 {_fmt(stats['wall_p95'])}s",
        "",
        "| book | landings | fetch p50 (s) | fetch p95 (s) |",
        "|---|---:|---:|---:|",
    ]
    for book, n, p50, p95 in stats["per_book"]:
        lines.append(f"| {book} | {n} | {_fmt(p50)} | {_fmt(p95)} |")
    lines += [
        "",
        "_Confirm-lane re-fetch latency is not persisted; its failures appear "
        "as `voided_no_fresh_books` in §6._",
    ]
    return "\n".join(lines)


def _render_funnel(d24: dict, d7: dict) -> str:
    if d24["unavailable"]:
        return d24["unavailable"]
    lines = [
        "| stage | 24h | 7d |",
        "|---|---:|---:|",
    ]
    for stage in ("seen", "in_scope", "quoted", "accepted", "confirmed",
                  "filled"):
        lines.append(f"| {stage} | {d24[stage]} | {d7[stage]} |")
    lines += [
        "",
        "**Path split (rfq-level, 7d)** — share of flow by pricing path:",
        "",
        "| path | rfqs | of which quoted |",
        "|---|---:|---:|",
    ]
    total = sum(p["rfqs"] for p in d7["paths"].values()) or 1
    for path in ("grid", "on_demand", "out_of_scope"):
        p = d7["paths"][path]
        share = 100.0 * p["rfqs"] / total
        lines.append(f"| {path} | {p['rfqs']} ({share:.0f}%) | {p['quoted']} |")
    lines += [
        "",
        "_Approximation: an on-demand RFQ whose fetch lands within its first "
        "discovery tick never logs `on_demand_pending` and counts as grid._",
        "",
        "**Decision / reason breakdown (7d):**",
        "",
        "| decision | reason | n |",
        "|---|---|---:|",
    ]
    for decision, reason, n in d7["decisions"]:
        lines.append(f"| {decision} | {reason or '—'} | {n} |")
    if not d7["decisions"]:
        lines.append("| _none_ | | |")
    return "\n".join(lines)

def _section_pnl(state_db: str) -> str:
    rows = _read(state_db,
                 "SELECT count(*) FROM fills WHERE contracts > 0")
    note = _unavailable_note(rows, "state")
    if note:
        return note
    fill_count = rows[0][0] if rows else 0
    if fill_count == 0:
        return "no fills yet."
    return f"{fill_count} fills recorded."


def build_report(state_db: str, research_db: str, market_db: str,
                 now: datetime | None = None) -> str:
    """Assemble the full markdown report. Never raises on empty/missing DBs."""
    now = now or datetime.now(timezone.utc)
    lines = [
        f"# State of the Maker — {now.strftime('%Y-%m-%d %H:%M UTC')}",
        "",
    ]
    for label, path in (("state", state_db), ("research", research_db),
                        ("market", market_db)):
        if not Path(path).exists():
            lines.append(f"- `{label}` DB not found: `{path}`")
    h24 = now - timedelta(hours=24)
    d7 = now - timedelta(days=7)
    lines += [
        "",
        "## 1. RFQ funnel",
        "",
        _render_funnel(funnel_counts(state_db, h24),
                       funnel_counts(state_db, d7)),
        "",
        "### 1b. On-demand pricing coverage (7d)",
        "",
        _render_on_demand(on_demand_stats(research_db, d7)),
        "",
        "## 5. Settlement P&L",
        "",
        _section_pnl(state_db),
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-db", default=str(config.DB_PATH))
    parser.add_argument("--research-db", default=str(config.RESEARCH_DB_PATH))
    parser.add_argument("--market-db", default=str(config.MARKET_DB))
    args = parser.parse_args(argv)
    print(build_report(args.state_db, args.research_db, args.market_db))


if __name__ == "__main__":
    main()
