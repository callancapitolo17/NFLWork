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


def universe_stats(market_db: str, since_naive: datetime,
                   now_naive: datetime) -> dict:
    """Quotable universe: combos passing consensus per day + book coverage.

    A combo is (game_id, combo, spread_line, total_line) — GROUP BY treats
    NULLs as equal, so ML x total rows (spread_line NULL by design) group
    correctly. "Passing" = >= MIN_AGREEING_BOOKS distinct books that day
    (an approximation of the live consensus gate, which also applies the
    freshness window and the agreement band at price time).
    """
    out = {"per_day": [], "per_book": [], "min_books":
           config.MIN_AGREEING_BOOKS, "unavailable": None}

    per_day = _read(
        market_db,
        "WITH day_combos AS ( "
        "  SELECT CAST(fetch_time AS DATE) AS day, game_id, combo, "
        "         spread_line, total_line, "
        "         count(DISTINCT bookmaker) AS books "
        "  FROM mlb_sgp_odds WHERE fetch_time >= ? "
        "  GROUP BY ALL) "
        "SELECT day, count(*) FILTER (books >= ?), count(*) "
        "FROM day_combos GROUP BY day ORDER BY day",
        [since_naive, config.MIN_AGREEING_BOOKS])
    if not _usable(per_day):
        out["unavailable"] = _unavailable_note(per_day, "market")
        return out
    out["per_day"] = [tuple(r) for r in per_day]

    per_book = _read(
        market_db,
        "SELECT bookmaker, "
        "       count(DISTINCT CAST(fetch_time AS DATE)), "
        "       count(DISTINCT (game_id, combo, spread_line, total_line)), "
        "       CAST(epoch(? - max(fetch_time)) / 60 AS INTEGER) "
        "FROM mlb_sgp_odds WHERE fetch_time >= ? "
        "GROUP BY bookmaker ORDER BY bookmaker",
        [now_naive, since_naive])
    if _usable(per_book):
        out["per_book"] = [tuple(r) for r in per_book]
    return out


def _render_universe(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    if not stats["per_day"]:
        return "_no SGP odds rows in window._"
    lines = [
        f"Combos passing consensus (>= {stats['min_books']} books) per day:",
        "",
        "| day | passing | total scraped |",
        "|---|---:|---:|",
    ]
    for day, passing, total in stats["per_day"]:
        lines.append(f"| {day} | {passing} | {total} |")
    lines += [
        "",
        "**Per-book participation (7d):**",
        "",
        "| book | days present | distinct combos | last row age (min) |",
        "|---|---:|---:|---:|",
    ]
    for book, days, combos, age_min in stats["per_book"]:
        lines.append(f"| {book} | {days} | {combos} | {age_min} |")
    lines += [
        "",
        "_Passing = distinct-book count per day; the live gate also applies "
        "freshness + agreement-band checks at price time._",
    ]
    return "\n".join(lines)


def _age_dist(ages: list[float]) -> dict:
    return {"n": len(ages), "p50": _percentile(ages, 0.50),
            "p90": _percentile(ages, 0.90), "p95": _percentile(ages, 0.95),
            "max": max(ages) if ages else None}


def staleness_stats(research_db: str, since: datetime) -> dict:
    """Book-data age at quote time + quote age at accept time (seconds).

    Quote-time age is a proxy: per-book scrape ages aren't persisted at
    quote time, so we measure each quote_priced event against the latest
    prior scrape_done (the SGP cycle that produced the data it priced on).
    Accept-time age is exact: confirm_singles_check.quote_age_sec.
    """
    import bisect
    import json

    out = {"quote_age": _age_dist([]), "quote_age_unknown": 0,
           "accept_age": _age_dist([]), "unavailable": None}

    scrapes = _read(research_db,
                    "SELECT ts FROM events WHERE event_type = 'scrape_done' "
                    "ORDER BY ts")
    if not _usable(scrapes):
        out["unavailable"] = _unavailable_note(scrapes, "research")
        return out
    quotes = _read(research_db,
                   "SELECT ts FROM events WHERE event_type = 'quote_priced' "
                   "AND ts >= ?", [since])
    confirms = _read(research_db,
                     "SELECT payload FROM events "
                     "WHERE event_type = 'confirm_singles_check' "
                     "AND ts >= ?", [since])
    if not _usable(quotes) or not _usable(confirms):
        out["unavailable"] = _unavailable_note(LOCKED, "research")
        return out

    # scrape_done is deliberately unwindowed: a quote just inside the window
    # may price off a scrape just before it.
    scrape_ts = [r[0] for r in scrapes]
    quote_ages = []
    for (quote_ts,) in quotes:
        idx = bisect.bisect_right(scrape_ts, quote_ts)
        if idx == 0:
            out["quote_age_unknown"] += 1
        else:
            quote_ages.append((quote_ts - scrape_ts[idx - 1]).total_seconds())
    out["quote_age"] = _age_dist(quote_ages)

    accept_ages = []
    for (payload_str,) in confirms:
        age = json.loads(payload_str).get("quote_age_sec")
        if age is not None:
            accept_ages.append(float(age))
    out["accept_age"] = _age_dist(accept_ages)
    return out


def _render_age_row(label: str, dist: dict) -> str:
    return (f"| {label} | {dist['n']} | {_fmt(dist['p50'])} | "
            f"{_fmt(dist['p90'])} | {_fmt(dist['p95'])} | "
            f"{_fmt(dist['max'])} |")


def _render_staleness(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    if stats["quote_age"]["n"] == 0 and stats["accept_age"]["n"] == 0:
        return "_no quotes in window._"
    lines = [
        "| age (s) | n | p50 | p90 | p95 | max |",
        "|---|---:|---:|---:|---:|---:|",
        _render_age_row("book data at quote", stats["quote_age"]),
        _render_age_row("quote at accept", stats["accept_age"]),
        "",
        "_Quote-time age = gap to the latest prior scrape cycle (per-book "
        "ages aren't persisted at quote time)._",
    ]
    if stats["quote_age_unknown"]:
        lines.append(f"_{stats['quote_age_unknown']} quotes had no prior "
                     "scrape_done event (bot warmup)._")
    return "\n".join(lines)


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
        "## 2. Quotable universe",
        "",
        _render_universe(universe_stats(
            market_db, (now - timedelta(days=7)).replace(tzinfo=None),
            now.replace(tzinfo=None))),
        "",
        "## 3. Staleness (7d)",
        "",
        _render_staleness(staleness_stats(research_db, d7)),
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
