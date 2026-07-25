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
    """RFQ funnel + path split for one window, derived from quote_decisions.

    quote_decisions is the only complete record of flow: seen_rfqs only
    stores out-of-scope RFQs and LIVE-quoted ones (main.py:1124 writes it
    inside the quote transaction), so in-scope-but-skipped and dry-run RFQs
    never appear there. "seen" = RFQs with >=1 decision in the window;
    an out-of-scope RFQ logs exactly one decision, at first sight.

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

    # One row per RFQ with its lifecycle flags; halted_high_void_rate rows
    # carry rfq_id NULL and are deliberately excluded.
    per_rfq = _read(
        state_db,
        "WITH per_rfq AS ( "
        "  SELECT rfq_id, "
        "         max(CASE WHEN coalesce(reason, '') LIKE 'out_of_scope%' "
        "             THEN 1 ELSE 0 END) AS oos, "
        "         max(CASE WHEN reason = 'on_demand_pending' "
        "             THEN 1 ELSE 0 END) AS on_demand, "
        f"        max(CASE WHEN decision IN {_sql_in(QUOTED_DECISIONS)} "
        "             THEN 1 ELSE 0 END) AS quoted, "
        f"        max(CASE WHEN decision IN {_sql_in(ACCEPT_DECISIONS)} "
        "             THEN 1 ELSE 0 END) AS accepted, "
        "         max(CASE WHEN decision = 'confirmed' "
        "             THEN 1 ELSE 0 END) AS confirmed "
        "  FROM quote_decisions "
        "  WHERE observed_at >= ? AND rfq_id IS NOT NULL "
        "  GROUP BY rfq_id) "
        "SELECT CASE WHEN oos = 1 THEN 'out_of_scope' "
        "            WHEN on_demand = 1 THEN 'on_demand' "
        "            ELSE 'grid' END AS path, "
        "       count(*), sum(quoted), sum(accepted), sum(confirmed) "
        "FROM per_rfq GROUP BY 1", [since])
    if not _usable(per_rfq):
        out["unavailable"] = _unavailable_note(per_rfq, "state")
        return out
    for path, rfqs, quoted, accepted, confirmed in per_rfq:
        out["paths"][path] = {"rfqs": rfqs, "quoted": quoted}
        out["seen"] += rfqs
        if path != "out_of_scope":
            out["in_scope"] += rfqs
        out["quoted"] += quoted
        out["accepted"] += accepted
        out["confirmed"] += confirmed

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


MARGIN_BUCKETS = ((0.01, "<1c"), (0.02, "1-2c"), (0.03, "2-3c"),
                  (0.04, "3-4c"), (0.05, "4-5c"), (None, ">=5c"))
FAIR_BANDS = ((0.10, "<.10"), (0.25, "[.10,.25)"), (0.50, "[.25,.50)"),
              (0.75, "[.50,.75)"), (0.90, "[.75,.90)"), (None, ">=.90"))


def _bucket(value: float, edges) -> str:
    for edge, label in edges:
        if edge is None or value < edge:
            return label
    raise AssertionError("unreachable")


def demand_stats(state_db: str, since: datetime) -> dict:
    """Quotes vs accepts by (margin bucket x fair-prob band) — feeds #26.

    Each quoted decision contributes two quote-sides: yes at fair
    blended_fair with margin fair - yes_bid, no at fair 1 - blended_fair
    with margin (1 - fair) - no_bid. An accept (fill) lands in the bucket
    of the side actually held, margin = side fair at quote - price paid.
    """
    out = {"cells": {}, "quote_sides": 0, "accepts": 0, "unavailable": None}

    def _cell(margin: float, fair: float) -> dict:
        key = (_bucket(margin, MARGIN_BUCKETS), _bucket(fair, FAIR_BANDS))
        return out["cells"].setdefault(key, {"quotes": 0, "accepts": 0})

    quotes = _read(
        state_db,
        f"SELECT blended_fair, yes_bid, no_bid FROM quote_decisions "
        f"WHERE decision IN {_sql_in(QUOTED_DECISIONS)} "
        "AND observed_at >= ? AND blended_fair IS NOT NULL", [since])
    if not _usable(quotes):
        out["unavailable"] = _unavailable_note(quotes, "state")
        return out
    for fair, yes_bid, no_bid in quotes:
        for side_fair, bid in ((fair, yes_bid), (1.0 - fair, no_bid)):
            if bid is None:
                continue
            _cell(side_fair - bid, side_fair)["quotes"] += 1
            out["quote_sides"] += 1

    fills = _read(state_db,
                  "SELECT side_held, price, blended_fair_at_quote FROM fills "
                  "WHERE contracts > 0 AND filled_at >= ? "
                  "AND blended_fair_at_quote IS NOT NULL "
                  "AND side_held IN ('yes', 'no')", [since])
    if _usable(fills):
        for side_held, price, fair_at_quote in fills:
            side_fair = (fair_at_quote if side_held == "yes"
                         else 1.0 - fair_at_quote)
            _cell(side_fair - price, side_fair)["accepts"] += 1
            out["accepts"] += 1
    return out


def _render_demand(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    if stats["quote_sides"] == 0:
        return "_no quotes in window._"
    band_labels = [label for _edge, label in FAIR_BANDS]
    lines = [
        f"{stats['quote_sides']} quote-sides, {stats['accepts']} accepts. "
        "Cells are `quotes/accepts` by margin (rows) x fair prob (cols):",
        "",
        "| margin \\ fair | " + " | ".join(band_labels) + " |",
        "|---|" + "---:|" * len(band_labels),
    ]
    for _edge, margin_label in MARGIN_BUCKETS:
        row = [f"| {margin_label} "]
        for band_label in band_labels:
            cell = stats["cells"].get((margin_label, band_label))
            row.append(f"| {cell['quotes']}/{cell['accepts']} " if cell
                       else "| ")
        lines.append("".join(row) + "|")
    lines += ["", "_Feeds #26 (margin experimentation ladder)._"]
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
        "**Decision / reason breakdown (7d)** — decision ROWS, not RFQs "
        "(an open in-scope RFQ is re-decided every discovery tick):",
        "",
        "| decision | reason | n |",
        "|---|---|---:|",
    ]
    for decision, reason, n in d7["decisions"]:
        lines.append(f"| {decision} | {reason or '—'} | {n} |")
    if not d7["decisions"]:
        lines.append("| _none_ | | |")
    return "\n".join(lines)

def pnl_stats(state_db: str) -> dict:
    """All-time settlement P&L (research_queries.sql §7-8 math, no markouts).

    Quoted margin per contract is side-adjusted (yes: fair - price;
    no: (1 - fair) - price). Realized below quoted beyond noise = adverse
    selection eating the margin — the measurement-phase headline.
    """
    out = {"fills": 0, "contracts": 0.0, "settled_fills": 0,
           "unsettled_fills": 0, "realized_pnl_total": None,
           "avg_quoted_margin_per_ct": None, "realized_per_ct": None,
           "by_result": [], "unavailable": None}

    totals = _read(
        state_db,
        "SELECT count(*), coalesce(sum(contracts), 0), "
        "       count(*) FILTER (realized_pnl IS NOT NULL) "
        "FROM fills WHERE contracts > 0")
    if not _usable(totals):
        out["unavailable"] = _unavailable_note(totals, "state")
        return out
    if totals:
        out["fills"], out["contracts"], out["settled_fills"] = totals[0]
        out["unsettled_fills"] = out["fills"] - out["settled_fills"]

    settled = _read(
        state_db,
        "SELECT sum(realized_pnl), sum(contracts), "
        "       avg(CASE WHEN side_held = 'yes' "
        "                THEN blended_fair_at_quote - price "
        "                ELSE (1 - blended_fair_at_quote) - price END) "
        "FROM fills WHERE contracts > 0 AND realized_pnl IS NOT NULL "
        "AND blended_fair_at_quote IS NOT NULL")
    if _usable(settled) and settled and settled[0][0] is not None:
        pnl_total, settled_contracts, avg_margin = settled[0]
        out["realized_pnl_total"] = pnl_total
        out["avg_quoted_margin_per_ct"] = avg_margin
        if settled_contracts:
            out["realized_per_ct"] = pnl_total / settled_contracts

    by_result = _read(
        state_db,
        "SELECT s.result, count(*), sum(f.realized_pnl) "
        "FROM fills f JOIN settlements s "
        "  ON f.combo_market_ticker = s.combo_market_ticker "
        "WHERE f.contracts > 0 AND f.realized_pnl IS NOT NULL "
        "GROUP BY s.result ORDER BY s.result")
    if _usable(by_result):
        out["by_result"] = [tuple(r) for r in by_result]
    return out


def _render_pnl(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    if stats["fills"] == 0:
        return "no fills yet."
    if stats["settled_fills"] == 0:
        return (f"{stats['fills']} fills recorded, none settled yet "
                "(settlement sweep may not have run since the #12 merge).")
    lines = [
        f"- fills: {stats['fills']} ({stats['contracts']:.0f} contracts), "
        f"settled {stats['settled_fills']}, "
        f"unsettled {stats['unsettled_fills']}",
        f"- realized P&L (settled): ${stats['realized_pnl_total']:.2f}",
        f"- quoted margin/ct {_fmt(stats['avg_quoted_margin_per_ct'], '.3f')} "
        f"vs realized/ct {_fmt(stats['realized_per_ct'], '.3f')} — "
        "realized persistently below quoted = adverse selection",
        "",
        "| result | fills | P&L ($) |",
        "|---|---:|---:|",
    ]
    for result, n, pnl in stats["by_result"]:
        lines.append(f"| {result} | {n} | {pnl:.2f} |")
    lines += ["", "_Markouts descoped (#13); settlement P&L is the "
              "adverse-selection signal. Attribution detail: "
              "research_queries.sql §7-8._"]
    return "\n".join(lines)


VOID_DECISIONS = ("voided_no_legs", "voided_singles_moved",
                  "voided_no_fresh_books", "voided_last_look")


def health_stats(state_db: str, research_db: str, since: datetime) -> dict:
    """Void rate, sweep cancels, halts, phantom fills, reconcile outcomes."""
    import json

    out = {"confirmed": 0, "voids_by_decision": {}, "void_rate": None,
           "sweep_cancels_by_reason": {}, "circuit_breakers": 0,
           "void_rate_halts": 0, "halt_transitions": {}, "phantom_fills": 0,
           "reconcile_outcomes": {}, "unavailable": None}

    decisions = _read(
        state_db,
        "SELECT decision, coalesce(reason, ''), count(*) "
        "FROM quote_decisions WHERE observed_at >= ? "
        f"AND (decision IN {_sql_in(VOID_DECISIONS)} "
        "     OR decision IN ('confirmed', 'sweep_cancel', "
        "                     'circuit_breaker', 'halted_high_void_rate')) "
        "GROUP BY 1, 2", [since])
    if not _usable(decisions):
        out["unavailable"] = _unavailable_note(decisions, "state")
        return out
    for decision, reason, n in decisions:
        if decision == "confirmed":
            out["confirmed"] += n
        elif decision in VOID_DECISIONS:
            out["voids_by_decision"][decision] = (
                out["voids_by_decision"].get(decision, 0) + n)
        elif decision == "sweep_cancel":
            out["sweep_cancels_by_reason"][reason] = (
                out["sweep_cancels_by_reason"].get(reason, 0) + n)
        elif decision == "circuit_breaker":
            out["circuit_breakers"] += n
        elif decision == "halted_high_void_rate":
            out["void_rate_halts"] += n
    voids = sum(out["voids_by_decision"].values())
    if voids + out["confirmed"] > 0:
        out["void_rate"] = voids / (voids + out["confirmed"])

    phantoms = _read(state_db,
                     "SELECT count(*) FROM fills "
                     "WHERE contracts = 0 AND reconciled AND filled_at >= ?",
                     [since])
    if _usable(phantoms) and phantoms:
        out["phantom_fills"] = phantoms[0][0]

    events = _read(research_db,
                   "SELECT event_type, payload FROM events "
                   "WHERE event_type IN ('halt_event', 'reconcile_done') "
                   "AND ts >= ?", [since])
    if _usable(events):
        for event_type, payload_str in events:
            payload = json.loads(payload_str)
            if event_type == "halt_event":
                key = payload.get("transition", "?")
                out["halt_transitions"][key] = (
                    out["halt_transitions"].get(key, 0) + 1)
            else:
                key = payload.get("outcome", "?")
                out["reconcile_outcomes"][key] = (
                    out["reconcile_outcomes"].get(key, 0) + 1)
    return out


def _render_counts(d: dict) -> str:
    return ", ".join(f"{k}: {v}" for k, v in sorted(d.items())) or "none"


def _render_health(stats: dict) -> str:
    if stats["unavailable"]:
        return stats["unavailable"]
    voids = sum(stats["voids_by_decision"].values())
    lines = [
        f"- accepts: {stats['confirmed'] + voids} → confirmed "
        f"{stats['confirmed']}, voided {voids} "
        f"(void rate {_fmt(stats['void_rate'], '.0%')})",
        f"- voids by reason: {_render_counts(stats['voids_by_decision'])}",
        f"- sweep cancels: {_render_counts(stats['sweep_cancels_by_reason'])}",
        f"- circuit breakers: {stats['circuit_breakers']}, "
        f"void-rate halts: {stats['void_rate_halts']} "
        f"(transitions: {_render_counts(stats['halt_transitions'])})",
        f"- phantom fills: {stats['phantom_fills']}, reconcile outcomes: "
        f"{_render_counts(stats['reconcile_outcomes'])}",
    ]
    return "\n".join(lines)


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
        "## 4. Demand curve (7d)",
        "",
        _render_demand(demand_stats(state_db, d7)),
        "",
        "## 5. Settlement P&L (all-time)",
        "",
        _render_pnl(pnl_stats(state_db)),
        "",
        "## 6. Health (7d)",
        "",
        _render_health(health_stats(state_db, research_db, d7)),
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
