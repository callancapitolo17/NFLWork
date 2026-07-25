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


# ---------------------------------------------------------------------------
# Sections (each returns a markdown fragment; none may raise)
# ---------------------------------------------------------------------------

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
    lines += [
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
