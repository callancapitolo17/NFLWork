"""Loader for the named queries in ``fetch_health_queries.sql`` (issue #37).

WHY A LOADER AND NOT COPIES
    The monitor dashboard has to answer "is this book dead?", which is
    exactly what ``current_failure_streak`` already answers — including the
    ``outcome IN ('ok','empty')`` predicate that is the trap of this whole
    feature ('empty' means the book ANSWERED with no markets). A second copy
    of that predicate in Python is a second place to get it wrong.

WHY ITS OWN MODULE
    ``sgp_health`` imports ``mlb_sgp._shared``. ``kalshi_mlb_monitor`` is a
    read-only dashboard that deliberately imports no bot code and runs as a
    separate process; this module imports nothing but ``pathlib`` so it can
    be used from there without dragging in the scraper package.
"""
from __future__ import annotations

from pathlib import Path

QUERY_FILE = Path(__file__).resolve().parent / "fetch_health_queries.sql"

_HEADER = "-- name:"


def _parse(text: str) -> dict[str, str]:
    blocks: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if line.startswith(_HEADER):
            current = line.split(":", 1)[1].strip()
            blocks[current] = []
        elif current is not None:
            blocks[current].append(line)
    return {name: "\n".join(lines).strip() for name, lines in blocks.items()}


def all_queries() -> dict[str, str]:
    """Every ``-- name: <id>`` block in the shipped .sql, id -> SQL."""
    return _parse(QUERY_FILE.read_text())


def named_query(name: str) -> str:
    """One query by its ``-- name:`` header. Raises KeyError if absent —
    a typo'd name must fail loudly, not return an empty string that reads
    downstream as "no unhealthy books"."""
    queries = all_queries()
    if name not in queries:
        raise KeyError(f"no query named {name!r} in {QUERY_FILE.name}; "
                       f"have {sorted(queries)}")
    return queries[name]
