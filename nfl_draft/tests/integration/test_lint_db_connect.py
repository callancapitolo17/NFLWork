"""Forbid raw duckdb.connect() outside allowlist — concurrency discipline.

Every new DuckDB callsite in ``nfl_draft/`` and ``kalshi_draft/`` must go
through the short-lived context-managed helpers in ``nfl_draft/lib/db.py``
(``write_connection`` / ``read_connection``). This keeps the lock-held
window bounded to milliseconds, which is what makes the cron writer + Dash
reader coexist safely.

The allowlist below documents the few files that legitimately need a raw
``duckdb.connect()`` — each with an inline comment explaining why.
"""
import re
from pathlib import Path

ALLOWLIST = {
    # Defines the context-managed helpers themselves.
    "nfl_draft/lib/db.py",
    # One-time migration opens two DBs concurrently (old + new). Not a hot path.
    "nfl_draft/lib/migrate_from_kalshi_draft.py",
    # Legacy writer shim; callers are fetcher.py / consensus.py / edge_detector.py
    # which run inside the cron process (they ARE the writer, not a reader).
    # New code must still go through write_connection / read_connection.
    "kalshi_draft/db.py",
}

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
PATTERN = re.compile(r"\bduckdb\.connect\b")


def test_no_raw_duckdb_connect_outside_allowlist():
    violations = []
    # Scan BOTH nfl_draft/ and kalshi_draft/ — we want this discipline across
    # the whole draft pipeline, not just the new package.
    roots = [REPO_ROOT / "nfl_draft", REPO_ROOT / "kalshi_draft"]
    for root in roots:
        if not root.exists():
            continue
        for py_file in root.rglob("*.py"):
            rel = py_file.relative_to(REPO_ROOT).as_posix()
            if rel in ALLOWLIST:
                continue
            if "/tests/" in rel or rel.endswith("/conftest.py"):
                continue
            text = py_file.read_text()
            if PATTERN.search(text):
                violations.append(rel)
    assert not violations, f"Raw duckdb.connect found in: {violations}"
