"""Forbid raw duckdb.connect() outside lib/db.py — concurrency discipline."""
import re
from pathlib import Path

ALLOWLIST = {
    "nfl_draft/lib/db.py",
    "nfl_draft/lib/migrate_from_kalshi_draft.py",  # legitimately needs two-DB open
}

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
PATTERN = re.compile(r"\bduckdb\.connect\b")


def test_no_raw_duckdb_connect_outside_allowlist():
    violations = []
    for py_file in (REPO_ROOT / "nfl_draft").rglob("*.py"):
        rel = py_file.relative_to(REPO_ROOT).as_posix()
        if rel in ALLOWLIST:
            continue
        if "/tests/" in rel or rel.endswith("/conftest.py"):
            continue
        text = py_file.read_text()
        if PATTERN.search(text):
            violations.append(rel)
    assert not violations, f"Raw duckdb.connect found in: {violations}"
