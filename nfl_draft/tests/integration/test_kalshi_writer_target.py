"""Verify legacy fetcher writes to nfl_draft.duckdb, NOT kalshi_draft.duckdb."""
import pytest
from pathlib import Path
from nfl_draft.lib import db as db_module


def test_fetcher_db_path_resolves_to_nfl_draft(monkeypatch):
    """Importing the legacy fetcher must use the new DB_PATH."""
    from kalshi_draft import db as legacy_db
    # The legacy db.py was repointed in Task 4
    assert "nfl_draft.duckdb" in str(legacy_db.DB_PATH)
    assert "kalshi_draft.duckdb" not in str(legacy_db.DB_PATH)
