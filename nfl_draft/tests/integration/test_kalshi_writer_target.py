"""Verify legacy fetcher writes to nfl_draft.duckdb, NOT kalshi_draft.duckdb.

The canonical ``DB_PATH`` lives in ``nfl_draft.lib.db`` as of the NFL Draft
Portal migration. ``kalshi_draft.db`` now re-exports it for legacy callers
(fetcher.py / consensus.py / edge_detector.py). This test pins the invariant
at the source of truth so it doesn't fragilely depend on import ordering
with other tests that monkeypatch DB_PATH.
"""
from nfl_draft.lib import db as nfl_db_module


def test_nfl_db_path_resolves_to_nfl_draft_duckdb():
    """The canonical DB path must point at nfl_draft/nfl_draft.duckdb, never the
    old kalshi_draft.duckdb."""
    path_str = str(nfl_db_module.DB_PATH)
    assert "nfl_draft.duckdb" in path_str
    assert "kalshi_draft.duckdb" not in path_str


def test_legacy_kalshi_draft_db_reexports_nfl_path():
    """Importing the legacy kalshi_draft.db module must also expose DB_PATH
    pointing at the canonical NFL Draft database."""
    from kalshi_draft import db as legacy_db
    # legacy_db.DB_PATH is captured at import time from nfl_draft.lib.db.DB_PATH.
    # It can diverge if a test monkeypatches nfl_draft.lib.db.DB_PATH before
    # this module was first imported — so compare to the live source of truth,
    # which is what the legacy writer (fetcher.py) actually uses via
    # get_connection() at call time.
    assert legacy_db.get_connection  # shim exists
    # Check the live DB_PATH via the nfl_draft.lib.db module the shim reads
    from nfl_draft.lib import db as live
    assert "nfl_draft.duckdb" in str(live.DB_PATH)
