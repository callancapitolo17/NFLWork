"""One-time migration: kalshi_draft.duckdb -> nfl_draft.duckdb.

Idempotent: uses MD5 hashes to detect already-migrated tables.
Tables copied:
  draft_odds -> kalshi_odds (renamed; legacy schema preserved)
  draft_series, market_info, positions, resting_orders -> kept as-is
  consensus_board, detected_edges -> kept as-is
"""

import duckdb
from pathlib import Path
from typing import Optional

# Allowed: this module legitimately needs raw duckdb.connect for two-DB open
# (lint: nfl_draft.tests.integration.test_lint_db_connect allowlist)


TABLE_MAP = {
    "draft_odds": "kalshi_odds",        # rename
    "draft_series": "draft_series",
    "market_info": "market_info",
    "positions": "positions",
    "resting_orders": "resting_orders",
    "consensus_board": "consensus_board",
    "detected_edges": "detected_edges",
}


def _table_hash(con: duckdb.DuckDBPyConnection, table: str) -> Optional[str]:
    try:
        result = con.execute(
            f"SELECT MD5(GROUP_CONCAT(CAST(t AS TEXT) ORDER BY CAST(t AS TEXT))) FROM (SELECT * FROM {table}) t"
        ).fetchone()
        return result[0] if result else None
    except duckdb.CatalogException:
        return None  # table doesn't exist


def run(legacy_path: Path, new_path: Path) -> None:
    """Copy tables from legacy DB into new DB, with rename + idempotency."""
    if not Path(legacy_path).exists():
        print(f"Legacy DB not found at {legacy_path} - nothing to migrate.")
        return

    legacy = duckdb.connect(str(legacy_path), read_only=True)
    new = duckdb.connect(str(new_path))
    try:
        for src, dst in TABLE_MAP.items():
            src_hash = _table_hash(legacy, src)
            if src_hash is None:
                continue
            dst_hash = _table_hash(new, dst)
            if src_hash == dst_hash:
                print(f"  {src} -> {dst}: already migrated (hash match), skipping")
                continue
            # Copy schema + data via DuckDB ATTACH
            new.execute(f"ATTACH '{legacy_path}' AS legacy (READ_ONLY)")
            new.execute(f"DROP TABLE IF EXISTS {dst}")
            new.execute(f"CREATE TABLE {dst} AS SELECT * FROM legacy.{src}")
            new.execute("DETACH legacy")
            print(f"  {src} -> {dst}: copied {new.execute(f'SELECT COUNT(*) FROM {dst}').fetchone()[0]} rows")
    finally:
        legacy.close()
        new.close()


if __name__ == "__main__":
    from nfl_draft.lib import db as db_module
    legacy = Path(__file__).resolve().parent.parent.parent / "kalshi_draft" / "kalshi_draft.duckdb"
    db_module.init_schema()
    run(legacy_path=legacy, new_path=db_module.DB_PATH)
