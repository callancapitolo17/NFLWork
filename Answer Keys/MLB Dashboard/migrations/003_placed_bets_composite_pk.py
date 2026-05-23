"""Promote placed_bets PRIMARY KEY from bet_hash to (bet_hash, account).

Backfills NULL or empty `account` values to 'Wagerzon' (the legacy
single-account default) before rebuilding the table with the new PK.

Idempotent — safe to re-run; second invocation is a no-op.

Run once against the live mlb_dashboard.duckdb:

    python "Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py" \\
        "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"

Restart the dashboard server after running — the running R process holds
a connection to the old schema and will not see the PK change without a
restart.
"""
from __future__ import annotations
import sys
import duckdb


def _pk_columns(con: duckdb.DuckDBPyConnection) -> set[str]:
    rows = con.execute("""
        SELECT UNNEST(constraint_column_names)
        FROM duckdb_constraints()
        WHERE table_name = 'placed_bets' AND constraint_type = 'PRIMARY KEY'
    """).fetchall()
    return {r[0] for r in rows}


def run(db_path: str) -> None:
    con = duckdb.connect(db_path)
    try:
        if _pk_columns(con) == {"bet_hash", "account"}:
            return  # already migrated

        con.execute("""
            UPDATE placed_bets
               SET account = 'Wagerzon'
             WHERE account IS NULL OR account = ''
        """)

        cols = con.execute("DESCRIBE placed_bets").fetchall()
        # DESCRIBE returns (column_name, column_type, null, key, default, extra)
        col_defs = ", ".join(f'"{c[0]}" {c[1]}' for c in cols)

        con.execute("DROP TABLE IF EXISTS placed_bets_new")
        con.execute(
            f"CREATE TABLE placed_bets_new ({col_defs}, "
            f"PRIMARY KEY (bet_hash, account))"
        )
        con.execute("INSERT INTO placed_bets_new SELECT * FROM placed_bets")
        con.execute("DROP TABLE placed_bets")
        con.execute("ALTER TABLE placed_bets_new RENAME TO placed_bets")
    finally:
        con.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path-to-mlb_dashboard.duckdb>",
              file=sys.stderr)
        sys.exit(1)
    run(sys.argv[1])
    print(f"Migration 003 applied to {sys.argv[1]}")
