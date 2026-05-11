"""One-shot migration: add WZ single-placer columns to placed_bets.

Mirrors 001_combined_parlay_columns.py. Idempotent — safe to re-run.

Run once against the live mlb_dashboard.duckdb:

    python "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py" \
        "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
"""
from __future__ import annotations
import sys
import duckdb


COLUMNS = [
    ("account",          "VARCHAR"),
    ("status",           "VARCHAR DEFAULT 'placed'"),
    ("ticket_number",    "VARCHAR"),
    ("error_msg",        "VARCHAR"),
    ("error_msg_key",    "VARCHAR"),
    ("wz_odds_at_place", "INTEGER"),
]


def run(db_path: str) -> None:
    con = duckdb.connect(db_path)
    try:
        existing = {row[0] for row in con.execute("DESCRIBE placed_bets").fetchall()}
        for name, ddl in COLUMNS:
            if name not in existing:
                con.execute(f"ALTER TABLE placed_bets ADD COLUMN {name} {ddl}")
    finally:
        con.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path-to-mlb_dashboard.duckdb>", file=sys.stderr)
        sys.exit(1)
    run(sys.argv[1])
    print(f"Migration 002 applied to {sys.argv[1]}")
