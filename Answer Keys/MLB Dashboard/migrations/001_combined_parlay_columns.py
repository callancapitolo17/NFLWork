"""One-shot migration: add combo tracking columns to placed_parlays.

Run once against the live mlb_dashboard.duckdb to make the dashboard combo-aware:

    python "Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py" \
        "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
"""
from __future__ import annotations
import sys
import duckdb


def run(db_path: str) -> None:
    con = duckdb.connect(db_path)
    try:
        existing = {row[0] for row in con.execute("DESCRIBE placed_parlays").fetchall()}

        if "is_combo" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN is_combo BOOLEAN DEFAULT FALSE")
        if "combo_leg_ids" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN combo_leg_ids VARCHAR")
        if "parent_combo_id" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN parent_combo_id INTEGER")
    finally:
        con.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path-to-mlb_dashboard.duckdb>", file=sys.stderr)
        sys.exit(2)
    run(sys.argv[1])
    print(f"Migration 001 applied to {sys.argv[1]}")
