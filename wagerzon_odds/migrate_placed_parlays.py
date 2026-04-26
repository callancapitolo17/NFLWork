# wagerzon_odds/migrate_placed_parlays.py
"""One-shot migration: ensures placed_parlays has all columns required by
parlay_placer + creates placement_orphans. Idempotent (safe to re-run)."""
import os
import duckdb

DB_PATH = os.path.join(
    os.path.dirname(__file__), "..", "Answer Keys", "MLB Dashboard",
    "mlb_dashboard.duckdb"
)

REQUIRED_COLS = {
    "parlay_hash":      "VARCHAR PRIMARY KEY",
    "status":           "VARCHAR NOT NULL",
    "combo":            "VARCHAR",
    "game_id":          "VARCHAR",
    "game_time":        "TIMESTAMP",
    "recommended_size": "DOUBLE",
    "expected_odds":    "INTEGER",
    "expected_win":     "DOUBLE",
    "actual_size":      "DOUBLE",
    "actual_win":       "DOUBLE",
    "ticket_number":    "VARCHAR",
    "idwt":             "BIGINT",
    "legs_json":        "VARCHAR",
    "error_msg":        "VARCHAR",
    "error_msg_key":    "VARCHAR",
    "placed_at":        "TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
    "updated_at":       "TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
}

def main():
    con = duckdb.connect(DB_PATH)
    try:
        # Create placed_parlays if missing
        con.execute(f"""
            CREATE TABLE IF NOT EXISTS placed_parlays (
                parlay_hash VARCHAR PRIMARY KEY,
                status      VARCHAR NOT NULL,
                placed_at   TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Add any missing columns
        existing = {r[1] for r in con.execute(
            "PRAGMA table_info(placed_parlays)"
        ).fetchall()}
        for col, ddl in REQUIRED_COLS.items():
            if col in existing:
                continue
            # PK / NOT NULL constraints can't be added via ALTER. If they're
            # missing on an existing table, that's fine (the create-table
            # above guarantees PK; NOT NULL is enforced application-side).
            ddl_alter = ddl.replace("PRIMARY KEY", "").replace("NOT NULL", "").strip()
            con.execute(f"ALTER TABLE placed_parlays ADD COLUMN IF NOT EXISTS {col} {ddl_alter}")
            print(f"  added column placed_parlays.{col}")

        # Create placement_orphans
        con.execute("""
            CREATE TABLE IF NOT EXISTS placement_orphans (
                idwt          BIGINT PRIMARY KEY,
                ticket_number VARCHAR,
                parlay_hash   VARCHAR,
                raw_response  VARCHAR,
                error         VARCHAR,
                created_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        print("placement_orphans ready")
    finally:
        con.close()

if __name__ == "__main__":
    main()
