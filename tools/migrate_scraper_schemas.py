"""Coordinated schema migration for all scraper DuckDBs.

Part of the MLB scrapers TZ-standardization refactor. Phase 1 audit
(see tools/TZ_AUDIT_FINDINGS.md) confirmed that BKM and Hoop88 record
game times in PT despite being treated as ET by R-side code, while
WZ is EDT and Bet105 is UTC. The fix replaces the per-scraper
(game_date VARCHAR, game_time VARCHAR) pair with a single
game_start_time TIMESTAMPTZ column in UTC.

Ephemeral scraper tables (dk_odds, fd_odds, bfa_odds, kalshi_odds,
wagerzon_odds, hoop88_odds, bookmaker_odds, bet105_odds): DROP TABLE
mlb_odds. The scraper recreates with the new schema on next cycle.

Persistent table (wagerzon_specials, lives alongside mlb_odds inside
wagerzon_odds/wagerzon.duckdb): ALTER ADD game_start_time TIMESTAMPTZ,
UPDATE from existing naive ET pair, leave game_date+game_time in place
for one cycle as a safety net (drop in a follow-up commit).

Idempotent. --rollback drops new column / restores old schema.

Usage:
    python tools/migrate_scraper_schemas.py            # apply migration
    python tools/migrate_scraper_schemas.py --rollback # undo migration

DO NOT RUN until the corresponding Phase 3 scraper code has landed —
otherwise the next scrape cycle will fail trying to write into a
table that no longer exists (it'll auto-recreate with the OLD schema,
which is harmless but wastes a cycle).
"""
import argparse
import duckdb
from pathlib import Path

# Resolve repo root relative to this file (tools/ is one level below root)
# so the script works regardless of which directory it's invoked from.
REPO_ROOT = Path(__file__).resolve().parent.parent

EPHEMERAL_DBS = {
    "dk_odds/dk.duckdb":               "dk",
    "fd_odds/fd.duckdb":               "fd",
    "bfa_odds/bfa.duckdb":             "bfa",
    "kalshi_odds/kalshi.duckdb":       "kalshi",
    "wagerzon_odds/wagerzon.duckdb":   "wagerzon",
    "hoop88_odds/hoop88.duckdb":       "hoop88",
    "bookmaker_odds/bookmaker.duckdb": "bookmaker",
    "bet105_odds/bet105.duckdb":       "bet105",
}
WZ_SPECIALS_DB = "wagerzon_odds/wagerzon.duckdb"


def migrate_ephemeral(rollback=False):
    """Drop mlb_odds in each ephemeral scraper DB.

    Rollback is a no-op: the table will be recreated by the next
    scraper cycle either way, and we can't restore the pre-migration
    schema from this side (the schema is owned by the scraper).
    """
    if rollback:
        print("[ephemeral] rollback is a no-op — scrapers own the schema; "
              "revert scraper code to restore the old shape")
        return
    for rel_path, label in EPHEMERAL_DBS.items():
        db_path = REPO_ROOT / rel_path
        if not db_path.exists():
            print(f"[{label}] DB missing at {db_path} — skip")
            continue
        con = duckdb.connect(str(db_path))
        try:
            con.execute("DROP TABLE IF EXISTS mlb_odds")
            print(f"[{label}] dropped mlb_odds (will be recreated by next scraper run)")
        finally:
            con.close()


def migrate_specials(rollback=False):
    """Add/backfill game_start_time on wagerzon_specials (persistent).

    Per scraper_specials.py convention, the existing game_date/game_time
    columns are naive Eastern Time. We interpret them as America/New_York
    when backfilling. If the table is missing or the source columns
    aren't present, we log and skip rather than crash.
    """
    db_path = REPO_ROOT / WZ_SPECIALS_DB
    if not db_path.exists():
        print(f"[wz_specials] DB missing at {db_path} — skip")
        return
    con = duckdb.connect(str(db_path))
    try:
        # Verify the target table exists before touching it.
        tables = [r[0] for r in con.execute("SHOW TABLES").fetchall()]
        if "wagerzon_specials" not in tables:
            print("[wz_specials] table wagerzon_specials not found — skip "
                  "(nothing to migrate; will be created by scraper)")
            return

        cols = [r[0] for r in con.execute("DESCRIBE wagerzon_specials").fetchall()]

        if rollback:
            if "game_start_time" in cols:
                con.execute("ALTER TABLE wagerzon_specials DROP COLUMN game_start_time")
                print("[wz_specials] rolled back: dropped game_start_time")
            else:
                print("[wz_specials] rollback no-op: game_start_time already absent")
            return

        # Forward migration: add the column if it isn't already there.
        if "game_start_time" not in cols:
            con.execute("ALTER TABLE wagerzon_specials ADD COLUMN game_start_time TIMESTAMPTZ")
            print("[wz_specials] added column game_start_time TIMESTAMPTZ")
        else:
            print("[wz_specials] game_start_time already present — will backfill any NULLs")

        # Only attempt the backfill if both source columns exist. Without
        # them there's nothing to convert and the UPDATE would error.
        if "game_date" in cols and "game_time" in cols:
            print("[wz_specials] backfill assumes existing rows are "
                  "America/New_York — see scraper_specials.py")
            con.execute("""
                UPDATE wagerzon_specials
                SET game_start_time = (CAST(game_date AS DATE) + CAST(game_time AS TIME))
                                       AT TIME ZONE 'America/New_York'
                WHERE game_start_time IS NULL
                  AND game_date IS NOT NULL AND game_time IS NOT NULL
            """)
            filled = con.execute(
                "SELECT count(*) FROM wagerzon_specials WHERE game_start_time IS NOT NULL"
            ).fetchone()[0]
            total = con.execute("SELECT count(*) FROM wagerzon_specials").fetchone()[0]
            print(f"[wz_specials] backfill complete: {filled}/{total} rows have game_start_time")
        else:
            missing = [c for c in ("game_date", "game_time") if c not in cols]
            print(f"[wz_specials] source columns missing {missing} — skipping backfill "
                  "(column added, but nothing to populate from)")

        # NOTE: Drop game_date/game_time in a follow-up commit after we confirm
        # post-merge stability. Keep them in place for one cycle as a safety net.
    finally:
        con.close()


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--rollback", action="store_true",
                   help="Undo the migration (drops game_start_time on wagerzon_specials)")
    args = p.parse_args()
    migrate_ephemeral(rollback=args.rollback)
    migrate_specials(rollback=args.rollback)


if __name__ == "__main__":
    main()
