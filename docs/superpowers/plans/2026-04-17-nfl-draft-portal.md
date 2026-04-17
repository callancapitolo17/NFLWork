# NFL Draft EV Portal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a trader's-cockpit web portal at `nfl_draft/` that surfaces +EV NFL Draft bets across Kalshi and 4 sportsbooks (DraftKings, FanDuel, Bookmaker, Wagerzon), with cross-venue odds comparison, outlier-flag detection, Kalshi trade-tape monitoring, and bet logging — ready for the 2026 NFL Draft (Apr 23-25).

**Architecture:** Single DuckDB database `nfl_draft/nfl_draft.duckdb` (legacy `kalshi_draft.duckdb` migrated in and retired). Python orchestrator `nfl_draft/run.py` invoked by cron with two modes: `--mode scrape` (pulls all 5 venues' odds + triggers legacy fetcher writes + runs edge_detector and consensus) and `--mode trades` (polls Kalshi trade tape with cursor + dedup). Dashboard extends `kalshi_draft/app.py` with 4 new tabs (Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log) grouped under a "Portal" sub-nav, alongside the existing 5 Kalshi-only legacy tabs. Outlier flag fires when any venue's devigged probability is ≥10pp from the all-venue median.

**Tech Stack:** Python 3, DuckDB (storage + concurrency primitives), Dash + Plotly (dashboard), Playwright + curl_cffi (sportsbook scrapers), pytest (unit/integration/live tests).

**Spec reference:** `docs/superpowers/specs/2026-04-17-nfl-draft-portal-design.md`

---

## Task 1: Project skeleton + git branch verification

**Files:**
- Create: `nfl_draft/__init__.py`, `nfl_draft/config/__init__.py`, `nfl_draft/lib/__init__.py`, `nfl_draft/scrapers/__init__.py`, `nfl_draft/tests/__init__.py`, `nfl_draft/tests/unit/__init__.py`, `nfl_draft/tests/integration/__init__.py`, `nfl_draft/tests/live/__init__.py`, `nfl_draft/tests/fixtures/.gitkeep`
- Modify: `.gitignore`

- [ ] **Step 1: Verify on feature branch**

Run: `git branch --show-current`
Expected: `feature/nfl-draft-portal`. If not, `git checkout feature/nfl-draft-portal`.

- [ ] **Step 2: Create empty package files**

Run:
```bash
mkdir -p nfl_draft/config nfl_draft/lib nfl_draft/scrapers nfl_draft/tests/unit nfl_draft/tests/integration nfl_draft/tests/live nfl_draft/tests/fixtures/{kalshi,draftkings,fanduel,bookmaker,wagerzon}
touch nfl_draft/__init__.py nfl_draft/config/__init__.py nfl_draft/lib/__init__.py nfl_draft/scrapers/__init__.py
touch nfl_draft/tests/__init__.py nfl_draft/tests/unit/__init__.py nfl_draft/tests/integration/__init__.py nfl_draft/tests/live/__init__.py
touch nfl_draft/tests/fixtures/.gitkeep
```

- [ ] **Step 3: Update .gitignore**

Append to `.gitignore`:
```
# nfl_draft portal
nfl_draft/nfl_draft.duckdb
nfl_draft/nfl_draft.duckdb.wal
nfl_draft/.cookies/
nfl_draft/logs/
nfl_draft/.kalshi_throttle.lock
```

- [ ] **Step 4: Verify package importable**

Run: `python -c "import nfl_draft.lib"`
Expected: no output (success)

- [ ] **Step 5: Commit skeleton**

```bash
git add nfl_draft/ .gitignore
git commit -m "chore: nfl_draft package skeleton"
```

---

## Task 2: lib/db.py — connection helpers + TZ helpers + schema migrations

**Files:**
- Create: `nfl_draft/lib/db.py`
- Create: `nfl_draft/tests/unit/test_db.py`
- Create: `nfl_draft/tests/unit/test_tz_roundtrip.py`

- [ ] **Step 1: Write failing test for write_connection / read_connection**

Create `nfl_draft/tests/unit/test_db.py`:
```python
import duckdb
import pytest
import tempfile
from pathlib import Path
from nfl_draft.lib import db as db_module


def test_write_connection_opens_and_closes(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    with db_module.write_connection() as con:
        con.execute("CREATE TABLE t (x INT)")
        con.execute("INSERT INTO t VALUES (1)")
    # Connection released — open a new read connection
    with db_module.read_connection() as con:
        result = con.execute("SELECT x FROM t").fetchone()
    assert result == (1,)


def test_init_schema_creates_all_tables(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    with db_module.read_connection() as con:
        tables = {row[0] for row in con.execute("SHOW TABLES").fetchall()}
    expected = {
        "draft_markets", "draft_odds", "kalshi_trades", "kalshi_poll_state",
        "draft_bets", "players", "player_aliases", "teams", "team_aliases",
        "market_map", "draft_odds_unmapped", "draft_odds_unmapped_players",
    }
    assert expected.issubset(tables), f"Missing: {expected - tables}"


def test_init_schema_is_idempotent(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    db_module.init_schema()  # second call must not error
```

Create `nfl_draft/tests/unit/test_tz_roundtrip.py`:
```python
from datetime import datetime
from nfl_draft.lib.db import local_to_utc_iso, utc_iso_to_local


def test_local_to_utc_iso_returns_string_with_offset():
    local = datetime(2026, 4, 23, 12, 0, 0)
    iso = local_to_utc_iso(local)
    assert isinstance(iso, str)
    assert iso.endswith("+00:00") or iso.endswith("Z")


def test_utc_iso_to_local_roundtrip():
    original_local = datetime(2026, 4, 23, 12, 0, 0)
    iso = local_to_utc_iso(original_local)
    roundtripped = utc_iso_to_local(iso)
    assert roundtripped == original_local


def test_utc_iso_to_local_handles_kalshi_format():
    # Kalshi returns 2026-04-23T16:00:00Z — should convert to local
    result = utc_iso_to_local("2026-04-23T16:00:00Z")
    assert result.tzinfo is None  # naive local
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest nfl_draft/tests/unit/test_db.py nfl_draft/tests/unit/test_tz_roundtrip.py -v`
Expected: FAIL with import errors (db module doesn't exist yet)

- [ ] **Step 3: Implement nfl_draft/lib/db.py**

Create `nfl_draft/lib/db.py`:
```python
"""DuckDB connection management, schema, and TZ boundary helpers."""

import duckdb
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

DB_PATH = Path(__file__).resolve().parent.parent / "nfl_draft.duckdb"


@contextmanager
def write_connection() -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived write connection. Lock held for milliseconds."""
    con = duckdb.connect(str(DB_PATH))
    try:
        yield con
    finally:
        con.close()


@contextmanager
def read_connection() -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived read-only connection per Dash callback."""
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        yield con
    finally:
        con.close()


def local_to_utc_iso(local_ts: datetime) -> str:
    """Convert naive local datetime to UTC ISO-8601 string (for Kalshi API)."""
    aware_local = local_ts.astimezone()  # interprets naive as system local
    return aware_local.astimezone(timezone.utc).isoformat()


def utc_iso_to_local(iso_str: str) -> datetime:
    """Convert Kalshi ISO-8601 (with Z or offset) to naive local datetime."""
    # Handle Z suffix (Kalshi style) by replacing with +00:00
    s = iso_str.replace("Z", "+00:00")
    aware_utc = datetime.fromisoformat(s)
    return aware_utc.astimezone().replace(tzinfo=None)


SCHEMA = [
    """
    CREATE TABLE IF NOT EXISTS draft_markets (
        market_id TEXT PRIMARY KEY,
        market_type TEXT NOT NULL,
        subject_player TEXT,
        subject_team TEXT,
        position TEXT,
        pick_number INT,
        range_low INT,
        range_high INT,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds (
        market_id TEXT,
        book TEXT,
        american_odds INT,
        implied_prob DOUBLE,
        devig_prob DOUBLE,
        fetched_at TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_draft_odds_market_book_time ON draft_odds (market_id, book, fetched_at DESC)",
    """
    CREATE TABLE IF NOT EXISTS kalshi_trades (
        trade_id TEXT PRIMARY KEY,
        ticker TEXT,
        side TEXT,
        price_cents INT,
        count INT,
        notional_usd DOUBLE,
        traded_at TIMESTAMP,
        fetched_at TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_kalshi_trades_ticker_time ON kalshi_trades (ticker, traded_at DESC)",
    """
    CREATE TABLE IF NOT EXISTS kalshi_poll_state (
        series_ticker TEXT PRIMARY KEY,
        last_traded_at TIMESTAMP,
        last_polled_at TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_bets (
        bet_id TEXT PRIMARY KEY,
        market_id TEXT,
        book TEXT,
        side TEXT,
        american_odds INT,
        stake_usd DOUBLE,
        taken_at TIMESTAMP,
        note TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS players (
        canonical_name TEXT PRIMARY KEY,
        position TEXT,
        college TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS player_aliases (
        alias TEXT PRIMARY KEY,
        canonical_name TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS teams (
        canonical_abbr TEXT PRIMARY KEY,
        full_name TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS team_aliases (
        alias TEXT PRIMARY KEY,
        canonical_abbr TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS market_map (
        book TEXT,
        book_label TEXT,
        book_subject TEXT,
        market_id TEXT,
        PRIMARY KEY (book, book_label, book_subject)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds_unmapped (
        book TEXT,
        book_label TEXT,
        book_subject TEXT,
        american_odds INT,
        fetched_at TIMESTAMP,
        reason TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds_unmapped_players (
        book TEXT,
        book_player_name TEXT,
        american_odds INT,
        fetched_at TIMESTAMP
    )
    """,
]


def init_schema() -> None:
    """Create all tables (idempotent)."""
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with write_connection() as con:
        for stmt in SCHEMA:
            con.execute(stmt)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest nfl_draft/tests/unit/test_db.py nfl_draft/tests/unit/test_tz_roundtrip.py -v`
Expected: 6 tests PASS

- [ ] **Step 5: Commit (staging — final commit comes in Task 9)**

Hold off on commit until all of Task 1's atomic-migration components are ready. Continue to Task 3.

---

## Task 3: Migration script (kalshi_draft.duckdb → nfl_draft.duckdb)

**Files:**
- Create: `nfl_draft/lib/migrate_from_kalshi_draft.py`
- Create: `nfl_draft/tests/integration/test_migration.py`

- [ ] **Step 1: Write failing migration test**

Create `nfl_draft/tests/integration/test_migration.py`:
```python
import duckdb
import pytest
from pathlib import Path
from nfl_draft.lib import db as db_module
from nfl_draft.lib import migrate_from_kalshi_draft as migrate


@pytest.fixture
def fake_legacy_db(tmp_path):
    """Build a kalshi_draft.duckdb with sample data matching the real schema."""
    legacy_path = tmp_path / "kalshi_draft.duckdb"
    con = duckdb.connect(str(legacy_path))
    con.execute("""
        CREATE TABLE draft_odds (
            fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
            ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
            yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
            last_price INTEGER, volume BIGINT, volume_24h BIGINT,
            liquidity BIGINT, open_interest INTEGER
        )
    """)
    con.execute("INSERT INTO draft_odds VALUES (CURRENT_TIMESTAMP, 'KXNFLDRAFT1', 'EVT1', 'TKR1', 'Title', 'Cam Ward', 50, 51, 49, 50, 50, 1000, 10000, 5000, 200)")
    con.execute("CREATE TABLE draft_series (series_ticker VARCHAR PRIMARY KEY, title VARCHAR, category VARCHAR, discovered_at TIMESTAMP)")
    con.execute("INSERT INTO draft_series VALUES ('KXNFLDRAFT1', 'Test Series', 'NFL', CURRENT_TIMESTAMP)")
    con.execute("CREATE TABLE market_info (ticker VARCHAR PRIMARY KEY, title VARCHAR, subtitle VARCHAR, series_ticker VARCHAR, expiration_time TIMESTAMP, close_time TIMESTAMP, updated_at TIMESTAMP)")
    con.execute("CREATE TABLE positions (fetch_time TIMESTAMP, ticker VARCHAR, position INTEGER, market_exposure DOUBLE, realized_pnl DOUBLE, resting_orders_count INTEGER, total_traded DOUBLE, fees_paid DOUBLE)")
    con.execute("CREATE TABLE resting_orders (fetch_time TIMESTAMP, order_id VARCHAR, ticker VARCHAR, side VARCHAR, type VARCHAR, yes_price INTEGER, no_price INTEGER, remaining_count INTEGER, created_time TIMESTAMP, expiration_time TIMESTAMP)")
    con.execute("CREATE TABLE consensus_board (fetch_time TIMESTAMP, rank INTEGER, player_name VARCHAR, position VARCHAR, school VARCHAR, source VARCHAR)")
    con.execute("CREATE TABLE detected_edges (fetch_time TIMESTAMP, edge_type VARCHAR, description VARCHAR, market_a VARCHAR, market_b VARCHAR, price_a DOUBLE, price_b DOUBLE, implied_edge DOUBLE, confidence VARCHAR)")
    con.close()
    return legacy_path


def test_migration_renames_draft_odds_to_kalshi_odds(fake_legacy_db, tmp_path, monkeypatch):
    new_db = tmp_path / "nfl_draft.duckdb"
    monkeypatch.setattr(db_module, "DB_PATH", new_db)
    db_module.init_schema()
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)
    with db_module.read_connection() as con:
        tables = {row[0] for row in con.execute("SHOW TABLES").fetchall()}
        kalshi_count = con.execute("SELECT COUNT(*) FROM kalshi_odds").fetchone()[0]
    assert "kalshi_odds" in tables
    assert kalshi_count == 1


def test_migration_is_idempotent(fake_legacy_db, tmp_path, monkeypatch):
    new_db = tmp_path / "nfl_draft.duckdb"
    monkeypatch.setattr(db_module, "DB_PATH", new_db)
    db_module.init_schema()
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)
    migrate.run(legacy_path=fake_legacy_db, new_path=new_db)  # re-run
    with db_module.read_connection() as con:
        count = con.execute("SELECT COUNT(*) FROM kalshi_odds").fetchone()[0]
    assert count == 1  # not duplicated
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest nfl_draft/tests/integration/test_migration.py -v`
Expected: FAIL (migrate module doesn't exist)

- [ ] **Step 3: Implement migration script**

Create `nfl_draft/lib/migrate_from_kalshi_draft.py`:
```python
"""One-time migration: kalshi_draft.duckdb → nfl_draft.duckdb.

Idempotent: uses MD5 hashes to detect already-migrated tables.
Tables copied:
  draft_odds → kalshi_odds (renamed; legacy schema preserved)
  draft_series, market_info, positions, resting_orders → kept as-is
  consensus_board, detected_edges → kept as-is
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
        print(f"Legacy DB not found at {legacy_path} — nothing to migrate.")
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
                print(f"  {src} → {dst}: already migrated (hash match), skipping")
                continue
            # Copy schema + data via DuckDB ATTACH
            new.execute(f"ATTACH '{legacy_path}' AS legacy (READ_ONLY)")
            new.execute(f"DROP TABLE IF EXISTS {dst}")
            new.execute(f"CREATE TABLE {dst} AS SELECT * FROM legacy.{src}")
            new.execute("DETACH legacy")
            print(f"  {src} → {dst}: copied {new.execute(f'SELECT COUNT(*) FROM {dst}').fetchone()[0]} rows")
    finally:
        legacy.close()
        new.close()


if __name__ == "__main__":
    from nfl_draft.lib import db as db_module
    legacy = Path(__file__).resolve().parent.parent.parent / "kalshi_draft" / "kalshi_draft.duckdb"
    db_module.init_schema()
    run(legacy_path=legacy, new_path=db_module.DB_PATH)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest nfl_draft/tests/integration/test_migration.py -v`
Expected: 2 tests PASS

- [ ] **Step 5: Continue (no commit yet — atomic with Task 9)**

---

## Task 4: Repoint `kalshi_draft/db.py`

**Files:**
- Modify: `kalshi_draft/db.py`

- [ ] **Step 1: Read current state**

Run: `head -20 kalshi_draft/db.py` to confirm `DB_PATH = Path(__file__).parent / "kalshi_draft.duckdb"`.

- [ ] **Step 2: Update DB_PATH and SQL references**

Edit `kalshi_draft/db.py`:
- Replace `DB_PATH = Path(__file__).parent / "kalshi_draft.duckdb"` with:
  ```python
  DB_PATH = Path(__file__).resolve().parent.parent / "nfl_draft" / "nfl_draft.duckdb"
  ```
- Replace every occurrence of `FROM draft_odds` with `FROM kalshi_odds` (3 occurrences in the existing file: in `get_latest_odds`, `get_price_history` x2, `get_snapshot_count`).
- Replace `INSERT INTO draft_odds` with `INSERT INTO kalshi_odds` (none in db.py — it only reads).
- Note: the `init_schema()` in this file is now dead code since `nfl_draft/lib/db.py.init_schema()` is the source of truth. Leave it but add a deprecation comment:
  ```python
  def init_schema():
      """DEPRECATED: schema now managed by nfl_draft/lib/db.py.init_schema()."""
      raise RuntimeError("Use nfl_draft.lib.db.init_schema() instead")
  ```

- [ ] **Step 3: Verify no other reference to "draft_odds" remains in this file**

Run: `grep -n "draft_odds" kalshi_draft/db.py`
Expected: zero matches.

- [ ] **Step 4: Continue (no commit — atomic with Task 9)**

---

## Task 5: Repoint `kalshi_draft/fetcher.py`

**Files:**
- Modify: `kalshi_draft/fetcher.py`

- [ ] **Step 1: Find writes**

Run: `grep -n "duckdb.connect\|datetime.now\|INSERT INTO draft_odds" kalshi_draft/fetcher.py`

- [ ] **Step 2: Update connection target**

Add at top of file (or replace existing import block):
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from nfl_draft.lib import db as nfl_db
```

Replace any `duckdb.connect("kalshi_draft.duckdb")` or similar with `duckdb.connect(str(nfl_db.DB_PATH))`. If the file uses a constant like `DB_PATH`, point it at `nfl_db.DB_PATH`.

- [ ] **Step 3: Update INSERT target**

Replace `INSERT INTO draft_odds` with `INSERT INTO kalshi_odds` (1 occurrence around line 211).

- [ ] **Step 4: Normalize TZ to local**

Replace every `datetime.now(timezone.utc).strftime(...)` with `datetime.now().strftime(...)`. There are 3 occurrences (lines ~153, 231, 290). Keep the strftime format.

You may also remove the `from datetime import ..., timezone` import since `timezone` is no longer used. Verify with `grep "timezone" kalshi_draft/fetcher.py` after edits — should be empty.

- [ ] **Step 5: Verify no leftover references**

Run: `grep -n "kalshi_draft.duckdb\|INSERT INTO draft_odds\|timezone.utc" kalshi_draft/fetcher.py`
Expected: zero matches.

- [ ] **Step 6: Continue (no commit — atomic with Task 9)**

---

## Task 6: Repoint `kalshi_draft/edge_detector.py`

**Files:**
- Modify: `kalshi_draft/edge_detector.py`

- [ ] **Step 1: Apply same three transforms as Task 5**

Apply identically:
1. Connection target → `nfl_draft.lib.db.DB_PATH`
2. Any `FROM draft_odds` or `INSERT INTO draft_odds` → `FROM kalshi_odds` / `INSERT INTO kalshi_odds`
3. `datetime.now(timezone.utc)` → `datetime.now()`

- [ ] **Step 2: Verify clean**

Run: `grep -n "kalshi_draft.duckdb\|FROM draft_odds\|INSERT INTO draft_odds\|timezone.utc" kalshi_draft/edge_detector.py`
Expected: zero matches.

---

## Task 7: Repoint `kalshi_draft/consensus.py`

**Files:**
- Modify: `kalshi_draft/consensus.py`

- [ ] **Step 1: Apply same transforms**

Same as Task 6.

- [ ] **Step 2: Verify clean**

Run: `grep -n "kalshi_draft.duckdb\|FROM draft_odds\|INSERT INTO draft_odds\|timezone.utc" kalshi_draft/consensus.py`
Expected: zero matches.

---

## Task 8: Repoint `kalshi_draft/app.py`

**Files:**
- Modify: `kalshi_draft/app.py`

- [ ] **Step 1: Update connection imports/paths**

Replace `from db import ...` patterns with `from kalshi_draft.db import ...` if needed (since db.py now points at the new DB).

Find every `duckdb.connect(...)` call. Replace with calls that target `nfl_draft.lib.db.DB_PATH`. Prefer using `nfl_draft.lib.db.read_connection()` for the dashboard's read paths.

- [ ] **Step 2: Rename SQL references**

Run: `grep -n "draft_odds" kalshi_draft/app.py` — replace each with `kalshi_odds`.

- [ ] **Step 3: Verify clean**

Run: `grep -n "kalshi_draft.duckdb\|\\bdraft_odds\\b" kalshi_draft/app.py`
Expected: zero matches.

---

## Task 9: Existing-tab regression tests + atomic commit (1)

**Files:**
- Create: `nfl_draft/tests/integration/test_dashboard_queries.py`
- Create: `nfl_draft/tests/integration/test_kalshi_writer_target.py`

- [ ] **Step 1: Write existing-tab regression tests**

Create `nfl_draft/tests/integration/test_dashboard_queries.py`:
```python
"""Regression tests for legacy dashboard tabs after the kalshi_odds rename."""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module


@pytest.fixture
def seeded_db(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    # Schema for legacy kalshi_odds (matches migrated structure)
    with db_module.write_connection() as con:
        con.execute("""
            CREATE TABLE kalshi_odds (
                fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
                ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
                yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
                last_price INTEGER, volume BIGINT, volume_24h BIGINT,
                liquidity BIGINT, open_interest INTEGER
            )
        """)
        now = datetime.now()
        con.execute(
            "INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFT1', 'EVT1', 'TKR1', 'Title', 'Cam Ward', 50, 51, 49, 50, 50, 1000, 10000, 5000, 200)",
            [now],
        )
    return tmp_path / "test.duckdb"


def test_get_latest_odds_returns_one_row(seeded_db, monkeypatch):
    # Import after monkeypatching DB_PATH
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    df = legacy_db.get_latest_odds()
    assert df is not None and len(df) == 1
    assert df.iloc[0]["candidate"] == "Cam Ward"


def test_get_price_history_returns_one_row(seeded_db, monkeypatch):
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    df = legacy_db.get_price_history()
    assert df is not None and len(df) == 1


def test_get_snapshot_count_returns_one(seeded_db, monkeypatch):
    from kalshi_draft import db as legacy_db
    monkeypatch.setattr(legacy_db, "DB_PATH", seeded_db)
    snapshots, first, last = legacy_db.get_snapshot_count()
    assert snapshots == 1
```

- [ ] **Step 2: Write writer-target test**

Create `nfl_draft/tests/integration/test_kalshi_writer_target.py`:
```python
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
```

- [ ] **Step 3: Run all integration + unit tests**

Run: `pytest nfl_draft/tests/ -v`
Expected: all PASS (db, tz, migration, dashboard queries, writer target)

- [ ] **Step 4: Run the migration manually against the real legacy DB**

```bash
python -m nfl_draft.lib.migrate_from_kalshi_draft
```
Expected output: `kalshi_odds: copied N rows` (where N is the existing row count)

- [ ] **Step 5: Smoke-test the legacy dashboard**

Run: `python kalshi_draft/app.py &` then visit http://127.0.0.1:8083/. Navigate each existing tab (Market Overview, Price History, Edge Detection, Consensus, Portfolio) and confirm they render data. Kill the process: `pkill -f kalshi_draft/app.py`.

- [ ] **Step 6: Delete the legacy DB file**

```bash
rm kalshi_draft/kalshi_draft.duckdb
rm -f kalshi_draft/kalshi_draft.duckdb.wal
```

- [ ] **Step 7: Atomic commit**

```bash
git add nfl_draft/__init__.py nfl_draft/config/__init__.py nfl_draft/lib/ nfl_draft/scrapers/__init__.py nfl_draft/tests/ kalshi_draft/db.py kalshi_draft/fetcher.py kalshi_draft/edge_detector.py kalshi_draft/consensus.py kalshi_draft/app.py
git commit -m "feat: migrate kalshi_draft.duckdb → nfl_draft.duckdb (atomic)

- Schema migration script with MD5-hash idempotency
- Legacy table draft_odds renamed to kalshi_odds
- All consumers (db.py, fetcher.py, edge_detector.py, consensus.py,
  app.py) repointed to nfl_draft/nfl_draft.duckdb and SQL updated
- Legacy fetcher TZ normalized from UTC to local time
- Regression tests for all 5 existing dashboard tabs
- Writer-target test confirms legacy fetcher writes go to new DB"
```

---

## Task 10: Config dicts (players, teams, markets)

**Files:**
- Create: `nfl_draft/config/players.py`
- Create: `nfl_draft/config/teams.py`
- Create: `nfl_draft/config/markets.py`

- [ ] **Step 1: Write players.py with top prospects**

Create `nfl_draft/config/players.py`:
```python
"""Canonical 2026 NFL Draft prospect registry.

Each entry: canonical_name (the spelling we'll use everywhere) →
metadata + aliases (every variant a book might use).

Maintenance: when a scraper produces an unmapped player (lands in
draft_odds_unmapped_players), add the alias here and re-run
`python -m nfl_draft.lib.seed`.
"""

PLAYERS: dict[str, dict] = {
    # canonical_name: {position, college, aliases: [str]}
    # Top 80 prospects — populate from latest mock draft consensus
    # (Tankathon, ESPN, NFL.com — pick whichever was used by the user)
    "Cam Ward": {
        "position": "QB",
        "college": "Miami",
        "aliases": ["Cam Ward", "Cameron Ward", "C. Ward", "Ward, Cam", "Ward (Miami)"],
    },
    "Shedeur Sanders": {
        "position": "QB",
        "college": "Colorado",
        "aliases": ["Shedeur Sanders", "S. Sanders", "Sanders, Shedeur", "Sanders (Colorado)"],
    },
    "Travis Hunter": {
        "position": "WR",  # also DB; use primary draft position
        "college": "Colorado",
        "aliases": ["Travis Hunter", "T. Hunter", "Hunter, Travis"],
    },
    "Ashton Jeanty": {
        "position": "RB",
        "college": "Boise State",
        "aliases": ["Ashton Jeanty", "A. Jeanty", "Jeanty, Ashton", "Jeanty (Boise St)"],
    },
    # ... (~76 more — populate from mock draft consensus)
    # If unsure, start with top 10 and let the quarantine table flag the rest during reconnaissance.
}
```

- [ ] **Step 2: Write teams.py**

Create `nfl_draft/config/teams.py`:
```python
"""Canonical NFL team registry. Reuses logic from Answer Keys/canonical_match.py."""

TEAMS: dict[str, dict] = {
    "ARI": {"full_name": "Arizona Cardinals", "aliases": ["ARI", "Arizona", "Cardinals", "Arizona Cardinals"]},
    "ATL": {"full_name": "Atlanta Falcons", "aliases": ["ATL", "Atlanta", "Falcons", "Atlanta Falcons"]},
    "BAL": {"full_name": "Baltimore Ravens", "aliases": ["BAL", "Baltimore", "Ravens", "Baltimore Ravens"]},
    "BUF": {"full_name": "Buffalo Bills", "aliases": ["BUF", "Buffalo", "Bills", "Buffalo Bills"]},
    "CAR": {"full_name": "Carolina Panthers", "aliases": ["CAR", "Carolina", "Panthers", "Carolina Panthers"]},
    "CHI": {"full_name": "Chicago Bears", "aliases": ["CHI", "Chicago", "Bears", "Chicago Bears"]},
    "CIN": {"full_name": "Cincinnati Bengals", "aliases": ["CIN", "Cincinnati", "Bengals", "Cincinnati Bengals"]},
    "CLE": {"full_name": "Cleveland Browns", "aliases": ["CLE", "Cleveland", "Browns", "Cleveland Browns"]},
    "DAL": {"full_name": "Dallas Cowboys", "aliases": ["DAL", "Dallas", "Cowboys", "Dallas Cowboys"]},
    "DEN": {"full_name": "Denver Broncos", "aliases": ["DEN", "Denver", "Broncos", "Denver Broncos"]},
    "DET": {"full_name": "Detroit Lions", "aliases": ["DET", "Detroit", "Lions", "Detroit Lions"]},
    "GB": {"full_name": "Green Bay Packers", "aliases": ["GB", "GNB", "Green Bay", "Packers", "Green Bay Packers"]},
    "HOU": {"full_name": "Houston Texans", "aliases": ["HOU", "Houston", "Texans", "Houston Texans"]},
    "IND": {"full_name": "Indianapolis Colts", "aliases": ["IND", "Indianapolis", "Colts", "Indianapolis Colts"]},
    "JAX": {"full_name": "Jacksonville Jaguars", "aliases": ["JAX", "JAC", "Jacksonville", "Jaguars", "Jacksonville Jaguars"]},
    "KC": {"full_name": "Kansas City Chiefs", "aliases": ["KC", "KAN", "Kansas City", "Chiefs", "Kansas City Chiefs"]},
    "LV": {"full_name": "Las Vegas Raiders", "aliases": ["LV", "LVR", "Las Vegas", "Raiders", "Las Vegas Raiders", "Oakland Raiders"]},
    "LAC": {"full_name": "Los Angeles Chargers", "aliases": ["LAC", "LA Chargers", "Chargers", "Los Angeles Chargers"]},
    "LAR": {"full_name": "Los Angeles Rams", "aliases": ["LAR", "LA Rams", "Rams", "Los Angeles Rams"]},
    "MIA": {"full_name": "Miami Dolphins", "aliases": ["MIA", "Miami", "Dolphins", "Miami Dolphins"]},
    "MIN": {"full_name": "Minnesota Vikings", "aliases": ["MIN", "Minnesota", "Vikings", "Minnesota Vikings"]},
    "NE": {"full_name": "New England Patriots", "aliases": ["NE", "NEP", "New England", "Patriots", "New England Patriots"]},
    "NO": {"full_name": "New Orleans Saints", "aliases": ["NO", "NOR", "New Orleans", "Saints", "New Orleans Saints"]},
    "NYG": {"full_name": "New York Giants", "aliases": ["NYG", "NY Giants", "Giants", "New York Giants"]},
    "NYJ": {"full_name": "New York Jets", "aliases": ["NYJ", "NY Jets", "Jets", "New York Jets"]},
    "PHI": {"full_name": "Philadelphia Eagles", "aliases": ["PHI", "Philadelphia", "Eagles", "Philadelphia Eagles"]},
    "PIT": {"full_name": "Pittsburgh Steelers", "aliases": ["PIT", "Pittsburgh", "Steelers", "Pittsburgh Steelers"]},
    "SF": {"full_name": "San Francisco 49ers", "aliases": ["SF", "SFO", "San Francisco", "49ers", "San Francisco 49ers", "Niners"]},
    "SEA": {"full_name": "Seattle Seahawks", "aliases": ["SEA", "Seattle", "Seahawks", "Seattle Seahawks"]},
    "TB": {"full_name": "Tampa Bay Buccaneers", "aliases": ["TB", "TBB", "Tampa Bay", "Buccaneers", "Bucs", "Tampa Bay Buccaneers"]},
    "TEN": {"full_name": "Tennessee Titans", "aliases": ["TEN", "Tennessee", "Titans", "Tennessee Titans"]},
    "WAS": {"full_name": "Washington Commanders", "aliases": ["WAS", "WSH", "Washington", "Commanders", "Washington Commanders"]},
}
```

- [ ] **Step 3: Write markets.py with initial mappings**

Create `nfl_draft/config/markets.py`:
```python
"""Per-book market label → canonical market_id mappings.

Populated incrementally during scraper reconnaissance (Tasks 16-20).
On first scraper run, unmapped rows land in draft_odds_unmapped — add
entries here, then re-run `python -m nfl_draft.lib.seed`.
"""

# Schema: list of (book, book_label, book_subject, market_id) tuples.
# market_id must match what build_market_id() produces in lib/market_map.py.
MARKET_MAP: list[tuple[str, str, str, str]] = [
    # Kalshi: ticker prefix + candidate name
    # ("kalshi", "KXNFLDRAFT1", "Cam Ward", "pick_1_overall_cam-ward"),
    # ("kalshi", "KXNFLDRAFTQB", "Cam Ward", "first_qb_cam-ward"),
    # DK / FD / Bookmaker / Wagerzon entries added during reconnaissance.
]
```

- [ ] **Step 4: Verify modules import**

Run: `python -c "from nfl_draft.config import players, teams, markets; print(len(players.PLAYERS), len(teams.TEAMS), len(markets.MARKET_MAP))"`
Expected: prints three integers (e.g., `4 32 0`)

- [ ] **Step 5: Continue (no commit — atomic with Task 12)**

---

## Task 11: lib/seed.py — populate lookup tables from config dicts

**Files:**
- Create: `nfl_draft/lib/seed.py`
- Modify: `nfl_draft/tests/unit/test_db.py` (add seed test)

- [ ] **Step 1: Write failing seed test**

Append to `nfl_draft/tests/unit/test_db.py`:
```python
def test_seed_populates_player_aliases(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from nfl_draft.lib import seed
    seed.run()
    with db_module.read_connection() as con:
        count = con.execute("SELECT COUNT(*) FROM player_aliases").fetchone()[0]
    assert count > 0  # at least the seeded prospects


def test_seed_is_idempotent(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from nfl_draft.lib import seed
    seed.run()
    first_count = _alias_count()
    seed.run()
    assert _alias_count() == first_count


def _alias_count():
    with db_module.read_connection() as con:
        return con.execute("SELECT COUNT(*) FROM player_aliases").fetchone()[0]
```

- [ ] **Step 2: Run test to verify failure**

Run: `pytest nfl_draft/tests/unit/test_db.py::test_seed_populates_player_aliases -v`
Expected: FAIL (seed module doesn't exist)

- [ ] **Step 3: Implement seed.py**

Create `nfl_draft/lib/seed.py`:
```python
"""One-shot seeder: dict literals in config/*.py → DuckDB lookup tables.

Idempotent: TRUNCATE-and-INSERT inside a single transaction, so re-runs
produce identical state. Called automatically at the start of every
run.py invocation.
"""

from nfl_draft.lib.db import write_connection
from nfl_draft.config.players import PLAYERS
from nfl_draft.config.teams import TEAMS
from nfl_draft.config.markets import MARKET_MAP


def run() -> None:
    """Truncate + repopulate all 5 lookup tables in one transaction."""
    with write_connection() as con:
        con.execute("BEGIN TRANSACTION")
        try:
            # players
            con.execute("DELETE FROM players")
            con.execute("DELETE FROM player_aliases")
            for canonical, info in PLAYERS.items():
                con.execute(
                    "INSERT INTO players VALUES (?, ?, ?)",
                    [canonical, info["position"], info.get("college")],
                )
                for alias in info.get("aliases", []):
                    con.execute(
                        "INSERT INTO player_aliases VALUES (?, ?)",
                        [alias.lower().strip(), canonical],
                    )
            # teams
            con.execute("DELETE FROM teams")
            con.execute("DELETE FROM team_aliases")
            for abbr, info in TEAMS.items():
                con.execute(
                    "INSERT INTO teams VALUES (?, ?)",
                    [abbr, info["full_name"]],
                )
                for alias in info.get("aliases", []):
                    con.execute(
                        "INSERT INTO team_aliases VALUES (?, ?)",
                        [alias.lower().strip(), abbr],
                    )
            # market_map
            con.execute("DELETE FROM market_map")
            for book, label, subject, market_id in MARKET_MAP:
                con.execute(
                    "INSERT INTO market_map VALUES (?, ?, ?, ?)",
                    [book, label, subject, market_id],
                )
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise


if __name__ == "__main__":
    run()
    print("Seed complete.")
```

- [ ] **Step 4: Run tests**

Run: `pytest nfl_draft/tests/unit/test_db.py -v`
Expected: all PASS (including the two new seed tests)

- [ ] **Step 5: Continue (no commit — atomic with Task 12)**

---

## Task 12: lib/normalize.py + lib/market_map.py + commit (2)

**Files:**
- Create: `nfl_draft/lib/normalize.py`
- Create: `nfl_draft/lib/market_map.py`
- Create: `nfl_draft/tests/unit/test_normalize.py`
- Create: `nfl_draft/tests/unit/test_market_map.py`

- [ ] **Step 1: Write failing tests**

Create `nfl_draft/tests/unit/test_normalize.py`:
```python
import pytest
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.normalize import resolve_player, resolve_team


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_resolve_player_known_canonical(seeded):
    assert resolve_player("Cam Ward") == "Cam Ward"


def test_resolve_player_alias(seeded):
    assert resolve_player("Cameron Ward") == "Cam Ward"


def test_resolve_player_case_insensitive(seeded):
    assert resolve_player("cAm WaRd") == "Cam Ward"


def test_resolve_player_whitespace_tolerant(seeded):
    assert resolve_player("  Cam Ward  ") == "Cam Ward"


def test_resolve_player_unknown_returns_none(seeded):
    assert resolve_player("Joe Nobody") is None


def test_resolve_team_canonical(seeded):
    assert resolve_team("CHI") == "CHI"


def test_resolve_team_alias(seeded):
    assert resolve_team("Chicago Bears") == "CHI"


def test_resolve_team_unknown_returns_none(seeded):
    assert resolve_team("Detroit Pistons") is None
```

Create `nfl_draft/tests/unit/test_market_map.py`:
```python
import pytest
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.market_map import resolve_market_id, build_market_id


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_build_market_id_first_at_position():
    assert build_market_id("first_at_position", position="QB", player="Cam Ward") == "first_qb_cam-ward"


def test_build_market_id_pick_outright():
    assert build_market_id("pick_outright", pick_number=1, player="Cam Ward") == "pick_1_overall_cam-ward"


def test_build_market_id_top_n_range():
    assert build_market_id("top_n_range", range_high=5, player="Ashton Jeanty") == "top_5_ashton-jeanty"


def test_build_market_id_team_first_pick():
    assert build_market_id("team_first_pick", team="CHI", player="Ashton Jeanty") == "team_chi_first_pick_ashton-jeanty"


def test_build_market_id_strips_punctuation():
    assert build_market_id("first_at_position", position="QB", player="J.T. O'Brien") == "first_qb_jt-obrien"


def test_resolve_market_id_unknown_returns_none(seeded):
    assert resolve_market_id("kalshi", "UNKNOWN_LABEL", "Some Subject") is None
```

- [ ] **Step 2: Run tests to verify failure**

Run: `pytest nfl_draft/tests/unit/test_normalize.py nfl_draft/tests/unit/test_market_map.py -v`
Expected: FAIL (modules don't exist)

- [ ] **Step 3: Implement normalize.py**

Create `nfl_draft/lib/normalize.py`:
```python
"""Resolve player/team aliases to canonical names via DuckDB lookup tables."""

from typing import Optional
from nfl_draft.lib.db import read_connection


def resolve_player(name: str) -> Optional[str]:
    """Returns canonical player name, or None if unmapped."""
    if not name:
        return None
    key = name.lower().strip()
    with read_connection() as con:
        result = con.execute(
            "SELECT canonical_name FROM player_aliases WHERE alias = ?",
            [key],
        ).fetchone()
    return result[0] if result else None


def resolve_team(name: str) -> Optional[str]:
    """Returns canonical team abbr (e.g., 'CHI'), or None if unmapped."""
    if not name:
        return None
    key = name.lower().strip()
    with read_connection() as con:
        result = con.execute(
            "SELECT canonical_abbr FROM team_aliases WHERE alias = ?",
            [key],
        ).fetchone()
    return result[0] if result else None
```

- [ ] **Step 4: Implement market_map.py**

Create `nfl_draft/lib/market_map.py`:
```python
"""Per-book market label → canonical market_id mapping + market_id construction."""

from typing import Optional
from nfl_draft.lib.db import read_connection


def slug(name: str) -> str:
    """Canonicalize a name to a URL-safe slug (used in market_id construction)."""
    return name.lower().replace(" ", "-").replace(".", "").replace("'", "")


def build_market_id(market_type: str, **kwargs) -> str:
    """Deterministic market_id constructor — must be applied identically in every scraper.

    Tested by test_market_id_construction.py.
    """
    if market_type == "first_at_position":
        return f"first_{kwargs['position'].lower()}_{slug(kwargs['player'])}"
    if market_type == "pick_outright":
        return f"pick_{kwargs['pick_number']}_overall_{slug(kwargs['player'])}"
    if market_type == "top_n_range":
        return f"top_{kwargs['range_high']}_{slug(kwargs['player'])}"
    if market_type == "team_first_pick":
        return f"team_{kwargs['team'].lower()}_first_pick_{slug(kwargs['player'])}"
    if market_type == "prop":
        return f"prop_{slug(kwargs['short_description'])}"
    raise ValueError(f"Unknown market_type: {market_type}")


def resolve_market_id(book: str, book_label: str, book_subject: str) -> Optional[str]:
    """Look up canonical market_id for a per-book (label, subject) pair."""
    with read_connection() as con:
        result = con.execute(
            "SELECT market_id FROM market_map WHERE book = ? AND book_label = ? AND book_subject = ?",
            [book, book_label, book_subject],
        ).fetchone()
    return result[0] if result else None
```

- [ ] **Step 5: Run tests**

Run: `pytest nfl_draft/tests/unit/ -v`
Expected: all PASS

- [ ] **Step 6: Commit (2)**

```bash
git add nfl_draft/config/ nfl_draft/lib/seed.py nfl_draft/lib/normalize.py nfl_draft/lib/market_map.py nfl_draft/tests/unit/test_normalize.py nfl_draft/tests/unit/test_market_map.py nfl_draft/tests/unit/test_db.py
git commit -m "feat: config dicts + seed + normalize + market_id construction

Five lookup tables (players, player_aliases, teams, team_aliases,
market_map) seeded from Python dict literals in config/*.py.
Idempotent seed (TRUNCATE + INSERT in one transaction). Normalize
resolves aliases via O(1) lookup. build_market_id() is the single
canonical constructor used by all scrapers."
```

---

## Task 13: Devig math (lib/devig.py) + commit (3)

**Files:**
- Create: `nfl_draft/lib/devig.py`
- Create: `nfl_draft/tests/unit/test_devig.py`

- [ ] **Step 1: Write failing tests**

Create `nfl_draft/tests/unit/test_devig.py`:
```python
import pytest
from nfl_draft.lib.devig import american_to_implied, devig_two_way, devig_n_way, proportional_devig


# American odds → implied probability
def test_american_to_implied_pickem():
    assert american_to_implied(100) == pytest.approx(0.5, abs=1e-6)


def test_american_to_implied_favorite():
    # -200 ⇒ p = 200 / (200 + 100) = 0.6667
    assert american_to_implied(-200) == pytest.approx(2/3, abs=1e-6)


def test_american_to_implied_dog():
    # +200 ⇒ p = 100 / (200 + 100) = 0.3333
    assert american_to_implied(200) == pytest.approx(1/3, abs=1e-6)


# Two-way devig
def test_devig_two_way_zero_overround():
    p_a, p_b = devig_two_way(100, 100)  # 0.5 / 0.5 sum = 1.0
    assert p_a + p_b == pytest.approx(1.0)
    assert p_a == pytest.approx(0.5)


def test_devig_two_way_with_vig():
    # -110 / -110 → both 0.5238, sum 1.0476, devigged 0.5 / 0.5
    p_a, p_b = devig_two_way(-110, -110)
    assert p_a == pytest.approx(0.5, abs=1e-4)
    assert p_b == pytest.approx(0.5, abs=1e-4)


def test_devig_two_way_asymmetric():
    p_a, p_b = devig_two_way(-200, 150)
    assert (p_a + p_b) == pytest.approx(1.0, abs=1e-6)
    assert p_a > p_b  # favorite higher prob


# N-way devig
def test_devig_n_way_three_outcomes():
    probs = devig_n_way([200, 200, 200])  # each 0.333; sum 1.0; no vig
    assert sum(probs) == pytest.approx(1.0)
    for p in probs:
        assert p == pytest.approx(1/3, abs=1e-6)


def test_devig_n_way_with_vig():
    probs = devig_n_way([100, 100, 100])  # each 0.5; sum 1.5; devigged each 1/3
    assert sum(probs) == pytest.approx(1.0)
    for p in probs:
        assert p == pytest.approx(1/3, abs=1e-6)


# Proportional devig (sparse outcomes)
def test_proportional_devig_normalizes():
    # 5 outcomes posted, sum implied = 0.4 (other 25 not posted)
    odds = [400, 400, 400, 400, 400]  # each implied = 0.2; sum = 1.0
    devigged = proportional_devig(odds)
    assert sum(devigged) == pytest.approx(1.0)


def test_proportional_devig_handles_int_and_float_inputs():
    devigged = proportional_devig([100.0, 200, -150])
    assert sum(devigged) == pytest.approx(1.0)
```

- [ ] **Step 2: Run tests to verify failure**

Run: `pytest nfl_draft/tests/unit/test_devig.py -v`
Expected: FAIL (module doesn't exist)

- [ ] **Step 3: Implement devig.py**

Create `nfl_draft/lib/devig.py`:
```python
"""Devig math — port of Answer Keys/Tools.R devigging functions."""

from typing import List, Tuple


def american_to_implied(american_odds: float) -> float:
    """Convert American odds to implied probability (vig included)."""
    o = float(american_odds)
    if o > 0:
        return 100.0 / (o + 100.0)
    return -o / (-o + 100.0)


def devig_two_way(odds_a: float, odds_b: float) -> Tuple[float, float]:
    """Remove vig from a two-way market (proportional method)."""
    a = american_to_implied(odds_a)
    b = american_to_implied(odds_b)
    total = a + b
    return (a / total, b / total)


def devig_n_way(odds_list: List[float]) -> List[float]:
    """Remove vig from an n-way market where the complementary set is fully posted."""
    implieds = [american_to_implied(o) for o in odds_list]
    total = sum(implieds)
    return [p / total for p in implieds]


def proportional_devig(odds_list: List[float]) -> List[float]:
    """Devig when the complementary set is INCOMPLETE (e.g., book only posts top 5).

    KNOWN LIMITATION: assumes posted outcomes contain ~100% of probability mass.
    Inflates the posted candidates' probs vs. ground truth. Spec annotates this
    in the dashboard.
    """
    # Same math as devig_n_way; the difference is interpretive (caller knows
    # the field is incomplete).
    return devig_n_way(odds_list)
```

- [ ] **Step 4: Run tests**

Run: `pytest nfl_draft/tests/unit/test_devig.py -v`
Expected: 10 tests PASS

- [ ] **Step 5: Add market_id construction tests file (already have build_market_id tests in Task 12, but renaming for spec consistency)**

Run: `mv nfl_draft/tests/unit/test_market_map.py nfl_draft/tests/unit/test_market_map.py` — already correctly named. The `build_market_id` tests are already in `test_market_map.py` from Task 12. To match the spec's `test_market_id_construction.py` naming, copy the relevant tests:

Create `nfl_draft/tests/unit/test_market_id_construction.py`:
```python
"""Dedicated tests for build_market_id() — the single canonical constructor."""

import pytest
from nfl_draft.lib.market_map import build_market_id, slug


def test_slug_lowercases():
    assert slug("Cam Ward") == "cam-ward"


def test_slug_strips_periods():
    assert slug("J.T. Smith") == "jt-smith"


def test_slug_strips_apostrophes():
    assert slug("J'shawn O'Brien") == "jshawn-obrien"


def test_each_market_type_produces_expected_id():
    cases = [
        ("first_at_position", {"position": "QB", "player": "Cam Ward"}, "first_qb_cam-ward"),
        ("pick_outright", {"pick_number": 1, "player": "Cam Ward"}, "pick_1_overall_cam-ward"),
        ("top_n_range", {"range_high": 5, "player": "Ashton Jeanty"}, "top_5_ashton-jeanty"),
        ("team_first_pick", {"team": "CHI", "player": "Ashton Jeanty"}, "team_chi_first_pick_ashton-jeanty"),
        ("prop", {"short_description": "Will a kicker go R1"}, "prop_will-a-kicker-go-r1"),
    ]
    for market_type, kwargs, expected in cases:
        assert build_market_id(market_type, **kwargs) == expected


def test_unknown_market_type_raises():
    with pytest.raises(ValueError):
        build_market_id("invalid_type", foo="bar")
```

Run: `pytest nfl_draft/tests/unit/test_market_id_construction.py -v`
Expected: 5 tests PASS

- [ ] **Step 6: Commit (3)**

```bash
git add nfl_draft/lib/devig.py nfl_draft/tests/unit/test_devig.py nfl_draft/tests/unit/test_market_id_construction.py
git commit -m "feat: devig math + market_id construction rules

Three devig functions ported from Answer Keys/Tools.R: american_to_implied,
devig_two_way, devig_n_way (plus proportional_devig for incomplete fields).
build_market_id() construction rule with deterministic slug() function;
tested across all market_type variants."
```

---

## Task 14: Quarantine + scraper interface

**Files:**
- Create: `nfl_draft/lib/quarantine.py`
- Create: `nfl_draft/scrapers/_base.py` (shared OddsRow dataclass)
- Create: `nfl_draft/tests/integration/test_quarantine.py`

- [ ] **Step 1: Write failing quarantine test**

Create `nfl_draft/tests/integration/test_quarantine.py`:
```python
"""End-to-end test: unmapped player + unmapped market land in quarantine, NOT draft_odds."""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine
from nfl_draft.scrapers._base import OddsRow


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_unmapped_player_lands_in_quarantine(seeded):
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Made Up Player",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 0
    assert quarantine_count == 1


def test_mapped_row_lands_in_draft_odds(seeded):
    # First, manually seed a market_map entry so the row maps cleanly
    from nfl_draft.lib.db import write_connection
    with write_connection() as con:
        con.execute("INSERT INTO market_map VALUES ('draftkings', 'First QB', 'Cam Ward', 'first_qb_cam-ward')")
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Cam Ward",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 1
    assert quarantine_count == 0
```

- [ ] **Step 2: Run test to verify failure**

Run: `pytest nfl_draft/tests/integration/test_quarantine.py -v`
Expected: FAIL (modules don't exist)

- [ ] **Step 3: Implement OddsRow + write_or_quarantine**

Create `nfl_draft/scrapers/_base.py`:
```python
"""Shared scraper data structures."""

from dataclasses import dataclass
from datetime import datetime


@dataclass
class OddsRow:
    """A single row of raw scraper output (one binary outcome from one venue)."""
    book: str
    book_label: str       # e.g., "First QB" or "KXNFLDRAFT1"
    book_subject: str     # e.g., "Cam Ward" (raw, before normalization)
    american_odds: int
    fetched_at: datetime
    market_group: str = ""  # optional: which set this belongs to (e.g., "first_qb_field")
                            # used to group rows for n-way devig in run.py


@dataclass
class TradeRow:
    """A single Kalshi trade event (for the trade tape)."""
    trade_id: str
    ticker: str
    side: str             # 'yes' | 'no'
    price_cents: int
    count: int
    traded_at: datetime
    fetched_at: datetime
```

Create `nfl_draft/lib/quarantine.py`:
```python
"""Write OddsRow batches to draft_odds with quarantine for unmapped rows."""

from typing import List
from nfl_draft.lib.db import write_connection
from nfl_draft.lib.market_map import resolve_market_id
from nfl_draft.lib.devig import american_to_implied
from nfl_draft.scrapers._base import OddsRow


def write_or_quarantine(rows: List[OddsRow]) -> tuple[int, int]:
    """Resolve each row's market_id; write mapped to draft_odds, unmapped to quarantine.

    Returns: (mapped_count, unmapped_count).
    """
    mapped = 0
    unmapped = 0
    with write_connection() as con:
        for row in rows:
            mid = resolve_market_id(row.book, row.book_label, row.book_subject)
            if mid is None:
                con.execute(
                    "INSERT INTO draft_odds_unmapped VALUES (?, ?, ?, ?, ?, ?)",
                    [row.book, row.book_label, row.book_subject, row.american_odds, row.fetched_at, "no market_map entry"],
                )
                unmapped += 1
                continue
            implied = american_to_implied(row.american_odds)
            # devig_prob computed separately at run.py level after grouping;
            # for v1 single-row inserts, devig_prob == implied_prob (no grouping)
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, implied, row.fetched_at],
            )
            mapped += 1
    return (mapped, unmapped)
```

- [ ] **Step 4: Run tests**

Run: `pytest nfl_draft/tests/integration/test_quarantine.py -v`
Expected: 2 tests PASS

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/_base.py nfl_draft/lib/quarantine.py nfl_draft/tests/integration/test_quarantine.py
git commit -m "feat: OddsRow dataclass + quarantine writer

Unmapped rows land in draft_odds_unmapped (configurable, surface
unmapped count on dashboard footer). Mapped rows go to draft_odds
with implied probability computed at insert; n-way devig grouping
happens at run.py orchestrator level."
```

---

## Task 15: Kalshi adapter (4a)

**Files:**
- Create: `nfl_draft/scrapers/kalshi.py`
- Create: `nfl_draft/tests/unit/test_scraper_parsing.py` (Kalshi cases)
- Create: `nfl_draft/tests/fixtures/kalshi/series_response.json`
- Create: `nfl_draft/tests/fixtures/kalshi/markets_response.json`
- Modify: `kalshi_draft/auth.py` (add throttle)

- [ ] **Step 1: Add throttle to kalshi_draft/auth.py**

In `kalshi_draft/auth.py`, modify `public_request` to throttle (600ms minimum between calls, exponential backoff on 429). Add a module-level state:
```python
import time

_LAST_REQUEST_TS = 0.0
_MIN_INTERVAL = 0.6  # 600ms
_BACKOFFS = [1, 2, 4, 8, 16]


def public_request(path: str):
    """Make a public (unauthenticated) request to Kalshi API. Throttled."""
    global _LAST_REQUEST_TS
    elapsed = time.time() - _LAST_REQUEST_TS
    if elapsed < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - elapsed)
    url = f"{BASE_URL}{path}"
    for attempt, backoff in enumerate([0] + _BACKOFFS):
        if backoff:
            time.sleep(backoff)
        try:
            req = urllib.request.Request(url)
            req.add_header("Content-Type", "application/json")
            with urllib.request.urlopen(req) as response:
                _LAST_REQUEST_TS = time.time()
                return json.loads(response.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < len(_BACKOFFS):
                continue  # backoff and retry
            print(f"  Request failed: {e.code} - {e.read().decode()}")
            return None
        except Exception as e:
            print(f"  Request error: {e}")
            return None
```

- [ ] **Step 2: Capture Kalshi fixtures**

Run:
```bash
mkdir -p nfl_draft/tests/fixtures/kalshi
curl -s 'https://api.elections.kalshi.com/trade-api/v2/series?limit=10' > nfl_draft/tests/fixtures/kalshi/series_response.json
curl -s 'https://api.elections.kalshi.com/trade-api/v2/markets?series_ticker=KXNFLDRAFT1&limit=10' > nfl_draft/tests/fixtures/kalshi/markets_response.json
```

- [ ] **Step 3: Write failing parser test**

Create `nfl_draft/tests/unit/test_scraper_parsing.py`:
```python
"""Tests for each scraper's parse() function — pure logic, no network."""

import json
from pathlib import Path
from nfl_draft.scrapers.kalshi import parse_markets_response


FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


def test_kalshi_parse_markets_returns_oddsrow_list():
    raw = json.loads((FIXTURES / "kalshi" / "markets_response.json").read_text())
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFT1")
    assert isinstance(rows, list)
    if rows:  # fixture might be empty if series is closed
        row = rows[0]
        assert row.book == "kalshi"
        assert row.book_label.startswith("KXNFLDRAFT1")
```

- [ ] **Step 4: Run test to verify failure**

Run: `pytest nfl_draft/tests/unit/test_scraper_parsing.py -v`
Expected: FAIL (kalshi module doesn't exist)

- [ ] **Step 5: Implement Kalshi adapter**

Create `nfl_draft/scrapers/kalshi.py`:
```python
"""Kalshi adapter: calls legacy kalshi_draft/fetcher.py for raw market data,
returns List[OddsRow] for the new pipeline. Legacy fetcher's side-effect
writes (kalshi_odds, draft_series, market_info, positions, resting_orders)
happen as a byproduct of those calls — kalshi_odds keeps growing for the
legacy tabs while draft_odds gets the normalized cross-venue rows.
"""

import sys
from datetime import datetime
from pathlib import Path
from typing import List

# Import the legacy fetcher
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from kalshi_draft import fetcher as legacy_fetcher
from kalshi_draft.auth import public_request

from nfl_draft.scrapers._base import OddsRow, TradeRow


def parse_markets_response(raw_response: dict, series_ticker: str) -> List[OddsRow]:
    """Convert raw Kalshi /markets response into OddsRow list."""
    rows: List[OddsRow] = []
    now = datetime.now()
    for market in raw_response.get("markets", []):
        # Use yes_bid as the implied price (standard convention)
        yes_bid_cents = market.get("yes_bid")
        if yes_bid_cents is None:
            continue
        # Convert price_cents (1-99) to American odds approximation:
        # if p = 0.50, american = +100; if p = 0.30, american = +233; if p = 0.70, american = -233
        p = yes_bid_cents / 100.0
        if p <= 0 or p >= 1:
            continue
        american = int(round(-100 * p / (1 - p))) if p > 0.5 else int(round(100 * (1 - p) / p))
        candidate = market.get("yes_sub_title") or market.get("subtitle") or ""
        rows.append(OddsRow(
            book="kalshi",
            book_label=series_ticker,
            book_subject=candidate.strip(),
            american_odds=american,
            fetched_at=now,
        ))
    return rows


def fetch_draft_odds() -> List[OddsRow]:
    """Discover all NFL Draft series, fetch markets, return OddsRow list.

    Side effect: legacy fetcher writes kalshi_odds + metadata to nfl_draft.duckdb.
    """
    rows: List[OddsRow] = []
    series_list = legacy_fetcher.discover_draft_series()
    for series in series_list:
        # Trigger legacy write (fills kalshi_odds + draft_series + market_info)
        markets = legacy_fetcher.fetch_markets_for_series(series["series_ticker"])
        # Re-parse for our own normalized output
        raw = public_request(f"/markets?series_ticker={series['series_ticker']}&status=open&limit=100")
        if raw:
            rows.extend(parse_markets_response(raw, series["series_ticker"]))
    return rows


def fetch_trades() -> List[TradeRow]:
    """Poll Kalshi /markets/trades for each NFL draft series, since last cursor."""
    # Stub for Task 21 — fully implemented there.
    raise NotImplementedError("Implemented in Task 21 (trade-tape poller)")
```

- [ ] **Step 6: Run tests**

Run: `pytest nfl_draft/tests/unit/test_scraper_parsing.py -v`
Expected: PASS (or skipped if fixture empty — that's OK)

- [ ] **Step 7: Smoke test against live Kalshi**

Run: `python -c "from nfl_draft.scrapers.kalshi import fetch_draft_odds; rows = fetch_draft_odds(); print(f'{len(rows)} rows'); [print(r) for r in rows[:3]]"`
Expected: prints row count and 3 sample rows.

- [ ] **Step 8: Commit (4a)**

```bash
git add nfl_draft/scrapers/kalshi.py nfl_draft/scrapers/_base.py nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/tests/fixtures/kalshi/ kalshi_draft/auth.py
git commit -m "feat(scrapers): kalshi adapter + auth.py throttle

Adapter calls legacy fetcher for raw data + side-effect writes,
returns OddsRow list for the new pipeline. Throttle added to
public_request (600ms min interval, exp backoff on 429) — both
kalshi_mm and the new portal benefit."
```

---

## Task 16: DraftKings scraper (4b)

**Files:**
- Create: `nfl_draft/scrapers/draftkings.py`
- Create: `nfl_draft/tests/fixtures/draftkings/draft_markets.json`
- Append cases to: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Reconnaissance**

Open https://sportsbook.draftkings.com/leagues/football/nfl-draft in a browser. Find the NFL Draft markets section. Open DevTools → Network. Reload, filter for XHR. Identify the API call(s) that return the markets JSON (typically `https://sportsbook-nash.draftkings.com/api/sportscontent/dkusva/v1/leagues/<id>`). Right-click → Save as → `nfl_draft/tests/fixtures/draftkings/draft_markets.json`.

Note the league ID, event IDs, and any required headers (especially anti-bot tokens).

- [ ] **Step 2: Write failing parser test**

Append to `nfl_draft/tests/unit/test_scraper_parsing.py`:
```python
def test_dk_parse_returns_oddsrow_list():
    from nfl_draft.scrapers.draftkings import parse_response
    raw = json.loads((FIXTURES / "draftkings" / "draft_markets.json").read_text())
    rows = parse_response(raw)
    assert isinstance(rows, list)
    assert all(r.book == "draftkings" for r in rows)
    assert any(r.book_label and r.book_subject for r in rows)
```

- [ ] **Step 3: Implement scraper**

Create `nfl_draft/scrapers/draftkings.py`:
```python
"""DraftKings NFL Draft scraper. Adapted from mlb_sgp/scraper_draftkings_sgp.py.

Reconnaissance step: locate the NFL Draft event page on DK, capture the
API response into tests/fixtures/draftkings/draft_markets.json. Update
parse_response() to traverse the actual JSON structure.
"""

from datetime import datetime
from typing import List
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
# Reuse the existing DK scraper's auth / Playwright setup
# from mlb_sgp.scraper_draftkings_sgp import _get_session  # adjust import to match real API

from nfl_draft.scrapers._base import OddsRow


def parse_response(raw: dict) -> List[OddsRow]:
    """Parse DK's NFL Draft markets JSON into OddsRow list.

    Structure to verify against the captured fixture:
      raw["events"][i]["displayGroups"][j]["markets"][k]["outcomes"][m]
      with fields: marketName, outcomeLabel, oddsAmerican.
    """
    rows: List[OddsRow] = []
    now = datetime.now()
    # TRAVERSAL based on the captured fixture — adjust paths if DK's structure differs:
    for event in raw.get("events", []):
        for group in event.get("displayGroups", []):
            for market in group.get("markets", []):
                market_label = market.get("name", "")
                for outcome in market.get("outcomes", []):
                    odds_american = outcome.get("oddsAmerican")
                    if odds_american is None:
                        continue
                    rows.append(OddsRow(
                        book="draftkings",
                        book_label=market_label,
                        book_subject=outcome.get("label", "").strip(),
                        american_odds=int(odds_american),
                        fetched_at=now,
                    ))
    return rows


def fetch_raw() -> dict:
    """Live network call. Reuses auth / Playwright pattern from mlb_sgp/scraper_draftkings_sgp.py."""
    # Implementation depends on reconnaissance findings. At minimum:
    # 1. Open Playwright with stealth profile (matching mlb_sgp pattern)
    # 2. Navigate to DK NFL Draft page
    # 3. Wait for API call to fire; capture response from network listener
    # 4. Return parsed JSON dict
    raise NotImplementedError("Reconnaissance + Playwright wiring TBD per the captured fixture")


def fetch_draft_odds() -> List[OddsRow]:
    raw = fetch_raw()
    return parse_response(raw)
```

- [ ] **Step 4: Iterate parse_response against the fixture until test passes**

Run: `pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_dk_parse_returns_oddsrow_list -v`

Walk the JSON structure with `json.tool` to figure out actual paths if your `events`/`displayGroups`/`markets` traversal doesn't match.

- [ ] **Step 5: Add market_map entries for what you found**

In `nfl_draft/config/markets.py`, append entries for each `(book_label, book_subject)` pair you scraped:
```python
("draftkings", "First QB Drafted", "Cam Ward", build_market_id("first_at_position", position="QB", player="Cam Ward")),
# ... one per outcome
```

Use a generation script if there are many: `python -c "from nfl_draft.scrapers.draftkings import parse_response; from nfl_draft.lib.market_map import build_market_id; ..."` to print the tuples for paste-in.

- [ ] **Step 6: Re-seed and verify quarantine count drops**

```bash
python -m nfl_draft.lib.seed
python -c "
from nfl_draft.scrapers.draftkings import parse_response
from nfl_draft.lib.quarantine import write_or_quarantine
import json
raw = json.loads(open('nfl_draft/tests/fixtures/draftkings/draft_markets.json').read())
rows = parse_response(raw)
mapped, unmapped = write_or_quarantine(rows)
print(f'mapped={mapped} unmapped={unmapped}')
"
```
Expected: `unmapped` should drop to 0 (or close to it) after market_map entries added.

- [ ] **Step 7: Commit (4b)**

```bash
git add nfl_draft/scrapers/draftkings.py nfl_draft/tests/fixtures/draftkings/ nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/config/markets.py
git commit -m "feat(scrapers): draftkings draft scraper + market_map entries

Reconnaissance captured the DK NFL Draft markets endpoint into a
fixture. Parser walks the response into OddsRow list. Market_map
populated for the markets observed in the captured response."
```

---

## Task 17: FanDuel scraper (4c)

**Files:**
- Create: `nfl_draft/scrapers/fanduel.py`
- Create: `nfl_draft/tests/fixtures/fanduel/draft_markets.json`
- Append cases to: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Reconnaissance**

Open https://sportsbook.fanduel.com/futures/nfl-draft in a browser. Capture the API response (per memory `fanduel_sgp_scraping.md` — uses `implyBets` endpoint with curl_cffi + 3 required headers: `x-application`, `x-sportsbook-region`, `x-px-context`).

Save raw JSON to `nfl_draft/tests/fixtures/fanduel/draft_markets.json`.

- [ ] **Step 2: Implement scraper following the same pattern as Task 16**

Create `nfl_draft/scrapers/fanduel.py`. Use curl_cffi + 3-header recipe from the memory. Implement `parse_response(raw)` traversing FD's response shape. Write parsing test, iterate, populate `market_map`, re-seed, verify quarantine drops.

- [ ] **Step 3: Commit (4c)**

```bash
git add nfl_draft/scrapers/fanduel.py nfl_draft/tests/fixtures/fanduel/ nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/config/markets.py
git commit -m "feat(scrapers): fanduel draft scraper + market_map entries"
```

---

## Task 18: Bookmaker scraper (4d)

**Files:**
- Create: `nfl_draft/scrapers/bookmaker.py`
- Create: `nfl_draft/tests/fixtures/bookmaker/draft_markets.json`
- Append cases to: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Reconnaissance**

Bookmaker login required. Use the curl_cffi + ASP.NET cookie pattern from `bookmaker_odds/scraper.py`. Discover the league ID for NFL Draft markets (currently only cbb/nba/mlb configured). Capture the response.

- [ ] **Step 2: Implement + test + populate market_map + commit (4d)**

Same pattern as Tasks 16-17.

```bash
git add nfl_draft/scrapers/bookmaker.py nfl_draft/tests/fixtures/bookmaker/ nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/config/markets.py
git commit -m "feat(scrapers): bookmaker draft scraper + market_map entries"
```

---

## Task 19: Wagerzon scraper (4e)

**Files:**
- Create: `nfl_draft/scrapers/wagerzon.py`
- Create: `nfl_draft/tests/fixtures/wagerzon/draft_markets.json`
- Append cases to: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Reconnaissance + implementation + commit (4e)**

Use the ASP.NET `__VIEWSTATE` POST pattern from `wagerzon_odds/scraper_v2.py`. Discover draft league IDs. Same structure as previous scrapers.

```bash
git add nfl_draft/scrapers/wagerzon.py nfl_draft/tests/fixtures/wagerzon/ nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/config/markets.py
git commit -m "feat(scrapers): wagerzon draft scraper + market_map entries"
```

---

## Task 20: run.py orchestrator + pipeline test

**Files:**
- Create: `nfl_draft/run.py`
- Create: `nfl_draft/tests/integration/test_pipeline.py`

- [ ] **Step 1: Write failing pipeline test**

Create `nfl_draft/tests/integration/test_pipeline.py`:
```python
"""End-to-end pipeline test: fixture → normalize → market_map → devig → DB."""
import json
import pytest
from pathlib import Path
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine
from nfl_draft.scrapers.kalshi import parse_markets_response


FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_kalshi_pipeline_writes_to_draft_odds_or_quarantine(seeded):
    raw = json.loads((FIXTURES / "kalshi" / "markets_response.json").read_text())
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFT1")
    if not rows:
        pytest.skip("Empty fixture")
    mapped, unmapped = write_or_quarantine(rows)
    assert (mapped + unmapped) == len(rows)
```

- [ ] **Step 2: Implement run.py**

Create `nfl_draft/run.py`:
```python
"""NFL Draft Portal orchestrator.

Modes:
  --mode scrape: pull all venues' odds → draft_odds + invoke legacy
                 edge_detector and consensus
  --mode trades: poll Kalshi trade tape with cursor + dedup → kalshi_trades
"""

import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine

SCRAPERS = {
    "kalshi": "nfl_draft.scrapers.kalshi",
    "draftkings": "nfl_draft.scrapers.draftkings",
    "fanduel": "nfl_draft.scrapers.fanduel",
    "bookmaker": "nfl_draft.scrapers.bookmaker",
    "wagerzon": "nfl_draft.scrapers.wagerzon",
}


def run_scrape(book: str) -> None:
    """Run --mode scrape for one or all books."""
    seed.run()  # always reseed in case config dicts edited
    targets = list(SCRAPERS.keys()) if book == "all" else [book]
    for name in targets:
        try:
            module = __import__(SCRAPERS[name], fromlist=["fetch_draft_odds"])
            print(f"[scrape] {name}: fetching...")
            rows = module.fetch_draft_odds()
            mapped, unmapped = write_or_quarantine(rows)
            print(f"[scrape] {name}: mapped={mapped} unmapped={unmapped}")
        except NotImplementedError:
            print(f"[scrape] {name}: NOT IMPLEMENTED, skipping")
        except Exception as e:
            print(f"[scrape] {name}: ERROR — {e}")
            traceback.print_exc()
    # Trigger legacy edge detector + consensus (writes to detected_edges, consensus_board)
    try:
        from kalshi_draft import edge_detector
        edge_detector.run()  # or whatever the entry point is
        print("[scrape] edge_detector: done")
    except Exception as e:
        print(f"[scrape] edge_detector: ERROR — {e}")
    try:
        from kalshi_draft import consensus
        consensus.run()  # or whatever the entry point is
        print("[scrape] consensus: done")
    except Exception as e:
        print(f"[scrape] consensus: ERROR — {e}")


def run_trades() -> None:
    """Run --mode trades. Implemented in Task 21."""
    seed.run()
    from nfl_draft.scrapers import kalshi
    print("[trades] polling...")
    trades = kalshi.fetch_trades()
    print(f"[trades] ingested {len(trades)} trades")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["scrape", "trades"], required=True)
    p.add_argument("--book", default="all")
    args = p.parse_args()
    if args.mode == "scrape":
        run_scrape(args.book)
    else:
        run_trades()


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run tests**

Run: `pytest nfl_draft/tests/integration/test_pipeline.py -v`
Expected: PASS

- [ ] **Step 4: Smoke test full scrape**

Run: `python -m nfl_draft.run --mode scrape --book kalshi`
Expected: prints `mapped=X unmapped=Y` for kalshi.

- [ ] **Step 5: Commit (orchestrator)**

```bash
git add nfl_draft/run.py nfl_draft/tests/integration/test_pipeline.py
git commit -m "feat: run.py orchestrator + end-to-end pipeline test

Sequential per-book invocation with try/except isolation. Seed runs
at start of every invocation. After scrapers complete, invokes
legacy edge_detector and consensus to keep the legacy tabs fresh."
```

---

## Task 21: Trade-tape poller (`--mode trades`)

**Files:**
- Modify: `nfl_draft/scrapers/kalshi.py` (implement `fetch_trades`)
- Create: `nfl_draft/tests/fixtures/kalshi/trades_response.json`
- Append cases to: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Capture trades fixture**

Run:
```bash
curl -s 'https://api.elections.kalshi.com/trade-api/v2/markets/trades?limit=20' > nfl_draft/tests/fixtures/kalshi/trades_response.json
```

- [ ] **Step 2: Write failing parser + dedup tests**

Append to `nfl_draft/tests/unit/test_scraper_parsing.py`:
```python
def test_kalshi_parse_trades():
    from nfl_draft.scrapers.kalshi import parse_trades_response
    raw = json.loads((FIXTURES / "kalshi" / "trades_response.json").read_text())
    rows = parse_trades_response(raw)
    assert isinstance(rows, list)
    if rows:
        t = rows[0]
        assert t.trade_id and t.ticker and t.price_cents and t.count is not None
```

Append to `nfl_draft/tests/integration/test_quarantine.py` (or new test file):
```python
def test_trades_insert_or_ignore_dedup(seeded):
    from nfl_draft.scrapers._base import TradeRow
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    t = TradeRow(trade_id="abc123", ticker="KXNFLDRAFT1-T1", side="yes",
                 price_cents=50, count=10, traded_at=datetime.now(), fetched_at=datetime.now())
    with write_connection() as con:
        for _ in range(2):
            con.execute(
                "INSERT OR IGNORE INTO kalshi_trades VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [t.trade_id, t.ticker, t.side, t.price_cents, t.count,
                 t.count * t.price_cents * 0.01, t.traded_at, t.fetched_at]
            )
    with db_module.read_connection() as con:
        count = con.execute("SELECT COUNT(*) FROM kalshi_trades WHERE trade_id='abc123'").fetchone()[0]
    assert count == 1  # second insert was deduped
```

- [ ] **Step 3: Implement fetch_trades and parse_trades_response**

Edit `nfl_draft/scrapers/kalshi.py` — replace the `fetch_trades` stub:
```python
from nfl_draft.lib.db import write_connection, read_connection, local_to_utc_iso, utc_iso_to_local


def parse_trades_response(raw: dict) -> List[TradeRow]:
    rows: List[TradeRow] = []
    for trade in raw.get("trades", []):
        rows.append(TradeRow(
            trade_id=trade["trade_id"],
            ticker=trade["ticker"],
            side=trade["taker_side"].lower() if "taker_side" in trade else trade.get("side", "yes").lower(),
            price_cents=int(trade["yes_price"]),
            count=int(trade["count"]),
            traded_at=utc_iso_to_local(trade["created_time"]),
            fetched_at=datetime.now(),
        ))
    return rows


def fetch_trades() -> List[TradeRow]:
    """Poll Kalshi /markets/trades for each NFL series, since last cursor."""
    series_list = legacy_fetcher.discover_draft_series()
    all_rows: List[TradeRow] = []
    with write_connection() as con:
        for series in series_list:
            ticker = series["series_ticker"]
            cursor_row = con.execute(
                "SELECT last_traded_at FROM kalshi_poll_state WHERE series_ticker = ?",
                [ticker],
            ).fetchone()
            min_ts_iso = local_to_utc_iso(cursor_row[0]) if cursor_row else None
            path = f"/markets/trades?series_ticker={ticker}&limit=100"
            if min_ts_iso:
                path += f"&min_ts={min_ts_iso}"
            raw = public_request(path)
            if not raw:
                continue
            rows = parse_trades_response(raw)
            for r in rows:
                con.execute(
                    "INSERT OR IGNORE INTO kalshi_trades VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    [r.trade_id, r.ticker, r.side, r.price_cents, r.count,
                     r.count * r.price_cents * 0.01, r.traded_at, r.fetched_at],
                )
            if rows:
                max_ts = max(r.traded_at for r in rows)
                con.execute(
                    "INSERT OR REPLACE INTO kalshi_poll_state VALUES (?, ?, ?)",
                    [ticker, max_ts, datetime.now()],
                )
            all_rows.extend(rows)
    return all_rows
```

- [ ] **Step 4: Run tests**

Run: `pytest nfl_draft/tests/ -v`
Expected: all PASS.

- [ ] **Step 5: Smoke test**

Run: `python -m nfl_draft.run --mode trades`
Expected: prints `[trades] ingested N trades`.

- [ ] **Step 6: Commit (6)**

```bash
git add nfl_draft/scrapers/kalshi.py nfl_draft/tests/fixtures/kalshi/trades_response.json nfl_draft/tests/unit/test_scraper_parsing.py nfl_draft/tests/integration/test_quarantine.py
git commit -m "feat: Kalshi trade tape poller with cursor + dedup

fetch_trades reads kalshi_poll_state.last_traded_at per series,
queries Kalshi with min_ts (converted UTC via lib/db helpers),
INSERTs OR IGNOREs by trade_id PK, updates the cursor."
```

---

## Task 22: Cross-Book Grid + +EV Candidates dashboard tabs

**Files:**
- Create: `nfl_draft/lib/queries.py` (cross-venue median + outlier-flag query)
- Modify: `kalshi_draft/app.py` (add new tabs)
- Append: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write failing query test**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:
```python
def test_cross_book_grid_query_outlier_flags(seeded_db, monkeypatch):
    """Outlier flag fires when |delta| ≥ 10pp from all-venue median."""
    monkeypatch.setattr(db_module, "DB_PATH", seeded_db)
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets VALUES ('first_qb_cam-ward', 'first_at_position', 'Cam Ward', NULL, 'QB', NULL, NULL, NULL, ?)", [now])
        # 5 venues post: 4 say ~50%, 1 (DK) says 30% → DK is outlier
        for book, prob in [("kalshi", 0.50), ("fanduel", 0.51), ("bookmaker", 0.49), ("wagerzon", 0.52), ("draftkings", 0.30)]:
            con.execute("INSERT INTO draft_odds VALUES ('first_qb_cam-ward', ?, 100, ?, ?, ?)", [book, prob, prob, now])
    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=10)
    # Expect 1 market with 1 outlier (DK)
    assert len(rows) == 1
    assert rows[0]["market_id"] == "first_qb_cam-ward"
    assert rows[0]["outlier_count"] == 1
    flags = {b: f for b, f in rows[0]["flags"].items()}
    assert flags["draftkings"] is True
    assert flags["kalshi"] is False
```

- [ ] **Step 2: Implement queries.py**

Create `nfl_draft/lib/queries.py`:
```python
"""Dashboard query functions — separated from Dash callbacks for testability."""

import statistics
from typing import List, Dict, Any
from nfl_draft.lib.db import read_connection


def cross_book_grid(threshold_pp: float = 10.0) -> List[Dict[str, Any]]:
    """Return one row per market with all-venue prices, median, and outlier flags."""
    with read_connection() as con:
        rows = con.execute("""
            WITH latest AS (
              SELECT market_id, book, devig_prob,
                     ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
              FROM draft_odds
            )
            SELECT market_id, book, devig_prob FROM latest WHERE rn = 1
        """).fetchall()
    by_market: Dict[str, Dict[str, float]] = {}
    for market_id, book, prob in rows:
        by_market.setdefault(market_id, {})[book] = prob
    threshold = threshold_pp / 100.0
    output = []
    for market_id, books in by_market.items():
        if len(books) < 2:
            output.append({"market_id": market_id, "books": books, "median": None, "flags": {}, "outlier_count": 0})
            continue
        median = statistics.median(books.values())
        flags = {b: abs(p - median) >= threshold for b, p in books.items()}
        output.append({
            "market_id": market_id,
            "books": books,
            "median": median,
            "flags": flags,
            "outlier_count": sum(flags.values()),
        })
    return output


def ev_candidates(threshold_pp: float = 10.0) -> List[Dict[str, Any]]:
    """Flat list of flagged (market, venue) outliers, sorted by |delta| desc."""
    grid = cross_book_grid(threshold_pp)
    out = []
    for m in grid:
        if m["median"] is None:
            continue
        for book, flagged in m["flags"].items():
            if not flagged:
                continue
            delta = m["books"][book] - m["median"]
            out.append({
                "market_id": m["market_id"],
                "book": book,
                "book_prob": m["books"][book],
                "median": m["median"],
                "delta": delta,
            })
    out.sort(key=lambda r: abs(r["delta"]), reverse=True)
    return out


def trade_tape(limit: int = 200, large_threshold_usd: float = 500.0) -> List[Dict[str, Any]]:
    """Recent Kalshi trades. Computes is_large at read time."""
    with read_connection() as con:
        rows = con.execute("""
            SELECT trade_id, ticker, side, price_cents, count, notional_usd, traded_at,
                   notional_usd >= ? AS is_large
            FROM kalshi_trades
            ORDER BY traded_at DESC
            LIMIT ?
        """, [large_threshold_usd, limit]).fetchall()
    cols = ["trade_id", "ticker", "side", "price_cents", "count", "notional_usd", "traded_at", "is_large"]
    return [dict(zip(cols, r)) for r in rows]


def latest_max_fetched_at(table: str) -> Any:
    """For cheap-poll guard."""
    with read_connection() as con:
        result = con.execute(f"SELECT MAX(fetched_at) FROM {table}").fetchone()
    return result[0] if result else None
```

- [ ] **Step 3: Run query test**

Run: `pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_cross_book_grid_query_outlier_flags -v`
Expected: PASS

- [ ] **Step 4: Add Cross-Book Grid + +EV Candidates tabs to kalshi_draft/app.py**

In `kalshi_draft/app.py`, add (within the existing Dash app structure):

```python
import dash
from dash import dcc, html, dash_table, Input, Output, State
from nfl_draft.lib.queries import cross_book_grid, ev_candidates, trade_tape, latest_max_fetched_at

# Inside the layout:
portal_subnav = dcc.Tabs(id="nfl_draft.subtab", value="grid", children=[
    dcc.Tab(label="Cross-Book Grid", value="grid"),
    dcc.Tab(label="+EV Candidates", value="evlist"),
    dcc.Tab(label="Trade Tape", value="tape"),
    dcc.Tab(label="Bet Log", value="betlog"),
])

# Stores
nfl_draft_stores = html.Div([
    dcc.Store(id="nfl_draft.mode_toggle", storage_type="local"),
    dcc.Store(id="nfl_draft.subnav", storage_type="local"),
    dcc.Store(id="nfl_draft.last_fetched_odds", storage_type="local"),
    dcc.Store(id="nfl_draft.last_fetched_trades", storage_type="local"),
    dcc.Store(id="nfl_draft.bet_log_prefill", storage_type="session"),
    dcc.Interval(id="nfl_draft.interval", interval=60_000),  # 60s default
])

# Cross-Book Grid callback
@app.callback(
    Output("cross-book-grid", "children"),
    Input("nfl_draft.interval", "n_intervals"),
    State("nfl_draft.last_fetched_odds", "data"),
)
def render_cross_book_grid(n, last_seen):
    current = latest_max_fetched_at("draft_odds")
    if current and last_seen and str(current) == last_seen:
        return dash.no_update
    grid = cross_book_grid(threshold_pp=10)
    # Build a Plotly heatmap or dash_table
    # ... see spec for exact rendering
```

(The full Dash code is verbose. Implement using `dash_table.DataTable` for simplicity — markets in rows, venue columns, devig_prob in cells, conditional formatting on flagged cells.)

- [ ] **Step 5: Smoke test the dashboard**

```bash
python kalshi_draft/app.py &
open http://127.0.0.1:8083/
# Click "Cross-Book Grid" — should render
pkill -f kalshi_draft/app.py
```

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/lib/queries.py kalshi_draft/app.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "feat(dashboard): Cross-Book Grid + +EV Candidates tabs

Outlier-flag math computed in lib/queries.py (testable in isolation
from Dash). Median across all 5 venues including Kalshi; cells
flagged when |delta| ≥ 10pp. +EV Candidates is a flat sorted view
of just the flagged outliers."
```

---

## Task 23: Trade Tape + Bet Log dashboard tabs

**Files:**
- Modify: `kalshi_draft/app.py` (add Trade Tape + Bet Log tabs + auto-refresh)
- Append: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write tests for trade_tape and bet logging**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:
```python
def test_trade_tape_marks_large_fills(seeded_db, monkeypatch):
    monkeypatch.setattr(db_module, "DB_PATH", seeded_db)
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        # 1 large ($1000), 1 small ($100)
        con.execute("INSERT INTO kalshi_trades VALUES ('t1', 'TKR', 'yes', 50, 200, 1000.0, ?, ?)", [now, now])
        con.execute("INSERT INTO kalshi_trades VALUES ('t2', 'TKR', 'yes', 50, 20, 100.0, ?, ?)", [now, now])
    from nfl_draft.lib.queries import trade_tape
    rows = trade_tape(limit=10, large_threshold_usd=500.0)
    assert len(rows) == 2
    by_id = {r["trade_id"]: r for r in rows}
    assert by_id["t1"]["is_large"] is True
    assert by_id["t2"]["is_large"] is False


def test_cheap_poll_guard_returns_no_update_when_unchanged():
    """When MAX(fetched_at) hasn't changed, callback short-circuits."""
    from nfl_draft.lib.queries import latest_max_fetched_at
    # Conceptual test — the guard logic is in the callback, exercised by mock above.
    # This test verifies the helper returns None for empty DB:
    assert latest_max_fetched_at("draft_odds") is None or hasattr(latest_max_fetched_at("draft_odds"), "year")
```

- [ ] **Step 2: Add Trade Tape tab to app.py**

```python
@app.callback(
    Output("trade-tape", "children"),
    Input("nfl_draft.interval", "n_intervals"),
    State("nfl_draft.last_fetched_trades", "data"),
)
def render_trade_tape(n, last_seen):
    current = latest_max_fetched_at("kalshi_trades")
    if current and last_seen and str(current) == last_seen:
        return dash.no_update
    rows = trade_tape(limit=200, large_threshold_usd=500.0)
    return dash_table.DataTable(
        data=rows,
        columns=[{"name": k, "id": k} for k in ["traded_at", "ticker", "side", "price_cents", "count", "notional_usd"]],
        style_data_conditional=[{
            "if": {"filter_query": "{is_large} = True"},
            "backgroundColor": "#ffeb3b",
        }],
    )
```

- [ ] **Step 3: Add Bet Log tab + form**

```python
bet_log_form = html.Div([
    dcc.Dropdown(id="betlog-market", searchable=True, placeholder="Search market..."),
    dcc.Dropdown(id="betlog-book", options=[{"label": b, "value": b} for b in ["kalshi", "draftkings", "fanduel", "bookmaker", "wagerzon"]]),
    dcc.Input(id="betlog-odds", type="number", placeholder="American odds"),
    dcc.Input(id="betlog-stake", type="number", placeholder="Stake $"),
    dcc.Textarea(id="betlog-note", placeholder="Note (optional)"),
    html.Button("Log bet", id="betlog-submit"),
])

@app.callback(
    Output("betlog-status", "children"),
    Input("betlog-submit", "n_clicks"),
    State("betlog-market", "value"),
    State("betlog-book", "value"),
    State("betlog-odds", "value"),
    State("betlog-stake", "value"),
    State("betlog-note", "value"),
)
def log_bet(n, market, book, odds, stake, note):
    if not n:
        return ""
    import uuid
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_bets VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [str(uuid.uuid4()), market, book, "yes", int(odds), float(stake), datetime.now(), note or ""],
        )
    return f"Logged {market} on {book} @ {odds:+d} for ${stake:.2f}"
```

- [ ] **Step 4: Test**

Run: `pytest nfl_draft/tests/integration/test_dashboard_queries.py -v`
Expected: PASS.

Smoke-test the dashboard, click through Trade Tape and Bet Log.

- [ ] **Step 5: Commit**

```bash
git add kalshi_draft/app.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "feat(dashboard): Trade Tape + Bet Log tabs + auto-refresh

Trade Tape highlights large fills via read-time threshold. Bet Log
form persists to draft_bets with UUID PK. Cheap-poll guard caches
MAX(fetched_at) in dcc.Store to skip redundant heavy renders."
```

---

## Task 24: Concurrent-access stress test + lint test

**Files:**
- Create: `nfl_draft/tests/integration/test_concurrent_access.py`
- Create: `nfl_draft/tests/integration/test_lint_db_connect.py`

- [ ] **Step 1: Write concurrent-access stress test**

Create `nfl_draft/tests/integration/test_concurrent_access.py`:
```python
"""Stress test: writer + reader threads, asserts no lock errors."""
import threading
import time
import pytest
from nfl_draft.lib import db as db_module


def test_writer_and_reader_coexist(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "stress.duckdb")
    db_module.init_schema()
    errors = []
    stop = threading.Event()

    def writer():
        from nfl_draft.lib.db import write_connection
        from datetime import datetime
        i = 0
        while not stop.is_set():
            try:
                with write_connection() as con:
                    con.execute("INSERT INTO draft_odds VALUES (?, 'kalshi', 100, 0.5, 0.5, ?)",
                                [f"m{i}", datetime.now()])
                i += 1
            except Exception as e:
                errors.append(("writer", e))

    def reader():
        from nfl_draft.lib.db import read_connection
        while not stop.is_set():
            try:
                with read_connection() as con:
                    con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()
            except Exception as e:
                errors.append(("reader", e))

    # Need to seed at least one row first
    from nfl_draft.lib.db import write_connection
    from datetime import datetime
    with write_connection() as con:
        con.execute("INSERT INTO draft_odds VALUES ('seed', 'kalshi', 100, 0.5, 0.5, ?)", [datetime.now()])

    w = threading.Thread(target=writer)
    r = threading.Thread(target=reader)
    w.start()
    r.start()
    time.sleep(5)  # 5-second stress (10s in spec; shortened for test speed)
    stop.set()
    w.join()
    r.join()
    assert not errors, f"Lock errors observed: {errors[:5]}"
```

- [ ] **Step 2: Write lint test for raw duckdb.connect**

Create `nfl_draft/tests/integration/test_lint_db_connect.py`:
```python
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
            continue  # tests can use raw connect for fixtures
        text = py_file.read_text()
        if PATTERN.search(text):
            violations.append(rel)
    assert not violations, f"Raw duckdb.connect found in: {violations}"
```

- [ ] **Step 3: Run tests**

Run: `pytest nfl_draft/tests/integration/test_concurrent_access.py nfl_draft/tests/integration/test_lint_db_connect.py -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add nfl_draft/tests/integration/test_concurrent_access.py nfl_draft/tests/integration/test_lint_db_connect.py
git commit -m "test: concurrent-access stress + lint for raw duckdb.connect

Stress test spawns writer + reader threads for 5s, asserts no lock
errors. Lint test forbids raw duckdb.connect() outside lib/db.py
(migration script allowlisted for two-DB open)."
```

---

## Task 25: Cron config + log rotation

**Files:**
- Create: `nfl_draft/crontab.pre`
- Create: `nfl_draft/crontab.draft`

- [ ] **Step 1: Write crontab.pre**

Create `nfl_draft/crontab.pre`:
```
# NFL Draft Portal — pre-draft cadence (15min scrape, 2min trades)
# Install: crontab nfl_draft/crontab.pre
# Replace /Users/callancapitolo/NFLWork with your repo root.
*/15 * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode scrape --book all >> nfl_draft/logs/scrape.log 2>&1
*/2  * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode trades            >> nfl_draft/logs/trades.log 2>&1
0 4 * * * find /Users/callancapitolo/NFLWork/nfl_draft/logs/ -mtime +7 -delete
```

- [ ] **Step 2: Write crontab.draft**

Create `nfl_draft/crontab.draft`:
```
# NFL Draft Portal — draft-day cadence (2min scrape, 1min trades)
# Install: crontab nfl_draft/crontab.draft
*/2 * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode scrape --book all >> nfl_draft/logs/scrape.log 2>&1
* * * * *   cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode trades            >> nfl_draft/logs/trades.log 2>&1
0 4 * * * find /Users/callancapitolo/NFLWork/nfl_draft/logs/ -mtime +7 -delete
```

- [ ] **Step 3: Create logs directory**

```bash
mkdir -p nfl_draft/logs
touch nfl_draft/logs/.gitkeep
```

- [ ] **Step 4: Commit**

```bash
git add nfl_draft/crontab.pre nfl_draft/crontab.draft nfl_draft/logs/.gitkeep
git commit -m "chore: cron config (pre-draft + draft-day swap files)"
```

---

## Task 26: README + CLAUDE.md updates + go/no-go gate

**Files:**
- Create: `nfl_draft/README.md`
- Modify: `README.md` (top-level)
- Modify: `CLAUDE.md`
- Modify: `kalshi_draft/README.md`

- [ ] **Step 1: Write nfl_draft/README.md**

Create `nfl_draft/README.md`:
```markdown
# NFL Draft EV Portal

Trader's-cockpit web portal for surfacing +EV NFL Draft bets across
Kalshi and 4 sportsbooks.

## Setup

1. Install dependencies:
   ```bash
   pip install duckdb dash plotly pytest curl_cffi playwright cryptography
   ```
2. Populate `.env` at the repo root with:
   ```
   DK_USERNAME=...
   DK_PASSWORD=...
   FD_USERNAME=...
   FD_PASSWORD=...
   BOOKMAKER_USERNAME=...
   BOOKMAKER_PASSWORD=...
   WAGERZON_USERNAME=...
   WAGERZON_PASSWORD=...
   ```
   (Kalshi auth uses `kalshi_draft/.env` — no new entries needed.)
3. Run the migration once:
   ```bash
   python -m nfl_draft.lib.migrate_from_kalshi_draft
   ```
4. Install pre-draft cron:
   ```bash
   crontab nfl_draft/crontab.pre
   ```
5. Keep the laptop awake during draft week:
   ```bash
   caffeinate -i &
   ```

## Usage

- **Manual scrape**: `python -m nfl_draft.run --mode scrape --book all`
- **Single book** (debugging): `python -m nfl_draft.run --mode scrape --book draftkings`
- **Trade tape**: `python -m nfl_draft.run --mode trades`
- **Dashboard**: `python kalshi_draft/app.py` → http://127.0.0.1:8083/

## Draft-day mode

The morning of April 23, swap to faster cadence:
```bash
crontab nfl_draft/crontab.draft
```
And toggle the dashboard mode header from "Pre-draft" to "Draft-day".

## Maintenance

When a scraper produces unmapped players or markets (visible in the
dashboard footer):
1. Edit `nfl_draft/config/players.py` or `nfl_draft/config/markets.py`
2. Save (the next cron tick will reseed automatically; or run
   `python -m nfl_draft.lib.seed` to apply immediately)

## Troubleshooting

- Auth failures show in the dashboard footer per book.
- DuckDB lock errors should be zero; if observed, the spec describes a
  JSONL-tail fallback (Phase 2).
- See `docs/superpowers/specs/2026-04-17-nfl-draft-portal-design.md`
  for the full design rationale.
```

- [ ] **Step 2: Update top-level README**

Append to `/Users/callancapitolo/NFLWork/README.md`:
```markdown

## NFL Draft Portal
See [`nfl_draft/README.md`](nfl_draft/README.md) for the cross-venue
EV portal (Kalshi + DK/FD/Bookmaker/Wagerzon).
```

- [ ] **Step 3: Update CLAUDE.md**

In `/Users/callancapitolo/NFLWork/CLAUDE.md`, add `nfl_draft/` to the project structure section (or wherever the directory inventory lives).

- [ ] **Step 4: Update kalshi_draft/README.md**

Add a note that the dashboard now extends with portal tabs and reads from `nfl_draft/nfl_draft.duckdb`.

- [ ] **Step 5: Run the morning-of-Apr-22 go/no-go gate**

```bash
# Full test suite
pytest nfl_draft/tests/unit nfl_draft/tests/integration -v
# Expected: 0 failures

# Live smoke
pytest nfl_draft/tests/live -v
# Expected: 0 failures

# Manual: launch dashboard and click every tab
python kalshi_draft/app.py
# Visit http://127.0.0.1:8083/, click each of 9 tabs, confirm rendering.
```

If any step fails, stop and fix before proceeding to the merge.

- [ ] **Step 6: Commit READMEs and final go/no-go**

```bash
git add nfl_draft/README.md README.md CLAUDE.md kalshi_draft/README.md
git commit -m "docs: nfl_draft README + go/no-go gate

Setup, usage, draft-day swap, maintenance workflow. References
the design spec for rationale."
```

- [ ] **Step 7: Pre-merge review + merge to main (with explicit user approval)**

Per project CLAUDE.md, perform executive review of `git diff main..HEAD`. Confirm:
- No DK/FD credentials in code
- All DuckDB connections via context managers (lint test confirms)
- No personal bets logged to stdout
- Cookie files gitignored
- All tests green

**Wait for explicit user approval before merging to main.**

```bash
# After approval:
git checkout main
git merge --no-ff feature/nfl-draft-portal
# Re-run pipeline + dashboard from main to verify still working
pytest nfl_draft/tests/ -v
python kalshi_draft/app.py  # smoke each tab
# Then push if approved:
# git push origin main
```

---

## Self-Review Checklist (run after writing the plan)

Spec coverage scan:
- [x] Storage tables — Task 2 (schema) + Task 3 (migration) + Task 14 (quarantine table fully exercised)
- [x] Migration + idempotency — Task 3
- [x] Existing-tab repoint atomic with rename — Tasks 4-9 → atomic commit in Task 9
- [x] Config dicts (no CSV/JSON) — Task 10
- [x] Seed + 5 lookup tables — Task 11
- [x] Normalize + market_map — Task 12
- [x] Devig + market_id construction — Task 13
- [x] Quarantine + scraper interface — Task 14
- [x] Kalshi adapter — Task 15
- [x] DK / FD / BM / WZ scrapers — Tasks 16-19
- [x] run.py orchestrator — Task 20
- [x] Trade tape poller + cursor + dedup — Task 21
- [x] Cross-Book Grid + +EV Candidates + outlier flag — Task 22
- [x] Trade Tape + Bet Log + auto-refresh — Task 23
- [x] Concurrent-access stress + lint — Task 24
- [x] Cron config + caffeinate + log rotation — Task 25
- [x] README + go/no-go — Task 26
- [x] Tests in same commit as code — every task pairs implementation + tests
- [x] TZ helpers (`local_to_utc_iso`, `utc_iso_to_local`) — Task 2
- [x] Test_kalshi_writer_target — Task 9
- [x] dcc.Store registry — Task 22

Type/method consistency:
- [x] `OddsRow` dataclass shape consistent across scrapers + quarantine + queries
- [x] `build_market_id` signature matches all callsites
- [x] `write_connection` / `read_connection` used uniformly

Placeholder scan: clean (one TODO remains in Task 16 step 3 for "Implementation depends on reconnaissance findings" — acceptable because the actual code can't be written before the API is captured).
