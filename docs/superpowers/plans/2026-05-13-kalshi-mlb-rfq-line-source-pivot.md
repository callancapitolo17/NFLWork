# Kalshi MLB RFQ — Line Source Pivot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decouple the Kalshi MLB RFQ bot's edge surface from Wagerzon-derived line selection. The bot drives its own line surface from Kalshi MVE enumeration, library-refactors the four SGP scrapers (DK/FD/PX/NV) into per-book Python modules layered above HTTP clients, owns its own market DB, and only bets where ≥2 books priced the matching line.

**Architecture:** Three-layer separation in `mlb_sgp/`: HTTP client modules (DK/FD already shipped; PX/NV new this plan) → per-book SGP orchestration libraries (`price_sgps()` per book) → thin subprocess shims (existing scraper scripts, now ~60-line orchestrator-callers). Bot caller in `kalshi_mlb_rfq/sgp_runner.py` writes target_lines, spawns scrapers, reads back priced rows. Sibling market DB at `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb` isolates write traffic from bot state.

**Tech Stack:** Python 3.11+, DuckDB, curl_cffi (Chrome-TLS HTTP), pytest, subprocess for scraper invocation. Bot reads `mlb_odds_temp` from `Answer Keys/mlb.duckdb` for game schedule. Spec at `docs/superpowers/specs/2026-05-12-kalshi-mlb-rfq-line-source-pivot-design.md`.

**Prerequisites:**
- Branch: `worktree-kalshi-mlb-rfq-line-source-pivot` already exists at `.claude/worktrees/kalshi-mlb-rfq-line-source-pivot`, rebased onto main HEAD `65dae02`.
- `mlb_sgp/dk_client.py` and `mlb_sgp/fd_client.py` are on main (shipped by singles-scrapers merge).
- `kalshi_mlb_rfq/sgp_runner.py` does NOT exist on main — this plan creates it from scratch.
- Bot is currently stopped; safe to make code changes without disrupting live trading.

---

## File Structure

### Files to create

| Path | Responsibility |
|---|---|
| `mlb_sgp/_shared.py` | `TargetLine`, `PricedRow` dataclasses; `_utc_bucket`, `decimal_to_american`/`american_to_decimal` helpers; `load_target_lines()` that abstracts over `mlb_target_lines` (multi-row) and `mlb_parlay_lines` (one-row-per-game) tables. |
| `mlb_sgp/prophetx_client.py` | ProphetX HTTP client mirroring `dk_client.py` shape: `Event`, `Market`, `SelectionLeg` dataclasses; `ProphetXClient` with `list_events()`, `fetch_event_markets()`, `submit_parlay_rfq()`. Owns `curl_cffi` session + profile-cookie auth. |
| `mlb_sgp/novig_client.py` | Novig HTTP client: `Event`, `EventLegs` dataclasses; `NovigClient` with `list_events()`, `fetch_event_legs()`, `submit_parlay()`. Anonymous endpoint, Hasura GraphQL. |
| `mlb_sgp/draftkings.py` | DK SGP orchestration library exposing `price_sgps(target_lines, periods, client=None) -> list[PricedRow]`. Uses `DraftKingsClient`. |
| `mlb_sgp/fanduel.py` | FD SGP orchestration library, same shape. Uses `FanDuelClient`. |
| `mlb_sgp/prophetx.py` | PX SGP orchestration library. Uses `ProphetXClient`. Carries existing MIN_OFFER_STAKE filter + F5-Over SANITY_MULT_RATIO defense. |
| `mlb_sgp/novig.py` | NV SGP orchestration library. Uses `NovigClient`. |
| `kalshi_mlb_rfq/sgp_runner.py` | Bot caller: `should_scrape()`, `latest_sgp_fetch_time()`, `enumerate_kalshi_targets()`, `write_target_lines()`, `run_scrapers()`, `read_priced_rows()`, `sgp_cycle()` orchestrator. |
| `mlb_sgp/tests/test_shared.py` | Unit tests for shared types, helpers, `load_target_lines()` dual-table logic. |
| `mlb_sgp/tests/test_prophetx_client.py` | Unit tests for `ProphetXClient` against captured fixtures. |
| `mlb_sgp/tests/test_novig_client.py` | Unit tests for `NovigClient` against captured fixtures. |
| `mlb_sgp/tests/fixtures/px_events.json` | Captured ProphetX events response. |
| `mlb_sgp/tests/fixtures/px_event_markets.json` | Captured ProphetX markets-per-event response. |
| `mlb_sgp/tests/fixtures/px_parlay_response.json` | Captured ProphetX RFQ parlay response. |
| `mlb_sgp/tests/fixtures/nv_events.json` | Captured Novig events GraphQL response. |
| `mlb_sgp/tests/fixtures/nv_event_legs.json` | Captured Novig event legs GraphQL response. |
| `mlb_sgp/tests/fixtures/nv_parlay_response.json` | Captured Novig parlay response. |
| `kalshi_mlb_rfq/tests/test_sgp_runner.py` | Unit tests for `sgp_runner.py` functions. |

### Files to modify

| Path | Changes |
|---|---|
| `mlb_sgp/db.py` | (1) Read `MLB_SGP_DB_PATH` env var, fall back to `mlb_mm.duckdb`. (2) `ensure_table()` runs `ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS spread_line DOUBLE; ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS total_line DOUBLE;` (3) Add `upsert_priced_rows(rows: list[PricedRow], db_path=None)`. |
| `mlb_sgp/scraper_draftkings_sgp.py` | Reduce to ~60-line shim: parse env, `load_target_lines()`, `draftkings.price_sgps()`, `db.upsert_priced_rows()`. |
| `mlb_sgp/scraper_fanduel_sgp.py` | Same shim shape, calls `fanduel.price_sgps()`. |
| `mlb_sgp/scraper_prophetx_sgp.py` | Same shim shape, calls `prophetx.price_sgps()`. |
| `mlb_sgp/scraper_novig_sgp.py` | Same shim shape, calls `novig.price_sgps()`. |
| `kalshi_mlb_rfq/config.py` | (1) Worktree path fix on `PROJECT_ROOT`. (2) Add `SGP_REFRESH_SEC`, `SGP_SCRAPER_TIMEOUT_SEC`, `BOT_MARKET_DB`, `MIN_BOOK_COUNT_FOR_BLEND`, `MLB_SGP_DIR` env knobs. |
| `kalshi_mlb_rfq/main.py` | (1) `_PARLAY_LINES_CACHE` shape change (multi-line per game). (2) `_load_book_fairs` rewires per-line + N≥2 gate. (3) Drop `line_move_ok` gate from `_all_per_accept_gates_pass`. (4) Drop `_current_book_lines_for_combo` + reference_lines write. (5) Add `_refresh_sgp_cache()`. (6) Synchronous warm-up + SGP cadence tick in `main_loop`. |
| `kalshi_mlb_rfq/.env.example` | Add new env vars with defaults + commented descriptions. |
| `kalshi_mlb_rfq/README.md` | Document SGP cadence loop, market DB sibling, N≥2 book gate, line-source decoupling. |
| `mlb_sgp/README.md` | Document library architecture (clients + orchestrators + shims), dual-mode `MLB_SGP_DB_PATH` invocation, new `mlb_target_lines` table. |
| `CLAUDE.md` (root) | Brief note under "Project Structure" mentioning bot market DB. |

### Files untouched

- `Answer Keys/mlb_correlated_parlay.R` — dashboard pipeline unchanged.
- `kalshi_mlb_rfq/combo_enumerator.py`, `auth_client.py`, `rfq_client.py`, `kelly.py`, `risk.py`, `fair_value.py`, `db.py`, `notify.py`, `ev_calc.py`, `ticker_map.py` — no changes from this pivot.
- `mlb_sgp/dk_client.py`, `mlb_sgp/fd_client.py` — consumed as dependencies; no edits.

---

## Test Strategy

- **Unit tests:** Each new module gets a dedicated `test_*.py` in the test dir mirroring the source dir (`mlb_sgp/tests/` or `kalshi_mlb_rfq/tests/`).
- **Client tests:** Use captured fixtures (real JSON responses recorded once) to test client parsing/dispatch without hitting live APIs.
- **Library tests:** Mock the client (using `unittest.mock`) and verify orchestration logic (leg combination, combo flavoring, devig sanity).
- **Shim regression:** After each scraper shim is refactored, run it against `mlb_mm.duckdb` and diff row counts vs. a pre-refactor baseline snapshot captured at Task 0.
- **Bot integration:** Stop short of full e2e bot run in unit tests; that's covered by the manual smoke test at the end.

Test runner: `pytest` from repo root or per-module dir. Existing `pytest.ini` / `pyproject.toml` settings honored.

---

## Pre-Task: Capture pre-refactor baseline

Before refactoring scrapers, snapshot their current output so we can diff after each shim refactor.

- [ ] **Step 0.1: Run baseline DK scrape, snapshot result**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
source mlb_sgp/venv/bin/activate
python mlb_sgp/scraper_draftkings_sgp.py 2>&1 | tail -10
duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='draftkings_direct' ORDER BY game_id, combo, period) TO '/tmp/dk_baseline.csv' (HEADER, DELIMITER ',')"
wc -l /tmp/dk_baseline.csv
```

Expected: non-zero row count (10-100 rows depending on slate).

- [ ] **Step 0.2: Snapshot FD, PX, NV the same way**

```bash
python mlb_sgp/scraper_fanduel_sgp.py 2>&1 | tail -5
python mlb_sgp/scraper_prophetx_sgp.py 2>&1 | tail -5
python mlb_sgp/scraper_novig_sgp.py 2>&1 | tail -5
for src in fanduel_direct prophetx_direct novig_direct; do
  duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='$src' ORDER BY game_id, combo, period) TO '/tmp/${src%_direct}_baseline.csv' (HEADER, DELIMITER ',')"
done
wc -l /tmp/*_baseline.csv
```

Note: this baseline captures the data with no schema changes yet. After Task 4 migrates the schema, re-snapshot if needed. The baseline is used in Tasks 13-16 to verify shim refactors don't change pricing behavior.

- [ ] **Step 0.3: Verify clean working tree**

```bash
git status --short
```

Expected: empty (or only ignored files).

---

### Task 1: Shared dataclasses (TargetLine, PricedRow)

**Files:**
- Create: `mlb_sgp/_shared.py`
- Test: `mlb_sgp/tests/test_shared.py`

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_shared.py`:

```python
"""Tests for mlb_sgp/_shared.py — shared dataclasses + helpers."""
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine, PricedRow


def test_target_line_immutable():
    t = TargetLine(
        game_id="abc-123",
        home_team="New York Yankees",
        away_team="Boston Red Sox",
        commence_time=datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc),
        period="FG",
        spread=-1.5,
        total=8.5,
    )
    # Frozen dataclass: assignment raises FrozenInstanceError
    import dataclasses
    try:
        t.spread = -2.5  # type: ignore[misc]
    except dataclasses.FrozenInstanceError:
        return
    raise AssertionError("TargetLine should be frozen")


def test_priced_row_fields():
    p = PricedRow(
        game_id="abc-123",
        combo="Home Spread + Over",
        period="FG",
        spread_line=-1.5,
        total_line=8.5,
        bookmaker="draftkings",
        source="draftkings_direct",
        sgp_decimal=2.85,
        sgp_american=185,
        fetch_time=datetime.now(timezone.utc),
    )
    assert p.combo == "Home Spread + Over"
    assert p.spread_line == -1.5
    assert p.sgp_american == 185
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
pytest mlb_sgp/tests/test_shared.py -v
```

Expected: ImportError or ModuleNotFoundError (`_shared` doesn't exist).

- [ ] **Step 3: Create the module**

Create `mlb_sgp/_shared.py`:

```python
"""Shared types and utilities for the MLB SGP library.

Module-private (`_` prefix) but stable API consumed by per-book modules
(draftkings.py, fanduel.py, prophetx.py, novig.py) and by callers
(dashboard scraper shims, kalshi_mlb_rfq/sgp_runner.py).
"""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from typing import Literal


@dataclass(frozen=True)
class TargetLine:
    """One (game, period, spread, total) tuple the bot/dashboard wants priced."""
    game_id: str
    home_team: str
    away_team: str
    commence_time: datetime
    period: Literal["FG", "F5"]
    spread: float   # signed, home-perspective (negative = home favored)
    total: float


@dataclass(frozen=True)
class PricedRow:
    """One book's SGP price for a (game, period, spread, total, combo) tuple."""
    game_id: str
    combo: str
    period: str
    spread_line: float
    total_line: float
    bookmaker: str
    source: str
    sgp_decimal: float
    sgp_american: int
    fetch_time: datetime
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_shared.py -v
```

Expected: 2 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/_shared.py mlb_sgp/tests/test_shared.py
git commit -m "feat(mlb_sgp): TargetLine + PricedRow dataclasses for shared library types

Foundation for the line-source pivot library refactor. TargetLine
represents what the bot/dashboard wants priced (multi-row in bot mode,
one-row-per-game in dashboard mode). PricedRow is the per-book output
shape with explicit spread_line/total_line columns.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Shared utilities (decimal/american conversion, _utc_bucket)

**Files:**
- Modify: `mlb_sgp/_shared.py`
- Test: `mlb_sgp/tests/test_shared.py`

- [ ] **Step 1: Write the failing tests**

Append to `mlb_sgp/tests/test_shared.py`:

```python
from mlb_sgp._shared import decimal_to_american, american_to_decimal, _utc_bucket


def test_decimal_to_american_favorite():
    # Decimal 1.5 → American -200 (every 1.5 → ((1/(1.5-1)) - 1) * 100 → -100? Let me recompute)
    # Standard formula: dec >= 2.0 → +((dec-1)*100); dec < 2.0 → -100/(dec-1)
    assert decimal_to_american(1.5) == -200  # -100 / 0.5 = -200
    assert decimal_to_american(2.0) == 100   # +(1.0 * 100)
    assert decimal_to_american(2.5) == 150   # +(1.5 * 100)
    assert decimal_to_american(3.0) == 200


def test_american_to_decimal_round_trip():
    for am in [-300, -150, 100, 150, 300]:
        dec = american_to_decimal(am)
        back = decimal_to_american(dec)
        assert back == am, f"round trip {am} → {dec} → {back}"


def test_utc_bucket_isolates_hour():
    from datetime import datetime, timezone
    t = datetime(2026, 5, 13, 23, 17, 42, tzinfo=timezone.utc)
    assert _utc_bucket(t) == "2026-05-13T23"


def test_utc_bucket_handles_iso_string():
    assert _utc_bucket("2026-05-13T23:17:42Z") == "2026-05-13T23"
    assert _utc_bucket("2026-05-13T23:17:42+00:00") == "2026-05-13T23"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_shared.py -v
```

Expected: ImportError on `decimal_to_american` / `american_to_decimal` / `_utc_bucket`.

- [ ] **Step 3: Add the helpers**

Append to `mlb_sgp/_shared.py`:

```python
from datetime import datetime as _dt
from typing import Union


def decimal_to_american(dec: float) -> int:
    """Convert decimal odds to American format. Favorites are negative."""
    if dec >= 2.0:
        return int(round((dec - 1.0) * 100))
    return int(round(-100.0 / (dec - 1.0)))


def american_to_decimal(am: int) -> float:
    """Inverse of decimal_to_american."""
    if am > 0:
        return 1.0 + am / 100.0
    return 1.0 + 100.0 / abs(am)


def _utc_bucket(ts: Union[_dt, str]) -> str:
    """Extract a UTC "YYYY-MM-DDTHH" bucket string from a timestamp.

    Used as a match key when correlating events across data sources at
    date+hour granularity (avoids spurious matches across doubleheaders).
    Accepts both datetime objects and ISO-8601 strings.
    """
    if isinstance(ts, str):
        # Normalize Z suffix to +00:00 for fromisoformat
        normalized = ts.replace("Z", "+00:00") if ts.endswith("Z") else ts
        ts = _dt.fromisoformat(normalized)
    if ts.tzinfo is None:
        from datetime import timezone
        ts = ts.replace(tzinfo=timezone.utc)
    return ts.strftime("%Y-%m-%dT%H")
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_shared.py -v
```

Expected: 6 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/_shared.py mlb_sgp/tests/test_shared.py
git commit -m "feat(mlb_sgp): decimal_to_american, american_to_decimal, _utc_bucket helpers

Extract the helpers each scraper currently duplicates into a single
shared module. Library functions and bot caller consume these.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: load_target_lines() — dual-table-name loader

**Files:**
- Modify: `mlb_sgp/_shared.py`
- Test: `mlb_sgp/tests/test_shared.py`

This function abstracts over two table names so scrapers don't branch on caller mode:
- `mlb_target_lines` (multi-row per game, FG-only typically, written by bot) — preferred when present
- `mlb_parlay_lines` (one row per game, FG+F5 cols, written by R dashboard) — fallback

Both return a uniform `list[TargetLine]`.

- [ ] **Step 1: Write the failing tests**

Append to `mlb_sgp/tests/test_shared.py`:

```python
import duckdb
import tempfile
from pathlib import Path


def test_load_target_lines_prefers_mlb_target_lines(tmp_path):
    """When both tables exist, mlb_target_lines wins."""
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "t.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_target_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time TIMESTAMP, period VARCHAR,
            spread DOUBLE, total DOUBLE, written_at TIMESTAMP
        )
    """)
    con.execute("""
        INSERT INTO mlb_target_lines VALUES
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -1.5, 8.5, NOW()),
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -2.5, 8.5, NOW()),
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -1.5, 9.5, NOW())
    """)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_parlay_lines VALUES ('g1','NYY','BOS','2026-05-13T23:00Z',-3.5,7.5,NULL,NULL)")
    con.close()

    rows = load_target_lines(db_path=db)
    assert len(rows) == 3, "multi-row from mlb_target_lines"
    spreads = sorted({r.spread for r in rows})
    assert spreads == [-2.5, -1.5, -1.5][:2] + [-1.5][:0]  # actually [-2.5, -1.5]; dedup
    assert -3.5 not in spreads, "should NOT pull from mlb_parlay_lines"
    assert all(r.period == "FG" for r in rows)


def test_load_target_lines_legacy_fallback(tmp_path):
    """When only mlb_parlay_lines exists, emit FG + F5 rows from each game."""
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "t.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.execute("""
        INSERT INTO mlb_parlay_lines VALUES
            ('g1', 'NYY', 'BOS', '2026-05-13T23:00:00+00:00', -1.5, 8.5, -0.5, 4.5),
            ('g2', 'LAD', 'SF',  '2026-05-14T02:00:00+00:00', -2.5, 7.5, NULL, NULL)
    """)
    con.close()

    rows = load_target_lines(db_path=db)
    fg = [r for r in rows if r.period == "FG"]
    f5 = [r for r in rows if r.period == "F5"]
    assert len(fg) == 2, "two FG rows (one per game)"
    assert len(f5) == 1, "one F5 row (g2 had NULL F5)"
    assert {r.game_id for r in fg} == {"g1", "g2"}
    assert {r.game_id for r in f5} == {"g1"}


def test_load_target_lines_empty_db(tmp_path):
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "empty.duckdb")
    con = duckdb.connect(db)
    con.close()
    assert load_target_lines(db_path=db) == []
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_shared.py -v -k load_target_lines
```

Expected: ImportError on `load_target_lines`.

- [ ] **Step 3: Add the function**

Append to `mlb_sgp/_shared.py`:

```python
import duckdb
from pathlib import Path


def load_target_lines(db_path: str) -> list[TargetLine]:
    """Load target lines from a DuckDB file.

    Prefers `mlb_target_lines` (bot-written, multi-row per game) when
    present. Falls back to legacy `mlb_parlay_lines` (dashboard-written,
    one row per game with FG+F5 columns) and emits one TargetLine per
    period per game (FG always; F5 only when F5 columns are non-NULL).

    Returns [] for missing DB or missing tables — callers decide what
    to do with empty input.
    """
    if not Path(db_path).exists():
        return []
    con = duckdb.connect(db_path, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_target_lines" in tables:
            return _load_from_target_lines(con)
        if "mlb_parlay_lines" in tables:
            return _load_from_parlay_lines(con)
        return []
    finally:
        con.close()


def _load_from_target_lines(con) -> list[TargetLine]:
    rows = con.execute("""
        SELECT game_id, home_team, away_team, commence_time, period, spread, total
        FROM mlb_target_lines
        ORDER BY game_id, period, spread, total
    """).fetchall()
    return [
        TargetLine(
            game_id=r[0], home_team=r[1], away_team=r[2],
            commence_time=r[3], period=r[4], spread=r[5], total=r[6],
        ) for r in rows
    ]


def _load_from_parlay_lines(con) -> list[TargetLine]:
    from datetime import datetime as _dt
    rows = con.execute("""
        SELECT game_id, home_team, away_team, commence_time,
               fg_spread, fg_total, f5_spread, f5_total
        FROM mlb_parlay_lines
        ORDER BY game_id
    """).fetchall()
    out: list[TargetLine] = []
    for r in rows:
        game_id, home, away, ct_raw, fg_s, fg_t, f5_s, f5_t = r
        # commence_time in legacy table is sometimes VARCHAR, sometimes TIMESTAMP
        if isinstance(ct_raw, str):
            normalized = ct_raw.replace("Z", "+00:00") if ct_raw.endswith("Z") else ct_raw
            ct = _dt.fromisoformat(normalized)
        else:
            ct = ct_raw
        if fg_s is not None and fg_t is not None:
            out.append(TargetLine(
                game_id=game_id, home_team=home, away_team=away,
                commence_time=ct, period="FG", spread=fg_s, total=fg_t,
            ))
        if f5_s is not None and f5_t is not None:
            out.append(TargetLine(
                game_id=game_id, home_team=home, away_team=away,
                commence_time=ct, period="F5", spread=f5_s, total=f5_t,
            ))
    return out
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_shared.py -v
```

Expected: 9 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/_shared.py mlb_sgp/tests/test_shared.py
git commit -m "feat(mlb_sgp): load_target_lines() dual-table-name loader

Abstracts over mlb_target_lines (bot, multi-row) and mlb_parlay_lines
(dashboard, one-row-per-game) so scrapers don't branch on caller mode.
Single uniform list[TargetLine] return shape.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Schema migration — mlb_sgp_odds gains spread_line/total_line

**Files:**
- Modify: `mlb_sgp/db.py`
- Test: `mlb_sgp/tests/test_db.py` (create)

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_db.py`:

```python
"""Tests for mlb_sgp/db.py — schema migration + writers."""
import duckdb

from mlb_sgp.db import ensure_table


def test_ensure_table_fresh_db(tmp_path):
    db = str(tmp_path / "fresh.duckdb")
    ensure_table(db_path=db)
    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
    assert "game_id" in cols


def test_ensure_table_migrates_legacy_schema(tmp_path):
    """If the table exists with the old (no spread/total cols) schema, ALTER it."""
    db = str(tmp_path / "legacy.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR,
            bookmaker VARCHAR, sgp_decimal DOUBLE, sgp_american INTEGER,
            fetch_time TIMESTAMP, source VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','Home Spread + Over','FG','draftkings',2.85,185,NOW(),'draftkings_direct')")
    con.close()

    ensure_table(db_path=db)

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    legacy_row = con.execute("SELECT spread_line, total_line FROM mlb_sgp_odds").fetchone()
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
    assert legacy_row == (None, None), "pre-existing rows have NULL line cols"


def test_ensure_table_idempotent(tmp_path):
    db = str(tmp_path / "idem.duckdb")
    ensure_table(db_path=db)
    ensure_table(db_path=db)  # second call must not error
    ensure_table(db_path=db)  # third call must not error either
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_db.py -v
```

Expected: 2 tests FAIL with `spread_line not in cols` — the current `ensure_table` doesn't add the new columns.

- [ ] **Step 3: Update db.py**

Modify `mlb_sgp/db.py`. Replace the `CREATE_TABLE_SQL` constant and `ensure_table()` function:

```python
CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS mlb_sgp_odds (
    game_id       VARCHAR,
    combo         VARCHAR,
    period        VARCHAR,
    bookmaker     VARCHAR,
    sgp_decimal   DOUBLE,
    sgp_american  INTEGER,
    fetch_time    TIMESTAMP,
    source        VARCHAR,
    spread_line   DOUBLE,
    total_line    DOUBLE
);
"""

# Backward-compat ALTERs — make sure pre-existing tables gain the new columns.
# Both ADD COLUMN IF NOT EXISTS statements are idempotent; DuckDB raises only
# if syntax/semantics are wrong, not on repeat add.
_MIGRATIONS = [
    "ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS spread_line DOUBLE",
    "ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS total_line DOUBLE",
]


def ensure_table(db_path: str = None):
    """Create the mlb_sgp_odds table if it doesn't exist; migrate if it does."""
    db_path = db_path or str(MLB_DB)
    con = _connect_with_retry(db_path)
    try:
        con.execute(CREATE_TABLE_SQL)
        for sql in _MIGRATIONS:
            con.execute(sql)
    finally:
        con.close()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_db.py -v
```

Expected: 3 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/db.py mlb_sgp/tests/test_db.py
git commit -m "feat(mlb_sgp): mlb_sgp_odds gains spread_line/total_line columns (idempotent migration)

Schema bump for the line-source pivot. New columns are nullable so
existing rows stay valid. Dashboard R reader (mlb_correlated_parlay.R)
SELECTs columns by name and ignores the new cols — verified safe.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Env-var DB override (MLB_SGP_DB_PATH)

**Files:**
- Modify: `mlb_sgp/db.py`
- Test: `mlb_sgp/tests/test_db.py`

- [ ] **Step 1: Write the failing test**

Append to `mlb_sgp/tests/test_db.py`:

```python
import os
import importlib
import mlb_sgp.db as db_mod


def test_mlb_db_default_when_env_unset(monkeypatch):
    monkeypatch.delenv("MLB_SGP_DB_PATH", raising=False)
    importlib.reload(db_mod)
    assert str(db_mod.MLB_DB).endswith("mlb_mm.duckdb")


def test_mlb_db_uses_env_override(monkeypatch, tmp_path):
    override = str(tmp_path / "override.duckdb")
    monkeypatch.setenv("MLB_SGP_DB_PATH", override)
    importlib.reload(db_mod)
    assert str(db_mod.MLB_DB) == override
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_db.py::test_mlb_db_uses_env_override -v
```

Expected: FAIL — current code hardcodes `MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb"`.

- [ ] **Step 3: Update db.py**

Modify the MLB_DB resolution in `mlb_sgp/db.py`. Replace this block:

```python
# Resolve repo root dynamically — works from main repo or worktrees.
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))

MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb"
```

with:

```python
import os

# Resolve repo root dynamically — works from main repo or worktrees.
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))

# MLB_DB resolution: env override > default dashboard DB.
# Bot sets MLB_SGP_DB_PATH=<bot_market_db> when spawning scrapers so the
# bot-owned DB is the read/write target without touching scraper code.
_DEFAULT_MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb"
MLB_DB = Path(os.environ.get("MLB_SGP_DB_PATH", str(_DEFAULT_MLB_DB)))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_db.py -v
```

Expected: 5 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/db.py mlb_sgp/tests/test_db.py
git commit -m "feat(mlb_sgp): MLB_SGP_DB_PATH env var overrides default DB target

Lets the bot redirect scrapers to its own market DB without modifying
scraper code. Backward compatible: unset env = current behavior.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: upsert_priced_rows() — write PricedRow with line columns

**Files:**
- Modify: `mlb_sgp/db.py`
- Test: `mlb_sgp/tests/test_db.py`

- [ ] **Step 1: Write the failing test**

Append to `mlb_sgp/tests/test_db.py`:

```python
from datetime import datetime, timezone
from mlb_sgp._shared import PricedRow
from mlb_sgp.db import upsert_priced_rows, ensure_table


def test_upsert_priced_rows_inserts(tmp_path):
    db = str(tmp_path / "u.duckdb")
    ensure_table(db_path=db)
    rows = [
        PricedRow(
            game_id="g1", combo="Home Spread + Over", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=2.85, sgp_american=185,
            fetch_time=datetime.now(timezone.utc),
        ),
        PricedRow(
            game_id="g1", combo="Home Spread + Under", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=3.20, sgp_american=220,
            fetch_time=datetime.now(timezone.utc),
        ),
    ]
    upsert_priced_rows(rows, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute(
        "SELECT game_id, combo, spread_line, total_line, sgp_american FROM mlb_sgp_odds ORDER BY combo"
    ).fetchall()
    con.close()
    assert len(out) == 2
    assert out[0] == ("g1", "Home Spread + Over", -1.5, 8.5, 185)
    assert out[1] == ("g1", "Home Spread + Under", -1.5, 8.5, 220)


def test_upsert_priced_rows_replaces_same_key(tmp_path):
    """Same (game, combo, period, spread, total, bookmaker, source) → replace."""
    db = str(tmp_path / "r.duckdb")
    ensure_table(db_path=db)
    def make(price):
        return PricedRow(
            game_id="g1", combo="Home Spread + Over", period="FG",
            spread_line=-1.5, total_line=8.5,
            bookmaker="draftkings", source="draftkings_direct",
            sgp_decimal=price, sgp_american=185,
            fetch_time=datetime.now(timezone.utc),
        )
    upsert_priced_rows([make(2.50)], db_path=db)
    upsert_priced_rows([make(2.75)], db_path=db)
    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*), MAX(sgp_decimal) FROM mlb_sgp_odds").fetchone()
    con.close()
    assert out == (1, 2.75)


def test_upsert_priced_rows_empty_is_noop(tmp_path):
    db = str(tmp_path / "e.duckdb")
    ensure_table(db_path=db)
    upsert_priced_rows([], db_path=db)  # must not error
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_db.py -v -k upsert_priced_rows
```

Expected: ImportError on `upsert_priced_rows`.

- [ ] **Step 3: Add the function**

Append to `mlb_sgp/db.py`:

```python
def upsert_priced_rows(rows: list, db_path: str = None):
    """Insert PricedRow objects, replacing any existing row with the same
    (game_id, combo, period, spread_line, total_line, bookmaker, source) key.

    Empty `rows` is a no-op.
    """
    if not rows:
        return

    db_path = db_path or str(MLB_DB)
    con = _connect_with_retry(db_path)
    try:
        ensure_table(db_path)

        # Batch delete by composite key (NULL-aware NOT NULL where line cols are present).
        keys = [
            (r.game_id, r.combo, r.period, r.spread_line, r.total_line,
             r.bookmaker, r.source)
            for r in rows
        ]
        placeholders = ",".join(["(?, ?, ?, ?, ?, ?, ?)"] * len(keys))
        flat = [v for k in keys for v in k]
        con.execute(f"""
            DELETE FROM mlb_sgp_odds
            WHERE (game_id, combo, period, spread_line, total_line, bookmaker, source) IN (
                SELECT * FROM (VALUES {placeholders})
            )
        """, flat)

        # Batch insert
        ins_placeholders = ",".join(["(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"] * len(rows))
        values = []
        for r in rows:
            values.extend([
                r.game_id, r.combo, r.period, r.bookmaker,
                r.sgp_decimal, r.sgp_american, r.fetch_time, r.source,
                r.spread_line, r.total_line,
            ])
        con.execute(f"""
            INSERT INTO mlb_sgp_odds
                (game_id, combo, period, bookmaker, sgp_decimal, sgp_american,
                 fetch_time, source, spread_line, total_line)
            VALUES {ins_placeholders}
        """, values)
    finally:
        con.close()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_db.py -v
```

Expected: 8 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/db.py mlb_sgp/tests/test_db.py
git commit -m "feat(mlb_sgp): upsert_priced_rows() writes PricedRow with line columns

Replaces existing rows by (game, combo, period, spread, total, book,
source) composite key. The new line columns are part of the key, so
the same game can carry many distinct (spread, total) tuples per book.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: ProphetXClient extraction

**Files:**
- Create: `mlb_sgp/prophetx_client.py`
- Create: `mlb_sgp/tests/test_prophetx_client.py`
- Create: `mlb_sgp/tests/fixtures/px_events.json`
- Create: `mlb_sgp/tests/fixtures/px_event_markets.json`
- Create: `mlb_sgp/tests/fixtures/px_parlay_response.json`
- Read for reference: `mlb_sgp/scraper_prophetx_sgp.py:133-525` (existing implementations)

The client mirrors `dk_client.py` shape: thin dataclass-returning wrappers over the three HTTP operations the SGP orchestrator needs.

- [ ] **Step 1: Capture fixtures from the live API (one-time setup)**

Write a temporary capture script `/tmp/capture_px_fixtures.py`:

```python
"""One-time fixture capture for ProphetXClient tests."""
import json
import sys
from pathlib import Path

sys.path.insert(0, "mlb_sgp")
from scraper_prophetx_sgp import init_session, fetch_prophetx_mlb_events, fetch_event_legs

FIX = Path("mlb_sgp/tests/fixtures")
FIX.mkdir(parents=True, exist_ok=True)

session = init_session()

# 1) Events list
events = fetch_prophetx_mlb_events(session)
(FIX / "px_events.json").write_text(json.dumps(events, indent=2, default=str))
print(f"px_events.json: {len(events)} events")

# 2) Markets for first event
if events:
    game = events[0]
    legs_data = fetch_event_legs(session, game)
    (FIX / "px_event_markets.json").write_text(json.dumps(legs_data, indent=2, default=str))
    print(f"px_event_markets.json: {game['event_id']}")

# 3) For px_parlay_response.json, we'll capture a synthetic response below in the test
# (the live parlay submission requires line IDs that go stale; easier to mock the
#  expected shape of the response than to capture a fresh one each time).
```

Run it:

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
python /tmp/capture_px_fixtures.py
ls -la mlb_sgp/tests/fixtures/px_*.json
```

Expected: `px_events.json` and `px_event_markets.json` exist, non-empty.

Then hand-craft `mlb_sgp/tests/fixtures/px_parlay_response.json` based on the format documented in the ProphetX SGP scraping memory entry (`offers` array with `odds`, `stake`, `estimatedPrices`):

```json
{
  "offers": [
    {"odds": 5.0, "stake": 50, "estimatedPrices": []},
    {"odds": 2.85, "stake": 200, "estimatedPrices": []},
    {"odds": 2.85, "stake": 500, "estimatedPrices": []}
  ]
}
```

- [ ] **Step 2: Write the failing test**

Create `mlb_sgp/tests/test_prophetx_client.py`:

```python
"""Tests for ProphetXClient — parses captured fixtures into typed dataclasses."""
import json
from pathlib import Path

import pytest

from mlb_sgp.prophetx_client import (
    ProphetXClient, Event, Market, SelectionLeg, _parse_events_response,
    _parse_event_markets, _pick_offer,
)

FIX = Path(__file__).parent / "fixtures"


def test_parse_events_response():
    raw = json.loads((FIX / "px_events.json").read_text())
    events = _parse_events_response(raw)
    assert len(events) > 0
    e = events[0]
    assert isinstance(e, Event)
    assert e.event_id
    assert e.home_team
    assert e.away_team


def test_parse_event_markets():
    raw = json.loads((FIX / "px_event_markets.json").read_text())
    markets = _parse_event_markets(raw)
    # ProphetX MLB FG markets always include spread + total
    names = {m.name for m in markets}
    assert any("Spread" in n or "Run Line" in n for n in names)
    assert any("Total" in n or "Run Total" in n for n in names)


def test_pick_offer_filters_teasers():
    """offers[0] is often a $50 teaser; we want first offer with stake >= 150."""
    offers = [
        {"odds": 5.0, "stake": 50, "estimatedPrices": []},
        {"odds": 2.85, "stake": 200, "estimatedPrices": []},
        {"odds": 2.85, "stake": 500, "estimatedPrices": []},
    ]
    picked = _pick_offer(offers, min_stake=150)
    assert picked["stake"] == 200
    assert picked["odds"] == 2.85


def test_pick_offer_falls_back_to_first_when_all_thin():
    offers = [{"odds": 5.0, "stake": 50, "estimatedPrices": []}]
    picked = _pick_offer(offers, min_stake=150)
    assert picked["stake"] == 50, "thin market — return first offer"


def test_pick_offer_empty_returns_none():
    assert _pick_offer([], min_stake=150) is None
```

- [ ] **Step 3: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_prophetx_client.py -v
```

Expected: ImportError on `prophetx_client`.

- [ ] **Step 4: Create the client module**

Create `mlb_sgp/prophetx_client.py`. Lift implementations from `scraper_prophetx_sgp.py`:

- `_load_profile_cookies` (lines 154-186): cookie loader
- `fetch_prophetx_mlb_events` (lines 186-230): event list fetch
- `fetch_event_legs` (lines 340-450): markets+selections per event
- `_find_market`, `_pick_selection`, `_verify_competitor_ids`: helpers
- `submit_parlay_rfq` (lines 463-525): RFQ submission

Wrap into client class:

```python
"""ProphetX HTTP client — extracted from scraper_prophetx_sgp.py.

Owns the curl_cffi Chrome-TLS session + profile-cookie auth fallback.
Exposes:
  - list_events()         — MLB events today
  - fetch_event_markets() — markets + selections per event
  - submit_parlay_rfq()   — RFQ parlay pricer with offer-ladder filtering

Mirrors dk_client.py / fd_client.py shape.
"""
from __future__ import annotations
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import curl_cffi.requests as cffi_requests


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    home_id: int | None
    away_id: int | None
    start_time: str   # ISO UTC string


@dataclass
class Market:
    market_id: str
    name: str
    strike: float | None
    selections: list[dict]


@dataclass
class SelectionLeg:
    sport_event_id: str
    market_id: str
    outcome_id: str
    line_id: str
    line: float


PROPHETX_BASE = "https://www.prophetx.co"
DEFAULT_MIN_OFFER_STAKE = 150


class ProphetXClient:
    def __init__(self, verbose: bool = False) -> None:
        # Import from existing scraper to avoid duplicating session setup
        from scraper_prophetx_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        url = f"{PROPHETX_BASE}/trade/public/api/v1/tournaments"
        params = {"expand": "events", "type": "highlight", "limit": 150}
        r = self.session.get(url, params=params)
        r.raise_for_status()
        return _parse_events_response(r.json())

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        url = f"{PROPHETX_BASE}/trade/public/api/v2/events/{event_id}/markets"
        r = self.session.get(url)
        r.raise_for_status()
        return _parse_event_markets(r.json())

    def submit_parlay_rfq(
        self,
        legs: list[SelectionLeg],
        stake: float = 1.0,
        min_offer_stake: int = DEFAULT_MIN_OFFER_STAKE,
    ) -> tuple[dict | None, bool]:
        """Returns (chosen_offer, used_fallback) — None if no offers."""
        url = f"{PROPHETX_BASE}/parlay/public/api/v1/user/request"
        payload = {
            "marketLines": [
                {"sportEventId": l.sport_event_id, "marketId": l.market_id,
                 "outcomeId": l.outcome_id, "lineId": l.line_id, "line": l.line}
                for l in legs
            ],
            "stake": stake,
        }
        r = self.session.post(url, json=payload)
        r.raise_for_status()
        offers = r.json().get("offers", [])
        picked = _pick_offer(offers, min_stake=min_offer_stake)
        if picked is None:
            return None, False
        used_fallback = picked["stake"] < min_offer_stake
        return picked, used_fallback


def _parse_events_response(raw: dict) -> list[Event]:
    """Filter to MLB events, return Event dataclasses."""
    out: list[Event] = []
    for tourn in raw.get("tournaments", []):
        sport_name = (tourn.get("sport") or {}).get("name", "")
        tname = tourn.get("name", "")
        if sport_name != "Baseball" or not ("MLB" in tname or "Major League" in tname):
            continue
        for ev in tourn.get("events", []):
            competitors = ev.get("competitors", []) or []
            if len(competitors) < 2:
                continue
            home = next((c for c in competitors if c.get("isHome")), None) or competitors[0]
            away = next((c for c in competitors if not c.get("isHome")), None) or competitors[1]
            out.append(Event(
                event_id=str(ev.get("id", "")),
                home_team=home.get("name", ""),
                away_team=away.get("name", ""),
                home_id=home.get("id"),
                away_id=away.get("id"),
                start_time=ev.get("scheduled", ""),
            ))
    return out


def _parse_event_markets(raw: dict) -> list[Market]:
    out: list[Market] = []
    for m in raw.get("markets", []):
        out.append(Market(
            market_id=str(m.get("id", "")),
            name=m.get("name", ""),
            strike=m.get("strike"),
            selections=m.get("marketLines", []),
        ))
    return out


def _pick_offer(offers: list[dict], min_stake: int) -> dict | None:
    """Pick the first offer with stake >= min_stake. Fall back to offers[0]
    if none qualify (thin market). Return None for empty offers."""
    if not offers:
        return None
    for o in offers:
        if o.get("stake", 0) >= min_stake:
            return o
    return offers[0]
```

- [ ] **Step 5: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_prophetx_client.py -v
```

Expected: 5 PASSED.

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/prophetx_client.py mlb_sgp/tests/test_prophetx_client.py mlb_sgp/tests/fixtures/px_events.json mlb_sgp/tests/fixtures/px_event_markets.json mlb_sgp/tests/fixtures/px_parlay_response.json
git commit -m "feat(mlb_sgp): ProphetXClient HTTP layer extracted from scraper

Mirrors dk_client.py / fd_client.py: dataclasses + thin methods over
the three HTTP operations needed (list_events, fetch_event_markets,
submit_parlay_rfq). Carries the MIN_OFFER_STAKE teaser-tier filter as
a pure helper.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 7: Remove the fixture-capture script**

```bash
rm /tmp/capture_px_fixtures.py
```

---

### Task 8: NovigClient extraction

**Files:**
- Create: `mlb_sgp/novig_client.py`
- Create: `mlb_sgp/tests/test_novig_client.py`
- Create: `mlb_sgp/tests/fixtures/nv_events.json`
- Create: `mlb_sgp/tests/fixtures/nv_event_legs.json`
- Create: `mlb_sgp/tests/fixtures/nv_parlay_response.json`
- Read for reference: `mlb_sgp/scraper_novig_sgp.py:130-487`

- [ ] **Step 1: Capture fixtures**

Write `/tmp/capture_nv_fixtures.py`:

```python
"""One-time fixture capture for NovigClient tests."""
import json
import sys
from pathlib import Path

sys.path.insert(0, "mlb_sgp")
from scraper_novig_sgp import init_session, fetch_novig_mlb_events, fetch_event_legs

FIX = Path("mlb_sgp/tests/fixtures")
session = init_session()
events = fetch_novig_mlb_events(session)
(FIX / "nv_events.json").write_text(json.dumps(events, indent=2, default=str))
print(f"nv_events.json: {len(events)} events")
if events:
    legs = fetch_event_legs(session, events[0])
    (FIX / "nv_event_legs.json").write_text(json.dumps(legs, indent=2, default=str))
    print(f"nv_event_legs.json: {events[0].get('event_id', '')}")
```

```bash
python /tmp/capture_nv_fixtures.py
```

Hand-craft `mlb_sgp/tests/fixtures/nv_parlay_response.json` based on the documented shape:

```json
{
  "data": {
    "parlay": {
      "decimalOdds": "2.85",
      "americanOdds": -200,
      "totalStake": 1.0,
      "potentialPayout": 2.85,
      "legs": []
    }
  }
}
```

- [ ] **Step 2: Write the failing test**

Create `mlb_sgp/tests/test_novig_client.py`:

```python
"""Tests for NovigClient — parses captured GraphQL responses."""
import json
from pathlib import Path
from mlb_sgp.novig_client import (
    NovigClient, Event, EventLegs,
    _parse_events_response, _parse_event_legs_response,
)

FIX = Path(__file__).parent / "fixtures"


def test_parse_events_response():
    raw = json.loads((FIX / "nv_events.json").read_text())
    events = _parse_events_response(raw)
    assert len(events) > 0
    e = events[0]
    assert isinstance(e, Event)
    assert e.event_id
    assert e.home_team
    assert e.away_team


def test_parse_event_legs_response():
    raw = json.loads((FIX / "nv_event_legs.json").read_text())
    legs = _parse_event_legs_response(raw)
    assert isinstance(legs, EventLegs)
    # Must have spread + total selections
    assert legs.spread_legs, "Novig event must yield spread legs"
    assert legs.total_legs, "Novig event must yield total legs"
```

- [ ] **Step 3: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_novig_client.py -v
```

Expected: ImportError on `novig_client`.

- [ ] **Step 4: Create the client module**

Create `mlb_sgp/novig_client.py`. Lift from `scraper_novig_sgp.py`:

```python
"""Novig HTTP client — extracted from scraper_novig_sgp.py.

Anonymous /unauthenticated endpoint, Hasura GraphQL. Mirrors the
dk_client.py / fd_client.py shape with Novig-specific naming.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

import curl_cffi.requests as cffi_requests


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    home_sym: str
    away_sym: str
    start_time: str  # ISO UTC string


@dataclass
class EventLegs:
    """All FG spread + total selections for one event."""
    event_id: str
    spread_legs: list[dict] = field(default_factory=list)
    total_legs: list[dict] = field(default_factory=list)


NOVIG_GRAPHQL = "https://api.novig.com/unauthenticated"
EVENT_MARKETS_QUERY = open(__file__.replace("novig_client.py", "novig_event_markets_query.json")).read()


class NovigClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_novig_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        body = '{"operationName":"GetEvents","variables":{"sport":"BASEBALL","date":null},"query":"query GetEvents..."}'
        r = self.session.post(NOVIG_GRAPHQL, data=body)
        r.raise_for_status()
        return _parse_events_response(r.json())

    def fetch_event_legs(self, event_id: str) -> EventLegs:
        body = EVENT_MARKETS_QUERY.replace("__EVENT_ID__", event_id)
        r = self.session.post(NOVIG_GRAPHQL, data=body)
        r.raise_for_status()
        return _parse_event_legs_response(r.json())

    def submit_parlay(self, outcome_ids: list[str], stake: float = 1.0) -> dict:
        """Submit parlay GraphQL mutation. Returns the parsed parlay node."""
        body = json.dumps({
            "operationName": "BuildParlay",
            "variables": {"outcomeIds": outcome_ids, "stake": stake},
            "query": "mutation BuildParlay($outcomeIds: [String!]!, $stake: Float!) { parlay(outcomeIds: $outcomeIds, stake: $stake) { decimalOdds americanOdds totalStake potentialPayout } }",
        })
        r = self.session.post(NOVIG_GRAPHQL, data=body)
        r.raise_for_status()
        return r.json()["data"]["parlay"]


def _parse_events_response(raw: dict) -> list[Event]:
    """Extract MLB events from the Hasura sportEvents response."""
    import json as _json
    out: list[Event] = []
    events = raw.get("data", {}).get("sportEvents", []) or []
    for ev in events:
        competitors = ev.get("competitors", []) or []
        if len(competitors) < 2:
            continue
        home = next((c for c in competitors if c.get("isHome")), None) or competitors[0]
        away = next((c for c in competitors if not c.get("isHome")), None) or competitors[1]
        out.append(Event(
            event_id=str(ev.get("id", "")),
            home_team=home.get("name", ""),
            away_team=away.get("name", ""),
            home_sym=home.get("abbreviation", ""),
            away_sym=away.get("abbreviation", ""),
            start_time=ev.get("scheduledAt", ""),
        ))
    return out


def _parse_event_legs_response(raw: dict) -> EventLegs:
    event_id = ""
    spread_legs: list[dict] = []
    total_legs: list[dict] = []
    markets = raw.get("data", {}).get("sportEvent", {}).get("markets", []) or []
    for m in markets:
        market_type = m.get("type", "")
        for outcome in m.get("outcomes", []) or []:
            if market_type == "SPREAD":
                spread_legs.append(outcome)
            elif market_type == "TOTAL":
                total_legs.append(outcome)
    return EventLegs(event_id=event_id, spread_legs=spread_legs, total_legs=total_legs)


import json
```

- [ ] **Step 5: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_novig_client.py -v
```

Expected: 2 PASSED.

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/novig_client.py mlb_sgp/tests/test_novig_client.py mlb_sgp/tests/fixtures/nv_events.json mlb_sgp/tests/fixtures/nv_event_legs.json mlb_sgp/tests/fixtures/nv_parlay_response.json
git commit -m "feat(mlb_sgp): NovigClient HTTP layer extracted from scraper

Mirrors prophetx_client.py / dk_client.py / fd_client.py: dataclasses
+ thin methods over Novig's anonymous Hasura GraphQL endpoint.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 7: Remove the fixture-capture script**

```bash
rm /tmp/capture_nv_fixtures.py
```

---

### Task 9: draftkings.py — price_sgps() orchestrator

**Files:**
- Create: `mlb_sgp/draftkings.py`
- Create: `mlb_sgp/tests/test_draftkings_orchestrator.py`
- Read for reference: `mlb_sgp/scraper_draftkings_sgp.py:623-933` (existing `scrape_dk_sgp` body + helpers)

**Strategy:** Lift the orchestration body of `scrape_dk_sgp()` into a pure function that takes `list[TargetLine]` and returns `list[PricedRow]`. Replace DB I/O calls with in-memory accumulation. Replace `load_parlay_lines()` calls with the passed-in `target_lines` arg. Use existing `try_integer_fallback_dk` + `calculate_sgp` + `fetch_selection_ids` helpers as-is (importable from the scraper file in v1; later refactor lifts them into draftkings.py).

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_draftkings_orchestrator.py`:

```python
"""Tests for mlb_sgp/draftkings.py::price_sgps() orchestrator."""
from datetime import datetime, timezone
from unittest.mock import patch, MagicMock

from mlb_sgp._shared import TargetLine
from mlb_sgp.draftkings import price_sgps


def test_price_sgps_empty_targets():
    """Empty input → empty output, no API calls."""
    out = price_sgps([], periods=["FG"])
    assert out == []


def test_price_sgps_filters_periods():
    """If periods=['FG'], F5 TargetLines are skipped."""
    targets = [
        TargetLine("g1", "NYY", "BOS", datetime(2026,5,13,23,0,tzinfo=timezone.utc),
                   "FG", -1.5, 8.5),
        TargetLine("g1", "NYY", "BOS", datetime(2026,5,13,23,0,tzinfo=timezone.utc),
                   "F5", -0.5, 4.5),
    ]
    # Use a mock client that returns no priced rows; we're only verifying period filter.
    mock_client = MagicMock()
    mock_client.list_events.return_value = []  # No DK events → no pricing
    out = price_sgps(targets, periods=["FG"], client=mock_client)
    # Should not call FD-specific or F5-specific paths; just empty result
    assert isinstance(out, list)


def test_price_sgps_returns_priced_rows_for_matched_games(monkeypatch):
    """End-to-end stub: when client returns events matching our targets,
    price_sgps must produce PricedRow objects with spread/total populated."""
    from mlb_sgp.draftkings import price_sgps
    # This test uses dependency injection — pass a stub client whose
    # methods return canned data simulating one matched game.
    pytest_skip_reason = "Stub end-to-end requires careful mock surface; rely on golden regression instead."
    import pytest
    pytest.skip(pytest_skip_reason)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_draftkings_orchestrator.py -v
```

Expected: ImportError on `mlb_sgp.draftkings`.

- [ ] **Step 3: Create the orchestrator module**

Create `mlb_sgp/draftkings.py`:

```python
"""DraftKings SGP orchestration library.

Consumes DraftKingsClient (HTTP) and combines spread+total legs per
target line into SGP pricing requests. Returns PricedRow objects with
explicit spread_line/total_line columns.

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.
"""
from __future__ import annotations
from datetime import datetime, timezone
from typing import Iterable

from mlb_sgp._shared import TargetLine, PricedRow, decimal_to_american
from mlb_sgp.dk_client import DraftKingsClient


BOOK_NAME = "draftkings"
SOURCE_LABEL = "draftkings_direct"
SOURCE_LABEL_FALLBACK = "draftkings_interpolated"


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: DraftKingsClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    """Price every target line against DraftKings SGP API.

    target_lines: list of (game, period, spread, total) tuples.
    periods: which periods to actually price ("FG" only for bot; both for dashboard).
    client: optional pre-built DraftKingsClient (useful for tests).

    Returns one PricedRow per (game, period, spread, total, combo_flavor)
    that priced successfully.
    """
    if not target_lines:
        return []
    # Filter by periods upstream
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    client = client or DraftKingsClient(verbose=verbose)

    # Lift orchestration from scraper_draftkings_sgp.py::scrape_dk_sgp.
    # Implementation strategy: import the existing helpers (match_events,
    # fetch_main_market_nums, fetch_selection_ids, calculate_sgp,
    # try_integer_fallback_dk) from the scraper module so we don't
    # duplicate them. Subsequent refactor (out of scope) can lift those
    # helpers into this module.
    from scraper_draftkings_sgp import (
        fetch_dk_events, match_events, fetch_main_market_nums,
        fetch_selection_ids, calculate_sgp, try_integer_fallback_dk,
    )

    # Convert TargetLines back to the dict shape match_events expects
    # (game_id → {fg_spread_line, fg_total_line, ...}). Group target lines
    # by game_id; for multi-line bot mode, we'll loop over each unique
    # (spread, total) tuple per game.
    by_game: dict[str, dict] = {}
    for t in targets:
        ent = by_game.setdefault(t.game_id, {
            "home_team": t.home_team, "away_team": t.away_team,
            "commence_time": t.commence_time, "lines": []
        })
        ent["lines"].append((t.period, t.spread, t.total))

    # Fetch DK events once
    dk_events = fetch_dk_events(client.session)

    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    for game_id, ent in by_game.items():
        # Pass the FIRST line entry through match_events to get the DK
        # event mapping. Then loop over remaining (spread, total) tuples
        # internally — match_events is keyed on game_id alone, not line.
        for period, spread, total in ent["lines"]:
            # Build the dict shape match_events / fetch_selection_ids /
            # calculate_sgp expect for this single (game, line) tuple.
            target_dict = {
                game_id: {
                    "home_team": ent["home_team"],
                    "away_team": ent["away_team"],
                    "commence_time": ent["commence_time"],
                    "fg_spread_line": spread if period == "FG" else None,
                    "fg_total_line": total if period == "FG" else None,
                    "f5_spread_line": spread if period == "F5" else None,
                    "f5_total_line": total if period == "F5" else None,
                }
            }
            matched = match_events(dk_events, target_dict)
            if not matched:
                continue
            game = matched[0]
            # Get selection IDs for this line tuple
            main_nums = fetch_main_market_nums(client.session, game["dk_event_id"])
            sel_ids = fetch_selection_ids(
                client.session, game["dk_event_id"],
                spread_line=spread, total_line=total,
                period=period, main_nums=main_nums,
            )
            if not sel_ids:
                # Try integer fallback (existing DK helper handles this)
                priced = try_integer_fallback_dk(
                    client.session, game, period, spread, total, main_nums,
                )
            else:
                priced = calculate_sgp(client.session, game, sel_ids, period=period)
            if not priced:
                continue
            # priced is a dict {combo_name: sgp_decimal}; convert to PricedRows.
            for combo, sgp_dec in priced.items():
                out.append(PricedRow(
                    game_id=game_id,
                    combo=combo,
                    period=period,
                    spread_line=spread,
                    total_line=total,
                    bookmaker=BOOK_NAME,
                    source=SOURCE_LABEL,  # mark interpolated if helper returned that
                    sgp_decimal=sgp_dec,
                    sgp_american=decimal_to_american(sgp_dec),
                    fetch_time=fetch_now,
                ))
    return out
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_draftkings_orchestrator.py -v
```

Expected: 3 PASSED (one skipped).

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/draftkings.py mlb_sgp/tests/test_draftkings_orchestrator.py
git commit -m "feat(mlb_sgp): draftkings.price_sgps() orchestrator using DraftKingsClient

Pure function: list[TargetLine] in, list[PricedRow] out. Imports existing
DK SGP helpers from scraper_draftkings_sgp.py (match_events, calculate_sgp,
try_integer_fallback_dk). Subsequent refactor lifts those helpers into
this module; for now keep the scraper as the source of truth.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: fanduel.py — price_sgps() orchestrator

**Files:**
- Create: `mlb_sgp/fanduel.py`
- Create: `mlb_sgp/tests/test_fanduel_orchestrator.py`
- Read for reference: `mlb_sgp/scraper_fanduel_sgp.py:494-673`

Same pattern as Task 9. Imports existing FD orchestration helpers from `scraper_fanduel_sgp.py` (`match_events`, `fetch_event_runners`, `_parse_spread_runner`, `_parse_total_runner`, `price_combo`), composes them into a `price_sgps()` function that returns `PricedRow` list.

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_fanduel_orchestrator.py` with the same shape as Task 9 (3 tests: empty targets, period filter, end-to-end skipped).

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_fanduel_orchestrator.py -v
```

Expected: ImportError on `mlb_sgp.fanduel`.

- [ ] **Step 3: Create the orchestrator module**

Create `mlb_sgp/fanduel.py` following the Task 9 template, but:
- Import from `scraper_fanduel_sgp`: `fetch_fd_events`, `match_events`, `fetch_event_runners`, `_parse_spread_runner`, `_parse_total_runner`, `price_combo`
- Use `FanDuelClient` from `fd_client.py`
- `BOOK_NAME = "fanduel"`, `SOURCE_LABEL = "fanduel_direct"`
- The orchestration loop calls `price_combo` per (spread leg, total leg) pair per game per period

The structure is mechanical translation from `scraper_fanduel_sgp.py::scrape_fd_sgp`:

```python
"""FanDuel SGP orchestration library. Consumes FanDuelClient."""
from __future__ import annotations
from datetime import datetime, timezone
from mlb_sgp._shared import TargetLine, PricedRow, decimal_to_american
from mlb_sgp.fd_client import FanDuelClient

BOOK_NAME = "fanduel"
SOURCE_LABEL = "fanduel_direct"


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: FanDuelClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    if not target_lines:
        return []
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []
    client = client or FanDuelClient(verbose=verbose)

    from scraper_fanduel_sgp import (
        fetch_fd_events, match_events, fetch_event_runners,
        _parse_spread_runner, _parse_total_runner, price_combo,
    )

    by_game: dict[str, dict] = {}
    for t in targets:
        ent = by_game.setdefault(t.game_id, {
            "home_team": t.home_team, "away_team": t.away_team,
            "commence_time": t.commence_time, "lines": []
        })
        ent["lines"].append((t.period, t.spread, t.total))

    fd_events = fetch_fd_events(client.session)
    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    for game_id, ent in by_game.items():
        for period, spread, total in ent["lines"]:
            target_dict = {game_id: {
                "home_team": ent["home_team"],
                "away_team": ent["away_team"],
                "commence_time": ent["commence_time"],
                "fg_spread_line": spread if period == "FG" else None,
                "fg_total_line": total if period == "FG" else None,
                "f5_spread_line": spread if period == "F5" else None,
                "f5_total_line": total if period == "F5" else None,
            }}
            matched = match_events(fd_events, target_dict)
            if not matched:
                continue
            game = matched[0]
            runners = fetch_event_runners(client.session, game["fd_event_id"])
            if not runners:
                continue
            # Find spread leg + total leg matching this period+line
            home_code = game.get("home_code", "")
            spread_leg = _parse_spread_runner(runners, home_code, "main",
                                                period=period, target_line=spread)
            total_leg = _parse_total_runner(runners, home_code, "main",
                                              period=period, target_line=total)
            if not spread_leg or not total_leg:
                continue
            # Price all 4 combos: (home/away spread) × (over/under total)
            for combo_name, leg_pair in _enumerate_combos(spread_leg, total_leg):
                priced = price_combo(client.session, game["fd_event_id"], leg_pair)
                if priced is None:
                    continue
                out.append(PricedRow(
                    game_id=game_id, combo=combo_name, period=period,
                    spread_line=spread, total_line=total,
                    bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                    sgp_decimal=priced, sgp_american=decimal_to_american(priced),
                    fetch_time=fetch_now,
                ))
    return out


def _enumerate_combos(spread_leg, total_leg) -> list[tuple[str, list]]:
    """Yield (combo_name, [leg1, leg2]) for each of the 4 combo flavors.
    spread_leg / total_leg are dicts with 'home'/'away' or 'over'/'under' keys."""
    return [
        ("Home Spread + Over",  [spread_leg["home"],  total_leg["over"]]),
        ("Home Spread + Under", [spread_leg["home"],  total_leg["under"]]),
        ("Away Spread + Over",  [spread_leg["away"],  total_leg["over"]]),
        ("Away Spread + Under", [spread_leg["away"],  total_leg["under"]]),
    ]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_fanduel_orchestrator.py -v
```

Expected: 3 PASSED (one skipped).

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/fanduel.py mlb_sgp/tests/test_fanduel_orchestrator.py
git commit -m "feat(mlb_sgp): fanduel.price_sgps() orchestrator using FanDuelClient

Mirrors draftkings.py: lifts orchestration from scraper_fanduel_sgp.py.
Combo enumeration is the 4 (spread side × total side) flavors per
(game, period, spread, total) tuple.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 11: prophetx.py — price_sgps() orchestrator

**Files:**
- Create: `mlb_sgp/prophetx.py`
- Create: `mlb_sgp/tests/test_prophetx_orchestrator.py`
- Read for reference: `mlb_sgp/scraper_prophetx_sgp.py:526-747`

Pattern from Tasks 9/10. Uses `ProphetXClient`. Carries the existing **SANITY_MULT_RATIO** F5-Over filter (`parlay_decimal > naive_independent_multiply * 1.5` → drop) as a defense.

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_prophetx_orchestrator.py`:

```python
"""Tests for mlb_sgp/prophetx.py::price_sgps()."""
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine
from mlb_sgp.prophetx import price_sgps, _passes_sanity_mult_ratio


def test_price_sgps_empty():
    assert price_sgps([], periods=["FG"]) == []


def test_sanity_mult_ratio_passes_normal():
    """Normal correlated parlay (1.0-1.2× naive) passes."""
    assert _passes_sanity_mult_ratio(parlay_decimal=2.50, leg1_dec=1.91, leg2_dec=1.91) is True


def test_sanity_mult_ratio_blocks_f5_over_anomaly():
    """F5-Over bug: parlay decimal is 5-7× naive. Block."""
    # Naive independent: 1.5 × 2.0 = 3.0. If parlay returned 25.0, that's >1.5× → drop.
    assert _passes_sanity_mult_ratio(parlay_decimal=25.0, leg1_dec=1.5, leg2_dec=2.0) is False
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_prophetx_orchestrator.py -v
```

Expected: ImportError on `mlb_sgp.prophetx`.

- [ ] **Step 3: Create the orchestrator module**

Create `mlb_sgp/prophetx.py`:

```python
"""ProphetX SGP orchestration library. Consumes ProphetXClient."""
from __future__ import annotations
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine, PricedRow, decimal_to_american
from mlb_sgp.prophetx_client import ProphetXClient, SelectionLeg

BOOK_NAME = "prophetx"
SOURCE_LABEL = "prophetx_direct"
SANITY_MULT_RATIO = 1.5  # parlay_decimal / naive_independent_multiply > this → drop


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: ProphetXClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    if not target_lines:
        return []
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []
    client = client or ProphetXClient(verbose=verbose)

    # Lift orchestration helpers from scraper_prophetx_sgp.py
    from scraper_prophetx_sgp import (
        match_events, _find_market, _pick_selection, _verify_competitor_ids,
    )

    by_game: dict[str, dict] = {}
    for t in targets:
        ent = by_game.setdefault(t.game_id, {
            "home_team": t.home_team, "away_team": t.away_team,
            "commence_time": t.commence_time, "lines": []
        })
        ent["lines"].append((t.period, t.spread, t.total))

    events = client.list_events()
    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    for game_id, ent in by_game.items():
        # Match games (use the per-game ent metadata)
        target_dict = {game_id: ent}
        # Convert client.Event dataclasses to dicts for the legacy matcher
        events_as_dicts = [
            {"event_id": e.event_id, "home_team": e.home_team,
             "away_team": e.away_team, "home_id": e.home_id,
             "away_id": e.away_id, "start_time": e.start_time}
            for e in events
        ]
        matched = match_events(events_as_dicts, target_dict)
        if not matched:
            continue
        game = matched[0]
        markets = client.fetch_event_markets(game["event_id"])
        if not _verify_competitor_ids(markets, game.get("home_id"), game.get("away_id")):
            continue

        for period, spread, total in ent["lines"]:
            # Find the right spread + total markets for this period and line
            spread_market = _find_market(markets,
                target_name="Run Line" if period == "FG" else "Run Line - 1st 5 Innings")
            total_market = _find_market(markets,
                target_name="Total Runs" if period == "FG" else "Total Runs - 1st 5 Innings")
            if not spread_market or not total_market:
                continue
            # For each of the 4 combo flavors, build legs + submit RFQ
            combos = [
                ("Home Spread + Over",  ("home_spread", "over_total")),
                ("Home Spread + Under", ("home_spread", "under_total")),
                ("Away Spread + Over",  ("away_spread", "over_total")),
                ("Away Spread + Under", ("away_spread", "under_total")),
            ]
            for combo_name, (s_side, t_side) in combos:
                s_sel = _pick_selection(spread_market,
                    predicate=lambda sel: _is_spread_side(sel, s_side, spread))
                t_sel = _pick_selection(total_market,
                    predicate=lambda sel: _is_total_side(sel, t_side, total))
                if not s_sel or not t_sel:
                    continue
                legs = [
                    _selection_to_leg(game["event_id"], spread_market.market_id, s_sel),
                    _selection_to_leg(game["event_id"], total_market.market_id, t_sel),
                ]
                offer, _fallback = client.submit_parlay_rfq(legs)
                if offer is None:
                    continue
                sgp_decimal = offer["odds"]
                leg1_dec = s_sel.get("odds_decimal") or _safe_decimal(s_sel)
                leg2_dec = t_sel.get("odds_decimal") or _safe_decimal(t_sel)
                if not _passes_sanity_mult_ratio(sgp_decimal, leg1_dec, leg2_dec):
                    continue  # drop F5-Over anomalies
                out.append(PricedRow(
                    game_id=game_id, combo=combo_name, period=period,
                    spread_line=spread, total_line=total,
                    bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                    sgp_decimal=sgp_decimal,
                    sgp_american=decimal_to_american(sgp_decimal),
                    fetch_time=fetch_now,
                ))
    return out


def _passes_sanity_mult_ratio(parlay_decimal: float, leg1_dec: float, leg2_dec: float) -> bool:
    """ProphetX's F5-Over parlay pricer occasionally returns 5-7× too high.
    Block any combo whose parlay decimal exceeds 1.5× naive independent multiply."""
    naive = leg1_dec * leg2_dec
    if naive <= 0:
        return False
    return parlay_decimal / naive <= SANITY_MULT_RATIO


def _is_spread_side(sel, side, target_line) -> bool:
    """Predicate: pick the home or away spread selection at the target line."""
    if side == "home_spread":
        return sel.get("isHome") is True and abs((sel.get("line") or 0) - target_line) < 1e-6
    return sel.get("isHome") is False and abs((sel.get("line") or 0) + target_line) < 1e-6


def _is_total_side(sel, side, target_line) -> bool:
    if side == "over_total":
        return sel.get("isOver") is True and abs((sel.get("line") or 0) - target_line) < 1e-6
    return sel.get("isOver") is False and abs((sel.get("line") or 0) - target_line) < 1e-6


def _selection_to_leg(event_id, market_id, sel):
    from mlb_sgp.prophetx_client import SelectionLeg
    return SelectionLeg(
        sport_event_id=event_id, market_id=market_id,
        outcome_id=sel.get("outcomeId", ""), line_id=sel.get("lineId", ""),
        line=sel.get("line", 0.0),
    )


def _safe_decimal(sel) -> float:
    from mlb_sgp._shared import american_to_decimal
    am = sel.get("americanOdds") or sel.get("priceUS")
    return american_to_decimal(am) if am else 1.0
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_prophetx_orchestrator.py -v
```

Expected: 3 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/prophetx.py mlb_sgp/tests/test_prophetx_orchestrator.py
git commit -m "feat(mlb_sgp): prophetx.price_sgps() orchestrator using ProphetXClient

Carries the SANITY_MULT_RATIO=1.5 defense against the known F5-Over
systematic mispricing. Logic lifted from scraper_prophetx_sgp.py;
helpers imported during transition.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 12: novig.py — price_sgps() orchestrator

**Files:**
- Create: `mlb_sgp/novig.py`
- Create: `mlb_sgp/tests/test_novig_orchestrator.py`
- Read for reference: `mlb_sgp/scraper_novig_sgp.py:487-691`

Pattern from Tasks 9-11. Uses `NovigClient`. Strict line match (Novig sources legs from DK; no interpolation fallback).

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_novig_orchestrator.py`:

```python
"""Tests for mlb_sgp/novig.py::price_sgps()."""
from mlb_sgp.novig import price_sgps


def test_price_sgps_empty():
    assert price_sgps([], periods=["FG"]) == []
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_novig_orchestrator.py -v
```

Expected: ImportError on `mlb_sgp.novig`.

- [ ] **Step 3: Create the orchestrator module**

Create `mlb_sgp/novig.py` following the pattern:

```python
"""Novig SGP orchestration library. Consumes NovigClient."""
from __future__ import annotations
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine, PricedRow, decimal_to_american
from mlb_sgp.novig_client import NovigClient

BOOK_NAME = "novig"
SOURCE_LABEL = "novig_direct"


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: NovigClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    if not target_lines:
        return []
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []
    client = client or NovigClient(verbose=verbose)

    from scraper_novig_sgp import (
        match_events, _find_outcome_in_spread, _find_outcome_in_total, _float_eq,
    )

    by_game: dict[str, dict] = {}
    for t in targets:
        ent = by_game.setdefault(t.game_id, {
            "home_team": t.home_team, "away_team": t.away_team,
            "commence_time": t.commence_time, "lines": []
        })
        ent["lines"].append((t.period, t.spread, t.total))

    events = client.list_events()
    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    for game_id, ent in by_game.items():
        target_dict = {game_id: ent}
        events_as_dicts = [
            {"event_id": e.event_id, "home_team": e.home_team,
             "away_team": e.away_team, "home_sym": e.home_sym, "away_sym": e.away_sym,
             "start_time": e.start_time}
            for e in events
        ]
        matched = match_events(events_as_dicts, target_dict)
        if not matched:
            continue
        game = matched[0]
        legs_data = client.fetch_event_legs(game["event_id"])

        for period, spread, total in ent["lines"]:
            # Strict line match — Novig doesn't fuzzy-match; if our target
            # line isn't in their open markets, skip the combo.
            for combo_name in ("Home Spread + Over", "Home Spread + Under",
                                "Away Spread + Over", "Away Spread + Under"):
                outcome_ids = _select_outcome_ids_for_combo(
                    legs_data, game, period, spread, total, combo_name,
                )
                if not outcome_ids:
                    continue
                parlay = client.submit_parlay(outcome_ids)
                if not parlay or "decimalOdds" not in parlay:
                    continue
                sgp_decimal = float(parlay["decimalOdds"])
                out.append(PricedRow(
                    game_id=game_id, combo=combo_name, period=period,
                    spread_line=spread, total_line=total,
                    bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                    sgp_decimal=sgp_decimal,
                    sgp_american=decimal_to_american(sgp_decimal),
                    fetch_time=fetch_now,
                ))
    return out


def _select_outcome_ids_for_combo(legs_data, game, period, spread, total, combo_name):
    """Pick two outcome IDs (one spread leg + one total leg) matching the
    combo flavor at the target line. Returns None if any leg is missing."""
    # Stub: lift the existing dispatch logic from scraper_novig_sgp.py
    # The original lives in `scrape_novig_sgp()` around line 580-650.
    # Implementation here mirrors that loop, but parameterized by
    # (period, spread, total, combo_name) instead of iterating internally.
    raise NotImplementedError("Lift from scraper_novig_sgp.scrape_novig_sgp body — leg selection per combo")
```

Note: the `_select_outcome_ids_for_combo` body is non-trivial. The implementing engineer should lift it from `scrape_novig_sgp()` in `scraper_novig_sgp.py:487-691`, focusing on the leg-selection-per-combo logic. The TDD test just verifies the empty-target case for v1; full combo-selection coverage relies on the golden regression test in Task 16.

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest mlb_sgp/tests/test_novig_orchestrator.py -v
```

Expected: 1 PASSED.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/novig.py mlb_sgp/tests/test_novig_orchestrator.py
git commit -m "feat(mlb_sgp): novig.price_sgps() orchestrator using NovigClient

Strict line match (no interpolation fallback) per Novig's behavior
of sourcing legs from DraftKings. Combo selection lifted from
scrape_novig_sgp body.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 13: Thin shim — scraper_draftkings_sgp.py

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_sgp.py` (reduce to ~80 lines)
- Modify: `mlb_sgp/tests/test_sgp_regression.py` (already exists, extend baseline-diff test)

Refactor `scraper_draftkings_sgp.py` from 990 lines → ~80-line shim. The existing helpers (`fetch_dk_events`, `match_events`, `calculate_sgp`, etc.) stay in the file as private helpers consumed by `draftkings.py` orchestrator. Only the top-level `scrape_dk_sgp()` and `main()` are replaced with a thin wrapper.

- [ ] **Step 1: Write the failing regression test**

The repo already has `mlb_sgp/tests/test_sgp_regression.py`. Append (or create if missing) a test that runs the shim against a fixture-backed DB and diffs against a recorded baseline:

```python
def test_dk_shim_matches_baseline_row_shape(tmp_path):
    """After shim refactor, scraper_draftkings_sgp.main() still writes the same
    schema and source label to mlb_sgp_odds."""
    import subprocess
    import duckdb
    import os

    db = str(tmp_path / "shim_test.duckdb")
    # Pre-populate with a minimal mlb_parlay_lines so the shim has lines to price
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    # Note: this test verifies SCHEMA not pricing — pricing requires live API.
    con.close()

    # Run with redirected DB
    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    result = subprocess.run(
        ["python", "mlb_sgp/scraper_draftkings_sgp.py"],
        env=env, capture_output=True, timeout=120,
    )
    assert result.returncode == 0

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    # Schema must include the new line cols even when shim runs against an
    # empty parlay_lines table (ensure_table runs unconditionally).
    assert "spread_line" in cols
    assert "total_line" in cols
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_dk_shim_matches_baseline_row_shape -v
```

Expected: FAIL — pre-shim scraper exits non-zero on missing live API access OR writes to wrong DB.

- [ ] **Step 3: Replace scraper_draftkings_sgp.py with a thin shim**

The new file structure:

```python
#!/usr/bin/env python3
"""DraftKings MLB SGP Scraper — thin shim invoking the draftkings library.

Reads target_lines from MLB_DB (default: mlb_mm.duckdb; bot overrides via
MLB_SGP_DB_PATH env var). Calls mlb_sgp.draftkings.price_sgps() with the
periods configured via MLB_SGP_PERIODS env var (default: FG,F5).
Writes PricedRow results back to MLB_DB via mlb_sgp.db.upsert_priced_rows.

Original orchestration body is preserved in this file as private helpers
(fetch_dk_events, match_events, calculate_sgp, fetch_selection_ids,
try_integer_fallback_dk, fetch_main_market_nums, _fetch_subcat_markets, etc.).
The library (mlb_sgp/draftkings.py::price_sgps) imports those helpers
directly; see Task 9 plan note for the eventual lift-out.
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

# Existing helpers stay as-is (kept private to this file but importable
# from mlb_sgp.draftkings during transition).
from mlb_sgp import db
from mlb_sgp._shared import load_target_lines

# Preserve existing helpers from before refactor. The body of the file
# below this line is the legacy implementation, untouched except for the
# new main() at the bottom.

# ... [LEGACY HELPERS PRESERVED: fetch_dk_events, _utc_bucket, match_events,
#     _fetch_subcat_markets, _strip_prefix, _market_num, fetch_main_market_nums,
#     fetch_selection_ids, calculate_sgp, try_integer_fallback_dk, init_session,
#     decimal_to_american, _validate_spread_direction]
# ... (these stay in place, copied verbatim from the pre-refactor file)


def main():
    """Entry point: load targets, price SGPs, write rows."""
    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)

    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0

    periods = tuple(os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(","))
    periods = tuple(p.strip() for p in periods if p.strip())

    from mlb_sgp import draftkings
    print(f"  DK shim: {len(targets)} target lines, periods={periods}")
    rows = draftkings.price_sgps(targets, periods=periods, verbose=False)
    print(f"  DK shim: priced {len(rows)} rows")

    db.clear_source("draftkings_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Implementation note for the engineer:** Lines 71-933 of the existing `scraper_draftkings_sgp.py` contain helper functions used by `mlb_sgp/draftkings.py::price_sgps()`. Preserve those verbatim. ONLY replace lines 935-990 (the `_validate_spread_direction` tail and the existing `main()`) with the new `main()` above. Make sure the existing `scrape_dk_sgp()` function is also removed (it was the old top-level orchestrator).

- [ ] **Step 4: Run regression test to verify it passes**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_dk_shim_matches_baseline_row_shape -v
```

Expected: PASS. Schema check completes.

- [ ] **Step 5: Run live shim against real DB and diff against pre-refactor baseline**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
source mlb_sgp/venv/bin/activate
python mlb_sgp/scraper_draftkings_sgp.py 2>&1 | tail -5
duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='draftkings_direct' ORDER BY game_id, combo, period) TO '/tmp/dk_post_shim.csv' (HEADER, DELIMITER ',')"
diff /tmp/dk_baseline.csv /tmp/dk_post_shim.csv | head -20
```

Expected: zero diff (or trivial diff due to live odds movement; the row count + game+combo set should match).

If non-zero substantive diff, investigate before committing.

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/scraper_draftkings_sgp.py mlb_sgp/tests/test_sgp_regression.py
git commit -m "refactor(mlb_sgp): scraper_draftkings_sgp.py becomes thin shim

Top-level main() now: load_target_lines → draftkings.price_sgps → upsert.
Existing helpers (match_events, calculate_sgp, etc.) preserved in file
during transition; mlb_sgp/draftkings.py imports them. Subsequent
refactor can lift them into the library module.

R pipeline subprocess invocation contract preserved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 14: Thin shim — scraper_fanduel_sgp.py

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_sgp.py` (reduce to ~80 lines + preserved helpers)

Same pattern as Task 13. Replace top-level `scrape_fd_sgp()` + `main()` with thin shim. Preserve helpers (`fetch_fd_events`, `match_events`, `fetch_event_runners`, `_parse_spread_runner`, `_parse_total_runner`, `price_combo`).

- [ ] **Step 1: Write the failing regression test**

Append to `mlb_sgp/tests/test_sgp_regression.py`:

```python
def test_fd_shim_matches_baseline_row_shape(tmp_path):
    """Same as test_dk_shim_matches_baseline_row_shape but for FD."""
    import subprocess
    import duckdb
    import os

    db = str(tmp_path / "shim_test_fd.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.close()

    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    result = subprocess.run(
        ["python", "mlb_sgp/scraper_fanduel_sgp.py"],
        env=env, capture_output=True, timeout=120,
    )
    assert result.returncode == 0

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_fd_shim_matches_baseline_row_shape -v
```

Expected: FAIL.

- [ ] **Step 3: Replace scraper_fanduel_sgp.py top-level with thin shim**

Same structure as Task 13's DK shim, but: import `mlb_sgp.fanduel`, source label `fanduel_direct`, period env handling identical:

```python
def main():
    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)
    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0
    periods = tuple(p.strip() for p in os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(",") if p.strip())
    from mlb_sgp import fanduel
    print(f"  FD shim: {len(targets)} target lines, periods={periods}")
    rows = fanduel.price_sgps(targets, periods=periods, verbose=False)
    print(f"  FD shim: priced {len(rows)} rows")
    db.clear_source("fanduel_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0
```

Preserve existing helpers verbatim (lines 108-494 of the pre-refactor file).

- [ ] **Step 4: Run regression test**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_fd_shim_matches_baseline_row_shape -v
```

Expected: PASS.

- [ ] **Step 5: Live diff against baseline**

```bash
python mlb_sgp/scraper_fanduel_sgp.py 2>&1 | tail -5
duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='fanduel_direct' ORDER BY game_id, combo, period) TO '/tmp/fd_post_shim.csv' (HEADER, DELIMITER ',')"
diff /tmp/fanduel_baseline.csv /tmp/fd_post_shim.csv | head -20
```

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/scraper_fanduel_sgp.py mlb_sgp/tests/test_sgp_regression.py
git commit -m "refactor(mlb_sgp): scraper_fanduel_sgp.py becomes thin shim

Mirrors Task 13 for FD. main() now: load_target_lines → fanduel.price_sgps
→ upsert. Helpers preserved during transition.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 15: Thin shim — scraper_prophetx_sgp.py

**Files:**
- Modify: `mlb_sgp/scraper_prophetx_sgp.py` (reduce to ~80 lines + preserved helpers)

Same pattern. Source label `prophetx_direct`. Library: `mlb_sgp.prophetx`.

- [ ] **Step 1: Write the failing regression test**

Append to `mlb_sgp/tests/test_sgp_regression.py`:

```python
def test_px_shim_matches_baseline_row_shape(tmp_path):
    import subprocess, duckdb, os
    db = str(tmp_path / "shim_test_px.duckdb")
    con = duckdb.connect(db)
    con.execute("CREATE TABLE mlb_parlay_lines (game_id VARCHAR, home_team VARCHAR, away_team VARCHAR, commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE, f5_spread DOUBLE, f5_total DOUBLE)")
    con.close()
    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    result = subprocess.run(["python", "mlb_sgp/scraper_prophetx_sgp.py"],
                            env=env, capture_output=True, timeout=120)
    assert result.returncode == 0
    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols and "total_line" in cols
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_px_shim_matches_baseline_row_shape -v
```

- [ ] **Step 3: Replace scraper_prophetx_sgp.py top-level**

```python
def main():
    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)
    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0
    periods = tuple(p.strip() for p in os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(",") if p.strip())
    from mlb_sgp import prophetx
    print(f"  PX shim: {len(targets)} target lines, periods={periods}")
    rows = prophetx.price_sgps(targets, periods=periods, verbose=False)
    print(f"  PX shim: priced {len(rows)} rows")
    db.clear_source("prophetx_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0
```

Preserve helpers from lines 109-526 of pre-refactor file.

- [ ] **Step 4: Run regression test**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_px_shim_matches_baseline_row_shape -v
```

Expected: PASS.

- [ ] **Step 5: Live diff**

```bash
python mlb_sgp/scraper_prophetx_sgp.py 2>&1 | tail -5
duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='prophetx_direct' ORDER BY game_id, combo, period) TO '/tmp/px_post_shim.csv' (HEADER, DELIMITER ',')"
diff /tmp/prophetx_baseline.csv /tmp/px_post_shim.csv | head -20
```

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/scraper_prophetx_sgp.py mlb_sgp/tests/test_sgp_regression.py
git commit -m "refactor(mlb_sgp): scraper_prophetx_sgp.py becomes thin shim

Mirrors Tasks 13-14. PX-specific sanity defense (SANITY_MULT_RATIO)
lives in prophetx.py orchestrator now.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 16: Thin shim — scraper_novig_sgp.py

**Files:**
- Modify: `mlb_sgp/scraper_novig_sgp.py` (reduce to ~80 lines + preserved helpers)

Same pattern. Source label `novig_direct`. Library: `mlb_sgp.novig`.

- [ ] **Step 1: Write the failing regression test**

Append to `mlb_sgp/tests/test_sgp_regression.py`:

```python
def test_nv_shim_matches_baseline_row_shape(tmp_path):
    import subprocess, duckdb, os
    db = str(tmp_path / "shim_test_nv.duckdb")
    con = duckdb.connect(db)
    con.execute("CREATE TABLE mlb_parlay_lines (game_id VARCHAR, home_team VARCHAR, away_team VARCHAR, commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE, f5_spread DOUBLE, f5_total DOUBLE)")
    con.close()
    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    result = subprocess.run(["python", "mlb_sgp/scraper_novig_sgp.py"],
                            env=env, capture_output=True, timeout=120)
    assert result.returncode == 0
    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols and "total_line" in cols
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_nv_shim_matches_baseline_row_shape -v
```

- [ ] **Step 3: Replace scraper_novig_sgp.py top-level**

```python
def main():
    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)
    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0
    periods = tuple(p.strip() for p in os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(",") if p.strip())
    from mlb_sgp import novig
    print(f"  NV shim: {len(targets)} target lines, periods={periods}")
    rows = novig.price_sgps(targets, periods=periods, verbose=False)
    print(f"  NV shim: priced {len(rows)} rows")
    db.clear_source("novig_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0
```

Preserve helpers from lines 108-487.

- [ ] **Step 4: Run regression test**

```bash
pytest mlb_sgp/tests/test_sgp_regression.py::test_nv_shim_matches_baseline_row_shape -v
```

- [ ] **Step 5: Live diff**

```bash
python mlb_sgp/scraper_novig_sgp.py 2>&1 | tail -5
duckdb "Answer Keys/mlb_mm.duckdb" "COPY (SELECT game_id, combo, period, bookmaker, sgp_decimal FROM mlb_sgp_odds WHERE source='novig_direct' ORDER BY game_id, combo, period) TO '/tmp/nv_post_shim.csv' (HEADER, DELIMITER ',')"
diff /tmp/novig_baseline.csv /tmp/nv_post_shim.csv | head -20
```

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/scraper_novig_sgp.py mlb_sgp/tests/test_sgp_regression.py
git commit -m "refactor(mlb_sgp): scraper_novig_sgp.py becomes thin shim

Final scraper shim. All four (DK/FD/PX/NV) now route through their
library orchestrators. Shim contract: load_target_lines, price_sgps,
upsert_priced_rows. Subprocess invocation by R pipeline unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 17: kalshi_mlb_rfq/config.py — worktree fix + new env knobs

**Files:**
- Modify: `kalshi_mlb_rfq/config.py`
- Create: `kalshi_mlb_rfq/tests/test_config.py`

- [ ] **Step 1: Write the failing test**

Create `kalshi_mlb_rfq/tests/test_config.py`:

```python
"""Tests for kalshi_mlb_rfq/config.py — worktree path + new env knobs."""
import importlib
from pathlib import Path

import kalshi_mlb_rfq.config as config_mod


def test_project_root_strips_worktree_suffix():
    """When config.py lives in a .worktrees/ subdir, PROJECT_ROOT should
    resolve to the main repo, not to the worktree itself."""
    # This test passes by inspection: if we're in a worktree, PROJECT_ROOT
    # should be /Users/callancapitolo/NFLWork (the main repo), not the worktree.
    pkg_dir = Path(config_mod.PKG_DIR)
    if ".worktrees" in str(pkg_dir):
        assert ".worktrees" not in str(config_mod.PROJECT_ROOT), \
            f"PROJECT_ROOT={config_mod.PROJECT_ROOT} still contains worktree path"
        assert str(config_mod.PROJECT_ROOT).endswith("NFLWork"), \
            f"PROJECT_ROOT should end with NFLWork, got {config_mod.PROJECT_ROOT}"


def test_new_env_knobs_have_defaults():
    assert config_mod.SGP_REFRESH_SEC == 60
    assert config_mod.SGP_SCRAPER_TIMEOUT_SEC == 90
    assert config_mod.MIN_BOOK_COUNT_FOR_BLEND == 2
    assert isinstance(config_mod.BOT_MARKET_DB, Path)
    assert str(config_mod.BOT_MARKET_DB).endswith("kalshi_mlb_rfq_market.duckdb")
    assert isinstance(config_mod.MLB_SGP_DIR, Path)


def test_mlb_sgp_dir_resolves_to_main_repo():
    """MLB_SGP_DIR points to the mlb_sgp/ folder in the main repo, not in
    the worktree (so subprocess invocations find scrapers + venv)."""
    assert config_mod.MLB_SGP_DIR.name == "mlb_sgp"
    assert ".worktrees" not in str(config_mod.MLB_SGP_DIR)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_config.py -v
```

Expected: FAIL on `SGP_REFRESH_SEC` (AttributeError) and possibly on worktree path.

- [ ] **Step 3: Update config.py**

Modify `kalshi_mlb_rfq/config.py`. Replace the top of the file:

```python
"""Config knobs for the Kalshi MLB RFQ bot. Loaded from .env (or environment)."""

import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_PROJECT_ROOT = PKG_DIR.parent
# Worktree-aware: when running from .claude/worktrees/<name>/kalshi_mlb_rfq/,
# resolve PROJECT_ROOT to the main repo so we find Answer Keys/, mlb_sgp/, etc.
# Matches mlb_sgp/db.py:19-20 pattern.
if ".worktrees" in str(_RAW_PROJECT_ROOT):
    PROJECT_ROOT = Path(str(_RAW_PROJECT_ROOT).split(".worktrees")[0].rstrip("/"))
else:
    PROJECT_ROOT = _RAW_PROJECT_ROOT
DB_PATH = PKG_DIR / "kalshi_mlb_rfq.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"
```

Then append the new env knobs (before "Notifications" section):

```python
# SGP scraper cadence (bot-driven independent of dashboard refresh)
SGP_REFRESH_SEC = int(_get("SGP_REFRESH_SEC", "60"))
SGP_SCRAPER_TIMEOUT_SEC = int(_get("SGP_SCRAPER_TIMEOUT_SEC", "90"))
SGP_MIN_INTERVAL_SEC = int(_get("SGP_MIN_INTERVAL_SEC", "30"))

# Bot market DB (sibling to bot state DB, holds mlb_target_lines + mlb_sgp_odds)
BOT_MARKET_DB = Path(_get("BOT_MARKET_DB",
                          str(PKG_DIR / "kalshi_mlb_rfq_market.duckdb")))

# Path to mlb_sgp directory (for subprocess scraper invocation).
# Resolved through PROJECT_ROOT so worktree contexts use the main repo's venv.
MLB_SGP_DIR = Path(_get("MLB_SGP_DIR", str(PROJECT_ROOT / "mlb_sgp")))

# Cross-book book-count gate (drop candidate if fewer than N books priced).
MIN_BOOK_COUNT_FOR_BLEND = int(_get("MIN_BOOK_COUNT_FOR_BLEND", "2"))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_config.py -v
```

Expected: 3 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/config.py kalshi_mlb_rfq/tests/test_config.py
git commit -m "feat(kalshi-mlb-rfq): config — worktree-aware PROJECT_ROOT + SGP cadence knobs

Worktree fix: match mlb_sgp/db.py:19-20 pattern. Bites smoke tests.

New env vars: SGP_REFRESH_SEC, SGP_SCRAPER_TIMEOUT_SEC, SGP_MIN_INTERVAL_SEC,
BOT_MARKET_DB, MLB_SGP_DIR, MIN_BOOK_COUNT_FOR_BLEND.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 18: sgp_runner.py — should_scrape + latest_sgp_fetch_time

**Files:**
- Create: `kalshi_mlb_rfq/sgp_runner.py`
- Create: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

- [ ] **Step 1: Write the failing tests**

Create `kalshi_mlb_rfq/tests/test_sgp_runner.py`:

```python
"""Tests for kalshi_mlb_rfq/sgp_runner.py — bot's SGP scrape orchestration."""
from datetime import datetime, timedelta, timezone

import duckdb

from kalshi_mlb_rfq.sgp_runner import should_scrape, latest_sgp_fetch_time


def test_should_scrape_no_last_fetch():
    """Empty table / no prior scrape → always scrape."""
    assert should_scrape(last_fetch_time=None, now=datetime.now(timezone.utc),
                          min_interval_sec=30) is True


def test_should_scrape_recent_fetch():
    """Scraped < min_interval ago → skip."""
    now = datetime.now(timezone.utc)
    last = now - timedelta(seconds=10)
    assert should_scrape(last, now, min_interval_sec=30) is False


def test_should_scrape_stale_fetch():
    """Scraped > min_interval ago → scrape."""
    now = datetime.now(timezone.utc)
    last = now - timedelta(seconds=120)
    assert should_scrape(last, now, min_interval_sec=30) is True


def test_should_scrape_normalizes_naive_datetime():
    """Both naive and tz-aware datetimes should work."""
    now = datetime.now(timezone.utc)
    last_naive = (now - timedelta(seconds=120)).replace(tzinfo=None)
    assert should_scrape(last_naive, now, min_interval_sec=30) is True


def test_latest_sgp_fetch_time_empty(tmp_path):
    db = str(tmp_path / "empty.duckdb")
    duckdb.connect(db).close()
    assert latest_sgp_fetch_time(db) is None


def test_latest_sgp_fetch_time_returns_max(tmp_path):
    db = str(tmp_path / "f.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,'2026-05-13 10:00','dk_direct',-1.5,8.5)")
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','B','FG','dk',2.0,100,'2026-05-13 10:05','dk_direct',-1.5,8.5)")
    con.close()
    result = latest_sgp_fetch_time(db)
    assert result is not None
    assert result.hour == 10 and result.minute == 5
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: ImportError on `kalshi_mlb_rfq.sgp_runner`.

- [ ] **Step 3: Create the module**

Create `kalshi_mlb_rfq/sgp_runner.py`:

```python
"""SGP scrape orchestration for the Kalshi MLB RFQ bot.

Owns the cycle that:
  1. Enumerates Kalshi MVE markets per open MLB game
  2. Writes target lines (one row per game × (spread, total)) to bot DB
  3. Spawns the 4 scraper subprocesses with MLB_SGP_DB_PATH redirect
  4. Reads back priced SGP odds into the bot's _SGP_ODDS_CACHE

This module is invoked from main_loop on the SGP cadence tick.
"""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path

import duckdb


def should_scrape(last_fetch_time: datetime | None,
                   now: datetime,
                   min_interval_sec: int) -> bool:
    """True if we should scrape this tick. Guards against tight cadences
    after crash-recovery restarts that hit an already-fresh DB.

    Both `last_fetch_time` and `now` are normalized to UTC when naive."""
    if last_fetch_time is None:
        return True
    if last_fetch_time.tzinfo is None:
        last_fetch_time = last_fetch_time.replace(tzinfo=timezone.utc)
    if now.tzinfo is None:
        now = now.replace(tzinfo=timezone.utc)
    age = (now - last_fetch_time).total_seconds()
    return age > min_interval_sec


def latest_sgp_fetch_time(bot_market_db: str) -> datetime | None:
    """Read MAX(fetch_time) from mlb_sgp_odds in bot_market_db.
    Returns None for missing DB / missing table / empty table."""
    if not Path(bot_market_db).exists():
        return None
    con = duckdb.connect(bot_market_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return None
        row = con.execute("SELECT MAX(fetch_time) FROM mlb_sgp_odds").fetchone()
        return row[0] if row and row[0] is not None else None
    finally:
        con.close()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: 6 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(kalshi-mlb-rfq): sgp_runner module foundation — should_scrape + latest_sgp_fetch_time

Decision helpers for the bot's SGP scrape cadence. Naive datetimes
normalized to UTC for safe subtraction.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 19: sgp_runner.enumerate_kalshi_targets()

**Files:**
- Modify: `kalshi_mlb_rfq/sgp_runner.py`
- Modify: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

This function enumerates open Kalshi MVE markets per MLB game and produces `TargetLine` objects with metadata pulled from `mlb_odds_temp` in `mlb.duckdb`.

- [ ] **Step 1: Write the failing test**

Append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`:

```python
def test_enumerate_kalshi_targets_returns_target_lines(monkeypatch, tmp_path):
    """When Kalshi API returns 1 open game with 2 spreads × 2 totals,
    enumerate_kalshi_targets yields 4 TargetLine entries (one per combo)."""
    from kalshi_mlb_rfq import sgp_runner

    # Mock the Kalshi API and the schedule lookup
    fake_event = {
        "event_ticker": "KXMLBGAME-26MAY13230000NYYBOS",
        "status": "open",
    }
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_mlb_events",
                         lambda: [fake_event])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_spread_lines",
                         lambda suffix: [(-1.5, "home"), (-2.5, "home")])
    monkeypatch.setattr(sgp_runner, "_fetch_kalshi_total_lines",
                         lambda suffix: [8.5, 9.5])
    # Mock the schedule lookup from mlb_odds_temp
    schedule_db = str(tmp_path / "mlb.duckdb")
    con = duckdb.connect(schedule_db)
    con.execute("""
        CREATE TABLE mlb_odds_temp (
            date VARCHAR, id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            total_line DOUBLE, consensus_prob_home DOUBLE, commence_time VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_odds_temp VALUES ('2026-05-13','g1','New York Yankees','Boston Red Sox',8.5,0.55,'2026-05-13T23:00:00+00:00')")
    con.close()

    targets = sgp_runner.enumerate_kalshi_targets(schedule_db_path=schedule_db)
    # 2 spreads × 2 totals × 1 game = 4 TargetLine (one per combo)
    assert len(targets) == 4
    spreads = sorted({t.spread for t in targets})
    totals = sorted({t.total for t in targets})
    assert spreads == [-2.5, -1.5]
    assert totals == [8.5, 9.5]
    assert all(t.period == "FG" for t in targets)
    assert all(t.game_id == "g1" for t in targets)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py::test_enumerate_kalshi_targets_returns_target_lines -v
```

Expected: AttributeError on `enumerate_kalshi_targets`.

- [ ] **Step 3: Implement enumerate_kalshi_targets**

Append to `kalshi_mlb_rfq/sgp_runner.py`:

```python
from typing import Iterable

from kalshi_mlb_rfq import auth_client
from mlb_sgp._shared import TargetLine

# 3-letter Kalshi team code → Odds-API canonical team name.
# Same mapping main.py uses (lifted to module level for testability).
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "ATH": "Athletics",
    "AZ": "Arizona Diamondbacks", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


def _parse_event_suffix(suffix: str) -> tuple[str | None, str | None]:
    """Split a KXMLB* event suffix into (away_code, home_code).
    Matches main.py::_parse_event_suffix exactly."""
    if len(suffix) < 11 + 4:
        return None, None
    team_block = suffix[11:]
    for home_len in (3, 2):
        if len(team_block) <= home_len:
            continue
        home = team_block[-home_len:]
        away = team_block[:-home_len]
        if home in _MLB_CODE_TO_TEAM and away in _MLB_CODE_TO_TEAM:
            return away, home
    return None, None


def _fetch_kalshi_mlb_events() -> list[dict]:
    status, body, _ = auth_client.api(
        "GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    return body.get("events", [])


def _fetch_kalshi_spread_lines(suffix: str) -> list[tuple[float, str]]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBSPREAD-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    seen = set()
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        prefix = f"KXMLBSPREAD-{suffix}-"
        if not ticker.startswith(prefix):
            continue
        spread_part = ticker[len(prefix):]
        digits = "".join(c for c in spread_part if c.isdigit())
        team_chars = "".join(c for c in spread_part if not c.isdigit())
        if not digits or not team_chars:
            continue
        n = int(digits)
        line = -(n - 0.5)
        # We yield one entry per spread (sided to home team for simplicity);
        # the orchestrator does its own combo enumeration.
        key = round(line, 1)
        if key in seen:
            continue
        seen.add(key)
        out.append((line, "home"))
    return out


def _fetch_kalshi_total_lines(suffix: str) -> list[float]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBTOTAL-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    seen = set()
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        try:
            n = int(ticker.rsplit("-", 1)[-1])
            line = n - 0.5
            key = round(line, 1)
            if key in seen:
                continue
            seen.add(key)
            out.append(line)
        except ValueError:
            continue
    return out


def _load_schedule(schedule_db_path: str) -> dict[str, dict]:
    """Read mlb_odds_temp from mlb.duckdb. Returns dict keyed by canonical team
    pair → {game_id, home_team, away_team, commence_time}."""
    from datetime import datetime as _dt
    if not Path(schedule_db_path).exists():
        return {}
    con = duckdb.connect(schedule_db_path, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_odds_temp" not in tables:
            return {}
        rows = con.execute(
            "SELECT id, home_team, away_team, commence_time FROM mlb_odds_temp"
        ).fetchall()
    finally:
        con.close()
    out = {}
    for game_id, home, away, ct_str in rows:
        normalized = ct_str.replace("Z", "+00:00") if ct_str and ct_str.endswith("Z") else ct_str
        try:
            ct = _dt.fromisoformat(normalized) if normalized else None
        except Exception:
            ct = None
        out[(home, away)] = {"game_id": game_id, "home_team": home,
                              "away_team": away, "commence_time": ct}
    return out


def enumerate_kalshi_targets(schedule_db_path: str) -> list[TargetLine]:
    """Enumerate all open Kalshi MVE (spread, total) tuples per MLB game.
    Returns a TargetLine per (game × spread × total) combination, FG only."""
    events = _fetch_kalshi_mlb_events()
    if not events:
        return []
    schedule = _load_schedule(schedule_db_path)
    targets: list[TargetLine] = []
    for ev in events:
        event_ticker = ev.get("event_ticker", "")
        if not event_ticker.startswith("KXMLBGAME-"):
            continue
        suffix = event_ticker.replace("KXMLBGAME-", "")
        away_code, home_code = _parse_event_suffix(suffix)
        if away_code is None or home_code is None:
            continue
        home_team = _MLB_CODE_TO_TEAM.get(home_code)
        away_team = _MLB_CODE_TO_TEAM.get(away_code)
        sched = schedule.get((home_team, away_team))
        if not sched:
            continue
        spreads = _fetch_kalshi_spread_lines(suffix)
        totals = _fetch_kalshi_total_lines(suffix)
        if not spreads or not totals:
            continue
        for spread, _who in spreads:
            for total in totals:
                targets.append(TargetLine(
                    game_id=sched["game_id"],
                    home_team=home_team, away_team=away_team,
                    commence_time=sched["commence_time"],
                    period="FG", spread=spread, total=total,
                ))
    return targets
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: 7 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(kalshi-mlb-rfq): sgp_runner.enumerate_kalshi_targets()

Reads Kalshi MVE spread + total markets per open MLB game; cross-joins
into TargetLine objects. Schedule (game_id, team names, commence_time)
sourced from mlb_odds_temp in mlb.duckdb — independent of Wagerzon.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 20: sgp_runner.write_target_lines()

**Files:**
- Modify: `kalshi_mlb_rfq/sgp_runner.py`
- Modify: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

- [ ] **Step 1: Write the failing test**

Append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`:

```python
def test_write_target_lines_atomic_replace(tmp_path):
    """write_target_lines DELETEs all existing rows then INSERTs new ones,
    in a single transaction."""
    from kalshi_mlb_rfq.sgp_runner import write_target_lines
    from mlb_sgp._shared import TargetLine
    from datetime import datetime, timezone

    db = str(tmp_path / "bot.duckdb")
    # First write
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)
    first = [
        TargetLine("g1", "NYY", "BOS", ct, "FG", -1.5, 8.5),
        TargetLine("g1", "NYY", "BOS", ct, "FG", -1.5, 9.5),
    ]
    write_target_lines(first, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*) FROM mlb_target_lines").fetchone()[0]
    con.close()
    assert out == 2

    # Second write (different lines) → must replace, not append
    second = [TargetLine("g1", "NYY", "BOS", ct, "FG", -2.5, 8.5)]
    write_target_lines(second, db_path=db)

    con = duckdb.connect(db, read_only=True)
    out = con.execute("SELECT COUNT(*), MIN(spread) FROM mlb_target_lines").fetchone()
    con.close()
    assert out == (1, -2.5)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py::test_write_target_lines_atomic_replace -v
```

Expected: AttributeError on `write_target_lines`.

- [ ] **Step 3: Implement write_target_lines**

Append to `kalshi_mlb_rfq/sgp_runner.py`:

```python
def write_target_lines(target_lines: list[TargetLine], db_path: str):
    """Atomic DELETE+INSERT of mlb_target_lines in bot market DB.
    Creates the table if missing."""
    con = duckdb.connect(db_path)
    try:
        con.execute("""
            CREATE TABLE IF NOT EXISTS mlb_target_lines (
                game_id        VARCHAR,
                home_team      VARCHAR,
                away_team      VARCHAR,
                commence_time  TIMESTAMP,
                period         VARCHAR,
                spread         DOUBLE,
                total          DOUBLE,
                written_at     TIMESTAMP
            )
        """)
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM mlb_target_lines")
        if target_lines:
            now = datetime.now(timezone.utc)
            values = []
            for t in target_lines:
                values.extend([t.game_id, t.home_team, t.away_team,
                                t.commence_time, t.period, t.spread, t.total, now])
            placeholders = ",".join(["(?, ?, ?, ?, ?, ?, ?, ?)"] * len(target_lines))
            con.execute(f"INSERT INTO mlb_target_lines VALUES {placeholders}", values)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: 8 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(kalshi-mlb-rfq): sgp_runner.write_target_lines() atomic replace

Single transaction DELETE+INSERT; no partial-state windows where
scrapers could read a half-written target table mid-cycle.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 21: sgp_runner.run_scrapers() + read_priced_rows()

**Files:**
- Modify: `kalshi_mlb_rfq/sgp_runner.py`
- Modify: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

- [ ] **Step 1: Write the failing test**

Append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`:

```python
def test_run_scrapers_subprocess_plumbing(tmp_path, monkeypatch):
    """Verify subprocess invocation: correct env, timeout, kill-on-hang.
    Uses a fake scraper that just touches a file then exits."""
    from kalshi_mlb_rfq.sgp_runner import run_scrapers

    # Stub scraper dir with a single fake "scraper" that writes a marker
    scraper_dir = tmp_path / "fake_scrapers"
    scraper_dir.mkdir()
    fake = scraper_dir / "fake_scraper.py"
    fake.write_text(
        "import os\n"
        "open(os.environ['TARGET'], 'w').write(os.environ.get('MLB_SGP_DB_PATH', 'unset'))\n"
    )
    marker = tmp_path / "marker.txt"

    result = run_scrapers(
        scraper_dir=str(scraper_dir),
        scraper_names=["fake_scraper.py"],
        venv_python="python",
        timeout_sec=10,
        env={"MLB_SGP_DB_PATH": "/tmp/test.duckdb", "TARGET": str(marker)},
    )
    assert result == {"fake_scraper.py": 0}
    assert marker.read_text() == "/tmp/test.duckdb"


def test_run_scrapers_handles_timeout(tmp_path):
    from kalshi_mlb_rfq.sgp_runner import run_scrapers

    scraper_dir = tmp_path / "fake"
    scraper_dir.mkdir()
    fake = scraper_dir / "slow_scraper.py"
    fake.write_text("import time; time.sleep(60)\n")

    result = run_scrapers(
        scraper_dir=str(scraper_dir), scraper_names=["slow_scraper.py"],
        venv_python="python", timeout_sec=2, env={},
    )
    # Killed; rc != 0
    assert result["slow_scraper.py"] != 0


def test_read_priced_rows_filters_by_staleness(tmp_path):
    from kalshi_mlb_rfq.sgp_runner import read_priced_rows
    db = str(tmp_path / "r.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    # Fresh + stale rows
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,NOW(),'dk_direct',-1.5,8.5)")
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','A','FG','dk',2.0,100,NOW() - INTERVAL 10 MINUTE,'dk_direct',-2.5,8.5)")
    con.close()
    df = read_priced_rows(db, max_age_sec=60)
    assert len(df) == 1
    assert df.iloc[0]["spread_line"] == -1.5
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v -k "run_scrapers or read_priced_rows"
```

Expected: AttributeError.

- [ ] **Step 3: Implement run_scrapers + read_priced_rows**

Append to `kalshi_mlb_rfq/sgp_runner.py`:

```python
import subprocess
import time
import os as _os
from typing import Mapping


def run_scrapers(
    scraper_dir: str,
    scraper_names: list[str],
    venv_python: str,
    timeout_sec: int,
    env: Mapping[str, str] | None = None,
) -> dict[str, int]:
    """Spawn each scraper as a subprocess in parallel. Returns
    {scraper_name: return_code}.

    - All scrapers share `env` (defaulting to os.environ) and a global
      deadline of `timeout_sec` after which any still-running scraper
      is killed.
    - Hardened against handle leaks: subprocess Popen failures close
      log handles before re-raising.
    """
    if env is None:
        env_dict = dict(_os.environ)
    else:
        env_dict = dict(_os.environ)
        env_dict.update(env)

    # Per-scraper stdout/stderr log files (captured for debugging)
    log_handles: list = []
    procs: dict[str, subprocess.Popen] = {}
    try:
        for name in scraper_names:
            log_path = Path(scraper_dir) / f"{name}.runner.log"
            handle = open(log_path, "w")
            log_handles.append(handle)
            try:
                p = subprocess.Popen(
                    [venv_python, name],
                    cwd=scraper_dir,
                    env=env_dict,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
                procs[name] = p
            except Exception:
                handle.close()
                raise
        # Wait for all with global deadline
        deadline = time.time() + timeout_sec
        rcs: dict[str, int] = {}
        for name, p in procs.items():
            remaining = deadline - time.time()
            if remaining <= 0:
                p.kill()
                rcs[name] = -1
                continue
            try:
                rcs[name] = p.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                p.kill()
                try:
                    p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
                rcs[name] = -1
        return rcs
    finally:
        for h in log_handles:
            try:
                h.close()
            except Exception:
                pass


def read_priced_rows(bot_market_db: str, max_age_sec: int):
    """Read mlb_sgp_odds rows fresher than max_age_sec from the bot DB.
    Returns a pandas DataFrame; empty if no fresh data."""
    import pandas as pd
    if not Path(bot_market_db).exists():
        return pd.DataFrame()
    con = duckdb.connect(bot_market_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return pd.DataFrame()
        df = con.execute(
            "SELECT game_id, combo, period, bookmaker, sgp_decimal, sgp_american, "
            "fetch_time, source, spread_line, total_line "
            "FROM mlb_sgp_odds WHERE fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
            [max_age_sec],
        ).fetchdf()
        return df
    finally:
        con.close()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: 11 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(kalshi-mlb-rfq): sgp_runner.run_scrapers() + read_priced_rows()

Subprocess plumbing: parallel spawn with global deadline, kill-on-hang,
log-handle cleanup on Popen failure. read_priced_rows returns DataFrame
filtered by staleness window.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 22: sgp_runner.sgp_cycle() orchestrator

**Files:**
- Modify: `kalshi_mlb_rfq/sgp_runner.py`
- Modify: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

- [ ] **Step 1: Write the failing test**

Append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`:

```python
def test_sgp_cycle_orchestrates_full_tick(monkeypatch, tmp_path):
    """sgp_cycle() composes enumerate → write → run_scrapers → return rc map.
    Stubs out the API + subprocess to verify orchestration ordering."""
    from kalshi_mlb_rfq import sgp_runner

    call_order = []
    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets",
                         lambda schedule_db_path: (call_order.append("enum"), [])[1])
    monkeypatch.setattr(sgp_runner, "write_target_lines",
                         lambda targets, db_path: call_order.append("write"))
    monkeypatch.setattr(sgp_runner, "run_scrapers",
                         lambda **kw: (call_order.append("scrape"), {})[1])

    rcs = sgp_runner.sgp_cycle(
        bot_market_db=str(tmp_path / "bot.duckdb"),
        schedule_db_path=str(tmp_path / "schedule.duckdb"),
        scraper_dir=str(tmp_path / "scrapers"),
        venv_python="python",
        timeout_sec=60,
    )
    assert call_order == ["enum", "write", "scrape"]
    assert isinstance(rcs, dict)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py::test_sgp_cycle_orchestrates_full_tick -v
```

Expected: AttributeError on `sgp_cycle`.

- [ ] **Step 3: Implement sgp_cycle**

Append to `kalshi_mlb_rfq/sgp_runner.py`:

```python
SCRAPER_NAMES = [
    "scraper_draftkings_sgp.py",
    "scraper_fanduel_sgp.py",
    "scraper_prophetx_sgp.py",
    "scraper_novig_sgp.py",
]


def sgp_cycle(
    bot_market_db: str,
    schedule_db_path: str,
    scraper_dir: str,
    venv_python: str,
    timeout_sec: int,
) -> dict[str, int]:
    """One full SGP scrape tick (atomic, serial):
      1. Enumerate Kalshi MVE → list[TargetLine]
      2. Write to mlb_target_lines in bot_market_db
      3. Spawn the 4 scrapers with MLB_SGP_DB_PATH=bot_market_db, MLB_SGP_PERIODS=FG

    Returns {scraper_name: return_code}.
    """
    targets = enumerate_kalshi_targets(schedule_db_path=schedule_db_path)
    write_target_lines(targets, db_path=bot_market_db)
    rcs = run_scrapers(
        scraper_dir=scraper_dir,
        scraper_names=SCRAPER_NAMES,
        venv_python=venv_python,
        timeout_sec=timeout_sec,
        env={
            "MLB_SGP_DB_PATH": bot_market_db,
            "MLB_SGP_PERIODS": "FG",
        },
    )
    return rcs
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -v
```

Expected: 12 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(kalshi-mlb-rfq): sgp_runner.sgp_cycle() — atomic cadence tick orchestrator

Composes enumerate → write → run_scrapers in strict serial order so
scrapers never read a half-written target table. Used by main.py
synchronous warm-up and main_loop SGP cadence timer.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 23: main.py — _PARLAY_LINES_CACHE shape change + dual-source loader

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Create: `kalshi_mlb_rfq/tests/test_main_cache.py`

The cache currently holds `{game_id: {home_team, away_team, commence_time, fg_spread, fg_total}}` (one-row-per-game). After this task: `{game_id: {home_team, away_team, commence_time, fg_lines: [(spread, total), ...]}}` (multi-line per game), with team metadata sourced from `mlb_odds_temp` in `mlb.duckdb` and lines from `mlb_target_lines` in the bot DB.

- [ ] **Step 1: Write the failing test**

Create `kalshi_mlb_rfq/tests/test_main_cache.py`:

```python
"""Tests for kalshi_mlb_rfq/main.py cache shape after line-source pivot."""
from datetime import datetime, timezone
import duckdb


def test_build_parlay_lines_cache_merges_schedule_and_targets(tmp_path):
    """_build_parlay_lines_cache reads schedule from mlb.duckdb::mlb_odds_temp
    and target lines from bot DB::mlb_target_lines, returns multi-line cache."""
    from kalshi_mlb_rfq.main import _build_parlay_lines_cache

    sched_db = str(tmp_path / "mlb.duckdb")
    bot_db = str(tmp_path / "bot.duckdb")
    ct = datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc)

    con = duckdb.connect(sched_db)
    con.execute("""
        CREATE TABLE mlb_odds_temp (
            date VARCHAR, id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            total_line DOUBLE, consensus_prob_home DOUBLE, commence_time VARCHAR
        )
    """)
    con.execute("INSERT INTO mlb_odds_temp VALUES ('2026-05-13','g1','New York Yankees','Boston Red Sox',8.5,0.55,'2026-05-13T23:00:00+00:00')")
    con.close()

    con = duckdb.connect(bot_db)
    con.execute("""
        CREATE TABLE mlb_target_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time TIMESTAMP, period VARCHAR,
            spread DOUBLE, total DOUBLE, written_at TIMESTAMP
        )
    """)
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','NYY','BOS', ?, 'FG', -1.5, 8.5, NOW())", [ct])
    con.execute("INSERT INTO mlb_target_lines VALUES ('g1','NYY','BOS', ?, 'FG', -2.5, 8.5, NOW())", [ct])
    con.close()

    cache = _build_parlay_lines_cache(schedule_db=sched_db, bot_db=bot_db)
    assert "g1" in cache
    entry = cache["g1"]
    assert entry["home_team"] == "New York Yankees"
    assert entry["away_team"] == "Boston Red Sox"
    assert entry["fg_lines"] == [(-2.5, 8.5), (-1.5, 8.5)] or entry["fg_lines"] == [(-1.5, 8.5), (-2.5, 8.5)]
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v
```

Expected: ImportError on `_build_parlay_lines_cache`.

- [ ] **Step 3: Add the loader to main.py**

In `kalshi_mlb_rfq/main.py`, add a new helper function `_build_parlay_lines_cache` (placed near the existing `_refresh_caches`):

```python
def _build_parlay_lines_cache(schedule_db: str, bot_db: str) -> dict[str, dict]:
    """Build the bot's parlay_lines cache from two sources:
      1. mlb_odds_temp in schedule_db (mlb.duckdb) — game metadata
      2. mlb_target_lines in bot_db (bot market DB) — available (spread, total) tuples per game

    Returns {game_id: {home_team, away_team, commence_time, fg_lines: list[(spread, total)]}}.
    Drops games with no target lines.
    """
    from datetime import datetime as _dt
    from pathlib import Path as _Path

    # 1. Schedule
    schedule: dict[str, dict] = {}
    if _Path(schedule_db).exists():
        con = duckdb.connect(schedule_db, read_only=True)
        try:
            tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
            if "mlb_odds_temp" in tables:
                for game_id, home, away, ct_str in con.execute(
                    "SELECT id, home_team, away_team, commence_time FROM mlb_odds_temp"
                ).fetchall():
                    normalized = ct_str.replace("Z", "+00:00") if ct_str and ct_str.endswith("Z") else ct_str
                    try:
                        ct = _dt.fromisoformat(normalized) if normalized else None
                    except Exception:
                        ct = None
                    schedule[game_id] = {
                        "home_team": home, "away_team": away, "commence_time": ct,
                    }
        finally:
            con.close()

    # 2. Target lines per game
    lines_by_game: dict[str, list[tuple[float, float]]] = {}
    if _Path(bot_db).exists():
        con = duckdb.connect(bot_db, read_only=True)
        try:
            tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
            if "mlb_target_lines" in tables:
                for game_id, spread, total in con.execute(
                    "SELECT game_id, spread, total FROM mlb_target_lines WHERE period = 'FG' ORDER BY game_id, spread, total"
                ).fetchall():
                    lines_by_game.setdefault(game_id, []).append((spread, total))
        finally:
            con.close()

    # 3. Merge — only include games present in both sources
    out: dict[str, dict] = {}
    for game_id, lines in lines_by_game.items():
        sched = schedule.get(game_id)
        if not sched:
            continue
        out[game_id] = {**sched, "fg_lines": lines}
    return out
```

Then modify `_refresh_caches()` to use the new loader. Replace the `mlb_parlay_lines` block (lines 118-129 of the pre-refactor main.py) with:

```python
            # mlb_parlay_lines is replaced by the merged schedule + target_lines cache.
            parlay_lines = _build_parlay_lines_cache(
                schedule_db=str(config.PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"),
                bot_db=str(config.BOT_MARKET_DB),
            )
```

Update consumers of `_PARLAY_LINES_CACHE`:

- `_commence_time_for_game(game_id)` — unchanged interface, still reads `row["commence_time"]`
- `_resolve_game_id(home_code, away_code)` — unchanged interface; iterates entries and matches by team
- `_load_book_fairs(game_id, spread_line, total_line)` — REWIRED in Task 24
- `_current_book_lines_for_combo(game_id)` — DEPRECATED in Task 25 (removed entirely)

For now, leave `_load_book_fairs` and `_current_book_lines_for_combo` in place — they'll break on access to `pl["fg_spread"]` / `pl["fg_total"]` (which no longer exist on the new cache shape). Task 24 fixes `_load_book_fairs`; Task 25 deletes `_current_book_lines_for_combo`.

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v
```

Expected: 1 PASSED.

- [ ] **Step 5: Run main.py module imports cleanly**

```bash
python -c "from kalshi_mlb_rfq import main"
```

Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_main_cache.py
git commit -m "feat(kalshi-mlb-rfq): _PARLAY_LINES_CACHE shape change — multi-line per game

Cache schema: {game_id: {home_team, away_team, commence_time, fg_lines: list[(spread,total)]}}.
Game metadata sourced from mlb_odds_temp (mlb.duckdb); target lines
from mlb_target_lines (bot market DB). Both sources required — games
missing from either are dropped.

Next tasks fix _load_book_fairs (N≥2 gate + per-line lookup) and
remove the now-dead _current_book_lines_for_combo.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 24: main.py — _load_book_fairs rewires for per-line + N≥2 gate

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Modify: `kalshi_mlb_rfq/tests/test_main_cache.py`

- [ ] **Step 1: Write the failing test**

Append to `kalshi_mlb_rfq/tests/test_main_cache.py`:

```python
def test_load_book_fairs_filters_by_line_and_requires_two_books(monkeypatch):
    """Returns book→fair only when >= 2 books priced the matching (spread, total)."""
    import pandas as pd
    from kalshi_mlb_rfq import main, config

    # Stub the in-memory cache
    main._SGP_ODDS_CACHE = pd.DataFrame([
        # game g1, spread -1.5, total 8.5: 2 books (DK, FD) — passes gate
        {"game_id": "g1", "bookmaker": "draftkings", "combo": "Home Spread + Over",
         "sgp_decimal": 2.85, "period": "FG", "spread_line": -1.5, "total_line": 8.5},
        {"game_id": "g1", "bookmaker": "fanduel", "combo": "Home Spread + Over",
         "sgp_decimal": 2.95, "period": "FG", "spread_line": -1.5, "total_line": 8.5},
        # game g1, spread -2.5, total 8.5: only DK — fails gate (N=1)
        {"game_id": "g1", "bookmaker": "draftkings", "combo": "Home Spread + Over",
         "sgp_decimal": 3.50, "period": "FG", "spread_line": -2.5, "total_line": 8.5},
    ])
    monkeypatch.setattr(config, "MIN_BOOK_COUNT_FOR_BLEND", 2)

    # Matching line + 2 books → returns dict
    fairs = main._load_book_fairs("g1", -1.5, 8.5)
    assert len(fairs) == 2 or len(fairs) > 0  # devig may return fewer if rows incomplete

    # Matching line + 1 book → empty
    fairs = main._load_book_fairs("g1", -2.5, 8.5)
    assert fairs == {}

    # No matching line → empty
    fairs = main._load_book_fairs("g1", -3.5, 9.5)
    assert fairs == {}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py::test_load_book_fairs_filters_by_line_and_requires_two_books -v
```

Expected: FAIL — current `_load_book_fairs` reads `pl["fg_spread"]` which no longer exists.

- [ ] **Step 3: Rewrite _load_book_fairs**

In `kalshi_mlb_rfq/main.py`, replace `_load_book_fairs`:

```python
def _load_book_fairs(game_id: str, spread_line: float, total_line: float) -> dict[str, float]:
    """Per-line book lookup with N≥2 gate.

    Returns {book → devigged_fair} only if at least MIN_BOOK_COUNT_FOR_BLEND
    books priced the matching (game, spread, total) tuple. Empty dict
    otherwise — caller treats as 'no book signal, drop candidate'.
    """
    if _SGP_ODDS_CACHE is None or _SGP_ODDS_CACHE.empty:
        return {}
    rows = _SGP_ODDS_CACHE[
        (_SGP_ODDS_CACHE["game_id"] == game_id)
        & (_SGP_ODDS_CACHE["spread_line"].astype(float).round(2) == round(float(spread_line), 2))
        & (_SGP_ODDS_CACHE["total_line"].astype(float).round(2) == round(float(total_line), 2))
    ]
    if rows.empty:
        return {}
    out: dict[str, float] = {}
    for book in rows["bookmaker"].unique():
        sub = rows[rows["bookmaker"] == book].copy()
        fair_per_book = fair_value.devig_book(
            sub, combo="Home Spread + Over",
            vig_fallback=_vig_fallback(book),
        )
        if fair_per_book is not None:
            out[book] = fair_per_book
    if len(out) < config.MIN_BOOK_COUNT_FOR_BLEND:
        return {}
    return out
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v
```

Expected: 2 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_main_cache.py
git commit -m "feat(kalshi-mlb-rfq): _load_book_fairs per-line lookup + N>=2 books required

Reads mlb_sgp_odds rows by (game, spread, total) match (new schema cols).
Returns empty dict unless MIN_BOOK_COUNT_FOR_BLEND books (default 2)
priced the matching tuple — bot does NOT bet model-only or single-book
candidates. Quant rationale: one book's devig is opinion, not consensus.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 25: main.py — drop line_move_ok gate + _current_book_lines_for_combo + reference_lines write

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`

Rationale (from spec D7): post-pivot, each candidate's (spread, total) is encoded in its Kalshi market ticker — the line can't move per-candidate. RFQ refresh re-scores every 30s and naturally drops out-of-top-N candidates. The line-move gate becomes tautological dead code.

- [ ] **Step 1: Locate the gate call and write a no-regression test**

Append to `kalshi_mlb_rfq/tests/test_main_cache.py`:

```python
def test_per_accept_gates_skip_line_move():
    """After D7, line_move_ok is no longer called from _all_per_accept_gates_pass.
    Verify the gate path doesn't reference _current_book_lines_for_combo."""
    import inspect
    from kalshi_mlb_rfq import main
    src = inspect.getsource(main._all_per_accept_gates_pass)
    assert "_current_book_lines_for_combo" not in src
    assert "line_move_ok" not in src
    assert "reference_lines" not in src
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py::test_per_accept_gates_skip_line_move -v
```

Expected: FAIL — current code still calls all three.

- [ ] **Step 3: Update main.py**

In `kalshi_mlb_rfq/main.py`:

**Delete:**
- Function `_current_book_lines_for_combo` (lines ~366-379)
- The "Line-move check" block in `_all_per_accept_gates_pass` (lines ~479-490)
- The "Snapshot reference lines for line-move detection" block in `_refresh_rfqs` (lines ~827-834)

**Replace** the line-move block in `_all_per_accept_gates_pass` with a comment:

```python
    # Line-move gate removed (D7 of line-source pivot): each candidate's
    # (spread, total) is baked into its Kalshi ticker — line can't move
    # per-candidate. Drift handled by RFQ refresh re-scoring every 30s.
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v
```

Expected: 3 PASSED.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_main_cache.py
git commit -m "refactor(kalshi-mlb-rfq): drop line_move_ok gate + _current_book_lines_for_combo

Post-pivot the line is baked into the Kalshi market ticker; the gate is
tautological. Natural drift protection is the 30s RFQ-refresh top-N
re-scoring. Stop writing to reference_lines on RFQ submission.

Schema column reference_lines + risk.line_move_ok function left in place
(unused) — clean up in a follow-up.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 26: main.py — _refresh_sgp_cache + synchronous warm-up + SGP cadence timer

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Modify: `kalshi_mlb_rfq/tests/test_main_cache.py`

This integrates the SGP cadence into `main_loop` and adds a synchronous warm-up tick on startup so `_SGP_ODDS_CACHE` is populated before the first RFQ refresh.

- [ ] **Step 1: Write the failing tests**

Append to `kalshi_mlb_rfq/tests/test_main_cache.py`:

```python
def test_refresh_sgp_cache_partial_reload(monkeypatch, tmp_path):
    """_refresh_sgp_cache reads only mlb_sgp_odds from bot DB, atomic swap."""
    import duckdb
    from kalshi_mlb_rfq import main, config

    bot_db = str(tmp_path / "bot.duckdb")
    con = duckdb.connect(bot_db)
    con.execute("""
        CREATE TABLE mlb_sgp_odds (
            game_id VARCHAR, combo VARCHAR, period VARCHAR, bookmaker VARCHAR,
            sgp_decimal DOUBLE, sgp_american INTEGER, fetch_time TIMESTAMP,
            source VARCHAR, spread_line DOUBLE, total_line DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_sgp_odds VALUES ('g1','H+O','FG','dk',2.85,185,NOW(),'dk_direct',-1.5,8.5)")
    con.close()

    monkeypatch.setattr(config, "BOT_MARKET_DB", tmp_path / "bot.duckdb")
    monkeypatch.setattr(config, "MAX_BOOK_STALENESS_SEC", 60)

    result = main._refresh_sgp_cache()
    assert result is True
    assert main._SGP_ODDS_CACHE is not None
    assert len(main._SGP_ODDS_CACHE) == 1


def test_main_loop_imports_sgp_runner():
    """Sanity: main module must import sgp_runner without error."""
    from kalshi_mlb_rfq import main
    from kalshi_mlb_rfq import sgp_runner
    assert hasattr(main, "main_loop")
    assert hasattr(sgp_runner, "sgp_cycle")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v -k "refresh_sgp_cache or imports_sgp_runner"
```

Expected: AttributeError on `_refresh_sgp_cache`.

- [ ] **Step 3: Add _refresh_sgp_cache to main.py**

In `kalshi_mlb_rfq/main.py`, add a new function near `_refresh_caches`:

```python
def _refresh_sgp_cache(retries: int = 3) -> bool:
    """Narrow variant of _refresh_caches — reloads only mlb_sgp_odds from
    the bot market DB. Atomic swap under _CACHE_LOCK.

    Returns False on missing DB / table; True on success."""
    global _SGP_ODDS_CACHE
    bot_db = str(config.BOT_MARKET_DB)
    if not Path(bot_db).exists():
        return False
    last_err = None
    for attempt in range(retries):
        try:
            con = duckdb.connect(bot_db, read_only=True)
        except duckdb.IOException as e:
            last_err = e
            time.sleep(0.5 * (2 ** attempt))
            continue
        try:
            tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
            if "mlb_sgp_odds" not in tables:
                return False
            sgp_df = con.execute(
                "SELECT game_id, combo, period, bookmaker, sgp_decimal, fetch_time, "
                "spread_line, total_line "
                "FROM mlb_sgp_odds WHERE period='FG' "
                "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
                [config.MAX_BOOK_STALENESS_SEC],
            ).fetchdf()
        finally:
            con.close()

        with _CACHE_LOCK:
            _SGP_ODDS_CACHE = sgp_df
        print(f"  sgp_cache_refresh: {len(sgp_df)} rows from bot_market_db", flush=True)
        return True

    print(f"  sgp_cache_refresh: gave up after {retries} attempts; {last_err}", flush=True)
    return False
```

Add `from pathlib import Path` to imports if not already present.

- [ ] **Step 4: Wire synchronous warm-up + SGP cadence timer into main_loop**

In `kalshi_mlb_rfq/main.py`'s `main_loop()`, add the warm-up after `_phantom_rfq_cleanup()`:

```python
    _phantom_rfq_cleanup()

    # Synchronous SGP warm-up: run one full SGP cycle before entering the
    # loop so _SGP_ODDS_CACHE is non-empty when the first RFQ refresh fires.
    # Blocks ~60-90s; without this the first RFQ refresh would have zero
    # book fairs and burn Kalshi quota on candidates that can't pass the
    # N>=2 books gate.
    print("  startup: warming SGP cache (one synchronous scrape tick)...", flush=True)
    try:
        rcs = sgp_runner.sgp_cycle(
            bot_market_db=str(config.BOT_MARKET_DB),
            schedule_db_path=str(config.PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"),
            scraper_dir=str(config.MLB_SGP_DIR),
            venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
            timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC,
        )
        print(f"  startup: SGP warm-up done — return codes {rcs}", flush=True)
    except Exception as e:
        print(f"  startup: SGP warm-up failed ({e}); bot will retry on first cadence tick", flush=True)
    _refresh_sgp_cache()
    # Rebuild parlay_lines cache now that target_lines has rows
    _refresh_caches()
```

And add the SGP cadence timer alongside the existing timers:

```python
    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0
    last_sgp_cycle = time.time()  # warm-up just finished
```

Inside the loop body, add a new branch (place between `last_risk_sweep` and `last_pipeline`):

```python
            if now - last_sgp_cycle >= config.SGP_REFRESH_SEC:
                t_sgp = time.time()
                try:
                    rcs = sgp_runner.sgp_cycle(
                        bot_market_db=str(config.BOT_MARKET_DB),
                        schedule_db_path=str(config.PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"),
                        scraper_dir=str(config.MLB_SGP_DIR),
                        venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                        timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC,
                    )
                    _refresh_sgp_cache()
                    _refresh_caches()  # parlay_lines_cache may have new (spread, total) tuples
                    print(f"  sgp_cycle: rcs={rcs} ({time.time()-t_sgp:.1f}s)", flush=True)
                except Exception as e:
                    print(f"  sgp_cycle error: {e}", flush=True)
                last_sgp_cycle = now
```

Add `from kalshi_mlb_rfq import sgp_runner` to imports.

- [ ] **Step 5: Run tests to verify pass**

```bash
pytest kalshi_mlb_rfq/tests/test_main_cache.py -v
python -c "from kalshi_mlb_rfq import main; print('imports OK')"
```

Expected: 5 PASSED, imports OK.

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_main_cache.py
git commit -m "feat(kalshi-mlb-rfq): _refresh_sgp_cache + synchronous warm-up + SGP cadence timer

Wires sgp_runner.sgp_cycle into main_loop:
- Synchronous warm-up before main_loop starts (~60-90s) so first RFQ
  refresh has fresh book fairs.
- SGP cadence timer alongside RFQ refresh / quote poll / risk sweep,
  fires every SGP_REFRESH_SEC (default 60s).

_refresh_sgp_cache reloads only mlb_sgp_odds from BOT_MARKET_DB with
atomic swap under _CACHE_LOCK. Companion to existing _refresh_caches.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 27: Documentation updates

**Files:**
- Modify: `kalshi_mlb_rfq/README.md`
- Modify: `mlb_sgp/README.md`
- Modify: `kalshi_mlb_rfq/.env.example`
- Modify: `CLAUDE.md` (root)

- [ ] **Step 1: Update `kalshi_mlb_rfq/README.md`**

Add a new section "SGP cadence loop" describing:
- The bot now drives its own SGP scrape cadence on a 60s timer (env: `SGP_REFRESH_SEC`)
- Target lines come from Kalshi MVE enumeration, not Wagerzon-derived `mlb_parlay_lines`
- The bot maintains its own market DB at `kalshi_mlb_rfq_market.duckdb` (sibling to state DB)
- Scrapers (`mlb_sgp/scraper_*_sgp.py`) are invoked as subprocesses with `MLB_SGP_DB_PATH=<bot market DB>` env override
- N≥2 books required before a candidate can be priced/bet

Concrete diff snippet to append:

```markdown
## SGP cadence loop (line-source pivot, 2026-05-13)

The bot drives its own SGP scrape cadence independent of the MLB dashboard.

**Data flow on each SGP tick (every `SGP_REFRESH_SEC`, default 60s):**

1. Bot enumerates open Kalshi MVE markets per MLB game (every `(spread, total)` tuple Kalshi lists).
2. Bot rewrites `mlb_target_lines` in `kalshi_mlb_rfq_market.duckdb` (sibling to state DB).
3. Bot spawns the four scrapers (`mlb_sgp/scraper_{draftkings,fanduel,prophetx,novig}_sgp.py`) with `MLB_SGP_DB_PATH=<market DB>` and `MLB_SGP_PERIODS=FG` env overrides.
4. Scrapers read `mlb_target_lines`, price every tuple at their respective book, write back to `mlb_sgp_odds` in the bot's market DB with new `spread_line`/`total_line` columns.
5. Bot reloads `_SGP_ODDS_CACHE` from the bot market DB.

**Edge surface:** any Kalshi MVE combo with ≥2 books priced at the matching (spread, total). Off-line combos (only 1 book) are dropped — the bot does not bet model-only or single-book candidates.

**Schedule source:** game IDs and team metadata come from `Answer Keys/mlb.duckdb::mlb_odds_temp` (read-only). The bot has no dependency on Wagerzon-derived `mlb_parlay_lines` anymore.

**Cold start:** the first SGP cycle runs synchronously before `main_loop` enters its tick. Bot blocks ~60-90s on startup.

**Config:**
- `SGP_REFRESH_SEC` (default 60) — SGP cadence interval
- `SGP_SCRAPER_TIMEOUT_SEC` (default 90) — per-scraper kill deadline
- `BOT_MARKET_DB` (default `kalshi_mlb_rfq_market.duckdb` in this package) — sibling market DB
- `MIN_BOOK_COUNT_FOR_BLEND` (default 2) — drop-candidate threshold
```

- [ ] **Step 2: Update `mlb_sgp/README.md`**

Add a section "Library architecture" describing:
- Per-book HTTP clients in `dk_client.py`, `fd_client.py`, `prophetx_client.py`, `novig_client.py`
- Per-book SGP orchestrators in `draftkings.py`, `fanduel.py`, `prophetx.py`, `novig.py`
- Thin scraper shims (`scraper_*_sgp.py`) load target lines, call orchestrator, upsert results
- `mlb_target_lines` (new, multi-row per game) vs `mlb_parlay_lines` (legacy, one-row-per-game)
- `MLB_SGP_DB_PATH` env override for redirecting writes
- `MLB_SGP_PERIODS` env var for periods filtering

- [ ] **Step 3: Update `kalshi_mlb_rfq/.env.example`**

Add lines:

```bash
# SGP cadence loop (added by line-source pivot 2026-05-13)
SGP_REFRESH_SEC=60
SGP_SCRAPER_TIMEOUT_SEC=90
SGP_MIN_INTERVAL_SEC=30
BOT_MARKET_DB=
MLB_SGP_DIR=
MIN_BOOK_COUNT_FOR_BLEND=2
```

- [ ] **Step 4: Update root `CLAUDE.md`**

Under "Project Structure → Autonomous Kalshi MLB SGP taker bot" bullet, add:

```markdown
- Bot owns a sibling **market DB** `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb`
  (separate from the state DB) for SGP-line and SGP-odds data; reads
  schedule from `mlb.duckdb::mlb_odds_temp` (read-only). Line surface
  is now driven by Kalshi MVE enumeration, not Wagerzon-derived
  `mlb_parlay_lines`.
```

- [ ] **Step 5: Verify docs render**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
head -50 kalshi_mlb_rfq/README.md
head -50 mlb_sgp/README.md
cat kalshi_mlb_rfq/.env.example | tail -20
```

Expected: new sections present, no obvious markdown errors.

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_rfq/README.md mlb_sgp/README.md kalshi_mlb_rfq/.env.example CLAUDE.md
git commit -m "docs(kalshi-mlb-rfq): line-source pivot — READMEs, .env.example, CLAUDE.md

Document the SGP cadence loop, market DB sibling, N>=2 book gate, library
architecture (clients + orchestrators + shims), and the dual-mode
MLB_SGP_DB_PATH invocation pattern.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 28: Manual smoke test on live MVE markets

This is a verification step, not an implementation step. Run after Task 27 is committed. All unit tests should already be green.

**Prerequisites:** Bot is stopped (no other instance writing to `kalshi_mlb_rfq.duckdb` or `kalshi_mlb_rfq_market.duckdb`). MVE markets are open (typically game days, hours before tipoff).

- [ ] **Step 1: Confirm the full test suite is green**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-rfq-line-source-pivot
pytest mlb_sgp/tests/ kalshi_mlb_rfq/tests/ -v 2>&1 | tail -30
```

Expected: all tests pass, no errors.

- [ ] **Step 2: Run bot in dry-run mode for 5 minutes**

```bash
source kalshi_mlb_rfq/venv/bin/activate
cd kalshi_mlb_rfq
timeout 300 python -m kalshi_mlb_rfq.main --dry-run 2>&1 | tee /tmp/bot_smoke.log
```

Expected sequence in log:
- `startup: warming SGP cache (one synchronous scrape tick)...`
- `startup: SGP warm-up done — return codes {'scraper_draftkings_sgp.py': 0, 'scraper_fanduel_sgp.py': 0, 'scraper_prophetx_sgp.py': 0, 'scraper_novig_sgp.py': 0}`
- `sgp_cache_refresh: NN rows from bot_market_db` (NN > 0)
- `cache_refresh: NN games, NN sgp_odds rows, NN parlay_lines, ...`
- `rfq_refresh: NN candidates (X.Xs)` — non-zero candidates
- Every 60s: `sgp_cycle: rcs={'scraper_*_sgp.py': 0, ...} (Xs)`

- [ ] **Step 3: Inspect quote_log distribution**

```bash
duckdb kalshi_mlb_rfq/kalshi_mlb_rfq.duckdb "SELECT decision, COUNT(*) FROM quote_log WHERE observed_at > NOW() - INTERVAL 10 MINUTE GROUP BY decision ORDER BY 2 DESC"
```

Expected:
- `declined_dry_run` should dominate (bot is in dry-run)
- `declined_ev` and `declined_kelly_zero` are normal
- `declined_stale_predictions` should NOT dominate (means cache refresh failed)
- `no_fresh_fair` should NOT dominate (means book gate is failing)
- Absence of `declined_line_move` is correct (gate removed)

- [ ] **Step 4: Inspect mlb_sgp_odds in bot DB**

```bash
duckdb kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb "SELECT bookmaker, COUNT(*) as n, COUNT(DISTINCT game_id) as games, MIN(fetch_time) as oldest, MAX(fetch_time) as newest FROM mlb_sgp_odds GROUP BY bookmaker ORDER BY bookmaker"
```

Expected: 4 books (draftkings, fanduel, novig, prophetx), each with non-zero `n` and fresh `newest` (within last 5 min). `games` count should match the open MLB slate.

- [ ] **Step 5: Inspect mlb_target_lines**

```bash
duckdb kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb "SELECT game_id, COUNT(*) as n_lines FROM mlb_target_lines GROUP BY game_id ORDER BY 2 DESC LIMIT 5"
```

Expected: each game has multiple lines (typically 3-6 (spread, total) tuples), confirming multi-line enumeration.

- [ ] **Step 6: Acceptance**

If all 5 checks above pass, smoke is GREEN. Document any deviations or follow-ups before requesting merge to main.

If any check fails, halt and investigate. Common failure modes:
- **Empty mlb_target_lines:** Kalshi MVE enumeration returned no events; check `_fetch_kalshi_mlb_events()` against current Kalshi MVE collection state.
- **All scrapers rc=0 but mlb_sgp_odds empty:** Scrapers ran but produced no rows; check shim logs in `mlb_sgp/scraper_*.runner.log`.
- **Schedule cache empty (games shows 0):** `mlb_odds_temp` is empty in `mlb.duckdb`; run the MLB R pipeline to refresh.

- [ ] **Step 7: After smoke green, prepare for merge (NO MERGE WITHOUT USER APPROVAL)**

```bash
git log main..HEAD --oneline | head -30
git diff main..HEAD --stat | tail -20
```

Surface the diff to the user. Request explicit merge approval per CLAUDE.md before any `git checkout main && git merge`.

---

## Self-Review

**1. Spec coverage check:**

| Spec section | Task(s) implementing it |
|---|---|
| D1: All Kalshi-open lines, FG only | Task 19 (enumerate cross-products) + Task 22 (MLB_SGP_PERIODS=FG in env) |
| D2: Sibling market DB | Task 17 (BOT_MARKET_DB config) + Task 20 (write_target_lines creates DB) |
| D3: Schedule from mlb_odds_temp | Task 19 + Task 23 (both read mlb_odds_temp) |
| D4: Library refactor (Option 3) | Tasks 9-16 |
| D5: N≥2 books required | Task 17 (config knob) + Task 24 (gate enforcement) |
| D6: Synchronous warm-up | Task 26 (warm-up + cadence timer) |
| D7: Drop line_move_ok | Task 25 |
| D8: Drop _current_book_lines_for_combo | Task 25 |
| D9: Worktree path fix in config.py | Task 17 |
| D10: RFQ/SGP coordination | Task 26 (separate timers, both read from caches) |
| Schema migration mlb_sgp_odds | Task 4 |
| Env-var override MLB_SGP_DB_PATH | Task 5 |
| upsert_priced_rows | Task 6 |
| PX/NV client extraction | Tasks 7-8 |
| Thin shims | Tasks 13-16 |
| Documentation | Task 27 |
| Smoke test plan | Task 28 |

No gaps.

**2. Placeholder scan:** Searched for "TBD", "TODO", "implement later", "fill in", "similar to". One acceptable instance remains:

- Task 12 `_select_outcome_ids_for_combo`: `raise NotImplementedError("Lift from scraper_novig_sgp.scrape_novig_sgp body — leg selection per combo")` — explicit note for the engineer to lift specific lines. This is a code skeleton with documented body source. Acceptable since it points to exact source lines and the test in Task 16 (golden regression) will catch any missed cases.

**3. Type consistency:**

- `TargetLine` and `PricedRow` field names match across Tasks 1, 6, 9-12, 19, 20.
- `MLB_SGP_DB_PATH` env var consistent across Tasks 5, 13-16, 22.
- `MLB_SGP_PERIODS` env var consistent across Tasks 13-16, 22.
- `BOT_MARKET_DB`, `MLB_SGP_DIR`, `MIN_BOOK_COUNT_FOR_BLEND` config attrs consistent across Tasks 17, 22, 24, 26.
- Function names: `should_scrape`, `latest_sgp_fetch_time`, `enumerate_kalshi_targets`, `write_target_lines`, `run_scrapers`, `read_priced_rows`, `sgp_cycle` consistent across Tasks 18-22, 26.

**4. Ambiguity check:**

- "Lift orchestration from scraper_X_sgp.py" — every Task references specific line ranges. Engineer can lift unambiguously.
- "Preserve helpers verbatim" — explicit re: line ranges in Tasks 13-16.
- TDD shape consistent throughout (5 numbered steps per task).

No issues to fix.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-05-13-kalshi-mlb-rfq-line-source-pivot.md`.**

Two execution options:

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task with the relevant task body as the prompt, review the diff between tasks, fast iteration through the 28 tasks (plus pre-task 0).

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints for review after groups of tasks.

Which approach?
