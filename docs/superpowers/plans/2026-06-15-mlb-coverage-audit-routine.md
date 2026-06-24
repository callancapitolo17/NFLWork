# MLB Scraper Coverage Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a deterministic, network-free Python script (`coverage/coverage_audit.py`) that reads every MLB book's local DuckDB and flags coverage regressions, staleness, row-count drops, and books missing from the odds screen — plus the agent playbook + docs that turn it into a daily Desktop scheduled task that wires fixes on worktrees.

**Architecture:** A single CLI script reads each book's `mlb_odds` table read-only, compares the latest run against the book's own trailing baseline (no cross-book normalization → no false positives), writes findings to `coverage/coverage.duckdb::coverage_gaps`, fires a macOS notification only when the gap set changes, and prints a JSON summary. A separate markdown playbook tells the daily agent how to consume that table and wire fixes (fuzzy scraped-vs-screen matching and live menu discovery live in the agent, not the script).

**Tech Stack:** Python 3, `duckdb` (1.4.x, already installed), `pytest` for tests, macOS `osascript` for notifications.

## Global Constraints

- All timestamp columns are `TIMESTAMP WITH TIME ZONE` (UTC). Never introduce naive timestamps. (CLAUDE.md rule #6)
- Never symlink `.duckdb` files; copy if needed. WAL lives next to the path. (CLAUDE.md rule #5)
- All DuckDB connections from the audit are **read-only** against book DBs; only `coverage/coverage.duckdb` is opened for writing. Every connection must be closed in a `finally`.
- No new data/CSV/temp files on disk — state lives in `coverage/coverage.duckdb`. (CLAUDE.md housekeeping)
- The script must never raise on a single book's failure (missing/locked DB, empty table) — it degrades to a reported gap and continues.
- Work happens on worktree branch `worktree-mlb-coverage-audit-routine`; merge to `main` only after review + tests.

## Reality notes (verified against live DBs, 2026-06-15)

- Each book DB lives at `<book>_odds/<book>.duckdb`, table `mlb_odds`, with columns
  `fetch_time TIMESTAMPTZ`, `market VARCHAR`, `period VARCHAR`, `sport_key`, plus odds columns.
- One scrape run writes one distinct `fetch_time` per book (verified: `COUNT(DISTINCT fetch_time)` over a day = 1 per run). So "latest run" = rows where `fetch_time = (SELECT MAX(fetch_time) FROM mlb_odds)`.
- The screen table is `Answer Keys/mlb_mm.duckdb::mlb_bets_book_prices`, columns include `bookmaker`, `market`, `period`, `fetch_time`.
- Only 5 `bookmaker` values ever appear on the screen: `bet105, bookmaker, draftkings, fanduel, wagerzon`. Hoop88/BFA/Kalshi are scraped but intentionally not pill-rendered.
- Label conventions differ per book (e.g. WZ `alternate_spreads_fg`/`fg` vs FD `alternate_spreads`/`FG`; FD `main` packs several markets in one row). This is why cross-source diffing is the agent's job, not the script's.

## File Structure

- `coverage/__init__.py` — marks package (empty).
- `coverage/registry.py` — `BOOK_REGISTRY` (per-book DB path, table, screen name, expected-on-screen flag, recon hint) + path resolution rooted at repo root.
- `coverage/db.py` — `connect_readonly(path)` and `connect_coverage()` helpers; graceful failure.
- `coverage/detectors.py` — pure functions: latest-run extraction, trailing baseline, regression/freshness/rowcount/screen-presence detectors. All take a DuckDB connection or plain data; no I/O orchestration.
- `coverage/audit.py` — orchestrator: loops books, runs detectors, writes `coverage_gaps`, diffs vs previous, notifies, prints JSON. Has `main()` + CLI.
- `coverage/coverage.duckdb` — output DB (gitignored, created at runtime).
- `coverage/AGENT_PLAYBOOK.md` — the daily agent's instructions.
- `coverage/README.md` — setup, tiers, schema, Desktop-task registration.
- `coverage/tests/` — pytest tests with synthetic DuckDB fixtures.

Detectors are split from orchestration so each detector is unit-testable against a tiny in-memory fixture without touching real DBs or the notification side effect.

---

### Task 1: Package scaffold + book registry

**Files:**
- Create: `coverage/__init__.py`
- Create: `coverage/registry.py`
- Test: `coverage/tests/test_registry.py`

**Interfaces:**
- Produces: `REPO_ROOT: Path`; `BOOK_REGISTRY: list[Book]` where `Book` is a dataclass with fields `name: str`, `db_path: Path` (absolute), `table: str`, `screen_name: str | None`, `expected_on_screen: bool`, `recon_hint: str | None`; `get_books() -> list[Book]`.

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_registry.py
from coverage.registry import get_books, BOOK_REGISTRY

def test_registry_has_nine_books():
    names = {b.name for b in get_books()}
    assert names == {
        "wagerzon", "wagerzon_specials", "hoop88", "bfa", "bookmaker",
        "bet105", "kalshi", "draftkings_singles", "fanduel_singles",
    }

def test_only_five_books_expected_on_screen():
    on_screen = {b.screen_name for b in get_books() if b.expected_on_screen}
    assert on_screen == {"wagerzon", "bookmaker", "bet105", "draftkings", "fanduel"}

def test_db_paths_are_absolute_and_named_correctly():
    by = {b.name: b for b in get_books()}
    assert by["fanduel_singles"].db_path.name == "fd.duckdb"
    assert by["draftkings_singles"].db_path.name == "dk.duckdb"
    assert by["wagerzon"].db_path.is_absolute()
    # wagerzon_specials reads a different table in the same DB
    assert by["wagerzon_specials"].db_path.name == "wagerzon.duckdb"
    assert by["wagerzon_specials"].table == "wagerzon_specials"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/mlb-coverage-audit-routine && python -m pytest coverage/tests/test_registry.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'coverage.registry'`

- [ ] **Step 3: Write minimal implementation**

```python
# coverage/__init__.py
# (empty — marks the package)
```

```python
# coverage/registry.py
from dataclasses import dataclass
from pathlib import Path

# repo root = three levels up from this file: coverage/registry.py -> coverage -> repo root
REPO_ROOT = Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class Book:
    name: str
    db_path: Path
    table: str
    screen_name: str | None      # bookmaker label in mlb_bets_book_prices, or None
    expected_on_screen: bool
    recon_hint: str | None       # how the agent refreshes auth for this book

def _db(rel: str) -> Path:
    return (REPO_ROOT / rel).resolve()

BOOK_REGISTRY: list[Book] = [
    Book("wagerzon", _db("wagerzon_odds/wagerzon.duckdb"), "mlb_odds", "wagerzon", True, "wagerzon_odds/recon_wagerzon.py"),
    Book("wagerzon_specials", _db("wagerzon_odds/wagerzon.duckdb"), "wagerzon_specials", None, False, "wagerzon_odds/recon_wagerzon.py"),
    Book("hoop88", _db("hoop88_odds/hoop88.duckdb"), "mlb_odds", None, False, "hoop88_odds/recon_hoop88.py"),
    Book("bfa", _db("bfa_odds/bfa.duckdb"), "mlb_odds", None, False, "bfa_odds (single API call)"),
    Book("bookmaker", _db("bookmaker_odds/bookmaker.duckdb"), "mlb_odds", "bookmaker", True, "bookmaker_odds/recon_bookmaker.py"),
    Book("bet105", _db("bet105_odds/bet105.duckdb"), "mlb_odds", "bet105", True, "bet105_odds/recon_bet105.py"),
    Book("kalshi", _db("kalshi_odds/kalshi.duckdb"), "mlb_odds", None, False, None),
    Book("draftkings_singles", _db("dk_odds/dk.duckdb"), "mlb_odds", "draftkings", True, None),
    Book("fanduel_singles", _db("fd_odds/fd.duckdb"), "mlb_odds", "fanduel", True, None),
]

def get_books() -> list[Book]:
    return list(BOOK_REGISTRY)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_registry.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/__init__.py coverage/registry.py coverage/tests/test_registry.py
git commit -m "feat(coverage): book registry with per-book DB paths + screen flags"
```

---

### Task 2: Read-only DB helpers

**Files:**
- Create: `coverage/db.py`
- Test: `coverage/tests/test_db.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `connect_readonly(path: Path) -> duckdb.DuckDBPyConnection | None` (returns `None` if the file is missing or the connection fails — e.g. locked/corrupt); `connect_coverage() -> duckdb.DuckDBPyConnection` (read-write, creating `coverage/coverage.duckdb` if absent).

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_db.py
from pathlib import Path
import duckdb
from coverage.db import connect_readonly, connect_coverage

def test_connect_readonly_missing_file_returns_none(tmp_path):
    assert connect_readonly(tmp_path / "nope.duckdb") is None

def test_connect_readonly_opens_existing(tmp_path):
    p = tmp_path / "x.duckdb"
    c = duckdb.connect(str(p)); c.execute("CREATE TABLE t(a INT)"); c.close()
    con = connect_readonly(p)
    assert con is not None
    assert con.execute("SELECT count(*) FROM t").fetchone()[0] == 0
    con.close()

def test_connect_coverage_creates_db(tmp_path, monkeypatch):
    import coverage.db as dbmod
    monkeypatch.setattr(dbmod, "COVERAGE_DB_PATH", tmp_path / "coverage.duckdb")
    con = connect_coverage()
    assert (tmp_path / "coverage.duckdb").exists()
    con.close()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_db.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'coverage.db'`

- [ ] **Step 3: Write minimal implementation**

```python
# coverage/db.py
from pathlib import Path
import duckdb
from coverage.registry import REPO_ROOT

COVERAGE_DB_PATH = (REPO_ROOT / "coverage" / "coverage.duckdb").resolve()

def connect_readonly(path: Path):
    """Open a book DB read-only. Returns None on any failure (missing, locked, corrupt)."""
    if not Path(path).exists():
        return None
    try:
        return duckdb.connect(str(path), read_only=True)
    except Exception:
        return None

def connect_coverage():
    """Open (creating if needed) the coverage output DB read-write."""
    COVERAGE_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(COVERAGE_DB_PATH))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_db.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/db.py coverage/tests/test_db.py
git commit -m "feat(coverage): read-only book DB + coverage DB connection helpers"
```

---

### Task 3: Latest-run + trailing-baseline extractors

**Files:**
- Create: `coverage/detectors.py`
- Test: `coverage/tests/test_extractors.py`
- Create: `coverage/tests/conftest.py` (fixture builder)

**Interfaces:**
- Consumes: a `duckdb` connection, a `table` name, a `days` int.
- Produces:
  - `latest_run(con, table) -> dict | None` → `{"fetch_time": datetime, "markets": set[tuple[str,str]], "row_count": int}` (None if table empty/missing).
  - `trailing_baseline(con, table, days) -> dict` → `{"run_count": int, "market_run_fraction": dict[tuple[str,str], float], "median_row_count": float | None}` where `market_run_fraction[(market,period)]` = fraction of the last `days` of runs in which that pair appeared (0..1).

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/conftest.py
import duckdb, pytest
from datetime import datetime, timezone, timedelta

@pytest.fixture
def book_db(tmp_path):
    """Build a tiny mlb_odds DB. runs = list of (fetch_time, [(market,period),...])."""
    def _build(runs, table="mlb_odds"):
        p = tmp_path / "book.duckdb"
        c = duckdb.connect(str(p))
        c.execute(f"CREATE TABLE {table}(fetch_time TIMESTAMPTZ, market VARCHAR, period VARCHAR)")
        for ft, pairs in runs:
            for m, per in pairs:
                c.execute(f"INSERT INTO {table} VALUES (?, ?, ?)", [ft, m, per])
        c.close()
        return duckdb.connect(str(p), read_only=True)
    return _build

def days_ago(n):
    return datetime.now(timezone.utc) - timedelta(days=n)
```

```python
# coverage/tests/test_extractors.py
from coverage.detectors import latest_run, trailing_baseline
from coverage.tests.conftest import days_ago

def test_latest_run_picks_max_fetch_time(book_db):
    con = book_db([
        (days_ago(2), [("spreads","FG"), ("totals","FG")]),
        (days_ago(0), [("spreads","FG")]),
    ])
    lr = latest_run(con, "mlb_odds")
    assert lr["markets"] == {("spreads","FG")}
    assert lr["row_count"] == 1

def test_latest_run_empty_table_returns_none(book_db):
    con = book_db([])
    assert latest_run(con, "mlb_odds") is None

def test_trailing_baseline_fractions(book_db):
    con = book_db([
        (days_ago(3), [("spreads","FG"), ("totals","FG")]),
        (days_ago(2), [("spreads","FG"), ("totals","FG")]),
        (days_ago(1), [("spreads","FG")]),               # totals missing once
    ])
    b = trailing_baseline(con, "mlb_odds", days=7)
    assert b["run_count"] == 3
    assert b["market_run_fraction"][("spreads","FG")] == 1.0
    assert abs(b["market_run_fraction"][("totals","FG")] - 2/3) < 1e-9
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_extractors.py -v`
Expected: FAIL with `ImportError: cannot import name 'latest_run'`

- [ ] **Step 3: Write minimal implementation**

```python
# coverage/detectors.py
def latest_run(con, table):
    row = con.execute(f"SELECT MAX(fetch_time) FROM {table}").fetchone()
    if row is None or row[0] is None:
        return None
    max_ft = row[0]
    pairs = con.execute(
        f"SELECT DISTINCT market, period FROM {table} WHERE fetch_time = ?", [max_ft]
    ).fetchall()
    count = con.execute(
        f"SELECT COUNT(*) FROM {table} WHERE fetch_time = ?", [max_ft]
    ).fetchone()[0]
    return {"fetch_time": max_ft, "markets": {(m, p) for m, p in pairs}, "row_count": count}

def trailing_baseline(con, table, days):
    runs = con.execute(
        f"SELECT DISTINCT fetch_time FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days'"
    ).fetchall()
    run_times = [r[0] for r in runs]
    run_count = len(run_times)
    if run_count == 0:
        return {"run_count": 0, "market_run_fraction": {}, "median_row_count": None}
    # fraction of runs each (market,period) appears in
    rows = con.execute(
        f"SELECT market, period, COUNT(DISTINCT fetch_time) AS runs_present FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days' GROUP BY 1,2"
    ).fetchall()
    fraction = {(m, p): present / run_count for m, p, present in rows}
    median_rc = con.execute(
        f"SELECT median(rc) FROM (SELECT fetch_time, COUNT(*) AS rc FROM {table} "
        f"WHERE fetch_time > now() - INTERVAL '{int(days)} days' GROUP BY fetch_time)"
    ).fetchone()[0]
    return {"run_count": run_count, "market_run_fraction": fraction, "median_row_count": median_rc}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_extractors.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/detectors.py coverage/tests/test_extractors.py coverage/tests/conftest.py
git commit -m "feat(coverage): latest-run + trailing-baseline extractors"
```

---

### Task 4: Regression / freshness / row-count detectors

**Files:**
- Modify: `coverage/detectors.py`
- Test: `coverage/tests/test_detectors.py`

**Interfaces:**
- Consumes: `latest_run`/`trailing_baseline` outputs from Task 3.
- Produces a `Gap` dataclass `{book, gap_type, severity, market, period, detail, metric_value, baseline_value}` and three functions returning `list[Gap]`:
  - `detect_regressions(book_name, latest, baseline, presence_threshold=0.8) -> list[Gap]` — pairs present in ≥`presence_threshold` of trailing runs but absent from the latest run. `gap_type="regression"`, `severity="alert"`.
  - `detect_freshness(book_name, latest, max_age_hours=26.0) -> list[Gap]` — latest run older than `max_age_hours` (or no latest at all). `gap_type="freshness"`. (26h default = a daily cadence missing one full day.)
  - `detect_rowcount(book_name, latest, baseline, min_ratio=0.5) -> list[Gap]` — latest row count below `min_ratio × median`. `gap_type="rowcount"`, `severity="warn"`.

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_detectors.py
from datetime import datetime, timezone, timedelta
from coverage.detectors import Gap, detect_regressions, detect_freshness, detect_rowcount

NOW = datetime.now(timezone.utc)

def test_regression_flags_disappeared_market():
    latest = {"fetch_time": NOW, "markets": {("spreads","FG")}, "row_count": 10}
    baseline = {"run_count": 5, "market_run_fraction": {
        ("spreads","FG"): 1.0, ("alternate_totals","FG"): 1.0}, "median_row_count": 10}
    gaps = detect_regressions("fanduel_singles", latest, baseline)
    assert len(gaps) == 1
    g = gaps[0]
    assert (g.market, g.period) == ("alternate_totals","FG")
    assert g.gap_type == "regression" and g.severity == "alert"

def test_regression_ignores_rarely_seen_market():
    latest = {"fetch_time": NOW, "markets": set(), "row_count": 0}
    baseline = {"run_count": 10, "market_run_fraction": {("odd_even","FG"): 0.2}, "median_row_count": 1}
    assert detect_regressions("wagerzon", latest, baseline) == []

def test_freshness_flags_stale_and_missing():
    stale = {"fetch_time": NOW - timedelta(hours=30), "markets": set(), "row_count": 0}
    assert detect_freshness("bfa", stale)[0].gap_type == "freshness"
    assert detect_freshness("bfa", None)[0].severity == "alert"

def test_rowcount_flags_drop():
    latest = {"fetch_time": NOW, "markets": set(), "row_count": 40}
    baseline = {"run_count": 7, "market_run_fraction": {}, "median_row_count": 600}
    g = detect_rowcount("fanduel_singles", latest, baseline)
    assert g[0].gap_type == "rowcount" and g[0].metric_value == 40
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_detectors.py -v`
Expected: FAIL with `ImportError: cannot import name 'Gap'`

- [ ] **Step 3: Write minimal implementation** (append to `coverage/detectors.py`)

```python
from dataclasses import dataclass
from datetime import datetime, timezone

@dataclass
class Gap:
    book: str
    gap_type: str          # regression | freshness | rowcount | book_absent
    severity: str          # alert | warn
    detail: str
    market: str | None = None
    period: str | None = None
    metric_value: float | None = None
    baseline_value: float | None = None

def detect_regressions(book_name, latest, baseline, presence_threshold=0.8):
    if latest is None or baseline["run_count"] == 0:
        return []
    gaps = []
    for (m, p), frac in baseline["market_run_fraction"].items():
        if frac >= presence_threshold and (m, p) not in latest["markets"]:
            gaps.append(Gap(
                book=book_name, gap_type="regression", severity="alert",
                market=m, period=p, metric_value=0.0, baseline_value=frac,
                detail=f"{m}/{p} present in {frac:.0%} of last runs but absent today "
                       f"(possible parse break or auth expiry — check recon)"))
    return gaps

def detect_freshness(book_name, latest, max_age_hours=26.0):
    if latest is None:
        return [Gap(book=book_name, gap_type="freshness", severity="alert",
                    detail="no rows at all in book DB (scraper never wrote / DB unreadable)")]
    age_h = (datetime.now(timezone.utc) - latest["fetch_time"]).total_seconds() / 3600.0
    if age_h > max_age_hours:
        return [Gap(book=book_name, gap_type="freshness", severity="alert",
                    metric_value=round(age_h, 1), baseline_value=max_age_hours,
                    detail=f"latest run is {age_h:.1f}h old (> {max_age_hours}h threshold)")]
    return []

def detect_rowcount(book_name, latest, baseline, min_ratio=0.5):
    med = baseline.get("median_row_count")
    if latest is None or not med:
        return []
    if latest["row_count"] < min_ratio * med:
        return [Gap(book=book_name, gap_type="rowcount", severity="warn",
                    metric_value=float(latest["row_count"]), baseline_value=float(med),
                    detail=f"row count {latest['row_count']} < {min_ratio:.0%} of trailing "
                           f"median {med:.0f} (possible partial parse break)")]
    return []
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_detectors.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/detectors.py coverage/tests/test_detectors.py
git commit -m "feat(coverage): regression/freshness/rowcount detectors"
```

---

### Task 5: Screen-presence detector

**Files:**
- Modify: `coverage/detectors.py`
- Test: `coverage/tests/test_screen_presence.py`

**Interfaces:**
- Consumes: a set of `bookmaker` labels present in the latest screen run; `BOOK_REGISTRY`.
- Produces: `screen_bookmakers(con) -> set[str]` (distinct `bookmaker` at the screen table's MAX `fetch_time`); `detect_screen_absence(books, screen_labels) -> list[Gap]` — any `expected_on_screen` book whose `screen_name` is absent. `gap_type="book_absent"`, `severity="alert"`.

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_screen_presence.py
import duckdb
from datetime import datetime, timezone, timedelta
from coverage.detectors import screen_bookmakers, detect_screen_absence
from coverage.registry import get_books

def test_screen_bookmakers_at_latest_only(tmp_path):
    p = tmp_path / "mm.duckdb"; c = duckdb.connect(str(p))
    c.execute("CREATE TABLE mlb_bets_book_prices(bookmaker VARCHAR, fetch_time TIMESTAMPTZ)")
    now = datetime.now(timezone.utc)
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('wagerzon', ?)", [now - timedelta(days=1)])
    c.execute("INSERT INTO mlb_bets_book_prices VALUES ('fanduel', ?)", [now])
    c.close()
    con = duckdb.connect(str(p), read_only=True)
    assert screen_bookmakers(con) == {"fanduel"}
    con.close()

def test_detect_screen_absence_flags_expected_missing():
    gaps = detect_screen_absence(get_books(), screen_labels={"wagerzon", "fanduel"})
    absent = {g.book for g in gaps}
    # bet105, bookmaker, draftkings expected but missing; hoop88/bfa/kalshi never expected
    assert "draftkings_singles" in absent and "bet105" in absent
    assert "hoop88" not in absent and "kalshi" not in absent
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_screen_presence.py -v`
Expected: FAIL with `ImportError: cannot import name 'screen_bookmakers'`

- [ ] **Step 3: Write minimal implementation** (append to `coverage/detectors.py`)

```python
def screen_bookmakers(con, table="mlb_bets_book_prices"):
    row = con.execute(f"SELECT MAX(fetch_time) FROM {table}").fetchone()
    if row is None or row[0] is None:
        return set()
    labels = con.execute(
        f"SELECT DISTINCT bookmaker FROM {table} WHERE fetch_time = ?", [row[0]]
    ).fetchall()
    return {l[0] for l in labels}

def detect_screen_absence(books, screen_labels):
    gaps = []
    for b in books:
        if b.expected_on_screen and b.screen_name not in screen_labels:
            gaps.append(Gap(
                book=b.name, gap_type="book_absent", severity="alert",
                detail=f"book '{b.screen_name}' expected on odds screen but absent "
                       f"from latest mlb_bets_book_prices run"))
    return gaps
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_screen_presence.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/detectors.py coverage/tests/test_screen_presence.py
git commit -m "feat(coverage): odds-screen presence detector"
```

---

### Task 6: Coverage-gaps storage + gap-set change diff

**Files:**
- Create: `coverage/store.py`
- Test: `coverage/tests/test_store.py`

**Interfaces:**
- Consumes: `connect_coverage` (Task 2), `Gap` (Task 4).
- Produces:
  - `ensure_schema(con)` — creates `coverage_gaps(audit_ts TIMESTAMPTZ, book, gap_type, severity, market, period, detail, metric_value DOUBLE, baseline_value DOUBLE)` if absent.
  - `gap_key(g) -> tuple` = `(book, gap_type, market, period)`.
  - `previous_gap_keys(con) -> set[tuple]` — keys from the most recent prior `audit_ts` (empty on first run).
  - `write_gaps(con, gaps, audit_ts) -> None` — inserts all gaps stamped with `audit_ts`.
  - `new_gap_keys(prev_keys, gaps) -> set[tuple]` — keys in `gaps` not in `prev_keys`.

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_store.py
import duckdb
from datetime import datetime, timezone, timedelta
from coverage.detectors import Gap
from coverage.store import ensure_schema, write_gaps, previous_gap_keys, new_gap_keys, gap_key

def _con(tmp_path):
    return duckdb.connect(str(tmp_path / "cov.duckdb"))

def test_write_and_previous_keys(tmp_path):
    con = _con(tmp_path); ensure_schema(con)
    t1 = datetime.now(timezone.utc) - timedelta(days=1)
    write_gaps(con, [Gap("bfa","freshness","alert","stale")], t1)
    assert previous_gap_keys(con) == {("bfa","freshness",None,None)}

def test_previous_keys_uses_only_latest_audit(tmp_path):
    con = _con(tmp_path); ensure_schema(con)
    t1 = datetime.now(timezone.utc) - timedelta(days=2)
    t2 = datetime.now(timezone.utc) - timedelta(days=1)
    write_gaps(con, [Gap("bfa","freshness","alert","old")], t1)
    write_gaps(con, [Gap("hoop88","rowcount","warn","drop")], t2)
    assert previous_gap_keys(con) == {("hoop88","rowcount",None,None)}

def test_new_gap_keys_detects_only_fresh():
    prev = {("bfa","freshness",None,None)}
    gaps = [Gap("bfa","freshness","alert","still"), Gap("fanduel_singles","regression","alert","x",market="totals",period="FG")]
    assert new_gap_keys(prev, gaps) == {("fanduel_singles","regression","totals","FG")}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_store.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'coverage.store'`

- [ ] **Step 3: Write minimal implementation**

```python
# coverage/store.py
def ensure_schema(con):
    con.execute("""
        CREATE TABLE IF NOT EXISTS coverage_gaps(
            audit_ts TIMESTAMPTZ, book VARCHAR, gap_type VARCHAR, severity VARCHAR,
            market VARCHAR, period VARCHAR, detail VARCHAR,
            metric_value DOUBLE, baseline_value DOUBLE)
    """)

def gap_key(g):
    return (g.book, g.gap_type, g.market, g.period)

def previous_gap_keys(con):
    row = con.execute("SELECT MAX(audit_ts) FROM coverage_gaps").fetchone()
    if row is None or row[0] is None:
        return set()
    keys = con.execute(
        "SELECT book, gap_type, market, period FROM coverage_gaps WHERE audit_ts = ?", [row[0]]
    ).fetchall()
    return {(b, t, m, p) for b, t, m, p in keys}

def write_gaps(con, gaps, audit_ts):
    for g in gaps:
        con.execute(
            "INSERT INTO coverage_gaps VALUES (?,?,?,?,?,?,?,?,?)",
            [audit_ts, g.book, g.gap_type, g.severity, g.market, g.period,
             g.detail, g.metric_value, g.baseline_value])

def new_gap_keys(prev_keys, gaps):
    return {gap_key(g) for g in gaps} - set(prev_keys)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_store.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add coverage/store.py coverage/tests/test_store.py
git commit -m "feat(coverage): coverage_gaps storage + new-gap diff"
```

---

### Task 7: Orchestrator, notification, JSON summary, CLI

**Files:**
- Create: `coverage/audit.py`
- Test: `coverage/tests/test_audit.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `notify(title, message, enabled=True) -> None` — fires `osascript -e 'display notification ...'`; swallows all errors; no-op when `enabled=False`.
  - `run_audit(audit_ts, notify_enabled=True) -> dict` — orchestrates all books + screen check, writes gaps, returns a JSON-serializable summary `{"audit_ts", "total_gaps", "new_gaps", "by_severity", "gaps":[...]}`; fires notification only when `new_gaps > 0`.
  - `main()` — parses `--no-notify` / `--json-only`, calls `run_audit`, prints summary JSON, exits 0.

- [ ] **Step 1: Write the failing test**

```python
# coverage/tests/test_audit.py
import json, subprocess
from datetime import datetime, timezone
import coverage.audit as audit

def test_notify_disabled_is_noop(monkeypatch):
    called = {"n": 0}
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.__setitem__("n", called["n"]+1))
    audit.notify("t", "m", enabled=False)
    assert called["n"] == 0

def test_notify_swallows_errors(monkeypatch):
    def boom(*a, **k): raise OSError("no osascript")
    monkeypatch.setattr(subprocess, "run", boom)
    audit.notify("t", "m", enabled=True)  # must not raise

def test_run_audit_returns_summary_and_writes(monkeypatch, tmp_path):
    # Point coverage DB at tmp and stub the per-book + screen scans to fixed gaps.
    import coverage.db as dbmod
    monkeypatch.setattr(dbmod, "COVERAGE_DB_PATH", tmp_path / "cov.duckdb")
    from coverage.detectors import Gap
    monkeypatch.setattr(audit, "scan_books", lambda: [Gap("bfa","freshness","alert","stale")])
    monkeypatch.setattr(audit, "scan_screen", lambda: [])
    monkeypatch.setattr(audit, "notify", lambda *a, **k: None)
    ts = datetime.now(timezone.utc)
    summary = audit.run_audit(audit_ts=ts, notify_enabled=False)
    assert summary["total_gaps"] == 1
    assert summary["new_gaps"] == 1
    assert json.dumps(summary)  # serializable
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest coverage/tests/test_audit.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'coverage.audit'`

- [ ] **Step 3: Write minimal implementation**

```python
# coverage/audit.py
import argparse, json, subprocess, sys
from dataclasses import asdict
from datetime import datetime, timezone
from coverage.registry import get_books, REPO_ROOT
from coverage.db import connect_readonly, connect_coverage
from coverage import detectors as D
from coverage.store import ensure_schema, previous_gap_keys, write_gaps, new_gap_keys

SCREEN_DB = (REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb").resolve()
BASELINE_DAYS = 7

def notify(title, message, enabled=True):
    if not enabled:
        return
    try:
        subprocess.run(
            ["osascript", "-e", f'display notification "{message}" with title "{title}"'],
            check=False, capture_output=True, timeout=10)
    except Exception:
        pass

def scan_books():
    gaps = []
    for b in get_books():
        con = connect_readonly(b.db_path)
        if con is None:
            gaps.append(D.Gap(book=b.name, gap_type="freshness", severity="alert",
                              detail=f"DB unreadable or missing at {b.db_path}"))
            continue
        try:
            latest = D.latest_run(con, b.table)
            baseline = D.trailing_baseline(con, b.table, BASELINE_DAYS)
            gaps += D.detect_freshness(b.name, latest)
            gaps += D.detect_regressions(b.name, latest, baseline)
            gaps += D.detect_rowcount(b.name, latest, baseline)
        finally:
            con.close()
    return gaps

def scan_screen():
    con = connect_readonly(SCREEN_DB)
    if con is None:
        return [D.Gap(book="(screen)", gap_type="book_absent", severity="alert",
                      detail=f"screen DB unreadable at {SCREEN_DB}")]
    try:
        labels = D.screen_bookmakers(con)
    finally:
        con.close()
    return D.detect_screen_absence(get_books(), labels)

def run_audit(audit_ts=None, notify_enabled=True):
    audit_ts = audit_ts or datetime.now(timezone.utc)
    gaps = scan_books() + scan_screen()
    con = connect_coverage()
    try:
        ensure_schema(con)
        prev = previous_gap_keys(con)
        write_gaps(con, gaps, audit_ts)
    finally:
        con.close()
    fresh = new_gap_keys(prev, gaps)
    by_sev = {}
    for g in gaps:
        by_sev[g.severity] = by_sev.get(g.severity, 0) + 1
    summary = {
        "audit_ts": audit_ts.isoformat(),
        "total_gaps": len(gaps),
        "new_gaps": len(fresh),
        "by_severity": by_sev,
        "gaps": [asdict(g) for g in gaps],
    }
    if fresh:
        notify("MLB Coverage Audit",
               f"{len(fresh)} new gap(s), {len(gaps)} total", enabled=notify_enabled)
    return summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-notify", action="store_true")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args()
    summary = run_audit(notify_enabled=not args.no_notify)
    print(json.dumps(summary, indent=None if args.json_only else 2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest coverage/tests/test_audit.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Full suite + live smoke run**

Run: `python -m pytest coverage/ -v`
Expected: PASS (all tasks' tests)
Run: `python -m coverage.audit --no-notify`
Expected: prints a JSON summary against the REAL book DBs; exits 0. Eyeball that book names appear and no traceback is raised. (This is a real read-only run; it writes one audit row to `coverage/coverage.duckdb`.)

- [ ] **Step 6: Commit**

```bash
git add coverage/audit.py coverage/tests/test_audit.py
git commit -m "feat(coverage): orchestrator + notification + JSON summary + CLI"
```

---

### Task 8: Agent playbook

**Files:**
- Create: `coverage/AGENT_PLAYBOOK.md`

This is the instruction set the daily Desktop scheduled task points at. No tests (it's a prompt document), but it must be concrete.

- [ ] **Step 1: Write the playbook**

````markdown
# MLB Coverage Audit — Daily Agent Playbook

You are the daily MLB scraper coverage agent. Your job: run the deterministic
audit, then fix what it found — safely, on worktrees, never touching `main`.

## 1. Run detection (deterministic, do not skip)
From the repo root:
```
python -m coverage.audit
```
Read the JSON summary it prints. The authoritative gap list is the latest
`audit_ts` in `coverage/coverage.duckdb::coverage_gaps`. If `total_gaps == 0`,
STOP — reply "no gaps" and end. Do not spend tokens hunting.

## 2. Triage each gap
- `gap_type=freshness` or DB-unreadable → likely auth expiry or a dead scraper.
  Run the book's recon (see `coverage/registry.py` `recon_hint`) or `/check-auth`.
  If auth is the cause and needs a manual browser login, you CANNOT fix it
  unattended — leave it for the human (see step 5).
- `gap_type=regression` → a market the book used to post stopped appearing.
  Usual cause: a parse break (a selector/regex/endpoint changed). This is the
  FD-paren-bug class. Fixable in code.
- `gap_type=rowcount` → partial parse break. Investigate like a regression.
- `gap_type=book_absent` → a book that should render on the odds screen is gone.
  Trace whether the scraper wrote rows (Stage B) but R dropped them (Stage C,
  the wiring layer in `Answer Keys/MLB Answer Key/odds_screen.R`).

## 3. Fuzzy scraped-vs-screen check (your job, not the script's)
For each `book_absent` / regression, decide if the market is truly missing or
just wired under a different label. Normalize across conventions before
concluding (period case `fg`→`FG`; market suffix `_fg/_f3/_f7`; FD `main`
packs spread+total+ml). The script deliberately does NOT do this — you do,
because it needs judgment.

## 4. Optional: discover never-seen markets (Tier 3, time-box to 10 min)
For the richest books (FanDuel, DraftKings), fetch the book's live menu and
compare to what the scraper requests. Surface markets the book posts that no
scraper captures. Do NOT rabbit-hole — one pass, then stop.

## 5. Wire fixes — ONE worktree per fixable gap
For each fixable gap:
1. `EnterWorktree` (or `git worktree add`) — never edit the shared checkout.
2. Make the minimal fix (new parser / regex / R market-type union).
3. Run `python tests/timezone_parity_test.py` AND the relevant scraper smoke
   test. To test a scraper DB safely, follow CLAUDE.md rule #5: never symlink a
   `.duckdb` into the worktree — run the scraper against the real DB path from
   the main checkout, or copy the DB in. Do not symlink.
4. Commit to the branch. **Do NOT merge to `main`.**
5. Leave a one-line summary of the branch + what it fixes.

## 6. Report
End with: how many gaps, how many fixed (with branch names), how many need the
human (auth / ambiguous / Tier-3 findings). That summary is what the human
wakes up to.
````

- [ ] **Step 2: Commit**

```bash
git add coverage/AGENT_PLAYBOOK.md
git commit -m "docs(coverage): daily agent playbook for wiring gaps on worktrees"
```

---

### Task 9: README, Desktop-task registration, repo docs, memory

**Files:**
- Create: `coverage/README.md`
- Modify: `CLAUDE.md` (root — Project Structure bullet)
- Modify: `.gitignore` (ensure `coverage/coverage.duckdb` ignored)
- Create: a memory entry (see step 4)

- [ ] **Step 1: Write `coverage/README.md`**

````markdown
# MLB Scraper Coverage Audit

Daily check that every MLB book you scrape is still producing the markets it
used to, is fresh, and (for the 5 pill-rendered books) reaches the odds screen.
Detection is a deterministic, read-only Python script; fixing is done by a daily
Claude Desktop scheduled task following `AGENT_PLAYBOOK.md`.

## What it checks (deterministic, network-free)
- **Regression** — a `(market, period)` the book posted in ≥80% of the last 7
  days of runs but not in today's run (parse break / auth expiry).
- **Freshness** — latest run older than 26h, or DB unreadable/empty.
- **Row-count** — today's run below 50% of the trailing median (partial break).
- **Screen presence** — an `expected_on_screen` book absent from the latest
  `mlb_bets_book_prices` run.

Fuzzy scraped-vs-screen label matching and live "never-seen market" discovery
are intentionally the agent's job at runtime — see `AGENT_PLAYBOOK.md`.

## Run it manually
```
python -m coverage.audit            # writes a row, notifies on new gaps, prints JSON
python -m coverage.audit --no-notify --json-only
```
Output table: `coverage/coverage.duckdb::coverage_gaps` (one row per gap per run).

## Register the daily Desktop scheduled task
This runs locally on your Mac (NOT a `/schedule` cloud routine — that can't
reach your DBs/auth/IP). In the **Claude Code Desktop app**:
1. Sidebar → **Routines** → **New routine** → choose **Local**.
2. Working folder: the repo root.
3. Schedule: daily ~9:00 AM Pacific.
4. Instructions: "Follow `coverage/AGENT_PLAYBOOK.md`."

Caveat: the Desktop app must be running for the task to fire (it does not need
an open chat session). If you ever want detection to run even with the app
closed, schedule `python -m coverage.audit` alone via launchd and keep only the
wiring as a Desktop task.

## Add a book or change thresholds
- New book → add a `Book(...)` row to `coverage/registry.py`.
- Thresholds (`BASELINE_DAYS`, `presence_threshold`, `max_age_hours`,
  `min_ratio`) are arguments/constants in `audit.py` / `detectors.py`.
````

- [ ] **Step 2: Update root `CLAUDE.md`** — add under "Project Structure":

```markdown
- **Scraper coverage audit** (`coverage/`) — daily deterministic check that each
  MLB book still posts the markets it used to (regression), is fresh, and reaches
  the odds screen. Read-only over per-book DuckDBs; writes `coverage/coverage.duckdb`.
  A Claude Desktop scheduled task follows `coverage/AGENT_PLAYBOOK.md` to wire
  fixes on worktrees (never auto-merges). See `coverage/README.md`.
```

- [ ] **Step 3: Ensure `.gitignore` covers the output DB**

Run: `grep -q "coverage/coverage.duckdb" .gitignore || echo "coverage/coverage.duckdb" >> .gitignore`
(The repo already ignores `*.duckdb`, but add the explicit path for clarity.)

- [ ] **Step 4: Write the memory entry**

Create `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_coverage_audit.md`:

```markdown
---
name: mlb-coverage-audit
description: Daily MLB scraper coverage-audit routine — deterministic detection + agent wiring; cloud routines can't reach local scrapers
metadata:
  type: project
---

Built 2026-06-15 (branch `worktree-mlb-coverage-audit-routine`). `coverage/`
holds a deterministic, read-only Python audit (`python -m coverage.audit`) that
flags per-book regression / freshness / rowcount / screen-absence into
`coverage/coverage.duckdb::coverage_gaps`, notifying via osascript only on NEW
gaps. A Claude **Desktop scheduled task** (NOT a `/schedule` cloud routine —
those run in Anthropic cloud with only a git clone, no local DB/auth/IP) follows
`coverage/AGENT_PLAYBOOK.md` to wire fixes on worktrees, never auto-merging.

Key facts: only 5 books render on the odds screen (wagerzon, bet105, bookmaker,
draftkings, fanduel); hoop88/bfa/kalshi feed pricing only. Cross-source
scraped-vs-screen label diffing is the agent's job (FD `main` packs multiple
markets; period case + `_fg` suffixes differ) — the script only does per-book
self-baseline so it never false-alarms. Related: [[timezone_standardization]],
[[fd_alt_total_paren_format]], [[fd_event_page_tab_coverage_gap]].
```

Then add to `MEMORY.md` under "Active Work":
```markdown
- [MLB coverage audit](mlb_coverage_audit.md) — daily deterministic scraper coverage check (`coverage/`) + Desktop-task agent wiring on worktrees; cloud routines can't reach local scrapers
```

- [ ] **Step 5: Commit**

```bash
git add coverage/README.md CLAUDE.md .gitignore
git commit -m "docs(coverage): README, Desktop-task registration, CLAUDE.md + memory"
```

---

## Worktree / Version Control

- All work on `worktree-mlb-coverage-audit-routine` (already active).
- New dir `coverage/` is additive; the only existing file modified is root
  `CLAUDE.md` and `.gitignore` — low merge-conflict risk.
- After all tasks: run `python -m pytest coverage/ -v` + one live
  `python -m coverage.audit --no-notify`, then request review before merge.
  Never merge to `main` without explicit approval.

## Documentation (in this branch)

- `coverage/README.md` (new), `coverage/AGENT_PLAYBOOK.md` (new), root
  `CLAUDE.md` bullet, memory entry. All in the same branch as the code.

## Self-Review Notes

- **Spec coverage:** Tier 1 regression → Task 4; freshness/rowcount ride-along →
  Task 4; screen presence (refinement) → Task 5; storage + change-only notify →
  Tasks 6–7; Desktop task + never-merge wiring → Tasks 8–9; tiers 2–3 → playbook.
- **Deviation from spec (flagged):** deterministic Stage-A−B / B−C cross-source
  diffing moved from the script to the agent playbook, because live data showed
  cross-book label packing makes it false-positive-prone. Spec updated to match.
- **Out of scope (unchanged):** cloud routines, fast refresh loop, auto-merge,
  the Unabated screen itself.
````
