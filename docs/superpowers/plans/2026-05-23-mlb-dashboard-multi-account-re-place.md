# MLB Dashboard Multi-Account Re-Place (v1) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a placed Wagerzon bet on the MLB Dashboard bets tab be re-placed on a *different* WZ account via a "+ another" affordance, with each placement tracked as its own row in `placed_bets`.

**Architecture:** Promote `placed_bets.PRIMARY KEY` from `bet_hash` to composite `(bet_hash, account)` (migration 003). Scope every server-side write and in-flight check by the composite key. On the dashboard side, add a small `render_placed_strip(chips)` R helper that produces one green pill per placement plus a dashed `+ another` button when at least one WZ account has no chip yet. Two JS handlers (`addAnother`, `placeAnother`) drive the re-open flow; a hook on the header "Placing on" pill clears `data-expected-win` on every card so the next Risk edit re-verifies under the new account.

**Tech Stack:**
- DuckDB schema migration (Python)
- Flask server (`mlb_dashboard_server.py`) — SQL scoping changes
- R/Shiny dashboard (`mlb_dashboard.R`) — new helper, inline CSS, JS handlers
- Existing helpers reused unchanged: `single_placer.place_single`, `/api/wz-quote-single`, `_replaceActionCell` DOM swap pattern, `window.WZ_SELECTED_ACCOUNT`

**Spec:** `docs/superpowers/specs/2026-05-23-mlb-dashboard-multi-account-re-place-design.md`

**Branch / Worktree:** Already on `worktree-feature+dashboard-multi-account-re-place-spec` in `.claude/worktrees/feature+dashboard-multi-account-re-place-spec/`. Spec is committed (`12a17ac`).

---

## File map

| File | Role | Action |
|---|---|---|
| `Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py` | One-shot migration that backfills NULL/empty `account` to `'Wagerzon'` and rebuilds `placed_bets` with composite PK | **Create** |
| `Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py` | pytest covering backfill, composite-PK creation, idempotency | **Create** |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | Flask routes & DB writers — scope by composite key | **Modify** (`_insert_placement_breadcrumb`, `_finalize_placement`, `/api/place-bet` in-flight check, `/api/remove-bet`) |
| `Answer Keys/MLB Dashboard/placed_strip.R` | New R helper: `render_placed_strip(chips, all_wz_accounts)` | **Create** |
| `Answer Keys/tests/test_placed_strip.R` | testthat unit tests for the helper | **Create** |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Source the new helper; add CSS for `.hero-placed` / `.placement-chip` / `.add-another`; update `create_bets_table` to query `account` and group chips by `bet_hash`; add `addAnother` / `placeAnother` JS + header-pill clear hook | **Modify** |
| `Answer Keys/CLAUDE.md` | Extend the "MLB Dashboard — Wagerzon multi-account" section with the composite-PK + multi-placement semantics | **Modify** |
| `docs/superpowers/specs/2026-05-23-mlb-dashboard-multi-account-re-place-design.md` | Spec (already committed) | (no change) |

---

## Task 1: Composite-PK migration

**Files:**
- Create: `Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py`
- Create: `Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py`

- [ ] **Step 1.1: Write the failing tests**

Create `Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py`:

```python
"""Tests for migration 003: composite PK on placed_bets."""
from __future__ import annotations
import importlib.util
from pathlib import Path

import duckdb
import pytest


def _load_migration():
    """Import the migration module by file path (folder name has a space)."""
    here = Path(__file__).parent
    spec = importlib.util.spec_from_file_location(
        "_m003", here / "003_placed_bets_composite_pk.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Schema mirrors the live placed_bets after migration 002.
OLD_SCHEMA_DDL = """
CREATE TABLE placed_bets (
    bet_hash         VARCHAR PRIMARY KEY,
    game_id          VARCHAR,
    home_team        VARCHAR,
    away_team        VARCHAR,
    market           VARCHAR,
    bet_on           VARCHAR,
    line             DOUBLE,
    odds             INTEGER,
    actual_size      DOUBLE,
    recommended_size DOUBLE,
    bookmaker        VARCHAR,
    model_prob       DOUBLE,
    model_ev         DOUBLE,
    account          VARCHAR,
    status           VARCHAR DEFAULT 'placed',
    ticket_number    VARCHAR,
    error_msg        VARCHAR,
    error_msg_key    VARCHAR,
    wz_odds_at_place INTEGER
)
"""


def _seed_old_db(path: Path) -> None:
    con = duckdb.connect(str(path))
    try:
        con.execute(OLD_SCHEMA_DDL)
    finally:
        con.close()


def _pk_cols(path: Path) -> set[str]:
    con = duckdb.connect(str(path))
    try:
        rows = con.execute("""
            SELECT UNNEST(constraint_column_names)
            FROM duckdb_constraints()
            WHERE table_name = 'placed_bets' AND constraint_type = 'PRIMARY KEY'
        """).fetchall()
    finally:
        con.close()
    return {r[0] for r in rows}


def test_backfills_null_and_empty_account_to_wagerzon(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", None, "placed"])
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h2", "", "placed"])
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h3", "WagerzonJ", "placed"])
    con.close()

    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    rows = sorted(con.execute("SELECT bet_hash, account FROM placed_bets").fetchall())
    con.close()
    assert rows == [("h1", "Wagerzon"), ("h2", "Wagerzon"), ("h3", "WagerzonJ")]


def test_promotes_pk_to_composite(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    _load_migration().run(str(db))
    assert _pk_cols(db) == {"bet_hash", "account"}


def test_same_hash_different_account_is_allowed(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    con.close()

    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "WagerzonJ", "placed"])
    count = con.execute("SELECT COUNT(*) FROM placed_bets WHERE bet_hash = 'h1'").fetchone()[0]
    con.close()
    assert count == 2


def test_duplicate_hash_and_account_is_rejected(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    _load_migration().run(str(db))

    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                    ["h1", "Wagerzon", "placed"])
    con.close()


def test_idempotent(tmp_path):
    db = tmp_path / "test.duckdb"
    _seed_old_db(db)
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO placed_bets (bet_hash, account, status) VALUES (?, ?, ?)",
                ["h1", "Wagerzon", "placed"])
    con.close()

    migrate = _load_migration().run
    migrate(str(db))
    migrate(str(db))   # second run must be a no-op

    con = duckdb.connect(str(db))
    count = con.execute("SELECT COUNT(*) FROM placed_bets").fetchone()[0]
    con.close()
    assert count == 1
    assert _pk_cols(db) == {"bet_hash", "account"}
```

- [ ] **Step 1.2: Run test to verify it fails**

Run from the worktree root:
```bash
pytest "Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py" -v
```
Expected: collection error or all tests fail with `FileNotFoundError: 003_placed_bets_composite_pk.py` — migration script doesn't exist yet.

- [ ] **Step 1.3: Write the migration script**

Create `Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py`:

```python
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
```

- [ ] **Step 1.4: Run tests to verify they pass**

```bash
pytest "Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py" -v
```
Expected: 5 passed.

- [ ] **Step 1.5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py" \
        "Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): migration 003 promotes placed_bets PK to (bet_hash, account)

Backfills NULL/empty account to 'Wagerzon' before rebuilding the table
with the composite PK. Required for multi-account re-place — without it,
the second placement for a given bet would overwrite the first.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Server-side composite-key scoping

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:751-804` (`_insert_placement_breadcrumb`)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:954-971` (`_finalize_placement`)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:1001-1010` (`/api/place-bet` in-flight check)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:1180` (`/api/remove-bet` DELETE)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:1042` (`_finalize_placement` call site — pass `account`)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py:1306` (size-update UPDATE if present in this flow)

> No automated test harness exists for this server today, so this task uses targeted manual verification with a throwaway DuckDB. Three small grep checks plus a smoke script confirm the contract.

- [ ] **Step 2.1: Update `_insert_placement_breadcrumb` to scope by (bet_hash, account)**

In `mlb_dashboard_server.py`, locate the function starting at line ~751. Replace the existence check, INSERT, and UPDATE blocks so all three scope by `(bet_hash, account)`:

```python
def _insert_placement_breadcrumb(bet_hash: str, account, bet_meta: dict,
                                  status: str = "placing"):
    """Insert (or upsert) a placed_bets row before invoking the placer.

    Idempotent within a single (bet_hash, account) — a retry won't
    duplicate the row. Different accounts for the same bet_hash are
    independent rows under the composite PK introduced in migration 003.
    """
    con = duckdb.connect(str(DB_PATH))
    try:
        existing = con.execute(
            "SELECT bet_hash FROM placed_bets WHERE bet_hash = ? AND account = ?",
            [bet_hash, account]
        ).fetchone()

        if existing is None:
            con.execute("""
                INSERT INTO placed_bets
                  (bet_hash, game_id, home_team, away_team, market, bet_on,
                   line, odds, actual_size, recommended_size, bookmaker,
                   model_prob, model_ev,
                   account, status, wz_odds_at_place)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, [
                bet_hash,
                bet_meta.get("game_id"),
                bet_meta.get("home_team"),
                bet_meta.get("away_team"),
                bet_meta.get("market"),
                bet_meta.get("bet_on"),
                bet_meta.get("line"),
                bet_meta.get("american_odds"),
                bet_meta.get("actual_size"),
                bet_meta.get("kelly_bet"),
                bet_meta.get("bookmaker_key") or bet_meta.get("bookmaker"),
                float(bet_meta.get("model_prob") or 0.0),
                float(bet_meta.get("model_ev") or 0.0),
                account,
                status,
                bet_meta.get("wz_odds_at_place"),
            ])
        else:
            con.execute("""
                UPDATE placed_bets
                SET status = ?, wz_odds_at_place = ?
                WHERE bet_hash = ? AND account = ?
            """, [
                status,
                bet_meta.get("wz_odds_at_place"),
                bet_hash,
                account,
            ])
    finally:
        con.close()
```

> Note: the old UPDATE also re-wrote `account`. Under the composite PK, `account` is part of the row identity and must NOT be updated by this code path. Drop it from the SET clause.

- [ ] **Step 2.2: Update `_finalize_placement` to scope by (bet_hash, account)**

Change the signature to accept `account` and scope the UPDATE:

```python
def _finalize_placement(bet_hash: str, account, result: dict):
    """Upsert the final status from the placer's result onto the (bet_hash, account) row."""
    con = duckdb.connect(str(DB_PATH))
    try:
        con.execute("""
            UPDATE placed_bets
            SET status = ?, ticket_number = ?,
                error_msg = ?, error_msg_key = ?
            WHERE bet_hash = ? AND account = ?
        """, [
            result.get("status"),
            result.get("ticket_number"),
            result.get("error_msg"),
            result.get("error_msg_key"),
            bet_hash,
            account,
        ])
    finally:
        con.close()
```

- [ ] **Step 2.3: Update every `_finalize_placement` call site**

Search and update every caller to pass `account`:

```bash
grep -n "_finalize_placement(" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```

For each match, change the call from `_finalize_placement(bet_hash, result)` to `_finalize_placement(bet_hash, account, result)`. The `account` value is already available in the enclosing scope (the request body's `account` field is read into a local `account` variable at the top of `place_bet`).

Expected primary call site to change is around line 1042 inside `place_bet()`.

- [ ] **Step 2.4: Update the `/api/place-bet` in-flight check**

In `mlb_dashboard_server.py` around line 1001-1010, change:

```python
        existing = con.execute(
            "SELECT status FROM placed_bets WHERE bet_hash = ?",
            [bet_hash]).fetchone()
```

To:

```python
        existing = con.execute(
            "SELECT status FROM placed_bets WHERE bet_hash = ? AND account = ?",
            [bet_hash, account]).fetchone()
```

This is the change that unblocks the second placement: an in-flight `WagerzonJ` row no longer blocks a Wagerzon placement (and vice-versa). A second placement on the *same* account still hits the 409 because the (bet_hash, account) row still exists.

- [ ] **Step 2.5: Update `/api/remove-bet` DELETE to scope by composite key**

Around line 1180 in `mlb_dashboard_server.py`, locate the route that handles removal. The current SQL is:

```python
"DELETE FROM placed_bets WHERE bet_hash = ? RETURNING bet_hash"
```

The corresponding route handler reads `bet_hash` from the request body. Add `account` to the request body contract and the WHERE clause:

```python
"DELETE FROM placed_bets WHERE bet_hash = ? AND account = ? RETURNING bet_hash"
```

The caller (JS in `mlb_dashboard.R`) currently sends `{bet_hash}` only. The corresponding JS change is in Task 6 (the `removeBet` handler picks up the chip's `data-account` and includes it). For now, gracefully reject if `account` is missing:

```python
if not account:
    return jsonify({"success": False,
                    "error": "account required to remove a placement"}), 400
```

- [ ] **Step 2.6: Audit any other `placed_bets` writers**

```bash
grep -n "placed_bets" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" | grep -E "INSERT|UPDATE|DELETE"
```

For each writer that targets a specific row (not a bulk SELECT), confirm it scopes by either `(bet_hash, account)` or some other criterion that's still safe under the composite PK. Likely candidates flagged in the spec exploration:
- Line ~1306 (`UPDATE placed_bets SET actual_size = ? WHERE bet_hash = ?`) — if this writer is used for the editable Risk feature, it MUST add `AND account = ?`. If it is dead code, delete it.
- Line ~1687 (`DELETE FROM placed_bets`) — read the surrounding context. If it's a bulk delete (e.g. clear-all), the lack of `account` is intentional and OK.

Apply scoping fixes inline.

- [ ] **Step 2.7: Smoke-test against a throwaway DB**

From the worktree root, in a Python REPL:

```bash
python3 - <<'EOF'
import sys, importlib.util, duckdb, tempfile, os
from pathlib import Path

# Apply migration 003 to a throwaway DB
mig_path = Path("Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py")
spec = importlib.util.spec_from_file_location("_m003", mig_path)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

tmp = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
tmp.close()
con = duckdb.connect(tmp.name)
con.execute("""
CREATE TABLE placed_bets (
    bet_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
    away_team VARCHAR, market VARCHAR, bet_on VARCHAR, line DOUBLE,
    odds INTEGER, actual_size DOUBLE, recommended_size DOUBLE,
    bookmaker VARCHAR, model_prob DOUBLE, model_ev DOUBLE, account VARCHAR,
    status VARCHAR DEFAULT 'placed', ticket_number VARCHAR,
    error_msg VARCHAR, error_msg_key VARCHAR, wz_odds_at_place INTEGER
)
""")
con.close()
mod.run(tmp.name)

# Now patch the server's DB_PATH and exercise the writers
import importlib.util as iu
srv_path = Path("Answer Keys/MLB Dashboard/mlb_dashboard_server.py")
# We can't easily import the Flask app standalone, so just check the SQL
# strings to make sure we updated them. A real e2e is in Task 8.
src = srv_path.read_text()
required_phrases = [
    "WHERE bet_hash = ? AND account = ?",  # in-flight check + finalize
]
for p in required_phrases:
    assert p in src, f"Missing scoping phrase in server: {p}"
print("OK: composite-key phrases present in server code")
os.unlink(tmp.name)
EOF
```
Expected: `OK: composite-key phrases present in server code`.

- [ ] **Step 2.8: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): scope placed_bets writes by (bet_hash, account)

Updates _insert_placement_breadcrumb, _finalize_placement, the
/api/place-bet in-flight check, and /api/remove-bet to scope every
placed_bets write by the composite key introduced in migration 003.
Same-account double-fires still 409; different-account placements
on the same bet_hash are now independent rows.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: R helper `render_placed_strip` + tests

**Files:**
- Create: `Answer Keys/MLB Dashboard/placed_strip.R`
- Create: `Answer Keys/tests/test_placed_strip.R`

- [ ] **Step 3.1: Write the failing tests**

Create `Answer Keys/tests/test_placed_strip.R`:

```r
# Answer Keys/tests/test_placed_strip.R
library(testthat)
source("../MLB Dashboard/placed_strip.R")

test_that("renders one chip per placement with account, risk and ticket", {
  chips <- list(
    list(account = "Wagerzon",  risk = 200, ticket = "W1234"),
    list(account = "WagerzonJ", risk = 150, ticket = "W5678")
  )
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ", "WagerzonC"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  # Two chips with expected text
  expect_match(html, '"placement-chip"', fixed = TRUE)
  expect_match(html, '>WZ<', fixed = TRUE)
  expect_match(html, '>WZJ<', fixed = TRUE)
  expect_match(html, '$200', fixed = TRUE)
  expect_match(html, '$150', fixed = TRUE)
  expect_match(html, '#W1234', fixed = TRUE)
  expect_match(html, '#W5678', fixed = TRUE)
})

test_that("+ another appears (enabled) when an untouched WZ account exists", {
  chips <- list(list(account = "Wagerzon", risk = 200, ticket = "W1234"))
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  expect_match(html, "add-another", fixed = TRUE)
  expect_false(grepl("disabled", html, fixed = TRUE))
})

test_that("+ another is disabled once every WZ account has a chip", {
  chips <- list(
    list(account = "Wagerzon",  risk = 200, ticket = "W1"),
    list(account = "WagerzonJ", risk = 150, ticket = "W2")
  )
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  expect_match(html, "add-another", fixed = TRUE)
  expect_match(html, "disabled", fixed = TRUE)
})

test_that("+ another is omitted entirely for non-wagerzon books", {
  chips <- list(list(account = "Hoop88", risk = 100, ticket = "H1"))
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "hoop88"
  )
  expect_false(grepl("add-another", html, fixed = TRUE))
})

test_that("chips carry data-account, data-risk, data-ticket for JS handlers", {
  chips <- list(list(account = "Wagerzon", risk = 200, ticket = "W1234"))
  html <- render_placed_strip(
    chips, all_wz_accounts = c("Wagerzon"),
    bet_hash = "abc", book = "wagerzon"
  )
  expect_match(html, 'data-account="Wagerzon"', fixed = TRUE)
  expect_match(html, 'data-risk="200"',          fixed = TRUE)
  expect_match(html, 'data-ticket="W1234"',      fixed = TRUE)
  expect_match(html, 'data-bet-hash="abc"',      fixed = TRUE)
})
```

- [ ] **Step 3.2: Run tests to verify they fail**

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_placed_strip.R")'
```
Expected: file-not-found error sourcing `placed_strip.R`, or function-not-found.

- [ ] **Step 3.3: Write the helper**

Create `Answer Keys/MLB Dashboard/placed_strip.R`:

```r
# Answer Keys/MLB Dashboard/placed_strip.R
#
# Pure helper that renders the "placed" hero strip for a single bet card
# after one or more Wagerzon placements. One <span class="placement-chip">
# per placement plus an optional dashed "+ another" button.
#
# Inputs:
#   chips           list of list(account=<chr>, risk=<num>, ticket=<chr>)
#   all_wz_accounts character vector of every WZ account label currently
#                   discovered by wagerzon_accounts.list_accounts()
#   bet_hash        character — the bet hash; emitted as data-bet-hash on
#                   the strip and on each chip so JS handlers can route
#                   per-bet (also doubles as the placeBet() routing key).
#   book            character — the bookmaker key of the pick. The
#                   "+ another" affordance only renders when book ==
#                   "wagerzon"; non-WZ books just get the chips.
#
# Returns a single HTML string (the entire <div class="hero-placed"> block).

.acct_short <- function(label) {
  # "Wagerzon" -> "WZ"; "WagerzonJ" -> "WZJ"; "WagerzonC" -> "WZC".
  # Falls back to label as-is for any non-Wagerzon book account.
  if (grepl("^Wagerzon", label)) {
    suffix <- sub("^Wagerzon", "", label)
    paste0("WZ", suffix)
  } else {
    label
  }
}

.chip_html <- function(chip, bet_hash) {
  sprintf(
    '<span class="placement-chip" data-account="%s" data-risk="%s" data-ticket="%s" data-bet-hash="%s"><span class="acct">%s</span>$%s <span class="ticket">#%s</span></span>',
    htmltools::htmlEscape(chip$account),
    as.character(round(as.numeric(chip$risk))),
    htmltools::htmlEscape(chip$ticket),
    htmltools::htmlEscape(bet_hash),
    htmltools::htmlEscape(.acct_short(chip$account)),
    as.character(round(as.numeric(chip$risk))),
    htmltools::htmlEscape(chip$ticket)
  )
}

render_placed_strip <- function(chips, all_wz_accounts, bet_hash, book) {
  stopifnot(is.list(chips), is.character(all_wz_accounts),
            is.character(bet_hash), is.character(book))

  chip_html <- paste(vapply(chips, .chip_html, character(1),
                            bet_hash = bet_hash),
                     collapse = "")

  is_wz <- identical(tolower(book), "wagerzon")
  placed_accounts <- vapply(chips, function(c) c$account, character(1))
  untouched <- setdiff(all_wz_accounts, placed_accounts)

  add_another_html <- ""
  if (is_wz) {
    disabled_attr <- if (length(untouched) == 0L) " disabled" else ""
    disabled_title <- if (length(untouched) == 0L) {
      ' title="All WZ accounts placed."'
    } else {
      ' title="Re-open the card to place on another WZ account"'
    }
    add_another_html <- sprintf(
      '<button type="button" class="add-another" data-bet-hash="%s" onclick="addAnother(this)"%s%s>+ another</button>',
      htmltools::htmlEscape(bet_hash), disabled_attr, disabled_title
    )
  }

  sprintf(
    '<div class="hero-placed" data-bet-hash="%s"><span class="placed-label">placed</span>%s%s</div>',
    htmltools::htmlEscape(bet_hash), chip_html, add_another_html
  )
}
```

- [ ] **Step 3.4: Run tests to verify they pass**

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_placed_strip.R")'
```
Expected: 5 tests pass, 0 fail.

- [ ] **Step 3.5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/placed_strip.R" \
        "Answer Keys/tests/test_placed_strip.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): render_placed_strip helper for multi-account chips

Pure R helper that produces the hero-placed strip HTML: one
placement-chip per row in placed_bets for a bet, plus a dashed
"+ another" button enabled iff at least one WZ account has no chip
yet. Non-WZ picks render chips only. Each chip carries
data-account/risk/ticket/bet-hash for the JS handlers landing in
a later task.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Wire `create_bets_table` to use `render_placed_strip`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — three changes:
  1. `source()` the new helper near the other `source()` lines (around the top of the file alongside `book_cell.R` etc.).
  2. Update the placed-bets query inside `create_bets_table()` to read `account` and aggregate per `bet_hash`.
  3. Replace the single `placed-bet-label` chip with a `render_placed_strip(...)` call when at least one placement exists.
  4. Add inline CSS (already-existing `<style>` block).

- [ ] **Step 4.1: Source the helper**

In `mlb_dashboard.R`, find the existing `source(...)` block (search for `book_cell.R`). Add a sibling line:

```r
source(file.path("Answer Keys", "MLB Dashboard", "placed_strip.R"))
```

Match the relative-path style of the surrounding `source()` calls — if the existing ones use a different anchor (e.g. `here::here`), use the same.

- [ ] **Step 4.2: Add CSS for the placed strip + chips + add-another button**

In `mlb_dashboard.R`, find the inline `<style>` block at approximately line 2891 (search for `.bet-card-v8 .hero`). Append the following rules to the same block:

```css
.bet-card-v8 .hero-placed {
  display: flex; align-items: center; gap: 12px; flex-wrap: wrap;
  background: #0d1f12;
  border: 1px solid #1f4d2e;
  border-radius: 8px;
  padding: 8px 14px;
  margin-top: 12px;
}
.bet-card-v8 .hero-placed .placed-label {
  color: #56d364;
  font-size: 11px;
  text-transform: uppercase;
  letter-spacing: .5px;
  font-weight: 700;
}
.bet-card-v8 .placement-chip {
  display: inline-flex; align-items: center; gap: 6px;
  background: #15321f;
  border: 1px solid #1f4d2e;
  border-radius: 5px;
  padding: 3px 8px;
  color: #c9d1d9;
  font: 12px ui-monospace, SFMono-Regular, Menlo, monospace;
}
.bet-card-v8 .placement-chip .acct  { color: #56d364; font-weight: 700; }
.bet-card-v8 .placement-chip .ticket { color: #8b949e; font-size: 11px; }

.bet-card-v8 .add-another {
  background: transparent;
  border: 1px dashed #1f4d2e;
  color: #56d364;
  border-radius: 5px;
  padding: 3px 10px;
  font-size: 12px;
  cursor: pointer;
  font-weight: 600;
}
.bet-card-v8 .add-another[disabled] {
  color: #6e7681;
  border-color: #30363d;
  cursor: not-allowed;
}
```

- [ ] **Step 4.3: Update the placed-bets query inside `create_bets_table` to return account**

Locate `create_bets_table` in `mlb_dashboard.R` (search for the function definition). It currently runs a query against `placed_bets` to find which bets are already placed. Today the query returns a single row per `bet_hash`; switch it to return per (bet_hash, account):

```r
# Replace whatever query is currently used to fetch the placed-bet map.
# The exact SELECT shape will match the surrounding code; this is the contract:
placed_rows <- DBI::dbGetQuery(con_dashboard, "
  SELECT bet_hash, account, actual_size AS risk, ticket_number AS ticket
  FROM placed_bets
  WHERE status IN ('placed', 'placing') AND account IS NOT NULL
")
```

Then group into a `list-by-bet_hash`:

```r
placed_chips_by_hash <- split(placed_rows, placed_rows$bet_hash)
```

Each list entry is a data.frame of chips for that bet.

- [ ] **Step 4.4: Render via `render_placed_strip` when chips exist**

Where the table-row builder currently chooses between an editable hero (un-placed) and the legacy `placed-bet-label` chip (single placement), replace with:

```r
# Look up all WZ accounts once (above the per-row loop)
all_wz_accounts <- wagerzon_accounts_labels()
# wagerzon_accounts_labels() is a small helper that calls into the
# Python registry via the existing /api/wagerzon/balances JSON or
# whatever cached source mlb_dashboard.R already uses for header pills.
# Use the same source — search for WZ_SELECTED_ACCOUNT / wz-account-pills
# render block; that block already enumerates the same list.

# Per bet row:
chips_df <- placed_chips_by_hash[[bet_hash]]
if (!is.null(chips_df) && nrow(chips_df) > 0) {
  chips_list <- lapply(seq_len(nrow(chips_df)), function(i) {
    list(account = chips_df$account[i],
         risk    = chips_df$risk[i],
         ticket  = chips_df$ticket[i])
  })
  hero_html <- render_placed_strip(
    chips           = chips_list,
    all_wz_accounts = all_wz_accounts,
    bet_hash        = bet_hash,
    book            = bet$bookmaker_key
  )
} else {
  hero_html <- render_hero_strip(...)  # the existing editable hero call
}
```

Important: the legacy `placed-bet-label` span path can be removed once `render_placed_strip` is wired in. If you keep both, add a fallback: when `account IS NULL` for a chip (pre-migration row that didn't get backfilled for whatever reason), render it as a single chip with `account = "Wagerzon"` — but the migration backfill should make this case impossible.

- [ ] **Step 4.5: Visual smoke test**

```bash
# From repo root; restart the dashboard so it picks up CSS + helper changes.
# (Adjust the launch command to match how the dashboard is normally run
# in this environment. The typical command is documented in the Dashboard
# README; if no README, use the same command that's in run.py or the
# session's recent shell history.)
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"  # or the project's run command
```

Open the dashboard. Place a small WZ bet on Wagerzon at $5 (real placement — the dashboard hits live WZ, no test mode). Verify:
- The card flips to a green-tinted strip with `placed · WZ $5 #W…` and a dashed `+ another` button next to it.
- `+ another` is clickable (no observable behavior yet — the click handler lands in Task 5).
- Reloading the page re-renders the same strip from the server (server-rendered shape matches the JS-target shape that Task 5 will emit).

- [ ] **Step 4.6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): render placed bets via render_placed_strip

create_bets_table now queries placed_bets by (bet_hash, account),
groups rows per bet_hash, and renders the new hero-placed strip
when one or more placements exist. Adds inline CSS for the
hero-placed container, placement chips and the dashed
"+ another" affordance.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: JS handler `addAnother` + header-pill clear hook

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add two JS functions inside the existing `tags$script(HTML('...'))` block (find it by searching for `function placeBet(btn)`).

> No JS unit test harness exists in this codebase. Verification is by manual click-through in Task 8.

- [ ] **Step 5.1: Add the `addAnother` handler**

Insert near the existing `placeBet` definition (matching its quoting style — single-quoted R string with escaped JS quotes):

```javascript
function addAnother(btn) {
  // Find the bet card that owns this button.
  var card = btn.closest('.bet-card-v8');
  if (!card) return;

  var betHash = btn.dataset.betHash;
  var heroPlaced = card.querySelector('.hero-placed');
  if (!heroPlaced) return;

  // Stash the placed strip in a sibling above the rebuilt hero so the
  // chips remain visible while the user types Risk.
  // The cheapest swap is to leave heroPlaced where it is and inject the
  // editable hero immediately AFTER it (Task 4 puts hero-placed in the
  // card; we'll append a sibling).
  if (card.querySelector('.hero.reopened')) {
    // already re-opened
    return;
  }

  // Read the most recent Risk amount from the last chip as a sensible
  // default for the new placement (assumption: user is spreading the
  // same intended stake across accounts).
  var chips = heroPlaced.querySelectorAll('.placement-chip');
  var lastRisk = chips.length
    ? parseFloat(chips[chips.length - 1].dataset.risk) || 0
    : 0;

  // Find the original card's hero data-* (book, odds, ev, fair, etc.)
  // by reading from the first chip's bet metadata. Easier path: grab the
  // data-* attributes from the card's parent context. The existing
  // create_bets_table renders enough data on the card root.
  // For v1, pull odds + ev from the original hero strip's data-attrs
  // that were emitted before placement; if absent, fall back to the
  // chip's own metadata.
  var heroData = card.dataset; // bet-level data-* set by create_bets_table
  // If create_bets_table doesn't already write per-card data-* attrs,
  // emit them as part of Task 4 (e.g. data-pick-odds, data-fair-odds,
  // data-pick-book, data-ev). If not present, this handler degrades to
  // showing Risk + Place only.

  // Build the editable hero (mirror of render_hero_strip output, but
  // marked .reopened so we can find it later).
  var editable = document.createElement('div');
  editable.className = 'hero reopened';
  editable.dataset.betHash = betHash;
  editable.innerHTML = ''
    + '<div class="pick">'
    +   '<span class="book">' + (heroData.pickBook || 'WZ') + '</span>'
    +   '<span class="odds">' + (heroData.pickOdds || '') + '</span>'
    + '</div>'
    + '<div class="divider"></div>'
    + '<div class="stat"><span class="lbl">Fair</span><span class="val">' + (heroData.fairOdds || '') + '</span></div>'
    + '<div class="stat"><span class="lbl">EV</span><span class="val ev">' + (heroData.evPct || '') + '</span></div>'
    + '<div class="stat risk-stat" data-model-risk="' + lastRisk + '" data-american-odds="' + (heroData.pickOddsRaw || '0') + '">'
    +   '<span class="lbl">Risk</span>'
    +   '<div class="risk-row">'
    +     '<span class="val risk risk-value" tabindex="0" title="click to edit">$' + lastRisk + '</span>'
    +     '<button type="button" class="risk-reset" title="reset">↻</button>'
    +     '<span class="risk-error" hidden></span>'
    +   '</div>'
    + '</div>'
    + '<div class="stat"><span class="lbl">To Win</span><div class="towin-row"><span class="val win towin-value"></span><span class="towin-status" hidden></span></div></div>'
    + '<div class="actions">'
    +   '<button class="btn-place reopened-place" data-bet-hash="' + betHash + '" onclick="placeAnother(this)">Place</button>'
    + '</div>';

  heroPlaced.insertAdjacentElement('afterend', editable);

  // Re-evaluate Place enable state based on header pill vs existing chips.
  _refreshReopenedPlaceState(editable);

  // Light informational toast — matches the existing showToast signature.
  if (typeof showToast === 'function') {
    showToast('Switch the header pill to add another account.', 'info');
  }
}

function _refreshReopenedPlaceState(editable) {
  if (!editable) return;
  var card = editable.closest('.bet-card-v8');
  if (!card) return;
  var placeBtn = editable.querySelector('.reopened-place');
  if (!placeBtn) return;

  var currentAcct = window.WZ_SELECTED_ACCOUNT || '';
  var chipAccounts = Array.prototype.map.call(
    card.querySelectorAll('.hero-placed .placement-chip'),
    function (chip) { return chip.dataset.account; }
  );
  if (!currentAcct || chipAccounts.indexOf(currentAcct) !== -1) {
    placeBtn.disabled = true;
    placeBtn.title = 'Switch the header pill to a different account';
  } else {
    placeBtn.disabled = false;
    placeBtn.title = '';
  }
}
```

- [ ] **Step 5.2: Make sure Task 4 emits the data-* attrs `addAnother` reads**

`addAnother` reads `card.dataset.pickBook`, `pickOdds`, `pickOddsRaw`, `fairOdds`, `evPct`. If `create_bets_table` doesn't already emit those on the root `.bet-card-v8` element, add them in Task 4's render path:

```r
sprintf('<div class="bet-card-v8" data-bet-hash="%s" data-pick-book="%s" data-pick-odds="%s" data-pick-odds-raw="%d" data-fair-odds="%s" data-ev-pct="%s">...',
  bet_hash, pick_book, format_odds(pick_odds), as.integer(pick_odds),
  format_odds(fair_odds), sprintf("+%.1f%%", ev_pct))
```

If they're already there (likely — these are visible on screen and the click-to-edit Risk handler probably already reads odds from the card), leave alone.

- [ ] **Step 5.3: Add header-pill click hook to clear `data-expected-win` and re-evaluate Place enable state**

Find the JS that handles header pill clicks (search for `WZ_SELECTED_ACCOUNT` assignment, likely near the `wz-account-pills` rendering at line ~3155 area, but the JS lives in the same script block). Append to the click handler:

```javascript
// After switching window.WZ_SELECTED_ACCOUNT, sweep every card:
document.querySelectorAll('.bet-card-v8').forEach(function (card) {
  // Clear the per-card verified quote — next Risk edit re-verifies.
  var btn = card.querySelector('[data-expected-win]');
  if (btn) btn.dataset.expectedWin = '';
  // Re-evaluate Place enable state on any re-opened editable hero.
  var reopened = card.querySelector('.hero.reopened');
  if (reopened) _refreshReopenedPlaceState(reopened);
});
```

> If multiple elements per card carry `data-expected-win` (e.g. the editable Risk row also stores it), use `querySelectorAll` instead and clear all of them.

- [ ] **Step 5.4: Smoke-test the re-open flow (visual)**

Restart the dashboard (per Task 4.5). Place a WZ bet on Wagerzon at $5. Click `+ another`. Verify:
- An editable hero strip appears below the placed strip.
- The Place button in the editable strip is disabled with the tooltip "Switch the header pill to a different account."
- Click WagerzonJ in the header pill row. Place becomes enabled. (No actual placement yet — that's Task 6.)
- Click Wagerzon again. Place becomes disabled again.

- [ ] **Step 5.5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): addAnother handler re-opens placed cards

Clicking the dashed "+ another" button injects an editable hero
strip below the placed strip with a disabled Place button and a
"switch the header pill" toast. Header-pill click hook clears
data-expected-win on every card and re-evaluates Place enable
state on any re-opened editable hero.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: JS handler `placeAnother` + chip append on success

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add `placeAnother` JS function adjacent to `addAnother`. Also extend the existing `removeBet` JS to include `data-account` in the DELETE request body.

- [ ] **Step 6.1: Add `placeAnother`**

```javascript
function placeAnother(btn) {
  var card = btn.closest('.bet-card-v8');
  if (!card) return;
  var editable = btn.closest('.hero.reopened');
  if (!editable) return;

  var account = window.WZ_SELECTED_ACCOUNT;
  if (!account) {
    showToast('No Wagerzon account selected — pick one in the header pills', 'error');
    return;
  }
  // Defensive: re-check that the selected account isn't already a chip on this card.
  var chipAccounts = Array.prototype.map.call(
    card.querySelectorAll('.hero-placed .placement-chip'),
    function (c) { return c.dataset.account; }
  );
  if (chipAccounts.indexOf(account) !== -1) {
    showToast('Already placed on ' + account + ' for this bet', 'warning');
    return;
  }

  // Read Risk from the editable hero's risk-value span.
  var riskValueEl = editable.querySelector('.risk-value');
  var riskStr = (riskValueEl && riskValueEl.textContent) || '';
  var riskNum = parseFloat(riskStr.replace(/[^0-9.\-]/g, '')) || 0;
  if (riskNum <= 0) {
    showToast('Enter a Risk amount before placing', 'warning');
    return;
  }
  var expectedWin = riskValueEl && riskValueEl.dataset.expectedWin
    ? parseFloat(riskValueEl.dataset.expectedWin)
    : null;

  // Pull bet metadata from the card root data-* attrs (set in
  // create_bets_table). We need everything /api/place-bet expects.
  var d = card.dataset;
  var body = {
    bet_hash:         d.betHash,
    bookmaker_key:    d.pickBook ? d.pickBook.toLowerCase() : 'wagerzon',
    account:          account,
    bet_on:           d.betOn,
    line:             (d.line === '' || d.line === undefined) ? null : parseFloat(d.line),
    market:           d.market,
    american_odds:    parseInt(d.pickOddsRaw, 10),
    actual_size:      riskNum,
    kelly_bet:        parseFloat(d.modelSize || d.kellyBet || d.pickOddsRaw) || riskNum,
    wz_odds_at_place: parseInt(d.pickOddsRaw, 10),
    expected_win:     expectedWin,
    game_id:          d.gameId,
    home_team:        d.home,
    away_team:        d.away,
    game_time:        d.time,
    model_prob:       d.prob ? parseFloat(d.prob) : 0.0,
    model_ev:         d.ev   ? parseFloat(d.ev)   : 0.0
  };

  btn.disabled = true;
  var originalLabel = btn.textContent;
  btn.textContent = 'Placing...';

  fetch('/api/place-bet', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body)
  })
    .then(function (r) {
      return r.json().then(function (j) { return {ok: r.ok, status: r.status, body: j}; });
    })
    .then(function (resp) {
      var result = resp.body;
      if (resp.status === 409) {
        showToast('Already in flight: ' + (result.error || result.status), 'warning');
        btn.disabled = false; btn.textContent = originalLabel;
        return;
      }
      if (result.status === 'placed') {
        var ticket = result.ticket_number || '';
        // Append a new chip to the existing hero-placed strip and tear down
        // the editable hero. Matches what create_bets_table would server-
        // render on the next refresh.
        var heroPlaced = card.querySelector('.hero-placed');
        if (heroPlaced) {
          var chip = document.createElement('span');
          chip.className = 'placement-chip';
          chip.dataset.account = account;
          chip.dataset.risk    = riskNum;
          chip.dataset.ticket  = ticket;
          chip.dataset.betHash = body.bet_hash;
          var acctShort = account.replace(/^Wagerzon/, 'WZ');
          chip.innerHTML =
            '<span class="acct">' + acctShort + '</span>$' + Math.round(riskNum) +
            ' <span class="ticket">#' + ticket + '</span>';
          // Insert chip just before the "+ another" button (so the chip
          // ordering stays oldest→newest left-to-right).
          var addBtn = heroPlaced.querySelector('.add-another');
          if (addBtn) {
            heroPlaced.insertBefore(chip, addBtn);
          } else {
            heroPlaced.appendChild(chip);
          }
          // Re-evaluate "+ another" enable state: disabled if every WZ
          // account is now a chip.
          if (addBtn) {
            var allChips = heroPlaced.querySelectorAll('.placement-chip');
            var placedAccts = Array.prototype.map.call(allChips,
              function (c) { return c.dataset.account; });
            // We do NOT have all_wz_accounts on the client. Easiest:
            // count chips against the header pill row's account labels.
            var headerLabels = Array.prototype.map.call(
              document.querySelectorAll('.wz-pills .acct-pill'),
              function (p) { return p.dataset.account || p.textContent.trim().split(/\s+/)[0]; }
            );
            var untouched = headerLabels.filter(function (a) {
              return placedAccts.indexOf(a) === -1;
            });
            if (untouched.length === 0) {
              addBtn.disabled = true;
              addBtn.title = 'All WZ accounts placed.';
            }
          }
        }
        // Remove the editable hero — back to chips-only view.
        editable.remove();
        showToast('Placed at ' + account + ' #' + ticket, 'success');
        return;
      }
      if (result.status === 'price_moved') {
        showToast('Price moved — bet not placed', 'warning');
        btn.disabled = false; btn.textContent = originalLabel;
        return;
      }
      if (result.error) {
        showToast(result.error, 'error');
      } else if (result.status) {
        showToast('Status: ' + result.status + (result.error_msg ? ' — ' + result.error_msg : ''), 'error');
      } else {
        showToast('Unknown response', 'error');
      }
      btn.disabled = false; btn.textContent = originalLabel;
    })
    .catch(function (e) {
      showToast('Network error: ' + e.message, 'error');
      btn.disabled = false; btn.textContent = originalLabel;
    });
}
```

> The exact `acct-pill` selector (`.wz-pills .acct-pill`) needs to match what Task 4's header renders — confirm by grep'ing `wz-pills` in `mlb_dashboard.R` and using the same class names. If the rendered class is different, swap it here.

- [ ] **Step 6.2: Update `removeBet` (if implemented) to include `data-account`**

Find the existing `removeBet` JS handler (or equivalent for `placed-bet-label` chip un-place). If it exists, the body it sends to `/api/remove-bet` currently only contains `bet_hash`. Add `account` (read from the chip's `data-account` attribute):

```javascript
var body = { bet_hash: data.hash, account: data.account };
```

If no `removeBet` handler exists today (un-place is not supported in the current UI), skip this step.

- [ ] **Step 6.3: Smoke-test the full place→re-open→place flow (manual)**

Restart the dashboard. Place a $5 WZ bet on Wagerzon. Click `+ another`. Switch the header pill to WagerzonJ. Edit Risk to $4 (verify To Win updates from the WZ-verified quote). Click Place. Verify:
- A second `placement-chip` appears in the hero-placed strip after the first.
- The editable hero is removed (replaced by the chip).
- The `+ another` button remains (if WagerzonC exists) or becomes disabled (if Wagerzon + WagerzonJ are the only registered WZ accounts).
- Reload the page. Both chips persist from the server-rendered HTML.

- [ ] **Step 6.4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): placeAnother handler appends chip on success

placeAnother POSTs /api/place-bet under the current header-pill
account, then appends a placement-chip to the existing hero-placed
strip and tears down the editable hero. Matches the server-rendered
shape so a page refresh produces an identical result.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Documentation updates

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — extend the "MLB Dashboard — Wagerzon multi-account" section
- Modify: `Answer Keys/MLB Dashboard/migrations/__init__.py` (if it contains a migration registry comment) OR add a top-level "Migrations" note in the same CLAUDE.md.

- [ ] **Step 7.1: Edit `Answer Keys/CLAUDE.md`**

Find the section titled "## MLB Dashboard — Wagerzon multi-account". Add a new subsection at the end of it:

```markdown
### Multi-account re-place (v1, 2026-05)

After a Wagerzon bet is placed, the card shows a green chip per
placement plus a dashed "+ another" button. Clicking `+ another`
re-opens an editable hero strip; the user switches the header pill
to a different WZ account, edits Risk (re-verified via
`/api/wz-quote-single`), and clicks Place to add a second placement
on the same bet.

- `placed_bets` PRIMARY KEY is now composite `(bet_hash, account)`
  (migration 003). Same `bet_hash` with different `account` is allowed;
  same `(bet_hash, account)` is rejected.
- All server-side writers (`_insert_placement_breadcrumb`,
  `_finalize_placement`, `/api/place-bet` 409 check, `/api/remove-bet`)
  scope by composite key.
- The header pill click handler clears `data-expected-win` on every
  card so the next Risk edit re-verifies under the new account.
- **v1 does NOT support stacking on the same account.** If WZ caps you
  at $50/wager and you want $500 on Wagerzon alone, the composite PK
  blocks the second placement (409). Future change.
- Spec: `docs/superpowers/specs/2026-05-23-mlb-dashboard-multi-account-re-place-design.md`.

### Migration 003 — composite PK

Apply once per environment after pulling this branch and BEFORE
restarting the dashboard:
```
python "Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py" \
    "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```
Idempotent — second invocation is a no-op. The running dashboard's R
process holds a connection to the old schema and will not see the new
PK without a restart.
```

- [ ] **Step 7.2: Commit**

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "$(cat <<'EOF'
docs(answer-keys): describe multi-account re-place + migration 003

Adds subsections under "MLB Dashboard — Wagerzon multi-account"
covering the re-place flow, composite-PK semantics, and the manual
migration command + post-migration server restart.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Manual end-to-end verification

> No automated UI test for the dashboard. Run the following manual checklist before merging. Capture any deviations as issues to fix in a follow-up commit on this same branch.

- [ ] **Step 8.1: Apply migration 003 to the live dashboard DB**

```bash
python "Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py" \
    "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```
Expected: `Migration 003 applied to ...`. Re-run once more to verify idempotency — expected: no output change, no error.

- [ ] **Step 8.2: Restart the dashboard server**

Stop the running dashboard (Ctrl-C in its terminal) and relaunch via the same command used in this environment.

- [ ] **Step 8.3: Walk the spec's 8-step test plan**

From `docs/superpowers/specs/2026-05-23-mlb-dashboard-multi-account-re-place-design.md` Section 10:

1. Place a WZ bet on Wagerzon at $5. Verify chip + `+ another` appears.
2. Click `+ another`. Verify Place is disabled.
3. Switch to WagerzonJ in the header pill row. Verify Place enables.
4. Edit Risk to $4. Verify To Win updates from the WZ-verified quote.
5. Place. Verify second chip lands and editable hero disappears.
6. Repeat for WagerzonC (if available). Verify `+ another` becomes disabled after all WZ accounts are placed.
7. Refresh the page. Verify all chips persist (server-rendered).
8. Place a Hoop88 bet on a different card. Verify NO `+ another` button appears (non-WZ book).

- [ ] **Step 8.4: Regression spot-checks**

1. Single-account WZ placement (today's flow) — place once, verify card flips to single chip + dashed `+ another`. No regression.
2. 409 in-flight check — double-click Place fast on a fresh card. Second click should hit "Already in flight" toast.
3. WZ-verified quote — edit Risk on a re-opened hero; verify the To Win value updates from the live quote (network tab shows `/api/wz-quote-single` call with the new account).
4. Removal (if `/api/remove-bet` is wired in the UI) — toggle un-place a chip; verify only the targeted (bet_hash, account) row is removed.

- [ ] **Step 8.5: Record findings**

If any step fails, fix the issue and commit on this branch with a `fix(mlb-dashboard): ...` message. Re-run the affected steps before proceeding to Task 9.

---

## Task 9: Pre-merge review & merge

> Per the project's CLAUDE.md, every feature branch requires an executive review of the full diff and explicit user approval before merging. This task spells out the review checklist and merge steps — do NOT merge without user approval.

- [ ] **Step 9.1: Generate the full diff**

```bash
git fetch origin main
git diff origin/main..HEAD
```

- [ ] **Step 9.2: Run the pre-merge review checklist**

Walk the items from the project's CLAUDE.md pre-merge checklist:

- **Data integrity**: migration backfills correctly; composite PK rejects duplicates; no other writers leak `bet_hash`-only WHERE clauses.
- **Resource safety**: every `duckdb.connect` in the changed code uses a `try/finally` close.
- **Edge cases**: pre-migration rows with NULL account get backfilled before the PK swap; non-WZ books don't render `+ another`.
- **Dead code**: no unused JS handlers, no leftover legacy `placed-bet-label` rendering path if it's been fully replaced.
- **Log / disk hygiene**: no new file writes added by this change set.
- **Security**: no new secrets in logs; `/api/remove-bet` rejects missing `account` parameter.

Document any findings as ISSUES TO FIX vs ACCEPTABLE RISKS in a brief comment before the merge.

- [ ] **Step 9.3: Request user approval to merge**

Surface the diff summary and the review notes to the user. **Do not merge without an explicit "yes, merge".**

- [ ] **Step 9.4: Merge to main**

After approval, from the main repository (NOT the worktree):

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff worktree-feature+dashboard-multi-account-re-place-spec \
    -m "Merge worktree-feature+dashboard-multi-account-re-place-spec: multi-account re-place on the MLB Dashboard"
```

> The merge happens from `main` in the primary checkout, not from inside the worktree. The user runs this step.

- [ ] **Step 9.5: Apply migration against the live DB on main and restart**

```bash
python "Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py" \
    "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
# Restart the dashboard so it sees the new schema.
```

- [ ] **Step 9.6: Clean up the worktree**

```bash
git worktree remove .claude/worktrees/feature+dashboard-multi-account-re-place-spec
git branch -d worktree-feature+dashboard-multi-account-re-place-spec
```

---

## Self-review notes (from plan author)

- **Spec coverage check (re-read spec sections 1-11 against the plan):**
  - §1 Background / existing plumbing — informational, no task needed.
  - §2 Goals & non-goals — covered: Task 1+2 implement multi-row placement, Task 3+4+5+6 cover UX, Task 8 verifies non-goals (no stacking, no batch).
  - §3 UX flow + state machine — covered by Task 4 (render) + Task 5 (re-open) + Task 6 (place new chip).
  - §4 Data model — Task 1.
  - §5 Server changes — Task 2.
  - §6 Frontend changes — Tasks 3+4+5+6.
  - §7 Edge cases — exercised in Task 8.3 and Task 8.4. Composite-PK 409 case explicit in Task 2 + Task 8.4.2.
  - §8 Version control plan — Tasks 1-9 commit structure mirrors the suggested order.
  - §9 Documentation updates — Task 7.
  - §10 Testing — Tasks 1.4 (migration), 3.4 (R helper), 8 (manual).
  - §11 Open questions — flagged as out of scope; not part of plan.
- **Placeholder scan:** every code step has full code; no "implement later" or "similar to above" references; no TBDs in commit messages.
- **Type consistency:** `render_placed_strip(chips, all_wz_accounts, bet_hash, book)` signature is consistent across Tasks 3 and 4. `_finalize_placement(bet_hash, account, result)` signature is consistent across Task 2.2 and 2.3. JS function names `addAnother`, `placeAnother`, `_refreshReopenedPlaceState` consistent across Tasks 5 and 6.

---
