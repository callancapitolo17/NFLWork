# MLB Dashboard Branch-Testing Isolation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a git worktree a self-contained place to fully test MLB dashboard changes before merging, where the worktree reads/writes its own seeded copy of the data and `main` stays byte-identical.

**Architecture:** The codebase is already ~90% script-relative. We add a repo-root resolver to `Tools.R` (derives root from the running script, falls back to `~/NFLWork`), route the nine `get_*_odds` readers through it, remove three hardcoded "jump back to main" overrides (`MLB.R` `setwd`, `mlb_dashboard.R` strip-hack, server strip-hack), and add a `seed_test_data.sh` that CoW-clones live DBs into the worktree. A guard test proves no DB path escapes the worktree (and that `main` is unchanged).

**Tech Stack:** R (testthat), Python (Flask), bash, DuckDB, macOS APFS `cp -c` (copy-on-write).

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `Answer Keys/Tools.R` | Repo-root resolver + 9 reader path defaults | Modify |
| `Answer Keys/MLB Answer Key/MLB.R` | Script-relative `setwd` + set repo root | Modify |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Remove worktree strip-hack | Modify |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | Remove strip-hack + configurable port | Modify |
| `Answer Keys/MLB Dashboard/seed_test_data.sh` | CoW-clone live DBs into the worktree | Create |
| `Answer Keys/MLB Dashboard/test_dashboard.sh` | Wrapper: seed → run slice → serve on test port | Create |
| `Answer Keys/tests/test_nflwork_root.R` | Guard test: routing + no-regression | Create |
| `.gitignore` | Ensure DuckDB WAL is ignored | Modify |

**Pre-flight:** Confirm you are in the worktree and on the feature branch.

- [ ] **Step 0: Verify branch + worktree**

Run: `git rev-parse --show-toplevel && git branch --show-current`
Expected: a path under `.claude/worktrees/` and branch `worktree-mlb-dashboard-worktree-testing`.

---

## Task 1: Repo-root resolver in Tools.R (TDD)

Add three functions so worktree-resident code can point at the worktree, with a `~/NFLWork` fallback that keeps every existing caller (CBB/NFL) unchanged. `derive_repo_root()` takes an optional `script_path` so it is unit-testable.

**Files:**
- Modify: `Answer Keys/Tools.R` (top of file, after the `library()` block)
- Test: `Answer Keys/tests/test_nflwork_root.R`

- [ ] **Step 1: Write the failing test**

Create `Answer Keys/tests/test_nflwork_root.R`:

```r
# Answer Keys/tests/test_nflwork_root.R
library(testthat)
source("../Tools.R")

test_that("nflwork_root falls back to ~/NFLWork when unset", {
  set_nflwork_root(NULL)
  expect_equal(nflwork_root(), normalizePath(path.expand("~/NFLWork"), mustWork = FALSE))
})

test_that("set_nflwork_root overrides the root", {
  tmp <- tempfile("root_"); dir.create(tmp)
  set_nflwork_root(tmp)
  expect_equal(nflwork_root(), normalizePath(tmp, mustWork = FALSE))
  set_nflwork_root(NULL)  # reset for other tests
})

test_that("derive_repo_root walks up to the dir containing 'Answer Keys'", {
  root <- tempfile("repo_"); dir.create(root)
  deep <- file.path(root, "Answer Keys", "MLB Answer Key")
  dir.create(deep, recursive = TRUE)
  script <- file.path(deep, "MLB.R")
  expect_equal(
    derive_repo_root(script_path = script),
    normalizePath(root, mustWork = FALSE)
  )
})

test_that("derive_repo_root falls back to ~/NFLWork when no marker found", {
  orphan <- file.path(tempfile("orphan_"))
  dir.create(orphan, recursive = TRUE)
  expect_equal(
    derive_repo_root(script_path = file.path(orphan, "x.R")),
    normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_nflwork_root.R")'`
Expected: FAIL — `could not find function "set_nflwork_root"`.

- [ ] **Step 3: Add the resolver to Tools.R**

In `Answer Keys/Tools.R`, immediately after the `library()`/`suppressPackageStartupMessages` block at the top, insert:

```r
# ---------------------------------------------------------------------------
# Repo-root resolution
# Lets worktree-resident code point at the worktree instead of the main repo.
# Entry scripts call set_nflwork_root(derive_repo_root()); any caller that does
# NOT opt in falls back to ~/NFLWork, so CBB/NFL behaviour is unchanged and
# `main` stays byte-identical.
# ---------------------------------------------------------------------------
.NFLWORK_ROOT <- NULL

set_nflwork_root <- function(root) {
  if (is.null(root)) {
    .NFLWORK_ROOT <<- NULL
  } else {
    .NFLWORK_ROOT <<- normalizePath(root, mustWork = FALSE)
  }
}

nflwork_root <- function() {
  if (!is.null(.NFLWORK_ROOT)) return(.NFLWORK_ROOT)
  normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)
}

# Derive the repo root from the running Rscript's --file= path (or an explicit
# script_path for testing) by walking up until we find a dir that contains an
# "Answer Keys" subdirectory. Falls back to ~/NFLWork if no marker is found.
derive_repo_root <- function(script_path = NULL) {
  if (is.null(script_path)) {
    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    if (length(fa) == 0) {
      return(normalizePath(path.expand("~/NFLWork"), mustWork = FALSE))
    }
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", fa[1]))
  }
  d <- normalizePath(dirname(script_path), mustWork = FALSE)
  while (d != dirname(d)) {
    if (dir.exists(file.path(d, "Answer Keys"))) {
      return(normalizePath(d, mustWork = FALSE))
    }
    d <- dirname(d)
  }
  normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_nflwork_root.R")'`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_nflwork_root.R"
git commit -m "feat(mlb): add repo-root resolver to Tools.R (fallback ~/NFLWork)"
```

---

## Task 2: Route the nine readers + team-dict through nflwork_root() (TDD)

Change each `get_*_odds` default and the team-dict path from a hardcoded `~/NFLWork/...` to `file.path(nflwork_root(), ...)`. R evaluates default args lazily at call time, so `nflwork_root()` resolves per call. The readers already `warning(... not found at %s ...)` with the resolved path when the DB is missing — we assert on that to prove routing.

**Files:**
- Modify: `Answer Keys/Tools.R` (lines for each reader default + `dict_db`)
- Test: `Answer Keys/tests/test_nflwork_root.R` (extend)

- [ ] **Step 1: Add the failing routing test**

Append to `Answer Keys/tests/test_nflwork_root.R`:

```r
test_that("get_*_odds defaults route through nflwork_root()", {
  # Inspect each function's db_path DEFAULT source — deterministic, no dependency
  # on DB existence or warning text. (nflwork_root() correctness is proven above.)
  book_dirs <- list(
    get_wagerzon_odds  = "wagerzon_odds",
    get_hoop88_odds    = "hoop88_odds",
    get_bfa_odds       = "bfa_odds",
    get_bookmaker_odds = "bookmaker_odds",
    get_dk_odds        = "dk_odds",
    get_fd_odds        = "fd_odds",
    get_bet105_odds    = "bet105_odds",
    get_kalshi_odds    = "kalshi_odds"
  )
  for (fn in names(book_dirs)) {
    default_src <- paste(deparse(formals(get(fn))$db_path), collapse = " ")
    expect_match(default_src, "nflwork_root()", fixed = TRUE, info = fn)
    expect_match(default_src, book_dirs[[fn]], fixed = TRUE, info = fn)
  }
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_nflwork_root.R")'`
Expected: FAIL — the defaults are still literal `"~/NFLWork/..."` strings, so their deparsed source does not contain `nflwork_root()`.

- [ ] **Step 3: Change each reader default in Tools.R**

Replace each hardcoded default (exact line content shown — search-and-replace each):

```r
# get_wagerzon_odds (~line 3066)
    db_path = "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "wagerzon_odds", "wagerzon.duckdb")

# get_hoop88_odds (~3220)
    db_path = "~/NFLWork/hoop88_odds/hoop88.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "hoop88_odds", "hoop88.duckdb")

# get_bfa_odds (~3362)
    db_path = "~/NFLWork/bfa_odds/bfa.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "bfa_odds", "bfa.duckdb")

# get_bookmaker_odds (~3498)
    db_path = "~/NFLWork/bookmaker_odds/bookmaker.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "bookmaker_odds", "bookmaker.duckdb")

# get_dk_odds (~3736)
    db_path = "~/NFLWork/dk_odds/dk.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "dk_odds", "dk.duckdb")

# get_fd_odds (~3862)
    db_path = "~/NFLWork/fd_odds/fd.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "fd_odds", "fd.duckdb")

# get_bet105_odds (~3986)
    db_path = "~/NFLWork/bet105_odds/bet105.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "bet105_odds", "bet105.duckdb")

# get_kalshi_odds (~4122)
    db_path = "~/NFLWork/kalshi_odds/kalshi.duckdb"
# becomes:
    db_path = file.path(nflwork_root(), "kalshi_odds", "kalshi.duckdb")
```

And the team-dict path (~line 4261):

```r
  dict_db <- sprintf("~/NFLWork/Answer Keys/%s.duckdb", sport)
# becomes:
  dict_db <- file.path(nflwork_root(), "Answer Keys", sprintf("%s.duckdb", sport))
```

Note: `get_wagerzon_betting_odds` (~line 4511) keeps `db_path = "~/NFLWork/wagerzon_odds/wagerzon.duckdb"` as a forwarded default — change it the same way to `file.path(nflwork_root(), "wagerzon_odds", "wagerzon.duckdb")` so the forwarded path also routes.

- [ ] **Step 4: Run test to verify it passes**

Run: `cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_nflwork_root.R")'`
Expected: PASS (all tests, including the routing test).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_nflwork_root.R"
git commit -m "feat(mlb): route Tools.R odds readers through nflwork_root()"
```

---

## Task 3: MLB.R — script-relative setwd + set repo root

Replace the hardcoded `setwd("~/NFLWork/Answer Keys")` with derivation from the script location, then opt this process into worktree-aware paths via `set_nflwork_root()`. On `main` the derived dir is identical (`~/NFLWork/Answer Keys`).

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (lines 5-22)

- [ ] **Step 1: Replace the setwd + add root opt-in**

In `Answer Keys/MLB Answer Key/MLB.R`, replace:

```r
.t_script_start <- Sys.time()
setwd("~/NFLWork/Answer Keys")
```

with:

```r
.t_script_start <- Sys.time()
# Resolve cwd from the script location so a worktree copy stays self-contained.
# MLB.R lives at <root>/Answer Keys/MLB Answer Key/MLB.R, so its Answer Keys dir
# is the parent of the script dir. On main this is identical to ~/NFLWork/Answer Keys.
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
.script_dir <- if (length(.file_arg) > 0) {
  .raw <- gsub("~\\+~", " ", sub("^--file=", "", .file_arg[1]))
  normalizePath(dirname(.raw), mustWork = FALSE)
} else {
  normalizePath("~/NFLWork/Answer Keys/MLB Answer Key", mustWork = FALSE)
}
setwd(dirname(.script_dir))  # <root>/Answer Keys
```

- [ ] **Step 2: Opt this process into worktree-aware reader paths**

Immediately after the existing `source("Tools.R")` line, add:

```r
set_nflwork_root(derive_repo_root())  # so get_*_odds read this repo's DBs
```

- [ ] **Step 3: Verify MLB.R parses and resolves correctly on main**

Run (from the worktree — this checks the derivation logic without a full pipeline run):
```bash
Rscript -e 'src <- readLines("Answer Keys/MLB Answer Key/MLB.R"); cat("setwd line:\n"); cat(grep("setwd", src, value=TRUE), sep="\n")'
```
Expected: shows `setwd(dirname(.script_dir))` and no remaining `setwd("~/NFLWork/Answer Keys")`.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "feat(mlb): MLB.R derives cwd + repo root from script location"
```

---

## Task 4: mlb_dashboard.R — remove the strip-hack

`DASHBOARD_DIR` is already derived from the script location, so `dirname(dirname(DASHBOARD_DIR))` already yields the correct root (the worktree in a worktree, `~/NFLWork` on main). Only the strip-hack forces it back to main.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the `NFLWORK_ROOT` block near line 6587)

- [ ] **Step 1: Delete the strip-hack**

Replace this block:

```r
# Resolve NFLWork root — if running from a worktree, use main repo for mlb.duckdb
NFLWORK_ROOT <- dirname(dirname(DASHBOARD_DIR))
if (grepl(".claude/worktrees", NFLWORK_ROOT, fixed = TRUE)) {
  NFLWORK_ROOT <- sub("/.claude/worktrees.*", "", NFLWORK_ROOT)
}
setwd(NFLWORK_ROOT)
```

with:

```r
# Resolve NFLWork root from the script location so a worktree render reads the
# worktree's own DBs. On main this is identical to ~/NFLWork.
NFLWORK_ROOT <- dirname(dirname(DASHBOARD_DIR))
setwd(NFLWORK_ROOT)
```

- [ ] **Step 2: Verify the hack is gone**

Run: `grep -n "claude/worktrees" "Answer Keys/MLB Dashboard/mlb_dashboard.R"`
Expected: no output (no matches).

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb): dashboard render reads worktree DBs (remove strip-hack)"
```

---

## Task 5: mlb_dashboard_server.py — remove strip-hack + configurable port

Remove the `_REPO_ROOT` strip-hack so closing-odds capture uses the natural script-relative `PROJECT_ROOT` (the worktree). Make the port overridable via env var so a test instance can run on 8093 alongside the live one on 8083 (default unchanged → `main` identical).

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (lines 138-140, 2730)

- [ ] **Step 1: Remove the strip-hack**

Replace:

```python
# DuckDB with pipeline bets — always in main repo (not worktree)
_REPO_ROOT = PROJECT_ROOT
if ".claude/worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".claude/worktrees")[0])
```

with:

```python
# Repo root for closing-odds capture — script-relative so a worktree test reads
# the worktree's own scraper DBs. On main this is identical to ~/NFLWork.
_REPO_ROOT = PROJECT_ROOT
```

- [ ] **Step 2: Make the port configurable**

Add `import os` if not already imported (check the top of the file first). Replace:

```python
    print(f"\nOpen in browser: http://localhost:8083")

    app.run(host="0.0.0.0", port=8083, debug=False)
```

with:

```python
    _port = int(os.environ.get("MLB_DASHBOARD_PORT", "8083"))
    print(f"\nOpen in browser: http://localhost:{_port}")

    app.run(host="0.0.0.0", port=_port, debug=False)
```

- [ ] **Step 3: Verify**

Run:
```bash
grep -n "claude/worktrees" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"; \
grep -n "MLB_DASHBOARD_PORT" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```
Expected: no `claude/worktrees` match; one `MLB_DASHBOARD_PORT` match.

- [ ] **Step 4: Confirm Python still imports**

Run: `cd "Answer Keys/MLB Dashboard" && python3 -c "import ast; ast.parse(open('mlb_dashboard_server.py').read()); print('OK')"`
Expected: `OK` (syntax valid).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "feat(mlb): server uses worktree root + configurable port"
```

---

## Task 6: seed_test_data.sh — CoW-clone live DBs into the worktree (TDD)

Seed the worktree with the DBs scrapers don't produce, plus current scraper DBs so `MLB.R`/render can run without re-scraping. Uses `cp -c` (APFS copy-on-write). Idempotent; `--refresh` re-clones.

**Files:**
- Create: `Answer Keys/MLB Dashboard/seed_test_data.sh`
- Test: inline bash assertions (Step 4)

- [ ] **Step 1: Write the script**

Create `Answer Keys/MLB Dashboard/seed_test_data.sh`:

```bash
#!/bin/bash
# Seed the current worktree with CoW clones of main's live DuckDB files so the
# worktree can run the MLB pipeline + dashboard in full isolation. Safe to run
# while the live dashboard / RFQ bot are running: cp is a byte-level read of
# main's files and never blocks or mutates them.
set -euo pipefail

MAIN_ROOT="$HOME/NFLWork"
WT_ROOT="$(git rev-parse --show-toplevel)"

if [ "$WT_ROOT" = "$MAIN_ROOT" ]; then
  echo "Refusing to seed: you are in the main repo, not a worktree." >&2
  exit 1
fi

# DBs to clone: <relative path under repo root>
DBS=(
  "Answer Keys/pbp.duckdb"
  "Answer Keys/mlb.duckdb"
  "Answer Keys/mlb_mm.duckdb"
  "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
  "dk_odds/dk.duckdb"
  "fd_odds/fd.duckdb"
  "wagerzon_odds/wagerzon.duckdb"
  "hoop88_odds/hoop88.duckdb"
  "bfa_odds/bfa.duckdb"
  "bookmaker_odds/bookmaker.duckdb"
  "bet105_odds/bet105.duckdb"
)

REFRESH=0
[ "${1:-}" = "--refresh" ] && REFRESH=1

for rel in "${DBS[@]}"; do
  src="$MAIN_ROOT/$rel"
  dst="$WT_ROOT/$rel"
  if [ ! -f "$src" ]; then
    echo "  skip (no source): $rel"
    continue
  fi
  if [ -f "$dst" ] && [ "$REFRESH" -eq 0 ]; then
    echo "  keep (exists):   $rel"
    continue
  fi
  mkdir -p "$(dirname "$dst")"
  # CoW clone on APFS; falls back to a normal copy on non-APFS volumes.
  cp -c "$src" "$dst" 2>/dev/null || cp "$src" "$dst"
  echo "  seeded:          $rel"
done

echo "Seed complete in: $WT_ROOT"
```

- [ ] **Step 2: Make it executable**

Run: `chmod +x "Answer Keys/MLB Dashboard/seed_test_data.sh"`

- [ ] **Step 3: Run it**

Run: `bash "Answer Keys/MLB Dashboard/seed_test_data.sh"`
Expected: a list of `seeded:` / `keep:` / `skip:` lines and `Seed complete in: <worktree path>`.

- [ ] **Step 4: Assert clones landed in the worktree, not main**

Run:
```bash
WT="$(git rev-parse --show-toplevel)"; \
test -f "$WT/Answer Keys/pbp.duckdb" && echo "pbp present in worktree: OK"; \
test -f "$WT/dk_odds/dk.duckdb" && echo "dk present in worktree: OK"; \
git status --porcelain | grep -E '\.duckdb' && echo "ERROR: duckdb tracked by git" || echo "duckdb gitignored: OK"
```
Expected: both `OK` presence lines and `duckdb gitignored: OK`.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/seed_test_data.sh"
git commit -m "feat(mlb): seed_test_data.sh CoW-clones live DBs into a worktree"
```

---

## Task 7: test_dashboard.sh — seed → run slice → serve on test port

A thin wrapper so testing is one command. Defaults to render-only (fastest); flags run more.

**Files:**
- Create: `Answer Keys/MLB Dashboard/test_dashboard.sh`

- [ ] **Step 1: Write the wrapper**

Create `Answer Keys/MLB Dashboard/test_dashboard.sh`:

```bash
#!/bin/bash
# Run the MLB dashboard from the current worktree against seeded, isolated data.
# Usage:
#   ./test_dashboard.sh                 # seed if needed, render-only, serve :8093
#   ./test_dashboard.sh --mlb           # re-run MLB.R against seeded DBs, then serve
#   ./test_dashboard.sh --pipeline      # full run.py mlb (live re-scrape), then serve
#   ./test_dashboard.sh --refresh ...   # re-clone seed data first
set -euo pipefail

DASH_DIR="$(cd "$(dirname "$0")" && pwd)"
WT_ROOT="$(git -C "$DASH_DIR" rev-parse --show-toplevel)"
PORT="${MLB_DASHBOARD_PORT:-8093}"
export MLB_DASHBOARD_PORT="$PORT"

MODE="render"
REFRESH=""
for arg in "$@"; do
  case "$arg" in
    --refresh)  REFRESH="--refresh" ;;
    --mlb)      MODE="mlb" ;;
    --pipeline) MODE="pipeline" ;;
    --render)   MODE="render" ;;
    *) echo "Unknown arg: $arg" >&2; exit 1 ;;
  esac
done

echo "== Seeding worktree data =="
bash "$DASH_DIR/seed_test_data.sh" $REFRESH

cd "$WT_ROOT"
case "$MODE" in
  pipeline)
    echo "== Full pipeline (run.py mlb) =="
    python3 "Answer Keys/run.py" mlb || echo "pipeline had errors, continuing"
    ;;
  mlb)
    echo "== MLB.R only =="
    Rscript "Answer Keys/MLB Answer Key/MLB.R" || echo "MLB.R had errors, continuing"
    ;;
esac

echo "== Rendering dashboard =="
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" || echo "render had errors, continuing"

echo "== Starting test server on :$PORT =="
(sleep 2 && open "http://localhost:$PORT") &
python3 "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```

- [ ] **Step 2: Make it executable + syntax-check**

Run:
```bash
chmod +x "Answer Keys/MLB Dashboard/test_dashboard.sh"; \
bash -n "Answer Keys/MLB Dashboard/test_dashboard.sh" && echo "syntax OK"
```
Expected: `syntax OK`.

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/test_dashboard.sh"
git commit -m "feat(mlb): test_dashboard.sh wrapper (seed + run slice on :8093)"
```

---

## Task 8: Integration verification (guard + no-regression + CBB smoke)

Prove the invariant both directions and confirm shared `Tools.R` didn't break CBB.

**Files:** none created — verification only.

- [ ] **Step 1: Worktree render produces a report from seeded data**

Run (from the worktree):
```bash
bash "Answer Keys/MLB Dashboard/seed_test_data.sh"; \
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"; \
WT="$(git rev-parse --show-toplevel)"; \
test -f "$WT/Answer Keys/MLB Dashboard/report.html" && echo "report rendered in worktree: OK"
```
Expected: `report rendered in worktree: OK`.

- [ ] **Step 2: Prove the readers route to the worktree, not main**

Run:
```bash
cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_nflwork_root.R")'
```
Expected: PASS — the routing test confirms every `get_*_odds` resolves under the set root.

- [ ] **Step 3: No-regression — main's paths unchanged**

Run (asserts default/fallback resolves to `~/NFLWork`, i.e. unset root behaves as today):
```bash
cd "Answer Keys/tests" && Rscript -e '
source("../Tools.R")
set_nflwork_root(NULL)
stopifnot(identical(nflwork_root(), normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)))
cat("fallback resolves to ~/NFLWork: OK\n")'
```
Expected: `fallback resolves to ~/NFLWork: OK`.

- [ ] **Step 4: CBB smoke (shared Tools.R)**

Run the existing Tools-dependent CBB test(s) to confirm the reader change didn't break other sports:
```bash
cd "Answer Keys/tests" && Rscript -e 'testthat::test_file("test_probit_devig.R")'
```
Expected: PASS. (If CBB has a dedicated odds-reader test, run that too.)

- [ ] **Step 5: Confirm no `.duckdb` is staged for commit anywhere**

Run: `git status --porcelain | grep -E '\.duckdb' && echo "ERROR" || echo "clean: OK"`
Expected: `clean: OK`.

---

## Task 9: Documentation

Update docs in the same branch so they merge with the code.

**Files:**
- Modify: `Answer Keys/CLAUDE.md`, `CLAUDE.md` (root), `.gitignore`
- Create: `Answer Keys/MLB Dashboard/README.md` (if absent)
- Modify: memory `verify_ui_features_by_rendering.md`

- [ ] **Step 1: .gitignore — ensure WAL is ignored**

Confirm `*.duckdb` is present (it is, line 44). Add a global WAL ignore after it if not already global:

```
*.duckdb.wal
```

Run: `grep -n "duckdb.wal" .gitignore || echo "add it"` and add the line under the existing `*.duckdb` entry if missing.

- [ ] **Step 2: Answer Keys/CLAUDE.md — DuckDB section**

In the `### DuckDB Databases` section, replace the line:

```
- Never symlink DuckDB files. Always copy if needed in worktrees.
```

with:

```
- Never symlink DuckDB files. To test from a worktree, run
  `Answer Keys/MLB Dashboard/seed_test_data.sh` to CoW-clone the live DBs into
  the worktree, then `test_dashboard.sh` to render on port 8093. Paths now
  derive from the running code's location (`Tools.R::nflwork_root()` /
  `derive_repo_root()`), so worktree code reads worktree DBs and `main` is
  unaffected.
```

- [ ] **Step 3: Root CLAUDE.md — worktree testing note**

In the **Branching workflow** / worktree section, add a bullet:

```
- **Testing the MLB dashboard from a worktree:** run
  `"Answer Keys/MLB Dashboard/seed_test_data.sh"` then
  `"Answer Keys/MLB Dashboard/test_dashboard.sh"` (serves on :8093). The worktree
  is self-contained — it reads/writes its own seeded DBs and never touches main's.
```

- [ ] **Step 4: Create Answer Keys/MLB Dashboard/README.md (if absent)**

If the file does not exist, create it:

```markdown
# MLB Dashboard

Port 8083 (live). See `Answer Keys/CLAUDE.md` for full architecture.

## Testing changes from a worktree

The worktree is self-contained: all DB paths derive from where the code lives,
so a worktree reads/writes its own seeded copies and never touches main's live
data.

1. `./seed_test_data.sh` — CoW-clone live DBs into this worktree (instant on
   APFS). `--refresh` re-clones.
2. `./test_dashboard.sh` — render against seeded data and serve on **:8093**.
   - `--mlb` re-runs `MLB.R` (re-price) first.
   - `--pipeline` runs the full `run.py mlb` (live re-scrape) first.
3. Seeded `*.duckdb` files are gitignored and deleted when the worktree is
   removed.

`MLB_DASHBOARD_PORT` overrides the server port (default 8083 live / 8093 test).
```

- [ ] **Step 5: Update memory**

Edit `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/verify_ui_features_by_rendering.md` to note: render from a worktree via `seed_test_data.sh` + `test_dashboard.sh` (self-contained worktree), not by copying DBs manually.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/CLAUDE.md" "CLAUDE.md" ".gitignore" "Answer Keys/MLB Dashboard/README.md"
git commit -m "docs(mlb): document self-contained worktree testing"
```

---

## Pre-Merge Checklist (executive-engineer review)

- [ ] `git diff main..HEAD` reviewed; the only behavior change on `main` is *additive* (new helpers, configurable port default unchanged, strip-hacks removed are no-ops on main).
- [ ] `grep -rn "claude/worktrees" "Answer Keys"` returns nothing (both strip-hacks gone).
- [ ] Guard test (Task 8 Steps 2-3) passes both directions.
- [ ] CBB smoke (Task 8 Step 4) passes.
- [ ] No `.duckdb` staged anywhere.
- [ ] Docs updated (Task 9).
- [ ] Get explicit user approval, then merge and clean up the worktree + branch.
```
