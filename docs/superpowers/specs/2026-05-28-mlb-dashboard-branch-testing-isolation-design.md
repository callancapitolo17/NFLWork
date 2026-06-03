# MLB Dashboard — Isolated Branch Testing Design

**Date:** 2026-05-28
**Branch:** `worktree-mlb-dashboard-worktree-testing`
**Status:** Design — awaiting user review

---

## Review Pack

**What we're building** — A way to test the MLB dashboard (and the pipeline
behind it) from inside a git worktree against a *private snapshot* of the live
data, so a branch test never reads from or writes to the main repo's live
DuckDB files. One new env var (`NFLWORK_DATA_ROOT`) tells every script where the
data lives; a new `test_dashboard.sh` harness seeds a worktree-local copy and
runs whatever slice of the pipeline you want, then renders the dashboard.

**Key decisions**

1. **One env var as the single source of truth for the data location**
   (`NFLWORK_DATA_ROOT`, default `~/NFLWork`). Rejected alternative: per-DB
   read/write path splitting. Why: a single DuckDB file (e.g. `dk.duckdb`) is
   *both* written by the scraper and read by `MLB.R` in the same run — they must
   resolve to the same place, so one root for all pipeline DBs is the only
   correct model.

2. **Copy-on-write (`cp -c`) clones to seed the snapshot, and we only clone DBs
   the test won't regenerate.** Rejected alternative: a full byte copy, and/or
   pausing the live processes before copying. Why: on APFS `cp -c` is instant
   and ~free; and the one DB with a torn-write risk (`mlb.duckdb`, held under the
   pipeline's long write lock) is *never cloned* because `MLB.R` rebuilds it from
   scratch. So no live process ever needs to pause.

3. **Production behaves byte-identically when the env var is unset.** Rejected
   alternative: making all paths script-relative (so the worktree is implicitly
   self-contained). Why: script-relative resolution would change how the *live*
   system resolves paths too — more regression risk on exactly the thing we
   don't want to disturb. An unset env var → the literal `~/NFLWork` we use
   today.

4. **A mandatory audit of every `.duckdb` / `dbConnect` site in the MLB code
   paths.** Why: the failure mode we are eliminating is "one path still points at
   main." A single missed site silently reads/writes the live data and defeats
   the whole feature, so enumerating every site *is* the core of the work, not a
   side task.

**Risks / push back here**

- **Scope of the refactor touches the live pipeline.** `MLB.R`, `Tools.R`, and
  `mlb_dashboard.R` are production files the live dashboard and (indirectly) your
  workflow depend on. **Decision: full refactor of all path sites now.** The
  render-path-only half-measure was rejected because it does not achieve
  isolation for data-layer changes — which is half of the stated requirement
  ("both render and data layer must work"). A partial routing would leave the
  exact alt-spread-style mismatch in place. Veto at the review gate if you want
  it staged anyway.
- **The snapshot can be stale relative to "right now."** A CoW clone captures the
  bytes at clone time. If the live pipeline writes `mlb_mm.duckdb` a second after
  you clone, your snapshot misses it. For testing this is fine (re-run the
  harness), but flag if you expected the test to always track the absolute latest
  live state — that's a different (and harder) goal.
- **`mlb_dashboard.duckdb` is 147MB of state (placed bets, settings, CLV).**
  **Decision: clone the real state DB** so the test dashboard renders against
  real settings, bankroll, placed bets, and CLV — which is what makes a render
  test trustworthy. Any write a test makes (e.g. placing a bet) lands in the
  `$SNAP` copy, never main. CoW makes the clone free. Veto at the review gate if
  you'd rather branch tests start from an empty state DB.

**Worth understanding**

- **Environment variables as configuration (vs. R's `options()`).** An env var is
  a key/value the OS hands to a process and *all of its child processes*
  automatically. That auto-inheritance is why it's the right tool here: `run.sh`
  spawns `run.py`, which spawns `MLB.R` and the scrapers — set the var once at the
  top and every child sees it, no plumbing. The R analogue you know is
  `Sys.getenv("NAME", "default")`, which is exactly like `options()$name %||%
  default` but readable across the process boundary.
- **Copy-on-write (CoW) clones.** `cp -c` on APFS doesn't duplicate the bytes — it
  creates a second directory entry pointing at the *same* physical blocks, and
  only the blocks that later get *written* are duplicated ("diverge"). So cloning
  a 2.4GB file is instant and uses no extra disk until something changes it.
  There's no R analogue; the closest mental model is R's copy-on-modify for
  vectors (`y <- x` doesn't copy until you write to `y`), but at the filesystem
  level.

---

## Problem

Testing the MLB dashboard from a git worktree against live data is unreliable
because the code resolves DuckDB paths **four different ways**, and most of them
pin to the main repo regardless of where the code runs:

| Site | How it resolves today | Effect from a worktree |
|---|---|---|
| `MLB.R:6` `setwd("~/NFLWork/Answer Keys")` | absolute home path | reads/writes **main's** `pbp/mlb/mlb_mm.duckdb` |
| `Tools.R` `get_dk_odds()` / `get_fd_odds()` defaults | absolute `~/NFLWork/...` | reads **main's** `dk.duckdb` / `fd.duckdb` |
| `mlb_dashboard.R:6585` `.claude/worktrees` strip hack | rewrites worktree path → main | dashboard reads **main's** DBs |
| `mlb_sgp/scraper_*_singles.py` | **script-relative** (`Path(__file__).parent.parent`) | scraper writes the **worktree's** `dk.duckdb` |

The inconsistency is the bug. Concretely, this is exactly why the alt-spread
scraper fix was painful to verify: the worktree's *scraper* wrote corrected rows
to `<worktree>/dk_odds/dk.duckdb`, while the worktree's `MLB.R` (via `Tools.R`)
read `~/NFLWork/dk_odds/dk.duckdb` — main's untouched copy. The fix was invisible
on the dashboard.

### What must work

Both kinds of change must be testable from a worktree, against live-shaped data,
without touching main:

- **Render-layer changes** (CSS / layout / JS / table rendering) — need live data
  flowing in to see the visual result.
- **Data-layer changes** (scraper / pricing / `odds_screen.R` / `MLB.R`) — need
  the *data the worktree code produces* reflected on the dashboard.

And the verification must run any slice: one scraper, just `MLB.R`, or the full
`run.py mlb` — because "how much do I need to re-run" varies per change.

---

## Design

### Mechanism: one env var + two tiny resolvers

A single environment variable names where the *data* lives:

```
NFLWORK_DATA_ROOT   (default: ~/NFLWork)
```

Two resolvers read it:

- **R** (`Tools.R`, sourced by everything):
  ```r
  data_path <- function(rel) {
    root <- Sys.getenv("NFLWORK_DATA_ROOT", unset = path.expand("~/NFLWork"))
    file.path(root, rel)
  }
  ```
- **Python** (a shared helper, e.g. `mlb_sgp/paths.py`):
  ```python
  def data_root() -> Path:
      return Path(os.environ.get("NFLWORK_DATA_ROOT", str(_DEFAULT_ROOT)))
  ```

When the var is unset (production, the live dashboard, the RFQ bot) everything
resolves to `~/NFLWork` exactly as today. **No production behavior change.**

### The key subtlety: code paths ≠ data paths

`MLB.R` currently does `setwd("~/NFLWork/Answer Keys")` and uses that one cwd for
**both** `source("Tools.R")` (code) and `dbConnect("pbp.duckdb")` (data). These
must be split:

- **Code** resolves relative to the running script, so the worktree runs the
  worktree's code. (`mlb_dashboard.R` already derives `DASHBOARD_DIR`; `MLB.R`
  will derive its own script dir the same way and `setwd()` there for `source()`.)
- **Data** routes through `data_path()`, so it points at the snapshot.

This split is the heart of the refactor and where the regression risk lives.

### Touch points (the complete routing list)

1. **`MLB.R`** — replace the hardcoded `setwd` with a script-relative `setwd` for
   sourcing; wrap every `dbConnect(dbdir = ...)` for `pbp.duckdb`, `mlb.duckdb`,
   `mlb_mm.duckdb`, `mlb_dashboard.duckdb` in `data_path()`.
2. **`Tools.R` — ALL nine odds readers, not just DK/FD.** `get_wagerzon_odds`,
   `get_hoop88_odds`, `get_bfa_odds`, `get_bookmaker_odds`, `get_dk_odds`,
   `get_fd_odds`, `get_bet105_odds`, `get_kalshi_odds`, and the team-dict
   reader (`dict_db`, ~line 4266) each have an absolute `~/NFLWork/<book>/<db>`
   default. **Every one** must default to `data_path("<book>/<db>.duckdb")`, or
   `MLB.R` silently reads main's scraper DBs even in the minimal `--mlb` case.
   (Found in review — the original spec named only two.)
3. **`mlb_dashboard.R`** — delete the `.claude/worktrees` strip hack; route its DB
   paths (`mlb_mm.duckdb`, `DB_PATH` → `mlb_dashboard.duckdb`) through
   `data_path()`.
4. **`mlb_dashboard_server.py` — MANDATORY, the largest write surface.** (Missed
   in the original spec; caught in review.) The harness runs the server so you
   can interact with the test dashboard, and the server has **30+
   `duckdb.connect()` calls, many read-write** (`/api/place-bet`,
   `/api/place-parlay`, `/refresh`, `/api/update-bet`, `/api/remove-bet`). It
   also carries its **own copy of the strip-hack** (`_REPO_ROOT`, lines 138-140)
   plus `BASE_DIR` (line 28) and `REPO_ROOT` (line 123) feeding `DB_PATH`,
   `MLB_MM_DB`, `DASHBOARD_DB`, and the per-book DBs. Delete the strip-hack and
   route every connection through a Python `data_root()`. Without this, clicking
   Place/Refresh in a branch test writes to **main's** `mlb_dashboard.duckdb` —
   defeating the entire feature.
5. **Python scrapers** (`scraper_draftkings_singles.py`,
   `scraper_fanduel_singles.py`) — output path = `data_root() / "dk_odds" /
   "dk.duckdb"` etc. (falling back to the current script-relative default when
   the var is unset). v1 routes only the DK/FD writers (see Scope below); the
   other books' DBs are cloned into the snapshot, not re-scraped.
6. **`mlb_correlated_parlay.R` / `mlb_triple_play.R`** — same `data_path()`
   routing for any `mlb.duckdb` / `mlb_mm.duckdb` connections (they run in
   `run.sh` between the pipeline and the render).
7. **`run.py`** — no direct DB opens (only code-relative script/sentinel paths);
   the env var propagates to its scraper/R children automatically. No change
   needed, confirmed in review.
8. **Audit gate (mandatory):** grep every `.duckdb` / `dbConnect` / `dbdir` /
   `duckdb.connect` reference in the MLB code paths and confirm each is routed.
   A missed site silently escapes to main. The audit list doubles as the
   pre-merge review checklist. The review above already expanded the list twice
   (Tools.R readers, the server) — treat the audit as authoritative over this
   prose.

### The harness: `test_dashboard.sh`

Lives in `Answer Keys/MLB Dashboard/`. Steps:

1. `SNAP="$(git rev-parse --show-toplevel)/.test_data"` — a gitignored dir inside
   the worktree. **Wipe and recreate it at the start of every run** (unless
   `--keep` is passed) so snapshots never accumulate or go stale. CoW makes
   re-seeding instant, so starting fresh is free.
2. **Seed it** with CoW clones (`cp -c`) of the DBs the run won't regenerate:
   - always: `pbp.duckdb` (read-only historical), `mlb_dashboard.duckdb` (state).
   - the scraper DBs (`dk_odds/dk.duckdb`, `fd_odds/fd.duckdb`, etc.) **only if**
     this run isn't re-scraping them.
   - `mlb.duckdb` / `mlb_mm.duckdb` are **not** cloned when the run regenerates
     them (`MLB.R` rebuilds from scratch).
3. `export NFLWORK_DATA_ROOT="$SNAP"`.
4. Run the requested slice via a flag:
   - `--render-only` — just `Rscript mlb_dashboard.R` (render-layer changes).
   - `--scraper dk|fd` — re-run a routed scraper into the snapshot, then render
     (the alt-spread case).
   - `--mlb` — run `MLB.R` (+ parlay/trifecta) against the cloned scraper DBs,
     then render (pricing / `odds_screen.R` changes). This is the realistic
     "full refresh" for most data-layer work.
5. Serve on a **non-8083 port** (e.g. 8093) so it never collides with the live
   dashboard, and `open` it. The harness must pass this port to
   `mlb_dashboard_server.py` (which currently assumes 8083).

**Cleanup.** Copied data never lingers: `$SNAP` is wiped and re-seeded at the
start of each run (step 1), `--keep` opts out when you want to inspect a snapshot
afterward, and the whole directory is deleted when the worktree is removed
(gitignored, lives inside the worktree). WAL files land in `$SNAP` too (DuckDB
writes the WAL next to the DB path, and `$SNAP` is a real directory, not a
forbidden symlink), so they are cleaned up by the same mechanisms.

Because seeding is a filesystem CoW copy (independent of DuckDB's advisory
locks), and because the only long-write-locked DB (`mlb.duckdb`) is never cloned,
**no live process — dashboard or RFQ bot — ever needs to pause.** The RFQ bot is
categorically unaffected: it writes only its own `kalshi_mlb_rfq*.duckdb` and
reads main's `mlb.duckdb`, which a branch test never writes.

### Scope: why no full `--pipeline` flag in v1 (review finding)

A true full `run.py mlb` re-scrapes **all seven books** live. For that to stay
isolated, every scraper writer (`wagerzon`, `hoop88`, `bfa`, `bookmaker`,
`bet105`, `kalshi`, in addition to `dk`/`fd`) would also need `data_root()`
routing — roughly tripling the Python surface and the regression risk, against a
case you rarely actually test (you're testing *your* change, not a 7-book live
refresh).

**v1 covers the three realistic test shapes** — render-layer (`--render-only`),
a single changed scraper (`--scraper dk|fd`, the alt-spread case), and re-pricing
against the latest cloned book data (`--mlb`). All nine *readers* are routed
(cheap, mechanical) so `--mlb` sees the snapshot; only the DK/FD *writers* are
routed. The other books' DBs are CoW-cloned into the snapshot and not re-run.

**Deferred to v2:** route the remaining scraper writers and add `--pipeline` for a
fully-live multi-book refresh. Flagged here so the limitation is explicit, not a
silent gap. **This is the main scope call to confirm.**

### Safety / isolation argument (why main is never at risk)

- **Reads of main are byte-level only**, at seed time, via `cp`. They cannot
  block or corrupt anything.
- **All writes go to `$SNAP`** inside the worktree (env var redirects them).
- **A torn clone only harms the throwaway copy.** If `cp -c` catches
  `mlb_mm.duckdb` or `mlb_dashboard.duckdb` mid-write, the *snapshot* may be
  slightly inconsistent; re-run the harness. Main is untouched.

### Testing the refactor itself

The regression we most fear is "a path site still points at main." Add a guard
test:

- Set `NFLWORK_DATA_ROOT` to a temp dir, run a dry render (and a minimal pipeline
  invocation), and assert that **every** DuckDB file opened resolves *under* that
  temp dir — failing loudly if any path escapes to `~/NFLWork`.
- Mechanism options: a thin `dbConnect` wrapper in test scope that records opened
  paths, or `lsof`/`fs_usage` capture around the run. Decide during planning.

---

## Version Control

- **Branch / worktree:** `worktree-mlb-dashboard-worktree-testing` (already
  created). All code changes happen here; merge to `main` only after review +
  user approval. Clean up the worktree + branch after merge.
- **Files created:** `Answer Keys/MLB Dashboard/test_dashboard.sh`,
  `mlb_sgp/paths.py` (Python resolver), the guard test, this spec, and an
  implementation plan.
- **Files modified:** `MLB.R`, `Tools.R` (all nine odds readers),
  `mlb_dashboard.R`, `mlb_dashboard_server.py` (delete its strip-hack + route all
  connections), `mlb_correlated_parlay.R`, `mlb_triple_play.R`, the DK/FD Python
  scrapers, `.gitignore` (add `.test_data/`). The Python `data_root()` resolver
  must live somewhere both the scrapers (`mlb_sgp/`) and the server
  (`Answer Keys/MLB Dashboard/`) can import — decide the shared location during
  planning.
- **Commit structure:** (1) R + Python resolvers + `.gitignore`; (2) route all
  path sites through the resolvers (the audit); (3) `test_dashboard.sh` harness;
  (4) guard test; (5) docs. Each commit independently runnable on main with the
  env var unset.
- **Pre-merge review:** executive-engineer review of `git diff main..HEAD` with
  the audit list as the checklist — confirm no path site escapes to main when the
  env var is set, and that production (env unset) is byte-identical.

## Documentation

Updated in the same merge to `main`:

- **`Answer Keys/CLAUDE.md`** — DuckDB section: document `NFLWORK_DATA_ROOT`, the
  code-vs-data path split, and replace the "Never symlink … always copy in
  worktrees" note with a pointer to the harness.
- **`CLAUDE.md` (root)** — branch-hygiene / worktree section: how to test the MLB
  dashboard from a worktree (`test_dashboard.sh`).
- **`Answer Keys/MLB Dashboard/README.md`** (create if absent) — harness usage,
  flags, and the non-8083 test port.
- Memory: update `verify_ui_features_by_rendering.md` to reference the harness.
