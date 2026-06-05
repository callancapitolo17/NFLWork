# MLB Dashboard — Isolated Branch Testing (Self-Contained Worktree)

**Date:** 2026-05-28
**Branch:** `worktree-mlb-dashboard-worktree-testing`
**Status:** Design — awaiting user review

---

## Review Pack

**What we're building** — A way to fully test MLB dashboard changes inside a git
worktree before merging to `main`, where the worktree reads and writes its *own*
copy of the data and never touches main's live DuckDB files. We do this by making
the worktree **self-contained**: the codebase is already ~90% script-relative
(every scraper and the server already resolve paths relative to where the code
lives), so we only remove four deliberate "jump back to main" overrides, then
seed the worktree's data once with instant copy-on-write clones.

**Key decisions**

1. **Make the worktree self-contained by deriving every data path from the
   running code's location — no env var, no separate environment.** Rejected
   alternative: a `NFLWORK_DATA_ROOT` env var pointing at a `.test_data/`
   snapshot. Why rejected: that fights the existing convention — it would force
   us to change *every scraper* to honor the env var, whereas the scrapers
   *already* write script-relative. Deriving-from-location changes **zero
   scrapers** and just removes four special-case overrides. Less code, lower
   risk, consistent with what the live system already does.

2. **The change is byte-identical on `main`.** Each fix derives the repo root
   from the code's own path, which on main *is* `~/NFLWork`; removing the
   strip-hacks is a no-op on main (the path never contains `.claude/worktrees`
   there). So production behavior does not change at all. This is the
   non-negotiable invariant, proven by a guard test (below).

3. **Seed the worktree with copy-on-write (`cp -c`) clones of the live DBs the
   scrapers don't produce.** Rejected alternative: a full byte copy, or pausing
   live processes first. Why: on APFS `cp -c` is instant and ~free even for the
   2.4GB `pbp.duckdb`, it's a byte-level read that can't block or corrupt
   anything, and no live process (dashboard or RFQ bot) ever needs to pause.

4. **A guard test that fails if any DB path escapes the worktree.** Why: the
   failure mode we're eliminating is "a path still points at main." The guard
   makes that impossible to ship silently — it's the proof behind decision #2.

**Risks / push back here**

- **Lost convenience: a worktree can no longer peek at main's *live* data for a
  zero-setup render.** The old strip-hack let a worktree dashboard read main's
  current DBs directly. After this, you seed the worktree first (one fast
  command). Given you want isolation, that's the right trade — flag if you relied
  on the live peek.
- **`Tools.R` is shared across CBB / NFL / MLB.** The reader-path change touches
  all sports. The fallback (below) keeps any caller that doesn't opt in behaving
  exactly as today, but CBB must be smoke-tested before merge per the repo rule.
- **Seeded state includes your real placed bets (`mlb_dashboard.duckdb`, 147MB).**
  We clone it so the test dashboard looks real; writes during a test land in the
  worktree copy, never main. Flag if you'd rather tests start from empty state.

**Worth understanding**

- **Script-relative path resolution (`__file__` / `parents[N]`).** A file can ask
  "where am I on disk?" and compute other paths from there. Python:
  `Path(__file__).resolve().parents[2]` walks up two directories from the current
  file. R has no `__file__`, but you can parse it out of `commandArgs()` (the
  flags `Rscript` was launched with). The payoff: code that lives in a worktree
  automatically points at the worktree, with no configuration. There's no direct
  R analogue you've used; the closest idea is `here::here()` if you've seen it.
- **Copy-on-write (CoW) clones.** `cp -c` on APFS doesn't duplicate bytes — it
  points a second directory entry at the *same* physical blocks, and only blocks
  later *written* diverge. So cloning a 2.4GB file is instant and uses no extra
  disk until something changes it. Closest mental model: R's copy-on-modify for
  vectors (`y <- x` doesn't copy until you write `y`), but at the filesystem
  level.

---

## Problem

Testing the MLB dashboard from a worktree is unreliable because of an
**asymmetry**: most of the code already resolves paths relative to where it
lives (so it's worktree-local), but four sites deliberately override that to jump
back to the main repo.

**Already worktree-relative (no change needed):**
- Every scraper: `DB_PATH = Path(__file__).parent / "<book>.duckdb"`
  (`wagerzon`, `hoop88`, `bfa`, `bookmaker`, `bet105`, `dk`, `fd`).
- `mlb_dashboard_server.py`: `REPO_ROOT = Path(__file__).resolve().parents[2]`
  drives `MLB_MM_DB`, `DASHBOARD_DB`, and the per-book paths.

**The four overrides that break isolation:**

| Site | What it does today |
|---|---|
| `MLB.R:6` `setwd("~/NFLWork/Answer Keys")` | absolute → reads/writes **main's** `pbp/mlb/mlb_mm.duckdb` |
| `Tools.R` — 9 `get_*_odds` readers + team-dict | absolute `~/NFLWork/<book>/` → reads **main's** scraper DBs |
| `mlb_dashboard.R:6585` strip-hack | rewrites worktree path → main |
| `mlb_dashboard_server.py:138-140` (`_REPO_ROOT`) strip-hack | rewrites worktree path → main (closing-odds capture) |

The asymmetry *is* the bug. Concretely, this is why the alt-spread scraper fix
was invisible: the worktree's *scraper* correctly wrote to
`<worktree>/dk_odds/dk.duckdb`, but the worktree's `MLB.R` (via `Tools.R`) read
`~/NFLWork/dk_odds/dk.duckdb` — main's untouched copy.

### What must work
Both render-layer changes (CSS / layout / JS) and data-layer changes (scraper /
pricing / `odds_screen.R` / `MLB.R`) must be testable from a worktree against
live-shaped data, running any slice — one scraper, just `MLB.R`, or the full
`run.py mlb` — all without touching main.

---

## Design

### Principle: make the worktree self-contained
Fix the four overrides so they derive the repo root from the running code's
location instead of hardcoding `~/NFLWork`. Everything else already follows that
convention, so the whole worktree — scrapers, pricing, dashboard, server — then
reads and writes within the worktree. Nothing points at main.

### The four fixes

1. **`MLB.R`** — replace `setwd("~/NFLWork/Answer Keys")` with a `setwd` to the
   *script's own* `Answer Keys` dir (derived from `commandArgs()`), and expose a
   `REPO_ROOT` the rest of the script (and `Tools.R`) can use for DB paths.
2. **`Tools.R`** — add a small `nflwork_root()` helper: returns the
   caller-provided root if set, else **falls back to `path.expand("~/NFLWork")`**.
   The nine `get_*_odds` defaults (+ team-dict) become
   `file.path(nflwork_root(), "<book>/<db>.duckdb")`. The fallback is the safety
   net: any caller (CBB, NFL) that doesn't opt in behaves *exactly* as today.
3. **`mlb_dashboard.R`** — delete the `.claude/worktrees` strip-hack; let
   `NFLWORK_ROOT` derive from `DASHBOARD_DIR` (already computed) so it stays in
   the worktree.
4. **`mlb_dashboard_server.py`** — delete the `_REPO_ROOT` strip-hack (lines
   138-140) so closing-odds capture uses the natural script-relative `REPO_ROOT`
   it already computes. The rest of the server is already worktree-relative.

### The core invariant: no effect on `main`
This is the requirement, so here is exactly why it holds:

- **Derivation yields the same paths on main.** `MLB.R` on main lives at
  `~/NFLWork/Answer Keys/MLB Answer Key/MLB.R`, so its derived root is
  `~/NFLWork` — identical to the literal it replaces. The server's `parents[2]`
  on main is already `~/NFLWork`. `Tools.R`'s fallback is `~/NFLWork`.
- **Removing the strip-hacks is a no-op on main.** They only fire when the path
  contains `.claude/worktrees`, which it never does on main.
- **Other sports are untouched.** The `Tools.R` fallback means CBB/NFL callers
  that don't set a root read `~/NFLWork` exactly as before.

So on main, every resolved DB path is byte-identical to today. The change is
*additive*: it gives worktree-resident code a way to point at itself, and does
nothing otherwise.

### Seeding the worktree: `seed_test_data.sh`
A fresh worktree has the code but not the data the scrapers don't produce. The
seed script CoW-clones the live DBs into the worktree's (script-relative)
locations:

- always: `Answer Keys/pbp.duckdb` (2.4GB historical, read-only) and
  `Answer Keys/MLB Dashboard/mlb_dashboard.duckdb` (state).
- the current scraper DBs (`dk_odds/dk.duckdb`, `fd_odds/fd.duckdb`,
  `wagerzon_odds/wagerzon.duckdb`, `hoop88_odds/...`, `bfa_odds/...`,
  `bookmaker_odds/...`, `bet105_odds/...`) so you can run `MLB.R` or render
  *without* re-scraping.
- `mlb.duckdb` / `mlb_mm.duckdb` are produced by `MLB.R`; clone them only for a
  render-only test that skips `MLB.R`.

CoW makes this instant and ~free. The clone is a byte-level read of main's files
— it cannot block or corrupt anything, so **no live process pauses**, and the
RFQ bot (which only writes its own DBs and reads main's `mlb.duckdb`) is
categorically unaffected.

### Running the test
Run the worktree's normal pipeline + dashboard, just on a **non-8083 port** (e.g.
8093) so it never collides with the live dashboard. Options, fastest to fullest:

- **render-only** — `Rscript mlb_dashboard.R` + serve (render-layer changes).
- **re-price** — run `MLB.R` against the cloned scraper DBs, then render
  (pricing / `odds_screen.R` changes). The realistic "full refresh" for most
  data-layer work.
- **one scraper** — re-run any single scraper (it writes worktree-local), then
  `MLB.R`, then render (the alt-spread case — works for *any* book, since all
  scrapers are already worktree-relative).
- **full pipeline** — the worktree's `run.sh` / `run.py mlb` end-to-end; every
  step is worktree-local, so a full live re-scrape stays isolated with no extra
  routing. (Live scrapers need auth/network and can be flaky — that's a property
  of the scrapers, not our isolation.)

A thin wrapper (`test_dashboard.sh`) can bundle "seed if needed → run chosen
slice → serve on 8093 → open", but the underlying mechanism is just the
self-contained worktree.

### Cleanup
Seeded DBs are gitignored (`*.duckdb`) and live inside the worktree, so they
vanish when the worktree is removed (`git worktree remove`) — your normal
post-merge step. `seed_test_data.sh --refresh` re-clones for a fresh snapshot.
DuckDB WAL files land next to each DB inside the worktree (real dirs, not
forbidden symlinks), so they're cleaned up the same way. Copied data never
lingers in main and never leaks back.

### Guard test (the proof behind the invariant)
A test that makes "a path escaped to main" impossible to ship:

- Run a worktree render (and a minimal `MLB.R` invocation) and assert **every**
  DuckDB file opened resolves *under the worktree root* — fail loudly if any path
  resolves under `~/NFLWork` (main).
- Separately, run the same on main and assert paths resolve under `~/NFLWork`
  exactly as before (the no-regression half).
- Mechanism (a thin `dbConnect`/`duckdb.connect` wrapper that records opened
  paths, or `lsof`/`fs_usage` capture) decided during planning.

---

## Version Control

- **Branch / worktree:** `worktree-mlb-dashboard-worktree-testing` (already
  created). All changes here; merge only after review + explicit user approval;
  clean up worktree + branch after merge.
- **Files created:** `Answer Keys/MLB Dashboard/seed_test_data.sh`, optional
  `test_dashboard.sh` wrapper, the guard test, this spec, an implementation plan.
- **Files modified:** `MLB.R` (script-relative `setwd` + `REPO_ROOT`), `Tools.R`
  (`nflwork_root()` helper + 9 reader defaults, with `~/NFLWork` fallback),
  `mlb_dashboard.R` (remove strip-hack), `mlb_dashboard_server.py` (remove
  `_REPO_ROOT` strip-hack), `.gitignore` if needed. **Scrapers: unchanged.**
- **Commit structure:** (1) `Tools.R` helper + reader defaults (fallback keeps
  main identical); (2) remove the three remaining overrides (`MLB.R`,
  `mlb_dashboard.R`, server); (3) `seed_test_data.sh` + optional wrapper;
  (4) guard test; (5) docs. Each commit leaves main byte-identical.
- **Pre-merge review:** executive-engineer review of `git diff main..HEAD`. The
  checklist is the four overrides + the guard test passing in *both* directions
  (worktree paths stay in worktree; main paths unchanged) + a CBB smoke test
  (shared `Tools.R`).

## Documentation

Updated in the same merge to `main`:

- **`Answer Keys/CLAUDE.md`** — DuckDB section: document the self-contained
  worktree model, `seed_test_data.sh`, and that paths now derive from code
  location; replace the "always copy in worktrees" note with a pointer to the
  seed script.
- **`CLAUDE.md` (root)** — worktree/branch-hygiene section: how to test the MLB
  dashboard from a worktree (seed → run on 8093).
- **`Answer Keys/MLB Dashboard/README.md`** (create if absent) — seed + run usage
  and the test port.
- Memory: update `verify_ui_features_by_rendering.md` to reference the seed-and-run
  flow.
