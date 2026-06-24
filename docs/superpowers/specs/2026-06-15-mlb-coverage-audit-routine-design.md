# MLB Scraper Coverage Audit + Auto-Wiring Routine — Design

## Review Pack

**What we're building** — A local Claude **Desktop scheduled task** that runs
once a day. A deterministic Python script first diffs, per book, what each book
*posts* vs. what your scraper *captures* vs. what reaches the odds screen, and
writes a gap report. The agent then reads that report and immediately works each
fixable gap on its own git worktree (writing the parser/regex/R fix, running
tests, committing to a branch), and notifies you to review + merge. It never
touches `main` unattended. North star: a complete Unabated-style odds screen.

**Key decisions**

1. **Desktop scheduled task, not a cloud routine.** Cloud routines (`/schedule`)
   execute on Anthropic infrastructure with only a fresh `git clone` — no access
   to your local DuckDBs, auth tokens, or residential IP, so they physically
   cannot run your scrapers. Desktop scheduled tasks run on your Mac with local
   file + network access and don't need an open session. (Rejected: cloud
   routine — can't reach the data; `CronCreate`/`/loop` — need a live session.)
2. **Deterministic detection, intelligent wiring.** The diff is done by a plain
   Python script, not the LLM. The agent's judgment is reserved for writing fixes.
   (Rejected: letting the agent compute the gap report — unreliable and wasteful
   for a `SELECT ... EXCEPT`.)
3. **Tiered detection so "all 9 books" is feasible day one.** Tier 1 (trailing
   self-baseline) + Tier 2 (raw-vs-parsed) ship deterministically for all books
   and catch every *regression*. Tier 3 (discover never-seen markets) is done by
   the agent live when it runs, the same way a scraper gets built — not a
   pre-built 9-book menu-crawler. (Rejected: building a full menu enumerator for
   all 9 books upfront — too much work for v1.)
4. **Worktree + test + notify; you merge.** Every fix lands on an isolated
   worktree branch with tests run; the merge gate is yours. Honors your standing
   "never auto-merge to main" rule. (Rejected: auto-merge on green tests — a
   subtly-wrong parser could feed bad odds into the screen before you see it.)

**Risks / push back here**

- **Once-a-day timing misses late-posting markets.** A ~9 AM Pacific run won't see
  alt ladders a book only posts an hour before first pitch (e.g. Bookmaker's
  multi-index grid "near start time"). Acceptable for a structural *coverage*
  audit, but it means the audit is not a completeness guarantee for late markets.
- **A green test is not a correct parser.** Worktree isolation + mandatory test
  pass + never-merge are the guardrails, but a fix that passes tests yet
  mis-parses can still sit on a branch looking ready. Your diff review is the real
  backstop — do not rubber-stamp.
- **Token cost.** A full agent runs daily. It must no-op cheaply when the script
  reports zero gaps; it should only spend real tokens when there is wiring to do.
- **Tier-3 scope creep.** "Discover never-seen markets" is open-ended. The agent
  must be bounded (probe a defined set of endpoints/tabs per book, time-boxed)
  or it will wander.

**Worth understanding**

- **Idempotent / deterministic vs. agentic work** — The detection script is
  *deterministic*: same inputs → same output, every run, no surprises (like a pure
  R function with no random seed). The wiring agent is *non-deterministic*: it
  reasons and can produce different code each run. The design deliberately puts
  the trustworthy, repeatable work in the script and the creative work in the
  agent. This is why the script — not the LLM — is the source of truth for "what's
  broken."
- **Set difference as a diagnostic** — The whole audit is three set differences
  (`A − B`, `B − C`). In R you'd reach for `setdiff()`; in SQL / DuckDB it's
  `EXCEPT`. Framing coverage as set math is what makes it cheap and reliable.

---

## Design Body

### 1. Goal

Catch — and immediately begin fixing — any MLB market that one of *your* scrapers
(non–Odds-API books) could be surfacing but isn't, whether because the scraper
never captures it or because it's captured and dropped downstream before the odds
screen. Today there is no recurring live scrape and failures are silent
(`run.py` swallows scraper exceptions via `asyncio.gather(return_exceptions=True)`
and exits cleanly on "0 games" auth expiry), so gaps go unnoticed for weeks.

### 2. Architecture

```
Desktop scheduled task fires daily (~9 AM Pacific)
        │
        ▼
[1] run  coverage_audit.py   ← plain Python, deterministic, no LLM
        │   • computes the gap report (three-stage diff, per book)
        │   • writes it to coverage/coverage.duckdb::coverage_gaps
        │   • fires a macOS notification ONLY if the gap set changed
        │   • exits 0 with a machine-readable summary the agent reads
        ▼
[2] agent reads coverage_gaps
        │   • for each FIXABLE gap → create a worktree, write the
        │     parser/regex/R-union fix, run tests, commit to a branch
        │   • for AMBIGUOUS gaps → log + leave for the human
        │   • if zero gaps → no-op, cheap exit
        ▼
[3] notify: "3 gaps found, 2 fixed on branches X/Y, 1 needs you"
```

The deterministic script is the source of truth for *what is broken*. The agent
is the arm that *writes code*. Each does only what it is good at.

### 3. Detection — the three-stage chain

For every book, a market can break at three points:

```
A: book POSTS it   →   B: scraper CAPTURES it   →   C: screen SHOWS it
        A − B = "scraper gap" (book has it, scraper never wrote it)
        B − C = "wiring gap"  (scraper wrote it, downstream dropped it)
```

- **Stage B** = distinct market types written to the book's own DuckDB
  (`<book>_odds/<book>.duckdb::mlb_odds`) in the latest run.
- **Stage C** = distinct `market × book` rows actually present in
  `mlb_mm.duckdb::mlb_bets_book_prices` (what the dashboard renders).
- **Stage A** is the expensive one and ships in tiers (below).

`B − C` is the FD-paren-bug class (scraped, regex dropped it). `A − B` is the
FD-event-page-tab class (book posts on a tab the scraper never fetched).

**Planning-time refinement (verified against live DBs 2026-06-15):** the
*deterministic script* does NOT compute the cross-source `B − C` / `A − B` diff,
because real data showed it is false-positive-prone — only 5 of 9 books render on
the screen at all (hoop88/bfa/kalshi feed pricing, not pills), and books pack
markets under different labels (FD's `main` row carries spread+total+ML; period
case and `_fg/_f3/_f7` suffixes differ per book). Instead the script does the
*reliable, network-free* work — per-book regression (Tier 1, same-book
self-baseline, no normalization needed), freshness, row-count, and
**screen-presence** (an `expected_on_screen` book missing from the screen). The
fuzzy cross-source label matching (`B − C`) and live Tier-2/3 discovery move to
the agent at runtime, where judgment lives. This still catches the historical
bugs: the FD paren break shows up as "FD stopped emitting `alternate_totals` it
posted all last week" via Tier 1.

#### Detection tiers (all 9 books covered from v1)

| Tier | "What the book posts" (Stage A) | Catches | Cost |
|------|----------------------------------|---------|------|
| **1 — v1, all 9** | what *this book itself* produced over the trailing N days (self-baseline) | regressions: silent scraper death, broke parses, dropped market types | ~free |
| **2 — v1 where easy** | raw markets present in the scraper's *own API response* before its filter | parse-drops on books that already fetch broadly (FD, DK) | cheap |
| **3 — agent-driven** | live probe of the book's full menu for *never-seen* markets | true expansion gaps (Unabated completeness) | agent does this when it runs, time-boxed |

Tiers 1–2 are deterministic and catch every regression in the problem statement.
Tier 3 (discovering markets never captured) is what a Claude agent does when it
builds a scraper, so the agent performs it live rather than us pre-building a
menu-crawler for 9 books.

#### Ride-along checks (free while each DB is open)

- **Freshness:** `MAX(fetch_time)` per book; flag if older than ~2× expected
  cadence (catches the live-probe "132 min stale" case).
- **Row-count sanity:** latest run's row count vs trailing median per book
  (DK ~580, bet105 ~300, hoop88 ~44, bfa ~55, bookmaker ~130, kalshi ~132,
  wagerzon ~266, fanduel ~330) — catches parse breaks that still write *some* rows.
- **Auth:** if a book's detection fails on auth, reuse the existing `recon_*`
  scripts to classify it as an `auth` gap — reported, never a crash. BFA checked
  first (expires fastest).

### 4. Wiring — worktree, test, notify, never merge

- Each fixable gap → its own **git worktree** (isolated from the live checkout and
  any parallel jobs).
- The agent writes the fix, then runs `tests/timezone_parity_test.py` plus the
  relevant scraper smoke test. A fix that fails tests is left on the branch
  flagged "tests failing," not merged.
- It commits to a branch and **stops at the merge gate.** The human reviews the
  diff and merges. No unattended writes to `main`.
- **DuckDB-in-worktree hazard (CLAUDE.md rule #5):** never symlink a `.duckdb`
  into the worktree (WAL loss on removal). Testing a scraper fix writes to the
  book's real DB; the implementation plan must define a safe rule here (e.g. run
  the scraper from the main checkout path, or copy the DB).

### 5. Scheduling, storage, alerts

- **Mechanism:** Desktop scheduled task — local, runs on your Mac with local
  file + network access, no active Claude *session* required. Created via the
  Claude Code **Desktop app → Routines → New routine → Local** (or by describing
  it in a Desktop session), **not** the `/schedule` CLI command (that makes a
  cloud routine, which can't reach your data). **Caveat:** the Desktop *app* must
  be running for the task to fire — unlike launchd, which needs nothing open. If
  the detection half ever needs to be bulletproof with the app closed, run
  `coverage_audit.py` alone on launchd and keep only the wiring as a Desktop task.
- **Timing:** ~9 AM Pacific (noon ET) — late enough that the day's slate is
  posted, before first pitch. One run/day. (See risk: late-posting markets.)
- **Storage:** a dedicated `coverage/coverage.duckdb` with its own write lock
  (mirrors the one-DB-per-concern pattern; never shares a lock with the pipeline).
  Tables: `coverage_gaps` (current + history), and the trailing baselines the
  detection needs.
- **Alerts:** macOS `osascript` notification, the same channel `run_fetch.sh`
  uses, fired only when the gap set *changes* — no daily noise when all is wired.

### 6. Out of scope (YAGNI)

- Cloud routines (cannot reach your machine).
- Fast/live refresh loops (this is a once-a-day coverage tool, not a freshness loop).
- Auto-merge to `main`.
- The Unabated odds screen itself (north star, a separate build the audit feeds).

---

## Version Control

- **Branch / worktree:** `worktree-mlb-coverage-audit-routine` (already created).
  Implementation proceeds here; merge to `main` only after review + tests.
- **New files (anticipated):**
  - `coverage/coverage_audit.py` — deterministic detection script.
  - `coverage/README.md` — setup, the detection tiers, how to read the report.
  - A Desktop-scheduled-task definition (registered via Claude Code, documented
    in the README; mirrors the launchd pattern in spirit but is the local-agent
    mechanism).
- **Modified files (anticipated):** none in v1 detection beyond additive new
  files; wiring fixes are produced on their own per-gap worktrees, not in this
  branch.

## Documentation

- New `coverage/README.md` (setup, tiers, report schema, how the scheduled task
  is registered, how to add a book/market to the audit).
- Update root `CLAUDE.md` Project Structure with a bullet for the coverage audit.
- Update `Answer Keys/CLAUDE.md` if the audit reads any pipeline tables not yet
  documented (it reads `mlb_bets_book_prices`, already documented).
- A memory entry summarizing the routine + the cloud-vs-desktop-routine finding.

## Open implementation questions (for the plan, not blocking this spec)

1. Trailing-baseline window `N` (days) for Tier 1 — start at 7?
2. Exact "fixable vs ambiguous" decision rule the agent uses before opening a
   worktree (to bound Tier-3 wandering and token spend).
3. The safe DuckDB-test rule for scraper fixes in worktrees.
4. Notification payload shape (how much detail in the macOS toast vs. the report).
