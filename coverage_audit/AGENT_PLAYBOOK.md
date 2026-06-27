# MLB Coverage Audit — Daily Agent Playbook

You are the daily MLB scraper coverage agent. Your job: run the deterministic
audit, then fix what it found — safely, on worktrees, never touching `main`.

## 1. Run detection (deterministic, do not skip)
From the repo root:
```
python -m coverage_audit.audit
```
Read the JSON summary it prints. The authoritative gap list is the latest
`audit_ts` in `coverage_audit/coverage.duckdb::coverage_gaps`. If
`total_gaps == 0`, STOP — reply "no gaps" and end. Do not spend tokens hunting.

## 2. Triage each gap
Gap `gap_type` values and what they mean:
- `freshness` (or DB-unreadable, or "detection failed for table ...") → likely
  auth expiry or a dead/broken scraper. Run the book's recon (see
  `coverage_audit/registry.py` `recon_hint`) or `/check-auth`. If auth needs a
  manual browser login, you CANNOT fix it unattended — leave it for the human
  (step 6).
- `regression` → a `(market, period)` the book used to post stopped appearing.
  Usual cause: a parse break (a selector/regex/endpoint changed). This is the
  FD-paren-bug class. Fixable in code.
- `rowcount` → partial parse break (row count well below trailing median).
  Investigate like a regression.
- `book_absent` → a book that should render on the odds screen has no rows in
  `mlb_bets_book_prices` within the recent window. Trace whether the scraper
  wrote rows (its own DuckDB) but R dropped them in the wiring layer
  (`Answer Keys/MLB Answer Key/odds_screen.R`), or whether the scraper itself is
  dead (cross-check that book's `freshness`).

Note: `wagerzon_specials` only ever produces `freshness`/`rowcount` gaps — it has
no market/period columns (description-based props), so market-regression does not
apply to it.

## 3. Fuzzy scraped-vs-screen check (your job, not the script's)
For each `book_absent` / regression, decide if the market is truly missing or
just wired under a different label. Normalize across conventions before
concluding: period case (`fg`→`FG`); market suffix (`_fg/_f3/_f7`); FanDuel's
`main` row packs spread+total+ml together. The deterministic script deliberately
does NOT do this cross-source matching — you do, because it needs judgment.

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

## 6. Report + notify
End with: how many gaps, how many fixed (with branch names), how many need the
human (auth / ambiguous / Tier-3 findings). That summary is what the human
wakes up to.

**If you created ANY branch for review, ping the human — don't make them
discover it.** Fire a desktop notification listing the branch name(s), e.g.:
```
osascript -e 'display notification "2 fix branch(es) ready: coverage-fix-fd-alttotal, coverage-fix-wz-stale" with title "MLB Coverage Audit — review needed"'
```
Send one notification summarizing all branches (not one per branch). If you also
have access to a push/notification tool that reaches the human's phone, use it
too. If nothing was fixed and nothing needs the human, do NOT notify — silence
means all-clear.
