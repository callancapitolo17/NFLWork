# NFL Draft Dashboard — 20-Minute Staleness Filter

**Date:** 2026-04-23
**Status:** Approved (pre-implementation)
**Scope:** Tighten the cross-venue dashboard's row-age filter from 2 hours to 20 minutes.

## Problem

The dashboard at `kalshi_draft/app.py` surfaces cross-venue draft odds read from
`nfl_draft/nfl_draft.duckdb`. Today, any row newer than 2 hours is eligible for
display (`nfl_draft/lib/queries.py:42`, `MAX_AGE_HOURS = 2`).

Two hours is too wide. In pre-draft mode the scrape cron runs every 15 minutes
(`nfl_draft/crontab.pre`); on draft day it runs every 2 minutes
(`nfl_draft/crontab.draft`). If a single venue silently stops posting, its
stale price can linger in the grid for up to 2 hours and pollute both the
Cross-Book Grid outlier flag and the +EV Candidates tab. The FanDuel regression
on 2026-04-19 — where FD's draft tab disappeared upstream — is the motivating
incident.

## Decision

Replace the hour-scoped constant with a minute-scoped one:

```
MAX_AGE_MINUTES = 20
```

Twenty minutes is 5 minutes above the pre-draft scrape cadence (15 min) — a
small cushion for late cron runs — and ~10× the draft-day cadence, so a truly
dead venue is hidden within ~20 minutes instead of ~2 hours.

### Alternatives considered

| Option | Pros | Cons | Why rejected |
|---|---|---|---|
| Mode-aware threshold (15 min draft / 30 min pre-draft) | Matches each cadence exactly | Two constants to maintain; requires wiring to the mode toggle | YAGNI — flat 20 min covers both cases with one knob |
| Strict 15 min flat | Tightest dead-venue detection | Equals pre-draft cron cadence exactly; any cron delay causes healthy rows to flicker in/out | Flicker risk outweighs the 5-min gain |
| Per-venue thresholds | Each venue tuned to its own cadence | Adds config surface; no evidence venues need different values today | YAGNI |

## Scope of Change

### Source code

**`nfl_draft/lib/queries.py`**

- Rename the module-level constant `MAX_AGE_HOURS = 2` (line 42) to
  `MAX_AGE_MINUTES = 20`.
- Rewrite the comment block above it (lines 36–42) to document the new
  reasoning (15-min pre-draft scrape cadence + 5-min cushion, ~10×
  draft-day cadence).
- Update both SQL call sites to use minute-scoped intervals:
  - `cross_book_grid` (line 121):
    `INTERVAL '{MAX_AGE_HOURS} hours'` → `INTERVAL '{MAX_AGE_MINUTES} minutes'`
  - `kalshi_tooltip_data` (line 371): same substitution.

### Tests

**`nfl_draft/tests/integration/test_dashboard_queries.py`**

- Line 766: change the import from `MAX_AGE_HOURS` to `MAX_AGE_MINUTES`.
- Line 767: change the stale-timestamp computation from
  `timedelta(hours=MAX_AGE_HOURS + 1)` to
  `timedelta(minutes=MAX_AGE_MINUTES + 5)` — still comfortably past the
  cutoff, just in units that match the new constant.
- Lines 389 and 758: update docstrings that name `MAX_AGE_HOURS` so they
  reference the renamed constant.

### Documentation

**`nfl_draft/README.md`** — the existing line about "the dashboard's
staleness filter hides FD rows while they're gone" should state the current
threshold (20 min) so future readers don't need to open `queries.py` to find
it.

### Explicitly out of scope

- `docs/superpowers/specs/2026-04-21-kalshi-cross-grid-pricing-design.md`
  and its plan reference the 2-hour value as a fact about the code at that
  time. These are historical artifacts; the source of truth going forward is
  the code and this spec. We do not rewrite prior design docs.
- The Dash auto-refresh interval (60s default, 15s in draft-day mode) is
  unchanged — this spec only concerns which rows the queries are willing to
  return, not how often the UI polls DuckDB.
- Scrape cron cadences (`crontab.pre`, `crontab.draft`) are unchanged.

## Behavior After the Change

- Rows older than 20 minutes are excluded from Cross-Book Grid, +EV
  Candidates, and Kalshi tooltip data.
- In healthy pre-draft state, a row lives ~15 min between refreshes; the
  20-min threshold provides a 5-min cushion for a late cron run.
- On draft day, a dead venue is filtered out within ~20 min instead of ~2 h.

## Known Tradeoff

If a pre-draft cron run hangs for more than ~5 minutes (e.g. slow Playwright
startup, a temporary network stall), every venue's row can fall past the
threshold simultaneously and the grid will render empty until the next run
completes. The prior 2-hour filter would have concealed that by continuing to
show the last scrape. We accept this: an empty grid is a clearer signal than
a grid full of stale prices, and the +EV Candidates tab is more dangerous
than the Cross-Book Grid when it's showing phantom edges.

## Acceptance Criteria

1. Constant `MAX_AGE_MINUTES = 20` exists in `nfl_draft/lib/queries.py`;
   `MAX_AGE_HOURS` no longer appears anywhere in the `nfl_draft/` tree.
2. Both SQL queries (`cross_book_grid`, `kalshi_tooltip_data`) use
   `INTERVAL '20 minutes'` (via the constant) as the row-age cutoff.
3. `pytest nfl_draft/tests/integration/test_dashboard_queries.py` passes,
   including the two regression tests that assert stale rows are excluded.
4. `nfl_draft/README.md` states the current staleness threshold.
5. Dashboard boots cleanly against a live DuckDB and the Cross-Book Grid +
   +EV Candidates + Kalshi tooltip tabs render without errors.

## Version Control

- Feature branch: `feature/staleness-20min`.
- Worktree created before any code change; removed after merge.
- Commits:
  1. Source + test changes in a single commit (they are interdependent).
  2. README update in the same commit or immediately after (must land on
     `main` in the same merge).
- Pre-merge review: full diff against `main`, then explicit user approval
  before merging.
