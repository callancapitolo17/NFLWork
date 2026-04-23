# Cross-Book Grid — Actual-Price Display with Devigged Median Anchor

**Date:** 2026-04-22
**Status:** Design approved, ready for implementation plan.
**Owner:** callan
**Related files:** `nfl_draft/lib/queries.py`, `nfl_draft/lib/devig.py`, `nfl_draft/lib/quarantine.py`, `kalshi_draft/app.py`.

## Problem

The Cross-Book Grid's docstring claims each cell is a venue's devigged
fair-value estimate. Investigation (2026-04-22) of the ingest path and live
DuckDB contents shows that **no sportsbook rows are actually devigged**:

- `nfl_draft/lib/quarantine.py:46-55` writes `devig_prob = implied_prob`
  whenever the scraper doesn't pre-compute a devig (no sportsbook scraper does).
- Verified live: for every book in `('draftkings', 'fanduel', 'bookmaker',
  'wagerzon', 'hoop88', 'betonline')`, `implied_prob == devig_prob` on 100% of
  rows. Only Kalshi splits the two values.

Consequences:

1. Each sportsbook cell on the grid is the raw vig-inclusive implied
   probability, not a fair estimate.
2. The "median" column is dragged upward by sportsbook vig — typically
   5–10pp above the true fair.
3. Kalshi (no vig) sits correctly near fair but reads artificially "low"
   vs the inflated median, producing false-flag bias against Kalshi.
4. Genuine +EV edges on sportsbooks (where the sportsbook price sits far
   from true fair) are invisible: every sportsbook sits at a similar
   vig-inflated level, so the grid's median-deviation flag misses them.

## Goal

Turn the Cross-Book Grid into an Unabated-style actual-price screen:

- Cells show the **bettable price** (American odds for sportsbooks, cents for
  Kalshi).
- The median anchor is the **true devigged fair** across venues.
- The outlier flag compares each venue's **raw take price** to the **median
  of fairs**, so a flagged cell is a direct +EV signal.

## Non-goals

- No change to the Kalshi scraper or to Kalshi's existing bid/ask
  tooltip behavior.
- No change to the Trade Tape, Bet Log, or +EV Candidates layout. +EV
  Candidates benefits automatically from the new `delta` semantics.
- No sharp-book anchoring (Kalshi-anchored fair / weighted-sharp consensus).
  Median-of-all-devigged-venues stays the anchor.
- No alternative devig methods (power, worst-case). Proportional devig only,
  with the known incomplete-coverage bias documented.

## Model

Three layers:

### Math layer (hidden)

Each `(book, outright_group)` bucket is devigged at ingest time via
`proportional_devig` (= normalize so implieds sum to 1.0). `outright_group`
is derived **authoritatively from the scraper's `market_type`**, not by
string-stripping `market_id` — `nfl_draft/lib/market_map.py:89` already has
`build_market_id(market_type, **kwargs)`; we add a symmetric
`build_outright_group(market_type, **kwargs)` so both come from the same
shared metadata:

| market_type                 | outright_group                            | example |
|-----------------------------|--------------------------------------------|---------|
| `pick_outright`             | `pick_{N}_overall`                         | `pick_2_overall`   |
| `first_at_position`         | `first_{position}`                         | `first_wr`         |
| `top_n_range`               | `top_{range_high}`                         | `top_10`           |
| `team_first_pick`           | `team_{team}_first_pick`                   | `team_washington_first_pick` |
| `team_first_pick_position`  | `{team}_first_pick_pos`                    | `arizona_first_pick_pos`     |
| `nth_at_position`           | `{nth}_{position}`                         | `2_wr`             |
| `draft_position_over_under` | `draft_position_ou_{player}_{line}`        | (2-way: over/under) |
| `matchup_before`            | `matchup_{a}_before_{b}`                   | (2-way: sorted)     |
| `prop`, `mr_irrelevant_position` | *no group* (pass through, no devig)  | —                  |

A book's `devig_prob` for a market is `implied_prob / sum(implieds in that
book's outright_group)`.

**Why not a string-strip heuristic on `market_id`?** The heuristic succeeds
for `pick_N_overall` and `first_X` but silently produces 1-row groups for
`team_*_first_pick_pos_*` and `nth_at_position` markets because the
candidate suffix is variable-length (e.g. `wide`, `defensive_line`,
`kicker_punter_long_snapper`). Using `market_type` makes the grouping
deterministic and coverage-complete.

### Flag layer

For every venue, the outlier comparison is:

```
|implied_prob_venue − median(devig_prob across all posting venues)| ≥ threshold_pp
```

Signed delta (`implied_prob − median_fair`) preserves direction:

- `delta > 0` → price is above fair → YES is overpriced → bet NO at this book
- `delta < 0` → price is below fair → YES is underpriced → bet YES at this book

Kalshi's `implied_prob` is already `yes_ask / 100` (the take price), so the
same formula applies uniformly. The existing Kalshi-specific flag branch in
`queries.py` is dropped.

### Display layer

| Cell         | Value shown                        | Format                  |
|--------------|------------------------------------|-------------------------|
| Sportsbook   | `american_odds` from `draft_odds`  | `+110`, `−150`          |
| Kalshi       | `implied_prob * 100`               | `42¢`                   |
| Median col   | `median(devig_prob)` across venues | `42.5%`                 |
| Outliers col | count of flagged venues            | integer                 |
| ⚑ suffix     | unchanged                          | unchanged               |

Red-background conditional styling in `TABLE_STYLE_DATA_CONDITIONAL` already
matches "cell contains ⚑", so it keeps working with no CSS change.

## Components

### 1. Outright group constructor — `nfl_draft/lib/market_map.py`

Add `build_outright_group(market_type: str, **kwargs) -> Optional[str]`
symmetric to `build_market_id`. Returns `None` for `market_type` in
`{"prop", "mr_irrelevant_position"}` (single-sided or non-outright). All
other market types return the group key per the table above.

Unit-tested in `tests/unit/test_market_map.py` alongside `build_market_id`
with a case per market_type.

### 2. `OddsRow` extension — `nfl_draft/scrapers/_base.py`

Add optional field:

```python
outright_group: Optional[str] = None
```

No schema change to `draft_odds` — the group is only needed during ingest,
not persisted.

### 3. Scrapers emit the group — each file in `nfl_draft/scrapers/`

Every call site that constructs an `OddsRow` also computes
`build_outright_group(market_type, **kwargs)` with the same kwargs it
already passes to `build_market_id`. This is a mechanical mirror of ~3–5
call sites per scraper (DK, FD, Bookmaker, Wagerzon, Hoop88, BetOnline,
Kalshi).

For the Kalshi scraper, `outright_group` is set but unused (Kalshi rows
pre-populate `devig_prob` from mid; the ingest devig pass skips them).
Still emit it for consistency — future changes may need it.

### 4. Devig at ingest — `nfl_draft/lib/quarantine.py`

Replace the per-row write loop in `write_or_quarantine` with:

1. **Bucket** every incoming `OddsRow` by `(row.book, row.outright_group)`.
   Rows with `outright_group is None` skip bucketing entirely.
2. For each bucket:
   - If any row in the bucket has `devig_prob` already set (Kalshi), leave
     every row's `devig_prob` as the scraper emitted. In practice a bucket
     is from one scraper so this is all-or-nothing; assert the invariant
     defensively and log a warning if mixed.
   - Else if the bucket has ≥2 rows, compute
     `proportional_devig([row.implied_prob for row in bucket])` and assign
     in order.
   - Else (1-row bucket, or rows with `outright_group is None`): set
     `devig_prob = implied_prob`. Matches current behavior — grid still
     renders, flag still works using raw.
3. Write rows to `draft_odds` with the populated `devig_prob`.

`implied_prob` semantics are unchanged: `american_to_implied(odds)` for
sportsbooks, `yes_ask / 100` for Kalshi.

### 5. Flag semantics — `nfl_draft/lib/queries.py`

`cross_book_grid` (lines 89–179) changes:

- SQL adds `american_odds` to the latest-snapshot CTE and the returned row.
- Per-market Python loop:
  - `display_by_book` is no longer computed (the callback will render cells
    directly from `american_odds` / `implied_prob`).
  - `median` is still `statistics.median(devig_prob across books with valid
    devig_prob)`.
  - Flag is `|implied_prob_book − median| ≥ threshold` for **every** book —
    no Kalshi special-case.
- Return shape extends to include `american_odds` and `implied_prob`
  per book alongside `devig_prob`, so the callback can format cells.

`ev_candidates` needs a one-line change: `delta = book["implied_prob"] −
median` instead of `book["devig_prob"] − median`. This makes the delta
directly the EV in percentage points.

Docstrings are rewritten to describe the new semantics.

### 6. Cell rendering — `kalshi_draft/app.py`

In `_update_crossbook` (line 1069 onward), replace the cell-formatting block
(lines 1109–1120):

```python
for venue in VENUES:
    record = m["books"].get(venue)
    flagged = m["flags"].get(venue, False)
    if record is None:
        row[venue] = ""
        continue
    if venue == "kalshi":
        cell = f"{round(record['implied_prob']*100)}¢"
    else:
        cell = f"{int(record['american_odds']):+d}"
    row[venue] = cell + (" ⚑" if flagged else "")
```

Header copy on `render_crossbook_grid` is rewritten to reflect the new
semantics: "Cells show each book's posted price. The median is the
cross-venue devigged fair. ⚑ fires when a book's price sits more than
`threshold` pp from that fair — i.e. a direct +EV signal."

## Edge cases

| Case | Behavior |
|------|----------|
| Outright group with 1 posted candidate at a book | `devig_prob = implied_prob`; cell still renders; flag uses raw. |
| Single-sided Yes/No prop (`top_10_*` with no NO companion) | 1-row bucket per book. Median across venues still meaningful. No false math. |
| Outright group with total implied < 1.0 | Proportional devig would inflate instead of deflate. Surface a warning log line; don't special-case. Indicates a scraping gap. |
| Market posted by only 1 venue | Unchanged from today (`queries.py:145-153`): no median, no flags. |
| Kalshi one-sided market (bid only, no ask) | `implied_prob` is NULL. Cell renders empty, no flag — same as today. |
| Sub-3pp threshold + sportsbook vig | Documented caveat in grid header: sportsbook cells sit ~1–3pp above fair by construction even with zero edge. |

## File-level change summary

| File                                         | Change                                                | Size          |
|----------------------------------------------|-------------------------------------------------------|---------------|
| `nfl_draft/lib/market_map.py`                | Add `build_outright_group`                            | ~15 lines     |
| `nfl_draft/scrapers/_base.py`                | Add `outright_group` field to `OddsRow`               | 1 line        |
| `nfl_draft/scrapers/{dk,fd,bm,wz,h88,bol,kalshi}.py` | Populate `outright_group` at each OddsRow emit | ~3 lines × 7 scrapers |
| `nfl_draft/lib/quarantine.py`                | Bucket + proportional_devig before write              | ~25 lines     |
| `nfl_draft/lib/queries.py`                   | Add `american_odds` to CTE + result; drop Kalshi flag branch; update `ev_candidates.delta` | ~15 lines |
| `kalshi_draft/app.py`                        | Cell formatter per venue; header copy                 | ~20 lines     |
| `nfl_draft/tests/unit/test_market_map.py`    | Add `build_outright_group` tests                      | new tests     |
| `nfl_draft/tests/unit/test_queries.py`       | Update flag assertions; new devig-at-ingest tests     | new tests     |
| `nfl_draft/README.md`                        | Rewrite "Cross-Book Grid" paragraph                   | 1 paragraph   |
| `nfl_draft/lib/backfill_devig.py`            | One-shot historical re-devig script (delete after run) | ~30 lines    |

## Testing

Unit:

- `build_outright_group` — table-driven tests mirroring
  `test_market_map.py::test_build_market_id` with a case per
  `market_type`, including the `prop`/`mr_irrelevant_position` → `None`
  cases.
- `write_or_quarantine` devig path — fixture batch with mixed groups
  (n-way, 1-row, Kalshi pre-devigged, `outright_group=None` prop) → assert
  per-row `devig_prob` matches hand-computed values.
- `cross_book_grid` — fixture DB with known fair / raw values per book →
  assert flag fires iff `|implied_prob − median_fair| ≥ threshold` for
  every book including Kalshi.
- `ev_candidates.delta` — assert it's `implied_prob − median`, not
  `devig_prob − median`.

Integration: run the real `nfl_draft/run.py --mode scrape --book draftkings`
against a temp DuckDB fixture, inspect that `devig_prob` is populated
correctly for an outright group.

Manual: start the dashboard, compare a row I know the current behavior of
(e.g. `first_wr_jordyn-tyson` per the brainstorm example) and confirm the
new grid shows raw prices + a lower median + a Bookmaker flag at ~+11pp
delta.

## Backfill

Historical `draft_odds` rows keep their wrong `devig_prob` unless rewritten.
A one-shot `nfl_draft/lib/backfill_devig.py` reads all existing rows, groups
by `(fetched_at batch, book, outright_group)`, re-runs `proportional_devig`,
and updates `devig_prob` in place. Run once from main after merge, then
delete the script (temp file per repo convention).

## Version control + worktree

- **Branch:** `feature/grid-devig-and-price-display`
- **Worktree path:** `/Users/callancapitolo/NFLWork-worktrees/grid-devig-display`
- **DuckDB for testing in worktree:** copy `nfl_draft/nfl_draft.duckdb` (never
  symlink — WAL loss on cleanup).
- **Commits (roughly):**
  1. `build_outright_group` + `OddsRow.outright_group` field + unit tests
  2. Each scraper populates `outright_group` on emitted rows
  3. Devig at ingest in `write_or_quarantine` + integration test
  4. `cross_book_grid` + `ev_candidates` flag semantics change + tests
  5. `app.py` cell renderer + header copy
  6. `nfl_draft/README.md` docs update
  7. Backfill script + run it + delete it
- **Merge flow:** test on worktree → re-run full `pytest nfl_draft/tests/` → pre-merge
  executive review → get explicit user approval → merge to `main` → remove worktree +
  delete branch.

## Documentation

Updated in the same merge to `main`:

- `nfl_draft/README.md` — Cross-Book Grid section (~1 paragraph rewrite).
- No `CLAUDE.md` update needed; this change doesn't alter project
  conventions, only fixes an implementation bug + UI formatting.

## Risks

- **Proportional devig known bias** on incomplete posted sets (`devig.py:29-36`
  already documents this). Numerically, books posting only top-5 of a 30-way
  outright will under-devig by ~10–15%. Not a blocker at the 10pp default
  threshold; revisit if we tighten the slider below 5pp.
- **Backfill correctness.** The script re-devigs *per scrape batch*, not per
  all-time. If we've ever written inconsistent `fetched_at` times within a
  single scrape, groups could split. Mitigation: bucket by
  `date_trunc('minute', fetched_at)` rather than exact timestamp.
- **Market-map drift.** Unmapped rows that land in `draft_odds_unmapped`
  are unchanged by this work. The grid never saw them.
