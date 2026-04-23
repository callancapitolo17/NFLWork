# Asymmetric Outlier Flags (Cross-Book Grid)

**Date:** 2026-04-23
**Status:** Approved for implementation
**Scope:** `nfl_draft/lib/queries.py`, `kalshi_draft/app.py` (UI copy),
associated tests and docs.

## Motivation

The Cross-Book Grid flags a venue's price on a given market as an outlier
whenever `|implied_prob - median_devig_prob| >= threshold_pp`. The flag is
direction-agnostic, so it fires in two cases:

1. **Book's YES is cheap vs consensus** (`implied_prob < median`) — a
   bettable +EV YES wager.
2. **Book's YES is expensive vs consensus** (`implied_prob > median`) —
   the book is offering bad odds. On a Kalshi binary you can still
   profit by buying NO. On a YES-only sportsbook futures market (DK,
   FD, Bookmaker, Wagerzon, Hoop88, BetOnline on NFL draft markets),
   there is no way to take the other side. The flag is noise.

The user has confirmed they only want flags for case 1 on non-Kalshi
venues. Kalshi stays symmetric because both YES and NO are purchasable.

## Design

### Flag rule change

In `nfl_draft/lib/queries.py::cross_book_grid`, replace:

```python
flags[book] = (take is not None and abs(take - median) >= threshold)
```

with:

```python
if book == "kalshi":
    flags[book] = (take is not None and abs(take - median) >= threshold)
else:
    flags[book] = (take is not None and (median - take) >= threshold)
```

That is the only logic change. Everything else — median computation,
`MAX_AGE_HOURS` staleness filter, `outlier_count`, signed delta
downstream — stays identical.

### Downstream effects

- `ev_candidates` delegates to `cross_book_grid`, so the filter
  propagates automatically. For non-Kalshi rows in the +EV list,
  `delta = book_prob - median` will always be negative going forward
  (always "bet YES"). Kalshi rows can still be either sign.
- The Cross-Book Grid's `⚑` glyph for non-Kalshi cells will only appear
  on bettable YES edges. `outlier_count` per market drops accordingly.

### UI copy

Update the blurb in `render_crossbook_grid` (`kalshi_draft/app.py`,
around line 834) to explain the asymmetric rule. Current copy claims
every flag is direct +EV, which will still be true — but the copy
should make clear *why* non-Kalshi venues only flag in one direction.

Proposed replacement:

> Each cell shows the bettable price — American odds for sportsbooks,
> cents for Kalshi. The median column is the cross-venue devigged fair
> (the true probability with vig stripped). A flagged cell (⚑) means
> there is a bettable +EV edge: for sportsbooks, the book's YES price
> is cheaper than the fair median by at least the threshold (bet YES);
> for Kalshi, the take price differs from fair in either direction by
> at least the threshold (bet YES if underpriced, NO if overpriced).

### Docstring update

Rewrite the `cross_book_grid` docstring so the rule is explicit:

- Kalshi flag: `abs(implied_prob - median) >= threshold`.
- All other books: `(median - implied_prob) >= threshold`.
- Rationale: sportsbook draft markets are YES-only; flagging an
  overpriced YES is noise the user cannot act on.

## Tests

### Unit (`nfl_draft/tests/unit/test_cross_book_grid.py`)

- **Update** `test_cross_book_grid_flag_uses_implied_vs_median_fair_all_books`:
  Bookmaker @ 0.535 (above median 0.420 by 11.5pp) is currently asserted
  to flag. Under the new rule, a non-Kalshi venue above the median does
  NOT flag. Change the assertion to `False` and rename the test to
  reflect asymmetric-rule semantics, OR add a sibling test and keep this
  one with adjusted setup values so the direction becomes negative
  (book below median).
- **Add** `test_non_kalshi_above_median_not_flagged`: non-Kalshi venue
  with `implied_prob` well above the median and delta >> threshold
  must have `flags[book] == False`.
- **Add** `test_non_kalshi_below_median_flagged`: non-Kalshi venue with
  `implied_prob` well below the median (>= threshold) must have
  `flags[book] == True`.
- **Keep** existing Kalshi tests; they already verify Kalshi is
  symmetric (one test exercises the below-median side, and both
  directions should continue to flag).
- **Add** `test_kalshi_above_median_still_flags`: Kalshi with
  `implied_prob` above median by >= threshold must still flag. This is
  a regression guard against accidentally applying the asymmetric rule
  to Kalshi.

### Integration (`nfl_draft/tests/integration/test_dashboard_queries.py`)

- Audit `test_cross_book_grid_query_outlier_flags` (around line 57) to
  confirm the asserted flag direction is "below median" for
  draftkings. If the fixture has DK above the median, either invert
  the fixture values or update the expected `outlier_count` and
  `flags` values.
- Confirm the Kalshi-specific integration tests
  (`test_cross_book_grid_kalshi_flag_uses_implied_prob_not_devig_prob`,
  `test_cross_book_grid_kalshi_flag_fires_on_big_take_edge`,
  `test_cross_book_grid_kalshi_flag_suppressed_when_no_take_price`)
  still pass unchanged; they already cover Kalshi-side semantics
  accurately.

## Documentation

- `nfl_draft/README.md` — the Cross-Book Grid section (rewritten
  recently for actual-price display) documents flag semantics. Add a
  paragraph describing the asymmetric rule and its motivation (YES-only
  sportsbook futures vs. two-sided Kalshi binaries).
- `nfl_draft/CLAUDE.md` — if it currently says anything about outlier
  flagging semantics, update it. If not, no change.

## Version control

- **Branch:** `feature/asymmetric-outlier-flags`
- **Worktree:** create a worktree off main before any edits
  (`/Users/callancapitolo/NFLWork-<branch>` or wherever the standard
  location is). Do not symlink the DuckDB file — copy it if the tests
  need a populated database; the unit tests use a fresh temp DB
  fixture, which is preferred.
- **Files touched:**
  - `nfl_draft/lib/queries.py` — flag rule + docstring.
  - `kalshi_draft/app.py` — grid blurb copy only.
  - `nfl_draft/tests/unit/test_cross_book_grid.py` — update + add tests.
  - `nfl_draft/tests/integration/test_dashboard_queries.py` — audit /
    adjust.
  - `nfl_draft/README.md` — doc paragraph.
- **Commit structure:** a single commit is acceptable given the small
  footprint. Commit message should reference the motivation (bettable
  edges only on non-Kalshi venues) and not just the diff.
- **Merge:** pre-merge executive review per CLAUDE.md, then explicit
  user approval before merging to `main`. Clean up worktree + branch
  after merge.

## Non-goals

- No change to Kalshi flag semantics.
- No config/per-book abstraction (Approach B rejected as YAGNI). If
  a future sportsbook ships two-sided NFL draft markets, we revisit.
- No change to the signed delta math in `ev_candidates`; consumers
  continue to see `book_prob - median`. For non-Kalshi rows this will
  always be negative after the change, which is the natural
  consequence of the new rule, not a separately-designed property.
- No dashboard-layout changes beyond the copy update.
