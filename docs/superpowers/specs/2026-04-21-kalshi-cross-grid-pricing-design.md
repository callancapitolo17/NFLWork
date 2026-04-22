# Kalshi Cross-Book Grid Pricing — Design

**Date:** 2026-04-21
**Scope:** Cross-Book Grid only (dashboard Tab 1). EV Candidates tab is out of scope.

## Problem

The Cross-Book Grid shows `pick_5_overall_carnell-tate` at **Kalshi = 2.0%** while Kalshi's own UI reports Tate at **4% chance** with `Buy 5¢ / Sell 2¢`. Same problem drives false outlier flags on `top_5` (8% vs median ~16%) and `top_10` (75% vs median ~89%).

**Root cause:** the scraper (`nfl_draft/scrapers/kalshi.py:53-72, 135-138`) uses **`yes_bid`** — the sell-side price — as Kalshi's implied probability. Every Kalshi cell in the grid is therefore the sell price, which is systematically below fair value by half the bid-ask spread.

Secondary bug discovered during investigation: `_write_legacy_kalshi_odds` (`kalshi.py:256-265`) uses `m.get("yes_bid", 0)` without handling Kalshi's current `yes_bid_dollars: "0.02"` response format. All five price fields (`yes_bid`, `yes_ask`, `no_bid`, `no_ask`, `last_price`) are therefore persisted as zeros in the `kalshi_odds` table, which also breaks the legacy Kalshi dashboard tabs (Price History / Edge Detection / Consensus).

## Principle

The Cross-Book Grid is a **fair-value comparison** across venues. Sportsbook cells show each book's fair-value estimate (vig-stripped implied probability). For internal consistency, Kalshi's cell should also be a fair-value estimate, not a tradeable price.

A tradeable price is relevant for **outlier flagging** — the question "can I actually exploit this disagreement?" is about what you'd pay to enter a position, not about the venue's opinion of fair. We can separate these two concerns cleanly by using Kalshi's mid for fair value and Kalshi's buy price for the flag.

## Design

### Per-venue semantics

| Schema column | Sportsbooks | Kalshi |
|---|---|---|
| `implied_prob` | Raw implied from American odds (includes vig) | **Buy price** (`yes_ask`) / 100 |
| `devig_prob` | Vig stripped → fair-value estimate | **Mid** = (`yes_bid` + `yes_ask`) / 200 |

`devig_prob` retains its meaning across venues ("this book's fair-value estimate"). `implied_prob` retains its meaning too ("raw take price — includes the venue's friction"). For sportsbooks the friction is vig; for Kalshi the friction is half the bid-ask spread. Both are already-present columns — no schema migration.

### Cell display and median

Grid renders **`devig_prob * 100`** for every venue including Kalshi — unchanged from current behavior. Median math continues to run on `devig_prob` across all posting venues. Kalshi's mid participates in the median on equal footing with each sportsbook's devigged probability.

### Outlier flag — per-venue logic

| Venue | Flag formula |
|---|---|
| DK, FD, BM, WZ, Hoop88 | `|devig_prob − median| ≥ threshold_pp` *(unchanged)* |
| Kalshi | **`|implied_prob − median| ≥ threshold_pp`** *(new)* |

The asymmetry reflects real structural differences: sportsbook raw prices are inflated by 5-8% vig, so raw-vs-fair flags would be all noise for them. Kalshi's spread is typically 1-3¢, so Kalshi's buy price is close enough to fair that `|buy − median|` yields clean edge signal.

### Resolution cascades

**`devig_prob` (fair value)** — populated by the scraper in this order:
1. `(yes_bid + yes_ask) / 200` — if both sides posted
2. `last_price / 100` — if only one side posted
3. The single posted side / 100 — if `last_price` is also absent
4. No row is emitted — the market has no signal at all

**`implied_prob` (take price)** — populated in this order:
1. `yes_ask / 100` — if posted
2. `last_price / 100` — fallback when no ask
3. Left `NULL` — Kalshi cell renders fair value but outlier flag is suppressed for that row

### Tooltip

Hovering a Kalshi cell shows:

```
Last trade:  4¢
Buy:  5¢     Sell:  2¢
```

Terminology is **Buy price** and **Sell price** throughout the UI — no `bid` / `ask` jargon. Tooltip values come from `kalshi_odds` via a join keyed on `market_map.book_label + book_subject`.

Staleness handling for `last_price`: use the value from the latest `kalshi_odds` snapshot as-is. The 2-hour `MAX_AGE_HOURS` filter at `nfl_draft/lib/queries.py:42` already drops whole rows when the scrape is stale, so if we render a Kalshi row at all, its snapshot is fresh — even if the embedded `last_price` field was recorded at an earlier trade. If this proves misleading in practice, v2 upgrades to joining `kalshi_trades` for an actual `MAX(traded_at)` gate.

## Implementation Surface

Three modules change. No schema migration.

### 1. `nfl_draft/scrapers/kalshi.py`

- Add `_extract_yes_ask_cents`, `_extract_no_bid_cents`, `_extract_no_ask_cents`, `_extract_last_price_cents` — mirrors of the existing `_extract_yes_bid_cents`, each tolerating both the legacy integer and new `*_dollars` response formats.
- Rewrite `parse_markets_response` to:
  - Compute `yes_bid_cents` and `yes_ask_cents`
  - Derive `fair_cents` per cascade above
  - Derive `take_cents` per cascade above
  - Emit `OddsRow` with `american_odds` reflecting `take_cents` (so `implied_prob` downstream resolves to the buy price), and carry `fair_cents` through so `devig_prob` resolves to the mid
  - This likely requires a small change to `scrapers/_base.py::OddsRow` or to how `normalize.py` populates `devig_prob` — see Open Question #1.
- Rewrite `_write_legacy_kalshi_odds` to use the `_extract_*_cents` helpers for all five price fields. This is the side-effect fix that restores the legacy `kalshi_odds` table to live data.

### 2. `nfl_draft/lib/queries.py::cross_book_grid`

Extend the CTE to pull `implied_prob` alongside `devig_prob`:

```sql
WITH latest AS (
  SELECT market_id, book, implied_prob, devig_prob,
         ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
  FROM draft_odds
  WHERE fetched_at > NOW() - INTERVAL 'MAX_AGE_HOURS hours'
)
SELECT market_id, book, implied_prob, devig_prob FROM latest WHERE rn = 1
```

In the aggregation loop, store both numbers per `(market_id, book)`. Flag logic becomes:

```python
for book, row in books.items():
    if book == "kalshi":
        take_price = row["implied_prob"]
        if take_price is None:
            flag = False  # no buy price and no last trade → no actionable signal
        else:
            flag = abs(take_price - median) >= threshold
    else:
        flag = abs(row["devig_prob"] - median) >= threshold
```

Median is still computed from `devig_prob` across books. If a Kalshi row is emitted with valid `devig_prob` but NULL `implied_prob` (one-sided market with no last trade), it still contributes to the median but the flag falls through to `False` — the cell appears in the grid without an outlier mark.

### 3. Dashboard render layer

Cell value continues to be `devig_prob * 100` — no visual change to the number when Kalshi's price matches its mid. Add a tooltip renderer that, for Kalshi cells only, pulls `{ticker, yes_bid, yes_ask, last_price}` for every Kalshi market visible on the grid (batched, not per-hover):

```sql
WITH latest_ko AS (
  SELECT ticker, candidate, yes_bid, yes_ask, last_price,
         ROW_NUMBER() OVER (PARTITION BY ticker ORDER BY fetch_time DESC) AS rn
  FROM kalshi_odds
)
SELECT mm.market_id, ko.ticker, ko.yes_bid, ko.yes_ask, ko.last_price
FROM market_map mm
LEFT JOIN latest_ko ko
  ON ko.ticker LIKE mm.book_label || '-%'
 AND ko.candidate = mm.book_subject
 AND ko.rn = 1
WHERE mm.book = 'kalshi'
```

Results are keyed by `market_id` and handed to the dashboard's row renderer as a lookup dict. Hovering a cell reads from the dict — no per-hover database call. The join key is `(book_label, book_subject)` because `market_map` stores `book_label` as a ticker prefix (e.g. `KXNFLDRAFTPICK-26-5`) and the full ticker is that prefix plus a candidate shortcode (`KXNFLDRAFTPICK-26-5-CTAT`). `LIKE prefix || '-%'` combined with candidate-name equality is unique per `(label, candidate)`.

The scrape code path that populates `market_map` and `kalshi_odds` already runs on every scrape tick, so this query is cheap (tens of ms at the full 500+ Kalshi markets) and the tooltip data is never more stale than the grid data.

## Expected Outcomes

**Tate example before/after** — `pick_5_overall_carnell-tate` with DK 7.7%, BM 6.2%, WZ 6.7% on the row:

| | Before | After |
|---|---|---|
| Kalshi cell value | 2.0% (sell price) | 3.5% (mid) |
| Kalshi tooltip | none | `Last 4¢ · Buy 5¢ · Sell 2¢` |
| Grid median (4 books: sort → average middle two) | 6.45% — median of [2.0, 6.2, 6.7, 7.7] | 6.45% — median of [3.5, 6.2, 6.7, 7.7] |
| Kalshi outlier flag at threshold 8pp | flagged — `|2.0 − 6.45| = 4.45pp` at threshold 5pp triggers, false positive | not flagged — buy price 5% vs median 6.45% = `1.45pp`, correct |

The median is unchanged in this case because Kalshi moves from 2% to 3.5%, both below the middle-two of the sorted list. In markets where Kalshi sits near the middle of the distribution, the move from bid to mid *will* shift the median — expect small drift on many grid rows after the switch.

**`top_10_carnell-tate`** (currently 75% flagged vs median 88.6%):
- Cell value moves from 75% (bid) to mid of buy/sell, likely ~80%
- Flag recomputes against buy price (probably ~85%) vs updated median — likely still in sync with consensus, no false flag

**Side effect:** legacy `kalshi_odds` table gets populated with real numbers for the first time since Kalshi migrated to the `*_dollars` response format. Legacy Kalshi dashboard tabs (Price History / Edge Detection / Consensus) recover.

## Non-Goals

- **EV Candidates tab** — out of scope for this spec. Its math logic and what Kalshi price it uses for +EV detection are a separate design discussion.
- **Kalshi No-side pricing** — the grid is Yes-only. `no_bid` and `no_ask` are persisted (so the legacy write fix covers them) but not rendered.
- **Bet sizing** — no changes to Kelly calc, stake suggestions, or liquidity gating.
- **Grid threshold defaults** — the default outlier threshold (`threshold_pp`) doesn't change. Users may want to tune it after seeing the new numbers, but that's a post-merge observation, not a design decision.

## Open Questions

1. **Where does `devig_prob = mid` get computed — in the scraper or in `lib/normalize.py`?** The current pipeline has `normalize.py` running over a set of `OddsRow`s per market to compute `devig_prob`. For Kalshi we want to **bypass** that group-normalization (it's vig-stripping logic written for sportsbooks) and pass the mid straight through. Two candidate paths:

   - (a) Add a nullable `precomputed_devig` field to `OddsRow`. `normalize.py` checks it — if present, emit as-is; if absent, run the existing per-group devig. Minimally invasive.
   - (b) Add a `per_book_normalizer` registry that dispatches by `book`. `normalize.py` becomes a thin dispatcher; Kalshi gets a passthrough normalizer. More principled but larger refactor.

   Recommend (a) for this spec; evaluate (b) later if more venues need custom treatment. Final choice is a plan-phase decision.

## Version Control

- **Branch:** `feature/kalshi-cross-grid-pricing` (already created)
- **Worktree:** brainstorm + design on `main`; implementation plan will use a worktree per `/feedback_worktree_planning`.
- **Files expected to change:**
  - `nfl_draft/scrapers/kalshi.py` (mods)
  - `nfl_draft/scrapers/_base.py` or `nfl_draft/lib/normalize.py` (mods — TBD per Open Question #1)
  - `nfl_draft/lib/queries.py` (mods to `cross_book_grid`)
  - Dashboard module rendering the grid (mods — tooltip rendering + per-book flag interpretation)
  - `nfl_draft/README.md` (docs — semantics of Kalshi cell + tooltip)
  - `CLAUDE.md` (docs — if Kalshi cell semantics is worth surfacing in project context)
  - Tests: extend `tests/unit/test_scraper_parsing.py` + `tests/integration/test_dashboard_queries.py`

## Documentation

- `nfl_draft/README.md` — new subsection under the Cross-Book Grid explanation documenting that Kalshi's cell is the mid of buy/sell, that outlier flags compare buy price to median, and that tooltip reveals raw buy/sell/last.
- No `CLAUDE.md` update required unless the semantic split (mid-for-display / buy-for-flag) becomes a cross-cutting pattern we'd reapply to another venue.

## Pre-merge Review Hooks

Per CLAUDE.md executive-engineer-review checklist, the plan phase will verify:
- No duplicate `OddsRow` emission when both bid and ask are present
- Off-season behavior: Kalshi returns no markets → scraper exits cleanly with zero rows, grid query returns `[]`, dashboard handles empty gracefully
- Stale-snapshot path: if scrape hasn't run for > 2h, `MAX_AGE_HOURS` filter drops the Kalshi row from the grid; tooltip query returns no rows → tooltip hidden rather than showing stale numbers
- Tooltip query performance at 500+ visible markets (whole Kalshi PICK series)
- Test coverage for the `yes_ask` zero / `last_price` zero / all-three-zero fallback chain
