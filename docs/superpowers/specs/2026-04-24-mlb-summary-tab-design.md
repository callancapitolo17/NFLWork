# MLB Summary Tab — Design

**Date:** 2026-04-24
**Status:** Draft, pending user review
**Owner:** Callan

## Problem

We place MLB correlated parlays (Spread + Total combos for full game and F5) on Wagerzon, sourced from the model in `Answer Keys/mlb_correlated_parlay.R`. Today there is no clean way to see how those bets are performing. The existing `CBB Summary` tab in the bet logging Google Sheet (built by `bet_logger/create_summary.py`) is a useful template but is also brittle — fragile description regex, a hidden helper sheet that occasionally throws `#SPILL!`, and a hardcoded season end date.

We want a new tab that tracks overall performance of MLB correlated parlays, auto-updates as bets land in `Sheet1`, and avoids the failure modes of the CBB tab.

## Goal

Add a single new tab — `MLB Summary` — to the bet logging Google Sheet (ID `1t9_7HmsrQAu34HI_gMIDPErC2kThz2MHBwwzZIYU6H4`) that auto-updates from `Sheet1` and shows aggregate performance of MLB correlated parlays placed on Wagerzon.

This is **Phase B**. Phase A (joining each placed bet back to `mlb_parlay_opportunities` in `mlb.duckdb` to attach the model's edge% and fair odds at placement) is out of scope here and will be a follow-on spec.

## Non-Goals

- No per-bet detail list. The user explicitly only wants aggregates.
- No combo-type breakdown (Home Fav+Over, Away Dog+Under, etc.). That requires the join to `mlb_parlay_opportunities` and is Phase A.
- No CLV, no edge-vs-realized-ROI, no model calibration. Phase A.
- No new scraper logic. The existing `bet_logger/scraper_wagerzon.py` already captures everything needed.
- No DuckDB read/write, no scheduled job, no ongoing Python script. The tab's contents are pure spreadsheet formulas. (A one-time Python setup script — covered in *Files Touched* — writes those formulas to the tab once; from then on, nothing runs.)

## Approach

**One new tab built entirely from spreadsheet formulas referencing `Sheet1`.** No hidden helper sheet. No Python script. No scheduled job.

Why this shape:
- The user wants the tab to update live as bets are added, with zero manual steps. Formulas give that.
- The filter for "is this row an MLB correlated parlay" is simple enough that formulas can express it cleanly without the regex gymnastics that broke the CBB tab.
- Phase A (joining to the model's recommendation table) will eventually need a Python script anyway. When that ships, it will write Phase-A-only columns to a separate region of the same tab; the Phase B formulas keep working untouched.

### Filter — what counts as a tracked bet

A `Sheet1` row qualifies as an MLB correlated parlay iff **all** of:

1. Column C (Sport) = `"MLB"`
2. Column E (Bet Type) = `"Parlay"`
3. Column D (Description) starts with `"PARLAY (2 TEAMS)"` — enforces exactly 2 legs and excludes 3+ leg parlays
4. Column D (Description) contains a spread token, matched by regex roughly `[+\-]\d+(\.\d+)?` near a team-name fragment (specific regex tuned during implementation against real Wagerzon descriptions)
5. Column D (Description) contains a total token, matched by regex roughly `(Over|Under|^O |^U )\s*\d+`

Filter rationale:
- The user confirmed: "any spread+total parlay, FG or F5, counts."
- 3+ leg parlays are excluded by the header check, since the Wagerzon scraper sets `HeaderDesc` based on leg count.
- 2-leg parlays that are *not* spread+total (e.g., ML+ML, total+total) naturally fail one of the regex checks and are excluded — consistent with what the user wants.

### FG vs F5 classification

Within the qualifying set, a bet is **F5** if Description matches `(?i)1st 5|F5|First 5|H1`, otherwise **FG**. Used only for the FG/F5 split block; everything else aggregates across both.

(The exact F5 regex will be confirmed against real Wagerzon F5 leg descriptions before implementation. If F5 markers turn out to vary, the regex gets widened — this is a small, contained risk.)

### Settlement treatment

The Wagerzon scraper already normalizes results to `win` / `loss` / `push` (with `cancelled`, `void`, `no action` mapping to `push`). Treatment in the tab:

- **Win**: counted in settled, contributes positive P&L (column J payout − column H stake), counts as W.
- **Loss**: counted in settled, contributes negative P&L (− column H stake), counts as L.
- **Push** (incl. void/cancelled/no action): counted in settled, P&L = $0, counts as P. Excluded from win-rate denominator.
- **Pending / unsettled**: not displayed (user does not want a pending bets metric).

## Tab Layout

Three blocks, top to bottom in the new `MLB Summary` tab.

### Block 1 — Headline stats (overall)

Single column of metrics:

| Metric | Formula sketch |
|---|---|
| Total bets placed | `SUMPRODUCT(filter_mask)` |
| Settled bets (W+L+P) | `SUMPRODUCT(filter_mask * settled_mask)` |
| Total wagered | `SUMPRODUCT(filter_mask * settled_mask * Sheet1!H2:H)` |
| Net P&L | `SUMPRODUCT(filter_mask * win_mask * (J − H)) − SUMPRODUCT(filter_mask * loss_mask * H)` |
| ROI | `Net P&L / Total wagered` |
| Win rate | `W / (W + L)` |
| Record | `"W-L-P"` (concatenation of the three counts) |
| Avg decimal odds | `AVERAGE` of column I across the qualifying set |

`filter_mask` is the inline product of the five filter conditions; it is repeated in each formula rather than precomputed in a helper sheet.

### Block 2 — FG vs F5 split

Two-row table with the same metrics as Block 1, one row labeled `FG` and one labeled `F5`. Same formulas as Block 1 with one extra factor in `filter_mask`: `NOT(REGEXMATCH(D, F5_regex))` for FG, `REGEXMATCH(D, F5_regex)` for F5.

### Block 3 — Trend over time

Two charts, with their backing data tables sitting below them on the same tab. Charts are the primary view; tables exist because charts need a cell range as a data source and because sometimes you want to look up an exact number.

**Chart 1 — Equity curve (cumulative P&L by day).** Line chart. X-axis: date. Y-axis: cumulative P&L. The most important chart for evaluating a betting strategy — its endpoint matches the headline Net P&L exactly. Backed by the daily P&L table's Date and Cumulative P&L columns.

**Chart 2 — Weekly P&L (bar chart).** Vertical bars, one per week, colored green for positive weeks and red for negative weeks. Lets the user spot streaks and identify which weeks dragged the bottom line. Backed by the weekly P&L table's Week-of and P&L columns.

**Daily P&L table** (placed below the charts; feeds the equity curve):

| Date | # Bets | Wagered | P&L | Cumulative P&L |
|---|---|---|---|---|

The Date column is dynamic:
```
=SORT(UNIQUE(FILTER(Sheet1!A:A, filter_mask_arrayformula)))
```
This auto-extends as bets on new dates appear. Other columns use `SUMPRODUCT` against the date in column A of the row. The chart's data range is set to a generous block (e.g., 1000 rows) so it picks up new dates automatically as the column extends.

**Weekly P&L table** (placed below the charts; feeds the bar chart): identical structure but grouped by week-start (using `A − WEEKDAY(A, 2) + 1` as the bucket key). Same dynamic auto-extension and same generous-range trick for the chart.

## Failure Modes & How They're Handled

- **Wagerzon changes its description format.** Filter would silently start excluding bets. Mitigation: Block 1's "Total bets placed" gives a sanity check the user can eyeball against their own memory. (Phase A's join will catch this much earlier — placed bets that don't match any recommendation will surface as unmatched.)
- **`#SPILL!` from `SORT(UNIQUE(FILTER(...)))`.** The CBB tab hit this because helper-sheet formulas overlapped with each other and with manually-edited rows. We avoid it by putting the trend tables far enough down the tab that they have unbounded vertical room, and by keeping Block 1 and Block 2 fixed-cell (no array spill).
- **A 2-leg total+total parlay slips through.** Possible if the total regex is too permissive. The header check (`PARLAY (2 TEAMS)`) bounds the false-positive rate. Worst case: the user adds a description-text exclusion later.
- **Wagerzon scraper voids a leg mid-game.** Already normalized to `push` upstream. Treated as $0 P&L.

## Edge Cases — Defaults Locked

- Push / void / cancelled / no action → settled, P&L = $0, excluded from win-rate denominator.
- 2-leg totals-only parlay → excluded by filter.
- 3+ leg parlay (even one with a spread+total inside) → excluded by the `PARLAY (2 TEAMS)` header check.

## Out of Scope (Phase A and Beyond)

These are intentionally not in this design and will be addressed in a follow-on spec:

- Joining each placed parlay to `mlb_parlay_opportunities` in `mlb.duckdb` via game + combo + line, attaching the model's fair odds, edge%, and combo classification at placement.
- Combo-type breakdowns (Home Fav+Over, Away Dog+Under, etc.).
- Realized vs expected ROI calibration.
- CLV tracking.
- Per-bet detail list.

## Files Touched / Created

- **New file:** `bet_logger/create_mlb_summary.py` — one-shot script that creates the `MLB Summary` tab, writes the formulas, and creates the two embedded charts (equity curve + weekly bars) via the Sheets API's `addChart` request, binding each chart to its backing data table. Run once at setup; from then on the tab and its charts are self-maintaining (charts auto-redraw whenever their data range changes). Modeled on the structure of `bet_logger/create_summary.py` but writes formulas only — no hidden helper sheet.
- **New tab in Google Sheet:** `MLB Summary` (sheet ID `1t9_7HmsrQAu34HI_gMIDPErC2kThz2MHBwwzZIYU6H4`).
- **No changes to** `scraper_wagerzon.py`, `sheets.py`, `create_summary.py`, `run_all_scrapers.sh`, or any DuckDB layer.

## Version Control

- **Branch:** `feature/mlb-summary-tab`
- **Worktree:** create at `../NFLWork-mlb-summary` before any code changes; remove after merge.
- **Commits:** single commit for the new script + a README update; one commit at most for any follow-up fixes from running the script against the live sheet.
- **Merge:** to `main` only after the script has been run against the live sheet and the user has eyeballed the resulting tab.

## Documentation

- **Update:** `bet_logger/README.md` — add a section describing the new `MLB Summary` tab, what bets it tracks, and how to run `create_mlb_summary.py` for initial setup.
- **No CLAUDE.md changes needed** — this doesn't introduce new architecture, just a new artifact in the existing bet-logger pattern.

## Open Questions / Risks

1. **Exact F5 leg-description format.** Need to confirm against real Wagerzon F5 bets before finalizing the F5 regex. If no F5 bet has been placed yet, the regex stays best-guess and the FG/F5 split goes 100% FG until the first F5 bet lands and we tune.
2. **Spread regex precision.** Wagerzon descriptions contain various numeric tokens (odds, lines, fractions like `½`). The spread-detection regex will need to be tuned against real bets, not designed in isolation. This is an implementation-time concern, not a design-time one.
