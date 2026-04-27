# MLB Combined Parlay — Design Spec

**Date:** 2026-04-26
**Status:** Approved (pending user spec review)
**Owner:** Callan

## Summary

Add a feature to the MLB Dashboard's Parlay tab that lets the user select two existing recommended parlay rows (each a 2-leg spread+total combo for a single game) and combine them into one 4-leg Wagerzon ticket. The dashboard prices the combined ticket via Wagerzon's `ConfirmWagerHelper` API, recommends a stake using Kelly on the joint distribution, and dynamically updates the source rows' Kelly recommendations to the **conditional residuals** so the user can optionally top up.

## Problem & Motivation

Wagerzon balance is the binding constraint on bet sizing. Today, placing both single-game parlays as separate WZ tickets ties up the full Kelly stake of each (e.g., $200 + $150 = $350). Combining them into a single 4-leg ticket lets the user achieve some exposure to both events for far less cash (e.g., $50 for the combined Kelly stake) — at the cost of compounded WZ vig and lower absolute EV.

The feature exists **purely as a cash-efficiency tool**. For independent cross-game legs, combining doesn't capture additional EV the singles wouldn't already provide; it lets the user express the same edges with less deposit.

## Scope

### In scope (v1)
- Combine **2 rows** from the existing Parlay tab into a 4-leg WZ ticket
- Live `ConfirmWagerHelper` pricing for the combined ticket
- Kelly recommendation for the combined stake using existing `parlay_bankroll` / `parlay_kelly_mult`
- **Conditional Kelly residuals** on source rows after the combo is placed
- Source rows stay active with reduced Kelly; user can place residual single bets if they want
- Persistence: combo + leg relationships survive page refreshes

### Out of scope (v1, considered for later)
- N>2 leg combinations
- Combining a parlay with a single from the Bets tab
- Combining two singles from the Bets tab
- Cash-budget-aware portfolio optimization (set "available WZ balance" → dashboard solves)
- Same-game two-row combinations (each game's spread+total is already a single Parlay row)
- "Top up to original Kelly" hint after singles are placed first
- Auto-place via WZ API (we still record intent in `placed_parlays`; user manually places at WZ)

## User Story

1. User is on the Parlay tab. Multiple +EV recommended parlays are listed.
2. User checks two rows (different games).
3. A slim **Combined Parlay banner** appears above the table, showing:
   - The two combined legs
   - Joint fair odds (= product of `fair_dec`) and exact WZ payout (from `ConfirmWagerHelper`)
   - Joint edge (%)
   - Recommended Kelly stake for the combined ticket
   - "Place Combined →" button
4. User clicks "Place Combined".
5. Combined ticket is recorded in `placed_parlays` with leg references.
6. Source rows stay visible. Their **Kelly column updates** to the conditional residual (e.g., $100 → $85). A small annotation reads `(reduced from $100 — combo #N placed)`.
7. User can place either source row's residual single, or stop.

## UI Changes

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard.R` (Parlay tab rendering, lines ~245–400)

### Parlay table
- New leftmost **checkbox column** ("Sel"), one checkbox per row
- Existing Kelly column gains conditional rendering: shows residual amount + annotation when a combo is placed against this row

### Combined Parlay banner
- Sticky div at the top of the Parlay tab, hidden by default
- Appears when exactly 2 rows are checked
- Maximum 2 rows can be checked at once; checking a 3rd row clears the oldest selection (FIFO)
- Layout (single line):
  ```
  ▍2 legs combined · stake $X · WZ pays $Y · edge +Z%   [Place combined →]
  ```
- States:
  - **Loading** — "Pricing combined ticket…" while waiting for `ConfirmWagerHelper`
  - **Ready** — full info as above
  - **Error** — "WZ rejected: <reason>" with no Place button
  - **Negative edge** — same layout, edge displayed in red (no special block)

### Edge-case interactions
- If either checked row's `game_id` matches the other → checkboxes refuse to allow both (tooltip: "Same-game combos are already a single row")
- If either source row is already in `placed_bets` or `placed_parlays` as a single → checkbox is disabled (tooltip: "Already placed; combine before placing singles")

## Sizing Math

### Combined stake (Kelly on the combined ticket)

Given the two rows:
- Joint fair probability: `p_joint = (1/fair_dec_A) * (1/fair_dec_B)` (independence across games)
- WZ payout decimal `d_wz`: returned by `ConfirmWagerHelper` (4 plays)
- Joint edge: `p_joint * d_wz - 1`
- Kelly fraction: `(p_joint * d_wz - 1) / (d_wz - 1)`
- Recommended stake: `kelly_fraction * parlay_bankroll * parlay_kelly_mult`

### Conditional residuals (after combo is placed)

Given the combo stake `s_C` is fixed, find `s_A` and `s_B` that maximize expected log-growth over the 4-outcome joint distribution (both win, A only, B only, neither). Implementation:

1. Generate a synthetic outcome matrix with one row per outcome:
   - Columns: `[s_A_payoff, s_B_payoff, s_C_payoff]`
   - 4 rows, with probabilities `p_A * p_B`, `p_A * (1-p_B)`, `(1-p_A) * p_B`, `(1-p_A) * (1-p_B)`
2. Pass to extended `parlay_multivariate_kelly` with `s_C` fixed
3. Solver returns optimal `s_A`, `s_B` given `s_C` constraint
4. Floor at $0; if residual ≤ minimum bet threshold, mark row as fully reduced

The math reuses the same machinery as the existing same-game multivariate Kelly. The new piece is the "fix one position, solve the rest" constraint — a small extension to the existing solver (~30 lines).

## Wagerzon Pricing

**File:** `wagerzon_odds/parlay_pricer.py`

- Existing `ConfirmWagerHelper` flow (4 plays per game) extends to 4-leg parlays naturally
- Plays for combined ticket: `[home/away_spread_A, over/under_total_A, home/away_spread_B, over/under_total_B]`
- Same `RiskWin="2"` preview-mode call
- Same `MAXPARLAYRISKEXCEED` fallback (try $100 if $10,000 fails)
- Returned exact payout becomes `d_wz` for the sizing math above

**Caching:** Server caches priced combos by `(row_a_id, row_b_id)` for the duration of the session, so checking back to a previously priced combo is instant.

## Data Model Changes

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard.duckdb`

### `placed_parlays` (existing table — extend)
Add columns:
- `is_combo` BOOLEAN (default FALSE) — true if this row is a 4-leg combined ticket
- `combo_leg_ids` TEXT — JSON array of `parlay_hash` values referenced by this combo (NULL when `is_combo = FALSE`)
- `parent_combo_id` INTEGER — for source rows, FK back to the parent combo row (NULL when not in a combo)

### Behavior
- Placing a combo: insert one row with `is_combo = TRUE`, `combo_leg_ids = [hash_a, hash_b]`, no `parent_combo_id`
- Source rows are **not** inserted into `placed_parlays` when only the combo is placed — they remain available to the dashboard. Their reduced Kelly is computed on the fly from `mlb_dashboard.duckdb` queries.
- If user later places one source row's residual single, that placement gets a normal `placed_parlays` row with `parent_combo_id` pointing back to the combo (so we know it's a top-up, not a duplicate)
- Unchecking the source rows in the UI **after** placing a combo does not remove the combo from `placed_parlays` — the residual Kelly stays in effect. To remove a combo, use the existing parlay-removal flow on the combo row in the placed bets summary.

## Code Changes (file-level)

| File | Change |
|---|---|
| `Answer Keys/Tools.R` | Add `compute_combined_parlay_pricing()` for joint fair odds + edge across two existing rows |
| `Answer Keys/mlb_correlated_parlay.R` | Extend `parlay_multivariate_kelly()` to accept a fixed-position vector for conditional Kelly |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Add checkbox column + Combined Parlay banner; conditional Kelly column rendering |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | New endpoints: `POST /api/price-combined-parlay`, `POST /api/place-combined-parlay`; modify residual Kelly computation in existing endpoints |
| `Answer Keys/MLB Dashboard/mlb_dashboard.duckdb` | Schema migration: add columns to `placed_parlays` |
| `wagerzon_odds/parlay_pricer.py` | Refactor `ConfirmWagerHelper` call to accept arbitrary play count (currently hardcoded to 2 games); add 4-leg helper |

## Error Handling

| Error | Behavior |
|---|---|
| `ConfirmWagerHelper` returns `MAXPARLAYRISKEXCEED` at $10k | Fall back to $100 stake (existing pattern) |
| `ConfirmWagerHelper` returns other error | Banner shows "WZ pricing unavailable: <code>" — no Place button |
| Network timeout on price call (>5s) | Banner shows "Pricing timed out — retry" with retry button |
| User checks 2 rows from same game | Checkboxes refuse second selection; tooltip explanation |
| User checks rows where one is already placed | Checkbox disabled; tooltip explanation |
| Conditional Kelly residual computes to negative | Floor at $0; row's Kelly column shows "0 (covered by combo)" |
| Joint edge is negative | Banner displays edge in red; Place button still enabled (per user decision: don't block) |
| Combo placement persistence fails | Roll back; show toast "Failed to record combo — try again" |

## Testing

### Unit tests (Tools.R)
- `compute_combined_parlay_pricing()`: joint fair odds matches product of `fair_dec`
- Conditional Kelly: synthetic 50/50 example with known closed-form answer (compare numerical result to expected within tolerance)
- Conditional Kelly: residual ≤ original Kelly always

### Integration tests
- End-to-end: insert mock parlay rows in `mlb_dashboard.duckdb`, hit `POST /api/price-combined-parlay`, assert response shape + math
- Place combo, fetch Parlay tab data, assert source rows show reduced Kelly

### Manual UI tests
- Check 2 rows → banner appears with correct numbers
- Check 3rd row → first selection clears (only 2 allowed at once)
- Place combo → source rows update Kelly, banner clears
- Refresh page → combo persists, source row reductions persist
- Try to combo same-game rows → blocked
- Try to combo a row that's already placed → blocked

## Documentation Updates

| File | Update |
|---|---|
| `Answer Keys/MLB Dashboard/README.md` | New section "Combined Parlay" — describe banner UX, conditional Kelly residuals, when to use (cash efficiency), when not to (unconstrained cash) |
| `CLAUDE.md` (repo root) | No changes needed — no architectural shifts |

Documentation goes in the same merge commit as the feature, per repo discipline.

## Branching & Worktree Plan

### Worktree
- Path: `/Users/callancapitolo/NFLWork-combined-parlay`
- Created via `/worktree` skill before any code changes
- DuckDB databases (`mlb.duckdb`, `mlb_dashboard.duckdb`) are **copied** into the worktree, not symlinked (per CLAUDE.md rule)
- Test from the worktree using a separate Flask port (e.g., 8084) to avoid colliding with running dashboard on 8083

### Branch
- Name: `feature/mlb-combined-parlay`
- Base: `main`

### Commit structure
1. **`feat(answer-keys): conditional Kelly solver`** — extend `parlay_multivariate_kelly` for fixed-position constraint + unit tests
2. **`feat(answer-keys): combined parlay pricing helper`** — `compute_combined_parlay_pricing` in Tools.R
3. **`feat(wagerzon): N-leg ConfirmWagerHelper`** — refactor `parlay_pricer.py` for arbitrary play count
4. **`feat(mlb-dashboard): combined parlay schema`** — DuckDB migration for `placed_parlays`
5. **`feat(mlb-dashboard): combined parlay backend`** — Flask endpoints
6. **`feat(mlb-dashboard): combined parlay UI`** — checkbox column + banner + residual rendering
7. **`docs(mlb-dashboard): document combined parlay`** — README updates

### Pre-merge review (per CLAUDE.md)
- Executive engineer review of `git diff main..HEAD`
- Specifically check: deduplication of combo placements, residual Kelly never negative, WZ rejection handling, schema migration is reversible, no stale combo references when source rows update
- Document findings as ISSUES TO FIX vs ACCEPTABLE RISKS
- Get explicit user approval before merging to `main`

### Cleanup
- After merge: `git worktree remove /Users/callancapitolo/NFLWork-combined-parlay`
- Delete branch: `git branch -d feature/mlb-combined-parlay`

## Open Questions / Followups (post-v1)

- Should we eventually surface a "WZ available cash" input and switch to budget-aware portfolio Kelly? (Option 3 from brainstorm — most powerful, requires UX add)
- Should we extend to N>2 legs? (Variance grows fast; probably diminishing returns)
- Should we extend to combine a parlay with a single from the Bets tab? (Cleanest first version stays within Parlay tab)
- Do we want CLV tracking to decompose combo outcomes into per-leg results? (Would help track which legs contributed to combo wins/losses over time)
