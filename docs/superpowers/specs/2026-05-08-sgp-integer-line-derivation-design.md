# SGP Integer-Line Derivation Design

**Date:** 2026-05-08
**Author:** callancapitolo + Claude
**Status:** Approved, ready for implementation plan
**Branch:** `feature/sgp-integer-line-derivation`
**Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation`

---

## 1. Problem

The MLB SGP comparison pipeline writes a "fair probability" per book per combo to `mlb_sgp_odds`. The R blender (`Answer Keys/mlb_correlated_parlay.R`) averages model + DK + FD + PX + NV fair probs into `blended_prob`, then computes edge vs Wagerzon's offered SGP price.

When Wagerzon posts an SGP on an integer total line (e.g., FG Over 8, F5 Over 4) and a comparison book's alt-total dictionary only contains half-points (7.5 / 8.5, 3.5 / 4.5), the scraper today does an exact-line dict lookup, misses, and **skips the entire game's period for that book**. All 4 combos for that period drop out.

Confirmed coverage holes:
- **FanDuel F5:** alts in 1.0 increments (2.5, 3.5, 4.5, 5.5...) — every WZ F5 integer total skips
- **Novig FG:** observed in production — `nv_null = 2` on the only integer-FG game in the current slate
- **DK / FD / PX FG:** observed by user ("happens a lot")

The skip is graceful (R blender handles missing books), but it loses signal: integer-line games run with 3-4 sources in the blend instead of 5, making edge calculations noisier and less robust on exactly the games where push concentration is most relevant.

## 2. Goal

When a comparison book lacks the exact integer line WZ is using, derive that book's implied joint fair probability for each of the 4 SGP combos at the integer line, using the book's two adjacent half-point alt SGP prices. Write the derived probabilities to `mlb_sgp_odds` with a tagged `source` so downstream consumers and post-hoc analysis can distinguish derived vs directly-priced rows.

Out of scope:
- Changes to the EV / Kelly-sizing math in the R blender (the derived `fair_prob` slots into the existing pipeline as if the book had quoted the line directly)
- Any change to half-point line handling (existing exact-match path is unchanged)
- Spread-side integer-line interpolation (not currently a coverage hole — FG run-line is universally ±1.5; F5 spreads can be integer but the user confirmed totals are the dominant gap)

## 3. The math

### 3.1 Setup

Each comparison book offers SGP pricing for 4 combos at each alt total it supports:

- Home Spread + Over
- Home Spread + Under
- Away Spread + Over
- Away Spread + Under

For an integer line `X` that the book doesn't quote, we use its alts at `X − 0.5` and `X + 0.5`.

After per-alt devigging (multiplicative method, dividing each combo's implied prob by the sum of the 4 combos at that alt), each devigged value represents a joint probability `P(spread ∧ over)` or `P(spread ∧ under)` from that book's perspective.

### 3.2 Joint push mass

For the Sox-cover side, the joint probability the bet pushes (Sox cover ∧ T = X) is recoverable two ways from the same data:

```
Δ_cover = devig_HomeOver_lo  − devig_HomeOver_hi   (from the Over derivation)
       = devig_HomeUnder_hi − devig_HomeUnder_lo   (from the Under derivation)
```

Same identity for the Don't-cover side. The total marginal push mass:

```
Δ_total = Δ_cover + Δ_uncover  =  P(T = X)
```

This is a **pure-joint** computation — we sum joint push masses across spread sides to recover the marginal `P(T = X)` without ever touching the book's singles totals market.

### 3.3 Fair probability formula

For each of the 4 combos at the integer line:

```
fair_prob_X = p_hi_combo / (1 − Δ_total)
```

Where `p_hi_combo` is the devigged "harder" alt probability for that combo:
- For Over combos: devigged at `X + 0.5` (T ≥ X+1)
- For Under combos: devigged at `X − 0.5` (T ≤ X−1)

### 3.4 Interpretation

The derived `fair_prob_X` is the conditional joint probability "this combo wins outright, given the bet doesn't push." The conditioning event is `T ≠ X`, which is a marginal property of total — books refund pushed legs regardless of which spread side wins, so this matches real-world settlement behavior.

Under perfect devig the 4 derived `fair_prob_X` values sum to exactly 1.0:
```
Σ fair_prob_X = (Σ p_hi_combo) / (1 − Δ_total)
              = (1 − Δ_total) / (1 − Δ_total)
              = 1
```
This is the cornerstone bounds check (§5.4).

### 3.5 Worked example

Yankees @ Red Sox, FG. WZ has Sox -1.5 + Over 8 SGP. DK has alts at Over 7.5 and Over 8.5.

Devigged DK SGP probabilities at the two alts:

| Combo | Over 7.5 (`p_lo`) | Over 8.5 (`p_hi`) | Joint push mass |
|---|---|---|---|
| Sox -1.5 + Over | 0.320 | 0.270 | 0.050 (Sox cover side) |
| Sox -1.5 + Under | 0.160 | 0.213 | 0.053 (sanity-check duplicate) |
| Yankees +1.5 + Over | 0.347 | 0.302 | 0.045 (DC side) |
| Yankees +1.5 + Under | 0.173 | 0.218 | 0.045 (sanity-check duplicate) |

`Δ_total = 0.050 + 0.045 = 0.095` (about a 9.5% chance the game lands on exactly 8 runs, joint with any spread outcome)

Derived fair probs at integer line 8:
```
Sox -1.5 + Over 8     = 0.270 / (1 − 0.095) = 0.2983
Sox -1.5 + Under 8    = 0.160 / (1 − 0.095) = 0.1768
Yankees +1.5 + Over 8 = 0.302 / (1 − 0.095) = 0.3337
Yankees +1.5 + Under 8 = 0.173 / (1 − 0.095) = 0.1912
                                       Sum  = 1.000  ✓
```

If WZ offers Sox -1.5 + Over 8 at +260 (decimal 3.60):
```
Edge = fair_prob × WZ_decimal − 1
     = 0.2983 × 3.60 − 1
     = +7.4%
```

Vs. naive alternatives:
- Skip game (today) → no DK in blend → noisier edge
- Use Over 7.5 prob (0.320) → false positive +15.2%
- Use Over 8.5 prob (0.270) → false negative −2.8%
- Naive midpoint (0.295) → close but biases toward overstating edge

## 4. Architecture

### 4.1 New shared module

`mlb_sgp/integer_line_derivation.py`:

```python
def is_integer_line(value: float) -> bool: ...

def derive_fair_probs(
    devigged_lo: dict[str, float],   # 4 combos at X-0.5
    devigged_hi: dict[str, float],   # 4 combos at X+0.5
) -> dict[str, float] | None:
    """
    Returns dict of 4 derived fair_probs keyed by combo name,
    or None if any bounds check fails.
    """

def validate_per_alt_vig(implied_probs_sum: float) -> bool: ...
def validate_push_mass_consistency(...) -> bool: ...
def validate_delta_total(delta_total: float) -> bool: ...
def validate_sum_to_one(fair_probs: list[float]) -> bool: ...
def validate_per_combo_bounds(fair_prob: float) -> bool: ...
```

The shared module owns the math and the bounds checks. Each scraper calls `derive_fair_probs` after issuing its own pricing calls and devigging.

### 4.2 Per-scraper changes

Four scrapers each get a fallback block (~30 LOC) inserted in the existing combo-matching loop, before the current "skip if exact line missing" branch:

- `mlb_sgp/scraper_draftkings_sgp.py`
- `mlb_sgp/scraper_fanduel_sgp.py`
- `mlb_sgp/scraper_novig_sgp.py`
- `mlb_sgp/scraper_prophetx_sgp.py`

Each scraper handles the scraper-specific work (calling its own pricing API, mapping selection IDs at the two alt levels) and delegates the math to the shared module.

### 4.3 R blender change

One-line WHERE-clause expansion in `Answer Keys/mlb_correlated_parlay.R` to pick up `*_interpolated` source rows alongside `*_direct`:

```r
# Before
WHERE source IN ('draftkings_direct', 'fanduel_direct', 'novig_direct', 'prophetx_direct')
# After
WHERE source LIKE 'draftkings_%' OR source LIKE 'fanduel_%' OR source LIKE 'novig_%' OR source LIKE 'prophetx_%'
```

No changes to devig logic, blending, or EV math.

## 5. Algorithm (per scraper, per game, per period)

1. **Existing exact-line lookup** runs first. If hit, no change — proceed to `calculateBets` / `implyBets` / equivalent as today.

2. **On miss**, evaluate fallback eligibility:
   - Is `total_line` an integer (`is_integer_line(total_line)`)?
   - Are both `(total_line − 0.5)` and `(total_line + 0.5)` neighbors present in the alt totals dict?
   - Are both home and away spread selection IDs present at this WZ-matched spread line?
   - If any condition fails → fall through to existing skip behavior.

3. **Fire 8 SGP pricing calls** (4 combos × 2 alts) using the scraper's existing pricing function. Reuse the same parallel-pool / retry logic the scraper uses today.

4. **Per-alt devig** of each set of 4 combos:
   ```
   for each alt level:
       implied_probs = [1/decimal for each combo]
       vig_sum = sum(implied_probs)
       devigged = {combo: p/vig_sum for combo, p in implied_probs}
   ```

5. **Bounds check (1):** per-alt vig sum ∈ `[1.05, 1.30]`. If violated for either alt, skip + log WARN.

6. **Bounds check (2):** push-mass cross-consistency:
   - `|Δ_cover_from_over − Δ_cover_from_under| / max(...)` < 0.10
   - `|Δ_uncover_from_over − Δ_uncover_from_under| / max(...)` < 0.10

   If violated, skip + log WARN.

7. **Compute** `Δ_total = Δ_cover + Δ_uncover` using the average of the two consistent estimates per side.

8. **Bounds check (3):** `Δ_total ∈ [0.03, 0.18]`. If violated, skip + log WARN.

9. **Apply formula** to derive 4 `fair_prob_X` values.

10. **Bounds check (4):** sum of 4 derived fair_probs ∈ `[0.97, 1.03]`. If violated, skip + log WARN.

11. **Bounds check (5):** each `fair_prob_X` ∈ `(0, 1)`. If any violated, skip + log WARN (this should never trip but is cheap defensive).

12. **Write 4 rows to `mlb_sgp_odds`:**
    - `sgp_decimal = 1 / fair_prob_X`
    - `sgp_american` = corresponding American odds
    - `source = '<book>_interpolated'`
    - all other columns identical to direct-priced rows

## 6. Bounds checks summary

| # | Check | Threshold | Failure mode caught |
|---|---|---|---|
| 1 | Per-alt vig sum | `[1.05, 1.30]` | Scrape corruption, partial loads, transient API issues |
| 2 | Push-mass cross-consistency | <10% relative diff per side | Asymmetric per-side vig errors |
| 3 | `Δ_total` range | `[0.03, 0.18]` | Degenerate inputs (stale alt, mispriced edge) |
| 4 | Sum-to-1 | `[0.97, 1.03]` | Catches anything (1)-(3) miss; mathematically guaranteed under perfect inputs |
| 5 | Per-combo `[0, 1]` | strict | Defensive; should never trip |

All thresholds are config constants (top of `integer_line_derivation.py`), tunable in a follow-up commit once production logs reveal real-world distributions.

On any failure: skip the game for that book (do **not** write any of the 4 rows), log a structured WARN with inputs and the violated check. The R blender treats the book as missing for that game — same graceful degradation as today's "WZ line not found" path.

## 7. Data flow

```
Answer Keys/mlb_correlated_parlay.R
    │ writes mlb_parlay_lines (WZ's lines per game)
    ▼
4 scrapers in parallel (DK, FD, PX, NV)
    │ existing path: exact line match → 1 SGP price → mlb_sgp_odds
    │
    │ NEW path: integer-line fallback
    │   - 8 SGP pricing calls (4 combos × 2 alts)
    │   - per-alt devig
    │   - 5 bounds checks
    │   - 4 derived fair_probs → mlb_sgp_odds with source='*_interpolated'
    ▼
mlb_sgp_odds table
    │ unified consumption — R blender treats _direct and _interpolated rows identically
    ▼
Answer Keys/mlb_correlated_parlay.R (continued)
    │ per-game devig (no-op for _interpolated rows since they already sum to 1.0)
    │ blend with model + other books → blended_prob
    │ compare to WZ price → edge_pct → Kelly sizing → mlb_parlay_opportunities
```

## 8. Testing

### 8.1 Unit tests

`mlb_sgp/tests/test_integer_line_derivation.py`:

- **Math correctness:**
  - Worked example from §3.5 (`p_lo=0.32, p_hi=0.27` etc → `fair_prob = 0.2983`)
  - Boundary: `Δ_total = 0` (no push) → `fair_prob = p_hi`
  - Boundary: `Δ_total = 0.10` (high push concentration) → still produces valid output
  - Sum-to-1 invariant holds across multiple synthetic input sets

- **Bounds-check triggers:**
  - Per-alt vig 1.40 (corrupt) → returns `None`, logs WARN
  - Push-mass disagreement 30% relative → returns `None`
  - `Δ_total = 0.25` (out of range) → returns `None`
  - `Δ_total = 0.005` (out of range) → returns `None`
  - Sum-to-1 drift to 0.92 (e.g., from synthetic asymmetric vig) → returns `None`
  - Negative `fair_prob` from synthetic input → returns `None`

- **Helper functions:**
  - `is_integer_line(8.0)` → True; `is_integer_line(8.5)` → False
  - Tolerance handling for floating-point inputs (`8.000001` should treat as 8)

### 8.2 Integration test

A scripted dry-run on a recent slate where one or more books skipped due to integer-line miss. Verify:
- The new path triggers
- All 5 bounds checks pass on real data
- The derived `fair_prob` lands in plausible range
- The R blender picks up the new rows and computes a non-null `dk_fair_prob` / etc.

Captured as a one-off script (`mlb_sgp/tests/test_integration_integer_line.py`), not added to the cron — used for pre-merge validation.

### 8.3 Production monitoring (post-deploy)

- Add a daily summary line to the scraper logs: count of `_interpolated` rows written + count of bounds-check skips by check type
- After 1-2 weeks of slate data, review:
  - Distribution of measured per-alt vigs per book — tighten thresholds if real data is narrower than `[1.05, 1.30]`
  - Distribution of `Δ_total` — confirm `[0.03, 0.18]` is appropriate
  - Realized P&L on `_interpolated` rows vs `_direct` rows on the same combos — if interpolated edge significantly underperforms direct, investigate (likely candidate: void-leg refinement using spread-singles correction)

## 9. Documentation

### 9.1 mlb_sgp/README.md updates

- Replace existing "Exact line matching only" section with new fallback behavior description
- Add a new section "Integer-Line Derivation" with:
  - The formula `fair_prob_X = p_hi / (1 − Δ_total)`
  - Worked example
  - Citation to Breeden-Litzenberger (1978) as the underlying principle, noting the sports-betting application is a derivation rather than a published standard
- Update FD-specific note about F5 integer totals (no longer skipped — falls through to interpolation when alts present)
- Add bounds-check thresholds table

### 9.2 No CLAUDE.md changes required

The change is internal to the SGP pipeline; no global conventions or workflow changes for the project as a whole.

## 10. Version control

### 10.1 Worktree lifecycle

- **Created:** `/Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation` on branch `feature/sgp-integer-line-derivation`
- **Pre-merge:** full `git diff main..HEAD` review per CLAUDE.md executive engineer checklist
- **Merge to main:** explicit user approval required
- **Post-merge cleanup:** `git worktree remove .worktrees/sgp-integer-derivation` + `git branch -d feature/sgp-integer-line-derivation`

### 10.2 Commit structure

Suggested ordering (each commit independently testable):

1. `mlb_sgp: add integer_line_derivation shared helper + tests`
   - New file: `mlb_sgp/integer_line_derivation.py`
   - New file: `mlb_sgp/tests/test_integer_line_derivation.py`
   - All 5 bounds checks with thresholds as module-level constants

2. `mlb_sgp: DK scraper integer-line fallback`
   - Modifies `scraper_draftkings_sgp.py`
   - Reuses existing `calculateBets` pricing + parallel pool

3. `mlb_sgp: FD scraper integer-line fallback`
   - Modifies `scraper_fanduel_sgp.py`
   - Reuses existing `implyBets` pricing

4. `mlb_sgp: NV scraper integer-line fallback`
   - Modifies `scraper_novig_sgp.py`

5. `mlb_sgp: PX scraper integer-line fallback`
   - Modifies `scraper_prophetx_sgp.py`

6. `mlb_correlated_parlay: pick up *_interpolated source rows`
   - One-line WHERE-clause update in R orchestrator

7. `mlb_sgp: README integer-line derivation documentation`
   - Updates `mlb_sgp/README.md`

### 10.3 Pre-merge review checklist (per CLAUDE.md)

- **Data integrity:** No duplicate writes; bounds-check failures result in zero rows written, not partial
- **Resource safety:** No new DB connections without `on.exit` cleanup; existing scraper retry logic reused
- **Edge cases:** Off-season behavior (no integer-line games on slate → no fallback triggers, zero impact); first run after merge (no `_interpolated` rows yet → R WHERE clause LIKE pattern handles cleanly)
- **Dead code:** No unused config constants; per-scraper changes are minimal
- **Log/disk hygiene:** WARN logs go to existing `mlb_sgp/logs/<scraper>.log` (already rotated); no new files
- **Security:** No secrets in logs; the structured WARN line includes only numerical inputs and the violated threshold

## 11. Open questions / future work

- **Void-leg refinement:** the current formula treats pushes as "bet didn't happen" (refund-style). Under WZ's confirmed void-leg semantics, spread-winning pushes pay singles-spread payout (positive EV for the customer). The current formula slightly understates true fair_prob (and thus understates edge). If post-deploy P&L analysis shows interpolated rows underperforming directs by >1pp, refinement: multiply the conditional formula by `(1 + Δ_total × spread_singles_correction)` using each book's published spread singles odds. Treated as a v2 enhancement, not blocking v1.

- **Spread-side integer derivation:** F5 spreads can be integer (e.g., -1) where books only offer half-points. Same formula applies symmetrically (spread joins with totals; integer-mass derivable from adjacent half-point alts). Not currently a coverage hole per user; revisit if F5 spread skip rate becomes meaningful.

- **PX/NV vig calibration:** thresholds shipped as `[1.05, 1.30]` universal; tighten per book once production logs show their typical ranges.

- **Backtest before scaling sizing:** v1 ships with bounds-checks as the only guardrail. Backtest against historical pushes deliberately deferred per user direction. Realized P&L on `_interpolated` rows is the de facto validation — review at 2-week mark.

---

## Citations

- **Breeden, D. T. & Litzenberger, R. H. (1978).** "Prices of State-Contingent Claims Implicit in Option Prices." *Journal of Business*, 51(4), 621-651. The continuous-market analog: extracting risk-neutral probability density from option prices via second derivative w.r.t. strike. Sports application here is the discrete version (difference between adjacent half-point alts).

- **Half-point pricing logic.** Sportsbook Review and Bookmakers Review half-point calculators apply push-rate adjustments to single-leg fair odds. Same skeleton — `fair_at_integer = f(fair_at_half_point, push_rate)` — extended here to joint SGP probabilities.

- **Multiplicative devig.** Joseph Buchdahl, *Squares and Sharps, Suckers and Sharks*; also Unabated and OddsJam no-vig calculators. Already used elsewhere in this pipeline.

- **No published source** for the specific application to two-leg correlated parlays at integer total lines. Derived from first principles in the design discussion (2026-05-08); see `git log` for design conversation reference.
