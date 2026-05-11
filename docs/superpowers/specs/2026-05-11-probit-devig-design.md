# Probit Devigging — Design Spec

**Date:** 2026-05-11
**Branch:** `feature/probit-devig`
**Worktree:** `.worktrees/probit-devig`

---

## 1. Motivation

The codebase devigs sportsbook odds in many places — MLB pipeline, CBB pipeline, NFL backtests, MLB dashboard, CBB dashboard, Kalshi MLB RFQ bot. Today every one of them uses **multiplicative devigging**:

```
p_i = p_raw_i / sum_j(p_raw_j)
```

This is the simplest devig method but it under-weights extremes (longshots and big favorites) relative to the canonical **probit (z-shift) method** common in the quantitative betting literature. The goal of this work is to:

1. Replace multiplicative with additive z-shift probit devigging across every active code path.
2. Regenerate the MLB and CBB historical sample pools so historical and live consensus use the same method.
3. Fix latent inconsistencies in how devig functions are defined and reused (4 shadow copies that silently override the canonical helper).
4. Add defensive input validation and a 3-way → 2-way fallback to harden the helpers.

The math change is small (a few hundredths of a probability in the middle of the distribution, larger at the tails) but the *consistency* improvement is significant: today live odds, historical samples, dashboards, and downstream pricers can disagree subtly because multiple shadow definitions exist. After this change, one canonical implementation per language (R, Python) drives every devig call.

---

## 2. The math — probit (additive z-shift)

For each market with raw vigged implied probabilities `p_raw_1, ..., p_raw_n`:

1. Clip each `p_raw_i` to `[eps, 1 - eps]` (eps = 1e-9) to avoid `qnorm(0) = -Inf` / `qnorm(1) = +Inf`.
2. Compute z-scores: `z_i = qnorm(p_clipped_i)` — R's inverse normal CDF (Python: `scipy.stats.norm.ppf`).
3. Find a scalar `c` such that `sum(pnorm(z_i + c)) = 1` via 1-D root-find (R: `uniroot`, Python: `scipy.optimize.brentq`), bracket `[-5, 5]`.
4. Return `pnorm(z_i + c)` for each side.

**Properties:**
- Output sums to exactly 1.0 (the root-find tolerance, 1e-9).
- Same algorithm handles 2-way (moneylines, spreads, totals), 3-way (soccer / 1H markets with ties), and n-way (SGP joint distribution cells).
- Reduces to multiplicative when sides are symmetric (-110/-110 → 50/50 either method).
- At the tails (e.g., -1000/+800), probit gives the favorite about 0.3% more credit than multiplicative.

**Worked example** (-150 / +130):

| Method | Favorite | Dog |
|---|---|---|
| Multiplicative (current) | 57.98% | 42.02% |
| Probit additive z-shift (new) | 58.28% | 41.72% |

---

## 3. API surface

### R — `Answer Keys/Tools.R`

Public functions keep their exact signatures and return shapes. Internal body swapped for probit.

```r
# Internal: shared n-way probit devig
.probit_devig_n <- function(p_raw, eps = 1e-9) {
  if (any(is.na(p_raw))) return(rep(NA_real_, length(p_raw)))
  p_clipped <- pmin(pmax(p_raw, eps), 1 - eps)
  z <- qnorm(p_clipped)
  f <- function(c) sum(pnorm(z + c)) - 1
  c_star <- uniroot(f, interval = c(-5, 5), tol = 1e-9)$root
  pnorm(z + c_star)
}

# Public: 2-way (signature unchanged, return data.frame(p1, p2) unchanged)
devig_american <- function(odd1, odd2) {
  # Defensive: 0 or NA inputs → NA outputs (was: silent miscalculation)
  bad <- is.na(odd1) | is.na(odd2) | odd1 == 0 | odd2 == 0
  p1_raw <- ifelse(odd1 > 0, 100/(odd1+100), -odd1/(-odd1+100))
  p2_raw <- ifelse(odd2 > 0, 100/(odd2+100), -odd2/(-odd2+100))
  result <- mapply(function(p1, p2, is_bad) {
    if (is_bad) return(c(NA_real_, NA_real_))
    .probit_devig_n(c(p1, p2))
  }, p1_raw, p2_raw, bad)
  data.frame(p1 = result[1, ], p2 = result[2, ])
}

# Public: 3-way with 2-way fallback when odd_tie is NA
devig_american_3way <- function(odd_home, odd_away, odd_tie) {
  if (all(is.na(odd_tie))) {
    res <- devig_american(odd_home, odd_away)
    return(data.frame(p_home = res$p1, p_away = res$p2, p_tie = NA_real_))
  }
  # ... full 3-way probit per row (analogous structure) ...
}
```

### Python — `kalshi_mlb_rfq/fair_value.py`

```python
from scipy.stats import norm
from scipy.optimize import brentq

def _probit_devig_n(p_raw: list[float], eps: float = 1e-9) -> list[float]:
    p_clipped = [min(max(p, eps), 1 - eps) for p in p_raw]
    z = [norm.ppf(p) for p in p_clipped]
    f = lambda c: sum(norm.cdf(zi + c) for zi in z) - 1
    c_star = brentq(f, -5, 5, xtol=1e-9)
    return [norm.cdf(zi + c_star) for zi in z]

def devig_book(book_rows, combo, vig_fallback=0.0):
    if book_rows.empty: return None
    target = book_rows.loc[book_rows["combo"] == combo]
    if target.empty: return None
    target_decimal = float(target["sgp_decimal"].iloc[0])

    if len(book_rows) >= 4:
        raw_probs = (1.0 / book_rows["sgp_decimal"]).tolist()
        devigged = _probit_devig_n(raw_probs)
        target_idx = book_rows.index.get_loc(target.index[0])
        return devigged[target_idx]

    # Single-side fallback: heuristic, no devig math possible with 1 cell
    return (1.0 / target_decimal) / (1.0 + vig_fallback)
```

Adds `scipy>=1.10` to `kalshi_mlb_rfq/requirements.txt` (scipy already in use elsewhere in the project, e.g. `clv_compute.py`).

---

## 4. Bug fixes (in the same change)

### 4.1 Shadow definitions — kill all 4

Today four files define their own `devig_american*` functions that silently override or sit alongside the canonical Tools.R helpers:

| Location | Action |
|---|---|
| `NFL Answer Key/NFLAnswerKey.R:9` (local `devig_american`, file does NOT source Tools.R; one commit ever — likely legacy after `NFLAnswerKey2.0.R`) | Delete local helper; add `source("Tools.R")`. If file confirmed legacy, recommend deletion in a follow-up commit. |
| `NFL Answer Key/All_Quarters_Backtest.R:20,79` (sources Tools.R but redefines after the source — shadow wins) | Delete local `devig_american_3way` and `devig_american_pair`. Rewrite callers: `devigged$home_prob → devigged$p_home`, `devigged$away_prob → devigged$p_away`. Replace `devig_american_pair(...)` calls with `devig_american(...)`. |
| `NFL Answer Key/Spreads_Totals_Backtest.R:23` (same pattern) | Same treatment as above. |
| `Tools.R:5561` (nested `devig_american_pair` inside an outer function — scope-shadows the global) | Delete the nested helper; use the global `devig_american` directly. |

### 4.2 Input validation

`devig_american(0, +200)` today silently returns `(0, 1)` — American odds of `0` are invalid but the function treats the side as a 0% prob. After the change: NA on either leg, or `odd == 0` on either leg, returns `(NA, NA)`. Same for 3-way.

### 4.3 3-way fallback

`devig_american_3way(odd_home, odd_away, NA)` today returns `(NA, NA, NA)`. After the change: when `odd_tie` is fully NA, the function transparently routes to 2-way devig and returns `(p_home, p_away, NA)`.

---

## 5. Historical pool regeneration

### 5.1 MLB

- Re-run `Answer Keys/MLB Answer Key/Consensus Betting History.R`.
- This rebuilds `mlb_betting_pbp` in `pbp.duckdb` with probit-devigged consensus probabilities.
- **Before running:** snapshot the old table — `CREATE TABLE mlb_betting_pbp_backup_2026_05_11 AS SELECT * FROM mlb_betting_pbp` — for rollback.

### 5.2 CBB

The CBB pipeline today does NOT have a pre-built consensus history. `cbb_betting_pbp` contains only outcomes; the `0.5` hardcode block in `CBB.R:51-63` stamps placeholder probability columns at runtime.

After the change:
- **Extend `Answer Keys/CBB Answer Key/build_betting_pbp.R`** to compute per-game per-book probit devig + sharp-weighted consensus from `cbb_closing_odds`, mirroring MLB's pattern (sharp weights via `SHARP_BOOKS` in `Tools.R`).
- **Write the consensus columns** (`consensus_devig_home_odds`, `consensus_devig_away_odds`, `consensus_devig_over_odds`, `consensus_devig_under_odds`) into `cbb_betting_pbp` alongside outcomes — same self-contained shape as MLB's `mlb_betting_pbp`.
- **Drop the `0.5` hardcode in `CBB.R:51-63`** — CBB.R reads consensus columns directly from `cbb_betting_pbp`.
- **Before running:** snapshot — `CREATE TABLE cbb_betting_pbp_backup_2026_05_11 AS SELECT * FROM cbb_betting_pbp`.
- **After the schema change**, rerun `build_betting_pbp.R` to repopulate.

### 5.3 NFL

`Consensus Betting History.R` (NFL) and `Historical Betting History 1999-2019.R` both source `Tools.R`. After the helper swap they'll automatically use probit. **Optionally** re-run them to refresh `nfl_betting_pbp` and `nfl_pre_20_betting_history`. Not strictly required for live pricing — only matters if NFL backtests are run before re-generation. Decision: defer NFL regen unless backtests need it; flag in PR description.

---

## 6. Tests

### 6.1 R unit tests — `Answer Keys/tests/test_probit_devig.R` (new file)

| Case | Expected |
|---|---|
| `devig_american(-110, -110)` | `(0.5, 0.5)` exactly (symmetric — probit ≡ multiplicative) |
| `devig_american(-200, +180)` | known probit values, frozen via high-precision compute |
| `devig_american(-1000, +800)` | tail case — probit visibly diverges from multiplicative |
| `devig_american(0, +100)` | `(NA, NA)` (invalid input) |
| `devig_american(NA, -110)` | `(NA, NA)` (NA propagation) |
| `devig_american_3way(+150, +200, +400)` | sums to 1.0 exactly |
| `devig_american_3way(-110, -110, NA)` | falls back to 2-way: `(0.5, 0.5, NA)` |

### 6.2 Python unit tests — `kalshi_mlb_rfq/tests/test_probit_devig.py` (new file)

| Case | Expected |
|---|---|
| `_probit_devig_n([0.55, 0.50])` | sums to 1.0 |
| `_probit_devig_n([0.27, 0.27, 0.27, 0.27])` | all four ≈ 0.25 (symmetric 4-cell) |
| `devig_book(synthetic_4_row_df, "Home Spread + Over")` | round-trip matches `_probit_devig_n` output for that index |
| `devig_book(single_row_df, "Home Spread + Over", vig_fallback=0.10)` | matches `(1/decimal) / 1.10` (fallback heuristic unchanged) |

### 6.3 Backtest validation (gates the merge)

- **MLB:** snapshot `mlb_betting_pbp`. Run new `Consensus Betting History.R`. Run `MLB_Backtest.R` and `MLB_ROI_Backtest.R`. Record ROI, win rate, CLV. Compare against pre-change baseline. Expected: small shifts (≤1% on most metrics). If anything moves >5%, investigate before merging.
- **CBB:** same procedure with `cbb_betting_pbp` + `CBB_Backtest.R`.
- **NFL:** if NFL regen deferred, just run `All_Quarters_Backtest.R` and `Spreads_Totals_Backtest.R` on the live (post-helper-swap) code and confirm they still execute. They'll use probit live but old multiplicative-baked historical pool — that's the deferred regen state.

### 6.4 Live sanity check (post-merge, before declaring done)

- Start MLB dashboard (port 8083), hover over a game, verify "Books (devigged fair %)" shows plausible probit values.
- Start CBB dashboard, same check.
- Start Kalshi MLB RFQ bot in dry-run, verify the book-fair logs show reasonable values for a known game.

---

## 7. Version control plan

**Worktree:** `.worktrees/probit-devig`
**Branch:** `feature/probit-devig`

**Commit sequence:**

1. `feat(devig): add probit z-shift helper, refactor Tools.R`
   - New `.probit_devig_n` internal helper
   - Swap `devig_american` and `devig_american_3way` bodies
   - Add input validation + 3-way fallback
   - Add `Answer Keys/tests/test_probit_devig.R`

2. `refactor(devig): consolidate shadow definitions`
   - Delete local devig helpers in `NFLAnswerKey.R`, `All_Quarters_Backtest.R`, `Spreads_Totals_Backtest.R`, `Tools.R:5561`
   - Update callers to use canonical Tools.R helpers
   - Verify backtest scripts still run

3. `feat(cbb): probit-devigged consensus in build_betting_pbp.R`
   - Extend `build_betting_pbp.R` to compute per-book devig + sharp-weighted consensus
   - Write `consensus_devig_*` columns into `cbb_betting_pbp`
   - Drop `0.5` hardcode in `CBB.R:51-63`
   - Note in commit body: requires `Rscript build_betting_pbp.R` to regenerate after merge

4. `feat(kalshi-mlb-rfq): probit devig for SGP joint distribution`
   - Port `_probit_devig_n` to Python
   - Swap `devig_book` body to use probit on the 4-cell grid
   - Single-side fallback unchanged
   - Add `scipy>=1.10` to `requirements.txt`
   - Add `kalshi_mlb_rfq/tests/test_probit_devig.py`

5. `docs: update CLAUDE.md + READMEs for probit devig`
   - `Answer Keys/CLAUDE.md`: update Consensus Architecture section (CBB no longer 0.5 hardcoded)
   - `Answer Keys/MLB Dashboard/README.md`: note probit method
   - `Answer Keys/CBB Dashboard/README.md`: same
   - `kalshi_mlb_rfq/README.md`: note probit method in fair-value layer
   - This repo's CLAUDE.md doesn't mention devig method, so no top-level update needed

**Files modified** (8):
- `Answer Keys/Tools.R`
- `Answer Keys/NFL Answer Key/NFLAnswerKey.R`
- `Answer Keys/NFL Answer Key/All_Quarters_Backtest.R`
- `Answer Keys/NFL Answer Key/Spreads_Totals_Backtest.R`
- `Answer Keys/CBB Answer Key/CBB.R`
- `Answer Keys/CBB Answer Key/build_betting_pbp.R`
- `kalshi_mlb_rfq/fair_value.py`
- `kalshi_mlb_rfq/requirements.txt`

**Files created** (3):
- `Answer Keys/tests/test_probit_devig.R`
- `kalshi_mlb_rfq/tests/test_probit_devig.py`
- This spec doc

**Files for docs updates** (4):
- `Answer Keys/CLAUDE.md`
- `Answer Keys/MLB Dashboard/README.md`
- `Answer Keys/CBB Dashboard/README.md`
- `kalshi_mlb_rfq/README.md`

---

## 8. Pre-merge review checklist

Per project rules (`Answer Keys/CLAUDE.md` pre-merge review):

- [ ] All unit tests pass (R + Python)
- [ ] MLB backtest before/after comparison reviewed — no metric moves >5%
- [ ] CBB backtest before/after comparison reviewed — no metric moves >5%
- [ ] `mlb_betting_pbp_backup_2026_05_11` snapshot exists
- [ ] `cbb_betting_pbp_backup_2026_05_11` snapshot exists
- [ ] MLB dashboard hover text visually validated
- [ ] CBB dashboard hover text visually validated
- [ ] Kalshi MLB RFQ bot dry-run shows plausible book fairs
- [ ] All 4 shadow definitions confirmed gone (`grep -rn "devig_american_pair\|local devig" "Answer Keys"`)
- [ ] No regressions in `Answer Keys/tests/test_answer_key.R`
- [ ] Documentation updates merged in the same PR
- [ ] User has given explicit approval to merge

---

## 9. Worktree lifecycle

1. **Create:** `git worktree add .worktrees/probit-devig -b feature/probit-devig` (already done)
2. **Work:** all edits happen inside `.worktrees/probit-devig/`
3. **Test:** runs from within the worktree (R scripts source local `Tools.R`, Python imports local `kalshi_mlb_rfq` modules)
4. **Backtest/regen:** done from `main` after merge (because historical pool tables live in `pbp.duckdb` at the repo root, not in the worktree)
5. **Merge:** after user approval, `git checkout main && git merge --no-ff feature/probit-devig`
6. **Regenerate historical pools** on `main`: `Rscript "MLB Answer Key/Consensus Betting History.R"` and `Rscript "CBB Answer Key/build_betting_pbp.R"`
7. **Cleanup:** `git worktree remove .worktrees/probit-devig && git branch -d feature/probit-devig`

---

## 10. Out of scope

- `nfl_draft/scrapers/_base.py:devig_prob` field — used for portal display only, not pricing-critical. Defer unless flagged.
- `Answer Keys/MLB Answer Key/MLB_ROI_Backtest.R` extensive single-side devig paths inside the body — uses Tools.R helpers, will pick up probit automatically.
- Top-level `NFLWork/CLAUDE.md` — doesn't reference devig method, no update needed.

---

## 11. Risk register

| Risk | Mitigation |
|---|---|
| Probit `uniroot` fails to bracket in extreme markets (e.g., overround > 30%) | `[-5, 5]` bracket covers any realistic vig; fallback to multiplicative if `uniroot` errors (with a logged warning) |
| Regenerated `mlb_betting_pbp` changes sample-matching outputs noticeably, affecting live edges | Backtest gate before merge; rollback table available |
| CBB.R reads renamed columns that aren't yet in `cbb_betting_pbp` (out-of-order rollout) | Commit order: schema change first (build_betting_pbp.R), then CBB.R reads. Same PR enforces ordering. |
| Kalshi RFQ bot reads `scipy` not installed on production venv | Add `scipy>=1.10` to `requirements.txt`; user must `pip install -r requirements.txt` post-merge |
| Hidden additional shadow devig copies I missed | Final grep pre-merge: `grep -rn "devig_american\|devig_book" --include="*.R" --include="*.py"` and review every hit |

---

## 12. Memory note (post-merge)

Save a project memory documenting the methodology change, since this affects every downstream backtest and is a non-obvious choice future-me may need to recall. Suggested location: `memory/probit_devig_methodology.md` indexed in `MEMORY.md`.
