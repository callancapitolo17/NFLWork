# MLB Combined Parlay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a feature to the MLB Dashboard's Parlay tab that lets the user select two existing recommended parlay rows, combine them into one 4-leg Wagerzon ticket priced via `ConfirmWagerHelper`, and dynamically reduce the source rows' Kelly recommendations to the conditional residuals after the combo is placed.

**Architecture:** Two halves. (1) Python in the Flask dashboard server handles combined-parlay pricing endpoints (joint fair odds, exact WZ payout, recommended Kelly stake) and persistence. (2) R in `create_parlays_table()` reads placed combos from `mlb_dashboard.duckdb` and applies conditional Kelly to source rows' Kelly column when the table HTML is rendered. The combined parlay banner UI is checkbox-driven on the Parlay tab.

**Tech Stack:** R (Shiny + reactable for the dashboard), Python (Flask backend, requests for WZ API), DuckDB (state), JavaScript (banner + checkbox interactions). R uses `testthat`; Python uses `pytest`.

---

## Performance Targets

Two perceptual latencies matter:

1. **Banner pricing on 2nd checkbox click** — currently bounded by the live `ConfirmWagerHelper` call (~1-3s). Mitigation: **server-side TTL cache** of priced combos keyed on `(hash_a, hash_b)`, TTL 60s. Cache hits return instantly (re-checking the same pair, undo/redo, etc.).
2. **Post-place "see updated Kelly" delay** — currently bounded by full-page reload + R re-render (potentially several seconds on a busy slate). Mitigation: **client-side partial DOM update** for the source rows' Kelly cells, with full reload as a fallback path.

Both are wired into the plan tasks below. WZ session reuse (`requests.Session()`) is already the existing pattern in `parlay_pricer.py` — preserve it.

Out of scope for v1 (defer if v1 still feels slow):
- Background pre-warming of top-N pairs (would need a background scheduler — overkill until usage justifies it)
- Client-side conditional Kelly (would need to port the optim() solver to JS — adds complexity for marginal gain since residuals are computed once on place)

---

## File Structure

### New files
- `Answer Keys/conditional_kelly.R` — `conditional_kelly_residuals()` for R-side residual math
- `Answer Keys/MLB Dashboard/combined_parlay.py` — Python helpers for joint pricing + Kelly stake (no R dependency)
- `Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py` — One-shot migration script for `placed_parlays` schema
- `Answer Keys/tests/test_conditional_kelly.R` — testthat tests for residual math
- `Answer Keys/MLB Dashboard/tests/test_combined_parlay.py` — pytest tests for combined-parlay helpers + endpoints

### Modified files
- `Answer Keys/Tools.R` — add `compute_combined_parlay_pricing()` (joint fair odds across two parlay rows)
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add Sel column + banner div + residual-aware Kelly column rendering in `create_parlays_table()`
- `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` — new `/api/price-combined-parlay` and `/api/place-combined-parlay` endpoints
- `wagerzon_odds/parlay_pricer.py` — refactor `get_parlay_price()` to allow legs from multiple games (pass per-leg `idgm`); add `get_combined_parlay_price()` wrapper
- `Answer Keys/MLB Dashboard/lib/dashboard.js` (or wherever the existing parlay JS lives — verify in Task 8) — checkbox handler, banner rendering, place-combined POST
- `Answer Keys/MLB Dashboard/README.md` — document the new feature

---

## Task Sequence and Dependencies

```
Task 1 (R: conditional Kelly math)
   └→ Task 3 (R: apply residuals in table) ─────────┐
                                                     │
Task 2 (R: joint pricing helper) ─────┐              │
                                       │              │
Task 4 (Py: cross-game WZ pricer) ─┐  │              │
                                    │  │              │
Task 5 (Py: schema migration) ──┐  │  │              │
                                 │  │  │              │
                                 ↓  ↓  ↓              ↓
                       Task 6 (Py: /api/price-combined)
                                 │
                                 ↓
                       Task 7 (Py: /api/place-combined)
                                 │
                                 ↓
                       Task 8 (R/JS: checkbox column)
                                 │
                                 ↓
                       Task 9 (HTML/JS: banner)
                                 │
                                 ↓
                       Task 10 (JS: Place Combined wiring)
                                 │
                                 ↓
                       Task 11 (docs: README)
```

---

## Task 1: R conditional Kelly residual math

**Files:**
- Create: `Answer Keys/conditional_kelly.R`
- Test: `Answer Keys/tests/test_conditional_kelly.R`

- [ ] **Step 1: Write the failing test**

Create `Answer Keys/tests/test_conditional_kelly.R`:

```r
library(testthat)
source("../conditional_kelly.R")

test_that("conditional Kelly returns non-negative residuals", {
  res <- conditional_kelly_residuals(
    p_a = 0.50, d_a = 2.10,
    p_b = 0.50, d_b = 2.10,
    s_combo = 50, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  expect_true(res$s_a >= 0)
  expect_true(res$s_b >= 0)
})

test_that("conditional Kelly residual <= unconstrained Kelly", {
  # Single A's full unconstrained Kelly: f = (p*d - 1) / (d - 1)
  full_kelly_a <- ((0.55 * 2.10) - 1) / 1.10 * 1000  # ~$45
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10,
    p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  expect_true(res$s_a <= full_kelly_a + 0.01)  # tolerance for numerical
  expect_true(res$s_b <= full_kelly_a + 0.01)
})

test_that("conditional Kelly with s_combo = 0 matches independent Kelly", {
  # Sanity: if no combo, residual should match standard Kelly per leg
  full_kelly_a <- ((0.55 * 2.10) - 1) / 1.10 * 1000
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10,
    p_b = 0.55, d_b = 2.10,
    s_combo = 0, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  # Allow ~$1 tolerance — numerical solver isn't exact
  expect_equal(res$s_a, full_kelly_a, tolerance = 1.0)
  expect_equal(res$s_b, full_kelly_a, tolerance = 1.0)
})

test_that("conditional Kelly with kelly_mult halves the residual", {
  full <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  half <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 0.5
  )
  expect_equal(half$s_a, full$s_a * 0.5, tolerance = 0.5)
})
```

- [ ] **Step 2: Run tests to verify they fail (function doesn't exist yet)**

Run from `Answer Keys/tests/`:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_conditional_kelly.R")'
```
Expected: FAIL with `could not find function "conditional_kelly_residuals"`.

- [ ] **Step 3: Implement `conditional_kelly_residuals()`**

Create `Answer Keys/conditional_kelly.R`:

```r
# conditional_kelly.R — Kelly residuals for the singles when a combined parlay
# is already placed. Used by the MLB Dashboard's Parlay tab to update the
# source rows' Kelly column after a combined ticket is locked in.
#
# Approach: for two independent legs A and B and a fixed combo stake s_combo,
# numerically maximize E[log(1 + total return)] over (s_a, s_b) using the
# 4-outcome joint distribution. Reuses the same Kelly-portfolio framing as
# parlay_multivariate_kelly() but with one position held fixed.

conditional_kelly_residuals <- function(p_a, d_a,
                                         p_b, d_b,
                                         s_combo, d_combo,
                                         bankroll, kelly_mult = 1.0) {
  # Per-bet returns: profit per $1 wagered if win; -$1 if lose
  b_a <- d_a - 1
  b_b <- d_b - 1
  b_c <- d_combo - 1

  # Combo as bankroll fraction (held fixed during optimization)
  f_c <- s_combo / bankroll

  # Joint probabilities under cross-game independence
  p_ab      <- p_a * p_b
  p_a_only  <- p_a * (1 - p_b)
  p_b_only  <- (1 - p_a) * p_b
  p_neither <- (1 - p_a) * (1 - p_b)

  # Negative expected log return — minimize this = maximize log growth
  neg_e_log <- function(f) {
    f_a <- f[1]; f_b <- f[2]
    r_both    <- f_a * b_a + f_b * b_b + f_c * b_c
    r_a_only  <- f_a * b_a - f_b           - f_c
    r_b_only  <- -f_a       + f_b * b_b    - f_c
    r_neither <- -f_a       - f_b           - f_c
    -(p_ab      * log1p(r_both)    +
      p_a_only  * log1p(r_a_only)  +
      p_b_only  * log1p(r_b_only)  +
      p_neither * log1p(r_neither))
  }

  result <- optim(
    par    = c(0.01, 0.01),
    fn     = neg_e_log,
    method = "L-BFGS-B",
    lower  = c(0, 0),
    upper  = c(0.5, 0.5)  # safe upper bound — Kelly never recommends > 50% bankroll
  )

  s_a <- result$par[1] * kelly_mult * bankroll
  s_b <- result$par[2] * kelly_mult * bankroll

  list(s_a = round(max(s_a, 0), 2),
       s_b = round(max(s_b, 0), 2))
}
```

- [ ] **Step 4: Run tests — expect all 4 to pass**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_conditional_kelly.R")'
```
Expected: PASS — 4 / 4.

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay
git add "Answer Keys/conditional_kelly.R" "Answer Keys/tests/test_conditional_kelly.R"
git commit -m "feat(answer-keys): conditional Kelly residual solver for combined parlays

Numerical L-BFGS-B optimizer over the 4-outcome joint distribution.
Given a fixed combo stake, returns the optimal residual stakes on each
single. Reuses the same Kelly-portfolio framing as parlay_multivariate_kelly
but with one position held fixed (per the design spec).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: R combined parlay pricing helper

**Files:**
- Modify: `Answer Keys/Tools.R` (add new function near `compute_parlay_fair_odds` at line 4924)
- Test: append to `Answer Keys/tests/test_conditional_kelly.R` (or new file)

- [ ] **Step 1: Write the failing test**

Append to `Answer Keys/tests/test_conditional_kelly.R`:

```r
source("../Tools.R")

test_that("compute_combined_parlay_pricing multiplies fair decimal odds", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10
  )
  expect_equal(res$joint_fair_dec, 4.32 * 4.85, tolerance = 0.001)
  expect_equal(res$joint_fair_prob, 1 / (4.32 * 4.85), tolerance = 0.0001)
})

test_that("compute_combined_parlay_pricing computes joint edge correctly", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10
  )
  expected_edge <- (1 / (4.32 * 4.85)) * 22.10 - 1
  expect_equal(res$joint_edge, expected_edge, tolerance = 0.0001)
})

test_that("compute_combined_parlay_pricing returns kelly stake using parlay_kelly_mult", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10,
    bankroll = 1000, kelly_mult = 0.5
  )
  expect_true(res$kelly_stake >= 0)
  expect_true(res$kelly_stake <= 500)  # never more than half bankroll at Kelly_mult = 0.5
})
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_conditional_kelly.R")'
```
Expected: 3 NEW failures (`could not find function "compute_combined_parlay_pricing"`).

- [ ] **Step 3: Implement `compute_combined_parlay_pricing()`**

Add to `Answer Keys/Tools.R` immediately after `compute_parlay_fair_odds()` (around line 4974):

```r
# compute_combined_parlay_pricing — joint pricing for a 4-leg combined parlay
# made from two existing same-game spread+total parlay rows. Cross-game legs
# are independent, so joint fair odds = product of fair decimal odds.
#
# Args:
#   fair_dec_a, fair_dec_b : fair decimal odds of each source parlay row
#   wz_dec                 : Wagerzon's exact 4-leg payout (from
#                            ConfirmWagerHelper) as decimal odds
#   bankroll, kelly_mult   : optional, for kelly_stake computation
#
# Returns: list(joint_fair_dec, joint_fair_prob, joint_edge, kelly_stake)
compute_combined_parlay_pricing <- function(fair_dec_a, fair_dec_b, wz_dec,
                                             bankroll = NULL, kelly_mult = NULL) {
  joint_fair_dec  <- fair_dec_a * fair_dec_b
  joint_fair_prob <- 1 / joint_fair_dec
  joint_edge      <- joint_fair_prob * wz_dec - 1

  kelly_stake_dollars <- if (!is.null(bankroll) && !is.null(kelly_mult)) {
    edge_fraction <- joint_edge / (wz_dec - 1)
    round(max(edge_fraction * kelly_mult * bankroll, 0), 2)
  } else NA_real_

  list(
    joint_fair_dec  = round(joint_fair_dec, 4),
    joint_fair_prob = round(joint_fair_prob, 6),
    joint_edge      = round(joint_edge, 4),
    kelly_stake     = kelly_stake_dollars
  )
}
```

- [ ] **Step 4: Run tests — expect all to pass**

Same command as Step 2. Expected: PASS — 7 / 7 (4 from Task 1 + 3 new).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_conditional_kelly.R"
git commit -m "feat(answer-keys): combined parlay pricing helper

Pure-arithmetic joint pricing for two cross-game parlays. Independence
across games means joint fair = product of decimals; the only book-side
input is Wagerzon's exact 4-leg payout from ConfirmWagerHelper.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: R apply residuals in `create_parlays_table()`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (around lines 249–438)

- [ ] **Step 1: Write the test (manual smoke test for now — table rendering is hard to unit test)**

We'll verify by inspecting `report.html` after running the dashboard pipeline with a mock combo placed. No automated test in this step — Task 7's integration test covers the API; this step is purely the rendering plumbing.

- [ ] **Step 2: Source the new helper at the top of `mlb_dashboard.R`**

Add near the existing `source()` calls at the top of `Answer Keys/MLB Dashboard/mlb_dashboard.R`:

```r
source(file.path("..", "conditional_kelly.R"))
```

- [ ] **Step 3: Add a helper that computes residuals for any source rows in placed combos**

Add this new function above `create_parlays_table()` in `mlb_dashboard.R`:

```r
# Returns parlay_opps with kelly_bet / size_display overridden to the conditional
# residual when the row appears in any combo in placed_parlays.
apply_combo_residuals <- function(parlay_opps, placed_parlays, parlay_bankroll, parlay_kelly_mult) {
  combos <- placed_parlays %>% filter(is_combo == TRUE)
  if (nrow(combos) == 0) return(parlay_opps)

  for (i in seq_len(nrow(combos))) {
    leg_ids   <- jsonlite::fromJSON(combos$combo_leg_ids[i])
    combo_dec <- combos$wz_odds[i] / 100 + 1   # convert American to decimal if needed
    if (combos$wz_odds[i] < 0) {
      combo_dec <- 100 / abs(combos$wz_odds[i]) + 1
    }
    s_combo <- combos$actual_size[i]

    row_a <- parlay_opps %>% filter(parlay_hash == leg_ids[1])
    row_b <- parlay_opps %>% filter(parlay_hash == leg_ids[2])
    if (nrow(row_a) == 0 || nrow(row_b) == 0) next

    res <- conditional_kelly_residuals(
      p_a     = 1 / row_a$fair_dec,
      d_a     = row_a$wz_dec,
      p_b     = 1 / row_b$fair_dec,
      d_b     = row_b$wz_dec,
      s_combo = s_combo,
      d_combo = combo_dec,
      bankroll   = parlay_bankroll,
      kelly_mult = parlay_kelly_mult
    )

    parlay_opps$kelly_bet[parlay_opps$parlay_hash == leg_ids[1]] <- res$s_a
    parlay_opps$kelly_bet[parlay_opps$parlay_hash == leg_ids[2]] <- res$s_b
    parlay_opps$combo_residual_note[parlay_opps$parlay_hash %in% leg_ids] <-
      sprintf("(residual after combo $%.0f)", s_combo)
  }
  parlay_opps
}
```

- [ ] **Step 4: Wire it into `create_parlays_table()`**

In `create_parlays_table()` (around line 249 of `mlb_dashboard.R`), at the top of the function, add a `combo_residual_note` column and call the helper. Modify the start of the function from:

```r
create_parlays_table <- function(parlay_opps, placed_parlays) {
  placed_hashes <- if (nrow(placed_parlays) > 0) placed_parlays$parlay_hash else character()
```

to:

```r
create_parlays_table <- function(parlay_opps, placed_parlays, parlay_bankroll = 100, parlay_kelly_mult = 0.25) {
  parlay_opps$combo_residual_note <- NA_character_
  parlay_opps <- apply_combo_residuals(parlay_opps, placed_parlays, parlay_bankroll, parlay_kelly_mult)
  placed_hashes <- if (nrow(placed_parlays) > 0) placed_parlays$parlay_hash %>% setdiff(placed_parlays$parlay_hash[placed_parlays$is_combo == TRUE]) else character()
```

- [ ] **Step 5: Update the caller to pass parlay_bankroll/kelly_mult**

Find the call site of `create_parlays_table()` in `mlb_dashboard.R` (use `grep -n "create_parlays_table" "Answer Keys/MLB Dashboard/mlb_dashboard.R"`) and update it to pass the existing `parlay_bankroll` and `parlay_kelly_mult` settings (already loaded earlier in the file — same source as the existing single-row Kelly).

- [ ] **Step 6: Update Size column rendering to show the residual annotation when present**

Find the `size_display` mutate (around line 280):

```r
size_display = sprintf("$%.0f", coalesce(as.numeric(exact_wager), kelly_bet)),
```

Replace with:

```r
size_display = ifelse(
  !is.na(combo_residual_note),
  sprintf("$%.0f<br><span class='combo-note'>%s</span>", kelly_bet, combo_residual_note),
  sprintf("$%.0f", coalesce(as.numeric(exact_wager), kelly_bet))
),
```

And mark the `size_display` column as `html = TRUE` in the `colDef` for it (find the `colDef(name = "Size", ...)` block and add `html = TRUE`).

- [ ] **Step 7: Smoke-test by running the dashboard pipeline**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay/Answer Keys/MLB Dashboard"
# Copy the live mlb_dashboard.duckdb so we test against real placed_parlays
cp /Users/callancapitolo/NFLWork/Answer\ Keys/MLB\ Dashboard/mlb_dashboard.duckdb ./
Rscript mlb_dashboard.R   # generates report.html
open report.html          # visual confirm: existing rows unchanged when no combos placed
```
Expected: report.html builds without error; existing rows render unchanged (no combos exist yet, so `apply_combo_residuals` is a no-op).

- [ ] **Step 8: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): apply conditional Kelly residuals to source rows

When a combined parlay row exists in placed_parlays, look up the two
source parlay_hashes and override their kelly_bet with the residual
computed by conditional_kelly_residuals(). Residual annotation shows
the combo stake that triggered the reduction.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Cross-game ConfirmWagerHelper in Wagerzon pricer

**Files:**
- Modify: `wagerzon_odds/parlay_pricer.py` (around lines 37–128)
- Test: `Answer Keys/MLB Dashboard/tests/test_combined_parlay.py` (new file)

- [ ] **Step 1: Write the failing test**

Create `Answer Keys/MLB Dashboard/tests/test_combined_parlay.py`:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "wagerzon_odds"))

import pytest
from unittest.mock import MagicMock, patch
import parlay_pricer


def test_get_combined_parlay_price_builds_sel_with_per_leg_idgm():
    """Cross-game parlay must encode each leg's own idgm in the sel string."""
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},   # game A home spread
        {"idgm": 100001, "play": 2, "points": "9.5", "odds": -105},   # game A over
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},  # game B away spread
        {"idgm": 100002, "play": 3, "points": "7.5", "odds": -110},   # game B under
    ]
    expected_sel_parts = [
        "1_100001_-1.5_110",
        "2_100001_9.5_-105",
        "0_100002_+1.5_-120",
        "3_100002_7.5_-110",
    ]

    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 21000, "Risk": 1000}],
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=1000)

    # Inspect the actual call
    call_args = fake_session.post.call_args
    sel_value = call_args.kwargs["data"]["sel"]
    for part in expected_sel_parts:
        assert part in sel_value, f"Missing leg in sel: {part}"
    assert result is not None
    assert result["win"] == 21000
    assert result["decimal"] == pytest.approx(22.0, rel=0.01)


def test_get_combined_parlay_price_returns_none_on_error():
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},
    ]
    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 0, "Risk": 0}],
            "ErrorMsgKey": "MAXPARLAYRISKEXCEED",
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=10000)
    assert result is None
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay
mkdir -p "Answer Keys/MLB Dashboard/tests"
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py" -v
```
Expected: FAIL with `AttributeError: module 'parlay_pricer' has no attribute 'get_combined_parlay_price'`.

- [ ] **Step 3: Implement `get_combined_parlay_price()` in `parlay_pricer.py`**

Add this new function in `wagerzon_odds/parlay_pricer.py` immediately after the existing `get_parlay_price()`:

```python
def get_combined_parlay_price(session: requests.Session, legs: list[dict],
                               amount: int = 10000) -> dict | None:
    """Call ConfirmWagerHelper for a cross-game parlay (legs from multiple games).

    Differs from get_parlay_price() in that each leg carries its own idgm.

    Args:
        session: Authenticated requests session
        legs: List of dicts with {idgm, play, points, odds} per leg
        amount: Bet amount for price query (defaults to 10000 for precision;
                falls back to 100 if MAXPARLAYRISKEXCEED — caller's responsibility)

    Returns:
        Dict with {win, decimal, american, amount} or None on error
    """
    # sel encoding: per-leg "play_idgm_points_odds"
    sel = ",".join(
        f"{l['play']}_{l['idgm']}_{l['points']}_{l['odds']}" for l in legs
    )
    detail_data = [
        {
            "Amount": str(amount),
            "RiskWin": "2",
            "TeaserPointsPurchased": 0,
            "IdGame": l["idgm"],
            "Play": l["play"],
            "Pitcher": 3,
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
        for l in legs
    ]

    try:
        resp = session.post(
            CONFIRM_URL,
            data={
                "IDWT": "0",
                "WT": "1",
                "amountType": "0",
                "open": "0",
                "sameAmount": "false",
                "sameAmountNumber": "0",
                "useFreePlayAmount": "false",
                "sel": sel,
                "detailData": json.dumps(detail_data),
            },
            timeout=15,
        )
        resp.raise_for_status()
        result = resp.json().get("result", {})

        details = result.get("details", [])
        if not details:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        win = details[0].get("Win", 0)
        risk = details[0].get("Risk", 0)
        if win <= 0 or risk <= 0:
            err = result.get("ErrorMsgKey") or result.get("ErrorMsg")
            if err:
                print(f"  API error: {err}")
            return None

        decimal_odds = 1 + win / amount
        american = round(win / amount * 100) if win > amount else round(-amount / win * 100)

        return {
            "win": win,
            "decimal": round(decimal_odds, 4),
            "american": american,
            "amount": amount,
        }
    except Exception as e:
        print(f"  Request error: {e}")
        return None
```

- [ ] **Step 4: Run tests — expect to pass**

```bash
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py" -v
```
Expected: PASS — 2 / 2.

- [ ] **Step 5: Commit**

```bash
git add wagerzon_odds/parlay_pricer.py "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py"
git commit -m "feat(wagerzon): cross-game ConfirmWagerHelper for combined parlays

Adds get_combined_parlay_price() — like get_parlay_price() but each leg
carries its own idgm. Required for combining two single-game parlays
into one 4-leg WZ ticket.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: DB schema migration for `placed_parlays`

**Files:**
- Create: `Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py`

- [ ] **Step 1: Write the failing test**

Append to `Answer Keys/MLB Dashboard/tests/test_combined_parlay.py`:

```python
import duckdb
import importlib.util


def _load_migration(path):
    spec = importlib.util.spec_from_file_location("migration", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_migration_adds_combo_columns(tmp_path):
    db = tmp_path / "test.duckdb"
    # Seed a minimal placed_parlays table matching the production schema
    con = duckdb.connect(str(db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY,
            game_id VARCHAR NOT NULL,
            home_team VARCHAR NOT NULL,
            away_team VARCHAR NOT NULL,
            game_time TIMESTAMP,
            combo VARCHAR NOT NULL,
            spread_line FLOAT,
            total_line FLOAT,
            fair_odds INTEGER,
            wz_odds INTEGER NOT NULL,
            edge_pct FLOAT,
            kelly_bet FLOAT NOT NULL,
            actual_size FLOAT,
            placed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            status VARCHAR DEFAULT 'pending'
        )
    """)
    con.close()

    migration_path = (
        Path(__file__).resolve().parents[1]
        / "migrations" / "001_combined_parlay_columns.py"
    )
    mod = _load_migration(migration_path)
    mod.run(str(db))

    con = duckdb.connect(str(db))
    cols = {row[0]: row[1] for row in con.execute("DESCRIBE placed_parlays").fetchall()}
    con.close()

    assert "is_combo" in cols
    assert "combo_leg_ids" in cols
    assert "parent_combo_id" in cols
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py::test_migration_adds_combo_columns" -v
```
Expected: FAIL — migration file doesn't exist.

- [ ] **Step 3: Implement the migration**

Create `Answer Keys/MLB Dashboard/migrations/__init__.py` (empty file) and `Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py`:

```python
"""One-shot migration: add combo tracking columns to placed_parlays.

Run once against the live mlb_dashboard.duckdb to make the dashboard combo-aware:

    python "Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py" \\
        "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
"""
from __future__ import annotations
import sys
import duckdb


def run(db_path: str) -> None:
    con = duckdb.connect(db_path)
    try:
        existing = {row[0] for row in con.execute("DESCRIBE placed_parlays").fetchall()}

        if "is_combo" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN is_combo BOOLEAN DEFAULT FALSE")
        if "combo_leg_ids" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN combo_leg_ids VARCHAR")
        if "parent_combo_id" not in existing:
            con.execute("ALTER TABLE placed_parlays ADD COLUMN parent_combo_id INTEGER")
    finally:
        con.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path-to-mlb_dashboard.duckdb>", file=sys.stderr)
        sys.exit(2)
    run(sys.argv[1])
    print(f"Migration 001 applied to {sys.argv[1]}")
```

- [ ] **Step 4: Run tests — expect to pass**

```bash
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py::test_migration_adds_combo_columns" -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/migrations/"
git commit -m "feat(mlb-dashboard): schema migration for combined parlay columns

Idempotent ALTER TABLE adding is_combo, combo_leg_ids, parent_combo_id.
Run once against the live mlb_dashboard.duckdb.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Backend endpoint — `POST /api/price-combined-parlay`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`
- Create: `Answer Keys/MLB Dashboard/combined_parlay.py` (Python helper for joint pricing)

- [ ] **Step 1: Write the failing test**

Append to `Answer Keys/MLB Dashboard/tests/test_combined_parlay.py`:

```python
def test_price_combined_parlay_endpoint(monkeypatch, tmp_path):
    # Set up isolated DB with two parlay opportunity rows + sizing settings
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    test_mlb_db = tmp_path / "mlb.duckdb"

    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            is_combo BOOLEAN DEFAULT FALSE, combo_leg_ids VARCHAR,
            parent_combo_id INTEGER, status VARCHAR DEFAULT 'pending'
        )
    """)
    con.execute("""
        CREATE TABLE sizing_settings (param VARCHAR PRIMARY KEY, value FLOAT)
    """)
    con.execute("INSERT INTO sizing_settings VALUES ('parlay_bankroll', 1000), ('parlay_kelly_mult', 0.5)")
    con.close()

    con = duckdb.connect(str(test_mlb_db))
    con.execute("""
        CREATE TABLE mlb_parlay_opportunities (
            parlay_hash VARCHAR, fair_dec FLOAT, wz_dec FLOAT, idgm INTEGER,
            spread_line FLOAT, total_line FLOAT, spread_price INTEGER,
            total_price INTEGER, combo VARCHAR
        )
    """)
    con.execute("""
        INSERT INTO mlb_parlay_opportunities VALUES
        ('hash_a', 4.32, 4.55, 100001, -1.5, 9.5, 110, -105, 'Home -1.5 + Over 9.5'),
        ('hash_b', 4.85, 5.10, 100002, -1.5, 7.5, 120, -110, 'Home -1.5 + Under 7.5')
    """)
    con.close()

    # Patch the server module's DB paths
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)
    monkeypatch.setattr(svr, "MLB_DB_PATH", test_mlb_db)

    # Stub out the WZ pricer so the test doesn't hit the real API
    def fake_get_combined_parlay_price(session, legs, amount=10000):
        return {"win": 18000, "decimal": 19.0, "american": 1800, "amount": 1000}
    monkeypatch.setattr(svr, "wz_get_combined_parlay_price", fake_get_combined_parlay_price)

    # Hit the endpoint
    client = svr.app.test_client()
    resp = client.post("/api/price-combined-parlay", json={
        "parlay_hash_a": "hash_a", "parlay_hash_b": "hash_b"
    })
    assert resp.status_code == 200, resp.get_data(as_text=True)
    body = resp.get_json()
    assert body["joint_fair_dec"] == pytest.approx(4.32 * 4.85, rel=0.001)
    assert body["wz_dec"] == 19.0
    assert "kelly_stake" in body
    assert "joint_edge" in body
```

- [ ] **Step 2: Run test — expect failure**

```bash
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py::test_price_combined_parlay_endpoint" -v
```
Expected: FAIL — endpoint missing or import errors.

- [ ] **Step 3: Create `combined_parlay.py` helper**

Create `Answer Keys/MLB Dashboard/combined_parlay.py`:

```python
"""Pure-Python helpers for combined parlay pricing in the MLB Dashboard.

Mirrors compute_combined_parlay_pricing() in Tools.R but in Python so the
Flask server doesn't need to call out to R for the live banner pricing.
"""
from __future__ import annotations


def joint_pricing(fair_dec_a: float, fair_dec_b: float, wz_dec: float,
                   bankroll: float | None = None,
                   kelly_mult: float | None = None) -> dict:
    """Joint pricing for two cross-game parlays combined into one ticket.

    Independence across games means joint fair = product of decimals.
    Returns joint fair odds, fair prob, edge, and (if bankroll/kelly_mult
    provided) the recommended Kelly stake in dollars.
    """
    joint_fair_dec = fair_dec_a * fair_dec_b
    joint_fair_prob = 1.0 / joint_fair_dec
    joint_edge = joint_fair_prob * wz_dec - 1.0

    kelly_stake = None
    if bankroll is not None and kelly_mult is not None:
        edge_fraction = joint_edge / (wz_dec - 1.0)
        kelly_stake = round(max(edge_fraction * kelly_mult * bankroll, 0.0), 2)

    return {
        "joint_fair_dec": round(joint_fair_dec, 4),
        "joint_fair_prob": round(joint_fair_prob, 6),
        "joint_edge": round(joint_edge, 4),
        "kelly_stake": kelly_stake,
    }
```

- [ ] **Step 4: Add endpoint to `mlb_dashboard_server.py`**

In `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`, near the top imports add:

```python
import time
from combined_parlay import joint_pricing
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "wagerzon_odds"))
from parlay_pricer import get_combined_parlay_price as wz_get_combined_parlay_price
```

Add a constant near `DB_PATH`:

```python
MLB_DB_PATH = Path(__file__).resolve().parents[1] / "mlb.duckdb"

# Combined parlay pricing cache: keyed on sorted (hash_a, hash_b) tuple.
# TTL prevents stale prices when WZ moves their lines. 60s is a safe default
# for the first version — the live banner on a 2nd-check click after a recent
# pricing returns instantly instead of waiting on a fresh ConfirmWagerHelper call.
_COMBO_PRICE_CACHE: dict[tuple[str, str], tuple[float, dict]] = {}
_COMBO_PRICE_TTL_SECONDS = 60
```

Add the new endpoint above the existing `/api/place-parlay` definition (around line 614):

```python
@app.route("/api/price-combined-parlay", methods=["POST"])
def price_combined_parlay():
    """Compute joint pricing + recommended Kelly stake for two parlay rows.

    Body: {"parlay_hash_a": "...", "parlay_hash_b": "..."}
    Returns: joint_fair_dec, joint_fair_prob, joint_edge, wz_dec, wz_win,
             kelly_stake, exact_payout (if WZ pricing succeeded).
    """
    data = request.json
    hash_a = data.get("parlay_hash_a")
    hash_b = data.get("parlay_hash_b")
    if not hash_a or not hash_b:
        return jsonify({"success": False, "error": "Missing parlay_hash_a or parlay_hash_b"}), 400
    if hash_a == hash_b:
        return jsonify({"success": False, "error": "Cannot combine a row with itself"}), 400

    # Pull the two source rows from mlb.duckdb
    mcon = duckdb.connect(str(MLB_DB_PATH))
    try:
        rows = mcon.execute("""
            SELECT parlay_hash, fair_dec, wz_dec, idgm,
                   spread_line, total_line, spread_price, total_price, combo, game_id
            FROM mlb_parlay_opportunities
            WHERE parlay_hash IN (?, ?)
        """, [hash_a, hash_b]).fetchall()
    finally:
        mcon.close()

    if len(rows) != 2:
        return jsonify({"success": False, "error": "One or both parlays not found"}), 404

    # Reject same-game combos (each game is already a single parlay row)
    if rows[0][9] == rows[1][9]:
        return jsonify({"success": False, "error": "Same-game combos not supported"}), 400

    by_hash = {r[0]: r for r in rows}
    a = by_hash[hash_a]; b = by_hash[hash_b]

    # Cache check (keyed on sorted hashes so order-of-selection doesn't matter)
    cache_key = tuple(sorted([hash_a, hash_b]))
    now = time.time()
    cached = _COMBO_PRICE_CACHE.get(cache_key)
    if cached and (now - cached[0]) < _COMBO_PRICE_TTL_SECONDS:
        return jsonify(cached[1])

    # Pull parlay sizing settings from dashboard DB
    dcon = duckdb.connect(str(DB_PATH))
    try:
        bankroll = dcon.execute(
            "SELECT value FROM sizing_settings WHERE param = 'parlay_bankroll'"
        ).fetchone()[0]
        kmult = dcon.execute(
            "SELECT value FROM sizing_settings WHERE param = 'parlay_kelly_mult'"
        ).fetchone()[0]
    finally:
        dcon.close()

    # Build legs for WZ pricing
    legs = [
        {"idgm": a[3], "play": _spread_play(a[6], a[4]), "points": str(a[4]), "odds": int(a[6])},
        {"idgm": a[3], "play": _total_play(a[7], a[5]),  "points": str(a[5]), "odds": int(a[7])},
        {"idgm": b[3], "play": _spread_play(b[6], b[4]), "points": str(b[4]), "odds": int(b[6])},
        {"idgm": b[3], "play": _total_play(b[7], b[5]),  "points": str(b[5]), "odds": int(b[7])},
    ]

    # Live WZ price (10000 first, fallback to 100)
    wz = wz_get_combined_parlay_price(_get_wz_session(), legs, amount=10000)
    if wz is None:
        wz = wz_get_combined_parlay_price(_get_wz_session(), legs, amount=100)
    if wz is None:
        return jsonify({"success": False, "error": "WZ pricing unavailable"}), 502

    pricing = joint_pricing(
        fair_dec_a=float(a[1]), fair_dec_b=float(b[1]), wz_dec=float(wz["decimal"]),
        bankroll=bankroll, kelly_mult=kmult,
    )

    response = {
        "success": True,
        "joint_fair_dec": pricing["joint_fair_dec"],
        "joint_fair_prob": pricing["joint_fair_prob"],
        "joint_edge": pricing["joint_edge"],
        "kelly_stake": pricing["kelly_stake"],
        "wz_dec": wz["decimal"],
        "wz_win": wz["win"],
        "amount": wz["amount"],
        "leg_a_combo": a[8],
        "leg_b_combo": b[8],
    }
    _COMBO_PRICE_CACHE[cache_key] = (now, response)
    return jsonify(response)


def _spread_play(odds: int, line: float) -> int:
    # Wagerzon play codes: 0=away spread, 1=home spread (heuristic from existing scraper)
    return 1 if line < 0 else 0


def _total_play(odds: int, line: float) -> int:
    # 2=over, 3=under — distinguished by sign convention; assume positive line = over
    # NOTE: Verify against parlay_pricer.py existing usage; may need to track over/under
    # explicitly in mlb_parlay_opportunities. Filed as TODO if heuristic breaks.
    return 2  # placeholder — task 6a below will refine using the combo string


# Helper for getting an authenticated WZ session (delegated to parlay_pricer.py)
from parlay_pricer import get_wz_session as _get_wz_session
```

- [ ] **Step 4a: Add `get_wz_session()` to `parlay_pricer.py` (one-time refactor)**

Find the existing CLI entrypoint in `wagerzon_odds/parlay_pricer.py` (likely under `if __name__ == "__main__":` or a `main()` function). It currently builds a `requests.Session()`, sets cookies/headers from a saved token file, and uses that session for all `get_parlay_price()` calls. Lift that initialization into a top-level reusable function:

```python
def get_wz_session() -> requests.Session:
    """Return an authenticated Wagerzon session.
    Reads token cache from the same path the CLI does. Raises RuntimeError
    if no valid token is available — caller should surface to user.
    """
    # ... lift existing token-load + session-init code here ...
```

Then update the CLI entrypoint to call `get_wz_session()` instead of doing the work inline. Verify the CLI still works:

```bash
python wagerzon_odds/parlay_pricer.py mlb --help
```

Expected: same output as before (no behavior change).

- [ ] **Step 4b: Resolve the over/under play-code heuristic**

The `_total_play()` placeholder above needs to use the actual combo string ("Home -1.5 + Over 9.5" vs "Home -1.5 + Under 7.5") to pick play=2 (over) or play=3 (under). Update `_total_play()` to take the combo string:

```python
def _total_play(combo: str) -> int:
    return 2 if " + Over " in combo else 3
```

And update the leg-building to pass the combo string. Same kind of fix for `_spread_play()` if the line-sign heuristic doesn't hold for "Away" spreads — confirm by inspecting an existing parlay row's combo field and `play` mapping in parlay_pricer.py.

- [ ] **Step 5: Run test — expect to pass**

```bash
pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py::test_price_combined_parlay_endpoint" -v
```
Expected: PASS.

- [ ] **Step 5a: Add a cache hit test**

Append to `tests/test_combined_parlay.py`:

```python
def test_price_combined_parlay_caches_within_ttl(monkeypatch, tmp_path):
    """Re-pricing the same pair within TTL must NOT call WZ a second time."""
    # Reuse the same setup as test_price_combined_parlay_endpoint
    # ... (same DB seeding + monkeypatching of svr.DB_PATH, MLB_DB_PATH) ...

    call_count = {"n": 0}
    def fake_wz(session, legs, amount=10000):
        call_count["n"] += 1
        return {"win": 18000, "decimal": 19.0, "american": 1800, "amount": 1000}
    monkeypatch.setattr(svr, "wz_get_combined_parlay_price", fake_wz)
    svr._COMBO_PRICE_CACHE.clear()  # ensure clean cache for the test

    client = svr.app.test_client()
    payload = {"parlay_hash_a": "hash_a", "parlay_hash_b": "hash_b"}
    r1 = client.post("/api/price-combined-parlay", json=payload)
    r2 = client.post("/api/price-combined-parlay", json=payload)
    assert r1.status_code == 200 and r2.status_code == 200
    assert call_count["n"] == 1, "WZ pricer should be called exactly once across two requests"
```

Run: `pytest "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py::test_price_combined_parlay_caches_within_ttl" -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/combined_parlay.py" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py"
git commit -m "feat(mlb-dashboard): /api/price-combined-parlay endpoint

Returns joint pricing + WZ exact payout + recommended Kelly stake for
two parlay rows. Rejects same-game combos and same-row pairs. Falls
back from amount=10000 to amount=100 on MAXPARLAYRISKEXCEED.

In-memory TTL cache (60s) keyed on sorted (hash_a, hash_b) tuple so
re-pricing the same pair returns instantly — banner pop-up on second
checkbox click stays under ~1ms after the first call.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Backend endpoint — `POST /api/place-combined-parlay`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_combined_parlay.py`:

```python
def test_place_combined_parlay_creates_combo_row(monkeypatch, tmp_path):
    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            actual_size FLOAT, status VARCHAR DEFAULT 'pending',
            placed_at TIMESTAMP, is_combo BOOLEAN DEFAULT FALSE,
            combo_leg_ids VARCHAR, parent_combo_id INTEGER
        )
    """)
    con.close()

    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)

    client = svr.app.test_client()
    resp = client.post("/api/place-combined-parlay", json={
        "combo_hash": "combo_hash_xyz",
        "parlay_hash_a": "hash_a",
        "parlay_hash_b": "hash_b",
        "wz_odds": 1810,
        "kelly_bet": 9.40,
        "actual_size": 9.40,
        "combo_label": "NYY @ BOS + LAD @ SD (4-leg)",
    })
    assert resp.status_code == 200

    con = duckdb.connect(str(test_dashboard_db))
    rows = con.execute("SELECT parlay_hash, is_combo, combo_leg_ids FROM placed_parlays").fetchall()
    con.close()
    assert len(rows) == 1
    assert rows[0][0] == "combo_hash_xyz"
    assert rows[0][1] is True
    import json as _json
    assert sorted(_json.loads(rows[0][2])) == ["hash_a", "hash_b"]
```

- [ ] **Step 2: Run test — expect failure**

Same `pytest` command. Expected: FAIL — endpoint missing.

- [ ] **Step 3: Add endpoint to `mlb_dashboard_server.py`**

Add this above (or below) `/api/place-parlay`:

```python
@app.route("/api/place-combined-parlay", methods=["POST"])
def place_combined_parlay():
    """Record a combined (cross-game) parlay placement.

    Inserts a single row into placed_parlays with is_combo=TRUE and
    combo_leg_ids = JSON list of the two source parlay_hashes.
    """
    data = request.json
    required = ["combo_hash", "parlay_hash_a", "parlay_hash_b",
                "wz_odds", "kelly_bet", "actual_size", "combo_label"]
    missing = [k for k in required if k not in data]
    if missing:
        return jsonify({"success": False, "error": f"Missing fields: {missing}"}), 400

    leg_ids_json = json.dumps([data["parlay_hash_a"], data["parlay_hash_b"]])

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            existing = con.execute(
                "SELECT parlay_hash, status FROM placed_parlays WHERE parlay_hash = ?",
                [data["combo_hash"]]
            ).fetchone()
            if existing and existing[1] == "pending":
                return jsonify({"success": False, "error": "Combo already placed"}), 409

            con.execute("""
                INSERT INTO placed_parlays (
                    parlay_hash, game_id, home_team, away_team, combo,
                    wz_odds, kelly_bet, actual_size, placed_at, status,
                    is_combo, combo_leg_ids
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 'pending', TRUE, ?)
            """, [
                data["combo_hash"],
                "COMBO",  # synthetic game_id for combos
                "(combined)", "(combined)",
                data["combo_label"],
                int(data["wz_odds"]),
                float(data["kelly_bet"]),
                float(data["actual_size"]),
                datetime.now().isoformat(),
                leg_ids_json,
            ])
        finally:
            con.close()
        return jsonify({"success": True, "message": "Combined parlay placed"})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500
```

- [ ] **Step 4: Run test — expect to pass**

Same `pytest` command. Expected: PASS — 4 / 4 (all tests in file).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" "Answer Keys/MLB Dashboard/tests/test_combined_parlay.py"
git commit -m "feat(mlb-dashboard): /api/place-combined-parlay endpoint

Inserts a single placed_parlays row with is_combo=TRUE and combo_leg_ids
referencing the two source parlay_hashes. Combo's own game_id is the
synthetic 'COMBO' to keep schema constraints happy.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Frontend — checkbox column on Parlay tab

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (`create_parlays_table()`, around line 249)

Note: this task is rendering plumbing; the JS-side selection state lives in Task 9 (it makes more sense for the checkbox HTML and the banner JS to land in the same change conceptually, but we split them for smaller commits).

- [ ] **Step 1: Add a hidden Sel column with checkbox HTML**

In the `columns = list(...)` block in `create_parlays_table()`, insert this colDef as the FIRST column:

```r
sel = colDef(
  name = "",
  minWidth = 30,
  align = "center",
  filterable = FALSE,
  sortable = FALSE,
  html = TRUE,
  cell = function(value, index) {
    row <- table_data[index, ]
    if (row$is_placed || isTRUE(row$is_combo)) {
      ''  # no checkbox if already placed
    } else {
      sprintf(
        '<input type="checkbox" class="combo-select" data-hash="%s" data-game-id="%s" onchange="onComboSelectChange(this)">',
        row$parlay_hash, row$game_id
      )
    }
  }
),
```

And add `sel = ""` as a column in the `table_data` (so reactable has something to render). Add at the start of the mutate chain:

```r
sel = "",
```

- [ ] **Step 2: Manual smoke check**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-combined-parlay/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R
open report.html
```
Expected: Parlay tab shows a checkbox in the leftmost column for each unplaced row. Clicking does nothing yet (no JS handler — that's Task 9).

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): add Sel checkbox column to Parlay tab

Checkboxes appear on every unplaced row. The onComboSelectChange handler
will be wired in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Frontend — Combined Parlay banner

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the HTML scaffolding around the Parlay tab)
- May also touch: any JS file the existing `placeParlay()` lives in (find it via `grep -r "placeParlay" "Answer Keys/MLB Dashboard"`)

- [ ] **Step 1: Inject the banner HTML above the Parlay table**

In `mlb_dashboard.R`, find the section that renders the Parlay tab (the call site of `create_parlays_table()`). Just before the table is rendered, add this HTML:

```r
HTML('
<div id="combined-parlay-banner" style="display:none; padding:10px 14px; background:#1f3a5f; border-left:3px solid #4a9eff; margin-bottom:10px; border-radius:4px;">
  <span id="combo-banner-status">Pricing combined ticket…</span>
  <button id="combo-place-btn" onclick="placeCombinedParlay()" style="display:none; float:right; background:#4a9eff; color:#fff; border:none; padding:6px 14px; border-radius:4px; cursor:pointer;">Place combined →</button>
</div>
')
```

- [ ] **Step 2: Add the JS in a `<script>` block at the bottom of the rendered HTML**

In the same file, near where the existing JS for `placeParlay` lives (find it with `grep -n "placeParlay" "Answer Keys/MLB Dashboard/mlb_dashboard.R"`), add:

```javascript
// Combined parlay state — at most 2 selected at a time, FIFO
window.comboSelections = [];
window.comboPricing = null;

function onComboSelectChange(checkbox) {
  const hash = checkbox.dataset.hash;
  const gameId = checkbox.dataset.gameId;
  if (checkbox.checked) {
    // FIFO: if 2 already selected, uncheck the oldest
    if (window.comboSelections.length >= 2) {
      const oldest = window.comboSelections.shift();
      const oldEl = document.querySelector(`.combo-select[data-hash="${oldest.hash}"]`);
      if (oldEl) oldEl.checked = false;
    }
    window.comboSelections.push({ hash, gameId });
  } else {
    window.comboSelections = window.comboSelections.filter(s => s.hash !== hash);
  }
  refreshComboBanner();
}

async function refreshComboBanner() {
  const banner = document.getElementById("combined-parlay-banner");
  const status = document.getElementById("combo-banner-status");
  const placeBtn = document.getElementById("combo-place-btn");

  if (window.comboSelections.length !== 2) {
    banner.style.display = "none";
    return;
  }
  // Same-game guard
  if (window.comboSelections[0].gameId === window.comboSelections[1].gameId) {
    banner.style.display = "block";
    status.textContent = "Same-game combos are already a single row — pick rows from different games.";
    placeBtn.style.display = "none";
    window.comboPricing = null;
    return;
  }

  banner.style.display = "block";
  status.textContent = "Pricing combined ticket…";
  placeBtn.style.display = "none";

  try {
    const res = await fetch("/api/price-combined-parlay", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        parlay_hash_a: window.comboSelections[0].hash,
        parlay_hash_b: window.comboSelections[1].hash,
      }),
    });
    const body = await res.json();
    if (!body.success) {
      status.textContent = "WZ rejected: " + (body.error || "unknown");
      return;
    }
    window.comboPricing = body;
    const edgeColor = body.joint_edge >= 0 ? "#5fd97a" : "#e57373";
    status.innerHTML = `<b>2 legs combined</b> · stake <b>$${body.kelly_stake}</b> · WZ pays <b>${body.wz_dec}</b> · edge <span style="color:${edgeColor}">+${(body.joint_edge * 100).toFixed(1)}%</span>`;
    placeBtn.style.display = "inline-block";
  } catch (err) {
    status.textContent = "Pricing failed — retry by re-checking a row.";
  }
}
```

- [ ] **Step 3: Manual smoke test**

Restart the Flask server (`mlb_dashboard_server.py`), refresh the dashboard, check 2 different-game rows. The banner should appear, show "Pricing…", then update with the joint odds and a Place Combined button (which doesn't do anything yet — that's Task 10).

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): combined parlay banner + selection JS

Banner div above the Parlay table, hidden by default. Checkbox onchange
triggers onComboSelectChange — FIFO state, same-game guard, fetch from
/api/price-combined-parlay and render joint odds, edge, recommended
stake. Place Combined button wired in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: Frontend — Place Combined button wiring

**Files:**
- Modify: same JS section in `mlb_dashboard.R` as Task 9

- [ ] **Step 1: Add `placeCombinedParlay()` JS function**

Append to the same `<script>` block as Task 9:

```javascript
async function placeCombinedParlay() {
  if (!window.comboPricing || window.comboSelections.length !== 2) return;
  const placeBtn = document.getElementById("combo-place-btn");
  placeBtn.disabled = true;
  placeBtn.textContent = "Placing…";

  // Synthetic combo_hash from the two leg hashes (deterministic)
  const sels = window.comboSelections.map(s => s.hash).sort();
  const comboHash = "combo_" + sels.join("_");

  try {
    const res = await fetch("/api/place-combined-parlay", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        combo_hash: comboHash,
        parlay_hash_a: sels[0],
        parlay_hash_b: sels[1],
        wz_odds: Math.round((window.comboPricing.wz_dec - 1) * 100),
        kelly_bet: window.comboPricing.kelly_stake,
        actual_size: window.comboPricing.kelly_stake,
        combo_label: window.comboPricing.leg_a_combo + " + " + window.comboPricing.leg_b_combo,
      }),
    });
    const body = await res.json();
    if (!body.success) {
      alert("Failed to place combined parlay: " + (body.error || "unknown"));
      placeBtn.disabled = false;
      placeBtn.textContent = "Place combined →";
      return;
    }
    // Reload page so R re-renders with the new combo's residual Kelly applied
    window.location.reload();
  } catch (err) {
    alert("Network error placing combined parlay");
    placeBtn.disabled = false;
    placeBtn.textContent = "Place combined →";
  }
}
```

- [ ] **Step 2: End-to-end manual test**

1. Refresh the dashboard
2. Check 2 different-game rows
3. Wait for banner to populate
4. Click "Place Combined"
5. Page reloads — verify:
   - The two source rows now show reduced Kelly (e.g., `$X` with `(residual after combo $Y)` annotation)
   - A new combo row appears in the placed bets section (the "(combined)" placeholder)
   - The two source rows' checkboxes are gone (since they're considered "active but reduced")
6. Refresh again — state persists

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): wire Place Combined Parlay button

POST /api/place-combined-parlay then reload so the R-side residual
Kelly application takes effect. Combo_hash is deterministic from sorted
leg hashes so re-clicking is idempotent.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 11: Documentation

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 1: Add a "Combined Parlay" section**

Append to `Answer Keys/MLB Dashboard/README.md`:

```markdown
## Combined Parlay (Wagerzon cash-efficiency)

The Parlay tab supports combining **two recommended parlay rows from different games** into a single 4-leg Wagerzon ticket. Useful when WZ balance is the binding constraint — the combined ticket needs far less cash than placing both source parlays separately.

### How to use
1. On the Parlay tab, check the leftmost checkbox on **two rows from different games**.
2. The "Combined Parlay" banner appears above the table with:
   - Joint fair odds (= product of the two `fair_dec`s, since cross-game legs are independent)
   - Wagerzon's exact 4-leg payout (live `ConfirmWagerHelper` call)
   - Joint edge and recommended Kelly stake (computed from your `parlay_bankroll` and `parlay_kelly_mult` settings)
3. Click **Place Combined →** to record the placement. The page reloads.
4. After reload, the two source rows' Kelly column shows the **conditional residual** — the optimal additional stake on each single given the combo is already placed. You can place those residuals as singles too if you want more exposure.

### When to use it
- Wagerzon balance is below the cost of placing both source parlays separately
- You're OK trading EV per dollar bet for cash efficiency

### When NOT to use it
- Wagerzon balance is plentiful — placing both as singles captures more EV
- Same-game combos: blocked, since each game's spread+total combo is already a single Parlay row

### Math
- Combined Kelly stake: `kelly_stake(joint_edge, wz_decimal_payout, parlay_bankroll, parlay_kelly_mult)`
- Conditional residual: numerical max-log-growth optimization over the 4-outcome joint distribution `(p_A * p_B, p_A * (1-p_B), (1-p_A) * p_B, (1-p_A) * (1-p_B))` with the combo stake fixed
- Implementation: `Answer Keys/conditional_kelly.R::conditional_kelly_residuals()` and `Answer Keys/MLB Dashboard/combined_parlay.py::joint_pricing()`

### One-time setup
After pulling this feature for the first time, run the schema migration once:
```bash
python "Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py" \
    "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```

### Out-of-scope (v1)
- N>2 leg combinations
- Combining a parlay with a single from the Bets tab
- Cash-budget-aware portfolio Kelly (set "available WZ balance" → dashboard solves)
```

- [ ] **Step 2: Commit**

```bash
git add "Answer Keys/MLB Dashboard/README.md"
git commit -m "docs(mlb-dashboard): document the combined parlay feature

How to use, when to use, when not to use, the math, and the one-time
schema migration step.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Pre-merge review checklist

Before merging `feature/mlb-combined-parlay` to `main`, run an executive engineer review of `git diff main..HEAD` (per CLAUDE.md). Specifically:

- [ ] **Combo deduplication** — re-clicking "Place Combined" with the same selections returns 409 (existing pattern), not a duplicate row
- [ ] **Same-game guard** — the API rejects same-game combos (Task 6 endpoint), AND the JS prevents the banner from showing the Place button for same-game (Task 9 frontend)
- [ ] **Already-placed guard** — checkboxes don't render for already-placed rows (Task 8); the API rejects if either source row is already placed as a single
- [ ] **Residual ≤ original** — the testthat suite from Task 1 already asserts this; verify still passing on full pipeline
- [ ] **Schema migration is idempotent** — running it twice doesn't error; existing rows get the default values (Task 5)
- [ ] **Negative residual handling** — `conditional_kelly_residuals()` floors at 0 (Task 1); table renders "$0 (residual after combo)" gracefully
- [ ] **WZ session leak** — the `_get_wz_session()` helper from Task 6b doesn't open a new session per request unbounded; verify reuse pattern matches existing `parlay_pricer.py`
- [ ] **DuckDB connections** — every `duckdb.connect()` is wrapped in `try/finally con.close()` (existing pattern)
- [ ] **No secrets in logs** — the WZ pricing code doesn't print authentication tokens

Document findings as ISSUES TO FIX vs ACCEPTABLE RISKS, fix all issues, then get explicit user approval to merge.

## Worktree cleanup

After merge:

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .worktrees/mlb-combined-parlay
git branch -d feature/mlb-combined-parlay
```

## Self-review notes

- **Spec coverage:** All in-scope items from the spec map to a task. Out-of-scope items (N>2 legs, cross-tab combos, budget Kelly) explicitly deferred per the spec.
- **Known sharp edges:**
  - Task 6 has TWO sub-steps (4a, 4b) flagged as needing local verification — over/under play-code mapping and WZ session helper. These can't be specified more precisely without inspecting the existing `parlay_pricer.py` CLI flow firsthand.
  - Task 3's residual computation in R reads `wz_odds` as American — needs verification that the schema's `wz_odds` column is American (it appears to be from the existing `place_parlay` endpoint's `int(data["wz_odds"])`). If it's already decimal, the conversion code in `apply_combo_residuals()` is a no-op or wrong; double-check before merge.
  - Task 9's JS selectors assume a single Parlay tab on the page; if the dashboard ever shows multiple parlay tables, the global `window.comboSelections` will collide.
- **Ambiguity check:** "Source rows stay active" in the spec means their checkboxes should still appear (so user can place residuals). Task 8 currently hides the checkbox after placement (`row$is_placed`). For combo legs the kelly is reduced but the row is NOT marked is_placed, so the checkbox stays — verify this works in the manual test in Task 10.
