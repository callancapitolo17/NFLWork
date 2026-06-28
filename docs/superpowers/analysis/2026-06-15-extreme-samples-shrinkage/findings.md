# Extreme Samples Overhaul — Findings (2026-06-15)

## The problem
The Elihu/Feustel matched-sample pricer estimates P(cover/over/win) as the
empirical frequency in ~N historically-similar games. For **extreme lines**
(sparse regions: total>p99≈167.5, |spread|>p99≈32), `mean_match` cannot hold
the target, so the shrink-loop **collapses `final_N`** to force a match. The
collapsed sample is small + biased → overconfident probability → phantom EV →
oversized Kelly bet. These are the "extreme samples."

## Research conclusion (statistically correct fix)
Estimating a probability at the edge of / outside the support of your data.
No estimator manufactures absent information. Correct response: **detect the
sparsity and defer to the market in proportion to how little local data you
have.** Two composable tools:
1. **Hard support guard** — when the local sample is too thin/biased to trust
   (sample collapse / line out of support), abstain or force p̂ = market.
2. **Beta-Binomial shrinkage toward the market** — p̂=(x+n₀·p_mkt)/(n+n₀)
   = w·(x/n)+(1−w)·p_mkt, w=n/(n+n₀); graceful degradation as n thins.
3. **Baker-McHale Kelly** — multiply fraction by α=edge²/(edge²+Var(p̂)).
Full citations + math in `docs/.../2026-06-15-...` research report (in chat).

## Validation protocol (anti-overfit)
- Temporal split: TRAIN = games ≤ 2025-01-11 (6,102 games), TEST = later
  (5,622). 11,724 CBB games, 160,064 H1 predictions. n₀/thresholds fit on
  TRAIN, evaluated frozen on TEST. Bets graded into CBB H1 closing odds
  (conservative: closing line is near-unbeatable, so this UNDER-states edge).
- Primary evidence = proper scores (log-loss/Brier) + calibration (ECE,
  reliability) + Diebold-Mariano; ROI is corroboration (high variance).

## Results (held-out TEST)

### 1. Hard support guard — the clean, robust win
`is_extreme = |spread|>p99 OR total>p99 OR final_N<0.5·N` (a-priori rule).
Decomposition shows toxicity is **sample collapse**, not line value:

| +EV>5% subset            |   n  |   ROI   | win% |
|--------------------------|-----:|--------:|-----:|
| sample collapsed (<0.5N) |  325 | −18.93% | 37.2 |
| NOT collapsed            | 4443 |  +3.92% | 51.2 |
| RAW (all, production)    | 4768 |  +2.16% | 50.3 |
| RAW − extremes (guard)   | 4341 |  +4.10% | 51.4 |

Abstaining on ~4% of games (the collapsed/out-of-support tail) ~doubles
held-out ROI (2.16→4.10%) and adds +72% PnL (+2069→+3551). `final_N` is
already returned by `run_answer_key_sample`; threshold is not ROI-tuned.

### 2. Beta-Binomial shrinkage — clean calibration win, ROI suggestive
- ECE: raw 0.0081 → shrunk 0.0041 (book 0.0033) — overconfidence halved.
- Extreme-subset log-loss: raw 0.6653 → shrunk 0.6464 (book 0.6457) — 96% of
  the gap to the (sharp) book closed.
- log-loss/Brier optimum stable TRAIN↔TEST (n₀≈300–500) — not overfit.
- ROI-vs-n₀ (TEST) has a broad plateau: moderate n₀≈30–75 lifts ROI
  2.36→3.5–4.5% AND raises PnL. BUT no simple TRAIN criterion uniquely
  selects that moderate n₀ (train-PnL picks n₀=0, train-ROI picks n₀=500),
  so the ROI lift is **suggestive, not cleanly OOS-proven**. n₀≈50 is
  justified on first principles (prior ≈ min_N: dense games barely move,
  w≈0.9; collapsed samples pulled hard, w≈0.55) + the calibration win.

### Honest summary
- **High confidence:** sample-collapse guard removes a genuinely toxic
  (−19% ROI) bet population; nearly doubles holdout ROI; least overfittable.
- **Medium confidence:** shrinkage toward market cleanly improves calibration
  / tames overconfidence (graceful degradation); ROI benefit real-direction
  but not provable from train selection alone.
- **Caveat:** test = bets into CBB H1 *closing* lines (pessimistic vs live
  soft-book edge). Effect should be re-confirmed on MLB / live CLV.

## Ablation — which piece earns its keep (TEST holdout, settled at real odds)

Reproducible block: `CBB_Backtest_Shrinkage_NoOverfit.R` "PRODUCTION-MIRROR
ABLATION" (160,064 preds / 11,724 games; bets = +EV>5% at CBB H1 closing odds;
guard = collapse-only `final_N<0.5N`, matching production exactly).

| Version | Bets | Staked | ROI | Profit |
|---|---:|---:|---:|---:|
| RAW (today) | 4,770 | $95,978 | 2.18% | +$2,092 |
| **#1 Guard** (skip collapsed samples) | 4,445 | $88,577 | 3.94% | +$3,493 |
| **#1 + #3 Baker-McHale (SHIPPED)** | 4,445 | $65,733 | **4.11%** | +$2,703 |

- **#1 (guard) is the dominant win** (+$2,092 → +$3,493 PnL; ROI 2.18→3.94%). It
  drops exactly **325 bets** — *precisely* the collapsed-sample population the
  decomposition flagged at −18.9% ROI. Surgical.
- **#2 (shrink) was DROPPED**: it *lowers* total profit (trims good bets too),
  barely moves ROI, and overlaps the existing SGP `blend_dk_with_model` (a fixed
  50/50 model↔DK blend). Not worth the complexity.
- **#3 (Baker-McHale Kelly), dimensionally corrected** trims stake on noisy
  edges (a 2.6pp edge on n=200 matched games is within 1 SE). It **raises ROI
  (3.94→4.11%) but lowers PnL ~23% ($3,493→$2,703)** by shrinking total stake
  26%. User chose to ship it (2026-06-28): max ROI / lowest variance over max
  PnL. (An earlier first-cut fed the *EV-edge* against a *probability* variance,
  under-shrinking by ~decimal_odds² — same direction, far milder; corrected to
  probability-space edge `prob − breakeven` so edge and variance share units.)

## Implementation (worktree `extreme-samples-shrinkage`) — shipped = #1 + #3

Centralized, one kill switch. In `Tools.R`:
- Config: `EXTREME_GUARD_ENABLE` (master kill switch), `EXTREME_GUARD_FLOOR=0.5`.
- `run_answer_key_sample()` also returns `target_N` + `low_confidence`
  (`final_N < FLOOR*N`) — additive, backward compatible.
- `bakermchale_alpha(edge, var)`; `build_sample_meta(samples)` (ends with
  `distinct(id)` so the downstream join can never fan a bet row into dups).
- `apply_extreme_samples_correction(bets, sample_meta)`: per bet, abstain
  (bet_size 0, ev NA → filtered) on sample collapse; else scale the existing
  Kelly stake by the dimensionally-correct Baker-McHale factor
  `alpha = pedge²/(pedge² + p(1-p)/n_eff)`, where `pedge = prob − breakeven` is
  the probability edge (breakeven recovered self-contained as `prob/(ev+1)`).
  Only removes/shrinks, never adds. On MLB, touches only `edge_source=="model"`
  rows (market/both pass through). Rows with no sample meta pass through.

`CBB.R`, `MLB.R`, `NFLAnswerKey2.0.R`: call it on `all_bets_combined` then
re-filter `ev >= EV_THRESHOLD`, before `adjust_kelly_for_correlation`.

### Verification (offline)
- Unit/integration: collapsed→abstain (size 0, ev NA), BM trims dense bets,
  no-meta passthrough, **schema preserved** (dashboards safe), no temp leak.
  `edge_source` gate verified (market/both unchanged). Corrected alpha unit-test:
  ≤ old alpha everywhere (strictly more conservative), in [0,1], breakeven
  recovery `prob/(ev+1)` == `1/decimal_odds` exact.
- TEST holdout (re-priced 2026-06-28, 160,064 preds): RAW ROI **2.18%** / +$2,092
  → guard **3.94%** / +$3,493 → **#1+#3 4.11%** / **+$2,703**. All scripts parse.

### Limitations / not covered
- **Live pipeline run + dashboard render not done offline** (needs Odds API +
  scrapers) — run `run.py cbb` / `mlb` once and eyeball the dashboard.
- **Kalshi bots NOT covered**: `kalshi_mlb_rfq` reads raw `mlb_game_samples`;
  `kalshi_mlb_mm` is book-consensus-only — they price independently of the
  corrected bet list.

## Harness (worktree only — measurement instruments, not production)
- `Answer Keys/CBB Answer Key/CBB_Backtest_Shrinkage_AB.R` — prices games, saves preds.
- `..._Analysis.R` — temporal split, n₀ fit, scores/calibration/DM/ROI sweep.
- `..._NoOverfit.R` — train-chosen-n₀ confirmation + EV-bucket/extreme breakdowns
  + the reproducible PRODUCTION-MIRROR ABLATION (the headline table above).
