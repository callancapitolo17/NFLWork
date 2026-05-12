# FG `-1.5` Favorite Correlated Parlay Leak — Findings

**Date:** 2026-05-11  
**Branch:** `analysis/fg-fav-leak-2026-05-11`  
**Data:** 627 MLB correlated parlays (2026-03-28 → 2026-05-10) from `bet_tracking.xlsx`, joined to dashboard `placed_parlays` for stored model `edge_pct` where available (~119 matched).

## TL;DR

The FG `-1.5` favorite + over leak (-$1,575, -11.6% ROI on 222 bets) is concentrated almost entirely in one sub-slice: **bets where the total line ≤ 7.0** (n=33, ROI -73%, bootstrap 95% CI [-97%, -38%] fully below zero). The other three total-line buckets (7-8, 8-9, 9+) are statistically indistinguishable from break-even.

The mechanism is **NOT** correlation factor overstatement (H6: model corr is mild, ~1.04 at low totals). It is **marginal-probability overestimation** of the over leg at low totals: the simulator predicts ~28% joint hit rate at total_line=6.5, but realized rate at total ≤7 is 9.4% — a 19pp gap. Most likely cause: the resampled distribution overstates how often games with low-total lines go over, or WZ's over line at low totals is sharper than the model believes.

**Recommended action:** Skip FG `-1.5` favorite + over parlays when the total line is ≤ 7.0. On the observed sample this would have eliminated -$1,721 of losses while leaving the other +EV combos untouched. Plus follow-up: investigate why the simulator overstates over probability at low totals.

## Baseline (Task 2)

| Slice | n | ROI | 95% CI ROI | Hit Rate | Break-even |
|---|---:|---:|---|---:|---:|
| FG-fav-1.5 (all) | 218 | -11.63% | [-35.08%, +14.04%] | 24.41% | 22.55% |
| FG-fav-1.5 (ex-bad-week) | 200 | -8.38% | [-32.99%, +18.89%] | 25.13% | 22.60% |
| FG-fav-1.5 + over | 217 | -11.95% | [-35.03%, +13.79%] | 24.06% | 22.53% |
| FG-fav-1.5 + under | 1 | +247.00% | [+247.00%, +247.00%] | 100.00% | 28.82% |
| FG-dog+1.5 | 74 | +8.85% | [-32.87%, +54.33%] | 33.80% | 28.57% |
| F5 half | 325 | +7.03% | [-14.97%, +30.05%] | 28.25% | 26.40% |

Note: hit rate is slightly **above** break-even for FG-fav-1.5 — the loss is in the odds distribution, not raw win frequency.

## H1 — Over vs Under within FG-fav-1.5

- `+over`: n=217, ROI -11.95%, hit_rate 24.06%, CI [-35.03%, +13.79%]
- `+under`: n=1 — **insufficient sample** (1 settled).
- 2-proportion test: skipped — under arm has insufficient sample (n_under < 5)

**Verdict:** Moot. The under arm of FG-fav-1.5 is essentially never bet (the model rarely flags it as +EV). The entire FG-fav-1.5 leak universe is the over arm.

## H2 — Calibration on placed_parlays subset

- n_settled = 30, expected wins = 7.06, observed wins = 7, gap = -0.06
- By model_prob quartile:
  - n=8: model_prob 21.56%, actual hit 12.50%, gap +9.1 pp
  - n=7: model_prob 22.52%, actual hit 28.57%, gap -6.1 pp
  - n=7: model_prob 23.90%, actual hit 14.29%, gap +9.6 pp
  - n=8: model_prob 26.05%, actual hit 37.50%, gap -11.4 pp

**Verdict:** On the small matched subset (n=30), the model is well-calibrated overall (7.06 expected vs 7 observed wins). Per-bin gaps are noisy but no monotonic over- or under-prediction pattern. The leak is NOT a calibration miss on probability.

## H3 — Total-line stratification (SMOKING GUN)

| total bucket | n | wins | losses | ROI | hit rate | break-even | 95% CI ROI |
|---|---:|---:|---:|---:|---:|---:|---|
| ≤7 | 33 | 3 | 29 | -72.84% | 9.38% | 22.11% | [-97.00%, -37.54%] |
| 7-8 | 119 | 33 | 83 | +1.01% | 28.45% | 22.64% | [-33.71%, +39.72%] |
| 8-9 | 48 | 12 | 36 | +13.95% | 25.00% | 22.38% | [-43.47%, +80.25%] |
| 9+ | 18 | 4 | 13 | -24.71% | 23.53% | 23.27% | [-77.46%, +51.01%] |

**Verdict:** The `≤7` bucket is catastrophic and its CI does not include zero — the leak is real, not variance. The other buckets are noisy break-even. This single bucket accounts for **more than 100%** of the FG-fav-1.5 loss.

## H4 — Effective WZ shave reality check

- n = 214, model assumes shave = 0.989
- Mean realized shave = 0.9848, median = 0.9909, p10 = 0.9876, p90 = 0.9961
- Gap (model - realized mean) = +0.0042

**Verdict:** Shave matches model assumption almost exactly (median 0.9909 vs assumed 0.989). NOT the leak.

## H5 — Push EV accounting

- n_pushes = 5, push_rate = 0.0229
- Counterfactual net if WZ collapsed pushes to single leg: $0.00

**Verdict:** Push rate is only 2.3%. WZ refunds full stake on push. Opportunity cost from push-collapse is essentially zero. NOT the leak.

## H6 — Sample correlation vs realized

Sample-derived correlation factor and predicted joint probability by total line (from `mlb_game_samples`):
| total_line | mean corr (home fav) | mean corr (away fav) | predicted joint (home) | predicted joint (away) |
|---:|---:|---:|---:|---:|
| 6.5 | 1.039 | 1.048 | 0.285 | 0.251 |
| 7.5 | 1.084 | 1.119 | 0.254 | 0.228 |
| 8.5 | 1.034 | 1.082 | 0.212 | 0.193 |
| 9.5 | 1.090 | 1.150 | 0.180 | 0.165 |

Realized joint hit rate by total bucket (from our actual bets):
| total bucket | n | realized joint hit rate |
|---|---:|---:|
| ≤7 | 32 | 9.38% |
| 7-8 | 115 | 27.83% |
| 8-9 | 48 | 25.00% |
| 9+ | 17 | 23.53% |

**Critical comparison:** at total_line=6.5 the model's predicted joint hit prob (samples) is ~28%; realized rate at total ≤7 (our bets) is only 9.4%. Gap: ~19 pp.

**Verdict:** Correlation factor (1.04-1.15 range) is mild and consistent — NOT overstated. The leak is in the **underlying marginal probability** the simulator assigns to OVER at low totals. Either the resampled distribution overstates how often low-scoring games go over their (low) totals, or WZ's over line at total ≤7 is sharper than the simulator expects.

## Ranking of Hypotheses by Evidence Strength

1. **H3 (total-line stratification): CONFIRMED.** CI [-97%, -38%] for the ≤7 bucket excludes zero. n=33 is modest but the effect size is large enough that this is signal, not noise. This is the only hypothesis with a CI fully on one side of zero.
2. **H6 (marginal probability mismatch at low totals): SUPPORTING.** Model-predicted joint at total=6.5 (~28%) is ~3× higher than realized (~9.4%) at total ≤7. Cannot quantify per-bet gap rigorously without re-running the pricer, but the directional evidence is consistent with H3 and identifies the mechanism (over-leg marginal, not correlation).
3. **H4, H5: ELIMINATED.** Shave matches model; push EV is negligible. Both ruled out as drivers.
4. **H1: MOOT.** Under arm has 1 bet; cannot test.
5. **H2: NEUTRAL.** Model calibration on n=30 is fine overall, but this subset over-represents games where the bet was placed (selection effect). Doesn't contradict H3/H6.

## Recommended Action

**Primary recommendation:** Skip FG `-1.5` favorite + over parlays when total_line ≤ 7.0.

Implementation: add a filter in `Answer Keys/mlb_correlated_parlay.R`'s combo-construction loop that drops Home Spread + Over / Away Spread + Over combos where `row$total_line <= 7.0` for FG. (Spread + Under and dog combos remain eligible.)

**Expected impact on the observed sample:** -$1,721 of losses avoided (the ≤7 bucket's P&L), while leaving the other three buckets — net +$659 — untouched. Net P&L improvement on this slice: +$1,721 if perfectly executed.

**Caveats / known limits:**
- n=33 in the ≤7 bucket is small. The point estimate -73% ROI is noisy; CI [-97%, -38%] is wide. The directional conclusion (sub-slice is unprofitable) is robust, but the precise magnitude is uncertain. Re-evaluate after another month of data.
- The mechanism (H6 marginal-probability mismatch) implies that the simulator's run-distribution at low-total games is wrong. A more thorough fix would investigate why and recalibrate the underlying samples — but the filter recommendation is the immediate +EV move.
- Generalization: it is possible other slices (F5, FG-dog) also have a low-total bleed but we didn't subset-test them. Worth a follow-up stratification.
- This analysis is conditional on the WZ shave assumption holding (H4 confirmed it does, within 0.4pp).

**Secondary recommendation:** Investigate the simulator's over-distribution at low-total games. Specifically:
- For games where MLB Odds API posts a total ≤ 7.0, does `mlb_game_samples` predict over% above the over leg's implied break-even prob from WZ's line?
- If yes, the upstream sample generation is biased toward higher-scoring outcomes than current MLB delivers. Possible fixes: (a) restrict historical samples to 2024+ for low-total games, (b) introduce a total-line conditional adjustment in the resampler, (c) skip parlay generation when the simulator's over-prob disagrees with the line by >X%.

## Stop Conditions Reached

Per the plan's stop conditions:
- H3 found one total-line bucket with CI fully below zero → **actionable filter triggered** (skip total ≤ 7).
- H4 + H5 ruled out as drivers (shave OK, push EV ~$0).
- H6 explains the mechanism (over-marginal mispriced at low totals) without contradicting H3.
- No hypothesis tested inconclusive across the board, so the "need 2x more bets" exit condition is not met.