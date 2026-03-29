# March Madness Simulator — 2027 Improvement Plan

Post-mortem from 2026 tournament. This document captures what worked, what didn't,
and a prioritized roadmap for next year's build.

## What Worked Well (2026)

- **Rcpp sim speed**: 50K sims in 1 second — no bottleneck on computation
- **Dynamic bracket**: C++ `rounds_won` approach correctly handles mid-tournament state
- **6-source composite**: Median of BPI/KenPom/Torvik/EvanMiya/TeamRankings/PowerRank is robust
- **Kalshi edge finder**: Real-time EV computation against live markets, correct fee model
- **Shiny dashboard**: 9 tabs, searchable, filterable, refreshable — good operational tool
- **ESPN API bracket fetch**: Reliable, no HTML scraping fragility

## Known Bugs to Fix

### 1. EvanMiya parser is brittle
Hard-coded 21-column pivot. When the Google Sheet format changes (it will), the parser silently fails.
**Fix**: Parse by column name, not position. Add validation: if Team column is numeric, skip source gracefully.

---

## Statistical Model Improvements

### A. Variable Game Variance (Priority: HIGH)

**Current**: Fixed `sd_margin = 11.2` for all games.

**Problem**: NCAA tournament games have different variance by round and matchup quality:
- 1-seed vs 16-seed: lower variance (outcome more certain)
- 8-seed vs 9-seed: higher variance (toss-up game)
- Elite 8+: best teams remaining, possibly lower variance

**Proposed**: Seed-dependent variance model:
```r
# Base SD adjusted by seed differential
seed_diff <- abs(seed_A - seed_B)
sd_game <- 11.2 - 0.15 * seed_diff  # Tighter spread for mismatches
sd_game <- max(sd_game, 8.0)         # Floor at 8 points
```

**Validation**: Backtest against 2015-2025 tournament margins. Compute actual SD by seed-pair
(1v16, 2v15, etc.) and fit the adjustment curve.

### B. Non-Normal Margin Distribution (Priority: MEDIUM)

**Current**: `rnorm(mean = diff, sd = 11.2)` — assumes margins are normally distributed.

**Reality**: Tournament margins have fatter tails than normal (more blowouts AND more upsets
than normal predicts). Historical 12-over-5 upset rate is ~35%, but our model may give ~28%.

**Proposed**: Use t-distribution with low degrees of freedom:
```r
actual_margin <- diff + sd_margin * rt(1, df = 7)  # df=7 gives fatter tails
```

**Validation**: Compare model upset rates to historical rates by seed matchup (1995-2025).
Target: model predictions within 2% of historical for each seed pair.

### C. Weighted Composite Rating (Priority: HIGH)

**Current**: Unweighted median of 6 sources.

**Problem**: Not all sources are equally predictive. KenPom has 20+ years of track record.
PowerRank is a newer, less-tested aggregator. Weighting them equally is suboptimal.

**Proposed**: Calibrate weights via historical prediction accuracy:
```r
# Calibrated weights (determine from 2015-2025 backtest)
weights <- c(KenPom = 1.5, BPI = 1.3, Torvik = 1.2,
             EvanMiya = 0.8, TeamRankings = 0.6, PowerRank = 0.5)
composite <- weighted.mean(ratings, weights, na.rm = TRUE)
```

**Calibration method**: For each historical tournament, compute each source's prediction accuracy
(log-loss or Brier score against actual outcomes). Weight inversely proportional to error.

### D. Closing Line as Prior (Priority: HIGH)

**Current**: Model uses only power ratings. Ignores market information entirely.

**Problem**: The market (Vegas closing lines, Kalshi prices) aggregates information from
thousands of bettors, sharp and recreational. Ignoring this is throwing away signal.

**Proposed**: Bayesian blend of model and market:
```r
# Market-implied probability from Kalshi/Pinnacle closing line
market_prob <- kalshi_price / 100

# Model probability from sim
model_prob <- sim_advancement_rate

# Blend: 40% model, 60% market (market is more efficient)
blended_prob <- 0.4 * model_prob + 0.6 * market_prob
```

The edge then comes from WHERE the model and market disagree — not from the raw model output.
This is the Nate Silver / 538 approach: your model adds value at the margins, not as a replacement.

### E. Strength of Schedule Adjustment (Priority: MEDIUM)

**Current**: Raw power ratings with no SOS adjustment.

**Problem**: A team from the Big 12 with a +15 rating played against much tougher competition
than a team from the Patriot League with the same rating. Tournament performance correlates
with regular-season schedule strength.

**Proposed**: Regress ratings toward conference strength:
```r
# Conference average rating
conf_avg <- ratings %>% group_by(conference) %>% summarise(avg = mean(composite_rating))

# Shrink toward conference average (Bayesian shrinkage)
shrink_factor <- 0.1  # 10% shrinkage
adjusted_rating <- composite_rating * (1 - shrink_factor) + conf_avg * shrink_factor
```

---

## Kalshi Edge Improvements

### F. Liquidity-Weighted Edge Scoring (Priority: HIGH)

**Current**: Edges sorted by EV% regardless of liquidity. A +150% EV edge on a 1¢ market
with no volume is useless.

**Proposed**: Expected profit metric:
```r
# Expected profit = EV% × executable size × liquidity score
liquidity_score <- min(1, volume / 10000) * (1 / (1 + spread/10))
expected_profit <- ev_pct * max_executable_size * liquidity_score
```

Sort by expected profit, not raw EV%.

### G. Model Confidence Filter (Priority: MEDIUM)

**Current**: Shows all edges regardless of model confidence.

**Problem**: Our model is weakest on:
- Extreme favorites/underdogs (calibration breaks down at tails)
- Partially completed rounds (remaining games affect count-based props)
- Teams with missing rating sources

**Proposed**: Add confidence column based on:
- Number of rating sources available for the team (6/6 = high, 3/6 = low)
- Distance from 50% fair value (closer to 50% = more confident)
- Spread between model and market (>20% divergence = flag for review)

### H. Historical Upset Calibration for Count Props (Priority: HIGH)

**Current**: Upset counts use bracket-position heuristics. R32+ upsets approximate
matchup pairs by seed rank within region.

**Problem**: The approximation can miscount when multiple seeds from the same bracket
half advance (e.g., if both the 3 and 6 seed advance from one half, pairing them
together is wrong — they couldn't have played each other).

**Proposed for 2027**: Track actual matchup graph in the C++ sim. Each game produces
a winner and a loser with known seeds. Count upsets directly from `winner_seed > loser_seed`.
This requires the C++ sim to output per-game results, not just per-team advancement.

---

## Infrastructure Improvements

### I. Backtest Framework (Priority: CRITICAL)

**Current**: No historical validation. We don't know if the model is calibrated.

**Proposed**: Build a backtest harness that:
1. Loads historical bracket + ratings for year Y (2015-2025)
2. Runs 50K simulations using pre-tournament ratings
3. Compares predicted probabilities to actual outcomes
4. Computes: Brier score, log-loss, calibration curve, upset accuracy by seed pair

**Deliverable**: A single script `backtest.R` that produces a calibration report:
```
Year | Brier Score | Log-Loss | R64 Upset Accuracy | Championship Correct?
2025 |    0.18     |   0.52   |       82%          |         No
2024 |    0.17     |   0.49   |       79%          |         Yes
...
```

### J. Rating Source Caching (Priority: LOW)

**Current**: Fetches all 6 rating sources on every dashboard restart (~30 seconds).

**Proposed**: Cache to DuckDB with TTL:
```r
# Only refetch if >1 hour old
cached <- dbGetQuery(con, "SELECT * FROM rating_cache WHERE fetched_at > NOW() - INTERVAL 1 HOUR")
if (nrow(cached) > 0) return(cached)
# Otherwise fetch fresh
```

### K. Automated Google Sheet Refresh (Priority: MEDIUM)

**Current**: KenPom/Torvik/EvanMiya require manual Google Sheet updates.

**Proposed**: Direct API or scraping:
- KenPom: Use `kenpomR` package or direct API if available
- Torvik: Already have `fetch_torvik.py` Playwright scraper as fallback
- EvanMiya: Investigate if evanmiya.com has a public API

---

## Parameter Sensitivity Analysis (Do Before 2027 Tournament)

Run the following experiments on 2025 tournament data:

| Experiment | Parameter | Range | Metric |
|------------|-----------|-------|--------|
| Game variance | `sd_margin` | 9.0 - 13.0 | Brier score |
| Rating update speed | `beta1` | 0.05 - 0.20 | Championship accuracy |
| Temporal decay | `beta2` | 0.0 - 0.10 | Late-round accuracy |
| Composite method | median vs weighted mean | N/A | Overall log-loss |
| Distribution | normal vs t(df=5) vs t(df=10) | N/A | Upset accuracy |
| Sim count | 10K vs 50K vs 100K | N/A | Estimate stability |

---

## 2027 Pre-Tournament Checklist

- [ ] Run backtests on 2020-2026 tournaments
- [ ] Calibrate source weights from backtest
- [ ] Implement seed-dependent variance
- [ ] Test t-distribution vs normal
- [ ] Build closing-line blending pipeline
- [ ] Add liquidity scoring to Kalshi edges
- [ ] Automate EvanMiya data fetch (replace manual Google Sheet)
- [ ] Add model confidence column to Kalshi Edges tab
- [ ] Cache rating fetches in DuckDB
- [ ] Run parameter sensitivity sweep
- [ ] Validate upset rates against historical base rates
- [ ] Update `TEAM_NAME_FIXES` for any new D1 programs
