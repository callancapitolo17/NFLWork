# March Madness Simulator — 2027 Improvement Plan

Post-mortem from 2026 tournament. Prioritized roadmap for building a
genuinely sharp tournament model.

## What Worked Well (2026)

- **Rcpp sim speed**: 50K sims in 1 second — computation is not the bottleneck
- **Dynamic bracket**: C++ `rounds_won` approach handles mid-tournament state correctly
- **6-source composite**: Median of BPI/KenPom/Torvik/EvanMiya/TeamRankings/PowerRank is robust to outliers
- **Kalshi edge finder**: Real-time EV after fees, live refresh, sorted by edge size
- **Shiny dashboard**: Team advancement, survivor, seed props, custom queries, Kalshi edges — all in one place
- **ESPN API bracket fetch**: Reliable, auto-detects tournament state

## What Didn't Work Well

- **No calibration**: We never validated whether the model's probabilities match reality. A model
  that says "Duke has a 25% championship chance" is useless if we don't know whether teams it
  rates at 25% actually win ~25% of the time historically.
- **Fixed game variance**: Every game used sd=11.2 regardless of matchup quality. 1v16 games
  are much more predictable than 8v9 games, but the model treated them identically.
- **No market integration**: The model used only power ratings, ignoring the single most
  informative signal: the betting market's closing line.
- **Upset undercounting**: The `abs(rating_change)` bug inflated favorites throughout the
  tournament. Fixed now, but we never validated upset rates against historical base rates.

---

## The Core Problem: Calibration

Everything below is secondary to this: **we don't know if the model is calibrated**.

Calibration means: when the model says a team has a 70% chance of advancing, teams in that
bucket should actually advance ~70% of the time across historical tournaments. Without this,
every other improvement is guesswork — we're tuning parameters against vibes, not data.

### Step 1: Build a Backtest Framework (DO THIS FIRST)

**Why this matters**: Every serious quant model — from Renaissance Technologies to Nate Silver's
538 — is built on backtesting. You can't improve what you can't measure. The backtest tells you:
- Is the model overconfident? (predicts 80% but teams only win 60%)
- Is it underconfident? (predicts 50% but teams win 70%)
- Where does it break? (good at R64, bad at F4?)
- Which parameters actually matter? (does sd=10 beat sd=12?)

**What to build**: A script `backtest.R` that:

```r
# For each historical tournament (2015-2026):
for (year in 2015:2026) {
  # 1. Load that year's pre-tournament power ratings
  #    (KenPom historical archive, BPI historical, etc.)
  ratings <- load_historical_ratings(year)

  # 2. Load that year's actual bracket and results
  bracket <- load_historical_bracket(year)
  actual_results <- load_actual_results(year)

  # 3. Run 50K sims using those ratings
  sim_probs <- run_simulation(ratings, bracket, n_sims = 50000)

  # 4. Compare: for each team, what did we predict vs what happened?
  comparison <- sim_probs %>%
    left_join(actual_results, by = "team") %>%
    mutate(
      predicted_r32 = Round_32,
      actual_r32 = as.integer(made_round_32),
      predicted_champ = Champion,
      actual_champ = as.integer(won_championship)
    )

  # 5. Score the model
  brier_score <- mean((comparison$predicted_champ - comparison$actual_champ)^2)
  log_loss <- -mean(actual * log(predicted) + (1-actual) * log(1-predicted))
}
```

**Data sources for historical ratings**:
- KenPom: kenpom.com archive (subscription) or `cbbdata::cbd_kenpom_ratings(year)`
- BPI: ESPN API supports historical years
- Torvik: barttorvik.com has full historical archive
- Brackets: ESPN API works for past years, or use `hoopR::espn_mbb_tournament(year)`

**Output**: A calibration table and plot:
```
Year | Brier | Log-Loss | R64 Upset Rate (pred/actual) | Champ in Top 5?
2026 | 0.18  | 0.52     | 28% / 34%                    | Yes
2025 | 0.17  | 0.49     | 30% / 31%                    | No
2024 | 0.19  | 0.55     | 25% / 38%                    | Yes
```

Plus a **calibration curve**: bucket predicted probabilities into deciles (0-10%, 10-20%, etc.)
and plot predicted vs actual win rate. A perfect model lies on the diagonal.

### Step 2: Use the Backtest to Guide Every Other Improvement

Once you have the backtest, EVERY parameter change becomes testable:
- "Does sd=10 beat sd=12?" → Run both, compare Brier scores across 10 years
- "Does weighted mean beat median?" → Run both, compare log-loss
- "Does the t-distribution improve upset accuracy?" → Run both, compare R64 upset rates

Without the backtest, you're flying blind. With it, you have a scientific framework.

---

## Statistical Model Improvements

### A. Variable Game Variance

**Current**: Fixed `sd_margin = 11.2` for all games.

**Why it's wrong**: The standard deviation of game margins is NOT constant. Empirically:
- 1 vs 16: actual SD ≈ 14-16 (blowouts vary widely)
- 8 vs 9: actual SD ≈ 10-11 (competitive games cluster tighter)
- Elite 8+: actual SD ≈ 9-10 (best teams, lower variance)

But there's a subtlety: the SD we use isn't the raw margin SD — it's the **residual** SD
(actual minus predicted). If our ratings are good, the residual SD should be lower for
well-predicted games. The question is whether residual SD varies by matchup type.

**How to measure this**:
```r
# From historical data (2015-2025):
historical_games <- load_all_tournament_games()

historical_games %>%
  mutate(
    predicted_margin = rating_fav - rating_dog,
    residual = actual_margin - predicted_margin,
    seed_diff = abs(seed_fav - seed_dog),
    round = tournament_round
  ) %>%
  group_by(seed_diff_bucket = cut(seed_diff, breaks = c(0, 3, 7, 11, 16))) %>%
  summarise(
    n_games = n(),
    residual_sd = sd(residual),
    upset_rate = mean(actual_margin < 0)
  )
```

If residual SD varies significantly by seed differential, implement:
```r
sd_game <- base_sd + slope * seed_diff  # fit from data
```

If it doesn't vary much (within ±1 point), keep the fixed SD — simpler is better.

**Key insight**: Don't assume the answer. Measure it first, then decide.

### B. Margin Distribution Shape

**Current**: Normal distribution (`rnorm`).

**Why it might be wrong**: Normal distributions underpredict extreme events. If the true
distribution has fatter tails (more blowouts and more upsets than normal predicts), we'll
systematically undercount upsets and overcount chalk.

**How to test**:
```r
# Compute historical residuals
residuals <- actual_margin - predicted_margin  # across all tournament games

# Test for normality
shapiro.test(residuals)  # p < 0.05 → not normal

# Compare tail behavior
# What fraction of residuals exceed 2 SD?
observed_tail <- mean(abs(residuals) > 2 * sd(residuals))
expected_normal <- 2 * pnorm(-2)  # 4.6%

# If observed_tail > 6%, consider t-distribution
# Fit degrees of freedom:
library(MASS)
fit <- fitdistr(residuals / sd(residuals), "t")
optimal_df <- fit$estimate["df"]
```

If `optimal_df < 15`, the t-distribution is materially better. If `df > 30`, stick with normal.

**Implementation** (in C++ for speed):
```cpp
// Replace: R::rnorm(diff, sd_margin)
// With:    diff + sd_margin * R::rt(optimal_df)
double margin = diff + sd_margin * R::rt(optimal_df);
```

### C. Weighted Composite Rating

**Current**: Unweighted median of 6 sources.

**Why the median is decent but not optimal**: The median is robust — if one source is garbage,
it gets ignored. But it also ignores the fact that some sources are consistently more predictive
than others. KenPom has been the gold standard for 20 years. PowerRank is a newer aggregator
with less track record.

**How to calibrate weights**:

For each historical year (2015-2025):
1. Get each source's pre-tournament ratings
2. For each tournament game, compute each source's predicted margin
3. Compute each source's prediction error (RMSE of residuals)
4. Weight = 1 / RMSE² (inverse variance weighting — the statistically optimal approach)

```r
# Per-source RMSE across all historical tournament games
source_rmse <- historical_games %>%
  group_by(source) %>%
  summarise(rmse = sqrt(mean(residual^2)))

# Inverse variance weights
source_rmse$weight <- 1 / source_rmse$rmse^2
source_rmse$weight <- source_rmse$weight / sum(source_rmse$weight)  # normalize

# Apply: weighted mean instead of median
composite <- weighted.mean(c(kenpom, bpi, torvik, evanmiya, teamrankings, powerrank),
                           w = source_rmse$weight, na.rm = TRUE)
```

**Important**: If a source is missing for a team (NA), redistribute its weight proportionally
to the remaining sources. Don't just drop it.

**Expected improvement**: Small but consistent. Probably 0.01-0.02 Brier score improvement.
The median is already 80% of the way there.

### D. Integrating Market Prices (The Biggest Edge)

**Current**: Model ignores betting markets entirely.

**Why this is the single biggest improvement**: The betting market is the most efficient
predictor in sports. Pinnacle's closing line has been shown to be more accurate than any
public model, including KenPom. The reason: it aggregates private information from thousands
of sophisticated bettors.

**But**: The market isn't perfect. It has biases:
- **Favorite-longshot bias**: Markets slightly overvalue favorites in futures
- **Public money bias**: Popular teams (Duke, Kansas) get inflated prices
- **Information lag**: Markets take hours to fully adjust to injury/lineup news

**Your edge comes from the GAP between your model and the market** — not from your model alone.

**Implementation approach (Bayesian blend)**:

```r
# Step 1: Convert Kalshi/Vegas prices to implied probabilities
market_prob <- devig(kalshi_yes_price, kalshi_no_price)

# Step 2: Your model's probability
model_prob <- sim_results %>%
  filter(team == "Duke") %>%
  summarise(prob = mean(Champion))

# Step 3: Blend with calibrated weight
# The weight should come from backtest: how much does adding your model
# improve predictions over market-only?
blend_weight <- 0.3  # Start here, calibrate from backtest
blended_prob <- blend_weight * model_prob + (1 - blend_weight) * market_prob

# Step 4: The edge is where blended_prob diverges from market_prob
edge <- blended_prob - market_prob
# If edge > fee_threshold → trade
```

**How to calibrate `blend_weight`**:
1. For each historical tournament, get pre-tournament market prices (Pinnacle or consensus)
2. Get your model's pre-tournament probabilities
3. Try blend_weight = 0.1, 0.2, 0.3, 0.4, 0.5
4. For each weight, compute Brier score of blended prediction vs actual outcome
5. Pick the weight that minimizes Brier score

**Expected finding**: blend_weight will probably be 0.2-0.4. If it's <0.1, your model
doesn't add much over the market. If it's >0.5, your model is genuinely better than
the market (unlikely but possible for derivative props).

**Where your model adds the MOST value**: Not on championship odds (market is very efficient
there) but on **derivative props** — seed sums, upset counts, highest seed to advance.
The market prices these lazily because they're low-volume. Your sim can compute them rigorously.

### E. Strength of Schedule / Conference Adjustment

**Current**: Raw power ratings, no adjustment.

**Why this matters**: Mid-major teams with inflated ratings (from beating weak opponents) tend
to underperform in the tournament. Power conference teams with "worse" ratings (from playing
strong opponents) tend to overperform.

**How to implement**:

```r
# Method 1: Conference prior (simple)
conf_strength <- ratings %>%
  group_by(conference) %>%
  summarise(conf_avg = mean(composite_rating))

adjusted_rating <- composite_rating +
  0.1 * (conf_strength$conf_avg[match(conference, conf_strength$conference)] - mean(conf_strength$conf_avg))

# Method 2: Bayesian shrinkage toward historical tournament performance
# Key insight: how do teams from this conference historically perform
# vs their seed expectations in the tournament?
# Big 12 teams historically OUTPERFORM their seed by +1.5 pts
# Mid-major teams historically UNDERPERFORM by -2.0 pts
# Apply this as an adjustment to composite_rating
```

**Validation**: Check if adjusted ratings predict tournament outcomes better than raw ratings
in backtests. If the improvement is <0.5% Brier score, skip it — not worth the complexity.

---

## Kalshi Edge Improvements

### F. Liquidity-Weighted Edge Scoring

**Current problem**: A +150% EV edge on a 1¢ market with $50 volume is worthless. You can't
execute it at size.

**Proposed**: Replace raw EV% ranking with expected profit:
```r
# Max you can trade at the displayed price ≈ min(your_size, ask_depth)
executable_size <- min(target_bet_size, volume * 0.05)  # 5% of volume as liquidity proxy

# Spread penalty: wide spreads mean you're crossing more to get filled
spread_penalty <- 1 / (1 + spread / 10)  # 1.0 for 0¢ spread, 0.5 for 10¢ spread

# Expected profit in dollars
expected_profit <- ev_pct * executable_size * spread_penalty

# Sort by this, not raw EV%
```

### G. Model Confidence Filter

Add a column to the Kalshi Edges table showing how confident we are in the fair value:

```r
confidence <- case_when(
  n_rating_sources >= 5 & abs(fair - 50) < 30 ~ "High",
  n_rating_sources >= 3 & abs(fair - 50) < 40 ~ "Medium",
  TRUE ~ "Low"
)
```

**Rationale**: Our model is most accurate near 50% (balanced matchups with lots of data).
At the extremes (95% or 5%), small errors in the rating get amplified into large probability
errors. A team rated +25 vs +23 might be 97% vs 95% — the 2-point rating error causes a
2% probability error. But for close matchups, the same 2-point error causes only a 5%
probability error.

### H. Track Per-Game Results in C++ for Accurate Props

**Current**: C++ sim only outputs which teams advance to each round. For upset counts, we
reconstruct matchup pairs from advancing teams — this is approximate and breaks for R32+.

**Proposed for 2027**: Modify C++ to output a per-game results matrix:
```cpp
// Additional output: game_results[sim * n_games + game_idx] = winner_seed | loser_seed
// Then in R: upset = (winner_seed > loser_seed)
```

This gives exact upset counts, exact seed matchup tracking, and enables new prop types
like "will a specific seed pair matchup occur?"

---

## Infrastructure

### I. Rating Source Reliability

Current fragility ranking (most to least fragile):
1. **EvanMiya** — Hard-coded pivot table parsing. Will break.
2. **TeamRankings** — HTML scraping. Will break if they redesign.
3. **PowerRank** — HTML scraping. Same risk.
4. **Torvik** — Google Sheet (manual paste) + Playwright fallback. Moderate risk.
5. **KenPom** — Google Sheet + cbbdata fallback. Moderate risk.
6. **BPI** — ESPN API. Most stable, least likely to break.

**For 2027**: Add graceful degradation. If a source fails, log a warning and continue with
remaining sources. The median of 4-5 sources is still robust. Don't crash the whole pipeline
because EvanMiya's sheet format changed.

### J. Rating Source Caching

Cache fetched ratings in DuckDB with a 1-hour TTL. This makes dashboard restarts instant
(~2 seconds instead of ~30 seconds for fresh fetches).

### K. Automate Data Entry

The manual Google Sheet workflow (copy KenPom data, paste into Sheet, refresh dashboard)
is error-prone and slow. For 2027:
- **KenPom**: Use `cbbdata::cbd_kenpom_ratings()` directly (no Sheet needed)
- **Torvik**: Use `fetch_torvik.py` Playwright scraper (already built, use as primary)
- **EvanMiya**: Scrape evanmiya.com directly or find API

---

## Execution Order for 2027

Do these in order. Each builds on the previous.

### Phase 1: Foundation (January, before tournament)
1. Build backtest framework (`backtest.R`)
2. Collect historical ratings data (2015-2026)
3. Run baseline backtest with current model
4. Measure: Brier score, log-loss, calibration curve, upset accuracy

### Phase 2: Model Tuning (February, use backtest to validate)
5. Test variable game variance — does seed-dependent SD improve Brier?
6. Test t-distribution — does it improve upset rates?
7. Calibrate composite weights — inverse variance weighting from historical RMSE
8. Run parameter sensitivity sweep (sd_margin, beta1, beta2)

### Phase 3: Market Integration (March, before tournament starts)
9. Collect pre-tournament Kalshi/Vegas prices
10. Implement Bayesian blend of model + market
11. Calibrate blend weight from backtest
12. Add liquidity scoring + confidence filter to Kalshi Edges tab

### Phase 4: Operations (Tournament time)
13. Update TEAM_NAME_FIXES for new D1 programs
14. Automate rating data fetch (reduce manual Google Sheet work)
15. Run dashboard, find edges, trade

---

## Known Bug to Fix

### EvanMiya parser is brittle
Hard-coded 21-column pivot. When the Google Sheet format changes, the parser silently fails.
**Fix**: Parse by column name, not position. Add validation: if Team column is numeric, skip
source gracefully.

---

## Key Metrics to Track Year-Over-Year

| Metric | 2026 Value | 2027 Target |
|--------|------------|-------------|
| Brier score (championship) | Unknown (no backtest) | < 0.15 |
| R64 upset rate accuracy | Unknown | Within 3% of historical |
| Calibration slope | Unknown | 0.9 - 1.1 (well-calibrated) |
| Kalshi edges traded | ~0 (display only) | 10+ with +EV |
| Rating sources available | 5-6 / 6 | 6 / 6 automated |
| Dashboard startup time | ~90 seconds | < 30 seconds (cached) |
