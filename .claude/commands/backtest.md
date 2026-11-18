# Run Backtest

Run a backtesting analysis. Parameters: $ARGUMENTS

Expected format: `<sport> [date_range] [market_types] [parameters]`
Examples:
- `/backtest cbb 2024-11-01 to 2025-03-01`
- `/backtest cbb --sample-pct 0.01,0.02,0.03` (CBB: one labelled run per value, then one report)
- `/backtest nfl spreads,totals`

## Steps

### 1. Identify the right backtest script
- CBB: `CBB Answer Key/CBB_Backtest.R` (`run --label X`, then `report --labels a,b`; sample-size sweeps are one run per `--sample-pct`, single-book runs use `--books`; see `CBB Answer Key/README.md` § Backtest)
- CBB with correlation: `CBB_Backtest_Correlation.R`
- NFL: `NFL Answer Key/Spreads_Totals_Backtest.R` or `All_Quarters_Backtest.R`

### 2. Configure parameters
- CBB: `--sample-pct`, `--edge` (min EV), `--books`, `--max-games`, `--code-dir` (a worktree, to compare commits). It reports flat-bet ROI, not Kelly; it has no date-range, CLV, drawdown or by-book option, so say so rather than inventing those numbers.
- NFL scripts: date range, market types, Kelly fraction, minimum EV, books, as each script exposes them

### 3. Run the backtest
Execute the R script with configured parameters.

### 4. Analyze results
Report:
- **Overall ROI** — net profit / total wagered
- **CLV** — average closing line value
- **Hit rate** — wins / total bets
- **By market type** — ROI breakdown per market
- **By book** — which books provide the most edge
- **Sample size** — total bets (need 1000+ to trust)
- **Drawdown** — worst peak-to-trough
- **Edge decay** — is ROI decreasing over time?

### 5. Statistical significance
- Is ROI > 0 at 95% confidence?
- Is CLV consistently positive?
- Does edge survive transaction costs (vig)?

## Rules
- Never trust results with < 500 bets
- Always report confidence intervals, not just point estimates
- Compare against a "bet everything" baseline to confirm edge is real
