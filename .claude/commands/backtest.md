# Run Backtest

Run a backtesting analysis. Parameters: $ARGUMENTS

Expected format: `<sport> [date_range] [market_types] [parameters]`
Examples:
- `/backtest cbb 2024-11-01 to 2025-03-01`
- `/backtest cbb --sweep-kelly 0.1,0.25,0.5`
- `/backtest nfl spreads,totals`

## Steps

### 1. Identify the right backtest script
- CBB: `CBB Answer Key/CBB_Backtest.R` or `CBB_Backtest_v2.R`
- CBB with correlation: `CBB_Backtest_Correlation.R`
- CBB parameter sweep: `CBB_Parameter_Sweep.R`
- NFL: `NFL Answer Key/Spreads_Totals_Backtest.R` or `All_Quarters_Backtest.R`

### 2. Configure parameters
- Date range (default: full available history)
- Market types to include
- Kelly fraction to test
- Minimum EV threshold
- Books to include/exclude

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
