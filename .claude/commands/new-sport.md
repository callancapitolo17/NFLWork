# New Sport Expansion

Expand the answer key system to support: $ARGUMENTS

## Process

### 1. Research
- Check what data is available from The Odds API for this sport (`SPORT_KEYS` in `run.py`)
- Identify which existing scrapers already support this sport (check `SCRAPER_CONFIGS` sports arrays)
- Determine what market types are available (spreads, totals, moneylines, 1H, team totals)

### 2. Create Answer Key
Use an existing answer key as the template:
- **For team sports**: Use `CBB Answer Key/CBB.R` as the guide (merged single-script pattern)
- **For player props**: Use `props.R` as the guide

Create `Answer Keys/<Sport> Answer Key/<SPORT>.R` with:
- Odds API data acquisition
- Offshore odds integration from DuckDB
- Sample generation (use Rcpp sampler from `src/`)
- Fair odds calculation via devigging
- EV computation against each book
- Kelly sizing with correlation adjustments
- HTML dashboard output via reactable

### 3. Create Dashboard
Create `Answer Keys/<Sport> Dashboard/`:
- `<sport>_dashboard.R` — reactable HTML generation
- `<sport>_dashboard_server.py` — Flask server (pick unused port)
  - Use `cbb_dashboard_server.py` as template
  - Include: bet placement, exposure tracking, book settings, filter settings

### 4. Wire into Pipeline
- Add R script config to `R_SCRIPTS` in `run.py`
- Add sport key to `SPORT_KEYS` if not present
- Enable relevant scrapers for this sport in `SCRAPER_CONFIGS`

### 5. Backtest
- Create `<Sport> Answer Key/<SPORT>_Backtest.R`
- Run against historical data to validate edge exists
- Require 1000+ game sample before trusting results
- Report: ROI, CLV, hit rate, by market type

### 6. Validate
- Run full pipeline: `cd "Answer Keys" && python run.py --sport <key>`
- Verify dashboard loads and shows odds from all books
- Confirm team name matching works
- Test bet placement flow

## Rules
- Use a worktree for this work
- Use `Tools.R` shared utilities — don't duplicate devig/Kelly functions
- Follow the same dark theme and layout as existing dashboards
- Port assignments: NFL=8081, CBB=8082, new sport=next available
