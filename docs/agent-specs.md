# Agent Specifications

Autonomous subagents that handle complex multi-step workflows with minimal human oversight.

---

## 1. Daily Pipeline Agent

**Problem**: The CBB pipeline (`run.py`) is partially automated but frequently breaks. Manual intervention required for token refreshes, missing data, timeout retries. User runs the pipeline and babysits it ~daily.

**What it does autonomously**:
1. Check auth tokens for all configured books (auto-refresh if possible)
2. Run `Answer Keys/run.py --sport cbb`
3. If a scraper fails:
   - Token expired → attempt refresh → retry once
   - Timeout → retry with longer timeout
   - Site structure changed → skip book, alert user
4. If R script fails:
   - Missing data → check which scraper failed, report
   - Team name mismatch → log the unmatched names
5. Verify dashboard updated (check `report.html` mtime)
6. Report summary: books scraped, games found, markets found, any issues

**Escalation**: Only alert the user for novel failures (site redesign, new error types). Common failures are handled automatically.

**Schedule**: Run via LaunchAgent/cron at configured time daily (during season).

**Implementation**: Claude Code task or a Python script with retry logic that wraps `run.py`.

---

## 2. Backtest Agent

**Problem**: Running backtests requires manually configuring parameters, running R scripts, and interpreting results. Parameter sweeps across date ranges, sports, and market types are tedious.

**What it does autonomously**:
1. Accept a backtest specification (sport, date range, markets, parameters to sweep)
2. Generate all parameter combinations
3. Run each backtest variation (parallelize where possible)
4. Collect results into a summary table
5. Identify the best-performing configuration
6. Generate a report with:
   - ROI, CLV, hit rate per configuration
   - Statistical significance (confidence intervals)
   - Drawdown analysis
   - Edge decay over time
   - Recommendation: which config to deploy

**Guard rails**:
- Reject results with < 500 bets as insufficient sample
- Flag potential overfitting if best config is an outlier
- Always report out-of-sample performance (train/test split)

**Implementation**: Agent that orchestrates multiple R script executions and aggregates results.

---

## 3. Market Expansion Agent

**Problem**: Adding new market types (team totals, alt lines, 1H, player props) requires changes in 4+ places: scraper config, R model, correlation matrix, dashboard filters. This has been done repeatedly and follows a pattern.

**What it does autonomously**:
1. Accept the new market type definition
2. Check which scrapers already provide this market type
3. Update scraper parsing to extract the new market
4. Add market to R answer key processing
5. Define correlation relationships with existing markets
6. Add to dashboard filter options
7. Run pipeline to verify end-to-end
8. Report: what changed, what to review

**Guard rails**:
- Always use a worktree
- Require user review before merge
- Validate that new market doesn't break existing markets

**Implementation**: Claude Code agent working in a worktree, following the existing patterns in each layer.

---

## 4. Bet Auto-Placer Agent (Existing Worktree)

**Problem**: Manually placing bets on soft books is slow. By the time you log in and navigate, the line may have moved.

**What it does autonomously**:
1. Monitor dashboard for new +EV opportunities above threshold
2. Apply Kelly sizing with correlation adjustments
3. Check current exposure to avoid overconcentration
4. Place bet via Playwright automation on the target book
5. Log the bet to Google Sheets
6. Verify the bet was placed at the expected odds (or better)

**Guard rails**:
- **Human-in-the-loop confirmation** for bets above configurable threshold
- Maximum single bet size cap
- Maximum daily exposure cap
- Never place on a game already at max exposure
- Pause if consecutive placement failures (site may be blocking)

**Status**: Worktree already exists at `.claude/worktrees/bet-auto-placer/`. Continue from existing work.

---

## Implementation Priority

| Agent | Impact | Complexity | Priority |
|-------|--------|------------|----------|
| Daily Pipeline | HIGH — eliminates daily babysitting | MEDIUM | 1 |
| Backtest | MEDIUM — saves hours per analysis | MEDIUM | 2 |
| Market Expansion | MEDIUM — saves ~2 sessions per market | HIGH | 3 |
| Bet Auto-Placer | HIGH — captures fleeting edge | HIGH | 4 |
