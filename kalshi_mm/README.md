# Kalshi CBB 1H Market Maker

Model-driven market maker for Kalshi's CBB 1st-half spread, total, and moneyline markets. Uses conditional Kelly criterion for position sizing.

## How It Works

1. **Pricing oracle**: CBB answer key (R) generates sample-based fair value probabilities for every 1H market
2. **Kelly sizing**: Computes optimal bet size per market based on edge magnitude, bankroll, and covariance with existing positions
3. **Quoting**: Posts resting YES/NO orders at 5% EV floor with orderbook-aware penny logic
4. **Taking**: Crosses the spread on +EV opportunities (5% after 7% taker fee)
5. **Monitoring**: Runs Bookmaker + Bet105 scrapers every 60s to detect line moves

## Setup

```bash
# 1. Copy env template and add your Kalshi credentials
cp .env.example .env
# Edit .env with your KALSHI_API_KEY_ID and KALSHI_PRIVATE_KEY_PATH

# 2. Install dependencies
pip install cryptography duckdb numpy

# 3. Run the CBB pipeline first to generate predictions + game samples
cd "../Answer Keys" && python run.py --sport cbb

# 4. Run the market maker
python main.py              # Live mode
python main.py --dry-run    # Compute quotes without placing orders
```

## Kelly Criterion Sizing

Position sizing is driven by the conditional Kelly criterion. No hard position limits — Kelly is the sole sizing authority.

### How it works

1. **R pipeline exports game simulation samples** (home_margin_h1, total_h1) to `cbb_game_samples` table in `cbb.duckdb`
2. **kelly.py loads samples** and evaluates binary outcomes for each bet type (spreads, totals, ML) across all simulations
3. **Covariance matrix** computed from simulation outcomes — captures correlation between bets on the same game (e.g., home spread and home ML are highly correlated; home spread and total are nearly independent)
4. **Conditional Kelly** adjusts new bet sizes based on existing positions: `f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)`
5. **Fallback**: If covariance matrix is ill-conditioned, scales single-bet Kelly by `1/sqrt(1 + (n-1)*avg_ρ)`

### Key parameters

| Setting | Default | Description |
|---------|---------|-------------|
| `BANKROLL` | 1000.0 | Total bankroll in dollars |
| `KELLY_FRACTION` | 0.25 | Fractional Kelly (0.25 = quarter Kelly) |
| `USE_KELLY_SIZING` | true | Set to false to revert to fixed CONTRACT_SIZE |

### Per-side sizing

Buying YES at 40c has different edge (and therefore different Kelly size) than selling YES at 60c. The bot computes separate `bid_size` and `ask_size` for every market.

### Maker vs taker sizing

The taker uses `kelly_size_for_take()` which differs from maker sizing in two ways:
1. **Execution price**: Uses the actual ask (YES) or 100-bid (NO) price, not the resting bid/ask
2. **Fee-adjusted**: Adds the 7% taker fee to the effective price before computing Kelly, producing smaller sizes that reflect true post-fee edge

### Conditional adjustment

When you already hold positions on a game, Kelly accounts for correlation:
- **Long home spread + quoting home ML**: Kelly reduces ML size (highly correlated — both pay when home wins big)
- **Long home spread + quoting total**: Kelly doesn't reduce size (nearly independent)
- **Long home spread + quoting away spread at different line**: Kelly may reduce or increase depending on the correlation between strike points

### No samples = no quotes

If game simulation samples aren't available (pipeline hasn't run, R export failed), Kelly returns size 0 and the bot doesn't quote that market. There is no fallback to fixed sizing — this is intentional. Trading without the covariance data is flying blind.

## Risk Philosophy

**Kelly is the only risk control for position sizing.** All hard position limits (per-market, per-event, per-game, directional, exposure) have been removed. The 5% EV floor prevents -EV bets; Kelly prevents oversizing.

What's still enforced:
- **5% minimum EV** on all maker quotes and taker fills (after fees)
- **Staleness check** — pulls all quotes if predictions are >10 min old
- **Tipoff pullback** — pulls quotes 30 min before game start
- **Line move detection** — monitors Bookmaker/Bet105 for sharp moves
- **Fair value range** — won't quote if fair prob is <10c or >90c (model unreliable at extremes)
- **Anti-penny-loop** — stops chasing if counterparty is walking our bid up or ask down
- **Minimum spread** — 4c minimum spread enforced

### Why no position limits?

The old regime used fixed 5-contract sizing with hard caps (5/ticker, 8/event, 40 directional, $250 exposure). These were belt-and-suspenders because fixed sizing has no concept of edge magnitude or correlation. Kelly replaces all of this:
- **Per-market limits** → Kelly sizes to 0 when there's no edge
- **Per-event limits** → Conditional Kelly reduces size for correlated positions
- **Directional limits** → Kelly reduces size when existing positions are directionally correlated
- **Exposure limits** → Quarter Kelly naturally limits exposure relative to bankroll

### SKEW_PER_CONTRACT (disabled)

The quoter has a `SKEW_PER_CONTRACT` parameter (default: 0) that shifts bid/ask quotes based on inventory to attract offsetting flow. This is a standard market-making technique for **pure market makers who have no directional view** — their goal is to capture the spread and stay flat.

This bot is different: it has a model that produces fair values, and it *wants* directional positions when the model says there's edge. Kelly already handles inventory management through conditional adjustment (reducing new bet sizes when correlated positions exist). Skew would fight Kelly — Kelly says "buy more, the edge is huge" while skew says "quote worse because you already hold some." Set to 0 so they don't contradict.

## Configuration

All settings in `.env` (see `.env.example`):

| Setting | Default | Description |
|---------|---------|-------------|
| `BANKROLL` | 1000.0 | Total bankroll in dollars |
| `KELLY_FRACTION` | 0.25 | Fractional Kelly multiplier |
| `USE_KELLY_SIZING` | true | Toggle Kelly on/off |
| `MIN_EV_PCT` | 0.05 | Minimum EV% to post a maker quote |
| `MIN_TAKE_EV_PCT` | 0.05 | Minimum EV% to take (after 7% fee) |
| `CONTRACT_SIZE` | 5 | Fallback size when Kelly is off |
| `TAKE_CONTRACT_SIZE` | 5 | Fallback taker size when Kelly is off |
| `SKEW_PER_CONTRACT` | 0 | Inventory skew per contract (disabled) |
| `MAX_MARKETS` | 50 | Max tickers to quote (API rate guard) |
| `MAX_EVENTS` | 30 | Max games to quote (API rate guard) |
| `MAX_STALENESS_SEC` | 600 | Pull quotes if predictions older than this |
| `LINE_MOVE_THRESHOLD` | 0.5 | Pull quotes if offshore line moves by this much |
| `ENABLED_MARKETS` | spreads,totals,moneyline | Market types to quote |
| `PIPELINE_REFRESH_SEC` | 600 | Auto-refresh predictions interval |

## Architecture

```
CBB.R → cbb_raw_predictions (DuckDB) → main.py → Kalshi API
      → cbb_game_samples (DuckDB) → kelly.py ↗
                                        ↑
                      Bookmaker/Bet105 scrapers (line-move detection)
```

### Files

| File | Purpose |
|------|---------|
| `main.py` | Orchestrator — quote cycle, fill polling, pipeline refresh |
| `kelly.py` | Conditional Kelly sizing — loads samples, computes covariance, sizes bets |
| `quoter.py` | Price computation — EV floor, penny logic, anti-penny-loop |
| `taker.py` | Crosses the spread on +EV opportunities |
| `risk.py` | Staleness, tipoff proximity, line move detection |
| `config.py` | All settings loaded from .env |
| `db.py` | DuckDB state — positions, fills, resting orders, quote log |
| `orders.py` | Kalshi API wrapper — place, amend, cancel, auth |

### Data

All state stored in `kalshi_mm.duckdb`:
- `positions` — net position per market (includes `fair_prob`, `line_value`, `contract_team` for Kelly)
- `fills` — fill audit trail
- `resting_orders` — our current orders on Kalshi
- `quote_log` — every quote decision (24h retention)
- `sessions` — session summaries
- `reference_lines` — line snapshots for move detection
- `take_log` — taker attempt log
