# Kalshi CBB 1H Market Maker

Low-stakes prototype for providing liquidity on Kalshi's CBB 1st-half spread markets.

## How It Works

1. **Pricing oracle**: The CBB answer key generates sample-based fair value probabilities for every 1H spread
2. **Quoting**: Posts resting YES/NO orders at `fair_value ± half_spread` on Kalshi
3. **Monitoring**: Runs Bookmaker + Bet105 scrapers every 60s to detect line moves; pulls quotes if lines shift >1 point
4. **Risk**: Position limits, exposure limits, staleness checks, and automatic kill switch

## Setup

```bash
# 1. Copy env template and add your Kalshi credentials
cp .env.example .env
# Edit .env with your KALSHI_API_KEY_ID and KALSHI_PRIVATE_KEY_PATH

# 2. Install dependencies (same as kalshi_draft)
pip install cryptography duckdb

# 3. Run the CBB pipeline first to generate predictions
cd "../Answer Keys" && python run.py --sport cbb

# 4. Run the market maker
python main.py              # Live mode
python main.py --dry-run    # Compute quotes without placing orders
```

## Configuration

All settings in `.env` (see `.env.example`):

| Setting | Default | Description |
|---------|---------|-------------|
| `HALF_SPREAD_CENTS` | 5 | Half-spread for quoting (total spread = 2x) |
| `CONTRACT_SIZE` | 2 | Contracts per side per market |
| `MAX_POSITION_PER_MARKET` | 5 | Max net contracts per ticker |
| `MAX_TOTAL_EXPOSURE` | 50.0 | Max dollars at risk |
| `MAX_MARKETS` | 10 | Max simultaneous markets to quote |
| `MAX_STALENESS_SEC` | 600 | Pull quotes if predictions older than this |
| `LINE_MOVE_THRESHOLD` | 1.0 | Pull quotes if offshore line moves by this much |

## Architecture

```
Answer Key (R) → cbb_raw_predictions (DuckDB) → Market Maker (Python) → Kalshi API
                                                        ↑
                                    Bookmaker/Bet105 scrapers (line-move detection)
```

## Data

All state stored in `kalshi_mm.duckdb`:
- `positions` — net position per market
- `fills` — fill audit trail
- `resting_orders` — our current orders
- `quote_log` — every quote decision for analysis
- `sessions` — session summaries
- `reference_lines` — line snapshots for move detection
