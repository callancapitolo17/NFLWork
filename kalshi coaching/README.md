# Kalshi NFL Coaching Odds Dashboard

Interactive dashboard for tracking NFL head coach hiring odds from Kalshi, with portfolio tracking and real-time refresh.

## Quick Start

```bash
./run.sh
```

This will:
1. Fetch latest odds from Kalshi API
2. Fetch your positions and orders (if authenticated)
3. Generate the interactive dashboard
4. Start local server at `http://127.0.0.1:8080`
5. Open the dashboard in your browser

Click the **Refresh Data** button anytime to re-fetch live data.

## Features

### Interactive Heatmap
- Candidates × Teams grid with probability colors
- Hover for details (bid/ask, liquidity, volume)
- Your positions marked with ★
- Zoom and pan support

### Portfolio Tracking
- **Open Positions** - All positions with exposure, P&L, entry price
- **Changes Since Last Run** - What you bought/sold with before/after
- **Resting Orders** - Pending orders with expiration times

### Modern UI
- Dark theme with glassmorphism effects
- Sidebar navigation for quick jumping
- Sortable, filterable, searchable tables
- Hover animations and smooth transitions

## Portfolio Setup (Optional)

To track your positions and orders:

1. Log in to [kalshi.com](https://kalshi.com) → **Settings** → **API Keys**
2. Create an API key and download the private key file
3. Copy `.env.example` to `.env`:

```bash
cp .env.example .env
```

4. Edit `.env` with your credentials:
```
KALSHI_API_KEY_ID=your-key-id
KALSHI_PRIVATE_KEY_PATH=/path/to/your/private-key.pem
```

5. Install cryptography:
```bash
pip3 install cryptography
```

## Usage Options

```bash
./run.sh            # Start server with refresh button (recommended)
./run.sh --static   # Just generate HTML, no server
```

## Files

| File | Description |
|------|-------------|
| `kalshi_coaching.py` | Fetches odds + portfolio from Kalshi API |
| `kalshi_coaching_display.R` | Generates interactive dashboard |
| `server.py` | Local server with /refresh endpoint |
| `kalshi_coaching.duckdb` | Historical data (odds, positions, orders) |
| `report.html` | Generated dashboard |
| `lib/` | JavaScript libraries for interactivity |
| `run.sh` | Main entry point |

## Requirements

**Python:**
- `duckdb`
- `flask`
- `cryptography` (for portfolio tracking)

**R:**
- `tidyverse`
- `duckdb`
- `reactable`
- `plotly`
- `htmltools`
- `htmlwidgets`

Install R packages:
```r
install.packages(c("reactable", "plotly", "htmltools", "htmlwidgets"))
```

## Data Stored

The DuckDB database tracks:
- **coaching_odds_v2** - Historical odds snapshots
- **positions** - Your position history
- **resting_orders** - Order history
- **market_info** - Cached market titles/subtitles

## API Notes

- Public API: No auth needed for odds
- Authenticated API: Required for positions/orders
- Rate limit: 20 req/sec (Basic tier)
