# Hoop88 Correlated Parlay Edge Finder

Find +EV opportunities by exploiting correlation mispricing between hoop88 and SGP odds from DraftKings/FanDuel/BetMGM.

## The Edge

Hoop88 prices parlays by multiplying individual leg odds (ignoring correlation). SGP engines account for correlation between spread and total, pricing correlated parlays tighter.

**Positive correlations to exploit:**
- **Underdog spread + Under** - Low-scoring games favor underdogs covering
- **Favorite spread + Over** - High-scoring games favor favorites covering

If hoop88 offers +260 on a parlay that SGP prices at +220, that's a +12.5% edge.

## Quick Start

```bash
cd hoop88_correlation
source venv/bin/activate

cd "/Users/callancapitolo/NFLWork/hoop88_correlation"
python find_edge.py --pikkit --leagues NFL --periods FG 

# Auto-fetch SGP odds from Pikkit (recommended)
python find_edge.py --pikkit

# With visible browsers for debugging
python find_edge.py --pikkit --visible
```

## Installation

```bash
cd hoop88_correlation
python -m venv venv
source venv/bin/activate
pip install playwright python-dotenv beautifulsoup4
playwright install chromium
```

Create `.env` with credentials:
```
HOOP88_URL=https://hoop88.com
HOOP88_USERNAME=your_username
HOOP88_PASSWORD=your_password
```

## Pikkit Setup (First Time)

Pikkit requires phone number + SMS login. Run once to save session:

```bash
python scraper_pikkit.py --visible
```

1. Browser opens to Pikkit
2. Login with your phone number
3. Complete SMS verification
4. Session saves automatically after 120 seconds

Re-run this command if the session expires.

## Usage

```bash
source venv/bin/activate

# Auto-fetch SGP odds from Pikkit (recommended)
python find_edge.py --pikkit

# With visible browsers
python find_edge.py --pikkit --visible

# Specific leagues and periods
python find_edge.py --pikkit --leagues NCAAF --periods FG 1H
python find_edge.py --pikkit --leagues NFL --periods FG 1H 1Q 2Q

# Custom bankroll and Kelly fraction
python find_edge.py --pikkit --bankroll 50000 --kelly-fraction 0.3

# Allow half-point line differences
python find_edge.py --pikkit --tolerance 0.5

# Use specific trusted books
python find_edge.py --pikkit --books draftkings fanduel novig

# Manual SGP entry (old behavior)
python find_edge.py

# Just view hoop88 odds without comparison
python find_edge.py --no-compare

# Fast mode (estimated odds, no bet slip interaction)
python find_edge.py --fast
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--pikkit` | Off | Auto-fetch SGP odds from Pikkit |
| `--leagues` | NFL NCAAF | Leagues to scrape |
| `--periods` | FG 1H 1Q | Periods: FG, 1H, 1Q, 2Q, 3Q, 4Q |
| `--no-compare` | Off | Skip SGP comparison (just show hoop88 odds) |
| `--fast` | Off | Use estimated odds (faster, less precise) |
| `--visible` | Off | Show browser windows |
| `--bankroll` | 26000 | Your bankroll for Kelly sizing |
| `--kelly-fraction` | 0.25 | Fraction of Kelly to use (25% = quarter Kelly) |
| `--books` | prophetx novig fanduel draftkings | Trusted books to filter SGP odds |
| `--multiplier` | 1.1 | Conservative multiplier for decimal odds |
| `--tolerance` | 0 | Line match tolerance (0 = exact match required) |

## Supported Markets

**Leagues:**
- NFL
- NCAAF (College Football)

**Periods:**
- FG (Full Game)
- 1H (1st Half)
- 1Q (1st Quarter)
- 2Q (2nd Quarter)
- 3Q (3rd Quarter)
- 4Q (4th Quarter)

## How It Works

1. **Scrape hoop88** - Logs in, navigates to each league/period, extracts spreads and totals via bet slip
2. **Scrape Pikkit** - Finds same games, builds parlays, extracts SGP odds from multiple books
3. **Filter books** - Excludes books with adjusted/altered legs (different lines than requested)
4. **Calculate edge** - `(hoop88 - Pikkit_adjusted) / Pikkit_adjusted × 100`
5. **Kelly sizing** - Recommends bet size based on edge and your bankroll

## Line Matching

The scraper ensures you're comparing apples to apples:

- **Exact line matching** (`--tolerance 0`, default): Only compares when Pikkit has the exact same spread and total as hoop88
- **Adjusted leg filtering**: Books that alter the line (e.g., changing total from 46 to 45.5) are automatically excluded
- **Half-point tolerance** (`--tolerance 0.5`): Allows comparing lines within half a point

When a book changes your requested line, Pikkit marks it with an orange dot. The scraper detects this via the `adjusted` field in the API and filters those books out.

## Example Output

```
======================================================================
CORRELATED PARLAY EDGE FINDER
======================================================================
Leagues: NFL
Periods: FG
SGP Source: Pikkit (auto)
Trusted Books: prophetx, novig, fanduel, draftkings
Multiplier: 1.1x
Line Tolerance: 0

Scraping hoop88 for ACTUAL parlay odds (interacting with bet slip)...
Login successful

============================================================
Scraping NFL FG
============================================================
Found 1 games

Processing: Seahawks @ Patriots
  Patriots +4.5 + Under 46.0: +260
  Seahawks -4.5 + Over 46.0: +260

======================================================================
FETCHING SGP ODDS FROM PIKKIT
======================================================================
Pikkit logged in!

============================================================
Seahawks @ Patriots
============================================================

  NFL FG
  ----------------------------------------

    Patriots +4.5 + Under 46.0
    Hoop88: +260
    Pikkit SGP: +232 (adj) [avg:+211 | Fanduel:+215, Draftkings:+210, Novig:+200, ProphetX:+220]
    --> EDGE: +8.2%  |  Kelly: $650

    Seahawks -4.5 + Over 46.0
    Hoop88: +260
    Pikkit SGP: +242 (adj) [avg:+220 | Fanduel:+225, Draftkings:+218, Novig:+210, ProphetX:+228]
    --> EDGE: +5.1%  |  Kelly: $380

======================================================================
RECOMMENDED BETS - SORTED BY EDGE
======================================================================
Bankroll: $26,000 | Kelly: 25%

  +8.2% EDGE -> Bet $650
    Seahawks @ Patriots (FG)
    Patriots +4.5 + Under 46.0
    Hoop88: +260 vs Pikkit: +232

  +5.1% EDGE -> Bet $380
    Seahawks @ Patriots (FG)
    Seahawks -4.5 + Over 46.0
    Hoop88: +260 vs Pikkit: +242

----------------------------------------------------------------------
  TOTAL RECOMMENDED ACTION: $1,030 (4.0% of bankroll)
```

**Note:** Books with adjusted legs (different lines than requested) are automatically filtered out. The output shows only books offering the exact lines you requested.

## Kelly Criterion

The tool uses Pikkit SGP odds as the "true" probability (since they account for correlation) and calculates optimal bet size:

```
Kelly % = (bp - q) / b
where:
  b = hoop88 decimal odds - 1
  p = true win probability (from SGP)
  q = 1 - p
```

Default is 25% Kelly (quarter Kelly) for safety. Full Kelly is too aggressive for most bankrolls.

## Files

| File | Purpose |
|------|---------|
| `find_edge.py` | Main tool - scrapes hoop88 + Pikkit, calculates edge and Kelly |
| `scraper_hoop88_odds.py` | Hoop88 scraping with bet slip interaction |
| `scraper_pikkit.py` | Pikkit SGP scraper (gets odds from multiple books) |
| `.env` | Credentials |

## Troubleshooting

**Pikkit session expired:**
```bash
python scraper_pikkit.py --visible
# Login manually, wait for session to save
```

**No games found on Pikkit:**
- Check if the sport filter is correct (Pikkit uses "NCAAFB" not "NCAAF")
- Verify games are listed on the Events page

**Hoop88 login failed:**
- Check credentials in `.env`
- Run with `--visible` to debug
