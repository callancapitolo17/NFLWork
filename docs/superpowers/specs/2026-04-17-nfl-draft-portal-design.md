# NFL Draft EV Portal — Design

**Date:** 2026-04-17
**Author:** Callan + Claude (brainstorming)
**Status:** Draft, awaiting user review
**Target ship date:** 2026-04-23 (Day 1 of 2026 NFL Draft)

---

## Goal

A trader's-cockpit web portal for surfacing +EV bets on 2026 NFL Draft markets across Kalshi and four sportsbooks (DraftKings, FanDuel, Bookmaker, Wagerzon). The user grinds mock-draft sites and other research themselves; the portal's job is to (a) put every venue's prices side-by-side so disagreement is obvious, (b) stream Kalshi trade flow so large bets are visible as they print, and (c) capture the bets the user takes for later analysis.

This is **Phase 1** of a multi-phase build. Phase 2 work (formal fair-value methodology, P&L view, market maker) is explicitly out of scope.

---

## Non-Goals (Phase 1)

- Formal fair-value / "true probability" oracle. v1 surfaces every venue's devigged probability and **flags outlier books that are ≥10pp from the cross-book median**; the user eyeballs edges. A methodology will be designed in Phase 2 after watching live behavior during the 2026 draft.
- Auto-betting or Kelly sizing recommendations.
- WebSocket / push streaming. Everything is poll-based in v1.
- Full P&L tab, CLV analysis, or backtesting view. v1 only *captures* bets; analysis is Phase 2.
- Kalshi market making for draft markets.
- Refactoring or extending non-draft scrapers (DK/FD/BM/WZ general-odds scrapers stay as-is).
- Refactoring `kalshi_draft/fetcher.py`'s discovery / API / portfolio code beyond the minimum needed to repoint its writes from `kalshi_draft.duckdb` to `nfl_draft/nfl_draft.duckdb`.

---

## Architecture

### Storage

A single DuckDB database at `nfl_draft/nfl_draft.duckdb` — the **only** draft database for the project. The existing `kalshi_draft/kalshi_draft.duckdb` is retired in v1: its data (a few months of Kalshi draft odds + portfolio history) is migrated into the new DB at the start of implementation, then the old file is deleted. Two DBs for the same logical concern (NFL draft data) is confusing; one is correct.

**Migration of existing data**: a Python script (`nfl_draft/lib/migrate_from_kalshi_draft.py`) copies tables from `kalshi_draft/kalshi_draft.duckdb` into `nfl_draft/nfl_draft.duckdb` with renames where needed:
- `draft_odds` (Kalshi-only) → `kalshi_odds_history` (preserved for back-reference; new portal writes go to the new `draft_odds` table described below)
- `draft_series`, `market_info` → kept as-is
- `positions`, `resting_orders` → kept as-is

The script is **idempotent**: it checks whether each destination table already exists with row count matching the source; if yes, it skips that table; if no, it copies. Re-running is safe and cheap (the existence check is O(1)). Tested via `test_migration.py` which runs the script twice in succession and asserts identical destination state.

After migration is verified, `kalshi_draft/kalshi_draft.duckdb` is removed and three things are updated together in the same commit as the migration script:

1. **Existing dashboard tabs** in `kalshi_draft/app.py` — point `duckdb.connect(...)` at `nfl_draft/nfl_draft.duckdb` AND rewrite every SQL query referencing `draft_odds` to reference `kalshi_odds_history` instead.
2. **Existing fetcher / portfolio writers** in `kalshi_draft/fetcher.py` (and any related modules) — repoint every `duckdb.connect(...)` write target from `kalshi_draft.duckdb` to `nfl_draft/nfl_draft.duckdb`. The fetcher's writes to `draft_series`, `market_info`, `positions`, `resting_orders` continue functioning, just to the new DB. This is a connection-target change, not a logic change.
3. **`kalshi_mm/`** — review for any direct DuckDB references to `kalshi_draft.duckdb` (the running CBB MM bot uses `kalshi_draft/auth.py` but shouldn't touch the draft DB; verify with grep). If found, repoint similarly.

Regression covered by:
- `test_dashboard_queries.py` — every existing tab returns expected row count from `kalshi_odds_history`.
- `test_kalshi_writer_target.py` — invokes the existing fetcher in a sandbox, confirms it writes to `nfl_draft.duckdb`, NOT `kalshi_draft.duckdb`.

**Migration verification**: idempotency check uses `MD5(GROUP_CONCAT(...))` hash of each table's content (not just row count), since two different datasets can have the same row count. The hash is computed inside DuckDB via aggregation; cheap and reliable.

**No local CSV / JSON config files** (per project rule "no temp files"). All static configuration — player aliases, team aliases, per-book market label mappings — lives as Python dict literals in `nfl_draft/config/*.py` modules (mirrors the existing `wagerzon_odds/config.py` and `bookmaker_odds/scraper.py` pattern). At runtime these dicts are loaded into in-memory lookups; for cross-process consistency they are also seeded into five DuckDB lookup tables (`players`, `player_aliases`, `teams`, `team_aliases`, `market_map`) by `python -m nfl_draft.lib.seed`. **Seed runs at the start of every `run.py` invocation** (idempotent, ~50ms) so any mid-week edit to a `config/*.py` file takes effect on the next cron tick automatically — no separate manual step.

**Modeling convention**: each `market_id` represents a single **binary outcome** (matches Kalshi's YES/NO contract model natively). Sportsbook outright menus (e.g., DK's "First QB Drafted" with multiple player options) are split into one binary contract per option at ingestion. So "Will Cam Ward be the first QB drafted" is `first_qb_ward`, "Will Shedeur Sanders" is `first_qb_sanders`, etc. This keeps the schema flat and the join-across-venues logic simple — every book's offering of a given outcome lands on the same `market_id`.

**Scraper module shape**: each `nfl_draft/scrapers/<book>.py` exposes `fetch_draft_odds() -> List[OddsRow]`. The Kalshi module additionally exposes `fetch_trades() -> List[TradeRow]` for the trade-tape feed.

**Time zones**: all timestamps are stored as **local time** (the user's machine timezone, presumably ET). DuckDB `TIMESTAMP` (not `TIMESTAMPTZ`). Rationale: the user reads the dashboard in local time, the cron fires in local time, and the draft itself is a single time-zoned event. Mixing UTC and local in storage creates silent drift bugs that are painful during draft week.

**TZ boundary discipline** (encoded in `lib/db.py` helpers):
- `local_to_utc_iso(local_ts) -> str` — converts a stored local TIMESTAMP to UTC ISO-8601 with offset, used when calling Kalshi's API (which expects ISO-8601). Example: `min_ts` parameter in trades poll.
- `utc_iso_to_local(iso_str) -> datetime` — converts Kalshi's ISO-8601 response timestamps to naïve local TIMESTAMP for insertion.
- All Kalshi I/O routes through these two functions; no raw datetime arithmetic at the API boundary. Tested via `test_tz_roundtrip.py` (round-trip stability for known values, DST transitions).
- Other books (DK/FD/BM/WZ) return local-time-context endpoints already; no conversion needed beyond storing what they return.

**Credentials & secrets**: each sportsbook scraper reads its login credentials from `.env` at the repo root (already in `.gitignore`). Required vars per book: `DK_USERNAME` / `DK_PASSWORD`, `FD_USERNAME` / `FD_PASSWORD`, `BOOKMAKER_USERNAME` / `BOOKMAKER_PASSWORD`, `WAGERZON_USERNAME` / `WAGERZON_PASSWORD`. Kalshi reuses the existing API-key pair from `kalshi_draft/.env` (no new secrets). Per-book session cookies persist to `nfl_draft/.cookies/<book>.json` (in `.gitignore`). On session expiry, the scraper attempts re-login once; if that fails, the run is logged as auth-failed and surfaces in the dashboard footer (per Error handling).

**Tables:**

- `draft_markets` — the canonical market catalog. **Each row is a single binary outcome** (matches Kalshi's YES/NO contract model — see Modeling convention above).
  - `market_id` TEXT PRIMARY KEY — built by deterministic construction rule (see below)
  - `market_type` TEXT (`first_at_position` | `pick_outright` | `top_n_range` | `team_first_pick` | `prop`)
  - `subject_player` TEXT (canonical player name; null when not applicable)
  - `subject_team` TEXT (canonical team abbr; null when not applicable)
  - `position` TEXT (`QB` | `RB` | `WR` | `TE` | `OL` | `DL` | `LB` | `DB` | `K` | null)
  - `pick_number` INT (for `pick_outright` markets, e.g., 1, 2, 3...; null otherwise)
  - `range_low`, `range_high` INT (for `top_n_range` markets; null otherwise)
  - `created_at` TIMESTAMP

  **`market_id` construction rule** (must be applied identically in every scraper adapter — encoded once in `lib/market_map.py` as `build_market_id(market_type, **kwargs)`):
  - `first_at_position`: `f"first_{position.lower()}_{slug(player)}"` → `first_qb_cam-ward`
  - `pick_outright`: `f"pick_{pick_number}_overall_{slug(player)}"` → `pick_1_overall_cam-ward`
  - `top_n_range`: `f"top_{range_high}_{slug(player)}"` → `top_5_ashton-jeanty`
  - `team_first_pick`: `f"team_{team.lower()}_first_pick_{slug(player)}"` → `team_chi_first_pick_ashton-jeanty`
  - `prop`: `f"prop_{slug(short_description)}"` → `prop_will-a-kicker-go-r1`

  Where `slug(name)` = `name.lower().replace(' ', '-').replace('.', '').replace("'", '')`. The single canonical builder eliminates the risk of two scrapers producing different `market_id`s for the same outcome. Tested by `test_market_id_construction.py`.

- `draft_odds` — append-only odds snapshots
  - `market_id` TEXT FK → `draft_markets.market_id`
  - `book` TEXT (`kalshi` | `draftkings` | `fanduel` | `bookmaker` | `wagerzon`)
  - `american_odds` INT (e.g., -110, +250)
  - `implied_prob` DOUBLE (raw, before devig)
  - `devig_prob` DOUBLE (after vig removal across complementary outcomes)
  - `fetched_at` TIMESTAMP
  - INDEX on (market_id, book, fetched_at DESC)

- `kalshi_trades` — Kalshi public trade tape
  - `trade_id` TEXT PRIMARY KEY (Kalshi's unique trade event ID; enables `INSERT OR IGNORE` for dedup)
  - `ticker` TEXT
  - `side` TEXT (`yes` | `no`)
  - `price_cents` INT (1-99)
  - `count` INT (number of contracts traded)
  - `notional_usd` DOUBLE (computed at insert: count × price_cents × 0.01)
  - `traded_at` TIMESTAMP
  - `fetched_at` TIMESTAMP
  - INDEX on (ticker, traded_at DESC)

  Note: `is_large` is **not stored** — it's computed at read time in the dashboard query (`notional_usd ≥ threshold`) so the threshold can be tuned without backfill. Threshold configurable via dashboard slider, default $500.

- `kalshi_poll_state` — incremental polling cursor for the trade tape
  - `series_ticker` TEXT PRIMARY KEY (e.g., `KXNFLDRAFT1`)
  - `last_traded_at` TIMESTAMP (the `traded_at` of the most recent trade we've ingested for this series)
  - `last_polled_at` TIMESTAMP

  Each `--mode trades` cycle reads `last_traded_at` per series, requests trades from Kalshi with `min_ts = last_traded_at`, paginates by cursor, INSERTs (with `INSERT OR IGNORE` for safety), updates `last_traded_at`. Avoids re-fetching the entire trade history every cycle.

- `draft_bets` — light bet logging
  - `bet_id` UUID PRIMARY KEY
  - `market_id` TEXT FK
  - `book` TEXT
  - `side` TEXT
  - `american_odds` INT
  - `stake_usd` DOUBLE
  - `taken_at` TIMESTAMP
  - `note` TEXT (optional free-form, e.g., "after Schefter Ward news")

- `players` — canonical name registry, seeded from `nfl_draft/config/players.py` (Python dict literal)
  - `canonical_name` TEXT PRIMARY KEY
  - `position` TEXT
  - `college` TEXT (for disambiguation when two prospects share a name)

- `player_aliases` — alias → canonical lookup, seeded from the same `players.py` (one row per alias)
  - `alias` TEXT PRIMARY KEY (case-folded at seed time, e.g., `'cam ward'`)
  - `canonical_name` TEXT FK → `players.canonical_name`

  Why a separate table: alias resolution is on the hot path (every scraped row). DuckDB doesn't index inside `TEXT[]` list columns, so `WHERE 'Cam Ward' = ANY(aliases)` is O(N) per lookup. A flat `player_aliases` table with a PK on `alias` is O(1). Same conceptual model, just normalized for performance.

- `teams` — canonical NFL team registry, seeded from `nfl_draft/config/teams.py` (Python dict literal). Built off the existing canonical mapping in `Answer Keys/canonical_match.py` — do not duplicate.

- `team_aliases` — alias → canonical lookup, same pattern as `player_aliases`.

- `market_map` — per-book label → canonical `market_id` mapping, seeded from `nfl_draft/config/markets.py` (Python dict literal)
  - `book` TEXT (`kalshi` | `draftkings` | `fanduel` | `bookmaker` | `wagerzon`)
  - `book_label` TEXT (the raw market title or ticker prefix from the book)
  - `book_subject` TEXT (raw player/team string from the book, before normalization)
  - `market_id` TEXT FK → `draft_markets.market_id`
  - PRIMARY KEY (book, book_label, book_subject)

The five lookup tables (`players`, `player_aliases`, `teams`, `team_aliases`, `market_map`) are populated by `python -m nfl_draft.lib.seed`, which reads the dict literals from `nfl_draft/config/*.py` and writes to DuckDB. The seed is **idempotent** (truncate-and-rewrite each lookup table inside a single transaction) and runs automatically at the start of every `run.py` invocation, so editing a config dict mid-week takes effect on the next cron tick. The dict literals are version-controlled (source code), so changes show up in `git diff` cleanly. No CSV or JSON config files anywhere.

Manual maintenance is unavoidable for rookie-name canonicalization — any auto-fuzzy-match will produce silent join errors that silently corrupt EV calculations. Better to fail loudly: scraper rows with unmapped players land in a quarantine table (see Error handling) and the dashboard footer surfaces the unmapped count, which is the user's daily prompt to extend the dict.

### Code layout

```
nfl_draft/
  __init__.py                  # makes nfl_draft an importable package
  config/
    __init__.py
    players.py                 # PLAYERS dict literal (canonical → metadata + aliases)
    teams.py                   # TEAMS dict literal (canonical → aliases)
    markets.py                 # MARKET_MAP dict literal (per-book label → canonical market_id)
  scrapers/
    __init__.py
    kalshi.py                  # extends existing kalshi_draft/fetcher.py for ALL series
    draftkings.py              # adapted from mlb_sgp/scraper_draftkings_sgp.py
    fanduel.py                 # adapted from mlb_sgp/scraper_fanduel_sgp.py
    bookmaker.py               # adapted from bookmaker_odds/scraper.py
    wagerzon.py                # adapted from wagerzon_odds/scraper_v2.py
  lib/
    __init__.py
    db.py                      # DuckDB connection + schema migrations
    devig.py                   # ported from Answer Keys/Tools.R (~30 LOC)
    normalize.py               # reads player_aliases/team_aliases → alias resolution
    market_map.py              # reads market_map table → per-book label → canonical market_id
    seed.py                    # reads config/*.py dicts, writes to DuckDB lookup tables (idempotent)
    migrate_from_kalshi_draft.py  # migrates kalshi_draft.duckdb → nfl_draft.duckdb (idempotent)
  tests/                       # see Testing section
    __init__.py
    ...
  run.py                       # orchestrator: --mode {pre-draft, draft-day, trades} --book all
  README.md
  .cookies/                    # gitignored — per-book session cookies
  nfl_draft.duckdb             # gitignored — the SOLE draft database
```

**Directory naming wart, accepted in Phase 1**: the `kalshi_draft/` directory now houses a *multi-venue* draft dashboard, not just Kalshi. Renaming would break `kalshi_mm/` which imports `kalshi_draft/auth.py` (4 hard references in `kalshi_mm/{orders,analyze_performance,config,dashboard}.py`) — the running CBB MM bot would go down. Phase 2 cleanup: rename `kalshi_draft/` → `kalshi_lib/` (it's really the Kalshi auth+fetcher library that `kalshi_mm` depends on), update kalshi_mm imports, and move `app.py` into `nfl_draft/`.

### Devig math

Port three functions from `Answer Keys/Tools.R` to `nfl_draft/lib/devig.py`:

- `american_to_implied(american_odds) -> float`
- `devig_two_way(odds_a, odds_b) -> (prob_a, prob_b)`
- `devig_n_way(odds_list) -> List[float]` (for first-QB / first-RB style multi-outcome markets where the complementary set is the full position pool)

For markets where the complementary set isn't fully posted (e.g., book only offers Yes on 5 of the 30+ first-QB candidates), use **proportional devig**: divide each implied prob by the sum of all posted implieds for that market group on that book.

**Known limitation**: proportional devig assumes the posted candidates contain ~100% of the probability mass. When a book only posts the top 5 of a 30-candidate field, the unposted 25 still carry positive probability — so the 5 posted, when normalized to sum to 1.0, are systematically inflated. Comparing these inflated devigged probs against Kalshi (which prices the full field) gives a directional bias: book candidates look "too high" relative to Kalshi by a roughly constant offset. **The Cross-Book Grid annotates each book's devigged prob with the count of posted candidates** (e.g., "DK: 38% (5/30)") so the user can mentally discount accordingly. A more accurate approach (Bayesian field-completion using the full Kalshi field as a prior) is queued for Phase 2.

### Scraper orchestration

`nfl_draft/run.py` is the single entry point. Two modes — they have **functionally different work**, not just different cadences:

```
python run.py --mode scrape --book all          # all books → draft_odds
python run.py --mode scrape --book kalshi       # single book, for debugging
python run.py --mode trades                     # Kalshi trade tape only → kalshi_trades
```

- Every invocation runs `lib/seed.py` first (idempotent, ~50ms) so dict edits take effect immediately.
- `--mode scrape`: each scraper module exposes `fetch_draft_odds() -> List[OddsRow]`. `run.py` calls them sequentially (parallelism deferred — scrape errors easier to debug serially), normalizes each result through `lib/normalize.py` and `lib/market_map.py`, devigs via `lib/devig.py`, writes to `draft_odds`.
- `--mode trades`: polls `/markets/trades` for each Kalshi series, requesting trades after `kalshi_poll_state.last_traded_at` for that series, paginating by cursor, and appending to `kalshi_trades` with `INSERT OR IGNORE` on `trade_id`. Updates `kalshi_poll_state.last_traded_at` to the max `traded_at` ingested per series.

Cadence is controlled by **cron**, not by mode flags.

**Cron** (macOS `cron` invoked under `caffeinate` for the draft window so the laptop doesn't sleep). Each entry must use absolute paths (cron has empty `PATH` / `PYTHONPATH`):

```
# Pre-draft (default — installed at v1 ship)
*/15 * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode scrape --book all >> nfl_draft/logs/scrape.log 2>&1
*/2  * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode trades            >> nfl_draft/logs/trades.log 2>&1

# Draft-day (swap manually morning of Apr 23 — see crontab.draft)
*/2 * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode scrape --book all >> nfl_draft/logs/scrape.log 2>&1
*    * * * * cd /Users/callancapitolo/NFLWork && /usr/bin/env python -m nfl_draft.run --mode trades            >> nfl_draft/logs/trades.log 2>&1
```

Two crontab files (`nfl_draft/crontab.pre`, `nfl_draft/crontab.draft`) shipped in the repo. Swap with `crontab nfl_draft/crontab.draft`. Stale-data thresholds in the dashboard footer (30 min pre-draft = 2 missed cycles; 5 min draft-day = 2.5 cycles) match these cadences.

`nfl_draft/logs/` is gitignored; rotation handled by `logrotate` config in README.

### Kalshi expansion

Today `kalshi_draft/fetcher.py` discovers via `/series` keyword filter ("nfl"+"draft") plus a hardcoded prefix fallback. v1 changes:

- Add all observed `KXNFLDRAFT*` prefixes to the fallback list (incl. `KXNFLDRAFTQB`, `KXNFLDRAFTWR`, `KXNFLDRAFTRB`, `KXNFLDRAFTPOS`, `KXNFLDRAFTTOP5`, `KXNFLDRAFTTOP10`, `KXNFLDRAFTPICK`, `KXNFLDRAFTTOP`).
- After discovery, log every series found that has open markets but no downstream consumer in the dashboard — these are markets we're storing but not yet rendering. Used to catch new Kalshi series (e.g., a `KXNFLDRAFTTRADE` for "first trade" if Kalshi adds one) without code changes.
- Map each Kalshi market to a canonical `market_id` via the `MARKET_MAP` dict in `nfl_draft/config/markets.py` (Kalshi-side keyed by ticker prefix + market subject extracted from the market title), seeded into the `market_map` table.

**Rate limiting**: Kalshi public endpoints allow ~100 req/min, authenticated ~10 req/min. First-run discovery (8+ series × paginated markets ≈ 100+ requests) sits near the public limit. The Kalshi scraper module enforces:
- 600ms minimum sleep between requests (caps at ~100 req/min)
- Exponential backoff on HTTP 429 (1s, 2s, 4s, 8s, 16s, then fail-loud)
- Single shared rate-limit token bucket across odds + trades subcommands so they don't compound when run in parallel

These are baked into a `kalshi_request()` wrapper in `nfl_draft/scrapers/kalshi.py`; no caller bypasses it. Implementation step 1 for the Kalshi adapter: read `kalshi_draft/auth.py`'s `public_request` — if it already throttles, the new wrapper just delegates to it; if not, add throttling at the existing function so both `kalshi_mm/` and the new portal benefit. Don't fork the throttle logic.

**Cross-process coordination**: cron runs `--mode scrape` and `--mode trades` as separate processes that both call Kalshi. A Python-local token bucket doesn't coordinate across processes. v1 accepts this: each process can hit ~50 req/min (half the throttle), total ≤ 100 req/min worst case — under the public limit. If observed 429s exceed 1/hour in pre-draft, escalate to a file-lock-based shared bucket (file at `nfl_draft/.kalshi_throttle.lock` with last-N-request timestamps; deferred to Phase 2).

### Per-book scraper notes

- **Kalshi**: existing `kalshi_draft/fetcher.py` is the base; `nfl_draft/scrapers/kalshi.py` is a thin adapter that calls it and remaps to the new schema.
- **DraftKings**: reuse the CDP click-and-capture + Akamai handling pattern from `mlb_sgp/scraper_draftkings_sgp.py`. Reconnaissance step required: locate the NFL Draft event page, capture the API call(s) that fetch the draft markets menu.
- **FanDuel**: reuse the curl_cffi + 3-header recipe (per memory `fanduel_sgp_scraping.md`). Reconnaissance: same as DK.
- **Bookmaker**: reuse the curl_cffi + ASP.NET cookie pattern from `bookmaker_odds/scraper.py`. Reconnaissance: discover the league ID for NFL Draft markets (currently only cbb/nba/mlb configured).
- **Wagerzon**: reuse the ASP.NET `__VIEWSTATE` POST pattern from `wagerzon_odds/scraper_v2.py`. Reconnaissance: same — find draft league IDs.

For each book, reconnaissance and parsing is the unknown — auth and request mechanics are solved.

### Dashboard tabs (extending `kalshi_draft/app.py`)

New tabs added to the existing Dash app:

1. **Cross-Book Grid** — Markets on rows, books on columns, devigged prob in cells. For each market with ≥ 2 books posting, compute the **median devigged probability** across all posting venues. For each (market, book) cell, compute `delta = book_prob - median_prob`. **Flag the cell** (colored highlight + inline `±Npp` indicator) when `abs(delta) ≥ threshold`. Threshold default: 10 percentage points; user-configurable via a slider in the tab header. The flagged cell is the "this book is way off the market" signal — that's where +EV likely lives.
   Final column shows the count of flagged books for that market (so a market with 3 outliers stands out from one with 1).
   Click a row to drill into a per-market detail view (line history if data exists, all bid/ask, take/log button).

2. **+EV Candidates** — Flat-listed view of every flagged (market, book) pair from the grid above, ranked by `abs(delta)` descending. One row per outlier — so a single market with 3 flagged books contributes 3 rows. Columns: market, position, book, book_prob, market_median, delta (signed: + means book is too high vs market, − means too low), implied edge direction. Filter by market type, position, book, threshold. This is the user's "scan the field for action" view.

   Why outlier-flag and not raw variance: when 4 books say 50% and one says 30%, variance and outlier-flag give the same signal (one book is off). But when 5 books each disagree by ~5pp, variance is high but no single book is the obvious target — outlier-flag (correctly) flags nothing. The flag is more directly actionable than the variance number.
3. **Trade Tape** — Streaming-style table of recent Kalshi trades (last 200 rows), large fills highlighted (`notional_usd ≥ threshold` evaluated at query time; threshold defaults to $500, configurable via tab slider). Filter by ticker / series.
4. **Bet Log** — Form to log a bet (market_id searchable-dropdown via `dcc.Dropdown(searchable=True)`, book dropdown, american_odds, stake, note). Below the form, table of past bets. The market_id dropdown is search-as-you-type because hard-scrolling 300+ markets is unusable.
   **Pre-fill from context**: rows in the +EV Candidates tab include a "Log this bet" button that navigates to the Bet Log tab with `market_id`, `book`, and `american_odds` pre-filled (passed via `dcc.Store(storage_type='session')`). This keeps the act-on-an-edge flow to two clicks.

**Tab organization** — adding 4 new tabs to the existing 5 makes 9 total. Tabs are grouped into two top-level sections via `dcc.Tabs` with sub-tabs (or two horizontal nav rows):
- **Portal** (new — the cockpit): Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log
- **Kalshi-only legacy** (existing): Market Overview, Price History, Edge Detection, Consensus, Portfolio

Default landing tab: Cross-Book Grid (the workhorse of the trading session). Sub-tab selection persisted in `dcc.Store(storage_type='local')`.

### Dashboard auto-refresh

All four new tabs (Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log) use a `dcc.Interval` component to auto-refresh without user F5. Cadence is mode-driven (matches the cron schedule):

- **Pre-draft mode**: 60 seconds (slightly faster than the 15-min scrape so the user sees fresh writes within a minute of arrival).
- **Draft-day mode**: 15 seconds (matches the 1-2 min scrape; user sees new prices nearly immediately).

A mode toggle in the dashboard header (`Pre-draft / Draft-day`) controls the interval. Default: pre-draft. The user manually flips to draft-day the morning of April 23 (same time they swap the cron file).

**Cheap-poll guard** (avoids redundant work when nothing has changed): the auto-refresh callback first queries `SELECT MAX(fetched_at) FROM draft_odds`, compares against the value cached in a `dcc.Store`, and skips the heavy grid re-render if unchanged. This means polling every 15s against a 90s scrape costs one cheap MAX query per tick (≈ 5ms), with the expensive grid-render firing only when there's new data. Same pattern for Trade Tape against `MAX(fetched_at) FROM kalshi_trades`. Eliminates the "14 of 15 polls return identical data" waste. Tested by `test_dashboard_queries.py::test_cheap_poll_guard_skips_render_when_unchanged`.

Each callback re-queries DuckDB using a short-lived read-only connection (per the concurrency pattern in Error Handling), so the refresh cost is bounded and lock-free.

**`dcc.Store` registry** (single source of truth, prevents id collisions):

| `id` | scope | type | purpose |
|---|---|---|---|
| `nfl_draft.mode_toggle` | local | string | "pre-draft" or "draft-day" mode |
| `nfl_draft.subnav` | local | string | "portal" or "kalshi-legacy" top-level section |
| `nfl_draft.last_fetched_odds` | local | ISO timestamp | cheap-poll guard for draft_odds |
| `nfl_draft.last_fetched_trades` | local | ISO timestamp | cheap-poll guard for kalshi_trades |
| `nfl_draft.bet_log_prefill` | session | object `{market_id, book, american_odds}` | "Log this bet" button payload from +EV Candidates |

All store ids prefixed `nfl_draft.` to avoid collision with the existing `kalshi_draft/app.py` Store ids. Type stays `local` for state that should survive reload; `session` for ephemeral cross-tab handoffs.

---

## Data flow

```
[ DK | FD | BM | WZ | Kalshi ]                      Cron @ 15min (pre) / 2min (draft)
        │
        ▼
   scrapers/<book>.py  →  raw odds rows
        │
        ▼
   lib/normalize.py       (resolve player/team aliases → canonical)
        │
        ▼
   lib/market_map.py      (per-book label → canonical market_id)
        │
        ▼
   lib/devig.py           (american → implied → devig per market group)
        │
        ▼
   nfl_draft.duckdb       draft_markets, draft_odds

[ Kalshi /markets/trades ]                          Cron @ 2min (pre) / 1min (draft)
        │
        ▼
   scrapers/kalshi.py (trades subcmd)
        │
        ▼
   kalshi_trades         (INSERT OR IGNORE on trade_id PK; is_large computed at read time)
   kalshi_poll_state     (per-series last_traded_at cursor for incremental polls)

[ Dash app — kalshi_draft/app.py extended ]
        │
        ├── Cross-Book Grid ←  reads draft_markets + latest draft_odds per (market, book); computes median + outlier flags
        ├── +EV Candidates  ←  same query, flat + filtered to flagged outliers, sorted by |delta|
        ├── Trade Tape      ←  reads kalshi_trades ORDER BY traded_at DESC LIMIT 200
        ├── Bet Log         ←  writes draft_bets; reads it for table view
        └── (existing tabs unchanged)
```

---

## Error handling

- **Per-scraper isolation**: a failure in one book's scraper must not prevent the others from running. `run.py` wraps each book in try/except, logs the error, continues.
- **Auth failure**: log and surface in dashboard footer (e.g., "DK auth expired — re-login needed"). Don't silently lose data.
- **New / unknown market**: when a scraper emits a market label that doesn't map to a canonical `market_id`, log a warning and write the row to a `draft_odds_unmapped` quarantine table. Daily review during draft week to extend `nfl_draft/config/markets.py` (then re-run `python -m nfl_draft.lib.seed`). Better than silently dropping.
- **Player alias miss**: same pattern — quarantine to `draft_odds_unmapped_players`, daily review to extend `nfl_draft/config/players.py` (then re-run `seed`). Surface unmapped count on the dashboard footer.
- **Stale data warning**: dashboard shows `fetched_at` per book. If any book is > 30 min stale during pre-draft mode (or > 5 min during draft mode), highlight in red.
- **DuckDB lock conflicts**: DuckDB (0.10+) allows multiple readers + one writer cross-process when WAL is enabled (default). Two design rules ensure clean concurrency without locks becoming user-visible:

  1. **Short-lived write connections**. `run.py` opens the DB, performs the bulk INSERT using DuckDB's `appender` API (fastest bulk-insert path), and closes the connection immediately. Total lock duration: milliseconds per cycle, even with several hundred rows.
  2. **Short-lived read connections**. Each Dash callback opens a fresh `read_only=True` connection, executes its query, returns the result, closes. No long-lived connection reuse across callbacks. Reads typically complete in <100ms.

  Both sides use `with duckdb.connect(...) as con:` to guarantee close-on-exit. If a read happens to land mid-write, it briefly blocks (typically <50ms) — invisible to the user.

  All connection management is encapsulated in `nfl_draft/lib/db.py` via two helpers:
  ```python
  @contextmanager
  def write_connection() -> Iterator[duckdb.DuckDBPyConnection]: ...

  @contextmanager
  def read_connection() -> Iterator[duckdb.DuckDBPyConnection]: ...
  ```
  Every caller uses these. No raw `duckdb.connect()` outside `db.py`. This is enforced via a unit test that greps the codebase for `duckdb.connect` outside of `lib/db.py` and fails if any unauthorized matches are found. **Allowlist**: `nfl_draft/lib/migrate_from_kalshi_draft.py` is permitted to call `duckdb.connect()` directly because it must open both the legacy and new DBs simultaneously (which the context-manager pattern doesn't support cleanly). The lint test reads its allowlist from a constant at the top of the test file; any future additions require an explicit code change.

  **Stress test**: an integration test (`test_concurrent_access.py`) spawns a writer thread doing INSERTs in a loop for 10 seconds while a reader thread does SELECTs in a loop. Asserts zero unhandled lock errors and that all writes are visible to the reader within 1 second. Catches regressions to the connection-management discipline.

  **Fallback** if contention nonetheless emerges under load (>1 lock error per minute observed): writes go to append-to-JSONL by `run.py`, with a separate batch worker draining JSONL → DuckDB every 30s. Adds 30s of staleness but eliminates contention entirely. Not implemented in v1; reserved for Phase 2.

---

## Testing

Tests are a first-class deliverable, not an afterthought. Two reasons:

1. **Rollout confidence**: when the dashboard is doing real-money calculations on draft day, "the math worked yesterday" isn't enough. Every commit must demonstrably preserve correctness.
2. **Auto-mode feedback loop**: when implementation runs autonomously (Claude iterating without per-step user review), tests are the *only* signal that a change worked. Without them, regressions land silently and compound. With them, each run produces a binary pass/fail with a precise failure location.

### Test organization

```
nfl_draft/tests/
  unit/                              # deterministic, fast (<1s total)
    test_devig.py                    # american→implied→devig math, all edge cases
    test_normalize.py                # player/team alias resolution
    test_market_map.py               # per-book label → canonical market_id
    test_db.py                       # schema migration idempotency, seed correctness
    test_scraper_parsing.py          # raw response → OddsRow conversion (pure logic)
  integration/                       # uses fixtures + temp DuckDB, no network (<30s total)
    test_pipeline.py                 # scraper rows → normalize → market_map → devig → DB
    test_migration.py                # kalshi_draft.duckdb → nfl_draft.duckdb roundtrip + idempotency
    test_dashboard_queries.py        # query functions return correct shape + values
    test_quarantine.py               # unmapped player/market lands in quarantine table
    test_concurrent_access.py        # writer + reader threads, asserts no lock errors (per Error handling)
    test_lint_db_connect.py          # greps for duckdb.connect() outside lib/db.py (allowlist enforced)
  live/                              # hits real APIs, run manually only (~2-5 min)
    test_scrapers_live.py            # each scraper produces ≥1 row from live endpoint
    test_dashboard_renders.py        # Dash app subprocess + GET each tab URL
  fixtures/
    kalshi/                          # captured API responses (snapshot once, replay forever)
    draftkings/
    fanduel/
    bookmaker/
    wagerzon/
  conftest.py                        # shared pytest fixtures: temp DB, seeded lookup tables
```

### What each tier covers

**Unit (`tests/unit/`)** — pure functions, no I/O, no network, no DB. Run after every code change.
- `test_devig.py`: every odds-conversion edge case — favorite/dog american odds, overround = 0% / 5% / 30%, two-way and n-way devig, proportional devig with sparse outcomes, integer/float input handling.
- `test_normalize.py`: alias resolution case-insensitive, whitespace-tolerant; known player + every alias → canonical; unknown player → `None`; conflicting aliases (same alias for two players) → raises.
- `test_market_map.py`: every per-book label → market_id; unknown label → `None`; collision detection (two book_labels mapping to same market_id is allowed; one book_label mapping to two market_ids raises).
- `test_db.py`: schema migration runs twice without error (idempotent); seed runs from PYTHON dict literal and produces expected rows; seed re-run is idempotent; lookup tables truncated-and-rewritten cleanly.
- `test_scraper_parsing.py`: each scraper's response-parsing function (separated from network calls) takes a fixture JSON/HTML and produces expected `OddsRow`s. This isolates parsing bugs from auth bugs.

**Integration (`tests/integration/`)** — uses a temp DuckDB and fixture data, no network. Run before every commit.
- `test_pipeline.py`: feed each book's fixture through the full pipeline; verify (a) correct row count in `draft_odds`, (b) devig sums to ~1.0 per market group, (c) quarantine count matches expected unmapped rows.
- `test_migration.py`: build a fake `kalshi_draft.duckdb` with known seed data, run migration script, verify all rows present in destination + correct table renames + idempotent on re-run.
- `test_dashboard_queries.py`: seed temp DB with known data, call each tab's query function (broken out from Dash callbacks), verify result shape matches dashboard expectations. **This is the test that catches "existing tabs break after DB repoint" — every existing tab gets a query test before its query is touched.** Specific cases include:
  - **Outlier-flag math**: feed 5 books' devigged probs, assert median is correct, assert exactly the books with `|delta| ≥ threshold` are flagged, assert delta sign is correct.
  - **Cheap-poll guard**: simulate two consecutive callback fires with no new data; assert the heavy render runs once, the second call short-circuits.
  - **Variance NULL handling**: market with only 1 book → variance is `None`/N/A, doesn't break sort.
- `test_quarantine.py`: feed scraper rows with deliberately unmapped player and unmapped market; verify they land in `draft_odds_unmapped` and `draft_odds_unmapped_players` and do NOT land in `draft_odds`.

**Live (`tests/live/`)** — hits real endpoints. Run manually before merge to main, and as the morning-of-April-22 go/no-go gate.
- `test_scrapers_live.py`: each scraper produces ≥ 1 row from the actual book; auth works; rate limits not exceeded.
- `test_dashboard_renders.py`: spin up the Dash app in a subprocess, GET each tab URL, assert 200 + key DOM elements present.

### Fixtures

For each book scraper, capture one full real response during reconnaissance (curl the endpoint, save to `tests/fixtures/<book>/<scenario>.json`). These fixtures are committed to git (no PII; just market data). Re-capture if the book changes its response shape.

Fixture rule: scrapers should be split into `fetch_raw()` (network, hard to test) + `parse(raw_response)` (pure, easy to test). Tests target `parse()`. Smoke tests target `fetch_raw()`.

### Test cadence

- **Per code change**: `pytest tests/unit` (~1s). No excuse to skip.
- **Per commit**: `pytest tests/unit tests/integration` (~30s). Pre-commit hook recommended (just runs pytest).
- **Per merge to main**: full `pytest` including live, plus the morning-of dashboard render check. Failure here = block merge.
- **Quarantine watch**: each scraper run logs unmapped counts; if counts grow during the day, that's a signal that the dict needs extending. Tests pin known-good aliases so dict edits don't accidentally break previously-working ones.

### What's NOT tested

- Visual rendering / CSS / Plotly chart appearance — manual eyeball.
- Actual cron scheduling — manual verification that the cron entries fire.
- Bet placement (out of scope entirely in Phase 1).

### Auto-mode contract

When implementation runs in auto mode, the loop is:

1. Make the change.
2. Run `pytest tests/unit tests/integration`.
3. If green: commit + move on. If red: fix, GOTO 2.

This requires that **every change be testable**. If a change can't be expressed as a test, it requires manual review before commit. Acceptable cases for manual review: dashboard CSS tweaks, README edits. Non-acceptable: anything that touches `lib/`, `scrapers/`, or query functions.

---

## Phasing

### Phase 1 — Ship by Apr 23 (this week)

In scope:
- DuckDB schema, all tables (`draft_markets`, `draft_odds`, `kalshi_trades`, `kalshi_poll_state`, `draft_bets`, `players`, `player_aliases`, `teams`, `team_aliases`, `market_map`, `draft_odds_unmapped`, `draft_odds_unmapped_players`)
- Migration of existing `kalshi_draft.duckdb` data into `nfl_draft.duckdb`; deletion of the old DB; existing dashboard tabs repointed AND queries rewritten to reference `kalshi_odds_history` (with regression query tests landing first)
- Player/team/market canonical config in `nfl_draft/config/*.py` (Python dict literals, populated for the top ~80 prospects), seeded into DuckDB lookup tables (idempotent seed)
- All 5 venue scrapers (Kalshi expansion + 4 books) writing to `draft_odds`, with rate limiting on Kalshi and credentials loaded from `.env`
- Kalshi trade tape polling → `kalshi_trades` with `kalshi_poll_state` cursor + `INSERT OR IGNORE` dedup
- Devig math ported to Python (with proportional-devig limitation annotated in the dashboard)
- Dashboard tabs grouped Portal vs Kalshi-only legacy: Cross-Book Grid (with outlier flags), +EV Candidates (with "Log this bet" button), Trade Tape, Bet Log (searchable dropdown)
- Auto-refresh with cheap-poll guard (MAX(fetched_at) check before re-render)
- Cron setup: pre-draft + draft-day swap
- All timestamps in local time (TIMESTAMP, not TIMESTAMPTZ)
- **Comprehensive test suite** (unit + integration + live; see Testing section). Tests are not optional — every `lib/`, `scrapers/`, and dashboard query function ships with tests. Pre-commit hook runs unit + integration on every commit.
- `requirements.txt` updates as needed
- README in `nfl_draft/`

Quality-over-speed posture: the deadline (Apr 23) is firm but speed is not the constraint. If a tradeoff appears between "ship faster with shortcuts" vs "ship a day later with full tests + correct error handling," the latter wins. Auto-mode iteration depends on tests; rolling out untested code on draft day is the failure mode to avoid.

### Phase 2 — Post-draft, weeks following

Out of scope for v1, planned for after:
- Formal fair-value oracle (sharp-weighted consensus, possibly with mock-draft prior)
- Real-time Kalshi WebSocket trade tape (replace polling)
- Line-movement charts per market
- P&L view + CLV analysis on `draft_bets`
- Kalshi MM adapted from `kalshi_mm/`
- Generalize the framework for NFL Draft 2027, NBA Draft, NHL Draft (rename `nfl_draft/` → `draft/`, parameterize sport)
- **Rename `kalshi_draft/` → `kalshi_lib/`** and move `app.py` into `nfl_draft/`. Update `kalshi_mm/` imports. Eliminates the directory-naming wart from Phase 1.

### Phase 3 — Out

- Auto-betting, Kelly sizing recommendations, alerting

---

## Version control

- **Branch**: `feature/nfl-draft-portal` (already created at brainstorm time)
- **Files created**: `nfl_draft/` directory tree (see Code layout); `docs/superpowers/specs/2026-04-17-nfl-draft-portal-design.md`
- **Files modified**: `kalshi_draft/app.py` (new tabs + queries rewritten from `draft_odds` → `kalshi_odds_history` + repointed at `nfl_draft/nfl_draft.duckdb`), `.gitignore` (add `nfl_draft/nfl_draft.duckdb`, `nfl_draft/.cookies/`), `README.md` (link to nfl_draft/README.md), `requirements.txt` if present (add any new pinned versions of `dash`, `plotly`, `duckdb`, `pytest`, `curl_cffi`, `playwright` not already declared)
- **Files deleted**: `kalshi_draft/kalshi_draft.duckdb` (after migration verified). Note: this file is already gitignored, so the deletion is filesystem-only — no commit needed for the removal itself.
- **Commit structure**: per-phase commits, each ships with its own tests in the same commit — (1) schema + migration script + migration tests, (2) seed + lookup tables + seed/normalize tests, (3) devig math + devig tests, (4) each scraper + that scraper's parsing tests + a fixture, (5) dashboard query functions + query tests + repoint existing tabs (with their pre-existing-tab query tests landing first as a regression baseline), (6) trade-tape poller + tests, (7) cron config + README. Review checkpoint at each scraper commit. No commit lands without tests for the code in it.
- **Worktree**: optional but recommended given other active work on `feature/cbb-consensus-calibration` and `feature/mlb-sgp-scrapers`. Use `/worktree` for isolation if any of those are likely to need attention during the 7 days.
- **Pre-merge review**: full executive review of `git diff main..HEAD` before merging — focus on (a) no DK/FD credentials in code, (b) DuckDB connections all closed, (c) no logging of personal bets to stdout, (d) cookie files gitignored.
- **Approval to merge**: explicit ask before `git merge` or any push.

## Documentation

Required updates in the same final merge:
- `nfl_draft/README.md` — setup (env vars, auth steps per book, cron install), usage (`run.py` flags), troubleshooting (auth failures, schema migrations).
- Top-level `README.md` — one-line link to `nfl_draft/README.md`.
- `CLAUDE.md` (project root) — add `nfl_draft/` to the project structure section so future Claude sessions discover it.
- `kalshi_draft/README.md` — note that the dashboard has been extended with draft-portal tabs and **repointed at `nfl_draft/nfl_draft.duckdb`**; the old `kalshi_draft.duckdb` has been retired (data migrated, file deleted).

---

## Open questions deferred to Phase 2

- How should fair-value be computed? (Sharp-weighted consensus? Mock-draft Bayesian prior? Hybrid with user override?)
- Is Kalshi a market or an oracle? (Affects EV direction — bet against Kalshi when books move, or bet against books when Kalshi prints heavy?)
- Threshold for "large bet" on Kalshi tape — default $500 notional in v1; tune from real data after the draft.

---

## Success criteria for Phase 1

By end of day April 22 (the day before draft Day 1):

1. All 5 venues populating `draft_odds` for at least the #1-overall and first-QB markets, refreshed in the last 15 minutes
2. Cross-Book Grid renders all populated markets with devigged probs per book, median across books, and outlier flags on cells ≥10pp from median
3. Kalshi Trade Tape shows the last hour of trades, large fills highlighted
4. Bet Log accepts and persists a test bet
5. Cron jobs running on the user's machine; draft-day swap rehearsed
6. **All tests green**: `pytest tests/unit tests/integration` returns 0 failures; `pytest tests/live` returns 0 failures on the morning-of run; existing dashboard tabs (Market Overview, Price History, Edge Detection, Consensus, Portfolio) verified to render correct data after the DB repoint via `test_dashboard_queries.py`.
