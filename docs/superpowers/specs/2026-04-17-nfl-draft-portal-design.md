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

- Formal fair-value / "true probability" oracle. v1 surfaces every venue's devigged probability and a variance metric; the user eyeballs edges. A methodology will be designed in Phase 2 after watching live behavior during the 2026 draft.
- Auto-betting or Kelly sizing recommendations.
- WebSocket / push streaming. Everything is poll-based in v1.
- Full P&L tab, CLV analysis, or backtesting view. v1 only *captures* bets; analysis is Phase 2.
- Kalshi market making for draft markets.
- Refactoring or extending non-draft scrapers (DK/FD/BM/WZ general-odds scrapers stay as-is).

---

## Architecture

### Storage

A new DuckDB database at `nfl_draft/nfl_draft.duckdb`. Kept separate from `kalshi_draft/kalshi_draft.duckdb` so this build can iterate (schema changes, wipes) without disturbing the existing Kalshi-only history or the Kalshi portfolio data.

**Modeling convention**: each `market_id` represents a single **binary outcome** (matches Kalshi's YES/NO contract model natively). Sportsbook outright menus (e.g., DK's "First QB Drafted" with multiple player options) are split into one binary contract per option at ingestion. So "Will Cam Ward be the first QB drafted" is `first_qb_ward`, "Will Shedeur Sanders" is `first_qb_sanders`, etc. This keeps the schema flat and the join-across-venues logic simple — every book's offering of a given outcome lands on the same `market_id`.

**Scraper module shape**: each `nfl_draft/scrapers/<book>.py` exposes `fetch_draft_odds() -> List[OddsRow]`. The Kalshi module additionally exposes `fetch_trades() -> List[TradeRow]` for the trade-tape feed.

**Tables:**

- `draft_markets` — the canonical market catalog
  - `market_id` TEXT PRIMARY KEY (e.g., `first_qb`, `pick_1_overall`, `team_chi_first_pick`, `top_5_jeanty`)
  - `market_type` TEXT (`first_at_position` | `pick_outright` | `top_n_range` | `team_first_pick` | `prop`)
  - `subject` TEXT (player canonical name, team canonical abbr, or null for type-level markets)
  - `position` TEXT (`QB` | `RB` | `WR` | `TE` | `OL` | `DL` | `LB` | `DB` | `K` | null)
  - `pick_number` INT (for `pick_outright` markets, e.g., 1, 2, 3...; null otherwise)
  - `range_low`, `range_high` INT (for `top_n_range` markets; null otherwise)
  - `created_at` TIMESTAMP

- `draft_odds` — append-only odds snapshots
  - `market_id` TEXT FK → `draft_markets.market_id`
  - `book` TEXT (`kalshi` | `draftkings` | `fanduel` | `bookmaker` | `wagerzon`)
  - `american_odds` INT (e.g., -110, +250)
  - `implied_prob` DOUBLE (raw, before devig)
  - `devig_prob` DOUBLE (after vig removal across complementary outcomes)
  - `fetched_at` TIMESTAMP
  - INDEX on (market_id, book, fetched_at DESC)

- `kalshi_trades` — Kalshi public trade tape
  - `ticker` TEXT
  - `side` TEXT (`yes` | `no`)
  - `price_cents` INT (1-99)
  - `count` INT (number of contracts traded)
  - `notional_usd` DOUBLE (computed: count × price_cents × 0.01 × $1)
  - `is_large` BOOLEAN (true if notional ≥ threshold; threshold configurable, default $500)
  - `traded_at` TIMESTAMP
  - `fetched_at` TIMESTAMP
  - INDEX on (ticker, traded_at DESC)

- `draft_bets` — light bet logging
  - `bet_id` UUID PRIMARY KEY
  - `market_id` TEXT FK
  - `book` TEXT
  - `side` TEXT
  - `american_odds` INT
  - `stake_usd` DOUBLE
  - `taken_at` TIMESTAMP
  - `note` TEXT (optional free-form, e.g., "after Schefter Ward news")

- `players_2026` — canonical name registry
  - `canonical_name` TEXT PRIMARY KEY
  - `position` TEXT
  - `college` TEXT (for disambiguation)
  - `aliases` JSON (list of strings, e.g., `["Cam Ward", "Cameron Ward", "C. Ward"]`)

  Loaded once from a hand-maintained `players_2026.csv` in `nfl_draft/config/`. Owned manually because rookie-name canonicalization is unavoidable cross-book and any auto-fuzzy-match will produce silent join errors. Same pattern for `teams_nfl` (already exists in canonical form via `Answer Keys/canonical_match.py`; reuse).

### Code layout

```
nfl_draft/
  config/
    players_2026.csv           # canonical name + aliases
    market_map.json            # per-book → canonical market_id mapping
  scrapers/
    kalshi.py                  # extends existing kalshi_draft/fetcher.py for ALL series
    draftkings.py              # adapted from mlb_sgp/scraper_draftkings_sgp.py
    fanduel.py                 # adapted from mlb_sgp/scraper_fanduel_sgp.py
    bookmaker.py               # adapted from bookmaker_odds/scraper.py
    wagerzon.py                # adapted from wagerzon_odds/scraper_v2.py
  lib/
    db.py                      # DuckDB connection + schema migrations
    devig.py                   # ported from Answer Keys/Tools.R (~30 LOC)
    normalize.py               # player/team alias resolution
    market_map.py              # per-book market label → canonical market_id
  run.py                       # orchestrator: --mode {pre-draft, draft-day} --book all
  app.py                       # extends kalshi_draft/app.py with new tabs (do NOT fork)
  README.md
  nfl_draft.duckdb             # gitignored
```

The dashboard **extends** `kalshi_draft/app.py` rather than living separately. New tabs read from the new `nfl_draft.duckdb`; existing tabs continue to read from `kalshi_draft.duckdb`. Two DBs is fine — Dash callbacks open whichever connection a tab needs. No schema sharing required.

### Devig math

Port three functions from `Answer Keys/Tools.R` to `nfl_draft/lib/devig.py`:

- `american_to_implied(american_odds) -> float`
- `devig_two_way(odds_a, odds_b) -> (prob_a, prob_b)`
- `devig_n_way(odds_list) -> List[float]` (for first-QB / first-RB style multi-outcome markets where the complementary set is the full position pool)

For markets where the complementary set isn't fully posted (e.g., book only offers Yes on 5 of the 30+ first-QB candidates), use **proportional devig**: divide each implied prob by the sum of all posted implieds for that market group on that book. Imperfect but standard practice.

### Scraper orchestration

`nfl_draft/run.py` is the single entry point.

```
python run.py --mode pre-draft --book all
python run.py --mode draft-day --book all
python run.py --mode draft-day --book kalshi   # one book at a time, for debugging
```

- Each scraper module exposes `fetch_draft_odds() -> List[OddsRow]`.
- `run.py` calls them sequentially (parallelism deferred — scrape errors easier to debug serially), normalizes each result through `lib/normalize.py` and `lib/market_map.py`, devigs via `lib/devig.py`, writes to `draft_odds`.
- Kalshi trade tape is its own subcommand: `python run.py --mode trades` polls `/markets/trades` since the last `fetched_at` and appends to `kalshi_trades`.

**Cron** (macOS launchd or simple `cron`):

- Pre-draft: `*/15 * * * *` for `--mode pre-draft --book all` (15-minute cycle)
- Pre-draft: `*/2 * * * *` for `--mode trades` (Kalshi tape every 2 min — public endpoint, low cost)
- Draft-day: manually swap to `*/2` for odds and `*/1` for trades the morning of April 23. Two crontab files (`crontab.pre`, `crontab.draft`) shipped in the repo, swap with one command.

### Kalshi expansion

Today `kalshi_draft/fetcher.py` discovers via `/series` keyword filter ("nfl"+"draft") plus a hardcoded prefix fallback. v1 changes:

- Add all observed `KXNFLDRAFT*` prefixes to the fallback list (incl. `KXNFLDRAFTQB`, `KXNFLDRAFTWR`, `KXNFLDRAFTRB`, `KXNFLDRAFTPOS`, `KXNFLDRAFTTOP5`, `KXNFLDRAFTTOP10`, `KXNFLDRAFTPICK`, `KXNFLDRAFTTOP`).
- After discovery, log every series found that has open markets but no downstream consumer in the dashboard — these are markets we're storing but not yet rendering. Used to catch new Kalshi series (e.g., a `KXNFLDRAFTTRADE` for "first trade" if Kalshi adds one) without code changes.
- Map each Kalshi market to a canonical `market_id` via `nfl_draft/config/market_map.json` (Kalshi-side keyed by ticker prefix + market subject extracted from the market title).

### Per-book scraper notes

- **Kalshi**: existing `kalshi_draft/fetcher.py` is the base; `nfl_draft/scrapers/kalshi.py` is a thin adapter that calls it and remaps to the new schema.
- **DraftKings**: reuse the CDP click-and-capture + Akamai handling pattern from `mlb_sgp/scraper_draftkings_sgp.py`. Reconnaissance step required: locate the NFL Draft event page, capture the API call(s) that fetch the draft markets menu.
- **FanDuel**: reuse the curl_cffi + 3-header recipe (per memory `fanduel_sgp_scraping.md`). Reconnaissance: same as DK.
- **Bookmaker**: reuse the curl_cffi + ASP.NET cookie pattern from `bookmaker_odds/scraper.py`. Reconnaissance: discover the league ID for NFL Draft markets (currently only cbb/nba/mlb configured).
- **Wagerzon**: reuse the ASP.NET `__VIEWSTATE` POST pattern from `wagerzon_odds/scraper_v2.py`. Reconnaissance: same — find draft league IDs.

For each book, reconnaissance and parsing is the unknown — auth and request mechanics are solved.

### Dashboard tabs (extending `kalshi_draft/app.py`)

New tabs added to the existing Dash app:

1. **Cross-Book Grid** — Markets on rows, books on columns, devigged prob in cells. Final column shows variance (max − min across books that posted the market). Click a row to drill into a per-market detail view (line history if data exists, all bid/ask, take/log button).
2. **+EV Candidates** — Same data, but flat-listed and ranked by variance descending. Filter by market type, position, book. This is the user's "scan the field" view.
3. **Trade Tape** — Streaming-style table of recent Kalshi trades (last 200 rows), large fills (≥ $500 notional, configurable) highlighted. Filter by ticker / series.
4. **Bet Log** — Form to log a bet (market_id dropdown, book dropdown, american_odds, stake, note). Below the form, table of past bets.

Existing tabs (Market Overview, Price History, Edge Detection, Consensus, Portfolio) stay unchanged in v1.

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
   kalshi_trades         (append-only, is_large flag set on insert)

[ Dash app — kalshi_draft/app.py extended ]
        │
        ├── Cross-Book Grid ←  reads draft_markets + latest draft_odds per (market, book)
        ├── +EV Candidates  ←  same query, flat + sorted by variance
        ├── Trade Tape      ←  reads kalshi_trades ORDER BY traded_at DESC LIMIT 200
        ├── Bet Log         ←  writes draft_bets; reads it for table view
        └── (existing tabs unchanged)
```

---

## Error handling

- **Per-scraper isolation**: a failure in one book's scraper must not prevent the others from running. `run.py` wraps each book in try/except, logs the error, continues.
- **Auth failure**: log and surface in dashboard footer (e.g., "DK auth expired — re-login needed"). Don't silently lose data.
- **New / unknown market**: when a scraper emits a market label that doesn't map to a canonical `market_id`, log a warning and write the row to a `draft_odds_unmapped` quarantine table. Daily review during draft week to extend `market_map.json`. Better than silently dropping.
- **Player alias miss**: same pattern — quarantine to `draft_odds_unmapped_players`, daily review to extend `players_2026.csv`. Surface unmapped count on the dashboard footer.
- **Stale data warning**: dashboard shows `fetched_at` per book. If any book is > 30 min stale during pre-draft mode (or > 5 min during draft mode), highlight in red.
- **DuckDB lock conflicts**: `run.py` and `app.py` both touch the DB. Use DuckDB's `read_only=true` for the dashboard (Dash callbacks open short-lived read-only connections) and `read_only=false` for `run.py`. If contention emerges, switch to a write-only worker that batches inserts.

---

## Testing

Phase 1 is time-boxed; testing is pragmatic, not exhaustive.

- **Unit tests**: only for `lib/devig.py` and `lib/normalize.py` — these are deterministic math/lookup, easy to test, and incorrect output silently corrupts every downstream calculation. ~10 tests total. Pytest.
- **Scraper smoke test**: each scraper module has a `__main__` block that runs once and prints rows. Used during reconnaissance and as the smoke check after auth changes.
- **End-to-end manual check**: the morning of April 22, run `--mode pre-draft --book all` end-to-end, verify all 5 venues populated `draft_odds` rows for at least the #1 overall market, verify the dashboard renders. This is the go/no-go gate.
- **No CI**: this is a personal project on a 7-day clock. Tests run locally on demand.

---

## Phasing

### Phase 1 — Ship by Apr 23 (this week)

In scope:
- DuckDB schema, all tables
- Player/team canonical config (manual, but populated for the top ~80 prospects)
- All 5 venue scrapers (Kalshi expansion + 4 books) writing to `draft_odds`
- Kalshi trade tape polling → `kalshi_trades`
- Devig math ported to Python
- Dashboard tabs: Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log
- Cron setup: pre-draft + draft-day swap
- README in `nfl_draft/`

Explicit risk: if reconnaissance on any book exceeds its 1-day budget, skip and ship without it. Better to ship 3 books on time than 4 books a day late.

### Phase 2 — Post-draft, weeks following

Out of scope for v1, planned for after:
- Formal fair-value oracle (sharp-weighted consensus, possibly with mock-draft prior)
- Real-time Kalshi WebSocket trade tape (replace polling)
- Line-movement charts per market
- P&L view + CLV analysis on `draft_bets`
- Kalshi MM adapted from `kalshi_mm/`
- Generalize the framework for NFL Draft 2027, NBA Draft, NHL Draft (rename `nfl_draft/` → `draft/`, parameterize sport)

### Phase 3 — Out

- Auto-betting, Kelly sizing recommendations, alerting

---

## Version control

- **Branch**: `feature/nfl-draft-portal` (already created at brainstorm time)
- **Files created**: `nfl_draft/` directory tree (see Code layout); `docs/superpowers/specs/2026-04-17-nfl-draft-portal-design.md`
- **Files modified**: `kalshi_draft/app.py` (new tabs), `.gitignore` (add `nfl_draft/nfl_draft.duckdb`, `nfl_draft/.cookies/`), `README.md` (link to nfl_draft/README.md)
- **Commit structure**: per-phase commits where possible — schema commit, then each scraper as its own commit, then dashboard, then cron config. Review checkpoint at each scraper before merging.
- **Worktree**: optional but recommended given other active work on `feature/cbb-consensus-calibration` and `feature/mlb-sgp-scrapers`. Use `/worktree` for isolation if any of those are likely to need attention during the 7 days.
- **Pre-merge review**: full executive review of `git diff main..HEAD` before merging — focus on (a) no DK/FD credentials in code, (b) DuckDB connections all closed, (c) no logging of personal bets to stdout, (d) cookie files gitignored.
- **Approval to merge**: explicit ask before `git merge` or any push.

## Documentation

Required updates in the same final merge:
- `nfl_draft/README.md` — setup (env vars, auth steps per book, cron install), usage (`run.py` flags), troubleshooting (auth failures, schema migrations).
- Top-level `README.md` — one-line link to `nfl_draft/README.md`.
- `CLAUDE.md` (project root) — add `nfl_draft/` to the project structure section so future Claude sessions discover it.
- `kalshi_draft/README.md` — note that the dashboard has been extended with draft-portal tabs reading from `nfl_draft.duckdb`.

---

## Open questions deferred to Phase 2

- How should fair-value be computed? (Sharp-weighted consensus? Mock-draft Bayesian prior? Hybrid with user override?)
- Is Kalshi a market or an oracle? (Affects EV direction — bet against Kalshi when books move, or bet against books when Kalshi prints heavy?)
- Threshold for "large bet" on Kalshi tape — default $500 notional in v1; tune from real data after the draft.

---

## Success criteria for Phase 1

By end of day April 22 (the day before draft Day 1):

1. All 5 venues populating `draft_odds` for at least the #1-overall and first-QB markets, refreshed in the last 15 minutes
2. Cross-Book Grid renders all populated markets with devigged probs and variance
3. Kalshi Trade Tape shows the last hour of trades, large fills highlighted
4. Bet Log accepts and persists a test bet
5. Cron jobs running on the user's machine; draft-day swap rehearsed
