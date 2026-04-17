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

A single DuckDB database at `nfl_draft/nfl_draft.duckdb` — the **only** draft database for the project. The existing `kalshi_draft/kalshi_draft.duckdb` is retired in v1: its data (a few months of Kalshi draft odds + portfolio history) is migrated into the new DB at the start of implementation, then the old file is deleted. Two DBs for the same logical concern (NFL draft data) is confusing; one is correct.

**Migration of existing data**: a one-time Python script (`nfl_draft/lib/migrate_from_kalshi_draft.py`) copies tables from `kalshi_draft/kalshi_draft.duckdb` into `nfl_draft/nfl_draft.duckdb` with renames where needed:
- `draft_odds` (Kalshi-only) → `kalshi_odds_history` (preserved for back-reference; new portal writes go to the new `draft_odds` table described below)
- `draft_series`, `market_info` → kept as-is
- `positions`, `resting_orders` → kept as-is

After migration is verified, `kalshi_draft/kalshi_draft.duckdb` is removed and the existing dashboard tabs in `kalshi_draft/app.py` are updated to point at `nfl_draft/nfl_draft.duckdb`.

**No local CSV / JSON config files** (per project rule "no temp files"). All static configuration — player aliases, team aliases, per-book market label mappings — lives as Python dict literals in `nfl_draft/config/*.py` modules (mirrors the existing `wagerzon_odds/config.py` and `bookmaker_odds/scraper.py` pattern). At runtime these dicts are loaded into in-memory lookups; for cross-process consistency they are also seeded into DuckDB lookup tables (`players`, `teams`, `market_map`) by a one-time `python -m nfl_draft.lib.seed` command and re-run whenever the source dict changes.

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

- `players` — canonical name registry, seeded from `nfl_draft/config/players.py` (Python dict literal)
  - `canonical_name` TEXT PRIMARY KEY
  - `position` TEXT
  - `college` TEXT (for disambiguation)
  - `aliases` TEXT[] (DuckDB list type, e.g., `['Cam Ward', 'Cameron Ward', 'C. Ward']`)

- `teams` — canonical NFL team registry, seeded from `nfl_draft/config/teams.py` (Python dict literal). Built off the existing canonical mapping in `Answer Keys/canonical_match.py` — do not duplicate.

- `market_map` — per-book label → canonical `market_id` mapping, seeded from `nfl_draft/config/markets.py` (Python dict literal)
  - `book` TEXT (`kalshi` | `draftkings` | `fanduel` | `bookmaker` | `wagerzon`)
  - `book_label` TEXT (the raw market title or ticker prefix from the book)
  - `book_subject` TEXT (raw player/team string from the book, before normalization)
  - `market_id` TEXT FK → `draft_markets.market_id`
  - PRIMARY KEY (book, book_label, book_subject)

The three lookup tables are populated by a one-time seed command (`python -m nfl_draft.lib.seed`) that reads the dict literals from `nfl_draft/config/*.py` and writes to DuckDB. Editing the source `.py` file + re-running `seed` is how aliases and mappings are maintained. The dict literals are version-controlled (source code), so changes show up in `git diff` cleanly. No CSV or JSON config files anywhere.

Manual maintenance is unavoidable for rookie-name canonicalization — any auto-fuzzy-match will produce silent join errors that silently corrupt EV calculations. Better to fail loudly: scraper rows with unmapped players land in a quarantine table (see Error handling) and the dashboard footer surfaces the unmapped count, which is the user's daily prompt to extend the dict.

### Code layout

```
nfl_draft/
  config/
    players.py                 # PLAYERS dict literal (canonical → metadata + aliases)
    teams.py                   # TEAMS dict literal (canonical → aliases)
    markets.py                 # MARKET_MAP dict literal (per-book label → canonical market_id)
  scrapers/
    kalshi.py                  # extends existing kalshi_draft/fetcher.py for ALL series
    draftkings.py              # adapted from mlb_sgp/scraper_draftkings_sgp.py
    fanduel.py                 # adapted from mlb_sgp/scraper_fanduel_sgp.py
    bookmaker.py               # adapted from bookmaker_odds/scraper.py
    wagerzon.py                # adapted from wagerzon_odds/scraper_v2.py
  lib/
    db.py                      # DuckDB connection + schema migrations
    devig.py                   # ported from Answer Keys/Tools.R (~30 LOC)
    normalize.py               # reads players/teams tables → alias resolution
    market_map.py              # reads market_map table → per-book label → canonical market_id
    seed.py                    # one-time: reads config/*.py dicts, writes to DuckDB lookup tables
    migrate_from_kalshi_draft.py  # one-time: migrates kalshi_draft.duckdb → nfl_draft.duckdb
  run.py                       # orchestrator: --mode {pre-draft, draft-day, trades} --book all
  README.md
  nfl_draft.duckdb             # gitignored — the SOLE draft database
```

The dashboard **extends** `kalshi_draft/app.py` rather than living separately, and is updated to point all tabs (existing + new) at `nfl_draft/nfl_draft.duckdb`. The old `kalshi_draft/kalshi_draft.duckdb` is deleted after migration.

**Directory naming wart, accepted in Phase 1**: the `kalshi_draft/` directory now houses a *multi-venue* draft dashboard, not just Kalshi. Renaming would break `kalshi_mm/` which imports `kalshi_draft/auth.py` (4 hard references in `kalshi_mm/{orders,analyze_performance,config,dashboard}.py`) — the running CBB MM bot would go down. Phase 2 cleanup: rename `kalshi_draft/` → `kalshi_lib/` (it's really the Kalshi auth+fetcher library that `kalshi_mm` depends on), update kalshi_mm imports, and move `app.py` into `nfl_draft/`.

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
- **New / unknown market**: when a scraper emits a market label that doesn't map to a canonical `market_id`, log a warning and write the row to a `draft_odds_unmapped` quarantine table. Daily review during draft week to extend `nfl_draft/config/markets.py` (then re-run `python -m nfl_draft.lib.seed`). Better than silently dropping.
- **Player alias miss**: same pattern — quarantine to `draft_odds_unmapped_players`, daily review to extend `nfl_draft/config/players.py` (then re-run `seed`). Surface unmapped count on the dashboard footer.
- **Stale data warning**: dashboard shows `fetched_at` per book. If any book is > 30 min stale during pre-draft mode (or > 5 min during draft mode), highlight in red.
- **DuckDB lock conflicts**: `run.py` and `app.py` both touch the DB. Use DuckDB's `read_only=true` for the dashboard (Dash callbacks open short-lived read-only connections) and `read_only=false` for `run.py`. If contention emerges, switch to a write-only worker that batches inserts.

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
    test_migration.py                # kalshi_draft.duckdb → nfl_draft.duckdb roundtrip
    test_dashboard_queries.py        # query functions return correct shape + values
    test_quarantine.py               # unmapped player/market lands in quarantine table
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
- `test_dashboard_queries.py`: seed temp DB with known data, call each tab's query function (broken out from Dash callbacks), verify result shape matches dashboard expectations. **This is the test that catches "existing tabs break after DB repoint" — every existing tab gets a query test before its query is touched.**
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
- DuckDB schema, all tables
- Migration of existing `kalshi_draft.duckdb` data into `nfl_draft.duckdb`; deletion of the old DB; existing dashboard tabs repointed (with query tests covering each tab before the repoint)
- Player/team/market canonical config in `nfl_draft/config/*.py` (Python dict literals, populated for the top ~80 prospects), seeded into DuckDB lookup tables
- All 5 venue scrapers (Kalshi expansion + 4 books) writing to `draft_odds`
- Kalshi trade tape polling → `kalshi_trades`
- Devig math ported to Python
- Dashboard tabs: Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log
- Cron setup: pre-draft + draft-day swap
- **Comprehensive test suite** (unit + integration + live; see Testing section). Tests are not optional — every `lib/`, `scrapers/`, and dashboard query function ships with tests. Pre-commit hook runs unit + integration on every commit.
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
- **Files modified**: `kalshi_draft/app.py` (new tabs + repointed at `nfl_draft/nfl_draft.duckdb`), `.gitignore` (add `nfl_draft/nfl_draft.duckdb`, `nfl_draft/.cookies/`), `README.md` (link to nfl_draft/README.md)
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
2. Cross-Book Grid renders all populated markets with devigged probs and variance
3. Kalshi Trade Tape shows the last hour of trades, large fills highlighted
4. Bet Log accepts and persists a test bet
5. Cron jobs running on the user's machine; draft-day swap rehearsed
6. **All tests green**: `pytest tests/unit tests/integration` returns 0 failures; `pytest tests/live` returns 0 failures on the morning-of run; existing dashboard tabs (Market Overview, Price History, Edge Detection, Consensus, Portfolio) verified to render correct data after the DB repoint via `test_dashboard_queries.py`.
