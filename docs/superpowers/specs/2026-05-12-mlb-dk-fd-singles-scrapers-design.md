# MLB DK/FD Single-Leg Scrapers — Design Spec

**Date:** 2026-05-12
**Branch:** `feature/mlb-dk-fd-singles-scrapers`
**Status:** Design — pending user approval before writing implementation plan

---

## Motivation

The MLB Dashboard bets-tab card layout (merged on `claude/mlb-sportsbook-comparison-98fM8`) renders per-book pills for every bet. DraftKings, FanDuel, and Pinnacle pills currently come from the Odds API (`prefetched_odds` JSON cache). The Odds API has structural gaps:

- Pinnacle is available only on the `eu` region.
- F-period markets (F3 / F7) are not reliably returned for DK / FD even when the books actually post them.
- Alt-line coverage is patchy.

DK and FD have publicly accessible REST APIs that already power the existing `mlb_sgp/scraper_*.py` files. This project replaces the Odds API as the source for DK and FD pill data with dedicated single-leg scrapers built on those same APIs. **Pinnacle keeps using the Odds API** — it has no public API and a Pinnacle-specific scraper is out of scope.

## Goals

1. Replace Odds API as the source for **DraftKings and FanDuel** single-leg pill data on the MLB Dashboard bets tab.
2. Cover **all market types the model bets on**: FG + F5 + F3 + F7, including main and alt lines for spreads and totals, plus h2h / moneyline where posted.
3. Run as part of the existing `run.py mlb` orchestrator, alongside the offshore scrapers, gated by the standard `.scrapers_done_mlb` sentinel.
4. Reuse existing DK/FD API code through a shared client abstraction; do not duplicate.
5. Leave existing SGP scrapers functionally unchanged.

## Non-goals

- Pinnacle scraping (kept on Odds API).
- Restructuring the SGP scrapers to drop their dependency on `mlb_parlay_lines` (separate workstream).
- Falling back to Odds API when a DK/FD scrape fails — failures degrade pills to empty, same as today's behavior when Odds API misses a book.
- Adding new market types (props, futures, alternate moneylines) — scope is limited to the markets `all_bets_combined` already produces.

---

## Architecture

### File / directory layout

```
mlb_sgp/
├── dk_client.py                       [NEW]       DraftKingsClient class
├── fd_client.py                       [NEW]       FanDuelClient class
├── scraper_draftkings_singles.py      [NEW]       Singles scrape → dk_odds/dk.duckdb
├── scraper_fanduel_singles.py         [NEW]       Singles scrape → fd_odds/fd.duckdb
├── scraper_draftkings_sgp.py          [REFACTOR]  Imports dk_client; behavior unchanged
├── scraper_fanduel_sgp.py             [REFACTOR]  Imports fd_client; behavior unchanged
├── scraper_novig_sgp.py               [unchanged]
├── scraper_prophetx_sgp.py            [unchanged]
├── scraper_pikkit_mlb.py              [unchanged]
└── README.md                          [EDIT]      +"Singles scrapers" section

dk_odds/                               [NEW DIR]   Mirrors offshore scraper layout
└── dk.duckdb                          [NEW]       Table: mlb_odds

fd_odds/                               [NEW DIR]
└── fd.duckdb                          [NEW]       Table: mlb_odds

Answer Keys/
├── run.py                             [EDIT]      +2 SCRAPER_CONFIGS entries
├── MLB Answer Key/MLB.R               [EDIT]      book_odds_by_book DK/FD lines
└── CLAUDE.md                          [EDIT]      Pipeline diagram updated
```

### Client classes — surface

`mlb_sgp/dk_client.py`:

```python
class DraftKingsClient:
    def __init__(self, verbose: bool = False) -> None:
        # Owns the curl_cffi Chrome-TLS session that bypasses Akamai.

    def list_events(self) -> list[Event]:
        # GET sportsbook-nash/.../leagueSubcategory/v1/markets
        # Returns one Event per MLB game today: (event_id, home, away, start_time).

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        # GET event/eventSubcategory/v1/markets
        # Returns market metadata: (market_id, name, subcategory).
        # Used to identify which markets are FG vs F5 vs F3 vs F7,
        # and which are main vs alt.

    def fetch_event_selections(self, event_id: str) -> list[Selection]:
        # GET parlays/v1/sgp/events/{id}
        # Returns ALL selections with their underlying prices:
        # (selection_id, market_id, name, line, american_odds).
        # Same call the SGP scraper already makes.
```

`mlb_sgp/fd_client.py` mirrors the surface: `list_events`, `fetch_event_runners`. Method names reflect FD terminology (runners) rather than DK's (selections).

### Singles scraper logic (skeleton — both books)

```python
def scrape_singles():
    client = DraftKingsClient()
    events = client.list_events()
    rows = []
    for event in events:
        try:
            selections = client.fetch_event_selections(event.id)
            rows.extend(parse_selections_to_wide_rows(event, selections))
        except Exception as e:
            log_warning(f"Skipping {event.id}: {e}")
            continue                          # per-game isolation
    write_to_duckdb(rows, "dk_odds/dk.duckdb", "mlb_odds")  # atomic write
```

### SGP scraper refactor

`scraper_draftkings_sgp.py` and `scraper_fanduel_sgp.py` import from the new client modules instead of inlining the curl_cffi session setup, event discovery, and selection-ID parsing. The SGP combo logic and `calculateBets` integration stay in the SGP files. **No behavior change** — regression-tested against a captured golden output.

---

## Data flow

### Output schema

Both singles scrapers write to a `mlb_odds` table in their per-book DuckDB, matching the offshore convention exactly (byte-identical to `wagerzon.duckdb::mlb_odds`):

```
fetch_time        TIMESTAMP
sport_key         VARCHAR     -- "baseball_mlb"
game_id           VARCHAR     -- DK/FD event id (informational; MLB.R joins by team names)
game_date         VARCHAR
game_time         VARCHAR
away_team         VARCHAR     -- canonical (DK/FD return full canonical names)
home_team         VARCHAR     -- canonical
market            VARCHAR     -- "main" | "alternate_spreads" | "alternate_totals"
period            VARCHAR     -- "FG" | "F5" | "F3" | "F7"
away_spread       FLOAT       -- NULL if not a spread row
away_spread_price INTEGER
home_spread       FLOAT
home_spread_price INTEGER
total             FLOAT       -- NULL if not a totals row
over_price        INTEGER
under_price       INTEGER
away_ml           INTEGER     -- main rows only
home_ml           INTEGER     -- main rows only
```

Volume estimate: ~30 games × 4 periods × (1 main + ~5 alt-spread + ~5 alt-total rows) ≈ ~1,300 rows per scrape per book. Trivial DB footprint.

### MLB.R integration

`Answer Keys/MLB Answer Key/MLB.R`, around lines 882–894:

```r
# BEFORE
book_odds_by_book <- list(
  wagerzon  = scraper_to_canonical(wagerzon_odds,  .game_id_lookup),
  hoop88    = scraper_to_canonical(hoop88_odds,    .game_id_lookup),
  bfa       = scraper_to_canonical(bfa_odds,       .game_id_lookup),
  bookmaker = scraper_to_canonical(bookmaker_odds, .game_id_lookup),
  bet105    = scraper_to_canonical(bet105_odds,    .game_id_lookup),
  draftkings = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "draftkings")),
  fanduel    = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "fanduel")),
  pinnacle   = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "pinnacle"))
)

# AFTER
book_odds_by_book <- list(
  wagerzon  = scraper_to_canonical(wagerzon_odds,  .game_id_lookup),
  hoop88    = scraper_to_canonical(hoop88_odds,    .game_id_lookup),
  bfa       = scraper_to_canonical(bfa_odds,       .game_id_lookup),
  bookmaker = scraper_to_canonical(bookmaker_odds, .game_id_lookup),
  bet105    = scraper_to_canonical(bet105_odds,    .game_id_lookup),
  draftkings = scraper_to_canonical(dk_odds,       .game_id_lookup),
  fanduel    = scraper_to_canonical(fd_odds,       .game_id_lookup),
  pinnacle   = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "pinnacle"))
)
```

`dk_odds` and `fd_odds` are loaded earlier with `dbGetQuery(con, "SELECT * FROM mlb_odds")` from `dk_odds/dk.duckdb` and `fd_odds/fd.duckdb` — same pattern as `wagerzon_odds`. `parse_prefetched_to_long()` still runs but only emits Pinnacle rows.

### Team-name resolution

DK and FD both return full team names like `"Kansas City Royals"` — same canonical format as the Odds API. The existing `resolve_offshore_teams()` in `Tools.R` already handles this; no new dictionary entries expected. **Confirm during implementation**, and add entries if any drift is found.

---

## Orchestration

### `run.py` changes

Two new entries appended to `SCRAPER_CONFIGS` in `Answer Keys/run.py`:

```python
"draftkings_singles": {
    "script": "../mlb_sgp/scraper_draftkings_singles.py",
    "sports": ["mlb"]
},
"fanduel_singles": {
    "script": "../mlb_sgp/scraper_fanduel_singles.py",
    "sports": ["mlb"]
},
```

The existing `run_scraper()` and sentinel logic at line 183 (`for name, scraper_config in SCRAPER_CONFIGS.items()`) picks these up automatically.

### One-cycle timeline

```
T+0min   run.py mlb starts
   ├──── Phase 1: sharp scrapers run first (bookmaker, bet105)
   ├──── Phase 2: parallel scrapers fire:
   │       • wagerzon, hoop88, bfa, kalshi
   │       • draftkings_singles  [NEW]
   │       • fanduel_singles     [NEW]
   ├──── Sentinel fires: .scrapers_done_mlb
T+~5m    MLB.R reads scraper outputs (dk.duckdb, fd.duckdb new)
T+~10m   MLB.R writes mlb_bets_combined + mlb_bets_book_prices
   ▼
[separate trigger, unchanged]
T>>10m   scraper_draftkings_sgp.py + scraper_fanduel_sgp.py run
         (reads mlb_parlay_lines, uses refactored clients internally)
```

---

## Failure modes & graceful degradation

| Failure | Effect | User-visible result |
|---|---|---|
| DK Akamai block (whole scrape fails) | `mlb_odds` empty for DK | DK pills show `—` (same as today when Odds API misses DK) |
| FD endpoint changes shape | Parsing fails, scraper exits non-zero | FD pills empty; sentinel still fires (verify against existing scraper behavior in impl) |
| Partial scrape (some games fail) | Only successful games have rows | Pills empty for missing games only |
| F3/F7 not posted by FD | Scraper writes only FG + F5 rows for FD | F3/F7 FD pills empty (same as today) |
| One DK selection has malformed price | Skip that row, continue | Single pill missing; rest of scrape unaffected |

Logging: `mlb_sgp/logs/dk_singles_YYYY-MM-DD.log`, `fd_singles_YYYY-MM-DD.log`. Existing log-rotation pattern in `mlb_sgp/logs/` continues.

---

## Open uncertainties (verify during implementation, not blocking design)

1. **Does the DK `parlays/v1/sgp/events/{id}` payload include single-leg prices, or just selection IDs?** If just IDs, the scraper adds one extra call per event to `event/eventSubcategory/v1/markets` (which returns prices for game lines). Either way the design holds.
2. **Does FD post F3/F7 markets?** If not, those rows are absent and FD F3/F7 pills stay empty (matching today).
3. **DK/FD team-name drift.** Confirm both return canonical names; add `resolve_offshore_teams()` entries if not.

---

## Testing strategy

- **Unit tests** for client classes — mock HTTP, assert parsers produce expected `Selection` / `Runner` objects from captured fixtures committed under `mlb_sgp/tests/fixtures/`.
- **Integration test** (network, opt-in) — run full singles scrape against live API once; assert `mlb_odds` is non-empty and schema-conforming.
- **MLB.R smoke** — run pipeline once after wiring; verify on a known game that bets-tab DK/FD pills show expected odds.
- **SGP regression test** — after `scraper_draftkings_sgp.py` and `scraper_fanduel_sgp.py` import from the new clients, confirm `mlb_sgp_odds` rows for a sample game match the unrefactored version (golden-file comparison).

---

## Documentation surface

Updates required before merge (exact wording in the implementation plan):

- `mlb_sgp/README.md` — new "Singles scrapers" section explaining files and their relationship to the SGP scrapers.
- `Answer Keys/CLAUDE.md` — Pipeline Flow diagram shows DK/FD as scrapers in the parallel-scrape phase.
- `Answer Keys/MLB Dashboard/PLAN_odds_screen.md` (or the project root MLB section) — DK/FD pill data now sourced from scrapers; Pinnacle still from Odds API.

---

## Out of scope (named for clarity)

- Pinnacle scraper.
- SGP scraper restructure to drop `mlb_parlay_lines` dependency.
- New market types (props, futures, alternate moneylines).
- Frontend changes — the dashboard is agnostic to the upstream source.
- Performance optimization beyond per-game isolation and atomic DB write.

---

## Phase 1 findings (2026-05-12)

- **DK selection price field:** `displayOdds.american` (string e.g. `"+260"`) plus `displayOdds.decimal` and a top-level numeric `trueOdds` (present in `parlays/v1/sgp/events/{id}` payload: **YES**). The parlays endpoint returns full `data.markets[].selections[]` with prices baked in, so the singles scraper is **one-call-per-event** — no need to also call `event/eventSubcategory/v1/markets`. The subcategory endpoint only returns `(market_id, name)` tuples (no selections), so it's not useful for pricing. Selection objects also expose `outcomeType` (`Over`/`Under`/team-side), `points` (line value as float), `name`, `status`, `isDisabled`, and the team-prefixed market `name` (e.g. `"Run Line"`, `"Total"`, `"Run Line - 1st 5 Innings"`).
- **DK F3 markets:** `ABSENT` — no markets in the parlays payload contain `"1st 3"` / `"F3"` / `"3rd Innings"`. DK does post per-inning `3rd Inning (3 Way)` style markets but no cumulative first-3-innings spread/total.
- **DK F5 markets:** `PRESENT` — exact `marketType.name` values: `"Run Line - 1st 5 Innings"`, `"Total Runs - 1st 5 Innings"`, `"Total Alternate - 1st 5 Innings"`, `"Team Total Runs - 1st 5 Innings"`.
- **DK F7 markets:** `PRESENT` — exact `marketType.name` values: `"Run Line - 1st 7 Innings"`, `"Total Runs - 1st 7 Innings"`, `"Total Alternate - 1st 7 Innings"`, `"Team Total Runs - 1st 7 Innings"`.
- **FD F3/F5/F7 markets:**
  - F3 spread/total: `ABSENT` — only `"First 3 Innings Result"` (3-way moneyline-style) exists; no run line, no total.
  - F5: `PRESENT` — `"First 5 Innings Run Line"`, `"First 5 Innings Alternate Run Lines"`, `"First 5 Innings Total Runs"`, `"First 5 Innings Alternate Total Runs"`, `"First 5 Innings Money Line"`.
  - F7 spread/total: `ABSENT` — only `"First 7 Innings Result"`; no spread, no total.
- **FD runner price field:** `runner.winRunnerOdds.americanDisplayOdds.americanOdds` (integer, e.g. `520`). Also available: `winRunnerOdds.americanDisplayOdds.americanOddsInt` (same value), `winRunnerOdds.trueOdds.decimalOdds.decimalOdds` (float). Runners also expose `handicap` (numeric line), `runnerName`, `runnerStatus`, `previousWinRunnerOdds` (history list).
- **DK team-name drift:** `LIST` — DK uses abbreviated city prefixes ("CLE Guardians", "LA Angels", "STL Cardinals") for **all 30 teams**; canonical (`mlb_consensus_temp`) uses full names ("Cleveland Guardians", "Los Angeles Angels", "St. Louis Cardinals"). Existing `resolve_offshore_teams()` in `Tools.R` does substring/word-overlap matching but DK's prefix form (e.g. "cleguardians" stripped) is not a substring of canonical "clevelandguardians" — needs explicit mapping or game-level word-overlap rescue. Note that the SGP scraper today already handles this via the `mlb_parlay_lines` join on team names from the R-side resolver, so DK is already being canonicalized somewhere — but the singles scraper must not skip that step.
- **FD team-name drift:** `NONE` — all 30 FD team names match canonical exactly (e.g. `"St. Louis Cardinals"`, `"Athletics"`, `"Cleveland Guardians"`).

**Implications for plan:**
- Tasks 4 + 5 client parsers should use field names: DK = `displayOdds.american` (parse `"+260"` → `260`) extracted from `data.markets[].selections[]` in the parlays endpoint response; FD = `runner.winRunnerOdds.americanDisplayOdds.americanOdds` (already an integer).
- Task 4 (DK singles parser) only needs **one HTTP call per event** (`parlays/v1/sgp/events/{id}`). The subcategory markets endpoint can be skipped for pricing — only useful if a future task needs market metadata lookup independent of selections.
- If FD doesn't post F3/F7 spread/total (confirmed), `scraper_fanduel_singles.py` emits only FG + F5 rows (matches today's SGP scraper behavior).
- DK posts F5 + F7 spread/total. The singles scraper for DK should emit rows for FG + F5 + F7 (F7 is a **new** period not currently produced by the SGP scraper — confirm this is desired before Task 9 wires DB writes).
- Team-name dictionary updates: `NEEDED` for DK. Two options: (a) extend `Answer Keys/Tools.R::resolve_offshore_teams()` substring rules to handle 2-3-char city prefixes ("CLE" → "Cleveland", "LA" → "Los Angeles"), or (b) build a DK-specific mapping table inside `scraper_draftkings_singles.py` keyed off the team-dict canonical list. Recommend (b) — keeps the resolver generic and the mapping co-located with the scraper that introduces the variant naming. FD needs no mapping changes.

---

## Pre-merge executive review (2026-05-12)

**Scope:** singles scrapers feature only. The odds-screen feature (merged in
via 6056fd5) was reviewed separately on its own branch.

**Reviewed commits (singles-only, 13 non-merge):** 1ad39bb, 389e4d9, f1c074b,
4f3f1b1, d4885ea, 0ec6541, b35d612, c51d6a0, 03c8639, 7ff70db, 8eb9b71,
5613ff6, 94c8757, f61433e.

**Unit-test verification:** `pytest tests/test_dk_singles_parser.py
tests/test_fd_singles_parser.py tests/test_dk_client.py
tests/test_fd_client.py` → 22/22 passed.

### Checklist results

**1. Data integrity**

- *No duplicate writes?* **YES.** Both `scraper_draftkings_singles.write_to_duckdb`
  and `scraper_fanduel_singles.write_to_duckdb` wrap the rewrite in a single
  `BEGIN TRANSACTION` / `DELETE FROM mlb_odds` / `executemany INSERT` /
  `COMMIT` with a try/except `ROLLBACK`. Partial-commit is impossible.
- *Proper deduplication?* **YES.** `mlb_odds` is fully rewritten each cycle
  (DELETE-then-INSERT inside a transaction). The parser's
  `parse_*_to_wide_rows` buckets by `(period, market_type, line)` so duplicate
  selections within a single scrape collapse to one row.
- *Incomplete records filtered?* **YES.**
  - DK: `fetch_event_selections` skips selections with `isDisabled` or empty
    `displayOdds.american`; `_parse_american_odds` failure → skip
    (dk_client.py:111–120).
  - FD: `fetch_event_runners` skips `runnerStatus != "ACTIVE"` (empty allowed),
    skips missing `selectionId`, skips missing american odds
    (fd_client.py:122–132).
  - Both parsers skip selections whose `market_id` isn't in `market_meta`
    (props/futures/team totals filtered out by `classify_market` →
    market never reaches `market_meta` map).
  - DK alt-spread-shaped selections whose name doesn't match either team are
    skipped (singles parser :191–193). Same for FD :183–185.

**2. Resource safety**

- *All DB connections use try/finally?* **YES.**
  - Python: both `write_to_duckdb` functions use `con = duckdb.connect(...)`
    followed by `try: ... finally: con.close()`
    (scraper_draftkings_singles.py:305–353; scraper_fanduel_singles.py:257–305).
  - R: `get_dk_odds` and `get_fd_odds` both use
    `on.exit(dbDisconnect(con, shutdown = TRUE))` immediately after
    `dbConnect` (Tools.R additions). Mirrors the pattern of the surrounding
    `get_bookmaker_odds`/`get_bet105_odds` helpers.
  - MLB.R loaders use `tryCatch(...)` to demote any scraper read to an empty
    tibble — but DB disconnection is owned by `get_dk_odds`/`get_fd_odds` via
    `on.exit`, so an exception inside the helper still closes the connection.
- *No lock-file leaks?* **YES.** Inspected `dk_odds/` and `fd_odds/` after
  recent scrapes — only `dk.duckdb`/`fd.duckdb` + `README.md`. No stale
  `.wal`/`.lock`/`*-journal` artifacts. (Note: `.gitignore` excludes
  `*.duckdb` and `*.duckdb-journal`, so even if WAL files transiently exist
  during a scrape they won't be committed.)

**3. Edge cases**

- *Off-season behavior?* **YES.**
  - DK: `client.list_events()` returns `[]` → loop body skipped → `all_rows`
    empty → `write_to_duckdb([])` still creates the table if missing and
    issues `DELETE FROM mlb_odds` (no-op on empty table). Resulting `mlb_odds`
    is empty.
  - FD: same flow.
  - MLB.R: `get_dk_odds("mlb")` finds an empty table → emits warning and
    returns `data.frame()` → `map_scraper_markets_mlb()` on an empty frame is
    a no-op → `scraper_to_canonical(dk_odds, ...)` returns NULL on empty input
    → `Filter(Negate(is.null), book_odds_by_book)` drops it → DK pills empty
    on the dashboard. Same for FD. **Graceful degradation matches the design
    spec's failure table.**
- *First-run with no existing DuckDB?* **YES.** `write_to_duckdb` runs
  `db_path.parent.mkdir(exist_ok=True)` and `CREATE TABLE IF NOT EXISTS
  mlb_odds (...)` before the DELETE. `duckdb.connect(...)` creates the file
  if missing. MLB.R-side `get_dk_odds` checks `file.exists(db_path)` and
  short-circuits to `data.frame()` with a warning if the DB hasn't been
  produced yet (Tools.R additions).
- *Timezone boundaries?* **YES.** `fetch_time = datetime.utcnow()` stored as
  naive UTC TIMESTAMP. Singles scrapers don't filter by tipoff — they snap
  whatever DK/FD post; downstream tipoff filtering is `MLB.R`'s
  responsibility via the `pt_start_time > Sys.time()` filter on
  `all_bets_combined`. No TZ ambiguity in the scraper itself.
- *DK Unicode-minus handling?* **YES.** `_parse_american_odds`
  (dk_client.py:140–155) normalizes `"−"` (U+2212) → `"-"` before `int(...)`.
  Covered by `test_dk_client.py`.

**4. Dead code**

- *No unused imports / flags?* **YES.** All imports in
  `scraper_draftkings_singles.py` (`argparse`, `re`, `datetime`, `Path`,
  `Any`, `duckdb`, `DraftKingsClient`/`Event`/`Selection`) are used.
  `scraper_fanduel_singles.py` similarly. Both `dk_client.py` and
  `fd_client.py` import only what they need (cleanup commit f1c074b removed
  an unused typing import from dk_client).
- *No dead branches?* **YES.** Reviewed `classify_market` (DK keyword-blacklist
  + period detection + market-type detection; every branch is reachable from
  the test fixtures) and FD's whitelist (10 entries — checked against
  `fd_events.json` fixture).
- **Acceptable duplication noted by implementer:** `get_dk_odds` and
  `get_fd_odds` in Tools.R are near-identical (~120 lines each, differ only in
  `bookmaker_key`, `dk_game_id`/`fd_game_id`, `db_path`). They mirror the
  established `get_bookmaker_odds` / `get_bet105_odds` pattern (also
  near-duplicates) so the redundancy is consistent with the codebase. **Not
  blocking — flagged as follow-up refactor.**

**5. Log/disk hygiene**

- *Log rotation?* **N/A.** Singles scrapers log to stdout via `print(...,
  flush=True)`. The `run.py` orchestrator captures stdout into its own
  per-cycle log; no new file-based logging introduced.
- *Unbounded file growth?* **NO** (i.e. growth is bounded). `dk.duckdb` and
  `fd.duckdb` are atomically rewritten each cycle (DELETE-then-INSERT inside
  a transaction). Observed size ~525 KB per DB; bounded by per-cycle row
  count (~1,300 rows estimated in design spec, 409–910 measured in commit
  messages).

**6. Security**

- *Secrets in logs?* **NO.** Every `print(...)` statement in the singles
  scrapers prints only event counts, event IDs, exception messages, or row
  counts — no cookies, headers, session tokens. `verbose` mode just adds
  per-event row counts and a parlays-enrich exception (also limited to the
  exception's `str(e)`).
- *API keys exposed?* **N/A.** Neither DK nor FD use authenticated REST
  endpoints. Both rely on curl_cffi Chrome-TLS impersonation; no API key,
  bearer token, or username/password is ever in the codebase or logs.
- *New .env / .pem / credentials files?* **NO.** Only new files are Python
  source, JSON fixtures (sanitized public DK/FD payloads), CSV golden
  baselines (4-column odds data, no PII), READMEs, R helpers, and CLAUDE.md
  edits. Inspected `git diff --stat 1ad39bb^..HEAD` — no secret files.

### Issues to fix before merge

- **None.**

### Acceptable risks (noted, not blocking)

- **Live-API fragility.** Both scrapers depend on DK's `parlays/v1/sgp/events`
  and FD's `event-page` endpoints. If DK changes its `displayOdds.american`
  shape (e.g. switches to a different Unicode minus, returns floats, or
  changes path), the parser silently emits zero rows for affected selections.
  Same for FD's `winRunnerOdds.americanDisplayOdds.americanOdds`. Mitigated by
  per-event isolation (one bad event doesn't tank the scrape) and by the
  dashboard's graceful-degrade-to-empty behavior — but operational
  monitoring of row counts on `dk.duckdb`/`fd.duckdb` is recommended.
- **`unmapped_teams` warning is print-only.** DK adds expansion teams rarely;
  if/when DK adds a 31st team or rebrands, the singles scraper will emit a
  WARNING line and ship rows with the un-canonicalized DK name. MLB.R-side
  `resolve_offshore_teams` will probably fail to match, and that game's DK
  pills will be empty. Self-healing once `DK_TEAM_MAP` is updated.
- **Parlays-enrich is best-effort.** DK singles scraper does an additional
  enrichment GET on the parlays endpoint after the markets subcategory call;
  failure is caught and only logged in `--verbose`. If this call goes
  permanently 4xx, F7 totals/alts silently disappear from the DK output. Not
  catastrophic (downstream just sees empty pills for those markets) but worth
  watching during the first few production cycles.

### Follow-up refactors (not blocking)

- **Tools.R `get_*_odds` helper duplication.** `get_dk_odds`, `get_fd_odds`,
  `get_bookmaker_odds`, `get_bet105_odds` share ~85% of their bodies. A
  single `get_scraper_odds(bookmaker_key, db_path, sport)` parameterized
  helper would shrink Tools.R by ~400 lines and centralize schema changes.
  Punted because the redundancy was already established before this PR.
- **Singles scraper main-loop duplication.** `scrape_singles` in
  `scraper_draftkings_singles.py` and `scraper_fanduel_singles.py` differ only
  in (a) the client class, (b) FD's lack of parlays-enrich, and (c) DK's
  team-name canonicalization step. Could collapse into a shared
  `_scrape_singles(client, classify_fn, canonicalize_fn=None)` helper. Punted
  to keep this PR minimal.
- **No integration test in CI.** The SGP regression test
  (`test_sgp_regression.py`) requires live DK/FD APIs and is `@pytest.mark.
  integration` — not run by default. A unit-test-only smoke that mocks the
  HTTP layer end-to-end (events list → markets → selections → write_to_duckdb
  → query) would catch parser regressions before they hit production.

### Verification commands run

- `pytest mlb_sgp/tests/test_dk_singles_parser.py
  mlb_sgp/tests/test_fd_singles_parser.py mlb_sgp/tests/test_dk_client.py
  mlb_sgp/tests/test_fd_client.py -q` → 22/22 passed (0.31s).
- `ls -la dk_odds/ fd_odds/` → only `*.duckdb` + `README.md`; no stale WAL.
- `git diff --stat 1ad39bb^..HEAD -- ':!Answer Keys/MLB Dashboard/'
  ':!docs/superpowers/specs/2026-05-11-mlb-odds-screen-*'` → 34 files,
  +69,206 −11 (bulk is fixtures: `dk_event_markets.json` 49k lines,
  `dk_league_response.json` 9k lines).

**Verdict: cleared for merge to `main` once user explicitly approves.**

---

## Post-smoke addendum (2026-05-13)

The first visual smoke session against real bets-tab data surfaced **seven**
latent bugs. Five were fixed in the feature branch (commits `8889b74`,
`70537ad`, `e09ddec`, `5566784`, `f3d6875`); the remaining two are fixed in
this follow-up session. Root-cause traces for all seven are documented in
`docs/superpowers/specs/2026-05-13-dashboard-card-layout-bugs.md`. Headline:

| # | Bug | Status | Commit / file |
|---|-----|--------|---------------|
| 1 | Mixed-type NA columns crash reactable (cents/placed_actual/fill_diff) | ✅ Fixed | `8889b74`, `70537ad` |
| 2 | Bets-tab book filter reads wrong cell text | ✅ Fixed | `e09ddec` |
| 3 | DK/FD canonical `market_name` missed `_1st_X_innings` suffix | ✅ Fixed | `5566784` |
| 4 | CSS `.cell-pickside` missing `!important` → vertical pill stacking | ✅ Fixed | `f3d6875` |
| 5 | `LINE_MATCH_TOLERANCE = 1.0` too strict for F-period totals display | ✅ Fixed (raised to 3.0) | `f3d6875` |
| 6 | DK F5/F7 multi-line markets bundled into single bucket, last-write-wins | ✅ Fixed (this session) | `scraper_draftkings_singles.py` |
| 7 | Spread-bet pill line tag rendered as `O-1.5` instead of `-1.5` | ✅ Fixed (this session) | `book_pill.R`, `mlb_dashboard.R` |

### Bug #6 root cause and fix

DraftKings posts each F5/F7 run-line and total as a **single market** that
contains the main line *and* all alt lines side by side. Example payload from
`parlays/v1/sgp/events/<id>` on 2026-05-13:

```
'Run Line - 1st 7 Innings'   id=2_84775087  selections=8  lines={-1.5,-0.5,0.5,1.5}
'Total Runs - 1st 7 Innings' id=3_84775084  selections=6  lines={5.5, 6.5, 7.5}
```

`classify_market` correctly labeled both as `(F7, "main")`, but
`parse_selections_to_wide_rows` keyed its bucket by `(period, market_type,
None)` for all `"main"` rows — so all 6 over/under selections (and all 8
spread selections) collapsed into one row per period, and **only the last
selection written survived**. The DB therefore stored 15 main F7 rows total
across 15 games (one line each) when DK was actually posting 3 total lines
and 2 spread lines per game.

Fix: pre-pass over the selection list counts distinct totals lines (Over/
Under) and distinct spread `|line|` values per `market_id`. Any
`"main"`-classified market carrying more than one distinct line is
reclassified per selection as `alternate_totals` / `alternate_spreads`, so
each line gets its own bucketed row. FG markets are unaffected — DK already
splits FG main from FG alt into separate markets (`Run Line` vs
`Run Line Alternate`, `Total` vs `Total Alternate`), so the pre-pass count
stays at 1 and the existing `"main"` coalescing into one ML+RL+Total row
per game is preserved.

**Verification:**
- Worktree `dk.duckdb` after re-scrape: 30 F5 alt-spread rows (was 0), 45 F5
  alt-total rows (was 0); same for F7 (was 15 total → now 75 split by line).
  FG row counts unchanged (150 alt-spread, 305 alt-total, 15 main).
- The `_1st_5_innings`/`_1st_7_innings` ML market (DK names it bare, e.g.
  `1st 5 Innings`) is still **not captured** because `classify_market` only
  matches names containing one of `run line` / `moneyline` / `total`. This is
  a pre-existing gap (predates the singles-scraper feature) — flagged as
  follow-up, not blocking.

### Bug #7 root cause and fix

`render_book_pill()` rendered the mismatched-line tag with a hard-coded
`prefix <- if (side == "under") "U" else "O"`, regardless of whether the
underlying bet was a totals market. Spread bets pass `side_word = "over"`
by default (the cascade in `create_bets_table` only checks `bet_on` for
"Over" / "Under" prefix), so a Yankees -1.5 pill showed `O-1.5` as the
mismatched-line tag.

Fix: added an `is_totals` parameter (default `TRUE` for backwards
compatibility with existing call sites) to `render_book_pill`. When
`is_totals = FALSE`, the mismatched-line tag uses a new `signed = TRUE`
mode on `.format_line_value()` to render `-1.5` / `+0.5` directly, with no
O/U prefix. Plumbed through `render_side_row` and both `create_bets_table`
call sites by computing `is_totals_market <- grepl("^totals",
table_data$market[i])`.

**Verification:** all 24 unit tests in `tests/test_book_pill.R` pass,
including 3 new tests for the spread-tag case:
- spread mismatch on favorite side → tag shows `-1.5`, no `O` prefix
- spread mismatch on dog side → tag shows `+1.5`, no `U` prefix
- spread exact-line → no tag rendered (unchanged from totals exact-line)

### Pre-merge checklist (re-run after follow-up)

- **Data integrity.** DK `dk.duckdb` now correctly stores N rows per
  (game, period, market_type) when DK bundles multiple lines in one market.
  Atomic DELETE-then-INSERT transaction is unchanged. FG behaviour and FD
  behaviour are bit-identical to pre-fix output (FD already produced per-line
  rows because FD splits markets cleanly by name).
- **Resource safety.** No new connections, no new lock-file paths. The
  parse_selections_to_wide_rows change is purely in-memory.
- **Edge cases.** New tests cover: bundled multi-line "main" markets (DK
  F5/F7), spread mismatch on both sides, spread exact-line. Pre-existing FG
  tests still pass.
- **Dead code.** None introduced. The pre-pass and the
  `effective_market_type` variable are both used on every iteration.
- **Backwards compatibility.** `render_book_pill`'s new `is_totals`
  defaults to `TRUE`, so any existing caller that did not pass it keeps
  emitting the O/U prefix. All test cases that don't pass `is_totals` still
  exercise the totals path.

### Known gaps (not blocking — pre-existing, separate workstream)

- **DK F-period moneyline (`1st N Innings` market) still uncaptured.**
  `classify_market` requires a market name containing `run line` /
  `moneyline` / `total`; DK names its F5/F7 ML markets just `1st 5 Innings`
  / `1st 7 Innings`. Result: F5/F7 ML bets have no DK pill on the dashboard.
  This predates the singles-scraper feature (the previous Odds-API source
  also lacked DK F-period ML). To fix: extend `classify_market` to detect a
  bare period name with 2 no-line selections as `(period, "main")` ML; add a
  fixture + test. Out of scope for this merge.
- **Worktree DB path mismatch during testing.** The Python scrapers write to
  `<repo>/dk_odds/dk.duckdb` and `<repo>/fd_odds/fd.duckdb` (path is relative
  to the script via `Path(__file__).parent.parent`), so when a scraper is run
  from a `git worktree`, the DBs land in the worktree. But MLB.R has
  `setwd("~/NFLWork/Answer Keys")` hardcoded (line 6 of `MLB.R`), and
  `get_dk_odds()` / `get_fd_odds()` default `db_path = "~/NFLWork/dk_odds/
  dk.duckdb"` — both forcing main repo paths. Testing in a worktree therefore
  requires copying the per-book DBs from `worktree/dk_odds/dk.duckdb` →
  `~/NFLWork/dk_odds/dk.duckdb` (and same for FD) before re-running the
  pipeline. The CLAUDE.md guidance ("NEVER symlink DuckDB databases — always
  copy") covers this; flagging here as a pre-merge testing footgun. Not a
  code defect.

**Updated verdict: cleared for merge to `main` after user explicitly
approves.** Two follow-ups noted but not blocking this PR.
