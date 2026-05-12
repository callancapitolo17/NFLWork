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
