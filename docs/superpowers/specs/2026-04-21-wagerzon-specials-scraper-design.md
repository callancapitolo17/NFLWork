# Wagerzon Specials Scraper + Generic Prop Pricer — Design

**Status:** Draft for review
**Date:** 2026-04-21
**Owner:** Callan
**Related:** `docs/superpowers/plans/2026-04-21-mlb-triple-play-pricer.md` (the manual-tribble pricer this replaces as the input)

## Goal

Automate the input side of the MLB triple-play / grand-slam pricer. Today's `mlb_triple_play.R` uses a hand-typed `tribble` of 8 lines — every day, the user re-types the posted odds. This design replaces that manual step with a scheduled Wagerzon scraper that captures **all MLB specials** into DuckDB, plus a generic description parser that turns verbatim prop labels (e.g. `"GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2½)"`) into structured leg specs the pricer can evaluate against dispersion-matched samples.

Out of scope for this design: pricing logic for prop types other than triple-play and grand slam; automated bet placement; dashboard integration; CLV tracking (the append-only schema supports it, but wiring is future work).

## Architecture

```
run.py mlb  (existing orchestrator, unchanged interface)
  ├── [parallel] wagerzon scraper (existing) → wagerzon_odds (existing table)
  ├── [parallel] wagerzon_specials scraper (NEW) → wagerzon_specials (NEW table)
  ├── [parallel] hoop88 / bfa / bookmaker / bet105 / kalshi scrapers (existing)
  ├── [parallel] MLB.R pipeline (unchanged)
  └── Sentinel: .scrapers_done_mlb

mlb_triple_play.R  (modified — no longer uses hardcoded tribble)
  ├── Reads mlb_game_samples (as today)
  ├── Reads wagerzon_specials WHERE prop_type IN ('TRIPLE-PLAY','GRAND-SLAM')
  │     AND scraped_at = most-recent, filtered to games in next 12h
  ├── parse_legs(description) → structured leg list per row
  ├── compute_prop_fair(samples, side, legs) → joint probability
  └── Prints fair vs book + edge per posted line
```

Authentication reuses the existing `wagerzon_odds/.wagerzon_profile` Playwright profile — same pattern as `recon_wagerzon.py`, `scraper_wagerzon.py`, and `parlay_pricer.py`. No new auth infrastructure.

## Storage schema

One new table `wagerzon_specials` in the existing `wagerzon_odds/wagerzon.duckdb`:

```sql
CREATE TABLE wagerzon_specials (
  scraped_at        TIMESTAMP,          -- when this snapshot was captured
  sport             VARCHAR,            -- 'mlb' for now; future 'nba', 'nfl', etc.
  league_id         INTEGER,            -- Wagerzon's lg param (4899 for MLB specials)
  game_id           INTEGER,            -- Wagerzon's idgm (internal game ID)
  home_team         VARCHAR,            -- Wagerzon's home team name (raw)
  away_team         VARCHAR,            -- Wagerzon's away team name (raw)
  game_date         DATE,
  game_time         TIMESTAMP,
  rotation_number   INTEGER,            -- the "bet number" (e.g. 757041)
  prop_type         VARCHAR,            -- Wagerzon's section header verbatim
                                        --   (e.g. 'TRIPLE-PLAY', 'GRAND-SLAM')
  description       VARCHAR,            -- row text verbatim
                                        --   (e.g. 'GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)')
  team              VARCHAR,            -- best-effort parse from description, else NULL
  line              DOUBLE,             -- points value, NULL for moneyline-style props
  odds              INTEGER             -- American odds
);

CREATE INDEX ON wagerzon_specials (scraped_at);
CREATE INDEX ON wagerzon_specials (sport, prop_type);
CREATE INDEX ON wagerzon_specials (rotation_number);
```

Design choices:
- **Append-only, keyed on `scraped_at`.** History preserved for CLV analysis; pricer queries `scraped_at = (SELECT MAX(...))`.
- **`prop_type` is Wagerzon's raw section header.** No classification logic in the scraper — Wagerzon already tagged it. New prop types become available by adding pricing logic (pricer-side), not scraper changes.
- **`description` is verbatim.** All prop-specific info (lines, periods, side) is derivable from this string at pricer-time.
- **Team names stay raw.** Downstream `resolve_offshore_teams()` (in `Tools.R`) handles canonical mapping, same pattern as every other Wagerzon consumer.

This is a different shape from the project's 18-column scraper schema (which is spread/total/h2h-specific). Acceptable because props don't fit that shape — trying to force them creates mostly-NULL columns and awkward conventions.

## Scraper

### Recon-first

Before writing production code, run `wagerzon_odds/recon_specials.py` — a small Playwright script modeled on the existing `recon_wagerzon.py`. It:
1. Launches Chromium with `.wagerzon_profile` (reuses existing auth)
2. Navigates to the URL the user manually opens when they view specials (captured via browser DevTools Network tab — user confirms the URL; likely `NewSchedule.aspx?WT=0&lg=4899` or similar)
3. Captures all network traffic for that URL
4. Dumps results to `wagerzon_odds/recon_specials.json`: URL structure, response format, section header HTML, row format

Output of recon determines scraper style:
- **If JSON endpoint exists** — scraper is a REST call in a `requests` session seeded with profile cookies, same pattern as `parlay_pricer.py`
- **If HTML-only** — scraper uses Playwright to render, then BeautifulSoup to parse

Either way, the `wagerzon_specials` table schema is unchanged.

### Production scraper

New file: `wagerzon_odds/scraper_specials.py`

Responsibilities:
1. Authenticate (via existing profile)
2. Fetch specials page for each configured league (`4899` for MLB; future: NBA/NFL specials leagues)
3. Parse each section: header text → `prop_type`, each row → `(rotation_number, description, line, odds)` + team best-effort parse
4. Join game info (home_team, away_team, game_date, game_time) — this is either on the same page or via a supplementary call to the existing schedule endpoint
5. Write rows to `wagerzon_specials` with a single `scraped_at` timestamp

CLI:
```bash
python3 wagerzon_odds/scraper_specials.py --sport mlb
```

Integrates into `run.py mlb` as a parallel scraper step alongside existing ones.

### Error handling

- **Auth expired** — detect login redirect, exit with code 2 and a clear message; user re-runs the existing login flow
- **No specials posted** (off-day, pre-season) — write zero rows, exit 0. Pricer gracefully handles empty
- **Page format changed** — scraper raises loudly rather than writing corrupt rows. We diagnose via the `debug_wagerzon_page.html` convention already used elsewhere
- **Partial data** — a row missing rotation_number or odds is skipped with a warning log, the rest of the scrape proceeds

## Description parser

New file: `Answer Keys/parse_legs.R`

`parse_legs(description)` returns a list of structured leg specs. Example outputs:

```r
parse_legs("GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)")
# → list(
#     list(type = "scores_first"),
#     list(type = "wins_period", period = "F5"),
#     list(type = "wins_period", period = "FG")
#   )

parse_legs("GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2½)")
# → list(
#     list(type = "scores_first"),
#     list(type = "wins_period", period = "F5"),
#     list(type = "wins_period", period = "FG"),
#     list(type = "team_total_under", line = 2.5)
#   )
```

Implementation:
1. Extract parenthetical: match `\((.+)\)`, capture the inside (e.g. `"SCR 1ST, 1H, GM & SCR U2½"`)
2. Tokenize: split on `,` and `&`, trim whitespace
3. Map each token via a **token registry** → leg spec

### Token registry

```r
# Each entry: regex → function(matches) returning leg spec list
TOKEN_REGISTRY <- list(
  # Period-lead tokens (F3, 1H=F5, F7, GM=FG)
  list(pattern = "^F3$",      spec = function(m) list(type="wins_period", period="F3")),
  list(pattern = "^1H$",      spec = function(m) list(type="wins_period", period="F5")),
  list(pattern = "^F7$",      spec = function(m) list(type="wins_period", period="F7")),
  list(pattern = "^GM$",      spec = function(m) list(type="wins_period", period="FG")),
  # First-to-score
  list(pattern = "^SCR 1ST$", spec = function(m) list(type="scores_first")),
  # Team-total under / over (matches "SCR U2½", "SCR U2.5", "SCR U3", etc.)
  list(pattern = "^SCR U(\\d+(?:[.½]\\d*)?)$",
       spec = function(m) list(type="team_total_under",
                               line = parse_fraction(m[[1]]))),
  list(pattern = "^SCR O(\\d+(?:[.½]\\d*)?)$",
       spec = function(m) list(type="team_total_over",
                               line = parse_fraction(m[[1]]))),
  # Reserved for future: opponent totals (opp_total_under/over). Not used in v1.
)
```

**`parse_fraction()`** handles Wagerzon's mixed notation — `"2½"`, `"2.5"`, `"3"` all map to the numeric line. Plus the ½ unicode character (U+00BD).

**Unknown tokens** raise a warning via `parse_legs()` and return `NULL` for the leg list. The pricer treats a NULL leg list as "cannot price" and emits `NA` for `fair_prob`. The row still appears in the output table with `fair_prob = NA` so we see the prop exists but can't price it yet — a clear signal to extend the registry.

### Semantic choices confirmed

- **`SCR U<N>` = listed team's team-total under N** (user confirmation 2026-04-21). The parser emits `team_total_under` — where "team" means the prop's subject team (the team named in the description, e.g., Giants in "GIANTS GRAND-SLAM"). The pricer's `side` parameter tells the evaluator which side's runs to compare.
- **Period naming** uses existing answer-key convention: `F3`, `F5`, `F7`, `FG` — matches the column suffixes in `mlb_game_samples`.

## Generic pricer

Replace `compute_triple_play_fair()` and `compute_grand_slam_fair()` with a single leg-driven function in `Answer Keys/mlb_triple_play.R`:

```r
compute_prop_fair <- function(samples, side, legs) {
  if (is.null(legs) || length(legs) == 0) return(NA_real_)
  samples <- samples[!is.na(samples$home_scored_first), ]
  if (nrow(samples) == 0) return(NA_real_)

  # Derive team-specific run totals once per sample set
  home_runs <- (samples$total_final_score + samples$home_margin) / 2
  away_runs <- (samples$total_final_score - samples$home_margin) / 2
  team_runs <- if (side == "home") home_runs else away_runs
  opp_runs  <- if (side == "home") away_runs else home_runs

  # AND-reduce across legs
  hits <- rep(TRUE, nrow(samples))
  for (leg in legs) {
    hits <- hits & eval_leg(leg, samples, side, team_runs, opp_runs)
  }
  mean(hits)
}

eval_leg <- function(leg, samples, side, team_runs, opp_runs) {
  switch(leg$type,
    scores_first     = (samples$home_scored_first == if (side == "home") 1L else 0L),
    wins_period      = {
      col <- switch(leg$period,
                    F3 = samples$home_margin_f3,
                    F5 = samples$home_margin_f5,
                    F7 = samples$home_margin_f7,
                    FG = samples$home_margin)
      if (side == "home") col > 0 else col < 0
    },
    team_total_under = team_runs <  leg$line,
    team_total_over  = team_runs >  leg$line,
    opp_total_under  = opp_runs  <  leg$line,  # reserved, unused in v1
    opp_total_over   = opp_runs  >  leg$line,  # reserved, unused in v1
    stop("Unknown leg type: ", leg$type)
  )
}
```

**Sample-column requirement**: `compute_prop_fair` needs `total_final_score` in the samples frame. The `mlb_game_samples` table already stores this; the pricer just needs to add it to the `SELECT` list (currently selects only `game_id, home_margin, home_margin_f5, home_scored_first`). Add `total_final_score, home_margin_f3, home_margin_f7` so all period columns are available.

### Pricer wiring

The main block of `mlb_triple_play.R` becomes data-driven:

```r
# Read live lines from the scraper table instead of a hardcoded tribble
wz_con <- dbConnect(duckdb(),
                    dbdir = "~/NFLWork/wagerzon_odds/wagerzon.duckdb",
                    read_only = TRUE)
specials <- dbGetQuery(wz_con, "
  SELECT rotation_number, home_team AS wz_home, away_team AS wz_away,
         game_time, prop_type, description, team, line, odds
  FROM wagerzon_specials
  WHERE sport = 'mlb'
    AND prop_type IN ('TRIPLE-PLAY', 'GRAND-SLAM')
    AND scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials WHERE sport='mlb')
    AND game_time > NOW() AND game_time < NOW() + INTERVAL 12 HOUR
")
dbDisconnect(wz_con)

specials <- resolve_offshore_teams(specials, sport = "mlb")  # Tools.R helper

todays_lines <- specials %>%
  mutate(
    side = case_when(team == wz_home ~ "home",
                     team == wz_away ~ "away",
                     TRUE             ~ NA_character_),
    home_team = wz_home, away_team = wz_away,
    book_odds = odds, target_team = team,
    legs = lapply(description, parse_legs)
  ) %>%
  filter(!is.na(side))

# Join to consensus for game_id (same as current pricer)
matched <- todays_lines %>%
  inner_join(consensus, by = c("home_team", "away_team"))

priced <- matched %>% rowwise() %>% mutate(
  game_samples = list(samples_df[samples_df$game_id == id, ]),
  fair_prob    = compute_prop_fair(game_samples, side, legs),
  fair_odds    = prob_to_american(fair_prob),
  book_prob    = american_to_prob(book_odds),
  edge_pct     = (fair_prob / book_prob - 1) * 100
)
```

Output table adds `prop_type` column and groups triple-plays and grand slams together, sorted by edge desc.

## Testing

- **`tests/test_parse_legs.R`** — unit tests for every known pattern:
  - `"GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)"` → 3-leg list
  - `"GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2½)"` → 4-leg list
  - `"GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2.5)"` → same as above (fraction parsing)
  - `"GIANTS MYSTERY (XYZ)"` → NULL + warning (unknown token)
- **`tests/test_eval_leg.R`** — synthetic samples × every leg type:
  - `scores_first` home-side: only rows where `home_scored_first == 1L` pass
  - `wins_period` with F5 home: only rows where `home_margin_f5 > 0` pass
  - `team_total_under` with line=2.5 home: only rows where `home_runs < 2.5` pass
  - Each of the above, away-side mirrors
- **`tests/test_compute_prop_fair.R`** — end-to-end synthetic:
  - Triple-play legs on a fixture where home hits 3/10 rows → expect 0.3
  - Grand slam legs on a fixture where home hits 1/10 rows → expect 0.1
  - Empty legs → NA_real_
  - Empty samples → NA_real_
- **`tests/test_wagerzon_specials_schema.py`** — after first real scraper run, assert `wagerzon_specials` has the expected columns, all rows have `scraped_at` not null, odds are integers, description is non-empty
- **Integration smoke test** — once recon + scraper run once, assert `mlb_triple_play.R` prices at least one row with non-NA `fair_prob`

## Rollout

Single plan, two commits staged in one merge to main:

**Commit set A — Scraper + storage (foundations)**
1. `recon_specials.py` (one-time artifact; not production path)
2. `scraper_specials.py` with CLI + `run.py mlb` integration
3. DuckDB schema migration (create `wagerzon_specials`)
4. First successful scrape captured to DuckDB

**Commit set B — Parser + pricer (consumer)**
1. `parse_legs.R` with token registry + tests
2. `compute_prop_fair()` + `eval_leg()` in `mlb_triple_play.R` + tests
3. `mlb_triple_play.R` main block rewired to read DB instead of tribble
4. `Answer Keys/CLAUDE.md` updated with the new data flow

Ship tomorrow when new specials post. Manual verification before first production run:
1. Eye-check that scraped rows match the page visually (row count, odds values, prop types present)
2. Eye-check that `parse_legs()` produces sensible leg lists for 5 sampled rows
3. Compare first-run fair prices against what we got today via the manual tribble for the same prop types
4. Only then trust the pipeline for betting decisions

## Open items (resolve during implementation)

1. **Exact specials URL** — confirm via recon. Current guess: `backend.wagerzon.com/wager/NewSchedule.aspx?WT=0&lg=4899`.
2. **JSON vs HTML response** — determines scraper style (see Recon-first section).
3. **Game info source** — whether the specials page carries `home_team/away_team/game_time` directly, or we need a supplementary call to the schedule endpoint for each `idgm`.
4. **Team name parsing** — most description lines start with the team name in all caps (`GIANTS TRIPLE-PLAY`). A simple regex `^([A-Z ]+?) (TRIPLE-PLAY|GRAND-SLAM)` should cover v1. If Wagerzon uses nicknames/city variants inconsistently, we add a dictionary.
5. **½ unicode handling** — confirm `parse_fraction("U2½")` correctly interprets as 2.5. R regex needs `perl = TRUE` for the unicode char.
6. **`run.py` scheduling cadence** — matches whatever cadence the existing Wagerzon scraper uses (currently 15-min periodic per session memory). Confirm and match.

## Version control

- Branch: `feature/wagerzon-specials-scraper`
- Worktree: `.worktrees/wagerzon-specials` (created via `/worktree`)
- DuckDB: copy `wagerzon_odds/wagerzon.duckdb` into worktree (never symlink, per project rule)
- Never merge to main without explicit user approval
- Cleanup: remove worktree + delete branch after merge

## Documentation updates

- `Answer Keys/CLAUDE.md` — extend the "Triple-Play Data Flow" subsection with the new book-line input path (scraper → DB → parser → pricer)
- `wagerzon_odds/README.md` — add `scraper_specials.py` to the file inventory with a one-line description
- `nfl_draft/README.md` — N/A (this is MLB scope)
