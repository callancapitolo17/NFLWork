# DraftKings SGP Blend for MLB Trifectas — Design

**Status:** Draft for review
**Date:** 2026-04-22
**Owner:** Callan
**Related:** `docs/superpowers/specs/2026-04-21-wagerzon-specials-scraper-design.md` (the input side); `docs/superpowers/plans/2026-04-21-prop-parser-and-pricer.md` (the just-shipped parser this builds on)

## Goal

Add DraftKings SGP fair odds as a blended supporting data point in the MLB triple-play / grand-slam pricer, mirroring how `mlb_correlated_parlay.R` blends DK + FD SGP fair odds with the historical answer-key fair. The pricer's published `fair_odds` becomes the **mean of available probabilities** (model + DK), which feeds edge calculation and (eventually) Kelly sizing. Triangulating against an independent book reduces model risk and surfaces dangerous divergences (e.g. our model says +400 while DK has +900 — that's a signal to investigate before betting).

This is a parser-driven design: each leg type from `parse_legs()` (the just-shipped `Answer Keys/parse_legs.R`) maps to a DK selection-ID resolver. Triple-plays (3 legs) and grand slams (4 legs) are both priced by the same generic flow. Adding a new leg type to the parser later requires adding a corresponding resolver — no scraper architecture changes.

Out of scope: FanDuel SGP integration for trifectas (the existing FD SGP scraper handles correlated parlays only); Kelly sizing changes; dashboard integration; CLV tracking on DK SGP prices (the schema supports it, but no consumer is built).

## Two-plan staging

Recon must precede production code (we don't know DK's selection ID format for the new leg types). This design specifies **Plan #1** in full and **Plan #2** in outline:

- **Plan #1** ships the schema, the recon script, NA-safe R blend scaffolding, and unit tests. Fully executable without DK live data — the pricer continues to produce model-only fair odds, but is now structurally ready to consume DK data the moment Plan #2 lands. ~6 commits.
- **Plan #2** is written after recon completes (you run `recon_dk_trifecta.py`, we inspect `recon_dk_trifecta.json`). It populates the resolvers and ships the production scraper. ~5 commits.

Same pattern as the Wagerzon specials work: recon-first, design with that knowledge.

## Architecture

```
                  ┌─────────────────────────────────────────┐
                  │  Answer Keys/parse_legs.R (existing)    │
                  │  parse_legs(description) → list of legs │
                  └──────────────────┬──────────────────────┘
                                     │ (leg specs drive both branches)
                  ┌──────────────────┴──────────────────┐
                  │                                     │
                  ▼                                     ▼
    ┌──────────────────────────────┐    ┌────────────────────────────────────┐
    │  HISTORICAL FAIR (existing)  │    │  DK SGP FAIR (NEW)                 │
    │  compute_prop_fair(samples,  │    │  Plan #2: Python scraper invoked   │
    │    side, legs)               │    │  via system2() from R              │
    │                              │    │                                    │
    │  Returns model_fair_prob     │    │  Returns dk_fair_prob (devigged)   │
    │  using the answer-key        │    │  via DK_SGP_VIG_DEFAULT = 1.10     │
    │  dispersion-matched samples  │    │  Plan #1: NA stub (DK col always   │
    │                              │    │  NA until Plan #2 lands)           │
    └──────────────┬───────────────┘    └─────────────────┬──────────────────┘
                   │                                      │
                   └──────────────┬───────────────────────┘
                                  ▼
                      ┌──────────────────────────┐
                      │  blend(model, dk) →      │
                      │  mean of available probs │
                      │  (skip NA, never fail    │
                      │   if either is missing)  │
                      └────────────┬─────────────┘
                                   ▼
                      ┌──────────────────────────────┐
                      │  Output columns:             │
                      │   target_team, prop_type,    │
                      │   side, n_samples,           │
                      │   model_odds, dk_odds,       │
                      │   fair_odds, book_odds,      │
                      │   edge_pct                   │
                      └──────────────────────────────┘
```

### Where each component lives

- `mlb_sgp/recon_dk_trifecta.py` — NEW. Plan #1. Manually-run Playwright/curl_cffi recon script. Captures `parlays/v1/sgp/events/{id}` and (if SGP combinable) `calculateBets` responses for known triple-play and grand-slam scenarios.
- `mlb_sgp/scraper_draftkings_trifecta.py` — NEW. Plan #2. Production scraper invoked from R via `system2()`. Reads a JSON request file, fetches DK SGP odds via `calculateBets`, writes rows to `mlb_trifecta_sgp_odds`.
- `mlb_sgp/dk_leg_resolvers.py` — NEW. Plan #2. Maps each `leg_type` to a selection-ID resolver function. Mirrors the structure of `parse_legs.R::TOKEN_REGISTRY`.
- `Answer Keys/mlb_triple_play.R` — MODIFIED. Plan #1 adds the blend scaffolding (NA-safe, works with no DK data). Plan #2 wires the `system2()` call.
- `Answer Keys/mlb.duckdb` — gains a new `mlb_trifecta_sgp_odds` table.
- `Answer Keys/CLAUDE.md` — extends "Triple-Play Data Flow" section to document the DK blend.

### Authentication and network reuse

DK SGP endpoints are public (no login). Plan #2's production scraper reuses the existing `curl_cffi` + `impersonate="chrome"` pattern from `mlb_sgp/scraper_draftkings_sgp.py`, including the same shared HTTP utilities if they exist (or just lift the request layer). No new auth profile needed.

The recon script in Plan #1 may use Playwright (for full network capture) or curl_cffi (if we trust the existing pattern works). I'll lean toward curl_cffi for simplicity since we already know that bypass works for `calculateBets`.

## Storage schema

New table `mlb_trifecta_sgp_odds` in `Answer Keys/mlb.duckdb`:

```sql
CREATE TABLE mlb_trifecta_sgp_odds (
  fetch_time     TIMESTAMP,        -- when this snapshot was captured
  game_id        INTEGER,          -- DK event ID (also the game_id we pass via R)
  prop_type      VARCHAR,          -- 'TRIPLE-PLAY' or 'GRAND-SLAM' (verbatim from Wagerzon)
  side           VARCHAR,          -- 'home' or 'away'
  legs_json      VARCHAR,          -- parsed leg list serialized:
                                   --   '[{"type":"scores_first"},{"type":"wins_period","period":"F5"},...]'
  selection_ids  VARCHAR,          -- comma-separated DK selection IDs used (debug)
  sgp_decimal    DOUBLE,           -- DK's combined SGP decimal odds, NULL if any leg missing
  sgp_american   INTEGER,          -- same in American format, NULL if missing
  source         VARCHAR           -- 'draftkings_direct' (matches existing convention)
);

CREATE INDEX ON mlb_trifecta_sgp_odds (fetch_time);
CREATE INDEX ON mlb_trifecta_sgp_odds (game_id, prop_type, side);
```

**Design choices:**

- **Append-only, keyed on `fetch_time`.** Same pattern as `mlb_sgp_odds` and `wagerzon_specials`. History preserved for CLV / debug; pricer queries `WHERE fetch_time = (SELECT MAX(fetch_time) FROM mlb_trifecta_sgp_odds)`.
- **Dedicated table, NOT reusing `mlb_sgp_odds`.** Reasoning: trifecta SGPs have only 2 observations per game (home + away) — fundamentally different vig methodology than correlated parlay SGPs (which fit per-game vig from 4+ observations). Mixing them in one table forces every consumer to filter by combo prefix and branch on vig logic. Separate tables keep concerns separate. Also: `legs_json` and `prop_type` are first-class metadata for trifectas but not for correlated parlays.
- **`legs_json` is the source of truth for leg structure.** Lets us replay the exact leg list later if the parser changes or DK changes a market spec. Critical for debugging.
- **`selection_ids` is for debugging only.** Not consumed by any pricer; included so that when DK returns surprising odds we can verify exactly which selections were used.
- **NULL `sgp_decimal` is valid** — happens when DK doesn't offer one of the legs (e.g., F3-winner doesn't exist on DK) or rejects the SGP combination. The R blend logic treats this as "DK data not available for this prop" and falls back to model-only fair.

## Plan #1 — Components

### 1. Schema migration

A small idempotent migration: when `mlb_triple_play.R` runs, ensure the table exists. Use `CREATE TABLE IF NOT EXISTS` so re-running doesn't error. Place near the top of the main block:

```r
# Ensure DK trifecta SGP table exists (idempotent)
con_w <- dbConnect(duckdb(), dbdir = MLB_DB)
dbExecute(con_w, "
  CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds (
    fetch_time     TIMESTAMP,
    game_id        INTEGER,
    prop_type      VARCHAR,
    side           VARCHAR,
    legs_json      VARCHAR,
    selection_ids  VARCHAR,
    sgp_decimal    DOUBLE,
    sgp_american   INTEGER,
    source         VARCHAR
  )
")
dbDisconnect(con_w)
```

### 2. Recon script `mlb_sgp/recon_dk_trifecta.py`

Plain Python script using `curl_cffi.requests.Session(impersonate="chrome")`. Takes one or two known DK event IDs (hardcoded for v1 — user updates before running) and:

1. Fetches `parlays/v1/sgp/events/{event_id}` for each game
2. Saves the full JSON response to `mlb_sgp/recon_dk_trifecta_<event_id>.json`
3. Tries a `calculateBets` POST with placeholder selection IDs that look like the candidates we expect (e.g. searches the response for "First to Score" market, picks the home team's selection ID)
4. If `calculateBets` returns clean odds, records the SGP decimal in the recon output

**Output**: a directory of recon JSON files we hand-inspect to populate Plan #2's resolvers.

The recon script does not depend on any other Plan #1 code — fully standalone — so it can be developed and reviewed in isolation.

### 3. R-side blend scaffolding in `mlb_triple_play.R`

After Task 6's existing parse_legs + compute_prop_fair pipeline, insert blend logic. Plan #1's blend reads `mlb_trifecta_sgp_odds` (which is empty until Plan #2 ships), so DK columns will be all NA in Plan #1. The blend correctly falls back to model_only.

```r
DK_SGP_VIG_DEFAULT <- 1.10  # same fallback as mlb_correlated_parlay.R

# Read latest DK trifecta SGP odds (empty until Plan #2 populates the table)
dk_sgp <- tryCatch({
  dbGetQuery(con, "
    SELECT game_id, prop_type, side, sgp_decimal
    FROM mlb_trifecta_sgp_odds
    WHERE source = 'draftkings_direct'
      AND fetch_time = (SELECT MAX(fetch_time) FROM mlb_trifecta_sgp_odds)
  ")
}, error = function(e) data.frame(
  game_id = integer(0), prop_type = character(0),
  side = character(0), sgp_decimal = double(0)
))

priced <- matched %>%
  rowwise() %>%
  mutate(
    game_samples    = list(samples_df[samples_df$game_id == id, ]),
    n_samples       = nrow(game_samples),
    legs            = list(parse_legs(description)),
    model_fair_prob = compute_prop_fair(game_samples, side, legs),

    dk_match        = list(dk_sgp[dk_sgp$game_id   == id        &
                                  dk_sgp$prop_type == prop_type &
                                  dk_sgp$side      == side, ]),
    dk_decimal      = if (nrow(dk_match) > 0) dk_match$sgp_decimal[1] else NA_real_,

    # Single seam — extracted helper. Tested in test_dk_blend.R.
    fair_prob       = blend_dk_with_model(model_fair_prob, dk_decimal,
                                          DK_SGP_VIG_DEFAULT),
    dk_fair_prob    = if (!is.na(dk_decimal) && dk_decimal > 0) {
                        (1 / dk_decimal) / DK_SGP_VIG_DEFAULT
                      } else NA_real_,

    model_odds      = prob_to_american(model_fair_prob),
    dk_odds         = prob_to_american(dk_fair_prob),
    fair_odds       = prob_to_american(fair_prob),
    book_prob       = american_to_prob(book_odds),
    edge_pct        = (fair_prob / book_prob - 1) * 100,

    prop_type       = {
      known_props <- "(TRIPLE-PLAY|GRAND-SLAM)"
      m <- regmatches(description,
                      regexec(paste0("\\b", known_props, "\\b"), description))[[1]]
      if (length(m) >= 2) m[[2]] else NA_character_
    }
  ) %>%
  ungroup() %>%
  select(target_team, prop_type, side, n_samples,
         model_odds, dk_odds, fair_odds, book_odds, edge_pct) %>%
  arrange(desc(edge_pct))
```

**Why this is safe to ship in Plan #1 even though DK data isn't flowing yet:**
- `dk_sgp` is fetched defensively with `tryCatch` — empty data frame on missing table or schema mismatch
- `dk_match` joins return zero rows when `dk_sgp` is empty → `dk_fair_prob` = NA
- `probs_to_blend` becomes a 1-element vector containing only `model_fair_prob` → `fair_prob` = `model_fair_prob`
- `model_odds == fair_odds` for every row in Plan #1
- `dk_odds` is NA for every row in Plan #1 — the column exists but is empty

When Plan #2 lands and the scraper actually populates the table, the same code starts producing real DK columns and blended fair odds with no R-side changes.

### 4. Plan #2 stub (R-side)

In Plan #1, the `system2()` call to invoke the scraper is **commented out with a TODO marker**. Plan #2's first task is to uncomment and verify it works:

```r
# Plan #2 will activate the scraper invocation here. For now, the table
# stays empty and the blend gracefully falls back to model-only fair.
#
# trifecta_input <- todays_lines %>% ... write JSON ...
# system2(SGP_VENV_PYTHON,
#         args = c(file.path(SGP_DIR, "scraper_draftkings_trifecta.py"),
#                  "--input", request_path, "--db", MLB_DB))
```

This makes Plan #2's R-side change an unambiguous "uncomment + remove the marker comment" — easy to review, no risk of subtle behavior change.

### 5. Unit tests for the blend logic

The blend helper `blend_dk_with_model(model_prob, dk_decimal, vig)` is defined alongside `compute_prop_fair` in `Answer Keys/parse_legs.R`:

```r
#' Blend a model fair probability with a DK SGP devigged probability.
#'
#' If dk_decimal is NA / 0 / negative, returns model_prob alone.
#' If model_prob is NA but DK is valid, returns DK alone.
#' If both NA, returns NA_real_.
#' Otherwise returns the mean of the two devigged probabilities.
#'
#' @param model_prob Numeric, the historical fair probability from compute_prop_fair.
#' @param dk_decimal Numeric, DK SGP decimal odds (>1.0 if valid).
#' @param vig Numeric, the SGP vig fallback (e.g. 1.10).
blend_dk_with_model <- function(model_prob, dk_decimal, vig) {
  dk_valid <- !is.na(dk_decimal) && dk_decimal > 0
  dk_prob  <- if (dk_valid) (1 / dk_decimal) / vig else NA_real_
  probs    <- c(model_prob, dk_prob)
  probs    <- probs[!is.na(probs)]
  if (length(probs) == 0) return(NA_real_)
  mean(probs)
}
```

New tests in `Answer Keys/tests/test_dk_blend.R`:

```r
test_that("blend returns model_only when DK is NA", {
  expect_equal(blend_dk_with_model(0.30, NA_real_, 1.10), 0.30)
})

test_that("blend averages model and DK devigged when both present", {
  # DK SGP +200 → decimal 3.0 → implied 0.333. Devigged at 1.10 → 0.303
  # Model 0.30. Blend = mean(0.30, 0.303) = 0.3015
  expect_equal(blend_dk_with_model(0.30, 3.0, 1.10), (0.30 + 0.30303) / 2,
               tolerance = 1e-4)
})

test_that("blend treats DK 0 / negative as missing", {
  expect_equal(blend_dk_with_model(0.30, 0,    1.10), 0.30)
  expect_equal(blend_dk_with_model(0.30, -1.0, 1.10), 0.30)
})

test_that("blend handles model NA gracefully (returns DK alone)", {
  # Edge case: if model returns NA (e.g., empty samples), DK is the only signal
  expect_equal(blend_dk_with_model(NA_real_, 3.0, 1.10), 0.30303,
               tolerance = 1e-4)
})

test_that("blend is NA when both inputs are NA", {
  expect_true(is.na(blend_dk_with_model(NA_real_, NA_real_, 1.10)))
})
```

The `blend_dk_with_model` function gets defined alongside `compute_prop_fair` in `Answer Keys/parse_legs.R` (or a new sibling `dk_blend.R` if we want stronger separation).

### 6. Documentation

`Answer Keys/CLAUDE.md` "Triple-Play Data Flow" subsection extends with the DK blend layer. New text appended after the existing diagram:

```markdown
**DK SGP blend (Plan #1 scaffolding, Plan #2 wires production):**
- mlb_trifecta_sgp_odds table holds DK SGP odds per (game_id, prop_type, side)
- Pricer reads latest fetch_time, devigs DK with DK_SGP_VIG_DEFAULT=1.10
- blend_dk_with_model(model, dk, vig) → mean of available probs
- fair_odds (the published number, used for edge calc) is the blend
- model_odds + dk_odds shown alongside for transparency
```

## Plan #2 — Outline (written after recon)

Plan #2 is roughly 5 tasks, all dependent on `recon_dk_trifecta.json` content:

1. `mlb_sgp/dk_leg_resolvers.py` — populate `LEG_RESOLVERS` dict with one resolver function per leg type, derived from recon findings. Each resolver takes `(leg_spec, side, event_state)` and returns a DK selection ID (or None if the market is missing).
2. `mlb_sgp/scraper_draftkings_trifecta.py` — production scraper. CLI: `--input <json>` (list of trifecta requests) and `--db <path>`. For each request: resolve legs, POST to `calculateBets`, write row to `mlb_trifecta_sgp_odds`.
3. Resolver unit tests against captured recon fixtures. One test per leg type, one round-trip test for each prop type.
4. R-side wiring — uncomment the `system2()` call in `mlb_triple_play.R` and remove the Plan #2 stub marker.
5. End-to-end integration smoke test on a known game.

Plan #2 will be written when recon completes. The scope and acceptable cost should match Plan #1 — single merge, ~5 commits.

## Testing summary

| Layer | Plan #1 | Plan #2 |
|---|---|---|
| Schema | `CREATE IF NOT EXISTS` runs idempotently | (no change) |
| Recon | Syntax check + import smoke | (no change — recon already ran) |
| Resolvers (Python) | (no resolvers yet) | Mock event JSON → expected selection IDs |
| Scraper (Python) | (no scraper yet) | Live `calculateBets` round-trip on one game |
| Blend (R) | 5 unit tests on `blend_dk_with_model` synthetic inputs | (covered by Plan #1) |
| Pricer end-to-end (R) | Existing 8-row regression preserved (model-only) | New regression: `dk_odds` populated for at least one row |

## Error handling — combined view

| Failure mode | Plan #1 behavior | Plan #2 behavior |
|---|---|---|
| Table missing or schema mismatch | `tryCatch` returns empty df → blend = model_only | (same) |
| DK row exists with NULL `sgp_decimal` | dk_fair_prob = NA → blend = model_only | (same — set by scraper) |
| DK row missing for a (game,prop,side) | dk_match has 0 rows → blend = model_only | (same) |
| DK scraper returns nonzero exit | (no scraper invocation in Plan #1) | R logs warning, blend = model_only |
| DK scraper throws | (no scraper invocation in Plan #1) | R catches via `tryCatch`, logs, blend = model_only |
| DK selection ID for one leg returns None (market missing) | (no resolvers in Plan #1) | Scraper writes NULL `sgp_decimal` for that prop, R blend = model_only |
| Akamai blocks the SGP endpoint | (no requests in Plan #1) | curl_cffi handles; if it ever fails, exit 2, R logs |

The pricer is **always** functional with model-only fair — DK is purely additive.

## Open items (resolve during Plan #2 implementation)

1. **DK selection ID format for new leg types** — first-to-score, F5 winner, F3/F7 winners (do they exist on DK at all?). Recon answers all of these.
2. **3-way vs 2-way F5 winner on DK** — affects whether tie probability needs to be modeled. If 2-way (no tie), DK pays differently than our F5-strict-lead semantics; need to either match DK's semantics or add a vig adjustment.
3. **Cross-market SGP restrictions** — DK is known to block alt-total + main-total combinations; whether these specific 3-leg or 4-leg combos are allowed is unknown until recon.
4. **Vig fallback accuracy** — `DK_SGP_VIG_DEFAULT = 1.10` is the correlated-parlay fallback; trifecta vig may differ systematically. Worth validating after Plan #2 collects a week of data: compute observed vig on a few known games, see how 1.10 holds up.
5. **DK event ID ↔ Wagerzon idgm mapping** — both reference "the same MLB game" but with different keys. Plan #2 needs a mapping (likely via team names + commence_time + canonical_match.py). Plan #1 punts this since no scraper invocation happens.

## Version control plan

**Plan #1:**
- Branch: `feature/dk-trifecta-blend-scaffolding`
- Worktree: `.worktrees/dk-trifecta-scaffolding`
- DuckDB: copy `mlb.duckdb` into worktree for verification (per project rule), remove after
- Commits: ~6 (one per component listed above)
- Documentation updated in same merge

**Plan #2:**
- Branch: `feature/dk-trifecta-scraper`
- Worktree: `.worktrees/dk-trifecta-scraper`
- DuckDB: same copy/cleanup pattern
- Commits: ~5

Each plan merges to main after explicit user approval — never auto-merge.

## Documentation updates

**Plan #1:**
- `Answer Keys/CLAUDE.md` — extend "Triple-Play Data Flow" with DK blend section
- `mlb_sgp/README.md` — add `recon_dk_trifecta.py` to file inventory (if README exists; create if not)

**Plan #2:**
- `Answer Keys/CLAUDE.md` — replace "Plan #1 scaffolding" hedge with the live data flow
- `mlb_sgp/README.md` — add `scraper_draftkings_trifecta.py` and `dk_leg_resolvers.py`
