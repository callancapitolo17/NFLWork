# MLB Dual-Source Edge Flagging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Flag a bet on the MLB dashboard bets tab when EITHER the model sees ≥2% EV OR a sportsbook's price is ≥2% out of line with the leave-one-out market consensus, with a MODEL / MARKET / BOTH badge on each card.

**Architecture:** A new isolated `find_market_edges()` module reads the canonical per-book odds `MLB.R` already builds, devigs each book's two sides, and flags any book whose price beats the median devigged fair of the *other* books (leave-one-out). Its output is merged into `mlb_bets_combined` via a shared canonical `bet_row_id`, so the existing book-prices expansion and dashboard render work unchanged except for a new badge.

**Tech Stack:** R (dplyr, tibble, stringr, digest), testthat, DuckDB. All work in worktree `worktree-mlb-market-edge-flagging`.

**Spec:** `docs/superpowers/specs/2026-06-05-mlb-market-edge-flagging-design.md`

---

## File Structure

- **Create** `Answer Keys/MLB Answer Key/market_edge.R` — the `find_market_edges()` module (pure; reads canonical odds, returns a long frame of flagged market edges).
- **Create** `Answer Keys/tests/test_market_edge.R` — unit tests for the module.
- **Create** `Answer Keys/tests/test_compute_bet_row_id.R` — unit tests for the shared hash helper.
- **Modify** `Answer Keys/MLB Answer Key/odds_screen.R` — add shared `compute_bet_row_id()` helper beside `.derive_market_type` / `.derive_period`.
- **Modify** `Answer Keys/MLB Answer Key/MLB.R` — use the shared helper for the model hash; move the `book_odds_by_book` block up; insert the market merge.
- **Modify** `Answer Keys/MLB Dashboard/mlb_dashboard.R` — edge badge, dual-EV in the hero strip, `.edge-badge` CSS.
- **Modify** `Answer Keys/CLAUDE.md` — document the dual-source flagging.

**Test command (run from the `Answer Keys/` directory):**
```bash
Rscript -e 'testthat::test_file("tests/test_<name>.R")'
```

---

## Task 1: Shared `compute_bet_row_id()` helper

**Why:** The model path and the market path must produce an *identical* `bet_row_id` for the same wager so they can be joined. The model currently hashes the verbose `market` string (`"totals_1st_5_innings"`), which the market path cannot reconstruct from canonical odds. This helper hashes the *canonical identity* both paths can compute: `(id, base_market_type, period, normalized_line, bet_on)`.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (add helper after `.derive_market_type`, ~line 55)
- Test: `Answer Keys/tests/test_compute_bet_row_id.R`

- [ ] **Step 1: Write the failing test**

Create `Answer Keys/tests/test_compute_bet_row_id.R`:

```r
# Answer Keys/tests/test_compute_bet_row_id.R
# compute_bet_row_id() must hash the CANONICAL bet identity so the model and
# market paths produce identical ids for the same wager.
# Run from "Answer Keys/":
#   Rscript -e 'testthat::test_file("tests/test_compute_bet_row_id.R")'

library(testthat)
source("../MLB Answer Key/odds_screen.R")

test_that("alternate_ and main collapse to the same id", {
  main <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  alt  <- compute_bet_row_id("g1", "alternate_totals", "F5", 8.5, "Over")
  expect_equal(main, alt)
})

test_that("line is normalized (8.5 == 8.50) and NA -> empty", {
  expect_equal(compute_bet_row_id("g1", "totals", "F5", 8.5,  "Over"),
               compute_bet_row_id("g1", "totals", "F5", 8.50, "Over"))
  # moneyline (NA line) is stable and distinct from a 0 line
  id_na <- compute_bet_row_id("g1", "h2h", "F5", NA, "NYY")
  expect_type(id_na, "character")
  expect_false(is.na(id_na))
})

test_that("different identity fields produce different ids", {
  base <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  expect_false(base == compute_bet_row_id("g2", "totals", "F5", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "spreads", "F5", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "FG", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "F5", 9.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "F5", 8.5, "Under"))
})

test_that("vectorized over inputs", {
  ids <- compute_bet_row_id(c("g1","g2"), c("totals","spreads"),
                            c("F5","F5"), c(8.5, -1.5), c("Over","NYY"))
  expect_length(ids, 2)
  expect_false(ids[1] == ids[2])
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compute_bet_row_id.R")'`
Expected: FAIL with `could not find function "compute_bet_row_id"`.

- [ ] **Step 3: Write minimal implementation**

In `Answer Keys/MLB Answer Key/odds_screen.R`, immediately after the `.derive_market_type` function (ends ~line 55) and before `.related_market_types`, insert:

```r
#' Canonical, book-agnostic bet identity hash shared by the model path
#' (MLB.R) and the market-edge path (market_edge.R). Hashes
#' (id, base_market_type, period, normalized_line, bet_on) so that the same
#' wager always maps to the same id regardless of the verbose market string
#' or the alternate_/main labeling. `market_type` may be the verbose model
#' string OR a canonical type; callers should pass the canonical type, but the
#' alternate_ prefix is stripped here defensively.
compute_bet_row_id <- function(id, market_type, period, line, bet_on) {
  base_mt   <- gsub("^alternate_", "", market_type)
  line_norm <- ifelse(is.na(line), "", as.character(round(as.numeric(line), 1)))
  key <- paste(id, base_mt, period, line_norm, bet_on, sep = "|")
  vapply(key, function(s) digest::digest(s, algo = "md5"), character(1),
         USE.NAMES = FALSE)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compute_bet_row_id.R")'`
Expected: PASS (all 4 tests).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_compute_bet_row_id.R"
git commit -m "feat(mlb): add shared compute_bet_row_id() canonical identity hash"
```

---

## Task 2: `find_market_edges()` module

**Why:** This is the new market signal. Pure function: in = canonical per-book odds list, out = long frame of flagged market edges in a schema that merges onto `mlb_bets_combined`.

**Files:**
- Create: `Answer Keys/MLB Answer Key/market_edge.R`
- Test: `Answer Keys/tests/test_market_edge.R`

- [ ] **Step 1: Write the failing test**

Create `Answer Keys/tests/test_market_edge.R`:

```r
# Answer Keys/tests/test_market_edge.R
# Unit tests for find_market_edges() — leave-one-out market-consensus edges.
# Run from "Answer Keys/":
#   Rscript -e 'testthat::test_file("tests/test_market_edge.R")'

library(testthat)
library(dplyr)
library(tibble)
source("../Tools.R")
source("../MLB Answer Key/odds_screen.R")
source("../MLB Answer Key/market_edge.R")

FT <- as.POSIXct("2026-06-05 18:00:00", tz = "UTC")
NOW <- as.POSIXct("2026-06-05 18:05:00", tz = "UTC")   # 5 min later

# canonical book row (matches scraper_to_canonical output shape)
crow <- function(game_id, market, period, side, line, american_odds, fetch_time = FT) {
  tibble(game_id = game_id, market = market, period = period, side = side,
         line = line, american_odds = as.integer(american_odds), fetch_time = fetch_time)
}

# A totals market at one book = Over + Under at the same line.
tot <- function(game_id, book_side_odds, line = 8.5, ft = FT) {
  bind_rows(
    crow(game_id, "totals", "F5", "Over",  line, book_side_odds[1], ft),
    crow(game_id, "totals", "F5", "Under", line, book_side_odds[2], ft)
  )
}

test_that("a soft book beating the others flags as a market edge", {
  # Pinnacle & BFA price Over ~ -110 (fair ~ 52.4%); Wagerzon offers Over +104.
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-108, -112)),
    wagerzon = tot("g1", c(+104, -124))
  )
  out <- find_market_edges(books, now = NOW)
  over <- out %>% filter(bet_on == "Over")
  expect_equal(nrow(over), 1)
  expect_equal(over$bookmaker_key, "wagerzon")
  expect_equal(over$edge_source, "market")
  expect_true(over$ev >= 0.02)
  expect_equal(over$market_type, "totals")
  expect_equal(over$period, "F5")
  expect_equal(over$line, 8.5)
  # leave-one-out: Wagerzon excluded from its own yardstick
  expect_equal(over$n_books, 3L)
})

test_that("prices in line with the market produce no flag", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-108, -112)),
    wagerzon = tot("g1", c(-109, -111))
  )
  out <- find_market_edges(books, now = NOW)
  expect_equal(nrow(out), 0)
})

test_that("two books (edge + 1 other) is enough to flag", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+104, -124))
  )
  out <- find_market_edges(books, now = NOW)
  expect_true(any(out$bet_on == "Over" & out$bookmaker_key == "wagerzon"))
})

test_that("a single book (no other to compare) is skipped", {
  books <- list(wagerzon = tot("g1", c(+150, -180)))
  out <- find_market_edges(books, now = NOW)
  expect_equal(nrow(out), 0)
})

test_that("a one-sided book is excluded from the consensus", {
  # wagerzon quotes only Over; consensus must come from the two full books,
  # and the run must not error.
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-108, -112)),
    wagerzon = crow("g1", "totals", "F5", "Over", 8.5, +104)
  )
  out <- find_market_edges(books, now = NOW)
  # wagerzon has no pair -> not a yardstick contributor and not flagged here
  expect_false(any(out$bookmaker_key == "wagerzon"))
})

test_that("stale quotes are dropped before devigging", {
  stale_ft <- as.POSIXct("2026-06-05 17:00:00", tz = "UTC")  # 65 min before NOW
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+104, -124), ft = stale_ft)
  )
  out <- find_market_edges(books, now = NOW, staleness_min = 30)
  expect_equal(nrow(out), 0)
})

test_that("bet_row_id matches the shared helper for the same wager", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+104, -124))
  )
  out <- find_market_edges(books, now = NOW) %>% filter(bet_on == "Over")
  expected <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  expect_equal(out$bet_row_id, expected)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_market_edge.R")'`
Expected: FAIL with `could not find function "find_market_edges"`.

- [ ] **Step 3: Write minimal implementation**

Create `Answer Keys/MLB Answer Key/market_edge.R`:

```r
# Answer Keys/MLB Answer Key/market_edge.R
# Market-consensus edge detection for the MLB bets tab.
#
# find_market_edges() reads the canonical per-book odds list MLB.R already
# builds (book_odds_by_book) and flags any book whose offered price beats the
# LEAVE-ONE-OUT median devigged fair of the OTHER books at the same wager.
# Pure: no DuckDB, no model dependency. Output rows merge onto mlb_bets_combined.
#
# Depends on (source first): Tools.R (devig_american, compute_ev, kelly_stake),
# odds_screen.R (compute_bet_row_id). Tested by tests/test_market_edge.R.

library(dplyr)
library(tibble)

# American odds -> raw implied probability (no devig).
.implied_prob <- function(o) ifelse(o > 0, 100 / (o + 100), -o / (-o + 100))

# Devig a single 2-outcome pair, returning the fair prob for each input side.
.pair_fairs <- function(o1, o2) {
  d <- devig_american(o1, o2)   # data.frame(p1, p2)
  c(d$p1, d$p2)
}

.empty_market_edges <- function() {
  tibble(
    bet_row_id = character(), id = character(), game_id = character(),
    market = character(), market_type = character(), period = character(),
    line = numeric(), bet_on = character(), bookmaker_key = character(),
    odds = integer(), prob = numeric(), ev = numeric(), n_books = integer(),
    bet_size = numeric(), edge_source = character(),
    model_ev = numeric(), market_ev = numeric(),
    home_team = character(), away_team = character(),
    pt_start_time = as.POSIXct(character(), tz = "UTC")
  )
}

#' @param book_odds_by_book named list of canonical frames
#'        (game_id, market, period, side, line, american_odds, fetch_time)
#' @param game_info optional tibble(game_id, home_team, away_team, pt_start_time)
#' @param threshold minimum market EV to flag (default 0.02)
#' @param min_others minimum OTHER books to form the leave-one-out yardstick
#' @param staleness_min drop quotes older than this many minutes before `now`
#' @param now reference time for staleness (default Sys.time())
#' @param bankroll,kelly_mult sizing inputs for kelly_stake()
find_market_edges <- function(book_odds_by_book,
                              game_info     = NULL,
                              threshold     = 0.02,
                              min_others    = 1,
                              staleness_min = 30,
                              now           = Sys.time(),
                              bankroll      = 100,
                              kelly_mult    = 0.25) {

  if (length(book_odds_by_book) == 0) return(.empty_market_edges())

  # 1. Stack books into one long frame tagged with bookmaker_key.
  long <- bind_rows(lapply(names(book_odds_by_book), function(bk) {
    df <- book_odds_by_book[[bk]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$bookmaker_key <- bk
    df
  }))
  if (is.null(long) || nrow(long) == 0) return(.empty_market_edges())

  # 2. Drop stale quotes.
  long <- long %>%
    filter(!is.na(fetch_time),
           as.numeric(difftime(now, fetch_time, units = "mins")) <= staleness_min)
  if (nrow(long) == 0) return(.empty_market_edges())

  # base market type collapses alt + main into one bucket (like
  # .related_market_types); line_key normalizes 8.5 vs 8.50.
  long <- long %>%
    mutate(base_mt  = gsub("^alternate_", "", market),
           line_key = ifelse(is.na(line), NA_real_, round(as.numeric(line), 1)),
           absline  = ifelse(is.na(line_key), -1, abs(line_key)))

  # 3. Per book, devig each 2-outcome pair. A pair = the two rows sharing
  #    (game, base_mt, period, book, |line|): Over+Under for totals, the two
  #    teams (+L/-L) for spreads, the two teams (line NA) for h2h. Books that
  #    quote only one side form no pair and are dropped.
  paired <- long %>%
    group_by(game_id, base_mt, period, bookmaker_key, absline) %>%
    filter(n() == 2) %>%
    mutate(fair = .pair_fairs(american_odds[1], american_odds[2])) %>%
    ungroup() %>%
    filter(!is.na(fair))
  if (nrow(paired) == 0) return(.empty_market_edges())

  # 4. Leave-one-out consensus: each book judged vs the median fair of the
  #    OTHER books at the same (game, base_mt, period, side, line).
  judged <- paired %>%
    group_by(game_id, base_mt, period, side, line_key) %>%
    mutate(
      n_books  = n(),
      loo_fair = vapply(seq_len(n()),
                        function(i) stats::median(fair[-i]), numeric(1))
    ) %>%
    ungroup() %>%
    filter(n_books - 1 >= min_others, !is.na(loo_fair))
  if (nrow(judged) == 0) return(.empty_market_edges())

  # 5. EV of each book's price vs its leave-one-out consensus; keep flags.
  flagged <- judged %>%
    mutate(market_ev = compute_ev(loo_fair, .implied_prob(american_odds))) %>%
    filter(market_ev >= threshold)
  if (nrow(flagged) == 0) return(.empty_market_edges())

  # Best book per (side, line), then best line per side (mirror model dedup).
  best <- flagged %>%
    group_by(game_id, base_mt, period, side, line_key) %>%
    slice_max(market_ev, n = 1, with_ties = FALSE) %>%
    group_by(game_id, base_mt, period, side) %>%
    slice_max(market_ev, n = 1, with_ties = FALSE) %>%
    ungroup()

  out <- best %>%
    transmute(
      id            = game_id,
      game_id       = game_id,
      market        = base_mt,        # canonical type; .derive_market_type() is idempotent on it
      market_type   = base_mt,
      period        = period,
      line          = line_key,
      bet_on        = side,
      bookmaker_key = bookmaker_key,
      odds          = as.integer(american_odds),
      prob          = loo_fair,
      ev            = market_ev,
      n_books       = as.integer(n_books),
      edge_source   = "market",
      model_ev      = NA_real_,
      market_ev     = market_ev,
      bet_size      = kelly_stake(market_ev, .implied_prob(american_odds),
                                  bankroll, kelly_mult)
    ) %>%
    mutate(bet_row_id = compute_bet_row_id(id, market_type, period, line, bet_on))

  # Attach game header info if provided.
  if (!is.null(game_info) && nrow(game_info) > 0) {
    out <- out %>% left_join(game_info, by = "game_id")
  } else {
    out <- out %>% mutate(home_team = NA_character_, away_team = NA_character_,
                          pt_start_time = as.POSIXct(NA, tz = "UTC"))
  }
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_market_edge.R")'`
Expected: PASS (all 7 tests). If the "soft book" EV assertion is off, print `over$ev` to confirm the devig math — Wagerzon Over +104 vs a ~52% consensus should be ≈ +6%.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/market_edge.R" "Answer Keys/tests/test_market_edge.R"
git commit -m "feat(mlb): add find_market_edges() leave-one-out market-consensus module"
```

---

## Task 3: Wire the market path into `MLB.R`

**Why:** Use the shared hash for the model bets, build `book_odds_by_book` before correlation-Kelly, compute market edges, and merge so `mlb_bets_combined` carries `edge_source` / `model_ev` / `market_ev`.

This task is integration (no unit test); verification is running the pipeline against copied live DBs. Make the edits, then run the verification step.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R`

- [ ] **Step 1: Source the new module**

Find where `MLB.R` sources `odds_screen.R` (search `source(` near the top). Add, on the next line:

```r
source("market_edge.R")
```

(Use the same relative-path style already used for `odds_screen.R`; if it sources via an absolute/`here`-style path, mirror that.)

- [ ] **Step 2: Replace the inline model hash with the shared helper**

In `MLB.R` lines 849–853, replace:

```r
  mutate(bet_row_id = vapply(
    paste(id, market, ifelse(is.na(line), "", as.character(line)), bet_on, sep = "|"),
    function(s) digest::digest(s, algo = "md5"),
    character(1)
  )) %>%
```

with:

```r
  mutate(bet_row_id = compute_bet_row_id(
    id, .derive_market_type(market), .derive_period(market), line, bet_on
  )) %>%
  mutate(edge_source = "model", model_ev = ev, market_ev = NA_real_) %>%
```

- [ ] **Step 3: Move the `book_odds_by_book` block above the correlation-Kelly call**

Cut lines ~905–936 (the block starting `# Game id lookup:` through `book_odds_by_book <- Filter(Negate(is.null), book_odds_by_book)`) and paste it immediately AFTER the model `all_bets_combined` pipeline ends (after line 860, before the `# Load placed bets` block at 862). Delete the original block from its old location (so it is not built twice).

- [ ] **Step 4: Insert the market merge right after the moved block**

Immediately after the pasted `book_odds_by_book <- Filter(...)` line, add:

```r
# ---- Market-consensus edges (model OR market flagging) ----
# Game header info so market-only cards have teams + tipoff and survive the
# future-game filter.
.market_game_info <- mlb_odds %>%
  transmute(game_id = id, home_team, away_team, pt_start_time) %>%
  distinct(game_id, .keep_all = TRUE)

market_bets <- find_market_edges(
  book_odds_by_book,
  game_info  = .market_game_info,
  threshold  = EV_THRESHOLD,
  min_others = 1,
  staleness_min = BOOK_STALENESS_CUTOFF_MIN,
  now        = Sys.time(),
  bankroll   = bankroll,
  kelly_mult = kelly_mult
) %>%
  filter(is.na(pt_start_time) | pt_start_time > Sys.time())

cat(sprintf("find_market_edges: %d market-edge candidate(s).\n", nrow(market_bets)))

# Merge model + market on the shared canonical bet_row_id.
.model_cols <- names(all_bets_combined)
market_new  <- market_bets %>% filter(!(bet_row_id %in% all_bets_combined$bet_row_id))
market_both <- market_bets %>% filter(bet_row_id %in% all_bets_combined$bet_row_id) %>%
  select(bet_row_id, market_ev_join = market_ev)

# (a) Promote model rows that ALSO have a market edge to "both".
all_bets_combined <- all_bets_combined %>%
  left_join(market_both, by = "bet_row_id") %>%
  mutate(
    edge_source = ifelse(!is.na(market_ev_join), "both", edge_source),
    market_ev   = ifelse(!is.na(market_ev_join), market_ev_join, market_ev)
  ) %>%
  select(-market_ev_join)

# (b) Append market-only rows, aligning columns to the model schema.
if (nrow(market_new) > 0) {
  miss <- setdiff(.model_cols, names(market_new))
  for (cc in miss) market_new[[cc]] <- NA
  market_new <- market_new %>% select(all_of(.model_cols))
  all_bets_combined <- bind_rows(all_bets_combined, market_new)
}

# Rank by the better of the two signals; drop any exact-id duplicate.
all_bets_combined <- all_bets_combined %>%
  mutate(ev = pmax(model_ev, market_ev, na.rm = TRUE)) %>%
  distinct(bet_row_id, .keep_all = TRUE) %>%
  arrange(desc(ev))
```

**Note:** the `mutate(ev = pmax(...))` only changes the *ranking* `ev`; `model_ev` and `market_ev` are preserved for display. The existing `adjust_kelly_for_correlation()` call (line ~881) now runs on this merged frame unchanged.

- [ ] **Step 5: Verify the new columns reach the written table**

The write at line 971 is `dbWriteTable(con_bets, "mlb_bets_combined", all_bets_combined)` — it writes all columns, so `edge_source` / `model_ev` / `market_ev` flow automatically. No change needed. Confirm by reading the source around 970–971 that it is still `dbWriteTable(... all_bets_combined)` with no explicit column list.

- [ ] **Step 6: Run the pipeline against copied live DBs and inspect**

Per CLAUDE.md, do NOT symlink DuckDBs — copy the live DBs the pipeline reads/writes into
the worktree, or (simplest and safest) run this verification from `main` after merge. The
pipeline that produces `mlb_bets_combined` is invoked via `python3 run.py mlb` (sharp
scrapers → rec scrapers → `MLB.R`); see `Answer Keys/MLB Dashboard/run.sh`. From the repo
root:

```bash
python3 run.py mlb        # runs the scrapers + MLB Answer Key/MLB.R end-to-end
```

Then inspect:

```bash
Rscript -e '
  library(DBI); library(duckdb)
  con <- dbConnect(duckdb(), "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
  df <- dbGetQuery(con, "SELECT edge_source, COUNT(*) n, ROUND(AVG(ev),4) avg_ev
                         FROM mlb_bets_combined GROUP BY edge_source")
  print(df)
  cat("\nSample market/both rows:\n")
  print(dbGetQuery(con, "SELECT bet_on, market, line, bookmaker_key, edge_source,
                           ROUND(model_ev,4) model_ev, ROUND(market_ev,4) market_ev
                         FROM mlb_bets_combined
                         WHERE edge_source IN (\x27market\x27,\x27both\x27) LIMIT 10"))
  bp <- dbGetQuery(con, "SELECT COUNT(*) n FROM mlb_bets_book_prices")
  cat(sprintf("\nbook_prices rows: %d\n", bp$n))
  dbDisconnect(con, shutdown = TRUE)
'
```

Expected: `edge_source` has `model` (and likely `market` / `both`) groups; market/both rows show a populated `market_ev`; `mlb_bets_book_prices` is non-empty. Confirm no error and that pre-existing model bets are still present (regression: model rows did not vanish).

- [ ] **Step 7: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "feat(mlb): merge market-consensus edges into mlb_bets_combined"
```

---

## Task 4: Dashboard badge + dual-EV display

**Why:** Surface `edge_source` as a MODEL / MARKET / BOTH badge and show both EVs on BOTH cards.

This task is verified by rendering (per CLAUDE.md verify-by-rendering rule), not unit tests.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R`

- [ ] **Step 1: Compute display EV fields in the table mutate**

In `create_bets_table()`, in the `mutate()` block at lines ~1595–1619, after the `ev_pct = ev * 100,` line add:

```r
    model_ev_pct   = ifelse(is.na(model_ev),  NA_real_, model_ev  * 100),
    market_ev_pct  = ifelse(is.na(market_ev), NA_real_, market_ev * 100),
    edge_source    = ifelse(is.na(edge_source), "model", edge_source),
```

(If a stale `mlb_bets_combined` predates the new columns, `SELECT *` won't supply them. Guard at the top of `create_bets_table()`, right after it receives `all_bets`: add `if (!"edge_source" %in% names(all_bets)) all_bets$edge_source <- "model"`, and likewise default `model_ev <- all_bets$ev` and `market_ev <- NA_real_` when absent.)

- [ ] **Step 2: Build the badge HTML per row and pass it down**

Inside the row loop (starts ~line 1658), where the per-row card is assembled, compute a badge string. Add near the top of the loop body (after `row <- table_data[i, ]`):

```r
  edge_src   <- if (is.na(row$edge_source)) "model" else row$edge_source
  badge_label <- c(model = "Model", market = "Market", both = "★ Both")[[edge_src]]
  edge_badge_html <- sprintf('<span class="edge-badge %s">%s</span>',
                             edge_src, badge_label)
```

Pass `edge_badge_html` into the `render_bet_card()` call as the correlation-badge argument is passed, OR append it to the existing `corr_badge_html` so both badges render in the `.bet-line`. Concretely, where `corr_badge_html` is built for this row, change the final value to include the edge badge:

```r
  corr_badge_html <- paste0(edge_badge_html, corr_badge_html)
```

(If `corr_badge_html` is constructed inside `render_bet_card`, instead add a new `edge_badge_html` parameter to `render_bet_card()` and inject it in the `.bet-line` sprintf at line ~1322–1332 next to `corr_badge_html`.)

- [ ] **Step 3: Dual-EV in the hero strip**

Change `render_hero_strip()` (lines 1249–1294) to accept and render both EVs. Update the signature and the EV stat block.

Replace the signature line:

```r
render_hero_strip <- function(pick_book, pick_odds, fair_odds,
                               ev_pct, risk_dollars, towin_dollars,
                               action_html) {
```

with:

```r
render_hero_strip <- function(pick_book, pick_odds, fair_odds,
                               ev_pct, risk_dollars, towin_dollars,
                               action_html,
                               model_ev_pct = NA_real_, market_ev_pct = NA_real_) {
```

Replace the single EV `<div class="stat">` (the line containing `<span class="lbl">EV</span><span class="val ev">%s</span>`) with a computed `ev_block` inserted via `%s`. Just before the `sprintf('<div class="hero">...`, add:

```r
  ev_block <- if (!is.na(model_ev_pct) && !is.na(market_ev_pct)) {
    sprintf(
      '<div class="stat"><span class="lbl">Model EV</span><span class="val ev">+%.1f%%</span></div>
       <div class="stat"><span class="lbl">Market EV</span><span class="val evm">+%.1f%%</span></div>',
      model_ev_pct, market_ev_pct)
  } else if (!is.na(market_ev_pct)) {
    sprintf('<div class="stat"><span class="lbl">Market EV</span><span class="val evm">+%.1f%%</span></div>',
            market_ev_pct)
  } else {
    sprintf('<div class="stat"><span class="lbl">EV</span><span class="val ev">%s</span></div>',
            ev_str)
  }
```

In the `sprintf` template, replace the substring
`<div class="stat"><span class="lbl">EV</span><span class="val ev">%s</span></div>`
with just `%s`, and in the argument list change the 4th argument from `ev_str` to
`ev_block` (same position, still one `%s`). Keep all other stat blocks and arguments
unchanged. `ev_str` is still computed at the top of the function and is used inside
`ev_block`'s fallback branch, so leave its definition in place.

Then, at the `render_hero_strip(...)` call site (search for it within the row loop), pass the two new args: `model_ev_pct = row$model_ev_pct, market_ev_pct = row$market_ev_pct`.

- [ ] **Step 4: Add `.edge-badge` and `.val.evm` CSS**

In the `tags$style(HTML('...'))` block, after the `.bet-card-v8 .hero .stat .val.fair { ... }` rule (~line 3072), append:

```css
.bet-card-v8 .hero .stat .val.evm { color: #58a6ff; }
.bet-card-v8 .edge-badge {
  font-size: 11px; font-weight: 700; letter-spacing: 0.06em;
  padding: 3px 10px; border-radius: 999px; text-transform: uppercase;
  margin-left: 8px;
}
.bet-card-v8 .edge-badge.model  { color: #3fb950; background: rgba(63,185,80,0.15);  border: 1px solid rgba(63,185,80,0.45); }
.bet-card-v8 .edge-badge.market { color: #58a6ff; background: rgba(56,139,253,0.15); border: 1px solid rgba(56,139,253,0.45); }
.bet-card-v8 .edge-badge.both   { color: #e3b341; background: rgba(227,179,65,0.15); border: 1px solid rgba(227,179,65,0.5); }
```

- [ ] **Step 5: Render and verify visually**

Render the dashboard against copied live DBs (`mlb_mm.duckdb` must contain the
new-schema `mlb_bets_combined` from Task 3). `mlb_dashboard.R` writes `report.html`
(`OUTPUT_PATH`, line 43):

```bash
cd "Answer Keys/MLB Dashboard" && Rscript mlb_dashboard.R
# then open report.html
```

Open the output and confirm:
- A model-only bet shows a green **Model** badge and a single EV.
- A market-only bet shows a blue **Market** badge and a single Market EV.
- A both bet shows a gold **★ Both** badge and two EVs (Model + Market).
- No card is broken / no `NA%` leaks into the hero.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): MODEL/MARKET/BOTH edge badge + dual-EV hero"
```

---

## Task 5: Documentation

**Files:**
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Document the dual-source flagging**

In `Answer Keys/CLAUDE.md`, in the bets-tab / "Odds screen" section, add a subsection describing:
- The two signals (model EV ≥ 2% OR market EV ≥ 2%) and the `edge_source` / `model_ev` / `market_ev` columns now on `mlb_bets_combined`.
- `find_market_edges()` in `market_edge.R`: leave-one-out consensus, `min_others = 1`, like-for-like matching, best-line-per-side collapse.
- The shared `compute_bet_row_id()` (canonical identity) in `odds_screen.R` used by both the model and market paths.

- [ ] **Step 2: Commit**

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb): document dual-source edge flagging on the bets tab"
```

---

## Final verification before merge

- [ ] All unit tests pass:
  ```bash
  cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compute_bet_row_id.R")' \
    && Rscript -e 'testthat::test_file("tests/test_market_edge.R")' \
    && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")' \
    && Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")'
  ```
- [ ] Pipeline runs clean on copied live DBs; `mlb_bets_combined` has `edge_source`/`model_ev`/`market_ev`; market & both rows present; `mlb_bets_book_prices` populated; model bets not lost (regression).
- [ ] Dashboard renders with correct badges + dual-EV (verify-by-rendering).
- [ ] Timezone gate (scraper-adjacent safety): `python tests/timezone_parity_test.py` still passes (no scraper changes expected, but the market path reads `fetch_time` — confirm no regression).
- [ ] Executive engineer review of `git diff main..HEAD` (use review-my-code skill): data integrity (no dup writes, future-game filter applied to market rows), resource safety (no new DB leaks), edge cases (off-season → 0 market rows → tab behaves as today), dead code, security.
- [ ] Get explicit user approval, merge to `main`, then remove worktree + branch.
```

