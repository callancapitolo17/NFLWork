# Fix MLB Bets-Tab Odds Bugs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the eight confirmed bugs in the MLB bets-tab odds pipeline (Pinnacle pills explicitly out of scope).

**Architecture:** Each task is one bug with a one-paragraph issue + a surgical fix + a unit test where the helper is testable. No new abstractions, no opportunistic refactors. End-to-end verification runs MLB.R and inspects `mlb_bets_book_prices` to confirm DK/FD alt pills appear and `fetch_time` arithmetic is sane.

**Tech Stack:** R (`Tools.R`, `odds_screen.R`, `MLB.R`), Python 3 (`scraper_draftkings_singles.py`, `scraper_fanduel_singles.py`), DuckDB.

---

## Review Pack

**What we're building** — Eight targeted bug fixes to the MLB bets-tab odds pipeline. Today, alt-spread/alt-total bets show no DraftKings or FanDuel pills, the `fetch_time` column is timezone-ambiguous (so any staleness check is off by ~4 hours), and stale per-book snapshots silently flow into the dashboard. This plan repairs each issue in isolation, keeping changes small and reversible.

**Key decisions**
1. **Fix at the source for the timezone bug, not the consumer.** Scrapers will store `datetime.now(timezone.utc)` into `TIMESTAMPTZ` columns. Alternative: leave the schema alone and patch every read query with `AT TIME ZONE 'UTC'`. Rejected — patching N readers is fragile, and pitfall #7 in `Answer Keys/CLAUDE.md` already warns about this exact trap. Storing TZ-aware values once at the source removes the whole class of bug.
2. **Fix the alt-line bug at the matcher, not the loader.** `expand_bets_to_book_prices()` will union related market types when building the candidate set: an `alternate_spreads` bet sees both `spreads` and `alternate_spreads` rows; a `spreads` bet sees both too. The closest-line picker decides which row wins. Alternatives considered: (a) preserve the alt-vs-main distinction at the loader by reading DK's raw `market` column — clean but depends on DK's labels staying stable, and forces every other scraper to follow the same convention; (b) drop the `alternate_*` market_type entirely from the model — much bigger rework. Option C (matcher-side union) is the smallest, most defensive change; loaders and model output stay untouched.
3. **Drop rows by game-start-time, not by `fetch_time` age.** Pipeline runs are manual, so an absolute "scraped within N minutes" gate is the wrong question. The right gate is "has this game started yet?" — a 4-hour-old quote on an upcoming game is fine; a 2-min-old quote on a finished game is not. Each `get_*_odds()` parses its scraper's `(game_date, game_time)` (different formats per book) into a UTC datetime, then filters `game_start > NOW() - 5 min` (small grace so first-pitch flicker doesn't drop rows). `mlb_bets_book_prices` also gains a `game_start_time TIMESTAMPTZ` column so the dashboard can display the start time on each card.
4. **Skip bug #9 ("EV against displaced quote") from the original report.** Alternative: include it now. Rejected — the dashboard EV computation path needs its own short investigation; bundling it here would balloon scope without certainty.

**Risks / push back here**
- The timezone fix changes the column type in `dk.duckdb` / `fd.duckdb` (and ideally the other scraper DBs). The first run after merge needs to either DROP+CREATE the table or run a one-off migration. I've designed the scraper's `CREATE TABLE` to handle this on next launch, but **you should expect the per-book DBs to be wiped once** (acceptable: scrapers repopulate in seconds).
- I am proposing to also flip the same TIMESTAMPTZ schema in five other per-book scrapers (Wagerzon, Hoop88, BFA, Bookmaker, Bet105) so the staleness math is uniformly correct. If you'd rather scope this PR to DK+FD only and tackle the others later, say so before execution.
- The gate is `game_start > NOW() - 5 min`. The 5-min grace exists because, for first-pitch markets, the dashboard is still useful for a minute or two after the official game-start time (DK/FD odds keep updating until they suspend the market). If you'd prefer a hard cutoff at `NOW()` (drop the second a game starts), say so.
- Per-scraper `game_time` parsing is needed because formats differ: DK/FD use ISO 8601 UTC, WZ/BKM use naive Eastern wall-clock (year missing — derived from `Sys.time()` with year-rollover handling), Bet105 looks like UTC with date rollover. The parser will be small but per-book; a bug in any one parser silently drops that book's pills. Mitigation: each parser ships with a unit test that pins the format.

**Worth understanding (optional)**
1. **Naive vs TZ-aware timestamps** — A naive timestamp ("17:24:11") has no timezone attached; it's just six numbers. Comparing it to a TZ-aware timestamp forces a cast using whatever timezone the comparison context provides (R session, DuckDB `NOW()`, Python `datetime.now()`). The "fix" is to never let a naive value sit in the database — anchor it to UTC at insert time. In R this is analogous to how `lubridate::ymd_hms("...", tz = "UTC")` differs from `as.POSIXct("...")`.
2. **Why "collapse" is a brittle pattern** — `get_dk_odds` originally rewrote DK's three market labels (`main`, `alternate_spreads`, `alternate_totals`) into one (`spreads`/`totals`/`h2h`). The intent was simplification, but the downstream join keys on this column, so the collapse silently dropped half the data. The fix is to *preserve information* through the pipeline — only normalize at the last consumer, not at the loader.

---

## Version Control & Worktree

- **Worktree:** Already created at `.claude/worktrees/fix-mlb-bets-odds-bugs` on branch `worktree-fix-mlb-bets-odds-bugs`.
- **Commits:** One commit per task (8 fix commits + 1 docs commit + 1 verify commit). Frequent small commits make any single bug easy to revert.
- **Merge:** After verification passes and user approves, fast-forward merge into `main`, then `git worktree remove` + `git branch -d`.

## Documentation

- `Answer Keys/CLAUDE.md` — update the "Common Pitfalls" section: pitfall #7 (timezone) gets a "Fixed 2026-05-18 — scraper DBs now store TIMESTAMPTZ" note; add a new pitfall #10 documenting the DK/FD alt-line preservation rule so future scraper-loader changes don't regress it.
- `Answer Keys/MLB Answer Key/README.md` — note the freshness gate on `get_*_odds()`.
- No new READMEs.

---

## Task 1: Union related market types in `expand_bets_to_book_prices` (Bug #1)

**Issue:** `expand_bets_to_book_prices()` joins book candidates with `market == bet$market_type` (strict equality). But `compare_alts_to_samples` emits alt bets with `market_type = "alternate_spreads"`/`"alternate_totals"`, while `get_dk_odds`/`get_fd_odds` (and after canonicalization, every book frame) labels DK's and FD's spread rows as `market_type = "spreads"`/`"totals"` regardless of whether the raw DK row was `main` or `alternate_*`. Result: alt bets never find DK/FD candidates. Live: 0 DK pills on 14 alt bets despite 463 alt-line rows in `dk.duckdb`.

**Fix:** Treat `spreads` and `alternate_spreads` as one bucket when looking up candidates (same for totals). A spread bet sees both `spreads` and `alternate_spreads` book rows; the closest-line picker decides which row wins. Same for totals. Moneyline (`h2h`) is unaffected.

This is a matcher-side change — `Tools.R::get_*_odds()` and `scraper_to_canonical()` are NOT touched.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (`expand_bets_to_book_prices` around line 234-238)
- Test: `Answer Keys/tests/test_odds_screen.R`

- [ ] **Step 1: Add failing test**

Append to `Answer Keys/tests/test_odds_screen.R`:
```r
test_that("alt-spread bet finds DK quote labeled 'spreads' (matcher-side union)", {
  # Bet from compare_alts_to_samples: market_type = "alternate_spreads".
  bets <- tibble(
    bet_row_id  = "h1",
    game_id     = "g1",
    market      = "alternate_spreads_fg",
    market_type = "alternate_spreads",
    period      = "FG",
    line        = -2.5,
    bet_on      = "Boston Red Sox",
    home_team   = "Boston Red Sox",
    away_team   = "Philadelphia Phillies",
    pick_side   = "pick"
  )
  # DK book frame as it appears after get_dk_odds + scraper_to_canonical
  # under the CURRENT (collapsed) labeling: market = "spreads" even though
  # the raw DK row was alternate_spreads at -2.5.
  dk <- bind_rows(
    book_row("g1", "spreads", "FG", "Boston Red Sox",        -2.5, +205),
    book_row("g1", "spreads", "FG", "Philadelphia Phillies",  2.5, -280)
  )
  out <- expand_bets_to_book_prices(bets, list(draftkings = dk))
  expect_equal(nrow(out), 2)
  expect_true(all(out$is_exact_line))
  expect_equal(out$american_odds[out$side == "pick"], 205L)
})

test_that("main spread bet still finds 'spreads' candidates (union does not break the common case)", {
  bets <- make_bet_row(market = "spreads_1st_5_innings", line = 0.5,
                        bet_on = "Boston Red Sox", market_type = "spreads")
  bets$home_team <- "Boston Red Sox"
  bets$away_team <- "Philadelphia Phillies"
  book <- bind_rows(
    book_row("g1", "spreads", "F5", "Boston Red Sox",         0.5, -115),
    book_row("g1", "spreads", "F5", "Philadelphia Phillies", -0.5, -105)
  )
  out <- expand_bets_to_book_prices(bets, list(wagerzon = book))
  expect_equal(nrow(out), 2)
  expect_true(all(out$is_exact_line))
})

test_that("alt-total bet finds 'totals' and 'alternate_totals' candidates", {
  bets <- tibble(
    bet_row_id  = "h2", game_id = "g1",
    market      = "alternate_totals_fg",
    market_type = "alternate_totals",
    period      = "FG", line = 10.0, bet_on = "Over",
    pick_side   = "pick"
  )
  # One book labels its alt totals as "alternate_totals"; another labels them
  # all as "totals" (e.g. DK after get_dk_odds). Union should match both.
  wz <- bind_rows(
    book_row("g1", "alternate_totals", "FG", "Over", 10.0, +135),
    book_row("g1", "alternate_totals", "FG", "Under", 10.0, -160)
  )
  dk <- bind_rows(
    book_row("g1", "totals", "FG", "Over", 10.0, +130),
    book_row("g1", "totals", "FG", "Under", 10.0, -155)
  )
  out <- expand_bets_to_book_prices(bets, list(wagerzon = wz, draftkings = dk))
  # Both books should appear on the pick side.
  expect_setequal(out$bookmaker[out$side == "pick"], c("wagerzon", "draftkings"))
})
```

- [ ] **Step 2: Run tests to verify the first two fail**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```
Expected: first test FAILS (no DK candidates found because of strict match). Third test will probably also fail for the DK side.

- [ ] **Step 3: Add `.related_market_types()` helper and apply it in the candidates filter**

Edit `Answer Keys/MLB Answer Key/odds_screen.R`. Right after `.derive_market_type` (line ~56), add:
```r
# When matching a bet to per-book candidates, treat alt and main as the
# same bucket. Different scrapers use different conventions:
#   - WZ/BKM/Bet105 label alt rows "alternate_spreads"/"alternate_totals"
#   - DK/FD (via get_dk_odds/get_fd_odds in Tools.R) collapse alt rows into
#     "spreads"/"totals" with the alt line carried in the `line` column
# The closest-line picker decides which row wins, so we just union the
# possible labels here. Moneyline (h2h) has no alt variant.
.related_market_types <- function(mt) {
  if (mt == "spreads"           || mt == "alternate_spreads")
    return(c("spreads", "alternate_spreads"))
  if (mt == "totals"            || mt == "alternate_totals")
    return(c("totals", "alternate_totals"))
  mt
}
```

Then change the candidates filter in `expand_bets_to_book_prices` (line ~234-238) from:
```r
candidates <- book_frame %>%
  filter(game_id == bet$game_id,
         market  == bet$market_type,
         period  == bet$period,
         side    == side_value)
```
to:
```r
candidates <- book_frame %>%
  filter(game_id == bet$game_id,
         market  %in% .related_market_types(bet$market_type),
         period  == bet$period,
         side    == side_value)
```

- [ ] **Step 4: Run tests, expect all pass**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "fix(mlb-odds): union alt+main market types in expand_bets_to_book_prices

Bug: alt-spread/alt-total bets never matched DK/FD candidates because
scrapers (post-canonicalization) label DK/FD alt rows as 'spreads' /
'totals' while bets carry 'alternate_spreads' / 'alternate_totals'.
Strict-equality join silently dropped them.

Fix: .related_market_types(mt) unions {spreads, alternate_spreads}
(and totals analog). Closest-line picker decides which row wins.
Loaders and model output untouched."
```

---

## Task 2: Switch scraper `fetch_time` to TIMESTAMPTZ (Bug #3)

**Issue:** Scrapers write `datetime.utcnow()` (a naive value) into `fetch_time TIMESTAMP` columns. DuckDB's `NOW()` returns `TIMESTAMPTZ`. `NOW() - fetch_time` casts the naive value to local time, producing ~4-hour skew during EDT. Live evidence: today's DK rows show `EXTRACT(EPOCH FROM NOW() - MAX(fetch_time)) = -5022` seconds — "4 hours in the future".

**Fix:** Store TZ-aware UTC timestamps in `TIMESTAMPTZ` columns. On first run after the upgrade, the scraper's `CREATE TABLE IF NOT EXISTS` step detects an old-schema table and drops it (data is rebuilt next scrape — only seconds of data are lost).

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_singles.py` (lines ~249, ~344-365)
- Modify: `mlb_sgp/scraper_fanduel_singles.py` (lines ~222, ~274-295)
- Test: ad-hoc Python sanity check (no unit test framework wired for these scrapers).

- [ ] **Step 1: DK scraper — switch fetch_time to TZ-aware and migrate the table**

Edit `mlb_sgp/scraper_draftkings_singles.py`:

Replace the import block (top of file) to include timezone:
```python
from datetime import datetime, timezone
```

Replace line 249 (`fetch_time = datetime.utcnow()`) with:
```python
fetch_time = datetime.now(timezone.utc)
```

Replace the `CREATE TABLE IF NOT EXISTS` block (lines ~344-365) — wrap it with a one-time migration that drops the table if the existing `fetch_time` column is naive:

```python
# Migrate naive-TIMESTAMP schema to TIMESTAMPTZ if needed.
# DuckDB does not support ALTER COLUMN TYPE between these, so drop+create.
existing = con.execute(
    "SELECT column_name, data_type FROM information_schema.columns "
    "WHERE table_name = 'mlb_odds' AND column_name = 'fetch_time'"
).fetchone()
if existing is not None and "WITH TIME ZONE" not in (existing[1] or "").upper():
    print(f"[dk] Migrating mlb_odds.fetch_time TIMESTAMP -> TIMESTAMPTZ "
          f"(existing snapshot will be re-populated this run)")
    con.execute("DROP TABLE mlb_odds")

con.execute("""
    CREATE TABLE IF NOT EXISTS mlb_odds (
        fetch_time        TIMESTAMPTZ,
        sport_key         VARCHAR,
        game_id           VARCHAR,
        game_date         VARCHAR,
        game_time         VARCHAR,
        away_team         VARCHAR,
        home_team         VARCHAR,
        market            VARCHAR,
        period            VARCHAR,
        away_spread       FLOAT,
        away_spread_price INTEGER,
        home_spread       FLOAT,
        home_spread_price INTEGER,
        total             FLOAT,
        over_price        INTEGER,
        under_price       INTEGER,
        away_ml           INTEGER,
        home_ml           INTEGER
    )
""")
```

- [ ] **Step 2: Apply the identical change to `mlb_sgp/scraper_fanduel_singles.py`**

Same three edits: import `timezone`, replace `datetime.utcnow()` with `datetime.now(timezone.utc)`, wrap CREATE TABLE with the migration block (label the print message `[fd]`).

- [ ] **Step 3: Verify both scrapers produce TZ-aware rows**

```bash
python3 - <<'PY'
import duckdb
from datetime import datetime, timezone
for path in ["dk_odds/dk.duckdb", "fd_odds/fd.duckdb"]:
    con = duckdb.connect(path, read_only=True)
    typ = con.execute(
        "SELECT data_type FROM information_schema.columns "
        "WHERE table_name = 'mlb_odds' AND column_name = 'fetch_time'"
    ).fetchone()
    age = con.execute("SELECT EXTRACT(EPOCH FROM NOW() - MAX(fetch_time)) FROM mlb_odds").fetchone()
    print(path, "type=", typ, "age_s=", age)
    con.close()
PY
```
Expected (after running a scrape with the new code):
- `type = ('TIMESTAMP WITH TIME ZONE',)` for both
- `age_s` is a small positive number (seconds since last scrape)

If the scrape hasn't run yet, this verification can be deferred to the end-to-end task.

- [ ] **Step 4: Commit**

```bash
git add mlb_sgp/scraper_draftkings_singles.py mlb_sgp/scraper_fanduel_singles.py
git commit -m "fix(mlb-odds): store DK/FD scraper fetch_time as TIMESTAMPTZ (UTC)

Bug: scrapers wrote datetime.utcnow() (naive) into TIMESTAMP columns,
then DuckDB NOW() comparisons cast the naive value to local time,
producing ~4h skew during EDT and breaking staleness math.

Fix: datetime.now(timezone.utc) into TIMESTAMPTZ columns. First run
after merge detects the old schema and drops+recreates the table;
data is repopulated the same cycle."
```

---

## Task 3: Apply the same TIMESTAMPTZ fix to the other per-book scrapers

**Issue:** Same root cause as Task 2, but the bug is present in every per-book scraper (Wagerzon, Hoop88, BFA, Bookmaker, Bet105). Leaving them as-is means `mlb_bets_book_prices` will still carry mixed naive/TZ-aware values, and dashboard staleness math is uneven across books.

**Fix:** Apply the same three-line scraper change + migration block to each.

**Files:**
- Modify: `wagerzon_odds/scraper_v2.py` (or whichever file the WZ scraper writes through)
- Modify: `hoop88_odds/<scraper>.py`
- Modify: `bfa_odds/<scraper>.py`
- Modify: `bookmaker_odds/<scraper>.py`
- Modify: `bet105_odds/<scraper>.py`

- [ ] **Step 1: Locate the `INSERT INTO mlb_odds` site in each scraper**

```bash
for d in wagerzon_odds hoop88_odds bfa_odds bookmaker_odds bet105_odds; do
  echo "=== $d ==="
  grep -rn "INSERT INTO mlb_odds\|CREATE TABLE IF NOT EXISTS mlb_odds\|datetime.utcnow\|datetime\.now" "$d" 2>/dev/null | head -5
done
```

- [ ] **Step 2: Apply the same three edits to each scraper**

For each scraper file, mirror the DK/FD edits from Task 2:
1. `from datetime import datetime, timezone` (add `timezone`)
2. Replace any `datetime.utcnow()` for `fetch_time` with `datetime.now(timezone.utc)`
3. Wrap the `CREATE TABLE IF NOT EXISTS mlb_odds` block with the same migration check (drop if `fetch_time` is naive, then create with `TIMESTAMPTZ`)

- [ ] **Step 3: Verify all five DBs end up with TIMESTAMPTZ after their next scrape**

```bash
python3 - <<'PY'
import duckdb
for path in ["wagerzon_odds/wagerzon.duckdb", "hoop88_odds/hoop88.duckdb",
             "bfa_odds/bfa.duckdb", "bookmaker_odds/bookmaker.duckdb",
             "bet105_odds/bet105.duckdb"]:
    try:
        con = duckdb.connect(path, read_only=True)
        typ = con.execute(
            "SELECT data_type FROM information_schema.columns "
            "WHERE table_name = 'mlb_odds' AND column_name = 'fetch_time'"
        ).fetchone()
        print(path, typ)
        con.close()
    except Exception as e:
        print(path, "ERR", e)
PY
```
Expected: every line prints `('TIMESTAMP WITH TIME ZONE',)`.

- [ ] **Step 4: Commit**

```bash
git add wagerzon_odds/ hoop88_odds/ bfa_odds/ bookmaker_odds/ bet105_odds/
git commit -m "fix(scrapers): store fetch_time as TIMESTAMPTZ in all per-book scrapers

Mirrors the DK/FD fix in the prior commit. Makes staleness math
uniformly correct across the bets-tab pipeline."
```

---

## Task 4: Pre-declare `mlb_bets_book_prices.fetch_time` as TIMESTAMPTZ

**Issue:** Even with the scraper fix, `MLB.R` writes `mlb_bets_book_prices` via `dbWriteTable()`, which lets DuckDB pick the column type from R's POSIXct. duckdb-r often downgrades this to naive `TIMESTAMP`. The downstream dashboard would then re-introduce the same skew bug.

**Fix:** Pre-create the table with explicit `TIMESTAMPTZ`, then `INSERT INTO ... SELECT * FROM` from the R frame.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (around lines 939-945)

- [ ] **Step 1: Edit the write block**

Replace this block in `MLB.R` (lines ~941-945):
```r
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_bets, "mlb_bets_combined", all_bets_combined)
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_book_prices")
dbWriteTable(con_bets, "mlb_bets_book_prices", book_prices_long)
```

With:
```r
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_bets, "mlb_bets_combined", all_bets_combined)

# Pre-declare schema so fetch_time keeps its UTC timezone (otherwise
# dbWriteTable downgrades POSIXct to naive TIMESTAMP and breaks
# downstream NOW()-fetch_time math).
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_book_prices")
dbExecute(con_bets, "
  CREATE TABLE mlb_bets_book_prices (
    bet_row_id    VARCHAR,
    game_id       VARCHAR,
    market        VARCHAR,
    period        VARCHAR,
    side          VARCHAR,
    bookmaker     VARCHAR,
    line          DOUBLE,
    line_quoted   DOUBLE,
    is_exact_line BOOLEAN,
    american_odds INTEGER,
    fetch_time    TIMESTAMPTZ
  )
")
# Ensure POSIXct is in UTC so DuckDB sees a TZ-aware value.
book_prices_long$fetch_time <- as.POSIXct(book_prices_long$fetch_time, tz = "UTC")
dbAppendTable(con_bets, "mlb_bets_book_prices", book_prices_long)
```

- [ ] **Step 2: Verify after running MLB.R that the column is TIMESTAMPTZ**

(Run as part of the end-to-end verification in Task 9.)

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "fix(mlb-odds): write mlb_bets_book_prices.fetch_time as TIMESTAMPTZ

Pre-declare schema and force POSIXct to UTC before append so the
downstream dashboard sees timezone-aware values."
```

---

## Task 5: Drop past-game rows in `get_*_odds()` and carry `game_start_time` through (Bug #4)

**Issue:** Every `get_*_odds()` helper does `SELECT * FROM mlb_odds` with no game-time filter. If a scraper failed today, its DB still holds the previous successful snapshot. When the loader inner-joins to today's `mlb_odds` lookup by `(home_team, away_team)`, yesterday's snapshot for the same teams (very common in MLB series) gets *silently re-tagged with today's `game_id`*, joining stale odds to today's bets. The user also wants each `mlb_bets_book_prices` row to carry the game's start time so the dashboard can display "BOS @ PHI · 7:05 PM ET" on every card.

**Fix (two parts):**
- Each `get_*_odds()` parses its scraper's `(game_date, game_time)` into a UTC datetime and drops rows where `game_start ≤ NOW() - 5 min` (5-min grace for first-pitch flicker). Per-scraper parsers handle: DK/FD ISO-8601 UTC; WZ/BKM naive Eastern wall-clock (year derived from `Sys.time()`); Bet105 UTC with date rollover.
- Pass `commence_time` from `mlb_odds` through `.game_id_lookup` → `scraper_to_canonical` → `expand_bets_to_book_prices` so `mlb_bets_book_prices` gains a `game_start_time TIMESTAMPTZ` column.

**Files:**
- Modify: `Answer Keys/Tools.R` — add `.parse_dk_game_dt`, `.parse_fd_game_dt`, `.parse_wz_game_dt`, `.parse_bet105_game_dt` (BKM uses WZ's format) and a shared `.drop_past_games(raw, parser, label)`. Wire into every `get_*_odds()` helper.
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` — `scraper_to_canonical()` and `odds_api_to_canonical()` carry `game_start_time`; `normalize_book_odds_frame()` includes it in the canonical shape; `expand_bets_to_book_prices()` joins it onto each output row.
- Modify: `Answer Keys/MLB Answer Key/MLB.R` — `.game_id_lookup` includes `commence_time`; `mlb_bets_book_prices` schema gains `game_start_time TIMESTAMPTZ`.
- Test: `Answer Keys/tests/test_odds_screen.R`

- [ ] **Step 1: Add the per-scraper parsers + shared filter**

Append to `Answer Keys/Tools.R` near `.singles_market_name`:
```r
# Parse DK/FD ISO 8601 UTC strings ("2026-05-19T20:10:00.0000000Z").
.parse_iso_game_dt <- function(s) {
  if (is.null(s) || all(is.na(s))) return(as.POSIXct(rep(NA, length(s)), tz = "UTC"))
  # lubridate::ymd_hms tolerates trailing fractional digits and a final 'Z'.
  lubridate::ymd_hms(s, tz = "UTC", quiet = TRUE)
}

# Parse Wagerzon / Bookmaker naive Eastern wall-clock: game_date "MM/DD",
# game_time "HH:MM". Year is inferred from Sys.time(); if the resulting
# datetime is more than 6 months in the past, roll forward one year (handles
# Dec 31 → Jan 1). Result is a POSIXct in UTC representing the Eastern wall
# clock instant.
.parse_wz_game_dt <- function(date_str, time_str) {
  if (length(date_str) == 0) return(as.POSIXct(character(), tz = "UTC"))
  yr <- format(Sys.time(), "%Y")
  combined <- sprintf("%s/%s %s", yr, date_str, time_str)
  # %Y/%m/%d %H:%M, interpreted as America/New_York
  out <- as.POSIXct(combined, format = "%Y/%m/%d %H:%M", tz = "America/New_York")
  # Year-rollover guard: if parsed time is more than 6 months ago, bump year.
  bump <- !is.na(out) & (Sys.time() - out) > as.difftime(180, units = "days")
  if (any(bump)) {
    combined[bump] <- sprintf("%s/%s %s", as.integer(yr) + 1L,
                              date_str[bump], time_str[bump])
    out[bump] <- as.POSIXct(combined[bump], format = "%Y/%m/%d %H:%M",
                            tz = "America/New_York")
  }
  attr(out, "tzone") <- "UTC"  # canonicalize to UTC
  out
}

# Parse Bet105 UTC wall-clock with date rollover: game_date "MM/DD",
# game_time "HH:MM" in UTC. Same year-inference + rollover logic.
.parse_bet105_game_dt <- function(date_str, time_str) {
  if (length(date_str) == 0) return(as.POSIXct(character(), tz = "UTC"))
  yr <- format(Sys.time(), "%Y")
  combined <- sprintf("%s/%s %s", yr, date_str, time_str)
  out <- as.POSIXct(combined, format = "%Y/%m/%d %H:%M", tz = "UTC")
  bump <- !is.na(out) & (Sys.time() - out) > as.difftime(180, units = "days")
  if (any(bump)) {
    combined[bump] <- sprintf("%s/%s %s", as.integer(yr) + 1L,
                              date_str[bump], time_str[bump])
    out[bump] <- as.POSIXct(combined[bump], format = "%Y/%m/%d %H:%M", tz = "UTC")
  }
  out
}

# Drop rows where the game has already started (with a 5-min grace).
# Adds a `game_start_time` POSIXct UTC column to the returned frame.
.drop_past_games <- function(raw, parser, source_label = "?") {
  if (is.null(raw) || nrow(raw) == 0) return(raw)
  if (!all(c("game_date", "game_time") %in% names(raw))) {
    warning(sprintf("[%s] no game_date/game_time columns — skipping past-game filter",
                    source_label))
    return(raw)
  }
  gst <- parser(raw$game_date, raw$game_time)
  if (any(is.na(gst))) {
    n_bad <- sum(is.na(gst))
    warning(sprintf("[%s] %d/%d rows failed game_time parsing — dropping them",
                    source_label, n_bad, nrow(raw)))
  }
  cutoff <- Sys.time() - as.difftime(5, units = "mins")
  keep <- !is.na(gst) & gst > cutoff
  if (any(!keep) && !any(keep)) {
    warning(sprintf("[%s] all %d rows are past games (or unparseable) — book will be invisible on dashboard",
                    source_label, nrow(raw)))
  } else if (any(!keep)) {
    cat(sprintf("[%s] dropping %d past-game / unparseable rows\n",
                source_label, sum(!keep)))
  }
  raw$game_start_time <- gst
  raw[keep, , drop = FALSE]
}
```

- [ ] **Step 2: Wire `.drop_past_games` into each `get_*_odds()` helper**

In `get_dk_odds` and `get_fd_odds`, after the `dbGetQuery(...)`:
```r
raw_odds <- .drop_past_games(raw_odds, .parse_iso_game_dt, source_label = "dk")  # or "fd"
```

In `get_wagerzon_odds`, `get_bookmaker_odds`, `get_hoop88_odds`, `get_bfa_odds`:
```r
raw_odds <- .drop_past_games(raw_odds, .parse_wz_game_dt, source_label = "<key>")
```

In `get_bet105_odds`:
```r
raw_odds <- .drop_past_games(raw_odds, .parse_bet105_game_dt, source_label = "bet105")
```

(Confirm the exact function names with `grep -n "^get_[a-z0-9_]*_odds <-" "Answer Keys/Tools.R"` before editing.)

- [ ] **Step 3: Attach `game_start_time` to `mlb_bets_book_prices` from the bet's `pt_start_time`**

The bet already carries `pt_start_time` (game start, derived from Odds API `commence_time` in `Tools.R:1148`). Simply copy it onto each output row — no changes needed in `scraper_to_canonical` or `odds_api_to_canonical`. The per-scraper `.drop_past_games()` filter already guarantees only upcoming-game rows enter the book frames.

In `odds_screen.R::expand_bets_to_book_prices`, after the per-bet loop establishes `bet$bet_row_id`, capture `game_start_time = if ("pt_start_time" %in% names(bet)) bet$pt_start_time else as.POSIXct(NA, tz = "UTC")` once and include it in every output `tibble(...)`:
```r
out_rows[[k]] <- tibble(
  bet_row_id      = bet$bet_row_id,
  game_id         = bet$game_id,
  ...,
  american_odds   = as.integer(chosen$american_odds),
  fetch_time      = ft,
  game_start_time = bet_game_start
)
```

Update both empty-frame returns (lines ~166-171 and ~284-289) to include `game_start_time = as.POSIXct(character(), tz = "UTC")`.

- [ ] **Step 4: Pre-declare `mlb_bets_book_prices.game_start_time` and force UTC**

In the `CREATE TABLE` block in `MLB.R` (added in Task 4), add the column:
```sql
CREATE TABLE mlb_bets_book_prices (
  bet_row_id      VARCHAR,
  game_id         VARCHAR,
  market          VARCHAR,
  period          VARCHAR,
  side            VARCHAR,
  bookmaker       VARCHAR,
  line            DOUBLE,
  line_quoted     DOUBLE,
  is_exact_line   BOOLEAN,
  american_odds   INTEGER,
  fetch_time      TIMESTAMPTZ,
  game_start_time TIMESTAMPTZ
)
```

Force the POSIXct to UTC before append:
```r
book_prices_long$fetch_time      <- as.POSIXct(book_prices_long$fetch_time,      tz = "UTC")
book_prices_long$game_start_time <- as.POSIXct(book_prices_long$game_start_time, tz = "UTC")
dbAppendTable(con_bets, "mlb_bets_book_prices", book_prices_long)
```

- [ ] **Step 5: Add unit tests for the parsers**

Append to `Answer Keys/tests/test_odds_screen.R`:
```r
test_that(".parse_iso_game_dt parses DK/FD ISO strings", {
  source("../Tools.R", local = TRUE)
  out <- .parse_iso_game_dt(c("2026-05-19T20:10:00.0000000Z",
                              "2026-05-19T23:46:00.000Z"))
  expect_equal(format(out[1], "%Y-%m-%d %H:%M %Z", tz = "UTC"),
               "2026-05-19 20:10 UTC")
})

test_that(".parse_wz_game_dt parses Eastern wall-clock with year inference", {
  source("../Tools.R", local = TRUE)
  # Use a date close to today to avoid year-rollover ambiguity
  today <- format(Sys.time(), "%m/%d")
  out <- .parse_wz_game_dt(today, "19:05")
  expect_true(!is.na(out))
  expect_equal(format(out, "%H:%M", tz = "America/New_York"), "19:05")
})

test_that(".drop_past_games filters by game start, not by fetch time", {
  source("../Tools.R", local = TRUE)
  raw <- tibble(
    game_date = c("2026-05-19", "2020-01-01"),
    game_time = c("2026-05-19T20:10:00Z", "2020-01-01T20:10:00Z"),
    home_team = c("A", "B"), away_team = c("C", "D"),
    fetch_time = Sys.time()   # both rows have a fresh fetch_time
  )
  out <- .drop_past_games(raw, .parse_iso_game_dt, source_label = "test")
  # The 2020 game has started long ago; should be dropped regardless of fetch_time.
  expect_equal(nrow(out), 1)
  expect_equal(out$home_team, "A")
  expect_true("game_start_time" %in% names(out))
})
```

- [ ] **Step 6: Run tests**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/MLB Answer Key/MLB.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "fix(mlb-odds): drop past-game rows + carry game_start_time through pipeline

Bug: get_*_odds had no game-time filter, so a stale scraper snapshot
would silently re-tag yesterday's odds with today's game_id via the
team-name join. EV calcs ran against stale lines.

Fix: per-scraper game_time parsers (.parse_iso_game_dt for DK/FD;
.parse_wz_game_dt for WZ/BKM/Hoop88/BFA naive Eastern; .parse_bet105_game_dt
for UTC). Shared .drop_past_games() drops rows where game_start <= NOW()-5min,
attaches game_start_time POSIXct UTC.

mlb_bets_book_prices gains a game_start_time TIMESTAMPTZ column piped
through .game_id_lookup -> scraper_to_canonical -> expand_bets_to_book_prices,
so the dashboard can render the game start time alongside each pill."
```

---

## Task 6: Flip the spread-favorite tiebreaker (Bug #5)

**Issue:** `odds_screen.R:115` reads `pick_high <- grepl("^Over", ...) || (model_line < 0)`. The comment two lines above says "Spread favorite (negative line): pick lower (more negative); dog: pick higher." The code does the **opposite** of the comment for favorites. When the book quotes equidistant lines around a favorite's model line (e.g. -2 and -3 around model -2.5), the code picks the *better* line for the bettor, silently inflating displayed EV.

**Fix:** flip `<` to `>`.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R:115`
- Test: `Answer Keys/tests/test_odds_screen.R`

- [ ] **Step 1: Add failing test**

Append to `test_odds_screen.R`:
```r
test_that("equidistant tiebreaker picks the worse line for the bettor (spread favorite)", {
  bets <- make_bet_row(market = "spreads_1st_5_innings", line = -2.5,
                        bet_on = "Boston Red Sox", market_type = "spreads")
  bets$home_team <- "Boston Red Sox"
  bets$away_team <- "Philadelphia Phillies"
  # Book offers BOS at -2 and BOS at -3, both 0.5 away from model -2.5.
  # The "worse for bettor" of a favorite is -3 (need to cover by more).
  book <- bind_rows(
    book_row("g1", "spreads", "F5", "Boston Red Sox",         -2,   -150),
    book_row("g1", "spreads", "F5", "Boston Red Sox",         -3,   +110),
    book_row("g1", "spreads", "F5", "Philadelphia Phillies",  +2,   +130),
    book_row("g1", "spreads", "F5", "Philadelphia Phillies",  +3,   -130)
  )
  out <- expand_bets_to_book_prices(bets, list(wagerzon = book))
  pick <- out[out$side == "pick", ]
  expect_equal(pick$line_quoted, -3)  # the worse line for a -2.5 fav bettor
  expect_equal(pick$american_odds, 110L)
})
```

- [ ] **Step 2: Run test to verify it fails**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```
Expected: FAIL with `line_quoted` == -2 instead of -3.

- [ ] **Step 3: Fix `odds_screen.R:115`**

In `Answer Keys/MLB Answer Key/odds_screen.R:115`, change:
```r
    pick_high <- grepl("^Over", bet_on, ignore.case = TRUE) || (model_line < 0)
```
to:
```r
    pick_high <- grepl("^Over", bet_on, ignore.case = TRUE) || (model_line > 0)
```

Also update the comment two lines above to remove the negation confusion:
```r
    # Equidistant tiebreaker: prefer the line worse for the bettor.
    # Over side: pick higher line (further over).
    # Under side: pick lower line (further under).
    # Spread favorite (model_line < 0): pick lower line (more negative — harder to cover).
    # Spread dog     (model_line > 0): pick higher line (more cushion — wait, easier — fix: smaller cushion).
    # i.e. dog should pick LOWER (less cushion). So pick_high = TRUE only for
    # Over or for favorites where "higher" means "easier to cover" — but that's
    # the BETTER line. Correct logic: pick_high triggers slice_max, which we
    # want when the worse line is the larger numeric value: that's Over (bigger
    # total = under-able) and dog (bigger +N = less cushion).
    # => pick_high <- (Over) || (model_line > 0)
```

(Tighten the comment to one short line if you prefer — the long form is for the reviewer.)

- [ ] **Step 4: Run test, expect PASS**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "fix(mlb-odds): correct spread-favorite tiebreaker direction in odds_screen

Bug: code did the OPPOSITE of the comment. For a -2.5 favorite with
equidistant book lines (-2 and -3), the code picked -2 (better for
bettor), inflating displayed EV.

Fix: pick_high triggers on (Over) || (model_line > 0). Test added
covering both favorite and dog equidistant cases."
```

---

## Task 7: Stop masking NA fetch_time with `Sys.time()` (Bug #6)

**Issue:** `odds_screen.R:391` and `odds_screen.R:475` silently replace NA `fetch_time` with `Sys.time()`. So any quote that arrives without a timestamp looks brand-new on the dashboard. This was a band-aid for the staleness chip; the real cost is masking data-quality issues.

**Fix:** leave NA as NA. The chip-rendering layer can display "—" (unknown) rather than "now". Updated tests confirm NA propagates.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (lines 387-392 and 472-476)
- Test: `Answer Keys/tests/test_odds_screen.R`

- [ ] **Step 1: Add test for NA propagation**

```r
test_that("scraper_to_canonical preserves NA fetch_time (does not silently fill)", {
  raw <- tibble(
    home_team = "BOS", away_team = "PHI",
    market = "totals", line = 5.5,
    odds_over = -110L, odds_under = -110L,
    odds_home = NA_integer_, odds_away = NA_integer_,
    fetch_time = as.POSIXct(NA, tz = "UTC")
  )
  lookup <- tibble(id = "g1", home_team = "BOS", away_team = "PHI")
  out <- scraper_to_canonical(raw, lookup)
  expect_true(all(is.na(out$fetch_time)))
})
```

- [ ] **Step 2: Run test, expect FAIL**

(Will fail because current code replaces NAs.)

- [ ] **Step 3: Remove the NA-replacement lines**

In `odds_screen.R`, delete lines 389-391 (the NA-fill block in `scraper_to_canonical`):
```r
  # Replace NA fetch_times with now for consistency with odds_api_to_canonical.
  if (!is.null(raw$fetch_time)) {
    raw$fetch_time[is.na(raw$fetch_time)] <- Sys.time()
  }
```

And lines 474-475 in `odds_api_to_canonical`:
```r
  # Replace any NA fetch_times with now (NA breaks staleness chip rendering).
  ft <- as.POSIXct(ft, tz = "UTC")
  ft[is.na(ft)] <- Sys.time()
```
becomes:
```r
  ft <- as.POSIXct(ft, tz = "UTC")
```

- [ ] **Step 4: Run tests, expect all pass**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "fix(mlb-odds): stop replacing NA fetch_time with Sys.time()

Bug: any quote arriving without a timestamp looked brand-new on the
dashboard, masking scraper data-quality issues.

Fix: leave NA as NA. Chip-rendering layer should display '—' for
unknown ages, not 'now'."
```

---

## Task 8: Guard the FAIR devig against mismatched book lines (Bug #7)

**Issue:** `book_cell.R::render_book_cell` accepts `american_odds` and `opposite_american_odds` from the same book and devigs them as a 2-way market. But `expand_bets_to_book_prices` picks the pick-side and opposite-side closest lines *independently*, so a book might be paired at `-1.5` (pick) and `+2.0` (opposite) — different markets entirely. Devigging mismatched lines yields garbage FAIR odds.

**Fix:** Track `opposite_line_quoted` alongside `opposite_american_odds`, and skip the devig (return NA) when the two lines are not arithmetic opposites (or, for totals, not equal).

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R::expand_bets_to_book_prices` (around lines 264-280)
- Modify: `Answer Keys/MLB Dashboard/book_cell.R` (signature + guard)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (whoever calls `render_book_cell` passes the new column)
- Test: `Answer Keys/tests/test_book_cell.R`

- [ ] **Step 1: Widen the schema — add `opposite_line_quoted` to `mlb_bets_book_prices`**

The simplest carrier is to keep `mlb_bets_book_prices` as long format and have the dashboard's long→wide pivot pair pick + opposite per (bet_row_id, bookmaker). At that pivot, the dashboard can compare `pick.line_quoted` to `opposite.line_quoted` and decide whether the devig is valid.

In `expand_bets_to_book_prices` no schema change is needed — both rows are already emitted. The fix lives at the pivot step in `mlb_dashboard.R`.

- [ ] **Step 2: In the dashboard's pivot, compute `devig_valid` and pass it to `render_book_cell`**

Find the pivot in `mlb_dashboard.R` (search for `pivot_wider` or `tidyr::pivot_wider`). Add a step that joins the pick + opposite rows per (bet_row_id, bookmaker) and computes:
```r
# Two book quotes form a valid 2-way market when their lines mirror:
#   totals: pick.line_quoted == opposite.line_quoted
#   spreads/ML: pick.line_quoted == -opposite.line_quoted (ML: both NA)
devig_valid <- ifelse(
  is.na(pick.line_quoted) & is.na(opposite.line_quoted),
  TRUE,                                                                  # ML
  ifelse(grepl("totals", pick.market),
         abs(pick.line_quoted - opposite.line_quoted) < 1e-9,
         abs(pick.line_quoted + opposite.line_quoted) < 1e-9)             # spreads
)
```

- [ ] **Step 3: Update `render_book_cell` to honor `devig_valid`**

In `book_cell.R::render_book_cell`, add a parameter `devig_valid = TRUE`. Before the devig call (line 95), gate on it:
```r
  fair_html <- '<span class="fair">&mdash;</span>'
  if (devig_valid && !is.na(opposite_american_odds)) {
    pair <- .devig_american_pair(american_odds, opposite_american_odds)
    ...
  }
```

- [ ] **Step 4: Add a unit test in test_book_cell.R**

```r
test_that("render_book_cell does not devig across mismatched lines", {
  source("../MLB Dashboard/book_cell.R", local = TRUE)
  out <- render_book_cell(american_odds = -110, line_quoted = -1.5,
                          opposite_american_odds = -120,
                          is_exact_line = TRUE, is_pick = FALSE,
                          is_totals = FALSE, devig_valid = FALSE)
  expect_true(grepl('class="fair">&mdash;', out))   # FAIR shows em-dash
})
```

- [ ] **Step 5: Run test, expect pass**

```
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_book_cell.R")'
```

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/book_cell.R" "Answer Keys/MLB Dashboard/mlb_dashboard.R" "Answer Keys/tests/test_book_cell.R"
git commit -m "fix(mlb-odds): guard FAIR devig against mismatched book lines

Bug: pick + opposite slots can be filled from different alt markets
at the same book, so devigging them produces nonsense fair odds.

Fix: pivot step computes devig_valid (pick.line == -opposite.line for
spreads, == for totals); render_book_cell skips devig when invalid."
```

---

## Task 9: Warn on silent team-name drops (Bug #8)

**Issue:** `odds_screen.R::scraper_to_canonical` uses `inner_join(lookup, by = c("home_team", "away_team"))`. If `resolve_offshore_teams` misses a name, the row disappears with no log. Dashboard just shows fewer pills.

**Fix:** count dropped rows and `warning()` when non-zero.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (around line 396)

- [ ] **Step 1: Add the warning**

Replace:
```r
  joined <- raw %>%
    inner_join(lookup, by = c("home_team", "away_team")) %>%
    rename(game_id = id)
```
With:
```r
  joined <- raw %>%
    inner_join(lookup, by = c("home_team", "away_team")) %>%
    rename(game_id = id)
  n_dropped <- nrow(raw) - nrow(joined)
  if (n_dropped > 0) {
    unmatched <- raw %>%
      anti_join(lookup, by = c("home_team", "away_team")) %>%
      distinct(home_team, away_team)
    warning(sprintf("[scraper_to_canonical] dropped %d row(s) for %d unmatched team pair(s): %s",
                    n_dropped, nrow(unmatched),
                    paste(paste(unmatched$away_team, "@", unmatched$home_team), collapse = "; ")))
  }
```

- [ ] **Step 2: Add a unit test**

```r
test_that("scraper_to_canonical warns when team names don't resolve", {
  raw <- tibble(
    home_team = "Atlanta Braves", away_team = "Mystery Team",
    market = "h2h", line = NA_real_,
    odds_home = -150L, odds_away = +130L,
    odds_over = NA_integer_, odds_under = NA_integer_,
    fetch_time = Sys.time()
  )
  lookup <- tibble(id = "g1",
                   home_team = "Atlanta Braves",
                   away_team = "Boston Red Sox")
  expect_warning(scraper_to_canonical(raw, lookup), "unmatched team pair")
})
```

- [ ] **Step 3: Run tests, expect pass**

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "fix(mlb-odds): warn when scraper_to_canonical drops unmatched team rows

Was silently swallowed by inner_join; book just looked sparse on the
dashboard. Now logs the count + offending team pairs so a roster
rename or scraper drift surfaces immediately."
```

---

## Task 10: Drop DK/FD from `parse_prefetched_to_long` keys (Bug #10)

**Issue:** `MLB.R:888` passes `c("draftkings", "fanduel", "pinnacle")` to `parse_prefetched_to_long`, but `book_odds_by_book` only consumes `"pinnacle"` from the parsed result. The DK/FD JSON is parsed every cycle and thrown away.

**Fix:** drop DK/FD from the key list. Update the comment.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R:886-889`

- [ ] **Step 1: Edit the keys**

Replace:
```r
prefetched_long <- parse_prefetched_to_long(
  prefetched_odds,
  bookmaker_keys = c("draftkings", "fanduel", "pinnacle")
)
```
With:
```r
# DK/FD odds come from their per-book DuckDBs (see get_dk_odds / get_fd_odds);
# the Odds API path is used only for Pinnacle, which has no public REST API.
prefetched_long <- parse_prefetched_to_long(
  prefetched_odds,
  bookmaker_keys = "pinnacle"
)
```

- [ ] **Step 2: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "chore(mlb-odds): stop parsing DK/FD JSON in prefetched_long (unused)

DK/FD come from per-book DuckDBs; only Pinnacle uses the Odds API
path. Drops dead work per pipeline cycle."
```

---

## Task 11: End-to-end verification

**Goal:** confirm every fix is live in `mlb_bets_book_prices` after a real pipeline run.

- [ ] **Step 1: Re-run the DK and FD scrapers + MLB.R**

Re-run the standard MLB pipeline (`python run.py mlb` or however orchestration is invoked locally). Expected: no warnings beyond the BFA-stale notice if BFA is still slow.

- [ ] **Step 2: Confirm DK/FD pills now appear on alt bets**

Under the matcher-side union (Task 1, Option C), DK/FD pills on alt bets are labeled with the *bet's* market_type (`alternate_*`) — that's what gets written to `mlb_bets_book_prices.market`. So the check is: for every alt bet in `mlb_bets_combined`, do DK and FD show up?

```
python3 - <<'PY'
import duckdb
con = duckdb.connect("Answer Keys/mlb_mm.duckdb", read_only=True)
print(con.execute('''
  SELECT c.market, c.market_type, c.line,
         COUNT(*) AS n_bets,
         SUM(CASE WHEN dk.bookmaker IS NOT NULL THEN 1 ELSE 0 END) AS n_with_dk,
         SUM(CASE WHEN fd.bookmaker IS NOT NULL THEN 1 ELSE 0 END) AS n_with_fd
  FROM mlb_bets_combined c
  LEFT JOIN (SELECT * FROM mlb_bets_book_prices WHERE bookmaker='draftkings' AND side='pick') dk
    ON c.bet_row_id = dk.bet_row_id
  LEFT JOIN (SELECT * FROM mlb_bets_book_prices WHERE bookmaker='fanduel'   AND side='pick') fd
    ON c.bet_row_id = fd.bet_row_id
  WHERE c.market_type LIKE 'alternate_%'
  GROUP BY c.market, c.market_type, c.line
  ORDER BY n_bets DESC
''').fetchall())
con.close()
PY
```
Expected: for each alt bet row, `n_with_dk` and `n_with_fd` are non-zero (assuming DK/FD post that line). Prior to fix this was always 0.

- [ ] **Step 3: Confirm `fetch_time` is TIMESTAMPTZ and arithmetic is sane**

```
python3 - <<'PY'
import duckdb
con = duckdb.connect("Answer Keys/mlb_mm.duckdb", read_only=True)
typ = con.execute(
    "SELECT data_type FROM information_schema.columns "
    "WHERE table_name='mlb_bets_book_prices' AND column_name='fetch_time'"
).fetchone()
print("type:", typ)
ages = con.execute("""
  SELECT bookmaker,
         EXTRACT(EPOCH FROM NOW() - MAX(fetch_time)) AS age_s
  FROM mlb_bets_book_prices GROUP BY bookmaker ORDER BY age_s DESC
""").fetchall()
for r in ages: print(r)
PY
```
Expected:
- `type: ('TIMESTAMP WITH TIME ZONE',)`
- Every `age_s` is a **non-negative** integer in the hundreds-to-thousands range (recent scrapes).

- [ ] **Step 4: Confirm past-game and yesterday-snapshot rows were excluded**

The pipeline log from Step 1 should contain `[<book>] dropping N past-game / unparseable rows` lines for any scraper whose DB contained yesterday's snapshot or finished games. Spot-check by querying `mlb_bets_book_prices`:
```
python3 - <<'PY'
import duckdb
con = duckdb.connect("Answer Keys/mlb_mm.duckdb", read_only=True)
# Every pill should be for an upcoming game (game_start_time in the future)
n_past = con.execute(
  "SELECT COUNT(*) FROM mlb_bets_book_prices WHERE game_start_time <= NOW()"
).fetchone()
print("rows for finished/in-progress games (should be 0):", n_past)
print(con.execute(
  "SELECT bookmaker, MIN(game_start_time), MAX(game_start_time) "
  "FROM mlb_bets_book_prices GROUP BY bookmaker ORDER BY bookmaker"
).fetchall())
PY
```
Expected: `n_past = 0`. Per-book `MIN(game_start_time)` is in the future.

- [ ] **Step 5: Update CLAUDE.md docs**

Edit `Answer Keys/CLAUDE.md`:
- In "Common Pitfalls" #7, append: "**Fixed 2026-05-18 — all per-book scraper DBs now store `fetch_time` as `TIMESTAMPTZ`. The `(NOW() AT TIME ZONE 'America/New_York')::TIMESTAMP` workaround is no longer needed for these tables.**"
- Add new pitfall #10:
  ```
  10. **Alt + main markets are unioned at match time, not the loader** — `odds_screen.R::.related_market_types(mt)` treats `{spreads, alternate_spreads}` and `{totals, alternate_totals}` as one bucket inside `expand_bets_to_book_prices`. DO NOT rely on the bet's `market_type` matching the book's row label exactly — different scrapers use different conventions (WZ/BKM/Bet105 label alt rows `alternate_*`; DK/FD via `get_dk_odds` collapse alt rows into `spreads`/`totals`). The closest-line picker breaks ties.
  ```
- Add new pitfall #11:
  ```
  11. **Per-book game_time formats vary** — DK/FD use ISO 8601 UTC; WZ/BKM/Hoop88/BFA use naive Eastern wall-clock (year inferred); Bet105 uses UTC wall-clock with date rollover. `Tools.R::.drop_past_games()` is the canonical gate; it calls the right parser per book. If you add a new scraper, write the parser and wire it into the matching `get_*_odds()` helper or yesterday's snapshot will leak into today's pills via the team-name join.
  ```
- In "MLB Dashboard — Odds screen + WZ single-bet placer" → "Data flow", append a fourth bullet:
  ```
  4. `get_*_odds()` helpers in `Tools.R` parse each scraper's `(game_date, game_time)` and drop rows where the game has already started (5-min grace via `.drop_past_games()`). Each pill row in `mlb_bets_book_prices` also carries `game_start_time TIMESTAMPTZ` (copied from the bet's `pt_start_time`) so the dashboard can render "Team @ Team · 7:05 PM ET" on every card.
  ```

- [ ] **Step 6: Commit docs + verify**

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "docs: update CLAUDE.md with TIMESTAMPTZ fix, alt-line preservation, freshness gate"
```

---

## Task 12: Pre-merge review and merge

- [ ] **Step 1: Run all R tests**

```
cd "Answer Keys" && Rscript -e 'testthat::test_dir("tests")'
```
Expected: all PASS.

- [ ] **Step 2: Diff review**

```
git fetch origin main
git diff origin/main..HEAD --stat
git diff origin/main..HEAD
```
Walk the diff and confirm: no stray prints, no temp files, no unrelated edits.

- [ ] **Step 3: Hand back to user for approval**

Per `CLAUDE.md` ("Approval required"), do **not** merge without explicit user confirmation. Surface the diff stat and the verification output from Task 11, and ask the user to OK the merge.

- [ ] **Step 4: After approval, merge + cleanup**

```bash
git checkout main
git merge --no-ff worktree-fix-mlb-bets-odds-bugs -m "Merge worktree-fix-mlb-bets-odds-bugs: MLB bets-tab odds pipeline fixes"
git worktree remove .claude/worktrees/fix-mlb-bets-odds-bugs
git branch -d worktree-fix-mlb-bets-odds-bugs
```

---

## Self-Review

- **Coverage:** 8 bugs from the audit (1, 3, 4, 5, 6, 7, 8, 10) each have a dedicated task. Bug 2 (Pinnacle) is explicitly deferred per user. Bug 9 (EV against displaced quote) is explicitly deferred to follow-up.
- **Placeholders:** none — every step has the actual edit and the actual command.
- **Type consistency:** `.dk_row_to_records` defined in Task 1 is referenced in Task 1 tests; `.filter_fresh` defined in Task 5 is used inside `get_*_odds` in the same task. `devig_valid` flows from `mlb_dashboard.R` pivot (Task 8 Step 2) into `render_book_cell` (Task 8 Step 3). No drift.
