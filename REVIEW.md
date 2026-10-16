# REVIEW.md

What a reviewer of this repo looks for. The user-level `pre-pr-review` skill runs the procedure (evidence rules, tags, count-first report, round cap) and hands this file to every reviewer agent.

## Diff base

Local `main`, which often runs ahead of `origin/main`: `git diff $(git merge-base main HEAD)..HEAD`.

## Checks to run

Run what the change touches and paste the last lines of output.

| Touched | Run |
|---|---|
| `unabated_ticket/` | `./unabated_ticket/check.sh` (ESLint, node tests, bets-service pytest; `npm install` in `unabated_ticket/` once per checkout) |
| Any other Python package with a `tests/` dir | `python3 -m pytest <package>/tests -q` |
| R under `Answer Keys/` | from `Answer Keys/`: `Rscript -e 'testthat::test_file("tests/<file>")'` for each test file that covers a changed function |
| A scraper, or anything writing `game_start_time` | `python3 tests/timezone_parity_test.py`. In a worktree it finds no per-book DBs and prints PASS having checked nothing: copy the book's DB in (never symlink), run the changed scraper, then run the test. A PASS with "DB missing" lines is not evidence. |

## High tier

These paths always get the High tier: full review, `/security-review`, and Cal reads the diff and the check output before merging.

- Order placement and quoting: `unabated_ticket/server/taker.js`, `unabated_ticket/server/runner.js`, any `unabated_ticket/extension/` module `runner.js` requires (they size the live taker's stake), `unabated_ticket/bets_service/taker_orders.py`, `unabated_ticket/bets_service/bfa_teaser.py`, `wagerzon_odds/single_placer.py`, `wagerzon_odds/parlay_placer.py`, the `/api/place-*` routes in `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`, `nfl_specials/`, `kalshi_mlb_mm/`, `kalshi_mlb_rfq/`, `kalshi_rfi/`, `kalshi_mm/`, `unabated_edge/maker/`, `kalshi_common/`.
- Bets and P&L: `unabated_ticket/bets_service/`, `bet_logger/`, `unabated_ticket/server/tracker/`, `unabated_ticket/server/closefair.js`.
- Schedulers and deploy: `*.plist`, `unabated_ticket/deploy/`, `Answer Keys/run.py`, any cron entry.

Low tier is Markdown and docs only. Everything else is Medium.

## Reviewer focuses

One per agent, 4 agents:

1. **Data**: the Scrapers, DuckDB writes, Dashboards, and Schedulers failures below.
2. **Money**: the Pricing, Order placement, Bets and P&L, and Extension failures below.
3. **Reuse and cut**: for each new function, file, field, constant, or setting, find an existing one in this repo or the library; then comments, docstrings, and code the change does not need.
4. **Tests and edge cases**: tests that assert one field instead of the full output; empty tables, no games today, first run, off-season, timezone boundaries; PR text or README out of sync with the code.

## What failure looks like, by area

**Scrapers** (`*_odds/`, `mlb_sgp/scraper_*`, anything that writes a per-book DuckDB). Failure: a wrong or stale row stored as fresh.
- `game_start_time` naive or not UTC (Bookmaker once stored Pacific times as naive text and the odds screen dropped about half its MLB rows).
- A game matched across sources on the wrong key: team names alone, a doubleheader, a date in the wrong timezone.
- Duplicate rows for one (game, market, book, fetch).
- A line stored on the wrong side or with the wrong sign.
- An empty or failed fetch (expired cookies, HTTP 200 with 0 rows) that leaves the old snapshot looking current.
- A partial scrape (cut off mid-run) written as if complete.
- A Playwright browser not closed on an error path.
- A request rate that trips a book's limit (ProphetX 403s at our volume, DraftKings' Akamai burst gate).
- Odds API `commence_time` is a character column in R; it must be parsed with `ymd_hms(..., tz = "UTC")` before any filter.

**Pricing and fair value** (`kalshi_common/fair_value*`, `unabated_edge/pricing*`, `nfl_specials/pricing.py`, `Answer Keys/Tools.R`, any devig or Kelly code). Failure: a silent wrong number.
- American, decimal, and implied probability mixed in one expression.
- Vig left in, or removed twice.
- A spread or total sign flipped between books.
- `NA`/`None` that turns into 0 or is dropped from a median.
- Kelly at full stake instead of the fractional stake the caller passes.
- One book counted twice in a consensus (Novig relays vendor prices; see `router.drop_relayed_novig`).
- A probability outside (0, 1) that nothing rejects.

**Order placement and quoting** (the High-tier order paths above). Failure: it bets when it should not, at the wrong size, or twice.
- A retry after an ambiguous reply (timeout, 5xx, lost response).
- An order that can rest when it was meant to be immediate-or-cancel.
- A bet on a started game, or on a decision or snapshot older than its freshness window (`TAKER_MAX_SNAPSHOT_AGE_MS` and `TAKER_MAX_HELD_AGE_MS` in `runner.js` are such windows).
- An exposure, daily, or per-game cap the new path does not pass through.
- An error reply turned into a silent `None` where the caller needs a halt or a logged reason (`kalshi_mlb_mm/quote_gateway.py` `submit_quote` still does this; the bot placed 0 quotes for two weeks).
- A reply status the code does not map, leaving a row open forever (a 404 on a quote once matched no status and held a slot under `MAX_OPEN_QUOTES` until quoting stopped).
- A stop file or dead-man switch (`TAKER_STOP_PATH`, `MAKER_LIVE_ACK`, `RFI_LIVE_ACK`) the new path does not check.

**Bets and P&L** (the High-tier bets paths above). Failure: a bet goes missing, doubles, or lands on the wrong day.
- A dedup key built differently on the two sides it compares (`bet_logger/sheets.py` keys a $0 bet as `''` on one side and `'0.00'` on the other, so it re-appends every run).
- A settlement that overwrites an open row's fill data.
- A source that returns 0 rows on HTTP 200 and is read as "no bets".
- An open or unsettled bet read as settled.
- A venue clock or format fixed from one captured sample (BFA's open-bets offset was set to +7 h off one February capture and was wrong by September).
- Settle day computed in UTC where the tracker uses Pacific (`server/tracker/trackerstats.js`).

**Extension and phone page** (`unabated_ticket/extension/`, `unabated_ticket/server/phone/`). Failure: the screen shows a number the math did not produce, or stops updating.
- A cached value rendered as live after its source stopped updating.
- A suggested stake above the liquidity cap or the Kelly stake. (The Ticket shows the uncapped stake under the cap on purpose.)
- An awaited request with no timeout behind a busy flag: one hung fetch stops every later poll (Bet105 stalled 17 hours this way).

**Dashboards and monitors** (`kalshi_mlb_monitor/`, `nfl_draft/`, `Answer Keys/` dashboard code outside the `/api/place-*` routes). Failure: stale shown as fresh, or a reader blocks the writer.
- A read that opens the live bot's DuckDB without the retry and poll guard.
- A freshness label computed from file mtime instead of `max(fetch_time)`.
- A re-render or re-query on every tick when nothing changed.

**DuckDB writes, anywhere.** Failure: data lost or duplicated.
- The docstring does not say which file and table, and whether it appends, upserts, or replaces.
- `connect()` inside a per-item loop (17 ms each; three incidents, the worst in `kalshi_mlb_mm` discovery).
- A connection not closed on the error path (`on.exit(dbDisconnect(...))` in R; `with` or `try/finally` in Python), or a WAL or lock file left behind on a crash.
- `SELECT *`.
- A `CREATE TABLE AS` that drops a primary key.

**Schedulers, launchd, cron, `run.py`.** Failure: it never runs, or runs twice.
- A `cwd` or `PYTHONPATH` change that breaks an import the way it is launched (DK/FD singles died at import because `run.py` runs them with `cwd=mlb_sgp`).
- An env var read in Python that only exists in `~/.Renviron`.
- A schedule in UTC where the slate is Pacific.

**Any change.**
- A flag, function, import, field, or config value nothing reads.
- A log or file that grows without rotation.
- A secret, API key, or session token in a log line or printed output.

## Skip

`docs/`, `*.log`, `**/lib/`, `.claude/worktrees/`, generated `report.html`.
