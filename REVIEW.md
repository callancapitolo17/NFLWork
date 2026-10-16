# REVIEW.md

What a reviewer of this repo looks for. The user-level `pre-pr-review` skill runs the procedure (evidence rules, tags, count-first report, round cap) and hands this file to every reviewer agent.

## Diff base

Local `main`, which often runs ahead of `origin/main`: `git diff $(git merge-base main HEAD)..HEAD`.

## Checks to run

Run what the change touches and paste the last lines of output. Name any touched file no row covers; silence is not a pass. Quote paths with spaces.

| Touched | Run |
|---|---|
| `unabated_ticket/` | `./unabated_ticket/check.sh` (ESLint, node tests, bets-service pytest; `npm install` in `unabated_ticket/` once per checkout) |
| A Python module | Its package's `tests/` dir, any `test_*.py` beside it (`wagerzon_odds/test_parlay_placer.py` sits outside `tests/`), the root `tests/kalshi_mlb_rfq/` for `kalshi_mlb_rfq` and `kalshi_common`, and every test file that imports the module: `grep -rl "<module>" --include='test_*.py' . \| grep -v worktrees`, then `python3 -m pytest <those paths> -q`. `"Answer Keys/MLB Dashboard/tests"` covers `single_placer` and the place routes. |
| R under `Answer Keys/` | from `Answer Keys/`: `Rscript -e 'testthat::test_file("tests/<file>")'` for each test that covers a changed function, or the run line in the test's header when it has one. testthat runs from `tests/`, so `Tools.R` never finds `src/sampling.cpp` and only the pure-R sampler is tested. |
| A `.plist` or crontab | `plutil -lint <file>`; `crontab` files by eye |
| A scraper, or anything writing `game_start_time` | `python3 tests/timezone_parity_test.py`, from a checkout that has `bet_logger/.env` and the book's cookie file (gitignored, so not in a worktree), after running the changed scraper. It reads only the `mlb_odds` table of the 8 books in its `SCRAPERS` dict; for any other writer, show the column's type and one written row. A PASS next to `DB missing`, `table missing`, or `0 game(s)` lines checked nothing. |

## High tier

These paths always get the High tier: full review, `/security-review`, and Cal reads the diff and the check output before merging.

- Order placement and quoting: `unabated_ticket/server/taker.js` and the `extension/` modules it reaches through `require`, directly or not (`kelly.js` and `condkelly.js` size the stake); `unabated_ticket/server/runner.js`; `unabated_ticket/server/teasers_payload.js`; `unabated_ticket/bets_service/taker_orders.py`; `unabated_ticket/bets_service/bfa_teaser.py`; `wagerzon_odds/single_placer.py`, `parlay_placer.py`, `parlay_pricer.py`, `wagerzon_balance.py`, `wagerzon_accounts.py`; `bet_placer/`; the `/api/place-*` and `/api/auto-place` routes and the helpers they call in `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` and `Answer Keys/CBB Dashboard/cbb_dashboard_server.py`; `nfl_specials/`; `kalshi_mlb_mm/`; `kalshi_mlb_rfq/`; `kalshi_rfi/`; `kalshi_mm/`; `unabated_edge/` (its caps and kill switch live in `config.py`, `risk.py`, `runner.py`, not only `maker/`); `kalshi_common/`.
- Bets and P&L: `unabated_ticket/bets_service/`, `bet_logger/`, `unabated_ticket/server/tracker/`, `unabated_ticket/server/closefair.js`.
- Schedulers, deploy, and the checks themselves: `*.plist`, `nfl_draft/crontab.*`, `unabated_ticket/deploy/`, `Answer Keys/run.py`, `coverage_audit/AGENT_PLAYBOOK.md`, `unabated_ticket/check.sh`.

Low tier is Markdown and docs only. Everything else is Medium.

## Reviewer focuses

One per agent, 4 agents:

1. **Data**: the Scrapers, DuckDB writes, Dashboards, and Schedulers failures below.
2. **Money**: the Pricing, Order placement, Bets and P&L, and Extension failures below.
3. **Reuse and cut**: for each new function, file, field, constant, or setting, find an existing one in this repo or the library; then comments, docstrings, and code the change does not need.
4. **Tests and edge cases**: tests that assert one field instead of the full output; empty tables, no games today, first run, off-season, timezone boundaries; PR text or README out of sync with the code.

## What failure looks like, by area

**Scrapers** (`*_odds/`, `mlb_sgp/scraper_*`, `nfl_draft/scrapers/`, `unabated_edge/runner.py`, anything that writes odds to a DuckDB). Failure: a wrong or stale row stored as fresh.
- `game_start_time` naive or not UTC (Bookmaker once stored Pacific times as naive text and the odds screen dropped about half its MLB rows).
- A team name the canonical join will not match, so the book's row silently drops.
- A game matched across sources on the wrong key: team names alone, a doubleheader, a date in the wrong timezone.
- A new `get_*_odds()` in `Tools.R` that does not pass its rows through `.drop_past_games()`, or a change that lets that gate return rows when `game_start_time` is missing. It must fail closed: on 2026-09-01 about 100 of 110 Wagerzon rows on the odds screen were a two-month-old slate.
- Duplicate rows for one (game, market, book, fetch).
- A line stored on the wrong side or with the wrong sign.
- An empty or failed fetch (expired cookies, HTTP 200 with 0 rows) that leaves the old snapshot looking current.
- A partial scrape (cut off mid-run) written as if complete.
- A snapshot written after kickoff where the newest snapshot is read as the close (`unabated_edge` gates `line_snapshots` on pre-kickoff).
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
- An exposure, daily, or per-game cap the new path does not pass through, or a cap's default raised in a config file.
- A side or stake built by a helper that the route calls (`_build_spec_from_row` picks the Wagerzon side; `compute_exact_payouts` writes `exact_wager`).
- An error reply turned into a silent `None` where the caller needs a halt or a logged reason (`kalshi_mlb_mm/quote_gateway.py` `submit_quote` still does this; the bot placed 0 quotes for two weeks).
- A reply status the code does not map, leaving a row open forever (a 404 on a quote once matched no status and held a slot under `MAX_OPEN_QUOTES` until quoting stopped).
- A stop file or dead-man switch (`TAKER_STOP_PATH`, `MAKER_LIVE_ACK`, `RFI_LIVE_ACK`, the `unabated_edge` kill file) the new path does not check.

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

**Dashboards and monitors** (`kalshi_mlb_monitor/`, `kalshi_draft/app.py`, `Answer Keys/` dashboard code outside the placement routes). Failure: stale shown as fresh, or a reader blocks the writer.
- A read that opens the live bot's DuckDB without the retry and poll guard.
- A freshness label computed from file mtime instead of the newest row timestamp (`max(fetch_time)`, `max(decision_ts)`, `max(scraped_at)`).
- A re-render or re-query on every tick when nothing changed.

**DuckDB writes, anywhere.** Failure: data lost or duplicated.
- `connect()`, or a query, per item inside a loop (a connect is 17 ms; three incidents, the worst in `kalshi_mlb_mm` discovery).
- A connection not closed on the error path (`on.exit(dbDisconnect(...))` in R; `with` or `try/finally` in Python), or a WAL or lock file left behind on a crash.
- `SELECT *` or a positional `row[n]` read: a new column shifts every later field.
- A `CREATE TABLE AS` that drops a primary key.

**Schedulers, launchd, cron, `run.py`, the coverage-audit scheduled task.** Failure: it never runs, or runs twice.
- An import edited in a module that is launched as a script: script mode puts the script's folder first on `sys.path`, so `from mlb_sgp._shared import` in `dk_client.py` killed both singles scrapers at import. Still broken on `main`; the DK and FD pills have been stale since 2026-08-17.
- An env var read in Python that only exists in `~/.Renviron`.
- A schedule in UTC where the slate is Pacific.

**Any change.**
- A flag, function, import, field, or config value nothing reads.
- A log or file that grows without rotation.
- A secret, API key, or session token in a log line or printed output.

## Skip

`docs/`, `*.log`, `**/lib/`, `.claude/worktrees/`, generated `report.html`.
