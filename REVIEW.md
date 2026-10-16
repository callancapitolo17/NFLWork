# REVIEW.md

What a reviewer of this repo looks for. The `pre-pr-review` skill (user-level, `~/.claude/skills/`) runs the procedure and hands this file to every reviewer agent.

## Checks to run

Run what the change touches, from the repo root, and paste the last lines of output.

| Touched | Run |
|---|---|
| A Python package with a `tests/` dir | `python -m pytest <package>/tests -q` |
| `unabated_ticket/` JS | `cd unabated_ticket && npm test` |
| A scraper, or anything writing `game_start_time` | `python tests/timezone_parity_test.py` |
| R under `Answer Keys/` | the script's own run, with its stdout |

## High-risk paths

Always the High tier (full review, `/security-review`, Cal reads the diff): the "Order placement and quoting" and "Bets service" paths below; `bet_logger/`; `.github/`, launchd plists, `deploy/`, any `run.py` or cron entry.

## Reviewer focuses

One per agent, 4 agents:

1. **Data**: the Scrapers, DuckDB writes, Dashboards, and Schedulers failures below.
2. **Money**: the Pricing, Order placement, Bets service, and Extension failures below.
3. **Reuse and cut**: for each new function, file, field, constant, or setting, find an existing one in this repo or the library; then comments, docstrings, and code the change does not need.
4. **Tests and edge cases**: tests that assert one field instead of the full output; empty tables, no games today, first run, off-season, timezone boundaries; PR text or README out of sync with the code.

## Rules for any finding

- Cite `file:line` and the input or state that produces the wrong result. A finding with no concrete failure is not a finding.
- Correctness only. No style, naming, or cleanup findings unless CLAUDE.md names the rule.
- At most 3 nits per review.
- A line older than the branch (`git blame` before the merge-base) is **Pre-existing**: report it, do not block on it.
- The summary starts with a count: "2 to fix, 1 pre-existing". When nothing needs fixing, say "No changes needed."

## What failure looks like, by area

**Scrapers** (`*_odds/`, `mlb_sgp/scraper_*`, anything that writes a per-book DuckDB). Failure: a wrong or stale row stored as fresh.
- `game_start_time` naive or not UTC (`tests/timezone_parity_test.py` is the gate).
- A team name the canonical join will not match.
- Duplicate rows for one (game, market, book, fetch).
- A line stored on the wrong side or with the wrong sign.
- An empty or failed fetch (expired cookies, HTTP 200 with 0 rows) that leaves the old snapshot looking current.
- Odds API `commence_time` is a character column in R; it must be parsed with `ymd_hms(..., tz = "UTC")` before any filter.

**Pricing and fair value** (`kalshi_common/fair_value*`, `unabated_edge/pricing*`, `nfl_specials/pricing.py`, `Answer Keys/Tools.R`, any devig or Kelly code). Failure: a silent wrong number.
- American, decimal, and implied probability mixed in one expression.
- Vig left in, or removed twice.
- A spread or total sign flipped between books.
- `NA`/`None` that turns into 0 or is dropped from a median.
- Kelly at full stake instead of the fractional stake the caller passes.
- One book counted twice in a consensus (Novig relays vendor prices; see `router.drop_relayed_novig`).
- A probability outside (0, 1) that nothing rejects.

**Order placement and quoting** (`unabated_ticket/server/taker.js`, `unabated_ticket/bets_service/taker_orders.py`, `unabated_ticket/bets_service/bfa_teaser.py`, `wagerzon_odds/single_placer.py`, `nfl_specials/`, `kalshi_mlb_mm/`, `kalshi_mlb_rfq/`, `unabated_edge/maker/`). Failure: it bets when it should not, at the wrong size, or twice.
- A retry after an ambiguous reply (timeout, 5xx, lost response).
- An order that can rest when it was meant to be immediate-or-cancel.
- A bet on a started game, or on a decision older than its freshness window.
- An exposure, daily, or per-game cap the new path does not pass through.
- An exception caught and logged where the caller expects a halt (`submit_quote` once swallowed every error and the bot placed 0 quotes for two weeks).
- A stop file or dead-man switch the new path does not check.

**Bets service and bet history** (`unabated_ticket/bets_service/`). Failure: a bet goes missing or doubles.
- A dedup key that changes when the venue reformats a field ($0 BetOnline bets never deduped).
- A settlement that overwrites an open row's fill data.
- A source that returns 0 rows on HTTP 200 and is read as "no bets".
- Settle day computed in UTC where the tracker uses Pacific.

**Extension and phone page** (`unabated_ticket/extension/`, `unabated_ticket/server/phone/`). Failure: the screen shows a number the math did not produce.
- A cached value rendered as live after its source stopped updating.
- A stake shown above the liquidity cap or the Kelly stake.

**Dashboards and monitors** (`kalshi_mlb_monitor/`, `nfl_draft/`, `Answer Keys/` dashboard code). Failure: stale shown as fresh, or a reader blocks the writer.
- A read that opens the live bot's DuckDB without the retry and poll guard.
- A freshness label computed from file mtime instead of `max(fetch_time)`.

**DuckDB writes, anywhere.** Failure: data lost or duplicated.
- The docstring does not say which file and table, and whether it appends, upserts, or replaces.
- `connect()` inside a per-item loop (17 ms each; three incidents).
- A connection not closed on the error path (`on.exit(dbDisconnect(...))` in R; `with` or `try/finally` in Python).
- `SELECT *`.
- A `CREATE TABLE AS` that drops a primary key.

**Schedulers, launchd, cron, `run.py`.** Failure: it never runs, or runs twice.
- A `cwd` or `PYTHONPATH` change that breaks an import the way it is launched (DK/FD singles died at import for a month because `run.py` ran them with `cwd=mlb_sgp`).
- An env var read in Python that only exists in `~/.Renviron`.
- A schedule in UTC where the slate is Pacific.

**Any change.**
- A flag, function, import, field, or config value nothing reads.
- A log or file that grows without rotation.
- A secret, API key, or session token in a log line or printed output.

## Always a human

Cal reads the diff and the test output before merging anything that places orders, sizes stakes, or writes a bets table. An agent's approval never replaces that read.

## Skip

`docs/`, `*.log`, `**/lib/`, `.claude/worktrees/`, generated `report.html`.

## Score the reviewer

`review_evals/cases.json` holds real bugs that reached `main`, each with the diff that introduced it and the finding a reviewer should have made. After any change to this file, the skill, or the reviewer model, run the cases blind (the skill's "Score the reviewer" section) and record with `python review_evals/score.py`. Keep the change only if the caught count does not drop.
