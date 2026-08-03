# MLB SGP Odds Scrapers

Fetch Same Game Parlay (SGP) odds from multiple books (DraftKings, FanDuel, ProphetX, Novig, BetMGM, Caesars) for MLB correlated parlay edge finding. All write to the same `mlb_sgp_odds` table so the scanner can shop combos across books. See "Additional SGP books" below for BetMGM/Caesars/bet365.

Each book prices two FG combo families, both stored as 4-cell devig-able grids:
**spread × total** (`"Home Spread + Over"` …, `spread_line` set) and
**moneyline × total** (`"Home ML + Over"` …, `spread_line` = **NULL** — there is
no spread leg). Both are correlation-priced by the book and devigged the same way
(`devig_book`, n-way probit). Moneyline is FG-only across all 6 books.

## Quick Start

```bash
cd mlb_sgp
source venv/bin/activate
python scraper_draftkings_sgp.py           # all games, ~30s
python scraper_fanduel_sgp.py              # all games, ~1s
python scraper_draftkings_sgp.py --verbose
python scraper_fanduel_sgp.py --verbose
```

Requirements: `pip install curl_cffi duckdb` and the MLB pipeline must have run (`mlb_parlay_opportunities` table populated).

## How It Works

Pure REST API — no browser, no Chrome, no clicking. Uses `curl_cffi` with Chrome TLS impersonation to bypass DraftKings' Akamai bot protection.

```
DK League API          → event list (team names, event IDs)
DK Event Markets API   → main market number (Run Line ID)
DK SGP Parlays API     → ALL selection IDs (2MB response, main + alt lines)
Match Wagerzon total   → exact Over/Under selection ID
DK calculateBets       → correlation-adjusted SGP trueOdds
                       → mlb_sgp_odds table in DuckDB
```

## Singles scrapers

Two scrapers fetch single-leg odds for the MLB Dashboard bets tab,
replacing the Odds API for DraftKings and FanDuel pill data:

- `scraper_draftkings_singles.py` → writes `../dk_odds/dk.duckdb::mlb_odds`
- `scraper_fanduel_singles.py`    → writes `../fd_odds/fd.duckdb::mlb_odds`

Both use the shared client classes (`dk_client.py`, `fd_client.py`) — same
curl_cffi sessions, same Akamai bypass, same event-discovery code path as
the SGP scrapers. No second auth path, no duplicate rate-limit budget.

Output schema matches the wagerzon offshore convention (18-column wide
`mlb_odds` table). MLB.R consumes via `get_dk_odds()` / `get_fd_odds()`
(in `Tools.R`) → `scraper_to_canonical()` → `book_odds_by_book`.

**FanDuel market selection (2026-05-23).** FD's `event-page` endpoint returns
a *different* market slice per `tab`, and no single tab has everything: the
**default** tab (`tab=`) carries F7/F3 + per-inning lines, while the
**same-game-parlay** tab carries FG-alts + all of F5. `scraper_fanduel_singles`
therefore fetches **both** tabs (`FD_TABS`) via `fd_client.fetch_event_page()`
and merges them (dedup by market/runner id). Market selection is a **keyword
classifier** (`classify_market`, mirroring DK) — period + market-type detection
with FD-specific junk exclusions (`parlay`, `listed`, `bands`, `tri-bet`,
`specials`, per-inning, 3-way `Result`) and event-team team-total filtering —
**not** the old exact-name whitelist, which silently missed any market FD posts
under an unlisted name (that bug is why F7 totals were absent). It now covers
FG/F5/F7/F3 main + alt wherever FD posts them and auto-picks-up new markets.

**Parse-failure tripwire (2026-06-10; logging since #36).** Every scrape logs
`fd_singles: alt parse: alternate_spreads N/M, alternate_totals N/M`
(parsed/seen) at INFO. If FD hands us alt runners and **none** parse
(`M > 0, N == 0`) it logs a WARNING — that's the signature of an FD
name-format change (the FG alt-total parens bug shipped silently in exactly
this mode). A day where FD posts no alts at all stays quiet. Grep the runner
log for `PARSE TRIPWIRE` when alt pills go missing from the dashboard; the
wording is shared with the SGP-path tripwires (#35) so one grep finds both.

### Coverage

| Period | DK | FD | Notes |
|---|---|---|---|
| FG main (spread + total + ML)        | ✅ | ✅ | All 30 teams |
| FG alternate spreads                 | ✅ | ✅ | When DK/FD posts them |
| FG alternate totals                  | ✅ | ✅ | When DK/FD posts them |
| F5 main (spread + total + ML)        | ✅ | ✅ | F5 ML may be ✗ at DK |
| F5 alternate spreads / totals        | ✅ | ✅ | |
| F7 main (spread + total)             | ✅ | ✅ | FD posts `First 7 Innings Run Line` + `Total Runs` — but only on the **default** event-page tab (see note below). |
| F7 alternate totals/spreads          | ✅ | ✗ | FD posts no F7 alts (vendor gap). |
| F3 main (spread + total)             | ✅ | ✅ | FD posts `First 3 Innings Run Line` + `Total Runs` (default tab). FD also posts `First 3 Innings Result` (3-way ML) which we skip. |
| F3 alternate totals/spreads          | ✅ | ✗ | FD posts no F3 alts (vendor gap). |
| F3/F5/F7 winner (for pick'em)        | ✅ | ✅ | DK: **bare period** name ("1st 3/5/7 Innings"), captured as `(period, "main")` ML rows. FD: 2-way "First 5 Innings Money Line" + 3-way "First N Innings Result" (F3/F5/F7) captured as `h2h_3way` with `tie_ml`. Both feed the dashboard pick'em DNB path (2026-05-27). |

### Pick'em (line-0 spread) → period winner (2026-05-27)

The dashboard treats a spread bet at **line 0** as a draw-no-bet and compares
it to each book's period **winner** (not its run line). That needs the 2-way
period ML captured here. DK's bare-period winner ("1st N Innings") is now
captured (row above) and flows via `get_dk_odds` → `h2h_1st_N_innings`. FD's
First-5 "Money Line" was already captured (matched by the `"money line"`
keyword in `classify_market`), so FD lights up on F5 pick'em cards but shows
a DNB-collapsible winner on F3/F7 too (below).

**FD 3-way "First N Innings Result" capture (2026-05-27).** FD posts a 3-way
"First N Innings Result" (Home/Tie/Away) for F3/F5/F7. `classify_market` now
maps it to `(period, "h2h_3way")` (single-inning and F2/4/6 Results stay
excluded), the parser writes the tie price to a new `tie_ml` column, and
`get_fd_odds` emits an `h2h_3way_1st_N_innings` row with `odds_tie`. The
dashboard's existing 3-way path (`scraper_to_canonical`'s `h2h_3way` shape +
`derive_pickem_american`) devigs it, drops the tie, and fills FD's F3/F5/F7
pick'em cells. So FD now shows on pick'em cards for all of F3/F5/F7 (F5 also
has its own 2-way Money Line).

### Run timing

Orchestrated by `Answer Keys/run.py mlb` in the parallel scrape phase
(pre-MLB.R), gated by `.scrapers_done_mlb`. The SGP scrapers continue to
run post-MLB.R (they depend on `mlb_parlay_lines`) on a separate trigger.

### Refactor — SGP scrapers now share clients

`scraper_draftkings_sgp.py` and `scraper_fanduel_sgp.py` were refactored
to import event discovery from `dk_client` / `fd_client`. SGP combo logic
and `calculateBets` / `implyBets` pricing stay in the SGP files.
Behavior-preserving — regression tests in `tests/test_sgp_regression.py`
compare current sgp_decimal output against captured golden baselines
within 0.20 decimal-odds tolerance.

### Team-name canonicalization

DK returns abbreviated city prefixes (`"CLE Guardians"`, `"LA Angels"`,
`"STL Cardinals"`). `scraper_draftkings_singles.py` maintains a
30-entry `DK_TEAM_MAP` that translates DK names → canonical (Odds API
format) before writing rows. FD names already match canonical — no
mapping needed.

## Additional SGP books (2026-06)

Three more books were added to widen the SGP consensus (a "top-down SGP
fair-odds" calculator). All follow the same shim/orchestrator/client trio as
Novig and write to `mlb_sgp_odds`. Both bots price them in-process via
`kalshi_common/sgp_service.py::SGPService` (registered in `DEFAULT_BOOKS` with
a `_run_book` dispatch branch, mirroring ProphetX/Novig); the legacy subprocess
path (`sgp_runner.py::SCRAPER_NAMES`) also lists them for rollback, and the
dashboard spawns the shims + blends them in `mlb_correlated_parlay.R`.

- **BetMGM** (`scraper_betmgm_sgp.py` + `betmgm.py` + `betmgm_client.py`) —
  ✅ GREEN, pure `curl_cffi`. Entain CDS API: harvest a static per-state
  `x-bwin-accessid` from `clientconfig` (unlock = the `x-bwin-sports-api: prod`
  header; use state `pa`), then `POST /cds-api/bettingoffer/picks` with legs
  sharing a `pickGroupId` returns the Angstrom correlated price. Source
  `betmgm_direct`. Verified live (4-corner overround ~1.12–1.18). No browser.
- **Caesars** (`scraper_caesars_sgp.py` + `caesars.py` + `caesars_client.py`
  + `caesars_waf.py` + `caesars_waf_node.js`) — token-broker + REST, **browser
  free**. `POST /sb/v2/bets/details` with `combinationSelections: []` returns
  the ZeroFlucs correlated price (`parlays[0].price.decimal`), logged-out.
  Every call needs an `aws-waf-token`, minted **without a browser** by running
  AWS WAF's real `challenge.js` under `node` (the `NetworkBandwidth` challenge —
  see `caesars_waf.py`; ~0.5s, cached ~4 min, validated before use, emits no
  rows on failure → never bad data). Requires `node` on PATH (no Playwright/
  Chromium). **Never pre-define `window.AwsWafIntegration` in the node
  harness** (issue #41): AWS's bundle installs itself only when that global is
  absent, so a placeholder makes `getToken()` return an empty string and the
  book silently drops to zero rows — that is exactly how Caesars died for a
  month. The SDK's auto-init rejection is logged and ignored instead, and every
  mint failure now logs at ERROR with node's stderr. `CaesarsClient.last_mint_sec`
  carries the last mint's wall-clock cost. Parsing is **exact-name** (`Run Line`/`Total Runs` + Alternate/F5
  variants) to exclude player props, team totals, and in-play (`... Live`)
  markets; events filtered to the **MLB** competition and **pregame** only;
  the away run-line leg carries its own (negated) line so all 4 corners price.
  Source `caesars_direct`. Verified 2026-07-30 (issue #41): mint 0.43–0.68s,
  tabs feed 200 / 178 KB JSON, sweep 8 rows in 2.2s with a 4-corner partition
  summing to 1.128 (12.8% hold), on-demand 2.62s. Target parallelism stays at
  3 — the WAF rate-limits aggressive hits.
- **bet365** — DEFERRED. Recon (`recon_bet365_*.py`) proved no fast path: odds
  live only on the `zap` WebSocket, which is Cloudflare-fingerprint-blocked for
  any non-browser client. Only a persistent live browser works; revisit if that
  tradeoff becomes acceptable.

## Fetch-health history — `sgp_fetch_health` (issue #38)

Every book fetch the bots make leaves one row behind, so "when did this book
last return a price?" is a SQL question instead of an archaeology project.
(The taker's feed died 2026-06-28 and went unnoticed for four weeks; nothing
recorded it.)

**Where:** the bot's own market DB — the same file, and therefore the same
write lock, as `mlb_sgp_odds`:

| bot | DB |
|---|---|
| maker | `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` |
| taker | `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb` |

The dashboard's CLI-shim path writes **nothing** — health recording is opt-in
via `SGPService(health_db_path=...)`, which only the bots pass.

**Schema** (`kalshi_common/sgp_health.py`):

| column | meaning |
|---|---|
| `fetched_at` | TIMESTAMPTZ, UTC |
| `book` | draftkings / fanduel / prophetx / novig / betmgm / caesars |
| `path` | `sweep` (periodic full-slate) or `on_demand` (live at-quote fetch) |
| `outcome` | `ok` / `empty` / `transport_error` / `timeout` / `error` |
| `rows_or_prices` | sweep: PricedRows written. on-demand: price calls that returned |
| `duration_sec` | the WHOLE logical fetch, including #34's retry sleeps |
| `error_class` | e.g. `BookTransportError:events:403`, `FutureTimeout` |
| `counters_json` | the #35 per-fetch counter snapshot |

`outcome` distinguishes four failure modes on purpose: `empty` is a book with
no markets, `transport_error` is a dead endpoint (#33), `timeout` is our own
deadline firing while the book was still working, and `error` is a bug in our
parse code. Collapsing them is the meta-bug this epic exists to remove. A book
skipped by `min_refresh_sec` writes **no row** — it was never fetched.

**Reading it** — `kalshi_common/fetch_health_queries.sql` ships four queries
(`outcome_mix_24h`, `current_failure_streak`, `parse_health_24h`,
`retention_footprint`):

```bash
duckdb kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb -readonly \
    < kalshi_common/fetch_health_queries.sql
```

⚠️ **Always group by `path`.** `targets_attempted` and `prices_returned` carry
different units on the two paths (sweep ~1:4, on-demand 1:1 — see the warning
block in `_shared.COUNTER_NAMES`). A ratio pooled across paths is wrong, not
merely noisy.

**Cost and safety.** Writes are buffered in memory and flushed once per bot
tick (`SGPService.flush_health()`); the on-demand path never touches the DB
inside a quote. A flush failure (locked DB) is swallowed and the buffer is
retained for the next attempt — a health write can never break pricing.
Retention is 30 days, pruned at flush time at most hourly. Measured footprint
(50k rows written through the real recorder, extrapolated): **~21 MB per 1M
rows**, and 1M rows is roughly a 30-day worst case at post-#53 quote volume
(~33k rows/day = 6 books x 1,440 sweeps + ~2k on-demand fetches x 6 books x 2).

## Run-time book-health alerting (issue #37)

The history above is the *memory*; this is the *alarm*. Both bots run
`kalshi_common/book_health.py::BookHealthAlerter` inside their tick loop, so a
book that dies mid-slate is loud within a few cycles instead of four weeks
later.

**Two rules**, both keyed on consecutive failed fetches — never on data age,
because an age rule would false-fire by design once #57 slows the sweep to
background structure-warming:

| Rule | Fires when | Message |
|---|---|---|
| A — book degraded | `BOOK_ALERT_STREAK` (default 3) consecutive failures on one (book, path) | `SGP book degraded: draftkings/sweep — 3 consecutive failed fetches (last: BookTransportError:events:403)` |
| B — consensus capacity | healthy books drop below the bot's own gate (maker `MIN_AGREEING_BOOKS`, taker `MIN_BOOK_COUNT_FOR_BLEND`) | `SGP consensus DARK: 1 healthy book(s), need 2 — quoting effectively blind` |

Each fires **once per incident** — a book failing 400 times in a row notifies
once — and re-arms only when the book answers again. Recovery is logged, not
notified. Channel: one `ERROR` line in `bot.log` plus a macOS notification
(`osascript`, same mechanism as `coverage_audit/`); it runs on a daemon thread
so a slow notification can never stall a loop that ticks 4x/sec.

⚠️ **`empty` is NOT a failure.** The predicate is `outcome NOT IN
('ok','empty')` — identical to `current_failure_streak` in the shipped SQL. A
book that answers with no markets is healthy: that is an off-day, a thin
slate, or (on-demand) "this book won't price that combo", which is the
*default* verdict for a plain `None`. Alerting on `!= 'ok'` would page on
every quiet morning and on every RFQ for a game a book doesn't list.

⚠️ **A skipped book is not a failing book.** `min_refresh_sec` skips write no
health row and make no observation; state simply persists. A gap means "we
didn't ask".

**Zero alerts while the bots are down.** State is in-memory and per-process,
so weeks of intentional downtime are silent by construction. It also means
streaks reset on restart — the first alert after a restart takes
`BOOK_ALERT_STREAK` fresh failures.

**Known gap — the silent parser.** A book whose parser breaks such that it
returns zero rows *without raising* records `empty` forever: Rule A stays
quiet and Rule B still counts it as healthy. That case is what #35's
`parse_failures` / `legs_resolved` counters are for (`parse_health_24h`), not
these alerts.

**Why in-memory counters and not a query over `sgp_fetch_health`:** the table
is a lossy view (flushes are swallowed and retried, so rows land up to 5s
late and later under a stuck lock), and on duckdb 1.4.4 an in-process
read-only connect while any read-write handle is open raises
`ConnectionException` — which the bots' narrow `IOException`/`CatalogException`
guards don't catch, so it would fall through to a fail-safe and **cost a
quote**. The monitor dashboard reads the table instead; it is a separate
process and hits the ordinary cross-process lock, which it already retries.

**Knobs** (both bots): `BOOK_ALERT_ENABLED` (default true),
`BOOK_ALERT_STREAK` (3), `BOOK_ALERT_PATHS` (`sweep,on_demand` — #57 sets
this to `on_demand` when every quote is live-priced; config, not code).

## DraftKings API Endpoints

| Endpoint | Auth | Purpose |
|----------|------|---------|
| `sportsbook-nash.../league/leagueSubcategory/v1/markets` | None | List MLB events |
| `sportsbook-nash.../event/eventSubcategory/v1/markets` | None | Main market IDs per game |
| `sportsbook-nash.../parlays/v1/sgp/events/{id}` | curl_cffi | **All selection IDs** (2MB response) |
| `gaming-us-nj.../en/api/wager/v1/calculateBets` | curl_cffi | **SGP pricing** (POST, returns trueOdds) — note the `/en/`, see below |
| `sportsbook-nash.../sgp/dkusnj/sportsdata/v2/sgp` | Full Akamai | SGP pricing (DK frontend only — **inaccessible** via REST) |

### Why curl_cffi?

DraftKings uses Akamai Bot Manager. Plain `requests` gets blocked by TLS fingerprinting. `curl_cffi` impersonates Chrome's TLS signature, which is enough to bypass Akamai on `calculateBets` and `parlays/v1/sgp/events`. The `sportsdata/v2/sgp` endpoint has stricter protection and is inaccessible — that's why we use `calculateBets` instead.

### The `/en/` prefix on calculateBets (issue #39)

DK reads on `sportsbook-nash` but **prices on `gaming-us-nj`**, and only the
price host is bot-protected. On **2026-06-24** an Akamai edge rule began
denying `POST /api/wager/v1/calculateBets` — DK's rows fell from ~18k/day to
zero while events and structure stayed perfectly green (200s, full selection
lists). That asymmetry is the signature of a price-host block, and it is why
`PriceCallTally` exists.

The rule matches that path **exactly**. DK's own betslip builds the URL as
`{locale}/api/wager/v1/calculateBets` (English maps to an empty prefix), and
the origin still routes the explicit `/en/` form, which the rule misses:

```
POST /api/wager/v1/calculateBets     -> 403 AkamaiGHost "Access Denied"
POST /en/api/wager/v1/calculateBets  -> 200 correlated SGP price
```

**The `/en/` is load-bearing — do not "tidy" it away.** It is pinned by
`tests/test_dk_price_host_403.py`, along with a guard that no module may
hardcode a second copy of the wager URL (the trifecta scraper used to, so a
fix here would have left that one 403ing).

This is *not* a TLS-fingerprint problem: every `impersonate=` target, DK's own
betslip headers, its anonymous enterprise JWT, and a real logged-out Chrome all
get the identical 403 on the bare path. If DK closes `/en/` too, the tell is
`error_class` `BookTransportError:price:403` in `sgp_fetch_health` — the next
move is re-reading `dkBetSlip.js` for the current route, **not** bumping
curl_cffi (which is shared by all 6 books).

**Things that DON'T work:** direct HTTP requests, page.evaluate(fetch()), cookie transfer from browser to requests, Playwright stealth plugins. All tested extensively.

## Selection ID Format

```
Spread: 0HC{market_num}{N|P}{line*100}_{suffix}
Total:  0OU{market_num}{O|U}{line*100}_{suffix}

Examples:
  0HC84191361N150_1   → Home team -1.5, market 84191361, suffix _1
  0HC84191361P150_3   → Away team +1.5, market 84191361, suffix _3
  0OU84203528O750_1   → Over 7.5, alt market 84203528, suffix _1
  0OU84203528U750_3   → Under 7.5, alt market 84203528, suffix _3
```

- **Suffix** (`_1` vs `_3`) varies per game — can't predict, must read from SGP parlays data
- **N** = negative spread (favorite), **P** = positive spread (underdog)
- **Market number** comes from market ID (e.g., `2_84191361` → `84191361`)
- Main market and alt market have **different** market numbers

## DK Market Structure

- **Main market** (subcategory 4519): Run Line (±1.5) + one Total (DK's main line) + Moneyline
- **Alt market**: Alt spreads (±1.0, ±2.5, etc.) + alt totals (every 0.5 from 5.0 to 13.0+)
- Both appear in the SGP parlays response
- **Game line markets** have BOTH spread AND total selections — inning/prop markets have only one. The scraper uses this to filter out non-game markets.

## Known DK Restrictions

### Cross-Market Blocking Near Main Line
DK won't combine the main run line with alt totals within ±0.5 of their main total. Example: if DK's main total is O/U 8.5, you can SGP with O7.5 or O9.5, but NOT O8.0 or O9.0. This is confirmed on DK's website too ("Sorry, your picks cannot be parlayed"). The scraper returns no price for these games rather than using a different total.

### Transient Rejections
Some games temporarily return `SelectionsCannotBeCombined` then work minutes later. The scraper retries once after a 2-second delay.

### SelectionClosed
Games near first pitch may close SGP pricing entirely.

## calculateBets Request/Response

### Request
```json
{
  "selections": [],
  "selectionsForYourBet": [
    {"id": "0HC84191361N150_1", "yourBetGroup": 0},
    {"id": "0OU84203528O750_1", "yourBetGroup": 0}
  ],
  "selectionsForCombinator": [],
  "selectionsForProgressiveParlay": [],
  "oddsStyle": "american"
}
```

### Response (success)
```json
{
  "selectionsForYourBet": [
    {"id": "0HC84191361N150_1", "trueOdds": 2.59, "displayOdds": "+159", "points": -1.5},
    {"id": "0OU84203528O750_1", "trueOdds": 1.87, "displayOdds": "−115", "points": 7.5}
  ],
  "bets": [
    {
      "type": "YourBet",
      "selectionsMapped": [{"id": "0HC84191361N150_1"}, {"id": "0OU84203528O750_1"}],
      "trueOdds": 4.0,
      "displayOdds": "+300"
    }
  ]
}
```

### Error responses
- **422 `SelectionsCannotBeCombined`** — DK won't combine these selections (cross-market restriction)
- **422 `SelectionClosed`** — Game near first pitch, SGP unavailable
- **200 with `combinabilityRestrictions`** — Selections recognized but can't be parlayed

## Output

Writes to `mlb_sgp_odds` table in `Answer Keys/mlb_mm.duckdb`:

| Column | Type | Description |
|--------|------|-------------|
| game_id | VARCHAR | Odds API event ID (joins to mlb_parlay_opportunities) |
| combo | VARCHAR | `"Home Spread + Over"` … or `"Home ML + Over"` … (moneyline family) |
| period | VARCHAR | "FG" |
| bookmaker | VARCHAR | "draftkings" |
| spread_line | DOUBLE | home-perspective spread; **NULL** for moneyline×total combos (no spread leg) |
| total_line | DOUBLE | the Over/Under line |
| sgp_decimal | DOUBLE | Decimal odds |
| sgp_american | INTEGER | American odds |
| fetch_time | TIMESTAMP | When scraped |
| source | VARCHAR | "draftkings_direct" |

> **NULL-safe dedup:** because `spread_line` is NULL for moneyline combos,
> `db.upsert_priced_rows` matches the composite key with `IS NOT DISTINCT FROM`
> (NULL-safe equality), not a plain tuple-IN. Never use a numeric sentinel for
> the absent leg — NULL keeps the grid/dashboard clean.

## SGP Scraping Playbook (for adding new books)

Lessons learned from building DK and FD scrapers. Use this when expanding to new sites (BetMGM, Caesars, ESPN BET, etc.).

### The 3 endpoints you need to find

Every book has these, just named differently:

| What | DK | FD | What to look for |
|---|---|---|---|
| **Event listing** | `league/.../v1/markets` | `scan/.../facet/.../search` | Returns today's games with event IDs + team names |
| **Market catalog** | `parlays/v1/sgp/events/{id}` | `sbapi/event-page?eventId=X&tab=same-game-parlay-` | Returns ALL selection IDs (main + alt) for one game |
| **SGP pricing** | `wager/v1/calculateBets` | `fixedodds/transactional/v1/implyBets` | POST with 2 selection IDs, returns correlated price |

### Recon process

1. Open the book in Chrome with DevTools → Network tab
2. Build a 2-leg SGP manually (spread + total)
3. Capture every network request during "add 2nd leg" — the SGP pricing call fires here
4. **Capture request HEADERS, not just URLs.** Missing headers was a multi-hour debugging session on FD. The API returned 200 with degraded data (singles only) instead of an error.
5. Try replaying the call with `curl_cffi` → if it works, you're done. If 400/403, check for missing headers or bot protection tokens.

### Common gotchas

**Silent degradation > explicit errors.** FD returns 200 with single-leg prices when headers are wrong, instead of 403. DK returns 422 with a clear error code. Assume the worst: always verify the response contains the SGP combined entry, not just a 200 status.

**Markets are split across tabs — no single tab is a superset (FD).** FD's `event-page?tab=` returns a *different* market slice per tab, and they only partially overlap. The `same-game-parlay-` tab returns ~156 markets (FG-alts + all F5) but **omits** the cumulative `First 7/3 Innings Run Line / Total Runs` and per-inning lines; the **default** tab (`tab=`) returns ~99 markets that **include** those F7/F3 lines but **omit** FG-alts and F5. Full coverage = the **union** of both tabs (`scraper_fanduel_singles.FD_TABS`), deduped by market/runner id. Lesson: don't assume one tab has everything — enumerate the tabs the event advertises and merge. (This is the bug that hid FD F7 totals: the singles scraper only read the SGP tab.)

**Selection ID formats vary wildly.** DK encodes market number + sign + line + suffix into strings like `0HC84191361N150_1`. FD uses plain integer `selectionId` tied to a `marketId`. Don't assume one format.

**Team disambiguation in alt spread markets.** Alt run line markets list BOTH teams at each line (e.g., "Reds -1.5" AND "Marlins -1.5"). You need to match runner team names against known home/away to avoid pricing the wrong team's spread. Key by `("home", line)` / `("away", line)`, not by handicap sign.

**Live games return adjusted handicaps.** A game in progress shows Run Line at +2.5/-2.5 instead of ±1.5. Filter events with `openDate < now` before scraping. Two events can exist for the same matchup (live + tomorrow's pre-game).

**Doubleheaders use hour-level matching.** Both scrapers dedupe by `(home, away, start_hour)` so two pre-game events for the same matchup (doubleheader) are preserved as separate events. The R staging table (`mlb_parlay_lines`) carries `commence_time` from the Odds API, and `match_events()` in each scraper compares the UTC hour of the DK/FD event against the stored `commence_time` hour to assign the correct `game_id`. Typical doubleheader slots (noon + 7 PM ET) are ~5 hours apart in UTC, so hour-level matching is sufficient.

**PerimeterX / bot protection tokens.** Some books require a `x-px-context` or similar token that's set by JavaScript on page load. `curl_cffi` impersonation alone isn't enough — you need the token too. These tokens are semi-persistent (days/weeks) but eventually rotate. Hardcode for v1, add auto-refresh later if needed.

### Vig measurement

For any book, compute vig from the 4 mutually-exclusive combos per game:

```
vig = sum(1/D for all 4 combos: HomeSpread+Over, HomeSpread+Under, AwaySpread+Over, AwaySpread+Under)
```

This sum exceeds 1.0 by the book's vig charge. Measured values:
- **DK:** stable ~1.125 (12.5% vig), consistent across time-to-game
- **FD:** bimodal — ~1.13 for games >21h out, ~1.22 for games <16h out

**Don't hardcode a single vig constant** when the book's vig varies. Use per-game devigging: divide each combo's implied prob by the per-game sum. Falls back to a constant when <4 combos are available.

**FG and F5 are separate partitions.** Never sum across both periods — that doubles the measured vig. Group by `(game_id, period)`.

### Line matching

The scraper first attempts exact-line matching: Wagerzon's exact spread and total lines must exist in the book's selection-ID dictionary, and that pair gets priced via `calculateBets` / `implyBets` / equivalent.

**Integer-line fallback.** When Wagerzon posts an SGP at an integer total line (e.g., FG Over 8, F5 Over 4) that the book doesn't quote directly, the scraper falls back to interpolating from the two adjacent half-point alts (X-0.5 and X+0.5). The integer-line derivation is implemented in `mlb_sgp/integer_line_derivation.py`. See "Integer-Line Derivation" below.

**No-fallback behavior.** If the book doesn't have either adjacent half-point alt, or any of the 4 bounds checks fails on the derived prices, the game's period is skipped for that book — same graceful degradation as today's "WZ line not found" path. The R blender treats the book as missing for that game.

### Integer-Line Derivation

When an integer total line `X` is missing, the scraper issues 8 SGP pricing calls (4 combos × 2 alts at `X − 0.5` and `X + 0.5`), per-alt-devigs each set, and derives joint fair probabilities at the integer line:

```
Δ_total      = (devig_HomeOver_lo − devig_HomeOver_hi)
             + (devig_AwayOver_lo − devig_AwayOver_hi)
             = P(T = X)   (the joint marginal push mass)

fair_prob_X  = devig_hi_combo / (1 − Δ_total)   for Over combos
             = devig_lo_combo / (1 − Δ_total)   for Under combos
```

Interpretation: the conditional joint probability "this combo wins outright, given the bet doesn't push." Books refund pushed legs regardless of spread side, so the conditioning matches real-world settlement. Pure-joint computation — no singles totals market accessed.

**Underlying principle:** Breeden & Litzenberger (1978) showed in options markets that adjacent strike prices implicitly carry the probability mass of the in-between outcome (the second derivative of call prices w.r.t. strike). This is the discrete sports analog: the difference between adjacent half-point alt SGPs gives us the implied joint integer-mass at `X`. Half-point pricing calculators (Sportsbook Review, Bookmakers Review, Unabated) apply the same skeleton in singles markets; we extend it to two-leg correlated parlays.

The derived rows are written with `source = '<book>_interpolated'` so post-hoc analysis can compare realized P&L on interpolated vs directly-priced rows.

**Bounds checks** (config constants in `integer_line_derivation.py`):

| # | Check | Threshold |
|---|---|---|
| 1 | Per-alt vig sum | `[1.05, 1.30]` |
| 2 | `Δ_total` plausibility | `[0.03, 0.18]` |
| 3 | Sum of 4 derived fair_probs | `[0.97, 1.03]` |
| 4 | Per-combo bounds | `(0, 1)` strict |

The push-mass cross-consistency check (initially shipped as bounds check #2) was removed after production data showed it tripping on real-world FD vig asymmetry (~10-15% disagreement is normal market noise, not data corruption). The two push-mass derivations are still computed and averaged for robustness.

On any failure: the game's period is skipped for that book, structured WARN logged with inputs and the violated check.

**Manual integration test:** `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m mlb_sgp.tests.test_integration_integer_line` — verifies a recent slate has interpolated rows with sum-to-one invariant holding. Exit 0=pass, 1=fail, 2=inconclusive (no integer-line games).

### Integration into the R scanner

New books write to `mlb_sgp_odds` with their own `bookmaker` and `source` values. The R scanner (`mlb_correlated_parlay.R`) reads all books, computes per-game vig for each, and blends their devigged fair probs equally with the model. Adding a third book requires:

1. Add `'newbook_direct'` to the `WHERE source IN (...)` query
2. Add a `NEWBOOK_SGP_VIG_DEFAULT` fallback constant
3. Add a `nb_vig_lookup` table (same pattern as DK/FD)
4. Add the `nb_row` / `nb_fair_prob` block in the per-combo loop (copy DK or FD block)
5. Optionally add `nb_fair_prob` to the output tibble

The blend automatically scales: `mean(model, dk, fd, nb)` when all present, falls back gracefully when books are missing.

## Observability: logging + per-fetch counters (issues #35, #36)

### Logging, not printing

Every module a bot imports uses `logger = logging.getLogger(__name__)`. The
library attaches **no handlers** — the host process decides where logs go (the
same principle as R's `message()` vs `cat()`). Both bots run these scrapers
in-process via `kalshi_common.sgp_service`, so a `print()` here would bypass
their handlers and levels entirely and never reach `bot.log`.

CLI entry points (`scraper_*_sgp.py`, the two `*_singles.py` scrapers) call
`logging.basicConfig(..., stream=sys.stdout)` inside their `__main__` block
only, so subprocess/dashboard runs still write to `mlb_sgp/logs/<scraper>.log`
exactly as before.

Level policy — chosen so one broken book cannot drown `bot.log`:

| event | level |
|---|---|
| per-game "N targets → M offered" | DEBUG |
| individual `SANITY-DROP` | DEBUG (the per-fetch count rides in the summary) |
| broad-`except` catch | WARNING the **first** time per stage per fetch, DEBUG after |
| per-fetch summary | INFO (sweep) / DEBUG (on-demand — it runs per RFQ per book) |
| parse tripwire | WARNING |
| auth / token-mint failure | WARNING–ERROR |

`mlb_sgp/tests/test_library_logging.py` enforces this with an AST scan: no
`print()` outside a `__main__` guard, a module logger in every bot-path file,
and no handler configuration at import time. Standalone diagnostics
(`recon_*.py`, `probe_concurrency.py`, `quick_recon.py`, the Pikkit modules)
are deliberately out of scope — no bot imports them.

### Per-fetch counters and the parse tripwire

`_shared.FetchCounters` records what one book's one fetch actually did. It is
the parse-side counterpart to `PriceCallTally` (#33): the tally answers *"is
the book's price endpoint reachable?"*, the counters answer *"does our parser
still understand what the book sends back?"*

| counter | sweep | on-demand |
|---|---|---|
| `events_seen` / `events_matched` | events listed / matched to our slate | 1 / 1 when the game is listed |
| `legs_attempted` / `legs_resolved` | targets checked against parsed structure / offered by the book | `len(legs)` / `len(resolved)` |
| `targets_attempted` / `prices_returned` | targets priced / priced rows produced | price calls / non-`None` decimals |
| `sanity_drops`, `parse_failures`, `transport_errors` | same on both paths | |

Counters sit **outside** `request_with_retry` (#34): one logical fetch is one
tick no matter how many attempts the wire call took. The per-stage breakdown
rides along as `stages=price_combo:3,resolve_legs:1`.

**Transport vs parse.** A broad `except` that could have caught something off
the wire calls `counters.record_exception(...)`, which files a
`BookTransportError` under `transport_errors` and everything else under
`parse_failures`. This matters because CZR/MGM swallow a price-stage transport
error inside the client and return `None`, while DK/FD/NV/PX price through
module-level helpers whose `request_with_retry(stage="price")` **raises** — so
without the split, a 403'd FanDuel price endpoint reported `parse_failures=4`
and fired a PARSE tripwire, sending a fixer to the parser instead of the
endpoint.

**Units differ by path — do not compare across them.** On the sweep,
`targets_attempted` counts target LINES (each fanning out to ~4 combo calls)
and `prices_returned` counts PricedRows; on-demand both count PRICE CALLS.
Each path's tripwire is right in its own units, but a cross-path ratio is not
meaningful. `_shared.COUNTER_NAMES` carries the same warning for #38.

Both `price_sgps` (sweep) and `SGPService.price_on_demand` scope a fetch with
the `fetch_counters(...)` context manager, so the summary line is emitted on
every exit path — including the several `return []` short-circuits and the way
out of a `BookTransportError`. Every fetch logs one line:

```
sgp fetch book=fanduel path=sweep events_seen=15 events_matched=15 \
  targets_attempted=42 legs_attempted=60 legs_resolved=42 prices_returned=168 \
  sanity_drops=0 parse_failures=0 transport_errors=0
```

**Tripwires** (WARNING, `PARSE TRIPWIRE book=… path=… tripped=…`):

| name | condition | what it means |
|---|---|---|
| `events_unmatched` | `events_seen > 0 and events_matched == 0` | book is up and listing games, none map to our slate (team-name / canonical-match drift) |
| `legs_unresolved` | `legs_attempted > 0 and legs_resolved == 0` | **the FanDuel alt-total case** — the book offered lines, our parser recognized none |
| `prices_empty` | `targets_attempted > 0 and prices_returned == 0 and transport_errors == 0` | legs resolved but nothing priced, and transport never complained |

Each compares a stage's OUTPUT against its INPUT, so a genuinely empty slate
(no input) never fires, and a transport failure is never dressed up as a parse
bug — that is #33's story and is already loud there.

On the on-demand path the counters also ride back on
`OnDemandBookResult.counters` (an immutable `FetchCountersSnapshot`), which is
what #38's per-book fetch-health history will persist.

`counters` is **keyword-only** everywhere it was threaded into a book module,
so a future positional argument can never silently bind to it (BetMGM's
`price_selection_set(client, refs, fixture_id)` is the live hazard).
`mlb_sgp/tests/test_counter_signatures.py` proves this with
`inspect.signature().bind()`.

## Concurrency & Logs

The four SGP scrapers (DK, FD, ProphetX, Novig) are launched in parallel by
`Answer Keys/mlb_correlated_parlay.R` via `parallel::mclapply` (4 forked R
workers, one per scraper). Wall-clock time of the SGP refresh block is the
*max* of the four scrape times rather than the sum.

**Per-scraper logs:** stdout + stderr from each scraper are captured to
`mlb_sgp/logs/<scraper>.log` (overwritten each run, gitignored). The
orchestrator also prints a one-line summary per scraper (elapsed seconds,
exit code, log path) plus an overall wall-clock line — check the R console
output to see which book is the slowest on a given run.

**DuckDB write contention:** All four scrapers write to `Answer Keys/mlb_mm.duckdb`,
which DuckDB only allows one writer to open at a time. `db.py` wraps every
write `connect()` call in `_connect_with_retry()` (exponential backoff +
jitter, up to 10 attempts) so transient lock collisions between scrapers are
invisible to callers. Read connections are not retried (DuckDB allows
unlimited concurrent readers). Note: the trifecta scraper (`scraper_draftkings_trifecta.py`) is separate — it writes `mlb_trifecta_sgp_odds` to `Answer Keys/mlb.duckdb`.

### Target-level parallelism (2026-06)

Each orchestrator fans target lines out on a thread pool; the 4-combo
pool nests inside (total in-flight requests = parallelism × 4):

| Book | Default | Env override | Probed ceiling (2026-06-16) |
|---|---|---|---|
| DraftKings | 8 | `MLB_SGP_DK_PARALLELISM` | no backoff to 32; throughput plateaus ~16 (sel-id fetch floor), so >16 is pointless |
| FanDuel | 4 | `MLB_SGP_FD_PARALLELISM` | not probed (not the long pole) |
| ProphetX | 6 | `MLB_SGP_PX_PARALLELISM` | no backoff (429/403) to 12; default 6 for ≤60s margin + modest RFQ burst |
| Novig | 4 | n/a (module constant) | not probed |

Live ramp (2026-06-16) found **zero** rate-limit / bot-detection signals at
any tested width — DK to 32, PX to 12. DK's constant "422 ×40" is its
legitimate cross-market combinability rule, not backoff. The practical
limits are diminishing returns (DK's fixed sel-id-fetch floor) and
RFQ-burst prudence (PX), not the books rejecting us.

`price_sgps(..., parallelism=N)` overrides both env and default. Ceilings
come from `probe_concurrency.py` — a **manual-only** ramp harness
(budgeted: ≤150 DK calls, ≤40 PX RFQs, cooldowns, post-run health
check). Run it in the morning (~7–8am PT, no games live):

    mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book dk

**The Kalshi bots no longer spawn these scrapers as subprocesses** — they
price in-process via `kalshi_common/sgp_service.py::SGPService`
(persistent sessions + TTL-cached structure fetches). The dashboard still
uses the CLI shims (`scraper_*_sgp.py`), which behave exactly as before.

---

## FanDuel Scraper

`scraper_fanduel_sgp.py` mirrors the DK scraper but is much simpler since FD's selection IDs are plain integers tied to marketIds (no DK-style `0HC...` decoding) and FD doesn't lock its pricing endpoint behind Akamai.

### How it works

```
FD scan API           → MLB events list (competitionId 11196870)
canonical_match       → resolve FD team names to internal game_ids
FD event-page API     → all SGP-eligible runners (main + alt, FG + F5)
                        via &tab=same-game-parlay- (156+ markets)
Match Wagerzon lines  → exact spread + total lookup
FD implyBets API      → POST each combo, parse isSGM=true entry
                      → mlb_sgp_odds (bookmaker='fanduel')
```

### Required headers

FD's API silently strips the SGP combination from `implyBets` responses (returning only single-leg prices) if any of these are missing:

| Header | Value | Notes |
|---|---|---|
| `x-application` | `FhMFpcPWXMeyZxOx` | FD's API key, also passed as `?_ak=` query param |
| `x-sportsbook-region` | `NJ` | Geo header. The `nj.` hostnames are FD's backend routing — works from any state, no VPN required. |
| `x-px-context` | `_pxvid=...;pxcts=...;` | PerimeterX visitor token. Hardcoded in the scraper; long-lived but rotates eventually. |

### Combo logic

Per game, 4 combos × 2 periods (FG + F5) = up to 8 prices. Matches Wagerzon's exact lines:

- Home Spread + Over / Under
- Away Spread + Over / Under
- F5 Home Spread + Over / Under
- F5 Away Spread + Over / Under

**Exact line matching only.** The scraper fetches main + alt spreads and totals from FD's SGP tab, then looks up the exact Wagerzon spread and total. If FD doesn't have the precise line, that game is skipped (no approximate matching). Alt spreads (e.g., ±2.5, ±3.5) are resolved by matching runner team names against the event's home/away teams.

### When the PerimeterX token expires

If FD starts returning 400s with empty bodies, or `implyBets` stops returning the `isSGM=true` entry, the `x-px-context` token has rotated. To refresh:

1. Open `sportsbook.fanduel.com/navigation/mlb` in Chrome with DevTools open
2. Find any `event-page` request in the Network tab
3. Copy the `x-px-context` request header
4. Paste into `FD_PX_CONTEXT` constant in `scraper_fanduel_sgp.py`

A future v2 may bootstrap this automatically via headless Chrome.

### FD-specific lessons learned

- **Event name format:** `"Away Team (P Pitcher) @ Home Team (P Pitcher)"` — strip the pitcher parens before team matching.
- **`nj.` hostnames are NOT geo-restricted.** `sib.nj.sportsbook.fanduel.com` works from California with no VPN. The `nj` prefix is FD's backend routing, not a geo gate.
- **`implyBets` returns 3 entries for a 2-leg combo:** 2 `SINGLE` (one per leg, `isSGM=false`) and 1 `DOUBLE` (`isSGM=true` — this is the SGP price). Parse `winAvgOdds.trueOdds.decimalOdds.decimalOdds` for the decimal value.
- **FD's SGP vig is bimodal by time-to-game:** ~13% for games >21h out, ~21% for games <16h out. The step change happens around 16-21h before first pitch — possibly tied to lineup posting windows.
- **Alt total runner names use parens:** `"Over (8.5)"`, `"Under (7.5)"`. Main total runners are just `"Over"` / `"Under"` with the line in the `handicap` field.
- **Alt spread runner names embed the team:** `"Cincinnati Reds +3.5"`. The `handicap` field is 0 for alts. Parse team name + signed line from the string, match team to home/away.
- **F5 totals at integer values (e.g., 5.0) trigger interpolation.** FD's F5 alt totals jump in 1.0 increments (2.5, 3.5, 4.5, 5.5...). When Wagerzon has F5 total 5.0, FD's exact-line lookup misses → the integer-line fallback kicks in, using FD's F5 alts at 4.5 and 5.5. See "Integer-Line Derivation" above.
- **2026-05-13:** Fixed alt-spread bucket key in `scraper_fanduel_singles.py` — was `abs(effective_line)` which collapsed opposite-direction same-magnitude lines (e.g. KC -2.5 and KC +2.5) into one row. Now buckets by signed home-team line so both directions persist. Verified event 35600618 went 7 → 14 alt-spread rows. See `docs/superpowers/research/2026-05-13-fd-recon-findings.md`.

## Library architecture (line-source pivot, 2026-05-13)

The SGP scrapers were refactored into a 3-layer stack so the same pricing code
can serve both the MLB dashboard (writes to `mlb_mm.duckdb`) and the Kalshi
MLB RFQ bot (writes to a sibling `kalshi_mlb_rfq_market.duckdb`).

- **Per-book HTTP clients** — pure transport, no domain logic:
  - `dk_client.py` — DK leagues / event-markets / parlays / `calculateBets`
  - `fd_client.py` — FD scan / event-page / `implyBets`
  - `prophetx_client.py` — ProphetX RFQ endpoint
  - `novig_client.py` — Novig anonymous `/unauthenticated` SGP endpoint
- **Per-book SGP orchestrators** — `price_sgps(targets) -> List[PricedRow]`:
  - `draftkings.py`, `fanduel.py`, `prophetx.py`, `novig.py`
  - Each loads its client, walks the target `(game_id, period, spread_line, total_line)` tuples, prices all 4 combos per tuple, devigs, returns `PricedRow`s.
- **Thin scraper shims** — `scraper_{book}_sgp.py`:
  - Load target lines via `_shared.load_target_lines()`
  - Call the book's `price_sgps()`
  - Upsert via `_shared.upsert_priced_rows()`
  - No pricing logic — just I/O glue.

**Target-line tables.** Two sources, same shape:
- `mlb_parlay_lines` in `mlb_mm.duckdb` — Wagerzon-derived; one target line per game/period; used by the dashboard pipeline.
- `mlb_target_lines` in the bot's sibling market DB — Kalshi MVE-derived; many target lines per game (every `(spread, total)` tuple Kalshi lists); used by the bot's SGP cadence loop.

`load_target_lines()` reads whichever table exists in the connected DB.

**Env overrides** (consumed by every scraper shim and `_shared.py`):
- `MLB_SGP_DB_PATH` — full path to the DuckDB the scraper should read targets from and write priced rows back to. Defaults to `Answer Keys/mlb_mm.duckdb`. The bot sets it to its sibling market DB so dashboard data is never touched.
- `MLB_SGP_PERIODS` — comma-separated list of periods to price (`FG`, `F5`, `F7`). Defaults to all. The bot sets `FG` only.

**Output schema.** `mlb_sgp_odds` now carries `spread_line` and `total_line` columns alongside `combo`, so multi-line target lines from Kalshi can be priced and dedupe-merged without collision (the dashboard pipeline writes a single line per game/period, the bot writes many).

## Failure semantics — "book down" vs "no markets" (issue #33, 2026-07-25)

**The rule: an empty result means the book has nothing to offer. A dead book
raises.** Before this, every client swallowed a 403 / DNS failure / WAF
challenge into `[]`, which `SGPService._book_done` scored as SUCCESS (resetting
the strike counter so the client reinit never fired) and `sgp_cycle` answered by
running `clear_source` — *deleting* the book's rows. DK (403), Novig (DNS) and
Caesars (WAF) each rotted for a month with zero alerts.

Two of those three diagnoses turned out to be wrong once the wire was actually
read: DK was an Akamai path-exact rule, not a stale fingerprint (issue #39),
and Novig's DNS was a transient blip — the host resolves and the book prices
(issue #40). Treat a recorded cause as a hypothesis until reproduced.

### Novig: 400 vs 403 at the price stage (issue #40)

Novig's anonymous parlay endpoint has two non-200s that mean opposite things:

| status | meaning | how often |
|---|---|---|
| `400 {"message": "Cannot price parlay"}` | the book declines **this combo** | common and normal |
| `403 <html>` | we have been **rate-limited** | after burst traffic |

Both used to collapse into a bare `BookTransportError:price`. They now carry
the status (`price:400` vs `price:403` in `sgp_fetch_health.error_class`), on
both the sweep and the on-demand path.

**Novig rate-limits bursts.** Two full sweeps back-to-back drove 100/100 price
calls to 403; the same combos recovered after a cooldown. The production
sweep's ~7-minute cadence is what has been hiding this. Live-at-RFQ pricing
(epic #53) is exactly the burst pattern that will trigger it — watch
`price:403`.

#### Novig leg resolution is line-keyed on BOTH paths (issue #64)

The sweep and the on-demand path now resolve legs the same way: one GraphQL
EventMarkets call per game, parsed by `build_line_structure` into per-**line**
buckets, then `_lookup_leg`'d per target at exactly the lines that target
asked for. `fetch_event_legs` is still the fetch, but only its raw `markets`
half is read — its single-line leg dict is no longer used by the sweep.

This closes a defect found while verifying #40 and fixed in #64. The sweep
used to collapse every TargetLine of a game into ONE `fg_spread_line` /
`fg_total_line` (the grouping loop overwrote, so the last target won), resolve
the leg UUIDs once per *game*, then price all ~33 of that game's targets with
those SAME legs while stamping each row with its own `spread_line` /
`total_line` — N rows carrying one price under N labels, which fed the maker's
cross-book median a wrong number presented as a real one. It was reproduced
offline (3 targets at totals 7.5 / 8.5 / 12.5 requested a single leg set and
produced 12 rows) and live (totals 6.5 / 8.5 / 9.5 all priced
17.5439 / 4.1667 / 4.8544 / 1.4925, hold 18.6%). The on-demand path was never
affected, which is more evidence for the epic's de-sweep rule.

Two consequences worth knowing:

* A target whose exact `(spread, total)` Novig does not offer on both sides
  now yields **no row**, where it previously got a plausible-looking wrong
  one. Integer totals still route to `try_integer_fallback_nv` (which was
  always per-target correct) and stay labelled `novig_interpolated`.
* **ML × total** is priced once per game — pricing it at every offered total
  would multiply wire calls into Novig's rate limit. Its over/under legs are
  looked up at `fg_total_line`, so the label names the line actually priced,
  but which rung of the ladder that is remains arbitrary.

### The contract

| Situation | Client returns / raises | Caller must |
|---|---|---|
| Non-200 (not 404), connection error, timeout, malformed JSON, failed auth gate | raise `BookTransportError` | **preserve** the book's existing rows + count a failure |
| Valid payload listing no markets | `[]` / empty structure | clear + rewrite the source (unchanged) |
| `404` on ONE event's structure fetch | `[]` / `None` for that event | skip that game, keep the cycle |

`BookTransportError(book, stage, status_code=…, cause=…, detail=…)` lives in
`_shared.py`. `stage` is `"auth"` (MGM accessid, CZR AWS-WAF token),
`"events"`, `"structure"`, or `"price"`. Helpers `check_response()` and
`json_or_raise()` do the raising so the six clients don't hand-roll it.

The 404 carve-out matters: without it, one postponed or delisted game would
abort a whole book's cycle.

### Orchestrators never catch it

`price_sgps` lets `BookTransportError` propagate. `SGPService._run_book_safe`
converts it to `None` and logs at ERROR with book + stage + status; `None` is
the fail-safe signal that makes `sgp_cycle` skip `clear_source`. The on-demand
path (`price_on_demand`) feeds the *same* `_book_done` strike path, so sweep and
on-demand share one recovery. A clean "this book won't price this combo" (event
miss, `None` structure, devig rejection) still counts **no** strike.

### `PriceCallTally` — the all-price-calls-failed verdict

A book can read fine and price nothing: DK READS from `sportsbook-nash` but
PRICES on `gaming-us-nj`, so the price host can be blocked on its own. Each
`price_sgps` snapshots the book's tally, records every price call, and calls
`verdict()` before returning — raising `stage="price"` when a cycle made **≥ 8**
price attempts (`PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT`, ~2 games of combos)
with **zero** successes. Below that threshold a legitimately thin one-game slate
would look like a dead book, so the tally stays silent.

Individual declined combos remain row-drops — they are partial, not systemic.

### Consequence for consumers

A book judged down keeps its **previous** rows rather than losing them. That is
safe because the maker's staleness gate refuses stale fairs (it does not trade
them) — absence of fresh rows becomes a no-quote, not a bad quote. If that gate
is ever removed, revisit this fail-safe.

### Shims

All six `scraper_*_sgp.py` shims now catch a failed `price_sgps`, print
`preserving last cycle's rows`, and exit 1 instead of clearing the source. The
ProphetX shim additionally **defers** `clear_source` until the first batch
prices without a transport error (it used to clear up-front, so a dead PX wiped
the source before we knew it was dead). An empty-but-healthy slate still clears:
the flag flips on a successful *call*, not on a non-zero row count.

## Retry & backoff (issue #34, 2026-07-27)

Making failures loud (#33) measurably **widened their blast radius**: a 5xx on
one event's structure fetch used to return `[]` and skip that one game; after
#33 it raises and the book loses the whole cycle. Retry shrinks that window
back down — a transient blip must not cost a book its fetch, nor (post-#54)
its vote on a live quote.

`_shared.request_with_retry(call, profile=…, book=…, stage=…)` wraps every
wire call in the six clients.

### The two profiles

| Profile | Attempts | Backoff | Max added latency | Used by |
|---|---|---|---|---|
| `RETRY_BACKGROUND` | 3 | 1.5s → 3.0s, ±25% jitter | **≤ 5.63s** (hard cap 10s) | sweep / structure-warming fetches (auth, events, structure) |
| `RETRY_LIVE` | 2 | 0.5s, ±25% jitter | **≤ 0.63s** (hard cap 1.0s) | on-demand pricing fetches, **and every price call on both paths** |

Those ceilings are per **wire call**, not per pricing operation, and are pinned
by `test_retry_backoff.py::test_worst_case_added_delay_stays_within_profile`.
`max_total_delay_sec` is a hard cap on cumulative sleep, so the live path can
never trade a lost book for a blown RFQ quote budget.

**Why price calls are always LIVE:** they fan out across a thread pool
(targets × combos at DK), so a 3-attempt backoff there would stall a degraded
cycle for minutes. One fast retry is the right trade at that fan-out.

### What is and isn't retried

| Outcome | Retried? |
|---|---|
| 429 (honors `Retry-After` when it FITS the profile's budget, else backs off) | ✅ |
| 5xx | ✅ |
| Connection reset / DNS / timeout / TLS | ✅ — classified via `OSError`, which both `curl_cffi`'s and `requests`' exception bases subclass |
| 200 with an unusable body, where a `body_ok` predicate is supplied | ✅ — Caesars' AWS-WAF challenge is HTTP 200 + HTML |
| **404** | ❌ — #33's per-event carve-out; one delisted game costs exactly one attempt |
| **4xx auth/permission** (401/403) | ❌ — a real verdict; DK's production 403 fails immediately |
| Programming errors (`TypeError` etc.) | ❌ — still raise `BookTransportError` on attempt 1 |

The helper **returns the last response** rather than judging it, so status
classification stays in `check_response()` / `json_or_raise()` and
`BookTransportError` keeps being constructed in exactly one place. It raises
only when every attempt threw, and stamps the attempt count into `detail`
(`failed after 3 attempt(s) [background profile]`) so a `bot.log` reader can
tell a transient blip from a book that is simply gone.

### Choosing the profile at a call site

Every structure/events/auth function takes `profile: RetryProfile =
RETRY_BACKGROUND`, so sweep call sites are unchanged. The **on-demand path
binds `RETRY_LIVE` at the `SGPService` seam**: `_dk_fetchers(st, profile)` /
`_fd_fetchers(st, profile)` bake the profile into the fetcher lambdas, and the
PX/NV/MGM/CZR hooks pass it at the client call. That keeps the profile out of
the six orchestrators' signatures entirely. Both bindings share one `TTLCache`,
so a warm entry written by either path serves the other.

This is an explicit argument, deliberately **not** a `contextvars` ambient:
`ThreadPoolExecutor` does not propagate context, so an ambient profile would
silently fall back to BACKGROUND inside every book's worker pool — precisely
the latency blowup the LIVE profile exists to prevent.

### Caesars behavior change

`caesars_client._get_json` used to hand-roll a 4×1.5s loop that retried
**everything**, including a 403 auth gate — burning ~6s before admitting the
book was dead. It now uses the shared helper: 4xx fails on the first attempt,
while the 200-with-HTML WAF challenge retry is preserved via `body_ok`.

### Testing note

`mlb_sgp/tests/conftest.py` installs an autouse fixture that **records** backoff
instead of sleeping (and pins jitter to its midpoint), so the suite never burns
wall-clock on retries and the schedule is directly assertable. A test that cares
takes `recorded_sleeps` as an argument and reads the durations that would have
been slept.

## Six-book verification + the pre-restart checklist (issue #42)

`mlb_sgp/verify_books.py` is the bot-restart gate. It answers, per book, on the
path a quote actually uses (`SGPService.price_on_demand`): **can this book
price a real combo right now, and how fast?**

```bash
python3 mlb_sgp/verify_books.py                        # all six books
python3 mlb_sgp/verify_books.py --json /tmp/gate.json  # machine-readable
python3 mlb_sgp/verify_books.py --books novig --pacing 12   # rate-limit safe
```

It is read-only: no DB is opened for write, no bot is touched.

### Why it separates four outcomes

A book returning no price is ambiguous, and collapsing that ambiguity is how
DK, Novig and Caesars each rotted for a month. The harness never collapses it:

| outcome | meaning | action |
|---|---|---|
| `priced` | resolved legs, returned a fair | none |
| `not_offered` | structure healthy, carries none of the probed lines | none — a coverage fact |
| `no_price` | legs resolved, book declined to price them | check the shape (see below) |
| `down` | event match or structure fetch failed | **a real failure** |

Coverage is probed through each book's own `resolve` hook — a pure lookup
against an already-fetched structure — so it costs zero extra wire calls.

### Line coverage differs wildly per book

Measured 2026-08-03, one game (WSH @ PHI, main total 9.0), out of 4 candidate
spreads and 7 candidate totals:

| book | spreads | totals | shape |
|---|---|---|---|
| DraftKings | 4 | 7 | full ladder, half AND whole numbers |
| ProphetX | 4 | 5 | wide ladder |
| Novig | 4 | 4 | half-point ladder only (3.5–13.5) |
| BetMGM | 4 | 4 | half-point ladder only |
| FanDuel | 4 | **1** | main total only |
| Caesars | **1** | **1** | main line only |

**No single rung is priced by all six books.** The books split on whole-number
vs half-point totals: at a 9.0 main line, FD and Caesars carry it and Novig and
BetMGM do not; at 8.5 the reverse. Expect a quorum of ~4 of 6 at any rung, not
6 of 6 — this is the number #55 must design against, and it is a property of
the books, not a bug.

### Shape capability matrix

| shape | DK | FD | PX | NV | MGM | CZR |
|---|---|---|---|---|---|---|
| spread × total | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| ML × total | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| spread × ML | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 3-leg (spread+total+ML) | ❌ | ❌ | ⚠️ | ⚠️ | ❌ | ⚠️ |

**No book prices spread × ML.** Verified on DK for all three side
combinations (home+home, home+away, away+away) — it is not a same-team
restriction, the pairing is simply not offered. Since an MLB same-game 3-leg
can only draw on spread / total / ML, every 3-leg necessarily contains that
pair, which is why the 3-leg column is mostly ❌.

⚠️ **The three books that do return a 3-leg number return a WRONG one — do not
quote 3-leg shapes.** A home −1.5 spread *implies* the home ML, so
`P(spread ∧ total ∧ ML)` must equal `P(spread ∧ total)`. It does not:

```
Novig  [home -1.5, over 7.5]             dec=3.1949
Novig  [over 7.5, home ML]               dec=2.3981
Novig  [home -1.5, over 7.5, home ML]    dec=2.3981   <-- spread leg DROPPED
```

Novig silently drops the spread leg and returns the 2-leg price under a 3-leg
label. Its partition confirms it: cells differing only in the spread's side are
byte-identical, and the implied probabilities sum to **2.364** instead of
~1.18. ProphetX is less dramatic but still wrong — the book itself prices the
3-leg identically to the 2-leg (3.25 both), yet our engine devigs them to
0.2975 vs 0.2885, because Route A cannot complete (the logically impossible
cells — home covers −1.5 *and* away wins — are correctly refused by the book)
and Route B's transfer does not preserve the implication.

### Cross-book agreement

#39 closed with "cross-book price equivalence unverified at n=1". Closed here:
the harness re-prices one common rung at every book that carries it. Measured
2026-08-03, `home -1.5 + over 7.5`, devigged fair:

| book | fair |
|---|---|
| BetMGM | 0.2550 |
| Novig | 0.2652 |
| DraftKings | 0.2661 |
| ProphetX | 0.2885 |

Range 0.2550–0.2885, **spread 13.1% of the low, median 0.2656**. ProphetX runs
consistently high — it was also the high book at 8.5 (0.2595 vs ~0.23). Four
books is the realistic quorum at a half-point rung; FanDuel and Caesars carry
only the main line and so cannot vote on it.

This is the strongest correctness signal available without a ground truth: an
independently-wrong book shows up here, whereas #20's dispersion gate cannot
see books that are jointly wrong.

### Latency baseline (feeds #50 and #55)

`price_on_demand`, seconds, measured 2026-08-03. **Cold** is a fresh
`SGPService` — event discovery + structure fetch + price, i.e. what the maker
pays on a game's first RFQ. Warm reuses the cached structure.

| book | shape | cold | warm p50 | warm p95 | structure fetch |
|---|---|---|---|---|---|
| DraftKings | spread × total | 2.77 | 0.76 | 3.55 | 1.37 |
| DraftKings | ML × total | 2.63 | 0.77 | 0.77 | |
| FanDuel | spread × total | 1.12 | 0.85 | 1.05 | 0.64 |
| FanDuel | ML × total | 0.90 | 0.72 | 0.74 | |
| ProphetX | spread × total | 2.42 | 1.38 | 1.43 | 0.51 |
| ProphetX | ML × total | 1.80 | 1.45 | 2.36 | |
| ProphetX | 3-leg | 2.78 | 1.76 | 1.86 | |
| Novig ¹ | spread × total | 6.69 | 3.89 | 9.06 | 0.33 |
| Novig ¹ | ML × total | 4.63 | 3.68 | 4.49 | |
| BetMGM | spread × total | 1.71 | 1.30 | 1.35 | 0.85 |
| BetMGM | ML × total | 1.70 | 1.27 | 1.35 | |
| Caesars ² | spread × total | 2.68 | 2.20 | 2.36 | 1.36 |
| Caesars ² | ML × total | 2.54 | 2.36 | 2.54 | |
| Caesars ² | 3-leg | 7.04 | 6.48 | 7.26 | |
| Caesars ² | AWS-WAF token mint | 0.67 | (cached 240s) | | |

¹ **Novig must be measured with `--pacing 12`.** At the default 1s pacing it
returned zero warm samples on three of four shapes — that was its own rate
limit (26 rapid price calls → 403, #40), not latency. Every shape here fans
out into a 2^N partition, so a full verification run is 60+ wire price calls
per book. Novig is the slowest book by a wide margin and its p95 sits at the
edge of #50's 8–10s budget.

² Caesars numbers are from the run before its WAF began throttling this IP —
see the checklist note below.

### Caesars WAF throttles under verification load

#41 closed with "not verified: behaviour under WAF rate-limiting". Now
observed: after a full verification run plus several diagnostics, Caesars began
returning `WAF token minted but NOT validated (attempt n/3)` and stayed that
way for **over 30 minutes**. The mint itself succeeds — the WAF rejects the
minted token, which is the signature of IP throttling rather than a broken
integration (the same code minted in 0.67s an hour earlier).

Operationally this is handled correctly: the client raises
`BookTransportError(stage="auth")`, and #33's contract preserves Caesars' prior
rows instead of clearing them. But it means **Caesars should be verified first
in a session, not last**, and a verification run should not be repeated
back-to-back.

### Pre-restart checklist

Run this before restarting the maker or taker. Self-contained; nothing here
starts, stops or restarts a bot.

1. **Confirm you are on a healthy slate.** Games must be upcoming, not live:
   `python3 mlb_sgp/verify_books.py` prints the game it picked and its start.
2. **Verify Caesars FIRST if you will run anything else** —
   `python3 mlb_sgp/verify_books.py --books caesars`. Its WAF throttles under
   load and takes >30 min to clear.
3. **Run the gate:** `python3 mlb_sgp/verify_books.py --json /tmp/gate.json`.
   Every book must report `priced` or `not_offered`. Any `down` blocks the
   restart — read `coverage.error` in the JSON for the stage.
4. **Re-measure Novig separately** with `--books novig --pacing 12`; at default
   pacing its warm samples are rate-limit noise.
5. **Check cross-book agreement.** The harness re-prices one common rung at
   every book carrying it. A spread wider than ~15% of the low means one book
   is mispricing — do not restart into that.
6. **Run the offline suite:**
   `python3 -m pytest mlb_sgp/tests/ kalshi_common/tests/ kalshi_mlb_mm/tests/ kalshi_mlb_rfq/tests/ -q`.
   Goldens (`test_golden_baselines.py`) must pass for every active book.
7. **Sweep check.** Every active book should write rows with its own
   `*_direct` source label and no `PARSE TRIPWIRE` lines. A book that raises
   keeps its previous rows — that is #33 working, not a failure to fix.
8. **Known limits to carry into the restart:** no book prices spread × ML;
   3-leg shapes return incoherent numbers and must not be quoted; expect a
   4-of-6 quorum at any given rung.

### Goldens

`mlb_sgp/tests/test_golden_baselines.py` replays frozen per-book payloads
through each real `price_sgps` and pins the resulting rows. Fully offline — a
companion test poisons every HTTP transport to prove it.

Recapture only when a book's call pattern changes on purpose, and only from a
slate where that book is healthy:

```bash
python3 mlb_sgp/tests/capture_goldens.py --books novig
```

Capturing mid-repair bakes a broken book's output in as the baseline. The
capture refuses to write a fixture for a book that produced zero rows.

Fixtures are gzipped (`<book>_golden.json.gz`). A raw capture is up to 4.5 MB
of market tree per book — 7.3 MB across the set — and every recapture would add
another blob; compressed they total ~220 KB.

**Caesars has no golden fixture yet.** Capturing one needs a single healthy
live pass, and its WAF was throttling this IP for the whole capture window (see
above). `test_every_book_has_a_golden` fails until it is captured — that is the
gate correctly reporting a gap, not a broken test. Fix it with:

```bash
python3 mlb_sgp/tests/capture_goldens.py --books caesars
```

This replaced two earlier "golden" tests for DK and FD that re-ran the real
scraper against the live book, wrote into the shared production DuckDB from a
test run, and diffed against a CSV keyed by `game_id` — so they went stale
within hours and could never pass on an ordinary day.

## Files

| File | Purpose |
|------|---------|
| `verify_books.py` | Six-book live verification gate (issue #42) — read-only |
| `tests/golden_replay.py` | Record/replay harness behind the golden baselines |
| `tests/capture_goldens.py` | Hand-run capture of frozen golden fixtures |
| `scraper_draftkings_sgp.py` | DK SGP scraper shim (calls `draftkings.price_sgps`) |
| `scraper_fanduel_sgp.py` | FD SGP scraper shim (calls `fanduel.price_sgps`) |
| `scraper_prophetx_sgp.py` | ProphetX SGP scraper shim |
| `scraper_novig_sgp.py` | Novig SGP scraper shim |
| `draftkings.py` / `fanduel.py` / `prophetx.py` / `novig.py` | Per-book orchestrators (`price_sgps`) |
| `dk_client.py` / `fd_client.py` / `prophetx_client.py` / `novig_client.py` | Per-book HTTP clients |
| `_shared.py` | `TargetLine` / `PricedRow` dataclasses, `load_target_lines`, `upsert_priced_rows`, decimal/american helpers, `BookTransportError` + `PriceCallTally` (see "Failure semantics"), `request_with_retry` + the two `RetryProfile`s (see "Retry & backoff") |
| `scraper_pikkit_mlb.py` | Pikkit MLB SGP scraper (fallback) |
| `pikkit_common.py` | Reusable Pikkit functions |
| `recon_draftkings_sgp.py` | DK network recon tool |
| `quick_recon.py` | Lightweight CDP recon (attaches to running Chrome) |
| `db.py` | DuckDB helpers for mlb_sgp_odds |

## Troubleshooting

**"No mlb_parlay_opportunities table"** — Run the MLB pipeline first (`cd "Answer Keys" && python run.py --sport mlb`).

**All games return "no price"** — curl_cffi session may have expired. The scraper auto-reinits after 3 consecutive failures, but if all games fail, try running again.

**`BookTransportError` in `bot.log`** — the book is genuinely unreachable, not
quiet. Read `stage=` to localize it: `auth` (MGM accessid / CZR WAF token),
`events` (the fixture listing), `structure` (a per-event market fetch), `price`
(every price call in the cycle failed — check whether that book prices on a
different host). The book's previous rows were preserved, and its client is torn
down and rebuilt after 3 consecutive failures. See "Failure semantics" above.

**"Total X not found in DK selection IDs"** — Wagerzon total doesn't exist on DK for this game. Rare — DK offers totals from 5.0 to 13.0+.

**"Spread ±1.5 not found"** — Game may have been removed from DK or SGP not yet available (too far from game time).
