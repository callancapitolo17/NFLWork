# FanDuel Full Market Coverage — Design

## Review Pack

**What we're building** — The FanDuel singles scraper currently sees only the
markets in FanDuel's same-game-parlay tab, so it silently misses First-7-Innings
and First-3-Innings totals/run-lines (and per-inning lines) that FanDuel *does*
post on its default tab. We're (1) making the scraper fetch **both** tabs and
merge them, and (2) replacing FanDuel's brittle exact-name whitelist with a
DraftKings-style keyword classifier so FanDuel picks up every market type we
already price — and auto-picks-up new ones FanDuel adds later. Goal: FanDuel
pill coverage on the bets tab matches DraftKings wherever FanDuel posts the
market.

**Key decisions**

1. **Merge both tabs, not switch.** Decision: fetch the default tab *and* the
   SGP tab and union the markets/runners. Rejected: switching to the default
   tab. Why: neither tab is a superset — FG-alts + all of F5 live *only* in the
   SGP tab; F7/F3 lines live *only* in the default tab. Verified live across
   both games.
2. **Keyword classifier over extended whitelist.** Decision: port DK's
   `classify_market` (period + market-type detection + exclusion list) to
   FanDuel. Rejected: just adding F7 names to the existing whitelist. Why: the
   whitelist is whack-a-mole — it caused *this* bug. The keyword approach gives
   DK parity and future-proofs against new FanDuel markets. (User explicitly
   asked for "all available FanDuel markets.")
3. **Validation gate before merge.** Decision: prove the new classifier accepts
   the right markets and rejects all junk by running it against the live
   FanDuel market list for several real games and eyeballing the accept/reject
   split. Why: FanDuel posts ~150 markets/event with lots of keyword-colliding
   junk (parlays, bands, "Listed" ML variants, team totals) — the original
   whitelist existed *specifically* to avoid this. Empirical validation is the
   safety net that lets us take the more powerful approach.
4. **Python-only change.** Decision: touch only `fd_client.py` +
   `scraper_fanduel_singles.py`. Why: `get_fd_odds` (Tools.R) already passes
   `period` through generically via the same `.singles_market_name` helper DK
   uses for its F7 rows, and the model already prices F3/F5/F7. New `period="F7"`
   rows flow through the existing plumbing with zero R changes.

**Risks / push back here**

- **Junk misclassification is the real risk.** A keyword classifier is more
  permissive than a whitelist. If FanDuel posts a market we don't anticipate
  that keyword-matches "total"/"run line"/"moneyline" but isn't a clean game
  line, it could write a bad pill. Mitigation is the validation gate (decision
  #3) — but it's validation against *today's* market set, so a genuinely novel
  FanDuel market name in the future could still slip in. Acceptable? My read:
  yes, because the same risk already exists for DraftKings and hasn't bitten,
  and a wrong pill is visible/low-stakes (it's display + sizing input, not an
  auto-bet). Push back if you'd rather keep a belt-and-suspenders allow-list of
  period prefixes on top of the keyword classifier.
- **Doubling tab fetches.** Coverage requires reading two payloads per event
  instead of one. ~26 events → ~52 HTTP fetches/cycle (same as today's
  redundant double-fetch if we refactor; see Architecture). Well under
  FanDuel's observed rate limits, but it is more traffic.
- **Scope boundary.** Per-inning individual markets (1st–9th inning run
  line/total) and 3-way "First N Innings Result" markets are explicitly **out**
  — no dashboard card or model prices them. Push back if you want those too
  (separate, larger effort: new card types + pricing).

**Worth understanding** (opt-in)

- **Whitelist vs. classifier** — this is the allow-list vs. rule-based-filter
  tradeoff. The whitelist is like an R named-vector lookup keyed on exact
  strings (`markets[name]`): unmissable precision, zero generalization. The
  keyword classifier is like a series of `grepl()` rules with an exclusion
  `filter()` first: it generalizes to unseen inputs but you must reason about
  false positives. DK chose rules; FanDuel chose lookup; we're converging
  FanDuel onto rules.
- **Tabs as server-side market filters** — FanDuel's `event-page?tab=` is
  effectively a server-side `WHERE category = ?` over the same underlying
  market set. Different tabs return overlapping-but-different slices; there is
  no `tab=all`. That's why coverage = union of slices, not a single call.

---

## Background

`mlb_sgp/scraper_fanduel_singles.py` writes one DuckDB table
(`fd_odds/fd.duckdb::mlb_odds`) of FanDuel single-leg lines, consumed by
`Tools.R::get_fd_odds()` → `scraper_to_canonical()` → the bets-tab pills.

The scraper gets its markets from `fd_client.FanDuelClient`, whose
`fetch_event_markets` / `fetch_event_runners` both call FanDuel's
`sbapi/event-page` endpoint with a hardcoded `&tab=same-game-parlay-`
(`fd_client.py:163`). That tab omits the First-7-Innings and First-3-Innings
cumulative line markets. Confirmed live (2026-05-23, event 35639118
White Sox @ Giants):

- `tab=same-game-parlay-` → 156 markets, includes FG main+alt, F5 main+alt, but
  **no** `First 7 Innings Total Runs` / `Run Line`, **no** F3 lines, **no**
  per-inning lines.
- `tab=` (default) → 99 markets, includes `First 7 Innings Total Runs`
  (Over −128 / Under +104), `First 7 Innings Run Line`, F3 lines, per-inning
  lines — but **no** FG-alts and **no** F5 markets.

39 markets exist in default-but-not-SGP; the two sets overlap on FG main only.

## Goals / Non-Goals

**Goals**
- FanDuel pill coverage on the bets tab matches DraftKings wherever FanDuel
  posts the market: FG main+alt, F5 main+alt, F7 main (total + run line), F3
  main (total + run line).
- Classifier auto-picks-up future FanDuel markets that map to a priced card.
- Zero regressions to existing FG/F5 FanDuel pills.

**Non-Goals**
- Per-inning individual markets (no card/model).
- 3-way "First N Innings Result" markets (would need tie-drop + re-devig).
- Any R / dashboard / model change.

## Architecture

Two layers, both in Python.

### Layer 1 — Acquisition (both-tab merge), `fd_client.py`

FanDuel's `event-page` payload contains **both** markets and runners in one
response. Today the client fetches that same URL **twice** (once in
`fetch_event_markets`, once in `fetch_event_runners`) — redundant.

Refactor to a single combined fetch per tab:

```
fetch_event_page(event_id, tab) -> (list[Market], list[Runner])
    # one HTTP GET, walk the payload once, return both
```

The scraper calls it once per tab over `TABS = ("", "same-game-parlay-")`,
then merges:
- markets: dedup by `market_id` (first-seen wins; same market in both tabs is
  identical)
- runners: dedup by `runner_id`

Net effect: **2 HTTP fetches/event** (one per tab) — same count as today's
redundant pair, but now covering both tabs. Existing
`fetch_event_markets` / `fetch_event_runners` are kept as thin wrappers
(`tab` defaulting to `same-game-parlay-`) so nothing else breaks; the SGP
scraper, which uses different functions, is untouched.

### Layer 2 — Classification, `scraper_fanduel_singles.py`

Replace `_FD_MARKET_WHITELIST` + `classify_market` with a FanDuel port of
DraftKings' keyword classifier (`scraper_draftkings_singles.py:71`):

```
classify_market(name, home_team, away_team) -> (period, market_type) | None
```

1. **Exclusions first** (return None). Port DK's list (`team total`,
   `player/prop/futures/race to/winning margin/correct score/odd/even/...`,
   single-inning regex `_SINGLE_INNING_RE`) **plus** FanDuel-specific junk that
   keyword-collides with real lines:
   - `parlay` — `… Run Line / Total Runs Parlay`, `Line / Total Parlay N`
   - `listed` — `Moneyline Away/Home/Both Listed` (pitcher-conditional; we keep
     plain `Moneyline`)
   - `(bands)` / `bands` — `Total Runs (Bands)`
   - `tri-bet`, `specials`
   - **Team totals via event teams**: if the name contains `home_team` or
     `away_team`, it's a per-team market (`<Team> Total Runs`,
     `<Team> Alt. Total Runs`) → exclude. (This is why the signature takes the
     teams — cleaner than DK's static prefix map since the client already has
     them.) **Backward-compat:** `home_team`/`away_team` default to `None`; the
     team-total check is skipped when absent, so existing single-arg callers
     and the current `test_classify_market_*` tests keep working unchanged.
2. **Period detection**: `first 7 innings`→F7, `first 5 innings`→F5,
   `first 3 innings`→F3, else FG. (FanDuel uses "First N Innings"; the
   single-inning regex already filtered singular "Inning".)
3. **Market-type detection** (alt before main): `alternate`+`run line`→
   alternate_spreads; `alternate`+`total`→alternate_totals; else `run line` /
   `moneyline` / `total`→main.

The existing `parse_runners_to_wide_rows` (alt-name regex, sign-symmetry guard,
paired-side guard, bucketing) is unchanged — it already handles any
`(period, market_type, line)` and is book-agnostic.

## Data flow

```
scrape_singles()
  for each event:
    markets, runners = {}                 # merged across tabs
    for tab in ("", "same-game-parlay-"):
      m, r = client.fetch_event_page(event_id, tab)
      merge m into markets (dedup market_id)
      merge r into runners (dedup runner_id)
    market_meta = { market_id: classify_market(name, home, away) ... if not None }
    rows = parse_runners_to_wide_rows(event, runners, market_meta, fetch_time)
  write_to_duckdb(all_rows)               # unchanged: atomic swap, empty-guard
```

`get_fd_odds()` then emits `period = row$period` for the new F7/F3 rows; the
bets tab renders them on the existing F7/F3 cards. No R change.

## Error handling

- Per-event isolation already wraps the fetch in try/except (one event's
  failure doesn't tank the scrape). A tab returning non-200 returns
  `([], [])` and we still get the other tab's markets — partial coverage beats
  none.
- Empty-scrape guard in `write_to_duckdb` is unchanged (leaves prior snapshot
  if zero rows).
- Sign-symmetry / paired-side guards in the parser are unchanged.

## Testing & validation

Verification is layered from most-decisive (hard pass/fail) to confirmatory.

**Test invocation (confirmed working):** the worktree has no venv; run with the
main repo's venv pointed at the worktree code, with `mlb_sgp/` on the path:
```
PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python \
    -m pytest mlb_sgp/tests/test_fanduel_singles.py -q
```
Baseline is currently green (2 passed) — that's the regression floor.

1. **Classifier unit test — gold standard, written TDD-first.** Extend
   `mlb_sgp/tests/test_fanduel_singles.py` with a table-driven test:
   a list of `(market_name, home, away) → expected (period, market_type) | None`
   covering **every class** — one real line per period/type that must be
   ACCEPTED, and one example per junk family that must be REJECTED
   (`… Total Runs Parlay`→None, `Total Runs (Bands)`→None,
   `Moneyline Away Listed`→None, `<away_team> Total Runs`→None,
   `7th Inning Total Runs`→None, `First 7 Innings Result`→None,
   `First 7 Innings Total Runs`→("F7","main"), `First 7 Innings Run Line`→
   ("F7","main"), etc.). Write these assertions **before** implementing the
   classifier; "done" = all green. This turns the old "eyeball it" into a hard
   pass/fail and locks the behavior against future edits.
2. **Full live-market audit (catches names the unit test didn't anticipate).**
   A throwaway script (job dir, not committed) runs the finished classifier
   against the *full* live market list (both tabs) for ≥3 real games and prints
   ACCEPTED (name → period/market_type) vs REJECTED. Confirm the accept set is
   exactly the regular FG/F5/F7/F3 lines and the reject set contains all junk.
3. **End-to-end ground truth.** Re-pull a fresh FanDuel F7 quote at test time
   (odds drift, so don't hardcode — during design the Giants F7 total read
   Over −128 / Under +104, F7 run line WSX −134 / Giants +110). Run the scraper
   for real, then query `fd.duckdb::mlb_odds` and confirm F7 + F3 rows now exist
   for both games with prices matching the fresh quote.
4. **Regression diff (no FG/F5 loss, no double-count).** Snapshot current
   `fd.duckdb` coverage as a baseline *before* the change (row counts grouped by
   `period × market`); re-snapshot after. FG/F5 counts must be ≥ baseline (we
   only add), and there must be zero duplicate `(period, market, total/spread)`
   rows (proves the both-tab dedup is correct).
5. **FD-vs-DK parity (the actual goal).** For each game, list cards where DK has
   a pill; FanDuel should now have one wherever FanDuel posts that market. The
   only remaining gaps should be the known structural ones — F7/F3 alts and
   F7/F3 2-way ML, which FanDuel genuinely does not post.
6. **Timezone parity gate (required by repo CLAUDE.md):** run
   `tests/timezone_parity_test.py` — it cross-checks each scraper's
   `game_start_time` against Odds API `commence_time`. Merging tabs must not
   perturb `game_start_time`.
7. **Visual confirmation:** the two cards from the original screenshot
   (Cardinals @ Reds Over 7.5, White Sox @ Giants Over 6.5) render FanDuel pills
   in the dashboard.

**Honest limits of this verification:** items 1–2 validate against *today's*
market set, so a genuinely novel FanDuel market name in the future could still
slip through the keyword classifier (this is the accepted risk in the Review
Pack, mitigated by the team-total + junk exclusions and the fact that a wrong
pill is a visible display/sizing input, not an auto-bet).

## Version control

- **Branch / worktree:** `worktree-fanduel-full-market-coverage` (already
  created at `.claude/worktrees/fanduel-full-market-coverage`). All work here;
  test from the worktree; merge to `main` only after the validation gate + a
  fresh scraper run pass and you approve. Clean up the worktree + branch after
  merge.
- **Files modified:**
  - `mlb_sgp/fd_client.py` — add `fetch_event_page(event_id, tab)`; make the
    two existing fetchers thin wrappers with a `tab` param.
  - `mlb_sgp/scraper_fanduel_singles.py` — both-tab merge loop; replace
    whitelist with keyword `classify_market`; thread `home_team`/`away_team`
    into classification.
- **Commits:** (1) fd_client both-tab fetch; (2) scraper classifier + merge;
  (3) docs + memory. Squash-friendly but kept separate for review clarity.

## Documentation

Updated in the same merge:
- `mlb_sgp/README.md` — FanDuel scraper now fetches both the default and SGP
  tabs and uses a keyword classifier (note the F7/F3 coverage and the junk
  exclusions).
- Memory `fd_event_page_tab_coverage_gap.md` — mark resolved, point to this
  spec and the merge commit.
- Verify nothing in `Answer Keys/CLAUDE.md` needs a note (no R change, so
  likely not — confirm during review).
