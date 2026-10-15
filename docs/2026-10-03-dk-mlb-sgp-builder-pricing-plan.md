# DK same-game MLB pricing over the SGP-builder GET — plan

Status: BUILT 2026-10-04 on `feature/dk-mlb-sgp-builder-pricing`; not merged.

**As built — one change from this plan.** `kalshi_rfi` prices 1-leg sets at
DK, and the builder's 1-leg price is not DK's straight price (1st-inning total
+125/-160 straight vs +120/-155 in the builder, 0.85 pp after devig). Declining
1-leg sets would also have tripped the `prices_empty` parser tripwire on every
RFI fetch. So `fetch_selection_ids` now returns DK's straight odds
(`straight_odds`) and `resolve_legs` puts them on each leg: 1-leg fairs devig
the straight two-way with zero SGP calls, and DK's Route B uses straight
singles. The "DK Route B" follow-up in section 3 is therefore done. Live check
2026-10-04: 16/16 sets priced through the real hooks, 0.59 pp mean from the
other books' median.

## 1. Findings (live, 2026-10-03, MLB postseason slate)

**Endpoint.** `GET sportsbook-nash.draftkings.com/sites/US-WV-SB/api/sportscontent/sgp/dkuswv/sportsdata/v2/sgp`
with header `X-SportId: 7` (baseball; `sportId` on both the league listing and
`parlays/v1/sgp/events/{id}` payloads). Same request shape as NFL
(`nfl_specials/dk_book.py`), curl_cffi `impersonate="chrome"`, no warm-up.

**Market ids.** DK's market id is `{type}_{mnum}`: `1_` moneyline, `2_` run
line, `3_` total. Selection ids already carry the mnum (`0ML{mnum}_1`,
`0HC{mnum}N150_1`, `0OU{mnum}O850_1`), so the market id of every selection
`resolve_legs` returns derives from its id. Verified on main, alternate, F5 and
1st-inning markets (e.g. `0OU86560036O50_1` -> `3_86560036`, "Runs - 1st Inning").

**Accuracy.** 21 leg sets across 3 games (2- and 3-leg; FG/F5; main and alt
lines; ML+total), each with >= 2 other books priced through the production
`SGPService.price_on_demand` (FD/Novig/MGM/CZR), DK partition from the GET,
both devigged by `fair_value.devig_partition`:

| metric | value |
|---|---|
| mean / median / max abs(DK − other-books median) | 0.63 / 0.52 / 1.69 pp |
| DK inside the other books' [min, max] | 14 / 21 (rest within 0.73 pp) |
| DK partitions with a missing cell | 0 / 21 |

The 1.4–1.7 pp gaps are all on MIL@SD, where Novig and MGM disagree with each
other by 2.5–5 pp; DK sits between them.

**Semantics (each probed directly).**
- Base + one candidate market returns the price of base + EACH outcome of that
  market: one GET prices both sides of the last leg.
- Empty base prices nothing (`trueOdds: null`) -> a 1-leg set cannot be priced.
- Base order is irrelevant (same cells either way).
- Two sides of one market in the base: DK silently drops one
  (`selectionsNotMapped: IncompatibleSelections`) and prices the rest — must
  be a decline, never a price.
- No candidate -> 400 ValidationError; unknown selection -> 400 InvalidSelection.

**Burst tolerance** (one variable per run, control first and last, stop on 403/429):

| run | shape | calls | result |
|---|---|---|---|
| A control | sequential, 0.5 s gap | 20 | 20 × 200 |
| B | sequential, no gap (4.5/s) | 60 | 60 × 200 |
| C | sequential, no gap, sustained | 89 | 88 × 200, 1 hang (20 s, 0 bytes) |
| D control | sequential, 0.5 s gap | 5 | 4 × 200, 1 hang (20 s) |
| E, F control | sequential, 0.5 s gap | 100 | 100 × 200 |
| G | sequential, no gap, sustained 73 s (4.1/s) | 300 | 300 × 200 |
| H | concurrency 2 (5.5/s) | 120 | 119 × 200, 1 hang (8 s timeout) |
| Z control | sequential, 0.5 s gap | 20 | 20 × 200 |

~750 calls, zero 403/429. Concurrency above 2 deliberately untested.

**Latency tail.** p50 0.2 s, p95 0.35 s, but ~1% of calls take 2–8 s and
~0.4% hang >= 8 s. Not a throttle: the slow requests replay in 0.2 s, and the
same request slow once is fast on the next call (DK-side cold computation).

**Demand.** Maker's last run (2026-09-13..20, `sgp_fetch_health`): DK
on-demand p50 8 / p95 39 / max 101 pricings per active minute. All 77k ended
`empty`: the maker still fires ~15k doomed `calculateBets` POSTs a day at DK's
Akamai-gated host because `draftkings` is in `SGP_BOOKS`.

## 2. Design

Scope: `mlb_sgp/draftkings.py::price_selection_set` (the on-demand path —
maker and `kalshi_rfi`). Nothing else changes behavior.

- **Request.** `price_selection_set(client, refs, dk_event_id, *, counters)`:
  base = `refs[:-1]`, `marketCandidates` = market id of `refs[-1]`; returns
  that selection's `trueOdds` from `compatibleMarkets`. Same sibling-module
  convention as MGM's `price_selection_set(client, refs, fixture_id)`.
- **Declines (None).** `len(refs) < 2` (no wire call); `selectionsNotMapped`
  non-empty or mapped != base; `refs[-1]` absent from `compatibleMarkets`;
  any non-200 that is not 403/429. 403/429 stays a counted
  `BookTransportError` (issue #102 classification kept).
- **Market id.** `dk_market_id_of(selection_id)`: `0ML` -> `1_`, `0HC` -> `2_`,
  `0OU` -> `3_` + mnum; any other prefix raises ValueError (counted parse
  failure, never a guessed market).
- **One price lane per client.** A lock on `DraftKingsClient` holds DK to one
  in-flight price call and >= 0.5 s between calls
  (`MIN_SECONDS_BETWEEN_PRICE_CALLS = 0.5`, DK's own widget throttle, same
  constant as `nfl_specials/dk_book.py`). This is the "cap DK's concurrency to
  1" condition the README set for any re-enable. A caller that cannot get the
  lane within 3 s declines (the RFQ just has one fewer book). Per-call timeout
  3 s, `RETRY_LIVE` unchanged (worst lane hold ~6.5 s on a hang).
- **Sibling-cell reuse.** The answer for (event, base, candidate market) is
  kept on the client for 2 s and looked up under the lane lock, so the two
  cells that differ only in the last leg cost one GET: 2 GETs per 2-leg
  partition (was 4 POSTs), 4 per 3-leg. 2 s spans one flight's cells, so a
  quote never prices off another RFQ's stale number.
- **Hook.** `sgp_service._book_on_demand_hooks` passes `event["dk_event_id"]`.
- Capacity at 0.5 s: 120 GETs/min, enough for p95 demand (~80–100/min); at
  peak minutes DK declines the overflow and the other books carry the quote.

## 3. Out of scope (follow-ups, not in this branch)

- **Sweep path** (`price_sgps` -> `calculate_sgp`, used by the taker's
  `refresh()` and the dashboard parlay tab) stays on `calculateBets` and keeps
  403ing. Same GET would fix it; separate ticket.
- ~~DK Route B~~ — built (see "As built" above).
- **Concurrency > 1.** Measured clean at 2 for 22 s only; not used.

## 4. Version control and worktree

- Worktree: this session's (`.claude/worktrees/zen-heisenberg-4a05f5`), branch
  `feature/dk-mlb-sgp-builder-pricing` off local main.
- Commit 1 — `mlb_sgp/draftkings.py` (`price_selection_set`, `dk_market_id_of`),
  `mlb_sgp/dk_client.py` (lane lock, last-call time, 2 s answer cache),
  `kalshi_common/sgp_service.py` (hook passes event id).
- Commit 2 — tests, no network: rewrite
  `mlb_sgp/tests/test_dk_on_demand_price_status.py` for the GET (403/429
  counted, 400 clean, dropped base -> None, missing candidate -> None, n=1 ->
  no call); new `mlb_sgp/tests/test_dk_sgp_builder_price.py` (market-id
  derivation incl. unknown prefix, sibling-cell reuse inside 2 s and refetch
  after, lane serializes two threads >= 0.5 s apart via injected clock);
  update `test_dk_price_host_403.py::test_on_demand_price_path_shares_the_same_url_constant`
  (the wager URL constant stays for the sweep/trifecta).
- Commit 3 — docs (section 5).
- Before merge: full `mlb_sgp/tests` + `kalshi_common` tests; live check from
  the worktree — ~20 leg sets through the real `SGPService.price_on_demand('draftkings', ...)`
  next to the other books on a live slate; pre-merge review of
  `git diff main..HEAD`; then merge only on your yes, restart the maker, remove
  the worktree and branch.

## 5. Documentation

- `mlb_sgp/README.md`: endpoint table (the SGP-builder row said
  "inaccessible" — wrong), § DraftKings price host (on-demand now prices via
  the GET; sport id 7; burst and tail numbers; lane and reuse), § "Decision:
  DK stays OFF" -> superseded for on-demand.
- `kalshi_mlb_mm/README.md`: DK back in the same-game consensus.
- `NFLWork/CLAUDE.md`: one clause on the maker bullet (DK prices same-game via
  its SGP-builder GET; `calculateBets` is placement-only).

## 6. Rollout and kill switch

- No new flag: DK is already in `SGP_BOOKS`; a maker restart picks it up.
- Watch `sgp_fetch_health` for book `draftkings`: `ok` share, `transport_errors`.
- Any 403/429 on the GET: drop DK with the existing `SGP_BOOKS` env and
  restart — before the rule can spread to DK's read endpoints, which feed the
  cross-game leg surface. Rollback is `git revert`.
