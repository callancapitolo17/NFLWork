# Novig MLB SGP Scraper — Plan

**Status:** Pre-build. Recon complete but has one capture gap (parlay response body) that must be closed before writing the scraper.
**Branch:** `feature/prophetx-novig-sgp-recon` (same branch as ProphetX work)
**Worktree:** `.worktrees/prophetx-novig-recon`
**Target merge:** `main`

---

## Goal

Add Novig as a fourth SGP price source alongside DK, FD, and ProphetX. Writes to `mlb_sgp_odds` with `source='novig_direct'`, blends in `mlb_correlated_parlay.R` with per-game devig (same pattern as the other three).

## Edge thesis — what Novig uniquely adds

Live probes on 2026-04-24 against `/nbx/v1/parlay/request/unauthenticated` revealed the actual pricing lineage:

- **Leg prices come from DraftKings.** Every leg in 6 sample calls returned `vendor: "DRAFTKINGS"`. Novig does not originate its own leg prices — it consumes DK's.
- **BUT — Novig applies its own correlation model on top.** Combined parlay prices are NOT naive independent multiplications. Observed samples:

    | Combo | Leg1 × Leg2 (naive) | Novig parlay | Ratio |
    |---|---|---|---|
    | DET -1.5 + Over 8.5 | 0.2364 | 0.266 | 1.125× (more generous) |
    | DET -0.5 + Over 4.5 | 0.2825 | 0.307 | 1.087× (more generous) |

    Both combos are fav + over (positive correlation expected) — Novig gives better odds than independent multiply, exactly the direction a correlation-aware model should push.

- **Therefore, what Novig adds to the blend = the delta between Novig's correlation model and DK's `calculateBets` correlation model**, both applied to the same DK leg prices. If the two correlation layers differ even modestly, Novig carries real signal DK doesn't.

- **What Novig is NOT:** a sharp benchmark, a different provider for leg prices, or a Pinnacle-equivalent. Adding Novig is betting that Novig's correlation algorithm diverges non-trivially from DK's.

- **Hidden risk:** Novig's legs depend on DK being available. If DK has a pricing outage, Novig likely returns stale prices or errors — scraper must handle missing/old data gracefully.

## Prior work / reference

- **Recon**: `mlb_sgp/recon_novig_sgp.py` + JSON dump at `mlb_sgp/recon_novig_sgp.json` (3292 lines, 1.7 MB, 6 phases).
- **Analysis**: confirmed parlay-quote endpoint `POST /nbx/v1/parlay/request/unauthenticated`, GraphQL market tree `POST /v1/graphql` op `EventMarkets_Query`, outcome-UUID selection IDs, no auth required.
- **Reference scraper**: `mlb_sgp/scraper_prophetx_sgp.py` (just built, same architectural pattern).

## Design decisions

### D1. Request format
**Decision:** `POST https://api.novig.us/nbx/v1/parlay/request/unauthenticated` with payload `{"outcomes":[{"id":"<uuid>"}, ...], "boostId": null}`. Each leg is a single outcome UUID — no marketId / lineId / line fields like ProphetX. Cleanest selection encoding we've seen across any book.

**Rationale:** captured verbatim from recon (3 live calls in `add_leg2`).

### D2. Response parsing — **NEEDS LIVE PROBE FIRST**
**Decision:** Before writing the scraper, `curl` the endpoint live once with a UUID pair already in the recon to capture the response body. The recon JSON's `response_preview` was empty for all three `201` responses — a capture gap, not a missing response.

**Rationale:** we literally don't know the response shape — American odds? decimal odds? implied probability (which is Novig's convention elsewhere: `"available":0.514`)? Single scalar or an offer ladder? Writing the parser without this is blind guessing.

**Mitigation if response is unexpected:** `available: <float>` is Novig's pattern elsewhere, so our best guess is the parlay endpoint returns a similar implied-probability scalar. Scraper will be written defensively — it tries multiple field paths (`odds`, `price`, `available`, `impliedProbability`) and logs the full response at first call.

### D3. Market discovery — GraphQL
**Decision:** `POST https://api.novig.us/v1/graphql` with op `EventMarkets_Query` and `{"variables":{"eventId":"<uuid>", ...}}`. Returns the full market tree including outcome UUIDs.

**Prerequisites before scraper can run:**
1. Capture the exact GraphQL query body (operation name, variables, query string). Recon shows the op name but query text is truncated in the preview — need to extract via a live run or by inspecting the `app.js` bundle.
2. Enumerate MLB events: likely via `LiveEventTicker_Query` GraphQL op scoped to `league: "MLB"`.

### D4. F5 support — UNKNOWN
**Decision:** Attempt to fetch F5 markets. If the market tree response has no F5 markets for MLB (as in the recon sample), log and continue with FG only. Novig's schema DOES include `MONEY_1H`/`SPREAD_1H`/`TOTAL_1H` for NBA — segment markets exist as a concept — but the MLB sample in recon had only player props and standard team markets.

**Open question:** whether Novig offers MLB F5 markets at all. Reconnect answered: *can't tell from dump* (10 KB response-preview truncation, and only player-prop markets appeared in the visible portion).

**Mitigation:** scraper treats F5 as optional — runs FG-only if F5 markets aren't found. Better-than-nothing coverage.

### D5. Pricing to use — PROBABILITY → AMERICAN CONVERSION
**Decision:** If response body confirms implied-probability format (e.g., `{"available": 0.514, ...}`), convert to decimal via `1 / p` and to American via standard formulas, then write to `mlb_sgp_odds` in the same (decimal, American) columns as the other books. No schema change.

### D6. Vig default
**Decision:** `NOVIG_SGP_VIG_DEFAULT = 1.10` initially — same starting value as DK and ProphetX.

**Rationale:** Novig's SGP is third-party-priced. 1.10 is the "assume retail-book-like vig" conservative prior. Scraper logs measured per-game vig (same `<-- high vs default` pattern); tune after a few runs.

### D7. Auth
**Decision:** Anonymous, no cookies required. Endpoint literally named `/unauthenticated`; recon confirms 201 response without any auth-required errors.

**Defensive check:** scraper still tolerates 401/403 with a clear error (per the ProphetX M2 pattern) in case Novig ever adds a soft rate-limit disguised as auth.

### D8. Selection ID persistence
**Decision:** Novig outcome UUIDs are considered **stable within a session but not across sessions or line moves.** Scraper fetches the market tree and submits the RFQ in quick succession (same run), same idiom as ProphetX's `lineId` handling.

### D9. Match to canonical game_ids
**Decision:** Same team-name + UTC-hour-bucket canonicalization as DK/FD/PX. Novig's event object exposes team display names that `resolve_team_names` should handle (e.g. "Detroit Tigers", "Cincinnati Reds").

## Architecture

```
mlb_correlated_parlay.R
    |
    |— invokes scraper_draftkings_sgp.py
    |— invokes scraper_fanduel_sgp.py
    |— invokes scraper_prophetx_sgp.py
    |— invokes scraper_novig_sgp.py          [NEW]
    |       |
    |       |— POST /v1/graphql op LiveEventTicker_Query → MLB events
    |       |— match to canonical game_ids
    |       |— for each game: POST /v1/graphql op EventMarkets_Query
    |       |— for each (game, period, combo): POST /nbx/v1/parlay/request/unauthenticated
    |       |— write rows to mlb_sgp_odds (source='novig_direct')
    |
    v reads mlb_sgp_odds WHERE source IN (…'novig_direct')
    |
    |— novig_vig_lookup, novig_fair_prob in blending (same idiom as PX)
```

## Files to create / modify

### New
- `mlb_sgp/scraper_novig_sgp.py` (~400 lines — probably slightly shorter than ProphetX because selection encoding is simpler)
- Memory: `novig_sgp_scraping.md` (after first live run)

### Modify
- `Answer Keys/mlb_correlated_parlay.R`:
  - Add `NOVIG_SGP_VIG_DEFAULT <- 1.10` near other vig defaults
  - Add `system2(..., "scraper_novig_sgp.py")` invocation
  - Add `'novig_direct'` to SGP load WHERE clause
  - Add `nv_sgp`, `nv_vig_lookup`, blending block, tibble column — mirrors the PX diff
- `MEMORY.md`: add `novig_sgp_scraping.md` pointer

## Pre-build probes (before writing scraper)

Because D2 (response body) and D3 (GraphQL query text) are unknowns, the scraper **cannot be written blindly**. Two short live probes first:

### Probe P1: parlay response body
```
curl -X POST "https://api.novig.us/nbx/v1/parlay/request/unauthenticated" \
  -H "Content-Type: application/json" \
  -d '{"outcomes":[{"id":"265536a0-206d-41c7-bd17-deb66c3af804"},{"id":"1ceae72d-4cef-4c78-a6f8-ef2bb7df42d0"}],"boostId":null}'
```
(UUIDs come from the recon's `add_leg2` phase — may be stale but still informative.)
**Success criteria:** a JSON body that reveals the pricing field (probability? odds?). Even a 400 "stale outcome" error tells us the shape of errors.

### Probe P2: EventMarkets_Query full query text
Grab from Novig's JS bundle (`https://novig.com/` and crawl `_next/static/chunks`) or from a fresh recon run with the response-body capture filter widened to include 201 + status>=200.

**Alternative:** rerun `recon_novig_sgp.py` with a small patch to `handle_response` so it captures ALL statuses, not just 200. One-line fix.

## Test plan

### T1. Probe P1 + P2 succeed → we can write the scraper
Document the response body shape in `novig_sgp_scraping.md`. Block scraper development until this is done.

### T2. Scraper standalone smoke
```
cd mlb_sgp && venv/bin/python3 scraper_novig_sgp.py --verbose
```
Expect priced MLB FG combos. F5 optional. No 401s. Per-game vig logs.

### T3. R pipeline end-to-end
`bash "Answer Keys/MLB Dashboard/run.sh"` — expect `Loaded N Novig SGP odds for blending` + `novig_fair_prob` column in the parlay tibble.

### T4. Regression check
For 20 sample games, compare blended_prob before/after Novig joins the blend. Any game shifting >10% is an outlier — investigate.

### T5. Cross-provider sanity
If DK/FD also source from OpticOdds (check via a dev poke / network inspection), note in memory that Novig ~= DK/FD for SGP pricing — reduces marginal edge of adding Novig.

## Documentation updates

- **`mlb_sgp/README.md`** (if exists): mention `scraper_novig_sgp.py` alongside the other three.
- **Memory note `novig_sgp_scraping.md`** (new): endpoints (GraphQL + REST parlay), outcome-UUID scheme, response body format, F5 availability for MLB, OpticOdds provider context, Sportico incident reference.
- **`MEMORY.md` index**: add pointer.

## Pre-merge review checklist

- [ ] Data integrity: `clear_source('novig_direct')` at run start; dedupe on `(game_id, combo, bookmaker, source)`.
- [ ] Resource safety: DuckDB connections in `try/finally`.
- [ ] Edge cases: off-season, empty market tree, UUIDs expired between fetch and RFQ.
- [ ] Dead code: none.
- [ ] Log hygiene: stdout only, R suppresses.
- [ ] Security: no secrets (no auth in the first place).
- [ ] **Loud failure:** zero-prices-written alert with actionable diagnosis (auth failure vs matching failure vs RFQ failure) — same pattern as ProphetX M2.

## Worktree lifecycle

Shared with ProphetX work — same worktree, same branch. Cleanup happens after BOTH scrapers are merged together OR if Novig is deferred, after ProphetX-only merge.

1. ✅ Worktree created, branch created
2. ✅ Recon scripts written + committed
3. ✅ ProphetX scraper built (pending live test + commit)
4. 🔲 Probe P1 + P2 (Novig response body + GraphQL query text)
5. 🔲 Scraper written
6. 🔲 R integration
7. 🔲 Live test T2 / T3
8. 🔲 Memory note + README update
9. 🔲 Pre-merge review
10. 🔲 User approves merge
11. 🔲 Merge to `main`
12. 🔲 Worktree cleanup

## Open questions / risks

- **Response body parsing (D2).** Hardest unknown. Close via probe P1.
- **GraphQL query text (D3).** Needed to send valid GraphQL requests. Close via P2 or recon re-run.
- **F5 market availability (D4).** Unknown; scraper degrades gracefully if absent.
- **Provider overlap with DK/FD.** If all three use OpticOdds, Novig adds little diversity to the blend. **Can't confirm without reverse-engineering competitors** — accept as residual risk.
- **Novig's void-parlay history.** The Dec 2025 incident voided legitimate wins. *Does not affect scraping* (we're reading prices, not placing bets) but worth remembering: these prices are sometimes wrong in ways Novig itself retroactively admits.
- **Geofence.** Unknown whether `api.novig.us` enforces geo at the API layer or only on bet placement. Scraping from Florida: probably fine. Scraping from a non-supported state: unknown.
- **Rate limits.** Unobserved. Start with `max_workers=3` (even lighter than ProphetX's 4).

## Not in scope (explicit)

- **Novig singles** (spread/total/ML P2P book). Genuinely sharp P2P liquidity but separate project; focus here is SGP to match the other three sources.
- **Player-prop parlays.** Novig's MLB coverage in the recon was player-prop-heavy; we could in theory RFQ player-prop parlays. Out of scope — the correlated parlay tab is scoped to spread+total combos.
- **OpticOdds direct integration.** Even if Novig is literally reselling OpticOdds, going direct would be a different project (B2B data deal).

## Success criteria

1. Probe P1 + P2 pass — we understand the parlay response shape and can issue valid GraphQL queries.
2. T2 passes: ≥1 MLB FG combo priced on a live run.
3. T3 passes: pipeline runs end-to-end with `novig_fair_prob` visible in the parlays tab.
4. Memory note + README updated.
5. Pre-merge review done.
6. User explicitly approves merge.

## Self-review (holes before presenting)

- [x] Edge thesis is honest (not sharp, ensemble diversity only, provider-overlap risk called out).
- [x] Blocking probes explicitly enumerated with concrete commands.
- [x] Worktree lifecycle covered, shared with ProphetX work acknowledged.
- [x] Failure modes identified: silent-auth, UUID staleness, provider overlap, geofence, void-parlay history.
- [x] Non-goals explicit.
- [x] All decisions include rationale and reversibility notes.
- [ ] **Potential hole:** I'm assuming `LiveEventTicker_Query` is the right event-listing endpoint. Recon confirms it exists and takes a league filter, but the exact variables/response shape is inside the 10 KB preview window — should verify in Probe P2.
- [ ] **Potential hole:** we haven't answered whether Novig MLB SGPs are even available TODAY (off-season risk — 2026-04-24 is mid-season so probably fine). User observed SGP behavior in recon; should spot-check in live run.
