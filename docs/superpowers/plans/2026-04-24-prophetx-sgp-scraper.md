# ProphetX MLB SGP Scraper — Plan

**Status:** Implementation written, pre-test. Awaiting end-to-end verification + merge approval.
**Branch:** `feature/prophetx-novig-sgp-recon`
**Worktree:** `.worktrees/prophetx-novig-recon`
**Target merge:** `main`

---

## Goal

Add a third SGP price source (ProphetX) to the MLB correlated-parlay tab, so every game's 4 combos (home/away spread × over/under, for both FG and F5) are priced by three books — DK, FD, ProphetX — instead of two. The `mlb_correlated_parlay.R` blending step (`mean(model, dk, fd, px)`) gets one more independent data point per combo, tightening the devigged fair-price estimate.

## Why now / edge thesis

- DK and FD are both retail books; their SGP prices share common correlation-model lineage and can drift in the same direction.
- ProphetX is a P2P exchange; its RFQ-priced SGPs come from a *different* pricing regime. A third price that is *not* correlated with DK/FD error is more valuable for the blend than a third retail book would be.
- ProphetX SGPs are NOT a true sharp benchmark (the RFQ is quoted by a third-party market maker, not peer-matched). But adding diversity to the blend still reduces variance of the fair-prob estimate — standard ensemble reasoning.
- Cost is low: the web app exposes a clean `/parlay/public/api/v1/user/request` endpoint that takes structured legs and returns priced offers. No Akamai, no cert pinning, single origin.

Explicit non-goal: building an internal SGP simulator. Deferred per user decision (see memory `sgp_simulator_plan.md`).

## Prior work that informs this

- **Recon** (`mlb_sgp/recon_prophetx_sgp.py`) captured a complete DET@CIN SGP session on 2026-04-24. Output at `mlb_sgp/recon_prophetx_sgp.json`.
- **Analysis** of the recon confirmed: single-origin API, structured selection IDs, explicit F5 markets, working `/parlay/public/api/v1/user/request` endpoint with a 4-tier offer ladder.
- **DK/FD templates**: `scraper_draftkings_sgp.py`, `scraper_fanduel_sgp.py`, and the `mlb_correlated_parlay.R` blending block are the reference patterns we mirror.

## Design decisions

### D1. Which offer in the RFQ ladder counts as "the price"?

**Decision:** `offers[0]` — the tightest available combined American odds.

**Rationale:** Matches DK/FD's one-price semantics so the R blending code stays uniform. The full ladder is retained in memory and surfaced via the `[N tiers]` log line so we can switch to stake-weighted or minimum-stake pricing later without re-fetching. Reversible.

**Alternatives considered:**
- Minimum-stake offer — "what a small retail user actually gets." Rejected for v1: ladder ordering and minimum stake floors aren't documented; using offers[0] is unambiguous.
- Stake-weighted average — more liquidity-aware. Rejected for v1: overkill before we have any data to judge it by.

### D2. Vig default when <4 combos present

**Decision:** `PROPHETX_SGP_VIG_DEFAULT = 1.10` (matches DK; instrument and tune).

**Rationale:** ProphetX is a P2P exchange that takes commission from winnings, not spread-in-price — in theory vig ≈ 1.00. In practice, the RFQ's *tightest* offer is a combined quote from a third-party market maker, which builds in some spread; empirically expect closer to retail book levels. Starting at 1.10 matches DK's default and is a more realistic prior than 1.00 — under-devigging (1.00) makes every PX fair-prob slightly over-confident, which pollutes the blend. 1.10 errs toward caution.

**Instrumentation:** Scraper logs `[FG vig: X.XXXX]` / `[F5 vig: X.XXXX]` per game and flags values materially different from the default (`> 1.12` or `< 1.08`). After 5-10 runs we re-tune the default from the observed distribution (median, or a bimodal split like FD if needed).

### D3. Auth model

**Decision:** Probe anonymously first; fall back to cookies from the Playwright `.prophetx_profile` if the endpoint returns 401/403.

**Rationale:** The recon ran inside a logged-in profile, so we can't distinguish "public endpoint" from "cookie-required endpoint" from that data alone. The path is under `/parlay/public/...` which *suggests* anonymous works, but we must verify. Cookie fallback costs nothing — a `sqlite3` read of the profile's Cookies DB in read-only mode.

**Failure mode:** If both paths 401, the scraper logs an error and writes zero rows. R blending degrades gracefully (uses `PROPHETX_SGP_VIG_DEFAULT` fallback and omits the px term from the blend).

### D4. F5 handling

**Decision:** Reuse the same RFQ endpoint; select F5 legs by explicit market name (`"1st-5th Inning Spread"`, `"1st-5th Inning Total Runs"`) with a small alias table for rename resilience.

**Rationale:** Recon confirmed F5 markets ship in the same `/trade/public/api/v2/events/{id}/markets` response as FG markets, with distinct marketIds and explicit names. No heuristic period classification needed (DK's F5/FG classifier was the hairiest part of that scraper — we skip it entirely).

**Open question:** whether the RFQ endpoint accepts pure-F5 parlays or mixed FG+F5 legs. Scraper only submits pure-period combos (FG+FG or F5+F5); if either is rejected we fall back per-period per-game gracefully.

### D5. Selection ID construction

**Decision:** Read `(marketId, outcomeId, lineId, line)` directly from the market tree and pass them verbatim to the RFQ endpoint.

**Rationale:** `lineId` is a hex fingerprint that changes when ProphetX moves a line. Fetching markets and submitting the RFQ in quick succession (same session, same run) minimizes the window where a stale `lineId` could 4xx. No need to cache or predict IDs.

### D6. No cross-market canonicalization

**Decision:** Don't implement DK's "canonical market set" logic. Use exactly the (marketId, outcomeId) we extracted from the markets response.

**Rationale:** ProphetX's markets response is not split across sub-markets the way DK's is. Each market name corresponds to one `marketId`; each line variant has a unique `lineId`. There's nothing to canonicalize.

## Architecture

```
mlb pipeline (R)
    |
    v  writes mlb_parlay_lines (game_id, home, away, fg_spread, fg_total, f5_spread, f5_total, commence_time)
    |
mlb_correlated_parlay.R
    |
    |— invokes scraper_draftkings_sgp.py  (unchanged)
    |— invokes scraper_fanduel_sgp.py     (unchanged)
    |— invokes scraper_prophetx_sgp.py    [NEW]
    |       |
    |       |— GET  /trade/public/api/v1/tournaments?expand=events  → MLB events
    |       |— match to canonical game_ids
    |       |— for each game: GET /trade/public/api/v2/events/{id}/markets
    |       |— for each (game, period, combo): POST /parlay/public/api/v1/user/request
    |       |— write rows to mlb_sgp_odds (source='prophetx_direct')
    |
    v  reads mlb_sgp_odds WHERE source IN ('draftkings_direct','fanduel_direct','prophetx_direct')
    |
    |— computes dk_vig_lookup, fd_vig_lookup, px_vig_lookup (per game+period)
    |— blends: mean(model_fair_prob, dk_fair_prob, fd_fair_prob, px_fair_prob)
    |— writes mlb_parlay_opportunities with px_fair_prob column
    |
    v
mlb_dashboard.R  (reads mlb_parlay_opportunities) — renders parlays tab
```

## Files changed

### Created
- `mlb_sgp/scraper_prophetx_sgp.py` (~720 lines, tracked in git)
- `docs/superpowers/plans/2026-04-24-prophetx-sgp-scraper.md` (this doc)

### Modified
- `Answer Keys/mlb_correlated_parlay.R`:
  - Added `PROPHETX_SGP_VIG_DEFAULT <- 1.10` (later tuned from initial 1.00 prior)
  - Added `system2(...)` ProphetX scraper invocation alongside DK/FD
  - SQL WHERE clause extended to include `'prophetx_direct'`
  - Added `px_sgp` subset + `Loaded %d PX SGP` log line; no-fresh-SGP check includes px
  - Added `px_vig_lookup` block (same pattern as dk_vig_lookup / fd_vig_lookup)
  - Added PX blending block (`px_fair_prob` computed, appended to `probs_to_blend`)
  - Added `px_fair_prob = round(px_fair_prob, 3)` column to output tibble
  (Line numbers omitted intentionally — they shift as the R file evolves.)
- `mlb_sgp/recon_novig_sgp.py` (unrelated fix): landing URL changed to `https://novig.com/events`
- `.gitignore`: added `.prophetx_profile/` and `.novig_profile/` (already committed separately)

### NOT modified (intentionally)
- `mlb_dashboard.R`: reactable renders whatever columns the tibble has. `px_fair_prob` will appear as an extra column by default. If hide/reorder is wanted, change in a follow-up.
- `mlb_sgp/db.py`: schema unchanged; `bookmaker='prophetx'` and `source='prophetx_direct'` are the only new values.

## Test plan

### T1. Scraper standalone smoke test
From worktree `mlb_sgp/`:
```
venv/bin/python3 scraper_prophetx_sgp.py --verbose
```
**Pass criteria:**
- No Python exceptions
- "N ProphetX MLB events" is > 0 during active MLB season
- "M matched games" ≥ 50% of events (allows some non-MLB leagues in the tournaments endpoint to slip past the name filter; loss of 50%+ would indicate a filter bug)
- At least one game shows 4 priced combos per period (8 total for games with both FG + F5 lines)
- Per-game vig log lines appear: `[FG vig: X.XXXX]` / `[F5 vig: X.XXXX]`
- No 401/403 on the RFQ endpoint (or if present, the cookie fallback loads cookies and retries succeed)

### T2. DuckDB write verification
```
duckdb "Answer Keys/mlb.duckdb" "SELECT bookmaker, period, COUNT(*) FROM mlb_sgp_odds WHERE source='prophetx_direct' GROUP BY 1,2"
```
Expect: counts roughly equal to matched_games × 4 combos per (period, game) where lines exist.

### T3. R pipeline end-to-end
Run the MLB pipeline normally:
```
cd "Answer Keys/MLB Dashboard" && bash run.sh
```
**Pass criteria:**
- R prints `Loaded N PX SGP odds for blending`
- Dashboard parlays tab loads without error
- A new `px_fair_prob` column is visible
- `blended_prob` values changed (indicating PX joined the blend)
- No games went from valid to NA (i.e. px didn't break any existing rows)

### T4. Empirical vig calibration
After T3 passes on ≥5 separate runs (at different times of day), extract measured per-game vigs from the scraper logs. If median observed vig is materially different from 1.10 (the current default), update `PROPHETX_SGP_VIG_DEFAULT` in a follow-up commit.

### T5. Regression check
Compare blended probs on 20 sample games before and after this change:
- For games where ProphetX has prices: blended prob should shift ≤ 3% (adding a third roughly-consistent source shouldn't swing much)
- For games where ProphetX has no prices: blended prob should be identical to before
- Any game whose blended prob shifts > 10% indicates ProphetX is an outlier — investigate before merging

## Documentation updates (required before merge)

- **`mlb_sgp/README.md`** (if exists, check): add ProphetX scraper entry alongside DK/FD; note the `--verbose` flag and the RFQ ladder semantics.
- **Memory note `prophetx_sgp_scraping.md`** (new): document the endpoint, auth story (cookie fallback), selection ID structure, and vig observations. Mirror the style of `dk_sgp_scraping.md` / `fanduel_sgp_scraping.md`.
- **`MEMORY.md` index**: add pointer to `prophetx_sgp_scraping.md` under Reference section.
- **No CLAUDE.md updates needed** (no new dev-tooling conventions introduced).

## Pre-merge review checklist (per CLAUDE.md)

- [ ] Data integrity: `clear_source('prophetx_direct')` at run start prevents stale rows. Dedupe logic in `upsert_sgp_odds` uses `(game_id, combo, bookmaker, source)` key.
- [ ] Resource safety: DuckDB connections wrapped in `try/finally`.
- [ ] Edge cases: off-season (no MLB events) → "No matches found" exit; first run (no cookies) → anonymous probe; doubleheaders → UTC-hour bucket dedup.
- [ ] Dead code: none introduced.
- [ ] Log hygiene: scraper prints to stdout only; R pipeline suppresses stdout when invoking (`stdout = FALSE`).
- [ ] Security: no secrets logged; cookie SQLite read is read-only.

## Worktree lifecycle

1. ✅ Feature branch `feature/prophetx-novig-sgp-recon` created
2. ✅ Worktree `.worktrees/prophetx-novig-recon/` created
3. ✅ Recon scripts written, committed
4. ✅ ProphetX scraper + R integration written (uncommitted)
5. 🔲 Test T1–T5 (requires live MLB pipeline run)
6. 🔲 Commit scraper + R changes + Novig URL fix
7. 🔲 Memory note + README update, committed
8. 🔲 Pre-merge review (executive diff audit)
9. 🔲 User approves merge
10. 🔲 Merge to `main`
11. 🔲 Remove worktree: `git worktree remove .worktrees/prophetx-novig-recon`
12. 🔲 Delete branch: `git branch -d feature/prophetx-novig-sgp-recon`

## Open questions / risks

- **Auth model not confirmed.** Recon was inside a logged-in profile. Anonymous-first probe will resolve this on first test run.
- **Event-listing filter.** We filter tournaments by `"MLB" in name OR "Major League" in name`. If ProphetX lists MLB under a different tournament name (e.g., just "Baseball"), zero matches. Fallback: fetch whatever "Baseball" events exist and let canonical team-name matching handle it. Adjust after T1 if needed.
- **`offers[0]` ordering.** We assume the API returns the ladder tightest-first. Recon data supports this but is a sample of one. If observed ordering is inconsistent, add an explicit `max(offer.odds)` sort.
- **RFQ rate limits.** Unknown. Starting at `max_workers=4` (DK uses 6). If we see 429s we'll back off.
- **ProphetX moving lines between market fetch and RFQ.** `lineId` is a fingerprint; if it changes mid-request the POST returns an error. Low risk since fetches happen seconds apart, but not zero.

## Not in scope (explicit)

- **Novig integration.** Novig's SGP recon is blocked pending a working URL (now `novig.com/events`); scraper TBD after recon returns.
- **Mobile-app scraping** (for either venue). Mitmproxy/appium path only if web turns out to lack the data we need.
- **Internal SGP simulator** (bivariate Poisson). Tracked in memory `sgp_simulator_plan.md`; deferred.
- **Dashboard UI changes** beyond the auto-rendered `px_fair_prob` column. Any hide/reorder is a follow-up.

## Success criteria for this PR

1. T1 passes on a live session (≥1 game priced, no exceptions).
2. T3 passes (pipeline end-to-end works, parlays tab loads, `px_fair_prob` visible).
3. Memory note + README updated.
4. Pre-merge review done.
5. User explicitly approves merge.
