# Kalshi World Cup Soccer Taker — Design Spec

**Date:** 2026-06-24
**Status:** Draft for review
**Author:** brainstorming session (Callan + Claude)

---

## Review Pack

**What we're building** — An autonomous Kalshi taker bot (`kalshi_soccer/`) that finds
+EV World Cup soccer bets by pricing fair probabilities from Unabated's *sharp* odds
(unblurred via your subscription) and hitting underpriced contracts on Kalshi's order
book. Same shape as your MLB RFQ taker, but truth comes from a sharp-book consensus
instead of a model, and it takes resting order-book offers instead of sending RFQs.

**Key decisions**

1. **Truth = Unabated `Sharp Book Price` + Circa, not Pinnacle.** Pinnacle is absent from
   Unabated's feed entirely; the sharp anchors that *are* there (Sharp Book Price id 7,
   Circa id 6/68) are blurred to placeholder `-110/-1.0` for anonymous users and unblur
   with your subscription. *Rejected:* scraping Pinnacle directly (geo/Cloudflare, weeks of
   work while the Cup is live).

2. **Data via Unabated's public JSON feed, authenticated — not scraping, not the $3k API.**
   The site runs on a public snapshot (`b_gameodds.json` on CloudFront) + a cursor-polled
   deltas endpoint. We poll those directly as structured JSON, sending your session token so
   the sharp anchor unblurs. *Rejected:* the $3,000/mo official API (real-time, but we don't
   need millisecond latency as a taker), and DOM scraping (fragile).

3. **Taker first; maker ("quoting") deferred to a measured Phase 3.** As a taker we only fire
   when our (slightly lagged) fair beats the live Kalshi price by enough to absorb the lag —
   self-selecting and defensive. The taker's own logs *measure* whether making is viable
   before we ever rest a quote. *Rejected:* quoting as v1 — our feed lags the sharps paying
   for the real-time API, so resting quotes invites adverse selection (your own MLB MM
   red-team flagged this as a strategy-killer).

4. **Market-type-agnostic engine, lit up incrementally.** One generic
   pull-sharp→devig→compare-to-Kalshi→net-of-fees core, with per-market-type adapters.
   Start with match moneyline, then totals, then spreads. *Rejected:* hardcoding 3-way moneyline.

5. **The draw (3-way) is the one open pricing question.** The public feed shows only 2 sides
   per match; Kalshi match markets are 3-way (home/draw/away). Phase 0 resolves whether the
   authenticated feed hands us the draw price directly (preferred) or we derive it from a
   Poisson goals model on Circa's total+spread.

**Risks / push back here**

- **Existential dependency on the authenticated token.** Without your subscription's unblur,
  there is no sharp anchor and no edge. If the token is hard to capture/refresh, the project
  stalls. This is Phase 0 and gates everything.
- **The market may already be efficient.** $300M+ has traded on Kalshi's World Cup with sharps
  actively arbing vs Pinnacle. Our edge could be thin and short-lived. The dry-run + CLV phase
  exists precisely to confirm there's real edge *before* risking capital — and to kill the
  project honestly if there isn't.
- **Feed latency.** We're seconds behind the real-time-API sharps. Acceptable for a taker on
  pre-kickoff mispricings; fatal for a maker. This bounds the strategy.
- **ToS / blocking.** Authenticated programmatic polling of Unabated's feed may bump their
  terms. Your subscription, your call. We poll politely to avoid blocks.

**Worth understanding** (opt-in)

- **Cursor-based delta polling.** Instead of re-downloading everything, you fetch a big
  snapshot once, then repeatedly ask "what changed since sequence N?" and advance N. Like
  `tail -f` on a log, or in R terms: read the full data frame once, then only `rbind` the new
  rows each tick. It's how the site stays live without re-pulling 2 MB every second.
- **Devigging a 3-way market.** A book's three prices imply probabilities that sum to >100%
  (the vig). Removing the vig fairly across 3 outcomes is the n-way case your probit method
  already handles via root-finding — the same `uniroot`-style solve you use for MLB.

---

## 1. Architecture

A new standalone process **`kalshi_soccer/`**, mirroring `kalshi_mlb_rfq/`, reusing
`kalshi_common/` (auth, fee-aware `ev_calc`, probit devig, leg typing helpers).

```
Unabated authenticated JSON feed ──► fair P(outcome)  ──┐
  (Sharp Book Price + Circa, unblurred)   via probit devig │
                                                           ├─► EV net of Kalshi fees
Kalshi order book (live asks) ─────────────────────────────┘        │
                                                                     ▼
                                          net EV > threshold → take resting offer → log
```

It's an **order-book taker**, not an RFQ bot: Kalshi single-match markets are standard
central limit order books, so we hit resting YES/NO offers when our fair says they're
underpriced. (RFQ exists for MLB only because SGP combos have no standing book.)

### Components (each independently testable)

| Module | Responsibility | Depends on |
|---|---|---|
| `feed.py` | Poll Unabated snapshot + cursor deltas; maintain in-memory line book; auth token mgmt | requests, token |
| `mapping.py` | Reconcile Unabated events/lines ↔ Kalshi contracts (team dict, line match, validation gate) | feed, kalshi market list |
| `pricing.py` | Devig sharp anchor → fair P per outcome (3-way via draw price or Poisson fallback) | `kalshi_common.fair_value` |
| `ev.py` | Net EV vs Kalshi ask, incl. fees | `kalshi_common.ev_calc` |
| `execution.py` | Gates, Kelly sizing, order placement, dry-run switch | `kalshi_common.auth_client` |
| `state.py` | DuckDB state/market/research writes | duckdb |
| `main.py` | Orchestration loop | all |

---

## 2. Data acquisition (confirmed by live spike, 2026-06-24)

**Feed endpoints (public JSON, not scraping):**

```
SNAPSHOT  GET content.unabated.com/markets/game-odds/b_gameodds.json   (~2.2 MB gz, CloudFront, CORS *)
DELTAS    GET api-k.unabated.com/api/markets/changes/query/<cursor>    (cursor-polled)
BOOTSTRAP GET api-k.unabated.com/api/markets/changes/query?full_refresh_ISO=<iso>
```

**Snapshot schema (confirmed):**

- Top-level: `gameOddsEvents`, `marketSources`, `teams`, `people`, `propsPeopleEvents`,
  `propsPeopleEventsFutures`, `startDateTimeISO`, `endDateTimeISO`.
- `gameOddsEvents` keyed `lg<leagueId>:pt<periodTypeId>:pregame`. **World Cup = `lg21`**
  (52 matches confirmed: Argentina, France, Portugal, Norway… on 6/22–6/23).
- Each event: `eventId`, `eventStart`, `eventTeams` {`"0"`:{id}, `"1"`:{id}},
  `gameOddsMarketSourcesLines`.
- Lines keyed `si<sideId>:ms<marketSourceId>:an<alternateNumber>` → `bt1` (moneyline) /
  `bt2` (spread) / `bt3` (total) → `{ americanPrice, points, isBlurred, marketSourceId,
  sequenceNumber, modifiedOn, … }`.
- `marketSources`: id→name. Key ids: **Sharp Book Price 7, Circa 6, Circa Sports 68**,
  FanDuel 2, DraftKings 1, BetOnline 9, Bookmaker 8, Novig 89, Prophet Exchange 66,
  Buckeye 59, **Kalshi 105**, Polymarket 107. **Pinnacle: absent.**
- `teams`: id→{name, abbreviation, leagueId, …}.

**Blur mechanism (confirmed):** Sharp anchors (7, 6, 68) are **100% blurred** in the
anonymous feed with placeholder `-110 / -1.0` values — useless. Soft books (FanDuel, DK,
etc.) are unblurred and real. **Your subscription unblurs the anchors** — this is the entire
value of the subscription and the precondition for the edge.

**Auth (PROVEN 2026-06-27):** Auth0 flow (`/api/auth/login`). The unblur key is the
**`unabated_at_prod` access-token cookie** set on `.unabated.com`, sent automatically to
`api-k.unabated.com` (credentialed CORS). It is a JWT with `permissions: ["role:premium"]`,
**valid 30 days** (`iat`→`exp`), refreshed from the longer-lived `unabated` session cookie.

**Unblur — PROVEN.** Replaying an authenticated `changes/query` returned **real** sharp-anchor
values vs. the anonymous placeholder:

| Book | Anonymous | Authenticated |
|---|---|---|
| Sharp Book Price (7) | `-110 / -1.0` | real, varied (`-353,-251,-127…`) ✓ |
| Circa (6) | `-110 / -1.0` | real (`-250,150,215…`, ±1.5) ✓ |
| Circa Sports (68) | `-110 / -1.0` | real ✓ |

**Authenticated feed mechanism (confirmed):** `GET api-k.unabated.com/api/markets/changes/query/<cursor>`
returns `{ latestTimestamp, resultCode, results:[ { marketLineChanges:[ { gameOdds:{ gameOddsEvents:
{ "lgNN:ptN:pregame|live": [ events ] } } } ] } ] }`. Cursor = `latestTimestamp`. `eventStart`
is **TZ-aware UTC** (`+00:00`). Bet types: `bt1`=moneyline (price only), `bt2`=spread/Asian
handicap, `bt3`=total, `bt4`=present on soccer (encoding TBD — see below).

### Phase 0 status

1. **Unblur mechanism — ✓ PROVEN** (token = `unabated_at_prod` cookie; real anchor values).
2. **Token lifecycle — ✓** 30-day `role:premium` JWT; refresh via `unabated` session cookie.
   Store in gitignored `.env`; build a refresh step (don't hardcode a static token).
3. **`eventStart` timezone — ✓** TZ-aware UTC in the authenticated feed.
4. **Draw price — OPEN** (the one remaining item). World Cup uses only si0/si1 (no draw side);
   the 3-way "Moneyline With Draw" the UI shows must live in a bet type (`bt4`?) — confirm via
   one capture against a *live* WC 3-way market. Fallback remains the Poisson derivation.

**Token capture method (implemented in spike):** decrypt Chrome cookies via Keychain *or*
copy one `changes/query` request as cURL from the logged-in browser (used here). The bot will
read the token from `.env` and refresh it.

---

## 3. Mapping layer (the hardest engineering)

Per `(match, market_type)`, reconcile Unabated ↔ Kalshi:

- **Match identity:** normalized team pair + kickoff time. World Cup is finite (48 teams) →
  a canonical country dictionary is enumerable (e.g. "Korea Republic" = "South Korea"). Build
  it once from the `teams` map ↔ Kalshi team naming.
- **Line identity (totals/spreads):** Kalshi line must equal Unabated's exactly; use alt
  lines (`an>0`) to find the match, else **skip** (strict-match discipline, like Novig).
- **Kalshi discovery:** enumerate Kalshi World Cup match series/markets via `kalshi_common`
  (ticker structure to confirm in Phase 1).
- **Validation gate (fail-closed):** kickoff within tolerance, all required sides present,
  devigged probabilities sane → else the candidate never becomes an order. A mismapped
  team/line is the highest-severity failure (confident fair on the wrong contract).

---

## 4. Pricing & EV engine

- **Anchor:** `Sharp Book Price` (id 7) primary; `Circa` (6/68) cross-check / blend.
- **Devig:** `kalshi_common` probit — closed-form 2-way (totals/spreads), n-way uniroot
  (3-way moneyline).
- **Draw handling:** if Phase 0 yields a real draw price → devig 3-way directly. Else →
  derive P(home/draw/away) from a Poisson goals model fit to Circa's total + spread
  (standard soccer approach; small, well-understood).
- **EV:** `net_EV = fair_P − (kalshi_ask/100) − kalshi_fee` via `kalshi_common.ev_calc`.
  Flag only if `net_EV > threshold`.
- **Threshold buffer:** conservative initial buffer for feed lag + fair-value uncertainty.
  Note the MLB finding that edge is overstated above ~7%; let CLV calibrate the threshold.

---

## 5. Execution, gates, sizing

- **Taker:** walk the Kalshi book, take resting offers up to the size where net EV stays
  positive.
- **Gates (reuse `kalshi_mlb_rfq` scaffold):**
  - Feed staleness (anchor tick age) — hard stop if too old.
  - Pre-kickoff cutoff.
  - Anchor line-move recompute (re-price if sharp moved since fair computed).
  - Exposure caps: per-match, per-team, global bankroll.
  - Fill-ratio halt (detect being picked off).
  - **`--dry-run`** flag: price + log orders without sending (bring-up safety).
- **Sizing:** fractional Kelly on net edge, capped. The 3 contracts in one match are mutually
  exclusive → bounded by a single per-match exposure cap, not sized independently.

---

## 6. Data model

Three sibling DuckDBs (separate write locks), mirroring `kalshi_mlb_rfq`:

- **State** `kalshi_soccer.duckdb` — positions, orders, cooldowns.
- **Market** `kalshi_soccer_market.duckdb` — Unabated line snapshots + Kalshi book snapshots.
- **Research** `kalshi_soccer_research.duckdb` — every candidate, per-anchor fair, gate
  decision, Kelly, fill. Buffered/batch-flushed; can never raise into the trading loop.

All timestamps `TIMESTAMPTZ` UTC (repo rule).

---

## 7. Observability & CLV

- **CLV monitor (the thing the MLB stack lacks):** log the anchor's *closing* line per match;
  measure whether we beat the close. This is how we *prove* edge before scaling — and the
  data that gates the Phase 3 maker decision.
- Operational logging via Python `logging` + rotating `bot.log` (no `print`).

---

## 8. Phasing

```
Phase 0  Auth token capture + confirm unblur + draw price + eventStart TZ   ← GATES EVERYTHING
Phase 1  Match moneyline, --dry-run: validate mapping + measure edge/CLV on live matches
Phase 2  Flip live execution, small per-match caps
Phase 3  Totals adapter → spreads adapter
Decision Maker/"quoting" mode — ONLY if Phase 1–2 CLV shows fat, persistent edge
Later    Tournament futures (worse edge; revisit)
```

---

## 9. Version control & worktree

- Worktree: `.claude/worktrees/kalshi-soccer-wc`, branch `worktree-kalshi-soccer-wc`
  (already created for this spec).
- Files created (implementation): `kalshi_soccer/{feed,mapping,pricing,ev,execution,state,
  main}.py`, `kalshi_soccer/README.md`, `kalshi_soccer/.env.example`, tests under
  `kalshi_soccer/tests/`.
- Files modified: `CLAUDE.md` (project-structure bullet), root `README.md` (tool list).
- Commits structured per phase. Pre-merge executive review of full diff. **No merge to `main`
  without explicit approval.** Clean up worktree + branch after merge.

## 10. Documentation

- New `kalshi_soccer/README.md`: setup, token capture, run/dry-run usage, troubleshooting.
- Update root `CLAUDE.md` project-structure list and `kalshi_common/` note if shared code grows.
- Docs land in the same merge as the code.

---

## Decisions (resolved with user)

1. **Token capture path — from the user's normal logged-in browser.** Phase 0 grabs the
   Unabated session token by manual copy from the browser tab the user is already logged into
   (no password). Note: the chrome-devtools automation profile is separate from the user's
   normal Chrome, so the token is copied out, not read from the automation browser.
2. **Bankroll = $1,000.** Drives Phase 2 exposure caps and Kelly sizing (fractional Kelly,
   per-match and global caps scaled to $1k).

## Still open (resolved at Phase 0 by what the feed yields)

- **Draw fallback** — if Phase 0 finds no clean draw price in the authenticated feed, derive
  P(home/draw/away) from a Poisson goals model on Circa's total+spread (spec default). If the
  feed *does* expose the draw, devig 3-way directly and skip the model.
