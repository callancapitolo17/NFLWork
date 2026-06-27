# Unabated-Anchored Kalshi Edge Engine — Design Spec (soccer first)

**Date:** 2026-06-24
**Status:** Draft for review
**Author:** brainstorming session (Callan + Claude)

---

## Review Pack

**What we're building** — A **sport-agnostic** engine (`unabated_edge/`) that turns your
Unabated subscription into +EV bets on Kalshi: it pulls Unabated's *sharp* odds (unblurred
via your subscription), devigs them into fair probabilities, compares to Kalshi contracts
net of fees, and takes the underpriced ones. Each sport is a small **adapter**
(`sports/<sport>.py`); we **build the core once and ship/validate soccer (World Cup) first**,
then add sports by writing adapters — no rewrite. Truth comes from a sharp-book consensus,
not a model; it takes resting order-book offers, not RFQs.

**Key decisions**

1. **Generic core + per-sport adapters; prioritize sports by edge *quality*, not count.**
   The edge differs by sport: for **soccer Unabated computes no line**, so devigging the raw
   board is *unique information* nobody else on Unabated has → durable edge; for **US sports
   Unabated hands every subscriber the same "Unabated Line"/EV** (Kalshi is even a book in the
   feed, id 105) → the edge there is execution *speed*, thin and fast-decaying. So the
   uncomputed boards (soccer) are the flagship; US sports come later and only if dry-run/CLV
   data proves a real speed edge. *Rejected:* a soccer-only silo (painful retrofit) and a
   build-all-8-sports-now platform (over-builds an unproven thesis while the Cup is live).

2. **Truth = Unabated `Sharp Book Price` + Circa, not Pinnacle.** Pinnacle is absent from
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

A **sport-agnostic core** `unabated_edge/` plus thin per-sport adapters under
`unabated_edge/sports/`, reusing `kalshi_common/` (auth, fee-aware `ev_calc`, probit devig).
Soccer ships first; adding a sport = adding one adapter file.

```
Unabated authenticated JSON feed ──► fair P(outcome)  ──┐
  (Sharp Book Price + Circa, unblurred)   via probit devig │
                                                           ├─► EV net of Kalshi fees
Kalshi order book (live asks) ─────────────────────────────┘        │
                                                                     ▼
                            net EV > threshold → take resting offer → log + PERSIST
```

It's an **order-book taker**, not an RFQ bot: Kalshi single-match markets are standard
central limit order books, so we hit resting YES/NO offers when our fair says they're
underpriced. (RFQ exists for MLB only because SGP combos have no standing book.)

### The core/adapter split

The **core** is everything sport-independent. A **sport adapter** supplies only what varies:
its Unabated league key (e.g. `lg21`), its market types and devig arity (3-way ML+draw for
soccer vs 2-way ML/runline/total for US sports), its team/name dictionary, and its mapping to
the venue's per-sport ticker structure. Adapters implement a small interface; the core never
hardcodes a sport.

| Layer | Module | Responsibility | Depends on |
|---|---|---|---|
| core | `feed.py` | Poll Unabated snapshot + cursor deltas; in-memory line book; **token mgmt + refresh** | requests, token |
| core | `pricing.py` | Devig sharp anchor → fair P (2-way & n-way; arity from adapter) | `kalshi_common.fair_value` |
| core | `ev.py` | Net EV vs Kalshi ask, incl. fees | `kalshi_common.ev_calc` |
| core | `mapping.py` | Generic Unabated↔venue reconciliation + fail-closed validation gate | adapter hooks |
| core | `storage.py` | **DuckDB: line history + research firehose + flagged edges + CLV capture** | duckdb |
| core | `execution.py` | Gates, Kelly sizing, order placement, dry-run switch (Plan 2) | `kalshi_common.auth_client` |
| core | `venues/kalshi.py` | Kalshi market discovery, order book, order placement | `kalshi_common.auth_client` |
| core | `runner.py` | Orchestration loop; iterates the registered sport adapters | all |
| adapter | `sports/soccer.py` | `lg21`, 3-way ML+draw, team dict, Kalshi WC tickers | core interfaces |
| adapter | `sports/base.py` | The `SportAdapter` interface every sport implements | — |

### Storage (the validation backbone — required from day one)

Operating the bot needs only the *current* lines in memory, but **proving** the edge needs the
*history* — and the closing line is gone the instant a match starts. So the core persists from
the first dry-run tick (mirrors `kalshi_mlb_rfq`'s 3-DB pattern):

1. **`line_snapshots`** — anchor + book odds per event over time, with a guaranteed snapshot at
   kickoff → enables **CLV** (did we beat the close?) and backtests.
2. **research firehose** — *every* candidate priced: fair, Kalshi ask, EV, gate decision,
   Kelly. Buffered + batch-flushed; never raises into the trading loop. Answers "is the edge
   real and persistent?" per sport — the data that gates each sport's go-live decision.
3. **`flagged_edges`** — the actionable +EV subset.

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

Three sibling DuckDBs (separate write locks), mirroring `kalshi_mlb_rfq`, owned by the core
and **shared across sports** — every table carries a `sport` column so one engine serves all
adapters without per-sport DB proliferation:

- **State** `unabated_edge.duckdb` — positions, orders, cooldowns. (Plan 2.)
- **Market** `unabated_edge_market.duckdb` — `line_snapshots` (Unabated anchor+book history,
  guaranteed kickoff snapshot) + Kalshi book snapshots + `flagged_edges`.
- **Research** `unabated_edge_research.duckdb` — every candidate, per-anchor fair, gate
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

Two dimensions: **depth** (dry-run → live → more market types) on the flagship sport, and
**breadth** (add sports as adapters). Breadth is gated by per-sport CLV evidence — we don't
turn on a sport until its dry-run data shows real edge.

```
Phase 0  Auth token + unblur + draw price + eventStart TZ            ← GATES EVERYTHING  [✓ unblur proven]
Phase 1  CORE + soccer adapter, --dry-run: validate mapping + persist line history + CLV
Phase 2  Flip live execution on soccer, small per-match caps
Phase 3  Soccer totals/spreads adapters; tournament futures (worse edge; revisit)
Breadth  Add the next sport as an adapter — prioritize UNCOMPUTED boards (best edge) over
         US sports (where Unabated hands everyone the EV). Each sport: dry-run → measure CLV
         → go live only if edge is real.
Decision Maker/"quoting" mode — ONLY if CLV shows fat, persistent edge.
```

This plan (Plan 1) delivers **Phase 1**: the generic core + soccer adapter running in dry-run
with full persistence. Live execution (Phase 2) and adapters/breadth are later plans.

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
