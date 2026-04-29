# Kalshi MLB RFQ Bot — Design Spec

**Author:** brainstorm session 2026-04-27
**Status:** Draft pending user review
**Branch / worktree:** `feat/kalshi-mlb-rfq-bot` at `.worktrees/kalshi-mlb-rfq`
**Scope:** Taker-side autonomous daemon that submits RFQs and auto-accepts +EV maker quotes for MLB game-line same-game-parlay combos on Kalshi via the cross-category multivariate-event collection.

---

## 1. Mission & Problem

Kalshi launched custom-combo same-game parlays in 2025, priced by institutional market makers via an RFQ (Request for Quote) flow. MLB game-level combos are RFQ-able today via the **`KXMVECROSSCATEGORY-R`** multivariate-event collection (verified live on 2026-04-27 — 181 MLB events eligible as legs, 4–5 maker quotes returned within seconds of submission).

Our existing infrastructure already prices these combos two ways:

- **Model fair value** — `mlb_game_samples` is a 10K-path Monte Carlo simulation per game, columns include `home_margin`, `total_final_score`, `home_margin_f5`, `total_f5`, `home_scored_first`. Joint probability of any combination of {winner, run line, total, RFI} legs is `mean(samples that hit all legs)`.
- **Sportsbook fair value** — DK / FD / ProphetX / Novig SGP scrapers (`mlb_sgp/`) write per-book correlated combo prices to `mlb_sgp_odds`. After devigging, each book contributes a fair-value estimate.

This spec turns those into a continuously-running daemon that:

1. Sweeps every game-line spread × total combo on Kalshi (Wide mode — see §4).
2. Computes a **blended fair value** from model + at-least-one-book.
3. RFQs the top-N candidates by edge against Kalshi's reference price.
4. Auto-accepts any maker quote that beats blended fair by `MIN_EV_PCT` after the real `⌈0.07 × P × (1−P)⌉` taker fee.
5. Sizes positions via **conditional Kelly** lifted from `kalshi_mm/kelly.py`, using `mlb_game_samples` to compute per-game covariance with existing positions.

The bot is fully autonomous (no human-in-the-loop on accepts), runs as a long-running process like the CBB MM bot, and ships with a complete safety scaffold (§9).

---

## 2. Scope

### 2.1 In scope

- New module **`kalshi_mlb_rfq/`** at the repo root, sibling to `kalshi_mm/`.
- Long-running daemon with the lifecycle pattern of `kalshi_mm/main.py` (SIGTERM cancel-on-shutdown, kill-switch file, phantom-RFQ cleanup at startup).
- Wide-mode combo enumerator: every spread × every total combo for every open MLB game (§4).
- 2-source fair-value gate: model always counts; require ≥1 sportsbook also pricing the combo before RFQing (§5).
- Conditional Kelly sizing via `mlb_game_samples` (§6) with accept serialization lock and per-combo cooldown.
- Real quadratic Kalshi taker-fee math (§7).
- Self-managed RFQ idempotency (recon proved `replace_existing` does not actually replace) — track our open RFQs in DuckDB, explicitly DELETE stale ones each cycle.
- **Continuous RFQ pipeline:** maintain up to `MAX_LIVE_RFQS = 80` in-flight RFQs, replenished from a priority queue ranked by edge magnitude (§4.3).
- Full safety scaffold (§9): per-accept gates (staleness, tipoff, line-move, sanity, exposure caps, cooldown, inverse-combo, positions-API health), accept serialization lock, post-accept fill reconciliation, kill switch, dry-run, fill notifications.
- Operational tooling: `python3 main.py [--dry-run]`, `bot.log` with rotation, `python3 dashboard.py` CLI status tool.
- Documentation: new `kalshi_mlb_rfq/README.md`; updates to root `README.md` and `CLAUDE.md`.

### 2.2 Out of scope (explicit non-goals; deferred for v2 or beyond)

- Maker-side RFQ responder (Phase 2 of the larger Kalshi-MLB plan).
- F5 markets (`KXMLBF5SPREAD`, `KXMLBF5TOTAL`) — not in cross-category MVE eligible-leg list.
- Player props (HR, K, hits, TB, HRR) — no fair-value source today.
- Live / in-progress games.
- Multi-account / subaccount routing.
- Any R / Shiny dashboard surfacing.
- Combos with more than 2 legs (we ship 2-leg spread × total only; expansion is enumerator-only and doesn't change the spec).
- Combos with winner or RFI legs (eligible on Kalshi but not produced by existing sportsbook scrapers, so they fail the 2-source gate).

### 2.3 Standalone process (no Wagerzon coupling)

The bot is fully independent of `mlb_parlay_opportunities` and the Wagerzon parlay placer. It writes only to its own `kalshi_mlb_rfq/kalshi_mlb_rfq.duckdb` and reads `mlb.duckdb` (samples + sgp_odds) read-only. Zero touchpoints with `mlb_dashboard.duckdb` or any Wagerzon code. No shared placement locks. Operator accepts the negligible risk that, in theory, the same canonical combo could be placed on Wagerzon (manually or via the Wagerzon placer) and Kalshi (auto by this bot) in close succession; in practice the lines and wide-mode coverage rarely collide.

---

## 3. Architecture

### 3.1 Module map

```
kalshi_mlb_rfq/                       (NEW — sibling to kalshi_mm/)
├── main.py              orchestrator daemon (long-running)
├── combo_enumerator.py  yield (game_id, leg_set) tuples per cycle
├── fair_value.py        model_fair + book_fair + blended_fair
├── ticker_map.py        game_id → Kalshi event/market tickers
├── rfq_client.py        mint combo, create RFQ, poll quotes, accept, DELETE
├── ev_calc.py           post-fee EV using ⌈0.07·P·(1−P)⌉ formula
├── kelly.py             conditional Kelly from mlb_game_samples (port from kalshi_mm)
├── risk.py              tipoff, line-move, staleness, exposure cap, kill switch
├── db.py                kalshi_mlb_rfq.duckdb
├── notify.py            fill / halt notifications
├── config.py            .env loader, defaults
├── dashboard.py         CLI status tool
└── README.md
```

### 3.2 Data flow on a single cycle

```
┌────────────────────────────────────────────────────────────────────┐
│ RFQ refresh loop @ 30s                                             │
│                                                                    │
│  1. risk.tipoff_cancel_all()    cancel RFQs for games < 5min away  │
│  2. risk.line_move_pull()       cancel RFQs whose legs moved       │
│  3. risk.exposure_cap_check()   halt if today's $$ ≥ cap           │
│  4. combo_enumerator.run()      yield candidates per open MLB game │
│  5. fair_value.score(combo)     model + book(s) → blended_fair     │
│  6. filters: 2-source gate, fair-prob bounds, sanity bounds        │
│  7. rank by |blended_fair − kalshi_reference_price|                │
│  8. take top MAX_LIVE_RFQS                                         │
│  9. diff vs db.live_rfqs():                                        │
│       - drop:   DELETE RFQs no longer in top-N                     │
│       - keep:   leave existing RFQs alone                          │
│       - add:    mint combo ticker + POST /communications/rfqs      │
│ 10. db.write(rfq_log)                                              │
└────────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────────┐
│ Quote poll loop @ 2s                                               │
│                                                                    │
│  for each live RFQ:                                                │
│    quotes = GET /communications/quotes?rfq_id=...&user=...         │
│    for each new quote:                                             │
│      fair = fair_value.score(rfq.combo)   # re-read; may have moved │
│      ev   = ev_calc.post_fee_ev(quote, fair, side)                 │
│      gates = risk.all_gates_pass(rfq, quote, fair)                 │
│      sanity = abs(quote_price - fair) <= MAX_QUOTE_DEVIATION       │
│      if ev >= MIN_EV_PCT and gates and sanity and not DRY_RUN:     │
│        kelly_size = kelly.conditional_size(rfq.combo, fair)        │
│        accept_quote(quote.id, contracts=kelly_size)                │
│        notify.fill(...)                                            │
│        db.write(fills)                                             │
└────────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────────┐
│ Pipeline refresh loop @ 600s (mirrors kalshi_mm)                   │
│                                                                    │
│  trigger MLB answer-key R pipeline:                                │
│    cd "Answer Keys" && Rscript MLB.R                               │
│  (re-populates mlb_game_samples + mlb_sgp_odds)                    │
└────────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────────┐
│ Heartbeat loop @ 60s                                               │
│  log alive line + write sessions row                               │
└────────────────────────────────────────────────────────────────────┘
```

### 3.3 Why a new module instead of extending `kalshi_mm/`

`kalshi_mm/` is a CBB 1H resting-order maker. Different sport, different mechanic (resting orders vs RFQ), different lifecycle constraints (CBB 1H is short-lived; MLB games are 3+ hours). They share auth (`kalshi_draft/auth.py`) but otherwise have no business logic in common. Two independent daemons running simultaneously — different DBs, different logs, same Kalshi account.

API rate limits are shared per-account; the 100-open-RFQ cap only affects the new bot since `kalshi_mm` doesn't use RFQs.

---

## 4. Combo enumeration (Wide mode)

### 4.1 Eligible legs per game

Per recon (2026-04-27), the cross-category MVE collection accepts these MLB single-leg market types:

| Series | Tickers per event | YES means |
|---|---|---|
| `KXMLBSPREAD-{event_suffix}` | `-{TEAM}{N}` for N ∈ {2, 3, 4, ...} | "{TEAM} wins by over (N−1).5 runs" |
| `KXMLBTOTAL-{event_suffix}` | `-{N}` for N ∈ {2, 3, ..., 12+} | "Total runs > (N−1).5" |
| `KXMLBGAME-{event_suffix}` | `-{TEAM}` | "{TEAM} wins" |
| `KXMLBRFI-{event_suffix}` | (none) | "Run scored in 1st inning" |

**v1 enumerator scope:** spread × total only. Winner and RFI legs are eligible on Kalshi but our SGP scrapers don't price them — combos with winner / RFI legs would fail the 2-source gate (§5.4) and be skipped anyway.

### 4.2 Enumeration algorithm

```
for each open MLB event:
  spread_legs = list every {KXMLBSPREAD-{event}-{TEAM}{N}, side} pair
                where Kalshi has the ticker (verified by GET /markets)
                # 2 teams × ~3 lines × 2 sides (YES, NO) = up to 12 spread legs

  total_legs  = list every {KXMLBTOTAL-{event}-{N}, side}
                where Kalshi has the ticker
                # ~10 lines × 2 sides = up to 20 total legs

  for s in spread_legs:
    for t in total_legs:
      yield combo(legs=[s, t], game_id=event)
```

Wide-mode candidate count: ~12 × 20 = **~240 combos / game**. With 15 games, that's **~3,600 candidates per cycle**. The 2-source gate (§5.4) and top-N prioritization (§4.3) bring this down to a tractable RFQ count.

### 4.3 Continuous RFQ pipeline

Kalshi caps live RFQs at 100 / account. We soft-cap at `MAX_LIVE_RFQS = 80` to leave headroom for restart races and edge cases.

Rather than rebuilding a static top-N every refresh cycle, the bot maintains a **continuous pipeline** of up to `MAX_LIVE_RFQS` in-flight RFQs, replenished from a priority queue:

1. After §5 gating, every surviving candidate enters a priority queue keyed on `|blended_fair − kalshi_reference_price|` (where `kalshi_reference_price = last_price_dollars` from `GET /markets/{combo_ticker}`; if absent or zero, fall back to `|blended_fair − 0.5|` as a contentious-combo heuristic).
2. As long as `count(live_rfqs.status='open') < MAX_LIVE_RFQS`, the bot pulls the next-highest-edge candidate from the queue and submits an RFQ.
3. RFQs naturally close as quotes are accepted, the per-combo cooldown expires unsuccessfully, the line-move check fires, or tipoff cancel sweeps. The next refresh re-enumerates and re-ranks; combos with edge larger than the lowest in-flight RFQ's submit-time edge can displace it.
4. Net effect: the top ~80 highest-edge candidates are *always* in flight asking for prices. With ~5–15s mean RFQ lifecycle, the pipeline cycles through hundreds of distinct candidates per hour. No edge sits idle waiting for a slot.

Diff logic each cycle:

- **add:** queue candidate not currently RFQ'd, slot available → mint ticker + create RFQ.
- **keep:** in-flight RFQ still in the top portion of the priority queue → leave alone.
- **drop:** in-flight RFQ now ranked below `MAX_LIVE_RFQS` worth of fresher higher-edge candidates → `DELETE /communications/rfqs/{id}`.

### 4.4 Combo ticker caching

Recon confirmed combos are **persistent** — submitting the same `selected_markets` to the lookup endpoint returns the same `market_ticker`. We cache `(canonical_leg_set_hash → combo_market_ticker)` in `kalshi_mlb_rfq.duckdb::combo_cache` to avoid re-calling the lookup endpoint each cycle.

`canonical_leg_set_hash = sha256(sorted([f"{m.market_ticker}|{m.side}" for m in selected_markets]))`

---

## 5. Pricing

### 5.1 Three sources

| Source | Computation | Coverage |
|---|---|---|
| **Model** (`mlb_game_samples`) | `model_fair = mean(samples where ALL legs hit)` | Universal — every game, every line, every combo expressible in the samples columns |
| **Books** (`mlb_sgp_odds`) | For each book ∈ {DK, FD, PX, Novig} that priced this exact combo: `book_implied = 1 / book_decimal_odds`. Devig per game: `book_fair = book_implied / per_game_vig_sum` | Partial — each book publishes a different subset of lines |
| **Blended** | `blended_fair = mean(model_fair, [book_fairs that exist])` | Used as the comparison vs Kalshi quote |

This mirrors what `mlb_correlated_parlay.R` already does for `mlb_parlay_opportunities`. We lift the math (per-book devigging, blend by simple mean) into Python; we do not call R from the daemon.

### 5.2 Joint-probability calculation from samples

For a 2-leg combo with leg A and leg B:

```python
def model_fair(samples_df, leg_a, leg_b):
    hit_a = leg_a.hit_function(samples_df)   # boolean column
    hit_b = leg_b.hit_function(samples_df)
    return float((hit_a & hit_b).mean())
```

Hit functions per leg type:

| Leg | hit function |
|---|---|
| `KXMLBSPREAD-{TEAM}{N}` YES | `home_margin >= N` (if TEAM == home) else `home_margin <= -N` |
| `KXMLBSPREAD-{TEAM}{N}` NO | inverse |
| `KXMLBTOTAL-{N}` YES | `total_final_score >= N` |
| `KXMLBTOTAL-{N}` NO | inverse |

(Winner / RFI hit functions specified for completeness in `fair_value.py` even though enumerator skips those combos in v1.)

### 5.3 Per-book devigging for `mlb_sgp_odds`

The MLB SGP scrapers each write rows like `(game_id, combo, period, bookmaker, sgp_decimal)`. For our 2-leg spread × total combo at *Wagerzon's* exact lines, the existing answer key already devigs. For Wide-mode lines (alt totals, alt spreads not on Wagerzon), we devig in-bot using the 4-combo per-game vig method:

```
For each (game_id, period, bookmaker, [spread_line × total_line]):
  vig = sum(1 / decimal_odds for all 4 sides at this (spread, total) pair)
  book_fair[combo] = (1 / decimal_odds[combo]) / vig
```

If a book has fewer than 4 sides at a given (spread, total) pair (e.g., DK pulled one side), we fall back to a per-book vig constant (defaults: DK 0.125, FD 0.18, PX 0.05, Novig 0.05) and devig as `(1 / decimal_odds) / (1 + vig_fallback)`.

### 5.4 Two-source gate

A combo proceeds to RFQ only if `count(non_null fair-value sources) ≥ 2`. Model always contributes one source. So in practice this means: `model_fair` is defined AND at least one of `{dk_fair, fd_fair, px_fair, novig_fair}` is also defined for this combo.

A combo with model only and zero books is **skipped** — too much model-error exposure on alt lines that no sharp book is even bothering to publish.

### 5.5 Fair-prob bounds

Model and combined model+book fair value must satisfy `0.05 ≤ blended_fair ≤ 0.95`. Outside that range the model is unreliable (very long shots, near-locks) and even small pricing errors blow up Kelly.

### 5.6 Sanity bounds at quote-evaluation time

When a quote arrives, refuse to accept if `|quote_implied_price − blended_fair| > MAX_QUOTE_DEVIATION` (default 0.15). Maker quotes that diverge that wildly from our blended fair are more likely a bug on our side than free money.

---

## 6. Sizing — conditional Kelly

Lifted from `kalshi_mm/kelly.py`. Same shape, same math.

### 6.1 Inputs

- `mlb_game_samples` — 10K paths per game.
- For each candidate combo: a YES-side sample-outcome vector (1 if combo hits in that path, else 0).
- Existing positions on the same game (from `kalshi_mlb_rfq.duckdb::positions`) — each represented by its own outcome vector and current contract count.

### 6.2 Conditional Kelly formula

For a new bet with mean outcome μ_new and existing position vector f_placed:

```
f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)
```

Where `Σ_nn` is the variance of the new bet's outcome, and `Σ_np` is the covariance vector between the new bet and existing positions on the same game. This matches `kalshi_mm/kelly.py::kelly_size_for_take()`.

If the covariance matrix is ill-conditioned, fall back to single-bet Kelly scaled by `1/sqrt(1 + (n-1)·avg_ρ)` — same fallback as `kalshi_mm`.

### 6.3 Fee-adjusted Kelly

For taker accepts on Kalshi, the effective execution price includes the fee. Kelly uses the fee-adjusted price:

```
yes_ask         = 1 - quote.no_bid_dollars
fee_per_contract = math.ceil(0.07 * yes_ask * (1 - yes_ask) * 100) / 100  # rounded UP to nearest cent
effective_price = yes_ask + fee_per_contract
b = (1 - effective_price) / effective_price   # Kelly's "b"
p = blended_fair                               # win probability
q = 1 - p
kelly_fraction = (b*p - q) / b
contracts      = max(0, floor(KELLY_FRACTION * kelly_fraction * BANKROLL / effective_price))
```

`KELLY_FRACTION` defaults to `0.25` (quarter-Kelly), matching the CBB MM bot.

### 6.4 No samples = no bet

If `mlb_game_samples` doesn't have rows for `game_id` (pipeline failure, off-season game, etc.), Kelly returns 0 and the bot does not RFQ that combo. There is no fixed-size fallback. Same intentional choice as CBB MM.

### 6.5 Correlation safeguards (concurrency + duplicate-side defenses)

Conditional Kelly handles correlation against positions *already in the database*. Three additional defenses cover scenarios Kelly alone doesn't:

**6.5.1 Accept serialization lock.** A `threading.Lock()` wraps the entire quote-acceptance code path (sanity gates → Kelly sizing → API call → DuckDB write). Each accept observes a fully-updated `positions` table including all prior accepts. Without the lock, two correlated quotes arriving within ~50ms could both Kelly-size against an empty positions snapshot and both fire at full size.

**6.5.2 Inverse-combo guard.** Before accepting a quote, compute the inverse leg-set hash by flipping every leg's side (YES↔NO). If `positions` has an open position keyed on the inverse, refuse the new accept and log `decision='declined_inverse_lock'` in `quote_log`. Pathological case: model error makes both A and ¬A look +EV; without the guard we'd hold a perfectly hedged, fee-bleeding pair.

**6.5.3 Per-game exposure cap.** Hard backstop above Kelly. Sum `fills.contracts × fills.price_dollars` over today's fills with the same `game_id`. If sum ≥ `MAX_GAME_EXPOSURE_PCT × BANKROLL` (default 10% of bankroll), refuse all further accepts on that game (RFQs continue — data only). Decision logged as `'declined_per_game_cap'`.

**6.5.4 Per-combo cooldown.** After accepting on a leg-set, suppress new RFQ submissions on that exact leg-set for `COMBO_COOLDOWN_SEC` (default 30s). Prevents the priority-queue replenishment from immediately re-RFQ'ing the same combo and accepting a second time.

---

## 7. EV / fee math

### 7.1 Kalshi taker fee (real formula)

```
fee_per_contract_dollars = ceil(0.07 * P * (1 - P) * 100) / 100
                         (rounded UP to the nearest cent, per Feb 2026 schedule)
```

Where P is the contract price in dollars (0.01 to 0.99). Implementation: `math.ceil(0.07 * p * (1-p) * 100) / 100`.

### 7.2 Post-fee EV per contract (YES side)

```
yes_ask  = 1 - quote.no_bid_dollars
fee      = fee_per_contract_dollars(yes_ask)
post_fee_ev_dollars  = blended_fair * (1 - yes_ask) - (1 - blended_fair) * yes_ask - fee
post_fee_ev_pct      = post_fee_ev_dollars / yes_ask
```

Accept iff `post_fee_ev_pct >= MIN_EV_PCT` (default 0.05). Mirror logic for NO side using `quote.yes_bid_dollars`.

### 7.3 Which side to accept

Each maker quote provides both `yes_bid_dollars` and `no_bid_dollars`. We compute post-fee EV for both sides (accepting YES = SELL YES at the maker's bid; accepting NO = SELL NO at the maker's bid). We accept whichever side has the higher post-fee EV, provided it clears the floor.

In our taker mode we are typically **buying** the combo (taking the YES side because our model thinks it's underpriced). So the dominant case is: `yes_ask = 1 - no_bid` is our buy price, and we accept the NO side (maker SELLS YES to us = maker accepts our YES purchase).

---

## 8. Data model

All tables in **`kalshi_mlb_rfq/kalshi_mlb_rfq.duckdb`**. Connection guard: every connection uses `try/finally con.close()`.

### 8.1 `combo_cache`

```sql
CREATE TABLE IF NOT EXISTS combo_cache (
  leg_set_hash    VARCHAR PRIMARY KEY,
  collection_ticker VARCHAR NOT NULL,
  combo_market_ticker VARCHAR NOT NULL,
  combo_event_ticker  VARCHAR NOT NULL,
  legs_json       VARCHAR NOT NULL,   -- canonicalized selected_markets array
  game_id         VARCHAR NOT NULL,
  created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_combo_cache_game ON combo_cache(game_id);
```

### 8.2 `live_rfqs`

```sql
CREATE TABLE IF NOT EXISTS live_rfqs (
  rfq_id          VARCHAR PRIMARY KEY,
  combo_market_ticker VARCHAR NOT NULL,
  leg_set_hash    VARCHAR NOT NULL,
  game_id         VARCHAR NOT NULL,
  blended_fair_at_submit DOUBLE,
  kalshi_ref_at_submit   DOUBLE,
  edge_at_submit  DOUBLE,
  status          VARCHAR NOT NULL,    -- 'open' | 'cancelled' | 'accepted'
  submitted_at    TIMESTAMP NOT NULL,
  closed_at       TIMESTAMP,
  cancellation_reason VARCHAR
);
```

### 8.3 `quote_log`

```sql
CREATE TABLE IF NOT EXISTS quote_log (
  quote_id        VARCHAR PRIMARY KEY,
  rfq_id          VARCHAR NOT NULL,
  combo_market_ticker VARCHAR,
  yes_bid_dollars DOUBLE,
  no_bid_dollars  DOUBLE,
  blended_fair_at_eval DOUBLE,
  post_fee_ev_pct DOUBLE,
  decision        VARCHAR NOT NULL,    -- 'accepted' | 'declined_ev' | 'declined_sanity'
                                       -- | 'declined_tipoff' | 'declined_line_move'
                                       -- | 'declined_dry_run' | 'declined_daily_cap'
                                       -- | 'declined_per_game_cap' | 'declined_kelly_zero'
                                       -- | 'declined_stale_predictions'
                                       -- | 'declined_inverse_lock' | 'declined_cooldown'
                                       -- | 'declined_positions_unhealthy' | 'declined_killswitch'
  reason_detail   VARCHAR,
  observed_at     TIMESTAMP NOT NULL
);
```

### 8.4 `fills`

```sql
CREATE TABLE IF NOT EXISTS fills (
  fill_id         VARCHAR PRIMARY KEY,
  quote_id        VARCHAR NOT NULL,
  rfq_id          VARCHAR NOT NULL,
  combo_market_ticker VARCHAR NOT NULL,
  game_id         VARCHAR NOT NULL,
  side            VARCHAR NOT NULL,     -- 'yes' | 'no'
  contracts       DOUBLE NOT NULL,
  price_dollars   DOUBLE NOT NULL,
  fee_dollars     DOUBLE NOT NULL,
  blended_fair_at_fill DOUBLE NOT NULL,
  expected_ev_dollars  DOUBLE NOT NULL,
  filled_at       TIMESTAMP NOT NULL,
  raw_response    VARCHAR
);
```

### 8.5 `positions`

```sql
CREATE TABLE IF NOT EXISTS positions (
  combo_market_ticker VARCHAR PRIMARY KEY,
  game_id         VARCHAR NOT NULL,
  net_contracts   DOUBLE NOT NULL,
  weighted_price  DOUBLE NOT NULL,
  legs_json       VARCHAR NOT NULL,
  updated_at      TIMESTAMP NOT NULL
);
```

### 8.6 `sessions`

```sql
CREATE TABLE IF NOT EXISTS sessions (
  session_id      VARCHAR PRIMARY KEY,
  started_at      TIMESTAMP NOT NULL,
  ended_at        TIMESTAMP,
  pid             INTEGER,
  bot_version     VARCHAR,
  dry_run         BOOLEAN NOT NULL,
  notes           VARCHAR
);
```

### 8.7 `combo_cooldown`

```sql
CREATE TABLE IF NOT EXISTS combo_cooldown (
  leg_set_hash    VARCHAR PRIMARY KEY,
  game_id         VARCHAR NOT NULL,
  cooled_until    TIMESTAMP NOT NULL,
  reason          VARCHAR                -- usually 'post_accept'
);
```

Populated on accept with `cooled_until = now() + COMBO_COOLDOWN_SEC`. Read at RFQ submission time and at quote evaluation time; both block on `cooled_until > now()`.

---

## 9. Safety scaffold

The bot is taker-only. Most defenses are **per-accept gates** rather than global halts: an RFQ sitting open with no quote action is harmless (no money committed), so we don't need to aggressively pull on every transient hiccup. Money commits only at `accept_quote`, and that's where the gauntlet runs.

### 9.1 Per-accept gate set (every accept must pass ALL)

| Gate | Trigger | Decision logged |
|---|---|---|
| Min EV after fee | `post_fee_ev_pct < MIN_EV_PCT` | `declined_ev` |
| Fair-value bounds | `blended_fair < MIN_FAIR_PROB or > MAX_FAIR_PROB` | `declined_kelly_zero` |
| Sanity bound | `|quote_implied − blended_fair| > MAX_QUOTE_DEVIATION` | `declined_sanity` |
| Prediction staleness | `now() − mlb_samples_meta.generated_at > MAX_PREDICTION_STALENESS_SEC` (1hr); also fail-safe on negative age (clock skew) | `declined_stale_predictions` |
| Tipoff window | `commence_time − now() <= TIPOFF_CANCEL_MIN` (5 min) | `declined_tipoff` |
| Line-move check | Any leg's underlying book line moved > `LINE_MOVE_THRESHOLD` since RFQ submission | `declined_line_move` |
| Per-game exposure cap | Today's `sum(contracts × price)` for this `game_id` ≥ `MAX_GAME_EXPOSURE_PCT × BANKROLL` | `declined_per_game_cap` |
| Daily exposure cap | Today's `sum(contracts × price)` across all games ≥ `DAILY_EXPOSURE_CAP_USD` | `declined_daily_cap` |
| Kill-switch off | `kalshi_mlb_rfq/.kill` exists | `declined_killswitch` |
| Inverse-combo not held | Open position on inverse leg-set hash | `declined_inverse_lock` |
| 2-source gate (re-verify) | Combo no longer has model + ≥1 book fair (e.g., book row aged out) | `declined_ev` (folded in) |
| Per-combo cooldown | `combo_cooldown.cooled_until > now()` for this leg-set | `declined_cooldown` |
| Positions API health | Last `GET /portfolio/positions` failed `POSITIONS_HEALTH_RETRIES` consecutive times | `declined_positions_unhealthy` |

### 9.2 Tipoff cancel (RFQ cleanup)

For every game with live RFQs, if `commence_time − now() <= TIPOFF_CANCEL_MIN`, the risk loop DELETEs every live RFQ for that game. Cleanup, not strictly safety — the per-accept gate above catches anything that slips through. Runs every `RISK_SWEEP_SEC` (10s).

### 9.3 Line-move pulls (RFQ cleanup)

`risk.py` snapshots reference lines at RFQ submission. The risk sweep DELETEs any RFQ whose underlying book line has moved > `LINE_MOVE_THRESHOLD`. Same accept-gate fires as a backstop if a quote arrives in the window between line move and RFQ delete.

### 9.4 Kill switch

File-based: `kalshi_mlb_rfq/.kill`. The risk loop checks every `RISK_SWEEP_SEC`. If present:
1. DELETE every live RFQ.
2. Block all `accept_quote` calls (kill-switch gate).
3. Mark `sessions` row halted; emit `notify.halt('kill_switch')`.
4. Bot stays alive in halted mode until the file is removed.

### 9.5 Accept serialization lock

A `threading.Lock()` wraps the entire quote-acceptance flow (gate evaluation → Kelly sizing → API call → DuckDB write of fill/position/cooldown). Two simultaneous quotes accept sequentially, second sees first's state. See §6.5.1.

### 9.6 Post-accept fill reconciliation

After `accept_quote` returns, the bot calls `GET /portfolio/positions` to confirm the actual contract count owned for `combo_market_ticker`. The reconciled count (not the accept response's reported count) is what gets written to `fills` and `positions`. Handles partial fills and edge cases where the accept response is out of sync with Kalshi's books.

### 9.7 Dry-run

`python3 main.py --dry-run` runs the full pipeline including RFQ creation, quote receipt, gate evaluation — but never calls `accept_quote`. Decisions write to `quote_log` with `decision='declined_dry_run'` (always emitted as the last decision). No money commits.

### 9.8 Fill notifications

Every fill and every halt emits `notify.fill(...)` / `notify.halt(...)`:

- **Always:** structured `[FILL]` / `[HALT]` line appended to `bot.log`.
- **Optional push:** if `NOTIFY_WEBHOOK_URL` is set, POST a JSON payload (Slack/Discord/Telegram-compatible incoming-webhook). **(TBD-1: operator picks channel at deploy; spec doesn't pin.)**

---

## 10. Configuration

All in `kalshi_mlb_rfq/.env` (gitignored, with `.env.example` checked in).

| Setting | Default | Description |
|---|---:|---|
| `KALSHI_API_KEY_ID` | — | Same as `kalshi_mm/.env` |
| `KALSHI_PRIVATE_KEY_PATH` | — | Same as `kalshi_mm/.env` |
| `KALSHI_USER_ID` | — | UUID for filtering `GET /communications/quotes`. Required (recon-derived). |
| `KALSHI_BASE_URL` | `https://api.elections.kalshi.com/trade-api/v2` | |
| `MVE_COLLECTION_TICKER` | `KXMVECROSSCATEGORY-R` | The cross-category collection. |
| `BANKROLL` | `1000.0` | Bankroll in dollars, for Kelly. |
| `KELLY_FRACTION` | `0.25` | Quarter Kelly. |
| `MIN_EV_PCT` | `0.05` | Min post-fee EV% to auto-accept. |
| `MAX_QUOTE_DEVIATION` | `0.15` | Sanity bound — reject quote if Kalshi price is more than this off blended fair. |
| `MIN_FAIR_PROB` | `0.05` | Lower bound on blended fair. |
| `MAX_FAIR_PROB` | `0.95` | Upper bound on blended fair. |
| `MAX_GAME_EXPOSURE_PCT` | `0.10` | Per-game cap as % of bankroll (default 10% = $100 on $1000). |
| `DAILY_EXPOSURE_CAP_USD` | `200.0` | Hard daily cap. |
| `LINE_MOVE_THRESHOLD` | `0.5` | Spread-units / runs of book-line movement that pulls RFQs. |
| `MAX_PREDICTION_STALENESS_SEC` | `3600` | Decline accepts if `mlb_samples_meta.generated_at` older than this (1hr — taker is forgiving). |
| `MAX_BOOK_STALENESS_SEC` | `60` | Ignore `mlb_sgp_odds` rows older than this when blending; if 0 books survive, 2-source gate fails. |
| `COMBO_COOLDOWN_SEC` | `30` | Suppress new RFQs / accepts on a leg-set after a fill. |
| `POSITIONS_HEALTH_RETRIES` | `2` | Consecutive `/portfolio/positions` failures before the positions-API health gate trips. |
| `RFQ_REFRESH_SEC` | `30` | Priority-queue refresh cadence. |
| `QUOTE_POLL_SEC` | `2` | Quote poll cadence. |
| `RISK_SWEEP_SEC` | `10` | Tipoff/kill-switch/exposure check cadence. |
| `PIPELINE_REFRESH_SEC` | `600` | MLB R pipeline rerun cadence. |
| `TIPOFF_CANCEL_MIN` | `5` | Minutes before first pitch to cancel RFQs and decline accepts. |
| `MAX_LIVE_RFQS` | `80` | Soft cap below Kalshi's 100/account hard cap. |
| `NOTIFY_WEBHOOK_URL` | — | Optional. POST target for fill/halt notifications. |
| `DK_VIG_FALLBACK` / `FD_VIG_FALLBACK` / `PX_VIG_FALLBACK` / `NOVIG_VIG_FALLBACK` | `0.125` / `0.18` / `0.05` / `0.05` | Fallback vig when <4 sides. |

---

## 11. Error handling

### 11.1 Kalshi API failures

| Failure | Behavior |
|---|---|
| 429 rate-limited | Back off exponentially (1s, 2s, 4s, 8s, 16s) up to 5 retries then bail this cycle. |
| 5xx | Same as 429 — back off + retry. |
| 401 / signature | Bail the cycle, log; next cycle re-derives signature from scratch. |
| 4xx other (e.g., `expired` on DELETE) | Log and continue — these are usually idempotent edge cases. |

### 11.2 Pipeline / data failures

| Failure | Behavior |
|---|---|
| `mlb_game_samples` missing rows for a game_id | Skip that game's combos; log once per cycle. |
| `mlb_sgp_odds` empty for a combo | Combo fails 2-source gate; skipped silently (expected for many alts). |
| MLB R pipeline rerun fails | Log loudly, keep RFQing on stale data. The accept-staleness gate (§9.1) automatically blocks accepts once data ages past `MAX_PREDICTION_STALENESS_SEC`. Send `notify.halt('pipeline_refresh_failed')` so the operator knows. No active halt or RFQ pull — taker doesn't need them. |

### 11.3 DuckDB write contention

Same retry-with-backoff pattern as `mlb_sgp/db.py::_connect_with_retry()` (10 attempts, exponential + jitter). Read connections never retried (DuckDB allows unlimited concurrent readers).

### 11.4 Crash recovery (startup reconciliation)

On startup, before the main loop:

1. Load `live_rfqs` rows where `status='open'`.
2. For each: `GET /communications/rfqs/{rfq_id}` to check Kalshi-side state.
   - If Kalshi says `status='open'`: keep tracking.
   - If Kalshi says `status='closed'`: mark our row closed too, log reconciliation.
3. List Kalshi-side open RFQs (`GET /communications/rfqs?status=open&creator_user_id=...`).
4. Any Kalshi-side RFQ NOT in our DB → cancel it (orphaned from a prior crashed session, same pattern as `kalshi_mm` startup phantom-cancel).

### 11.5 Idempotency & uniqueness

- Combo ticker minting: cache prevents re-calling lookup endpoint; safe to call repeatedly anyway (returns same ticker).
- RFQ creation: `live_rfqs` rows prevent duplicate submission for an active leg-set hash.
- Quote acceptance: each quote_id can be accepted at most once; the daemon writes `quote_log` row first (with `decision='accepted'`) then calls accept; on subsequent crashes the recovery path doesn't re-attempt accepts for already-logged quote_ids.

---

## 12. Version control & docs

### 12.1 Branch & worktree (already created, per CLAUDE.md)

- Branch: `feat/kalshi-mlb-rfq-bot`
- Worktree: `/Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq`

### 12.2 Files created

| File | Purpose |
|---|---|
| `kalshi_mlb_rfq/main.py` | orchestrator |
| `kalshi_mlb_rfq/combo_enumerator.py` | per-cycle candidate generation |
| `kalshi_mlb_rfq/fair_value.py` | model + book + blend |
| `kalshi_mlb_rfq/ticker_map.py` | game_id → Kalshi tickers |
| `kalshi_mlb_rfq/rfq_client.py` | RFQ + quote API wrappers |
| `kalshi_mlb_rfq/ev_calc.py` | post-fee EV |
| `kalshi_mlb_rfq/kelly.py` | conditional Kelly (port from `kalshi_mm/kelly.py`) |
| `kalshi_mlb_rfq/risk.py` | tipoff, line-move, exposure, kill switch |
| `kalshi_mlb_rfq/db.py` | DuckDB helpers + schema migrations |
| `kalshi_mlb_rfq/notify.py` | log + webhook |
| `kalshi_mlb_rfq/config.py` | .env + defaults |
| `kalshi_mlb_rfq/dashboard.py` | CLI status tool |
| `kalshi_mlb_rfq/.env.example` | template |
| `kalshi_mlb_rfq/README.md` | setup, run, monitor |
| `mlb_sgp/recon_kalshi_mlb_rfq.py` | recon script (already drafted; revised before commit) |

### 12.3 Files modified

| File | Change |
|---|---|
| `README.md` (root) | Mention `kalshi_mlb_rfq/`. |
| `CLAUDE.md` (root) | Mention the new bot in the "Project Structure" section. |
| `.gitignore` | Add `kalshi_mlb_rfq/.env`, `kalshi_mlb_rfq/*.duckdb`, `kalshi_mlb_rfq/*.log`, `kalshi_mlb_rfq/.kill`. |

No code outside `kalshi_mlb_rfq/` is modified. The bot is fully standalone (per §2.3).

### 12.4 Commit structure

Commits, in order:
1. `feat(kalshi-mlb-rfq): initial scaffolding (config, db, auth)`
2. `feat(kalshi-mlb-rfq): rfq_client + ticker_map + recon harness`
3. `feat(kalshi-mlb-rfq): combo_enumerator + fair_value (model + books + blend)`
4. `feat(kalshi-mlb-rfq): ev_calc + kelly (port from kalshi_mm)`
5. `feat(kalshi-mlb-rfq): risk (per-accept gates + RFQ-cleanup sweeps)`
6. `feat(kalshi-mlb-rfq): main orchestrator + accept lock + fill reconciliation + dashboard CLI`
7. `docs: kalshi_mlb_rfq README + root README + CLAUDE.md updates`

### 12.5 Pre-merge review checklist (per CLAUDE.md)

- Data integrity: `live_rfqs` reconciliation handles every Kalshi-side state; combo cache invalidation on schema change.
- Resource safety: every DuckDB connection in try/finally; SIGTERM cancels every live RFQ; bot.log rotation in place.
- Edge cases: off-season behavior (zero open MLB events → bot idles cleanly); first-run with empty DB; doubleheaders (hour-precision matching, same as DK/FD scrapers); time-zone boundaries in `TIPOFF_CANCEL_MIN`.
- Dead code: any prototype combo dimensions not actually enumerated in v1 should be removed before merge.
- Log/disk hygiene: bot.log rotation; combo_cache pruning for closed games (keep last 7 days).
- Security: no creds in logs (existing kalshi_mm logging conventions); `.env` not tracked.

### 12.6 Worktree lifecycle (per CLAUDE.md)

1. Worktree already created: `git worktree add -b feat/kalshi-mlb-rfq-bot .worktrees/kalshi-mlb-rfq main` — done.
2. All implementation happens inside the worktree.
3. After implementation + tests + dry-run validation: present diff to operator (`git diff main..HEAD`); operator approves merge.
4. Merge to `main` via fast-forward or merge commit (operator's choice).
5. Cleanup: `git worktree remove .worktrees/kalshi-mlb-rfq` and `git branch -d feat/kalshi-mlb-rfq-bot`.
6. Never leave the worktree behind after merge.

### 12.7 Documentation updates (per CLAUDE.md)

- `kalshi_mlb_rfq/README.md`: setup, .env vars, start/stop commands, kill-switch usage, dry-run, dashboard CLI, troubleshooting.
- Root `README.md`: one-line entry for the new bot under Project Structure.
- Root `CLAUDE.md`: one-paragraph entry under Project Structure.
- Memory: navigator/learnings entry capturing recon findings (combo lookup endpoint path, replace_existing not actually replacing, KXMLBSPREAD/TOTAL ticker structure).

---

## 13. Testing strategy

### 13.1 Unit tests

- `ev_calc`: known fee values for P ∈ {0.05, 0.20, 0.50, 0.80, 0.95}; post-fee EV calculations.
- `combo_enumerator`: given a synthetic event with N spread legs and M total legs, yields exactly N×M combos.
- `fair_value`: joint probability for a known combo against a hand-built samples dataframe.
- `kelly`: covariance computation and conditional formula match a reference implementation; falls back gracefully on ill-conditioned matrix.
- `ticker_map`: round-trip game_id → Kalshi event ticker → game_id matches; doubleheaders disambiguate by hour.

### 13.2 Integration tests (recon-style, against real Kalshi)

- Reuse `mlb_sgp/recon_kalshi_mlb_rfq.py` (already drafted) as a smoke test. Exit code 0 only if all phases succeed and zero RFQs are left open at exit.

### 13.3 Dry-run validation

Before live deploy: run `python3 main.py --dry-run` for ≥6 hours. Verify:
- RFQs created and DELETED cleanly each cycle.
- `quote_log` shows `decision='declined_dry_run'` for every quote (no accidental accepts).
- `live_rfqs` has zero rows with `status='open'` at any point >30s after cycle boundary (i.e., explicit DELETE works).
- No phantom RFQs on the account (`GET /communications/rfqs?status=open&creator_user_id=...` returns 0 ours).

### 13.4 Live deploy validation

First live session: `MAX_LIVE_RFQS=10`, `KELLY_FRACTION=0.05` (one-twentieth Kelly), `DAILY_EXPOSURE_CAP_USD=20`. Watch closely for the first hour. Tune up only after a clean session.

### 13.5 Halt conditions to monitor

In `bot.log`, every `[HALT]` line should explain: kill-switch / exposure cap / pipeline failure / SIGTERM. No silent halts.

---

## 14. Open questions / TBDs

| ID | Question | Resolution path |
|---|---|---|
| **TBD-1** | Notification channel for `NOTIFY_WEBHOOK_URL` — Slack incoming webhook vs Discord vs Telegram bot vs email-to-SMS gateway. | Operator-picked at deploy time; spec doesn't pin. README documents Slack-incoming-webhook as the recommended default. |
| **TBD-2** | RFQ TTL on Kalshi side (how long an idle RFQ stays open before auto-expiry). | Discover empirically during dry-run. Doesn't block v1 — we explicitly DELETE on our cadence. |
| **TBD-3** | Kalshi API rate limits (no `X-RateLimit-*` headers returned). | Discover empirically; back-off on 429 is in place. |
| **TBD-4** | Whether to extend `mlb_parlay_lines` / answer-key pipeline to produce winner+total and total+RFI combos so they pass the 2-source gate. | Out of v1 scope; tracked as a follow-up after v1 ships. |
| **TBD-5** | Kalshi `accept_quote` price atomicity: does it execute at the price we evaluated, or whatever the quote currently is? Determines whether we need pre-accept quote confirmation. | Verify during dry-run validation by accepting a stale quote on a small test combo and inspecting the fill price. If atomic-at-current-price, add pre-accept confirm step. |

---

## 15. Success criteria

The v1 ships when all of:

1. ≥6 hour dry-run completes with zero phantom RFQs and zero unexplained errors.
2. Live deploy with `KELLY_FRACTION=0.05` and `MAX_LIVE_RFQS=10` runs for ≥48 hours and accepts at least one fill that resolves correctly.
3. Operator can stop the bot cleanly (SIGTERM cancels every live RFQ; account is clean after exit).
4. `kalshi_mlb_rfq/README.md`, root `README.md`, and `CLAUDE.md` are updated and reflect the deployed state.
5. Pre-merge review checklist (§12.5) passes with no unresolved items.

After v1 is profitable for ≥2 weeks of live trading, work on the maker-side responder (Phase 2) begins.
