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
- Conditional Kelly sizing via `mlb_game_samples` (§6).
- Real quadratic Kalshi taker-fee math (§7).
- Self-managed RFQ idempotency (recon proved `replace_existing` does not actually replace) — track our open RFQs in DuckDB, explicitly DELETE stale ones each cycle.
- Top-N prioritization with soft cap below Kalshi's 100 RFQ/account hard limit.
- Full safety scaffold (§9): tipoff cancel, line-move detection, sanity bounds, daily exposure cap, kill switch, dry-run, placement lock spanning Wagerzon + Kalshi, fill notifications.
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

### 2.3 Wagerzon coordination

The bot is independent of `mlb_parlay_opportunities` and the Wagerzon parlay placer. To prevent double-placement of the same canonical combo across Wagerzon and Kalshi, both placers consult a shared `placement_locks` table in `mlb_dashboard.duckdb` keyed on a `combo_key` (game_id + sorted-leg-descriptors hash). See §9.6.

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

### 4.3 Top-N prioritization

Kalshi caps live RFQs at 100 / account. We soft-cap at `MAX_LIVE_RFQS = 60` to leave headroom.

After computing `blended_fair` for every surviving candidate (§5):

1. Rank candidates by `|blended_fair − kalshi_reference_price|` where `kalshi_reference_price` is `last_price_dollars` from `GET /markets/{combo_ticker}`. If `last_price_dollars` is `0.0000` or absent, fall back to `|blended_fair − 0.5|` (heuristic: contentious combos are interesting).
2. Take top `MAX_LIVE_RFQS` by rank.
3. Diff against `db.live_rfqs`:
   - **add:** combos newly in top-N → mint ticker + create RFQ.
   - **keep:** combos still in top-N with a live RFQ → leave alone.
   - **drop:** combos no longer in top-N → `DELETE /communications/rfqs/{id}`.

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
                                       -- | 'declined_dry_run' | 'declined_exposure_cap'
                                       -- | 'declined_lock' | 'declined_kelly_zero'
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

### 8.7 `placement_locks` (in `mlb_dashboard.duckdb`, NOT the bot's DB)

```sql
-- in mlb_dashboard.duckdb (shared with Wagerzon parlay placer)
CREATE TABLE IF NOT EXISTS placement_locks (
  combo_key       VARCHAR PRIMARY KEY,   -- canonical hash of (game_id + sorted leg descriptors)
  venue           VARCHAR NOT NULL,      -- 'wagerzon' | 'kalshi'
  status          VARCHAR NOT NULL,      -- 'locked' | 'placed' | 'released'
  placed_by_pid   INTEGER NOT NULL,
  acquired_at     TIMESTAMP NOT NULL,
  released_at     TIMESTAMP
);
```

`combo_key` is computed identically by both placers. Acquire = `INSERT ... ON CONFLICT DO NOTHING` then check `last_insert`. Release = `UPDATE status='released'` after fill or rejection. See §9.6.

---

## 9. Safety scaffold

All non-negotiable for v1.

### 9.1 Tipoff cancel

For every game in `live_rfqs`: if `commence_time - now() <= TIPOFF_CANCEL_MIN` (default 5 min), DELETE every live RFQ for that game and refuse to create new RFQs for that game until next pipeline refresh (when the game will have rolled out of the open-events list).

Runs every 10s in `risk.py::tipoff_sweep()` independent of the main 30s RFQ refresh.

### 9.2 Line-move detection

Mirrors `kalshi_mm/risk.py::line_move_check`. Maintains a `reference_lines` snapshot at RFQ submission time. If any leg's underlying book line (from DK/FD/PX/Novig single-leg odds — separate scrapers, already running) moves by more than `LINE_MOVE_THRESHOLD` (default 0.5 spread units / 0.5 total runs), DELETE the RFQ and refuse to re-submit until next pipeline refresh.

### 9.3 Sanity bounds (already covered §5.6)

Refuse to accept if `|quote_implied_price − blended_fair| > MAX_QUOTE_DEVIATION`.

### 9.4 Daily exposure cap

Sum `fills.contracts × fills.price_dollars` for fills with `filled_at >= today_start_et`. If sum >= `DAILY_EXPOSURE_CAP_USD` (default $200), halt — refuse new RFQ creation and refuse all quote acceptances. The bot stays alive in halted mode; halt clears at next ET midnight.

### 9.5 Kill switch

File-based: `kalshi_mlb_rfq/.kill`. The risk loop checks for it every 10s. If present:
1. DELETE every live RFQ.
2. Mark sessions row halted.
3. Stop accepting quotes.
4. Stay alive until file is removed.

### 9.6 Placement lock (Wagerzon ↔ Kalshi)

Both the Wagerzon parlay placer and this bot acquire a row in `mlb_dashboard.duckdb::placement_locks` keyed on `combo_key` BEFORE placing on their respective venue. If acquisition fails (row already exists with another venue), the placer skips that combo.

`combo_key = sha256(game_id + ":" + sorted(leg_descriptors))` where `leg_descriptors` are canonical strings like `"spread:home:-1.5"`, `"total:over:8.5"`. Both placers compute it identically.

The Wagerzon placer needs a small change: add `combo_key` computation and `placement_locks` interaction. Out of this spec's primary scope but in scope as a coordinated change in the same merge. Tracked in §12.

### 9.7 Dry-run

`python3 main.py --dry-run` runs the full pipeline including RFQ creation, quote receipt, EV evaluation — but never calls `accept_quote`. Decisions are written to `quote_log` with `decision='declined_dry_run'`. Use to validate behavior without trading real money.

### 9.8 Fill notifications

Every fill and every halt sends a notification on `notify.fill(...)` / `notify.halt(...)`. v1 supports two channels, configured via env:

- **Always:** append to `bot.log` with a structured `[FILL]` / `[HALT]` line.
- **Optional push:** if `NOTIFY_WEBHOOK_URL` is set, POST a JSON payload to that URL. The webhook URL is meant to be a Slack/Discord/Telegram-compatible incoming-webhook endpoint chosen by the operator at deploy time. **(TBD: operator picks the channel; spec doesn't pin it.)**

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
| `DAILY_EXPOSURE_CAP_USD` | `200.0` | Hard daily cap. |
| `LINE_MOVE_THRESHOLD` | `0.5` | Spread-units / runs of book-line movement that pulls RFQs. |
| `RFQ_REFRESH_SEC` | `30` | RFQ refresh cadence. |
| `QUOTE_POLL_SEC` | `2` | Quote poll cadence. |
| `RISK_SWEEP_SEC` | `10` | Tipoff/kill-switch/exposure check cadence. |
| `PIPELINE_REFRESH_SEC` | `600` | MLB R pipeline rerun cadence. |
| `TIPOFF_CANCEL_MIN` | `5` | Minutes before first pitch to cancel. |
| `MAX_LIVE_RFQS` | `60` | Soft cap below Kalshi's 100/account hard cap. |
| `MIN_FAIR_PROB` | `0.05` | Lower bound on blended fair. |
| `MAX_FAIR_PROB` | `0.95` | Upper bound on blended fair. |
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
| MLB R pipeline rerun fails | Log loudly, keep running on stale data. Operator notification via `notify.halt('pipeline_refresh_failed')`. |

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

### 12.3 Files modified (out of primary scope but coordinated in same merge)

| File | Change |
|---|---|
| `bet_placer/parlay_placer.py` (Wagerzon) | Add `combo_key` computation + `placement_locks` acquire/release. |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | Same — coordinated `placement_locks` interaction. |
| `mlb_dashboard.duckdb::placement_locks` (schema) | New table; bot's `db.py` runs migration on startup if not present. |
| `README.md` (root) | Mention `kalshi_mlb_rfq/`. |
| `CLAUDE.md` (root) | Mention the new bot in the "Project Structure" section. |
| `.gitignore` | Add `kalshi_mlb_rfq/.env`, `kalshi_mlb_rfq/*.duckdb`, `kalshi_mlb_rfq/*.log`, `kalshi_mlb_rfq/.kill`. |

### 12.4 Commit structure

Commits, in order:
1. `feat(kalshi-mlb-rfq): initial scaffolding (config, db, auth)`
2. `feat(kalshi-mlb-rfq): rfq_client + ticker_map + recon harness`
3. `feat(kalshi-mlb-rfq): combo_enumerator + fair_value (model + books + blend)`
4. `feat(kalshi-mlb-rfq): ev_calc + kelly (port from kalshi_mm)`
5. `feat(kalshi-mlb-rfq): risk (tipoff, line-move, exposure, kill switch)`
6. `feat(kalshi-mlb-rfq): main orchestrator + dashboard CLI`
7. `feat(placement-locks): shared lock table for wagerzon ↔ kalshi-mlb-rfq`
8. `docs: kalshi_mlb_rfq README + root README + CLAUDE.md updates`

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

---

## 15. Success criteria

The v1 ships when all of:

1. ≥6 hour dry-run completes with zero phantom RFQs and zero unexplained errors.
2. Live deploy with `KELLY_FRACTION=0.05` and `MAX_LIVE_RFQS=10` runs for ≥48 hours and accepts at least one fill that resolves correctly.
3. Operator can stop the bot cleanly (SIGTERM cancels every live RFQ; account is clean after exit).
4. `kalshi_mlb_rfq/README.md`, root `README.md`, and `CLAUDE.md` are updated and reflect the deployed state.
5. Pre-merge review checklist (§12.5) passes with no unresolved items.

After v1 is profitable for ≥2 weeks of live trading, work on the maker-side responder (Phase 2) begins.
