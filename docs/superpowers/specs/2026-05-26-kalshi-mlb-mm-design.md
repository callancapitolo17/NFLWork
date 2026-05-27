# Kalshi MLB RFQ **Maker** Bot — Design Spec

**Date:** 2026-05-26
**Status:** Draft for review
**Author:** brainstorming session (Claude + Callan)
**Sibling-of:** `kalshi_mlb_rfq/` (the existing *taker* bot, which stays running untouched)

---

## Review Pack

**What we're building.** A second, independent Kalshi bot that plays the *opposite* role from the existing taker: instead of creating RFQs and accepting other makers' quotes, it *listens* for RFQs other users post, prices the ones it understands, and **provides** two-sided quotes — earning the spread when someone crosses it. v1 is deliberately small and heavily instrumented: its job is to *measure* whether the maker edge survives adverse selection on real fills before we scale.

**Key decisions.**
1. **New sibling package `kalshi_mlb_mm/`, fully independent process** — not a mode of the taker. Rejected "add a `--maker` flag to the taker" because the entire taker hot-path (RFQ creation, accept, walk-diagnostics) is taker-specific; coupling them risks destabilizing a live, profitable bot.
2. **Shared pure-function code extracted into `kalshi_common/`** (`fair_value`, `ev_calc`, `auth_client`, leg-typing, `sgp_runner`). Rejected "copy the files" (divergence debt) and "maker imports from `kalshi_mlb_rfq`" (couples maker to taker package). Cost: a one-time mechanical edit to the taker's imports + a re-test/restart.
3. **Fixed-margin pricing at 5% ROI per side**, priced `bid = side_prob / 1.05 − maker_fee` (NOT a flat percentage-point shave — that gives incoherent ROI across price levels). Rejected competitive-undercut (winner's curse) and inventory-skew (needs an edge we haven't measured) for v1.
4. **REST polling behind a transport seam** (`RFQSource` / `QuoteGateway` interfaces). Rejected WebSocket for v1: a fixed-margin bot doesn't race, and the 30s confirm window makes 2s polling completely safe. The seam makes a later WS swap a one-adapter change.
5. **`BANKROLL` is the master risk dial; everything scales off it.** Starting at **$500**, daily cap stays at the chosen 75% — but in dollars that's a small, bounded validation budget. Raise one number to scale, no ratio gymnastics.

**Risks / push back here.**
- **Daily cap = 75% of bankroll** is high for an *unproven* edge; mitigated by a small ($500) bankroll, the ~5% gross spread filling rarely, and per-game/per-fill caps. The residual is a *correlated burst* (one counterparty sweeping many quotes, or many same-direction combos hitting at once) — the validation logging is built to surface this early. **Push back if you want a lower starting cap.**
- **The edge is unproven and structurally adverse-selected.** Unlike the taker (which only trades when *it* sees +EV), the maker quotes continuously and the counterparty chooses when/which side to hit — preferentially when our fair is wrong. The whole v1 is a measurement of whether 5% quoted margin nets positive *realized* edge. It might not.
- **Maker-fee assumption.** We assume RFQ quote fills are charged the **maker** fee (25% of taker). Strongly implied by the fee schedule but **must be verified on the first real fill**.
- **Soft dependency on the R pipeline** for fresh samples (shared, read-only). If samples go stale the bot declines — safe, but it means the maker isn't *fully* self-contained.

**Worth understanding** (opt-in concepts).
1. **The transport seam = dependency inversion.** The pricing brain talks to `RFQSource`/`QuoteGateway` interfaces, never raw HTTP — like separating a `fetch_data()` from an `analyse()` function in R so you can swap `read_csv` for a DB query without touching the analysis. It's what lets us pick the simple transport now without locking out the fast one later.
2. **Adverse selection = selection bias, applied to order flow.** Our *quoted* edge (5%) assumes our fair is right. But the people who cross our quote are a *biased sample* — disproportionately those who know our fair is wrong on that combo. So *realized* edge ≠ quoted edge. Same logic as non-random sampling skewing an estimate in stats.
3. **EV vs ROI (dollars vs percent).** EV in dollars (expected profit `p − c`) is absolute; ROI = `(p − c)/c` divides it by stake. They diverge in dollars but coincide in percent: "5% EV" *is* "5% ROI." We price to a constant 5% ROI so our margin is uniform across all price levels.

---

## 1. Goal & success criteria

**Goal:** stand up a working, independent maker bot that quotes the MLB spread×total combos we already model, with real (tiny) money, and produces a **validation dataset** that answers one question: *does our quoted 5% margin net a positive realized edge after adverse selection?*

**v1 is "successful" when:** the bot reliably discovers in-scope RFQs, quotes them at the correct 5%-ROI price, confirms fills inside the window with last-look protection, and records every quote/fill/settlement so we can compute realized edge. **Not** "is profitable" — that's the *output* we're measuring, not a build criterion.

**Explicitly a measurement phase**, mirroring how the taker shipped at $1 diagnostic sizing.

## 2. Scope

**In scope (v1):**
- Quote **only** 2-leg combos that are exactly one `KXMLBSPREAD-` + one `KXMLBTOTAL-` leg for the same game (the combos our `fair_value` engine already prices). Decline everything else.
- Two-sided fixed-margin quotes at 5% ROI per side.
- Real submit + confirm with last-look, on a small bankroll.

**Out of scope (future, noted so we don't design them out):**
- Other combo structures (moneyline, 3-leg, team totals) — expand later behind the same `scope` gate.
- WebSocket transport — swap behind `RFQSource` when logs justify it.
- Edge-scaled (Kelly) max-fill sizing — turn on once realized edge is calibrated.
- A shared SGP producer service (one scraper feeding both bots + dashboard) — extract when double-scrape load justifies it.

## 3. Architecture

**Process model.** Standalone daemon, separate OS process, **zero runtime dependency on the taker**. Single main loop driving timed sub-loops:

| Loop | Cadence | Job |
|---|---|---|
| Discovery + quote | ~2s | Fetch open RFQs → keep in-scope → price → submit/refresh quote |
| Confirm | ~2s | Poll our open quotes → on `accepted`, last-look → confirm or let-void → record |
| Risk sweep | ~10s | Kill-switch, tipoff pull, exposure recompute |
| SGP scrape | 60s | Own scraper cadence → own market DB |
| Samples refresh | 600s | Read `mlb_mm.duckdb::mlb_game_samples` read-only (does NOT run the R pipeline) |

**The transport seam:**
- `RFQSource.poll() -> list[RFQEvent]` — v1 `RestRFQSource` hits `GET /communications/rfqs?status=open`; a future `WsRFQSource` is a drop-in.
- `QuoteGateway` — `submit_quote` / `confirm` / `cancel` / `get_competitors`. Pricing/risk code only ever calls these two interfaces.

**Module map — `kalshi_mlb_mm/`:**

| Module | Responsibility |
|---|---|
| `main.py` | The loops + wiring |
| `rfq_source.py` | `RFQSource` interface + `RestRFQSource` (poll open RFQs, decode `mve_selected_legs`) |
| `quote_gateway.py` | `QuoteGateway` interface + `RestQuoteGateway` |
| `scope.py` | In-scope test: exactly one spread + one total leg, same game, resolvable `game_id` |
| `pricing.py` | 5%-ROI margin math → `(yes_bid, no_bid)`, grid rounding, sum<1 guard |
| `risk.py` | Size gate, exposure caps, last-look, tipoff, kill-switch |
| `db.py` | Maker schema + migrations |
| `config.py`, `notify.py` | Knobs + alerts (mirror taker patterns) |
| `sgp_runner` wrapper | Own SGP scrape cadence → `kalshi_mlb_mm_market.duckdb` |

**Shared `kalshi_common/` (extracted from the taker):** `fair_value`, `ev_calc` (+ new `maker_fee_per_contract`), `auth_client`, leg-typing helpers (`_leg_dict_to_typed`, `_parse_event_suffix`, `_MLB_CODE_TO_TEAM`), `sgp_runner`. The taker's imports are updated to point at `kalshi_common`; taker tests re-run and the taker restarted as part of the extraction step.

## 4. Data flow (one RFQ's life)

**Discovery + quote (~2s):**
1. `RFQSource.poll()` → `GET /communications/rfqs?status=open` (no creator filter → market-wide list; taker's phantom-cleanup confirmed this).
2. Dedup against `seen_rfqs`.
3. `scope.in_scope(rfq)` → `GET /markets/{market_ticker}`, read `mve_selected_legs` (verified to return the exact `{market_ticker, event_ticker, side}` shape the leg-typer consumes). Keep only 2-leg spread+total, same game; cache the verdict in `seen_rfqs`.
4. Fair value: samples + per-line book fairs → `fair_value.blend` (2-source gate). No fresh fair → skip.
5. `pricing.quote(fair)` → `(yes_bid, no_bid)` (§6), grid-rounded, `sum<1` asserted.
6. Risk pre-checks: size gate, exposure caps, kill-switch, tipoff, staleness. Any fail → skip.
7. `QuoteGateway.submit_quote(...)` → `POST /communications/quotes` → store in `live_quotes`.
8. On re-seeing an open RFQ we already quote: recompute fair; replace the quote only if price moved beyond a small hysteresis band (Kalshi auto-cancels the prior quote on the same market).

**Confirm (~2s):**
1. Poll status of each open quote.
2. On `accepted`: **last-look** — recompute fair *now*; confirm only if the filled side is still ≥ `(price + fee)`, not stale, not past tipoff, and fair hasn't drifted past `FAIR_DRIFT_TOLERANCE`. Otherwise **let it void** (no confirm). This is the core adverse-selection defense the 30s window buys us.
3. `QuoteGateway.confirm(quote_id)` → `PUT /communications/quotes/{id}/confirm`, well inside 30s.
4. Reconcile real fill via `/portfolio/positions` (reuse `get_position_contracts`); write `fills` + `positions` + exposure.

**Side note:** a two-sided quote means the *requester* picks the side, so we hold whichever side they didn't take. Because we quoted both sides at 5% ROI inside fair, we hold margin **whichever side fills** — last-look only defends against fair moving between quote and accept.

## 5. Combo decomposition (verified)

`GET /markets/{combo_ticker}` returns `mve_selected_legs`, an array of `{event_ticker, market_ticker, side}` that matches the taker's stored `legs_json` field-for-field. So an incoming RFQ → one GET → legs → existing `_leg_dict_to_typed` → existing `model_fair`/`blend`. No new parsing code; scope-detection is `len(legs)==2` with one `KXMLBSPREAD-` + one `KXMLBTOTAL-`.

## 6. Pricing & risk math

**Quote (fixed margin = 5% ROI per side):**
```
TARGET_ROI = 0.05
For each side, p = side win-prob (fair for YES, 1−fair for NO):
    raw_bid = p / (1 + TARGET_ROI)            # = p / 1.05  → (p − cost)/cost = 5%
    bid     = raw_bid − maker_fee_per_contract(raw_bid)   # one cheap iteration; fee ~0.4¢
    bid     = round_down_to_grid(bid)         # deci_cent, $0.001 step; round DOWN = conservative

yes_bid = price_for(fair)
no_bid  = price_for(1 − fair)
assert 0 < yes_bid, 0 < no_bid, yes_bid + no_bid < 1   # sum ≈ 0.94, always valid
```
- The 5% is a **quoted/expected ROI** (≡ `ev_pct` in the codebase). What we *realize* is what v1 measures.
- Gross spread ≈ 5–6%, so fills are less rare than a 10% spread → faster data, but the daily cap is more "live."

**Maker fee:** `maker_fee_per_contract = 0.25 × taker_fee` (taker = `ceil(0.07·P·(1−P)·100)/100`). Plumbed into `kalshi_common.ev_calc`. **Verify on first real fill.**

**Risk knobs (master dial = `BANKROLL`):**

| Knob | v1 value | Purpose |
|---|---|---|
| `BANKROLL` | **$500** | Master risk dial; raise to scale everything |
| `DAILY_EXPOSURE_CAP` | 75% → **$375** | Hard daily stop (user-chosen; high vs unproven edge, mitigated by small bankroll) |
| `MAX_GAME_EXPOSURE_PCT` | 10% → $50 | Per-game cap (reuse taker helper) |
| `MAX_RFQ_CONTRACTS` | 5 (≤ ~$5/fill) | Only quote RFQs with requested size ≤ this (the only "sizing" lever a maker has) |
| `MAX_OPEN_QUOTES` | 25 | Bounds outstanding risk surface + polling load; well under Kalshi's 100 |
| `TARGET_ROI` | 0.05 | The margin |
| `FAIR_DRIFT_TOLERANCE` | 0.02 | Last-look voids confirm if fair drifted >2¢ against filled side |
| `MAX_PREDICTION_STALENESS_SEC` | 600 | Decline if samples stale (reuse) |
| `TIPOFF_CANCEL_MIN` | 5 | Pull quotes near first pitch (reuse) |
| `QUOTE_HYSTERESIS` | 0.005 | Don't replace a resting quote unless fair moved >½¢ |

**On "dynamic / edge-based sizing":** deferred to v2. A maker cannot choose its fill size (quotes are full-RFQ-size); its only lever is *which RFQs to fill* (the size gate). Kelly sizing needs a *measured* edge — which is exactly what v1 produces. v1 uses a flat size gate; v2 replaces it with an edge-scaled max-fill ceiling once realized edge is calibrated. The seam keeps that a one-function change.

## 7. State / schema (`kalshi_mlb_mm.duckdb`)

| Table | Purpose |
|---|---|
| `seen_rfqs` | RFQ id, market_ticker, in-scope verdict, decoded legs, game_id, first_seen_at — dedup + scope cache |
| `live_quotes` | quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, fair_at_quote, status, submitted_at, closed_at |
| `quote_decisions` | Every decision (quoted / skipped / voided-at-last-look) + reason + diagnostics (fair drift, competitor count) — maker analog of the taker's `quote_log` |
| `fills` | fill_id, quote_id, rfq_id, combo_market_ticker, game_id, side_held, contracts, price, fee, fair_at_quote, fair_at_confirm, realized_pnl (post-settlement) — **the validation dataset** |
| `positions` | Net per (combo_market_ticker, side), reused pattern |
| `settlements` | Sweep reconciling fills against market results → populates `realized_pnl` |
| `sessions` | Run bookkeeping (reuse) |

## 8. Error handling & safety

- **Crash = safe:** a quote unconfirmed at crash time voids automatically — no surprise position.
- **Clean shutdown:** SIGTERM cancels all open quotes (free). Loops use short sleeps + frequent `_running` checks so shutdown isn't starved — directly addresses the taker's known SIGTERM-starvation pain.
- **Kill switch** (`.kill` file): stop quoting, stay alive.
- **429 back-off** on rate limits; **DuckDB retry-with-backoff** (reuse pattern).
- **Dedup**: never double-quote an RFQ; respect Kalshi's replace-on-resubmit semantics.
- **`--dry-run`:** a **pre-launch wiring smoke-test only** (decode RFQs, detect in-scope, compute prices, no exchange writes). Explicitly **not** the validation mechanism — live-small is.

## 9. Testing

- **Unit:** `pricing.quote` (5%-ROI math, fee iteration, grid round-down, sum<1), `scope.in_scope` (from real `mve_selected_legs` fixtures already captured), last-look void conditions, risk caps, `maker_fee_per_contract`.
- **Reuse** the existing `fair_value` tests (move with the module into `kalshi_common`).
- **Dry-run integration** against live open RFQs: confirms it flags in-scope combos and computes sane prices without submitting.
- **First-fill verification:** confirm the actual fee charged matches the maker-fee assumption.
- **Taker regression:** after the `kalshi_common` extraction, re-run the taker's full test suite + a dry-run before restarting it.

## 10. Open items to verify during implementation

1. **Maker fee on RFQ fills** — confirm 25%-of-taker on the first real fill.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes` for our own quotes (`status`, accepted-side, size). Recon probe before wiring the confirm loop.
3. **`get_competitors`** — confirm we can read competing quotes on an RFQ (needed only for future competitive pricing; stub in v1).

## 11. Delivery (version control)

- **Worktree:** `worktree-kalshi-mlb-mm-maker-bot` (already created).
- **Commits, roughly:** (1) extract `kalshi_common/` + repoint taker imports + re-test; (2) `kalshi_mlb_mm/` skeleton + schema + config; (3) `RestRFQSource` + `scope`; (4) `pricing` + `RestQuoteGateway`; (5) main loop + last-look confirm; (6) risk caps + safety; (7) tests + README.
- **Docs:** new `kalshi_mlb_mm/README.md`; update root `CLAUDE.md` project-structure list; note the `kalshi_common/` extraction in the taker's README.
- **Detailed task breakdown:** produced next via the writing-plans skill.
