# World Cup Totals Market Maker — Design Spec

Date: 2026-07-10
Branch: `worktree-wc-maker-depth-capture` (builds on the depth-capture commit `8299cce`)
Status: awaiting user approval

---

## Review Pack

**What we're building.** A market maker for Kalshi's World Cup regulation-time
totals (`KXWCTOTAL`) that rests limit orders around the fair value the
`unabated_edge` engine already computes from the devigged Unabated sharp
anchor. It earns the bid/ask spread from retail flow instead of needing the
market to be wrong (the taker dry-run proved it isn't: 1-cent spreads sitting
on the anchor, every taker EV negative). It goes live immediately at
ledger-bounded size and learns from real fills — measurement-first, same
philosophy as the MLB maker v1.

**Key decisions.**

1. **Live-first, no shadow phase** (user decision, revised from shadow-first).
   Real fills measure queue position and adverse selection directly; inferred
   shadow fills only bound them. Rejected: shadow through the semis — burns
   half the remaining tournament for weaker data. A ShadowGateway still ships
   for testing future leagues.
2. **Maker lives inside `unabated_edge` as a module, same process** (Option A;
   rejected: standalone daemon à la `kalshi_mlb_mm`). Fair value, book
   snapshot, and quote decision come from the same 5s tick — no duplicated
   polling, and every quote is attributable to the exact market state that
   produced it. Extracting to a daemon later is mechanical.
3. **Risk = exact worst-case ledger on the goal grid, not per-market caps.**
   All of a match's rungs settle on one random variable (total goals), so
   portfolio P&L is evaluated exactly over g=0..10 and capped per match.
   Offsets (Under 2.5 + Over 1.5) and middles price correctly for free.
   Rejected: per-outcome caps — the known correlated-inventory defect in the
   taker path.
4. **Inventory rides through settlement** (user decision). Flattening at
   kickoff pays back the same 1-2c spread the maker just earned. The match
   cap IS the agreed amount at risk on any single game.
5. **All anchor-quoted rungs, two-sided, dynamically pruned.** Every rung the
   anchor ladder prices gets two-sided quotes; the ledger shrinks/pulls the
   side that would deepen the worst case. Alt-line rungs quote wider and
   smaller (provenance gate). Rejected: main-rung-only — joins the deepest
   1c queue, the least attractive spot for a maker.

**Risks / push back here.**

- **Aggressive risk caps (user-set):** 40% of bankroll worst-case per match,
  75% global, 30% per resting order. Two bad settlements = −80% of bankroll.
  This is a deliberate pay-for-information posture; all are single env vars
  to walk back after the first match day.
- **3-minute pull window (user-set, revised from 15):** we quote straight
  through the lineup-news window (~60-75 min pre-kick). Defenses are
  requote-on-anchor-move within one 5s tick and the fill-burst tripwire.
  If the first match day shows lineup-window pick-offs, widen it.
- **Alt-rung trust is unvalidated until the capture data says otherwise.**
  Alt rungs quote 1.5x wider and 0.5x smaller, but whether they should quote
  at all is an empirical question the Jul 11 capture answers.
- **Timeline:** live for the Jul 11 semis requires build + tests + your merge
  approval within ~1 day. Scoped v1 beats perfect design after the final
  (Jul 19), but the gate is real: no merge, no order without your explicit go.

**Worth understanding** (opt-in, anchored to R):

1. *The goal-grid ledger* is `sapply(0:10, function(g) sum(pnl_of_each_fill(g)))`
   then `min()` — risk as an exact function over a small outcome space,
   instead of Greeks/heuristics. Possible only because every instrument in
   the match settles on one integer.
2. *Interface + swappable implementations* (`QuoteGateway` → Shadow or Live):
   like writing one R function that takes a `write_fn` argument — callers
   don't change when the implementation does. This is how "live" stays a
   one-line, explicitly-acked config change instead of a code path fork.

---

## 1. Context

- Engine: `unabated_edge/` (merged to main `1408ee1`) polls Unabated's v2
  per-league file every 5s, deviggs the sharp anchor's total ladder per match
  (`sports/soccer.py::_anchor_ladder` → `{line: {p_over, book, alt,
  overround}}`), pairs events to Kalshi, and (as of `8299cce`, this branch)
  captures full Kalshi orderbook depth (`book_snapshots`) and the executed
  trades tape (`kalshi_trades`) per rung per tick, pre-kickoff only.
- Taker findings (measured, 2026-07-09): Kalshi crowd quotes 1c spreads with
  $10k+ top-of-book depth on main rungs, centered on the anchor; ~216
  trades/hour of retail flow clustered on the 2.5 rung. No taker edge exists;
  maker edge = spread capture from that flow.
- MLB maker (`kalshi_mlb_mm/`, live) supplies the patterns: fixed ROI margin,
  measurement-first v1, QuoteGateway/RFQSource interfaces, sibling-DB
  discipline, fill ledgers. This maker rests orders on standard order books
  (not RFQs).

## 2. Architecture

New module `unabated_edge/maker/`, invoked from `runner.run_tick` after
pricing, consuming what the tick already fetched:

```
run_tick (every 5s, per paired pre-kickoff match):
  anchor ladder (devigged)  ─┐
  kalshi books (get_book)   ─┼─→ maker.engine.on_tick(...)
  trades tape (30s)         ─┘        │
                                      ├─ ledger.py: worst-case over goal grid
                                      ├─ desired quotes = fair ± margin, capped
                                      ├─ diff vs resting orders
                                      ▼
                              QuoteGateway
                              ├─ LiveGateway   → Kalshi REST orders (v1 default)
                              └─ ShadowGateway → DB only (testing/new leagues)
```

Files:

| file | job |
|---|---|
| `maker/engine.py` | per-match quote decisions: margins, never-cross guard, alt gating, requote/pull triggers, cap enforcement |
| `maker/ledger.py` | pure math: fills → pnl(g) over g=0..10 → worst case; hypothetical-fill check |
| `maker/gateway.py` | `QuoteGateway` ABC; `LiveGateway` (place/cancel via `kalshi_common.auth_client`); `ShadowGateway` (record only) |
| `maker/state.py` | resting-order book, fills/positions reconciliation, maker-DB writes |

Sport-agnostic: the engine consumes any adapter's ladder; nothing in
`maker/` is soccer-specific. Adding a league to the maker = the adapter that
already exists for the taker.

## 3. Quoting logic

Per rung with devigged fair `f` (prob of Over), all prices in cents:

```
margin m   = max(maker_fee_per_contract(p) + PICKOFF_BUFFER_CENTS,
                 ROI_MARGIN * p)                    # p = relevant price
yes_bid    = floor_to_cent(f − m)                   # our bid on Over
no_bid     = floor_to_cent((1 − f) − m)             # our bid on Under
```

Guards, in order:

1. **Never cross:** each bid is capped 1 tick below the opposing best ask
   (resting at/above it would execute as a taker). If the cap pushes the bid
   below `f − MAX_MARGIN_CENTS` (default 5c), skip the rung — the crowd is
   tighter than we're willing to be. This automatically keeps us out of deep
   1c main-rung queues and concentrates quoting where spreads are wide.
2. **Alt rungs second-class:** `alt=True` rungs use `m × ALT_MARGIN_MULT`
   and `size × ALT_SIZE_MULT`, and only quote when the rung's `overround` is
   in `[ALT_OVERROUND_MIN=1.01, ALT_OVERROUND_MAX=1.15]` — a real two-way
   market's vig band; outside it the "quote" is likely a derived or broken
   number. Initial bounds; recalibrate from the Jul 11 capture.
3. **Sanity:** bids in [1c, 99c]; the two sides of one rung must not overlap.

Sizing: order size = `MAX_QUOTE_PCT × BANKROLL` of contract cost, then shrunk
so a full hypothetical fill keeps the match worst case within its cap and the
global cap (see §4). In practice the ledger, not the quote pct, binds.

**Requote / pull triggers** (checked per tick, per match):

| trigger | action |
|---|---|
| any rung's fair moved ≥ 1c since our quote | cancel/replace that rung this tick |
| feed `modifiedOn` age > `MAX_STALENESS_SEC` (20s) or fetch failed | pull match |
| within `QUOTE_PULL_MIN` (3 min) of kickoff | pull match; inventory rides |
| kill file / tripwire (§5) | pull everything |

## 4. Risk layer

**Ledger** (pure function): fills = `(line, side, contracts, price)`;
`pnl(g) = Σ` YES-on-Over-L: `n × (1{g>L} − p)`; NO-on-Over-L:
`n × (1{g≤L} − p)`; `worst = min over g∈0..10 of pnl(g)` (g=10 ≡ "10+";
P&L is constant above the top rung — unit test proves it). Every candidate
quote is admitted only if `worst(after hypothetical full fill)` stays within
caps; otherwise size shrinks to fit or the quote is dropped. Opposite-
direction fills raise the worst case and re-open room — two-sided quoting
prunes itself to one-sided as inventory skews, with no separate mode.

**Cap stack** (percent of `BANKROLL`, user-set 2026-07-10):

| cap | param | default | at $1,000 |
|---|---|---|---|
| per resting order | `MAX_QUOTE_PCT` | 0.30 | $300 (ledger binds first) |
| per match worst case | `MATCH_CAP_PCT` | 0.40 | $400 |
| global Σ worst cases | `GLOBAL_CAP_PCT` | 0.75 | $750 |
| daily realized loss halt | `DAILY_LOSS_HALT_PCT` | 0.40 | $400 |

Two-match day: global cap shrinks the second match's budget (400+400>750 →
$350). Daily halt = one match cap: one blown settlement stops the day.

**Positions source of truth:** Kalshi's `/portfolio/fills` and
`/portfolio/positions`, polled each tick and reconciled into `state.py` —
never trust only our own order bookkeeping (MLB lesson, 2026-03-22 session).
Unreconciled divergence → tripwire.

## 5. Tripwires (each: pull all quotes + loud log)

- **Crossed/impossible book** on any rung (`yes_ask + no_ask < 1 − 2·fee`):
  data error, pull match.
- **Fill burst**: more than `FILL_BURST_N` (default 3) fills on one match in
  60s → informed flow we don't see; pull match, `COOLOFF_MIN` (10 min).
- **Anchor regression**: ladder vanished / overround spike / blur pattern →
  pull match (feed-integrity failure from the adversarial review).
- **Reconciliation mismatch** between local fills and Kalshi positions →
  pull everything (never quote on state we can't trust).
- **Daily halt** (§4). Kill file `.kill` works as today.

Heartbeat gains: `quotes_live`, `fills_recent`, `worst_case_total`,
`halted` — "broken vs quiet" stays distinguishable.

## 6. Data model — new sibling DB `unabated_edge_maker.duckdb`

| table | row = | columns (abridged) |
|---|---|---|
| `maker_quotes` | every quote decision | ts, sport, event_id, ticker, side, action(rest/replace/cancel/skip), price, size, fair, margin, alt, reason, order_id |
| `maker_fills` | every real fill | ts, order_id, ticker, side, price, contracts, fee, ledger_worst_after |
| `ledger_snapshots` | match per tick | ts, event_id, worst_case, pnl_grid JSON, side_pulled |
| `maker_pnl` | match settled | event_id, spread_captured, fees, settlement_pnl, adverse_marks JSON |

Adverse-selection marks: each fill marked against anchor fair at +60s and
+5min (did fair move through us?) and settlement. The v1 verdict is
**quoted margin vs realized adverse-selection cost per rung class**
(main vs alt) — the same number the MLB maker v1 was built to measure.

## 7. Order lifecycle (LiveGateway)

- Place: `POST /portfolio/orders` (limit, side yes/no, price cents,
  `client_order_id` = deterministic hash of (ticker, side, tick, price) for
  idempotency on retry).
- Cancel: `DELETE /portfolio/orders/{order_id}`; replace = cancel+place (no
  amend on REST).
- Fills: `GET /portfolio/fills` cursor-polled each tick.
- All via `kalshi_common.auth_client` (already configured in `venues/kalshi.init`).
- `MAKER_MODE=live` refuses to start unless `MAKER_LIVE_ACK=1` is also set
  (live-by-typo impossible). `MAKER_MODE=shadow` records intended quotes to
  `maker_quotes` and infers fills from the tape (strict/queue bounds) — used
  for future leagues, not v1's path.

## 8. Config additions (all `.env`-overridable)

`ROI_MARGIN=0.03`, `PICKOFF_BUFFER_CENTS=1`, `MAX_MARGIN_CENTS=5`,
`ALT_MARGIN_MULT=1.5`, `ALT_SIZE_MULT=0.5`, `ALT_OVERROUND_MIN=1.01`,
`ALT_OVERROUND_MAX=1.15`, `QUOTE_PULL_MIN=3`,
`MAX_QUOTE_PCT=0.30`, `MATCH_CAP_PCT=0.40`, `GLOBAL_CAP_PCT=0.75`,
`DAILY_LOSS_HALT_PCT=0.40`, `FILL_BURST_N=3`, `COOLOFF_MIN=10`,
`MAKER_MODE=live`, `MAKER_LIVE_ACK` (unset). Existing reused:
`MAX_STALENESS_SEC=20` (now enforced), `TRADES_POLL_SEC=30`, `BANKROLL=1000`.

## 9. Testing

- `ledger.py` heaviest: hand-computed grids — single fill, same-direction
  stack, the Under-2.5/Over-1.5 offset (worst −$5, middle at g=2),
  cap-shrinking to fit, grid-truncation constancy above top rung.
- `engine.py`: margin floor (fee+buffer vs ROI), never-cross cap, skip-rung
  when crowd tighter than MAX_MARGIN, alt gating, requote-on-move, staleness
  pull, 3-min pull window, global-cap shrink on second match.
- `gateway.py`: Live place/cancel payloads against a mocked auth_client
  (idempotent client_order_id), MAKER_LIVE_ACK refusal, Shadow record path.
- `state.py`: fills reconciliation, divergence tripwire.
- Existing 55 tests stay green; runner tick test extended for the maker hook.

## 10. Rollout

```
today:   build + tests → executive review of diff → USER GATE: approve merge
         to main + first live start (one approval moment)
Jul 11:  Argentina–Switzerland, Norway–England — first live pre-kick windows
Jul 12:  measurement report: fills, spread captured, adverse marks, P&L per
         rung class → tune (margins, pull window, alt policy)
Jul 14:  France–Spain semi with tuned params
Jul 19:  final — sized on a week of measured data
```

No merge to main, no live order, without explicit user approval at the gate.

## 11. Version control, worktree, documentation

- **Worktree/branch:** all work on `worktree-wc-maker-depth-capture` (this
  branch already carries capture commit `8299cce`; the runner currently
  executes from this worktree and its DuckDBs hold live capture data — do
  NOT remove the worktree until merged AND data folded into main's DBs).
- **Commits:** spec (this file) → maker module + tests → README/CLAUDE.md
  docs, each its own commit; pre-merge executive review of
  `git diff main..HEAD` per repo rules.
- **Documentation (same merge):** `unabated_edge/README.md` maker section
  (architecture, config table, runbook, kill/halt semantics); root
  `CLAUDE.md` unabated_edge bullet gains the maker; memory update after
  merge.
- **Post-merge:** restart runner from main cwd (restart-gotcha memory),
  verify heartbeat + first quotes, fold worktree capture DBs into main's,
  then remove worktree + branch.
