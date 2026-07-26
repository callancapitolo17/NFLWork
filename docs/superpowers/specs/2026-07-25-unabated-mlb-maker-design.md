# Generalized Unabated-Anchored Maker — MLB First

**Date:** 2026-07-25
**Status:** Draft for review
**Predecessor:** `2026-07-10-wc-totals-maker-design.md` (WC totals maker), `2026-07-18-wc-touchjoin-pricing-design.md` (touch-join pricing)

---

## Review Pack

**What we're building** — Generalize the World Cup Kalshi maker
(`unabated_edge/maker/`) into a sport-agnostic framework and onboard MLB run
totals as the first sport. Unabated sharp anchors in, devigged fair out,
touch-join/fair±margin quotes on Kalshi MLB totals ladders. Go live at tiny
size with safety caps engaged; passive book capture runs alongside as
measurement context, not as a go-live gate.

**Key decisions**

1. **Generalize through onboarding** (rejected: framework-first build of
   issues #8/#9/#1; rejected: thin `soccer.py` fork). Generic pieces get
   built because MLB needs them, so nothing is speculative and MLB ships
   soonest.
2. **Totals only, interfaces shaped for spreads/ML later** (rejected:
   totals+spreads together). Smallest step that still fixes the
   architecture; spreads add new ledger/mapping surface better done after
   MLB proves the pattern.
3. **Keep the WC trading strategy unchanged** — touch-join when fair clears
   the crowd's touch by ≥1¢ net of maker fee, fair±margin fallback when the
   book is wide, per-rung whichever is strictly better, queue hysteresis.
   It is regime-adaptive, which is what an unknown book requires. No new
   pricing logic in v1.
4. **Capture is a dashboard, not a gate** (user call, 2026-07-25; rejected:
   scan-gated go-live). We quote live regardless of what the crowd looks
   like and use the captured tape to interpret results. Trial-and-error is
   the chosen learning mode.
5. **Generic interval ledger replaces the goal grid.** Worst-case exposure
   is evaluated over the intervals between rung strikes — sport-agnostic by
   construction; soccer migrates onto the same code with behavior unchanged.

**Risks / push back here**

- **MLB's book may be WC-like** (1¢ crowd, deep queues). If so, live
  quoting yields ~zero fills and the fills we do get skew adverse. You
  accepted this to keep momentum — the capture data will tell us quickly,
  but only after we're exposed. Push back if you'd rather re-add the gate.
- **Leashes ON is a step back from your WC-final config.** v1 runs with
  `MAKER_MAX_CONTRACTS` small (2–5) and a real `HARD_STOP_DOLLARS` (~$50
  mark-to-anchor quoting halt). You ran the WC final uncapped; confirm you
  accept caps for a new sport's first live run.
- **Third uncoordinated Kalshi MLB bot.** RFQ taker + SGP combo maker +
  this singles maker share game exposure with no cross-bot ledger.
  Accepted for v1 (measured, not prevented) — flag if you want a shared cap
  instead.
- **MLB news risk is faster than soccer's.** Pitcher scratches and lineup
  news move totals sharply pre-game. Our only defense is the kickoff-aware
  anchor-staleness gate; event-window halts (like `kalshi_mlb_mm` issue
  #24's design) are out of scope for v1.

**Worth understanding**

1. **Piecewise-constant payoff** — a totals ladder's P&L only changes when
   the total crosses a rung strike, so "check every outcome" reduces to
   "check each interval between strikes." Like R's `cut()`: infinitely many
   scores, finitely many buckets, one evaluation per bucket.
2. **Adverse selection** — a maker's fills are not a random sample of flow;
   you get filled most when the other side knows something (the anchor just
   moved and you haven't repriced). Spread earned must exceed this cost, and
   it's why "zero fills" and "only bad fills" are the two failure modes.
3. **Queue priority** — Kalshi fills price-time priority; canceling and
   re-quoting sends you to the back of the line. This is why the WC final
   produced 5 correct quotes and 0 fills, and why quote churn is a
   first-class cost, not a cosmetic one.

---

## 1. Goals and non-goals

**Goals**

- One maker framework where onboarding sport N+1 = one adapter file + one
  config block + one registry line.
- MLB run totals quoted live at tiny size during the 2026 season (runway
  through October; NFL/CFB onboard from this framework in September).
- Realized measurement the WC never produced: fills, markouts vs anchor,
  quoted-margin-vs-adverse-selection, per pricing mode (`touch_join` vs
  `quote`).

**Non-goals (v1)**

- Spreads and moneylines (interfaces accommodate them; no implementation).
- Proprietary modeling of any kind — fair is always the devigged anchor
  (standing constraint on this engine).
- Cross-bot MLB exposure coordination; event-window news halts; WebSocket
  feeds (issue #5); crowd-keyed breaker (issue #1) unless MLB shows the
  WC's queue-reset pattern.

## 2. Architecture

```
Unabated v2 per-league feed ──▶ feed.py (poll, parse)          [unchanged]
                                  │
                          sports/registry.py                    [+1 line]
                          ├── soccer.py                         [unchanged]
                          └── mlb.py            ◀── NEW: bt3 totals → KX rungs
                                  │  Candidate(ticker, side, fair_prob, label)
                                  ▼
                          maker/engine.py                       [unchanged]
                          (touch-join + fair±margin, hysteresis)
                                  │
              ┌───────────────────┼──────────────────┐
              ▼                   ▼                  ▼
        maker/ledger.py     maker/gateway.py    maker/state.py
        REWORKED: generic   [unchanged: v2     [unchanged: fills,
        interval worst-case  orders live/shadow] caps, kill file]
```

One runner process serves all enabled sports; each adapter polls its own
league file on the existing `V2_POLL_SEC` cadence. DBs unchanged:
`unabated_edge_market.duckdb` (snapshots/book/trades),
`unabated_edge_research.duckdb` (firehose), maker state in the existing
store.

## 3. Components

### 3.1 Task 0 — build inputs (live probing, ~hours, read-only)

Not a study; these are facts the adapter cannot be written without.
Deliverable: constants recorded in the sport config + a short note in this
spec's implementation plan.

- Unabated MLB league id in the v2 feed; confirm sharp anchors quote bt3
  run totals; confirm alt-ladder presence and `modifiedOn` cadence vs first
  pitch; confirm side convention (expect si0=Over/si1=Under as in soccer —
  verify, don't assume).
- Kalshi MLB total series: ticker (pattern suggests `KXMLBTOTAL` — verify),
  event/market structure, `floor_strike` convention (expect X.5 rungs),
  `orderbook_fp` dollar-string format parity.
- Event pairing: Kalshi event title team-pair vs Unabated event teams; MLB
  team-name normalization (reuse patterns from existing MLB canonical-match
  code where sensible, but the adapter stays self-contained).

### 3.2 Generic interval ledger (replaces goal grid)

`maker/ledger.py` currently enumerates goals 0..10. Replacement: given all
rungs the maker touches for one game (resting quotes counted as
hypothetically filled, as today), collect their strikes
`s_1 < s_2 < ... < s_k`, form intervals
`(-inf, s_1), (s_1, s_2), ..., (s_k, inf)`, evaluate portfolio P&L once per
interval, take the minimum. Exact for any totals ladder in any sport — no
per-sport grid, no truncation error (the goal grid's 0..10 cap becomes
irrelevant).

Soccer's existing ledger tests migrate to the new implementation and must
produce identical worst-case numbers — that parity is the proof of no
behavior change.

### 3.3 `sports/mlb.py` adapter

Same contract as soccer (`price_event()` → candidates, `event_teams()` for
pairing): devig anchor bt3 over/under per line, main-wins-over-alt
setdefault, match Kalshi rung by strike, over = buy YES / under = buy NO,
rung provenance {book, alt, overround} to the research firehose, fail
closed everywhere.

### 3.4 Per-sport config

Soccer constants (league id, series prefix, staleness windows) move into a
per-sport config block consumed by adapter + registry. Env override pattern
unchanged (repo-root `.env` — verify location before launch; flagged in the
WC post-mortem).

### 3.5 Capture-as-dashboard

The existing Kalshi book/trades capture points at MLB totals markets from
day one of the build. Standing readout (a query file, not a service):
per-rung spread distribution, depth, % of trades within ±1¢ of mid, queue
turnover, anchor coverage. Interprets live results; gates nothing.

**Disk guard:** WC capture died at 98% disk and the market DB grew
~675MB/day. Before starting: check free space; cap MLB capture (orderbook
top-of-book + trades only, retention pruned) so a 15-game slate doesn't
recreate the problem.

## 4. Go-live plan

1. **First-tick shadow check** (same session as launch, ~minutes, user
   call 2026-07-25 — no standalone shadow day): run the pipeline with the
   shadow gateway for the first poll cycles; verify pairing count, rung
   matches, candidate counts, heartbeat. The order path is unchanged
   WC-live-proven code, so shadow only needs to prove the new MLB wiring.
2. **Flip live, same slate**: `MAKER_MODE=live MAKER_LIVE_ACK=1`,
   `MAKER_MAX_CONTRACTS=2–5`, `HARD_STOP_DOLLARS≈50`, kickoff-aware
   `ANCHOR_STALE_SEC` tuned to MLB anchor cadence from Task 0, kill file
   armed, %-of-bankroll caps as in WC config.
3. **Review after ~1 week of slates**: fills by mode, markouts, worst
   rungs; decide scale-up, strategy changes (issue #1 breaker?), or stop.
   The multi-sport expansion epic is written here, ranked on real fill
   data (re-triage of unabated_edge issues #1–#10 plus deferred spec
   items: spreads/ML, shared MLB exposure cap, news halts, NFL).

Restart discipline: never mid-game while holding inventory (fills ledger
rebuild caveat carries over from WC README).

## 5. Testing

- Full existing suite stays green (131 at WC merge).
- Interval ledger: soccer-parity tests (identical worst case vs goal grid)
  + MLB-rung cases + property test (adding a losing interval never
  improves worst case).
- MLB adapter: fixture-based feed→candidate tests (mirroring soccer's),
  including pairing normalization and fail-closed paths.
- Live probes for gateway are already covered; no order-path changes.

## 6. Version control, worktree, documentation

- **Worktree:** `worktree-unabated-mlb-maker` (created 2026-07-25, this
  branch). All build work happens here; merge to `main` only after tests +
  pre-merge review + explicit user approval; clean up worktree + branch
  after merge.
- **Files:** new `unabated_edge/sports/mlb.py`, reworked
  `maker/ledger.py`, per-sport config, capture pointing + readout queries,
  tests. No changes to gateway/state/engine expected.
- **Commits:** Task-0 findings (config constants) → ledger generalization
  (with parity tests) → MLB adapter → capture/readout → docs.
- **Docs (same merge):** `unabated_edge/README.md` (multi-sport onboarding
  steps, MLB section, capture readout usage), root `CLAUDE.md` blurb
  (unabated_edge entry gains MLB + generic-ledger note).
