# Kalshi MLB Bots Monitor — Design Spec

_2026-06-18 · branch `worktree-kalshi-mlb-monitor`_

## Review Pack

**What we're building** — A standalone, auto-refreshing local web dashboard that
monitors **both** Kalshi MLB bots on one screen: the live **maker**
(`kalshi_mlb_mm`) and the **taker** (`kalshi_mlb_rfq`). It answers four questions —
how many RFQs, how many fills, what's outstanding, and *why am I not getting
filled* — and surfaces the strategic signal already buried in the bots' DuckDBs.
It is **read-only**: it never touches the bots or places orders.

**Key decisions**

1. **Standalone Dash app, not an extension of `kalshi_draft`.** Reuses the proven
   `kalshi_draft/app.py` patterns (dark theme, `dcc.Interval` refresh, read-only
   connections) but lives in its own dir/port. Rejected extending the draft portal:
   different domain, would couple unrelated tools. Rejected Streamlit: not in the stack.
2. **Direct read-only DuckDB reads, no snapshotting.** Empirically verified the
   monitor can read the maker DB while the live bot holds the write lock (DuckDB
   1.4.4 serves the last checkpoint). Rejected copy-on-refresh snapshots: unnecessary
   complexity once direct reads proved safe. Transient checkpoint-lock errors are
   caught and the last good render is kept.
3. **Per-bot "adapter" config abstracts schema differences.** The two bots have
   different column names (maker `fills.side_held/price/fee` vs taker
   `fills.side/price_dollars/fee_dollars`; taker positions carry `legs_json`, maker
   don't). One adapter dict per bot maps logical fields → real columns so the render
   code stays generic. Rejected duplicating every query twice.
4. **Reason vocabularies are read from the data, not hardcoded.** Decision/reason
   charts `GROUP BY` whatever strings exist. Grounding in the live DBs already
   revealed codes the source-reading missed (`declined_sanity`, `size_gate_dollars`),
   so hardcoding a list would silently drop categories.
5. **Stale-aware, not live-assuming.** Each bot panel shows "last activity Xm ago";
   the taker (down since Jun 3) gets an explicit STALE banner so 15-day-old history
   isn't read as current.

**Risks / push back here**

- **Taker is dormant.** Its tab is historical until you restart the bot. If you'd
  rather the monitor only cover what's live, say so — but you asked for both, and the
  taker's 1,103-fill history is worth seeing.
- **The maker's story is brutal and the dashboard makes it loud:** 296,980 RFQs seen,
  99.996% out-of-scope, 13 ever addressable, **0 fills.** The dashboard's job is to
  show this honestly, not flatter it. If you want it framed differently (e.g. "scope
  coverage" as an opportunity metric vs a failure metric), that's a wording call.
- **Settlement P&L is unavailable.** Maker `fills.realized_pnl` is null until games
  settle, and the maker has 0 fills anyway. P&L views show "pending"/realized-only.
  No settlement backfill in v1.

**Worth understanding** (opt-in)

- **The funnel is the core mental model.** Like an R `dplyr` pipeline that `filter()`s
  rows away at each step, an RFQ flows Seen → In-scope → Quoted → Filled, and each
  stage drops candidates. "Why am I not getting filled" = *where* the funnel collapses.
  For the maker it collapses at stage 1 (scope); for the taker, at the quote-eval gates.
- **Lock-safe reads.** DuckDB allows one writer OR many readers per file across
  processes — but here a separate read-only process can still read a writer-held file
  because it sees the last committed checkpoint. We wrap reads so a rare lock error
  degrades gracefully (keep last render) instead of crashing the page. Like a `tryCatch`
  in R that returns the previous value on error.

---

## 1. Architecture

```
kalshi_mlb_monitor/
  app.py        # Dash layout + callbacks (entrypoint)
  queries.py    # lock-safe read-only DuckDB helpers + per-bot adapters
  bots.py       # BOT adapter config (db paths, column maps, reason groupings)
  run.sh        # launcher (port 8092, env-overridable)
  README.md
```

- **Stack:** Python Dash + Plotly + DuckDB. Mirrors `kalshi_draft/app.py`.
- **Read-only:** every DB connection opens with `read_only=True`. The monitor never
  writes to any bot DB and imports no bot code (zero coupling, zero risk to the live
  $500 maker).
- **Refresh:** `dcc.Interval` default 30 s (UI toggle 10 s / 30 s / off). Each callback
  opens a short-lived read-only connection, queries, closes. A `QueryLocked` sentinel is
  returned on transient lock errors → callback raises `PreventUpdate` (keeps last render).
- **Port:** default `8092`, override via `KALSHI_MLB_MONITOR_PORT`.
- **Theme:** reuse the `kalshi_draft` GitHub-dark `COLORS` palette and card/table styles.

### 1.1 Bot adapter (`bots.py`)

A dict per bot describing where data lives and how columns map:

```python
BOTS = {
  "maker": {
    "label": "Maker (kalshi_mlb_mm)",
    "state_db":    ".../kalshi_mlb_mm/kalshi_mlb_mm.duckdb",
    "research_db": ".../kalshi_mlb_mm/kalshi_mlb_mm_research.duckdb",
    "rfq_table":   "seen_rfqs",      # one row per RFQ observed
    "rfq_ts":      "first_seen_at",
    "decision_table": "quote_decisions",
    "decision_ts": "observed_at",
    "fills": {"side": "side_held", "price": "price", "fee": "fee", "ts": "filled_at"},
    "has_position_legs": False,
    "kind": "maker",
  },
  "taker": {
    "label": "Taker (kalshi_mlb_rfq)",
    "state_db":    ".../kalshi_mlb_rfq/kalshi_mlb_rfq.duckdb",
    "research_db": ".../kalshi_mlb_rfq/kalshi_mlb_rfq_research.duckdb",
    "rfq_table":   "live_rfqs",       # one row per RFQ the bot sent
    "rfq_ts":      "submitted_at",
    "decision_table": "quote_log",
    "decision_ts": "observed_at",
    "fills": {"side": "side", "price": "price_dollars", "fee": "fee_dollars", "ts": "filled_at"},
    "has_position_legs": True,
    "kind": "taker",
  },
}
```

DB paths derive from a repo root computed from the running file location (so it works
in a worktree and on `main`), mirroring the repo's existing root-derivation pattern.

---

## 2. Data model (verified against live DBs)

**Maker** `kalshi_mlb_mm.duckdb`: `seen_rfqs`, `live_quotes`, `quote_decisions`,
`fills`, `positions`, `sessions`, `combo_cooldown`. Current state: 296,993 RFQs seen
(13 in-scope), 297,929 decisions, 13 quotes (all cancelled), **0 fills, 0 positions**,
not dry-run.

**Taker** `kalshi_mlb_rfq.duckdb`: `live_rfqs`, `quote_log`, `fills`, `positions`,
`combo_cache`, `combo_cooldown`, `sessions`, `reference_lines`. Current state: 32,771
RFQs, 126,316 quote-log rows, **1,103 fills / $539 staked**, 708 positions; last fill
Jun 2; 160 orphaned `open` RFQs (bot down).

Verified reason/decision vocabularies (live `GROUP BY`):

- **Maker `quote_decisions.reason`:** `out_of_scope` (296,980), `no_fair` (671),
  `size_gate` (237), `size_gate_dollars` (22), `no_game` (5), `per_combo_cap` (1).
  `decision` ∈ {`skipped`, `quoted`}.
- **Taker `quote_log.decision`:** `declined_stale_predictions` (53k), `declined_ev`
  (50k), `halted_low_fill_ratio` (14k), `failed_quote_walked` (2,308),
  `declined_sanity` (1,952), `declined_tipoff`, `accepted` (1,103),
  `declined_inverse_lock`, `declined_kelly_zero`, `declined_side_mismatch`,
  `declined_per_game_cap`, `declined_dry_run`.

Queries never hardcode these — they aggregate whatever is present.

---

## 3. Layout & tabs

Header: title · **bot-selector toggle** `[ Maker ● live | Taker ]` · status strip
(running?, last-activity age, dry-run flag, current session start, STALE banner) ·
refresh control (interval toggle + manual refresh).

### Tab 1 — Overview
KPI cards + the **funnel**.
- **Maker KPIs:** RFQs seen, in-scope count + %, quotes sent, fills, fill rate,
  open quotes, open positions, exposure $. Funnel bar: Seen → In-scope → Quoted → Filled.
- **Taker KPIs:** RFQs sent, quotes evaluated, accepted (fills), fill rate, halt-rate,
  walk-rate, open positions, exposure $, last fill age. Funnel: Sent → Evaluated →
  Accepted, with halt/walk siphons annotated.

### Tab 2 — Why Not Filled ★ (the core ask)
- Horizontal bar of decision/reason counts over a selectable window (1h / 24h / 7d / all).
- Stacked time-series of decisions by reason (so you see *when* a reason spikes —
  e.g. taker `declined_stale_predictions` bursts).
- Plain-language legend mapping each code → what it means (sourced from the README
  defense hierarchy / gate list).
- Maker-specific callout: scope rejection rate + the 13 in-scope RFQs detail (each one's
  quote + why it cancelled).

### Tab 3 — Fills & P&L
- Recent fills table (time, game, combo, side, contracts, price, fee, fair-at-quote).
- Cumulative fills + cumulative stake over time.
- Maker empty-state: "0 fills yet — see **Why Not Filled**."
- P&L: realized where available (`realized_pnl`), else "pending settlement."

### Tab 4 — Positions & Exposure
- Open positions table (combo, side, net contracts, weighted price, exposure $, age),
  taker shows leg detail from `legs_json`.
- Exposure by game + caps utilization (daily / per-game / per-combo vs configured caps).
- Open quotes (maker `live_quotes`) / open RFQs (taker `live_rfqs status='open'`),
  with **orphaned-open** flag when the bot is down but rows are still `open`.

### Tab 5 — Adverse Selection / Fill Quality
- **Maker:** quoted margin (5% target) vs realized fair drift
  (`fair_at_confirm − blended_fair_at_quote`) per fill; void-rate / circuit-breaker /
  last-look event counts from research `events`. (Mostly empty-state today → shows the
  measurement scaffold so it lights up on first fills.)
- **Taker:** `halted_low_fill_ratio` timeline, `failed_quote_walked` rate +
  `best_competitor_no_bid_dollars` (how close we lost), `edge_at_submit` distribution
  for accepted vs walked.

---

## 4. Freshness & error handling

- **Freshness badge** per bot = `now − max(decision_ts)`. Green < 2 min, amber < 30 min,
  red otherwise. Taker shows STALE + "bot down since `max(sessions.ended_at)`".
- **Lock degradation:** read helper returns `QueryLocked` on DuckDB lock/IO error;
  callbacks `raise PreventUpdate` so the page keeps the last good values.
- **Empty states:** every table/chart renders a friendly placeholder when 0 rows
  (maker fills/positions today) rather than erroring.
- **No write path, no auth, no secrets** read or logged.

---

## 5. Out of scope (YAGNI for v1)

Order placement / any write · authentication · settlement P&L backfill · alerting/paging
· historical replay beyond what's in the DBs · cross-bot correlation analytics.

---

## 6. Version control / docs / testing

- **Worktree:** `worktree-kalshi-mlb-monitor` (already created). New dir
  `kalshi_mlb_monitor/` only; **zero edits to bot code**.
- **Commits:** (1) queries+bots+app skeleton, (2) tabs, (3) launcher+README+CLAUDE.md line.
- **Docs:** new `kalshi_mlb_monitor/README.md` (setup, launch, what each tab shows,
  troubleshooting); one line in root `CLAUDE.md` project structure; memory note.
- **Testing (render-based, per repo policy):** launch against live maker + historical
  taker DBs; verify each of the 5 tabs loads for both bots; verify maker empty-states,
  taker stale banner, funnel math (counts reconcile to totals), lock-safety while the
  maker bot writes concurrently. Manual inspection via browser/screenshot before merge.
- **Merge:** executive-engineer diff review → explicit user approval → merge → remove
  worktree + branch.
