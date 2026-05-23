# RFQ Research Logging — Design Spec

*Date: 2026-05-23 · Branch: `worktree-rfq-research-logging` · Component: `kalshi_mlb_rfq/`*

---

## Review Pack

**What we're building** — A second observability layer for the Kalshi MLB RFQ *taker* bot whose single goal is to let us **fully understand the RFQ process end to end**. Today the bot keeps a solid audit trail of what it *acted on* (`quote_log`, `fills`), but everything it *considered and discarded* — rejected candidates, per-book fair values, gate internals, Kelly math — is computed in memory and thrown away. We add a separate "research firehose" that records the full lifecycle of every candidate combo so any one can be replayed step by step, plus we upgrade the bot's primitive `print()` logging to the standard Python `logging` module with rotation.

**Key decisions**

1. **Separate research DuckDB (`kalshi_mlb_rfq_research.duckdb`), written by the bot only, batched once per tick** — *rejected:* (a) adding tables to the trading DB, which would worsen the single-writer lock contention the bot already fights; (b) JSONL files, which add a second format; (c) a threaded async writer, which adds concurrency machinery we don't need. *Why:* a separate file has its own lock (zero contention with trading) yet keeps the familiar DuckDB SQL surface, and per-tick batching makes writes cheap.
2. **"Full movie" capture — one row per candidate per tick — with a `RESEARCH_CANDIDATE_SAMPLING` throttle knob** — *rejected:* log-only-on-change (leaner but loses tick-by-tick movement) and aggregate-counts-only (too coarse to replay). *Why:* seeing fairs/edges move minute-by-minute *is* the understanding we're after; ~290k rows/day is a non-issue on the batched separate DB; the knob is the safety valve.
3. **One wide `events` table with a flexible JSON `payload` column** — *rejected:* a table per event type. *Why:* adding a new event type or field needs no migration, and analysis is uniform.
4. **No automatic pruning in v1** — the `prune_research()` helper is written but left **unscheduled**. *Why:* user's call; growth is slow (~2.5 GB/month) so there's months of runway, and enabling retention later is a zero-code change (point a cron line at the helper).
5. **Outcome / settlement / CLV capture is explicitly out of scope** — *rejected* for now. *Why:* the goal is understanding the RFQ *process*, not building a track record or product. The same firehose can be extended with deferred outcome events later if that goal changes.

**Risks / push back here**

- **No prune = unbounded growth.** Accepted as a deliberate choice. It's slow (~2.5 GB/month) and `bot.log` is independently rotation-capped, so nothing runs away in the near term — but it is technically unbounded until the helper is scheduled.
- **The lifecycle events require touching hot trading code.** To log per-book fairs, gate headroom, and Kelly internals, helper functions (`fair_value._load_book_fairs`, `risk.*`, `kelly.*`) must *return* data they currently discard. This is behavior-preserving refactoring of code on the trading path — the risk is a subtle behavior change. Mitigation: refactors preserve existing return values, emit happens at the call sites, and tests pin the existing outputs.
- **Capture volume.** ~290k rows/day at full movie. Fine for 2b, but if you'd rather start lean, say so and we default the sampling knob below 1.0.

**Worth understanding** (opt-in, anchored to R)

- **Buffer-then-flush batching.** `emit()` appends event dicts to an in-memory list during a tick; `flush()` writes the whole list in one transaction at the end. This is the R habit of building up a data frame and writing it once with `dbWriteTable(...)` instead of looping single-row `INSERT`s — one expensive operation instead of hundreds of cheap-but-contended ones.
- **DuckDB's single-writer file lock.** A `.duckdb` file allows only one read-write connection at a time, like opening a file in exclusive-write mode. The current bot opens/closes that lock many times per tick, which is why `db.py` has a retry-backoff loop. A *separate* file has a *separate* lock — so research writes can't collide with trading writes.
- **Python `logging` vs `print()`.** A long-running process wants log *levels* (so you can filter `WARNING`-and-up in prod, flip to `DEBUG` when debugging, no code change) and a `RotatingFileHandler` (so the file can't grow forever). `print()` gives you neither.

---

## 1. Goal & scope

**Goal:** make the full RFQ lifecycle of every candidate combo reconstructable, so we can answer questions the current logging can't — "what did we pass on and why," "what did each book think vs. our model," "how close was that gate to firing," "what did Kelly actually want before flooring."

**In scope:** capture of the enumeration → pricing → selection → submission → quote → evaluation → sizing → accept/walk → position-update lifecycle; replacement of `print()` with structured rotating logging; an analysis/query surface.

**Out of scope (v1):** outcome/settlement capture, CLV computation, any product/data-service layer, a dashboard, automatic pruning, Parquet archival.

---

## 2. Architecture — "Approach 2b"

Three sibling DuckDB files in `kalshi_mlb_rfq/`:

| DB | Role | Writer | Lock pressure |
|----|------|--------|---------------|
| `kalshi_mlb_rfq.duckdb` (state) | what we acted on: `quote_log`, `fills`, `live_rfqs`, `positions`, `sessions`, `combo_cooldown` | bot (frequent, per-op) | high — existing bottleneck |
| `kalshi_mlb_rfq_market.duckdb` (market) | SGP lines + odds cache | bot + scraper subprocess | moderate |
| **`kalshi_mlb_rfq_research.duckdb` (research) — NEW** | the firehose: full lifecycle events | bot only, batched once per tick | **isolated — its own lock** |

The research DB is never read or written by anything except the bot's own `flush()` and ad-hoc analysis (read-only `ATTACH`). Because it is a distinct file, its write lock is independent of the trading DB's, so research logging cannot add contention to the trading path.

---

## 3. The research store — `research.py` + schema

New module `kalshi_mlb_rfq/research.py`. Public surface:

```python
set_session(session_id: str) -> None
    # Called once at startup so every event is stamped with the session.

emit(event_type: str, *, game_id=None, combo_ticker=None, rfq_id=None,
     quote_id=None, **payload) -> None
    # Append one event dict to the in-memory buffer. O(1). NEVER raises.

flush() -> None
    # Write the whole buffer to the research DB in ONE transaction, then clear it.
    # Wrapped in try/except: any failure is logged (operational logger) and the
    # buffer is dropped/retained per the cap rule. NEVER raises into the loop.

init_research_db() -> None
    # Idempotently create the research DB + events table + indexes. Called at startup.

prune_research(days: int = 90) -> int
    # DELETE FROM events WHERE ts < now - <days>. Returns rows deleted.
    # WRITTEN BUT NOT SCHEDULED in v1 — available for a future cron line.
```

**Buffer & safety rules:**
- The buffer is a module-level `list`. `emit()` only appends — it physically cannot throw into the hot path.
- **Cap:** if the buffer exceeds `RESEARCH_BUFFER_MAX` (default 50,000) without a successful flush, drop the oldest events. A persistent flush failure degrades to "we lose research data," never "the bot OOMs or stalls."
- `flush()` opens its own connection to the research DB (mirroring the `db.py::connect()` context-manager pattern with retry-backoff), writes, and closes it. It is the **only** path that holds the research DB write lock.
- `flush()` is called once per main-loop iteration (after the tick's work) **and once more in the shutdown handler** so the final tick's buffer isn't lost on a clean exit.

**Schema:**

```sql
CREATE TABLE IF NOT EXISTS events (
    event_id     VARCHAR PRIMARY KEY,    -- uuid4
    session_id   VARCHAR,                -- joins to state.sessions
    event_type   VARCHAR NOT NULL,       -- 'candidate_evaluated', 'gate_evaluated', ...
    ts           TIMESTAMPTZ NOT NULL,   -- UTC (repo TZ convention)
    game_id      VARCHAR,
    combo_ticker VARCHAR,
    rfq_id       VARCHAR,
    quote_id     VARCHAR,                -- links to state.quote_log when applicable
    payload      JSON                    -- event-specific fields; no migration to add one
);
CREATE INDEX IF NOT EXISTS idx_events_type_ts ON events(event_type, ts);
CREATE INDEX IF NOT EXISTS idx_events_game    ON events(game_id);
CREATE INDEX IF NOT EXISTS idx_events_combo   ON events(combo_ticker);
```

---

## 4. Event catalog — by RFQ lifecycle stage

Research events **complement** the state DB; they never duplicate it. Linkage keys: `leg_set_hash`/`combo_ticker` → `rfq_id` → `quote_id`.

| # | Lifecycle stage | Event | Payload (the discarded data we now keep) | Today |
|---|---|---|---|---|
| 1–3 | Enumerate → Price → Select | `candidate_evaluated` (per candidate per tick) | legs (spread/total/side), `model_fair`, `book_fairs` {dk, fd, pinnacle, …} + `n_books`, `blended_fair`, `kalshi_ref`, `edge`, `rank`, `outcome` ∈ {`submitted`, `rejected_no_mapping`, `rejected_no_book_data`, `rejected_fair_oob`, `rejected_below_topN`} | silent (gaps A+B) |
| 4 | Submit RFQ | `rfq_submit_failed` | add-failure exception text (success already in `live_rfqs`) | `print()` only |
| 5 | Receive quote | `quote_priced` (keyed `quote_id`) | the components behind `quote_log.blended_fair_at_eval`: `model_fair` + each book's fair at eval time | only the blend kept |
| 6 | Evaluate / gates | `gate_evaluated` (per quote) | each gate pass/fail + headroom: per-game spend vs cap, daily spend vs cap, fill-ratio window state, cooldown seconds left, staleness age, tipoff minutes | only final label in `quote_log` |
| 6.5 | Size (Kelly) | `kelly_sized` | base Kelly fraction, correlation adjustment vs existing positions, target fraction, target vs floored contracts, bankroll | silent (gap D) |
| 7 | Accept / walk | `walk_diagnosed` | parsed walk reason (`quote_expired` / `rfq_closed` / other) | raw error body only |
| 8 | Update position | `position_snapshot` (on each fill) | ticker, side, net_contracts & weighted_price **before and after** | mutated in place (gap D) |

**Capture rate:** `candidate_evaluated` defaults to one row per candidate per tick ("full movie"). `RESEARCH_CANDIDATE_SAMPLING` (default `1.0` = log all) can throttle it later without code changes. All other events are low-volume (only fire on real quotes/fills) and are always logged.

**Replay example** (the goal made literal): for one `leg_set_hash` you can read, in time order — *enumerated 7:02 → model 0.41, DK 0.39, FD 0.40, blended 0.40, edge +3% → ranked #4, submitted RFQ X → MM quoted 0.37 at 7:03 → EV +8%, all gates pass except fill-ratio (2/50 headroom) → Kelly wanted 1.2, floored to 1 → accepted → position 0 → 1.*

---

## 5. Operational logging — gap C

Replace `print(..., flush=True)` throughout with Python's `logging`, configured once in a new `kalshi_mlb_rfq/log_setup.py` and imported at startup.

- **Levels:** existing prints map to `DEBUG` (verbose internals), `INFO` (cache refresh, heartbeat, RFQ counts), `WARNING` (recoverable errors: `rfq_refresh error`, `quote_poll error`, **research `flush()` failures**), `ERROR`/`CRITICAL` (halts, math-invariant breaks).
- **Handlers:** a `StreamHandler` to stdout (so systemd/supervisor still captures everything) **plus** a `RotatingFileHandler` on `bot.log` — `maxBytes` ≈ 50 MB, `backupCount` = 5 → bounded at ~300 MB total. Both numbers in `config.py`.
- **Format:** `%(asctime)s %(levelname)s %(name)s: %(message)s`.
- **`notify.py`** keeps its `[FILL]`/`[HALT]` webhook alerting; its file writes route through the rotating logger.

This is independent of the research firehose and can ship first (Phase 0).

---

## 6. Analysis surface

The research DB is a normal DuckDB file, so analysis is plain SQL via read-only `ATTACH` + cross-DB join:

```sql
ATTACH 'kalshi_mlb_rfq.duckdb'          AS state    (READ_ONLY);
ATTACH 'kalshi_mlb_rfq_research.duckdb' AS research (READ_ONLY);

-- JSON fields via ->> :
SELECT ts,
       payload->>'model_fair'  AS model_fair,
       payload->>'blended_fair' AS blended_fair,
       payload->>'edge'        AS edge
FROM research.events
WHERE event_type = 'candidate_evaluated' AND combo_ticker = ?;
```

Ship **three canned example queries** (documented in README, optionally as DuckDB views):
1. **RFQ journey** — all events for one `leg_set_hash`/`combo_ticker` in time order (full replay).
2. **Missed edges** — `candidate_evaluated` where `outcome <> 'submitted'` but `edge` above a threshold.
3. **Gate breakdown** — counts of which gate rejected accepts over a date range.

**Retention:** none automatic in v1. `prune_research(days=90)` exists, unscheduled.

---

## 7. Error handling & safety (the non-negotiable)

**Research logging must never affect trading.**
- `emit()` only appends to a list — cannot raise into the loop.
- `flush()` catches every exception, logs a `WARNING` via the operational logger, and continues. A locked/full/erroring research DB never propagates to trading.
- Buffer cap (`RESEARCH_BUFFER_MAX`) bounds memory if flush keeps failing.
- The research connection is fully separate from the trading/market connections; `flush()` closes it on every call (no lingering lock).
- `flush()` is quick (one batched write) so it cannot starve the shutdown/SIGTERM path the way the long `rfq_refresh` block can (see restart-gotchas note).

---

## 8. Testing

The bot trades live, so tests must run **without** Kalshi:
- `research.py` units: `emit` appends; `flush` writes one batch and clears the buffer; **a forced DB error inside `flush` is swallowed and never raises** (the critical safety test); buffer cap drops oldest when exceeded; `init_research_db` is idempotent; `prune_research` deletes only rows older than N days.
- Refactor-safety: `fair_value._load_book_fairs` (and `kelly`/`risk` helpers) still return their original values after being extended to expose internals.
- Smoke test: init research DB on a temp path → emit a few events of each type → flush → query them back and assert shape.

No new scraper is added, so `tests/timezone_parity_test.py` does not apply — but the `ts` column follows the TIMESTAMPTZ-UTC convention regardless.

---

## 9. Phased rollout

Each phase is independently shippable and testable:

- **Phase 0 — Operational logging (C).** `print` → `logging` + rotation. Self-contained; nothing depends on it.
- **Phase 1 — Research plumbing.** `research.py` (`emit`/`flush`/`init_research_db`/`prune_research`) + schema + wire `flush()` into the loop and shutdown. Prove it writes safely and never blocks trading before emitting real events.
- **Phase 2 — Big unlock (A+B).** `candidate_evaluated`. Requires `_load_book_fairs` to return the per-book breakdown.
- **Phase 3 — Rest of lifecycle.** `gate_evaluated`, `kelly_sized`, `position_snapshot`, `quote_priced`, `walk_diagnosed`, `rfq_submit_failed`.
- **Phase 4 — Analysis.** Example queries + views + README docs.

---

## 10. Version control

- All work on branch `worktree-rfq-research-logging` (this worktree).
- **New files:** `kalshi_mlb_rfq/research.py`, `kalshi_mlb_rfq/log_setup.py`, tests under the bot's test path.
- **Modified files:** `main.py` (wire `flush` + `emit` calls, replace prints), `config.py` (new knobs: `RESEARCH_DB_PATH`, `RESEARCH_CANDIDATE_SAMPLING`, `RESEARCH_BUFFER_MAX`, `LOG_MAX_BYTES`, `LOG_BACKUP_COUNT`), `fair_value.py` (return per-book fairs), `kelly.py` (return sizing diagnostics), `risk.py` (return gate headroom), `notify.py` (route through logger).
- **gitignore:** confirm `bot.log` and `*.duckdb` (incl. the new research DB) are ignored. Add `bot.log*` (rotated suffixes) if not covered.
- Commits structured per phase. The live bot is **not** run from the worktree (per restart-gotchas: it must launch from the main repo cwd); validation is via unit/smoke tests, and live restart happens after merge to `main`.
- Pre-merge: executive-engineer review of the full diff; explicit user approval before merge.

---

## 11. Documentation

Updated in the same merge to `main`:
- `kalshi_mlb_rfq/README.md` — new "Observability / research logging" section: the three DBs, event catalog, the `ATTACH` query pattern, config knobs, the (unscheduled) prune helper.
- Repo `CLAUDE.md` — extend the `kalshi_mlb_rfq` bullet to mention the third sibling research DB and structured operational logging.
- Memory — add a note on the research-logging architecture once shipped.

---

## 12. Open questions / future

- Turn on `prune_research()` (cron) once the DB gets chunky, or add Parquet archival for long-run history.
- If the goal later expands from "understand the process" to "prove edge / build a data product," add deferred outcome + closing-line events keyed to `combo_ticker`/`quote_id`.
- A read-only analysis dashboard over the research DB (not v1).
