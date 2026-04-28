# MLB Trifecta Dashboard Tab — Design

**Status:** Draft for review
**Date:** 2026-04-27
**Owner:** Callan
**Related:**
- `docs/superpowers/specs/2026-04-22-dk-trifecta-blend-design.md` (the pricing side — Plans #1 & #2 already merged)
- `Answer Keys/mlb_correlated_parlay.R` + `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the patterns this plan mirrors)

## Goal

Add a "Trifectas" tab to the MLB +EV Dashboard, next to "Parlays", that shows priced TRIPLE-PLAY and GRAND-SLAM specials with per-row manual-log capability. Inherits the bets-tab UI conventions (color-coded EV, filtering, in-place Action button) and the parlay-tab plumbing conventions (pricer writes to its own DuckDB table on each refresh, dedicated `placed_trifectas` table in dashboard DB for dedup, dedicated `/api/place-trifecta` endpoint).

Auto-placement is explicitly out of scope for this plan — we ship a manual-log path first, the same way the bets tab supports both manual and auto and the parlay tab has both code paths. A follow-up plan can layer auto-placement once Wagerzon's Specials placement endpoint is reverse-engineered.

## Sequencing relative to other plans

Both prerequisites are already on `main`:
- **Plan #1** (`feature/dk-trifecta-blend-scaffolding`, merged at `919ef1b`): blend scaffolding, `mlb_trifecta_sgp_odds` schema, `blend_dk_with_model` helper, R-side blend wiring.
- **Plan #2** (`feature/dk-trifecta-scraper`, merged at `cdf1354`): live DK SGP scraper invoked via `system2()` from `mlb_triple_play.R`, leg resolvers, vig 1.25.

This plan ("Plan #3" informally) executes immediately after the spec is approved. No external blockers.

## Non-goals

- **Auto-placement** to Wagerzon. The existing `parlay_placer` hits Wagerzon's parlay endpoint; Specials use a different endpoint (per the existing scraper recon notes). A separate plan handles that.
- **Hot-swap fragment endpoint.** Trifectas are bounded at ~14 rows (max 30 MLB games × 2 sides × {TRIPLE-PLAY, GRAND-SLAM} where Wagerzon posts ~7 props/day). Full table re-render on next refresh is fine; we use bets-tab-style in-place button updates for placement instead.
- **Cross-game trifectas** ("ALL WIN" multi-team props). The pricer already filters these out via `parse_legs()` returning `NULL` for unsupported leg types.
- **CLV tracking on trifecta book odds.** The `placed_trifectas` schema is forward-compatible (snapshots `book_odds` at log time) but no consumer is built in this plan.

## Architecture

```
                                                 ┌─────────────────────────┐
                              already on main →  │ wagerzon_specials       │  (lines posted by WZ)
                                                 └────────────┬────────────┘
                                                              │
                              already on main →  ┌────────────┴─────────────┐
                                                 │ mlb_trifecta_sgp_odds    │  (DK fair, blended at vig 1.25)
                                                 └────────────┬─────────────┘
                                                              │
                       MODIFIED by this plan →    ┌───────────┴────────────┐
                                                  │ Answer Keys/           │
                                                  │   mlb_triple_play.R    │
                                                  │   • adds Kelly         │
                                                  │   • adds trifecta_hash │
                                                  │   • adds dbWriteTable  │
                                                  │   • reads sizing rows  │
                                                  └───────────┬────────────┘
                                                              │ writes
                                NEW   by this plan →  ┌───────┴───────────────────┐
                                                      │ mlb_trifecta_opportunities│  (in mlb.duckdb)
                                                      └───────────┬───────────────┘
                                                                  │ read on each dashboard refresh
                       MODIFIED by this plan →   ┌─────────────── ┴────────────────────┐
                                                 │ MLB Dashboard/mlb_dashboard.R       │
                                                 │   • new load_trifecta_opps()        │
                                                 │   • new create_trifectas_table()    │
                                                 │   • new "Trifectas" tab in HTML     │
                                                 └───────────┬─────────────────────────┘
                                                             │ /api/place-trifecta
                                NEW    by this plan → ┌──────┴───────────┐
                                                      │ placed_trifectas │  (in mlb_dashboard.duckdb)
                                                      └──────────────────┘
                       MODIFIED by this plan →   ┌────────────────────────────────────────────┐
                                                 │ MLB Dashboard/mlb_dashboard_server.py      │
                                                 │   • parallel Popen for triple_play.R       │
                                                 │   • init: placed_trifectas + sizing rows   │
                                                 │   • /api/place-trifecta                    │
                                                 │   • /api/remove-trifecta                   │
                                                 └────────────────────────────────────────────┘
```

### Refresh orchestration (efficiency win)

Today the refresh runs `mlb_correlated_parlay.R` sequentially as Step 2. `mlb_triple_play.R` shares no inputs/outputs with the parlay flow — it reads `wagerzon_specials`, `mlb_consensus_temp`, `mlb_game_samples` (all populated by Step 1) and writes its own table. **Run them in parallel via `subprocess.Popen`:**

```
Step 1   (run.py mlb)              [scrapers + R prepare]
   │
Step 1.5 (parlay_pricer.py mlb)
   │
Step 2  ┬─ Popen: mlb_correlated_parlay.R ─┐
        └─ Popen: mlb_triple_play.R       ─┤  ← NEW. independent, runs concurrently
                                           │
Step 2.5 (parlay_pricer.py mlb --exact-payouts)  ← waits on Step 2
   │
Step 3   (mlb_dashboard.R)         [render HTML]
```

Cost: ~10 lines in `run_pipeline()`. Win: zero added refresh latency — the trifecta DK fetch (5–30s typical) overlaps with correlated-parlay R work.

## Data layer

### `Answer Keys/mlb_triple_play.R` — additions to end of main block

Mirrors the end-of-script pattern in `mlb_correlated_parlay.R`:

1. **Read trifecta sizing settings** from `mlb_dashboard.duckdb::sizing_settings`:
   ```r
   trifecta_bankroll  <- read_setting("trifecta_bankroll",  default = 100)
   trifecta_kelly_mult <- read_setting("trifecta_kelly_mult", default = 0.10)
   trifecta_min_edge  <- read_setting("trifecta_min_edge",  default = 0.05)
   ```
   Dedicated settings (not shared with parlays) because trifecta risk profile is meaningfully different: 3–4 legs at +400 to +800, single-book DK blend, vig pinned at 1.25. `read_setting()` is a small helper that opens the DB read-only with `tryCatch` fallback.

2. **Compute Kelly per row** via `independent_kelly()` (already used by `mlb_correlated_parlay.R` for single-leg parlays). Filter to rows where `edge_pct >= trifecta_min_edge * 100`; non-qualifying rows get `kelly_bet = 0`.

3. **Add `trifecta_hash`** = `sha256(paste(game_id, target_team, prop_type, side, sep="|"))`. Stable across reruns; same hash dedups against `placed_trifectas`.

4. **Add `game` and `game_time`** columns from the `mlb_consensus_temp` join (already in scope — used for the home/away resolution).

5. **`dbWriteTable` step** (matches `mlb_correlated_parlay.R`'s pattern exactly):
   ```r
   write_con <- dbConnect(duckdb(), dbdir = MLB_DB)
   tryCatch({
     dbExecute(write_con, "DROP TABLE IF EXISTS mlb_trifecta_opportunities")
     dbWriteTable(write_con, "mlb_trifecta_opportunities", priced)
     cat(sprintf("Wrote %d trifecta opportunities to %s.\n", nrow(priced), MLB_DB))
   }, error = function(e) {
     cat(sprintf("Warning: Failed to write trifectas to DB: %s\n", e$message))
   })
   if (!is.null(write_con)) dbDisconnect(write_con)
   ```
   Drop+rewrite each run (same as parlay opps) — opportunities are derived data, no need for append history.

6. **Keep the existing console print.** Useful for terminal debugging; doesn't conflict with DB write.

### New table `mlb_trifecta_opportunities` (in `Answer Keys/mlb.duckdb`)

```sql
-- Recreated by mlb_triple_play.R on every dashboard refresh.
CREATE TABLE mlb_trifecta_opportunities (
  trifecta_hash  VARCHAR,    -- sha256(game_id|target_team|prop_type|side)
  game_id        VARCHAR,
  game           VARCHAR,    -- "Away @ Home" display string
  game_time      TIMESTAMP,  -- drives upcoming-only filter on read
  target_team    VARCHAR,    -- canonical team name
  prop_type      VARCHAR,    -- 'TRIPLE-PLAY' or 'GRAND-SLAM'
  side           VARCHAR,    -- 'home' or 'away'
  description    VARCHAR,    -- raw Wagerzon description (e.g. "RANGERS — SCR 1ST, 1H & GM")
  n_samples      INTEGER,
  model_odds     INTEGER,    -- American
  dk_odds        INTEGER,    -- American, NULL when DK didn't price the SGP
  fair_odds      INTEGER,    -- blended American
  book_odds      INTEGER,    -- Wagerzon posted price
  edge_pct       DOUBLE,     -- (fair_prob / book_prob - 1) * 100
  kelly_bet      DOUBLE      -- suggested stake at current bankroll/Kelly mult
);
```

No indexes needed — table never exceeds ~14 rows.

### New table `placed_trifectas` (in `MLB Dashboard/mlb_dashboard.duckdb`)

Created idempotently by `init_db()` in `mlb_dashboard_server.py`:

```sql
CREATE TABLE IF NOT EXISTS placed_trifectas (
  trifecta_hash  TEXT PRIMARY KEY,  -- same hash as in mlb_trifecta_opportunities
  placed_at      TIMESTAMP,
  game_id        TEXT,
  game           TEXT,
  game_time      TIMESTAMP,
  target_team    TEXT,
  prop_type      TEXT,
  side           TEXT,
  description    TEXT,
  book_odds      INTEGER,
  fair_odds      INTEGER,
  edge_pct       DOUBLE,
  kelly_bet      DOUBLE,
  actual_wager   DOUBLE,    -- defaults to round(kelly_bet) on manual log
  status         TEXT       -- 'placed' for manual log; reserved for 'pending'/'price_moved'/etc. for future auto-place
);
```

**Why a dedicated table and not extending `placed_parlays`:** schema cleanliness. `placed_parlays` has ~20 columns specific to combo + spread/total/correlation flow. Trifectas don't share most of those. Forcing them in means NULL columns + branching reads forever. Dedicated table mirrors how parlays themselves are decoupled from the bets tab.

### Sizing settings rows (added to `init_db()` in server)

```python
# After the existing parlay_bankroll / parlay_kelly_mult / parlay_min_edge inserts:
con.execute("INSERT OR IGNORE INTO sizing_settings (param, value) VALUES ('trifecta_bankroll',  100)")
con.execute("INSERT OR IGNORE INTO sizing_settings (param, value) VALUES ('trifecta_kelly_mult', 0.10)")
con.execute("INSERT OR IGNORE INTO sizing_settings (param, value) VALUES ('trifecta_min_edge',  0.05)")
```

Defaults rationale: 10% Kelly mult (half the parlay default of 25%) starts conservative until the 1.25 vig is validated against real fills. 5% min-edge (vs. 0% for parlays) demands more cushion given trifectas' higher payout volatility.

## UI / interaction layer

### Dashboard R (`Answer Keys/MLB Dashboard/mlb_dashboard.R`)

**New helpers** (mirroring `parlay_opps` load):
- `load_trifecta_opps(mlb_db)` — reads `mlb_trifecta_opportunities` filtered to `game_time IS NULL OR game_time > NOW()`. Defensive `tryCatch` returns empty tibble if table missing.
- `load_placed_trifectas(dash_db)` — reads `placed_trifectas` filtered the same way. Defensive `tryCatch`.

**New table renderer** `create_trifectas_table(trifecta_opps, placed_trifectas, trifecta_bankroll, trifecta_kelly_mult)` — reactable mirroring `create_bets_table()`'s ergonomics (NOT `create_parlays_table()`'s):
- `searchable = TRUE`, `filterable = TRUE`, `striped = TRUE`, `compact = TRUE`, `highlight = TRUE`
- `defaultPageSize = 25` (covers all rows; pagination chrome hidden in practice)
- Color-coded `edge_pct` cell, same gradient as bets tab EV: ≥15% bold green, ≥10% mid-green, ≥5% light-green, otherwise default
- Hidden helper columns (`trifecta_hash`, `game_time`, `game_id`) — used by JS via `data-*` attributes on the Action button
- Visible columns:

  | Column      | Source              | Notes                                          |
  |-------------|---------------------|------------------------------------------------|
  | Game        | `game`              | filterable                                     |
  | Time        | `game_time`         | client-side `toLocaleString` like bets tab     |
  | Pick        | `target_team` + `side` | bold team name + grey side hint             |
  | Prop        | `prop_type`         | filterable, narrow column                      |
  | Description | `description`       | raw WZ string, e.g. "SCR 1ST, 1H & GM"        |
  | Model       | `model_odds`        | monospace                                      |
  | DK          | `dk_odds`           | monospace, dash if NULL                        |
  | Fair        | `fair_odds`         | monospace, bold (this is the published number) |
  | Book        | `book_odds`         | monospace                                      |
  | Edge        | `edge_pct`          | color-coded gradient                           |
  | Stake       | `kelly_bet`         | rounded to nearest dollar                      |
  | Action      | conditional         | Place / Placed button                          |

**Action column conditional** (mirrors bets tab `placeBet` / `removeBet` pattern, NOT the parlay-tab fragment pattern):

```r
# Pseudocode for cell renderer
if (trifecta_hash %in% placed_trifectas$trifecta_hash) {
  # button class "btn-placed", onclick="removeTrifecta(this)"
} else if (kelly_bet > 0) {
  # button class "btn-place",  onclick="placeTrifecta(this)"
  # All data needed for the POST goes in data-* attrs
} else {
  # empty cell — row stays visible for context (model_odds, dk_odds, fair_odds)
  # but no Place button when edge_pct < trifecta_min_edge (kelly_bet == 0)
}
```

Below-threshold rows stay in the table (so the user can see what was priced and rejected) but show no Place button — same idea as bets tab filtering out below-threshold bets, but trifectas have so few rows that hiding them costs more than it gains.

The button has `data-trifecta-hash`, `data-actual-wager`, `data-game`, etc. as attributes. Clicking it fires a `fetch('/api/place-trifecta', ...)` and on success swaps the button class + onclick in place. No DOM rebuild, no fragment endpoint, no full-table refresh.

**Tab insertion in the rendered HTML.** The existing tabs are hand-rolled HTML in `mlb_dashboard.R`. Add:
- `<button data-tab="trifectas">Trifectas</button>` to the tab strip, positioned right after the "Parlays" button so it appears next-to-correlated-parlays as requested
- `<div id="tab-trifectas" class="tab-pane">…content with reactable widget…</div>`

The existing tab-switching JS (already handles `data-tab="bets"`, `data-tab="parlays"`, etc.) picks up the new tab automatically — no JS changes needed for tab navigation.

### Server (`Answer Keys/MLB Dashboard/mlb_dashboard_server.py`)

**`init_db()` additions:**
- `CREATE TABLE IF NOT EXISTS placed_trifectas` (DDL above)
- Three `INSERT OR IGNORE` rows into `sizing_settings`

**Refresh orchestration change** in `run_pipeline()` Step 2: replace the single `subprocess.run(...)` for `mlb_correlated_parlay.R` with parallel `subprocess.Popen` for both R scripts, then wait on both:

```python
# Step 2: Find correlated parlay opportunities AND price trifectas (parallel, non-fatal)
print("Finding parlay opportunities + pricing trifectas...")
parlay_proc = subprocess.Popen(
    ["Rscript", str(answer_keys_dir / "mlb_correlated_parlay.R")],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    cwd=str(nfl_work_dir)
)
trifecta_proc = subprocess.Popen(
    ["Rscript", str(answer_keys_dir / "mlb_triple_play.R")],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    cwd=str(nfl_work_dir)
)
parlay_out, parlay_err = parlay_proc.communicate()
trifecta_out, trifecta_err = trifecta_proc.communicate()
if parlay_proc.returncode != 0:
    print(f"Parlay finder warning: {parlay_err.decode('utf-8', errors='replace')[-300:]}")
if trifecta_proc.returncode != 0:
    print(f"Trifecta pricer warning: {trifecta_err.decode('utf-8', errors='replace')[-300:]}")
```

**New endpoints** (~50 lines total, mirroring `/api/place-bet` + `/api/remove-bet`):

- `POST /api/place-trifecta` — body: `{trifecta_hash, actual_wager}`. Server-side fetches the full opportunity row from `mlb_trifecta_opportunities` by hash and INSERTs into `placed_trifectas` (so client only needs to send hash + wager). Idempotent on `trifecta_hash` PK — second click is a no-op (`INSERT OR IGNORE`). Returns `{ok: true}`.
- `POST /api/remove-trifecta` — body: `{trifecta_hash}`. DELETEs from `placed_trifectas`. Returns `{ok: true}`.

No `/api/trifecta-table-fragment` endpoint — explicitly skipped per efficiency design.

### Client-side JS (in the rendered HTML)

Two small functions added to the existing `<script>` block in `mlb_dashboard.R`:

```javascript
async function placeTrifecta(btn) {
  const hash = btn.dataset.trifectaHash;
  const wager = parseFloat(btn.dataset.actualWager);
  const r = await fetch('/api/place-trifecta', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({trifecta_hash: hash, actual_wager: wager})
  });
  if (r.ok) {
    btn.textContent = 'Placed';
    btn.className = 'btn-placed';
    btn.onclick = () => removeTrifecta(btn);
  }
}

async function removeTrifecta(btn) {
  const hash = btn.dataset.trifectaHash;
  const r = await fetch('/api/remove-trifecta', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({trifecta_hash: hash})
  });
  if (r.ok) {
    btn.textContent = 'Place';
    btn.className = 'btn-place';
    btn.onclick = () => placeTrifecta(btn);
  }
}
```

## Efficiency summary

| Concern                       | Cost                                                       | Mitigation                                              |
|-------------------------------|------------------------------------------------------------|---------------------------------------------------------|
| Refresh latency               | DK scraper takes 5–30s                                     | Run parallel with correlated parlay R via `Popen`       |
| Place click → UI update       | Parlay tab does full table re-render (fragment endpoint)   | In-place button class swap; no fragment endpoint        |
| DuckDB lock contention        | Parlay R + triple_play R both write to `mlb.duckdb` concurrently | Each writes its own table; DuckDB serializes write transactions transparently. Both writes are `DROP + dbWriteTable` of small tables (~14 trifecta rows, ~50 parlay rows). Microsecond block at most. |
| HTML render time              | ~14 trifecta rows + existing tabs                          | Trivial — reactable handles 14 rows in microseconds     |
| Settings panel weight         | Three new sliders                                          | Same row pattern as parlay sliders; no new UI primitives|
| Endpoint round-trip on place  | Single INSERT, single hash lookup                          | Server reads opportunity row by indexed hash, single statement to placed_trifectas |

## Error handling matrix

| Failure mode                                   | Behavior                                                                                                |
|------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| `mlb_trifecta_opportunities` table missing     | `load_trifecta_opps()` returns empty tibble; tab renders "No trifectas priced yet" placeholder          |
| `placed_trifectas` table missing               | `load_placed_trifectas()` returns empty tibble; every row shows Place button                            |
| `mlb_triple_play.R` exits nonzero on refresh   | Server logs warning (non-fatal, mirrors parlay step). Last successful trifecta opportunities still served. |
| DK scraper inside pricer fails                 | Already handled by Plan #2 — pricer continues with model-only blend; `dk_odds` is NULL                  |
| Wagerzon Specials empty (no posted props)      | Pricer prints "No priceable specials found" and exits 0; opportunities table gets dropped + empty       |
| `/api/place-trifecta` called with unknown hash | Server returns 404; client leaves button unchanged                                                      |
| Double-click on Place                          | Second INSERT is no-op via `INSERT OR IGNORE` on PK; client's button is already in `Placed` state       |
| User edits `actual_wager` post-placement       | Out of scope; future plan can add an Edit/Adjust button mirroring bets tab's `updateBet` pattern        |

## Testing summary

| Layer       | Test                                                                                                |
|-------------|-----------------------------------------------------------------------------------------------------|
| R unit      | `test_triple_play.R` already covers the pricer math; extend it with one test confirming `trifecta_hash` is stable across reruns and `kelly_bet` reads the new sizing rows. |
| R end-to-end| Run `Rscript mlb_triple_play.R` against a copied `mlb.duckdb`; confirm `mlb_trifecta_opportunities` table exists with expected columns + row count matches the console print. |
| Server unit | One pytest covering `/api/place-trifecta` (idempotent on duplicate hash) and `/api/remove-trifecta`. |
| Manual UX   | Refresh the dashboard, switch to Trifectas tab, sort by Edge, click Place on top row, verify button toggles to Placed, refresh page and confirm placement persists. |
| Regression  | Confirm Bets and Parlays tabs unchanged after merging — placement flows for both untouched.         |

## Version control plan

- **Branch:** `feature/mlb-trifecta-dashboard-tab` (already created)
- **Worktree:** `.worktrees/mlb-trifecta-dashboard` (already created)
- **Source:** branched from `main` at `51e721b` (post Plan #1 + Plan #2 merges)
- **Approximate commits** (one per logical task):
  1. `feat(mlb): write priced trifectas to mlb_trifecta_opportunities` — pricer dbWriteTable + trifecta_hash + Kelly + sizing read
  2. `feat(dashboard): create placed_trifectas table + sizing rows on init` — `init_db()` migration
  3. `feat(dashboard): /api/place-trifecta + /api/remove-trifecta endpoints` — server endpoints
  4. `perf(dashboard): run mlb_correlated_parlay.R + mlb_triple_play.R in parallel` — Popen orchestration change
  5. `feat(dashboard): Trifectas tab — load, render, place/remove JS` — dashboard R + JS
  6. `docs: trifecta dashboard tab in Answer Keys CLAUDE.md + MLB Dashboard README`
- **DuckDB handling:** copy `mlb.duckdb` and `mlb_dashboard.duckdb` into worktree only for end-to-end verification; remove before commit (per project rule — never check in DBs, never symlink them).
- **Pre-merge:** `git diff main..HEAD --stat`, run pre-merge review checklist (below), get explicit user approval before merging to `main`.
- **Cleanup post-merge:** `git worktree remove .worktrees/mlb-trifecta-dashboard && git branch -d feature/mlb-trifecta-dashboard-tab`.

## Documentation plan

- **`Answer Keys/CLAUDE.md`** — extend the "Triple-Play Data Flow" section with a note that the pricer now writes `mlb_trifecta_opportunities` and the dashboard renders a Trifectas tab. ~5 lines.
- **`Answer Keys/MLB Dashboard/README.md`** — add a "Trifectas tab" subsection: data source (`mlb_trifecta_opportunities`), refresh wiring (parallel with parlay), placement flow (manual log into `placed_trifectas`), tab navigation, sizing settings (`trifecta_*` rows). ~20 lines.
- Updates land in the same commit as the corresponding code change (per project doc-discipline rule).

## Pre-merge review checklist

Run before requesting merge approval (per project CLAUDE.md):

- **Data integrity:** `mlb_trifecta_opportunities` is `DROP + dbWriteTable` per refresh — no duplicate rows. `placed_trifectas` PK on `trifecta_hash` prevents double-log. Hash includes only stable identifiers (game_id, team, prop_type, side) — does NOT include `book_odds` or `fair_odds`, so the same posted prop dedups across price moves.
- **Resource safety:** Every `dbConnect` paired with `dbDisconnect` (via `on.exit` or explicit close). No `DROP TABLE` runs on a connection still held by another process — write happens after pricer's other DB reads complete.
- **Edge cases verified:** off-day (zero specials posted), first run (no `placed_trifectas` table yet), DK scraper failure (model-only fair, NULL `dk_odds` column rendered as dash), placed game now in past (filtered out by `game_time` check on read).
- **Dead code:** No unused imports/helpers introduced. `load_placed_trifectas` is referenced by the new render path; `placeTrifecta` / `removeTrifecta` JS are referenced by the inline button onclick.
- **Log/disk hygiene:** Pricer's existing console output unchanged; no new long-lived log files introduced.
- **Security:** `/api/place-trifecta` and `/api/remove-trifecta` validate `trifecta_hash` is present and non-empty; the hash is opaque to clients (sha256). No user-supplied SQL.
- **Concurrency:** Parallel `Popen` of two writers to `mlb.duckdb` is safe (DuckDB serializes write transactions). Refresh is already mutex'd via `_refresh_lock`.

## Open items (none blocking)

1. **Edit/Adjust placed wager.** Bets tab has `updateBet` for partial fills. Trifectas don't get partial fills typically, but the user may want to log a different stake than `kelly_bet` recommends. Default behavior in this plan: `actual_wager = round(kelly_bet)` on initial place; an Edit button is a follow-up.
2. **Auto-placement to Wagerzon Specials.** Future plan; out of scope here. Will need its own recon of the Wagerzon Specials placement endpoint.
3. **Bankroll percentage sizing vs. flat Kelly.** Current `kelly_bet` is dollars at the slider's bankroll value. Some users prefer "always X% Kelly of current liquidity"; future enhancement.
4. **CLV tracking.** `placed_trifectas` snapshots `book_odds` + `fair_odds` at log time — schema-ready, no consumer yet.
