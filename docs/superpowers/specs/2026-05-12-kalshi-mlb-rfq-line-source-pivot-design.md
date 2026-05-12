# Kalshi MLB RFQ — Line Source Pivot Design

**Date:** 2026-05-12
**Status:** Design — awaiting plan
**Branch:** `worktree-kalshi-mlb-rfq-line-source-pivot`
**Predecessors:**
- `docs/2026-05-12-kalshi-mlb-rfq-architecture-pivot.md` — chat summary that motivated this design
- `docs/superpowers/specs/2026-04-27-kalshi-mlb-rfq-bot-design.md` — original bot spec

---

## TL;DR

Decouple the Kalshi MLB RFQ bot from Wagerzon-derived line selection. Bot drives its own line surface from Kalshi cross-category MVE enumeration. Library-refactor the four SGP scrapers (DK/FD/PX/NV) into per-book Python modules. Bot owns a sibling market DB. Bot only bets where ≥2 books priced the same line.

## Goal

Move the bot's edge surface from "combos Wagerzon happens to post" to "every Kalshi MVE combo with cross-book consensus." Eliminate the dashboard's `mlb_parlay_lines` as a hidden gating layer on bot behavior.

---

## Context

The four MLB SGP scrapers (`mlb_sgp/scraper_{draftkings,fanduel,prophetx,novig}_sgp.py`) all read `mlb_parlay_lines` from `MLB_DB` at startup to know which `(game_id, fg_spread, fg_total)` tuples to price. `mlb_parlay_lines` is written by `Answer Keys/mlb_correlated_parlay.R:407-410` from `wz_fg_matched` — Wagerzon's posted FG spread+total combinations.

So the bot's `_SGP_ODDS_CACHE` only covers combos Wagerzon happens to post. Smoke testing (2026-05-12) confirmed: bot-owned DB with empty `mlb_parlay_lines` → all four scrapers exit rc=0 having scraped nothing → bot's downstream pipeline silently degrades.

A prior plan tried to ship an SGP cadence loop over the same line source; the cadence work is orthogonal but pointless without solving line selection first. That plan was discarded; design fragments from it are reusable but re-derived cleanly here.

---

## Architecture

Three layers, clean separation:

```
┌─────────────────────────────────────────────────┐
│ Layer 1: Per-book library modules               │
│   mlb_sgp/draftkings.py                         │
│   mlb_sgp/fanduel.py                            │
│   mlb_sgp/prophetx.py                           │
│   mlb_sgp/novig.py                              │
│                                                  │
│   Each exposes: price_sgps(target_lines,        │
│                            periods,             │
│                            session=None)        │
│                 → list[PricedRow]               │
│                                                  │
│   Pure: no DB I/O, no env vars, no subprocess.  │
│   Only book-API logic.                          │
└─────────────────────────────────────────────────┘
                       ▲
                       │ imports
        ┌──────────────┴──────────────┐
        │                             │
┌───────────────────┐         ┌───────────────────┐
│ Layer 2a:         │         │ Layer 2b:         │
│ Dashboard caller  │         │ Bot caller        │
│                   │         │                   │
│ mlb_sgp/          │         │ kalshi_mlb_rfq/   │
│   scraper_X_sgp.py│         │   sgp_runner.py   │
│ (thin shims,      │         │                   │
│ subprocessed by R)│         │ Reads             │
│                   │         │ mlb_target_lines, │
│ Reads             │         │ calls libs FG-    │
│ mlb_parlay_lines, │         │ only, writes new- │
│ calls libs FG+F5, │         │ schema rows to    │
│ writes mlb_mm DB  │         │ bot DB            │
└───────────────────┘         └───────────────────┘
```

**Data flow on a bot SGP cadence tick (atomic, serial):**

1. Bot enumerates open Kalshi MVE markets per game (KXMLBSPREAD-*, KXMLBTOTAL-*).
2. Bot rewrites `mlb_target_lines` table in `kalshi_mlb_rfq_market.duckdb` (single transaction; DELETE + INSERT).
3. Bot spawns the four scraper subprocesses with `MLB_SGP_DB_PATH=kalshi_mlb_rfq_market.duckdb`, `MLB_SGP_PERIODS=FG`, global deadline, kill-on-hang.
4. Scrapers read `mlb_target_lines`, price every `(game, FG, spread, total)` tuple via their library module, write `mlb_sgp_odds` rows tagged with `spread_line` and `total_line`.
5. Bot reloads `_SGP_ODDS_CACHE` and `_PARLAY_LINES_CACHE` from the bot DB.

No concurrent writers to `mlb_target_lines`. SGP cadence and RFQ refresh share enumeration data via the cache, not via redundant Kalshi API calls.

**Dashboard data flow unchanged.** The dashboard's `mlb_correlated_parlay.R` still writes its single-line-per-game `mlb_parlay_lines` to `mlb_mm.duckdb`, subprocess-invokes `scraper_X_sgp.py` (now thinner), and reads `mlb_sgp_odds` back. The only dashboard-visible delta: `mlb_sgp_odds` rows now carry `spread_line`/`total_line` columns (nullable for pre-existing rows; populated for new rows).

---

## Library API

```python
# mlb_sgp/_shared.py
from dataclasses import dataclass
from datetime import datetime
from typing import Literal

@dataclass(frozen=True)
class TargetLine:
    game_id: str
    home_team: str
    away_team: str
    commence_time: datetime
    period: Literal["FG", "F5"]
    spread: float   # signed, home-perspective (e.g. -1.5 = home favored)
    total: float

@dataclass(frozen=True)
class PricedRow:
    game_id: str
    combo: str              # "Home Spread + Over", "Away Spread + Under", etc.
    period: str             # "FG" or "F5"
    spread_line: float
    total_line: float
    bookmaker: str          # "draftkings", "fanduel", "prophetx", "novig"
    source: str             # "draftkings_direct", etc.
    sgp_decimal: float
    sgp_american: int
    fetch_time: datetime

# mlb_sgp/draftkings.py  (and 3 siblings, same signature)
def price_sgps(
    target_lines: list[TargetLine],
    periods: list[str] = ("FG",),
    session: Session | None = None,
    verbose: bool = False,
) -> list[PricedRow]: ...
```

- `target_lines` is multi-row by design. Dashboard passes `[FG, F5] × N games`; bot passes `[FG, FG, ...]` (many FG lines per game).
- `periods` filters which subset of `target_lines` actually get priced. Bot sets `("FG",)`; dashboard sets `("FG", "F5")`.
- `session` is optional so callers can reuse a session across invocations (rare, but useful for testing).
- No DB I/O inside the library. Caller is responsible for loading inputs and writing outputs.

---

## File map (deltas only)

### New

| Path | Purpose |
|---|---|
| `mlb_sgp/_shared.py` | `TargetLine`, `PricedRow`, `_utc_bucket`, `decimal_to_american`/inverse, vig-sanity helpers, `load_target_lines()` (handles both table names) |
| `mlb_sgp/draftkings.py` | Pure DK book-API module exposing `price_sgps()`. Auth/session, market fetching, selection-ID resolution, batch SGP pricing call. |
| `mlb_sgp/fanduel.py` | Same for FanDuel (`implyBets` endpoint, sid encoding). |
| `mlb_sgp/prophetx.py` | Same for ProphetX (RFQ endpoint, profile-cookie auth, MIN_OFFER_STAKE filter, F5-Over SANITY_MULT_RATIO defense). |
| `mlb_sgp/novig.py` | Same for Novig (anonymous GraphQL, Hasura market tree). |
| `kalshi_mlb_rfq/sgp_runner.py` | Bot caller: `enumerate_kalshi_targets()`, `write_target_lines()`, `run_scrapers(env=...)`, `read_priced_rows()`. Subsumes the existing `should_scrape` helper (already committed). |
| `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb` | Bot's market data DB (sibling to `kalshi_mlb_rfq.duckdb`). Holds `mlb_target_lines` + `mlb_sgp_odds`. Created at startup if missing. |
| Tests: `mlb_sgp/tests/test_shared.py`, `kalshi_mlb_rfq/tests/test_sgp_runner.py` (extend existing) | Unit tests for shared utilities and bot caller. |

### Rewritten thin (preserved subprocess invocation contract)

| Path | New shape |
|---|---|
| `mlb_sgp/scraper_draftkings_sgp.py` | ~60-80 line shim: parse env, open MLB_DB, `load_target_lines()`, `draftkings.price_sgps(target_lines, periods=["FG","F5"])`, `db.upsert_priced_rows()`, exit. R pipeline invocation unchanged. |
| `mlb_sgp/scraper_fanduel_sgp.py` | Same shape, calls `fanduel.price_sgps`. |
| `mlb_sgp/scraper_prophetx_sgp.py` | Same shape, calls `prophetx.price_sgps`. |
| `mlb_sgp/scraper_novig_sgp.py` | Same shape, calls `novig.price_sgps`. |

### Modified

| Path | Changes |
|---|---|
| `mlb_sgp/db.py` | `MLB_DB` env-var override (`MLB_SGP_DB_PATH`); `ensure_table()` runs `CREATE TABLE IF NOT EXISTS` with new schema + `ALTER TABLE ... ADD COLUMN IF NOT EXISTS spread_line DOUBLE` / `total_line DOUBLE` (idempotent migration); new `upsert_priced_rows(rows: list[PricedRow])` function. |
| `kalshi_mlb_rfq/config.py` | Worktree path fix: `.worktrees` detection matching `mlb_sgp/db.py:19-20`; new env knobs (`SGP_REFRESH_SEC`, `SGP_SCRAPER_TIMEOUT_SEC`, `BOT_MARKET_DB`, `MIN_BOOK_COUNT_FOR_BLEND`). Remove `MIN_EV_PCT_MODEL_ONLY` (decision below). |
| `kalshi_mlb_rfq/main.py` | `_PARLAY_LINES_CACHE` shape change (multi-line per game); `_load_book_fairs` rewires to per-line lookup + N≥2 hard gate; `_current_book_lines_for_combo` removed (dead code post-pivot); `line_move_ok` per-accept gate removed; `_refresh_sgp_cache()` partial reload; SGP cadence timer added to `main_loop`; synchronous warm-up before entering loop. |

### Untouched

- `Answer Keys/mlb_correlated_parlay.R` — keeps writing legacy `mlb_parlay_lines` to `mlb_mm.duckdb`. Dashboard pipeline unchanged.
- `kalshi_mlb_rfq/combo_enumerator.py` — Kalshi MVE enumeration already in place.
- `kalshi_mlb_rfq/auth_client.py`, `rfq_client.py`, `kelly.py`, `risk.py` (mostly), `fair_value.py`, `db.py` — no changes from this pivot.

---

## Schema changes

### `mlb_sgp_odds` — two new nullable columns

```sql
CREATE TABLE IF NOT EXISTS mlb_sgp_odds (
    game_id       VARCHAR,
    combo         VARCHAR,
    period        VARCHAR,
    bookmaker     VARCHAR,
    sgp_decimal   DOUBLE,
    sgp_american  INTEGER,
    fetch_time    TIMESTAMP,
    source        VARCHAR,
    spread_line   DOUBLE,    -- NEW (nullable for legacy rows)
    total_line    DOUBLE     -- NEW (nullable for legacy rows)
);
```

Migration runs via `mlb_sgp/db.py::ensure_table()` on first scraper invocation after deployment:

```sql
ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS spread_line DOUBLE;
ALTER TABLE mlb_sgp_odds ADD COLUMN IF NOT EXISTS total_line DOUBLE;
```

**Backward compatibility verified:** `Answer Keys/mlb_correlated_parlay.R:521-525` reads `mlb_sgp_odds` with `SELECT game_id, combo, sgp_decimal, sgp_american, source` — named columns, not positional. New columns are ignored. No R code changes required.

### `mlb_target_lines` — new (bot DB only)

```sql
CREATE TABLE IF NOT EXISTS mlb_target_lines (
    game_id        VARCHAR,
    home_team      VARCHAR,
    away_team      VARCHAR,
    commence_time  TIMESTAMP,
    period         VARCHAR,    -- "FG" only for bot; future "F5" support for dashboard would work here
    spread         DOUBLE,
    total          DOUBLE,
    written_at     TIMESTAMP
);
```

Multi-row per game. Rewritten on each SGP cadence tick (DELETE all + INSERT).

### `kalshi_mlb_rfq_market.duckdb` — new DB file

Sibling to `kalshi_mlb_rfq/kalshi_mlb_rfq.duckdb`. Holds only `mlb_target_lines` and `mlb_sgp_odds`. No state tables (live_rfqs, fills, etc. — those stay in the state DB).

---

## Decisions

### D1 — All Kalshi-open lines, FG only

Bot enumerates every open Kalshi MVE `(spread, total)` tuple per game and writes one `mlb_target_lines` row per tuple. No canonical-line heuristic. Scrapers iterate every row. **F5 is out of scope** for the bot — Kalshi doesn't list F5 MVE markets, so F5 isn't quotable. Dashboard continues to price FG+F5 in its own pipeline.

### D2 — Sibling market DB

`kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb` sits next to the state DB. Scrapers write to it; bot reads from it. Avoids lock contention with frequent state writes (RFQs, fills, quote_log) in `kalshi_mlb_rfq.duckdb`.

### D3 — Schedule source

Bot reads `mlb_odds_temp` directly from `Answer Keys/mlb.duckdb` for `(game_id, home_team, away_team, commence_time)` mapping. Read-only. Zero R-pipeline changes. The bot's `_PARLAY_LINES_CACHE` is populated from two sources merged: `mlb_odds_temp` (schedule) + bot DB `mlb_target_lines` (lines).

### D4 — Library refactor (Option 3)

Each book becomes an importable Python module with no DB I/O. Existing `scraper_X_sgp.py` files become thin subprocess shims that call the library. Dashboard's R subprocess invocation contract is preserved.

**v1 benefit: code organization + testability.** Runtime speed gain (in-process scraping from bot, avoiding subprocess + DuckDB lock from child processes) is a v2 follow-up and explicitly NOT shipped here. The library refactor is what enables v2.

### D5 — N≥2 books required to bet (hardened)

Bot only bets candidates where ≥2 books priced the matching `(spread, total)`. If fewer than 2 books are available, the candidate is dropped before quote evaluation.

- `_load_book_fairs` returns `{book: fair}` dict only if `len(books) >= 2`, else returns `{}`.
- When `book_fairs` is empty, `fair_value.blend` returns None or model-only — and the candidate is filtered upstream in `_enumerate_and_score_all_games`.
- No `MIN_EV_PCT_MODEL_ONLY` config knob. No model-only bets.

Quant rationale: a single book's devigged fair is one opinion, not a consensus. The bot's edge comes from cross-book cross-model agreement, not from being smarter than one book.

### D6 — Cold start: synchronous SGP warm-up

On bot startup, before entering `main_loop`:

1. `db.init_database()`
2. `_refresh_caches()` (samples + schedule from `mlb_odds_temp`)
3. `_phantom_rfq_cleanup()`
4. **`_sgp_cycle()` synchronous** — runs target_line enumeration, scraper bundle, cache reload. Blocks ~60-90s.
5. Enter `main_loop`.

First RFQ refresh happens with a warm `_SGP_ODDS_CACHE`. No main-loop branches for "not-yet-warmed."

### D7 — Drop `line_move_ok` gate

The line-move per-accept gate becomes tautological post-pivot. Each candidate's `(spread, total)` is encoded in its Kalshi market ticker — the line can't move per-candidate. If Kalshi closes one market and opens another, that's a new candidate, not a line move.

Drift protection is already provided naturally: RFQ refresh re-scores every 30s and cancels out-of-top-N candidates. Blended fair drift → candidate drops out of top-N → RFQ cancelled.

**Action:** Remove `line_move_ok` call from `_all_per_accept_gates_pass`. Stop writing to `reference_lines` at RFQ submission. Schema can stay (unused). `risk.py::line_move_ok` function can stay for potential future use or be removed in a follow-up cleanup.

### D8 — `_current_book_lines_for_combo` removed

Dead code post-pivot. Was used by `line_move_ok` (D7) and by the reference_lines snapshot. Both removed.

### D9 — Worktree path fix in `config.py`

`kalshi_mlb_rfq/config.py:7`'s `PROJECT_ROOT = PKG_DIR.parent` doesn't account for worktrees. Match the `mlb_sgp/db.py:19-20` pattern:

```python
_REPO_ROOT = PKG_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))
PROJECT_ROOT = _REPO_ROOT
```

Bundled into this work because smoke tests need it.

### D10 — Coordination between RFQ refresh and SGP cadence

Two timers in `main_loop`. They share Kalshi MVE enumeration data via a cached snapshot so we don't re-enumerate twice:

- SGP cadence tick (every `SGP_REFRESH_SEC`, default 60s): re-enumerates Kalshi MVE → updates a module-level cache → writes `mlb_target_lines` → spawns scrapers → waits → reloads `_SGP_ODDS_CACHE`.
- RFQ refresh tick (every `RFQ_REFRESH_SEC`, default 30s): reads from the cached MVE enumeration (no API calls), scores candidates against `_SGP_ODDS_CACHE`, manages live RFQs.

If RFQ refresh fires before the first SGP cadence completes, it has no MVE data → skip. Post-warmup, this never happens.

---

## Cache & gate integration

### `_PARLAY_LINES_CACHE` — shape change

Today: `{game_id: {home_team, away_team, commence_time, fg_spread, fg_total}}`

After: `{game_id: {home_team, away_team, commence_time, fg_lines: [(spread, total), ...]}}`

Loader is two-step:
1. Read `(game_id → home_team, away_team, commence_time)` from `mlb_odds_temp` in `mlb.duckdb`.
2. Read `(game_id → [(spread, total)])` from `mlb_target_lines` in bot DB, group by game_id.
3. Merge into the cache dict. Drop games with no target lines (Kalshi has no open MVE).

### `_load_book_fairs(game_id, spread, total)` — rewires

```python
def _load_book_fairs(game_id, spread_line, total_line):
    """Per-line lookup of book-blended fair. N≥2 books required."""
    if _SGP_ODDS_CACHE is None or _SGP_ODDS_CACHE.empty:
        return {}
    rows = _SGP_ODDS_CACHE[
        (_SGP_ODDS_CACHE["game_id"] == game_id) &
        (_SGP_ODDS_CACHE["spread_line"] == spread_line) &
        (_SGP_ODDS_CACHE["total_line"] == total_line)
    ]
    if rows.empty:
        return {}
    out = {}
    for book in rows["bookmaker"].unique():
        sub = rows[rows["bookmaker"] == book]
        fair_per_book = fair_value.devig_book(
            sub, combo="Home Spread + Over",  # combo selection from candidate's legs
            vig_fallback=_vig_fallback(book),
        )
        if fair_per_book is not None:
            out[book] = fair_per_book
    if len(out) < config.MIN_BOOK_COUNT_FOR_BLEND:  # default 2
        return {}
    return out
```

### Per-accept gates after pivot

| Gate | Status |
|---|---|
| `staleness_ok` (prediction) | unchanged |
| `tipoff_ok` | unchanged |
| `line_move_ok` | **REMOVED** (D7) |
| `fair_in_bounds` | unchanged |
| `kill_switch_ok` | unchanged |
| `cooldown_ok` | unchanged |
| `inverse_combo_ok` | unchanged |
| `per_game_cap_ok` | unchanged |
| `daily_cap_ok` | unchanged |
| `fill_ratio_ok` | unchanged |
| Positions-API health | unchanged |
| **Book count ≥ 2** | **NEW upstream filter** (drops candidate before per-accept gates) |

---

## Env vars / config

| Var | Default | Purpose |
|---|---|---|
| `MLB_SGP_DB_PATH` | unset (= `mlb_mm.duckdb`) | Scraper write target override. Bot sets this to bot DB path. |
| `MLB_SGP_PERIODS` | `FG,F5` (dashboard) | Periods to price. Bot sets `FG`. |
| `SGP_REFRESH_SEC` | `60` | Bot SGP cadence tick interval. |
| `SGP_SCRAPER_TIMEOUT_SEC` | `90` | Per-scraper deadline. Killed and skipped if exceeded. |
| `BOT_MARKET_DB` | `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb` | Bot market DB path. |
| `MIN_BOOK_COUNT_FOR_BLEND` | `2` | Min books required for a candidate to proceed. |

---

## Non-goals

- **F5 quoting for the bot.** Kalshi doesn't have F5 MVE markets; no point pricing F5 from the bot side. Dashboard pipeline continues to price F5 unchanged.
- **In-process scraper calls from bot.** v1 is subprocess (process isolation; preserves the existing pattern, scraper crash doesn't kill bot). The library refactor enables this for v2.
- **Off-line book coverage betting.** N≥2 books required hard gate. Lines where only one book has priced are not bet.
- **Per-game cooldown layer.** Combo-level cooldown stays as-is. Per-game cooldown is a potential follow-up if observed needed.
- **`MAX_LIVE_RFQS` bump.** Stays at 80; revisit if consistently top-N-bound after the wider candidate pool ships.
- **Dashboard line selection changes.** `mlb_correlated_parlay.R` keeps its Wagerzon-matched line picker. Untouched.
- **Per-book unit tests with fixtures.** Existing scrapers have no unit tests. Library refactor introduces tests for shared utilities and bot caller; per-book modules tested by smoke run + dashboard regression diff.

---

## Test plan

1. **Shared library unit tests** — `mlb_sgp/tests/test_shared.py`: `TargetLine`/`PricedRow` round-trips, `load_target_lines()` dual-table-name logic with both schemas, F5 NULL handling preserves FG-only emission, vig helpers.
2. **Bot caller unit tests** — `kalshi_mlb_rfq/tests/test_sgp_runner.py` (extend existing): `enumerate_kalshi_targets()` mocks Kalshi API and verifies multi-row output, `write_target_lines()` is atomic (DELETE + INSERT in one transaction), `run_scrapers(env=...)` subprocess plumbing handles timeout/kill/log-handle cleanup, `read_priced_rows()` filters by staleness and returns expected shape.
3. **`main.py` integration tests** — synchronous warm-up path runs without error; `_load_book_fairs` returns `{}` when book count < 2; `_PARLAY_LINES_CACHE` shape change doesn't break `_resolve_game_id`/`_commence_time_for_game`; SGP cadence timer fires on the right interval.
4. **Dashboard regression** — run thin-shim `scraper_draftkings_sgp.py` (and 3 siblings) against `mlb_mm.duckdb` with the pre-refactor `mlb_parlay_lines` data. Diff row count + per-game pricing vs. a pre-refactor baseline snapshot. Acceptance: ≤2% row-count drift, no semantic differences in priced combos.
5. **Bot dry-run smoke** — boot bot with `--dry-run`, observe:
   - `mlb_target_lines` populated with multi-row Kalshi MVE combos per game
   - All 4 scrapers write non-zero `mlb_sgp_odds` rows to bot DB within `SGP_SCRAPER_TIMEOUT_SEC`
   - `_SGP_ODDS_CACHE` non-empty after first cycle
   - `quote_log` decisions: no `declined_stale_predictions` dominance, no `no_fresh_fair` dominance
   - Candidates without N≥2 books are filtered (verifiable via decision distribution)
6. **30-min dry-run on live MVE** — verify per-book consensus on the wider line surface (off-canonical lines now have book fairs where ≥2 books priced), sanity-check fill candidates would have hit (post-fee EV ≥ MIN_EV_PCT).

---

## Open follow-ups (NOT in this scope)

| Item | Trigger |
|---|---|
| In-process scraping from bot | If subprocess overhead bottlenecks SGP cadence below 60s. |
| Per-game cooldown layer | If observed candidate flurries on single games drive bad fills. |
| `MIN_EV_PCT` recalibration | After 1-2 weeks of fill data, tune based on actual EV realization. |
| `MAX_LIVE_RFQS` bump | If consistently top-N-bound after wider candidate pool ships. |
| `risk.line_move_ok` deletion | If `reference_lines` table is also dropped, clean up the dead function. |
| `f5_spread`/`f5_total` columns in `mlb_parlay_lines` → migrate to `mlb_target_lines` schema | If/when dashboard line selection is also moved to a multi-row model. |

---

## Version control

- **Branch:** `worktree-kalshi-mlb-rfq-line-source-pivot` (worktree). Active during implementation.
- **Files created/modified:** See "File map" section above.
- **Commit structure:** One commit per task in the implementation plan. Each commit passes its own tests before the next begins.
- **Merge prerequisite:** Full smoke test on live MVE markets (Test plan #5 + #6) before requesting merge approval.
- **Pre-merge review:** Executive engineer review of full diff per CLAUDE.md (data integrity, resource safety, edge cases, log/disk hygiene).
- **Approval gate:** Explicit user approval required before merging to main.

## Worktree lifecycle

- **Create:** Done — `EnterWorktree` invoked at design-writing time.
- **Implementation:** Conducted in this worktree, branch off main.
- **Cleanup:** After merge, `git worktree remove .claude/worktrees/kalshi-mlb-rfq-line-source-pivot` + `git branch -d worktree-kalshi-mlb-rfq-line-source-pivot`.

## Documentation updates

Required (same merge as code):

- `kalshi_mlb_rfq/README.md` — document the new SGP cadence loop, env vars, market DB sibling file, N≥2 book gate, line-source decoupling from Wagerzon.
- `mlb_sgp/README.md` (or create if missing) — document the library architecture (per-book modules + thin scraper shims), the dual-mode `MLB_SGP_DB_PATH` invocation pattern, the new `mlb_target_lines` table.
- `.env.example` in `kalshi_mlb_rfq/` — add new env vars.
- `CLAUDE.md` (root) — small note on the bot's market DB sibling under "Project Structure".
- Memory updates after merge: update `kalshi_mlb_rfq_line_source.md` to reflect resolved state; update `kalshi_mm.md` if relevant; update active-work pointers in `MEMORY.md`.
