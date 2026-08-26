# Issue #96 — leg surface: schema + ingest loop

**Status:** proposal, not approved. 2026-08-25.
**Branch:** `feature/mm-leg-surface` (worktree `/Users/callancapitolo/nflwork-surface`).
**Parent:** epic #94. Depends on #95 (done). Blocks #98, #99.

---

## 0. What #96 ships, and what it deliberately does not

Ships: a **standalone ingest loop** that refreshes single-leg book fairs for
the slate and writes them to a surface, plus the in-memory store #98 will
read.

Does **not** ship: any change to `main.py`, `router.py`, `on_demand.py`, or
the RFQ path. The maker is untouched. #98 adds the wiring (one construction +
one loop arm); #99 adds the age gate.

Consequence: `git diff main..HEAD` for #96 contains **zero lines** inside the
quote path. That is how "must not touch the RFQ path" gets verified, not by
inspection.

---

## 1. Identity — the one decision everything else hangs off

**The surface is keyed on the Kalshi game code (event-ticker suffix), not on
the Odds API `game_id`, and never on team names.**

Today's maker resolves a game like this
(`main._resolve_game_for_legs_uncached`):

```sql
SELECT game_id FROM mlb_target_lines
WHERE home_team=? AND away_team=? LIMIT 1
```

Team names only. `sgp_runner.enumerate_kalshi_targets` builds that table with
`schedule.get((home_team, away_team))` — also team names only. A doubleheader
collapses onto one row and `LIMIT 1` picks whichever. This is the exact
failure #95 hit with FanDuel's two `PHI @ SEA` rows: **a wrong number, not a
decline.**

The Kalshi suffix (`26AUG252138CLELAA`) already encodes date, ET first pitch
and both team codes, so it is unique per doubleheader game by construction.
It is also what `CanonicalLeg.game_id` carries, so #98's lookup becomes a
direct key hit with no translation layer.

The surface therefore owns its own slate discovery (§3) rather than reading
`mlb_target_lines`. Every row stores `game_start_time TIMESTAMPTZ` (UTC,
parsed from the suffix) alongside the key, and **every** book-side match — DK
and FD singles rows especially — is made on canonical teams **plus** start
time within tolerance.

> The team-only resolution in `mlb_target_lines` is a real pre-existing
> doubleheader bug in the live maker. Out of scope here; worth its own issue.

---

## 2. Schema

Own sibling DuckDB: **`kalshi_mlb_mm/kalshi_mlb_mm_surface.duckdb`**.

Separate write lock, per the issue's "do not contend with the maker's
existing write lock". Precedent: the maker already runs state / market /
research siblings. It does **not** go in `kalshi_mlb_mm_market.duckdb`, which
the pricing path reads.

### 2.1 `mlb_leg_surface` — current state, one row per leg per book

```sql
CREATE TABLE IF NOT EXISTS mlb_leg_surface (
    book             VARCHAR   NOT NULL,   -- 'fanduel' | 'betmgm' | ...
    game_id          VARCHAR   NOT NULL,   -- Kalshi suffix, e.g. 26AUG252138CLELAA
    game_start_time  TIMESTAMPTZ NOT NULL, -- UTC first pitch, parsed from the suffix
    period           VARCHAR   NOT NULL,   -- 'FG' | 'F5' | 'I1'
    market_type      VARCHAR   NOT NULL,   -- 'ml' | 'spread' | 'total'
    line             DOUBLE,               -- signed home-perspective; NULL for ml
    side             VARCHAR   NOT NULL,   -- 'home'/'away' | 'over'/'under'
    fair_prob        DOUBLE    NOT NULL,   -- devigged P(this side)
    -- raw prices retained for debugging (issue #96)
    raw_decimal      DOUBLE    NOT NULL,   -- this side, as posted
    raw_decimal_opp  DOUBLE    NOT NULL,   -- the other side, as posted
    raw_overround    DOUBLE    NOT NULL,   -- 1/raw + 1/raw_opp, pre-devig
    route            VARCHAR   NOT NULL,   -- 'structure' | 'singles'
    built_at         TIMESTAMPTZ NOT NULL  -- when the BOOK DATA was fetched
);
```

Notes that matter:

- **`built_at` is the fetch time of the underlying book payload**, not the
  time the row was devigged or written. A structure fetch that serves 72 legs
  stamps all 72 with one `built_at`. Otherwise #99 would read a stale price as
  fresh.
- **`line` is NULL for moneyline** — NULL, not a sentinel, matching the
  established `spread_line = NULL` convention for ml×total rows.
- **No PRIMARY KEY.** DuckDB PKs reject NULL, and `line` is legitimately NULL.
  Uniqueness is structural instead: the in-memory store is a dict on the key,
  and each flush is `DELETE WHERE book = ?` + insert for that book only. A
  duplicate key cannot be produced. Acceptance verifies it with a
  `GROUP BY … HAVING COUNT(*) > 1` query rather than trusting a constraint.
- **Both sides of a rung are stored** (2 rows per rung, `fair_prob` summing to
  1). #98 then never has to know which side the book "posted".
- `game_start_time` is denormalized onto every row so the monitor and the
  research queries never need a join to answer "how far from first pitch was
  this".

### 2.2 `surface_refresh_log` — one row per (book, pass)

This is how the acceptance criterion "`built_at` freshness stays within its
target cadence" gets measured, and how exclusions become **countable**.

```sql
CREATE TABLE IF NOT EXISTS surface_refresh_log (
    book              VARCHAR NOT NULL,
    route             VARCHAR NOT NULL,      -- 'structure' | 'singles'
    started_at        TIMESTAMPTZ NOT NULL,
    duration_sec      DOUBLE NOT NULL,
    games_attempted   INTEGER NOT NULL,
    games_priced      INTEGER NOT NULL,
    rungs_priced      INTEGER NOT NULL,      -- devig succeeded
    legs_written      INTEGER NOT NULL,      -- == 2 * rungs_priced
    -- fixed exclusion vocabulary, one column each; plain-SQL countable
    n_crossed         INTEGER NOT NULL,      -- implied sum < 1.0
    n_overround       INTEGER NOT NULL,      -- sum outside [1.005, 1.20]
    n_one_sided       INTEGER NOT NULL,      -- book posts only one side
    n_unresolved      INTEGER NOT NULL,      -- book does not offer that rung
    n_game_unmatched  INTEGER NOT NULL,      -- book does not list the game
    n_game_ambiguous  INTEGER NOT NULL,      -- 2+ book games matched (fail closed)
    error_class       VARCHAR                -- transport failure, else NULL
);
```

The vocabulary is **fixed and small**, so a reason is a column, not JSON.
Rows are pruned to `SURFACE_LOG_RETENTION_HOURS` (default 24) on a slow tick.

**Why exclusions are counted, not row-stored:** FanDuel alone carries ~51
one-sided deep-alt rungs per slate (#95). At a 20s cadence that is ~9k
exclusion rows/hour of pure noise. Counts answer every question the issue
asks; a DEBUG log line carries the first N examples per pass for when a
count looks wrong.

### 2.3 The quote path reads memory, not this table

`LegSurface` is an in-process dict keyed
`(book, game_id, period, market_type, line, side)` → `SurfaceRow`. The DuckDB
tables are the **durable mirror** for research, the monitor, and #96's own
acceptance — never the quote path's read.

This is not a nicety. The memory note is explicit: a DuckDB open costs ~17ms
and has caused three separate incidents when done per-item in a hot loop. The
epic's definition of done is "a cross-game RFQ prices with **zero** outbound
book requests"; paying 17ms of file I/O per RFQ instead would trade one
bottleneck for another. Flush cadence to disk: `SURFACE_DB_FLUSH_SEC` = 30.

---

## 3. Slate discovery (slow arm, Kalshi-only, zero book requests)

Every `SURFACE_SLATE_REFRESH_SEC` (default 300, matching
`TARGET_LINE_REFRESH_SEC`):

1. `GET /events?series_ticker=KXMLBGAME&status=open&limit=50`
2. Keep games whose suffix-derived first pitch is
   `SURFACE_GAME_MIN_MINUTES` (5, matching `TIPOFF_CANCEL_MIN`) to
   `SURFACE_GAME_MAX_HOURS` (12) away. **Measured 2026-08-25: 48 events are
   open at once (~3 days out).** Ingesting all of them triples cost for games
   where #95 showed books post main lines only. 12h covers the day's slate.
3. Per surviving game, enumerate the strike ladder from Kalshi:

   | series | markets/game (measured, 1 game, 2026-08-25) | legs |
   |---|---|---|
   | `KXMLBGAME` | 2 | 2 (ml home/away) |
   | `KXMLBSPREAD` | 8 | 16 |
   | `KXMLBTOTAL` | 13 | 26 |
   | `KXMLBF5` | 3 (incl. TIE) | 4 — TIE is unpriceable, dropped |
   | `KXMLBF5SPREAD` | 4 | 8 |
   | `KXMLBF5TOTAL` | 7 | 14 |
   | `KXMLBRFI` | 1 | 2 |
   | **total** | **38** | **~72 legs = ~36 rungs** |

   6 Kalshi calls per game per slate refresh — ~90 calls per 5 min on a
   15-game slate. Kalshi API only; **no book traffic.**
4. Convert each market to `CanonicalLeg`s with the **existing**
   `legset.parse_leg` re-encodings, so the surface can never disagree with the
   quote path about what a leg means: F5-winner → ±0.5 F5 spread (#86), RFI →
   I1 total at 0.5 (#87), spread signs home-perspective (#70).
5. `parse_suffix_start_utc` (today in `kalshi_rfi/discovery.py`) is promoted to
   `kalshi_common/leg_types.py`; `kalshi_rfi` re-exports it. A maker module
   importing `kalshi_rfi` would be wrong-shaped coupling.

The result is a `SurfaceSlate`: `{game_id: SurfaceGame(suffix, teams,
start_utc, legs)}`. Held in memory, mirrored to `mlb_surface_slate` for
debugging.

---

## 4. Ingest workers — per book, per route

Each book runs an **independent worker thread with its own cadence and its own
clock.** A slow book's pass never delays a fast book's; each writes only its
own keys and its own `built_at`.

### 4.1 Route assignment is explicit and exclusive

Exactly one route is authoritative per `(book, market_type, period)`, so a key
can never be written by two sources:

| book | route | why |
|---|---|---|
| betmgm | structure, all | #95: 21/21 leg-instances |
| novig | structure, all | #95: 20/21 |
| caesars | structure, all | #95: 8/21 — alive for singles; #90 is a *combo* stat, and its I1 mapping landed on main 2026-08-25 (7ead5ce) |
| fanduel | structure for ml/spread + **all** I1; **singles** for `(total, FG)` and `(total, F5)` | #95 fact 3: FD's SGP structure carries exactly ONE total line per period. Its singles scraper has the full ladder (312 FG-total rows/slate) |
| draftkings | **singles**, all | #95 fact 4: 21/21 `no_structure_odds` |
| prophetx | **disabled by default** | 403 at `events` on the first request of a session, 21/21. Budget zero legs |

Note the FD row: the split is keyed on `(market_type, period)`, not
`market_type` alone. Keying on `market_type` would send FD's **I1** legs
(which are `total` legs at 0.5) to the singles route, and neither singles
scraper emits I1 at all — FD would silently vanish from the thinnest row on
the surface.

### 4.2 Structure route — one flight per (book, game), ~36 devigs

Uses `SGPService(single_leg_structure_fair=True, structure_ttl_sec=<cadence>)`
— **a second, private instance**, so the maker's on-demand engine and its
420s structure cache are untouched.

Per game: `match_event` → `build_structure` (one wire fetch; `built_at` is
stamped here) → then, for each of ~36 rungs, a purely local
`resolve(structure, [leg], home, away)` → `(single_decimal,
opposite_decimal)` → devig. `resolve_legs` is all-or-nothing across the legs
handed to it, so legs are resolved **one at a time**: one missing rung must
not unprice the other 35.

This is exactly the #95 latency finding — one cold fetch per game (median
0.31–0.42s), then ~0.00s per additional leg, zero SGP price calls.

`GameRef.commence_time` is the **suffix-derived** first pitch. The books'
`match_events` helpers already bucket on UTC date+hour, so a correct
commence_time makes doubleheaders resolve correctly rather than needing to be
dropped.

### 4.3 Singles route — one whole-slate scrape per pass

`mlb_sgp/scraper_draftkings_singles.py` and `scraper_fanduel_singles.py` are
refactored:

```python
def collect_singles_rows(verbose=False) -> list[dict]:   # NEW: pure scrape
def scrape_singles(verbose=False) -> int:                # collect + write, unchanged
```

The surface calls `collect_singles_rows`. **It never calls `scrape_singles`.**

Deliberate, per the brief's question. Both alternatives are worse:
- *Read `dk_odds/dk.duckdb` / `fd_odds/fd.duckdb`* — those have no cron and
  were last written 2026-08-17. A stale snapshot as a maker's pricing input is
  the #54 failure again, with worse numbers.
- *Call `scrape_singles`* — each call does a full `CREATE OR REPLACE` of the
  production `mlb_odds` table that MLB.R and the dashboard read. The maker's
  ingest cadence would silently become the dashboard's refresh rate. An
  outward side effect the maker has no business owning.

Matching book rows to Kalshi games: canonical team names (both scrapers emit
Odds-API canonical forms) **plus** `|book_start − kalshi_start| ≤
SURFACE_START_TOLERANCE_MIN` (30). Zero matches → `n_game_unmatched`; **two or
more matches → `n_game_ambiguous`, fail closed, price nothing for that game at
that book.** Prices are American integers → `mlb_sgp._shared.american_to_decimal`.

Period vocabulary maps straight through (`FG`/`F5`); `F3`/`F7` rows are
ignored; there are no I1 rows to ignore.

### 4.4 Devig — two-way on the line the book posted

New shared helper in `kalshi_common/fair_value.py` (additive; nothing existing
changes behaviour):

```python
def two_way_fair(dec_chosen, dec_opposite, *, band_min, band_max):
    """(fair, None) or (None, reason) with reason in
    {'crossed', 'overround'}. Gate runs BEFORE devig."""
```

Order of operations, per the issue and the #73 precedent
(`unabated_edge/pricing.overround_reject`):

1. `overround = 1/dec_chosen + 1/dec_opposite`
2. `overround < 1.0` → **`crossed`**, excluded. One side refreshed while the
   other was stale; arithmetically not a two-way market.
3. `overround` outside `[1.005, 1.20]` → **`overround`**, excluded.
4. Otherwise `devig_two_way` (exact probit) → store both sides.

One side missing → **`one_sided`**, excluded, never haircut. The Route-B
vig-fallback haircut (`ON_DEMAND_VIG_FALLBACK`) is deliberately **not** reused
here: it is calibrated for cancelling a combo's compounded margin, and applying
it to a lone leg produces a confident wrong number where a decline costs
nothing — #95 showed 4–5 books per leg, so one book dropping out never
approaches quorum.

`devig_partition` is *not* used. It exists for 2^N cells with a per-leg
overround budget; a two-way rung wants the tighter, explicit `[1.005, 1.20]`
envelope the issue specifies.

**Absence ≠ failure.** #95 fact 3: at T−15h most books post main lines only
(DK went ~14 spread lines/game → exactly 1; MGM's F5 and YRFI markets were
absent entirely for next-day games). `n_unresolved` is therefore an expected,
time-of-day-dependent count, never a book-health signal, and it does not feed
`BookHealthAlerter`. Only transport failures do.

---

## 5. Cadence — and the traffic number that has to be settled first

**This is the one place I think the epic's stated numbers do not survive
contact.** The issue quotes "~6 requests per cycle, roughly 0.4–0.7 req/sec".
That reads a whole-book pull as one request. It is not — it is one fetch *per
game*.

Structure books, 15-game slate, ~1–2 HTTP calls per game-structure:

| cadence | book-HTTP req/sec (4 structure books) | req/day |
|---|---|---|
| 5s | ~12–24 | 1.0–2.1 M |
| 15s | ~4–8 | 350–690 k |
| **20s** | **~3–6** | **260–520 k** |
| 60s | ~1–2 | 86–170 k |

For scale: the congested on-demand path served **219,724 requests on
2026-08-19** — and ProphetX 403s on sight. A 5s cadence would be 5–10× that
volume, permanently. The epic's own framing ("compare to the 15,000+/hour the
on-demand path burned") is only true at a slow cadence.

**Proposed defaults**

| book | route | cadence | pass cost (#95 median) |
|---|---|---|---|
| betmgm | structure | 20s | 2.58s |
| novig | structure | 20s | 2.15s |
| caesars | structure | 20s | 2.82s |
| fanduel | structure | 20s | 3.78s |
| fanduel | singles (totals) | 45s | 13.2s |
| draftkings | singles | 60s | 28.4s |

Plus a per-book request-rate ceiling (`SURFACE_MAX_REQ_PER_SEC_PER_BOOK`,
default 2.0) that stretches a pass rather than firing it, so a misconfigured
cadence cannot become a self-inflicted 403.

**Two consequences #99 has to own, stated here so they are not a surprise:**

1. At 20s, structure rows are typically 0–20s old — inside a 30s
   `SURFACE_MAX_AGE_SEC`, but with no margin. Any pass that overruns pushes
   rows past the gate.
2. **DraftKings cannot meet a 30s gate at all.** Its slate scrape is 28.4s
   before cadence, so DK rows are routinely 30–90s old. Under a uniform 30s
   gate DK contributes ~nothing and the surface is effectively FD/MGM/NV/CZR.
   That still clears `MIN_AGREEING_BOOKS=2` comfortably (#95: 3–4 books on
   every FG/F5 leg structure-only) — so this is survivable, but it should be
   a decision, not a discovery.

My read: the age gate wants to be **90–120s**, not 30s. The design doc's own
§6 argues it — quotes already rest 30s (median, expired) to 1,654s
(cancelled), freshness at *fill* is what matters, and that is guarded by the
10s constituent-jump breaker plus #99's pre-quote Kalshi freshness veto, which
is book-independent and free. A 30s book-age gate mostly buys the right to
make 4× the book requests. The user has already decided 30s; I have built the
knob per-book so #99 can settle it with data instead of re-litigating it.

---

## 6. Config (all `kalshi_mlb_mm/config.py`, env-overridable)

```
SURFACE_ENABLED=false            # #96 ships dark; #98 turns it on
SURFACE_DB                       = kalshi_mlb_mm_surface.duckdb
SURFACE_BOOKS_STRUCTURE=fanduel,betmgm,novig,caesars
SURFACE_BOOKS_SINGLES=draftkings,fanduel
SURFACE_SLATE_REFRESH_SEC=300
SURFACE_GAME_MAX_HOURS=12
SURFACE_GAME_MIN_MINUTES=5       # == TIPOFF_CANCEL_MIN
SURFACE_CADENCE_DEFAULT_SEC=20
SURFACE_CADENCE_DRAFTKINGS_SEC=60
SURFACE_CADENCE_FANDUEL_SINGLES_SEC=45
SURFACE_MAX_REQ_PER_SEC_PER_BOOK=2.0
SURFACE_OVERROUND_MIN=1.005
SURFACE_OVERROUND_MAX=1.20
SURFACE_START_TOLERANCE_MIN=30
SURFACE_DB_FLUSH_SEC=30
SURFACE_LOG_RETENTION_HOURS=24
```

## 7. Files

```
kalshi_mlb_mm/leg_surface/__init__.py    LegSurface store + SurfaceRow
kalshi_mlb_mm/leg_surface/slate.py       Kalshi slate + leg enumeration
kalshi_mlb_mm/leg_surface/structure.py   structure-route worker
kalshi_mlb_mm/leg_surface/singles.py     singles-route worker + game matching
kalshi_mlb_mm/leg_surface/store.py       DuckDB schema, flush, prune
kalshi_mlb_mm/leg_surface/runner.py      workers + cadence; SurfaceIngest
kalshi_mlb_mm/leg_surface/__main__.py    standalone runner (acceptance)
kalshi_mlb_mm/tests/test_leg_surface_*.py

MODIFIED (additive only):
kalshi_common/fair_value.py              + two_way_fair()
kalshi_common/leg_types.py               + parse_suffix_start_utc()
kalshi_rfi/discovery.py                  re-export the moved parser
mlb_sgp/scraper_draftkings_singles.py    split collect_singles_rows / scrape_singles
mlb_sgp/scraper_fanduel_singles.py       same split
kalshi_mlb_mm/config.py                  + SURFACE_* knobs
```

`SurfaceIngest` is an object with `.tick()` (like `OnDemandEngine`), so #98
wires it into `main_loop` with one construction and one loop arm. #96 does not
add that arm.

## 8. Testing

- Unit: devig gate (crossed / overround / one-sided / pass) against
  hand-built decimals; leg enumeration for all 7 Kalshi series incl. the F5
  TIE drop and the RFI→I1 re-encoding; the singles game-matcher against the
  real FD two-`PHI @ SEA` case from #95 (must pick the right one, and must
  fail closed when two land inside tolerance).
- Integration (live, ~1h, per the acceptance criterion): run
  `python -m kalshi_mlb_mm.leg_surface` standalone; then verify with SQL —
  per-book `built_at` p50/p95 age vs target cadence, `GROUP BY key HAVING
  COUNT(*) > 1` returns 0 rows, rung counts per game sane vs the ~36 measured,
  exclusion counts non-zero and attributable.
- Isolation proof: `git diff main..HEAD -- kalshi_mlb_mm/main.py
  kalshi_mlb_mm/router.py kalshi_mlb_mm/on_demand.py` is empty.
- The maker bot is **not running** (died 2026-08-19) and this spec does not
  start it.

## 9. Version control

- Worktree `/Users/callancapitolo/nflwork-surface`, branch
  `feature/mm-leg-surface`, merged up to `main` @ 7ead5ce.
- Commits: (1) shared helpers — `two_way_fair`, `parse_suffix_start_utc`,
  singles-scraper split, all with tests; (2) slate + schema + store;
  (3) structure and singles workers + runner; (4) docs.
- Never symlink a `.duckdb` into the worktree — copy.
- Pre-merge review on the full diff against the CLAUDE.md checklist, then
  explicit user approval before merging. Clean up worktree + branch after.

## 10. Documentation (same merge)

- `kalshi_mlb_mm/README.md` — surface tables, the two routes, per-book
  cadence, the DK freshness caveat.
- `CLAUDE.md` — the maker blurb currently says live-only on-demand is the ONLY
  pricing path.
- `kalshi_mlb_mm/research_queries.sql` — freshness-by-book and
  exclusion-by-reason queries.

## 11. Flagged, not fixed here

1. **DraftKings `calculateBets` is 403 for every set size**, n=1 and n=2, on
   both documented paths, verified live 2026-08-25. Read endpoints are green.
   This is #39 recurring with its `/en/` workaround closed, and it is a **live
   outage of the maker's same-game path today**, independent of this epic.
   Deserves its own issue. The repo's docs — `CLAUDE.md`,
   `kalshi_common/sgp_service.py` comments, and the #95 comment — attribute it
   to DK "refusing 1-selection sets", which is the wrong cause;
   `kalshi_rfi/README.md` on main is the only place that states it correctly.
2. **Team-only game resolution** in `mlb_target_lines` /
   `_resolve_game_for_legs` — a live doubleheader wrong-number risk in the
   maker. The surface routes around it; the maker still has it.
3. **ProphetX 403** (#91) — upside only, not a blocker.
