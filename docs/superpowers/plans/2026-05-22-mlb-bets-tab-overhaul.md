# MLB Bets Tab Overhaul + Timezone Standardization

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

## Review Pack

**What we're building** — All 8 MLB-relevant scrapers stop writing `(game_date VARCHAR, game_time VARCHAR)` and start writing one `game_start_time TIMESTAMPTZ` column in UTC. The R side deletes 4 per-book parser functions. Along the way we fix the BKM-PT-not-ET bug that silently drops ~half of Bookmaker's MLB rows today, and add coverage / freshness / dedup defenses to the bets-tab pipeline on top.

**Key decisions**

1. **Merge the TZ refactor and the bets-tab fixes into one branch** rather than shipping them sequentially. Rejected: two PRs. Why: the BKM-PT bug is fully fixed only when both ship, and a single cutover with one revert path is cheaper than two coordinated ones.
2. **Drop the `(game_date, game_time)` pair entirely** rather than keeping it in parallel with `game_start_time`. Rejected: dual columns for backward compatibility. Why: timezone is metadata that belongs *inside* the timestamp; keeping both forms preserves the bug surface forever.
3. **Empirically audit each "claimed naive" scraper before writing migration code** (Phase 1) rather than trusting the existing labels. Rejected: trust the code comments. Why: BKM was labeled ET but actually wrote PT — that's hard evidence the labels are unreliable for any of WZ/Hoop88/Bet105.
4. **Treat the per-scraper DuckDB tables as ephemeral** (drop + recreate) rather than ALTER-in-place. Rejected: in-place migrations everywhere. Why: scraper tables fully repopulate every 1-5 minutes, so drop+recreate loses at most one cycle of stale data. Only `wagerzon_specials` is persistent and gets ALTER treatment.
5. **F3 markets added to BOTH DK and FD** (reverses the prior "F3 not in scope" decision). Rejected: leave F3 hidden. Why: your direction earlier in this session — you want F3 pills across all books, not just Wagerzon/Bookmaker.

**Risks / push back here**

- **Scope is ~3x the original plan.** Realistic 4-8 hours attended. If you'd rather start with a 1-line BKM hotfix in `Tools.R::.parse_wz_game_dt` (have it interpret BKM as PT) and do the rest later, the dashboard recovers tonight at near-zero risk.
- **Bots stopped during cutover.** Kalshi MM + RFQ should be down during Phases 3-4. Probably an off-hours window — flag if that's a problem.
- **`wagerzon_specials` is the only persistent table.** ALTER+UPDATE goes wrong → unrecoverable without a snapshot. Worth `cp wagerzon.duckdb wagerzon.bak.duckdb` before running the migration.
- **Audit may surface MORE surprises.** If WZ or Hoop88 also writes a wrong TZ, plan still works but each new offset is one more place to verify post-merge.
- **Out-of-scope by choice:** NFL Draft portal + Kalshi RFQ market DB have their own timestamp surfaces. Leaving them untouched is a judgment call — if you'd rather sweep them in too, scope grows again.

**Worth understanding** (3 portable concepts, opt-in)

1. **TIMESTAMPTZ vs naive TIMESTAMP** — In R the analogue is `POSIXct` with a non-empty `tz=` attribute vs an empty one. A TZ-aware timestamp carries its own UTC offset, so comparing it to `Sys.time()` or DuckDB `NOW()` is unambiguous. A naive timestamp is just digits — you have to *know* the assumed TZ to interpret it, and if your assumption is wrong (BKM!) you get silent corruption with no error. This refactor's whole premise is "stop relying on assumptions; put the TZ in the value."
2. **Atomic write via staging** — Today's DK/FD pattern is `BEGIN; DELETE; INSERT; COMMIT`. If the Python process is killed between DELETE and INSERT completing, the table can end up empty. The fix is `CREATE OR REPLACE TABLE x AS SELECT * FROM staging` in one statement. R analogue: rather than mutating a data frame in place (where a crash leaves it half-mutated), build the replacement as a separate object and assign it at the end. Reader sees either old or new, never halfway.
3. **`anti_join` for coverage diagnostics** — You know `inner_join` from R. `anti_join(raw, lookup)` is its complement: rows in `raw` that DON'T appear in `lookup`. It's the cleanest way to surface "what's about to be silently dropped" before a destructive join — exactly the pattern in Task 13. This pairs with the principle of "make silent failures loud."

---

**Goal:** Two merged workstreams in one branch:
1. **TZ standardization** — Replace per-scraper `(game_date VARCHAR, game_time VARCHAR)` with one canonical `game_start_time TIMESTAMPTZ` (UTC) column across all 8 MLB-relevant scrapers. Delete 4 per-book parser functions in R. Fix the **BKM-PT-not-ET bug** that silently drops ~half of BKM's MLB rows.
2. **Bets-tab overhaul** — Defensive logging, freshness gates, and bug fixes across the matching layer and dashboard render, on top of the new TZ-clean foundation.

**Architecture:** One canonical timestamp, written once, parsed never. The R side loses its 4 `.parse_*_game_dt` functions. The 18-column scraper schema becomes 17 columns. CLAUDE.md loses Pitfalls #7 and #11. Bets-tab defensive checks operate on the new TZ-clean substrate.

**Tech Stack:** Python 3.11+ (duckdb, zoneinfo, regex), R 4.x (dplyr, duckdb, tibble, stringr, lubridate), DuckDB, MLB Odds API.

**Scope explicitly excluded:**
- NFL Draft portal scrapers (separate pipeline, separate DuckDB; document but don't change).
- Kalshi RFQ bot internal market DB (read-only on mlb.duckdb; not touched).
- Historical PBP DuckDBs (cbb_2021..2026, pbp.duckdb — own conventions, already correct).
- `LINE_MATCH_TOLERANCE` tightening (deliberate 2026-05-13 widening).
- Main/alt market union "fix" (audit hallucination — actual filter is strict equality).

---

## Worktree lifecycle

- **Branch:** `worktree-mlb-bets-tab-overhaul` (already created via `EnterWorktree`).
- **Path:** `/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-overhaul`
- **Operating mode during migration:** Per Phase 6 verification, the Kalshi MM bot and Kalshi MLB RFQ bot should be **stopped** during the cutover window (verify with `ps aux | grep -E "mm_bot|kalshi_mlb_rfq"`). Dashboards can keep running; they'll re-read on next refresh.
- **At completion:** End-to-end smoke + dashboard verification + parity test green + diff review surfaced to user. On user approval: merge to `main`, `git worktree remove`, `git branch -d worktree-mlb-bets-tab-overhaul`. Bots restart against new schema.

---

## Phase map (high-level)

| Phase | What | Tasks | Goal |
|---|---|---|---|
| P1 | TZ audit | 1 | Empirically confirm each Group B scraper's actual TZ. Source of truth for Phase 3. |
| P2 | Migration tooling | 2 | One-shot DB schema migrator (drop+recreate ephemeral, ALTER for persistent). |
| P3 | Per-scraper migrations | 9 | 8 scrapers + wagerzon_specials. Each commit: schema + atomic write + book-specific fixes. |
| P4 | R-side collapse | 1 | Delete parsers, simplify `.drop_past_games`, prune dashboard server defensives. |
| P5 | Bets-tab defenses | 4 | Coverage logging, fetch_time gate, pivot dedup, bet_row_id collision detection. |
| P6 | Parity harness | 1 | `tests/timezone_parity_test.py` as the regression gate. |
| P7 | Verification | 1 | End-to-end MLB + CBB pipeline smoke + dashboard render. |
| P8 | Docs | 1 | CLAUDE.md / READMEs / memory updates. |

Total commits: ~17–20. Estimated worktree lifetime: 4–8 hours of attended execution.

---

## Pre-mortem (known dangers)

Before kicking off, the realistic risk inventory:

1. **Big-bang merge means the parity harness is load-bearing.** If `timezone_parity_test.py` misses a scenario, a bug ships. Mitigation: harness covers every (scraper × MLB-game) pair for the next 24h and is run *twice* — once pre-refactor (where some scrapers MUST fail, proving the harness detects offset) and once post-refactor (where all must pass).
2. **Scope creep risk into NFL Draft / Kalshi RFQ.** Both contain timestamp surfaces. Out of scope here; documented in CLAUDE.md as "follow this convention if rewritten."
3. **Migration timing during a live trading day.** Running mid-day could leave the pipeline broken for tens of minutes. Recommendation: execute Phases 3–4 with bots stopped; dashboard read-only during cutover.
4. **wagerzon_specials is the only persistent table.** Loss of rows is unrecoverable. Migration uses `ALTER ADD` first, `UPDATE` to populate, `ALTER DROP` last — atomic at each step.
5. **R sessions holding cached parser closures.** Any long-running R process must be restarted post-merge or it'll call deleted functions. Mitigation: pipeline orchestrator spawns fresh R processes per run; no long-lived R.
6. **BKM-PT bug discovered today implies the "claimed-ET" tags on WZ, Hoop88, Bet105 are unverified.** Mitigation: Phase 1 audit is mandatory before Phase 3 migrations.

---

## File map

| File | Phases that touch it |
|---|---|
| `tests/timezone_audit.py` | P1 (new) |
| `tests/timezone_parity_test.py` | P6 (new) |
| `tools/migrate_scraper_schemas.py` | P2 (new) |
| `mlb_sgp/scraper_draftkings_singles.py` | P3a |
| `mlb_sgp/scraper_fanduel_singles.py` | P3b |
| `mlb_sgp/test_fd_alt_parsing.py` | P3b (new) |
| `bfa_odds/scraper*.py` | P3c |
| `kalshi_odds/scraper*.py` | P3d |
| `wagerzon_odds/scraper.py` + `scraper_specials.py` | P3e, P3i |
| `hoop88_odds/scraper*.py` | P3f |
| `bookmaker_odds/scraper.py` | P3g (the headline fix) |
| `bet105_odds/scraper*.py` | P3h |
| `Answer Keys/Tools.R` | P4 |
| `Answer Keys/MLB Answer Key/MLB.R` | P5 (callers + freshness gate) |
| `Answer Keys/MLB Answer Key/odds_screen.R` | P5 (coverage logging) |
| `Answer Keys/MLB Answer Key/mlb_triple_play.R` | P4 (drop ICU workaround) |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | P5 (pivot dedup + collision) |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | P4 (drop UTC-defensive branches) |
| `Answer Keys/CBB Dashboard/cbb_dashboard_server.py` | P4 |
| `Answer Keys/CLAUDE.md`, scraper READMEs, memory | P8 |

---

# PHASE 1 — TZ Audit Foundation

### Task 1: `tests/timezone_audit.py` — empirically verify Group B scrapers

**Why:** BKM was claimed ET in code comments and CLAUDE.md but actually writes PT. We must not trust the other claimed-naive scrapers (WZ, Hoop88, Bet105) without verification.

**Files:**
- Create: `tests/timezone_audit.py`

- [ ] **1.1:** Read the existing `Answer Keys/Tools.R::.parse_wz_game_dt`, `.parse_bet105_game_dt`, `.parse_bfa_game_dt`, `.parse_iso_game_dt` (lines ~3618-3666) to document what each scraper CLAIMS:
  ```bash
  sed -n '3610,3700p' "Answer Keys/Tools.R"
  ```
- [ ] **1.2:** Create `tests/timezone_audit.py`:
  ```python
  """Empirically verify each scraper's recorded game_time TZ against Odds API.

  Picks the next 24h of MLB games. For each game appearing in scraper_db.mlb_odds
  and Odds API, computes Δ = scraper_recorded_dt − odds_api_commence_time
  (assuming naive scraper writes are UTC). The modal Δ across games reveals
  the actual TZ offset.
  """
  import duckdb, json, os
  from datetime import datetime, timezone, timedelta
  from pathlib import Path
  from urllib.request import urlopen, Request

  # Path → scraper label
  SCRAPERS = {
      "wagerzon_odds/wagerzon.duckdb":   "wagerzon",
      "hoop88_odds/hoop88.duckdb":       "hoop88",
      "bookmaker_odds/bookmaker.duckdb": "bookmaker",
      "bet105_odds/bet105.duckdb":       "bet105",
  }
  ODDS_API_KEY = os.environ.get("ODDS_API_KEY") or open(os.path.expanduser("~/.Renviron")).read().split("ODDS_API_KEY=")[1].split()[0].strip('"')

  def fetch_odds_api_games():
      url = f"https://api.the-odds-api.com/v4/sports/baseball_mlb/odds?apiKey={ODDS_API_KEY}&regions=us&markets=h2h&oddsFormat=american"
      return json.loads(urlopen(Request(url)).read())

  def parse_naive_as_utc(date_str, time_str):
      """Treat scraper-written (game_date, game_time) as naive UTC for diff baseline."""
      return datetime.strptime(f"{date_str} {time_str}", "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)

  def audit():
      api_games = {(g["home_team"], g["away_team"]): datetime.fromisoformat(g["commence_time"].replace("Z","+00:00"))
                   for g in fetch_odds_api_games()}
      for db_path, label in SCRAPERS.items():
          full_path = Path(db_path)
          if not full_path.exists():
              print(f"[{label}] missing DB at {db_path} — skip"); continue
          con = duckdb.connect(str(full_path), read_only=True)
          rows = con.execute("SELECT DISTINCT home_team, away_team, game_date, game_time FROM mlb_odds WHERE game_time IS NOT NULL").fetchall()
          deltas = []
          for ht, at, gd, gt in rows:
              api_dt = api_games.get((ht, at))
              if api_dt is None: continue
              try: scraper_dt = parse_naive_as_utc(gd, gt)
              except ValueError: continue
              delta_h = (scraper_dt - api_dt).total_seconds() / 3600
              deltas.append((ht, at, round(delta_h, 2)))
          print(f"\n=== {label} ({len(deltas)} matched games) ===")
          if not deltas: print("  no matches"); continue
          for ht, at, d in deltas: print(f"  {ht} vs {at}: Δ = {d:+.2f}h")
          modal = max(set(d for _,_,d in deltas), key=[d for _,_,d in deltas].count)
          print(f"  → MODAL Δ: {modal:+.2f}h  ⇒  {tz_guess(modal)}")
          con.close()

  def tz_guess(delta_h):
      """Given Δ (naive-treated-as-UTC minus actual UTC), guess source TZ."""
      mapping = {0.0: "UTC", 5.0: "America/New_York (EST)", 4.0: "America/New_York (EDT)",
                 8.0: "America/Los_Angeles (PST)", 7.0: "America/Los_Angeles (PDT)"}
      return mapping.get(round(delta_h, 1), f"unknown offset {delta_h:+.1f}h")

  if __name__ == "__main__": audit()
  ```
- [ ] **1.3:** Run it:
  ```bash
  python tests/timezone_audit.py 2>&1 | tee /tmp/tz_audit.log
  ```
- [ ] **1.4:** Record the modal Δ per scraper in a "P1 audit findings" comment at the top of `tools/migrate_scraper_schemas.py` (created in Phase 2). Expected based on prior chat: BKM ≈ +7h (PDT), WZ/Hoop88 ≈ +4h (EDT), Bet105 ≈ 0h (UTC). **If any value surprises, stop and investigate.**
- [ ] **1.5:** Commit.
  ```bash
  git add tests/timezone_audit.py
  git commit -m "test: timezone audit script for naive-claimed scrapers

  Empirically computes Δ between each scraper's recorded game_time and
  Odds API commence_time to reveal the actual source TZ. Output is the
  source of truth for Phase 3 per-scraper conversions. Stays in repo
  as a regression check."
  ```

---

# PHASE 2 — Migration tooling

### Task 2: `tools/migrate_scraper_schemas.py`

**Why:** All 8 scraper DuckDBs need their `mlb_odds` table replaced with the new 17-column schema. Per-scraper tables are ephemeral (fully repopulated every cycle), so drop+recreate is safe. wagerzon_specials is persistent and uses ALTER instead. One coordinated script avoids per-DB hand-mangling.

**Files:**
- Create: `tools/migrate_scraper_schemas.py`

- [ ] **2.1:** Skeleton:
  ```python
  """Coordinated schema migration for all scraper DuckDBs.

  Ephemeral scraper tables (dk_odds, fd_odds, bfa_odds, kalshi_odds,
  wagerzon_odds, hoop88_odds, bookmaker_odds, bet105_odds): DROP TABLE
  mlb_odds. The scraper recreates with the new schema on next cycle.

  Persistent table (wagerzon_specials): ALTER ADD game_start_time
  TIMESTAMPTZ, UPDATE from existing naive ET pair, DROP game_date+game_time.

  Idempotent. --rollback drops new column / restores old schema.
  """
  import argparse, duckdb
  from pathlib import Path

  EPHEMERAL_DBS = {
      "dk_odds/dk.duckdb":             "dk",
      "fd_odds/fd.duckdb":             "fd",
      "bfa_odds/bfa.duckdb":           "bfa",
      "kalshi_odds/kalshi.duckdb":     "kalshi",
      "wagerzon_odds/wagerzon.duckdb": "wagerzon",
      "hoop88_odds/hoop88.duckdb":     "hoop88",
      "bookmaker_odds/bookmaker.duckdb":"bookmaker",
      "bet105_odds/bet105.duckdb":     "bet105",
  }
  WZ_SPECIALS_DB = "wagerzon_odds/wagerzon.duckdb"

  def migrate_ephemeral(rollback=False):
      for db_path, label in EPHEMERAL_DBS.items():
          if not Path(db_path).exists():
              print(f"[{label}] DB missing — skip"); continue
          con = duckdb.connect(db_path)
          try:
              con.execute("DROP TABLE IF EXISTS mlb_odds")
              print(f"[{label}] dropped mlb_odds (will be recreated by next scraper run)")
          finally:
              con.close()

  def migrate_specials(rollback=False):
      con = duckdb.connect(WZ_SPECIALS_DB)
      try:
          cols = [r[0] for r in con.execute("DESCRIBE wagerzon_specials").fetchall()]
          if rollback:
              if "game_start_time" in cols:
                  con.execute("ALTER TABLE wagerzon_specials DROP COLUMN game_start_time")
                  print("[wz_specials] rolled back: dropped game_start_time")
              return
          if "game_start_time" not in cols:
              con.execute("ALTER TABLE wagerzon_specials ADD COLUMN game_start_time TIMESTAMPTZ")
              con.execute("""
                UPDATE wagerzon_specials
                SET game_start_time = (CAST(game_date AS DATE) + CAST(game_time AS TIME))
                                       AT TIME ZONE 'America/New_York'
                WHERE game_start_time IS NULL
                  AND game_date IS NOT NULL AND game_time IS NOT NULL
              """)
              print(f"[wz_specials] added game_start_time; rows updated: "
                    f"{con.execute('SELECT count(*) FROM wagerzon_specials WHERE game_start_time IS NOT NULL').fetchone()[0]}")
          # NOTE: Drop game_date/game_time in a follow-up commit after we confirm
          # post-merge stability. Keep them in place for one cycle as a safety net.
      finally:
          con.close()

  def main():
      p = argparse.ArgumentParser()
      p.add_argument("--rollback", action="store_true")
      args = p.parse_args()
      migrate_ephemeral(rollback=args.rollback)
      migrate_specials(rollback=args.rollback)

  if __name__ == "__main__": main()
  ```
- [ ] **2.2:** Dry-run mental check — script is idempotent (DROP IF EXISTS, ADD COLUMN only if absent). Safe to re-run.
- [ ] **2.3:** Do NOT execute yet. Phase 3 scraper code must land first or the next scrape cycle will fail.
- [ ] **2.4:** Commit.
  ```bash
  git add tools/migrate_scraper_schemas.py
  git commit -m "tool: migrate_scraper_schemas.py for TZ standardization cutover

  Idempotent: DROPs ephemeral mlb_odds tables (scrapers recreate them)
  and adds game_start_time TIMESTAMPTZ to wagerzon_specials with backfill
  from existing naive-ET pair. --rollback for safe revert.

  Run AFTER scraper code in Phase 3 lands, BEFORE next scrape cycle."
  ```

---

# PHASE 3 — Per-scraper migrations

Each task in Phase 3 follows the same template: rewrite schema + write path, atomic write, scraper-specific fixes, smoke. One commit per scraper.

### Task 3 (P3a): DK scraper — schema + atomic + F3 + exact-prefix team match + sign validation

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_singles.py`

- [ ] **3.1:** Update schema in `write_to_duckdb`. Replace the 18-column CREATE TABLE (lines 342-365) with 17-column:
  ```sql
  CREATE TABLE IF NOT EXISTS mlb_odds (
      fetch_time         TIMESTAMPTZ,
      sport_key          VARCHAR,
      game_id            VARCHAR,
      game_start_time    TIMESTAMPTZ,
      away_team          VARCHAR,
      home_team          VARCHAR,
      market             VARCHAR,
      period             VARCHAR,
      away_spread        FLOAT,
      away_spread_price  INTEGER,
      home_spread        FLOAT,
      home_spread_price  INTEGER,
      total              FLOAT,
      over_price         INTEGER,
      under_price        INTEGER,
      away_ml            INTEGER,
      home_ml            INTEGER
  )
  ```
- [ ] **3.2:** Update row-emission in `parse_selections_to_wide_rows`. Replace `"game_date": fetch_time.strftime(...)` and `"game_time": event.start_time` with:
  ```python
  "game_start_time": event.start_time,  # already a UTC datetime from dk_client
  ```
- [ ] **3.3:** Replace `DELETE FROM mlb_odds` + INSERT pattern with atomic:
  ```python
  if rows:
      con.execute("CREATE OR REPLACE TEMP TABLE mlb_odds_new AS SELECT * FROM mlb_odds LIMIT 0")
      tuples = [tuple(r.get(c) for c in cols) for r in rows]
      placeholders = ", ".join(["?"] * len(cols))
      con.executemany(f"INSERT INTO mlb_odds_new VALUES ({placeholders})", tuples)
      con.execute("CREATE OR REPLACE TABLE mlb_odds AS SELECT * FROM mlb_odds_new")
  # If rows is empty, leave prior snapshot untouched (don't nuke).
  ```
  Update `cols` list to drop `game_date`/`game_time`, add `game_start_time`.
- [ ] **3.4:** Remove F3 skip block (lines 95-97); add F3 to period detection:
  ```python
  if "1st 3 innings" in n or "first 3 innings" in n:    period = "F3"
  elif "1st 5 innings" in n or "first 5 innings" in n:  period = "F5"
  elif "1st 7 innings" in n or "first 7 innings" in n:  period = "F7"
  else:                                                  period = "FG"
  ```
- [ ] **3.5:** Exact-prefix team match replaces substring (~lines 227-235):
  ```python
  sel_name = sel.name.strip()
  if sel_name.startswith(event.home_team + " ") or sel_name == event.home_team:
      side = "home"
  elif sel_name.startswith(event.away_team + " ") or sel_name == event.away_team:
      side = "away"
  else:
      continue
  ```
- [ ] **3.6:** Spread sign validation before emit:
  ```python
  if row.get("home_spread") is not None and row.get("away_spread") is not None:
      if abs(row["home_spread"] + row["away_spread"]) > 1e-9:
          print(f"[dk_singles] WARN: asymmetric spread h={row['home_spread']} a={row['away_spread']}", flush=True)
          continue
  ```
- [ ] **3.7:** Live-data smoke:
  ```bash
  python tools/migrate_scraper_schemas.py   # drops dk_odds/mlb_odds
  python mlb_sgp/scraper_draftkings_singles.py mlb --verbose 2>&1 | tail -10
  python -c "
  import duckdb
  c = duckdb.connect('dk_odds/dk.duckdb', read_only=True)
  print('schema:', c.execute('DESCRIBE mlb_odds').fetchall())
  print('F3 rows:', c.execute(\"SELECT count(*) FROM mlb_odds WHERE period='F3'\").fetchone()[0])
  print('game_start_time type:', c.execute('SELECT typeof(game_start_time) FROM mlb_odds LIMIT 1').fetchone())
  "
  ```
  Confirm: schema is 17 cols, F3 rows > 0, game_start_time is TIMESTAMPTZ.
- [ ] **3.8:** Commit.
  ```bash
  git add mlb_sgp/scraper_draftkings_singles.py
  git commit -m "feat(dk-singles): TZ standardization + F3 + atomic write + hardening

  Schema migration: (game_date, game_time) → game_start_time TIMESTAMPTZ.
  Atomic write via CREATE OR REPLACE FROM staging (empty scrapes no
  longer nuke the table). F3 markets now emitted (was 'not in scope').
  Spread side detection: exact-prefix match (not substring). Spread sign
  symmetry enforced. dk_odds/mlb_odds requires drop+recreate; coordinate
  via tools/migrate_scraper_schemas.py."
  ```

---

### Task 4 (P3b): FD scraper — schema + atomic + F3 + alt regex + sign validation

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_singles.py`
- Create: `mlb_sgp/test_fd_alt_parsing.py`

- [ ] **4.1:** Same schema migration as Task 3 (write_to_duckdb at lines 270-318). Same atomic CREATE OR REPLACE.
- [ ] **4.2:** Same row-emission update (game_start_time = event.start_time, UTC).
- [ ] **4.3:** Add F3 whitelist entries (lines 61-76):
  ```python
  # F3 main
  "First 3 Innings Run Line":              ("F3", "main"),
  "First 3 Innings Total Runs":            ("F3", "main"),
  "First 3 Innings Money Line":            ("F3", "main"),
  # F3 alt
  "First 3 Innings Alternate Run Lines":   ("F3", "alternate_spreads"),
  "First 3 Innings Alternate Total Runs":  ("F3", "alternate_totals"),
  ```
  Replace placeholder names with exact strings from a live FD response if they differ.
- [ ] **4.4:** Add unicode/whitespace tolerance helper near regex defs:
  ```python
  _UNICODE_MINUS = "−–—"
  def _normalize_alt_name(raw: str) -> str:
      out = raw
      for ch in _UNICODE_MINUS:
          out = out.replace(ch, "-")
      return " ".join(out.split())
  ```
  Use it before each regex match in the alt branches. Skip (continue) on parse failure, don't fall through.
- [ ] **4.5:** Spread sign validation (same pattern as Task 3.6).
- [ ] **4.6:** Unit test at `mlb_sgp/test_fd_alt_parsing.py`:
  ```python
  from scraper_fanduel_singles import _FD_ALT_SPREAD_RE, _FD_ALT_TOTAL_RE, _normalize_alt_name

  def test_ascii_alt_spread():
      assert _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox +1.5")).group("line") == "+1.5"

  def test_unicode_minus():
      assert _FD_ALT_SPREAD_RE.match(_normalize_alt_name("Boston Red Sox −1.5")).group("line") == "-1.5"

  def test_extra_whitespace_total():
      assert _FD_ALT_TOTAL_RE.match(_normalize_alt_name("Over   (8.5)")).group("line") == "8.5"
  ```
  Run: `cd mlb_sgp && python -m pytest test_fd_alt_parsing.py -v`
- [ ] **4.7:** Smoke (same as Task 3.7, swap fd for dk).
- [ ] **4.8:** Commit.
  ```bash
  git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/test_fd_alt_parsing.py
  git commit -m "feat(fd-singles): TZ standardization + F3 + atomic + hardening"
  ```

---

### Task 5 (P3c): BFA scraper — schema migration (Group A)

**Files:**
- Modify: `bfa_odds/scraper*.py` (locate write path; structure should mirror DK/FD)

- [ ] **5.1:** Read BFA scraper to find write_to_duckdb / equivalent. Grep:
  ```bash
  grep -n "CREATE TABLE\|game_date\|game_time\|DELETE FROM mlb_odds" bfa_odds/scraper*.py
  ```
- [ ] **5.2:** BFA internally parses UTC `YYYY-MM-DD` + `HH:MM:SS` then writes back as strings. Skip the split: keep the parsed `datetime` (UTC-aware) and write it as `game_start_time TIMESTAMPTZ`.
- [ ] **5.3:** Atomic write (CREATE OR REPLACE pattern) — same as DK/FD if applicable.
- [ ] **5.4:** Smoke + verify `game_start_time` is TIMESTAMPTZ UTC.
- [ ] **5.5:** Commit.
  ```bash
  git add bfa_odds/scraper*.py
  git commit -m "feat(bfa): TZ standardization — write game_start_time TIMESTAMPTZ UTC"
  ```

---

### Task 6 (P3d): Kalshi scraper — schema migration (Group A)

**Files:**
- Modify: `kalshi_odds/scraper*.py`

- [ ] **6.1:** Locate Kalshi MLB scraper write path. Kalshi parses ISO UTC then re-stringifies to `%m/%d` + `%H:%M`; skip the round-trip.
- [ ] **6.2:** Write parsed datetime as `game_start_time TIMESTAMPTZ`.
- [ ] **6.3:** Verify kalshi_mlb_rfq bot read-only contract isn't broken — confirm bot reads `game_start_time` (or `commence_time` from mlb_odds_temp, not the Kalshi scraper DB directly):
  ```bash
  grep -rn "kalshi.duckdb\|kalshi_odds" /Users/callancapitolo/NFLWork/kalshi_mlb_rfq/
  ```
- [ ] **6.4:** Smoke + commit.

---

### Task 7 (P3e): WZ scraper — schema migration + TZ conversion (Group B)

**Files:**
- Modify: `wagerzon_odds/scraper.py`

- [ ] **7.1:** Apply Phase 1 audit's discovered WZ TZ (expected: `America/New_York`). Replace raw string slicing with:
  ```python
  from zoneinfo import ZoneInfo
  WZ_TZ = ZoneInfo("America/New_York")  # confirmed by P1 audit
  game_start_time = datetime.strptime(f"{gmdt} {gmtm}", "%Y%m%d %H:%M") \
                            .replace(tzinfo=WZ_TZ) \
                            .astimezone(timezone.utc)
  ```
- [ ] **7.2:** Schema migration to 17 cols (drop game_date/game_time, add game_start_time TIMESTAMPTZ).
- [ ] **7.3:** Atomic write.
- [ ] **7.4:** Smoke + commit.

---

### Task 8 (P3f): Hoop88 scraper — schema migration + TZ conversion (Group B)

Same pattern as Task 7, applying Phase 1's Hoop88 TZ finding.

---

### Task 9 (P3g): BKM scraper — schema migration + **PT→UTC fix**

**Why this is THE headline fix:** BKM was claimed ET; actually PT. ~half of BKM's MLB rows silently dropped by `.drop_past_games()` today.

**Files:**
- Modify: `bookmaker_odds/scraper.py`

- [ ] **9.1:** At the row-construction site (line ~427 `gmtm[:5]`):
  ```python
  from zoneinfo import ZoneInfo
  BKM_TZ = ZoneInfo("America/Los_Angeles")  # confirmed by P1 audit (was claimed ET, actually PT)
  # gmdt is the date string from BKM's API; gmtm is HH:MM:SS local
  bkm_local = datetime.strptime(f"{gmdt} {gmtm[:5]}", "%Y-%m-%d %H:%M").replace(tzinfo=BKM_TZ)
  game_start_time = bkm_local.astimezone(timezone.utc)
  ```
- [ ] **9.2:** Schema migration to 17 cols.
- [ ] **9.3:** Atomic write.
- [ ] **9.4:** Live verify: after smoke, confirm BKM's `game_start_time` for Mets@Nationals (or any current MLB game) matches Odds API `commence_time` to within seconds:
  ```bash
  python -c "
  import duckdb
  c = duckdb.connect('bookmaker_odds/bookmaker.duckdb', read_only=True)
  rows = c.execute(\"SELECT home_team, away_team, game_start_time FROM mlb_odds LIMIT 5\").fetchall()
  for r in rows: print(r)
  "
  ```
- [ ] **9.5:** Commit.
  ```bash
  git add bookmaker_odds/scraper.py
  git commit -m "fix(bookmaker): TZ standardization + PT→UTC conversion

  BKM's API returns Pacific Time naive strings; previously parsed as ET,
  dropping ~half of MLB rows via .drop_past_games(). Now writes
  game_start_time TIMESTAMPTZ (UTC) via explicit America/Los_Angeles
  conversion. This is THE bets-tab bug that made BKM pills empty on
  every game starting within ~3h of fetch."
  ```

---

### Task 10 (P3h): Bet105 scraper — schema migration + TZ verification

**Files:**
- Modify: `bet105_odds/scraper*.py`

- [ ] **10.1:** Apply Phase 1 audit's Bet105 finding. If audit confirms UTC, conversion is identity; just write `game_start_time TIMESTAMPTZ` (UTC). If audit shows offset, convert via discovered TZ.
- [ ] **10.2:** Schema migration + atomic write.
- [ ] **10.3:** Smoke + commit.

---

### Task 11 (P3i): wagerzon_specials in-place migration

**Why:** Persistent table — can't drop. Migration runs via `tools/migrate_scraper_schemas.py` (already coded in Task 2). This task updates the scraper that *writes* to it.

**Files:**
- Modify: `wagerzon_odds/scraper_specials.py`

- [ ] **11.1:** Update the row construction to write `game_start_time` (UTC TIMESTAMPTZ via WZ_TZ from P3e) into wagerzon_specials.
- [ ] **11.2:** Run the migration:
  ```bash
  python tools/migrate_scraper_schemas.py
  python -c "
  import duckdb
  c = duckdb.connect('wagerzon_odds/wagerzon.duckdb', read_only=True)
  print(c.execute(\"SELECT count(*), count(game_start_time) FROM wagerzon_specials\").fetchone())
  "
  ```
- [ ] **11.3:** Verify mlb_triple_play.R can still read wagerzon_specials with the new column (uses it via the ICU workaround today; Phase 4 will simplify).
- [ ] **11.4:** Commit.
  ```bash
  git add wagerzon_odds/scraper_specials.py
  git commit -m "feat(wz-specials): write game_start_time TIMESTAMPTZ UTC"
  ```

---

# PHASE 4 — R-side collapse

### Task 12: Delete 4 parsers + rewrite `.drop_past_games()` + simplify get_*_odds + drop ICU workaround

**Files:**
- Modify: `Answer Keys/Tools.R` (lines ~3608-3830 for get_*_odds, ~3618-3666 for parsers)
- Modify: `Answer Keys/MLB Answer Key/mlb_triple_play.R`
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`
- Modify: `Answer Keys/CBB Dashboard/cbb_dashboard_server.py`

- [ ] **12.1:** Delete the 4 parser functions. Grep first to verify call sites:
  ```bash
  grep -n "parse_wz_game_dt\|parse_bet105_game_dt\|parse_bfa_game_dt\|parse_iso_game_dt" "Answer Keys/Tools.R" "Answer Keys/MLB Answer Key/MLB.R" "Answer Keys/CBB Answer Key/CBB.R" 2>&1
  ```
- [ ] **12.2:** Rewrite `.drop_past_games()` in `Answer Keys/Tools.R`:
  ```r
  .drop_past_games <- function(raw, source_label = "scraper") {
    if (is.null(raw) || nrow(raw) == 0) return(raw)
    if (!"game_start_time" %in% names(raw)) {
      warning(sprintf("[%s] no game_start_time column — passing through", source_label))
      return(raw)
    }
    cutoff <- Sys.time() - 300  # 5 min grace
    before <- nrow(raw)
    out <- raw %>% filter(is.na(game_start_time) | game_start_time > cutoff)
    dropped <- before - nrow(out)
    if (dropped > 0) message(sprintf("[%s] dropped %d past-game rows", source_label, dropped))
    out
  }
  ```
- [ ] **12.3:** Update the 8 `get_*_odds` functions to drop the parser argument from `.drop_past_games(raw, .parse_*_game_dt, source_label = ...)`. Each becomes `.drop_past_games(raw, source_label = "dk")` etc.
- [ ] **12.4:** Update all SELECT statements that referenced `game_date`/`game_time` to use `game_start_time`:
  ```bash
  grep -n "game_date\|game_time" "Answer Keys/Tools.R" | head -40
  ```
  Audit each hit; replace where it's from a scraper table.
- [ ] **12.5:** `mlb_triple_play.R` — drop ICU workaround (lines ~98-105 install/load) and replace the `(NOW() AT TIME ZONE 'America/New_York')::TIMESTAMP` join (lines ~134-138) with plain `game_start_time > NOW()`.
- [ ] **12.6:** `mlb_dashboard_server.py` (lines 596-657) — delete the `if tzinfo is None: assume UTC` branches. TIMESTAMPTZ always returns aware.
- [ ] **12.7:** `cbb_dashboard_server.py` (lines 308-369) — same.
- [ ] **12.8:** Smoke: run MLB pipeline + dashboard. Confirm no errors about missing columns or parser functions.
- [ ] **12.9:** Commit.
  ```bash
  git add "Answer Keys/Tools.R" "Answer Keys/MLB Answer Key/mlb_triple_play.R" \
          "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" \
          "Answer Keys/CBB Dashboard/cbb_dashboard_server.py"
  git commit -m "refactor(R): delete parser zoo, simplify .drop_past_games

  All scraper-produced timestamps are now TIMESTAMPTZ UTC. Deleted
  .parse_wz_game_dt, .parse_bet105_game_dt, .parse_bfa_game_dt,
  .parse_iso_game_dt. .drop_past_games() filters on game_start_time
  directly. Dropped ICU workaround in mlb_triple_play.R. Dropped
  defensive UTC-assumption in dashboard servers."
  ```

---

# PHASE 5 — Bets-tab defenses

### Task 13: Coverage logging in `scraper_to_canonical`

**Why:** Silent inner_join drops in `odds_screen.R:397`. Even after TZ fix, unmapped team-name pairs disappear without a peep.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R::scraper_to_canonical`
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (callers — pass book_name)

- [ ] **13.1:** Add `book_name = NULL` parameter to `scraper_to_canonical()`.
- [ ] **13.2:** Anti-join + log before the inner_join:
  ```r
  dropped <- anti_join(raw, lookup, by = c("home_team", "away_team"))
  if (nrow(dropped) > 0) {
    dropped_pairs <- dropped %>% distinct(home_team, away_team) %>%
      mutate(pair = paste(home_team, "vs", away_team)) %>% pull(pair)
    message(sprintf("[scraper_to_canonical] %s: %d row(s) across %d game(s) dropped (no game_id match): %s",
                    book_name %||% "<book>", nrow(dropped), length(dropped_pairs),
                    paste(dropped_pairs, collapse = "; ")))
  }
  ```
- [ ] **13.3:** Update every MLB.R caller to pass `book_name = "dk"` etc.
- [ ] **13.4:** Smoke + commit.

---

### Task 14: `fetch_time` freshness gate in MLB.R

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R`

- [ ] **14.1:** Add top-of-script:
  ```r
  BOOK_STALENESS_CUTOFF_MIN <- 30
  ```
- [ ] **14.2:** Helper after per-book loads:
  ```r
  .drop_stale_book_rows <- function(df, book_label) {
    if (is.null(df) || nrow(df) == 0) return(df)
    if (!"fetch_time" %in% names(df)) return(df)
    cutoff <- Sys.time() - lubridate::minutes(BOOK_STALENESS_CUTOFF_MIN)
    fresh <- df %>% filter(is.na(fetch_time) | fetch_time >= cutoff)
    if ((dropped <- nrow(df) - nrow(fresh)) > 0)
      message(sprintf("[%s] dropped %d row(s) older than %dmin", book_label, dropped, BOOK_STALENESS_CUTOFF_MIN))
    fresh
  }
  dk_odds <- .drop_stale_book_rows(dk_odds, "dk")
  fd_odds <- .drop_stale_book_rows(fd_odds, "fd")
  # …repeat per book that exposes fetch_time
  ```
- [ ] **14.3:** Smoke + commit.

---

### Task 15: Pivot dedup guard

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:5970-5995`

- [ ] **15.1:** Insert before the pivot:
  ```r
  dup_check <- book_prices_long %>% count(bet_row_id, side, bookmaker, name = ".n") %>% filter(.n > 1)
  if (nrow(dup_check) > 0) {
    warning(sprintf("[bets-tab] %d duplicate (bet_row_id, side, bookmaker) — pivot keeps last only",
                    nrow(dup_check)))
    book_prices_long <- book_prices_long %>%
      distinct(bet_row_id, side, bookmaker, .keep_all = TRUE)
  }
  ```
- [ ] **15.2:** Smoke + commit.

---

### Task 16: `bet_row_id` collision detection

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (same region as Task 15)

- [ ] **16.1:** After loading `all_bets`:
  ```r
  if (nrow(all_bets) > 0) {
    bet_id_dups <- all_bets %>% count(bet_row_id, name = ".n") %>% filter(.n > 1)
    if (nrow(bet_id_dups) > 0)
      warning(sprintf("[bets-tab] %d bet_row_id collision(s) in all_bets", nrow(bet_id_dups)))
  }
  ```
- [ ] **16.2:** Smoke + commit (can combine with Task 15 into a single bets-tab-defenses commit if no other changes between).

---

# PHASE 6 — Parity harness

### Task 17: `tests/timezone_parity_test.py`

**Why:** This is the regression gate. Must pass post-refactor; must FAIL on `main` pre-refactor (proves the harness detects offsets).

**Files:**
- Create: `tests/timezone_parity_test.py`

- [ ] **17.1:** Implementation — runs against each scraper DB and Odds API:
  ```python
  """Regression gate: every scraper writes game_start_time matching Odds API."""
  import duckdb, json, os, sys
  from datetime import datetime, timezone
  from pathlib import Path
  from urllib.request import urlopen

  TOLERANCE_S = 60
  SCRAPERS = ["dk_odds/dk.duckdb", "fd_odds/fd.duckdb", "bfa_odds/bfa.duckdb",
              "kalshi_odds/kalshi.duckdb", "wagerzon_odds/wagerzon.duckdb",
              "hoop88_odds/hoop88.duckdb", "bookmaker_odds/bookmaker.duckdb",
              "bet105_odds/bet105.duckdb"]

  def main():
      key = os.environ["ODDS_API_KEY"]
      url = f"https://api.the-odds-api.com/v4/sports/baseball_mlb/odds?apiKey={key}&regions=us&markets=h2h"
      api = {(g["home_team"], g["away_team"]): datetime.fromisoformat(g["commence_time"].replace("Z","+00:00"))
             for g in json.loads(urlopen(url).read())}
      failures = []
      for db in SCRAPERS:
          if not Path(db).exists(): continue
          c = duckdb.connect(db, read_only=True)
          for ht, at, gst in c.execute("SELECT DISTINCT home_team, away_team, game_start_time FROM mlb_odds").fetchall():
              api_dt = api.get((ht, at))
              if api_dt is None: continue
              delta = abs((gst - api_dt).total_seconds())
              if delta > TOLERANCE_S:
                  failures.append((db, ht, at, delta))
          c.close()
      if failures:
          for db,ht,at,d in failures: print(f"FAIL {db}: {ht} vs {at} Δ={d:.0f}s", file=sys.stderr)
          sys.exit(1)
      print(f"PASS: all scrapers within {TOLERANCE_S}s of Odds API")

  if __name__ == "__main__": main()
  ```
- [ ] **17.2:** Run twice: first with all scrapers re-run on new code (`python tests/timezone_parity_test.py` → expect PASS), then mentally verify it would have FAILED on `main` (no need to actually run; the BKM-PT bug was the failure case).
- [ ] **17.3:** Commit.
  ```bash
  git add tests/timezone_parity_test.py
  git commit -m "test: timezone_parity_test.py regression gate

  Cross-references each scraper's game_start_time against Odds API
  commence_time within 60s tolerance. Run before any scraper-touching
  merge."
  ```

---

# PHASE 7 — Verification

### Task 18: End-to-end pipeline + dashboard smoke

- [ ] **18.1:** Stop bots if running:
  ```bash
  ps aux | grep -E "mm_bot|kalshi_mlb_rfq" | grep -v grep
  # If any: pkill -TERM -f mm_bot   # (avoid -9 unless necessary per CLAUDE.md)
  ```
- [ ] **18.2:** Full MLB pipeline:
  ```bash
  python "Answer Keys/run.py" mlb 2>&1 | tee /tmp/mlb_overhaul_smoke.log
  ```
  Inspect log for:
  - `[scraper_to_canonical] <book>: ... dropped` — record baseline per-book coverage.
  - `[<book>] dropped N row(s) older than 30min` — expect zero on healthy day.
  - `[dk_singles] WARN`, `[fd_singles] WARN` — expect zero.
  - `[bets-tab] duplicate / collision` — expect zero.
- [ ] **18.3:** Dashboard render:
  ```bash
  cd "Answer Keys/MLB Dashboard" && bash run.sh &
  sleep 5 && curl -s http://localhost:8083/ | grep -c "bet-card"
  ```
  Open in browser. Check explicitly:
  - **BKM pills are now present on Mets@Nationals-style cards** (the smoking-gun verification).
  - F3 pills exist on DK and FD where the model has F3 bets.
  - No cards have all-empty pills where one should exist.
- [ ] **18.4:** CBB pipeline smoke (lighter):
  ```bash
  python "Answer Keys/run.py" cbb 2>&1 | tail -20
  ```
- [ ] **18.5:** Parity test (Task 17):
  ```bash
  python tests/timezone_parity_test.py
  ```
  Must report PASS.
- [ ] **18.6:** If everything green, no commit needed (this is verification). If smoke reveals issues, loop back to Phase 3/4 to fix.

---

# PHASE 8 — Documentation

### Task 19: CLAUDE.md, scraper READMEs, memory

**Files:**
- Modify: `Answer Keys/CLAUDE.md` (Pitfalls #7 and #11)
- Modify: `/Users/callancapitolo/NFLWork/CLAUDE.md` (root, Housekeeping)
- Modify: 8 scraper READMEs (`*_odds/README.md`)
- Create: `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/timezone_standardization.md`

- [ ] **19.1:** `Answer Keys/CLAUDE.md` — delete Pitfall #7 (DuckDB naive TIMESTAMP vs NOW()) and Pitfall #11 (per-book game_time formats). Replace with:
  ```
  7. **Timestamp standard (2026-05-22)** — All scraper-produced
     timestamps are TIMESTAMPTZ in UTC. Use `column > NOW()` directly.
     `.drop_past_games()` filters game_start_time. Historical PBP
     timestamps (mlb_betting_pbp, cbb_betting_pbp) follow their own
     internal convention and are not produced by the in-scope scrapers.
  ```
- [ ] **19.2:** Root `CLAUDE.md` Housekeeping — add one bullet:
  ```
  6. **All new scrapers must write game_start_time TIMESTAMPTZ in UTC.**
     Do not introduce naive timestamp columns. See
     `tests/timezone_parity_test.py` for the regression gate.
  ```
- [ ] **19.3:** Each of 8 scraper READMEs — update Schema section to show `game_start_time TIMESTAMPTZ` instead of `game_date VARCHAR, game_time VARCHAR`.
- [ ] **19.4:** Create memory file `timezone_standardization.md`:
  ```markdown
  ---
  name: timezone-standardization
  description: Unified TIMESTAMPTZ UTC convention shipped 2026-05-22; replaces 4 R parsers and 3 per-book TZ conventions
  metadata:
    type: reference
  ---

  All 8 scrapers (DK, FD, BFA, Kalshi, WZ, Hoop88, BKM, Bet105) write
  `game_start_time TIMESTAMPTZ` in UTC. The 4 per-book parsers in
  Tools.R (`.parse_wz_game_dt`, `.parse_bet105_game_dt`, `.parse_bfa_game_dt`,
  `.parse_iso_game_dt`) were deleted. `.drop_past_games()` filters
  `game_start_time` directly.

  **Notable root cause:** BKM was tagged ET in code but actually wrote
  Pacific Time, silently dropping ~half of MLB rows via .drop_past_games
  pre-refactor. Confirmed via tests/timezone_audit.py.

  **Regression gate:** tests/timezone_parity_test.py. Must pass before
  any scraper-touching merge.

  Pre-refactor pitfalls (deleted from CLAUDE.md): #7 (naive TIMESTAMP vs
  NOW()) and #11 (per-book format zoo).
  ```
- [ ] **19.5:** Add pointer to memory MEMORY.md index:
  ```
  - [Timezone standardization](timezone_standardization.md) — Unified TIMESTAMPTZ UTC convention, regression gate, BKM-PT bug retro
  ```
- [ ] **19.6:** Commit.
  ```bash
  git add "Answer Keys/CLAUDE.md" CLAUDE.md \
          dk_odds/README.md fd_odds/README.md bfa_odds/README.md \
          kalshi_odds/README.md wagerzon_odds/README.md hoop88_odds/README.md \
          bookmaker_odds/README.md bet105_odds/README.md
  git commit -m "docs: timezone standardization — collapse 2 pitfalls into 1 standard

  CLAUDE.md Pitfalls #7 and #11 replaced with one rule. Root CLAUDE.md
  enforces the convention for new scrapers. Per-scraper READMEs updated.
  Memory entry added at timezone_standardization.md."
  ```

---

# Pre-merge executive engineer review checklist

Run before requesting merge approval:

- [ ] `python tests/timezone_parity_test.py` → PASS
- [ ] All 8 scraper `mlb_odds` tables have `game_start_time TIMESTAMPTZ` schema
- [ ] `wagerzon_specials` has `game_start_time TIMESTAMPTZ` populated for all rows
- [ ] `Answer Keys/Tools.R` has zero `.parse_*_game_dt` references
- [ ] `mlb_triple_play.R` has zero `INSTALL icu` / `LOAD icu` references
- [ ] MLB pipeline E2E green (Task 18)
- [ ] CBB pipeline E2E green
- [ ] Dashboard shows BKM pills on Mets@Nationals-style cards (the smoking-gun fix)
- [ ] No DuckDB connections left open (`on.exit(dbDisconnect(..., shutdown=TRUE))` in every R reader)
- [ ] `git diff main..HEAD --stat` shows only files in File map above
- [ ] CLAUDE.md Pitfalls #7 + #11 deleted; new unified entry present

---

# Out-of-scope follow-ups (tracked, not done here)

- NFL Draft portal scrapers (separate pipeline)
- Kalshi RFQ bot internal market DB timestamp surfaces
- Wagerzon-side equivalents of DK Task 5 (exact-prefix team match) — WZ uses different parsing, separate audit
- Per-pill staleness chips in dashboard UI
- Final DROP of wagerzon_specials.game_date/game_time after one stability cycle (kept as safety net)
- Centralized `mlb_odds` schema contract (8 independent DDLs today)
