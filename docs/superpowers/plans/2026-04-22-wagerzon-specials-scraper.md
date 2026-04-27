# Wagerzon Specials Scraper + Pricer DB Rewire Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the hardcoded daily tribble in `Answer Keys/mlb_triple_play.R` with a JSON-API scraper that captures all MLB specials into `wagerzon_odds/wagerzon.duckdb` → the pricer reads from DuckDB and prices every triple-play / grand-slam automatically.

**Architecture:** New Python module `wagerzon_odds/scraper_specials.py` reuses the existing `scraper_v2.login()` REST session to GET `NewScheduleHelper.aspx?WT=0&lg=4899`, parses the `result.listLeagues[0][0].Games` array, writes one row per prop into a new `wagerzon_specials` table. Pricer's main block replaces `todays_lines <- tribble(...)` with a DuckDB read filtered to `prop_type IN ('TRIPLE-PLAY','GRAND-SLAM')` for games starting in the next 12 hours.

**Tech Stack:** Python 3 (requests, duckdb — both already vendored in `wagerzon_odds/venv`), R (dplyr, duckdb, DBI). No new dependencies.

---

## Recon-confirmed schema (input)

```json
{
  "htm": "DODGERS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U5½)",
  "hnum": 757088,
  "idgm": 5635981,
  "gmdt": "20260427",  "gmtm": "19:40:00",
  "GameLines": [{ "odds": "190", "oddsh": "+190" }]
}
```

URL: `https://backend.wagerzon.com/wager/NewScheduleHelper.aspx?WT=0&lg=4899`
Confirmed today: 144 games in MLB - SPECIALS, of which 14 are priceable triple-plays + grand-slams (8 TP + 6 GS).

---

## File Structure

**Created:**
- `wagerzon_odds/scraper_specials.py` — fetcher + parser + DuckDB writer + CLI. Single file ~120 lines.
- `wagerzon_odds/tests/test_scraper_specials.py` — pytest with one fixture-based parser test (uses the recon JSON we already have on disk).

**Modified:**
- `Answer Keys/mlb_triple_play.R` — main block: replace tribble with DuckDB read.
- `Answer Keys/CLAUDE.md` — extend "Triple-Play Data Flow" with scraper input layer.
- `wagerzon_odds/README.md` — file inventory entry.

---

## Worktree & Version Control Plan

- Branch: `feature/wagerzon-specials-scraper` (already exists, already has the recon-fix commit `2074893`)
- Worktree: `.worktrees/wagerzon-specials-scraper` (already exists)
- Commits (one per task):
  1. Task 1 → `feat(wagerzon): parse_specials_json fixture-tested parser`
  2. Task 2 → `feat(wagerzon): scraper_specials writes wagerzon_specials table`
  3. Task 3 → `feat(wagerzon): scraper_specials CLI + main`
  4. Task 4 → `refactor(mlb): pricer reads wagerzon_specials from DuckDB`
  5. Task 5 → `docs: scraper_specials in CLAUDE.md + README`
- DuckDB files: scraper writes to `wagerzon_odds/wagerzon.duckdb` directly (this DB lives in the worktree's wagerzon_odds dir; the existing `wagerzon.duckdb` on main is a different file). For the regression check in Task 4, copy `Answer Keys/mlb.duckdb` from main into the worktree (read-only by pricer; remove after verification).
- Never merge to main without explicit user approval.

---

## Documentation Plan

Task 5 updates two files:
1. `Answer Keys/CLAUDE.md` — extend the "Triple-Play Data Flow" subsection with the scraper input layer (how `wagerzon_specials` feeds the pricer).
2. `wagerzon_odds/README.md` — add `scraper_specials.py` to the file inventory with one-line description.

---

## Pre-Merge Review Checklist

- **Data integrity**: scraper writes one row per (game_id, scraped_at). No duplicates within a snapshot. Table is append-only — pricer queries `MAX(scraped_at)` for current.
- **Resource safety**: `requests.Session` from `scraper_v2.login()` already manages cookies; one `duckdb.connect()` per scraper run with explicit `.close()`.
- **Graceful failure**: if login fails, scraper exits non-zero. If JSON shape is unexpected, scraper raises with the unexpected key name (no silent corruption).
- **No PII / secrets in logs**: scraper logs request URL + status + row count, never headers or response bodies.
- **Pricer regression**: today's 14 priceable props produce identical fair_odds when read from DB vs. when read from the manual tribble we just used (recorded in conversation: Marlins GS +627, Marlins TP +596, Rangers TP +373, etc.).
- **Idempotency**: re-running the scraper twice in 5 minutes produces 2 distinct snapshots (different `scraped_at` timestamps). Pricer always uses the latest.
- **No DuckDB symlinks** in the worktree (per project rule). Copy when needed, remove when done.

---

## Task 1: Pure parser — parse_specials_json

**Files:**
- Create: `wagerzon_odds/scraper_specials.py` (parser function only; no I/O yet)
- Create: `wagerzon_odds/tests/test_scraper_specials.py`

### - [ ] Step 1: Write failing tests

Create `wagerzon_odds/tests/test_scraper_specials.py`:

```python
"""Tests for wagerzon_odds.scraper_specials"""
import json
from pathlib import Path

# pytest discovers this via conftest path; for now just import directly
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from scraper_specials import parse_specials_json, extract_team_from_htm


# Fixture: today's recon JSON saved to disk during reconnaissance
RECON_FULL = Path(__file__).parent.parent / "recon_specials_full.json"


def test_extract_team_triple_play():
    assert extract_team_from_htm("DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)") == "DODGERS"

def test_extract_team_grand_slam_multiword():
    assert extract_team_from_htm("WHITE SOX GRAND-SLAM (SCR 1ST, 1H, GM & SCR U4½)") == "WHITE SOX"

def test_extract_team_no_prop_type_returns_none():
    assert extract_team_from_htm("MARLINS, ANGELS & YANKEES ALL WIN") is None

def test_extract_team_game_level_g_s_returns_none():
    # CHC-SDG G-S is game-level (cross-team) — no single subject team
    assert extract_team_from_htm("CHC-SDG G-S (GM O7½, 1H U4½, HITS U15½, 1ST HR 2R)") is None


def test_parse_specials_json_today_fixture():
    """Run against today's saved JSON. Should yield the 14 priceable rows."""
    if not RECON_FULL.exists():
        import pytest
        pytest.skip(f"{RECON_FULL} not present; run scraper_specials.py once first")
    raw = json.loads(RECON_FULL.read_text())
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    # Expect at least 14 priceable props (8 TP + 6 single-team GS)
    priceable = [r for r in rows if r["prop_type"] in ("TRIPLE-PLAY", "GRAND-SLAM") and r["team"]]
    assert len(priceable) >= 14, f"Got {len(priceable)} priceable, expected >=14"
    # Required fields present
    for r in priceable:
        assert r["scraped_at"] is not None
        assert r["sport"] == "mlb"
        assert r["league_id"] == 4899
        assert isinstance(r["odds"], int)
        assert r["description"]
        assert r["prop_type"] in ("TRIPLE-PLAY", "GRAND-SLAM")
        assert r["team"]


def test_odds_parsed_as_signed_int():
    raw = {"result": {"listLeagues": [[{"Description": "MLB - SPECIALS", "Games": [
        {"htm": "DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
         "hnum": 757041, "idgm": 5635934,
         "gmdt": "20260427", "gmtm": "19:40:00",
         "GameLines": [{"odds": "110", "oddsh": "-110"}]},
    ]}]]}}
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    dodgers = [r for r in rows if r["team"] == "DODGERS"][0]
    assert dodgers["odds"] == -110  # signed


def test_no_specials_section_returns_empty():
    raw = {"result": {"listLeagues": [[{"Description": "MLB - REGULAR", "Games": [
        {"htm": "Yankees", "hnum": 1, "idgm": 1,
         "gmdt": "20260427", "gmtm": "19:40:00",
         "GameLines": [{"odds": "100", "oddsh": "+100"}]},
    ]}]]}}
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    assert rows == []
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
wagerzon_odds/venv/bin/pip install pytest >/dev/null 2>&1 || true
wagerzon_odds/venv/bin/python -m pytest wagerzon_odds/tests/test_scraper_specials.py -v
```
Expected: failure with `ModuleNotFoundError: No module named 'scraper_specials'`.

### - [ ] Step 3: Implement parser

Create `wagerzon_odds/scraper_specials.py`:

```python
"""
Wagerzon MLB Specials Scraper

Fetches all MLB specials (lg=4899) from the NewScheduleHelper JSON endpoint
and writes one row per prop to wagerzon.duckdb / wagerzon_specials table.

The pricer (Answer Keys/mlb_triple_play.R) reads from this table.

Usage:
    python3 scraper_specials.py            # one-shot scrape, writes snapshot
    python3 scraper_specials.py --dry-run  # fetch + parse, don't write
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb
import requests

# Reuse the existing login + base-URL logic
sys.path.insert(0, str(Path(__file__).parent))
from scraper_v2 import login  # noqa: E402
from config import WAGERZON_BASE_URL  # noqa: E402

DB_PATH = Path(__file__).parent / "wagerzon.duckdb"
SPECIALS_URL_TPL = WAGERZON_BASE_URL + "/wager/NewSchedule" + "Helper.aspx?WT=0&lg={lg}"

# Regex used to extract the subject team. Anchors on known prop-type tokens
# so multi-word team names (WHITE SOX, RED SOX, BLUE JAYS) work. Extend the
# alternation when new prop types are added (e.g. DOUBLE-PLAY, MEGA).
PROP_TYPE_RE = re.compile(r"^(.+?)\s+(TRIPLE-PLAY|GRAND-SLAM)\b")
SECTION_RE = re.compile(r"\b(TRIPLE-PLAY|GRAND-SLAM)\b")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("scraper_specials")


def extract_team_from_htm(htm: str) -> Optional[str]:
    """Return the subject team (e.g. 'DODGERS', 'WHITE SOX') from a prop
    description, or None if the prop is cross-game / multi-team / unsupported.
    """
    if not htm:
        return None
    m = PROP_TYPE_RE.match(htm)
    if not m:
        return None
    team = m.group(1).strip()
    # Cross-game props ("CHC-SDG G-S", "MIA-LAD TRIPLE-PLAY") have a hyphenated
    # abbreviation in the team slot — those are game-level, not single-team.
    if "-" in team:
        return None
    return team


def extract_prop_type(htm: str) -> Optional[str]:
    m = SECTION_RE.search(htm or "")
    return m.group(1) if m else None


def parse_specials_json(payload: dict, sport: str, league_id: int) -> list[dict]:
    """Parse the NewScheduleHelper response into row dicts.

    Returns a list of dicts ready for DuckDB insert. Drops rows that aren't
    priceable single-team triple-plays or grand-slams (cross-game props,
    multi-team parlays, HR/hits combos, etc.). Filtering happens here so the
    DB stays focused; raw cross-game props can be added later by relaxing the
    filter without changing the schema.
    """
    scraped_at = datetime.now(timezone.utc).replace(microsecond=0)
    rows: list[dict] = []

    leagues_outer = payload.get("result", {}).get("listLeagues", [])
    if not leagues_outer or not leagues_outer[0]:
        return rows
    for league in leagues_outer[0]:
        # We only want MLB - SPECIALS leagues; other section headers (e.g.
        # MLB - REGULAR) are scraped by other tools.
        desc = (league.get("Description") or "").upper()
        if "SPECIALS" not in desc:
            continue
        for g in league.get("Games", []):
            htm = g.get("htm") or ""
            prop_type = extract_prop_type(htm)
            if prop_type not in ("TRIPLE-PLAY", "GRAND-SLAM"):
                continue
            team = extract_team_from_htm(htm)
            if team is None:
                continue  # cross-game props skipped
            game_lines = g.get("GameLines") or [{}]
            odds_str = (game_lines[0] or {}).get("odds") or ""
            oddsh    = (game_lines[0] or {}).get("oddsh") or ""
            try:
                # oddsh has explicit sign, e.g. "+190" or "-110"; odds is unsigned
                odds = int(oddsh) if oddsh else int(odds_str)
            except (TypeError, ValueError):
                log.warning("Could not parse odds for hnum=%s: oddsh=%r odds=%r",
                            g.get("hnum"), oddsh, odds_str)
                continue
            gmdt = g.get("gmdt") or ""
            gmtm = g.get("gmtm") or ""
            game_date = datetime.strptime(gmdt, "%Y%m%d").date() if gmdt else None
            try:
                # Wagerzon's gmtm appears to be local Eastern time — the scraper
                # stores it verbatim (no tz conversion). Pricer's 12-hour filter
                # uses commence_time from mlb_consensus_temp instead, which is
                # authoritative UTC.
                game_time = datetime.strptime(f"{gmdt} {gmtm}", "%Y%m%d %H:%M:%S")
            except ValueError:
                game_time = None
            rows.append({
                "scraped_at":      scraped_at,
                "sport":           sport,
                "league_id":       league_id,
                "game_id":         g.get("idgm"),
                "rotation_number": g.get("hnum"),
                "prop_type":       prop_type,
                "description":     htm,
                "team":            team,
                "line":            None,  # not applicable for SCR/period legs
                "odds":            odds,
                "game_date":       game_date,
                "game_time":       game_time,
            })
    return rows


def fetch_specials_json(lg: int) -> dict:
    s = requests.Session()
    login(s)
    url = SPECIALS_URL_TPL.format(lg=lg)
    r = s.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def write_rows(con: duckdb.DuckDBPyConnection, rows: list[dict]) -> None:
    """Idempotent: ensures table exists, then bulk inserts."""
    con.execute("""
        CREATE TABLE IF NOT EXISTS wagerzon_specials (
            scraped_at      TIMESTAMP,
            sport           VARCHAR,
            league_id       INTEGER,
            game_id         INTEGER,
            rotation_number INTEGER,
            prop_type       VARCHAR,
            description     VARCHAR,
            team            VARCHAR,
            line            DOUBLE,
            odds            INTEGER,
            game_date       DATE,
            game_time       TIMESTAMP
        )
    """)
    if not rows:
        log.info("No rows to write.")
        return
    con.executemany(
        "INSERT INTO wagerzon_specials VALUES "
        "(?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            (r["scraped_at"], r["sport"], r["league_id"], r["game_id"],
             r["rotation_number"], r["prop_type"], r["description"],
             r["team"], r["line"], r["odds"], r["game_date"], r["game_time"])
            for r in rows
        ],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lg", type=int, default=4899, help="Wagerzon league id (4899=MLB-SPECIALS)")
    ap.add_argument("--sport", default="mlb")
    ap.add_argument("--dry-run", action="store_true", help="Fetch + parse, don't write to DB")
    args = ap.parse_args()

    log.info("Fetching specials lg=%d", args.lg)
    payload = fetch_specials_json(args.lg)
    rows = parse_specials_json(payload, sport=args.sport, league_id=args.lg)
    log.info("Parsed %d priceable rows (TRIPLE-PLAY + GRAND-SLAM, single-team)", len(rows))

    if args.dry_run:
        for r in rows[:5]:
            print(r)
        log.info("Dry-run: not writing.")
        return 0

    con = duckdb.connect(str(DB_PATH))
    try:
        write_rows(con, rows)
        log.info("Wrote %d rows to %s/wagerzon_specials", len(rows), DB_PATH.name)
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
wagerzon_odds/venv/bin/python -m pytest wagerzon_odds/tests/test_scraper_specials.py -v
```
Expected: 6 tests pass (5 always, 1 skipped if recon_specials_full.json absent — but it IS present from earlier today).

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
git add wagerzon_odds/scraper_specials.py wagerzon_odds/tests/test_scraper_specials.py
git commit -m "feat(wagerzon): parse_specials_json fixture-tested parser"
```

---

## Task 2: Live scraper run + DB schema

**Files:**
- Modify: none (uses Task 1's `write_rows`)

### - [ ] Step 1: Run the scraper end-to-end

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
wagerzon_odds/venv/bin/python wagerzon_odds/scraper_specials.py
```
Expected log:
```
[INFO] Fetching specials lg=4899
[INFO] Logged in successfully
[INFO] Parsed 14 priceable rows (TRIPLE-PLAY + GRAND-SLAM, single-team)
[INFO] Wrote 14 rows to wagerzon.duckdb/wagerzon_specials
```

### - [ ] Step 2: Verify schema + content

```bash
wagerzon_odds/venv/bin/python -c "
import duckdb
con = duckdb.connect('wagerzon_odds/wagerzon.duckdb', read_only=True)
print('Schema:')
for c in con.execute(\"PRAGMA table_info('wagerzon_specials')\").fetchall():
    print(f'  {c[1]:<16} {c[2]}')
print()
print('Today rows:')
df = con.execute('''
  SELECT prop_type, team, odds, description
  FROM wagerzon_specials
  WHERE scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials)
  ORDER BY prop_type, team
''').fetchdf()
print(df.to_string(index=False))
print(f'\nTotal: {len(df)} rows')
"
```
Expected: 12-column schema, 14 rows (8 TRIPLE-PLAY + 6 GRAND-SLAM), all with team and odds populated.

### - [ ] Step 3: Re-run scraper, verify second snapshot is appended (idempotent)

```bash
wagerzon_odds/venv/bin/python wagerzon_odds/scraper_specials.py
wagerzon_odds/venv/bin/python -c "
import duckdb
con = duckdb.connect('wagerzon_odds/wagerzon.duckdb', read_only=True)
n_snaps = con.execute('SELECT COUNT(DISTINCT scraped_at) FROM wagerzon_specials').fetchone()[0]
n_rows  = con.execute('SELECT COUNT(*) FROM wagerzon_specials').fetchone()[0]
print(f'snapshots={n_snaps}  total_rows={n_rows}  per_snap={n_rows//n_snaps}')
"
```
Expected: `snapshots=2 total_rows=28 per_snap=14`. Confirms append-only + per-snapshot grouping work.

### - [ ] Step 4: Commit

(No code changes — committing Task 2 as a marker that the schema + first run were verified.)

```bash
git commit --allow-empty -m "feat(wagerzon): scraper_specials writes wagerzon_specials table

Verified: 14 rows per snapshot (8 TRIPLE-PLAY + 6 single-team GRAND-SLAM),
schema matches design spec, append-only behavior confirmed via 2 sequential
snapshots producing 28 total rows."
```

---

## Task 3: scraper_specials CLI integration

**Files:**
- Modify: `wagerzon_odds/scraper_specials.py` — add `--cron` mode (silent on success)

### - [ ] Step 1: Add a quiet / cron flag

Edit `wagerzon_odds/scraper_specials.py`. Find the `main()` function. Replace it with:

```python
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lg", type=int, default=4899,
                    help="Wagerzon league id (4899=MLB-SPECIALS)")
    ap.add_argument("--sport", default="mlb")
    ap.add_argument("--dry-run", action="store_true",
                    help="Fetch + parse, don't write to DB")
    ap.add_argument("--cron", action="store_true",
                    help="Quiet mode: log only WARN+, suitable for scheduled runs")
    args = ap.parse_args()

    if args.cron:
        logging.getLogger().setLevel(logging.WARNING)

    try:
        log.info("Fetching specials lg=%d", args.lg)
        payload = fetch_specials_json(args.lg)
        rows = parse_specials_json(payload, sport=args.sport, league_id=args.lg)
        log.info("Parsed %d priceable rows", len(rows))
    except Exception:
        log.exception("Scrape failed")
        return 2

    if args.dry_run:
        for r in rows[:5]:
            print(r)
        log.info("Dry-run: not writing.")
        return 0

    if not rows:
        log.warning("No priceable specials found — not writing snapshot")
        return 1

    con = duckdb.connect(str(DB_PATH))
    try:
        write_rows(con, rows)
        log.info("Wrote %d rows to %s/wagerzon_specials", len(rows), DB_PATH.name)
    finally:
        con.close()
    return 0
```

### - [ ] Step 2: Smoke-test the new flags

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
wagerzon_odds/venv/bin/python wagerzon_odds/scraper_specials.py --dry-run
wagerzon_odds/venv/bin/python wagerzon_odds/scraper_specials.py --cron
```
Expected: dry-run prints 5 sample rows + "Dry-run: not writing"; cron mode runs silently and exits 0.

### - [ ] Step 3: Commit

```bash
git add wagerzon_odds/scraper_specials.py
git commit -m "feat(wagerzon): scraper_specials CLI + main"
```

---

## Task 4: Pricer rewire to read DB

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — replace tribble with DuckDB read

### - [ ] Step 1: Replace the todays_lines tribble block

Edit `Answer Keys/mlb_triple_play.R`. Find the `todays_lines <- tribble(...)` block (~lines 46-59). Replace the ENTIRE tribble + comments with:

```r
  # Read today's posted specials from the wagerzon_specials scraper output.
  # Uses the most recent snapshot. Pricer is robust to empty / off-day cases:
  # if zero rows are posted, prints a clear message and exits.
  WZ_DB <- "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
  wz_con <- dbConnect(duckdb(), dbdir = path.expand(WZ_DB), read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(wz_con), error = function(e) NULL), add = TRUE)
  specials <- dbGetQuery(wz_con, "
    SELECT team, prop_type, description, odds AS book_odds
    FROM wagerzon_specials
    WHERE sport = 'mlb'
      AND prop_type IN ('TRIPLE-PLAY', 'GRAND-SLAM')
      AND scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials WHERE sport = 'mlb')
  ")
  dbDisconnect(wz_con)

  if (nrow(specials) == 0) {
    cat("No priceable specials found in wagerzon_specials. Run scraper_specials.py first.\n")
    quit(status = 0)
  }

  # Translate Wagerzon team names (UPPER) to Odds API canonical names so the
  # consensus join works. Map covers the 30 MLB teams.
  WZ_TO_CANONICAL <- c(
    "ANGELS"        = "Los Angeles Angels",
    "ASTROS"        = "Houston Astros",
    "ATHLETICS"     = "Athletics",
    "BLUE JAYS"     = "Toronto Blue Jays",
    "BRAVES"        = "Atlanta Braves",
    "BREWERS"       = "Milwaukee Brewers",
    "CARDINALS"     = "St. Louis Cardinals",
    "CUBS"          = "Chicago Cubs",
    "DBACKS"        = "Arizona Diamondbacks",
    "DIAMONDBACKS"  = "Arizona Diamondbacks",
    "DODGERS"       = "Los Angeles Dodgers",
    "GIANTS"        = "San Francisco Giants",
    "GUARDIANS"     = "Cleveland Guardians",
    "MARINERS"      = "Seattle Mariners",
    "MARLINS"       = "Miami Marlins",
    "METS"          = "New York Mets",
    "NATIONALS"     = "Washington Nationals",
    "ORIOLES"       = "Baltimore Orioles",
    "PADRES"        = "San Diego Padres",
    "PHILLIES"      = "Philadelphia Phillies",
    "PIRATES"       = "Pittsburgh Pirates",
    "RANGERS"       = "Texas Rangers",
    "RAYS"          = "Tampa Bay Rays",
    "RED SOX"       = "Boston Red Sox",
    "REDS"          = "Cincinnati Reds",
    "ROCKIES"       = "Colorado Rockies",
    "ROYALS"        = "Kansas City Royals",
    "TIGERS"        = "Detroit Tigers",
    "TWINS"         = "Minnesota Twins",
    "WHITE SOX"     = "Chicago White Sox",
    "YANKEES"       = "New York Yankees"
  )

  specials$canonical_team <- WZ_TO_CANONICAL[specials$team]
  unmapped <- specials[is.na(specials$canonical_team), ]
  if (nrow(unmapped) > 0) {
    warning(sprintf("Dropped %d specials with unmapped team names: %s",
                    nrow(unmapped),
                    paste(unique(unmapped$team), collapse = ", ")))
    specials <- specials[!is.na(specials$canonical_team), ]
  }

  # Resolve home/away by joining to mlb_consensus_temp (same approach used
  # for the prior tribble path). For each canonical_team, find the consensus
  # row where it's home OR away.
  con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)

  samples_df <- dbGetQuery(con,
    "SELECT game_id, home_margin, total_final_score,
            home_margin_f3, home_margin_f5, home_margin_f7,
            home_scored_first
     FROM mlb_game_samples")
  consensus  <- dbGetQuery(con,
    "SELECT id, home_team, away_team, commence_time FROM mlb_consensus_temp")

  if (!"home_scored_first" %in% names(samples_df)) {
    stop("mlb_game_samples is missing home_scored_first. Re-run MLB.R to regenerate.")
  }

  # 12-hour filter (same as before): mlb_consensus_temp carries multiple days
  consensus <- consensus %>%
    filter(commence_time > Sys.time() &
           commence_time < Sys.time() + 12 * 3600) %>%
    select(-commence_time)

  todays_lines <- specials %>%
    inner_join(consensus, by = c("canonical_team" = "home_team"),
               relationship = "many-to-many", keep = TRUE) %>%
    mutate(side = "home", target_team = team,
           home_team = canonical_team) %>%
    select(home_team, away_team, target_team, side, book_odds, description, id) %>%
    bind_rows(
      specials %>%
        inner_join(consensus, by = c("canonical_team" = "away_team"),
                   relationship = "many-to-many", keep = TRUE) %>%
        mutate(side = "away", target_team = team,
               away_team = canonical_team) %>%
        select(home_team, away_team, target_team, side, book_odds, description, id)
    )

  matched <- todays_lines  # consensus join already done; keep var name compatible
  n_matched <- nrow(matched)
  n_posted  <- nrow(specials)
  if (n_matched < n_posted) {
    warning(sprintf("Matched %d/%d posted specials. Some teams may have no game tonight.",
                    n_matched, n_posted))
  }
```

CRITICAL: this REPLACES the prior tribble + the `inner_join(consensus, ...)` block that followed it. Find both the tribble block AND the existing `# Join today's lines to their game_id` block AND the existing `dbGetQuery samples + consensus` block — all three are now part of this single replacement. Keep everything from `# Price each line` onward unchanged.

To make this concrete: in the existing file, the replacement spans roughly lines 46-95 (tribble + the original DB-read + consensus filter + inner_join + n_matched warn). The replacement body above includes all of that re-implemented to read from `wagerzon_specials` instead of the tribble.

### - [ ] Step 2: Copy mlb.duckdb into worktree for verification

```bash
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper/Answer Keys/mlb.duckdb"
```

### - [ ] Step 3: Run pricer end-to-end

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper/Answer Keys"
Rscript mlb_triple_play.R 2>&1 | tail -25
```
Expected: 14-row table with fair_odds populated. Compare to today's manual run:

| Team | Prop | Today's fair (manual) | Acceptable from DB |
|---|---|---:|---:|
| Marlins GS | away | +627 | ±5 |
| Marlins TP | away | +596 | ±5 |
| Rangers TP | home | +373 | ±5 |
| Angels GS | away | +292 | ±5 |
| White Sox TP | home | +297 | ±5 |
| Cubs TP | away | +223 | ±5 |
| Padres TP | home | +251 | ±5 |
| Yankees TP | away | +149 | ±3 |
| Angels TP | away | +220 | ±5 |
| Dodgers TP | home | +133 | ±3 |
| Rangers GS | home | +1555 | ±20 |
| White Sox GS | home | +1389 | ±20 |
| Dodgers GS | home | +439 | ±10 |
| Yankees GS | away | +805 | ±15 |

Identical values are expected — same samples, same legs, same compute_prop_fair. Any deviation > the bands above indicates a regression. Stop and investigate.

### - [ ] Step 4: Remove the copied mlb.duckdb

```bash
rm "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper/Answer Keys/mlb.duckdb"
```

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-specials-scraper
git add "Answer Keys/mlb_triple_play.R"
git commit -m "refactor(mlb): pricer reads wagerzon_specials from DuckDB"
```

---

## Task 5: Documentation

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — extend Triple-Play Data Flow with scraper input
- Modify: `wagerzon_odds/README.md` — file inventory entry

### - [ ] Step 1: Update Answer Keys/CLAUDE.md

In `Answer Keys/CLAUDE.md`, find the existing "Triple-Play Data Flow" subsection. Find the `parse_legs.R (generic prop parser)` block. Insert a new block ABOVE it (so the data flows top-down: scraper → parser → pricer):

```markdown

wagerzon_odds/scraper_specials.py (NEW)
  ├── Authenticated GET to NewScheduleHelper.aspx?WT=0&lg=4899
  ├── Filters to single-team TRIPLE-PLAY + GRAND-SLAM (cross-game props skipped)
  └── Appends one row per prop to wagerzon_specials in wagerzon.duckdb,
       keyed by scraped_at for snapshot history
```

Then update the `mlb_triple_play.R` description block in the same section. Replace its current content:

```
mlb_triple_play.R (standalone pricer)
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + scored_first)
  ├── For each tribble row: parse_legs(description) → compute_prop_fair(...)
  └── Prints fair odds vs book + edge per posted line, grouped by prop_type
```

with:

```
mlb_triple_play.R (standalone pricer)
  ├── Reads wagerzon_specials (latest scraped_at) for posted lines
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + scored_first)
  ├── Maps Wagerzon team names → Odds API canonical via WZ_TO_CANONICAL dict
  ├── Joins to consensus_temp for game_id + side (home/away)
  ├── For each row: parse_legs(description) → compute_prop_fair(...)
  └── Prints fair odds vs book + edge per posted line, grouped by prop_type
```

Add one new bullet to the existing bullet list:
- Adding new MLB teams (none today, but if expansion happens) requires updating `WZ_TO_CANONICAL` in the pricer.

### - [ ] Step 2: Update wagerzon_odds/README.md

In `wagerzon_odds/README.md`, find the file inventory section (or the Markets/Sports section, whichever lists scripts). Add this line under existing scripts:

```
- `scraper_specials.py` — captures all MLB - SPECIALS props (TRIPLE-PLAY, GRAND-SLAM)
  via the NewScheduleHelper JSON endpoint into `wagerzon_specials` table.
  Used by `Answer Keys/mlb_triple_play.R` for daily prop pricing.
```

### - [ ] Step 3: Commit

```bash
git add "Answer Keys/CLAUDE.md" wagerzon_odds/README.md
git commit -m "docs: scraper_specials in CLAUDE.md + README"
```

---

## Self-Review

**Spec coverage** (from `docs/superpowers/specs/2026-04-21-wagerzon-specials-scraper-design.md`):

| Spec section | Covered by |
|---|---|
| `wagerzon_specials` table schema | Task 1 (write_rows CREATE TABLE IF NOT EXISTS) |
| Authenticated REST scraper using existing scraper_v2.login() | Task 1 (fetch_specials_json) |
| Append-only snapshots keyed on scraped_at | Task 1 + Task 2 verify |
| Filter to triple-plays + grand-slams | Task 1 (parse_specials_json filter) |
| Pricer reads DB, drops tribble | Task 4 |
| Team-name canonical mapping | Task 4 (WZ_TO_CANONICAL dict) |
| Documentation updates | Task 5 |

**Placeholder scan**: no TBDs, no "implement later". Each step has either complete code or a concrete bash command with expected output.

**Type consistency**:
- `parse_specials_json` returns `list[dict]`; `write_rows` expects `list[dict]` — match.
- `odds` is `int` everywhere (signed in Python, INTEGER in DuckDB).
- `team` is uppercase Wagerzon name (str); R pricer translates via `WZ_TO_CANONICAL`.
- `prop_type` is one of `"TRIPLE-PLAY"` or `"GRAND-SLAM"` exactly (case-sensitive); R query uses same casing.
- `scraped_at` is timezone-aware UTC datetime; DuckDB stores as TIMESTAMP.

No issues found.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-22-wagerzon-specials-scraper.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — fresh subagent per task, review between tasks.

**2. Inline Execution** — execute tasks in this session with checkpoints.

**Which approach?**
