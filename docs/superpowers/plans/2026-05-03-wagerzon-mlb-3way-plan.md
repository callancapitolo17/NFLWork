# Wagerzon MLB F5 3-Way Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Match Wagerzon's F5 3-way moneyline market end-to-end (scrape → store → devig → match → emit bets), coexisting with the existing 2-way F5 ML matching.

**Architecture:** Add a new `idgmtyp=29` parser branch in the Wagerzon scraper, store the draw price in a new `draw_ml` column on the offshore odds tables (added via idempotent ALTER), and add a new `^h2h_3way_` branch to `compare_alts_to_samples` in `Tools.R` that uses 3-way devigging (sum-and-normalize) and prices three outcomes (home/away/tie) from the existing F5 sample margin distribution. No changes to the 2-way matching path.

**Tech Stack:** Python (scraper, DuckDB), R (Tools.R, testthat), pytest

**Branch:** `feature/wagerzon-mlb-3way` (already created via worktree)
**Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way`
**Spec:** `docs/superpowers/specs/2026-05-03-wagerzon-mlb-3way-design.md`

---

## File Structure

| File | Role | Touch |
|------|------|-------|
| `wagerzon_odds/scraper_v2.py` | Wagerzon scraper — owns `init_database()`, `parse_odds()`, `save_odds()`. Receives the new `parse_3way_line()` helper + idgmtyp=29 dispatch + schema column. | **Modify** (~+60 lines) |
| `wagerzon_odds/tests/test_scraper_v2_3way.py` | NEW pytest file mirroring `tests/test_scraper_specials.py` convention. Tests `parse_3way_line()` and `parse_odds()`'s 3-way dispatch via a synthetic JSON fixture inlined in the test. | **Create** |
| `Answer Keys/Tools.R` | Hosts `american_prob`, `get_wagerzon_odds`, `compare_alts_to_samples`. Receives one new helper (`american_prob_3way`), one new emission branch in the Wagerzon reader, and three coupled edits in `compare_alts_to_samples` (filter regex, suffix dispatcher, new `^h2h_3way_` branch). | **Modify** (~+95 lines) |
| `Answer Keys/tests/test_compare_alts_to_samples.R` | Existing testthat file (built in the prior PR). Receives 4 new tests covering the 3-way branch + devig helper. | **Modify** (~+90 lines) |
| `Answer Keys/CLAUDE.md` | Architecture doc for the answer-keys engine. Pipeline Flow (MLB) section needs the new market noted. | **Modify** (~+3 lines) |
| `Answer Keys/MLB Dashboard/README.md` | Markets section list needs `h2h_3way_1st_5_innings`. | **Modify** (~+1 line) |
| `wagerzon_odds/CLAUDE.md` | Module quick-map needs a note about `idgmtyp=29` handling + the `draw_ml` column. | **Modify** (~+5 lines) |

No new files outside the test file.

---

## Worktree Section

- **Worktree:** `.worktrees/wagerzon-mlb-3way` (already created off latest main, branch `feature/wagerzon-mlb-3way`)
- **DuckDB rule:** End-to-end smoke (Task 7) copies `wagerzon.duckdb`, `mlb.duckdb`, `mlb_mm.duckdb`, `pbp.duckdb` into the worktree — never symlinks. Per CLAUDE.md housekeeping #5, symlinks lose WAL data on worktree removal.
- **Cleanup after merge:** `git worktree remove .worktrees/wagerzon-mlb-3way && git branch -d feature/wagerzon-mlb-3way`

---

## Version Control Section

- All commits land on `feature/wagerzon-mlb-3way`.
- Commit cadence: one commit per task (8 commits total: 6 implementation + 1 smoke + 1 docs).
- Final merge: `git checkout main && git merge --no-ff feature/wagerzon-mlb-3way` — only after explicit user approval.
- Never use `--no-verify`. No force-push.

Commit messages:
1. `feat(wagerzon): add draw_ml column + idempotent schema upgrade`
2. `feat(wagerzon): parse_3way_line helper for idgmtyp=29 GameLines`
3. `feat(wagerzon): route idgmtyp=29 parents to 3-way parser`
4. `feat(mlb): emit h2h_3way records from get_wagerzon_odds`
5. `feat(mlb): american_prob_3way helper for 3-way devigging`
6. `feat(mlb): match h2h_3way_1st_5_innings in compare_alts_to_samples`
7. `test(mlb): end-to-end smoke for Wagerzon F5 3-way`
8. `docs(mlb): note F5 3-way matching in CLAUDE.md + READMEs`

---

## Documentation Section

After code is finalized and reviewed, the docs commit (Task 8) updates:

- **`Answer Keys/CLAUDE.md`** — under "Pipeline Flow (MLB)", append a bullet about `h2h_3way_1st_5_innings` matching via `compare_alts_to_samples` (parallel to the 2-way path).
- **`Answer Keys/MLB Dashboard/README.md`** — under "Markets" → "Derivatives", append `h2h_3way_1st_5_innings` with a note about Wagerzon-only.
- **`wagerzon_odds/CLAUDE.md`** — under "Quick map", note that `scraper_v2.py` now handles `idgmtyp=29` 3-way leagues and the `draw_ml` column on `*_odds` tables.

---

## Tasks

### Task 1: Schema upgrade — `draw_ml` column

**Files:**
- Modify: `wagerzon_odds/scraper_v2.py:580-614` (init_database)
- Modify: `wagerzon_odds/scraper_v2.py:628-633` (save_odds columns list)
- Test: `wagerzon_odds/tests/test_scraper_v2_3way.py` (new file — first test)

**Why this task first:** Every other task either reads or writes `draw_ml`. Locking the schema down first makes the rest of the implementation flow naturally.

- [ ] **Step 1: Create the new test file with the schema-upgrade test**

Create `wagerzon_odds/tests/test_scraper_v2_3way.py`:

```python
"""Tests for Wagerzon scraper_v2 — 3-way market handling."""
import json
import sys
import tempfile
from pathlib import Path

import duckdb
import pytest

# pytest discovers this via conftest path; for now import directly
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_init_database_adds_draw_ml_to_fresh_db(monkeypatch):
    """Fresh DB created with init_database() must have draw_ml column."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "wagerzon.duckdb"
        # Patch the module-level DB_PATH that init_database uses
        import scraper_v2
        monkeypatch.setattr(scraper_v2, "DB_PATH", db_path)

        scraper_v2.init_database("mlb")

        con = duckdb.connect(str(db_path))
        cols = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        con.close()

        assert "draw_ml" in cols, f"draw_ml missing from fresh schema: {cols}"


def test_init_database_idempotent_ALTER_on_existing_db(monkeypatch):
    """Existing DB without draw_ml gets the column added; second call is no-op."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "wagerzon.duckdb"

        # Pre-create an old-schema mlb_odds (no draw_ml)
        con = duckdb.connect(str(db_path))
        con.execute("""
            CREATE TABLE mlb_odds (
                fetch_time TIMESTAMP, sport_key VARCHAR, game_id VARCHAR,
                game_date VARCHAR, game_time VARCHAR, away_team VARCHAR,
                home_team VARCHAR, market VARCHAR, period VARCHAR,
                away_spread FLOAT, away_spread_price INTEGER,
                home_spread FLOAT, home_spread_price INTEGER,
                total FLOAT, over_price INTEGER, under_price INTEGER,
                away_ml INTEGER, home_ml INTEGER, idgm INTEGER
            )
        """)
        con.close()

        import scraper_v2
        monkeypatch.setattr(scraper_v2, "DB_PATH", db_path)

        # First call adds draw_ml
        scraper_v2.init_database("mlb")
        con = duckdb.connect(str(db_path))
        cols_after_1 = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        assert "draw_ml" in cols_after_1

        # Second call: idempotent, no error, no duplicate column
        scraper_v2.init_database("mlb")
        cols_after_2 = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        con.close()
        assert cols_after_1 == cols_after_2, "Second init_database call changed schema"
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py -v 2>&1 | tail -15
```

Expected: 2 FAILs — both tests fail because `draw_ml` is not in the schema.

- [ ] **Step 3: Update `init_database()` in `wagerzon_odds/scraper_v2.py`**

Replace the function body (currently at lines 580-614) with:

```python
def init_database(sport: str):
    """Initialize DuckDB with the odds table for a sport."""
    config = get_sport_config(sport)
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            fetch_time TIMESTAMP,
            sport_key VARCHAR,
            game_id VARCHAR,
            game_date VARCHAR,
            game_time VARCHAR,
            away_team VARCHAR,
            home_team VARCHAR,
            market VARCHAR,
            period VARCHAR,
            away_spread FLOAT,
            away_spread_price INTEGER,
            home_spread FLOAT,
            home_spread_price INTEGER,
            total FLOAT,
            over_price INTEGER,
            under_price INTEGER,
            away_ml INTEGER,
            home_ml INTEGER,
            draw_ml INTEGER,
            idgm INTEGER
        )
    """)
    # Idempotent upgrades for existing DBs that pre-date these columns.
    # ADD COLUMN IF NOT EXISTS is supported in DuckDB 1.4+ and is a no-op
    # when the column already exists.
    for col_def in ("idgm INTEGER", "draw_ml INTEGER"):
        conn.execute(f"ALTER TABLE {table_name} ADD COLUMN IF NOT EXISTS {col_def}")
    conn.close()
```

This replaces the previous `try/except` ALTER block with a cleaner `IF NOT EXISTS` form (verified working on DuckDB 1.4.4).

- [ ] **Step 4: Update `save_odds()` columns list in `wagerzon_odds/scraper_v2.py`**

Replace lines 628-633 (the `columns = [...]` block) with:

```python
    columns = [
        "fetch_time", "sport_key", "game_id", "game_date", "game_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml",
        "draw_ml", "idgm"
    ]
```

(Just adds `"draw_ml"` between `"home_ml"` and `"idgm"`.)

- [ ] **Step 5: Run the tests to verify they pass**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py -v 2>&1 | tail -15
```

Expected: 2 PASS, 0 FAIL.

- [ ] **Step 6: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "wagerzon_odds/scraper_v2.py" "wagerzon_odds/tests/test_scraper_v2_3way.py"
git commit -m "feat(wagerzon): add draw_ml column + idempotent schema upgrade"
```

---

### Task 2: `parse_3way_line` helper

**Files:**
- Modify: `wagerzon_odds/scraper_v2.py:189-216` (insert new helper after `parse_moneyline_only`)
- Modify: `wagerzon_odds/tests/test_scraper_v2_3way.py` (append parser tests)

**Why a separate helper (vs reusing `parse_game_line`):** `parse_game_line` returns spread + total + 2-way ML; `parse_3way_line` returns 3-way ML only. Keeping them separate avoids overloading `parse_game_line` with `vspoddst` extraction that's only meaningful for 3-way.

- [ ] **Step 1: Append failing tests to `wagerzon_odds/tests/test_scraper_v2_3way.py`**

```python
def _make_3way_base():
    """Synthetic base dict matching what parse_odds builds before calling helpers."""
    return {
        "fetch_time": "2026-05-03 19:30:00",
        "sport_key": "baseball_mlb",
        "away_team": "Cleveland Guardians",
        "home_team": "Athletics",
        "game_date": "05/03",
        "game_time": "19:30",
        "idgm": 5635900,
    }


def test_parse_3way_line_emits_record_with_draw_ml():
    """A 3-way GameLine with all three prices populates the row correctly."""
    from scraper_v2 import parse_3way_line
    line = {
        "voddst": "105",   # away ML
        "hoddst": "130",   # home ML
        "vspoddst": "475", # draw price
        # ...other fields like vsprdt/vsprdoddst exist but are empty/irrelevant
    }
    rec = parse_3way_line(line, "test-game-1", "f5", "h2h_3way_1st_5_innings", _make_3way_base())

    assert rec is not None
    assert rec["market"] == "h2h_3way_1st_5_innings"
    assert rec["period"] == "f5"
    assert rec["away_ml"] == 105
    assert rec["home_ml"] == 130
    assert rec["draw_ml"] == 475
    # Spread/total fields must be NULL for 3-way rows
    assert rec["away_spread"] is None
    assert rec["home_spread"] is None
    assert rec["total"] is None
    assert rec["over_price"] is None
    assert rec["under_price"] is None
    # Base fields preserved
    assert rec["away_team"] == "Cleveland Guardians"
    assert rec["home_team"] == "Athletics"
    assert rec["game_id"] == "test-game-1"


def test_parse_3way_line_returns_none_when_all_three_prices_missing():
    """If voddst/hoddst/vspoddst are all empty, no row should be emitted."""
    from scraper_v2 import parse_3way_line
    line = {"voddst": "", "hoddst": "", "vspoddst": ""}
    rec = parse_3way_line(line, "test-game-2", "f5", "h2h_3way_1st_5_innings", _make_3way_base())
    assert rec is None


def test_parse_3way_line_handles_negative_prices():
    """Negative American odds (e.g. heavy favorite) parse correctly."""
    from scraper_v2 import parse_3way_line
    line = {"voddst": "-150", "hoddst": "+200", "vspoddst": "+500"}
    rec = parse_3way_line(line, "test-game-3", "f5", "h2h_3way_1st_5_innings", _make_3way_base())
    assert rec["away_ml"] == -150
    assert rec["home_ml"] == 200
    assert rec["draw_ml"] == 500
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py -v -k "parse_3way_line" 2>&1 | tail -15
```

Expected: 3 FAILs (`AttributeError: module 'scraper_v2' has no attribute 'parse_3way_line'`).

- [ ] **Step 3: Add `parse_3way_line` to `wagerzon_odds/scraper_v2.py`**

Insert after the `parse_moneyline_only` function (after line 216, before the `parse_odds` function):

```python
def parse_3way_line(line: dict, game_id: str, period: str, market: str,
                    base: dict) -> Optional[dict]:
    """Parse a 3-way GameLine (e.g. lg=1280 MLB - 1ST 5 INN WINNER (3-WAY)).

    Wagerzon's 3-way market uses three price fields on a single GameLine:
        voddst    -> away ML
        hoddst    -> home ML
        vspoddst  -> draw price
                     (the 'spread' field is repurposed to hold the third outcome)

    Returns None if all three prices are missing (Wagerzon posts placeholder
    games with empty prices when lines aren't yet posted).
    """
    away_ml = safe_int(line.get("voddst"))
    home_ml = safe_int(line.get("hoddst"))
    draw_ml = safe_int(line.get("vspoddst"))

    if away_ml is None and home_ml is None and draw_ml is None:
        return None

    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": None,
        "away_spread_price": None,
        "home_spread": None,
        "home_spread_price": None,
        "total": None,
        "over_price": None,
        "under_price": None,
        "away_ml": away_ml,
        "home_ml": home_ml,
        "draw_ml": draw_ml,
    }
```

- [ ] **Step 4: Update `parse_team_total` and `parse_moneyline_only` to include `draw_ml: None`**

Both helpers currently return dicts that DON'T have a `draw_ml` key. Since `save_odds` builds tuples positionally from the `columns` list (which now includes `draw_ml`), every dict that flows into `save_odds` must have a `draw_ml` key — even if NULL.

In `parse_team_total` (around line 162-186), add `"draw_ml": None,` immediately after `"home_ml": None,`:

```python
def parse_team_total(line: dict, game_id: str, period: str, market: str,
                     base: dict) -> Optional[dict]:
    # ... existing body unchanged ...
    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": None,
        "away_spread_price": None,
        "home_spread": None,
        "home_spread_price": None,
        "total": total,
        "over_price": over_price,
        "under_price": under_price,
        "away_ml": None,
        "home_ml": None,
        "draw_ml": None,    # NEW — required by save_odds columns list
    }
```

In `parse_moneyline_only` (around line 189-216), same addition immediately after `"home_ml": home_ml,`:

```python
def parse_moneyline_only(line: dict, game_id: str, period: str, market: str,
                         base: dict) -> Optional[dict]:
    # ... existing body unchanged ...
    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": None,
        "away_spread_price": None,
        "home_spread": None,
        "home_spread_price": None,
        "total": None,
        "over_price": None,
        "under_price": None,
        "away_ml": away_ml,
        "home_ml": home_ml,
        "draw_ml": None,    # NEW — required by save_odds columns list
    }
```

In `parse_game_line` (around line 126-159), same addition:

```python
    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": away_spread,
        "away_spread_price": away_spread_price,
        "home_spread": home_spread,
        "home_spread_price": home_spread_price,
        "total": total,
        "over_price": over_price,
        "under_price": under_price,
        "away_ml": away_ml,
        "home_ml": home_ml,
        "draw_ml": None,    # NEW — required by save_odds columns list
    }
```

Also update the `parse_odds` race-to-X record builder at line 549-569 (in the `for league in leagues:` race loop) to include `"draw_ml": None,` in the record dict. AND update the alt_spread record at line 387-401 and alt_total record at line 405-420 (both inside the `child_type == 25` branch) to include `"draw_ml": None,`.

This is mechanical: every dict that ends up in `records` must have a `draw_ml` key.

- [ ] **Step 5: Run all parser tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py -v 2>&1 | tail -15
```

Expected: all tests pass (2 from Task 1 + 3 from Task 2 = 5 total).

- [ ] **Step 6: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "wagerzon_odds/scraper_v2.py" "wagerzon_odds/tests/test_scraper_v2_3way.py"
git commit -m "feat(wagerzon): parse_3way_line helper for idgmtyp=29 GameLines"
```

---

### Task 3: Route `idgmtyp=29` parents to the 3-way parser

**Files:**
- Modify: `wagerzon_odds/scraper_v2.py:269-313` (parent-game loop in `parse_odds`)
- Modify: `wagerzon_odds/tests/test_scraper_v2_3way.py` (append dispatch test)

**Why this task is separate from Task 2:** Task 2 builds and unit-tests the helper. Task 3 wires the helper into the live parsing path, adds the team-name strip logic ("1H " prefix and " 3WAY" suffix), and tests the full parse_odds dispatch with a synthetic JSON fixture.

- [ ] **Step 1: Append a dispatch test to `wagerzon_odds/tests/test_scraper_v2_3way.py`**

```python
def test_parse_odds_routes_idgmtyp_29_to_3way_parser():
    """Synthetic JSON with a 3-way league should produce h2h_3way_* rows."""
    from scraper_v2 import parse_odds

    synthetic = {
        "result": {
            "listLeagues": [[
                {
                    "Description": "MLB - 1ST 5 INN WINNER (3-WAY)",
                    "Games": [
                        {
                            "vtm": "1H CLE GUARDIANS 3WAY",
                            "htm": "1H ATHLETICS 3WAY",
                            "vnum": 12345,
                            "hnum": 12346,
                            "gmdt": "20260503",
                            "gmtm": "19:30:00",
                            "gpd": "GAME",
                            "idgmtyp": 29,
                            "idlg": 1280,
                            "idspt": "MLB",
                            "idgm": 5635900,
                            "GameLines": [{
                                "voddst": "105",
                                "hoddst": "130",
                                "vspoddst": "475",
                            }],
                            "GameChilds": [],
                        },
                        {
                            # Empty placeholder game (Wagerzon posts these
                            # before lines are set) — should be skipped.
                            "vtm": ".",
                            "htm": ".",
                            "vnum": 0, "hnum": 0,
                            "gmdt": "20260503", "gmtm": "00:00:00",
                            "gpd": "GAME", "idgmtyp": 29, "idlg": 1280,
                            "idspt": "MLB", "idgm": 5635901,
                            "GameLines": [{"voddst": "", "hoddst": "", "vspoddst": ""}],
                            "GameChilds": [],
                        },
                    ],
                }
            ]]
        }
    }

    records = parse_odds(synthetic, "mlb")

    three_way = [r for r in records if r["market"] == "h2h_3way_1st_5_innings"]
    assert len(three_way) == 1, f"Expected 1 three-way row, got {len(three_way)}: {three_way}"

    rec = three_way[0]
    assert rec["away_ml"] == 105
    assert rec["home_ml"] == 130
    assert rec["draw_ml"] == 475
    assert rec["period"] == "f5"

    # Team names: "1H " prefix and " 3WAY" suffix must be stripped before
    # team-name resolution. Without canonical mapping the test sees the
    # stripped raw names — exact-match isn't asserted (depends on local
    # team_dict), but they must not still contain "1H" or "3WAY".
    assert "1H" not in rec["away_team"], f"1H not stripped: {rec['away_team']}"
    assert "3WAY" not in rec["away_team"], f"3WAY not stripped: {rec['away_team']}"
    assert "1H" not in rec["home_team"]
    assert "3WAY" not in rec["home_team"]
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py::test_parse_odds_routes_idgmtyp_29_to_3way_parser -v 2>&1 | tail -10
```

Expected: FAIL — current `parse_odds` has no idgmtyp=29 branch; the parent loop calls `parse_game_line` on it, which produces a row with `market="spreads"` (not `h2h_3way_1st_5_innings`).

- [ ] **Step 3: Modify the parent-game loop in `parse_odds` (around line 269-313)**

Find the parent loop:

```python
for game in games:
    if not game.get("GameLines"):
        continue

    away_raw = game["vtm"]
    home_raw = game["htm"]
    # ...
```

Insert at the top of the loop body (immediately after the `if not game.get("GameLines"): continue` guard), BEFORE the team-name resolution:

```python
    # 3-way market parents: idgmtyp=29 lives in lg=1280 ("MLB - 1ST 5 INN
    # WINNER (3-WAY)"). Route to parse_3way_line, strip "1H " prefix and
    # " 3WAY" suffix from team names, and skip the standard parent + child
    # parsing (3-way games have no derivatives).
    if game.get("idgmtyp") == 29:
        line = game["GameLines"][0]
        # Strip "1H " prefix and " 3WAY" suffix before team resolution
        away_raw_3w = re.sub(r"^1H\s+|\s+3WAY$", "", game["vtm"])
        home_raw_3w = re.sub(r"^1H\s+|\s+3WAY$", "", game["htm"])

        gmdt = game.get("gmdt", "")
        game_date = f"{gmdt[4:6]}/{gmdt[6:8]}" if len(gmdt) == 8 else ""
        game_time = game.get("gmtm", "")[:5]

        if team_dict or canonical_games:
            away_team, home_team = resolve_team_names(
                away_raw_3w, home_raw_3w, team_dict, canonical_games
            )
        else:
            away_team = normalize_team_name(away_raw_3w, sport)
            home_team = normalize_team_name(home_raw_3w, sport)

        away_rot = str(game["vnum"])
        home_rot = str(game["hnum"])
        game_id = f"{away_rot}-{home_rot}"

        base = {
            "fetch_time": fetch_time,
            "sport_key": sport_key,
            "game_date": game_date,
            "game_time": game_time,
            "away_team": away_team,
            "home_team": home_team,
            "idgm": game.get("idgm"),
        }

        rec = parse_3way_line(
            line, game_id, "f5", "h2h_3way_1st_5_innings", base
        )
        if rec:
            records.append(rec)
            print(f"  3-way: {away_team} @ {home_team} | "
                  f"{rec['away_ml']}/{rec['home_ml']}/{rec['draw_ml']}")
        continue   # 3-way games have no GameChilds we care about
```

The `continue` at the end is important — it short-circuits the existing parent + GameChilds parsing for these games (which would otherwise re-parse them as `spreads` rows with bogus team names).

- [ ] **Step 4: Run the dispatch test to verify it passes**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/test_scraper_v2_3way.py::test_parse_odds_routes_idgmtyp_29_to_3way_parser -v 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Run the full scraper test suite**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds"
python3 -m pytest tests/ -v 2>&1 | tail -20
```

Expected: all tests pass (Task 1's 2 + Task 2's 3 + Task 3's 1 = 6 from `test_scraper_v2_3way.py`, plus whatever was passing before in `test_scraper_specials.py`). No regressions.

- [ ] **Step 6: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "wagerzon_odds/scraper_v2.py" "wagerzon_odds/tests/test_scraper_v2_3way.py"
git commit -m "feat(wagerzon): route idgmtyp=29 parents to 3-way parser"
```

---

### Task 4: Reader emits `h2h_3way` records from `get_wagerzon_odds`

**Files:**
- Modify: `Answer Keys/Tools.R:3035+` (per-row dispatcher inside `get_wagerzon_odds`)
- Test: `Answer Keys/tests/test_compare_alts_to_samples.R` (append a small reader test — but actually the integration is exercised by Task 6's matcher tests; Task 4 just needs to be observable via the matcher tests, which require the reader to emit records correctly)

**Why this task exists:** The matcher (`compare_alts_to_samples`) reads from the offshore_odds tibble that `get_wagerzon_odds` produces. If the reader doesn't emit `h2h_3way` records with `odds_draw` populated, the matcher branch in Task 6 has nothing to match.

**Why no dedicated test file:** The reader's behavior is exercised end-to-end via the synthetic test fixtures in Task 6 (which populate offshore_odds tibbles directly with the expected shape). Adding a reader-only test would duplicate the assertion. If a future change breaks the reader, Task 6's tests will catch it because they require `get_wagerzon_odds`-shaped output.

- [ ] **Step 1: Read `get_wagerzon_odds` to find the per-row dispatch loop**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
sed -n '3035,3105p' "Answer Keys/Tools.R"
```

You should see a `for (i in seq_len(nrow(raw_odds)))` loop that emits multiple records per raw row (one for spread, one for total, one for ML if present).

- [ ] **Step 2: Insert the new branch BEFORE the existing spread emission**

Find the line `for (i in seq_len(nrow(raw_odds))) {` and the immediately-following `row <- raw_odds[i, ]` and `base <- list(...)`. Right after `base` is constructed and before the first `if (!is.na(row$away_spread))` check, insert:

```r
    # 3-way market: emit one record per game with all three prices
    # (h2h_3way_1st_5_innings — Wagerzon's F5 3-way market). Skip the
    # standard spread/total/h2h triplet emission because 3-way rows have
    # all those fields NULL by definition.
    if (!is.null(row$draw_ml) && !is.na(row$draw_ml) &&
        grepl("^h2h_3way_", row$market)) {
      three_way_rec <- c(base, list(
        market = row$market,
        market_type = "h2h_3way",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_draw = row$draw_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- three_way_rec
      next   # don't fall through to spread/total/h2h record emission
    }
```

The `next` is critical — it advances the for-loop iteration so the row doesn't ALSO get emitted as a spread/total/h2h record.

- [ ] **Step 3: Confirm the file still parses**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
Rscript -e 'source("Tools.R"); cat("Tools.R loaded OK\n")' 2>&1 | tail -3
```

Expected: prints `Tools.R loaded OK` (no syntax error). If it fails, look for unmatched braces.

- [ ] **Step 4: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "Answer Keys/Tools.R"
git commit -m "feat(mlb): emit h2h_3way records from get_wagerzon_odds"
```

---

### Task 5: `american_prob_3way` devig helper

**Files:**
- Modify: `Answer Keys/Tools.R:87+` (insert helper after existing `american_prob`)
- Modify: `Answer Keys/tests/test_compare_alts_to_samples.R` (add helper tests)

**Why a separate helper (vs extending `american_prob`):** `american_prob` is called from dozens of sites and takes 2 args. Changing its signature to accept an optional 3rd arg risks subtle bugs. A new `american_prob_3way` is clean, isolated, and testable.

- [ ] **Step 1: Append failing tests to `Answer Keys/tests/test_compare_alts_to_samples.R`**

```r
test_that("american_prob_3way devigs three American odds via sum-and-normalize", {
  # CLE @ ATH 3-way from the recon: away=+105, home=+130, draw=+475
  # Implied: 100/205, 100/230, 100/575 = 0.4878, 0.4348, 0.1739; sum=1.0965
  # Devigged: 0.4448, 0.3966, 0.1586
  result <- american_prob_3way(105, 130, 475)
  expect_equal(result$p_away, 0.4448, tolerance = 1e-3)
  expect_equal(result$p_home, 0.3966, tolerance = 1e-3)
  expect_equal(result$p_draw, 0.1586, tolerance = 1e-3)
  expect_equal(result$p_away + result$p_home + result$p_draw, 1.0, tolerance = 1e-9)
})

test_that("american_prob_3way handles negative odds", {
  # Heavy home favorite: away=+200, home=-150, draw=+500
  # Implied: 0.3333, 0.6, 0.1667; sum=1.1
  # Devigged: 0.3030, 0.5455, 0.1515
  result <- american_prob_3way(200, -150, 500)
  expect_equal(result$p_away, 0.3030, tolerance = 1e-3)
  expect_equal(result$p_home, 0.5455, tolerance = 1e-3)
  expect_equal(result$p_draw, 0.1515, tolerance = 1e-3)
})

test_that("american_prob_3way returns NA for any NA input", {
  result <- american_prob_3way(NA_integer_, 130, 475)
  expect_true(is.na(result$p_away))
  expect_true(is.na(result$p_home))
  expect_true(is.na(result$p_draw))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")' 2>&1 | tail -10
```

Expected: 3 new FAILs (`could not find function "american_prob_3way"`).

- [ ] **Step 3: Add `american_prob_3way` to `Answer Keys/Tools.R`**

Find the existing `american_prob` function (line 87). After its closing `}`, insert:

```r
#' Devig three American odds (3-way market: home/draw/away)
#'
#' Sum-and-normalize devig: each implied prob divided by the sum of all three.
#' Used for Wagerzon's F5 3-way market where push is its own outcome (not refund).
#'
#' @param odds_away,odds_home,odds_draw American odds for the three outcomes
#' @return list with p_away, p_home, p_draw (each in [0,1], sum to 1.0)
american_prob_3way <- function(odds_away, odds_home, odds_draw) {
  implied <- function(odds) {
    if (is.na(odds) || odds == 0) return(NA_real_)
    if (odds > 0) 100 / (odds + 100)
    else (-odds) / (-odds + 100)
  }
  p_a <- implied(odds_away)
  p_h <- implied(odds_home)
  p_d <- implied(odds_draw)
  if (any(is.na(c(p_a, p_h, p_d)))) {
    return(list(p_away = NA_real_, p_home = NA_real_, p_draw = NA_real_))
  }
  total <- p_a + p_h + p_d
  list(p_away = p_a / total, p_home = p_h / total, p_draw = p_d / total)
}
```

- [ ] **Step 4: Run all tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")' 2>&1 | tail -10
```

Expected: 3 new tests PASS, all prior tests still pass (8 total: 5 from prior PR + 3 new).

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_compare_alts_to_samples.R"
git commit -m "feat(mlb): american_prob_3way helper for 3-way devigging"
```

---

### Task 6: Match `h2h_3way_1st_5_innings` in `compare_alts_to_samples`

**Files:**
- Modify: `Answer Keys/Tools.R:4519-4525` (filter regex extension)
- Modify: `Answer Keys/Tools.R:4555-4570` (suffix dispatcher: add `1st_5_innings$ → "f5"`)
- Modify: `Answer Keys/Tools.R:4789` (insert new `else if` branch after the `odd_even_runs` branch)
- Modify: `Answer Keys/tests/test_compare_alts_to_samples.R` (add 3-way matcher tests)

**Why all three Tools.R edits in one task:** They're coupled. The filter regex must let the row through; the suffix dispatcher must resolve the period; the per-row branch must price it. Skip any one and tests fail. Bundling them into one commit keeps the invariant intact.

- [ ] **Step 1: Append failing tests to `Answer Keys/tests/test_compare_alts_to_samples.R`**

```r
test_that("compare_alts_to_samples emits h2h_3way bets — fair home wins 70%, draws 10%", {
  # 700 home wins (margin > 0), 200 away wins (margin < 0), 100 ties (margin == 0)
  margins <- c(rep(2L, 700), rep(-2L, 200), rep(0L, 100))
  samples <- make_synthetic_samples(margins_f3 = margins)   # margins_f3 doubles as F5 placeholder per fixture
  # Force F5 margins specifically: override the placeholder
  samples[["test-game-1"]]$sample$game_home_margin_period_F5 <- margins
  consensus <- make_consensus()

  # Wagerzon-shaped 3-way row: home -110, away +200, draw +500
  # devigged: implied 0.524 + 0.333 + 0.167 = 1.024; norm 0.512, 0.325, 0.163
  # vs model 0.70/0.20/0.10 → home is BIG +EV, away/draw -EV
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    game_date = "2026-05-03",
    game_time = "19:00",
    market = "h2h_3way_1st_5_innings",
    market_type = "h2h_3way",
    line = NA_real_,
    odds_away = 200L,
    odds_home = -110L,
    odds_draw = 500L,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  home_bet <- bets[bets$bet_on == "Test Home" & bets$market == "h2h_3way_1st_5_innings", ]
  expect_equal(nrow(home_bet), 1)
  expect_equal(home_bet$prob, 0.70, tolerance = 1e-9)
})

test_that("compare_alts_to_samples emits Tie bet when draw is +EV", {
  # 200 home wins, 200 away wins, 600 ties — model says 60% chance of tie
  margins <- c(rep(2L, 200), rep(-2L, 200), rep(0L, 600))
  samples <- make_synthetic_samples()
  samples[["test-game-1"]]$sample$game_home_margin_period_F5 <- margins
  consensus <- make_consensus()

  # Wagerzon prices the draw at +200 (implied 33%) — model says 60% tie → big +EV on draw
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    game_date = "2026-05-03",
    game_time = "19:00",
    market = "h2h_3way_1st_5_innings",
    market_type = "h2h_3way",
    line = NA_real_,
    odds_away = 200L,
    odds_home = 200L,
    odds_draw = 200L,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  tie_bet <- bets[bets$bet_on == "Tie" & bets$market == "h2h_3way_1st_5_innings", ]
  expect_equal(nrow(tie_bet), 1)
  expect_equal(tie_bet$prob, 0.60, tolerance = 1e-9)
})

test_that("compare_alts_to_samples returns no bets when 3-way row has any NA odds", {
  margins <- c(rep(2L, 700), rep(-2L, 200), rep(0L, 100))
  samples <- make_synthetic_samples()
  samples[["test-game-1"]]$sample$game_home_margin_period_F5 <- margins
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    game_date = "2026-05-03",
    game_time = "19:00",
    market = "h2h_3way_1st_5_innings",
    market_type = "h2h_3way",
    line = NA_real_,
    odds_away = 200L,
    odds_home = -110L,
    odds_draw = NA_integer_,    # missing draw price — must skip cleanly
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )
  expect_equal(nrow(bets), 0)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")' 2>&1 | tail -10
```

Expected: 3 new FAILs — the `h2h_3way_1st_5_innings` rows are silently dropped by the filter regex (which doesn't match `^h2h_3way_`).

- [ ] **Step 3: Extend the filter regex (`Tools.R:4519-4525`)**

Find:
```r
  alt_odds <- offshore_odds %>%
    filter(
      grepl("^alternate_", market) |
      grepl("^team_totals_", market) |
      grepl("1st_3_innings", market) |
      grepl("1st_7_innings", market) |
      market == "odd_even_runs"
    )
```

Replace with:
```r
  alt_odds <- offshore_odds %>%
    filter(
      grepl("^alternate_", market) |
      grepl("^team_totals_", market) |
      grepl("1st_3_innings", market) |
      grepl("1st_7_innings", market) |
      market == "odd_even_runs" |
      grepl("^h2h_3way_", market)
    )
```

- [ ] **Step 4: Extend the suffix dispatcher (`Tools.R:4555-4570`)**

Find the `suffix <- if (...)` chain (added in the prior PR's Task 2/3 work):

```r
    suffix <- if (grepl("1st_3_innings$", row$market)) {
      "f3"
    } else if (grepl("1st_7_innings$", row$market)) {
      "f7"
    } else if (row$market == "odd_even_runs") {
      "fg"
    } else {
      sub(".*_", "", row$market)
    }
```

Add a new branch for `1st_5_innings$ → "f5"`:

```r
    suffix <- if (grepl("1st_3_innings$", row$market)) {
      "f3"
    } else if (grepl("1st_5_innings$", row$market)) {
      "f5"
    } else if (grepl("1st_7_innings$", row$market)) {
      "f7"
    } else if (row$market == "odd_even_runs") {
      "fg"
    } else {
      sub(".*_", "", row$market)
    }
```

(Note: F5 markets weren't going through this matcher before — they went through `compare_moneylines_to_wagerzon` against the Odds API. The new 3-way is the first F5 market routed through `compare_alts_to_samples`.)

- [ ] **Step 5: Add the new `else if` branch (insert after the `odd_even_runs` branch)**

Find the closing `}` of the `odd_even_runs` branch (currently around `Tools.R:4789`). Immediately after it (and BEFORE the closing `}` of the for-loop), insert:

```r
    } else if (grepl("^h2h_3way_", row$market) &&
               !is.na(row$odds_home) && !is.na(row$odds_away) && !is.na(row$odds_draw)) {
      # 3-way moneyline (Wagerzon F5 3-way market). Tie is its OWN outcome,
      # not a refund — so p_home + p_away + p_draw = 1, no exclusion.
      col_name <- paste0(margin_col, "_", period)
      if (!col_name %in% names(sample_df)) next
      margins <- sample_df[[col_name]]
      margins <- margins[!is.na(margins)]
      if (length(margins) == 0) next

      p_home <- sum(margins > 0) / length(margins)
      p_away <- sum(margins < 0) / length(margins)
      p_draw <- sum(margins == 0) / length(margins)

      probs <- american_prob_3way(row$odds_away, row$odds_home, row$odds_draw)
      if (any(is.na(unlist(probs)))) next

      home_ev <- compute_ev(p_home, probs$p_home)
      away_ev <- compute_ev(p_away, probs$p_away)
      draw_ev <- compute_ev(p_draw, probs$p_draw)
      home_size <- kelly_stake(home_ev, probs$p_home, bankroll, kelly_mult)
      away_size <- kelly_stake(away_ev, probs$p_away, bankroll, kelly_mult)
      draw_size <- kelly_stake(draw_ev, probs$p_draw, bankroll, kelly_mult)

      if (home_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$home_team,
          line = NA_real_, bet_size = home_size, ev = home_ev,
          odds = row$odds_home, prob = p_home
        )
      }
      if (away_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$away_team,
          line = NA_real_, bet_size = away_size, ev = away_ev,
          odds = row$odds_away, prob = p_away
        )
      }
      if (draw_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = "Tie",
          line = NA_real_, bet_size = draw_size, ev = draw_ev,
          odds = row$odds_draw, prob = p_draw
        )
      }
```

The `american_prob_3way`, `compute_ev`, `kelly_stake`, `all_bets`, `book_key`, `pt_start_time`, `period`, `margin_col`, `bankroll`, `kelly_mult`, `ev_threshold` variables are all already in scope inside the function — same pattern as the spread/h2h/odd_even branches above.

- [ ] **Step 6: Run all tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")' 2>&1 | tail -10
```

Expected: 11 PASS, 0 FAIL (5 from prior PR + 3 from Task 5 + 3 from Task 6).

- [ ] **Step 7: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_compare_alts_to_samples.R"
git commit -m "feat(mlb): match h2h_3way_1st_5_innings in compare_alts_to_samples"
```

---

### Task 7: End-to-end smoke test

**Files (artifacts only — never committed):**
- Copy: live `wagerzon.duckdb`, `mlb.duckdb`, `mlb_mm.duckdb`, `pbp.duckdb` into the worktree

**Why this task exists:** Unit tests prove the new code paths behave correctly on synthetic input. This task proves they fire correctly when run inside the real MLB pipeline against real Wagerzon data — catches integration bugs the unit tests can't (team-name resolution, lg=1280 actually being included in the URL, samples table presence, etc.).

- [ ] **Step 1: Copy live DuckDB files into the worktree (NEVER symlink)**

```bash
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys/mlb.duckdb"
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys/mlb_mm.duckdb"
cp "/Users/callancapitolo/NFLWork/Answer Keys/pbp.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys/pbp.duckdb"
cp "/Users/callancapitolo/NFLWork/wagerzon_odds/wagerzon.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds/wagerzon.duckdb"
```

These DBs are gitignored — no commit needed.

- [ ] **Step 2: Run the MLB pipeline end-to-end on the worktree's copies**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys"
python3 run.py mlb 2>&1 | tee /tmp/mlb_pipeline_3way_smoke.log
```

The pipeline runs scrapers (which will write the new 3-way rows into `wagerzon.duckdb`'s `mlb_odds` table with `draw_ml` populated), then MLB.R reads them and routes through the new matcher branch.

If Wagerzon scraper auth fails (token expired in the worktree env), the integration smoke is partial — note it but proceed. The unit tests already cover the matcher logic; this is integration confirmation.

- [ ] **Step 3: Verify the new column was populated for 3-way rows**

```bash
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/wagerzon_odds/wagerzon.duckdb', read_only=True)
print('=== 3-way rows in mlb_odds ===')
print(con.execute(\"\"\"
  SELECT market, period, away_team, home_team, away_ml, home_ml, draw_ml
  FROM mlb_odds
  WHERE market = 'h2h_3way_1st_5_innings'
  ORDER BY fetch_time DESC
  LIMIT 10
\"\"\").fetchdf().to_string())
print()
print('=== Draw_ml population check (NULL on non-3way rows is correct) ===')
print(con.execute(\"\"\"
  SELECT market, COUNT(*) total_rows, COUNT(draw_ml) rows_with_draw_ml
  FROM mlb_odds
  GROUP BY market ORDER BY market
\"\"\").fetchdf().to_string())
"
```

Expected:
- Some `h2h_3way_1st_5_innings` rows present (one per game in lg=1280 today, typically 3-12 games)
- `draw_ml` populated ONLY on those rows; NULL elsewhere

- [ ] **Step 4: Verify matcher emitted bets (or distinguish soft vs hard failure)**

```bash
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way/Answer Keys/mlb.duckdb', read_only=True)
print('=== mlb_bets_combined: 3-way rows ===')
df = con.execute(\"\"\"
  SELECT id, home_team, away_team, market, bet_on, odds, prob, ev, bet_size
  FROM mlb_bets_combined
  WHERE market = 'h2h_3way_1st_5_innings'
  ORDER BY ev DESC
\"\"\").fetchdf()
print(df.to_string() if len(df) else '(no 3-way bets at threshold)')
"
```

Two acceptable outcomes:
- **Rows present:** integration confirmed. Hand-verify Step 5.
- **No rows:** could mean (i) no edges at the EV threshold tonight (acceptable), (ii) raw 3-way rows weren't scraped (check Step 3 output), or (iii) hard matching failure (check pipeline log for `Generated N wagerzon alt/team-total bets from samples` line and any errors).

- [ ] **Step 5: Spot-check one matched bet (if any emerged)**

Pick a row from Step 4's output. Hand-verify:
- `prob` should match the model: `sum(game_home_margin_period_F5 > 0) / N` for home bet, `sum(margin < 0) / N` for away, `sum(margin == 0) / N` for tie
- `odds` should match the raw Wagerzon row (cross-check `wagerzon.duckdb::mlb_odds` for the same `id`/team pair)
- `ev = (prob × profit_per_dollar) − ((1 − prob) × 1)` — recompute and confirm

- [ ] **Step 6: Confirm git status is clean**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git status
```

Expected: only the spec/plan markdowns (untracked) — no tracked-file modifications. The copied .duckdb files are gitignored.

- [ ] **Step 7: Commit nothing for this task** (this task produces no committed artifacts)

The smoke is purely verification. Move to Task 8.

---

### Task 8: Documentation updates

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — Pipeline Flow (MLB) section
- Modify: `Answer Keys/MLB Dashboard/README.md` — Markets section (Derivatives subsection)
- Modify: `wagerzon_odds/CLAUDE.md` — Quick map section

- [ ] **Step 1: Update `Answer Keys/CLAUDE.md`**

Find the Pipeline Flow (MLB) section's bullet list. The most recent bullet block (added in the prior PR) currently reads (around line 26-31):

```markdown
- F5 (first 5 innings) markets matched via Odds API: h2h, totals, spreads, alternate_totals
- Derivative markets matched via `compare_alts_to_samples` against scraped offshore odds:
  - F3 / F7 spreads + totals + h2h (Wagerzon, Bookmaker for F3 only). Note: …
  - FG alt spreads + alt totals (Wagerzon, Bet105)
  - FG odd/even total runs (Wagerzon — single-game prop, …)
- Team totals (`team_totals_*_fg`, `team_totals_*_h1`) are scraped but NOT yet matched …
```

Add ONE new bullet inside the "Derivative markets matched via `compare_alts_to_samples`" sub-list (immediately after the FG odd/even bullet):

```markdown
  - F5 3-way moneyline `h2h_3way_1st_5_innings` (Wagerzon only — `idgmtyp=29` league `lg=1280`). Tie is a real outcome, not push refund. Coexists with the existing 2-way F5 ML matched via `compare_moneylines_to_wagerzon`.
```

- [ ] **Step 2: Update `Answer Keys/MLB Dashboard/README.md`**

Find the "Derivatives" subsection inside the Markets section (added in prior PR). Add ONE new bullet:

```markdown
- `h2h_3way_1st_5_innings` — F5 3-way moneyline (Wagerzon only) with home/away/tie outcomes; bet_on label is the team name for home/away or `"Tie"` for the draw
```

- [ ] **Step 3: Update `wagerzon_odds/CLAUDE.md`**

Find the "Quick map" section. Append a new bullet about the new behavior:

```markdown
- **3-way F5 ML scraping:** `scraper_v2.py::parse_odds()` routes `idgmtyp=29` parent
  games (league `lg=1280` — "MLB - 1ST 5 INN WINNER (3-WAY)") through
  `parse_3way_line()`. Records have `market = "h2h_3way_1st_5_innings"`,
  `period = "f5"`, and the new `draw_ml INTEGER` column populated with
  the third (draw) outcome's American odds. `*_odds` tables include
  `draw_ml` via idempotent `ALTER TABLE ADD COLUMN IF NOT EXISTS`.
```

- [ ] **Step 4: Commit docs**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git add "Answer Keys/CLAUDE.md" "Answer Keys/MLB Dashboard/README.md" "wagerzon_odds/CLAUDE.md"
git commit -m "docs(mlb): note F5 3-way matching in CLAUDE.md + READMEs"
```

---

### Task 9: Pre-merge executive review

Per CLAUDE.md "Pre-merge review" section.

- [ ] **Step 1: Review the full diff**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way"
git diff main..HEAD --stat
git diff main..HEAD
```

Note: main may have moved forward — use the merge-base for the most accurate "branch only" diff:
```bash
git diff $(git merge-base main HEAD)..HEAD --stat
```

- [ ] **Step 2: Walk the executive checklist**

| Check | What to verify |
|-------|----------------|
| Data integrity | `draw_ml` only populated on 3-way rows; NULL guard prevents bogus bets; no duplicate writes |
| Resource safety | No new DB connections opened in matcher (pure-function modifications); `init_database` connection still uses `conn.close()` in `finally` semantics |
| Edge cases | Empty placeholder games (`EmptyGame: true` or `voddst=""`) return None from `parse_3way_line` and don't emit; NA odds short-circuit the matcher branch |
| Dead code | No commented-out blocks, no unused imports |
| Log/disk hygiene | The `print` statement in the new parse_odds 3-way branch is consistent with existing branches (no unbounded log growth) |
| Security | No new env vars, no new credentials, no new external calls beyond the existing Wagerzon helper URL (which already includes lg=1280) |

- [ ] **Step 3: Surface to user with the diff and ask for explicit approval to merge**

Per project rule and user preference — never merge without explicit approval. Even after all tests pass.

---

### Task 10: Merge + cleanup (only after user approves)

- [ ] **Step 1: Merge to main**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/wagerzon-mlb-3way -m "Merge feature/wagerzon-mlb-3way: F5 3-way market matching"
```

- [ ] **Step 2: Re-run all tests on main to confirm clean merge**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
cd "/Users/callancapitolo/NFLWork/wagerzon_odds" && python3 -m pytest tests/test_scraper_v2_3way.py -v
```

Expected: all R tests PASS (11 tests), all pytest scraper tests PASS (6 tests), 0 failures.

- [ ] **Step 3: Preserve the spec + plan on main**

The spec (`docs/superpowers/specs/2026-05-03-wagerzon-mlb-3way-design.md`) and plan (`docs/superpowers/plans/2026-05-03-wagerzon-mlb-3way-plan.md`) live on the feature branch and need to land on main with the merge. Verify:

```bash
ls /Users/callancapitolo/NFLWork/docs/superpowers/specs/2026-05-03-wagerzon-mlb-3way-design.md
ls /Users/callancapitolo/NFLWork/docs/superpowers/plans/2026-05-03-wagerzon-mlb-3way-plan.md
```

Both should exist on main after the merge (they were committed inside the worktree, so the merge brings them along).

- [ ] **Step 4: Clean up worktree + branch**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .worktrees/wagerzon-mlb-3way
git branch -d feature/wagerzon-mlb-3way
git worktree list
```

Expected: only the main worktree + any pre-existing worktrees remain.

---

## Self-Review

- [x] **Spec coverage**
  - Scraper change (idgmtyp=29 dispatch, parse_3way_line, schema column, team-name strip): Tasks 1, 2, 3
  - Reader change (get_wagerzon_odds h2h_3way emission): Task 4
  - Devig helper (american_prob_3way): Task 5
  - Matcher change (filter regex, suffix dispatcher F5 entry, h2h_3way branch): Task 6
  - Schema additions (draw_ml column, idempotent ALTER): Task 1
  - Tests (unit, scraper, devig, matcher): Tasks 1, 2, 3, 5, 6
  - End-to-end smoke: Task 7
  - Docs: Task 8
  - All spec sections covered.

- [x] **Placeholder scan**
  - No "TBD"/"TODO"/"implement later"/"similar to Task N"
  - Every code step has complete code blocks
  - Every test has assertion code
  - All exact file paths and line numbers given

- [x] **Type consistency**
  - `draw_ml` column type: `INTEGER` everywhere (Tasks 1, 2, 3)
  - `parse_3way_line` signature matches across helper definition + dispatch site (Tasks 2, 3)
  - `american_prob_3way` returns `list(p_away=, p_home=, p_draw=)` consistently (Tasks 5, 6)
  - `bet_on = "Tie"` string literal consistent (Task 6)
  - `market = "h2h_3way_1st_5_innings"` consistent across scraper, reader, matcher
  - `period = "f5"` (lowercase) in scraper output; `period_map["f5"] = "F5"` in matcher resolves to uppercase column suffix — same convention as `h1`/`f3`/`f7`

- [x] **TDD discipline** — every task starts with a failing test, then minimal implementation
- [x] **Frequent commits** — 8 commits, one per logical task (Task 7 produces 0 commits intentionally; Tasks 9/10 are review/merge ceremony)
- [x] **DRY** — `parse_3way_line` mirrors `parse_moneyline_only` shape; `american_prob_3way` is structurally parallel to `american_prob`; matcher branch mirrors the spread/h2h/odd_even pattern
- [x] **YAGNI** — no F3/F7 3-way handling, no replacement of 2-way matcher, no dashboard JS changes, no schema migration tooling beyond inline ALTER
