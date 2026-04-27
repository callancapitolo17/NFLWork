# DK Trifecta Blend Scaffolding (Plan #1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship the schema, NA-safe blend helper, and recon script that prepare `Answer Keys/mlb_triple_play.R` to consume DraftKings SGP fair odds as a blended supporting data point — without yet wiring DK live data through (that's Plan #2 after recon).

**Architecture:** Three independent additions: (1) a pure helper `blend_dk_with_model()` in `Answer Keys/parse_legs.R` with unit tests, (2) a Playwright/curl_cffi recon script `mlb_sgp/recon_dk_trifecta.py` ready to run when DK trifecta data is needed, (3) schema migration + R-side blend wiring + Plan #2 stub in `Answer Keys/mlb_triple_play.R`. The first two have no dependencies; the third must merge AFTER `feature/wagerzon-specials-scraper` lands on main since both modify `mlb_triple_play.R` in overlapping regions. Plan executes the first two tasks today, then halts at the Wagerzon roadblock.

**Tech Stack:** R (testthat for unit tests, base R for the helper), Python 3 (curl_cffi for recon — same pattern as `mlb_sgp/scraper_draftkings_sgp.py`). DuckDB for the new `mlb_trifecta_sgp_odds` table (created lazily on first pricer run via `CREATE TABLE IF NOT EXISTS`).

---

## File Structure

**Created:**
- `Answer Keys/tests/test_dk_blend.R` — testthat file. 5 unit tests for `blend_dk_with_model`. ~50 lines.
- `mlb_sgp/recon_dk_trifecta.py` — recon script. Hardcoded DK event IDs at top, fetches `parlays/v1/sgp/events/{id}`, dumps to `recon_dk_trifecta_<id>.json`. ~100 lines. Manual run only; not in any pipeline.

**Modified:**
- `Answer Keys/parse_legs.R` — append `blend_dk_with_model(model_prob, dk_decimal, vig)` function with roxygen docs. ~15 lines added.
- `Answer Keys/mlb_triple_play.R` — schema migration (`CREATE TABLE IF NOT EXISTS`), DK SGP read with `tryCatch`, blend wiring in the rowwise mutate, output table gains `model_odds` + `dk_odds` columns. **BLOCKED** until `feature/wagerzon-specials-scraper` merges to main (overlapping edits). ~40 lines added.
- `Answer Keys/CLAUDE.md` — extend "Triple-Play Data Flow" subsection with the DK blend layer. **BLOCKED** until Wagerzon's CLAUDE.md edits land (same subsection). ~12 lines added.

---

## Sequencing & roadblock

| # | Task | Status |
|---:|---|---|
| 1 | Blend helper + unit tests (`parse_legs.R` + new test file) | **READY TODAY** |
| 2 | Recon script (`mlb_sgp/recon_dk_trifecta.py`) | **READY TODAY** |
| 3 | Schema migration + R blend wiring (`mlb_triple_play.R`) | **BLOCKED** by `feature/wagerzon-specials-scraper` merge |
| 4 | Plan #2 stub comment (`mlb_triple_play.R`) | **BLOCKED** (same file) |
| 5 | Documentation (`Answer Keys/CLAUDE.md`) | **BLOCKED** (same subsection as Wagerzon) |

When the Wagerzon plan lands on main, resume from Task 3.

---

## Worktree & Version Control Plan

- **Branch:** `feature/dk-trifecta-blend-scaffolding`
- **Worktree:** `.worktrees/dk-trifecta-scaffolding`
- **Setup (before Task 1):**
  ```bash
  cd /Users/callancapitolo/NFLWork
  git worktree add .worktrees/dk-trifecta-scaffolding -b feature/dk-trifecta-blend-scaffolding main
  cd .worktrees/dk-trifecta-scaffolding
  git branch  # confirm feature/dk-trifecta-blend-scaffolding
  ```
- **Commits** (one per task as it completes):
  1. Task 1 → `feat(mlb): add blend_dk_with_model helper for trifecta DK SGP integration`
  2. Task 2 → `feat(dk): add recon_dk_trifecta.py for DK SGP selection-id reconnaissance`
  3. Task 3 → `feat(mlb): wire DK trifecta blend into mlb_triple_play.R`
  4. Task 4 → `chore(mlb): document Plan #2 stub for DK scraper invocation`
  5. Task 5 → `docs(mlb): document DK trifecta blend in Triple-Play Data Flow`
- **DuckDB:** none of these tasks need a DuckDB copy. Task 3's `CREATE TABLE IF NOT EXISTS` runs against any DB on first pricer invocation.
- **Cleanup after merge:** `git worktree remove .worktrees/dk-trifecta-scaffolding && git branch -d feature/dk-trifecta-blend-scaffolding`
- Never merge to main without explicit user approval.

---

## Documentation Plan

Task 5 extends the existing "Triple-Play Data Flow" subsection in `Answer Keys/CLAUDE.md` (introduced 2026-04-21, modified 2026-04-22 by the Wagerzon plan) with the DK blend layer. The exact insertion point depends on what the Wagerzon plan's Task 5 wrote — Task 5 of THIS plan reads the post-Wagerzon CLAUDE.md and inserts the DK section after the `parse_legs.R` block.

No `mlb_sgp/README.md` update in this plan — the recon script is an artifact, not a permanent component. When Plan #2 ships the production scraper, that plan adds the README entry alongside it.

---

## Pre-Merge Review Checklist

Run `git diff main..HEAD` in the worktree before merging:

- **Behavior preservation:** running `mlb_triple_play.R` with no rows in `mlb_trifecta_sgp_odds` (the expected state at first run) MUST produce identical fair odds to today's committed values. The blend reduces to model-only when DK is unavailable. Verify by running the pricer in the worktree and comparing fair_odds row-by-row to the most recent main run.
- **Schema isolation:** `CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds` only adds the table, never modifies existing tables.
- **Resource safety:** new `dbConnect` for the trifecta read uses `tryCatch` for missing-table fallback; existing `on.exit(dbDisconnect(con))` covers the connection.
- **Dead code:** `blend_dk_with_model` is referenced by the rowwise mutate; no orphan function.
- **Type discipline:** `model_prob` is numeric; `dk_decimal` is numeric or NA; return is numeric or NA_real_. All paths covered by tests.
- **Pricer NA handling:** when `dk_match` returns 0 rows, `dk_decimal` becomes NA_real_, `blend_dk_with_model(model, NA, vig)` returns model_prob — verified by Task 1 tests.
- **No DuckDB writes anywhere new:** the schema migration is a `CREATE TABLE IF NOT EXISTS`, not an INSERT. The pricer remains read-only against `mlb_trifecta_sgp_odds`.
- **No secrets in logs:** recon script logs URLs and HTTP status only.
- **Plan #2 stub is unambiguous:** the commented-out `system2()` call has a `# Plan #2 will activate the scraper invocation here` marker so future-Plan-#2 has a clear single-edit-point.

---

## Task 1: Blend helper + unit tests

**Status: READY TODAY** — no dependency on Wagerzon plan.

**Files:**
- Modify: `Answer Keys/parse_legs.R` (append the helper function + roxygen docs)
- Create: `Answer Keys/tests/test_dk_blend.R`

### - [ ] Step 1: Write failing tests

Create `Answer Keys/tests/test_dk_blend.R`:

```r
# Answer Keys/tests/test_dk_blend.R
library(testthat)
source("../parse_legs.R")

test_that("blend returns model_prob when DK decimal is NA", {
  expect_equal(blend_dk_with_model(0.30, NA_real_, 1.10), 0.30)
})

test_that("blend averages model and DK devigged when both present", {
  # DK SGP +200 → decimal 3.0 → implied 0.3333... Devigged at vig=1.10 → 0.30303
  # Model 0.30. Expected blend = mean(0.30, 0.30303) ≈ 0.30152
  expect_equal(blend_dk_with_model(0.30, 3.0, 1.10),
               (0.30 + (1/3.0) / 1.10) / 2,
               tolerance = 1e-9)
})

test_that("blend treats DK decimal of 0 or negative as missing", {
  expect_equal(blend_dk_with_model(0.30, 0,    1.10), 0.30)
  expect_equal(blend_dk_with_model(0.30, -1.0, 1.10), 0.30)
})

test_that("blend returns DK alone when model is NA", {
  # Edge case: if compute_prop_fair returned NA (e.g. zero valid samples),
  # DK is the only available signal.
  expect_equal(blend_dk_with_model(NA_real_, 3.0, 1.10),
               (1/3.0) / 1.10,
               tolerance = 1e-9)
})

test_that("blend is NA when both inputs are NA / missing", {
  expect_true(is.na(blend_dk_with_model(NA_real_, NA_real_, 1.10)))
  expect_true(is.na(blend_dk_with_model(NA_real_, 0,        1.10)))
})
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_dk_blend.R")'
```
Expected: 5 test_that blocks fail with `could not find function "blend_dk_with_model"`.

### - [ ] Step 3: Implement the helper

Append to `Answer Keys/parse_legs.R` (at the end of the file, after `compute_prop_fair`):

```r

#' Blend a model fair probability with a DK SGP devigged probability.
#'
#' If `dk_decimal` is NA / 0 / negative, treats DK as missing and returns
#' `model_prob` alone.
#' If `model_prob` is NA but DK is valid, returns DK's devigged probability
#' alone.
#' If both are missing, returns NA_real_.
#' Otherwise returns the mean of the two devigged probabilities — the
#' blending strategy used by mlb_correlated_parlay.R for DK + FD SGPs.
#'
#' @param model_prob numeric, the historical fair probability from
#'   compute_prop_fair (may be NA for unpriceable props).
#' @param dk_decimal numeric, DK SGP decimal odds (>1.0 if valid; NA / 0
#'   / negative when DK couldn't price the prop).
#' @param vig numeric, the SGP vig fallback (default convention is 1.10,
#'   matching DK_SGP_VIG_DEFAULT in mlb_correlated_parlay.R).
#' @return numeric scalar — the blended probability — or NA_real_ if no
#'   inputs are available.
blend_dk_with_model <- function(model_prob, dk_decimal, vig) {
  dk_valid <- !is.na(dk_decimal) && dk_decimal > 0
  dk_prob  <- if (dk_valid) (1 / dk_decimal) / vig else NA_real_
  probs    <- c(model_prob, dk_prob)
  probs    <- probs[!is.na(probs)]
  if (length(probs) == 0) return(NA_real_)
  mean(probs)
}
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_dk_blend.R")'
```
Expected: 5 PASS, 0 FAIL.

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
git add "Answer Keys/parse_legs.R" "Answer Keys/tests/test_dk_blend.R"
git commit -m "feat(mlb): add blend_dk_with_model helper for trifecta DK SGP integration"
```

---

## Task 2: Recon script

**Status: READY TODAY** — no dependency on Wagerzon plan.

**Files:**
- Create: `mlb_sgp/recon_dk_trifecta.py`

### - [ ] Step 1: Create the recon script

Create `mlb_sgp/recon_dk_trifecta.py`:

```python
"""
DraftKings Trifecta Recon Script

Captures DK's `parlays/v1/sgp/events/{event_id}` response for one or two
known MLB games, plus an attempted `calculateBets` POST with candidate
selection IDs. Output guides the Plan #2 leg-resolver implementation:
- "First Team to Score" market name + selection ID format
- "1st 5 Innings Winner" market — 2-way ML or 3-way?
- "Team Total Over/Under" market for the GRAND-SLAM 4th leg
- Whether DK accepts these legs as a single SGP

Usage:
    1. Find one or two DK MLB game URLs (e.g. https://sportsbook.draftkings.com/event/<event_id>)
    2. Update DK_EVENT_IDS below
    3. Run: python3 mlb_sgp/recon_dk_trifecta.py
    4. Inspect the resulting recon_dk_trifecta_<event_id>.json files
    5. Plan #2 reads these to populate the leg resolvers

Manual one-shot artifact. Not part of any pipeline.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Optional

try:
    from curl_cffi import requests
except ImportError:
    print("ERROR: curl_cffi not installed. The DK SGP scraper venv has it.", file=sys.stderr)
    print("       cd mlb_sgp && python -m venv venv && venv/bin/pip install curl_cffi", file=sys.stderr)
    sys.exit(1)

# ===== CONFIG — UPDATE THESE BEFORE RUNNING =====
# DK event IDs for ~2 MLB games we want to characterize. Find these by:
# 1. Open https://sportsbook.draftkings.com/leagues/baseball/mlb in a browser
# 2. Click into a game, look at the URL: /event/<event_id>
# 3. Pick games that have triple-play / grand-slam comparable markets
#    (any MLB game should work — we only need the selection IDs)
DK_EVENT_IDS: list[int] = [
    # 32109876,  # example — replace with real IDs before running
    # 32109877,
]

# DK API endpoints (public, no auth required)
SGP_EVENTS_URL_TPL = "https://sportsbook-nash.draftkings.com/parlays/v1/sgp/events/{event_id}"

OUT_DIR = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("recon_dk_trifecta")


def fetch_event_selections(event_id: int) -> Optional[dict]:
    """Fetch the full SGP event payload (~2MB JSON of all selection IDs)."""
    url = SGP_EVENTS_URL_TPL.format(event_id=event_id)
    log.info("GET %s", url)
    try:
        r = requests.get(url, impersonate="chrome", timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log.exception("Failed to fetch event %d: %s", event_id, e)
        return None


def summarize_markets(payload: dict) -> dict:
    """Walk the payload and pull out market-name / selection-id summaries
    relevant to triple-play / grand-slam legs:
      - "First Team to Score" or similar
      - "1st 5 Innings Winner"
      - "Game Moneyline" (h2h)
      - Team Total Over/Under
    Returns a dict suitable for Plan #2's resolver design."""
    summary = {
        "first_to_score": [],
        "first_5_innings_winner": [],
        "moneyline": [],
        "team_totals": [],
        "other_markets_found": [],
    }
    # The DK SGP payload structure varies; this is a permissive walker
    # that records anything that looks like a market.
    def walk(node, path=""):
        if isinstance(node, dict):
            name = (node.get("marketName") or node.get("name") or "").strip()
            if name:
                lname = name.lower()
                bucket = None
                if "first" in lname and "score" in lname:
                    bucket = "first_to_score"
                elif "1st 5" in lname or "first 5" in lname:
                    bucket = "first_5_innings_winner"
                elif "moneyline" in lname or lname in ("h2h", "game"):
                    bucket = "moneyline"
                elif "team total" in lname:
                    bucket = "team_totals"
                else:
                    summary["other_markets_found"].append({"path": path, "name": name})
                if bucket:
                    selections = node.get("selections") or node.get("outcomes") or []
                    summary[bucket].append({
                        "path": path,
                        "name": name,
                        "selections": [
                            {"name": s.get("label") or s.get("name"),
                             "id": s.get("id") or s.get("selectionId")}
                            for s in selections
                            if isinstance(s, dict)
                        ],
                    })
            for k, v in node.items():
                walk(v, f"{path}.{k}" if path else k)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")
    walk(payload)
    return summary


def main() -> int:
    if not DK_EVENT_IDS:
        log.error("DK_EVENT_IDS is empty. Edit %s and add 1-2 event ids first.", __file__)
        return 1
    for event_id in DK_EVENT_IDS:
        log.info("=== Reconning event %d ===", event_id)
        payload = fetch_event_selections(event_id)
        if payload is None:
            log.warning("Skipping %d (fetch failed)", event_id)
            continue
        out_full = OUT_DIR / f"recon_dk_trifecta_{event_id}.json"
        out_full.write_text(json.dumps(payload, indent=2, default=str))
        log.info("Saved full payload to %s", out_full.name)

        summary = summarize_markets(payload)
        out_summary = OUT_DIR / f"recon_dk_trifecta_{event_id}_summary.json"
        out_summary.write_text(json.dumps(summary, indent=2, default=str))
        n_first = len(summary["first_to_score"])
        n_5inn  = len(summary["first_5_innings_winner"])
        n_ml    = len(summary["moneyline"])
        n_tt    = len(summary["team_totals"])
        log.info("Summary: first_to_score=%d  first_5_innings=%d  moneyline=%d  team_totals=%d",
                 n_first, n_5inn, n_ml, n_tt)
        log.info("Wrote market summary to %s", out_summary.name)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

### - [ ] Step 2: Syntax check

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
python3 -c "import ast; ast.parse(open('mlb_sgp/recon_dk_trifecta.py').read()); print('SYNTAX OK')"
```
Expected: `SYNTAX OK`.

### - [ ] Step 3: Import smoke test (no execution)

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
python3 -c "
import sys
sys.path.insert(0, 'mlb_sgp')
# Avoid import-time curl_cffi dep for this smoke test
import ast
src = open('mlb_sgp/recon_dk_trifecta.py').read()
tree = ast.parse(src)
fn_names = [n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]
print('Functions:', fn_names)
assert 'fetch_event_selections' in fn_names
assert 'summarize_markets' in fn_names
assert 'main' in fn_names
print('IMPORT-LEVEL CHECK OK')
"
```
Expected: prints `Functions: [...]` and `IMPORT-LEVEL CHECK OK`. We avoid actually importing the module here so it's runnable on machines without `curl_cffi`.

### - [ ] Step 4: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
git add mlb_sgp/recon_dk_trifecta.py
git commit -m "feat(dk): add recon_dk_trifecta.py for DK SGP selection-id reconnaissance"
```

---

## ROADBLOCK — pause here

**Tasks 3, 4, and 5 modify `Answer Keys/mlb_triple_play.R` and `Answer Keys/CLAUDE.md`. Both files are also modified by `feature/wagerzon-specials-scraper` (plan committed `a80ab5c`).** Trying to land Tasks 3-5 before that branch merges would create unnecessary conflicts.

**Resume conditions** (all must be true):
1. `feature/wagerzon-specials-scraper` is merged to `main`
2. The post-merge `Answer Keys/mlb_triple_play.R` reads `wagerzon_specials` and produces the 14-row triple-play + grand-slam table
3. `Answer Keys/CLAUDE.md` has the `wagerzon_odds/scraper_specials.py` block in its Triple-Play Data Flow subsection

When resuming: rebase this branch onto the new `main`, then continue from Task 3.

---

## Task 3: Schema migration + R-side blend wiring

**Status: BLOCKED until Wagerzon plan merges to main and this branch is rebased.**

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R`

### - [ ] Step 1: Add the schema migration near the top of the main block

In `Answer Keys/mlb_triple_play.R`, inside the `if (!interactive() && sys.nframe() == 0L)` guard, immediately after the `setwd("~/NFLWork/Answer Keys")` line and the `source("parse_legs.R")` line, add a new schema-migration block:

```r
  # Ensure DK trifecta SGP table exists (idempotent; created lazily on first
  # pricer run). The Plan #2 scraper writes here; Plan #1 reads it (empty).
  MLB_DB <- "mlb.duckdb"
  con_mig <- dbConnect(duckdb(), dbdir = MLB_DB)
  dbExecute(con_mig, "
    CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds (
      fetch_time     TIMESTAMP,
      game_id        VARCHAR,
      prop_type      VARCHAR,
      side           VARCHAR,
      legs_json      VARCHAR,
      selection_ids  VARCHAR,
      sgp_decimal    DOUBLE,
      sgp_american   INTEGER,
      source         VARCHAR
    )
  ")
  dbDisconnect(con_mig)
```

(`MLB_DB` may already be defined elsewhere in the post-Wagerzon main block. If so, reuse the existing assignment and don't duplicate it.)

### - [ ] Step 2: Add the DK SGP read after the samples query

In the same main block, after the post-Wagerzon `samples_df <- dbGetQuery(con, ...)` and `consensus <- dbGetQuery(con, ...)` calls, add the DK SGP read:

```r
  # Read latest DK trifecta SGP odds (table empty until Plan #2 populates it)
  dk_sgp <- tryCatch({
    dbGetQuery(con, "
      SELECT game_id, prop_type, side, sgp_decimal
      FROM mlb_trifecta_sgp_odds
      WHERE source = 'draftkings_direct'
        AND fetch_time = (SELECT MAX(fetch_time) FROM mlb_trifecta_sgp_odds)
    ")
  }, error = function(e) data.frame(
    game_id     = character(0),
    prop_type   = character(0),
    side        = character(0),
    sgp_decimal = double(0)
  ))
```

### - [ ] Step 3: Insert blend logic into the rowwise mutate

Find the existing `priced <- matched %>% rowwise() %>% mutate(...)` block. Replace the existing `mutate(...)` body with the blended version below. The four key changes vs. the post-Wagerzon-merge version: `dk_match` lookup, `dk_decimal` extraction, `blend_dk_with_model()` call replacing the direct `compute_prop_fair` assignment to `fair_prob`, and three new output columns:

```r
  DK_SGP_VIG_DEFAULT <- 1.10  # matches mlb_correlated_parlay.R convention

  priced <- matched %>%
    rowwise() %>%
    mutate(
      game_samples    = list(samples_df[samples_df$game_id == id, ]),
      n_samples       = nrow(game_samples),
      legs            = list(parse_legs(description)),
      model_fair_prob = compute_prop_fair(game_samples, side, legs),

      dk_match        = list(dk_sgp[dk_sgp$game_id   == id        &
                                    dk_sgp$prop_type == prop_type &
                                    dk_sgp$side      == side, ]),
      dk_decimal      = if (nrow(dk_match) > 0) dk_match$sgp_decimal[1] else NA_real_,
      dk_fair_prob    = if (!is.na(dk_decimal) && dk_decimal > 0) {
                          (1 / dk_decimal) / DK_SGP_VIG_DEFAULT
                        } else NA_real_,

      fair_prob       = blend_dk_with_model(model_fair_prob, dk_decimal,
                                            DK_SGP_VIG_DEFAULT),

      model_odds      = prob_to_american(model_fair_prob),
      dk_odds         = prob_to_american(dk_fair_prob),
      fair_odds       = prob_to_american(fair_prob),
      book_prob       = american_to_prob(book_odds),
      edge_pct        = (fair_prob / book_prob - 1) * 100,

      prop_type       = {
        known_props <- "(TRIPLE-PLAY|GRAND-SLAM)"
        m <- regmatches(description,
                        regexec(paste0("\\b", known_props, "\\b"), description))[[1]]
        if (length(m) >= 2) m[[2]] else NA_character_
      }
    ) %>%
    ungroup() %>%
    select(target_team, prop_type, side, n_samples,
           model_odds, dk_odds, fair_odds, book_odds, edge_pct) %>%
    arrange(desc(edge_pct))
```

### - [ ] Step 4: Verify regression — fair odds unchanged when DK table is empty

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cd "Answer Keys" && Rscript mlb_triple_play.R 2>&1 | tail -25
```

Expected: 14-row table prints with `fair_odds` populated. `model_odds` column equals `fair_odds` for every row (no DK data → blend reduces to model only). `dk_odds` is blank/NA for every row. The fair_odds values must match the most recent committed run from main (within ±1 American odds for rounding noise) for every row.

If any row's fair_odds differs by more than ±1 from the prior committed run, STOP. The blend logic must be a strict superset of the existing `compute_prop_fair` output when DK is absent.

### - [ ] Step 5: Remove the copied DuckDB

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
rm "Answer Keys/mlb.duckdb"
```

### - [ ] Step 6: Run R unit tests to confirm no regression

```bash
cd "Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_dk_blend.R"); test_file("test_parse_legs.R"); test_file("test_triple_play.R")'
```
Expected: all pass (Task 1's 5 + parse_legs's 32 + triple_play's existing tests).

### - [ ] Step 7: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
git add "Answer Keys/mlb_triple_play.R"
git commit -m "feat(mlb): wire DK trifecta blend into mlb_triple_play.R"
```

---

## Task 4: Plan #2 stub comment

**Status: BLOCKED until Task 3 lands.**

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R`

### - [ ] Step 1: Add Plan #2 stub marker

In `Answer Keys/mlb_triple_play.R`, immediately AFTER the schema-migration block (added in Task 3 Step 1) and BEFORE the existing samples-and-consensus reads, insert a clearly-marked stub block:

```r
  # ===== Plan #2 will activate the DK trifecta SGP scraper here =====
  # When Plan #2 lands, replace the comment block below with a writable
  # request file + system2() call to mlb_sgp/scraper_draftkings_trifecta.py.
  # The blend logic below already handles the populated table correctly;
  # no changes to the rowwise mutate are needed when activating.
  #
  # trifecta_input <- todays_lines %>% ... write JSON ...
  # system2(SGP_VENV_PYTHON,
  #         args = c(file.path(SGP_DIR, "scraper_draftkings_trifecta.py"),
  #                  "--input", request_path, "--db", MLB_DB))
  # ===================================================================
```

### - [ ] Step 2: Verify nothing else changed

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("SYNTAX OK\n")'
```
Expected: `SYNTAX OK`.

### - [ ] Step 3: Run pricer end-to-end again to confirm comment-only change

```bash
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cd "Answer Keys" && Rscript mlb_triple_play.R 2>&1 | tail -25
rm mlb.duckdb
```
Expected: same 14-row table as Task 3 Step 4, with identical fair_odds.

### - [ ] Step 4: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
git add "Answer Keys/mlb_triple_play.R"
git commit -m "chore(mlb): document Plan #2 stub for DK scraper invocation"
```

---

## Task 5: Documentation — Triple-Play Data Flow gains DK blend

**Status: BLOCKED until Task 4 lands AND the post-Wagerzon CLAUDE.md is on main.**

**Files:**
- Modify: `Answer Keys/CLAUDE.md`

### - [ ] Step 1: Insert the DK blend layer into the Triple-Play Data Flow subsection

After the Wagerzon plan lands, the subsection has a `wagerzon_odds/scraper_specials.py` block at the top, followed by `parse_legs.R`, `mlb_triple_play.R`, and bullets. Find the `mlb_triple_play.R (standalone pricer)` block. Replace it with the version below (adds DK blend description), and append the new bullets:

Find the existing `mlb_triple_play.R (standalone pricer)` block (post-Wagerzon version):

```
mlb_triple_play.R (standalone pricer)
  ├── Reads wagerzon_specials (latest scraped_at) for posted lines
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + scored_first)
  ├── Maps Wagerzon team names → Odds API canonical via WZ_TO_CANONICAL dict
  ├── Joins to consensus_temp for game_id + side (home/away)
  ├── For each row: parse_legs(description) → compute_prop_fair(...)
  └── Prints fair odds vs book + edge per posted line, grouped by prop_type
```

Replace with:

```
mlb_triple_play.R (standalone pricer)
  ├── Reads wagerzon_specials (latest scraped_at) for posted lines
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + scored_first)
  ├── Reads mlb_trifecta_sgp_odds (latest fetch_time) for DK SGP fair odds
  ├── Maps Wagerzon team names → Odds API canonical via WZ_TO_CANONICAL dict
  ├── Joins to consensus_temp for game_id + side (home/away)
  ├── For each row:
  │     parse_legs(description) → compute_prop_fair(...) → model_fair_prob
  │     dk_sgp lookup → devig with DK_SGP_VIG_DEFAULT=1.10 → dk_fair_prob
  │     blend_dk_with_model(model, dk, vig) → fair_prob (the published fair)
  └── Prints model_odds, dk_odds, fair_odds (blended), book_odds, edge_pct
```

Then add three bullets to the existing bullet list, after the `WZ_TO_CANONICAL` bullet:

```markdown
- DK trifecta SGP blend mirrors mlb_correlated_parlay.R: per-prop devig with `DK_SGP_VIG_DEFAULT = 1.10` (only 2 obs/game so per-game vig fitting isn't possible).
- Blend = mean of available probs (model + DK). When DK is unavailable (table empty, missing leg, scraper failure), blend reduces to model-only.
- Plan #2 (post-recon) populates `mlb_trifecta_sgp_odds` via `mlb_sgp/scraper_draftkings_trifecta.py`. Plan #1 ships the schema + blend scaffolding; the table stays empty until Plan #2 lands.
```

### - [ ] Step 2: Verify the file reads correctly

```bash
grep -A 60 "Triple-Play Data Flow" "Answer Keys/CLAUDE.md" | head -80
```
Expected: the updated section with all four data-flow layers (Wagerzon → parse_legs → mlb_triple_play → DK blend) visible, plus the three new bullets.

### - [ ] Step 3: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scaffolding
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb): document DK trifecta blend in Triple-Play Data Flow"
```

---

## Self-Review

**Spec coverage** (from `docs/superpowers/specs/2026-04-22-dk-trifecta-blend-design.md`):

| Spec section | Covered by |
|---|---|
| Schema migration (`mlb_trifecta_sgp_odds`) | Task 3 Step 1 (CREATE TABLE IF NOT EXISTS) |
| `blend_dk_with_model` helper | Task 1 |
| Recon script (`mlb_sgp/recon_dk_trifecta.py`) | Task 2 |
| R-side blend scaffolding (NA-safe) | Task 3 Steps 2-3 |
| Plan #2 stub comment | Task 4 |
| Output columns: `model_odds`, `dk_odds`, `fair_odds` | Task 3 Step 3 |
| Documentation update | Task 5 |
| Unit tests | Task 1 (5 tests on blend_dk_with_model) |
| Sequencing relative to Wagerzon plan | Roadblock between Tasks 2 and 3 |

**Placeholder scan:** No "TBD" / "TODO" / "implement later". Each step has either complete code or a concrete shell command with expected output.

**Type consistency:**
- `blend_dk_with_model(model_prob, dk_decimal, vig)` returns numeric or NA_real_; consumers (`prob_to_american`) handle both.
- `mlb_trifecta_sgp_odds.game_id` is VARCHAR; matches `mlb_consensus_temp.id` (the join key on the R side).
- `dk_decimal` is numeric or NA_real_ from start to finish through the rowwise mutate.
- `DK_SGP_VIG_DEFAULT = 1.10` matches `mlb_correlated_parlay.R` convention.
- Column names `model_odds`, `dk_odds`, `fair_odds` consistent across blend wiring (Task 3) and documentation (Task 5).

**Behavior preservation:** Task 3 Step 4's regression check confirms blended fair odds match model-only fair odds when DK table is empty. The blend reduces to `mean(c(model_prob))` = `model_prob`.

No issues found.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-22-dk-trifecta-blend-scaffolding.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session with checkpoints.

**Today's executable tasks: 1 and 2.** Task 3 onward halts at the roadblock until `feature/wagerzon-specials-scraper` merges.

**Which approach?**
